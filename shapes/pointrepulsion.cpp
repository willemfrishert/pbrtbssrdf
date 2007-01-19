#include "pointrepulsion.h"
#include "core/sampling.h"
#include "transform.h"
#include "TriangleMesh.h"

#ifndef round
#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))
#endif

PointRepulsion::PointRepulsion( int aNumberOfTriangles, int aNumberOfVertices, int* aVertexIndex, 
							   Point* aVertices, vector<Reference<Shape> >& aTriangleList  ) 
: iNumberOfTriangles( aNumberOfTriangles )
, iNumberOfVertices( aNumberOfVertices )
, iVertexIndex( aVertexIndex )
, iVertices( aVertices)
, iTotalAreaOfSurface( 0.0f )
{	
	// setup the triangle structure (find edge neighbors)
	SetupTriangleUseSets( aTriangleList );

	// compute the total area of the object
	vector<TriangleUseSet>::iterator triangleIter = this->iTriangles.begin();

	while ( triangleIter != this->iTriangles.end() ) 
	{
		iTotalAreaOfSurface += (*triangleIter).GetArea();
		printf("accumilated area = %f\n", iTotalAreaOfSurface);
		triangleIter++;
	}

	ComputePartialAreaSum( iPartialAreaSums );

	// add stratified sample points to the triangle mesh using a
	// list of partial sums of the polygon areas
	//SetupSamplePoints();

	////// compute the repulsive range r
	//iRepulsiveRadius = 2*sqrt( (iTotalAreaOfSurface/iNumberOfSamplePoints) );

	//// relaxation of the random points 
	//for(int k = 0; k<aNumberOfIterations; k++)
	//{
	//	ComputeRepulsiveForces();
	//	ComputeNewPositions();
	//}

	// testing binary search
	//iNumberOfTriangles = 4;
	//vector<float> areas;
	//areas.push_back( 0.1 );
	//areas.push_back( 0.3 );
	//areas.push_back( 0.6 );
	//areas.push_back( 1.0 );
	//int i = SearchForTriangle(0.05, areas, iNumberOfTriangles/2 );
	//i = SearchForTriangle(0.2, areas, iNumberOfTriangles/2 );
	//i = SearchForTriangle(0.5, areas, iNumberOfTriangles/2 );
	//i = SearchForTriangle(0.7, areas, iNumberOfTriangles/2 );
	//i = SearchForTriangle(1.0, areas, iNumberOfTriangles/2 );
}

PointRepulsion::~PointRepulsion()
{
	for ( int i=0, maxSamplePoints = iSamplePoints.size(); i<maxSamplePoints; i++)
	{
		delete iSamplePoints.at(i);
		iSamplePoints.at(i) = NULL;
	}
	iSamplePoints.clear();

	for ( int i=0, maxSamplePointContainers = iSamplePointContainer.size(); i<maxSamplePointContainers; i++)
	{
		delete iSamplePointContainer.at(i);
		iSamplePointContainer.at(i) = NULL;
	}
	iSamplePointContainer.clear();

	iTriangles.clear();
};

/*!
 * \brief
 * Sets up the triangle usesets. Store the neighboring triangles and the rotation 
 * matrices needed to make the triangles coplanar
 * 
 */
void PointRepulsion::SetupTriangleUseSets(vector<Reference<Shape> >& aTriangleList)
{
	int* vertices1 = NULL;

	//set up all the triangle use sets
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		Triangle* triangle = PointRepulsion::Cast(aTriangleList.at(i));
		triangle->GetVertexIndices( &vertices1 );

		TriangleUseSet triangleUseSet( aTriangleList.at(i), (iVertices+vertices1[0]), (iVertices+vertices1[1]), (iVertices+vertices1[2]) );
		this->iTriangles.push_back( triangleUseSet );
	}

	// loop through all the indices to find triangle neighbors
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		vector<Neighbor*> neighbor;
		this->iTriangles.at(i).GetEdgeNeighbors( neighbor );
		if ( neighbor.size() > 2)
		{
			continue;
		}

		vertices1 = &iVertexIndex[3*i];

		// go through the triangles i+1 (triangles that are 
		// before 1 have already registered as being neighbors).
		for (int j=i+1; j< iNumberOfTriangles; j++)
		{
			Point* vertexNeighbor[2] = { NULL };

			int vertexNeighborCount = 0;

			int* vertices2 =&iVertexIndex[3*j];

			// these indices contain the point numbers counting from 0 to 2
			// and are used to calculate the edge numbers.
			// the arrays needs to be kept in order of 0 to 2
			int indices[2] = {0,0};
			int indicesNeighbor[2] = {0,0};

			// contain edge number (either 0, 1 or 2)
			u_int edgePosition = 0;
			u_int edgePositionNeighbor = 0;

			// loop through the two triangles and try to match vertex indices.
			for (int k=0; k<3;k++)
			{
				for(int l=0; l<3;l++)
				{
					if (vertices1[k] == vertices2[l])
					{
						vertexNeighbor[vertexNeighborCount] = (iVertices+vertices1[k]);
						indices[vertexNeighborCount] = k;
						indicesNeighbor[vertexNeighborCount] = l;

						vertexNeighborCount++;
						continue;
					}
				}

				// once we found two neighbors, we have to find out to which edge they are connected and jump out of this loop
				if ( (NULL != vertexNeighbor[0]) && (NULL != vertexNeighbor[1]) )
				{
					// calculate the edge position for this triangle
					if ((indices[0] == 0) && (indices[1] == 1))
					{
						edgePosition = 0;
					}
					else if ((indices[0] == 1) && (indices[1] == 2))
					{
						edgePosition = 1;
					}
					else 
					{
						edgePosition = 2;
					}

					// make sure the numbers are in order
					// (for the neighbor triangle edge, this might not be the case)
					int minVal = min(indicesNeighbor[0], indicesNeighbor[1]);
					int maxVal = max(indicesNeighbor[0], indicesNeighbor[1]);

					indicesNeighbor[0] = minVal;
					indicesNeighbor[1] = maxVal;
					if ((indicesNeighbor[0] == 0) && (indicesNeighbor[1] == 1))
					{
						edgePositionNeighbor = 0;
					}
					else if ((indicesNeighbor[0] == 1) && (indicesNeighbor[1] == 2))
					{
						edgePositionNeighbor = 1;
					}
					else 
					{
						edgePositionNeighbor = 2;
					}
					break;
				}
			}


			// if neighbor is found, calculate rotation rotation angle, matrices and register neighbors
			if ( (NULL != vertexNeighbor[0]) && (NULL != vertexNeighbor[1]) )
			{
				Reference<Matrix4x4> alignToAxis = new Matrix4x4;
				Reference<Matrix4x4> alignToAxisInv = new Matrix4x4;
				float rotationAngle;

				CalculateTransformationMatrices( iTriangles.at(i), iTriangles.at(j), *vertexNeighbor[0], *vertexNeighbor[1], alignToAxis, alignToAxisInv, rotationAngle );
				Neighbor* neighborj = new Neighbor( &iTriangles.at(j), vertexNeighbor, alignToAxis, alignToAxisInv, rotationAngle );
				this->iTriangles.at(i).AddEdgeNeightbor( neighborj, edgePosition );

				// swap edge points for neighbor
				Point* tempPoint = vertexNeighbor[0];
				vertexNeighbor[0] = vertexNeighbor[1];
				vertexNeighbor[1] = tempPoint;

				Reference<Matrix4x4> tempMatrix = alignToAxis;
				alignToAxis = alignToAxisInv;
				alignToAxisInv = tempMatrix;

				Neighbor* neighbori = new Neighbor( &iTriangles.at(i), vertexNeighbor, alignToAxis, alignToAxisInv, rotationAngle );
				this->iTriangles.at(j).AddEdgeNeightbor( neighbori, edgePositionNeighbor );
			}
		}
	}
}

void PointRepulsion::CalculateTransformationMatrices( TriangleUseSet &aTriangle, TriangleUseSet& aNeighborTriangle,
													 Point& P0, Point& P1, Reference<Matrix4x4>& aArbitraryRotation,
													 Reference<Matrix4x4>& aArbitraryRotationInv, float& aRotationAngle )
{
	// Calculate the angle of rotation.
	// This is the angle between the surface of this triangle and the neighbor.
	Normal normal = aTriangle.GetNormal();
	Normal neighborNormal = aNeighborTriangle.GetNormal();

	float dotProduct = Dot( normal, neighborNormal);

	// clamping the dot product
	if (dotProduct>1.0f) dotProduct = 1.0f;
	if (dotProduct<-1.0f) dotProduct = -1.0f;

	aRotationAngle = Degrees( acos( dotProduct ) );
	if (aRotationAngle < 0.03f)
	{
		return;
	}

	// Calculate the edge vector using the normals. This will make sure the vector points in the right direction
	Vector V = Normalize( Cross( neighborNormal, Vector(normal) ) );

	//rotate around arbitrary vector V
	Rotate(aRotationAngle, V, aArbitraryRotation, aArbitraryRotationInv);

#ifdef DEBUG_POINTREPULSION
	Vector mass = Vector(P0);
	Transform testTransform = Translate( -mass )*Transform(aArbitraryRotation, aArbitraryRotationInv)*Translate( mass );
	Normal testNormal = Normalize( testTransform( neighborNormal ) );
	float testValue = Dot( testNormal, normal );
	if (!(testValue > 0.9995))
	{
		printf("Testvalue = %f\n", testValue);
	}
	assert( testValue > 0.9995 );
#endif
}

/*!
 * \brief
 * Add stratified sample points to the triangle mesh using a
 * list of partial sums of the polygon areas
 *
 * returns number of sample points created on the surface of the object
 */
int PointRepulsion::SetupSamplePoints( float aMeanFreePath )
{
	float numberOfSamplePoints = GetTotalSurfaceArea()/(M_PI*aMeanFreePath*aMeanFreePath);

	int xSamples = static_cast<int>( ceil( sqrt( numberOfSamplePoints ) ) );
	int ySamples = static_cast<int>( round( sqrt( numberOfSamplePoints ) ) );
//	int ySamples = xSamples;

	int iNumberOfSamplePoints = xSamples*ySamples;
	float* samples = new float[iNumberOfSamplePoints*2];

	StratifiedSample2D(samples, xSamples, ySamples, true);

	printf("Area = %f\nNumber Of Sample Points = %d\n", GetTotalSurfaceArea(), iNumberOfSamplePoints);
	for (int i=0; i<iNumberOfSamplePoints*2; i+=2)
	{
		//if (i%1000000 == 0)
		//{
		//	printf("Added %d sample points\n", i);
		//}
		ComputeSamplePointPosition(iPartialAreaSums , samples[i], samples[i+1]);
	}

	//// create testing sample points
	//Point* samplePoint = new Point( 1.0f, 0.1f, -0.5f);

	//iSamplePoints.push_back( samplePoint );

	//// constructor of PRSamplePointContainer will add a reference of samplePoint to the triangle the point is in
	//SamplePointContainer* samplePointContainer = new SamplePointContainer( samplePoint, &iTriangles.at(5) );
	//this->iSamplePointContainer.push_back( samplePointContainer );

	//samplePoint = new Point( 1.0f, 0.9f, -0.5f);

	//iSamplePoints.push_back( samplePoint );

	//// constructor of PRSamplePointContainer will add a reference of samplePoint to the triangle the point is in
	//samplePointContainer = new SamplePointContainer( samplePoint, &iTriangles.at(4) );
	//this->iSamplePointContainer.push_back( samplePointContainer );

	delete[] samples;


	// compute the repulsive range r. This is based on the total area surface of the object
	// and the amount of sample points (which has just been determined).
	iRepulsiveRadius = 2*sqrt( (iTotalAreaOfSurface/iNumberOfSamplePoints) );

	return iNumberOfSamplePoints;
}

void PointRepulsion::ComputePartialAreaSum( vector<float>& aArea )
{
	vector<TriangleUseSet>::iterator triangleIter = this->iTriangles.begin();
	int i = 0;

	float areaSum = 0;

	// compute relative areas of the sub-triangles of polygon
	for( ;triangleIter != iTriangles.end(); triangleIter++ )
	{
		areaSum += (*triangleIter).GetArea();
		aArea.push_back( areaSum );
		i++;
	}

	// normalize areas so that the sum of all sub-triangles is one
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		aArea[i] /= iTotalAreaOfSurface;
	}
}


/*!
 * \brief
 * Calculates the positions of sample points 
 * 
 * \param aAreas
 * An array containing the partial area sums
 * 
 * \param s
 * random position between 0 and 1
 * 
 * \param t
 * random position between 0 and 1
 * 
 */
void PointRepulsion::ComputeSamplePointPosition( vector<float>& aPartialAreaSums, float s, float t )
{ 
	int i;
	float areaSum = 0;
	float a,b,c;

	// use 's' to pick one sub-triangle, weighted by relative
	// area of triangles
	//for (i = 0; i < iNumberOfTriangles; i++)
	//{
	//	areaSum += aAreas[i];
	//	if (areaSum >= s)
	//	{
	//		break;
	//	}
	//}

	i = SearchForTriangle(s, aPartialAreaSums, 0, iNumberOfTriangles/2, iNumberOfTriangles );
	float normalizedTriangleArea = iTriangles.at(i).GetArea()/iTotalAreaOfSurface;

	/* map 's' into the interval [0,1] */
	s = (s - aPartialAreaSums[i] + normalizedTriangleArea) / normalizedTriangleArea;

	/* map (s,t) to a point in that sub-triangle */
	t = sqrt(t);

	a = 1 - t;
	b = (1 - s) * t;
	c = s * t;

	Point* samplePoint = new Point;
	Point* vertex0 = iTriangles.at(i).iVertices[0];
	Point* vertex1 = iTriangles.at(i).iVertices[1];
	Point* vertex2 = iTriangles.at(i).iVertices[2];

	samplePoint->x = a * vertex0->x + b * vertex1->x + c * vertex2->x;
	samplePoint->y = a * vertex0->y + b * vertex1->y + c * vertex2->y;
	samplePoint->z = a * vertex0->z + b * vertex1->z + c * vertex2->z;

	iSamplePoints.push_back( samplePoint );

	// constructor of PRSamplePointContainer will add a reference of samplePoint to the triangle the point is in
	SamplePointContainer* samplePointContainer = new SamplePointContainer( samplePoint, &iTriangles.at(i) );
	this->iSamplePointContainer.push_back( samplePointContainer );
}

void PointRepulsion::ComputeRepulsiveForces( const float& aForceScale )
{
	vector<TriangleUseSet>::iterator triangleIterator = iTriangles.begin();

	while ( triangleIterator != iTriangles.end() )
	{
		Reference<Matrix4x4> transform = new Matrix4x4;
		vector<bool> triangleMapped(iNumberOfTriangles, false);

		MapSamplePointsToPlane(  *triangleIterator, *triangleIterator, transform, triangleMapped, aForceScale );

		triangleIterator++;
	}
}

/*!
 * \brief
 * Rotating sample points to make them coplanar with the main triangle
 * 
 * \param aCurrentTriangle
 * The current triangle that is being evaluated
 * 
 * \param aMainTriangle
 * Description of parameter aMainTriangle.
 *
 * \param aEdgeRotationMatrix
 * the composed matrix which will rotate the sample points of the 
 * current triangle to the plane of the main triangle
 * 
 * \param aTriangleMapped
 * A list of triangles that have already been evaluated. This prevents triangles from being re-evaluated
 * 
 */
void PointRepulsion::MapSamplePointsToPlane(TriangleUseSet& aCurrentTriangle, 
											TriangleUseSet& aMainTriangle,
											Reference<Matrix4x4> aEdgeRotationMatrix,
											vector<bool>& aTriangleMapped,
											const float& aForceScale )
{
	// repel the sample points on the sample point Triangle
	ComputeRepulsiveForces(aCurrentTriangle, aMainTriangle, aEdgeRotationMatrix , aForceScale);

	// mark this triangle as "done". The points on this triangle have repelled the sample points on this triangle
	aTriangleMapped[aCurrentTriangle.GetTriangleId()] = true;

	vector<Neighbor*> neighbors;
	aCurrentTriangle.GetEdgeNeighbors( neighbors );
	vector<Neighbor*>::iterator neighborIter = neighbors.begin();
	vector<Neighbor*>::iterator neighborIterEnd = neighbors.end();
	while( neighborIter != neighborIterEnd )
	{
		//pointers to neighboring edge triangles of the current triangle that is being evaluated
		TriangleUseSet* neighborTriangle = (*neighborIter)->iEdgeNeighbor;

		// if the triangle has already been evaluated,
		// then it means that the sample points have already been rotated
		// and their repulsive forces have already been computed and stored.
		if (false == aTriangleMapped[neighborTriangle->GetTriangleId()])
		{	
			// if distance between the center of the neighboring triangle and 
			//   center of main triangle > radius repelling force
			// no need to continue checking neighbors
			Vector range = aMainTriangle.GetCentroid()-neighborTriangle->GetCentroid();
			if (range.Length() < iRepulsiveRadius)
			{
	
				// the edgeRotationMatrix concatenates the rotation and translation matrices which is used
				// to rotate the neigbors triangle to make it coplanar with the main triangle.
				// The code simply looks like:
				// edgeRotationMatrix = aEdgeRotationMatrix*TranslateToAxisInv*ArbitraryRotation*TranslationToAxis
				Reference<Matrix4x4> edgeRotationMatrix = Matrix4x4::Mul( aEdgeRotationMatrix, 
															Matrix4x4::Mul( (*neighborIter)->iTranslateToAxis,
																Matrix4x4::Mul( (*neighborIter)->iArbitraryRotation,
																		(*neighborIter)->iTranslateToAxisInv ) ) );
#ifdef DEBUG_POINTREPULSION
				Transform transform = Transform(edgeRotationMatrix, edgeRotationMatrix->Transpose());
				Vector mass = Vector( *(*neighborIter)->iP0 );
				Normal neighborNormal = neighborTriangle->GetNormal();
				Normal testNormal = Normalize( transform( neighborNormal ) );
				float testValue = Dot( testNormal, aMainTriangle.GetNormal() );
				assert( testValue > 0.9995 );

				Point P0 = transform( *((*neighborTriangle).iVertices[0]) );
				Point P1 = transform( *((*neighborTriangle).iVertices[1]) );
				Point P2 = transform( *((*neighborTriangle).iVertices[2]) );

				assert( VeryVerySmallDistancePointToPlane(aMainTriangle.GetNormal(), *aMainTriangle.iVertices[0], P0) );
				assert( VeryVerySmallDistancePointToPlane(aMainTriangle.GetNormal(), *aMainTriangle.iVertices[0], P1) );
				assert( VeryVerySmallDistancePointToPlane(aMainTriangle.GetNormal(), *aMainTriangle.iVertices[0], P2) );
#endif
				// start evaluating the neighbor triangle
				MapSamplePointsToPlane(*neighborTriangle, aMainTriangle, edgeRotationMatrix, aTriangleMapped, aForceScale );
			}
		}

		neighborIter++;
	}
}

/*!
 * \brief
 * Calculate the repulsive force each sample point laying on ACurrentTriangle 
 * exerts on each sample point laying on aMainTriangle.
 * 
 * \param aCurrentTriangle
 * the triangle containing the points which will be evaluated. The sample points on this
 * triangle will exert forces
 * 
 * \param aMainTriangle
 * Triangle containing sample points receiving repulsive forces.
 * 
 * \param aEdgeRotationMatrix
 * Matrix used to rotate the sample points laying on the aCurrentTriangle to the plane of
 * in which aMainTriangle lies.
 *
 */
void PointRepulsion::ComputeRepulsiveForces(TriangleUseSet& aCurrentTriangle,
											TriangleUseSet& aMainTriangle,
											Reference<Matrix4x4> aEdgeRotationMatrix,
											const float& aForceScale)
{
	list<SamplePointContainer* >::iterator mainTrianglePointIter = aMainTriangle.iSamplePoints.begin();
	list<SamplePointContainer* >::iterator currentTrianglePointIter;

	while ( mainTrianglePointIter != aMainTriangle.iSamplePoints.end() )
	{	
		currentTrianglePointIter = aCurrentTriangle.iSamplePoints.begin();

		while( currentTrianglePointIter != aCurrentTriangle.iSamplePoints.end())
		{
			// check if the sample point and repelling point are not the same
			if ( (*mainTrianglePointIter) != (*currentTrianglePointIter) )
			{

				Transform transform(aEdgeRotationMatrix);
				Point currentSamplePoint = transform( *((*currentTrianglePointIter)->GetSamplePoint()) );
				Point mainSamplePoint = *((*mainTrianglePointIter)->GetSamplePoint());

#ifdef DEBUG_POINTREPULSION
				assert( VeryVerySmallDistancePointToPlane(aMainTriangle.GetNormal(), mainSamplePoint, currentSamplePoint) );
#endif

				// calculate the direction of the force
				Vector v(mainSamplePoint - currentSamplePoint);

				// calculate the repulsive force. Falls off linearly with distance between the two points.
				float distance = v.Length();
				float repulsiveForce = 1 - (distance/iRepulsiveRadius);
				if (repulsiveForce > 0.0f)
				{
					(*mainTrianglePointIter)->AddForceVector( Normalize(v)*repulsiveForce*aForceScale );
				}
			}

			currentTrianglePointIter++;
		}

		mainTrianglePointIter++;
	}
}

void PointRepulsion::ComputeNewPositions()
{
	// iterate through all sample points
	vector<SamplePointContainer* >::iterator samplePointIter = iSamplePointContainer.begin();

	for(; samplePointIter != iSamplePointContainer.end(); samplePointIter++ )
	{
		SamplePointContainer* container = *samplePointIter;
		Point* p = container->GetSamplePoint();
		Point pPrime = *p+container->GetForceVector();
		TriangleUseSet* useSet = container->GetTriangle();
	
		ComputeNewPositions(pPrime, useSet, p, container );

		container->ResetForceVector();
	}
}

/*!
 * \brief
 * Moves the point using the computed force values
 */
void PointRepulsion::ComputeNewPositions(Point& aPPrime, TriangleUseSet* aUseSet, Point* p, SamplePointContainer* container)
{
	Neighbor* nextNeighbor = NULL;
	Point p0, p1, edgePoint;

	// find out if the "new" position of the sample point lies inside this triangle or outside.
	// if it lies outside hint of the direction of the new point are provided through:
	// the neighboring triangle where the point might lay in
	// the edge (p0 and p1) which is intersected
	// the point of intersection with the edge
	bool inside = PointInsideTriangle( aPPrime, *p, *aUseSet, &nextNeighbor, p0, p1, edgePoint );

#ifdef DEBUG_POINTREPULSION
	VeryVerySmallDistancePointToPlane(aUseSet->GetNormal(), *aUseSet->iVertices[0], edgePoint);
#endif
	// if the point is inside the triangle, point p becomes p' 
	if (inside)
	{
		*p = aPPrime;
		container->SetTriangle( aUseSet );
#ifdef DEBUG_POINTREPULSION
		assert( VeryVerySmallDistancePointToPlane(aUseSet->GetNormal(), *aUseSet->iVertices[0], *p) );
#endif
	}
	else // the point lies outside the triangle.
		 // We should by now have the neighbor p' might lay in, the intersection of the vector p->p' and with the edge
	{
		// The point is pushed over an edge so we compute the intersection point
		// shortening the begin point
		*p = edgePoint;

		if (NULL != nextNeighbor)
		{
			// calculate a displacement vector from the point on the edge to pPrime
			// and rotate the displacement vector so it's aligned with the neighboring triangle
			Vector displacement( aPPrime-*p );

			// rotate the displacement vector onto the triangle
			Transform transform = Transform( nextNeighbor->iArbitraryRotationInv, nextNeighbor->iArbitraryRotation );
			displacement = transform(displacement);

			// calculate a new pPrime using the point and a vector and 
			aPPrime = *p+displacement;

#ifdef DEBUG_POINTREPULSION
			Normal normal = aUseSet->GetNormal();
			Normal testNormal = Normalize( transform( normal ) );
			float testValue = Dot( testNormal, normal );
#endif

			ComputeNewPositions( aPPrime, nextNeighbor->iEdgeNeighbor, p, container);
		}
		else
		{
			container->SetTriangle( aUseSet );
#ifdef DEBUG_POINTREPULSION
			assert( VeryVerySmallDistancePointToPlane(aUseSet->GetNormal(), *aUseSet->iVertices[0], *p) );
#endif
		}
	}
}


/*!
 * \brief
 * Returns true if the point falls inside the triangle, otherwise it will be false.
 * In case it's outside, it stores the neighboring triangle in aNeighbor, the edge Points and the edge intersection point.
 * 
 */
bool PointRepulsion::PointInsideTriangle( const Point& aPoint, const Point& aStartPoint, const TriangleUseSet& aTriangle, 
										 Neighbor** aNeighbor, Point& aP0, Point& aP1, Point& aEdgePoint )
{
	Point a = *aTriangle.iVertices[0];
	Point b = *aTriangle.iVertices[1];
	Point c = *aTriangle.iVertices[2];

	Normal normal = aTriangle.GetNormal();

	vector<Neighbor*> neighbors;
	aTriangle.GetAllEdgeNeighbors(neighbors);

	Vector ab = b - a;
	Normalize(ab);
	Vector ap = aPoint - a;
	Normalize(ap);
	float f = Dot( Cross( ab, ap ), normal );

	// the point falls outside the triangle. Possibly it intersects this edge.
	if (f < 0)
	{
		// figure out if it intersects this edge of the triangle
		if ( LineSegLineSegIntersection(a, b, aStartPoint, aPoint, aEdgePoint) )
		{
			aP0 = *(aTriangle.iVertices[0]);
			aP1 = *(aTriangle.iVertices[1]);
			*aNeighbor = neighbors[0];
			return false;
		}
	}

	Vector bc = c - b;
	Normalize(bc);
	Vector bp = aPoint - b;
	Normalize(bp);
	float s = Dot( Cross( bc, bp ), normal );

	// the point falls outside the triangle.  Possibly it intersects this edge.
	if ( s < 0 )
	{
		// figure out if it intersects this edge of the triangle
		if ( LineSegLineSegIntersection(b, c, aStartPoint, aPoint, aEdgePoint) )
		{
			aP0 = *(aTriangle.iVertices[1]);
			aP1 = *(aTriangle.iVertices[2]);
			*aNeighbor = neighbors[1];
			return false;
		}
	}

	Vector ca = a - c;
	Normalize(ca);
	Vector cp = aPoint - c;
	Normalize(cp);
	float t = Dot( Cross( ca, cp ), normal );

	// the point falls outside the triangle. Possibly it intersects this edge.
	if ( t < 0 )
	{
		// figure out if it intersects this edge of the triangle
		if ( LineSegLineSegIntersection(a, c, aStartPoint, aPoint, aEdgePoint) )
		{
			aP0 = *(aTriangle.iVertices[2]);
			aP1 = *(aTriangle.iVertices[0]);
			*aNeighbor = neighbors[2];
			return false;
		}
	}

	aNeighbor = NULL;
	return true;
}

bool PointRepulsion::LineSegLineSegIntersection(const Point& x1, const Point& x2, const Point& x3, const Point& x4, Point& p )
{
	p = LineLineIntersection(x1, x2, x3, x4);

	return PointBetweenTwoPoints(x1, x2, p);
}

bool PointRepulsion::PointBetweenTwoPoints(const Point& x1, const Point x2, const Point& testPoint)
{
	float minX = min(x1.x, x2.x);
	float maxX = max(x1.x, x2.x);
	float minY = min(x1.y, x2.y);
	float maxY = max(x1.y, x2.y);
	float minZ = min(x1.z, x2.z);
	float maxZ = max(x1.z, x2.z);

	return ( (minX <= testPoint.x) && (testPoint.x <= maxX) &&
			 (minY <= testPoint.y) && (testPoint.y <= maxY) &&
			 (minZ <= testPoint.z) && (testPoint.z <= maxZ) );
}

Point PointRepulsion::LineLineIntersection( const Point& x1,const Point& x2,const Point& x3,const Point& x4  )
{
	Vector a(x2-x1);
	Vector b(x4-x3);
	Vector c(x3-x1);

	Vector crossAB = Cross(a,b);
	float numerator = Dot( Cross(c,b), crossAB );
	float denominator = crossAB.LengthSquared();

	return Point( x1+a*( numerator / denominator ) );
}


/**
* @description ugly function that does a dynamic cast on an object. 
* Hammer method used at full strength ;) in PBRT
*
* @param aTriangle
* @return the object or NULL if is not of the return type
*/
inline Triangle* PointRepulsion::Cast(Reference<Shape>& aTriangle)
{
	Shape* triangle = aTriangle.operator ->();
	return dynamic_cast<Triangle*> (triangle);
}

void PointRepulsion::FillUniformSamplePointStructure( vector<UniformPoint>& container )
{
	vector<SamplePointContainer *>::iterator it= iSamplePointContainer.begin();

	for (;it != iSamplePointContainer.end(); it++)
	{
		Point* point =(*it)->GetSamplePoint();
		TriangleUseSet* triangleUseSet = (*it)->GetTriangle();
		Reference<Shape> triangle = triangleUseSet->GetTriangle();

		UniformPoint uniformPoint;
		uniformPoint.p = *point;
		uniformPoint.n = triangleUseSet->GetNormal();
		uniformPoint.triangle= triangleUseSet->GetTriangle();

		container.push_back(uniformPoint);
	};
}

float PointRepulsion::GetTotalSurfaceArea()
{
	return iTotalAreaOfSurface;
}

int PointRepulsion::SearchForTriangle( const float aS, const vector<float>& aAreas, int min, int i, int max  )
{
	// hit the end of the vector of partial area sums
	if ( (i+1) >= iNumberOfTriangles)
	{
		return iNumberOfTriangles-1;
	}
	
	// hit the beginning of the vector of partial area sums
	if ( i <= 0)
	{
		return i;
	}

	// if current partial area sum is bigger than s...
	if ( aS < aAreas[i]	)
	{
		// check if the partial area sum on the left side is smaller than or equal to s
		if (aS >= aAreas[i-1] )
		{
			return i;
		}
		else 
		{
			max = i;
			i = static_cast<int>( min+(max-min)/2.0f );
			// the partial area sum on the left side is not smaller or equal to s, check futher on the left side
			return SearchForTriangle(aS, aAreas, min, i, max  );
		}
	}
	else
	{
		min = i;
		i = static_cast<int>( min+(max-min)/2.0f );
		//this partial area sum is smaller or equal to s, check on the right side
		return SearchForTriangle(aS, aAreas, min, i, max );
	}
}

#ifdef DEBUG_POINTREPULSION

bool PointRepulsion::VeryVerySmallDistancePointToPlane(Normal& aNormal, Point& aPoint, Point& aTestPoint )
{
	float a = aNormal.x;
	float b = aNormal.y;
	float c = aNormal.z;
	float d = -(Dot(aNormal, Vector(aPoint) ));

	float distance = ( abs( a*aTestPoint.x + b*aTestPoint.y + c*aTestPoint.z + d ) ) / ( aNormal.Length() );

	bool status = ( distance<0.0001 );

	if (status == false)
	{
		int i = 0;
	}
	return status;
}

#endif
#include "pointrepulsion.h"
#include "core/sampling.h"
#include "transform.h"
#include "TriangleMesh.h"

#include <iostream>
#include <fstream>
using namespace std;

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
	printf("Number Of Triangles: %d\n", iNumberOfTriangles);

	// Setup the triangle meta data:
	// - Total Surface Area
	// - Normalized Partial Surface Area
	// - Structure (find edge neighbors)
	SetupTriangleUseSets( aTriangleList );
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
	CreateTriangleUseSets(aTriangleList);
	LinkTriangleUseSets();
}

void PointRepulsion::CreateTriangleUseSets(vector<Reference<Shape> >& aTriangleList)
{
	ProgressReporter progress(iNumberOfTriangles, "Creating Triangle Usesets");
	int* vertices = NULL;

	//set up all the triangle use sets
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		Triangle* triangle = TriangleMesh::Cast(aTriangleList.at(i));
		triangle->GetVertexIndices( &vertices );

		TriangleUseSet triangleUseSet( aTriangleList.at(i), (iVertices+vertices[0]), (iVertices+vertices[1]), (iVertices+vertices[2]) );
		this->iTriangles.push_back( triangleUseSet );

		// calculate the total surface area;
		iTotalAreaOfSurface += triangleUseSet.GetArea();

		// compute relative areas of the sub-triangles of polygon
		iPartialAreaSums.push_back(iTotalAreaOfSurface);
		progress.Update();
	}
	progress.Done();
}

void PointRepulsion::LinkTriangleUseSets()
{
	ProgressReporter progress(iNumberOfTriangles, "Linking Triangle Usesets");

	int* vertices1 = NULL;

	// loop through all the indices to find triangle neighbors
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		// normalize areas so that the sum of all sub-triangles is one
		iPartialAreaSums[i] /= iTotalAreaOfSurface;

		vector<TriangleEdge*> neighbor;
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
				TriangleEdge* neighborj = new TriangleEdge( &iTriangles.at(j), vertexNeighbor, alignToAxis, alignToAxisInv, rotationAngle );
				this->iTriangles.at(i).AddEdgeNeightbor( neighborj, edgePosition );

				// swap edge points for neighbor
				Point* tempPoint = vertexNeighbor[0];
				vertexNeighbor[0] = vertexNeighbor[1];
				vertexNeighbor[1] = tempPoint;

				Reference<Matrix4x4> tempMatrix = alignToAxis;
				alignToAxis = alignToAxisInv;
				alignToAxisInv = tempMatrix;

				TriangleEdge* neighbori = new TriangleEdge( &iTriangles.at(i), vertexNeighbor, alignToAxis, alignToAxisInv, rotationAngle );
				this->iTriangles.at(j).AddEdgeNeightbor( neighbori, edgePositionNeighbor );
			}
		}

		progress.Update();
	}

	progress.Done();
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
	float numberOfSamplePoints = (GetTotalSurfaceArea()/(M_PI*aMeanFreePath*aMeanFreePath));
	printf("Number of sample points: %f\n", numberOfSamplePoints);
	int xSamples = static_cast<int>( ceil( sqrt( numberOfSamplePoints ) ) );
	int ySamples = static_cast<int>( round( sqrt( numberOfSamplePoints ) ) );
//	int ySamples = xSamples;

	//int iNumberOfSamplePoints = xSamples*ySamples;
	//float* samples = new float[iNumberOfSamplePoints*2];
	//StratifiedSample2D(samples, xSamples, ySamples, false);

	float* samples = new float[numberOfSamplePoints*2];
//	LatinHypercube(samples, numberOfSamplePoints, 2);
	LDShuffleScrambled2D(numberOfSamplePoints, 1, samples);

	int iNumberOfSamplePoints = numberOfSamplePoints;


	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
	//delete [] samples;

	//iNumberOfSamplePoints = 3;
	//samples = new float[iNumberOfSamplePoints*2];
	//samples[0] = 0.35f;
	//samples[1] = 0.6f;

	//samples[2] = 0.45f;
	//samples[3] = 0.3f;

	//samples[4] = 0.8f;
	//samples[5] = 0.6f;

	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
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

#ifdef POINTREPULSION_PRINTASCII
	fstream outfile("stratified2D.txt", fstream::out);
	outfile << "1" << endl;
	outfile << iNumberOfSamplePoints << endl;

	for (int i=0; i<iNumberOfSamplePoints*2; i+=2)
	{
		outfile << samples[i] << " " << samples[i+1] << " " << 0.0 << " ";
	}
	outfile.close();
#endif
	//vector<Point3d> points;
	//ParseFile( "stratified2D.txt", points );

	printf("Area = %f\nNumber Of Sample Points = %d\n", GetTotalSurfaceArea(), iNumberOfSamplePoints);
	ProgressReporter progress(iNumberOfTriangles, "Throwing sample points on the triangles");
	for (int i=0; i<iNumberOfSamplePoints*2; i+=2)
	{
		ComputeSamplePointPosition(iPartialAreaSums , samples[i], samples[i+1]);
		progress.Update();
	}
	progress.Done();

	delete[] samples;


	// compute the repulsive range r. This is based on the total area surface of the object
	// and the amount of sample points (which has just been determined).
//	iRepulsiveRadius = sqrt(iTotalAreaOfSurface/(M_PI*iNumberOfSamplePoints));
	iRepulsiveRadius = 2*sqrt( (iTotalAreaOfSurface/iNumberOfSamplePoints) );

	printf("iRepulsiveRadius = %f\n", iRepulsiveRadius);
	return iNumberOfSamplePoints;
}


/*!
 * \brief
 * Calculates the positions of sample points from 2d plane to 3d mesh
 *  
 */
void PointRepulsion::ComputeSamplePointPosition( vector<float>& aPartialAreaSums, float s, float t )
{ 
	int i;
	float areaSum = 0;
	float a,b,c;

	i = SearchForTriangle(s, aPartialAreaSums, 0, iNumberOfTriangles-1 );
	float normalizedTriangleArea = iTriangles.at(i).GetArea()/iTotalAreaOfSurface;

	// map 's' into the interval [0,1]
	s = (s - aPartialAreaSums[i] + normalizedTriangleArea) / normalizedTriangleArea;

	// map (s,t) to a point in that sub-triangle
	t = sqrt(t);

	a = 1 - t;
	b = (1 - s) * t;
	c = s * t;

	Point* samplePoint = new Point;
	Point* vertex0 = iTriangles.at(i).iVertices[0];
	Point* vertex1 = iTriangles.at(i).iVertices[1];
	Point* vertex2 = iTriangles.at(i).iVertices[2];

	// calculate the sample point position using the vertices of the triangle the weights connected to it.
	samplePoint->x = a * vertex0->x + b * vertex1->x + c * vertex2->x;
	samplePoint->y = a * vertex0->y + b * vertex1->y + c * vertex2->y;
	samplePoint->z = a * vertex0->z + b * vertex1->z + c * vertex2->z;

	iSamplePoints.push_back( samplePoint );

	// constructor of PRSamplePointContainer will add a reference of samplePoint to the triangle the point is in
	SamplePointContainer* samplePointContainer = new SamplePointContainer( samplePoint, &iTriangles.at(i) );
	this->iSamplePointContainer.push_back( samplePointContainer );
}

void PointRepulsion::ComputeRepulsiveForces( const float& aForceScale, const string& aProcessString )
{
	ProgressReporter progress(iNumberOfTriangles, aProcessString );
	vector<TriangleUseSet>::iterator triangleIterator = iTriangles.begin();

	list< pair<TriangleUseSet*, Reference<Matrix4x4> > > triangleQueue;
	while ( triangleIterator != iTriangles.end() )
	{

		Reference<Matrix4x4> transform = new Matrix4x4;
		vector<bool> triangleMapped(iNumberOfTriangles, false);

		triangleMapped[(*triangleIterator).GetTriangleId()] = true;

		triangleQueue.push_back( make_pair(&(*triangleIterator), transform) );

		while( !triangleQueue.empty() )
		{
			MapSamplePointsToPlane( triangleQueue, *triangleIterator, triangleMapped, aForceScale );
		}

		triangleQueue.clear();
		triangleIterator++;
		progress.Update();
	}
	progress.Done();
}
/*!
 * \brief
 * Rotating sample points to make them coplanar with the aMain triangle
 * 
 * \param aTriangleQueue
 * Description of parameter aTriangleQueue.
 * 
 */
void PointRepulsion::MapSamplePointsToPlane(list< pair<TriangleUseSet*, Reference<Matrix4x4> > >& aTriangleQueue,
											TriangleUseSet& aMainTriangle, vector<bool>& aTriangleMapped, 
											const float& aForceScale )
{
	pair<TriangleUseSet*, Reference<Matrix4x4> > trianglePair = aTriangleQueue.front();
	aTriangleQueue.pop_front();

	TriangleUseSet* evaluatedTriangle = trianglePair.first;
	Reference<Matrix4x4> aEdgeRotationMatrix = trianglePair.second;

	ComputeRepulsiveForces(*evaluatedTriangle, aMainTriangle, aEdgeRotationMatrix , aForceScale);

	vector<TriangleEdge*> neighbors;
	evaluatedTriangle->GetEdgeNeighbors( neighbors );

	vector<TriangleEdge*>::iterator neighborIter = neighbors.begin();
	vector<TriangleEdge*>::iterator neighborIterEnd = neighbors.end();

	for (; neighborIter != neighborIterEnd; neighborIter++)
	{
		//pointers to neighboring edge triangles of the current triangle that is being evaluated
		TriangleUseSet* neighborTriangle = (*neighborIter)->iEdgeNeighbor;

		// if the triangle has already been evaluated,
		// then it means that the sample points have already been rotated
		// and their repulsive forces have already been computed and stored.
		if (false == aTriangleMapped[neighborTriangle->GetTriangleId()])
		{
			Reference<Matrix4x4> edgeRotationMatrix = Matrix4x4::Mul( aEdgeRotationMatrix, 
														Matrix4x4::Mul( (*neighborIter)->iTranslateToAxis,
															Matrix4x4::Mul( (*neighborIter)->iArbitraryRotation,
																			(*neighborIter)->iTranslateToAxisInv ) ) );
			Transform transform(edgeRotationMatrix);

			float minimalDistance = ComputeSmallestDistanceBetweenTriangles(*neighborTriangle, aMainTriangle, transform);

			if (minimalDistance < iRepulsiveRadius)
			{
				// we add the triangle use set to the list to be evaluated and
				// mark it as evaluated so it won't be added again
				aTriangleQueue.push_back( make_pair(neighborTriangle, edgeRotationMatrix) );
			}
		}
		aTriangleMapped[neighborTriangle->GetTriangleId()] = true;
	}
}

float PointRepulsion::ComputeSmallestDistanceBetweenTriangles(TriangleUseSet& aEvaluatedTriangle, TriangleUseSet& aMainTriangle, Transform& aTransform)
{
	Point evalVertices[3];
	evalVertices[0] = aTransform(*aEvaluatedTriangle.iVertices[0]);
	evalVertices[1] = aTransform(*aEvaluatedTriangle.iVertices[1]);
	evalVertices[2] = aTransform(*aEvaluatedTriangle.iVertices[2]);

	Point* vertices[3];
	vertices[0] = aMainTriangle.iVertices[0];
	vertices[1] = aMainTriangle.iVertices[1];
	vertices[2] = aMainTriangle.iVertices[2];

	float minimal = INFINITY;

	Vector v;
	float length = 0.0f;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			v = evalVertices[i]-*(vertices[j]);
			length = v.Length();
			minimal = min(length, minimal);
		}
	}

	return minimal;
};

/*!
 * \brief
 * Calculate the repulsive force each sample point laying on ACurrentTriangle 
 * exerts on each sample point laying on aMainTriangle.
 * 
 * \param aEvaluatedTriangle
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
void PointRepulsion::ComputeRepulsiveForces(TriangleUseSet& aEvaluatedTriangle,
											TriangleUseSet& aMainTriangle,
											Reference<Matrix4x4> aEdgeRotationMatrix,
											const float& aForceScale)
{
	//bool pointIsRepelled = false;

	if (aEvaluatedTriangle.iSamplePoints.empty())
	{
		return;
	}

	list<SamplePointContainer* >::iterator mainTrianglePointIter = aMainTriangle.iSamplePoints.begin();
	list<SamplePointContainer* >::iterator currentTrianglePointIter;

	for (; mainTrianglePointIter != aMainTriangle.iSamplePoints.end();	mainTrianglePointIter++ )
	{	
		currentTrianglePointIter = aEvaluatedTriangle.iSamplePoints.begin();
		Point mainSamplePoint = *((*mainTrianglePointIter)->GetSamplePoint());

		for(; currentTrianglePointIter != aEvaluatedTriangle.iSamplePoints.end(); currentTrianglePointIter++)
		{
			// check if the sample point and repelling point are not the same
			if ( (*mainTrianglePointIter) != (*currentTrianglePointIter) )
			{

				Transform transform(aEdgeRotationMatrix);
				Point evaluatedSamplePoint = transform( *((*currentTrianglePointIter)->GetSamplePoint()) );

				Point tempPoint = evaluatedSamplePoint;
				//PushPointToPlane(aMainTriangle.GetNormal(), mainSamplePoint, tempPoint);
				//if ( !PointCloseToPlane(aMainTriangle.GetNormal(), mainSamplePoint, tempPoint) )
				//{
				//	int i = 0;
				//}

				// make sure the point is exactly on the plane
//				PushPointToPlane(aMainTriangle.GetNormal(), mainSamplePoint, evaluatedSamplePoint);

				//if (!PointCloseToPlane(aMainTriangle.GetNormal(), mainSamplePoint, evaluatedSamplePoint))
				//{
				//	PushPointToPlane(aMainTriangle.GetNormal(), mainSamplePoint, evaluatedSamplePoint);
				//}
				// calculate the direction of the force
				Vector v(mainSamplePoint - evaluatedSamplePoint);

				// calculate the repulsive force. Falls off linearly with distance between the two points.
				float distance = v.Length();
				float repulsiveForce = iRepulsiveRadius-distance;
				if (repulsiveForce > 0.0f)
				{
					(*mainTrianglePointIter)->AddForceVector( Normalize(v)*repulsiveForce*aForceScale );
				}
			}
		}
	}
}

void PointRepulsion::ComputeNewPositions( const string& aProcessString )
{
	ProgressReporter progress(iSamplePoints.size(), aProcessString);

	//static float averageOld = 0.0f;
	//float average = 0.0f;

	// iterate through all sample points
	vector<SamplePointContainer* >::iterator samplePointIter = iSamplePointContainer.begin();
	for(; samplePointIter != iSamplePointContainer.end(); samplePointIter++)
	{
		SamplePointContainer* container = *samplePointIter;
		Vector force = container->GetForceVector();

//		average += force.Length();

		Point* p = container->GetSamplePoint();
		Point pPrime = *p + force;

		bool reposition = true;
		while (reposition)
		{
			reposition = ComputeNewPositions( pPrime, p, container );
		}

		container->ResetForceVector();
		progress.Update();
	}
	progress.Done();

	//average = average/iSamplePointContainer.size();
	//printf("\nAverage force vector length = %f\n", average );
	//printf("Difference previous - current = %f\n", averageOld - average );
	//averageOld = average;
}

/*!
 * \brief
 * Moves the point using the computed force values
 */
bool PointRepulsion::ComputeNewPositions(Point& aPPrime, Point* p, SamplePointContainer* container)
{
	static int iterationCounter = 0;
	TriangleEdge* nextNeighbor = NULL;
	Point p0, p1, edgePoint;
	TriangleUseSet* useSet = container->GetTriangle();

	// find out if the "new" position of the sample point lies inside this triangle or outside.
	// if it lies outside, the direction of the new point are provided through:
	// the neighboring triangle where the point might lay in
	// the edge (p0 and p1) which is intersected
	// the point of intersection with the edge
	bool inside = PointInsideTriangle( aPPrime, *p, *useSet, &nextNeighbor, p0, p1, edgePoint );

//	PushPointToPlane(aUseSet->GetNormal(), *aUseSet->iVertices[0], edgePoint);
	// if the point is inside the triangle, point p becomes p' 
	if (inside)
	{
		*p = aPPrime;

		PushPointToPlane(useSet->GetNormal(), *useSet->iVertices[0], *p);

		iterationCounter = 0;
	}
	else // the point lies outside the triangle.
		 // We should by now have the neighbor p' might lay in, the intersection of the vector p->p' and with the edge
	{
		iterationCounter++;

		if (iterationCounter > 100)
		{
			printf("Number of iterations: %d\n", iterationCounter);
		}
		// The point is pushed over an edge so we compute the intersection point
		// shortening the begin point
		*p = edgePoint;
		
		if (NULL != nextNeighbor)
		{
			// move the point to a new triangle;
			container->SetTriangle( nextNeighbor->iEdgeNeighbor );

			// calculate a displacement vector from the point on the edge to pPrime
			// and rotate the displacement vector so it's aligned with the neighboring triangle
			Vector displacement( aPPrime-*p );

			// rotate the displacement vector onto the triangle
			Transform transform = Transform( nextNeighbor->iArbitraryRotationInv, nextNeighbor->iArbitraryRotation );
			displacement = transform(displacement);

			// calculate a new pPrime using the point and a vector and 
			aPPrime = *p+displacement;

#ifdef DEBUG_POINTREPULSION
			Normal normal = useSet->GetNormal();
			Normal testNormal = Normalize( transform( normal ) );
			float testValue = Dot( testNormal, normal );
#endif
		}
		else // there's no neighbor, we're done with our repulsion for this point
		{
			inside = true;
		}
	}


	return !inside;
}


/*!
 * \brief
 * Returns true if the point falls inside the triangle, otherwise it will be false.
 * In case it's outside, it stores the neighboring triangle in aNeighbor, the edge Points and the edge intersection point.
 * 
 */
bool PointRepulsion::PointInsideTriangle( const Point& aPoint, const Point& aStartPoint, const TriangleUseSet& aTriangle, 
										 TriangleEdge** aNeighbor, Point& aP0, Point& aP1, Point& aEdgePoint )
{
	Point a = *aTriangle.iVertices[0];
	Point b = *aTriangle.iVertices[1];
	Point c = *aTriangle.iVertices[2];

	Normal normal = aTriangle.GetNormal();

	vector<TriangleEdge*> neighbors;
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


/*!
 * \brief
 * return true when there is an intersection between two line segments and returns the point of intersection
 *	else returns false
 * 
 */
bool PointRepulsion::LineSegLineSegIntersection(const Point& x1, const Point& x2, const Point& x3, const Point& x4, Point& p )
{
	p = LineLineIntersection(x1, x2, x3, x4);

	return PointBetweenTwoPoints(x1, x2, p);
}

/*!
 * \brief
 * Checks whether a point falls on the line between two other points
 * 
 */
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

/*!
 * \brief
 * Calculate the intersection point of two lines. The method assumes an intersection will occur at all times.
 * 
 */
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

int PointRepulsion::SearchForTriangle( const float aS, const vector<float>& aAreas, int lower, int upper  )
{
	int numberOfElements = (upper - lower);

	// hit the beginning of the vector of partial area sums
	if ( numberOfElements <= 1)
	{
		return (aS>aAreas[lower]) ? upper: lower;
	}

	int middle = lower+(numberOfElements)/2;

	if ( aS < aAreas[middle] )
	{
		// the partial area sum on the left side is not smaller or equal to s, check futher on the left side
		return SearchForTriangle(aS, aAreas, lower, middle  );
	}
	else
	{
		//this partial area sum is smaller or equal to s, check on the right side
		return SearchForTriangle(aS, aAreas, middle, upper );
	}
}

void PointRepulsion::PushPointToPlane( const Normal& aNormal, const Point& aPoint, Point& aTestPoint  )
{
	float a = aNormal.x;
	float b = aNormal.y;
	float c = aNormal.z;
	float d = -(Dot(aNormal, Vector(aPoint) ));

	float distance = ( a*aTestPoint.x + b*aTestPoint.y + c*aTestPoint.z + d ) / ( aNormal.Length() );

	if (distance < -0.0000001f)
	{
		aTestPoint += Vector(aNormal)*distance;
	}

	if (distance > 0.0000001f)
	{
		aTestPoint -= Vector(aNormal)*distance;
	}
}

bool PointRepulsion::PointCloseToPlane( const Normal& aNormal, const Point& aPoint, Point& aTestPoint  )
{
	float a = aNormal.x;
	float b = aNormal.y;
	float c = aNormal.z;
	float d = -(Dot(aNormal, Vector(aPoint) ));

	float distance = ( abs( a*aTestPoint.x + b*aTestPoint.y + c*aTestPoint.z + d ) ) / ( aNormal.Length() );

	return (!(distance>0.0001));
}
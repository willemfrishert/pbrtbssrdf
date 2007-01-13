#include "pointrepulsion.h"
#include "core/sampling.h"
#include "transform.h"

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

PointRepulsion::PointRepulsion(int aNumberOfTriangles, int aNumberOfVertices, int* aVertexIndex,
							   Point* aVertices, int aNumberOfSamplePoints )
: iNumberOfTriangles( aNumberOfTriangles )
, iNumberOfVertices( aNumberOfVertices )
, iVertexIndex( aVertexIndex )
, iVertices( aVertices)
, iNumberOfSamplePoints( aNumberOfSamplePoints )
{
	// setup the triangle structure (find edge neighbors)
	SetupTriangles();

	// add stratified sample points to the triangle mesh using a
	// list of partial sums of the polygon areas
	// (Also calculates a total area sum)
	SetupSamplePoints();

	// compute the repulsive range r
	ComputeRepulsiveRange();

	// relaxation of the random points 
	for(int k = 0; k<40; k++)
	{
		RepelSamplePoints();
	}
}

PointRepulsion::~PointRepulsion()
{
	iSamplePointContainer.clear();

	iTriangles.clear();
};

void PointRepulsion::SetupTriangles()
{
	int* t;

	//set up all the triangle use sets
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		t = &iVertexIndex[3*i];

		TriangleUseSet triangle( (iVertices+t[0]), (iVertices+t[1]), (iVertices+t[2]) );
		this->iTriangles.push_back( triangle );
	}

	// loop through all the indices to find triangle neighbors
	for (int i=0; i<iNumberOfTriangles;i++)
	{
		if ( this->iTriangles.at(i).GetEdgeNeighbors().size() > 2)
		{
			continue;
		}

		t = &iVertexIndex[3*i];

		// go through the triangles i+1 (triangles that are 
		// before 1 have already registered as being neighbors).
		for (int j=i+1; j< iNumberOfTriangles; j++)
		{
			Point* vertexNeighbor[2] = { NULL };

			int vertexNeighborCount = 0;

			int* t2 =&iVertexIndex[3*j];

			// set the indices to see if they match
			for (int k=0; k<3;k++)
			{
				if (t[k] == t2[0])
				{
					vertexNeighbor[vertexNeighborCount] = (iVertices+t[k]);
					vertexNeighborCount++;
					continue;
				}
				if (t[k] == t2[1])
				{
					vertexNeighbor[vertexNeighborCount] = (iVertices+t[k]);
					vertexNeighborCount++;
					continue;
				}
				if (t[k] == t2[2])
				{
					vertexNeighbor[vertexNeighborCount] = (iVertices+t[k]);
					vertexNeighborCount++;
					continue;
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
				this->iTriangles.at(i).AddEdgeNeightbor( neighborj );

				// swap edge points for neighbor
				Point* tempPoint = vertexNeighbor[0];
				vertexNeighbor[0] = vertexNeighbor[1];
				vertexNeighbor[1] = tempPoint;

				Reference<Matrix4x4> tempMatrix = alignToAxis;
				alignToAxis = alignToAxisInv;
				alignToAxisInv = tempMatrix;

				Neighbor* neighbori = new Neighbor( &iTriangles.at(i), vertexNeighbor, alignToAxis, alignToAxisInv, rotationAngle );
				this->iTriangles.at(j).AddEdgeNeightbor( neighbori );
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
	Normal normal = aTriangle.ComputeNormal();
	Normal neighborNormal = aNeighborTriangle.ComputeNormal();

	aRotationAngle = Degrees( acos( Dot( normal, neighborNormal) ) );
	if (aRotationAngle == 0.0f)
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
	assert( testValue > 0.9995 );
#endif
}


void PointRepulsion::ComputeRepulsiveRange()
{
	iRepulsiveRadius = 2*sqrt( (iTotalAreaOfSurface/iNumberOfSamplePoints) );
}

void PointRepulsion::SetupSamplePoints()
{
	iNumberOfSamplePoints = 4;
	int xSamples = static_cast<int>( ceil( sqrt( static_cast<float>(iNumberOfSamplePoints) ) ) );
	int ySamples = xSamples;

	iNumberOfSamplePoints = xSamples*ySamples;
	float* samples = new float[iNumberOfSamplePoints*2];

	StratifiedSample2D(samples, xSamples, ySamples, true);

	float* area = new float[iNumberOfTriangles];
	ComputePartialAreaSum(area);

	for (int i=0; i<xSamples*ySamples*2; i+=2)
	{
		ComputeSamplePointPosition(area, samples[i], samples[i+1]);
	}

	delete[] samples;
	delete[] area;
	return;
}

void PointRepulsion::ComputePartialAreaSum( float* aArea )
{
	iTotalAreaOfSurface = 0;

	// compute relative areas of the sub-triangles of polygon
	for (int i = 0; i < iNumberOfTriangles; i++) 
	{
		aArea[i] = this->iTriangles.at(i).ComputeArea();
		iTotalAreaOfSurface += aArea[i];
	}

	// normalize areas so that the sum of all sub-triangles is one
	for (int i = 0; i < iNumberOfTriangles; i++)
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
void PointRepulsion::ComputeSamplePointPosition(float* aAreas, float s, float t)
{ 
	int i;
	float areaSum = 0;
	float a,b,c;

	/* use 's' to pick one sub-triangle, weighted by relative */
	/* area of triangles */
	for (i = 0; i < iNumberOfTriangles; i++)
	{
		areaSum += aAreas[i];
		if (areaSum >= s)
		{
			break;
		}
	}

	/* map 's' into the interval [0,1] */
	s = (s - areaSum + aAreas[i]) / aAreas[i];

	/* map (s,t) to a point in that sub-triangle */
	t = sqrt(t);

	a = 1 - t;
	b = (1 - s) * t;
	c = s * t;

	Point samplePoint;
	Point* vertex0 = iTriangles.at(i).iVertices[0];
	Point* vertex1 = iTriangles.at(i).iVertices[1];
	Point* vertex2 = iTriangles.at(i).iVertices[2];

	samplePoint.x = a * vertex0->x + b * vertex1->x + c * vertex2->x;
	samplePoint.y = a * vertex0->y + b * vertex1->y + c * vertex2->y;
	samplePoint.z = a * vertex0->z + b * vertex1->z + c * vertex2->z;

	iSamplePoints.push_back( samplePoint );

	// constructor of PRSamplePointContainer will add a reference of samplePoint to the triangle the point is in
	SamplePointContainer samplePointContainer( &samplePoint, &iTriangles.at(i) );
	this->iSamplePointContainer.push_back( samplePointContainer );
}

void PointRepulsion::RepelSamplePoints()
{
	vector<TriangleUseSet>::iterator triangleIterator = iTriangles.begin();

	while ( triangleIterator != iTriangles.end() )
	{
		Reference<Matrix4x4> transform = new Matrix4x4;
		vector<bool> triangleMapped(iNumberOfTriangles, false);

		MapSamplePointsToPlane( transform, triangleMapped, *triangleIterator, *triangleIterator );

		triangleIterator++;
	}
}

void PointRepulsion::MapSamplePointsToPlane(Reference<Matrix4x4> aTransform, 
											vector<bool>& aTriangleMapped,
											TriangleUseSet& aTriangle,
											TriangleUseSet& aSamplePointTriangle )
{
	// repel the sample points on the sample point Triangle
	ComputeRepulsiveForces(aTransform, aTriangle, aSamplePointTriangle);

	// mark this triangle as "done". The points on this triangle have repelled the sample points on this triangle
	aTriangleMapped[aTriangle.GetTriangleId()] = true;

	vector<Neighbor*> neighbors = aTriangle.GetEdgeNeighbors();
	vector<Neighbor*>::iterator neighborIter = neighbors.begin();
	vector<Neighbor*>::iterator endNeighborIter = neighbors.end();

	while( neighborIter != neighbors.end())
	{
		TriangleUseSet* neighborTriangle = (*neighborIter)->iEdgeNeighbor;

		// if distance between the center of this triangle and center of sample point triangle > radius repelling force
		// no need to continue checking neighbors
		Vector range = aTriangle.GetCentroid()-neighborTriangle->GetCentroid();
		if (range.Length() < iRepulsiveRadius)
		{
			if (false == aTriangleMapped[neighborTriangle->GetTriangleId()])
			{		
				// transform = TranslateToAxisInv*ArbitraryRotation*TranslationToAxis
				Reference<Matrix4x4> transform = Matrix4x4::Mul( (*neighborIter)->iTranslateToAxisInv , 
																	Matrix4x4::Mul( (*neighborIter)->iArbitraryRotation,
																		(*neighborIter)->iTranslateToAxis ) );

				MapSamplePointsToPlane(transform, aTriangleMapped, *neighborTriangle, aSamplePointTriangle);
			}
		}

		neighborIter++;
	}
}

void PointRepulsion::ComputeRepulsiveForces(Reference<Matrix4x4> aTransform,
											TriangleUseSet& aTriangle,
											TriangleUseSet& aSamplePointTriangle)
{
	list<SamplePointContainer* >::iterator sampleTrianglePointIter = aSamplePointTriangle.iSamplePoints.begin();
	list<SamplePointContainer* >::iterator trianglePointIter;

	while ( sampleTrianglePointIter != aSamplePointTriangle.iSamplePoints.end() )
	{	
		trianglePointIter = aTriangle.iSamplePoints.begin();

		while( trianglePointIter != aTriangle.iSamplePoints.end())
		{
			// check if the sample point and repelling point are not the same
			if ( (*sampleTrianglePointIter) != (*trianglePointIter) )
			{

				//Point point = aTransform( *((*trianglePointIter)->GetSamplePoint()) );
				//Point samplePoint = *((*sampleTrianglePointIter)->GetSamplePoint());
				//// calculate force using distance between points and repulsive radius
				//Vector v = point - samplePoint;
				//float distance = v.Length();

				//float force = 1 - (distance/iRepulsiveRadius);
				//if (force > 0.0f)
				//{
				//	(*sampleTrianglePointIter)->AddForceVector(v*force);
				//}
			}

			trianglePointIter++;
		}

		sampleTrianglePointIter++;
	}
}


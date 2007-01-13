#ifndef PBRT_POINTREPULSION_H
#define PBRT_POINTREPULSION_H

#include "geometry.h"
#include "triangleuseset.h"
#include "samplepointcontainer.h"

#include <vector>

using namespace std;

class PointRepulsion
{
public:
	PointRepulsion( int aNumberOfTriangles, int aNumberOfVertices, int* aVertexIndex,
					Point* aVertices, int aNumberOfSamplePoints );
	~PointRepulsion();
private:
	void SetupTriangles();
	void CalculateTransformationMatrices(TriangleUseSet &aTriangle, TriangleUseSet& aNeighborTriangle,
										 Point& P0, Point& P1, Reference<Matrix4x4>& aArbitraryRotation,
										 Reference<Matrix4x4>& aArbitraryRotationInv, float& aRotationAngle);
	void ComputeRepulsiveRange();

	void SetupSamplePoints();
	void ComputePartialAreaSum(float* aArea);
	void ComputeSamplePointPosition(float* aAreas, float s, float t);

	void RepelSamplePoints();
	void MapSamplePointsToPlane(Reference<Matrix4x4> aTransform, 
								vector<bool>& aTriangleMapped,
								TriangleUseSet& aTriangle,
								TriangleUseSet& aSamplePointTriangle);
	void ComputeRepulsiveForces(Reference<Matrix4x4> aTransform, TriangleUseSet& aTriangle, TriangleUseSet& aSamplePointTriangle);



	int iNumberOfTriangles;
	int iNumberOfVertices;
	int *iVertexIndex;
	Point *iVertices;
	int iNumberOfSamplePoints;
	float iRepulsiveRadius;
	float iTotalAreaOfSurface;

	vector<Point> iSamplePoints;
	vector<SamplePointContainer> iSamplePointContainer;
	vector<TriangleUseSet> iTriangles;
};


#endif
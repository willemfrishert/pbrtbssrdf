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

	void SetupSamplePoints();
	void ComputePartialAreaSum(float* aArea);
	void ComputeSamplePointPosition(float* aAreas, float s, float t);

	void ComputeRepulsiveForces();
	void MapSamplePointsToPlane(TriangleUseSet& aCurrentTriangle, TriangleUseSet& aMainTriangle, Reference<Matrix4x4> aEdgeRotationMatrix, vector<bool>& aTriangleMapped);
	void PointRepulsion::ComputeRepulsiveForces(TriangleUseSet& aCurrentTriangle,
												TriangleUseSet& aMainTriangle,
												Reference<Matrix4x4> aEdgeRotationMatrix);
	bool PointInsideTriangle(const Point& aPoint, const Point& aStartPoint, const TriangleUseSet& aTriangle, Neighbor** aNeighbor, Point& aP0, Point& aP1, Point& aEdgePoint);
	Point LineLineIntersection(const Point& x1,const Point& x2,const Point& x3,const Point& x4 );
	void ComputeNewPositions();
	void ComputeNewPositions(Point& aPPrime, const TriangleUseSet* aUseSet, Reference<Matrix4x4> aEdgeRotationMatrixInv, Point* p);
	bool LineSegLineSegIntersection(const Point& x1, const Point& x2, const Point& x3, const Point& x4, Point& p );
	bool PointBetweenTwoPoints(const Point& x1, const Point x2, const Point& testPoint);


#ifdef DEBUG_POINTREPULSION
	bool VeryVerySmallDistancePointToPlane(Normal& aNormal, Point& aPoint, Point& aTestPoint );

#endif

	int iNumberOfTriangles;
	int iNumberOfVertices;
	int *iVertexIndex;
	Point *iVertices;
	int iNumberOfSamplePoints;
	float iRepulsiveRadius;
	float iTotalAreaOfSurface;

	// scaling value to scale the samplepoint force vector
	float iK;

	vector<Point* > iSamplePoints;
	vector<SamplePointContainer *> iSamplePointContainer;
	vector<TriangleUseSet> iTriangles;
};


#endif
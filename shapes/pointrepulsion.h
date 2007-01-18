#ifndef PBRT_POINTREPULSION_H
#define PBRT_POINTREPULSION_H

#include "geometry.h"
#include "triangleuseset.h"
#include "samplepointcontainer.h"
#include "shape.h"
#include <vector>

using namespace std;

class Triangle;

class PointRepulsion
{
public:
	PointRepulsion(int aNumberOfTriangles, int aNumberOfVertices, int* aVertexIndex, Point* aVertices,vector<Reference<Shape> >& aTriangleList );
	~PointRepulsion();
	int SetupSamplePoints( float aMeanFreePath );
	void ComputeRepulsiveForces( const float& aForceScale );
	void ComputeNewPositions();
	void FillUniformSamplePointStructure( vector<UniformPoint>& container );
	float GetTotalSurfaceArea();

private:
	void SetupTriangleUseSets(vector<Reference<Shape> >& aTriangleList);
	void CalculateTransformationMatrices(TriangleUseSet &aTriangle, TriangleUseSet& aNeighborTriangle,
										 Point& P0, Point& P1, Reference<Matrix4x4>& aArbitraryRotation,
										 Reference<Matrix4x4>& aArbitraryRotationInv, float& aRotationAngle);


	void ComputePartialAreaSum(float* aArea);
	void ComputeSamplePointPosition(float* aAreas, float s, float t);

	void MapSamplePointsToPlane(TriangleUseSet& aCurrentTriangle, TriangleUseSet& aMainTriangle, Reference<Matrix4x4> aEdgeRotationMatrix, vector<bool>& aTriangleMapped, const float& aForceScale);
	void PointRepulsion::ComputeRepulsiveForces(TriangleUseSet& aCurrentTriangle,
												TriangleUseSet& aMainTriangle,
												Reference<Matrix4x4> aEdgeRotationMatrix,
												const float& aForceScale);
	bool PointInsideTriangle(const Point& aPoint, const Point& aStartPoint, const TriangleUseSet& aTriangle, Neighbor** aNeighbor, Point& aP0, Point& aP1, Point& aEdgePoint);
	Point LineLineIntersection(const Point& x1,const Point& x2,const Point& x3,const Point& x4 );

	void ComputeNewPositions(Point& aPPrime, TriangleUseSet* aUseSet, Point* p, SamplePointContainer* container);
	bool LineSegLineSegIntersection(const Point& x1, const Point& x2, const Point& x3, const Point& x4, Point& p );
	bool PointBetweenTwoPoints(const Point& x1, const Point x2, const Point& testPoint);
	static Triangle* Cast(Reference<Shape>& aTriangle);

#ifdef DEBUG_POINTREPULSION
	bool VeryVerySmallDistancePointToPlane(Normal& aNormal, Point& aPoint, Point& aTestPoint );
#endif

	int iNumberOfTriangles;
	int iNumberOfVertices;
	int *iVertexIndex;
	Point *iVertices;
	float iRepulsiveRadius;
	float iTotalAreaOfSurface;

	vector<Point* > iSamplePoints;
	vector<SamplePointContainer *> iSamplePointContainer;
	vector<TriangleUseSet> iTriangles;
};


#endif
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
	void ComputeRepulsiveForces( const float& aForceScale, const string& aProcessString);
	void ComputeNewPositions( const string& aProcessString );
	void FillUniformSamplePointStructure( vector<UniformPoint>& container );
	float GetTotalSurfaceArea();
	static float CreateTriangleUseSets(vector<Reference<Shape> >& aTriangleList, vector<UniformPoint>& aSamplePoints, 
		Point* iVertices);

private:
	void SetupTriangleUseSets(vector<Reference<Shape> >& aTriangleList);
	void CreateTriangleUseSets(vector<Reference<Shape> >& aTriangleList);
	void LinkTriangleUseSets();
	//void CalculateTransformationMatrices(TriangleUseSet &aTriangle, TriangleUseSet& aNeighborTriangle,
	//									 Point& P0, Point& P1, Reference<Matrix4x4>& aArbitraryRotation,
	//									 Reference<Matrix4x4>& aArbitraryRotationInv, float& aRotationAngle);
	void CalculateTransformationMatrices( TriangleUseSet &aTriangle, TriangleUseSet& aNeighborTriangle,
										  Point& P0, Point& P1,
										  Reference<Matrix4x4>& aArbitraryRotation,
										  Reference<Matrix4x4>& aArbitraryRotationInv,
										  Reference<Matrix4x4>& aTranslateToAxis,
										  Reference<Matrix4x4>& aTranslateToAxisInv );

	bool ComputeNewPositions(Point& aPPrime, Point* p, SamplePointContainer* container);
	void ComputeSamplePointPosition(vector<float>& aPartialAreaSums, float s, float t);
	void MapSamplePointsToPlane(list< pair<TriangleUseSet*, Reference<Matrix4x4> > >& aTriangleQueue, TriangleUseSet& aMainTriangle,
								vector<bool>& aTriangleMapped, const float& aForceScale);
	void ComputeRepulsiveForces(TriangleUseSet& aCurrentTriangle,
								TriangleUseSet& aMainTriangle,
								Reference<Matrix4x4> aEdgeRotationMatrix,
								const float& aForceScale);

	bool PointInsideTriangle(const Point& aPoint, const Point& aStartPoint, const TriangleUseSet& aTriangle, TriangleEdge** aNeighbor, Point& aP0, Point& aP1, Point& aEdgePoint);
	Point LineLineIntersection(const Point& x1,const Point& x2,const Point& x3,const Point& x4 );

	bool LineSegLineSegIntersection(const Point& x1, const Point& x2, const Point& x3, const Point& x4, Point& p );
	bool PointBetweenTwoPoints(const Point& x1, const Point x2, const Point& testPoint);
	int SearchForTriangle(const float aS, const vector<float>& aAreas, int lower, int upper );
	void PushPointToPlane(const Normal& aNormal, const Point& aPoint, Point& aTestPoint );
	bool PointCloseToPlane( const Normal& aNormal, const Point& aPoint, Point& aTestPoint );
	float ComputeSmallestDistanceBetweenTriangles(TriangleUseSet& aEvaluatedTriangle, TriangleUseSet& aMainTriangle, Transform& aTransform);

	int iNumberOfTriangles;
	int iNumberOfVertices;
	int *iVertexIndex;
	Point *iVertices;
	float iRepulsiveRadius;
	float iTotalAreaOfSurface;

	vector<Point* > iSamplePoints;
	vector<SamplePointContainer *> iSamplePointContainer;
	vector<TriangleUseSet> iTriangles;

	vector<float> iPartialAreaSums;
};


#endif
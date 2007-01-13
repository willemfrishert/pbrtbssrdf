#ifndef PBRT_TRIANGLE_H
#define PBRT_TRIANGLE_H

#include "geometry.h"
#include "transform.h"
#include <list>
using namespace std;

class TriangleUseSet;
class SamplePointContainer;

struct Neighbor
{
	Neighbor(TriangleUseSet* aEdgeNeighbor,
			 Point* aEdgeNeighborsVertex[2],
			 Reference<Matrix4x4> aArbitraryRotation,
			 Reference<Matrix4x4> aArbitraryRotationInv,
			 float aRotationAngle)
	: iEdgeNeighbor(aEdgeNeighbor)
	, iArbitraryRotation(aArbitraryRotation)
	, iArbitraryRotationInv(aArbitraryRotationInv)
	, iRotationAngle(aRotationAngle)
	{
		iP0 = aEdgeNeighborsVertex[0];
		iP1 = aEdgeNeighborsVertex[1];

		iTranslateToAxis = new Matrix4x4(1, 0, 0, iP0->x,
										 0, 1, 0, iP0->y,
										 0, 0, 1, iP0->z,
										 0, 0, 0,      1);

		iTranslateToAxisInv = new Matrix4x4(1, 0, 0, -iP0->x,
											0, 1, 0, -iP0->y,
											0, 0, 1, -iP0->z,
											0, 0, 0,       1);
	};

	// pointer to the triangle next door
	TriangleUseSet* iEdgeNeighbor;
	
	// Points making up the shared edge between two triangles
	Point* iP0;
	Point* iP1;
	
	// Matrices to rotate around an arbiraty axis
	Reference<Matrix4x4> iTranslateToAxis;
	Reference<Matrix4x4> iTranslateToAxisInv;

	Reference<Matrix4x4> iArbitraryRotation;
	Reference<Matrix4x4> iArbitraryRotationInv;

	// angle between two triangles
	float iRotationAngle;
};

class TriangleUseSet
{
	//METHODS
public:
	TriangleUseSet(Point* aV1, Point* aV2, Point* aV3);
	~TriangleUseSet();

	void AddSamplePoint(SamplePointContainer* aSamplePoint);
	void DeleteSamplePoint(SamplePointContainer* aSamplePoint);
	float ComputeArea();
	Normal ComputeNormal();
	void AddEdgeNeightbor(Neighbor* aEdgeNeighbor);
	vector<Neighbor*> GetEdgeNeighbors();
	inline int GetTriangleId();
	Vector GetCentroid();

	Point* iVertices[3];
	list<SamplePointContainer* > iSamplePoints;

private:
	// VARIABLES
	static int triangleCounter;

	int triangleId;
	vector<Neighbor* > iEdgeNeighbors;
	Vector iCentroid;
	Normal iUnNormalizedNormal;
};


inline 
int TriangleUseSet::GetTriangleId()
{
	return triangleId;
}

#endif
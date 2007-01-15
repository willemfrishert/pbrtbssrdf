#include "triangleuseset.h"
#include "samplepointcontainer.h"

int TriangleUseSet::triangleCounter=0;

TriangleUseSet::TriangleUseSet( Point* aP1, Point* aP2, Point* aP3 )
{
	triangleId = triangleCounter;
	triangleCounter++;

	iVertices[0] = aP1;
	iVertices[1] = aP2;
	iVertices[2] = aP3;

	iEdgeNeighbors[0] = NULL;
	iEdgeNeighbors[1] = NULL;
	iEdgeNeighbors[2] = NULL;

	iCentroid.x = (iVertices[0]->x+iVertices[1]->x+iVertices[2]->x)*0.33f;
	iCentroid.y = (iVertices[0]->y+iVertices[1]->y+iVertices[2]->y)*0.33f;
	iCentroid.z = (iVertices[0]->z+iVertices[1]->z+iVertices[2]->z)*0.33f;

	iUnNormalizedNormal = Normal( Cross( *aP2-*aP1, *aP3-*aP1 ) );
}

TriangleUseSet::~TriangleUseSet()
{
	iVertices[0] = NULL;
	iVertices[1] = NULL;
	iVertices[2] = NULL;

	iEdgeNeighbors[0] = NULL;
	iEdgeNeighbors[1] = NULL;
	iEdgeNeighbors[2] = NULL;
}


void TriangleUseSet::GetEdgeNeighbors( vector<Neighbor*>& aNeighbor ) const
{
	for (int i=0;i<3;i++)
	{
		if ( iEdgeNeighbors[i] != NULL )
		{
			aNeighbor.push_back(iEdgeNeighbors[i]);
		}
	}
}

void TriangleUseSet::GetAllEdgeNeighbors( vector<Neighbor*>& aNeighbor ) const
{
	for (int i=0;i<3;i++)
	{
		aNeighbor.push_back(iEdgeNeighbors[i]);
	}
}

void TriangleUseSet::AddEdgeNeightbor( Neighbor* aEdgeNeighbor, u_int aPosition )
{
	assert(aPosition < 3);
	iEdgeNeighbors[aPosition] = aEdgeNeighbor;
}

void TriangleUseSet::AddSamplePoint( SamplePointContainer* aSamplePoint )
{
	this->iSamplePoints.push_back( aSamplePoint );
}

void TriangleUseSet::DeleteSamplePoint( SamplePointContainer* aSamplePoint )
{
	this->iSamplePoints.remove( aSamplePoint );
}

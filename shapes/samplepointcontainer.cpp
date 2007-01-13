#include "samplepointcontainer.h"

SamplePointContainer::SamplePointContainer( Point* aPoint, TriangleUseSet* aTriangle )
: iPoint(aPoint)
, iTriangle( aTriangle )
{
	iTriangle->AddSamplePoint( this );
}

SamplePointContainer::~SamplePointContainer()
{
	iPoint = NULL;
	iTriangle = NULL;
}

void SamplePointContainer::AddForceVector( Vector& aForceVector )
{
	iS += aForceVector;
}

void SamplePointContainer::ResetForceVector()
{
	iS = Vector();
}

Vector SamplePointContainer::GetForceVector()
{
	return iS;
}
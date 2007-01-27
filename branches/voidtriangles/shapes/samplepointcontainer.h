#ifndef PBRT_SAMPLEPOINT_H
#define PBRT_SAMPLEPOINT_H

#include "geometry.h"
#include "triangleuseset.h"

class SamplePointContainer
{
public:
	SamplePointContainer(Point* aPoint, TriangleUseSet* aTriangle);
	~SamplePointContainer();
	Point* GetSamplePoint();
	TriangleUseSet* GetTriangle();
	void SetTriangle( TriangleUseSet* aTriangle );
	void AddForceVector( Vector& aForceVector);
	void ResetForceVector();
	Vector GetForceVector();
public:
protected:
private:
	Point* iPoint;
	TriangleUseSet* iTriangle;
	Vector iS;
};


inline
Point* SamplePointContainer::GetSamplePoint()
{
	return iPoint;
}


inline
TriangleUseSet* SamplePointContainer::GetTriangle()
{
	return iTriangle;
}

inline
void SamplePointContainer::SetTriangle( TriangleUseSet* aTriangle )
{
	iTriangle->DeleteSamplePoint( this );

	iTriangle = aTriangle;



	iTriangle->AddSamplePoint( this );
}

#endif
#include<iostream>
#include<cmath>
#include <limits>
#include<cstdlib>

using namespace std;

#define RealType double

struct SInteriorSpace
{
	RealType fX;///< ksi
	RealType fY;///< theta
	RealType fZ;///< dtheta / dksi

	void Add(const RealType fA, const RealType fB, const RealType fC);

};

ostream &operator<<( ostream &output, const SInteriorSpace &Point )
{
	output << scientific;
	output << Point.fX<<" "<<Point.fY<<" "<<Point.fZ;
	return output;
}

void SInteriorSpace::Add(const RealType fA, const RealType fB, const RealType fC)
{
	fX+=fA;
	fY+=fB;
	fZ+=fC;
}



class CPolytropeShooting
{
protected:
	RealType m_YPrime(const SInteriorSpace Point) const;
	RealType m_ZPrime(const SInteriorSpace Point) const;

	void m_RKMethod();
	SInteriorSpace m_ApproximateSolution(const RealType fX);

	SInteriorSpace m_PropagateRKMethod(const SInteriorSpace PreviousPoint, const RealType fDeltaX);
	RealType m_ApproximateSolutionY(const RealType fX);
	RealType m_ApproximateSolutionZ(const RealType fX);


public:
	RealType fPolytropicIndex;

	void ComputeInterior();

	CPolytropeShooting();
	CPolytropeShooting(const RealType fIndex);
	~CPolytropeShooting(){};

};

CPolytropeShooting::CPolytropeShooting()
{
	CPolytropeShooting(1.);
}

CPolytropeShooting::CPolytropeShooting(const RealType fIndex)
{
	fPolytropicIndex=fIndex;
}

SInteriorSpace CPolytropeShooting::m_ApproximateSolution(const RealType fX)
{
	SInteriorSpace Result;
	Result.fX=fX;
	Result.fY=m_ApproximateSolutionY(fX);
	Result.fZ=m_ApproximateSolutionY(fX);
	return Result;
}

RealType CPolytropeShooting::m_ApproximateSolutionY(const RealType fX)
{
	RealType fXX (fX*fX);
	return ( 1. + fXX *(-1./6. + fXX*( fPolytropicIndex/120. - fXX*fPolytropicIndex*(8*fPolytropicIndex-5)/15120. ))  );
}

RealType CPolytropeShooting::m_ApproximateSolutionZ(const RealType fX)
{
	RealType fXX (fX*fX);
	return ( fX *(-1./3. + fXX*( fPolytropicIndex/300. - fXX*fPolytropicIndex*(8*fPolytropicIndex-5)/2520. ))  );
}

RealType CPolytropeShooting::m_YPrime(const SInteriorSpace Point) const
{
	return (Point.fZ);
}

RealType CPolytropeShooting::m_ZPrime(const SInteriorSpace Point) const
{
	return (-pow(Point.fY,fPolytropicIndex) - Point.fZ*2./Point.fX);
}

SInteriorSpace CPolytropeShooting::m_PropagateRKMethod(const SInteriorSpace PreviousPoint, const RealType fDeltaX)
{
	SInteriorSpace WorkPoint;

	WorkPoint = PreviousPoint;

	RealType fK1( fDeltaX * m_YPrime( WorkPoint ) );
	RealType fL1( fDeltaX * m_ZPrime( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaX/2., fK1/2., fL1/2. );

	RealType fK2( fDeltaX * m_YPrime( WorkPoint ) );
	RealType fL2( fDeltaX * m_ZPrime( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaX/2., fK2/2., fL2/2. );


	RealType fK3( fDeltaX * m_YPrime( WorkPoint ) );
	RealType fL3( fDeltaX * m_ZPrime( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaX, fK3, fL3 );

	RealType fK4( fDeltaX * m_YPrime( WorkPoint ) );
	RealType fL4( fDeltaX * m_ZPrime( WorkPoint ) );


	SInteriorSpace Result;

	Result.fX = PreviousPoint.fX + fDeltaX;
	Result.fY = PreviousPoint.fY + fK1/6. + fK2/3. + fK3/3. +fK4/6.;
	Result.fZ = PreviousPoint.fZ + fL1/6. + fL2/3. + fL3/3. +fL4/6.;

	return Result;
}


void CPolytropeShooting::ComputeInterior()
{
	//const RealType fAlmostZero ( 10.*std::numeric_limits<RealType>::epsilon() );

	const RealType fAlmostZero ( 1e-5);
	const RealType fStep(1.e-4);
	SInteriorSpace PointOfInterest;
	PointOfInterest = m_ApproximateSolution(fAlmostZero);
	cout<<PointOfInterest<<endl;
	while (PointOfInterest.fY>fAlmostZero)
	{
		PointOfInterest = m_PropagateRKMethod(PointOfInterest,fStep);
		cout<<PointOfInterest<<endl;
	}
}


int main(int argc, char** argv)
{
	if(argc != 2)
	{
		cout<<"Wrong arguments"<<endl;
		cout<<"./exe fPolitropicIndex"<<endl;
		return 1;
	}

	CPolytropeShooting Poly(atof(argv[1]));
	Poly.ComputeInterior();


	return 0;
}

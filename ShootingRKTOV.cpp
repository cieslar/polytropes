#include<iostream>
#include<cmath>
#include <limits>
#include<cstdlib>

using namespace std;

#define RealType double

struct SInteriorSpaceTOV
{
	RealType fR;
	RealType fM;
	RealType fMProper;
	RealType fP;
	RealType fPhi;

	RealType fDMDR;
	RealType fDMProperDR;
	RealType fDPDR;
	RealType fDPhiDR;

	void Zero();
};


void SInteriorSpaceTOV::Zero()
{
	fR=0.;
	fM=0.;
	fMProper=0.;
	fP=0.;
	fPhi=0.;

	fDMDR=0.;
	fDMProperDR=0.;
	fDPDR=0.;
	fDPhiDR=0.;
}


SInteriorSpaceTOV &operator+=(SInteriorSpaceTOV &Left, const SInteriorSpaceTOV &Right)
{
	Left.fR+=Right.fR;
	Left.fM+=Right.fM;
	Left.fMProper+=Right.fMProper;
	Left.fP+=Right.fP;
	Left.fPhi+=Right.fPhi;

	Left.fDMDR+=Right.fDMDR;
	Left.fDMProperDR+=Right.fDMProperDR;
	Left.fDPDR+=Right.fDPDR;
	Left.fDPhiDR+=Right.fDPhiDR;

}

SInteriorSpaceTOV operator+(const SInteriorSpaceTOV &Left, const SInteriorSpaceTOV &Right)
{
	SInteriorSpaceTOV Result(Left);
	Result+=Right;
	return Result;
}


SInteriorSpaceTOV operator/(const SInteriorSpaceTOV Point, const RealType fNumber)
{
	SInteriorSpaceTOV Result(Point);
	Result.fR/=fNumber;
	Result.fM/=fNumber;
	Result.fMProper/=fNumber;
	Result.fP/=fNumber;
	Result.fPhi/=fNumber;

	Result.fDMDR/=fNumber;
	Result.fDMProperDR/=fNumber;
	Result.fDPDR/=fNumber;
	Result.fDPhiDR/=fNumber;
	return Result;
}

ostream &operator<<( ostream &output, const SInteriorSpaceTOV &Point )
{
	output << scientific;
	output << Point.fR<<" "<<Point.fM<<" "<<Point.fMProper<<" "<<Point.fP<<" "<<Point.fPhi<<" "<<Point.fDMDR<<" "<<Point.fDMProperDR<<" "<<Point.fDPDR<<" "<<Point.fDPhiDR;
	return output;
}

//how to impreve the name of this function???
RealType MultiDimRelativeDifference(const SInteriorSpaceTOV A, const SInteriorSpaceTOV B)
{
	return sqrt(pow((B.fM-A.fM)/B.fM, 2.) 
                  + pow((B.fP-A.fP)/B.fP, 2.) 
                  + pow((B.fPhi-A.fPhi)/B.fPhi, 2.) 
                  + pow((B.fDMDR-A.fDMDR)/B.fDMDR, 2.) 
                  + pow((B.fDPDR-A.fDPDR)/B.fDPDR, 2.) 
                  + pow((B.fDPhiDR-A.fDPhiDR)/B.fDPhiDR, 2.));
}

bool isnan(const SInteriorSpaceTOV A)
{
	return ( isnan(A.fR) || isnan(A.fP) || isnan(A.fM) || isnan(A.fMProper) || isnan(A.fPhi) || isnan(A.fDMDR) || isnan(A.fDMProperDR) || isnan(A.fDPhiDR) || isnan(A.fDPDR) );
}

class CPolytropeShooting
{
protected:
	RealType m_MDerivative(const SInteriorSpaceTOV Point) const;
	RealType m_MProperDerivative(const SInteriorSpaceTOV Point) const;
	RealType m_PDerivative(const SInteriorSpaceTOV Point) const;
	RealType m_PhiDerivative(const SInteriorSpaceTOV Point) const;

	SInteriorSpaceTOV m_ApproximateSolutionInZero(const RealType fEpsilon) const; 
	SInteriorSpaceTOV m_PropagateRKMethod(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const;

public:
	RealType fPolytropicGamma;
	RealType fRhoCentre;

	void ComputeInterior();
	void ComputeInteriorAdaptive();

	RealType ComputePhiConstant(const SInteriorSpaceTOV Point) const;

	CPolytropeShooting(const RealType fGamma, const RealType fRhoC);
	~CPolytropeShooting(){};

};


CPolytropeShooting::CPolytropeShooting(const RealType fGamma, const RealType fRhoC)
{
	fPolytropicGamma=fGamma;
	fRhoCentre=fRhoC;
}

SInteriorSpaceTOV CPolytropeShooting::m_ApproximateSolutionInZero(const RealType fEpsilon) const
{
	SInteriorSpaceTOV Result;

	Result.fR=fEpsilon;
	Result.fM=0.;
	Result.fP=0.;
	Result.fPhi=1.;

	Result.fDMDR=m_MDerivative(Result);
	Result.fDMProperDR=m_MProperDerivative(Result);
	Result.fDPDR=m_PDerivative(Result);
	Result.fDPhiDR=m_PhiDerivative(Result);

	return Result;
}


RealType CPolytropeShooting::m_MDerivative(const SInteriorSpaceTOV Point) const
{
	return 0.;
}

RealType CPolytropeShooting::m_MProperDerivative(const SInteriorSpaceTOV Point) const
{
	return 0.;
}

RealType CPolytropeShooting::m_PDerivative(const SInteriorSpaceTOV Point) const
{
	return 0.;
}

RealType CPolytropeShooting::m_PhiDerivative(const SInteriorSpaceTOV Point) const
{
	return 0.;
}

SInteriorSpaceTOV CPolytropeShooting::m_PropagateRKMethod(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const
{
	SInteriorSpaceTOV WorkPoint;
	SInteriorSpaceTOV DWorkPoint1;
	SInteriorSpaceTOV DWorkPoint2;
	SInteriorSpaceTOV DWorkPoint3;
	SInteriorSpaceTOV DWorkPoint4;

	WorkPoint = PreviousPoint;

	DWorkPoint1.Zero();
	DWorkPoint1.fR       = fDeltaR;
	DWorkPoint1.fM       = fDeltaR * m_MDerivative( WorkPoint );
	DWorkPoint1.fMProper = fDeltaR * m_MProperDerivative( WorkPoint );
	DWorkPoint1.fP       = fDeltaR * m_PDerivative( WorkPoint );
	DWorkPoint1.fPhi     = fDeltaR * m_PhiDerivative( WorkPoint );

	WorkPoint = PreviousPoint;
	WorkPoint += DWorkPoint1/2.; 
//	WorkPoint.Add(fDeltaR/2., fM1/2., fP1/2., fPhi1/2. );

	DWorkPoint2.Zero();
	DWorkPoint2.fR       = fDeltaR;
	DWorkPoint2.fM       = fDeltaR * m_MDerivative( WorkPoint );
	DWorkPoint2.fMProper = fDeltaR * m_MProperDerivative( WorkPoint );
	DWorkPoint2.fP       = fDeltaR * m_PDerivative( WorkPoint );
	DWorkPoint2.fPhi     = fDeltaR * m_PhiDerivative( WorkPoint );

	WorkPoint = PreviousPoint;
	WorkPoint += DWorkPoint2/2.; 

	DWorkPoint3.Zero();
	DWorkPoint3.fR       = fDeltaR;
	DWorkPoint3.fM       = fDeltaR * m_MDerivative( WorkPoint );
	DWorkPoint3.fMProper = fDeltaR * m_MProperDerivative( WorkPoint );
	DWorkPoint3.fP       = fDeltaR * m_PDerivative( WorkPoint );
	DWorkPoint3.fPhi     = fDeltaR * m_PhiDerivative( WorkPoint );

	WorkPoint = PreviousPoint;
	WorkPoint += DWorkPoint3; 

	DWorkPoint4.Zero();
	DWorkPoint4.fR       = fDeltaR;
	DWorkPoint4.fM       = fDeltaR * m_MDerivative( WorkPoint );
	DWorkPoint4.fMProper = fDeltaR * m_MProperDerivative( WorkPoint );
	DWorkPoint4.fP       = fDeltaR * m_PDerivative( WorkPoint );
	DWorkPoint4.fPhi     = fDeltaR * m_PhiDerivative( WorkPoint );

	SInteriorSpaceTOV Result;

	Result = PreviousPoint + DWorkPoint1/6. + DWorkPoint2/3. + DWorkPoint3/3. + DWorkPoint4/6.;
	//Correcting fR:
	Result.fR = PreviousPoint.fR + fDeltaR;

	Result.fDMDR =  m_MDerivative( Result );
	Result.fDMProperDR =  m_MProperDerivative( Result );
	Result.fDPDR =  m_PDerivative( Result );
	Result.fDPhiDR = m_PhiDerivative( Result );

	return Result;
}
	

void CPolytropeShooting::ComputeInterior()
{
	//const RealType fAlmostZero ( 10.*std::numeric_limits<RealType>::epsilon() );

	const RealType fAlmostZero ( 1e-5);
	const RealType fStep(1.e-4);
	SInteriorSpaceTOV PointOfInterest;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<PointOfInterest<<endl;
	while (PointOfInterest.fP>fAlmostZero)
	{
		PointOfInterest = m_PropagateRKMethod(PointOfInterest,fStep);
		cout<<PointOfInterest<<endl;
	}
}

void CPolytropeShooting::ComputeInteriorAdaptive()
{
	const RealType fAlmostZero ( 1e-10);
	const RealType fEnough(1.e-12);
	RealType fStep(1.e-12);

	SInteriorSpaceTOV PointOfInterest, Small, Big;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<PointOfInterest<<" "<<fStep<<endl;
	while (PointOfInterest.fP>fAlmostZero)
	{
		do
		{
			fStep /= 2.;
			Small = m_PropagateRKMethod( m_PropagateRKMethod(PointOfInterest,fStep/2.) , fStep/2.);
			Big   = m_PropagateRKMethod( PointOfInterest,fStep);
			//cout<<" "<<Small<<" "<<fStep/2.<<endl;
			//cout<<" "<<Big<<" "<<fStep<<endl;
			//cout<<" "<<MultiDimRelativeDifference(Small,Big) <<" "<< fEnough<<endl;
		}
		while(  (MultiDimRelativeDifference(Small,Big) > fEnough) || isnan(Big)  ||  Big.fP < 0. );
		PointOfInterest=Small;
		cout<<PointOfInterest<<" "<<fStep<<endl;
		fStep *= 4.;
	}
}
#if 0
void CPolytropeShooting::ComputeInteriorAdaptiveBisection()
{
	const RealType fAlmostZero ( 1e-5);
	const RealType fEnough(1.e-7);
	RealType fStep(1.e-12);

	SInteriorSpaceTOV PointOfInterest, Small, Big;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
	while (PointOfInterest.fPhi>fAlmostZero)
	{
		do
		{
			fStep /= 2.;
			Small = m_PropagateRKMethod( m_PropagateRKMethod(PointOfInterest,fStep/2.) , fStep/2.);
			Big   = m_PropagateRKMethod( PointOfInterest,fStep);
			//cout<<" "<<Small<<" "<<m_MassTilda(Small)<<" "<<fStep/2.<<endl;
			//cout<<" "<<Big<<" "<<m_MassTilda(Big)<<" "<<fStep<<endl;
			//cout<<" "<<MultiDimRelativeDifference(Small,Big) <<" "<< fEnough<<endl;
		}
		while( (MultiDimRelativeDifference(Small,Big) > fEnough) || isnan(Big) || Big.fPhi < 0. );
		PointOfInterest=Small;
		cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
		fStep *= 4.;
	}
	fStep /= 8.;
	//Find the such fR so the fPhi<0	
	SInteriorSpaceTOV Left, Right;
	Left=PointOfInterest; //>0
	//cout<<Left<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
	while(PointOfInterest.fPhi>0.)
	{
		PointOfInterest=m_PropagateRKMethod( PointOfInterest,fStep);
	}
	Right = PointOfInterest;
	//cout<<Right<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;

	const RealType fZero(1e-15);
	//Do the bisection
	while(Left.fPhi>fZero)
	{
		fStep=(Right.fR-Left.fR)/2.;
		PointOfInterest=m_PropagateRKMethod( Left,fStep);
		if(PointOfInterest.fPhi>0.) Left=PointOfInterest;
		else Right=PointOfInterest;

	}
	cout<<Left<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
}


#endif
int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout<<"Wrong arguments"<<endl;
		cout<<"./PolyTOVShooting fPolytropicGamma fRhoCentre"<<endl;
		return 1;
	}

	CPolytropeShooting Poly(atof(argv[1]),atof(argv[2]));
	Poly.ComputeInteriorAdaptive();
	//Poly.ComputeInterior();


	return 0;
}

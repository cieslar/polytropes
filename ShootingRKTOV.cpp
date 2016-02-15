#include<iostream>
#include<cmath>
#include<limits>
#include<cstdlib>
#include<Debug.hpp>

using namespace std;

#define RealType double

//Can we estimate the needed precision of constants?
const RealType gfG(6.67408e-8);///<Gravitational constant in [cm^3 g^-1 s^-2]
const RealType gfC(2.99792458e10);///<Speed of light in [cm s-1] 
const RealType gfBaryonMass(1.66e-24);///<Baryon mass in [g]
//const RealType gfBaryonMass(1.);///<Baryon mass in [g]

RealType DeKappify(const RealType fLoreneKappa, const RealType fGamma)
{
	const RealType cfRhoNuc(1.66e14);///<[g cm^-3]
	//const RealType cfNNuc(1.e38);///<[cm^-3]
	//Translates Lorene units to cgs
//	return gfBaryonMass*fLoreneKappa*cfRhoNuc*gfC*gfC/pow(cfNNuc,fGamma);

	return fLoreneKappa*gfC*gfC/pow(cfRhoNuc,fGamma-1.);

	//Changes the equations from mass density to barion density to be consistent with Lorene.
	//return fLoreneKappa*cfRhoNuc*gfC*gfC/pow(cfNNuc,fGamma);
}

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

SInteriorSpaceTOV operator*(const SInteriorSpaceTOV Point, const RealType fNumber)
{
	SInteriorSpaceTOV Result(Point);

	Result.fR *= fNumber;
	Result.fM *= fNumber;
	Result.fMProper *= fNumber;
	Result.fP *= fNumber;
	Result.fPhi *= fNumber;

	Result.fDMDR *= fNumber;
	Result.fDMProperDR *= fNumber;
	Result.fDPDR *= fNumber;
	Result.fDPhiDR *= fNumber;

	return Result;
}

SInteriorSpaceTOV operator*(const RealType fNumber, const SInteriorSpaceTOV Point)
{
	return Point*fNumber;
}
SInteriorSpaceTOV operator/(const SInteriorSpaceTOV Point, const RealType fNumber)
{
	SInteriorSpaceTOV Result(Point);

	Result.fR /= fNumber;
	Result.fM /= fNumber;
	Result.fMProper /= fNumber;
	Result.fP /= fNumber;
	Result.fPhi /= fNumber;

	Result.fDMDR /= fNumber;
	Result.fDMProperDR /= fNumber;
	Result.fDPDR /= fNumber;
	Result.fDPhiDR /= fNumber;

	return Result;
}
#include <sstream>
std::string DescInteriorSpaceTOV()
{
	ostringstream output;
//	output <<"#R[km] M[MSun] MProper[MSun] P[dyn*cm^-2] Phi[holerawie] DMDR DMProperDR DPDR DPhiDR"<<endl;
	output <<"#R[cm] M[g] MProper[g] P[dyn*cm^-2] Phi[holerawie] DMDR DMProperDR DPDR DPhiDR"<<endl;
	return output.str();
	
}
ostream &operator<<( ostream &output, const SInteriorSpaceTOV &Point )
{
	output << scientific;
	const RealType fMSun(1.9891e33);
	//output << Point.fR/1e5<<" "<<Point.fM/fMSun<<" "<<Point.fMProper/fMSun<<" "<<Point.fP<<" "<<Point.fPhi<<" "<<Point.fDMDR<<" "<<Point.fDMProperDR<<" "<<Point.fDPDR<<" "<<Point.fDPhiDR;
	output << Point.fR<<" "<<Point.fM<<" "<<Point.fMProper<<" "<<Point.fP<<" "<<Point.fPhi<<" "<<Point.fDMDR<<" "<<Point.fDMProperDR<<" "<<Point.fDPDR<<" "<<Point.fDPhiDR;
	return output;
}

//how to improve the name of this function???
//Computes the realtive multidim cartesian difference for values changed by integration method
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
	SInteriorSpaceTOV m_Propagate(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const;
	SInteriorSpaceTOV m_PropagateRK2Method(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const;

	RealType m_fPolytropicK2Gamma;
public:
	RealType fPolytropicGamma;
	RealType fPolytropicK;
	RealType fRhoCentre;

	void ComputeInterior();
	void ComputeInteriorAdaptive();

	RealType ComputePhiConstant(const SInteriorSpaceTOV Point) const;

	CPolytropeShooting(const RealType fGamma, const RealType fK, const RealType fRhoC);
	~CPolytropeShooting(){};

};




CPolytropeShooting::CPolytropeShooting(const RealType fGamma, const RealType fK, const RealType fRhoC)
{
	fPolytropicGamma=fGamma;
	fPolytropicK=fK;
	fRhoCentre=fRhoC;

	m_fPolytropicK2Gamma=pow(fPolytropicK,fPolytropicGamma);
}

SInteriorSpaceTOV CPolytropeShooting::m_ApproximateSolutionInZero(const RealType fEpsilon) const
{
	SInteriorSpaceTOV Result;

	Result.fR=fEpsilon;
	Result.fM=fRhoCentre*fEpsilon*fEpsilon*fEpsilon*(4./3.)*M_PI;
	Result.fP=fPolytropicK*pow(fRhoCentre,fPolytropicGamma);
	Result.fPhi=1.;

	Result.fDMDR=m_MDerivative(Result);
	Result.fDMProperDR=m_MProperDerivative(Result);

	//Result.fDPDR=0.;
	//Result.fDPhiDR=0.;
	Result.fDPDR=m_PDerivative(Result);
	Result.fDPhiDR=m_PhiDerivative(Result);
	
	return Result;
}


RealType CPolytropeShooting::m_MDerivative(const SInteriorSpaceTOV Point) const
{
	return ( 4. * M_PI * m_fPolytropicK2Gamma * Point.fR * Point.fR ) /
               pow(Point.fP, fPolytropicGamma);
}

//TODO check if the metric is ok (after Shapiro)
RealType CPolytropeShooting::m_MProperDerivative(const SInteriorSpaceTOV Point) const
{
	return m_MDerivative(Point) * pow( (1.- (2. * gfG * Point.fM ) / ( Point.fR * gfC * gfC) ), -0.5 );
}

RealType CPolytropeShooting::m_PDerivative(const SInteriorSpaceTOV Point) const
{
	return ( -gfG * Point.fM * m_fPolytropicK2Gamma / (pow(Point.fP,fPolytropicGamma) * (Point.fR*Point.fR)) ) *
               ( 1. + (pow(Point.fP,fPolytropicGamma+1.) / (gfC*gfC*m_fPolytropicK2Gamma) ) ) * 
               ( 1. + (4. * M_PI * Point.fR * Point.fR * Point.fR * Point.fP / ( Point.fM * gfC * gfC ) ) ) / 
               ( 1. - ( ( 2. * gfG * Point.fM ) / ( Point.fR * gfC * gfC ) ) );
}


//The metric
RealType CPolytropeShooting::m_PhiDerivative(const SInteriorSpaceTOV Point) const
{
	return (-2./(Point.fP + gfC*gfC * m_fPolytropicK2Gamma/pow(Point.fP,fPolytropicGamma))) * m_PDerivative(Point);
}


SInteriorSpaceTOV CPolytropeShooting::m_Propagate(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const
{
	SInteriorSpaceTOV DResult;

	DResult.Zero();

	DResult.fR       = fDeltaR;
	DResult.fM       = fDeltaR * m_MDerivative( PreviousPoint );
	DResult.fMProper = fDeltaR * m_MProperDerivative( PreviousPoint );
	DResult.fP       = fDeltaR * m_PDerivative( PreviousPoint );
	DResult.fPhi     = fDeltaR * m_PhiDerivative( PreviousPoint );

	return PreviousPoint+DResult;
}

SInteriorSpaceTOV CPolytropeShooting::m_PropagateRK2Method(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaR) const
{
	SInteriorSpaceTOV K1, K2;

	K1.Zero();
	K1.fR       = fDeltaR;
	K1.fM       = fDeltaR * m_MDerivative( PreviousPoint );
	K1.fMProper = fDeltaR * m_MProperDerivative( PreviousPoint );
	K1.fP       = fDeltaR * m_PDerivative( PreviousPoint );
	K1.fPhi     = fDeltaR * m_PhiDerivative( PreviousPoint );

	K1 = 0.5*K1 + PreviousPoint;

	K2.Zero();
	K2.fR       = fDeltaR;
	K2.fM       = fDeltaR * m_MDerivative( K1 );
	K2.fMProper = fDeltaR * m_MProperDerivative( K1 );
	K2.fP       = fDeltaR * m_PDerivative( K1 );
	K2.fPhi     = fDeltaR * m_PhiDerivative( K1 );

	return (PreviousPoint+K2);

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
	const RealType fStep(1.e-1);
	SInteriorSpaceTOV PointOfInterest;
	PointOfInterest = m_ApproximateSolutionInZero(0.01);
	cout<<DescInteriorSpaceTOV();
	cout<<PointOfInterest<<endl;
	while (PointOfInterest.fP>fAlmostZero)
	{
//		PointOfInterest = m_PropagateRKMethod(PointOfInterest,fStep);
		PointOfInterest = m_Propagate(PointOfInterest,fStep);
		cout<<PointOfInterest<<endl;
	}
}

void CPolytropeShooting::ComputeInteriorAdaptive()
{
	const RealType fAlmostZero ( 1.e-1);
	const RealType fEnough(1.e-5);
	RealType fStep(1.e-2);

	SInteriorSpaceTOV PointOfInterest, Small, Big;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<DescInteriorSpaceTOV();
	cout<<PointOfInterest<<" "<<fStep<<endl;
	while (PointOfInterest.fP>fAlmostZero)
	{
		do
		{
			fStep /= 2.;
			Small = m_PropagateRKMethod( m_PropagateRKMethod(PointOfInterest,fStep/2.) , fStep/2.);
			Big   = m_PropagateRKMethod( PointOfInterest,fStep);
			//Small = m_Propagate( m_PropagateRKMethod(PointOfInterest,fStep/2.) , fStep/2.);
			//Big   = m_Propagate( PointOfInterest,fStep);
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
	if(argc != 4)
	{
		cout<<"Wrong number of arguments."<<endl;
		cout<<"Correct usage (values in CGS):"<<endl;
		cout<<"./PolyTOVShooting fPolytropicGamma fPolytropicKappa(LoreneUnits) fRhoCentre"<<endl;
		return 1;
	}

	RealType fGamma(atof(argv[1]));
	RealType fKappa(atof(argv[2]));
	RealType fRhoCentre(atof(argv[3]));


	CPolytropeShooting Poly( fGamma, DeKappify(fKappa, fGamma), fRhoCentre);
	Poly.ComputeInteriorAdaptive();
	//Poly.ComputeInterior();


	return 0;
}

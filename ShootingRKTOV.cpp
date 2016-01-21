#include<iostream>
#include<cmath>
#include <limits>
#include<cstdlib>

using namespace std;

#define RealType double

struct SInteriorSpaceTOV
{
	RealType fXi;
	RealType fV;
	RealType fTheta;

	RealType fDVDXi;
	RealType fDThetaDXi;

	void Add(const RealType fDXi, const RealType fDV, const RealType fDTheta);

};

ostream &operator<<( ostream &output, const SInteriorSpaceTOV &Point )
{
	output << scientific;
	output << Point.fXi<<" "<<Point.fV<<" "<<Point.fTheta<<" "<<Point.fDVDXi<<" "<<Point.fDThetaDXi;
	return output;
}

void SInteriorSpaceTOV::Add(const RealType fDXi, const RealType fDV, const RealType fDTheta)
{
	fXi+=fDXi;
	fV+=fDV;
	fTheta+=fDTheta;
}

//how to impreve the name of this function???
RealType MultiDimRelativeDifference(const SInteriorSpaceTOV A, const SInteriorSpaceTOV B)
{
	return sqrt(pow((B.fXi-A.fXi)/B.fXi, 2.) + pow((B.fV-A.fV)/B.fV, 2.) + pow((B.fTheta-A.fTheta)/B.fTheta, 2.) + pow((B.fDVDXi-A.fDVDXi)/B.fDVDXi, 2.) + pow((B.fDThetaDXi-A.fDThetaDXi)/B.fDThetaDXi, 2.));
}

bool isnan(const SInteriorSpaceTOV A)
{
	return ( isnan(A.fXi) || isnan(A.fV) || isnan(A.fTheta) || isnan(A.fDVDXi) || isnan(A.fDThetaDXi) );
}

class CPolytropeShooting
{
protected:
	RealType m_VDerivative(const SInteriorSpaceTOV Point) const;
	RealType m_ThetaDerivative(const SInteriorSpaceTOV Point) const;

	SInteriorSpaceTOV m_ApproximateSolutionInZero(const RealType fEpsilon) const; 

	SInteriorSpaceTOV m_PropagateRKMethod(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaXi) const;

	RealType m_MassTilda(const SInteriorSpaceTOV Point) const;
public:
	RealType fPolytropicIndex;
	RealType fSigma;

	void ComputeInterior();
	void ComputeInteriorAdaptive();

	CPolytropeShooting();
	CPolytropeShooting(const RealType fN, const RealType fSig);
	~CPolytropeShooting(){};

};

CPolytropeShooting::CPolytropeShooting()
{
	CPolytropeShooting(3.,0.5);///testing
}

CPolytropeShooting::CPolytropeShooting(const RealType fN, const RealType fSig)
{
	fPolytropicIndex=fN;
	fSigma=fSig;
}

SInteriorSpaceTOV CPolytropeShooting::m_ApproximateSolutionInZero(const RealType fEpsilon) const
{
	SInteriorSpaceTOV Result;
	Result.fXi=fEpsilon;
	Result.fV=0.;
	Result.fTheta=1.;
	Result.fDVDXi=m_VDerivative(Result);
	Result.fDThetaDXi=m_ThetaDerivative(Result);
	return Result;
}


RealType CPolytropeShooting::m_VDerivative(const SInteriorSpaceTOV Point) const
{
	return Point.fXi*Point.fXi*pow(Point.fTheta,fPolytropicIndex);
}

RealType CPolytropeShooting::m_ThetaDerivative(const SInteriorSpaceTOV Point) const
{
	return -(1./(Point.fXi*Point.fXi))*(Point.fV + fSigma*pow(Point.fTheta,fPolytropicIndex+1)*Point.fXi*Point.fXi*Point.fXi  )*(1.+fSigma*Point.fTheta) / (1.-2.*fSigma*(fPolytropicIndex+1.)*Point.fV/Point.fXi ); 
}

SInteriorSpaceTOV CPolytropeShooting::m_PropagateRKMethod(const SInteriorSpaceTOV PreviousPoint, const RealType fDeltaXi) const
{
	SInteriorSpaceTOV WorkPoint;

	WorkPoint = PreviousPoint;

	RealType fK1( fDeltaXi * m_VDerivative( WorkPoint ) );
	RealType fL1( fDeltaXi * m_ThetaDerivative( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaXi/2., fK1/2., fL1/2. );

	RealType fK2( fDeltaXi * m_VDerivative( WorkPoint ) );
	RealType fL2( fDeltaXi * m_ThetaDerivative( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaXi/2., fK2/2., fL2/2. );


	RealType fK3( fDeltaXi * m_VDerivative( WorkPoint ) );
	RealType fL3( fDeltaXi * m_ThetaDerivative( WorkPoint ) );

	WorkPoint = PreviousPoint;
	WorkPoint.Add(fDeltaXi, fK3, fL3 );

	RealType fK4( fDeltaXi * m_VDerivative( WorkPoint ) );
	RealType fL4( fDeltaXi * m_ThetaDerivative( WorkPoint ) );


	SInteriorSpaceTOV Result;

	Result.fXi = PreviousPoint.fXi + fDeltaXi;
	Result.fV = PreviousPoint.fV + fK1/6. + fK2/3. + fK3/3. +fK4/6.;
	Result.fTheta = PreviousPoint.fTheta + fL1/6. + fL2/3. + fL3/3. +fL4/6.;
	Result.fDVDXi =  m_VDerivative( Result );
	Result.fDThetaDXi = m_ThetaDerivative( Result );


	return Result;
}
	
RealType CPolytropeShooting::m_MassTilda(const SInteriorSpaceTOV Point) const
{
	return pow(fSigma,(3.-fPolytropicIndex)/2.)*Point.fV;
}


void CPolytropeShooting::ComputeInterior()
{
	//const RealType fAlmostZero ( 10.*std::numeric_limits<RealType>::epsilon() );

	const RealType fAlmostZero ( 1e-5);
	const RealType fStep(1.e-4);
	SInteriorSpaceTOV PointOfInterest;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<endl;
	while (PointOfInterest.fTheta>fAlmostZero)
	{
		PointOfInterest = m_PropagateRKMethod(PointOfInterest,fStep);
		cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<endl;
	}
}


void CPolytropeShooting::ComputeInteriorAdaptive()
{
	const RealType fAlmostZero ( 1e-4);
	const RealType fEnough(1.e-1);
	RealType fStep(1.e-4);

	SInteriorSpaceTOV PointOfInterest, Small, Big;
	PointOfInterest = m_ApproximateSolutionInZero(fAlmostZero);
	cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
	while (PointOfInterest.fTheta>fAlmostZero)
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
		while( (MultiDimRelativeDifference(Small,Big) > fEnough) || isnan(Big) );
		PointOfInterest=Small;
		cout<<PointOfInterest<<" "<<m_MassTilda(PointOfInterest)<<" "<<fStep<<endl;
		fStep *= 4.;
	}
}


int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout<<"Wrong arguments"<<endl;
		cout<<"./PolyTOVShooting fPolytropicIndex fSigma"<<endl;
		return 1;
	}

	CPolytropeShooting Poly(atof(argv[1]),atof(argv[2]));
	Poly.ComputeInteriorAdaptive();
	//Poly.ComputeInterior();


	return 0;
}

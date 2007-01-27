#include "pbrt.h"
#include "transport.h"
#include "exoctree.h"
#include "bssrdfmaterial.h"

BSSRDFMaterial::BSSRDFMaterial( const Spectrum& aSigmaPrimeS, const Spectrum& aSigmaA, float aEta, float aScaleFactor ) 
: sigmaPrimeS( aSigmaPrimeS )
, sigmaA( aSigmaA )
, eta( aEta )
{
	sigmaA		*= aScaleFactor;
	sigmaPrimeS	*= aScaleFactor;

	sigmaPrimeT = (sigmaPrimeS + sigmaA);
	sigmaTr		= (sigmaA * sigmaPrimeT * 3).Sqrt();

	// the first, etai is the air's index of refraction, and the 
	// second one, the material's index of refraction
	// NOTE: file translucent.cpp line 60
	fresnel = new FresnelDielectric(1.0f, eta);

	// it's used the Luminance value of the Spectrum to compute the Mean Free Path
	//float sigmaTrLuminance = sigmaTr.y();
	float sigmaPrimeTLuminance = sigmaPrimeT.y();
	assert( sigmaPrimeTLuminance );

	lu			= (1 / sigmaPrimeTLuminance) * 3.0f;

	Fdr			= -(1.440f / (eta * eta)) + (0.710f / eta) + 0.668f + 0.0636f * eta;
	A			= (1 + Fdr) / (1 - Fdr);
	zv			= this->lu * (1 + 4/3 * A);

	printf("\n");
	printf("sigmaPrimeS: "); sigmaPrimeS.Print(stdout); printf("\n");
	printf("sigmaA: "); sigmaA.Print(stdout); printf("\n");
	printf("sigmaPrimeT: "); sigmaPrimeT.Print(stdout); printf("\n");
	printf("sigmaTR: "); sigmaTr.Print(stdout); printf("\n");
	printf("lu: %f", lu); printf("\n\n");
}

BSSRDFMaterial::~BSSRDFMaterial()
{
	delete fresnel;
}

BSDF* BSSRDFMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const
{
	DifferentialGeometry dgs;
	
	dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);

	float c[3] = {1.0f, 1.0f, 1.0f};
	Spectrum reflectance( c );
	MicrofacetDistribution* distribution = new Blinn(1.0f);
	
	// only add a SpecularReflection BSDF, the transmission is already 
	// contemplated on our implementation, right? :P
	//bsdf->Add(BSDF_ALLOC(SpecularReflection)(Spectrum(0.5f), fresnel));
	//bsdf->Add(BSDF_ALLOC(Microfacet(reflectance, fresnel, distribution)));
	bsdf->Add(BSDF_ALLOC(Lambertian)(reflectance));

	return bsdf;
}

extern "C" DLLEXPORT BSSRDFMaterial * CreateMaterial(const Transform &xform,
											   const TextureParams &mp) 
{
	Spectrum sigmaPrimeS = mp.FindSpectrum("sigmaPrimeS", Spectrum(1.f));
	Spectrum sigmaA = mp.FindSpectrum("sigmaA", Spectrum(1.f));
	float eta = mp.FindFloat("eta", 0.f);
	float scaleFactor = mp.FindFloat("scale", 1.f);
	
	return new BSSRDFMaterial(sigmaPrimeS, sigmaA, eta, scaleFactor);
}

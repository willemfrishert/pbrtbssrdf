#include "pbrt.h"
#include "transport.h"
#include "exoctree.h"
#include "bssrdfmaterial.h"

BSSRDFMaterial::BSSRDFMaterial( const Spectrum& sigmaPrimeS, const Spectrum& sigmaA, float eta ) : sigmaPrimeS( sigmaPrimeS )
, sigmaA( sigmaA )
, eta( eta )
{
	sigmaPrimeT = sigmaPrimeS + sigmaA;
	//// the first, etai is the air's index of refraction, and the 
	//// second one, the material's index of refraction
	//fresnel = new FresnelDielectric(1.0, eta);

	// the first, etai is the air's index of refraction, and the 
	// second one, the material's index of refraction
	// NOTE: file translucent.cpp line 60
	fresnel = new FresnelDielectric(1.0f, eta);

	// it's used the Luminance value of the Spectrum to compute the Mean Free Path
	//float sigmaTrLuminance = sigmaTr.y();
	//assert( sigmaTrLuminance );
	float sigmaPrimeTLuminance = sigmaPrimeT.y();
	assert( sigmaPrimeTLuminance );

	lu	= (1 / sigmaPrimeTLuminance);
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

	float c[3] = {0.8f, 0, 0};
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
	
	return new BSSRDFMaterial(sigmaPrimeS, sigmaA, eta);
}

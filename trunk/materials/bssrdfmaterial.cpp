#include "pbrt.h"
#include "transport.h"
#include "exoctree.h"
#include "bssrdfmaterial.h"

BSSRDFMaterial::BSSRDFMaterial(const Spectrum &sigmaPrimeS, const Spectrum &sigmaA, float lu, float eta)
: sigmaPrimeS( sigmaPrimeS )
, sigmaA( sigmaA )
, lu( lu )
, eta( eta )
{
	sigmaPrimeT = sigmaPrimeS + sigmaA;
	// the first, etai is the air's index of refraction, and the 
	// second one, the material's index of refraction
	fresnel = new FresnelDielectric(1.0, eta);
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
	
	// only add a SpecularReflection BSDF, the transmission is already 
	// contemplated on our implementation, right? :P
	bsdf->Add(BSDF_ALLOC(SpecularReflection)(Spectrum(0.5f), fresnel));

	return bsdf;
}

extern "C" DLLEXPORT BSSRDFMaterial * CreateMaterial(const Transform &xform,
											   const TextureParams &mp) 
{
	Spectrum sigmaPrimeS = mp.FindSpectrum("sigmaPrimeS", Spectrum(1.f));
	Spectrum sigmaA = mp.FindSpectrum("sigmaA", Spectrum(1.f));
	float eta = mp.FindFloat("eta", 0.f);
	float lu = mp.FindFloat("lu", 0.f);
	
	return new BSSRDFMaterial(sigmaPrimeS, sigmaA, lu, eta);
}

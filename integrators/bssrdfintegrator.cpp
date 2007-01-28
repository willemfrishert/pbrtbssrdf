/*
 * --> Based on the file Photonmap.cpp <--
 *
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// photonmap.cpp*
#include "pbrt.h"
#include "transport.h"
#include "scene.h"
#include "mc.h"
#include "kdtree.h"
#include "sampling.h"
#include "trianglemesh.h"
#include "exoctree.h"
#include "bssrdfmaterial.h"
#include "bssrdfintegrator.h"

using namespace std;

float BSSRDFIntegrator::maxOmega = 0.0f;
float BSSRDFIntegrator::minOmega = 0.0f;

struct Photon 
{
	// Photon Constructor
	Photon(const Point &pp, const Spectrum &wt, const Vector &w)
		: p(pp), alpha(wt), wi(w) {
	}
	Photon() { }
	
	// attributes
	Point p;

	Spectrum alpha;
	Vector wi;
};
struct PhotonProcess 
{
	// PhotonProcess Public Methods
	PhotonProcess(u_int mp, const Point &p);
	void operator()(const Photon &photon, float dist2, float &maxDistSquared) const;
	
	// attributes
	const Point &p;
	ClosePhoton *photons;
	u_int nLookup;
	mutable u_int foundPhotons;
};
struct ClosePhoton 
{
	ClosePhoton(const Photon *p = NULL, float md2 = INFINITY) 
	{
		photon = p;
		distanceSquared = md2;
	}
	bool operator<(const ClosePhoton &p2) const 
	{
		return distanceSquared == p2.distanceSquared ? (photon < p2.photon) :
			distanceSquared < p2.distanceSquared;
	}

	// attributes
	const Photon *photon;
	float distanceSquared;
};




// Photonmap Method Definitions
BSSRDFIntegrator::BSSRDFIntegrator(int ncaus, int ndir, int nind,
								   int nl,	int mdepth, float mdist, bool fg,
								   int gs, bool dp, int lightSamples, float epsilon, float nPointFactor) 
{
	nCausticPhotons = ncaus;
	nIndirectPhotons = nind;
	nDirectPhotons = ndir;
	nLookup = nl;
	maxDistSquared = mdist * mdist;
	maxSpecularDepth = mdepth;
	causticMap = directMap = indirectMap = NULL;
	specularDepth = 0;
	finalGather = fg;
	gatherSamples = gs;
	directWithPhotons = dp;
	this->lightSamples = lightSamples;
	this->epsilon = epsilon;
	this->maxOmega = 0.0f;
	this->minOmega = 0.0f;
	this->nPointFactor = nPointFactor;
}

BSSRDFIntegrator::~BSSRDFIntegrator() 
{
	delete causticMap;
	delete directMap;
	delete indirectMap;
}
void BSSRDFIntegrator::RequestSamples(Sample *sample,
									  const Scene *scene) 
{
	// Allocate and request samples for sampling all lights
	u_int nLights = scene->lights.size();
	lightSampleOffset = new int[nLights];
	bsdfSampleOffset = new int[nLights];
	bsdfComponentOffset = new int[nLights];

	for (u_int i = 0; i < nLights; ++i) 
	{
		const Light *light = scene->lights[i];
		int lightSamples =
			scene->sampler->RoundSize(light->nSamples);
		lightSampleOffset[i] = sample->Add2D(lightSamples);
		bsdfSampleOffset[i] = sample->Add2D(lightSamples);
		bsdfComponentOffset[i] = sample->Add1D(lightSamples);
	}

	lightNumOffset = -1;
	if (finalGather) 
	{
		gatherSamples = scene->sampler->RoundSize(gatherSamples);
		gatherSampleOffset = sample->Add2D(gatherSamples);
		gatherComponentOffset = sample->Add1D(gatherSamples);
	}

	// store a pointer to the Sample object
	this->mSample = sample;
}
/************************************************************************/
/* NEW FUNCTIONS														*/
/************************************************************************/
/**
 * @description since the samples are only generated within the call of Li() 
 * we use Stratified Sampling to generate samples to compute the 
 * Irradiance values for the Diffusion Approximation
 */
void BSSRDFIntegrator::GenerateStratifiedSamples()
{
	// Generate stratified samples for integrators
	for (u_int i = 0; i < this->mSample->n1D.size(); ++i)
	{
		for (u_int j = 0; j < this->mSample->n1D[i]; ++j)
		{
			this->mSample->oneD[i][j] = RandomFloat();
		}
	}
	
	for (u_int i = 0; i < this->mSample->n2D.size(); ++i)
	{
		for (u_int j = 0; j < 2*this->mSample->n2D[i]; ++j)
		{
			this->mSample->twoD[i][j] = RandomFloat();
		}
	}
}

/**
 * @param primitive
 * @return true if the material contained on this primitive is a translucent, i.e.,
 * a BSSRDF one.
 */
inline bool BSSRDFIntegrator::TranslucentMaterial( Reference<Primitive>& primitive, 
												   GeometricPrimitive** geoPrim, BSSRDFMaterial** material )
{
	*material = NULL;
	*geoPrim = BSSRDFIntegrator::Cast( primitive );
	if ( *geoPrim )
	{
		*material = BSSRDFIntegrator::Cast( (*geoPrim)->getMaterial() );
	}

	return ( *material != NULL );
}

inline bool BSSRDFIntegrator::TranslucentMaterial(Reference<Primitive> &primitive)
{
	GeometricPrimitive* prim;
	BSSRDFMaterial* material;
	
	return TranslucentMaterial(primitive, &prim, &material);
}

/**
 * @description ugly function that does a dynamic cast on an object. 
 * Hammer method used at full strength ;) in PBRT
 *
 * @param primitive
 * @return the object or NULL if is not of the return type
 */
inline GeometricPrimitive* BSSRDFIntegrator::Cast(Reference<Primitive>& primitive)
{
	Primitive* prim = primitive.operator ->();
	return dynamic_cast<GeometricPrimitive*> (prim);
}

/**
 * @description ugly function that does a dynamic cast on an object. 
 * Hammer method used at full strength ;) in PBRT
 *
 * @param material
 * @return the object or NULL if is not of the return type
 */
inline BSSRDFMaterial* BSSRDFIntegrator::Cast(Reference<Material>& material)
{
	Material* mat = material.operator ->();
	return dynamic_cast<BSSRDFMaterial*> (mat);
}

/**
* @description ugly function that does a dynamic cast on an object. 
* Hammer method used at full strength ;) in PBRT
*
* @param material
* @return the object or NULL if is not of the return type
*/
//inline const BSSRDFMaterial* BSSRDFIntegrator::Cast(const Reference<Material>& material)
//{
//	const Material* mat = material.operator ->();
//	return dynamic_cast<const BSSRDFMaterial*> (mat);
//}

/**
 * @param r
 * @param scene
 * @return 
 */
Spectrum BSSRDFIntegrator::ComputeRadianceEstimateAlongRay( RayDifferential& r, const Scene *scene ) const
{
	Spectrum L(0.0f);

	// Declare common path integration variables
	Spectrum pathThroughput = 1.0f;
	RayDifferential ray(r);
	bool specularBounce = false;

	for (int pathLength = 0; ; ++pathLength) 
	{
		// Find next vertex/point of path
		Intersection isect;
		if (!scene->Intersect(ray, &isect))
			break;

		// ### if first time store first ray ?? length ??
		if (pathLength == 0)
			r.maxt = ray.maxt;

		// ### if there is some participating media: absorption, scattering, etc...
		pathThroughput *= scene->Transmittance(ray);

		// ### Possibly add emitted light at path vertex: only if light source
		if (specularBounce)
			L += pathThroughput * isect.Le(-ray.d);

		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);

		// Sample illumination from lights to find path contribution
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;
		Vector wo = -ray.d;

		// mSample is in fact a pointer to the sample
		// NOTE: the sample is only passed to the Integrator on the Li method
		L += pathThroughput *
			UniformSampleOneLight(scene, p, n, wo, bsdf, this->mSample);

		// *********************** ################## ********** TMP ######
		int maxIndirectDepth = 5;

		// ### if reached max depth stop: NO RUSSIAN ROULETTE????
		if (pathLength + 1 == maxIndirectDepth) break;

		// Sample BSDF to get new path direction
		// Get random numbers for sampling new direction, \mono{bs1}, \mono{bs2}, and \mono{bcs}
		float bs1 = RandomFloat(), bs2 = RandomFloat(), bcs = RandomFloat();
		Vector wi;
		float pdf;
		BxDFType flags;
		
		// ### used for BxDF's with perfect specular features: book page 335
		Spectrum f = bsdf->Sample_f(wo, &wi, bs1, bs2, bcs,
			&pdf, BSDF_ALL, &flags);
		if (f.Black() || pdf == 0.)
			break;
		
		specularBounce = (flags & BSDF_SPECULAR) != 0;

		// ### Monte carlo evaluation (f * cos()) / pdf
		pathThroughput *= f * AbsDot(wi, n) / pdf;
		ray = RayDifferential(p, wi);

		// ### Possibly terminate the path: Russian Roulette
		if (pathLength > 3) 
		{
			float continueProbability = 0.5f;
			if (RandomFloat() > continueProbability)
				break;
			// ### ?????????????? ###
			pathThroughput /= continueProbability;
		}
	}

	return L;
}

/**
 * @description Computes the initial BSDF for a given point on a given shape using
 * a snippet from Triangle::Intersect().
 *
 * @param p
 * @param n
 * @param shape
 * @param material
 * @return 
 */
BSDF* BSSRDFIntegrator::ComputeInitialBSDF( const Point& p, const Normal& n, Reference<Shape>& shape, Reference<Material>& material ) const
{
	DifferentialGeometry dg;
	shape->GetDifferentialGeometry(p, &dg);

	// *** IF WE SUPPOSE THAT THE OBJECT DOESN'T CONTAIN ANY s OR n values
	return material->GetBSDF(dg, dg);
}

/**
 * @description Computes the irradiance for a BSSRDF material at a certain point
 *
 * @param p
 * @param n
 * @param pArea
 * @param triangle
 * @param material
 * @param scene
 */
void BSSRDFIntegrator::ComputeIrradiance( const Point &p, const Normal &n, float pArea, 
										 Reference<Shape>& triangle, Reference<Material>& material, 
										 const Scene *scene )
{
	Spectrum E;

	// Compute irradiance at current point
	u_int scramble[2] = { RandomUInt(), RandomUInt() };
	
	// *********************** ############## ********************/
	// COMPUTING THE TEMPORARY BSDF
	BSDF *bsdf = ComputeInitialBSDF(p, n, triangle, material );

	// Sample the hemisphere, computing the radiance for n sample 
	// directions in order to compute the irradiance
	for (int i = 0; i < this->lightSamples; ++i) 
	{
		float u[2];
		Sample02(i, scramble, u);

		// ### random direction over the hemisphere: initial ray
		Vector w = CosineSampleHemisphere(u[0], u[1]);

		// ### ray from the point on w direction
		//RayDifferential r(p, bsdf->LocalToWorld(w));

		w = bsdf->LocalToWorld(w);

		// ### flip ray if going on the wrong direction
		//if (Dot(r.d, n) < 0) r.d = -r.d;
		if (Dot(w, n) < 0) w = -w;

		// generate new sample values for the Direct and Indirect illumination computation
		GenerateStratifiedSamples();

		// ### Do path tracing to compute radiance along ray for estimate
		/************************************************************************/
		/* TODO: UNCOMMENT THIS LINE                                            */
		/************************************************************************/
		//E += ComputeRadianceEstimateAlongRay(r, scene);

		Spectrum directLight = UniformSampleAllLights(scene, p, n, w, bsdf, this->mSample, 
			this->lightSampleOffset, this->bsdfSampleOffset,
			this->bsdfComponentOffset);

		E += directLight;
	}

	// ### Finally, estimate the irradiance based on cosine-weighted distribution 
	// ### of directions: see Monte Carlo basic Estimator
	E *= M_PI / float( this->lightSamples );

	// Compute bounding box of irradiance sample: seen as a point mass 
	// because the hierarchical evaluation will take into account all the samples
	BBox sampleExtent( p );

	/************************************************************************/
	/* TODO: CHECK THIS PART TO CLEAN/OPTIMIZE: PROC/SAMPLE                 */
	/************************************************************************/

	BSSRDFMaterial* bssrdfMaterial = BSSRDFIntegrator::Cast( material );
	IrradBSSRDFSample irradSample(E, p, pArea);
	IrradBSSRDFProcess irradProcess( bssrdfMaterial, this->epsilon );

	// just to make sure :)
	assert( bssrdfMaterial );

	bssrdfMaterial->bssrdfIrradianceValues->Add(irradSample, sampleExtent, irradProcess);
}

void BSSRDFIntegrator::Postprocess(const Scene *scene)
{
	//vector<UniformPoint> container;

	//// Fetch all the BSSRDF-type objects from the scene
	//vector< Reference<GeometricPrimitive> > bssrdfObjects;
	//FindBSSRDFObjects(scene, bssrdfObjects);

	//vector< Reference<GeometricPrimitive> >::iterator primitiveIt = bssrdfObjects.begin();
	//BSSRDFMaterial* bssrdfMaterial = NULL;

	//FILE* f;
	//fopen_s(&f , "final_octree.out", "w");

	//for (int i = 0; primitiveIt != bssrdfObjects.end(); primitiveIt++, i++)
	//{
	//	Reference<GeometricPrimitive> primitive = *primitiveIt;
	//	Reference<Shape> mesh = primitive->getShape();
	//	bssrdfMaterial = BSSRDFIntegrator::Cast( primitive->getMaterial() );

	//	break;
	//}

	//if (bssrdfMaterial)
	//{
	//	bssrdfMaterial->bssrdfIrradianceValues->Print(f);
	//}

	//fclose(f);

	printf(">>>> Max Omega value: %f\n", maxOmega);
	printf(">>>> Min Omega value: %f\n\n", minOmega);
}

/**
 * @description Computes the irradiance values for all the BSSRDF objects in the scene
 *
 * @param scene
 */
void BSSRDFIntegrator::ComputeBSSRDFIrradianceValues( const Scene *scene )
{
	vector<UniformPoint> container;

	// Fetch all the BSSRDF-type objects from the scene
	vector< Reference<GeometricPrimitive> > bssrdfObjects;
	FindBSSRDFObjects(scene, bssrdfObjects);

	vector< Reference<GeometricPrimitive> >::iterator primitiveIt = bssrdfObjects.begin();
	BSSRDFMaterial* bssrdfMaterial = NULL;
	float pointArea;

	//FILE* f;
	//fopen_s(&f , "octree.out", "w");

	for (int i = 0; primitiveIt != bssrdfObjects.end(); primitiveIt++, i++)
	{
		Reference<GeometricPrimitive> primitive = *primitiveIt;
		Reference<Shape> mesh = primitive->getShape();
		bssrdfMaterial = BSSRDFIntegrator::Cast( primitive->getMaterial() );
		
		float lu = bssrdfMaterial->lu.y();
		mesh->GetUniformPointSamples(container, pointArea, lu / this->nPointFactor);

		/************************************************************************/
		/* TODO: RETRIEVE THE NUMBER OF SAMPLE POINTS!!!                        */
		/************************************************************************/
		int nPoints = container.size();
		int depth = static_cast<int>(ceil( log(static_cast<float>(nPoints)) / log(8.0f) ));
		printf("Number of Sample points: << %d >>\n", nPoints);
		printf("Optimal Tree Depth computed: << %d >>\n", depth);

		// Compute scene's BB, and extend it a little more 
		// (due to floating-point errors in scene intersections - p. 765): 
		// like in the IrradianceCache class; and finally build the octree
		BBox wb = mesh->WorldBound();
		Vector delta = 0.01f * (wb.pMax - wb.pMin);
		wb.pMin -= delta;
		wb.pMax += delta;

		bssrdfMaterial->bssrdfIrradianceValues = new ExOctree<IrradBSSRDFSample, IrradBSSRDFProcess>(wb, depth);

		vector<UniformPoint>::iterator pointIt = container.begin();

		ProgressReporter progress(container.size(), "Precomputing irradiance values");

		for (; pointIt != container.end(); pointIt++)
		{
			UniformPoint uniformPoint = *pointIt;
			ComputeIrradiance(uniformPoint.p, uniformPoint.n, pointArea, uniformPoint.triangle, 
				primitive->getMaterial(), scene);
			
			progress.Update();
		}
		progress.Done();

	}

	//if (bssrdfMaterial)
	//{
	//	bssrdfMaterial->bssrdfIrradianceValues->Print(f);

	//	/*IrradBSSRDFProcess irradProcess( bssrdfMaterial );
	//	bssrdfMaterial->bssrdfIrradianceValues->Lookup(Point(0.069f, 0.05f, 0.0f), irradProcess);*/
	//}

	//fclose(f);
}

/**
 * @param scene
 * @param container
 */
void BSSRDFIntegrator::FindBSSRDFObjects( const Scene *scene, vector< Reference<GeometricPrimitive> >& container )
{
	vector< Reference< Primitive > > primitives;
	
	// Get all the scene primitives by a "primitive" way :P
	//
	// *** NOTE: check the KD-tree Refine method, because it has been implemented on purpose ***
	scene->aggregate->Refine( primitives );

	vector< Reference< Primitive > >::iterator it = primitives.begin();

	// iterate through all the scene unrefined primitives
	for (; it != primitives.end(); it++)
	{
		Reference< Primitive > primitiveRef = *it;

		GeometricPrimitive* geoPrimitive;
		BSSRDFMaterial* material;

		// if truly a geoPrimitive AND the material 
		// has BSSRDF characteristics
		if ( BSSRDFIntegrator::TranslucentMaterial(primitiveRef, &geoPrimitive, &material) )
		{
			container.push_back( geoPrimitive );
		}
	}
}

/************************************************************************/
/*                                                                      */
/************************************************************************/

void BSSRDFIntegrator::Preprocess(const Scene *scene) 
{
	printf("Starting BSSRDF Integrator Preprocess\n");
	fflush(stdout);
	
	if (scene->lights.size() == 0) return;

	//
	// Precompute the BSSRDF irradiance values for all this kind of objects
	//
	ComputeBSSRDFIrradianceValues( scene );

	//
	// Start the Photon Mapping process
	//
	ProgressReporter progress(nCausticPhotons+nDirectPhotons+
		nIndirectPhotons, "Shooting photons");
	vector<Photon> causticPhotons;
	vector<Photon> directPhotons;
	vector<Photon> indirectPhotons;
	causticPhotons.reserve(nCausticPhotons);
	directPhotons.reserve(nDirectPhotons);
	indirectPhotons.reserve(nIndirectPhotons);

	// Initialize photon shooting statistics
	static StatsCounter nshot("Photon Map",
		"Number of photons shot from lights");
	bool causticDone = (nCausticPhotons == 0);
	bool directDone = (nDirectPhotons == 0);
	bool indirectDone = (nIndirectPhotons == 0);

	//DebugBreak();

	// ### Main while loop: shooting photons for caustics and photonmap ###
	while (!causticDone || !directDone || !indirectDone)
	{
		++nshot;
		
		// Give up if we're not storing enough photons
		if (nshot > 500000 &&
			(unsuccessful(nCausticPhotons,
			causticPhotons.size(),
			static_cast<int>( nshot )) ||
			unsuccessful(nDirectPhotons,
			directPhotons.size(),
			static_cast<int>( nshot )) ||
			unsuccessful(nIndirectPhotons,
			indirectPhotons.size(),
			static_cast<int>( nshot )))) 
		{
			Error("Unable to store enough photons.  Giving up.\n");
			return;
		}

		// Trace a photon path and store contribution
		// Choose 4D sample values for photon
		float u[4];
		u[0] = (float)RadicalInverse((int)nshot+1, 2);
		u[1] = (float)RadicalInverse((int)nshot+1, 3);
		u[2] = (float)RadicalInverse((int)nshot+1, 5);
		u[3] = (float)RadicalInverse((int)nshot+1, 7);

		// Choose light to shoot photon from
		int nLights = int(scene->lights.size());
		int lightNum =
			min(Floor2Int(nLights * (float)RadicalInverse((int)nshot+1, 11)),
			nLights-1);
		Light *light = scene->lights[lightNum];
		float lightPdf = 1.f / nLights;

		// Generate _photonRay_ from light source and initialize _alpha_
		RayDifferential photonRay;
		float pdf;
		Spectrum alpha =
			light->Sample_L(scene, u[0], u[1], u[2], u[3],
			&photonRay, &pdf);
		if (pdf == 0.f || alpha.Black()) continue;
		alpha /= pdf * lightPdf;

		if (!alpha.Black()) 
		{
			//
			// Follow photon path through scene and record intersections
			//
			bool specularPath = false;
			Intersection photonIsect;
			int nIntersections = 0;
			while (scene->Intersect(photonRay, &photonIsect)) 
			{
				++nIntersections;
				// Handle photon/surface intersection
				alpha *= scene->Transmittance(photonRay);
				Vector wo = -photonRay.d;
				BSDF *photonBSDF = photonIsect.GetBSDF(photonRay);
				BxDFType specularType = BxDFType(BSDF_REFLECTION |
					BSDF_TRANSMISSION | BSDF_SPECULAR);

				// not so good... but :)
				Reference<Primitive> primitiveRef = const_cast<Primitive*>( photonIsect.primitive );

				// if the hit surface has more components than the reflective 
				// ones, i.e., diffuse surface
				bool hasNonSpecular = (photonBSDF->NumComponents() >
					photonBSDF->NumComponents(specularType));
				
				if (hasNonSpecular && 
					! BSSRDFIntegrator::TranslucentMaterial( primitiveRef ) )
				{
					// Deposit photon at surface
					Photon photon(photonIsect.dg.p, alpha, wo);
					if (nIntersections == 1) 
					{
						// Process direct lighting photon intersection
						if (!directDone) 
						{
							directPhotons.push_back(photon);
							if (directPhotons.size() == nDirectPhotons) 
							{
								directDone = true;
								nDirectPaths = (int)nshot;
								directMap =
									new KdTree<Photon,
									PhotonProcess>(directPhotons);
							}
							progress.Update(); // NOBOOK
						}
					}
					else if (specularPath) 
					{
						// Process caustic photon intersection
						if (!causticDone) 
						{
							causticPhotons.push_back(photon);
							if (causticPhotons.size() == nCausticPhotons) 
							{
								causticDone = true;
								nCausticPaths = (int)nshot;
								causticMap =
									new KdTree<Photon,
									PhotonProcess>(causticPhotons);
							}
							progress.Update();
						}
					}
					else {
						// Process indirect lighting photon intersection
						if (!indirectDone) 
						{
							indirectPhotons.push_back(photon);
							if (indirectPhotons.size() == nIndirectPhotons) 
							{
								indirectDone = true;
								nIndirectPaths = (int)nshot;
								indirectMap =
									new KdTree<Photon,
									PhotonProcess>(indirectPhotons);
							}
							progress.Update();
						}
					}
				}

				// Sample new photon ray direction
				Vector wi;
				float pdf;
				BxDFType flags;

				// Get random numbers for sampling outgoing photon direction
				float u1, u2, u3;
				if (nIntersections == 1) 
				{
					u1 = (float)RadicalInverse((int)nshot+1, 13);
					u2 = (float)RadicalInverse((int)nshot+1, 17);
					u3 = (float)RadicalInverse((int)nshot+1, 19);
				}
				else 
				{
					u1 = RandomFloat();
					u2 = RandomFloat();
					u3 = RandomFloat();
				}
				Spectrum fr = photonBSDF->Sample_f(wo, &wi, u1, u2, u3,
					&pdf, BSDF_ALL, &flags);
				if (fr.Black() || pdf == 0.f)
					break;
				specularPath = (nIntersections == 1 || specularPath) &&
					((flags & BSDF_SPECULAR) != 0);
				alpha *= fr * AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
				photonRay = RayDifferential(photonIsect.dg.p, wi);

				// Possibly terminate photon path
				if (nIntersections > 3) 
				{
					float continueProbability = .5f;
					if (RandomFloat() > continueProbability)
						break;
					alpha /= continueProbability;
				}
			}
		}
		BSDF::FreeAll();
	}
	progress.Done(); // NOBOOK

	printf("Finished BSSRDF Integrator Preprocess\n");
	fflush(stdout);
}
Spectrum BSSRDFIntegrator::Li(const Scene *scene,
							  const RayDifferential &ray, const Sample *sample,
							  float *alpha) const 
{
	// Compute reflected radiance with photon map
	Spectrum L(0.);
	Intersection isect;

	if (scene->Intersect(ray, &isect)) 
	{
		// Evaluate BSDF at hit point
		BSDF *bsdf = isect.GetBSDF(ray);
		const Point &p = bsdf->dgShading.p;
		const Normal &n = bsdf->dgShading.nn;

		if (alpha) *alpha = 1.;

		/************************************************************************/
		/* THIS IS WHERE WE START :) TODO: process intersect with BSSRDF        */
		/************************************************************************/
		GeometricPrimitive* primitive;
		BSSRDFMaterial* material;

		// not so good... but :)
		Reference<Primitive> primitiveRef = const_cast<Primitive*>( isect.primitive );

		// if it is a Translucent Material (BSSRDF)
		if ( BSSRDFIntegrator::TranslucentMaterial( primitiveRef, &primitive, &material ) )
		{
			IrradBSSRDFProcess lookupProcess( material, this->epsilon );
			material->bssrdfIrradianceValues->Lookup(isect.dg.p, lookupProcess);

			BSSRDFIntegrator::maxOmega = max(BSSRDFIntegrator::maxOmega, lookupProcess.maxOmega);
			BSSRDFIntegrator::minOmega = min(BSSRDFIntegrator::minOmega, lookupProcess.minOmega);
			
			//
			// Return the computed Lo for the object: the cosi variable is cosine 
			// of the angle between the normal and incident ray
			//
			Vector wo = -ray.d;
			float cosi = Dot(wo, isect.dg.nn);
			
			Spectrum Lo = lookupProcess.Lo( cosi );

			return Lo;
		}
		

		/************************************************************************/
		/*                                                                      */
		/************************************************************************/
		Vector wo = -ray.d;

		// Compute emitted light if ray hit an area light source
		L += isect.Le(wo);

		// Compute direct lighting for photon map integrator
		if (directWithPhotons)
			L += LPhoton(directMap, nDirectPaths, nLookup,
			bsdf, isect, wo, maxDistSquared);
		else
			L += UniformSampleAllLights(scene, p, n,
			wo, bsdf, sample,
			lightSampleOffset, bsdfSampleOffset,
			bsdfComponentOffset);

		// Compute indirect lighting for photon map integrator
		L += LPhoton(causticMap, nCausticPaths, nLookup, bsdf,
			isect, wo, maxDistSquared);
		if (finalGather) 
		{
			// Do one-bounce final gather for photon map
			Spectrum Li(0.);
			for (int i = 0; i < gatherSamples; ++i) 
			{
				// Sample random direction for final gather ray
				Vector wi;
				float u1 = sample->twoD[gatherSampleOffset][2*i];
				float u2 = sample->twoD[gatherSampleOffset][2*i+1];
				float u3 = sample->oneD[gatherComponentOffset][i];
				float pdf;
				Spectrum fr = bsdf->Sample_f(wo, &wi, u1, u2, u3,
					&pdf, BxDFType(BSDF_ALL & (~BSDF_SPECULAR)));
				if (fr.Black() || pdf == 0.f) continue;
				RayDifferential bounceRay(p, wi);
				static StatsCounter gatherRays("Photon Map", // NOBOOK
					"Final gather rays traced"); // NOBOOK
				++gatherRays; // NOBOOK
				Intersection gatherIsect;
				if (scene->Intersect(bounceRay, &gatherIsect)) 
				{
					// Compute exitant radiance at final gather intersection
					BSDF *gatherBSDF = gatherIsect.GetBSDF(bounceRay);
					Vector bounceWo = -bounceRay.d;
					Spectrum Lindir =
						LPhoton(directMap, nDirectPaths, nLookup,
						gatherBSDF, gatherIsect, bounceWo, maxDistSquared) +
						LPhoton(indirectMap, nIndirectPaths, nLookup,
						gatherBSDF, gatherIsect, bounceWo, maxDistSquared) +
						LPhoton(causticMap, nCausticPaths, nLookup,
						gatherBSDF, gatherIsect, bounceWo, maxDistSquared);
					Lindir *= scene->Transmittance(bounceRay);
					Li += fr * Lindir * AbsDot(wi, n) / pdf;
				}
			}
			L += Li / float(gatherSamples);
		}
		else
			L += LPhoton(indirectMap, nIndirectPaths, nLookup,
			bsdf, isect, wo, maxDistSquared);
		if (specularDepth++ < maxSpecularDepth) 
		{
			Vector wi;

			// Trace rays for specular reflection and refraction
			Spectrum f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_REFLECTION | BSDF_SPECULAR));

			if (!f.Black()) 
			{
				// Compute ray differential _rd_ for specular reflection
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;
				// Compute differential reflected directions
				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy +
					bsdf->dgShading.dndv * bsdf->dgShading.dvdy;
				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);
				rd.rx.d = wi -
					dwodx + 2 * Vector(Dot(wo, n) * dndx +
					dDNdx * n);
				rd.ry.d = wi -
					dwody + 2 * Vector(Dot(wo, n) * dndy +
					dDNdy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
			f = bsdf->Sample_f(wo, &wi,
				BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));

			if (!f.Black()) 
			{
				// Compute ray differential _rd_ for specular transmission
				RayDifferential rd(p, wi);
				rd.hasDifferentials = true;
				rd.rx.o = p + isect.dg.dpdx;
				rd.ry.o = p + isect.dg.dpdy;

				float eta = bsdf->eta;
				Vector w = -wo;
				if (Dot(wo, n) < 0) eta = 1.f / eta;

				Normal dndx = bsdf->dgShading.dndu * bsdf->dgShading.dudx + bsdf->dgShading.dndv * bsdf->dgShading.dvdx;
				Normal dndy = bsdf->dgShading.dndu * bsdf->dgShading.dudy + bsdf->dgShading.dndv * bsdf->dgShading.dvdy;

				Vector dwodx = -ray.rx.d - wo, dwody = -ray.ry.d - wo;
				float dDNdx = Dot(dwodx, n) + Dot(wo, dndx);
				float dDNdy = Dot(dwody, n) + Dot(wo, dndy);

				float mu = eta * Dot(w, n) - Dot(wi, n);
				float dmudx = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdx;
				float dmudy = (eta - (eta*eta*Dot(w,n))/Dot(wi, n)) * dDNdy;

				rd.rx.d = wi + eta * dwodx - Vector(mu * dndx + dmudx * n);
				rd.ry.d = wi + eta * dwody - Vector(mu * dndy + dmudy * n);
				L += scene->Li(rd, sample) * f * AbsDot(wi, n);
			}
		}
		--specularDepth;
	}
	else 
	{
		// Handle ray with no intersection
		if (alpha) *alpha = 0.;
		for (u_int i = 0; i < scene->lights.size(); ++i)
			L += scene->lights[i]->Le(ray);
		if (alpha && !L.Black()) *alpha = 1.;
		return L;
	}
	return L;
}
Spectrum BSSRDFIntegrator::LPhoton(
								   KdTree<Photon, PhotonProcess> *map,
								   int nPaths, int nLookup, BSDF *bsdf,
								   const Intersection &isect, const Vector &wo,
								   float maxDistSquared) 
{
	Spectrum L(0.);
	if (!map) return L;
	BxDFType nonSpecular = BxDFType(BSDF_REFLECTION |
		BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY);
	if (bsdf->NumComponents(nonSpecular) == 0)
		return L;
	static StatsCounter lookups("Photon Map", "Total lookups"); // NOBOOK
	// Initialize _PhotonProcess_ object, _proc_, for photon map lookups
	PhotonProcess proc(nLookup, isect.dg.p);
	proc.photons =
		(ClosePhoton *)alloca(nLookup * sizeof(ClosePhoton));
	// Do photon map lookup
	++lookups;  // NOBOOK
	map->Lookup(isect.dg.p, proc, maxDistSquared);
	// Accumulate light from nearby photons
	static StatsRatio foundRate("Photon Map", "Photons found per lookup"); // NOBOOK
	foundRate.Add(proc.foundPhotons, 1); // NOBOOK
	float scale = 1.f / (float(nPaths) * maxDistSquared * M_PI);
	// Estimate reflected light from photons
	ClosePhoton *photons = proc.photons;
	int nFound = proc.foundPhotons;
	Normal Nf = Dot(wo, bsdf->dgShading.nn) < 0 ? -bsdf->dgShading.nn :
		bsdf->dgShading.nn;

	if (bsdf->NumComponents(BxDFType(BSDF_REFLECTION |
		BSDF_TRANSMISSION | BSDF_GLOSSY)) > 0) 
	{
		// Compute exitant radiance from photons for glossy surface
		for (int i = 0; i < nFound; ++i) {
			BxDFType flag = Dot(Nf, photons[i].photon->wi) > 0.f ?
BSDF_ALL_REFLECTION : BSDF_ALL_TRANSMISSION;
			L += bsdf->f(wo, photons[i].photon->wi,	flag) *
				(scale * photons[i].photon->alpha);
		}
	}
	else 
	{
		// Compute exitant radiance from photons for diffuse surface
		Spectrum Lr(0.), Lt(0.);
		for (int i = 0; i < nFound; ++i)
			if (Dot(Nf, photons[i].photon->wi) > 0.f)
				Lr += photons[i].photon->alpha;
			else
				Lt += photons[i].photon->alpha;
		L += (scale * INV_PI) * (Lr * bsdf->rho(wo, BSDF_ALL_REFLECTION) +
			Lt * bsdf->rho(wo, BSDF_ALL_TRANSMISSION));
	}
	return L;
}
PhotonProcess::PhotonProcess(u_int mp, const Point &P)
: p(P) 
{
	photons = 0;
	nLookup = mp;
	foundPhotons = 0;
}
void PhotonProcess::operator()(const Photon &photon,
							   float distSquared, float &maxDistSquared) const 
{
	static StatsPercentage discarded("Photon Map", "Discarded photons"); // NOBOOK
	discarded.Add(0, 1); // NOBOOK

	if (foundPhotons < nLookup) 
	{
		// Add photon to unordered array of photons
		photons[foundPhotons++] = ClosePhoton(&photon, distSquared);
		if (foundPhotons == nLookup) 
		{
			std::make_heap(&photons[0], &photons[nLookup]);
			maxDistSquared = photons[0].distanceSquared;
		}
	}
	else 
	{
		// Remove most distant photon from heap and add new photon
		discarded.Add(1, 0); // NOBOOK
		std::pop_heap(&photons[0], &photons[nLookup]);
		photons[nLookup-1] = ClosePhoton(&photon, distSquared);
		std::push_heap(&photons[0], &photons[nLookup]);
		maxDistSquared = photons[0].distanceSquared;
	}
}
extern "C" DLLEXPORT SurfaceIntegrator *CreateSurfaceIntegrator(const ParamSet &params) 
{
	int nCaustic		= params.FindOneInt("causticphotons", 20000);
	int nDirect			= params.FindOneInt("directphotons", 100000);
	int nIndirect		= params.FindOneInt("indirectphotons", 100000);
	int nUsed			= params.FindOneInt("nused", 50);
	int maxDepth		= params.FindOneInt("maxdepth", 5);
	bool finalGather	= params.FindOneBool("finalgather", true);
	bool directPhotons	= params.FindOneBool("directwithphotons", false);
	int gatherSamples	= params.FindOneInt("finalgathersamples", 32);
	float maxDist		= params.FindOneFloat("maxdist", .1f);
	int lightSamples	= params.FindOneInt("lightsamples", 128);
	float epsilon		= params.FindOneFloat("epsilon", 0.05f);
	float nPointFactor	= params.FindOneFloat("npointfactor", 3.0f);
	
	return new BSSRDFIntegrator(nCaustic, nDirect, nIndirect,
		nUsed, maxDepth, maxDist, finalGather, gatherSamples,
		directPhotons, lightSamples, epsilon, nPointFactor);
}
// Photonmap Local Declarations
struct Photon;
struct ClosePhoton;
struct PhotonProcess;

class TriangleMesh;
class Triangle;

class BSSRDFIntegrator : public SurfaceIntegrator {
public:
	// BSSRDFIntegrator Public Methods
	BSSRDFIntegrator(int ncaus, int ndir, int nindir, int nLookup, int mdepth,
		float maxdist, bool finalGather, int gatherSamples,
		bool directWithPhotons, int lightSamples, float epsilon, float nPointFactor);

	~BSSRDFIntegrator();

	Spectrum Li(const Scene *scene, const RayDifferential &ray,
		const Sample *sample, float *alpha) const;

	void RequestSamples(Sample *sample, const Scene *scene);

	void Preprocess(const Scene *);

public:
	static float maxOmega;
	static float minOmega;

private:
	// BSSRDFIntegrator Private Methods

	static inline bool unsuccessful(int needed, int found, int shot) {
		return (found < needed &&
			(found == 0 || found < shot / 1024));
	}

	static Spectrum LPhoton(KdTree<Photon, PhotonProcess> *map,
		int nPaths, int nLookup, BSDF *bsdf, const Intersection &isect,
		const Vector &w, float maxDistSquared);

	void ComputeBSSRDFIrradianceValues(const Scene *scene);
	
	void FindBSSRDFObjects(const Scene *scene, vector< Reference<GeometricPrimitive> >& container);
	
	void ComputeIrradiance(const Point &p, const Normal &n, float pArea, Reference<Shape>& triangle, Reference<Material>& material, const Scene *scene);
	
	Spectrum ComputeRadianceEstimateAlongRay(RayDifferential& r, const Scene *scene) const;
	
	BSDF* ComputeInitialBSDF(const Point& p, const Normal& n, Reference<Shape>& shape, Reference<Material>& material) const;

	inline static GeometricPrimitive* Cast(Reference<Primitive>& primitive);

	inline static BSSRDFMaterial* Cast(Reference<Material>& material);

	//inline static const BSSRDFMaterial* Cast(const Reference<Material>& material);

	inline static bool TranslucentMaterial(Reference<Primitive>& primitive, GeometricPrimitive** geoPrim, BSSRDFMaterial** material);

	inline static bool TranslucentMaterial(Reference<Primitive>& primitive);

	void GenerateStratifiedSamples();

	virtual void Postprocess(const Scene *scene);

	// BSSRDFIntegrator Private Data
	u_int nCausticPhotons, nIndirectPhotons, nDirectPhotons;
	u_int nLookup;
	mutable int specularDepth;
	int maxSpecularDepth;
	float maxDistSquared;
	bool directWithPhotons, finalGather;
	int gatherSamples;

	// Declare sample parameters for light source sampling
	int *lightSampleOffset, lightNumOffset;
	int *bsdfSampleOffset, *bsdfComponentOffset;
	int gatherSampleOffset, gatherComponentOffset;
	int nCausticPaths, nDirectPaths, nIndirectPaths;
	mutable KdTree<Photon, PhotonProcess> *causticMap;
	mutable KdTree<Photon, PhotonProcess> *directMap;
	mutable KdTree<Photon, PhotonProcess> *indirectMap;

	// Just going to store a pointer to the sample 
	// so I can use it on the Preprocess
	Sample* mSample;

	/**
	 * @description number of samples used to estimate 
	 * direct lighting for irradiance caching
	 */
	int lightSamples;

	/**
	 * @description maximum solid angle spanned by the points in a voxel
	 */
	float epsilon;

	/**
	 * @description ...nevermind; just to help computing the number of points
	 */
	float nPointFactor;
};
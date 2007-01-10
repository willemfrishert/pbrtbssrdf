// Photonmap Local Declarations
struct Photon;
struct ClosePhoton;
struct PhotonProcess;
struct IrradBSSRDFProcess;
struct IrradBSSRDFSample;

class BSSRDFIntegrator : public SurfaceIntegrator {
public:
	// BSSRDFIntegrator Public Methods
	BSSRDFIntegrator(int ncaus, int ndir, int nindir, int nLookup, int mdepth,
		float maxdist, bool finalGather, int gatherSamples,
		bool directWithPhotons);

	~BSSRDFIntegrator();

	Spectrum Li(const Scene *scene, const RayDifferential &ray,
		const Sample *sample, float *alpha) const;

	void RequestSamples(Sample *sample, const Scene *scene);

	void Preprocess(const Scene *);


private:
	// BSSRDFIntegrator Private Methods

	void storeIrradianceSamples(const Shape& bssrdfShape);

	static inline bool unsuccessful(int needed, int found, int shot) {
		return (found < needed &&
			(found == 0 || found < shot / 1024));
	}

	static Spectrum LPhoton(KdTree<Photon, PhotonProcess> *map,
		int nPaths, int nLookup, BSDF *bsdf, const Intersection &isect,
		const Vector &w, float maxDistSquared);

	void SampleBSSRDFIrradianceValues(const Scene *scene);
	
	void FindBSSRDFObjects(const Scene *scene, vector< Reference<GeometricPrimitive> >& container);
	
	void ComputeIrradiance(const Point &p, const Normal &n, float pArea, Reference<GeometricPrimitive>& primitive, const Scene *scene);
	
	Spectrum ComputeRadianceEstimateAlongRay(RayDifferential& r, const Scene *scene) const;
	
	BSDF* ComputeInitialBSDF(const Point& p, const Normal& n, Reference<GeometricPrimitive>& shape) const;

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

	/************************************************************************/
	/* OUR NEW STUFF                                                        */
	/************************************************************************/
	mutable Octree<IrradBSSRDFSample, IrradBSSRDFProcess>* bssrdfIrradianceValues;
};
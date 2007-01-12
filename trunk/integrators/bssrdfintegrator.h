// Photonmap Local Declarations
struct Photon;
struct ClosePhoton;
struct PhotonProcess;
//struct IrradBSSRDFProcess;
//struct IrradBSSRDFSample;

/************************************************************************/
/*                                                                      */
/************************************************************************/

struct IrradBSSRDFSample
{
	// IrradBSSRDFSample Constructor
	IrradBSSRDFSample()
		: Av( 0.0f )
	{ }

	IrradBSSRDFSample(const Spectrum &e, const Point &P, const float& a)
		: Ev(e), Pv(P), Av(a)
	{ }


	// The total irradiance on the node
	Spectrum Ev;

	// The total area represented by the point
	float Av;

	// The average location of the points, weighted by the irradiance
	Point Pv;
};

struct IrradBSSRDFProcess
{
	// IrradBSSRDFProcess Public Methods
	IrradBSSRDFProcess(float eps = 0.5f)
		: epsilon( eps )
	{
	}

	/**
	 * @param P
	 * @param sample
	 * @param childData
	 * @return true if the voxel shouldn't be subdivided, i.e., the lookup recursion should stop
	 */
	bool operator()(const Point &P, const vector<IrradBSSRDFSample> &sample, 
		const IrradBSSRDFSample* childData) const
	{
		// if it is not a leaf node
		if ( childData != NULL )
		{
			this->Av = childData->Av;
			this->Pv = childData->Pv;
			this->Ev = childData->Ev;
		} 
		else
		{
			// sum all the quantities on the leaf node
			vector<IrradBSSRDFSample>::const_iterator sampleIt = sample.begin();

			for (; sampleIt != sample.end(); sampleIt++)
			{
				this->Av = childData->Av;
				this->Ev = childData->Ev;
			}
		}
		

		

		Vector x( P );

		// As described on the paper: deltaW = Av / ||x - Pv||^2
		float delta = (x - Pv).LengthSquared();
		float omega = Av / delta;

		return ( omega > epsilon );
	}

	void addChildNode(ExOctNode<IrradBSSRDFSample> *node, const IrradBSSRDFSample &dataItem, 
		const BBox &dataBound, const BBox &nodeBound) const
	{
		/************************************************************************/
		/* HOW TO PROCEED IF IT IS A LEAF NODE???????                           */
		/************************************************************************/
		if ( ! node->childData )
		{
			node->childData = new IrradBSSRDFSample();
		}
		
		// Get both volumes and see how much of the node's volume the sample spans
		float nodeVolume	= nodeBound.Volume();
		float dataVolume	= dataBound.Volume();
		float ratio = 0.0f;
		
		// If not zero, or close to it
		if (nodeVolume > 0.0000000001f)
		{
			ratio = dataVolume / nodeVolume;
		}

		// Add the Irradiance, weighted by the volume it spans on the node
		node->childData->Ev += dataItem.Ev * ratio;
		
		// Add the point's area, weighted by the volume it spans on the node
		node->childData->Av += dataItem.Av * ratio;

		u_int t = node->childLeaves;
		float invT = 1 / static_cast<float>(t);

		// recursive average computation:
		// Add the average position weighted by the 'weighted irradiance'
		// It's assumed that the Spectrum's irradiance intensity is the Luminance value
		Point weightedPv = dataItem.Pv * dataItem.Ev.y() * ratio;
		
		node->childData->Pv = ((t - 1) * invT) * node->childData->Pv + ( invT * weightedPv);
	}

	/**
	 * @description maximum solid angle allowed to subdivide the voxels
	 */
	float epsilon;

	// The total irradiance on the node
	Spectrum Ev;

	// The total area represented by the point
	float Av;

	// The average location of the points, weighted by the irradiance
	Vector Pv;
};



class TriangleMesh;
class Triangle;

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

	static inline bool unsuccessful(int needed, int found, int shot) {
		return (found < needed &&
			(found == 0 || found < shot / 1024));
	}

	static Spectrum LPhoton(KdTree<Photon, PhotonProcess> *map,
		int nPaths, int nLookup, BSDF *bsdf, const Intersection &isect,
		const Vector &w, float maxDistSquared);

	void ComputeBSSRDFIrradianceValues(const Scene *scene);
	
	void FindBSSRDFObjects(const Scene *scene, vector< Reference<GeometricPrimitive> >& container);
	
	void ComputeIrradiance(const Point &p, const Normal &n, float pArea, Reference<Triangle>& triangle, Reference<Material>& material, const Scene *scene);
	
	Spectrum ComputeRadianceEstimateAlongRay(RayDifferential& r, const Scene *scene) const;
	
	BSDF* ComputeInitialBSDF(const Point& p, const Normal& n, Reference<Triangle>& shape, Reference<Material>& material) const;

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
	mutable ExOctree<IrradBSSRDFSample, IrradBSSRDFProcess>* bssrdfIrradianceValues;
};
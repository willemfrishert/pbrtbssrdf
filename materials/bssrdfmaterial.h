/**
 * @description The representation of a BSSRDF Material
 * 
 * @file bssrdfmaterial.h
 * @author João Pedro Jorge
 */
class Material;
class BSSRDFMaterial;

#include "exoctree.h"

/**
 * @description represents an irradiance sample
 */
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

/**
 * @description represents a lookup/add object. It's basically
 * the object that does the main decisions on the insertion and
 * lookup of the octree.
 */
struct IrradBSSRDFProcess
{
	// IrradBSSRDFProcess Public Methods
	IrradBSSRDFProcess(float eps = 0.0001f)
		: epsilon( eps )
	{
	}

	Spectrum Lo()
	{
		Spectrum maravilhaDoMundo;

		return maravilhaDoMundo;
	}

	void evaluate(const Point &P, const vector<IrradBSSRDFSample> &samples)
	{
		//if ( ! subdivide )
		//{
		//	// evaluate directly: get the node information about its children
		//	evaluateNode(P, samples[ 0 ]);

		//	return false;
		//}
	}

	/**
	* @param P
	* @param sample
	* @param childData
	* @return true if the voxel should be subdivided, i.e., the lookup recursion should continue
	*/
	bool subdivide(const Point &P, const vector<IrradBSSRDFSample> &samples)
	{
		// Get the node averaged values
		float Av = samples[ 0 ].Av;
		Vector Pv( samples[ 0 ].Pv );

		// As described on the paper: deltaW = Av / ||x - Pv||^2
		Vector x( P );
		float delta = (x - Pv).LengthSquared();
		float omega = Av / delta;

		// check if "voxel is small enough" to subdivide
		return omega > epsilon;
	}

	/**
	* @description adds information to an intermediate node regarding 
	*	a sample (future leaf) being inserted.
	*
	* @param node
	* @param dataItem
	* @param dataBound
	* @param nodeBound
	*/
	void addChildNode(ExOctNode<IrradBSSRDFSample> *node, const IrradBSSRDFSample &dataItem, 
		const BBox &dataBound, const BBox &nodeBound) const
	{
		IrradBSSRDFSample sample;

		// if not empty, it means it already has a sample 
		if ( ! node->data.empty() )
		{
			sample = node->data[ 0 ];
		}

		// Add up the point's area Av
		sample.Av += dataItem.Av;

		// Add up the Irradiance Ev
		sample.Ev += dataItem.Ev;

		u_int t = node->childLeaves;
		float invT = 1 / static_cast<float>(t);

		// recursive average computation:
		// Add the average position weighted by the 'weighted irradiance'
		// It's assumed that the Spectrum's irradiance intensity is the Luminance value
		Point weightedPv = dataItem.Pv * dataItem.Ev.y();

		sample.Pv = ((t - 1) * invT) * sample.Pv + ( invT * weightedPv);

		// if empty, push it into the vector
		if ( node->data.empty() )
		{
			node->data.push_back( sample );
		}
		else
		{
			node->data[ 0 ] = sample;
		}
	}

	/**
	* @description maximum solid angle allowed to subdivide the voxels
	*/
	float epsilon;

	// The total radiant exitance on the node
	Spectrum Mo;

	/**
	 * @description the material containing the parameters to compute the Radiance
	 */
	BSSRDFMaterial* material;
};

// BSSRDF Material Class Declarations
class BSSRDFMaterial : public Material {

	// methods
public:
	// BSSRDF Public Methods
	BSSRDFMaterial(const Spectrum& sigmaPrimeS, const Spectrum& sigmaA, float lu, float eta);

	virtual ~BSSRDFMaterial();

	BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading) const;

	// attributes
public:

	ExOctree<IrradBSSRDFSample, IrradBSSRDFProcess>* bssrdfIrradianceValues;

	// attributes
private:
	
	/**
	 * @description Reduced Scattering Coefficient
	 */
	Spectrum sigmaPrimeS;

	/**
	 * @description Reduced Scattering Coefficient
	 */
	Spectrum sigmaPrimeT;

	/**
	 * @description Absorption Coefficient
	 */
	Spectrum sigmaA;

	/**
	 * @description Mean Free Path
	 */
	float lu;

	/**
	 * @description Relative Index of Refraction
	 */
	float eta;

	/**
	 * @description Fresnel Dielectric term
	 * NOTE: the Evaluate() method returns the fresnel reflectance
	 * so, to get the transmittance, it should be done: Ft = 1 - Fr
	 */
	Fresnel* fresnel;
};
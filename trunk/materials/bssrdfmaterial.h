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
	IrradBSSRDFProcess(BSSRDFMaterial* material, float eps = 0.0001f);

	Spectrum Lo(float w);

	void evaluate(const Point &P, const vector<IrradBSSRDFSample> &samples);

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
	 * @description Effective Transport extinction coefficient
	 */
	Spectrum sigmaTr;

	/**
	 * @description Mean Free Path
	 */
	float lu;

	/**
	 * @description Relative Index of Refraction
	 */
	float eta;

	/**
	 * @description Distance to the real dipole light source
	 */
	float zr;

	/**
	 * @description Distance to the virtual dipole light source
	 */
	float zv;

	/**
	 * @description Diffuse Fresnel term
	 */
	float Fdr;

	/**
	 * @description ...
	 */
	float A;
	
	static const float inv4PI;

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

void IrradBSSRDFProcess::evaluate(const Point &P, const vector<IrradBSSRDFSample> &samples)
{
	vector<IrradBSSRDFSample>::const_iterator sampleIt = samples.begin();
	for (; sampleIt != samples.end(); sampleIt++)
	{
		IrradBSSRDFSample smp = *sampleIt;
		float r = (smp.Pv - P).Length();
		float rSqrt = r * r;
		float dr = sqrt(rSqrt + zr * zr);
		float dv = sqrt(rSqrt + zv * zv);
		
		Spectrum term1 = (sigmaTr + 1/dr) * zr;
		Spectrum term2 = Exp(-sigmaTr * dr) / (dr * dr);
		Spectrum term3 = (sigmaTr + 1/dv) * zv;
		Spectrum term4 = Exp(-sigmaTr * dv) / (dv * dv);

		Spectrum Rd = inv4PI * term1 * term2 * term3 * term4;

		/************************************************************************/
		/* TODO: Fdt can also be precomputed!!!!                                */
		/************************************************************************/
		Spectrum Fdt = 1 - Fdr;
		Mo += Fdt * Rd * smp.Ev * smp.Av;
	}

}

// static variable
float inv4PI = 1 / (4 * M_PI);

IrradBSSRDFProcess::IrradBSSRDFProcess(BSSRDFMaterial* material, float eps)
: epsilon( eps )
, material( material )
{
	/************************************************************************/
	/* TODO: Should I pass these computations to the MAterial???            */
	/************************************************************************/
	sigmaA		= this->material->sigmaA;
	sigmaPrimeT = this->material->sigmaPrimeT;
	sigmaPrimeS = this->material->sigmaPrimeS;
	sigmaTr		= (sigmaA * sigmaPrimeT * 3).Sqrt();
	zr			= this->material->lu;
	eta			= this->material->eta;
	Fdr			= -(1.440f / eta * eta) + (0.710f / eta) + 0.668f + 0.0636f * eta;
	A			= (1 + Fdr) / (1 - Fdr);
	zv			= this->material->lu * (1 + 4/3 * A);
	
}

Spectrum IrradBSSRDFProcess::Lo(float w)
{
	Spectrum Lo;

	// Fresnel transmittance: 1 - Fresnel reflectance
	Spectrum Ft = -this->material->fresnel->Evaluate(w) + 1;

	return (Ft / Fdr) * (Mo / M_PI);
}

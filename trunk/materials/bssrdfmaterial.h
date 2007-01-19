/**
 * @description The representation of a BSSRDF Material
 * 
 * @file bssrdfmaterial.h
 * @author Jo�o Pedro Jorge
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

	void Print(FILE* f)
	{
		fprintf(f, "Ev("); Ev.Print( f ); fprintf(f, ")"); 
		fprintf(f, " Pv(%.4f, %.4f, %.4f)", Pv.x, Pv.y, Pv.z);
	}

	static inline void CheckAverageLocation(IrradBSSRDFSample& sample);

	void operator+=(const IrradBSSRDFSample& sample)
	{
		this->Ev += sample.Ev;
		this->Av += sample.Av;
		this->Pv += sample.Pv;
	}


	// The total irradiance on the node
	Spectrum Ev;

	// The total area represented by the point
	float Av;

	// The average location of the points, weighted by the irradiance
	Point Pv;

	// True if the location is already averaged or it's a leaf node
	bool averaged;
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

	void evaluate(const Point &P, vector<IrradBSSRDFSample> &samples);

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
	 * @description adds information to an intermediate node about 
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
		IrradBSSRDFSample sample = dataItem;

		// Add the position weighted by the 'weighted irradiance'
		// It's assumed that the Spectrum's irradiance intensity is the Luminance value
		// It will be averaged when visited
		Point weightedPv = dataItem.Pv * dataItem.Ev.y();
		sample.Pv = weightedPv;

		// if empty, push it into the vector
		// it means it's the first time information is being added to the intermediate node
		if ( node->data.empty() )
		{
			sample.averaged = false;
			node->data.push_back( sample );
		}
		else // otherwise, just add it up
		{
			node->data[ 0 ] += sample;
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

// static variable
const float IrradBSSRDFProcess::inv4PI = 1 / (4 * M_PI);

// BSSRDF Material Class Declarations
class BSSRDFMaterial : public Material {

	// methods
public:
	// BSSRDF Public Methods
	BSSRDFMaterial(const Spectrum& sigmaPrimeS, const Spectrum& sigmaA, float eta);

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
	 * @description Relative Index of Refraction
	 */
	float eta;

	/**
	 * @description Fresnel Dielectric term
	 * NOTE: the Evaluate() method returns the fresnel reflectance
	 * so, to get the transmittance, it should be done: Ft = 1 - Fr
	 */
	Fresnel* fresnel;

	/**
	 * @description Mean Free Path
	 */
	float lu;
};

inline void IrradBSSRDFSample::CheckAverageLocation(IrradBSSRDFSample& sample)
{
	if ( ! sample.averaged )
	{
		sample.Pv /= sample.Ev.y();
		sample.averaged = true;
	}
}

void IrradBSSRDFProcess::evaluate(const Point &P, vector<IrradBSSRDFSample> &samples)
{
	vector<IrradBSSRDFSample>::iterator sampleIt = samples.begin();
	for (; sampleIt != samples.end(); sampleIt++)
	{
		IrradBSSRDFSample& smp = *sampleIt;

		// check if average location has already been computed
		IrradBSSRDFSample::CheckAverageLocation( smp );

		//float r = (smp.Pv - P).Length();
		float rSqr = (smp.Pv - P).LengthSquared();
		float dr = sqrt(rSqr + zr * zr);
		float dv = sqrt(rSqr + zv * zv);
		
		Spectrum C1 = (sigmaTr + 1/dr) * zr;
		Spectrum term2 = Exp(-sigmaTr * dr) / (dr * dr);
		Spectrum C2 = (sigmaTr + 1/dv) * zv;
		Spectrum term4 = Exp(-sigmaTr * dv) / (dv * dv);

		Spectrum Rd = inv4PI * C1 * term2 * C2 * term4;

		/************************************************************************/
		/* TODO: Fdt can also be precomputed!!!!                                */
		/************************************************************************/
		Spectrum Fdt = 1 - Fdr;
		//printf("Mo: Fdt("); Fdt.Print(stdout); printf(") "); 
		//printf("E("); smp.Ev.Print(stdout); printf(") ");
		//printf("Rd("); Rd.Print(stdout); printf(")-- ");
		//printf("C1("); C1.Print(stdout); printf(") -- ");
		//printf("term2("); term2.Print(stdout); printf(") -- ");
		//printf("C2("); C2.Print(stdout); printf(") -- ");
		//printf("term4("); term4.Print(stdout); printf(") -- ");
		//printf("MULT("); (C1 * term2 * C2 * term4).Print(stdout); printf(") -- ");
		
		Mo += Fdt * Rd * smp.Ev * smp.Av;
		/************************************************************************/
		/* TODO: REMOVE THE NEXT LINE. DEBUG ONLY                             */
		/************************************************************************/
		int i = 0;
		//printf("Mo("); Mo.Print(stdout); printf(")\n\n");
		//printf("*");
	}

}

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
	
	lu			= this->material->lu;
	zr			= this->lu;
	eta			= this->material->eta;
	Fdr			= -(1.440f / (eta * eta)) + (0.710f / eta) + 0.668f + 0.0636f * eta;
	A			= (1 + Fdr) / (1 - Fdr);
	zv			= this->lu * (1 + 4/3 * A);
	
}

Spectrum IrradBSSRDFProcess::Lo(float w)
{
	//FILE* f;
	//fopen_s(&f , "lo.out", "a");
	
	// Fresnel transmittance: 1 - Fresnel reflectance
	Spectrum Ft = -this->material->fresnel->Evaluate(w) + 1;

	Spectrum Lo = /*(Ft / Fdr) * */(Mo / M_PI);

	//fprintf(f, "Mo("); Mo.Print(f); fprintf(f, ")  ");
	//fprintf(f, "Lo("); Lo.Print(f); fprintf(f, ")\n");

	//fclose( f );

	return Lo;
}

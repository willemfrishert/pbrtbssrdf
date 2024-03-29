/*
* pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
*
* All Rights Reserved.
* For educational use only; commercial use expressly forbidden.
* NO WARRANTY, express or implied, for this software.
* (See file License.txt for complete license)
*/
#include "shape.h"
#include "paramset.h"

#include "core/sampling.h"

class PointRepulsion;
class Triangle;


// ###################### TriangleMesh class #############################
class TriangleMesh : public Shape 
{
public:
	// TriangleMesh Public Methods
	TriangleMesh(const Transform &o2w, bool ro, int ntris, int nverts, const int *vptr, const Point *P, const Normal *N, const Vector *S, const float *uv, const int aNumberOfIterations, const float aForceScalar, const string& aLoadFromFile, const string& aSaveToFile);
	~TriangleMesh();

	//	virtual void GetUniformPointSamples(vector<Point>& container) const;
	virtual void GetUniformPointSamples(vector<UniformPoint>& container, float& pointArea, float meanFreePath) const;

	BBox ObjectBound() const;
	BBox WorldBound() const;
	bool CanIntersect() const { return false; }

	void Refine(vector<Reference<Shape> > &refined) const;

	static Triangle* Cast(Reference<Shape>& aTriangle);
	void SaveSamplePoints(const vector<UniformPoint>& aContainer, const float& aPointArea) const;
	bool LoadFromFile(vector<Reference<Shape> > aTriangleList, vector<UniformPoint>& aContainer, float& aPointArea) const;
	void PrintSamplePointsToFile( vector<UniformPoint>& container, string aFileName ) const;


	friend class Triangle;
	template <class T> friend class VertexTexture;

protected:
	// TriangleMesh Data
	int ntris;			// Number of Triangles
	int nverts;			// Number of Vertices
	int *vertexIndex;	// Triangles
	Point *p;			// Vertices
	Normal *n;			// Normals
	Vector *s;
	float *uvs;

	int iNumberOfIterations;  // Number of iterations the point repulsion algorithm will take
	float iForceScale;		  // Scalar number
	string iLoadDataFile;	  // file to load sample points, normals and triangles from
	string iSaveDataFile;	  // file to save computed samplepoints, normals and triangles to
};

// ###################### Triangle class #############################
class Triangle : public Shape 
{
public:

	// Triangle Public Methods
	Triangle(const Transform &o2w, bool ro, TriangleMesh *m, int n)
		: Shape(o2w, ro) 
		, iID(n)
	{
		mesh = m;
		v = &mesh->vertexIndex[3*n];
		// Update created triangles stats
		static StatsCounter trisMade("Geometry", "Triangles created");
		++trisMade;
	}

	BBox ObjectBound() const;
	BBox WorldBound() const;

	bool Intersect(const Ray &ray, float *tHit,
		DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;

	void GetUVs(float uv[3][2]) const;

	float Area() const;
	
	int Id() const
	{
		return iID;
	}

	virtual void GetDifferentialGeometry(const Point& p, DifferentialGeometry *dg);

	virtual void GetShadingGeometry(const Transform &obj2world,const DifferentialGeometry &dg,DifferentialGeometry *dgShading) const 
	{
		if (!mesh->n && !mesh->s) 
		{
			*dgShading = dg;
			return;
		}

		// Initialize _Triangle_ shading geometry with _n_ and _s_
		// Compute barycentric coordinates for point
		float b[3];
		// Initialize _A_ and _C_ matrices for barycentric
		float uv[3][2];
		GetUVs(uv);
		float A[2][2] ={ { uv[1][0] - uv[0][0], uv[2][0] - uv[0][0] },
		{ uv[1][1] - uv[0][1], uv[2][1] - uv[0][1] } };
		float C[2] = { dg.u - uv[0][0], dg.v - uv[0][1] };

		if (!SolveLinearSystem2x2(A, C, &b[1])) 
		{
			// Handle degenerate parametric mapping
			b[0] = b[1] = b[2] = 1.f/3.f;
		}
		else
		{
			b[0] = 1.f - b[1] - b[2];
		}

		// Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
		Normal ns;
		Vector ss, ts;
		if (mesh->n)
		{
			ns = Normalize(obj2world(b[0] * mesh->n[v[0]] +	b[1] * mesh->n[v[1]] + b[2] * mesh->n[v[2]]));
		}
		else 
		{
			ns = dg.nn;
		}

		if (mesh->s)
		{
			ss = Normalize(obj2world(b[0] * mesh->s[v[0]] + b[1] * mesh->s[v[1]] + b[2] * mesh->s[v[2]]));
		}

		else 
		{
			ss = Normalize(dg.dpdu);
		}

		ts = Normalize(Cross(ss, ns));
		ss = Cross(ts, ns);
		Vector dndu, dndv;
		if (mesh->n) 
		{
			// Compute \dndu and \dndv for triangle shading geometry
			float uvs[3][2];
			GetUVs(uvs);
			// Compute deltas for triangle partial derivatives of normal
			float du1 = uvs[0][0] - uvs[2][0];
			float du2 = uvs[1][0] - uvs[2][0];
			float dv1 = uvs[0][1] - uvs[2][1];
			float dv2 = uvs[1][1] - uvs[2][1];
			Vector dn1 = Vector(mesh->n[v[0]] - mesh->n[v[2]]);
			Vector dn2 = Vector(mesh->n[v[1]] - mesh->n[v[2]]);
			float determinant = du1 * dv2 - dv1 * du2;
			if (determinant == 0)
			{
				dndu = dndv = Vector(0,0,0);
			}
			else 
			{
				float invdet = 1.f / determinant;
				dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
				dndv = (-du2 * dn1 + du1 * dn2) * invdet;
			}
		}
		else
		{
			dndu = dndv = Vector(0,0,0);
		}

		*dgShading = DifferentialGeometry(dg.p, ss, ts,	dndu, dndv, dg.u, dg.v, dg.shape);
		dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx; // NOBOOK
		dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy; // NOBOOK
		dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy; // NOBOOK
	}

	Point Sample(float u1, float u2, Normal *Ns) const;

	void GetVertexIndices(int** aIndices) const
	{
		*aIndices = v;
	}

	// methods
private:

	/************************************************************************/
	/* TODO: USE THIS METHOD ON THE INTERSECT AND INTERSECTP METHODS !!!!!! */
	/************************************************************************/

	inline void GetPartialDerivatives(Vector& dpdu, Vector& dpdv, const float uvs[3][2], 
		const Point& p1, const Point& p2, const Point& p3) const
	{
		Vector e1 = p2 - p1;
		Vector e2 = p3 - p1;

		// Compute deltas for triangle partial derivatives
		float du1 = uvs[0][0] - uvs[2][0];
		float du2 = uvs[1][0] - uvs[2][0];
		float dv1 = uvs[0][1] - uvs[2][1];
		float dv2 = uvs[1][1] - uvs[2][1];
		Vector dp1 = p1 - p3, dp2 = p2 - p3;
		float determinant = du1 * dv2 - dv1 * du2;

		if (determinant == 0.f) {
			// Handle zero determinant for triangle partial derivative matrix
			CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
		}
		else {
			float invdet = 1.f / determinant;
			dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
			dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
		}
	}

	// attributes
private:
	// Triangle Data
	Reference<TriangleMesh> mesh;
	int *v;
	int iID;
};
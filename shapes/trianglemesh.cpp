// trianglemesh.cpp*
#include "trianglemesh.h"
#include "pointrepulsion.h"

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
// TriangleMesh Method Definitions
// ###################### TriangleMesh Method Definitions #############################
TriangleMesh::TriangleMesh( const Transform &o2w, bool ro, int nt, int nv,
						   const int *vi, const Point *P, const Normal *N,
						   const Vector *S, const float *uv, 
						   const int aNumberOfIterations, const float aForceScalar ) 
: Shape(o2w, ro) 
, iNumberOfIterations(aNumberOfIterations)
, iForceScale(aForceScalar)
{
	ntris = nt;
	nverts = nv;
	vertexIndex = new int[3 * ntris];
	memcpy(vertexIndex, vi, 3 * ntris * sizeof(int));
	// Copy _uv_, _N_, and _S_ vertex data, if present
	if (uv) {
		uvs = new float[2*nverts];
		memcpy(uvs, uv, 2*nverts*sizeof(float));
	}
	else uvs = NULL;
	p = new Point[nverts];
	if (N) {
		n = new Normal[nverts];
		memcpy(n, N, nverts*sizeof(Normal));
	}
	else n = NULL;
	if (S) {
		s = new Vector[nverts];
		memcpy(s, S, nverts*sizeof(Vector));
	}
	else s = NULL;

	// Transform mesh vertices to world space
	for (int i  = 0; i < nverts; ++i)
		p[i] = ObjectToWorld(P[i]);

	//printf("Undo transformation to world space. Give me a copy of local space to debug :)\n");
	//printf("Don't forget to remove this code");
	//for (int i  = 0; i < nverts; ++i)
	//	p[i] = P[i];

	//vector<UniformPoint > container;
	//float i = 0;

	//this->GetUniformPointSamples(container, i, 2);
}

TriangleMesh::~TriangleMesh()
{
	delete[] vertexIndex;
	delete[] p;
	delete[] s;
	delete[] n;
	delete[] uvs;
}

BBox TriangleMesh::ObjectBound() const 
{
	BBox bobj;
	for (int i = 0; i < nverts; i++)
	{
		bobj = Union(bobj, WorldToObject(p[i]));
	}

	return bobj;
}

BBox TriangleMesh::WorldBound() const 
{
	BBox worldBounds;
	for (int i = 0; i < nverts; i++)
	{
		worldBounds = Union(worldBounds, p[i]);
	}

	return worldBounds;
}

void TriangleMesh::Refine(vector<Reference<Shape> > &refined) const 
{
	for (int i = 0; i < ntris; ++i)
	{
		refined.push_back(new Triangle(ObjectToWorld,
		                               reverseOrientation,
                                       (TriangleMesh *)this,
									   i));
	}
}

void TriangleMesh::GetUniformPointSamples( vector<UniformPoint>& container, float& pointArea, float meanFreePath ) const
{
	vector<Reference<Shape> > triangleList;	// container for the triangles
	Refine( triangleList );

	PointRepulsion pr( ntris, nverts, vertexIndex, p, triangleList);
	

	// create sample points based on the mean free path and the area of the object.
	// the area is already computed inside
	int numberOfSamplePoints = pr.SetupSamplePoints( meanFreePath );

	// save initial points
	//fstream outfile("mapped3D.txt", fstream::out);
	//outfile << "1" << endl;
	//outfile << numberOfSamplePoints << endl;
	//outfile.close();

	//pr.FillUniformSamplePointStructure( container );
	//PrintSamplePointsToFile( container, "mapped3D.txt" );
	//container.clear();



	//// setup header for turk points
	//outfile.open("turk3D.txt", fstream::out);
	//outfile << iNumberOfIterations << endl;
	//outfile << numberOfSamplePoints << endl;
	//outfile.close();

	ostringstream ostr;
	// relaxation of the random points 

	for(int k = 0; k < iNumberOfIterations; k++)
	{
		ostr.str("");
		ostr << "Computing repulsive forces " << k+1 << "/" << iNumberOfIterations;
		pr.ComputeRepulsiveForces( iForceScale, ostr.str() );
		ostr.str("");
		ostr << "Computing new point positions " << k+1 << "/" << iNumberOfIterations;
		pr.ComputeNewPositions( ostr.str() );

	////	//// collect the information
	////	//pr.FillUniformSamplePointStructure( container );
	////	//PrintSamplePointsToFile( container, "turk3D.txt" );
	////	//container.clear();

	}
	// collect the information
	pr.FillUniformSamplePointStructure( container );

	pointArea = PointRepulsion::CreateTriangleUseSets(triangleList, container, p);

	pointArea = pr.GetTotalSurfaceArea()/numberOfSamplePoints;
}


void TriangleMesh::PrintSamplePointsToFile( vector<UniformPoint>& container, string aFileName ) const
{
	fstream outfile;
	outfile.open(aFileName.c_str(), fstream::out | fstream::app);

	int maximum=container.size();
	for (int i=0; i<maximum; i++)
	{
		outfile << container.at(i).p.x << " " << container.at(i).p.y << " " << container.at(i).p.z;
		if ((i+1)<maximum)
		{
			outfile << " ";
		}
	}
	outfile << endl;
	outfile.close();
}

// ###################### Triangle Method Definitions #############################

BBox Triangle::ObjectBound() const 
{
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	return Union(BBox(WorldToObject(p1), WorldToObject(p2)),
		WorldToObject(p3));
}

BBox Triangle::WorldBound() const 
{
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	return Union(BBox(p1, p2), p3);
}

bool Triangle::Intersect(const Ray &ray, float *tHit,DifferentialGeometry *dg) const 
{
	// Initialize triangle intersection statistics
	static StatsPercentage triangleHits("Geometry", "Triangle Ray Intersections");

	// Update triangle tests count
	triangleHits.Add(0, 1);

	// Compute $\VEC{s}_1$
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];

	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);
	if (divisor == 0.)
	{
		return false;
	}
		
	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector d = ray.o - p1;
	float b1 = Dot(d, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
	{
		return false;
	}
	
	// Compute second barycentric coordinate
	Vector s2 = Cross(d, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
	{
		return false;
	}
	
	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
	{
		return false;
	}
		
	triangleHits.Add(1, 0); //NOBOOK

	// Fill in _DifferentialGeometry_ from triangle hit
	// Compute triangle partial derivatives
	Vector dpdu, dpdv;
	float uvs[3][2];
	GetUVs(uvs);

	// Compute deltas for triangle partial derivatives
	float du1 = uvs[0][0] - uvs[2][0];
	float du2 = uvs[1][0] - uvs[2][0];
	float dv1 = uvs[0][1] - uvs[2][1];
	float dv2 = uvs[1][1] - uvs[2][1];
	Vector dp1 = p1 - p3, dp2 = p2 - p3;
	float determinant = du1 * dv2 - dv1 * du2;
	if (determinant == 0.f) 
	{
		// Handle zero determinant for triangle partial derivative matrix
		CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
	}
	else
	{
		float invdet = 1.f / determinant;
		dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
		dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
	}

	// Interpolate $(u,v)$ triangle parametric coordinates
	float b0 = 1 - b1 - b2;
	float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
	float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
	*dg = DifferentialGeometry(ray(t), dpdu, dpdv, Vector(0,0,0), Vector(0,0,0), tu, tv, this);
	*tHit = t;

	return true;
}
bool Triangle::IntersectP(const Ray &ray) const 
{
	// Initialize triangle intersection statistics
	static StatsPercentage triangleHits("Geometry", "Triangle Ray Intersections");

	// Update triangle tests count
	triangleHits.Add(0, 1);

	// Compute $\VEC{s}_1$
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];

	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);
	if (divisor == 0.)
	{
		return false;
	}

	float invDivisor = 1.f / divisor;

	// Compute first barycentric coordinate
	Vector d = ray.o - p1;
	float b1 = Dot(d, s1) * invDivisor;
	if (b1 < 0. || b1 > 1.)
	{
		return false;
	}

	// Compute second barycentric coordinate
	Vector s2 = Cross(d, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
	{
		return false;
	}
		
	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
	{
		return false;
	}
		
	triangleHits.Add(1, 0); //NOBOOK

	return true;
}

/**
 * @param p
 * @param dg
 */
void Triangle::GetDifferentialGeometry( const Point& p, DifferentialGeometry *dg )
{
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	
	// Fill in _DifferentialGeometry_ from triangle hit
	// Compute triangle partial derivatives
	Vector dpdu, dpdv;
	float uvs[3][2];
	GetUVs(uvs);
	
	GetPartialDerivatives(dpdu, dpdv, uvs, p1, p2, p3);

	float A = p1.x - p3.x;
	float B = p2.x - p3.x;
	float C = p3.x - p.x;
	float F_I = p3.y - p.y + p3.z - p.z;
	float E_H = p2.y - p3.y + p2.z - p3.z;
	float D_G = p1.y - p3.y + p1.z - p3.z;

	float b1 = (B * F_I - C * E_H) / A * E_H - B * D_G;
	float b2 = (A * F_I - C * D_G) / B * D_G - A * E_H;

	// Interpolate $(u,v)$ triangle parametric coordinates
	float b0 = 1 - b1 - b2;
	float tu = b0*uvs[0][0] + b1*uvs[1][0] + b2*uvs[2][0];
	float tv = b0*uvs[0][1] + b1*uvs[1][1] + b2*uvs[2][1];
	*dg = DifferentialGeometry(p, dpdu, dpdv,
		Vector(0,0,0), Vector(0,0,0),
		tu, tv, this);
}

void Triangle::GetUVs(float uv[3][2]) const 
{
	if (mesh->uvs) 
	{
		uv[0][0] = mesh->uvs[2*v[0]];
		uv[0][1] = mesh->uvs[2*v[0]+1];
		uv[1][0] = mesh->uvs[2*v[1]];
		uv[1][1] = mesh->uvs[2*v[1]+1];
		uv[2][0] = mesh->uvs[2*v[2]];
		uv[2][1] = mesh->uvs[2*v[2]+1];
	} 
	else 
	{
		uv[0][0] = 0.; uv[0][1] = 0.;
		uv[1][0] = 1.; uv[1][1] = 0.;
		uv[2][0] = 1.; uv[2][1] = 1.;
	}
}

float Triangle::Area() const 
{
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];

	return 0.5f * Cross(p2-p1, p3-p1).Length();
}

Point Triangle::Sample(float u1, float u2, Normal *Ns) const
{
	float b1, b2;
	UniformSampleTriangle(u1, u2, &b1, &b2);

	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point &p1 = mesh->p[v[0]];
	const Point &p2 = mesh->p[v[1]];
	const Point &p3 = mesh->p[v[2]];
	Point p = b1 * p1 + b2 * p2 + (1.f - b1 - b2) * p3;

	Normal n = Normal(Cross(p2-p1, p3-p1));
	*Ns = Normalize(n);

	if (reverseOrientation) 
	{
 	 *Ns *= -1.f;
	}

	return p;
}

// ###################### external C function #############################

extern "C" DLLEXPORT Shape *CreateShape(const Transform &o2w, bool reverseOrientation, const ParamSet &params)
{
	int nvi, npi, nuvi, nsi, nni;
	const int *vi = params.FindInt("indices", &nvi);
	const Point *P = params.FindPoint("P", &npi);

	const float *uvs = params.FindFloat("uv", &nuvi);
	if (!uvs) 
	{
 	 uvs = params.FindFloat("st", &nuvi);
	}

	// XXX should complain if uvs aren't an array of 2...
	if (!vi || !P) 
	{
		return NULL;
	}
	
	const Vector *S = params.FindVector("S", &nsi);
	if (S && nsi != npi)
	{
		Error("Number of \"S\"s for triangle mesh must match \"P\"s");
		S = NULL;
	}

	const Normal *N = params.FindNormal("N", &nni);
	if (N && nni != npi)
	{
		Error("Number of \"N\"s for triangle mesh must match \"P\"s");
		N = NULL;
	}

	if (uvs && N)
	{
		// if there are normals, check for bad uv's that
		// give degenerate mappings; discard them if so
		const int *vp = vi;

		for (int i = 0; i < nvi; i += 3, vp += 3)
		{
			float area = .5f * Cross(P[vp[0]]-P[vp[1]], P[vp[2]]-P[vp[1]]).Length();
			if (area < 1e-7)
			{
				continue; // ignore degenerate tris.
			}
			
			if ((uvs[2*vp[0]] == uvs[2*vp[1]] &&
				uvs[2*vp[0]+1] == uvs[2*vp[1]+1]) ||
				(uvs[2*vp[1]] == uvs[2*vp[2]] &&
				uvs[2*vp[1]+1] == uvs[2*vp[2]+1]) ||
				(uvs[2*vp[2]] == uvs[2*vp[0]] &&
				uvs[2*vp[2]+1] == uvs[2*vp[0]+1])) 
			{
				Warning("Degenerate uv coordinates in triangle mesh.  Discarding all uvs.");
				uvs = NULL;
				break;
			}
		}
	}
	for (int i = 0; i < nvi; ++i)
	{
		if (vi[i] >= npi)
		{
			Error("trianglemesh has out of-bounds vertex index %d (%d \"P\" values were given",
				vi[i], npi);
			return NULL;
		}
	}

	int numberOfIterations = params.FindOneInt("pointrepulsioniterate", 5);
	float forceScale = params.FindOneFloat("forcescale", 0.1f);

	return new TriangleMesh(o2w, reverseOrientation, nvi/3, npi, vi, P,	N,
							S, uvs, numberOfIterations, forceScale);
}

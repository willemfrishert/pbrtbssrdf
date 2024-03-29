
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// paraboloid.cpp*
#include "shape.h"
// Paraboloid Declarations
class Paraboloid: public Shape {
public:
	// Paraboloid Public Methods
	Paraboloid(const Transform &o2w, bool ro, float rad,
	           float z0, float z1, float tm );
	BBox ObjectBound() const;
	bool Intersect(const Ray &ray, float *tHit,
	                  DifferentialGeometry *dg) const;
	bool IntersectP(const Ray &ray) const;
	float Area() const;
protected:
	// Paraboloid Data
	float radius;
	float zmin, zmax;
	float phiMax;
};
// Paraboloid Method Definitions
Paraboloid::Paraboloid(const Transform &o2w, bool ro,
                       float rad, float z0, float z1,
					   float tm)
	: Shape(o2w, ro) {
	radius = rad;
	zmin = min(z0,z1);
	zmax = max(z0,z1);
	phiMax = Radians( Clamp( tm, 0.0f, 360.0f ) );
}
BBox Paraboloid::ObjectBound() const {
	Point p1 = Point( -radius, -radius, zmin );
	Point p2 = Point(  radius,  radius, zmax );
	return BBox( p1, p2 );
}
bool Paraboloid::Intersect(const Ray &r, float *tHit,
		DifferentialGeometry *dg) const {
	float phi;
	Point phit;
	// Transform _Ray_ to object space
	Ray ray;
	WorldToObject(r, &ray);
	// Compute quadratic paraboloid coefficients
	float k = zmax/(radius*radius);
	float A =   k*(ray.d.x * ray.d.x + ray.d.y * ray.d.y);
	float B = 2*k*(ray.d.x * ray.o.x + ray.d.y * ray.o.y) -
	            ray.d.z;
	float C =   k*(ray.o.x * ray.o.x + ray.o.y * ray.o.y) -
	            ray.o.z;
	// Solve quadratic equation for _t_ values
	float t0, t1;
	if (!Quadratic(A, B, C, &t0, &t1))
		return false;
	// Compute intersection distance along ray
	if (t0 > ray.maxt || t1 < ray.mint)
		return false;
	float thit = t0;
	if (t0 < ray.mint) {
		thit = t1;
		if (thit > ray.maxt) return false;
	}
	// Compute paraboloid inverse mapping
	phit = ray(thit);
	phi = atan2f(phit.y, phit.x);
	if (phi < 0.) phi += 2.f*M_PI;
	// Test paraboloid intersection against clipping parameters
	if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
		if (thit == t1) return false;
		thit = t1;
		if (t1 > ray.maxt) return false;
		// Compute paraboloid inverse mapping
		phit = ray(thit);
		phi = atan2f(phit.y, phit.x);
		if (phi < 0.) phi += 2.f*M_PI;
		if (phit.z < zmin || phit.z > zmax || phi > phiMax)
			return false;
	}
	// Find parametric representation of paraboloid hit
	float u = phi / phiMax;
	float v = (phit.z-zmin) / (zmax-zmin);
	// Compute parabaloid \dpdu and \dpdv
	Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0.);
	Vector dpdv = (zmax - zmin) *
		Vector(phit.x / (2.f * phit.z), phit.y / (2.f * phit.z), 1.);
	// Compute parabaloid \dndu and \dndv
	Vector d2Pduu = -phiMax * phiMax *
	                Vector(phit.x, phit.y, 0);
	Vector d2Pduv = (zmax - zmin) * phiMax *
	                Vector(-phit.y / (2.f * phit.z),
					       phit.x / (2.f * phit.z),
						   0);
	Vector d2Pdvv = -(zmax - zmin) * (zmax - zmin) *
	                Vector(phit.x/(4.f*phit.z*phit.z),
					       phit.y/(4.f*phit.z*phit.z),
						   0.);
	// Compute coefficients for fundamental forms
	float E = Dot(dpdu, dpdu);
	float F = Dot(dpdu, dpdv);
	float G = Dot(dpdv, dpdv);
	Vector N = Normalize(Cross(dpdu, dpdv));
	float e = Dot(N, d2Pduu);
	float f = Dot(N, d2Pduv);
	float g = Dot(N, d2Pdvv);
	// Compute \dndu and \dndv from fundamental form coefficients
	float invEGF2 = 1.f / (E*G - F*F);
	Vector dndu = (f*F - e*G) * invEGF2 * dpdu +
		(e*F - f*E) * invEGF2 * dpdv;
	Vector dndv = (g*F - f*G) * invEGF2 * dpdu +
		(f*F - g*E) * invEGF2 * dpdv;
	// Initialize _DifferentialGeometry_ from parametric information
	*dg = DifferentialGeometry(ObjectToWorld(phit),
	                           ObjectToWorld(dpdu),
							   ObjectToWorld(dpdv),
	                           ObjectToWorld(dndu),
							   ObjectToWorld(dndv),
	                           u, v, this);
	// Update _tHit_ for quadric intersection
	*tHit = thit;
	return true;
}

bool Paraboloid::IntersectP(const Ray &r) const {
	float phi;
	Point phit;
	// Transform _Ray_ to object space
	Ray ray;
	WorldToObject(r, &ray);
	// Compute quadratic paraboloid coefficients
	float k = zmax/(radius*radius);
	float A =   k*(ray.d.x * ray.d.x + ray.d.y * ray.d.y);
	float B = 2*k*(ray.d.x * ray.o.x + ray.d.y * ray.o.y) -
	            ray.d.z;
	float C =   k*(ray.o.x * ray.o.x + ray.o.y * ray.o.y) -
	            ray.o.z;
	// Solve quadratic equation for _t_ values
	float t0, t1;
	if (!Quadratic(A, B, C, &t0, &t1))
		return false;
	// Compute intersection distance along ray
	if (t0 > ray.maxt || t1 < ray.mint)
		return false;
	float thit = t0;
	if (t0 < ray.mint) {
		thit = t1;
		if (thit > ray.maxt) return false;
	}
	// Compute paraboloid inverse mapping
	phit = ray(thit);
	phi = atan2f(phit.y, phit.x);
	if (phi < 0.) phi += 2.f*M_PI;
	// Test paraboloid intersection against clipping parameters
	if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
		if (thit == t1) return false;
		thit = t1;
		if (t1 > ray.maxt) return false;
		// Compute paraboloid inverse mapping
		phit = ray(thit);
		phi = atan2f(phit.y, phit.x);
		if (phi < 0.) phi += 2.f*M_PI;
		if (phit.z < zmin || phit.z > zmax || phi > phiMax)
			return false;
	}
	return true;
}
float Paraboloid::Area() const {
	return phiMax/12.0f *
		(powf(1+4*zmin, 1.5f) - powf(1+4*zmax, 1.5f));
}
extern "C" DLLEXPORT Shape *CreateShape(const Transform &o2w,
		bool reverseOrientation, const ParamSet &params) {
	float radius = params.FindOneFloat( "radius", 1 );
	float zmin = params.FindOneFloat( "zmin", 0 );
	float zmax = params.FindOneFloat( "zmax", 1 );
	float phimax = params.FindOneFloat( "phimax", 360 );
	return new Paraboloid(o2w, reverseOrientation, radius, zmin, zmax, phimax);
}

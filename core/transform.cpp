
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// transform.cpp*
#include "transform.h"
#include "shape.h"
// Transform Method Definitions
ostream &operator<<(ostream &os, const Transform &t) {
	t.m->Print(os);
	return os;
}
COREDLL Transform Translate(const Vector &delta) {
	Matrix4x4 *m, *minv;
	m = new Matrix4x4(1, 0, 0, delta.x,
                      0, 1, 0, delta.y,
                      0, 0, 1, delta.z,
                      0, 0, 0,       1);
	minv = new Matrix4x4(1, 0, 0, -delta.x,
                         0, 1, 0, -delta.y,
                         0, 0, 1, -delta.z,
                         0, 0, 0,        1);
	return Transform(m, minv);
}
COREDLL Transform Scale(float x, float y, float z) {
	Matrix4x4 *m, *minv;
	m = new Matrix4x4(x, 0, 0, 0,
                      0, y, 0, 0,
                      0, 0, z, 0,
                      0, 0, 0, 1);
	minv = new Matrix4x4(1.f/x,     0,     0, 0,
                             0, 1.f/y,     0, 0,
                             0,     0, 1.f/z, 0,
                             0,     0,     0, 1);
	return Transform(m, minv);
}
Transform RotateX(float angle) {
	float sin_t = sinf(Radians(angle));
	float cos_t = cosf(Radians(angle));
	Matrix4x4 *m = new Matrix4x4(1,     0,      0, 0,
                                 0, cos_t, -sin_t, 0,
                                 0, sin_t,  cos_t, 0,
                                 0,     0,      0, 1);
	return Transform(m, m->Transpose());
}
Transform RotateY(float angle) {
	float sin_t = sinf(Radians(angle));
	float cos_t = cosf(Radians(angle));
	Matrix4x4 *m = new Matrix4x4( cos_t,   0, sin_t, 0,
                                      0,   1,     0, 0,
                                 -sin_t,   0, cos_t, 0,
                                      0,   0,     0, 1);
	return Transform(m, m->Transpose());
}

Transform RotateZ(float angle) {
	float sin_t = sinf(Radians(angle));
	float cos_t = cosf(Radians(angle));
	Matrix4x4 *m = new Matrix4x4(cos_t, -sin_t, 0, 0,
                                 sin_t,  cos_t, 0, 0,
                                 0,      0, 1, 0,
                                 0,      0, 0, 1);
	return Transform(m, m->Transpose());
}
Transform Rotate(float angle, const Vector &axis) {
	Vector a = Normalize(axis);
	float s = sinf(Radians(angle));
	float c = cosf(Radians(angle));
	float m[4][4];

	m[0][0] = a.x * a.x + (1.f - a.x * a.x) * c;
	m[0][1] = a.x * a.y * (1.f - c) - a.z * s;
	m[0][2] = a.x * a.z * (1.f - c) + a.y * s;
	m[0][3] = 0;

	m[1][0] = a.x * a.y * (1.f - c) + a.z * s;
	m[1][1] = a.y * a.y + (1.f - a.y * a.y) * c;
	m[1][2] = a.y * a.z * (1.f - c) - a.x * s;
	m[1][3] = 0;

	m[2][0] = a.x * a.z * (1.f - c) - a.y * s;
	m[2][1] = a.y * a.z * (1.f - c) + a.x * s;
	m[2][2] = a.z * a.z + (1.f - a.z * a.z) * c;
	m[2][3] = 0;

	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 1;

	Matrix4x4 *mat = new Matrix4x4(m);
	return Transform(mat, mat->Transpose());
}

void Rotate(float angle, const Vector &axis, Reference<Matrix4x4>& aM, Reference<Matrix4x4>& aMinv)
{
	Vector a = Normalize(axis);
	float s = sinf(Radians(angle));
	float c = cosf(Radians(angle));
	float m[4][4];

	m[0][0] = a.x * a.x + (1.f - a.x * a.x) * c;
	m[0][1] = a.x * a.y * (1.f - c) - a.z * s;
	m[0][2] = a.x * a.z * (1.f - c) + a.y * s;
	m[0][3] = 0;

	m[1][0] = a.x * a.y * (1.f - c) + a.z * s;
	m[1][1] = a.y * a.y + (1.f - a.y * a.y) * c;
	m[1][2] = a.y * a.z * (1.f - c) - a.x * s;
	m[1][3] = 0;

	m[2][0] = a.x * a.z * (1.f - c) - a.y * s;
	m[2][1] = a.y * a.z * (1.f - c) + a.x * s;
	m[2][2] = a.z * a.z + (1.f - a.z * a.z) * c;
	m[2][3] = 0;

	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 1;

	aM = new Matrix4x4(m);
	aMinv = aM->Transpose();
}

Transform LookAt(const Point &pos, const Point &look, const Vector &up) {
	float m[4][4];
	// Initialize fourth column of viewing matrix
	m[0][3] = pos.x;
	m[1][3] = pos.y;
	m[2][3] = pos.z;
	m[3][3] = 1;
	// Initialize first three columns of viewing matrix
	Vector dir = Normalize(look - pos);
	Vector right = Normalize(Cross(dir, up));
	Vector newUp = Cross(right, dir);
	m[0][0] = right.x;
	m[1][0] = right.y;
	m[2][0] = right.z;
	m[3][0] = 0.;
	m[0][1] = newUp.x;
	m[1][1] = newUp.y;
	m[2][1] = newUp.z;
	m[3][1] = 0.;
	m[0][2] = dir.x;
	m[1][2] = dir.y;
	m[2][2] = dir.z;
	m[3][2] = 0.;
	Matrix4x4 *camToWorld = new Matrix4x4(m);
	return Transform(camToWorld->Inverse(), camToWorld);
}
bool Transform::HasScale() const {
#if 0
	float det = fabsf(m->m[0][0] * (m->m[1][1] * m->m[2][2] - m->m[1][2] * m->m[2][1])) -
		(m->m[0][1] * (m->m[1][0] * m->m[2][2] - m->m[1][2] * m->m[2][0])) +
		(m->m[0][2] * (m->m[1][0] * m->m[2][1] - m->m[1][1] * m->m[2][0]));
	return (det < .999f || det > 1.001f);
#endif
	return false;
}
BBox Transform::operator()(const BBox &b) const {
	const Transform &M = *this;
	BBox ret(       M(Point(b.pMin.x, b.pMin.y, b.pMin.z)));
	ret = Union(ret,M(Point(b.pMax.x, b.pMin.y, b.pMin.z)));
	ret = Union(ret,M(Point(b.pMin.x, b.pMax.y, b.pMin.z)));
	ret = Union(ret,M(Point(b.pMin.x, b.pMin.y, b.pMax.z)));
	ret = Union(ret,M(Point(b.pMin.x, b.pMax.y, b.pMax.z)));
	ret = Union(ret,M(Point(b.pMax.x, b.pMax.y, b.pMin.z)));
	ret = Union(ret,M(Point(b.pMax.x, b.pMin.y, b.pMax.z)));
	ret = Union(ret,M(Point(b.pMax.x, b.pMax.y, b.pMax.z)));
	return ret;
}
Transform Transform::operator*(const Transform &t2) const {
	Reference<Matrix4x4> m1 = Matrix4x4::Mul(m, t2.m);
	Reference<Matrix4x4> m2 = Matrix4x4::Mul(t2.mInv, mInv);
	return Transform(m1, m2);
}
bool Transform::SwapsHandedness() const {
	float det = ((m->m[0][0] *
                  (m->m[1][1] * m->m[2][2] -
                   m->m[1][2] * m->m[2][1])) -
                 (m->m[0][1] *
                  (m->m[1][0] * m->m[2][2] -
                   m->m[1][2] * m->m[2][0])) +
                 (m->m[0][2] *
                  (m->m[1][0] * m->m[2][1] -
                   m->m[1][1] * m->m[2][0])));
	return det < 0.f;
}
Transform COREDLL Orthographic(float znear, float zfar) {
	return Scale(1.f, 1.f, 1.f / (zfar-znear)) *
		Translate(Vector(0.f, 0.f, -znear));
}
COREDLL
Transform Perspective(float fov, float n, float f) {
	// Perform projective divide
	float inv_denom = 1.f/(f-n);
	Matrix4x4 *persp =
	    new Matrix4x4(1, 0,       0,          0,
	                  0, 1,       0,          0,
	                  0, 0, f*inv_denom, -f*n*inv_denom,
	                  0, 0,       1,          0);
	// Scale to canonical viewing volume
	float invTanAng = 1.f / tanf(Radians(fov) / 2.f);
	return Scale(invTanAng, invTanAng, 1) *
	       Transform(persp);
}

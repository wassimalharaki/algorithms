/**
 * 10:41:25 11/24/24
 * quaternion
 */
// ./algorithms/geometry/OrganizedTemplates/quaternion.cpp
#include <bits/stdc++.h>

using namespace std;

// To begin, we introduce the notation for the class and its elements, which we will use.
typedef double ftype;
typedef complex<ftype> point;
typedef complex<point> quater;
#define qA real()
#define qB imag()
#define qs qA.real()
#define qx qA.imag()
#define qy qB.real()
#define qz qB.imag()

const ftype pi = acos(-1.);
const ftype eps = 1e-9;

// Multiplication via formulas in complex terms.
quater operator * (quater a, quater b)
{
    return {{a.qA * b.qA - a.qB * conj(b.qB)},
            {a.qA * b.qB + a.qB * conj(b.qA)}};
}

// Conjugated element
quater conj(quater a)
{
    return {conj(a.qA), - a.qB};
}

// Quaternionic norm.
ftype norm(quater a)
{
    return (a * conj(a)).qs;
}

// Length.
ftype abs(quater a)
{
    return sqrt(norm(a));
}

// Dividing via inverse element.
quater operator / (quater a, quater b)
{
    return a * conj(b) / point(norm(b));
}

// Quaternion from vector coordinates.
quater vec(ftype x, ftype y, ftype z)
{
    return {{0, x}, {y, z}};
}

// Basis vectors.
const quater ex = vec(1, 0, 0);
const quater ey = vec(0, 1, 0);
const quater ez = vec(0, 0, 1);

// Vector part of quaternion.
quater vec(quater a)
{
    return a -= a.qs;
}

// Dot product.
ftype dot(quater a, quater b)
{
    return -(a * b).qs;
}

// Cross product.
quater cross(quater a, quater b)
{
    return vec(a * b);
}

// Triple (mixed) product.
ftype mix(quater a, quater b, quater c)
{
    return dot(a, cross(b, c));
}

// Conjugation of vector a by quaternion g.
quater conj(quater a, quater g)
{
    return g * a / g;
}

// Quaternion representing rotation around i on phi counter-clockwise.
quater rotation(quater i, ftype phi)
{
    return point(cos(phi / 2)) + i * point(sin(phi / 2));
}

// Rotate a around i counter-clockwise.
quater rotate(quater a, quater i, ftype phi)
{
    return conj(a, rotation(i, phi));
}

// Check if two vectors are equal.
bool cmp(const quater &a, const quater &b)
{
    return abs(a - b) < eps;
}

// Any vector which is orthogonal to v.
quater ortho(quater v)
{
    if(abs(v.qy) > eps)
        return vec(v.qy, -v.qx, 0);
    else
        return vec(v.qz, 0, -v.qx);
}

// Quaternion representing rotation from a to b via minimal angle.
quater min_rotation(quater a, quater b)
{
    a /= abs(a);
    b /= abs(b);
    if(cmp(a, -b)) // Degenerate case :(
        return rotation(ortho(b), pi);
    else
        return conj(a * (a + b));
}

// Angle of rotation from [-pi; pi] defined by quaternion.
ftype get_angle(quater a)
{
    a /= abs(a);
    return remainder(2 * acos(a.qs), 2 * pi);
}

// Quaternion representing rotation which turns nx in x, ny in y and [nx, ny] in z.
quater basis_rotation(quater nx, quater ny)
{
    nx /= abs(nx);
    ny /= abs(ny);
    quater a = min_rotation(nx, ex);
    ny = conj(ny, a);
    quater b = min_rotation(ny, ey);
    if(cmp(ny, -ey))
        b = rotation(ex, pi);
    return b * a;
}

// Copyright (c) 2015 burningmime
// 
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
// 
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgement in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.

module;

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include "glm_wrapper.hpp"

module bezierfit;

using namespace bezierfit;

CubicBezier::CubicBezier(const VECTOR& p0, const VECTOR& p1, const VECTOR& p2, const VECTOR& p3)
	: p0(p0), p1(p1), p2(p2), p3(p3)
{
}

VECTOR CubicBezier::Sample(FLOAT t) const
{
	FLOAT ti = 1.0 - t;
	FLOAT t0 = ti * ti * ti;
	FLOAT t1 = 3.0 * ti * ti * t;
	FLOAT t2 = 3.0 * ti * t * t;
	FLOAT t3 = t * t * t;
	return (t0 * p0) + (t1 * p1) + (t2 * p2) + (t3 * p3);
}

VECTOR CubicBezier::Derivative(FLOAT t) const
{
	FLOAT ti = 1.0 - t;
	FLOAT tp0 = 3.0 * ti * ti;
	FLOAT tp1 = 6.0 * t * ti;
	FLOAT tp2 = 3.0 * t * t;
	return (tp0 * (p1 - p0)) + (tp1 * (p2 - p1)) + (tp2 * (p3 - p2));
}

VECTOR CubicBezier::Tangent(FLOAT t) const
{
	return VectorHelper::Normalize(Derivative(t));
}

std::string CubicBezier::ToString() const
{
	std::ostringstream oss;
	oss << "CubicBezier: (<" << std::fixed << std::setprecision(3) << p0.x << ", " << p0.y << "> <"
		<< p1.x << ", " << p1.y << "> <" << p2.x << ", " << p2.y << "> <" << p3.x << ", " << p3.y << ">)";
	return oss.str();
}

// Equality members
bool CubicBezier::operator==(const CubicBezier& other) const
{
	return p0 == other.p0 && p1 == other.p1 && p2 == other.p2 && p3 == other.p3;
}

bool CubicBezier::operator!=(const CubicBezier& other) const
{
	return !(*this == other);
}

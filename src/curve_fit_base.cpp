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

#include "curve_fit.hpp"
#include "vector_helper.hpp"

using namespace bezierfit;

VECTOR CurveFitBase::GetLeftTangent(int last)
{
	int count = _pts.size();
	FLOAT totalLen = _arclen[count - 1];
	VECTOR p0 = _pts[0];
	VECTOR tanL = glm::normalize(_pts[1] - p0);
	VECTOR total = tanL;
	FLOAT weightTotal = 1;
	last = std::min(END_TANGENT_N_PTS, last - 1);
	for (int i = 2; i <= last; i++)
	{
		FLOAT ti = 1 - (_arclen[i] / totalLen);
		FLOAT weight = ti * ti * ti;
		VECTOR v = glm::normalize(_pts[i] - p0);
		total += v * weight;
		weightTotal += weight;
	}
	if (glm::length(total) > EPSILON)
		tanL = glm::normalize(total / weightTotal);
	return tanL;
}

VECTOR CurveFitBase::GetRightTangent(int first)
{
	int count = _pts.size();
	FLOAT totalLen = _arclen[count - 1];
	VECTOR p3 = _pts[count - 1];
	VECTOR tanR = glm::normalize(_pts[count - 2] - p3);
	VECTOR total = tanR;
	FLOAT weightTotal = 1;
	first = std::max(count - (END_TANGENT_N_PTS + 1), first + 1);
	for (int i = count - 3; i >= first; i--)
	{
		FLOAT t = _arclen[i] / totalLen;
		FLOAT weight = t * t * t;
		VECTOR v = glm::normalize(_pts[i] - p3);
		total += v * weight;
		weightTotal += weight;
	}
	if (glm::length(total) > EPSILON)
		tanR = glm::normalize(total / weightTotal);
	return tanR;
}

VECTOR CurveFitBase::GetCenterTangent(int first, int last, int split)
{
	int count = _pts.size();
	FLOAT splitLen = _arclen[split];
	VECTOR pSplit = _pts[split];

	// left side
	FLOAT firstLen = _arclen[first];
	FLOAT partLen = splitLen - firstLen;
	VECTOR total = VECTOR(0);
	FLOAT weightTotal = 0;
	for (int i = std::max(first, split - MID_TANGENT_N_PTS); i < split; i++)
	{
		FLOAT t = (_arclen[i] - firstLen) / partLen;
		FLOAT weight = t * t * t;
		VECTOR v = glm::normalize(_pts[i] - pSplit);
		total += v * weight;
		weightTotal += weight;
	}
	VECTOR tanL = glm::length(total) > EPSILON && weightTotal > EPSILON ?
		glm::normalize(total / weightTotal) :
		glm::normalize(_pts[split - 1] - pSplit);

	// right side
	partLen = _arclen[last] - splitLen;
	int rMax = std::min(last, split + MID_TANGENT_N_PTS);
	total = VECTOR(0);
	weightTotal = 0;
	for (int i = split + 1; i <= rMax; i++)
	{
		FLOAT ti = 1 - ((_arclen[i] - splitLen) / partLen);
		FLOAT weight = ti * ti * ti;
		VECTOR v = glm::normalize(pSplit - _pts[i]);
		total += v * weight;
		weightTotal += weight;
	}
	VECTOR tanR = glm::length(total) > EPSILON && weightTotal > EPSILON ?
		glm::normalize(total / weightTotal) :
		glm::normalize(pSplit - _pts[split + 1]);

	total = tanL + tanR;

	if (glm::length2(total) < EPSILON)
	{
		tanL = glm::normalize(_pts[split - 1] - pSplit);
		tanR = glm::normalize(pSplit - _pts[split + 1]);
		total = tanL + tanR;
		return glm::length2(total) < EPSILON ? tanL : glm::normalize(total / 2.0f);
	}
	else
	{
		return glm::normalize(total / 2.0f);
	}
}

void CurveFitBase::InitializeArcLengths()
{
	int count = _pts.size();
	_arclen.clear();
	_arclen.push_back(0);
	FLOAT clen = 0;
	VECTOR pp = _pts[0];
	for (int i = 1; i < count; i++)
	{
		VECTOR np = _pts[i];
		clen += glm::distance(pp, np);
		_arclen.push_back(clen);
		pp = np;
	}
}

void CurveFitBase::ArcLengthParamaterize(int first, int last)
{
	int count = _pts.size();
	_u.clear();
	FLOAT diff = _arclen[last] - _arclen[first];
	FLOAT start = _arclen[first];
	int nPts = last - first;
	_u.push_back(0);
	for (int i = 1; i < nPts; i++)
		_u.push_back((_arclen[first + i] - start) / diff);
	_u.push_back(1);
}

/// <summary>
 /// Generates a bezier curve for the segment using a least-squares approximation.
 /// </summary>
CubicBezier CurveFitBase::GenerateBezier(int first, int last, VECTOR tanL, VECTOR tanR)
{
	std::vector<VECTOR>& pts = _pts;
	std::vector<FLOAT>& u = _u;
	int nPts = last - first + 1;
	VECTOR p0 = pts[first], p3 = pts[last]; // first and last points of curve are actual points on data
	FLOAT c00 = 0, c01 = 0, c11 = 0, x0 = 0, x1 = 0; // matrix members -- both C[0,1] and C[1,0] are the same, stored in c01
	for (int i = 1; i < nPts; i++)
	{
		// Calculate cubic bezier multipliers
		FLOAT t = u[i];
		FLOAT ti = 1 - t;
		FLOAT t0 = ti * ti * ti;
		FLOAT t1 = 3 * ti * ti * t;
		FLOAT t2 = 3 * ti * t * t;
		FLOAT t3 = t * t * t;

		// For X matrix; moving this up here since profiling shows it's better up here (maybe a0/a1 not in registers vs only v not in regs)
		VECTOR s = (p0 * t0) + (p0 * t1) + (p3 * t2) + (p3 * t3); // NOTE: this would be Q(t) if p1=p0 and p2=p3
		VECTOR v = pts[first + i] - s;

		// C matrix
		VECTOR a0 = tanL * t1;
		VECTOR a1 = tanR * t2;
		c00 += VectorHelper::Dot(a0, a0);
		c01 += VectorHelper::Dot(a0, a1);
		c11 += VectorHelper::Dot(a1, a1);

		// X matrix
		x0 += VectorHelper::Dot(a0, v);
		x1 += VectorHelper::Dot(a1, v);
	}

	// determinants of X and C matrices
	FLOAT det_C0_C1 = c00 * c11 - c01 * c01;
	FLOAT det_C0_X = c00 * x1 - c01 * x0;
	FLOAT det_X_C1 = x0 * c11 - x1 * c01;
	FLOAT alphaL = det_X_C1 / det_C0_C1;
	FLOAT alphaR = det_C0_X / det_C0_C1;

	// if alpha is negative, zero, or very small (or we can't trust it since C matrix is small), fall back to Wu/Barsky heuristic
	FLOAT linDist = VectorHelper::Distance(p0, p3);
	FLOAT epsilon2 = EPSILON * linDist;
	if (std::abs(det_C0_C1) < EPSILON || alphaL < epsilon2 || alphaR < epsilon2)
	{
		FLOAT alpha = linDist / 3;
		VECTOR p1 = (tanL * alpha) + p0;
		VECTOR p2 = (tanR * alpha) + p3;
		return CubicBezier(p0, p1, p2, p3);
	}
	else
	{
		VECTOR p1 = (tanL * alphaL) + p0;
		VECTOR p2 = (tanR * alphaR) + p3;
		return CubicBezier(p0, p1, p2, p3);
	}
}

/// <summary>
 /// Attempts to find a slightly better parameterization for u on the given curve.
 /// </summary>
void CurveFitBase::Reparameterize(int first, int last, CubicBezier curve)
{
	std::vector<VECTOR>& pts = _pts;
	std::vector<FLOAT>& u = _u;
	int nPts = last - first;
	for (int i = 1; i < nPts; i++)
	{
		VECTOR p = pts[first + i];
		FLOAT t = u[i];
		FLOAT ti = 1 - t;

		// Control vertices for Q'
		VECTOR qp0 = (curve.p1 - curve.p0) * 3.f;
		VECTOR qp1 = (curve.p2 - curve.p1) * 3.f;
		VECTOR qp2 = (curve.p3 - curve.p2) * 3.f;

		// Control vertices for Q''
		VECTOR qpp0 = (qp1 - qp0) * 2.f;
		VECTOR qpp1 = (qp2 - qp1) * 2.f;

		// Evaluate Q(t), Q'(t), and Q''(t)
		VECTOR p0 = curve.Sample(t);
		VECTOR p1 = ((ti * ti) * qp0) + ((2 * ti * t) * qp1) + ((t * t) * qp2);
		VECTOR p2 = (ti * qpp0) + (t * qpp1);

		// these are the actual fitting calculations using http://en.wikipedia.org/wiki/Newton%27s_method
		// We can't just use .X and .Y because Unity uses lower-case "x" and "y".
		FLOAT num = ((VectorHelper::GetX(p0) - VectorHelper::GetX(p)) * VectorHelper::GetX(p1)) + ((VectorHelper::GetY(p0) - VectorHelper::GetY(p)) * VectorHelper::GetY(p1));
		FLOAT den = (VectorHelper::GetX(p1) * VectorHelper::GetX(p1)) + (VectorHelper::GetY(p1) * VectorHelper::GetY(p1)) + ((VectorHelper::GetX(p0) - VectorHelper::GetX(p)) * VectorHelper::GetX(p2)) + ((VectorHelper::GetY(p0) - VectorHelper::GetY(p)) * VectorHelper::GetY(p2));
		FLOAT newU = t - num / den;
		if (std::abs(den) > EPSILON && newU >= 0 && newU <= 1)
			u[i] = newU;
	}
}

/// <summary>
/// Computes the maximum squared distance from a point to the curve using the current parameterization.
/// </summary>
FLOAT CurveFitBase::FindMaxSquaredError(int first, int last, CubicBezier curve, int& split)
{
	std::vector<VECTOR>& pts = _pts;
	std::vector<FLOAT>& u = _u;
	int s = (last - first + 1) / 2;
	int nPts = last - first + 1;
	FLOAT max = 0;
	for (int i = 1; i < nPts; i++)
	{
		VECTOR v0 = pts[first + i];
		VECTOR v1 = curve.Sample(u[i]);
		FLOAT d = VectorHelper::DistanceSquared(v0, v1);
		if (d > max)
		{
			max = d;
			s = i;
		}
	}

	// split at the point of maximum error
	split = s + first;
	if (split <= first)
		split = first + 1;
	if (split >= last)
		split = last - 1;

	return max;
}

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

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

module bezierfit;

import :curve_fit;
import :curve_preprocess;

using namespace bezierfit;

std::vector<VECTOR> bezierfit::reduce(std::vector<VECTOR> points, FLOAT error)
{
	return CurvePreprocess::RdpReduce(points, error);
}

std::pair<VECTOR, VECTOR> bezierfit::calc_four_point_cubic_bezier(const VECTOR& p0, const VECTOR& p1, const VECTOR& p2, const VECTOR& p3)
{
	// See https://apoorvaj.io/cubic-bezier-through-four-points/
	constexpr auto alpha = 0.5f;
	auto d1 = powf(glm::distance(p1, p0), alpha);
	auto d2 = powf(glm::distance(p2, p1), alpha);
	auto d3 = powf(glm::distance(p3, p2), alpha);

	auto a = d1 * d1;
	auto b = d2 * d2;
	auto c = (2.f * d1 * d1) + (3 * d1 * d2) + (d2 * d2);
	auto d = 3.f * d1 * (d1 + d2);
	VECTOR t1{
		(a * p2.x - b * p0.x + c * p1.x) / d,
		(a * p2.y - b * p0.y + c * p1.y) / d
	};

	a = d3 * d3;
	b = d2 * d2;
	c = (2 * d3 * d3) + (3 * d3 * d2) + (d2 * d2);
	d = 3 * d3 * (d3 + d2);
	auto t2 = VECTOR{
		(a * p1.x - b * p3.x + c * p2.x) / d,
		(a * p1.y - b * p3.y + c * p2.y) / d
	};
	return { t1, t2 };
}

std::vector<std::array<VECTOR, 4>> bezierfit::fit(std::vector<VECTOR> data, FLOAT maxError)
{
	if (data.empty())
		return {};
	auto reduced = CurvePreprocess::RdpReduce(data, 0.03f);

	CurveFit curveFit{};
	auto bezierCurves = curveFit.Fit(reduced, maxError);
	std::vector<std::array<VECTOR, 4>> result;
	result.resize(bezierCurves.size());
	for (int i = 0; i < bezierCurves.size(); ++i)
	{
		auto& bc = bezierCurves[i];
		result[i] = { bc.p0,bc.p1,bc.p2,bc.p3 };
	}
	return result;
}

bool CurveFitBase::FitCurve(int first, int last, VECTOR tanL, VECTOR tanR, CubicBezier& curve, int& split)
{
	int nPts = last - first + 1;
	if (nPts < 2)
	{
		throw std::invalid_argument("INTERNAL ERROR: Should always have at least 2 points here");
	}
	else if (nPts == 2)
	{
		// if we only have 2 points left, estimate the curve using Wu/Barsky
		VECTOR p0 = _pts[first];
		VECTOR p3 = _pts[last];
		float alpha = glm::distance(p0, p3) / 3;
		VECTOR p1 = (tanL * alpha) + p0;
		VECTOR p2 = (tanR * alpha) + p3;
		curve = CubicBezier(p0, p1, p2, p3);
		split = 0;
		return true;
	}
	else
	{
		split = 0;
		ArcLengthParamaterize(first, last); // initially start u with a simple chord-length paramaterization
		for (int i = 0; i < MAX_ITERS + 1; i++)
		{
			if (i != 0)
				Reparameterize(first, last, curve); // use Newton's method to find better parameters (except on the first run, since we don't have a curve yet)
			curve = GenerateBezier(first, last, tanL, tanR); // generate the curve itself
			float error = FindMaxSquaredError(first, last, curve, split); // calculate error and get split point (point of max error)
			if (error < _squaredError)
				return true; // if we're within error tolerance, awesome!
		}
		return false;
	}
}

// Initialize the static member variable NO_CURVES.
const std::vector<CubicBezier> CurveFit::NO_CURVES;

std::vector< CubicBezier> CurveFit::Fit(std::vector<VECTOR> points, FLOAT maxError)
{
	if (maxError < EPSILON)
		throw std::invalid_argument("maxError cannot be negative/zero/less than epsilon value");
	if (points.size() < 2)
		return NO_CURVES; // need at least 2 points to do anything

	CurveFit instance;
	instance._pts = points;
	instance.InitializeArcLengths();
	instance._squaredError = maxError * maxError;

	// Find tangents at ends
	int last = points.size() - 1;
	VECTOR tanL = instance.GetLeftTangent(last);
	VECTOR tanR = instance.GetRightTangent(0);

	// do the actual fit
	instance.FitRecursive(0, last, tanL, tanR);
	return instance._result;
}

void CurveFit::FitRecursive(int first, int last, VECTOR tanL, VECTOR tanR)
{
	int split;
	CubicBezier curve;
	if (FitCurve(first, last, tanL, tanR, curve, split))
	{
		_result.push_back(curve);
	}
	else
	{
		// If we get here, fitting failed, so we need to recurse
		// first, get mid tangent
		VECTOR tanM1 = GetCenterTangent(first, last, split);
		VECTOR tanM2 = -tanM1;

		// our end tangents might be based on points outside the new curve (this is possible for mid tangents too
		// but since we need to maintain C1 continuity, it's too late to do anything about it)
		if (first == 0 && split < END_TANGENT_N_PTS)
			tanL = GetLeftTangent(split);
		if (last == _pts.size() - 1 && split > (_pts.size() - (END_TANGENT_N_PTS + 1)))
			tanR = GetRightTangent(split);

		// do actual recursion
		FitRecursive(first, split, tanL, tanM1);
		FitRecursive(split, last, tanM2, tanR);
	}
}


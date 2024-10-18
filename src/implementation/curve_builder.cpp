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
#include <vector>
#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include "glm_wrapper.hpp"

module bezierfit;

import :curve_builder;

using namespace bezierfit;

const CurveBuilder::AddPointResult bezierfit::CurveBuilder::AddPointResult::NO_CHANGE {};

bool  CurveBuilder::AddPointResult::WasChanged() const { return data != 0; }
int CurveBuilder::AddPointResult::FirstChangedIndex() const { return std::abs(data) - 1; }
bool CurveBuilder::AddPointResult::WasAdded() const { return data < 0; }

CurveBuilder::AddPointResult::AddPointResult(int firstChangedIndex, bool curveAdded)
	: data((firstChangedIndex + 1)* (curveAdded ? -1 : 1))
{
	assert(firstChangedIndex >= 0 && firstChangedIndex != std::numeric_limits<int>::max());
}

CurveBuilder::CurveBuilder(FLOAT linDist, FLOAT error)
	:_linDist(linDist), _totalLength(0.0f), _first(0), _tanL(VECTOR{ 0.0f, 0.0f })
{
	_squaredError = error * error;
}

CurveBuilder::AddPointResult CurveBuilder::AddPoint(const VECTOR& p)
{
	VECTOR prev = _prev;
	std::vector<VECTOR>& pts = _pts;
	int count = static_cast<int>(pts.size());
	if (count != 0)
	{
		FLOAT td = VectorHelper::Distance(prev, p);
		FLOAT md = _linDist;
		if (td > md)
		{
			int first = std::numeric_limits<int>::max();
			bool add = false;
			FLOAT rd = td - md;
			VECTOR dir = VectorHelper::Normalize(p - prev);
			do
			{
				VECTOR np = prev + dir * md;
				AddPointResult res = AddInternal(np);
				first = std::min(first, res.FirstChangedIndex());
				add |= res.WasAdded();
				prev = np;
				rd -= md;
			} while (rd > md);
			_prev = prev;
			return AddPointResult(first, add);
		}
		return AddPointResult::NO_CHANGE;
	}
	else
	{
		_prev = p;
		_pts.push_back(p);
		_arclen.push_back(0.0f);
		return AddPointResult::NO_CHANGE;
	}
}

const std::vector<CubicBezier>& CurveBuilder::Curves() const { return _result; }

void CurveBuilder::Clear()
{
	_result.clear();
	_pts.clear();
	_arclen.clear();
	_u.clear();
	_totalLength = 0.0f;
	_first = 0;
	_tanL = VECTOR{ 0.0f, 0.0f };
	_prev = VECTOR{ 0.0f, 0.0f };
}

CurveBuilder::AddPointResult CurveBuilder::AddInternal(const VECTOR& np)
{
	std::vector<VECTOR>& pts = _pts;
	int last = static_cast<int>(pts.size());
	assert(last != 0);

	pts.push_back(np);
	_arclen.push_back(_totalLength = _totalLength + _linDist);

	if (last == 1)
	{
		assert(_result.empty());
		VECTOR p0 = pts[0];
		VECTOR tanL = VectorHelper::Normalize(np - p0);
		VECTOR tanR = -tanL;
		_tanL = tanL;
		FLOAT alpha = _linDist / 3;
		VECTOR p1 = tanL * alpha + p0;
		VECTOR p2 = tanR * alpha + np;
		_result.push_back(CubicBezier(p0, p1, p2, np));
		return AddPointResult(0, true);
	}
	else
	{
		int lastCurve = static_cast<int>(_result.size()) - 1;
		int first = _first;

		VECTOR tanL = lastCurve == 0 ? GetLeftTangent(last) : _tanL;
		VECTOR tanR = GetRightTangent(first);

		// Try fitting with the new point
		int split;
		CubicBezier curve;
		if (FitCurve(first, last, tanL, tanR, curve, split))
		{
			_result[lastCurve] = curve;
			return AddPointResult(lastCurve, false);
		}
		else
		{
			// Need to split
			VECTOR tanM1 = GetCenterTangent(first, last, split);
			VECTOR tanM2 = -tanM1;

			if (first == 0 && split < END_TANGENT_N_PTS)
				tanL = GetLeftTangent(split);

			// Do a final pass on the first half of the curve
			int unused;
			FitCurve(first, split, tanL, tanM1, curve, unused);
			_result[lastCurve] = curve;

			// Prepare to fit the second half
			FitCurve(split, last, tanM2, tanR, curve, unused);
			_result.push_back(curve);
			_first = split;
			_tanL = tanM2;

			return AddPointResult(lastCurve, true);
		}
	}
}

bool CurveBuilder::FitCurve(int first, int last, const VECTOR& tanL, const VECTOR& tanR, CubicBezier& curve, int& split)
{
	std::vector<VECTOR> pts = _pts;
	int nPts = last - first + 1;
	if (nPts < 2)
	{
		throw new std::logic_error("INTERNAL ERROR: Should always have at least 2 points here");
	}
	else if (nPts == 2)
	{
		// if we only have 2 points left, estimate the curve using Wu/Barsky
		VECTOR p0 = pts[first];
		VECTOR p3 = pts[last];
		FLOAT alpha = VectorHelper::Distance(p0, p3) / 3;
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
		curve = CubicBezier{};
		for (int i = 0; i < MAX_ITERS + 1; i++)
		{
			if (i != 0) Reparameterize(first, last, curve);                                  // use newton's method to find better parameters (except on first run, since we don't have a curve yet)
			curve = GenerateBezier(first, last, tanL, tanR);                                // generate the curve itself
			FLOAT error = FindMaxSquaredError(first, last, curve, split);               // calculate error and get split point (point of max error)
			if (error < _squaredError)  return true;                                         // if we're within error tolerance, awesome!
		}
		return false;
	}
}

std::vector<CubicBezier>::const_iterator CurveBuilder::begin() const { return _result.cbegin(); }
std::vector<CubicBezier>::const_iterator CurveBuilder::end() const { return _result.cend(); }

std::vector<CubicBezier>::iterator CurveBuilder::begin() { return _result.begin(); }
std::vector<CubicBezier>::iterator CurveBuilder::end() { return _result.end(); }

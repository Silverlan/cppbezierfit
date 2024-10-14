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

#include <stdexcept>
#include <string>
#include <cassert>
#include "glm_wrapper.hpp"

module bezierfit;

using namespace bezierfit;

Spline::Spline(int samplesPerCurve) : _samplesPerCurve(samplesPerCurve)
{
	if (_samplesPerCurve < MIN_SAMPLES_PER_CURVE || _samplesPerCurve > MAX_SAMPLES_PER_CURVE)
		throw std::invalid_argument("samplesPerCurve must be between " + std::to_string(MIN_SAMPLES_PER_CURVE) + " and " + std::to_string(MAX_SAMPLES_PER_CURVE));
	_curves.resize(16);
	_arclen.resize(16 * samplesPerCurve);
}

Spline::Spline(const std::vector<CubicBezier>& curves, int samplesPerCurve) : _samplesPerCurve(samplesPerCurve)
{
	if (curves.empty())
		throw std::invalid_argument("curves cannot be empty");
	if (_samplesPerCurve < MIN_SAMPLES_PER_CURVE || _samplesPerCurve > MAX_SAMPLES_PER_CURVE)
		throw std::invalid_argument("samplesPerCurve must be between " + std::to_string(MIN_SAMPLES_PER_CURVE) + " and " + std::to_string(MAX_SAMPLES_PER_CURVE));
	_curves.resize(curves.size());
	_arclen.resize(_curves.size() * samplesPerCurve);
	for (auto& curve : curves)
		Add(curve);
}

void Spline::Add(const CubicBezier& curve)
{
	if (_curves.size() > 0 && !VectorHelper::EqualsOrClose(_curves[_curves.size() - 1].p3, curve.p0))
		throw std::invalid_argument("The new curve at index " + std::to_string(_curves.size()) + " does not connect with the previous curve at index " + std::to_string(_curves.size() - 1));
	_curves.push_back(curve);
	for (int i = 0; i < _samplesPerCurve; i++) // expand the array since updateArcLengths expects these values to be there
		_arclen.push_back(0);
	UpdateArcLengths(_curves.size() - 1);
}

void Spline::Update(int index, const CubicBezier& curve)
{
	if (index < 0)
		throw std::out_of_range("Negative index");
	if (index >= _curves.size())
		throw std::out_of_range("Curve index " + std::to_string(index) + " is out of range (there are " + std::to_string(_curves.size()) + " curves in the spline)");
	if (index > 0 && !VectorHelper::EqualsOrClose(_curves[index - 1].p3, curve.p0))
		throw std::invalid_argument("The updated curve at index " + std::to_string(index) + " does not connect with the previous curve at index " + std::to_string(index - 1));
	if (index < _curves.size() - 1 && !VectorHelper::EqualsOrClose(_curves[index + 1].p0, curve.p3))
		throw std::invalid_argument("The updated curve at index " + std::to_string(index) + " does not connect with the next curve at index " + std::to_string(index + 1));

	_curves[index] = curve;
	for (int i = index; i < _curves.size(); i++)
		UpdateArcLengths(i);
}

void Spline::Clear()
{
	_curves.clear();
	_arclen.clear();
}

FLOAT Spline::Length() const
{
	std::vector<FLOAT> arclen = _arclen;
	int count = arclen.size();
	return count == 0 ? 0 : arclen[count - 1];
}

const std::vector<CubicBezier>& Spline::Curves() const
{
	return _curves;
}

glm::vec2 Spline::Sample(FLOAT u) const
{
	SamplePos pos = GetSamplePosition(u);
	return _curves[pos.Index].Sample(pos.Time);
}

typename Spline::SamplePos Spline::GetSamplePosition(FLOAT u) const
{
	if (_curves.empty())
		throw std::invalid_argument("No curves have been added to the spline");
	if (u < 0)
		return SamplePos(0, 0);
	if (u > 1)
		return SamplePos(_curves.size() - 1, 1);

	const std::vector<FLOAT>& arclen = _arclen;
	FLOAT total = Length(); // Assuming Length() method is implemented
	FLOAT target = u * total;
	assert(target >= 0);

	// Binary search to find largest value <= target
	int index = 0;
	int low = 0;
	int high = static_cast<int>(arclen.size()) - 1;
	FLOAT found = std::numeric_limits<FLOAT>::quiet_NaN();
	while (low < high)
	{
		index = (low + high) / 2;
		found = arclen[index];
		if (found < target)
			low = index + 1;
		else
			high = index;
	}

	// this should be a rather rare scenario: we're past the end, but this wasn't picked up by the test for u >= 1
	if (index >= static_cast<int>(arclen.size()) - 1)
		return SamplePos(_curves.size() - 1, 1);

	// this can happen because the binary search can give us either index or index + 1
	if (found > target)
		index--;

	if (index < 0)
	{
		// We're at the beginning of the spline
		FLOAT max = arclen[0];
		assert(target <= max + EPSILON); // Assuming EPSILON is defined
		FLOAT part = target / max;
		FLOAT t = part / _samplesPerCurve;
		return SamplePos(0, t);
	}
	else
	{
		// interpolate between two values to see where the index would be if continuous values
		FLOAT min = arclen[index];
		FLOAT max = arclen[index + 1];
		assert(target >= min - EPSILON && target <= max + EPSILON); // Assuming EPSILON is defined
		FLOAT part = target < min ? 0 : target > max ? 1 : (target - min) / (max - min);
		FLOAT t = (((index + 1) % static_cast<int>(_samplesPerCurve)) + part) / _samplesPerCurve;
		int curveIndex = (index + 1) / static_cast<int>(_samplesPerCurve);
		return SamplePos(curveIndex, t);
	}
}

void Spline::UpdateArcLengths(int iCurve)
{
	assert(iCurve >= 0 && iCurve < _curves.size());

	CubicBezier curve = _curves[iCurve];
	int nSamples = static_cast<int>(_samplesPerCurve);
	std::vector<FLOAT>& arclen = _arclen;
	FLOAT clen = iCurve > 0 ? arclen[iCurve * nSamples - 1] : 0;
	VECTOR pp = curve.Sample(0); // Assuming t = 0 for the starting point
	assert(arclen.size() >= ((iCurve + 1) * nSamples));
	for (int iPoint = 0; iPoint < nSamples; iPoint++)
	{
		int idx = (iCurve * nSamples) + iPoint;
		FLOAT t = static_cast<FLOAT>(iPoint + 1) / static_cast<FLOAT>(nSamples);
		VECTOR np = curve.Sample(t);
		FLOAT d = glm::distance(np, pp);
		clen += d;
		arclen[idx] = clen;
		pp = np;
	}
}

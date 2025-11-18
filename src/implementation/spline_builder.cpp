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

#include <cassert>

module bezierfit;

import :spline_builder;

using namespace bezierfit;

SplineBuilder::SplineBuilder(FLOAT pointDistance, FLOAT error, int samplesPerCurve)
	: _builder(static_cast<float>(pointDistance), static_cast<float>(error)), _spline(samplesPerCurve)
{
}

bool SplineBuilder::Add(const glm::vec2& p)
{
	// Add point to CurveBuilder and check if the spline was modified
	CurveBuilder::AddPointResult res = _builder.AddPoint(p);
	if (!res.WasChanged())
		return false;

	// Update spline
	const std::vector<CubicBezier>& curves = _builder.Curves();
	if (res.WasAdded() && curves.size() == 1)
	{
		// First curve
		assert(_spline.Curves().empty());
		_spline.Add(curves[0]);
	}
	else if (res.WasAdded())
	{
		// Split
		_spline.Update(_spline.Curves().size() - 1, curves[res.FirstChangedIndex()]);
		for (int i = res.FirstChangedIndex() + 1; i < curves.size(); i++)
			_spline.Add(curves[i]);
	}
	else
	{
		// Last curve updated
        assert(res.FirstChangedIndex() == curves.size() - 1);
		_spline.Update(_spline.Curves().size() - 1, curves[curves.size() - 1]);
	}

	return true;
}

glm::vec2 SplineBuilder::Sample(FLOAT u) const
{
	return _spline.Sample(static_cast<float>(u));
}

glm::vec2 SplineBuilder::Tangent(FLOAT u) const
{
	Spline::SamplePos pos = _spline.GetSamplePosition(u);
	return _spline.Curves()[pos.Index].Tangent(pos.Time);
}

void SplineBuilder::Clear()
{
	_builder.Clear();
	_spline.Clear();
}

const std::vector<CubicBezier>& SplineBuilder::Curves() const
{
	return _spline.Curves();
}

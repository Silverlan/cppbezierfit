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

#include <vector>

export module bezierfit:spline_builder;

import :cubic_bezier;
import :curve_builder;
import :spline;

namespace bezierfit {
	class SplineBuilder
	{
	public:
		SplineBuilder(FLOAT pointDistance, FLOAT error, int samplesPerCurve);

		bool Add(const glm::vec2& p);
		glm::vec2 Sample(FLOAT u) const;
		glm::vec2 Tangent(FLOAT u) const;
		void Clear();
		const std::vector<CubicBezier>& Curves() const;

	private:
		CurveBuilder _builder;
		Spline _spline;
	};
};

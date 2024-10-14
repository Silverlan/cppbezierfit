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
#include "glm_wrapper.hpp"

export module bezierfit:spline;

import :cubic_bezier;

namespace bezierfit
{
	class Spline
	{
	public:
		static const int MIN_SAMPLES_PER_CURVE = 8;
		static const int MAX_SAMPLES_PER_CURVE = 1024;
		static const FLOAT EPSILON;

		struct SamplePos
		{
			int Index;
			FLOAT Time;

			SamplePos(int curveIndex, FLOAT t) : Index(curveIndex), Time(t) {}
		};

		Spline(int samplesPerCurve);
		Spline(const std::vector<CubicBezier>& curves, int samplesPerCurve);

		void Add(const CubicBezier& curve);
		void Update(int index, const CubicBezier& curve);
		void Clear();
		FLOAT Length() const;
		const std::vector<CubicBezier>& Curves() const;
		glm::vec2 Sample(FLOAT u) const;
		SamplePos GetSamplePosition(FLOAT u) const;

	private:
		void UpdateArcLengths(int iCurve);

		std::vector<CubicBezier> _curves;
		std::vector<FLOAT> _arclen;
		int _samplesPerCurve;
	};
}

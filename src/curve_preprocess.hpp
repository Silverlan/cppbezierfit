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

#ifndef __BEZIERFIT_CURVE_PREPROCESS_HPP__
#define __BEZIERFIT_CURVE_PREPROCESS_HPP__

#include "bezier_fit.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <glm/gtc/epsilon.hpp>

namespace bezierfit
{
	class CurvePreprocess
	{
	public:
		static constexpr FLOAT EPSILON = 0.000001; // Change the epsilon value as needed for FLOAT type

		static std::vector<glm::vec2> Linearize(const std::vector<glm::vec2>& src, FLOAT md);

		static std::vector<glm::vec2> RemoveDuplicates(const std::vector<glm::vec2>& pts);

		static std::vector<glm::vec2> RdpReduce(const std::vector<glm::vec2>& pts, FLOAT error);

	private:
		static void RdpRecursive(const std::vector<glm::vec2>& pts, FLOAT error, int first, int last, std::vector<int>& keepIndex);

		static FLOAT PerpendicularDistance(const glm::vec2& a, const glm::vec2& b, FLOAT abDist, FLOAT aCrossB, const glm::vec2& p);
	};
}

#endif

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

#ifndef __BEZIERFIT_VECTOR_HELPER_HPP__
#define __BEZIERFIT_VECTOR_HELPER_HPP__

#include "bezier_fit.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

namespace bezierfit {
	class VectorHelper
	{
	public:
		static const FLOAT EPSILON;
		static FLOAT Distance(const glm::vec2& a, const glm::vec2& b);
		static FLOAT DistanceSquared(const glm::vec2& a, const glm::vec2& b);
		static FLOAT Dot(const glm::vec2& a, const glm::vec2& b);
		static glm::vec2 Normalize(const glm::vec2& v);
		static FLOAT Length(const glm::vec2& v);
		static FLOAT LengthSquared(const glm::vec2& v);
		static glm::vec2 Lerp(const glm::vec2& a, const glm::vec2& b, FLOAT amount);
		static FLOAT GetX(const glm::vec2& v);
		static FLOAT GetY(const glm::vec2& v);
		static bool EqualsOrClose(const glm::vec2& v1, const glm::vec2& v2);
	};
};

#endif

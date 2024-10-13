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

module bezierfit;

import "interface/glm_wrapper.hpp";

using namespace bezierfit;

const FLOAT VectorHelper::EPSILON = 1.2e-12f;

FLOAT VectorHelper::Distance(const glm::vec2& a, const glm::vec2& b)
{
	return glm::distance(a, b);
}

FLOAT VectorHelper::DistanceSquared(const glm::vec2& a, const glm::vec2& b)
{
	return glm::distance2(a, b);
}

FLOAT VectorHelper::Dot(const glm::vec2& a, const glm::vec2& b)
{
	return glm::dot(a, b);
}

glm::vec2 VectorHelper::Normalize(const glm::vec2& v)
{
	return glm::normalize(v);
}

FLOAT VectorHelper::Length(const glm::vec2& v)
{
	return glm::length(v);
}

FLOAT VectorHelper::LengthSquared(const glm::vec2& v)
{
	return glm::length2(v);
}

glm::vec2 VectorHelper::Lerp(const glm::vec2& a, const glm::vec2& b, FLOAT amount)
{
	return glm::mix(a, b, amount);
}

FLOAT VectorHelper::GetX(const glm::vec2& v)
{
	return v.x;
}

FLOAT VectorHelper::GetY(const glm::vec2& v)
{
	return v.y;
}

bool VectorHelper::EqualsOrClose(const glm::vec2& v1, const glm::vec2& v2)
{
	return DistanceSquared(v1, v2) < EPSILON;
}


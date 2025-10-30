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


export module bezierfit:cubic_bezier;

import :vector_helper;

namespace bezierfit
{
	class CubicBezier
	{
	public:
		// Control points
		VECTOR p0;
		VECTOR p1;
		VECTOR p2;
		VECTOR p3;

		CubicBezier() = default;
		CubicBezier(const VECTOR& p0, const VECTOR& p1, const VECTOR& p2, const VECTOR& p3);
		CubicBezier& operator=(const CubicBezier& other) {
			p0 = other.p0;
			p1 = other.p1;
			p2 = other.p2;
			p3 = other.p3;
			return *this;
		}

		VECTOR Sample(FLOAT t) const;

		VECTOR Derivative(FLOAT t) const;

		VECTOR Tangent(FLOAT t) const;

		std::string ToString() const;

		// Equality members
		bool operator==(const CubicBezier& other) const;

		bool operator!=(const CubicBezier& other) const;
	};
}

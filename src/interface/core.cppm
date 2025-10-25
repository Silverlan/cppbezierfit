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
#include <array>

export module bezierfit:core;

export import glm;

export namespace bezierfit {
	using VECTOR = glm::vec2;
	using FLOAT = float;

	std::vector<VECTOR> reduce(std::vector<VECTOR> points, FLOAT error = 0.03f);
	std::vector<std::array<VECTOR, 4>> fit(std::vector<VECTOR> points, FLOAT maxError);
	std::pair<VECTOR, VECTOR> calc_four_point_cubic_bezier(const VECTOR &v0, const VECTOR &v1, const VECTOR &v2, const VECTOR &v3);
};

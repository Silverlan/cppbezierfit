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

#ifndef __BEZIER_FIT_HPP__
#define __BEZIER_FIT_HPP__

#include <glm/glm.hpp>
#include <glm/vec2.hpp>
#include <vector>
#include <array>

namespace bezierfit {
	using VECTOR = glm::vec2;
	using FLOAT = float;

	std::vector<VECTOR> reduce(std::vector<VECTOR> points, FLOAT error = 0.03f);
	std::vector<std::array<VECTOR, 4>> fit(std::vector<VECTOR> points, FLOAT maxError);
};

#endif

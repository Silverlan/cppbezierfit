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


export module bezierfit:curve_preprocess;

import :core;

namespace bezierfit
{
	class CurvePreprocess
	{
	public:
		static constexpr FLOAT EPSILON = 0.000001; // Change the epsilon value as needed for FLOAT type

		static std::vector<VECTOR> Linearize(const std::vector<VECTOR>& src, FLOAT md);

		static std::vector<VECTOR> RemoveDuplicates(const std::vector<VECTOR>& pts);

		static std::vector<VECTOR> RdpReduce(const std::vector<VECTOR>& pointList, float epsilon);

	private:
		static FLOAT PerpendicularDistance(const VECTOR& p, const VECTOR& lineP1, const VECTOR& lineP2);
	};
}

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

#include <stdexcept>
#include <vector>
#include <algorithm>
#include "glm_wrapper.hpp"

module bezierfit;

import :curve_preprocess;

using namespace bezierfit;
std::vector<VECTOR> CurvePreprocess::Linearize(const std::vector<VECTOR>& src, FLOAT md)
{
	if (src.empty())
		throw std::invalid_argument("src cannot be empty");
	if (md <= EPSILON)
		throw std::invalid_argument("md must be greater than epsilon");

	std::vector<VECTOR> dst;
	if (src.size() > 0)
	{
		VECTOR pp = src[0];
		dst.push_back(pp);
		FLOAT cd = 0;
		for (size_t ip = 1; ip < src.size(); ip++)
		{
			VECTOR p0 = src[ip - 1];
			VECTOR p1 = src[ip];
			FLOAT td = glm::distance(p0, p1);
			if (cd + td > md)
			{
				FLOAT pd = md - cd;
				dst.push_back(glm::mix(p0, p1, pd / td));
				FLOAT rd = td - pd;
				while (rd > md)
				{
					rd -= md;
					VECTOR np = glm::mix(p0, p1, (td - rd) / td);
					if (!glm::all(glm::epsilonEqual(np, pp, EPSILON)))
					{
						dst.push_back(np);
						pp = np;
					}
				}
				cd = rd;
			}
			else
			{
				cd += td;
			}
		}
		// last point
		VECTOR lp = src.back();
		if (!glm::all(glm::epsilonEqual(pp, lp, EPSILON)))
			dst.push_back(lp);
	}
	return dst;
}

std::vector<VECTOR> CurvePreprocess::RemoveDuplicates(const std::vector<VECTOR>& pts)
{
	if (pts.size() < 2)
		return pts;

	std::vector<VECTOR> dst;
	dst.reserve(pts.size());
	dst.push_back(pts[0]);
	for (size_t i = 1; i < pts.size(); i++)
	{
		VECTOR cur = pts[i];
		VECTOR prev = dst.back();
		if (!glm::all(glm::epsilonEqual(prev, cur, EPSILON)))
		{
			dst.push_back(cur);
		}
	}
	return dst;
}

std::vector<VECTOR> CurvePreprocess::RdpReduce(const std::vector<VECTOR>& pointList, float epsilon)
{
	std::vector<VECTOR> resultList;
	resultList.reserve(pointList.size() /2);

	// Find the point with the maximum distance
	float dmax = 0;
	int index = 0;
	for (int i = 1; i < pointList.size() - 1; ++i) {
		float d = PerpendicularDistance(pointList[i], pointList[0], pointList[pointList.size() - 1]);
		if (d > dmax) {
			index = i;
			dmax = d;
		}
	}
	// If max distance is greater than epsilon, recursively simplify
	if (dmax > epsilon) {
		std::vector<VECTOR> pre_part, next_part;
		pre_part.reserve(index +1);
		for (int i = 0; i <= index; ++i)
			pre_part.push_back(pointList[i]);
		next_part.reserve(pointList.size() -index);
		for (int i = index; i < pointList.size(); ++i)
			next_part.push_back(pointList[i]);
		// Recursive call
		std::vector<VECTOR> resultList1 = RdpReduce(pre_part, epsilon);
		std::vector<VECTOR> resultList2 = RdpReduce(next_part, epsilon);

		// combine
		resultList.insert(resultList.end(), resultList1.begin(), resultList1.end());
		resultList.insert(resultList.end(), resultList2.begin() + 1, resultList2.end());
	}
	else {
		resultList.push_back(pointList[0]);
		resultList.push_back(pointList[pointList.size() - 1]);
	}
	if (resultList.size() == resultList.capacity())
		resultList.reserve(pointList.size());

	return resultList;
}

FLOAT CurvePreprocess::PerpendicularDistance(const VECTOR& p, const VECTOR& lineP1, const VECTOR& lineP2)
{
	VECTOR vec1 = VECTOR(p.x - lineP1.x, p.y - lineP1.y);
	VECTOR vec2 = VECTOR(lineP2.x - lineP1.x, lineP2.y - lineP1.y);
	float d_vec2 = sqrt(vec2.x * vec2.x + vec2.y * vec2.y);
	float cross_product = vec1.x * vec2.y - vec2.x * vec1.y;
	float d = fabs(cross_product / d_vec2);
	return d;
}

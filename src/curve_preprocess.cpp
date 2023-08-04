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

#include "curve_preprocess.hpp"
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <glm/gtc/epsilon.hpp>

using namespace bezierfit;
std::vector<glm::vec2> CurvePreprocess::Linearize(const std::vector<glm::vec2>& src, FLOAT md)
{
	if (src.empty())
		throw std::invalid_argument("src cannot be empty");
	if (md <= EPSILON)
		throw std::invalid_argument("md must be greater than epsilon");

	std::vector<glm::vec2> dst;
	if (src.size() > 0)
	{
		glm::vec2 pp = src[0];
		dst.push_back(pp);
		FLOAT cd = 0;
		for (size_t ip = 1; ip < src.size(); ip++)
		{
			glm::vec2 p0 = src[ip - 1];
			glm::vec2 p1 = src[ip];
			FLOAT td = glm::distance(p0, p1);
			if (cd + td > md)
			{
				FLOAT pd = md - cd;
				dst.push_back(glm::mix(p0, p1, pd / td));
				FLOAT rd = td - pd;
				while (rd > md)
				{
					rd -= md;
					glm::vec2 np = glm::mix(p0, p1, (td - rd) / td);
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
		glm::vec2 lp = src.back();
		if (!glm::all(glm::epsilonEqual(pp, lp, EPSILON)))
			dst.push_back(lp);
	}
	return dst;
}

std::vector<glm::vec2> CurvePreprocess::RemoveDuplicates(const std::vector<glm::vec2>& pts)
{
	if (pts.size() < 2)
		return pts;

	std::vector<glm::vec2> dst;
	dst.reserve(pts.size());
	dst.push_back(pts[0]);
	for (size_t i = 1; i < pts.size(); i++)
	{
		glm::vec2 cur = pts[i];
		glm::vec2 prev = dst.back();
		if (!glm::all(glm::epsilonEqual(prev, cur, EPSILON)))
		{
			dst.push_back(cur);
		}
	}
	return dst;
}

std::vector<glm::vec2> CurvePreprocess::RdpReduce(const std::vector<glm::vec2>& pts, FLOAT error)
{
	if (pts.empty())
		throw std::invalid_argument("pts cannot be empty");
	std::vector<glm::vec2> uniquePts = RemoveDuplicates(pts);
	if (uniquePts.size() < 3)
		return uniquePts;

	std::vector<int> keepIndex;
	keepIndex.reserve(std::max(uniquePts.size() / 2, static_cast<size_t>(16)));
	keepIndex.push_back(0);
	keepIndex.push_back(uniquePts.size() - 1);
	RdpRecursive(uniquePts, error, 0, uniquePts.size() - 1, keepIndex);
	std::sort(keepIndex.begin(), keepIndex.end());

	std::vector<glm::vec2> res;
	res.reserve(keepIndex.size());
	for (int idx : keepIndex)
		res.push_back(uniquePts[idx]);
	return res;
}

void CurvePreprocess::RdpRecursive(const std::vector<glm::vec2>& pts, FLOAT error, int first, int last, std::vector<int>& keepIndex)
{
	int nPts = last - first + 1;
	if (nPts < 3)
		return;

	glm::vec2 a = pts[first];
	glm::vec2 b = pts[last];
	FLOAT abDist = glm::distance(a, b);
	FLOAT aCrossB = glm::cross(glm::vec3(a, 0), glm::vec3(b, 0)).z;
	FLOAT maxDist = error;
	int split = 0;
	for (int i = first + 1; i < last - 1; i++)
	{
		glm::vec2 p = pts[i];
		FLOAT pDist = PerpendicularDistance(a, b, abDist, aCrossB, p);
		if (pDist > maxDist)
		{
			maxDist = pDist;
			split = i;
		}
	}

	if (split != 0)
	{
		keepIndex.push_back(split);
		RdpRecursive(pts, error, first, split, keepIndex);
		RdpRecursive(pts, error, split, last, keepIndex);
	}
}

FLOAT CurvePreprocess::PerpendicularDistance(const glm::vec2& a, const glm::vec2& b, FLOAT abDist, FLOAT aCrossB, const glm::vec2& p)
{
	FLOAT area = std::abs(aCrossB +
		a.x * b.y + p.x * a.y -
		p.x * b.y - a.x * p.y);
	FLOAT height = area / abDist;
	return height;
}

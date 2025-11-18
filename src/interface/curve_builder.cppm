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


export module bezierfit:curve_builder;

import :curve_fit;

namespace bezierfit
{
	class CurveBuilder : public CurveFitBase
	{
	public:
		struct AddPointResult
		{
			static const AddPointResult NO_CHANGE;

			bool WasChanged() const;
			int FirstChangedIndex() const;
			bool WasAdded() const;

			AddPointResult() = default;
			AddPointResult(int firstChangedIndex, bool curveAdded);

		private:
			int data = 0;
		};

		explicit CurveBuilder(FLOAT linDist, FLOAT error);

		AddPointResult AddPoint(const VECTOR& p);

		const std::vector<CubicBezier>& Curves() const;

		void Clear();

	private:
		FLOAT _linDist;
		VECTOR _prev;
		VECTOR _tanL;
		FLOAT _totalLength;
		int _first;
		std::vector<CubicBezier> _result;

		AddPointResult AddInternal(const VECTOR& np);

		bool FitCurve(int first, int last, const VECTOR& tanL, const VECTOR& tanR, CubicBezier& curve, int& split);

		std::vector<CubicBezier>::const_iterator begin() const;
		std::vector<CubicBezier>::const_iterator end() const;

		std::vector<CubicBezier>::iterator begin();
		std::vector<CubicBezier>::iterator end();
	};
};

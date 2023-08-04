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

#ifndef __BEZIERFIT_CURVE_FIT_HPP__
#define __BEZIERFIT_CURVE_FIT_HPP__

#include "cubic_bezier.hpp"
#include <vector>

namespace bezierfit {
	const FLOAT EPSILON = std::numeric_limits<FLOAT>::epsilon();
	const int MAX_ITERS = 4;
	const int END_TANGENT_N_PTS = 8;
	const int MID_TANGENT_N_PTS = 4;

	class CurveFitBase
	{
	protected:
		std::vector<VECTOR> _pts;
		std::vector<FLOAT> _arclen;
		std::vector<FLOAT> _u;
		FLOAT _squaredError;

		VECTOR GetLeftTangent(int last);

		VECTOR GetRightTangent(int first);

		VECTOR GetCenterTangent(int first, int last, int split);

		void InitializeArcLengths();

		void ArcLengthParamaterize(int first, int last);

		/// <summary>
		 /// Generates a bezier curve for the segment using a least-squares approximation.
		 /// </summary>
		CubicBezier GenerateBezier(int first, int last, VECTOR tanL, VECTOR tanR);

		/// <summary>
		 /// Attempts to find a slightly better parameterization for u on the given curve.
		 /// </summary>
		void Reparameterize(int first, int last, CubicBezier curve);

		/// <summary>
		/// Computes the maximum squared distance from a point to the curve using the current parameterization.
		/// </summary>
		FLOAT FindMaxSquaredError(int first, int last, CubicBezier curve, int& split);

		/// <summary>
		/// Tries to fit single Bezier curve to the points in [first ... last]. Destroys anything in <see cref="_u"/> in the process.
		/// Assumes there are at least two points to fit.
		/// </summary>
		/// <param name="first">Index of first point to consider.</param>
		/// <param name="last">Index of last point to consider (inclusive).</param>
		/// <param name="tanL">Tangent at the start of the curve ("left").</param>
		/// <param name="tanR">Tangent on the end of the curve ("right").</param>
		/// <param name="curve">The fitted curve.</param>
		/// <param name="split">Point at which to split if this method returns false.</param>
		/// <returns>true if the fit was within error tolerance, false if the curve should be split. Even if this returns false, curve will contain
		/// a curve that somewhat fits the points; it's just outside error tolerance.</returns>
		bool FitCurve(int first, int last, VECTOR tanL, VECTOR tanR, CubicBezier& curve, int& split);

	};

	class CurveFit : public CurveFitBase
	{
	public:
		std::vector< CubicBezier> Fit(std::vector<VECTOR> points, FLOAT maxError);
	private:
		// Curves we've found so far.
		std::vector<CubicBezier> _result;

		// Shared zero-curve array.
		static const std::vector<CubicBezier> NO_CURVES;

		/// <summary>
		/// Main fit function that attempts to fit a segment of curve and recurses if unable to.
		/// </summary>
		void FitRecursive(int first, int last, VECTOR tanL, VECTOR tanR);

		// Other functions and variables go here...
	};
};

#endif

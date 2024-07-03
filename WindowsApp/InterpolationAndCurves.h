#pragma once

#include <cmath>
#include "LinearAlgebra.h"
#include "Calculus.h"
#include "ControlPoints.h"


template<typename V>
static inline V TLinearInterpolate(float t, V a, V b, float t0 = 0, float t1 = 1)
{
	/*if (fabs(a - b) < FLT_EPSILON)
	{
		return a;
	}*/

	//map t to [0,1]
	t = (t - t0) / (t1 - t0);

	return a * (1 - t) + b * t;
}

template<typename V>
static inline V TBarycentricInterpolate(V v0, V v1, V v2, float3& b)
{
	// Use the formula v = v0 * s + v1 * t + v2 * (1 - s - t)
	return v0 * b.x + v1 * b.y + v2 * b.z;
}

static inline float ConvertParameterToSmoothstep(float t)
{
	return 3 * pow(t, 2) - 2 * pow(t, 3);
}


static inline float CircularFadeIn(float x) {
	float y = 1 - sqrt(1 - x * x);
	return y;
}


static inline float CircularFadeOut(float x) {
	float y = sqrt(1 - pow(1 - x, 2));
	return y;
}


///Custom Fade In/Out through a point (x,y) in [0,1] range
static inline float QuadraticThroughAGivenPoint(float x, float2 point)
{
	float A = (1 - point.y) / (1 - point.x) - (point.y / point.x);
	float B = (A * (point.x * point.x) - point.y) / point.x;
	float y = A * (x * x) - B * (x);
	return y;
}

///Smooth Nearest-Neighbor
static inline float DoubleCircleSigmoidFade(float x, float a = 0.5f)
{
	float y = 0;
	if (x <= a) {
		y = a - sqrt(a * a - x * x);
	}
	else {
		y = a + sqrt(pow(1 - a, 2) - pow(x - 1, 2));
	}
	return y;
}

/*
* * * * * * * * * * * * * * * * * * * * * * * Oscillation Interpolation* * * * * * * * * * * * * * * * * * * * * * * * * 
*/

template<typename V>
static inline V TQuadraticInterpolate(float x, V y0, V y1, V y2, float x0 = 0, float x1 = 0.5f, float x2 = 1.f)
{
	return TLinearInterpolate(x, TLinearInterpolate(x, y0, y1, x0, x1), TLinearInterpolate(x, y1, y2, x1, x2), x0, x2);
}

template<typename V>
static inline V TCubicInterpolate(float t, V a, V b, V c, V d, float t0 = 0, float t1 = 0.33f, float t2 = 0.66f, float t3 = 1.f)
{
	auto quad1 = TQuadraticInterpolate(t, a, b, c, t0, t1, t2);
	auto quad2 = TQuadraticInterpolate(t, b, c, d, t1, t2, t3);
	return TLinearInterpolate(t, quad1, quad2, t0, t3);
}

template<typename V>
V TAitkenInterpolate(RDynamicVector<V> outputSet, const RDynamicVector<float>& inputSet, float t)
{
	int n = outputSet.size();
	int j = 0;
	while (n > 1)
	{
		n--;
		j++;
		for (int i = 0; i < n; i++)
		{
			V a = outputSet[i];
			V b = outputSet[i + 1];
			float t0 = inputSet[i];
			float t1 = inputSet[i + j];
			outputSet[i] = TLinearInterpolate(t, a, b, t0, t1);
		}
	}
	return outputSet[0];
}


//Given a set of inputs, current Input-Output index, and newInput, generate Basis (weight) function
inline float GetLagrangeBasis(const RDynamicVector<float>& inputSet, float i, float t)
{
	int n = inputSet.size();
	const float ti = inputSet[i];

	float res = 1;//Basis Polynomial

	//Compute Basis Polynomial
	for (int j = 0; j < n; j++)
	{
		if (i != j)
		{
			const float tj = inputSet[j];
			//basis polynomial: 
			//returns value closer to 1 when t is closer to ti
			//and value closer to 0 when t is closer to tj - disabling the node completely
			float basisPolynomialPart = (t - tj) / (ti - tj);
			res *= basisPolynomialPart;
		}
	}
	return res;
}
/*Oscillating function Interpolation/Approximation*/
template <typename V>
inline V PolynomialFunctionInterpolationLagrange(float newInput, const RDynamicVector<V>& outputSet, const RDynamicVector<float>& inputSet)
{
	//Sum of all terms of created n-degree polynomial
	V newOutput;
	std::memset(&newOutput, 0, sizeof(V));

	int n = outputSet.size();
	for (int i = 0; i < n; i++)
	{
		//Cur output weight
		float wi = GetLagrangeBasis(inputSet, i, newInput);//Get weight of output i for current t

		//Scale cur output
		V yi = outputSet[i];
		newOutput = newOutput + (yi * wi);
	}
	return newOutput;
}

/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * Bezier Interpolation* * * * * * * * * * * * * * * * * * * * * 
*/
template<typename V>
static inline V TQuadraticInterpolateBezier(float x, V y0, V y1, V y2, float x0 = 0, float x1 = 1.f)
{
	return TLinearInterpolate(x, TLinearInterpolate(x, y0, y1, x0, x1), TLinearInterpolate(x, y1, y2, x0, x1), x0, x1);
}

template<typename V>
static inline V TCubicInterpolateBezier(float t, V a, V b, V c, V d, float t0 = 0, float t1 = 1.f)
{
	auto quad1 = TQuadraticInterpolateBezier(t, a, b, c, t0, t1);
	auto quad2 = TQuadraticInterpolateBezier(t, b, c, d, t0, t1);
	return TLinearInterpolate(t, quad1, quad2, t0, t1);
}


/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * Curves* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*/


/*
* Polynomial Parametric
*/
///Specify Coefficients to describe a curve
static inline float GetPointOnPolynomialCurveParametric(float4 coefficients, float t)
{
	//c0 + c1 * t + c2 * pow(t, 2) + c3 * pow(t, 3);
	const float4 Parameters{ 1, t, (float)pow(t, 2), (float)pow(t, 3) };
	return VectorDot(coefficients, Parameters);
}

static inline float2 GetPointOnPolynomialCurveParametric(float4 xCoefficients, float4 yCoefficients, float t)
{
	const float4 Parameters{ 1, t, (float)pow(t, 2), (float)pow(t, 3) };
	return float2{ VectorDot(xCoefficients, Parameters), VectorDot(yCoefficients, Parameters) };
}



/*
* Hermite
*/

static inline float3 GetPointOnHermiteCurveParametric(float t, c_float3& p0, c_float3& v0, c_float3& p1, c_float3& v1)
{
	//Hermit Basis
	float H0 = 1 - 3 * pow(t, 2) + 2 * pow(t, 3);
	float H1 = t - 2 * pow(t, 2) + pow(t, 3);
	float H2 = -pow(t, 2) + pow(t, 3);
	float H3 = 3 * pow(t, 2) - 2 * pow(t, 3);

	return VectorMul(p0, H0) + VectorMul(v0, H1) + VectorMul(v1, H2) + VectorMul(p1, H3);
}

template <typename V>
static inline V TGetPointOnHermiteCurveParametric(float t, const V& y0, const V& d0, const V& y1, const V& d1)
{
	//Hermit Basis
	float H0 = 1 - 3 * pow(t, 2) + 2 * pow(t, 3);
	float H1 = t - 2 * pow(t, 2) + pow(t, 3);
	float H2 = 3 * pow(t, 2) - 2 * pow(t, 3);
	float H3 = -pow(t, 2) + pow(t, 3); 

	return (y0 * H0) + (d0 * H1) + (y1 * H2) + (d1 * H3);
}

template <typename V>
static inline V TGetPointOnQuinticHermiteCurveParametric(float t, const V& p0, const V& d0, const V& d20, const V& p1, const V& d1, const V& d21)
{
	//this is a logical extension of the cubic hermite spline's basis functions
	//that has one basis function for t0 position, one for t1 position
	//one for t0 tangent (1st derivative of position), and one for t1 tangent
	//this adds 2 more basis functions, one for t0 curvature (2nd derivative) and t1 curvature
	auto basis00 = (1 - t) * (1 - t) * (1 - t) * (t * (6 * t + 3) + 1);
	auto basis10 = t * (1 - t) * (1 - t) * (1 - t) * (3 * t + 1);
	auto basis20 = 0.5f * (1 - t) * (1 - t) * (1 - t) * t * t;
	auto basis01 = t * t * t * (t * (6 * t - 15) + 10);
	auto basis11 = t * t * t * (1 - t) * (t * 3 - 4);
	auto basis21 = 0.5f * (1 - t) * (1 - t) * t * t * t;

	return
		p0 * basis00 +
		d0 * basis10 +
		d20 * basis20 +

		d21 * basis21 +
		d1 * basis11 +
		p1 * basis01     ;
}


static inline void ExtractSubsectionOfHermiteCurve(float3& p0, float3& v0, float3& p1, float3& v1, float sStart, float sEnd)
{
	auto InterpolatePositionVelocityHermite = [](float3 p0, float3 v0, float3 p1, float3 v1, float t, float3& outp, float3& outv)
	{
		auto hermite = [&](float x)->float3
		{
			return GetPointOnHermiteCurveParametric(x, p0, v0, p1, v1);
		};

		//Hermit Basis
		outp = hermite(t);
		outv = TGetDerivative<float3>(hermite, t, 0.00001f);
	};

	float3 newStartP;
	float3 newStartV;
	InterpolatePositionVelocityHermite(p0, v0, p1, v1, sStart, newStartP, newStartV);

	float3 newEndP;
	float3 newEndV;
	InterpolatePositionVelocityHermite(p0, v0, p1, v1, sEnd, newEndP, newEndV);
	
	float newLength = sEnd - sStart;

	p0 = newStartP;
	v0 = newStartV * newLength;
	p1 = newEndP;
	v1 = newEndV * newLength;

}

/*
* Bezier
*/

template<typename V>
static inline V TCasteljauInterpolate(const RDynamicVector<V>& controlSet, float t, float t0 = 0, float t1 = 1)
{
	static RDynamicVector<V> tempStorage;
	tempStorage.resize(controlSet.size());
	tempStorage = controlSet;
	int n = controlSet.size();
	while (n > 1)
	{
		n--;
		for (int i = 0; i < n; i++)
		{
			V p1 = tempStorage[i];
			V p2 = tempStorage[i + 1];
			tempStorage[i] = TLinearInterpolate(t, p1, p2, t0, t1);
		}
	}

	return tempStorage[0];

}


static inline void ConvertCubicBezierToHermite(float3& b0, float3& b1, float3& b2, float3& b3)
{
	float3 p0 = b0;
	float3 v0 = (b1 - b0) * 3;
	float3 v1 = (b3 - b2) * 3;
	float3 p1 = b3;

	b0 = p0;
	b1 = v0;
	b2 = p1;
	b3 = v1;
}

static inline void ConvertHermiteToCubicBezier(float3& p0, float3& v0, float3& p1, float3& v1)
{
	float3 b0 = p0;
	float3 b1 = p0 + (v0 / 3.0f);
	float3 b2 = p1 - (v1 / 3.0f);
	float3 b3 = p1;

	p0 = b0;
	v0 = b1;
	p1 = b2;
	v1 = b3;
}

static inline void ExtractSubsectionOfBezierCurve(float3& b0, float3& b1, float3& b2, float3& b3, float sStart, float sEnd)
{
	//Convert to Hermite
	ConvertCubicBezierToHermite(b0, b1, b2, b3);
	//Extract Hermite
	ExtractSubsectionOfHermiteCurve(b0, b1, b2, b3, sStart, sEnd);
	//Convert back to Bezier
	ConvertHermiteToCubicBezier(b0, b1, b2, b3);
}

///************************************************ Cubic Interpolation ****************************************************

static inline float3 InterpolateBSplineSegment(float t, c_float3& p0, c_float3& p1, c_float3& p2, c_float3& p3)
{
	return (
		p0 * ((1 - t) * (1 - t) * (1 - t)) +
		p1 * (t * t * 3 * (t - 2) + 4) +
		p2 * (t * (t * (-3 * t + 3) + 3) + 1) +
		p3 * (t * t * t)
		) / 6.f;
}

static inline float3 GetPointOnBezierCurveParamteric(float t, c_float3& b0, c_float3& b1, c_float3& b2, c_float3& b3)
{
	//Bernstein Basis
	float B0 = pow(1 - t, 3);
	float B1 = 3 * t * pow(1 - t, 2);
	float B2 = 3 * t * t * (1 - t);
	float B3 = pow(t, 3);

	return VectorMul(b0, B0) + VectorMul(b1, B1) + VectorMul(b2, B2) + VectorMul(b3, B3);
}

template <typename V>
static inline V TGetPointOnCatmullRomCurveParametric(float t, const V& pPrev, const V& p0, const V& p1, const V& pNext)
{
	//Catmull Rom Basis
	float H0 = 0 - 1 * t + 2 * pow(t, 2) - 1 * pow(t, 3);
	float H1 = 2 + 0 * t - 5 * pow(t, 2) + 3 * pow(t, 3);
	float H2 = 0 + 1 * t + 4 * pow(t, 2) - 3 * pow(t, 3);
	float H3 = 1 - (H0 + H1 + H2);

	return (pPrev * H0) + (p0 * H1) + (p1 * H2) + (pNext * H3);
}


///************************************************ Curve Designers ****************************************************

class RLinearCurve
{

};

class RCurveDesigner
{
public:
	enum class CurveType
	{
		Bezier,
		Lagrange
	};

	RCurveDesigner()
	{
	}

	void RenderUI()
	{
		
	}

	void Generate()
	{
		OutputSet.clear();


		check(false);///TODO:Feed OutputSet into Control Points, not vice versa
		/*for (auto& point : RSpatialControlSet::Get().ControlPoints)
		{
			OutputSet.push_back(point.Position);
		}*/

		//Generate Input Set
		float singleSegmentSize = 1.f / (OutputSet.size() - 1);
		for (int i = 0; i < OutputSet.size(); i++)
		{
			InputSet.push_back(singleSegmentSize * i);
		}
	}

	float3 GetPointOnACurve(float t, bool bBezier = false)
	{
		///Casteljau
		if (bBezier)
		{
			return TCasteljauInterpolate(OutputSet, t);
		}

		///Aitken
		return TAitkenInterpolate(OutputSet, InputSet, t);
	}


	RDynamicVector<float3> OutputSet;

	RDynamicVector<float> InputSet;


};


///Spline is a piecewise function
///In Uniform spline we don't modify Input t value
///In Non-Uniform Spline we can modify Input & Weights of the Basis Functions
struct RSplineDesigner
{
	RDynamicVector<float3> ControlPoints;
	RDynamicVector<float> InputSet;
	RDynamicVector<float3> DerivativeSet;
	RDynamicVector<float3> SecondDerivativeSet;

	virtual void Generate(float alpha = 0.f) {}
	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal) = 0;
	virtual float3 GetPointOnSpline(float tGlobal)
	{
		uint segmentIndex;
		float tLocal;
		ComputeSegmentIndexAndTLocal(tGlobal, InputSet.empty(), segmentIndex, tLocal);

		return GetPointOnSpline(segmentIndex, tLocal);
	}
	float3 GetPointOnSplineNormalized(float t)
	{
		t = MapToRange(t, 0, 1, GetTValueFirst(), GetTValueLast());
		return GetPointOnSpline(t);
	}

	virtual void RenderUI()
	{}

	virtual uint GetNumSegments()
	{
		return ControlPoints.size() - 1;
	}


	///******************* Non Uniform Spline **********************

	///Non-Uniform Splines should have a generated set of t values per knot
	void ComputeSegmentIndexAndTLocal(float tGlobal, bool bUniformSpline, uint& outSegmentIndex, float& outTLocal)
	{
		if (bUniformSpline)
		{
			outSegmentIndex = int(tGlobal);
			outTLocal = tGlobal - outSegmentIndex;
		}
		else
		{
			outSegmentIndex = GetNonUniformSegmentIndex(tGlobal);

			float t0 = InputSet[outSegmentIndex];
			float t1 = InputSet[outSegmentIndex + 1];
			///Map tGlobal to [0,1]
			outTLocal = (tGlobal - t0) / (t1 - t0);
		}

		check(outSegmentIndex < GetNumSegments());
	}


	uint GetNonUniformSegmentIndex(float tGlobal)
	{
		check(!InputSet.empty());
		int numPoints = InputSet.size();
		for (int i = 1; i < numPoints; i++)
		{
			if (tGlobal <= InputSet[i])
			{
				return i - 1;
			}
		}
		check(tGlobal <= InputSet.back());
	};

	float GetTValueFirst()
	{
		return InputSet.size() > 0 ? InputSet.front() : 0.f;
	}

	float GetTValueLast()
	{
		return (InputSet.size() > 0 ? InputSet.back() : float(GetNumSegments())) - 0.001f;
	}

	static inline void GenerateInputSetForNonUniformSplineKnots(const RDynamicVector<float3>& knots, float alpha, RDynamicVector<float>& outtValues)
	{
		///alpha 0.0 for Default Catmull-Rom
		///alpha 0.5 for Centripetal Catmull-Rom
		///alpha 1.0 for Chordal
		auto ComputeTDiff = [alpha](float3 p1, float3 p2) -> float
		{
			float distanceSq = VectorGetLengthSquared(p1 - p2);
			//if these points are right on top of each other, don't bother with the power calculation
			if (distanceSq < .0001)
			{
				return 0;
			}
			else
			{
				return std::pow(distanceSq, alpha * 0.5);
			}
		};

		size_t numPoints = knots.size();

		outtValues.resize(numPoints);

		outtValues[0] = 0;

		for (size_t i = 1; i < numPoints; i++)
		{
			outtValues[i] = outtValues[i - 1] + ComputeTDiff(knots[i - 1], knots[i]);
		}
	}


};


///Derivatives specified by user
struct RHermiteSplineDesigner : public RSplineDesigner
{
	virtual void Generate(float alpha)
	{
		
	}

	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal)
	{
		if (segmentIndex > (ControlPoints.size() - 2))
		{
			return float3{ 0,0,0 };
		}

		auto& p0 = ControlPoints[segmentIndex];
		auto& v0 = DerivativeSet[segmentIndex];

		auto& p1 = ControlPoints[segmentIndex + 1];
		auto& v1 = DerivativeSet[segmentIndex + 1];

		return GetPointOnHermiteCurveParametric(tLocal, p0, v0, p1, v1);
	}

};


///Knot has 2 Derivatives assigned
struct RTCBSplineDesigner : public RSplineDesigner
{

public:

	struct TCBSplineSettings
	{
		float tension = 0;
		float continuity = 0;
		float bias = 0;
	} TCBSettings;

	RTCBSplineDesigner()
	{
	}

	virtual void RenderUI() override
	{
		if (ImGui::Begin("RHermiteSplineDesigner"))
		{
			ImGui::SliderFloat("tension", &TCBSettings.tension, -2, 2);
			ImGui::SliderFloat("continuity", &TCBSettings.continuity, -2, 2);
			ImGui::SliderFloat("bias", &TCBSettings.bias, -2, 2);
		}
		ImGui::End();
	}


	virtual void Generate(float alpha)
	{
		DerivativeSet = GenerateDerivativesForKnotsTCB(ControlPoints);
	}

	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal)
	{
		uint firstKnotIndex = segmentIndex;
		uint secondKnotIndex = segmentIndex + 1;
		auto& p0 = ControlPoints[firstKnotIndex];
		auto& p1 = ControlPoints[secondKnotIndex];
		///2 derivatives per knot
		auto& fisrstKnotvOut = DerivativeSet[(firstKnotIndex * 2) + 1];
		auto& secondKnotvIn = DerivativeSet[secondKnotIndex * 2];

		return GetPointOnHermiteCurveParametric(tLocal, p0, fisrstKnotvOut, p1, secondKnotvIn);
	}

private:

	RDynamicVector<float3> GenerateDerivativesForKnotsTCB(RDynamicVector<float3>& knots)
	{
		int n = knots.size();
		check(n > 2, "Can't generate Catmull-Rom spline with 2 points");

		float t = TCBSettings.tension;
		float c = TCBSettings.continuity;
		float b = TCBSettings.bias;

		///2 derivatives per knot since we allow C0 continuity
		RDynamicVector<float3> derivatives;

		for (int i = 0; i < n; i++)
		{
			float3 pPrev = knots[std::max(i - 1, 0)];
			float3 pCurrent = knots[i];
			float3 pNext = knots[std::min(i + 1, n-1)];
			
			float3 d0 = (pCurrent - pPrev);
			float3 d1 = (pNext - pCurrent);

			///TCB extension of Catmull-Rom
			float3 vIn = d0 * ((1 - t)*(1 + b)*(1 - c))*0.5f + d1 * ((1 - t)*(1 - b)*(1 + c))*0.5f;
			float3 vOut = d0 * ((1 - t)*(1 + b)*(1 + c))*0.5f + d1 * ((1 - t)*(1 - b)*(1 - c))*0.5f;

			derivatives.push_back(vIn);
			derivatives.push_back(vOut);

		}

		check((derivatives.size() / 2) == knots.size());

		return derivatives;
	}


};

///4 control points per segment (C1 max)
struct RBezierSplineDesigner : public RSplineDesigner
{

	virtual uint GetNumSegments() override
	{
		return (ControlPoints.size() - 1) / 3;
	}

	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal)
	{
		int pointIndex = segmentIndex * 3;

		if (pointIndex + 3 > (ControlPoints.size() - 1))
		{
			return float3{ 0,0,0 };
		}

		auto& b0 = ControlPoints[pointIndex];
		auto& b1 = ControlPoints[pointIndex + 1];
		auto& b2 = ControlPoints[pointIndex + 2];
		auto& b3 = ControlPoints[pointIndex + 3];

		return GetPointOnBezierCurveParamteric(tLocal, b0, b1, b2, b3);
	}
};



struct RNonUniformBSplineDesigner : public RSplineDesigner
{

public:
	virtual void RenderUI()
	{
		if (ImGui::Begin("BSpline"))
		{
			ImGui::DragFloat("tFirst", &InputSet.front(), 0.01, InputSet.front() - 1.f, InputSet[2]);
			auto numInputs = InputSet.size();

			for (int i = 2; i < numInputs - 2; i++)
			{
				char buffer[5];
				sprintf_s(buffer, "t%i", i);
				ImGui::DragFloat(buffer, &InputSet[i], 0.01, GetTValue(i - 1), GetTValue(i + 1));
			}
			ImGui::DragFloat("tLast", &InputSet.back(), 0.01, InputSet[numInputs - 2], InputSet.back() + 1.f);

		}
		ImGui::End();
	}

	virtual void Generate(float alpha = 0.f)
	{
		///Generate Uniform Input Set 
		auto numKnots = ControlPoints.size();

		if (InputSet.empty())
		{
			InputSet.push_back(0);
		}
		if (InputSet.size() < numKnots)
		{
			while (InputSet.size() < numKnots)
			{
				InputSet.push_back(InputSet.back() + 1.f);
			}
		}
	}

	virtual float3 GetPointOnSpline(float tGlobal)
	{
		auto numPoints = ControlPoints.size();

		float3 result{ 0,0,0 };

		for (int i = 0; i < numPoints; i++)
		{
			result = result + GetControlPoint(i) * GenerateWeightForPointCoxBoor(i, 3 /*numPoints*/, tGlobal);
		}
		
		return result;
	}


	virtual float3 GetPointOnSpline(uint segmentIndex, float tGlobal)
	{
		check(segmentIndex >= GetFirstSegmentIndex(), "BSpline Segment Index starts from 1");
		float3 result = float3{ 0,0,0 };
		for (int k = 0; k < 4; k++)
		{
			float3 point = GetControlPoint(segmentIndex + k - 1);
			float weight = GenerateWeightForPointCoxBoor(segmentIndex + k - 1, 3, tGlobal);
			result = result + point * weight;
		}
		return result;
	}

	///Places where BSpline Curves join together
	float3 GetBSplineKnot(int segmentIndex)
	{
		float t = GetTValue(segmentIndex + 1);
		return GetPointOnSpline(segmentIndex, t);
	}

private:

	virtual uint GetNumSegments()
	{
		return ControlPoints.size() - 2;
	}

	int GetFirstSegmentIndex()
	{
		return 1;
	}

	int GetSegmentIndex(float tGlobal)
	{
		for (int i = GetFirstSegmentIndex(); i < GetNumSegments(); i++)
		{
			if (tGlobal < GetTValue(i + 1))
			{
				return i;
			}
		}
	}

	float3 GetControlPoint(int index)
	{
		return ControlPoints[index];
	}

	float GetTValue(int index)
	{
		if (InputSet.empty())
		{
			return index;
		}
		if (index < 2)
		{
			return InputSet.front();
		}
		if (index > (InputSet.size() - 3))
		{
			return InputSet.back();
		}
		return InputSet[index];
	}


	float GenerateWeightForPointCoxBoor(int i, int k, float tGlobal)
	{
		if (k == 0)
		{
			if (tGlobal >= GetTValue(i - 2) && tGlobal < GetTValue(i - 1))
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
		auto s1 = (GetTValue(i + k - 2) - GetTValue(i - 2));
		auto s2 = (GetTValue(i + k - 1) - GetTValue(i - 1));
		if (s1 != 0)
		{
			s1 = GenerateWeightForPointCoxBoor(i, k - 1, tGlobal) / s1;
		}
		if (s2 != 0)
		{
			s2 = GenerateWeightForPointCoxBoor(i + 1, k - 1, tGlobal) / s2;
		}
		return (tGlobal - GetTValue(i - 2)) * s1
			+ (GetTValue(i + k - 1) - tGlobal) * s2;
	}


	///-------------------------------------- Visualize --------------------------------------

	void GenerateColorsForVisualizers()
	{
		auto numKnots = ControlPoints.size();

		if (ControlPointsColors.size() < numKnots)
		{
			while (ControlPointsColors.size() < numKnots)
			{
				float3 color{ RandomNumberGenerator::Get().NextFloat(), RandomNumberGenerator::Get().NextFloat(), RandomNumberGenerator::Get().NextFloat() };

				ControlPointsColors.push_back(color);
			}
		}

		int i = 0;
		for (auto& point : RSpatialControlSet::Get().ControlPoints)
		{
			//if (!SpatialControlSet.IsPointCurrentllyInteractedWith(point))
			{
				point.Color = ControlPointsColors[i];
			}
			i++;
		}
	}

public:
	///Main
	void VisualizeSpline(float tGlobal, RGeometry& visualizer)
	{

		//tGlobal = WrapValue(tGlobal, GetTValueFirst(), GetTValueLast());
		tGlobal = WrapValuef(tGlobal, GetTValueFirst(), GetTValueLast());

		///Generate Colors
		GenerateColorsForVisualizers();

		///Render Spline Itself
		visualizer.SetColor(EColor::Red);
		VisualizeSplineCurve(visualizer);

		///Render Spline KnotsBufferCPU
		visualizer.SetColor(EColor::White);
		VisualizeSplineKnots(visualizer);

		///Render Spline Weights
		visualizer.SetColor(EColor::White);
		VisualizeWeights(float2{ 5,1 }, tGlobal, visualizer);

		///Point On A Spline 
		VisualizePointOnASpline(visualizer, tGlobal);
	}

	void VisualizeWeights(float2 visualizerOffset, float tGlobal, RGeometry& visualizer)
	{
		auto numPoints = ControlPoints.size();

		float t0 = GetTValueFirst();
		float t1 = GetTValueLast();

		///Add Line
		visualizer.SetColor(EColor::Grey);
		float3 LineStart = float3{ visualizerOffset.x + t0, visualizerOffset.y, 0 };
		float3 LineEnd= float3{ visualizerOffset.x + t1, visualizerOffset.y, 0 };
		visualizer.AddLine(LineStart, LineEnd);

		///For each point
		for (int i = 0; i < numPoints; i++)
		{
			std::function func = [&](float tGlobal)
			{
				float weight = GenerateWeightForPointCoxBoor(i, 3, tGlobal);
				return float3{ tGlobal, weight, 0 } ;
			};

			const RDynamicVector<float3>& samples = SampleFunction(func, t0, t1);
			RDynamicVector<float3> samplesCopy;

			auto numSamples = samples.size();
			for (int i = 0; i < numSamples; i++)
			{
				if (samples[i].y > 0.001)
				{
					samplesCopy.push_back(samples[i] + float3{ visualizerOffset, 0 });
				}
			}

			visualizer.SetColorRaw(GetPointColor(i));
			if (samplesCopy.size() > 3)
			{
				visualizer.AddCurvePolyboard(samplesCopy);
			}

		}

		///Visualize Point Moving Along Weights Line
		tGlobal = std::fmod(tGlobal, t1);
		float3 pos{ visualizerOffset.x + tGlobal, visualizerOffset.y, 0 };
		visualizer.SetColor(EColor::White);
		visualizer.AddCircleCameraFacing(pos);

	}

	void VisualizePointOnASpline(RGeometry& visualizer, float tGlobal)
	{
		///Visualize Point Moving Along Spline
		float t1 = GetTValueLast();
		tGlobal = std::fmod(tGlobal, t1);
		float3 pos = GetPointOnSpline(tGlobal);
		float3 color = GetColorOnSpline(tGlobal);

		float radius = 0.1;
		
		visualizer.SetColor(EColor::White);
		visualizer.AddCircleCameraFacing(pos, radius * GoldenRatio);

		visualizer.SetColorRaw(color);
		visualizer.AddCircleCameraFacing(pos, radius);
	}


	RDynamicVector<float3> ControlPointsColors;

	float3 GetPointColor(int index)
	{
		return ControlPointsColors[index];
	}
	
	float3 GetColorOnSpline(float tGlobal)
	{
		auto numPoints = ControlPointsColors.size();

		float3 result{ 0,0,0 };
		for (int i = 0; i < numPoints; i++)
		{
			result = result + GetPointColor(i) * GenerateWeightForPointCoxBoor(i, 3, tGlobal);
		}
		return result;
	}

	void VisualizeSplineCurve(RGeometry& visualizer)
	{
		std::function func = [&](float tGlobal)
		{
			return GetPointOnSpline(tGlobal);
		};

		float t0 = GetTValueFirst();
		float tLast = GetTValueLast();

		RDynamicVector<float3> samples = SampleFunction(func, t0, tLast);

		///Each Point Gets a Color, when interpolating between points, interpolate between colors as well
		std::function funcColor = [&](float tGlobal)
		{
			return GetColorOnSpline(tGlobal);
		};

		RDynamicVector<float3> colorSamples = SampleFunction(funcColor, t0, tLast);

		visualizer.AddCurvePolyboard(samples, 1, colorSamples);
	}


	void VisualizeSplineKnots(RGeometry& visualizer)
	{
		auto numSegments = GetNumSegments();
		for (int i = GetFirstSegmentIndex(); i < numSegments; i++)
		{
			float3 point = GetBSplineKnot(i);
			visualizer.AddCircle(point, 0.1);
		}
	}

};



struct RNonUniformCatmullRomSplineDesigner : RSplineDesigner
{

public:
	virtual void Generate(float alpha = 0.f)
	{
		check(ControlPoints.size() > 0);

		if (alpha != 0.f)
		{
			///Generate Non-Uniform Input Set
			GenerateInputSetForNonUniformSplineKnots(ControlPoints, alpha, InputSet);
		}
		
		GenerateDerivativesNonUniform();

		//SecondDerivativeSet = GenerateDerivativesForNonUniformSpline(DerivativeSet, InputSet);
	}

	using RSplineDesigner::GetPointOnSpline;

	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal)
	{

		float3 y0 = ControlPoints[segmentIndex];
		float3 y1 = ControlPoints[segmentIndex + 1];

		float3 d0 = DerivativeSet[segmentIndex];
		float3 d1 = DerivativeSet[segmentIndex + 1];

		bool bQuintic = !SecondDerivativeSet.empty();
		if (bQuintic)
		{
			float3 d20 = SecondDerivativeSet[segmentIndex];
			float3 d21 = SecondDerivativeSet[segmentIndex + 1];

			return TGetPointOnQuinticHermiteCurveParametric(tLocal, y0, d0, d20, y1, d1, d21);
		}
		
		return TGetPointOnHermiteCurveParametric(tLocal, y0, d0, y1, d1);
		
	}

private:
	float3 GenerateDerivativeForNonUniformSplineKnot(c_float3& pPrev, c_float3& pCurrent, c_float3& pNext, float tPrev, float tCurrent, float tNext)
	{
		return
			pPrev * (tCurrent - tNext) / ((tNext - tPrev) * (tCurrent - tPrev))
			+ pNext * (tCurrent - tPrev) / ((tNext - tPrev) * (tNext - tCurrent))
			- pCurrent * ((tCurrent - tPrev) - (tNext - tCurrent)) / ((tNext - tCurrent) * (tCurrent - tPrev));
	};

	inline float3 GenerateDerivativeForKnot(c_float3& pPrev, c_float3& pNext)
	{
		return (pNext - pPrev) * 0.5;
	};

	void GenerateDerivativesNonUniform()
	{
		const uint numKnots = ControlPoints.size();

		DerivativeSet.resize(numKnots);

		DerivativeSet[0] = (float3{});//start

		if (InputSet.empty())
		{
			for (int i = 1; i < numKnots - 1; i++)
			{
				DerivativeSet[i] = GenerateDerivativeForKnot(ControlPoints[i - 1], ControlPoints[i + 1]);
			}
		}
		else
		{
			for (int i = 1; i < numKnots - 1; i++)
			{
				DerivativeSet[i] = GenerateDerivativeForNonUniformSplineKnot(ControlPoints[i - 1], ControlPoints[i], ControlPoints[i + 1], InputSet[i - 1], InputSet[i], InputSet[i + 1]);
			}
		}
		
		DerivativeSet[numKnots-1] = float3{};//end

	}

};


struct RNonUniformLinearSplineDesigner : public RSplineDesigner
{
	virtual void Generate(float alpha)
	{
		if (alpha != 0.f)
		{
			GenerateInputSetForNonUniformSplineKnots(ControlPoints, alpha, InputSet);
		}
	}

	virtual float3 GetPointOnSpline(uint segmentIndex, float tLocal)
	{
		float3 p0 = ControlPoints[segmentIndex];
		float3 p1 = ControlPoints[segmentIndex + 1];

		return TLinearInterpolate(tLocal, p0, p1);
	}

	///t value is stored in ControlPoint.x
	float GetPointOnSpline1D(float tGlobal)
	{
		///.y is t value
		WrapValuef(tGlobal, ControlPoints.front().x, ControlPoints.back().x - 0.001);

		for (int i = 1; i < ControlPoints.size(); i++)
		{
			if (tGlobal < ControlPoints[i].x)
			{
				return TLinearInterpolate(tGlobal, ControlPoints[i - 1].y, ControlPoints[i].y, ControlPoints[i - 1].x, ControlPoints[i].x);
			}
		}
	}

};


//--------------------------
static inline void CalcCircleFrom3Points(
	float2 pt1,
	float2 pt2,
	float2 pt3,
	float2& outCenter, float& outRadius)
{

	float2 d1 = pt2 - pt1;
	float2 d2 = pt3 - pt2;

	float epsilon = 0.000001;

	// IsPerpendicular() assure that xDelta(s) are not zero
	float aSlope = d1.y / d1.x;
	float bSlope = d2.y / d2.x;

	// calc center
	outCenter.x = (
		aSlope * bSlope * (pt1.y - pt3.y) +
		bSlope * (pt1.x + pt2.x) -
		aSlope * (pt2.x + pt3.x))
		/ (2 * (bSlope - aSlope));
	outCenter.y = -1 * (outCenter.x - (pt1.x + pt2.x) / 2) / aSlope + (pt1.y + pt2.y) / 2;
	outRadius = sqrt(pow(outCenter.x - pt1.x, 2) + pow(outCenter.y - pt1.y, 2));
}







///=================================================================================================================
/// 
///************************************* Multivariate Interpolation ************************************************
/// 
///=================================================================================================================




template<typename V>
static inline V TBiLinearInterpolate(float2 uv, V a, V b, V c, V d, float2 min = float2{ 0,0 }, float2 max = float2{ 1,1 })
{
	///Interpolate horizontal a-b & c-d 
	V e = TLinearInterpolate(uv.x, a, b, min.x, max.x);
	V f = TLinearInterpolate(uv.x, c, d, min.x, max.x);

	///Interpolate e|f vertically
	return TLinearInterpolate(uv.y, e, f, min.y, max.y);
}



struct RSurface
{

	RVector2D<float3> ControlPoints2D;
	RVector2D<float3> OutputSet2D;

	RVector2D<float3> FirstPartialDerivativeX;
	RVector2D<float3> FirstPartialDerivativeY;
	RVector2D<float3> Normals;

	//Control Points
	uint2 NumPoints2D{4,4};
	float3 StartPos;
	float Spacing = 1.f;
	bool bGenerateControlPoints = true;


	float SamplingPeriodLength = 0.01f;

	void Main(c_float3& color, RGeometry& visualizer)
	{
		if (bGenerateControlPoints)
		{
			GenerateControlPoints();
			bGenerateControlPoints = false;
		}

		GenerateSamples();

		Visualize(color, visualizer);

	}

	void GenerateControlPoints()
	{
		ControlPoints2D.Resize(NumPoints2D);

		for (uint y = 0; y < NumPoints2D.y; y++)
		{
			for (uint x = 0; x < NumPoints2D.x; x++)
			{
				uint2 coord{ x,y };
				c_float3 pos{ StartPos.x + float(x) * Spacing, 0,  StartPos.y + float(y) * Spacing };
				ControlPoints2D[coord] = pos;
			}
		}

		//Assign Spatial Control
		auto& spatialControlPoints = RSpatialControlSet::Get().ControlPoints;
		spatialControlPoints.clear();

		for (uint y = 0; y < NumPoints2D.y; y++)
		{
			for (uint x = 0; x < NumPoints2D.x; x++)
			{
				uint2 coord{ x,y };
				spatialControlPoints.emplace_back(&ControlPoints2D[coord]);
			}
		}
	}

	void GenerateSamples()
	{
		//Generate Samples
		{
			
			//For each stack
			
			uint2 numSamples;
			numSamples.x = ((1) / SamplingPeriodLength) + 1;
			numSamples.y = ((1) / SamplingPeriodLength) + 1;

			OutputSet2D.Resize(numSamples);

			float2 uv;




			//TODO:Inheritance or Enum for Surface Type
			static RNonUniformCatmullRomSplineDesigner Spline;
			//static RNonUniformBSplineDesigner Spline;
			//static RBezierSplineDesigner Spline;
			

			//We are using Spline interface here since the concept is the same, we are just interpolating multiple Splines


			//Interpolate in horizontal direction first, then interpolate generated splines (size depends in the OutputSet precision) vertically
			{
				static RVector2D<float3> horizontalInterpolationResultsPerStack;
				horizontalInterpolationResultsPerStack.Resize(uint2(OutputSet2D.Size.x, ControlPoints2D.Size.y));

				//For each stack or rows
				Spline.ControlPoints.resize(ControlPoints2D.Size.x);

				for (uint stackIndex = 0; stackIndex < ControlPoints2D.Size.y; stackIndex++)
				{
					for (uint pointX = 0; pointX < ControlPoints2D.Size.x; pointX++)
					{
						Spline.ControlPoints[pointX] = ControlPoints2D[uint2(pointX, stackIndex)];
					}

					Spline.Generate();

					for (uint x = 0; x < OutputSet2D.Size.x; x++)
					{
						uv.x = 0 + float(x) * SamplingPeriodLength;

						horizontalInterpolationResultsPerStack[uint2(x, stackIndex)] = Spline.GetPointOnSplineNormalized(uv.x);

					}
				}

				//For each generated sample on X axis
				Spline.ControlPoints.resize(ControlPoints2D.Size.y);

				for (uint x = 0; x < OutputSet2D.Size.x; x++)
				{
					for (uint pointY = 0; pointY < ControlPoints2D.Size.y; pointY++)
					{
						Spline.ControlPoints[pointY] = horizontalInterpolationResultsPerStack[uint2(x, pointY)];
					}

					Spline.Generate();

					for (uint y = 0; y < OutputSet2D.Size.y; y++)
					{
						uv.y = 0 + float(y) * SamplingPeriodLength;

						auto finalPoint = Spline.GetPointOnSplineNormalized(uv.y);

						OutputSet2D[uint2{x, y}] = finalPoint;
					}
				}
			}
		}

		//End

		GenerateDerivatives();
	}


	void GenerateDerivatives()
	{

		FirstPartialDerivativeX.Resize(OutputSet2D.Size - uint2(1, 1));
		FirstPartialDerivativeY.Resize(OutputSet2D.Size - uint2(1, 1));
		Normals.Resize(OutputSet2D.Size);

		for (uint y = 0; y < OutputSet2D.Size.y - 1; y++)
		{
			for (uint x = 0; x < OutputSet2D.Size.x - 1; x++)
			{
				uint2 coord{ x,y };
				c_float3& point = OutputSet2D[coord];
				c_float3& pointNextX = OutputSet2D[coord + uint2{ 1, 0 }];
				c_float3& pointNextY = OutputSet2D[coord + uint2{ 0, 1 }];

				c_float3 derivativeX = (pointNextX - point) / SamplingPeriodLength;
				c_float3 derivativeY = (pointNextY - point) / SamplingPeriodLength;

				FirstPartialDerivativeX[coord] = derivativeX;
				FirstPartialDerivativeY[coord] = derivativeY;

				Normals[coord] = VectorNormalize(VectorCross((derivativeY), (derivativeX)));

			}

		}

	}



	void Visualize(c_float3& color, RGeometry& visualizer)
	{

		bool bRenderNormals = true;

		if (bRenderNormals)
		{
			visualizer.EnableNormals();
		}

		for (uint u = 0; u < OutputSet2D.Size.x - 1; u++)
		{
			for (uint v = 0; v < OutputSet2D.Size.y - 1; v++)
			{
				uint2 coord{ u, v };

				c_float3& point = OutputSet2D[coord];
				c_float3& pointRight = OutputSet2D[coord + uint2{1, 0}];
				c_float3& pointUp = OutputSet2D[coord + uint2{ 0, 1 }];
				c_float3& pointRightUp = OutputSet2D[coord + uint2{ 1, 1 }];

				visualizer.SetColorRaw(color);

				const bool bWire = false;
				if (bWire)
				{
					visualizer.AddLine(point, pointRight);
					visualizer.AddLine(point, pointUp);
				}
				else
				{
					visualizer.AddQuadToBuffer(point, pointUp, pointRightUp, pointRight);
				}

				///Normals
				if (bRenderNormals)
				{
					c_float3& normal = Normals[coord];
					c_float3& normalRight = Normals[coord + uint2{ 1, 0 }];
					c_float3& normalUp = Normals[coord + uint2{ 0, 1 }];
					c_float3& normalightUp = Normals[coord + uint2{ 1, 1 }];

					visualizer.AddQuadNormalsToBuffer(normal, normalUp, normalightUp, normalRight);
				}
				

				bool bVisualizeDerivatives = false;
				if (bVisualizeDerivatives)
				{
					
					c_float3& derivativeX = FirstPartialDerivativeX[coord];
					c_float3& derivativeY = FirstPartialDerivativeY[coord];

					GlobalFloat(DerivativeVisualizerScale, 0.01f, EGlobalUIType::Drag, 0.001, 10, 0.01);

					//visualizer.SetColor(EColor::Red);
					//visualizer.AddArrowPolyboard(point, point + (/*VectorNormalize*/(derivativeX) * DerivativeVisualizerScale.Value));

					//visualizer.SetColor(EColor::Green);
					//visualizer.AddArrowPolyboard(point, point + (/*VectorNormalize*/(derivativeY)*DerivativeVisualizerScale.Value));

					float3 normal = VectorCross(VectorNormalize(derivativeY), VectorNormalize(derivativeX));
					visualizer.SetColor(EColor::YellowCanary);
					visualizer.AddArrowPolyboard(point, point + (VectorNormalize(normal)*DerivativeVisualizerScale.Value));

				}
				
			}
		}


		if (bRenderNormals)
		{
			visualizer.DIsableNormals();
		}

	}



};





#pragma once

#include "LinearAlgebra.h"
#include "GenerateSimpleGeometry.h"
//#include "InterpolationAndCurves.h"

static inline float GetIntegral(float a, float b, std::function<float(float)> integrand, float delta = 0.001f)
{
	uint numSamples = fabs(b - a) / delta;
	float res = 0;
	for (uint sampleId = 0; sampleId < numSamples; sampleId++)
	{
		float t = a + sampleId * delta;
		res = res + (integrand(t) * delta);
	}
	return res;
}

template <typename V>
static inline V TGetDerivative(std::function<V(float)> func, float input, float delta = 0.001f)
{
	float dInput = delta;
	V Output1 = func(input);
	V Output2 = func(input + dInput);
	V dOutput = (Output2 - Output1);
	return dOutput / dInput;
}

template <typename VOutput, typename VInput>
static inline VOutput TGetPartialDerivative(std::function<VOutput(VInput)> func, const VInput& input, const VInput& inputDirection, float delta = 0.001f)
{
    float dInput = delta;
    VOutput Output1 = func(input);
    VOutput Output2 = func(input + (inputDirection * dInput));
    VOutput dOutput = (Output2 - Output1);
    return dOutput / dInput;
}



class SplineLibraryCalculus {
private:
    SplineLibraryCalculus() = default;

public:
    //use the gauss-legendre quadrature algorithm to numerically integrate f from a to b
    //as of this writing, hardcoded to use 13 points
    template<class IntegrandType, class Function, typename floating_t>
    inline static IntegrandType GaussLegendreQuadratureIntegral(Function f, floating_t a, floating_t b)
    {
        const size_t NUM_POINTS = 13;

        //these are precomputed :( It would be cool to compute these at compile time, but apparently
        //it's not easy to compute the points/weights just given the number of points.
        //it involves computing every root of a polynomial. which can obviously be done, but not in a reasonable amount of code
        std::array<floating_t, NUM_POINTS> quadraturePoints = {
            floating_t(0.0000000000000000),
            floating_t(-0.2304583159551348),
            floating_t(0.2304583159551348),
            floating_t(-0.4484927510364469),
            floating_t(0.4484927510364469),
            floating_t(-0.6423493394403402),
            floating_t(0.6423493394403402),
            floating_t(-0.8015780907333099),
            floating_t(0.8015780907333099),
            floating_t(-0.9175983992229779),
            floating_t(0.9175983992229779),
            floating_t(-0.9841830547185881),
            floating_t(0.9841830547185881)
        };

        std::array<floating_t, NUM_POINTS> quadratureWeights = {
            floating_t(0.2325515532308739),
            floating_t(0.2262831802628972),
            floating_t(0.2262831802628972),
            floating_t(0.2078160475368885),
            floating_t(0.2078160475368885),
            floating_t(0.1781459807619457),
            floating_t(0.1781459807619457),
            floating_t(0.1388735102197872),
            floating_t(0.1388735102197872),
            floating_t(0.0921214998377285),
            floating_t(0.0921214998377285),
            floating_t(0.0404840047653159),
            floating_t(0.0404840047653159)
        };

        floating_t halfDiff = (b - a) / 2;
        floating_t halfSum = (a + b) / 2;

        IntegrandType sum{};
        for (size_t i = 0; i < NUM_POINTS; i++)
        {
            sum += quadratureWeights[i] * f(halfDiff * quadraturePoints[i] + halfSum);
        }
        return halfDiff * sum;
    }
};



extern RGlobalFloat GlobalSamplingPrecision;

///TODO: Variable Sampling precision or LOD Precision
template <typename VOutput>
static inline const RDynamicVector<VOutput>& SampleFunction(std::function<VOutput(float)> func, float tmin = 0, float tmax = 1, float samplingPrecision = 0.01f)
{
    samplingPrecision = GlobalSamplingPrecision.Value;

    int numSamples = ((tmax - tmin) / samplingPrecision) + 1;

    static RDynamicVector<VOutput> Samples;
    Samples.resize(numSamples);
    for (int i = 0; i < numSamples; i++)
    {
        float t = tmin + float(i) * samplingPrecision;
        Samples[i] = func(t);
    }

    return Samples;

}



template <typename VOutput>
static inline const RVector2DLinear<VOutput>& SampleFunction2D(std::function<VOutput(float2)> func, float2 tmin = float2{ 0,0 }, float2 tmax = float2{ 1,1 }, float samplingPrecision = 0.01f)
{
    samplingPrecision = GlobalSamplingPrecision.Value;

    uint2 numSamples;
    numSamples.x = ((tmax.x - tmin.x) / samplingPrecision) + 1;
    numSamples.y = ((tmax.y - tmin.y) / samplingPrecision) + 1;

    static RVector2DLinear<VOutput> Samples2D;

    Samples2D.resize(numSamples.x);
    for (auto& vec : Samples2D)
    {
        vec.resize(numSamples.y);
    }

    float2 uv;

    ///Scanline order starting from bottom
    //for (int v = numSamples - 1; v >= 0; v--)
    for (int v = 0; v < numSamples.y; v++)
    {
        for (int u = 0; u < numSamples.x; u++)
        {
            uv.x = tmin.x + float(u) * samplingPrecision;
            uv.y = tmin.y + float(v) * samplingPrecision;

            Samples2D[v][u] = func(uv);
        }
    }

    return Samples2D;

}


///=========================================================================================================
///                                             VECTOR FUNCTIONS
///=========================================================================================================

///Single Input Vector Output
class RParametricFunctionVisualizer
{

public:
    RDynamicVector<float> InputSet;
    RDynamicVector<float3> OutputSet;

    ///Integrals
    RDynamicVector<float3> FirstIntergralSet;
    RDynamicVector<float3> SecondIntegralSet;

    ///Derivatives
    RDynamicVector<float3> FirstDerivativeSet;
    RDynamicVector<float3> SecondDerivativeSet;

    ///Arc Length
    RDynamicVector<float3> SegmentLength;

    RDynamicVector<float3> TangentDirectionSet;///First Derivative Normalized Direction
    RDynamicVector<float3> NormalDirectionSet;///Direction of the Change of Direction of the first derivative
    RDynamicVector<float3> BinormalDirectionSet;

    RDynamicVector<float3> OutputPlusTangentSet;///Visualization helper

    ///Sampler Settings
    float2 DomainMinMax{ 0, 1 };
    uint NumSamples = 128;
    float SamplingPeriodLength = 0.1;

    ///Visualizer Settings
    bool bVisualizeDerivative = false;
    bool bVisualizeSecondDerivative = false;
    bool bVisualizeCurvature = false;
    bool bVisualizeTorsion = false;
    bool bVisualizeMagnitudeWithColor = true;
    float VisualizedDerivativeMagnitudeScale = 0.1f;
    int DerivativeVisualizationPrecision = 1;///Visualize every nth derivative

    float PointOnAFuncPosition = 0.5f;
    bool bDrawOscullatingCircle = false;
    bool bVisualizeFrenetBasis = true;

public:
    void VisualizeFunction(std::function<float3(float)> func, RGeometry& visualizer)
    {
        RenderUI();

        SampleFunction(func);

        static const float3 RandomColor = GenerateRandomColor();
        visualizer.SetColorRaw(RandomColor);
        VisualizeSet(OutputSet, visualizer);

        if (bVisualizeDerivative)
        {
            VisualizeDerivative(OutputSet, FirstDerivativeSet, visualizer, false);
        }
        if (bVisualizeSecondDerivative)
        {
           // VisualizeDerivative(OutputPlusTangentSet, SecondDerivativeSet, visualizer, false);
            VisualizeDerivative(OutputSet, SecondDerivativeSet, visualizer, false);
        }
        if (bVisualizeCurvature)
        {
            visualizer.SetColor(EColor::YellowCanary);
            VisualizeCurvature(visualizer);
        }
        if (bVisualizeTorsion)
        {
            visualizer.SetColor(EColor::Purple);
            VisualizeTorsion(visualizer);
        }

        VisualizePointOnAFunc(PointOnAFuncPosition, visualizer, true, bDrawOscullatingCircle);
        
    }

    ///Integrate a function over domain defined by this curve
    float LineIntegral(std::function<float(float3)> func)
    {
       // float3 result{ 0,0,0 };
        float result{ 0 };
        ///We use this curve as a domain to vector function integrand
        auto numSamples = OutputSet.size()-1;

        ///Integrate
        for (int i = 0; i < numSamples; i++)
        {
            float3 x = OutputSet[i];
            float ds = VectorGetLength(OutputSet[i+1] - OutputSet[i]);

            ///F(x) * ds
            result = result + func(x) * ds;
            
        }


        if (ImGui::Begin("Parametric Function Visualizer"))
        {
            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nLineIntegralResult");
            //auto direction = VectorNormalize(result);
            //ImGui::Text("Direction: %f, %f, %f", direction.x, direction.y, direction.z);
            //ImGui::Text("Length: %f", VectorGetLength(result));
            ImGui::Text("Result: %f", (result));

        }
        ImGui::End();


        return result;

    }

    ///Circulation / Flow
    float LineIntegralOverVectorField(std::function<float3(float3)> fieldfunc)
    {
        ///Field passes through a surface

        float result{ 0 };
        auto numSamples = OutputSet.size() - 1;

        ///Integrate
        for (int i = 0; i < numSamples; i++)
        {
            ///(F dot T) ds == F dot dr
            float3 x = OutputSet[i];
            float3 F = fieldfunc(x);
            //float3 T = TangentDirectionSet[i];
            float3 dr = (OutputSet[i + 1] - OutputSet[i]);

            //result = result + VectorDot(F, T) * ds;
            result = result + VectorDot(F, dr);

        }

        if (ImGui::Begin("Parametric Function Visualizer"))
        {
            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\Line Integral Over Vector Field");
            ImGui::Text("Result: %f", result);

        }
        ImGui::End();

        return result;
    }


private:

    void RenderUI()
    {
        if (ImGui::Begin("Parametric Function Visualizer"))
        {

            uint logBase2Size;
            logBase2Size = RMath::Log2(NumSamples);
            if (ImGui::SliderInt("NumSamplesPow2", (int*)&logBase2Size, 1, 11))
            {
                NumSamples = (uint)pow(2u, logBase2Size);
            }

            ImGui::DragFloat("DomainMin", &DomainMinMax.x, 0.1, -100, 0);
            ImGui::DragFloat("DomainMax", &DomainMinMax.y, 0.1, 1, 100);

            ImGui::Text("dt: %f", SamplingPeriodLength);

            float arcLength = GetArcLength();
            ImGui::Text("ArcLength: %f", arcLength);

            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nPointOnAFunc");
            ImGui::SliderFloat("Position", &PointOnAFuncPosition, 0, 1);
            ImGui::Checkbox("FrenetBasis", &bVisualizeFrenetBasis);
            if(OutputSet.size() > 0)
            {
                int sampleIndex = int(PointOnAFuncPosition * (SecondDerivativeSet.size() - 1));
                ImGui::Text("Curvature At Point: %f", GetCurvature(sampleIndex));
                ImGui::Text("Torsion At Point: %f", GetTorsion(sampleIndex));
            }
            ImGui::Checkbox("OscullatingCircle", &bDrawOscullatingCircle);

            static bool bAnimatePoint = false;
            static float animationSpeed = 0.1f;
            ImGui::Checkbox("AnimatePoint", &bAnimatePoint);

            if (bAnimatePoint)
            {
                ImGui::SliderFloat("AnimationSpeed", &animationSpeed, 0.01, 1.f);

                auto dt = Time::Get().GetDeltaTimeSec();

                PointOnAFuncPosition += (dt * animationSpeed);
                if (PointOnAFuncPosition > 1.f)
                {
                    float fracPart = frac(PointOnAFuncPosition);
                    PointOnAFuncPosition = 0.f + fracPart;
                }
            }

            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nAdditional");
            ImGui::Checkbox("Derivative", &bVisualizeDerivative);
            ImGui::Checkbox("SecondDerivative", &bVisualizeSecondDerivative);
            ImGui::Checkbox("Curvature", &bVisualizeCurvature);
            ImGui::Checkbox("Torsion", &bVisualizeTorsion);

            ImGui::Checkbox("bVisualizeMagnitudeWithColor", &bVisualizeMagnitudeWithColor);
            if (ImGui::DragInt("DerivativeVisualizationPrecision", &DerivativeVisualizationPrecision))
            {
                DerivativeVisualizationPrecision = std::max(DerivativeVisualizationPrecision, 1);
            }
            ImGui::SliderFloat("VisualizedDerivativeMagnitudeScale", &VisualizedDerivativeMagnitudeScale, 0.01, 1.f);



        }
        ImGui::End();
    }



    ///TODO: Adaptive Sampling Precision: less precise in areas with const or high derivative
    void SampleFunction(std::function<float3(float)> func)
    {
        float domainLength = DomainMinMax.y - DomainMinMax.x;

        SamplingPeriodLength = domainLength / (NumSamples - 1);

        OutputSet.resize(NumSamples);
        InputSet.resize(NumSamples);
        for (int i = 0; i < NumSamples; i++)
        {
            float t = DomainMinMax.x + float(i) * SamplingPeriodLength;
            OutputSet[i] = func(t);
            InputSet[i] = t;
        }

        ///First derivative
        GenerateDerivative(OutputSet, FirstDerivativeSet, SamplingPeriodLength);
        GenerateDerivative(OutputSet, TangentDirectionSet, SamplingPeriodLength, true);
        GenerateOutputPlusTangentSet();
        ///Second derivative
        GenerateDerivative(FirstDerivativeSet, SecondDerivativeSet, SamplingPeriodLength);
        ///Normal
        GenerateDerivative(TangentDirectionSet, NormalDirectionSet, SamplingPeriodLength);

        ///Binormal
        uint numSamples = NormalDirectionSet.size();
        BinormalDirectionSet.resize(numSamples);
        for (int i = 0; i < numSamples; i++)
        {
            BinormalDirectionSet[i] = VectorNormalize(VectorCross(TangentDirectionSet[i], VectorNormalize(NormalDirectionSet[i])));
        }
    }
    
    ///(how quickly tangent changes direction per segment length)
    float GetCurvature(uint sampleIndex)
    {
        ///dT / ds
        auto dT = VectorGetLength(NormalDirectionSet[sampleIndex] /** SamplingPeriodLength*/);
        auto ds = VectorGetLength(FirstDerivativeSet[sampleIndex] /** SamplingPeriodLength*/);
        float curvature = dT / ds;
        return curvature;
    }

    float GetTorsion(uint sampleIndex)
    {
        ///dB / ds
        c_float3& Output1 = BinormalDirectionSet[sampleIndex];
        c_float3& Output2 = BinormalDirectionSet[sampleIndex + 1];
        float dotProduct = VectorDot(Output2, Output1);
        float3 dB = (Output2 - Output1);
        auto ds = VectorGetLength(FirstDerivativeSet[sampleIndex] * SamplingPeriodLength);
        return VectorGetLength(dB / ds);
    }

    float GetArcLength()
    {
        ///Sum of the lengths of all segments
        auto numSamples = OutputSet.size();
        float sum = 0;

        ///Integrate magnitudes
        for (int i = 1; i < numSamples; i++)
        {
            float ds = VectorGetLength(OutputSet[i] - OutputSet[i - 1]);
            //float ds = VectorGetLength(FirstDerivativeSet[i] * SamplingPeriodLength);
            sum = sum + ds;
        }
        return sum;
    }


    ///Point On A Func with Tangent and Curvature
    void VisualizePointOnAFunc(float tFrom0To1, RGeometry& visualizer, bool bNormalizeDerivatives = true, bool bVisualizeCurvature = false)
    {
        ///t is parametric
        float t = std::clamp(tFrom0To1, 0.f, 1.f);

        ///Deduce Sample Index
        int sampleIndex = int(t * (SecondDerivativeSet.size() - 1));

        float3 pointPos = OutputSet[sampleIndex];

        auto firstDerivative = FirstDerivativeSet[sampleIndex];
        auto secondDerivative = SecondDerivativeSet[sampleIndex];
        auto tangent = TangentDirectionSet[sampleIndex];
        auto normal = NormalDirectionSet[sampleIndex];

        ///Draw Point itself
        visualizer.SetColor(EColor::White);
        visualizer.AddCircleCameraFacing(pointPos, SamplingPeriodLength);

        ///Draw curvature
        ///Magnitude of (Rate of change of tangent direction / Segment Length (Length of first derivative) )
        /// High curvature - high rate of rate of change of first derivative direction over distance
        if (bVisualizeCurvature)
        {
            float curvature = VectorGetLength(normal) / VectorGetLength(firstDerivative);

            float osculatingCircleRadius = 1.f / curvature;
            float3 osculatingCircleCenter = pointPos + VectorMul(VectorNormalize(normal), osculatingCircleRadius);

            visualizer.SetColor(EColor::Purple);
            visualizer.AddCircle(osculatingCircleCenter, osculatingCircleRadius);
            //visualizer.AddCircleCameraFacing(osculatingCircleCenter, osculatingCircleRadius);

        }

        ///Visualize Velocity & Acceleration
        if(!bVisualizeFrenetBasis)
        {
            visualizer.SetColor(EColor::LightBlue);
            auto velocityArrow = VectorMul(firstDerivative, VisualizedDerivativeMagnitudeScale);
            visualizer.AddArrowPolyboard(pointPos, pointPos + velocityArrow);

            visualizer.SetColor(EColor::Green);
            auto accelerationArrow = VectorMul(secondDerivative, VisualizedDerivativeMagnitudeScale);
            visualizer.AddArrowPolyboard(pointPos + velocityArrow, pointPos + velocityArrow + accelerationArrow);
        }
        else
        {
            ///Visualize Normal
            {
                visualizer.SetColor(EColor::Yellow);
                visualizer.AddArrowPolyboard(pointPos, pointPos + (normal));
            }

            ///Visualize tangent
            {
                visualizer.SetColor(EColor::Red);
                visualizer.AddArrowPolyboard(pointPos, pointPos + tangent);
            }

            ///Visualize Frenet Frame Axis
            {
                auto frenetAxis = VectorCross(tangent, (normal));
                visualizer.SetColor(EColor::Green);
                visualizer.AddArrowPolyboard(pointPos, pointPos + frenetAxis);
            }
        }

    }

    void VisualizeCurvature(RGeometry& visualizer)
    {
        const auto numPoints = NormalDirectionSet.size();

        static RDynamicVector<float3> outputPlusCurvature;
        outputPlusCurvature.resize(numPoints);

        for (uint i = 0; i < numPoints; i++)
        {
            float curvature = GetCurvature(i);
            float3 curvatureDirection = VectorMul(VectorNormalize(SecondDerivativeSet[i]), curvature * VisualizedDerivativeMagnitudeScale);
            //float3 curvatureDirection = VectorMul(VectorNormalize(NormalDirectionSet[i]), curvature * VisualizedDerivativeMagnitudeScale);

            const auto& output = OutputSet[i];
            visualizer.AddLine(output, output + curvatureDirection, 0.5);

            outputPlusCurvature[i] = output + curvatureDirection;
        }

        visualizer.AddCurvePolyboard(outputPlusCurvature);

    }

    void VisualizeTorsion(RGeometry& visualizer)
    {
        const auto numPoints = BinormalDirectionSet.size();

        static RDynamicVector<float3> outputPlusTorsion;
        outputPlusTorsion.resize(numPoints-1);

        for (uint i = 0; i < numPoints - 1; i++)
        {
            float torsion = GetTorsion(i);
            float3 torsionDirection = VectorNormalize(SecondDerivativeSet[i]) * torsion;
            //float3 torsionDirection = VectorNormalize(NormalDirectionSet[i]) * torsion;

            const auto& output = OutputSet[i];
            visualizer.AddLine(output, output + torsionDirection, 0.5);

            outputPlusTorsion[i] = output + torsionDirection;
        }

        visualizer.AddCurvePolyboard(outputPlusTorsion);

    }


    ///Derivative Direction & Magnitude Visualizer (with Color)
    void VisualizeDerivative(const RDynamicVector<float3>& outputSet, const RDynamicVector<float3> derivativeSet, RGeometry& visualizer, bool bDrawLine = false);


    ///
    /// Statics
    ///
    static void GenerateIntegral(const RDynamicVector<float3>& derivativeSet, RDynamicVector<float3>& outIntegralSet, float dInput, float3 constantOfIntegration = float3{})
    {
        outIntegralSet.resize(derivativeSet.size() + 1);

        outIntegralSet[0] = constantOfIntegration;

        ///Each sample is a sum of previous ones
        auto sum = constantOfIntegration;
        for (int i = 1; i < outIntegralSet.size(); i++)
        {
            sum = sum + derivativeSet[i - 1] * dInput;
            outIntegralSet[i] = sum;
        }
    }

    static void GenerateDerivative(const RDynamicVector<float3>& outputSet, RDynamicVector<float3>& outDerivativeSet, float dInput, bool bNormalize = false)
    {
        uint numSamples = outputSet.size();
        check(numSamples > 1);

        outDerivativeSet.resize(numSamples - 1);

        for (int i = 0; i < numSamples - 1; i++)
        {
            c_float3& Output1 = outputSet[i];
            c_float3& Output2 = outputSet[i + 1];
            c_float3& dOutput = (Output2 - Output1);
            if (bNormalize)
            {
                outDerivativeSet[i] = VectorNormalize(outDerivativeSet[i]);/// dv/ds
            }
            else
            {
                outDerivativeSet[i] = dOutput / dInput;/// dv/dt
            }
        }
    }

    static inline void VisualizeSet(const RDynamicVector<float3>& set, RGeometry& visualizer)
    {
        visualizer.AddCurvePolyboard(set);
    }

 private:

     ///Generate Direction Curve For 2nd Derivative Visualization
     void GenerateOutputPlusTangentSet()
     {
         int numPoints = FirstDerivativeSet.size();

         OutputPlusTangentSet.resize(numPoints);

         for (int i = 0; i < numPoints; i++)
         {
             const float3 derivative = FirstDerivativeSet[i];

             ///Get final Output and Derivative Direction
             float3 output;
             float3 derivativeDirection;
             bool bNormalizeDerivative = bVisualizeMagnitudeWithColor;

             output = OutputSet[i];
             derivativeDirection = bNormalizeDerivative ? VectorNormalize(derivative) : VectorMul(derivative, VisualizedDerivativeMagnitudeScale);

             ///Cache tangents endpoints to generate curve
             OutputPlusTangentSet[i] = output + derivativeDirection;

         }
     }

};



///Single Input Vector Output
class RParametricFunctionVisualizer2D
{

public:
    RVector2D<float2> InputSet;
    RVector2D<float3> OutputSet;


    ///Derivatives
    struct PartialDerivatives
    {
        RVector2D<float3> PartialX;
        RVector2D<float3> PartialY;
        void Resize(uint2 size)
        {
            PartialX.Resize(size);
            PartialY.Resize(size);
        }
    };
    PartialDerivatives FirstDerivativeSet;
    PartialDerivatives SecondDerivativeSet;

    PartialDerivatives TangentSet;
    RVector2D<float3> NormalSet;

    RVector2D<float> PatchArea;

    ///Sampler Settings
    float2 DomainMin{ 0, 0 };
    float2 DomainMax{ 1, 1 };
    uint NumSamples = 32;
    float2 SamplingPeriodLength;
    bool bInvertNormals = false;

    ///Visualizer Settings
    float VisualizedDerivativeMagnitudeScale = 0.1;
    bool bWireframe = false;
    bool bVisualizeDerivatives = false;
    bool bPerformShading = false;


    float2 PointOnAFuncPosition{0.5,0.5};
    bool bDrawOscullatingCircle = false;
    bool bVisualizeFrenetBasis = true;

public:
    void VisualizeFunction(std::function<float3(float2)> func, RGeometry& visualizer)
    {
        RenderUI();

        SampleFunction(func);

        static const float3 RandomColor = GenerateRandomColor();
        VisualizeAsSurface(RandomColor, visualizer);

        VisualizePointOnAFunc(visualizer);

    }

    ///Integrate a function over domain defined by this surface
    float3 SurfaceIntegral(std::function<float3(float3)> func)
    {

        float3 res{0,0,0};

        ///We use this Surface as a domain to vector function integrand
        for (uint u = 0; u < PatchArea.Size.x; u++)
        {
            for (uint v = 0; v < PatchArea.Size.y; v++)
            {
                uint2 coord{ u, v };

                ///F(x) * dA
                c_float3& x = OutputSet[coord];
                float dA = PatchArea[coord];

                res = res + func(x) * dA;

            }
        }


        if (ImGui::Begin("Parametric Function Visualizer"))
        {
            ImGui::TextColored(ImVec4(1,1,0,1),"\nSurfaceIntegralResult");
            auto direction = VectorNormalize(res);
            ImGui::Text("Direction: %f, %f, %f", direction.x, direction.y, direction.z);
            ImGui::Text("Length: %f", VectorGetLength(res));

        }
        ImGui::End();


        return res;

    }




    void ComputeFlux(std::function<float3(float3)> fieldfunc)
    {
        ///Field passes through a surface

        float result = 0;

        for (uint u = 0; u < PatchArea.Size.x; u++)
        {
            for (uint v = 0; v < PatchArea.Size.y; v++)
            {
                uint2 coord{ u, v };

                ///(E dot n) dA

                c_float3& x = OutputSet[coord];
                float3 E = fieldfunc(x);
                float3 n = NormalSet[coord];
                float dA = PatchArea[coord];

                result = result + VectorDot(E, n) * dA;

            }
        }

        if (ImGui::Begin("Parametric Function Visualizer"))
        {
            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nFluxResult");
            ImGui::Text("Flux: %f", result);

        }
        ImGui::End();

    }



//private:

    void RenderUI()
    {
        if (ImGui::Begin("Parametric Function Visualizer"))
        {

            uint logBase2Size;
            logBase2Size = RMath::Log2(NumSamples);
            if (ImGui::SliderInt("NumSamplesPow2", (int*)&logBase2Size, 1, 11))
            {
                NumSamples = (uint)pow(2u, logBase2Size);
            }

            ImGui::DragFloat2("DomainMin", &DomainMin.x, 0.1, -100, 0);
            ImGui::DragFloat2("DomainMax", &DomainMax.x, 0.1, 1, 100);

            ImGui::Text("du: %f, dv: %f", SamplingPeriodLength.x, SamplingPeriodLength.y);

            if (OutputSet.Size.x > 0)
            {
                float surfaceArea = GetSurfaceArea();
                ImGui::Text("Surface Area: %f", surfaceArea);
            }
            
            ImGui::Checkbox("InvertNormals", &bInvertNormals);
            
            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nPointOnAFunc");
            ImGui::SliderFloat2("Position", &PointOnAFuncPosition.x, 0, 1);
            if (OutputSet.Size.x > 0)
            {
                float2 uv = PointOnAFuncPosition;
                uint2 sampleIndex{ uint(uv.x * (TangentSet.PartialX.Size.x - 1)), uint(uv.y * (TangentSet.PartialX.Size.y - 1)) };
                float3 normal = NormalSet[sampleIndex];
                ImGui::Text("Normal: %f, %f, %f", normal.x, normal.y, normal.z);
            }
            ImGui::Checkbox("OscullatingCircle", &bDrawOscullatingCircle);

            ImGui::TextColored(ImVec4(1, 1, 0, 1), "\nAdditional");
            ImGui::Checkbox("Wireframe", &bWireframe);
            ImGui::Checkbox("VisualizeDerivatives", &bVisualizeDerivatives);
            ImGui::Checkbox("PerformShading", &bPerformShading);

            ImGui::SliderFloat("VisualizedDerivativeMagnitudeScale", &VisualizedDerivativeMagnitudeScale, 0.01, 1.f);

        }
        ImGui::End();
    }

    void VisualizeAsSurface(float3 color, RGeometry& visualizer)
    {
        for (uint u = 0; u < OutputSet.Size.x - 1; u++)
        {
            for (uint v = 0; v < OutputSet.Size.y - 1; v++)
            {
                uint2 coord{ u, v };

                c_float3& point = OutputSet[coord];
                c_float3& pointRight = OutputSet[coord + uint2{ 1, 0 }];
                c_float3& pointUp = OutputSet[coord + uint2{ 0, 1 }];
                c_float3& pointRightUp = OutputSet[coord + uint2{ 1, 1 }];

                visualizer.SetColorRaw(color);

                if (bWireframe)
                {
                    visualizer.AddLine(point, pointRight);
                    visualizer.AddLine(point, pointUp);
                }
                else
                {
                    visualizer.AddQuadToBuffer(point, pointUp, pointRightUp, pointRight);
                }

                ///Normals
                if (bPerformShading)
                {
                    c_float3& normal = NormalSet[coord];
                    c_float3& normalRight = NormalSet[coord + uint2{ 1, 0 }];
                    c_float3& normalUp = NormalSet[coord + uint2{ 0, 1 }];
                    c_float3& normalightUp = NormalSet[coord + uint2{ 1, 1 }];

                    visualizer.AddQuadNormalsToBuffer(normal, normalUp, normalightUp, normalRight);
                }

                if (bVisualizeDerivatives)
                {
                    c_float3& normal = NormalSet[coord];
                    c_float3& tangentX = TangentSet.PartialX[coord];
                    c_float3& tangentY = TangentSet.PartialY[coord];

                    visualizer.SetColor(EColor::BlueEgyptian);
                    visualizer.AddArrowPolyboard(point, point + (VectorNormalize(tangentX) * VisualizedDerivativeMagnitudeScale));

                    visualizer.SetColor(EColor::YellowCanary);
                    visualizer.AddArrowPolyboard(point, point + (VectorNormalize(tangentY) * VisualizedDerivativeMagnitudeScale));

                    visualizer.SetColor(EColor::RedOrange);
                    visualizer.AddArrowPolyboard(point, point + (VectorNormalize(normal) * VisualizedDerivativeMagnitudeScale));

                }

            }
        }

    }

    float GetSurfaceArea()
    {
        ///Sum of the areas of all patches
        float sum = 0;

        for (uint u = 0; u < PatchArea.Size.x; u++)
        {
            for (uint v = 0; v < PatchArea.Size.y; v++)
            {
                uint2 coord{ u, v };
                float dA = PatchArea[coord];
                sum = sum + dA;
            }
        }
        return sum;
    }

    float GetPatchArea(uint2 coord)
    {
        ///Cross product returns an area of parallelogram
        float3 sideA = FirstDerivativeSet.PartialX[coord] * SamplingPeriodLength.x;
        float3 sideB = FirstDerivativeSet.PartialY[coord] * SamplingPeriodLength.y;
        float dA = VectorGetLength(VectorCross(sideA, sideB));
        return dA;
    }


    ///Point On A Func with Tangent and Curvature
    void VisualizePointOnAFunc(RGeometry& visualizer)
    {
        ///t is parametric
        float2 uv = PointOnAFuncPosition;

        ///Deduce Sample Index
        uint2 sampleIndex{ uint(uv.x * (TangentSet.PartialX.Size.x - 1)), uint(uv.y * (TangentSet.PartialX.Size.y - 1)) };

        float3 point = OutputSet[sampleIndex];

        c_float3& tangentX = FirstDerivativeSet.PartialX[sampleIndex];
        c_float3& tangentY = FirstDerivativeSet.PartialY[sampleIndex];
        c_float3& normal = VectorCross(tangentY, tangentX);

        visualizer.SetColor(EColor::Blue);
        visualizer.AddArrowPolyboard(point, point + (VectorNormalize(tangentX) * VisualizedDerivativeMagnitudeScale));

        visualizer.SetColor(EColor::Yellow);
        visualizer.AddArrowPolyboard(point, point + (VectorNormalize(tangentY) * VisualizedDerivativeMagnitudeScale));

        visualizer.SetColor(EColor::Red);
        visualizer.AddArrowPolyboard(point, point + (VectorNormalize(normal) * VisualizedDerivativeMagnitudeScale));

        ///Draw Point itself
        visualizer.SetColor(EColor::Red);
        visualizer.AddCircleCameraFacing(point, SamplingPeriodLength.x);

    }


    ///TODO: Adaptive Sampling Precision: less precise in areas with const or high derivative
    void SampleFunction(std::function<float3(float2)> func)
    {
        float2 DomainLength = DomainMax - DomainMin;

        SamplingPeriodLength.x = DomainLength.x / (NumSamples - 1);
        SamplingPeriodLength.y = DomainLength.y / (NumSamples - 1);

        OutputSet.Resize(uint2{ NumSamples, NumSamples });
        InputSet.Resize(uint2{ NumSamples, NumSamples });

        for (uint y = 0; y < NumSamples; y++)
        {
            for (uint x = 0; x < NumSamples; x++)
            {
                uint2 sampleCoord{ x,y };
                float2 uv;
                uv.x = DomainMin.x + float(x) * SamplingPeriodLength.x;
                uv.y = DomainMin.y + float(y) * SamplingPeriodLength.y;

                OutputSet[sampleCoord] = func(uv);
                InputSet[sampleCoord] = uv;
            }
        }

        GenerateDerivatives();

    }

    void GenerateDerivatives()
    {
        GeneratePartialDerivatives(OutputSet, FirstDerivativeSet, SamplingPeriodLength, false);
        GeneratePartialDerivatives(OutputSet, TangentSet, SamplingPeriodLength, true);
        
        ///Generate Normal
        uint2 numSamples = TangentSet.PartialX.Size;
        NormalSet.Resize(numSamples);
        for (uint y = 0; y < numSamples.y; y++)
        {
            for (uint x = 0; x < numSamples.x; x++)
            {
                uint2 coord{ x,y };

                if (bInvertNormals)
                {
                    NormalSet[coord] = VectorNormalize(VectorCross(TangentSet.PartialX[coord], TangentSet.PartialY[coord]));
                }
                else
                {
                    NormalSet[coord] = VectorNormalize(VectorCross(TangentSet.PartialY[coord], TangentSet.PartialX[coord]));

                }
            }
        }

        ///Compute Area of each Patch
        {
            PatchArea.Resize(numSamples);
            for (uint y = 0; y < numSamples.y; y++)
            {
                for (uint x = 0; x < numSamples.x; x++)
                {
                    uint2 coord{ x,y };
                    PatchArea[coord] = GetPatchArea(coord);
                }
            }
        }
    }

    static void GeneratePartialDerivatives(const RVector2D<float3>& outputSet, PartialDerivatives& outDerivativeSet, float2 dInput, bool bNormalize = false)
    {
        uint2 numSamples = outputSet.Size - uint2{1,1};
        check(numSamples.x > 1);

        outDerivativeSet.Resize(numSamples);

        for (uint y = 0; y < numSamples.y; y++)
        {
            for (uint x = 0; x < numSamples.x; x++)
            {
                uint2 coord{ x,y };

                c_float3& vec = outputSet[coord];
                c_float3& vecNextX = outputSet[coord + uint2{ 1,0 }];
                c_float3& vecNextY = outputSet[coord + uint2{ 0,1 }];

                c_float3 partialX = (vecNextX - vec) / dInput.x;
                c_float3 partialY = (vecNextY - vec) / dInput.y;

                outDerivativeSet.PartialX[coord] = bNormalize ? VectorNormalize(partialX) : partialX;
                outDerivativeSet.PartialY[coord] = bNormalize ? VectorNormalize(partialY) : partialY;

            }
        }
    }



};






///=========================================================================================================
///                                             SCALAR FUNCTIONS
///=========================================================================================================

///2D Input Scalar Output
class RScalarFunctionVisualizer2D
{

public:
    RVector2D<float2> InputSet;
    RVector2D<float> OutputSet;

    ///Integrals
    RVector2D<float> FirstIntergralSet;
    RVector2D<float> SecondIntegralSet;
    RVector2D<float> VolumesList;///Each Sample represents a Prism with volume width*length*height = f(x,y) * dx * dy
    float DefiniteIntegral;

    ///Derivatives
    struct PartialDerivative
    {
        float DerivativeX;
        float DerivativeY;
    };

    RVector2D<float2> GradientSet;
    RVector2D<float> LaplacianSet;///Laplacian
    RVector2D<float2> NormalSet;

    ///Separate Derivatives Sampling
    bool bComputeDerivativesFromFunc = false;
    float SamplingPeriodLengthDerivative = 0.001;

    float2 SamplingPeriodLength{ 0.1,0.1 };
    uint2 NumSamples{ 128, 128 };

    ///Visualizer Settings
    bool bVisualizeMagnitudeWithColor = true;
    float VisualizedDerivativeMagnitudeScale = 0.1f;
    int DerivativeVisualizationPrecision = 1;///Visualize every nth derivative
    bool bVisualizeDomain = true;
    bool bVisualizePointOnSurface = false;
    bool bVisualizeGradientField = false;
    bool bVisualizeOutput = true;
    bool bVisualizeGradientDivergence = false;

    float2 PointOnASurfaceUV;

    ///Function
    std::function<float(float2)> Function;
    float2 InputMin;
    float2 InputMax{ 1,1 };


    ///
    float GradientMin, GradientMax;
    float LaplacianMin, LaplacianMax;

    void Main(std::function<float(float2)> func, RGeometry& visualizer)
    {
        DrawUI();

        SampleFunction(func);

        Visualize(visualizer);

    }

    void Visualize(RGeometry& visualizer)
    {
        if (bVisualizeOutput)
        {
            static const float3 RandomColor = GenerateRandomColor();
            VisualizeOutput(RandomColor, visualizer);
        }

        if (bVisualizePointOnSurface)
        {
            VisualizePointOnASurface(visualizer);
        }

        if (bVisualizeDomain)
        {
            VisualizeDomain(visualizer);
        }

        if (bVisualizeGradientField)
        {
            VisualizeGradientField(visualizer);
        }
    }

//private:

    void SampleFunction(std::function<float(float2)>& func)
    {

        Function = func;

        float2 DomainLength = InputMax - InputMin;

        SamplingPeriodLength.x = DomainLength.x / (NumSamples.x - 1);
        SamplingPeriodLength.y = DomainLength.y / (NumSamples.y - 1);

        ///Generate 2 more samples for derivatives

        OutputSet.Resize(NumSamples);
        InputSet.Resize(NumSamples);

        NormalSet.Resize(OutputSet.Size);
        GradientSet.Resize(OutputSet.Size);

        LaplacianSet.Resize(NumSamples);

        for (uint y = 0; y < NumSamples.y; y++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                uint2 coord{ x,y };
                float2 uv;
                uv.x = InputMin.x + float(x) * SamplingPeriodLength.x;
                uv.y = InputMin.y + float(y) * SamplingPeriodLength.y;

                OutputSet[coord] = func(uv);

                InputSet[coord] = uv;

                if (bComputeDerivativesFromFunc)
                {
                    GenerateDerivativesFromFunc(uv, coord);
                }
            }
        }

        if (!bComputeDerivativesFromFunc)
        {
            GenerateDerivatives();
        }

        GenerateIntegrals();

    }


    void DrawUI()
    {
        if (ImGui::Begin("RFunctionVisualizationGenerator"))
        {

            ImGui::Checkbox("Sample Derivatives From Func", &bComputeDerivativesFromFunc);
            ImGui::Checkbox("Zero Gradient At Boundary", &bZeroGradientAtBoundary);

            static bool bfloat2Domain = true;
            ImGui::Checkbox("SquareDomain", &bfloat2Domain);
            if (!bfloat2Domain)
            {
                ImGui::DragFloat2("DomainMin", &InputMin.x, 0.1, -100, 0);
                ImGui::DragFloat2("DomainMax", &InputMax.x, 0.1, 1, 100);
                uint2 logBase2Size;
                logBase2Size.x = RMath::Log2(NumSamples.x);
                logBase2Size.y = RMath::Log2(NumSamples.y);
                if (ImGui::SliderInt2("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
                {
                    NumSamples = uint2{ (uint)pow(2u, logBase2Size.x), (uint)pow(2u, logBase2Size.y) };
                }
            }
            else
            {
                uint logBase2Size;
                logBase2Size = RMath::Log2(NumSamples.x);
                if (ImGui::SliderInt("NumSamplesPow2", (int*)&logBase2Size, 1, 11))
                {
                    NumSamples = uint2{ (uint)pow(2u, logBase2Size), (uint)pow(2u, logBase2Size) };
                }
                if (ImGui::DragFloat("DomainMin", &InputMin.x, 0.1, -100, 0))
                {
                    InputMin.y = InputMin.x;
                }
                if (ImGui::DragFloat("DomainMax", &InputMax.x, 0.1, 1, 100))
                {
                    InputMax.y = InputMax.x;
                }
            }


            ImGui::Text("SamplingPeriodLength: %f, %f", SamplingPeriodLength.x, SamplingPeriodLength.y);
            

            ImGui::NewLine();
            ImGui::Checkbox("Render Output Surface", &bVisualizeOutput);
            ImGui::Checkbox("Visualize Output Magnitude", &bVisualizeMagnitudeWithColor);
            ImGui::Checkbox("VisualizeDomain", &bVisualizeDomain);
            ImGui::Checkbox("Visualize Gradient Field", &bVisualizeGradientField);
            if (bVisualizeGradientField)
            {
                ImGui::Checkbox("Visualize Gradient Divergence", &bVisualizeGradientDivergence);
            }
            ImGui::Text("Gradient MinMax: %f, %f", GradientMin, GradientMax);
            ImGui::Text("Laplacian MinMax: %f, %f", LaplacianMin, LaplacianMax);

            ImGui::TextColored(ImVec4{ 1,1,0,1 }, "\nPoint On A Surface:");
            ImGui::Checkbox("PointOnSurface", &bVisualizePointOnSurface);
            if (bVisualizePointOnSurface)
            {
                ImGui::DragFloat("PointOnASurfaceX", &PointOnASurfaceUV.x, SamplingPeriodLength.x, InputMin.x, InputMax.x);
                ImGui::DragFloat("PointOnASurfaceY", &PointOnASurfaceUV.y, SamplingPeriodLength.y, InputMin.y, InputMax.y);

                int2 sampleIndex;
                sampleIndex.x = (PointOnASurfaceUV.x - InputMin.x) / SamplingPeriodLength.x;
                sampleIndex.y = (PointOnASurfaceUV.y - InputMin.y) / SamplingPeriodLength.y;

                float output = OutputSet[sampleIndex];
                ImGui::Text("Cur Output: %f", output);

                float2 gradient = GradientSet[sampleIndex];

                ImGui::Text("Gradient: %f, %f", gradient.x, gradient.y);

                float gradDiv = LaplacianSet[sampleIndex];
                ImGui::Text("Gradient Divergence: %f", gradDiv);
            }
            
            ImGui::NewLine();
            ImGui::Text("Integral (whole domain): %f", DefiniteIntegral);


        }
        ImGui::End();
    }


    void VisualizeGradientField(RGeometry& visualizer);


    void VisualizeDomain(RGeometry& visualizer)
    {
        visualizer.SetColor(EColor::BlueEgyptian);
        float2 center = (InputMax + InputMin) / 2;
        float2 extent = (InputMax - InputMin) / 2;
        visualizer.AddPlaneLineListXZ(float3{ center.x, 0, center.y }, extent);
    }

    bool bZeroGradientAtBoundary = true;

    void GenerateDerivatives()
    {
        for (uint y = 0; y < GradientSet.Size.y; y++)
        {
            for (uint x = 0; x < GradientSet.Size.x; x++)
            {
                uint2 coord{ x,y };
                const float& pointCur = OutputSet[coord];

                float pointNextX = pointCur;
                float pointNextY = pointCur;

                if ((coord.x + 1) < NumSamples.x)
                {
                    pointNextX = OutputSet[coord + uint2{ 1, 0 }];
                }
                if ((coord.y + 1) < NumSamples.x)
                {
                    pointNextY = OutputSet[coord + uint2{ 0, 1 }];
                }

                const float derivativeX = (pointNextX - pointCur) / SamplingPeriodLength.x;
                const float derivativeY = (pointNextY - pointCur) / SamplingPeriodLength.y;

                GradientSet[coord] = float2{ derivativeX, derivativeY };
                NormalSet[coord] = float2{ -derivativeX, -derivativeY };

                if (bZeroGradientAtBoundary)
                {
                    if ((y == (GradientSet.Size.y - 1)) || (x == (GradientSet.Size.x - 1)))
                    {
                        GradientSet[coord] = float2{ 0,0 };
                    }
                }


            }

        }


        ///Compute Second Derivative
        for (uint y = 0; y < LaplacianSet.Size.y; y++)
        {
            for (uint x = 0; x < LaplacianSet.Size.x; x++)
            {
                uint2 coord{ x,y };

                c_float2& gradientCur = GradientSet[coord];

                float2 gradientPrevX = gradientCur;
                float2 gradientPrevY = gradientCur;

                gradientPrevX = float2{ 0,0 };
                gradientPrevY = float2{ 0,0 };

                if ((int(coord.x) - 1) >= 0)
                {
                    gradientPrevX = GradientSet[coord - uint2{ 1,0 }];
                }
                if ((int(coord.y) - 1) >= 0)
                {
                    gradientPrevY = GradientSet[coord - uint2{ 0,1 }];
                }

                c_float2 partialX = (gradientCur - gradientPrevX) / SamplingPeriodLength.x;
                c_float2 partialY = (gradientCur - gradientPrevY) / SamplingPeriodLength.y;

                LaplacianSet[coord] = partialX.x + partialY.y;


            }

        }

    }

    void GenerateDerivativesFromFunc(float2 uv, uint2 sampleIndex)
    {
        float dInput = SamplingPeriodLengthDerivative;

        const float& point = Function(uv);
        const float& pointNextX = Function(uv + float2{ dInput, 0 });
        const float& pointNextY = Function(uv + float2{ 0, dInput });
        const float& pointPrevX = Function(uv - float2{ dInput, 0 });
        const float& pointPrevY = Function(uv - float2{ 0, dInput });

        const float derivativeX = (pointNextX - point) / dInput;
        const float derivativeY = (pointNextY - point) / dInput;

        GradientSet[sampleIndex] = float2{ derivativeX, derivativeY };
        NormalSet[sampleIndex] = float2{ -derivativeX, -derivativeY };

        LaplacianSet[sampleIndex] = (pointPrevX + pointNextX + pointPrevY + pointNextY - (4 * point)) / pow(dInput, 2.f);

    }

    float2 GetGradientFromFunc(float2 uv)
    {
        float dInput = SamplingPeriodLengthDerivative;

        const float& point = Function(uv);
        const float& pointNextX = Function(uv + float2{dInput, 0});
        const float& pointNextY = Function(uv + float2{ 0, dInput });

        const float derivativeX = (pointNextX - point) / dInput;
        const float derivativeY = (pointNextY - point) / dInput;

        return float2{ derivativeX, derivativeY };
    }

    float GetLaplacianFromFunc(float2 uv)
    {
        float dInput = SamplingPeriodLengthDerivative;

        const float& point = Function(uv);
        const float& pointNextX = Function(uv + float2{ dInput, 0 });
        const float& pointNextY = Function(uv + float2{ 0, dInput });
        const float& pointPrevX = Function(uv - float2{ dInput, 0 });
        const float& pointPrevY = Function(uv - float2{ 0, dInput });

        return (pointPrevX + pointNextX + pointPrevY + pointNextY - (4 * point)) / pow(dInput, 2.f);
    }

    void GenerateIntegrals()
    {
        FirstIntergralSet.Resize(OutputSet.Size);

        float constant = 0.f;

       // FirstIntergralSet[uint2{ 0,0 }] = constant;

        VolumesList.Resize(NumSamples);

        float Sum = 0;

        for (int x = 0; x < NumSamples.x; x++)
        {
            for (int y = 0; y < NumSamples.y; y++)
            {
                int2 coord{ x, y };
                ///f(x,y) * dx * dy
                float volume = OutputSet[coord] * SamplingPeriodLength.x * SamplingPeriodLength.y;
                volume = isnan(volume) ? 0 : abs(volume);
                VolumesList[coord] = volume;
                Sum = Sum + volume;
            }
        }

        DefiniteIntegral = Sum;

    }

    void VisualizeOutput(float3 color, RGeometry& visualizer);

    void VisualizePointOnASurface(RGeometry& visualizer)
    {
        int2 sampleIndex;
        sampleIndex.x = (PointOnASurfaceUV.x - InputMin.x) / SamplingPeriodLength.x;
        sampleIndex.y = (PointOnASurfaceUV.y - InputMin.y) / SamplingPeriodLength.y;

        float3 pointPos;
        float2 input = InputSet[sampleIndex];
        float output = OutputSet[sampleIndex];
        pointPos = float3{ input.x , output , input.y };

        {
            if (bVisualizeDomain)
            {
                c_float3 inputPos{ input.x , 0 , input.y };
                visualizer.SetColor(EColor::LightBlue);
                visualizer.AddCircleCameraFacing(inputPos, SamplingPeriodLength.x);
            }
            
            visualizer.SetColor(EColor::White);
            visualizer.AddCircleCameraFacing(pointPos, SamplingPeriodLength.x);

        }

        ///Gradient & Normal
        float2 g = GradientSet[sampleIndex];
        float3 gradient{ g.x, 0, g.y };
        visualizer.SetColor(EColor::White);
        visualizer.AddArrowPolyboard(pointPos, pointPos + VectorNormalize(gradient));

        float3 normal{ -g.x, 1, -g.y };
        visualizer.SetColor(EColor::Purple);
        visualizer.AddArrowPolyboard(pointPos, pointPos + VectorNormalize(normal));
        
    }


    ///******************************* Polar Coordinates ******************************************

    float PolarDomainAngleMin = 0;
    float PolarDomainAngleMax = PiX2;
    float PolarDomainRadiusMin = 1;
    float PolarDomainRadiusMax = 2;

    uint PolarDomainAngleNumSamples = 32;
    uint PolarDomainRadiusNumSamples = 8;

    float PolarDomainAngleDelta;
    float PolarDomainRadiusDelta;

    RVector2D<float2> PolarInputSet;

public:
    void VisualizeFunctionPolarCoord(std::function<float(float2)> func, RGeometry& visualizer)
    {

        DrawUIPolarCoord();
        
        float angularLength = PolarDomainAngleMax - PolarDomainAngleMin;
        float RadiusLength = PolarDomainRadiusMax - PolarDomainRadiusMin;

        PolarDomainAngleDelta = angularLength / (PolarDomainAngleNumSamples - 1);
        PolarDomainRadiusDelta = RadiusLength / (PolarDomainRadiusNumSamples - 1);

        NumSamples = uint2{ PolarDomainAngleNumSamples, PolarDomainRadiusNumSamples };

        OutputSet.Resize(NumSamples);
        InputSet.Resize(NumSamples);
        PolarInputSet.Resize(NumSamples);

        for (uint y = 0; y < NumSamples.y; y++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                float theta = PolarDomainAngleMin + float(x) * PolarDomainAngleDelta;
                float r = PolarDomainRadiusMin + float(y) * PolarDomainRadiusDelta;

                OutputSet[uint2{ x,y }] = func(float2{theta, r});

                InputSet[uint2{ x,y }] = MapPolarToCartesian(theta, r);
                PolarInputSet[uint2{ x,y }] = float2{ theta, r };
            }
        }

        GenerateIntegralsPolar();


        static const float3 RandomColor = GenerateRandomColor();
        VisualizeOutput(RandomColor, visualizer);
    }

private:
    void GenerateIntegralsPolar()
    {
        FirstIntergralSet.Resize(OutputSet.Size);

        float constant = 0.f;

        // FirstIntergralSet[uint2{ 0,0 }] = constant;

        VolumesList.Resize(NumSamples);

        float Sum = 0;

        for (int y = 0; y < NumSamples.y; y++)
        {
            for (int x = 0; x < NumSamples.x; x++)
            {
                int2 coord{ x, y };
                ///f(x,y) * r * dr * dTheta
                float r = PolarInputSet[coord].y;
                float volume = OutputSet[coord] * r * PolarDomainRadiusDelta * PolarDomainAngleDelta;
                volume = isnan(volume) ? 0 : abs(volume);
                VolumesList[coord] = volume;
                Sum = Sum + volume;

                if (PolarDomainRadiusMin == 0.f && y == 0 )
                {
                    break;
                }
            }
        }

        DefiniteIntegral = Sum;

    }

    void DrawUIPolarCoord()
    {
        if (ImGui::Begin("RFunctionVisualizationGenerator"))
        {

            ImGui::Checkbox("bVisualizeMagnitudeWithColor", &bVisualizeMagnitudeWithColor);

            uint2 logBase2Size;
            logBase2Size.x = RMath::Log2(PolarDomainAngleNumSamples);
            logBase2Size.y = RMath::Log2(PolarDomainRadiusNumSamples);
            if (ImGui::SliderInt2("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
            {
                PolarDomainAngleNumSamples = (uint)pow(2u, logBase2Size.x);
                PolarDomainRadiusNumSamples = (uint)pow(2u, logBase2Size.y);
            }

            ImGui::Text("SamplingPeriodLength: %f, %f", SamplingPeriodLength.x, SamplingPeriodLength.y);

            ImGui::SliderAngle("AngleMin", &PolarDomainAngleMin, 0, (PolarDomainAngleMax / Pi * 180.f));
            ImGui::SliderAngle("AngleMax", &PolarDomainAngleMax, (PolarDomainAngleMin / Pi * 180.f) + 10);

            ImGui::DragFloat("RadiusMin", &PolarDomainRadiusMin, 0.1);
            ImGui::DragFloat("RadiusMax", &PolarDomainRadiusMax, 0.1, PolarDomainRadiusMin + 0.1, 50);


            ImGui::Text("Integral (whole domain): %f", DefiniteIntegral);


        }
        ImGui::End();
    }


};


///3D Input Scalar Output
class RScalarFunctionVisualizer3D
{
public:

    RVector3D<float3> InputGrid;
    RVector3D<float> OutputGrid;

    float3 SamplingPeriodLength{};
    uint3 NumSamples{ 4, 4, 4 };
    float3 DomainMin{-1,-1,-1};
    float3 DomainMax{ 1,1,1 };

    float DefiniteIntegral;///Integral of current domain


    ///Derivatives
    struct PartialDerivatives
    {
        float PartialX;
        float PartialY;
        float PartialZ;
    };
    RVector3D<PartialDerivatives> FirstDerivativeSet;
    RVector3D<PartialDerivatives> SecondDerivativeSet;

    RVector3D<float3> GradientSet;
    RVector3D<float> LaplacianSet;///Laplacian
    RVector3D<float3> NormalSet;

    ///Point on a Func
    float3 PointOnAFunc{ 0.5,0.5,0.5 };

    ///Visualizer Settings
    int DerivativeVisualizationPrecision = 1;///Visualize every nth derivative
    bool bVisualizeDomain = true;
    bool bVisualizePointOnSurface = true;
    bool bVisualizeGradientField = false;
    bool bVisualizeOutput = true;
    bool bVisualizeGradientDivergence = false;


    void VisualizeFunction(std::function<float(float3)> func, RGeometry& visualizer)
    {

        DrawUI();

        SampleFunction(func);
        
        if (bVisualizeOutput)
        {
            Visualize(visualizer);

        }

        VisualizePointOnAFunc(visualizer);

        if (bVisualizeGradientField)
        {
            VisualizeGradientField(visualizer);
        }
    }


    void SampleFunction(std::function<float(float3)>& func)
    {
        float3 DomainLength = DomainMax - DomainMin;

        SamplingPeriodLength.x = DomainLength.x / (NumSamples.x - 1);
        SamplingPeriodLength.y = DomainLength.y / (NumSamples.y - 1);
        SamplingPeriodLength.z = DomainLength.z / (NumSamples.z - 1);

        InputGrid.Resize(NumSamples);
        OutputGrid.Resize(NumSamples);

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    uint3 sampleCoord{ x,y,z };
                    float3 uvw;
                    uvw.x = DomainMin.x + float(x) * SamplingPeriodLength.x;
                    uvw.y = DomainMin.y + float(y) * SamplingPeriodLength.y;
                    uvw.z = DomainMin.z + float(z) * SamplingPeriodLength.z;

                    OutputGrid[sampleCoord] = func(uvw);

                    InputGrid[sampleCoord] = uvw;
                }
            }
        }


        GenerateDerivatives();


        GenerateIntegrals();

    }

    void GenerateDerivatives()
    {
        FirstDerivativeSet.Resize(OutputGrid.Size);
        GradientSet.Resize(OutputGrid.Size);
        SecondDerivativeSet.Resize(OutputGrid.Size);
        LaplacianSet.Resize(OutputGrid.Size);

        for (uint z = 0; z < NumSamples.z - 1; z++)
        {
            for (uint y = 0; y < NumSamples.y - 1; y++)
            {
                for (uint x = 0; x < NumSamples.x - 1; x++)
                {
                    uint3 coord{ x,y,z };

                    float point = OutputGrid[coord];
                    float pointNextX = OutputGrid[coord + uint3{ 1,0,0 }];
                    float pointNextY = OutputGrid[coord + uint3{ 0,1,0 }];
                    float pointNextZ = OutputGrid[coord + uint3{ 0,0,1 }];

                    float partialX = (pointNextX - point) / SamplingPeriodLength.x;
                    float partialY = (pointNextY - point) / SamplingPeriodLength.y;
                    float partialZ = (pointNextZ - point) / SamplingPeriodLength.z;

                    FirstDerivativeSet[coord].PartialX = partialX;
                    FirstDerivativeSet[coord].PartialY = partialY;
                    FirstDerivativeSet[coord].PartialZ = partialZ;

                    GradientSet[coord] = float3{ partialX , partialY, partialZ };
                }
            }
        }

        for (uint z = 0; z < NumSamples.z - 2; z++)
        {
            for (uint y = 0; y < NumSamples.y - 2; y++)
            {
                for (uint x = 0; x < NumSamples.x - 2; x++)
                {
                    uint3 coord{ x,y,z };

                    c_float3& gradient = GradientSet[coord];
                    c_float3& gradientNextX = GradientSet[coord + uint3{ 1,0,0 }];
                    c_float3& gradientNextY = GradientSet[coord + uint3{ 0,1,0 }];
                    c_float3& gradientNextZ = GradientSet[coord + uint3{ 0,0,1 }];

                    c_float3 partialX = (gradientNextX - gradient) / SamplingPeriodLength.x;
                    c_float3 partialY = (gradientNextY - gradient) / SamplingPeriodLength.y;
                    c_float3 partialZ = (gradientNextZ - gradient) / SamplingPeriodLength.z;

                    SecondDerivativeSet[coord].PartialX = partialX.x;
                    SecondDerivativeSet[coord].PartialY = partialY.y;
                    SecondDerivativeSet[coord].PartialZ = partialZ.z;

                    LaplacianSet[coord] = partialX.x + partialY.y + partialZ.z;
                }
            }
        }


        /*
        ///Copy last elements
        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                uint y = NumSamples.y - 1;

                uint3 coord{ x,y,z };

                FirstDerivativeSet.PartialX[coord] = FirstDerivativeSet.PartialX[coord - uint3{ 0,1,0 }];
                FirstDerivativeSet.PartialY[coord] = FirstDerivativeSet.PartialY[coord - uint3{ 0,1,0 }];
                Divergence[coord] = Divergence[coord - uint3{ 0,1,0 }];
                Curl[coord] = Curl[coord - uint3{ 0,1,0 }];
            }
            for (uint y = 0; y < NumSamples.x; y++)
            {
                uint x = NumSamples.x - 1;

                uint3 coord{ x,y,z };

                FirstDerivativeSet.PartialX[coord] = FirstDerivativeSet.PartialX[coord - uint3{ 1,0,0 }];
                FirstDerivativeSet.PartialY[coord] = FirstDerivativeSet.PartialY[coord - uint3{ 1,0,0 }];
                Divergence[coord] = Divergence[coord - uint3{ 1,0,0 }];
                Curl[coord] = Curl[coord - uint3{ 1,0,0 }];
            }
        }
        */
    }

    void VisualizePointOnAFunc(RGeometry& visualizer)
    {
        float3 fsampleIndex = MapToRange(PointOnAFunc, float3{ 0,0,0 }, float3{ 1,1,1 }, float3{ 0,0,0 }, float3{ (float)NumSamples.x-1, (float)NumSamples.y-1, (float)NumSamples.z-1 });
        uint3 sampleIndex{ (uint)fsampleIndex.x, (uint)fsampleIndex.y, (uint)fsampleIndex.z };

        float3 pointPos;
        float3 input = InputGrid[sampleIndex];
        float output = OutputGrid[sampleIndex];

        visualizer.SetColor(EColor::LightBlue);
        visualizer.AddCircleCameraFacing(input, SamplingPeriodLength.x * 0.1);

    }

private:

    void DrawUI()
    {
        if (ImGui::Begin("RFunctionVisualizationGenerator3D"))
        {
            
            static bool bfloat2Domain = true;
            ImGui::Checkbox("CubeDomain", &bfloat2Domain);
            if (!bfloat2Domain)
            {
                ImGui::DragFloat3("DomainMin", &DomainMin.x, 0.1, -100, 0);
                ImGui::DragFloat3("DomainMax", &DomainMax.x, 0.1, 1, 100);

                uint3 logBase2Size;
                logBase2Size.x = RMath::Log2(NumSamples.x);
                logBase2Size.y = RMath::Log2(NumSamples.y);
                logBase2Size.z = RMath::Log2(NumSamples.z);
                if (ImGui::SliderInt3("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
                {
                    NumSamples = uint3{ (uint)pow(2u, logBase2Size.x), (uint)pow(2u, logBase2Size.y), (uint)pow(2u, logBase2Size.z) };
                }
            }
            else
            {
                if (ImGui::DragFloat("DomainMin", &DomainMin.x, 0.1, -100, 0))
                {
                    DomainMin.y = DomainMin.z = DomainMin.x;
                }
                if (ImGui::DragFloat("DomainMax", &DomainMax.x, 0.1, 1, 100))
                {
                    DomainMax.y = DomainMax.z = DomainMax.x;
                }

                uint3 logBase2Size;
                logBase2Size.x = RMath::Log2(NumSamples.x);
                if (ImGui::SliderInt("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
                {
                    NumSamples = uint3{ (uint)pow(2u, logBase2Size.x), (uint)pow(2u, logBase2Size.x), (uint)pow(2u, logBase2Size.x) };
                }
            }

            ImGui::Text("SamplingPeriodLength: %f, %f, %f", SamplingPeriodLength.x, SamplingPeriodLength.y, SamplingPeriodLength.z);
            ImGui::Text("NumSamplesAll: %i", NumSamples.x * NumSamples.y * NumSamples.z);

            ImGui::NewLine();
            ImGui::Checkbox("Render Output Surface", &bVisualizeOutput);
            ImGui::Checkbox("Visualize Gradient Field", &bVisualizeGradientField);
            if (bVisualizeGradientField)
            {
                ImGui::Checkbox("Visualize Gradient Divergence", &bVisualizeGradientDivergence);
            }


            ImGui::TextColored(ImVec4{ 1,1,0,1 }, "\nPoint On A Surface:");
            ImGui::Checkbox("PointOnSurface", &bVisualizePointOnSurface);
            if (bVisualizePointOnSurface)
            {

                ImGui::DragFloat3("PointOnFunc", &PointOnAFunc.x, 0.01, 0, 1);

                float3 fsampleIndex = MapToRange(PointOnAFunc, float3{ 0,0,0 }, float3{ 1,1,1 }, float3{ 0,0,0 }, float3{ (float)NumSamples.x - 1, (float)NumSamples.y - 1, (float)NumSamples.z - 1 });
                uint3 sampleIndex{ (uint)fsampleIndex.x, (uint)fsampleIndex.y, (uint)fsampleIndex.z };

                if ((sampleIndex.x < OutputGrid.Size.x) && (sampleIndex.y < OutputGrid.Size.y) && (sampleIndex.z < OutputGrid.Size.z))
                {
                    float output = OutputGrid[sampleIndex];
                    ImGui::Text("Cur Output: %f", output);

                    float3 gradient = GradientSet[sampleIndex];

                    ImGui::Text("Gradient: %f, %f", gradient.x, gradient.y, gradient.z);

                    float gradDiv = LaplacianSet[sampleIndex];
                    ImGui::Text("Gradient Divergence: %f", gradDiv);
                }
                
            }


            ImGui::Text("Integral (whole domain): %f", DefiniteIntegral);

        }
        ImGui::End();
    }


    void GenerateIntegrals()
    {
        float Sum = 0;

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    uint3 coord{x, y, z};
                    ///f(x,y,z) * dx * dy * dz
                    float volume = OutputGrid[coord] * SamplingPeriodLength.x * SamplingPeriodLength.y * SamplingPeriodLength.z;
                    volume = abs(volume);
                    Sum = Sum + volume;
                }
            }
        }

        DefiniteIntegral = Sum;

    }

    void Visualize(RGeometry& visualizer)
    {
        float outputMax = OutputGrid[uint3{ 0,0,0 }];
        float outputMin = outputMax;

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                for (uint y = 0; y < NumSamples.y; y++)
                {
                    uint3 coord{ x,y,z };

                    ///Render Input Grid

                    float curOutput = OutputGrid[coord];

                    outputMax = std::max(outputMax, curOutput);
                    outputMin = std::min(outputMin, curOutput);
                }
            }
        }

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                for (uint y = 0; y < NumSamples.y; y++)
                {
                    uint3 coord{ x,y,z };

                    ///Render Input Grid

                    float3 curInput = InputGrid[coord];

                    float curOutput = OutputGrid[coord];

                    float finColor = MapTo01Range(curOutput, outputMin, outputMax);

                    if (bUnitSphere)
                    {
                        visualizer.SetColor(EColor::White);
                        //visualizer.AddCircleCameraFacing(VectorNormalize(curInput) * curOutput * SphericalDomainRadiusMax, 0.1);
                        visualizer.AddCircleCameraFacing(curInput * curOutput, 0.1);
                    }
                    else
                    {
                        visualizer.SetColorRaw(float3{ finColor ,finColor ,finColor });
                        visualizer.AddCircleCameraFacing(curInput, SamplingPeriodLength.x / 5);
                    }

                }
            }
        }

    }


    void VisualizeGradientField(RGeometry& visualizer);


    ///************************************************* Spherical Coordinates ***************************

    float2 SphericalDomainAnglesMin{ 0,-PiDiv2 };
    float2 SphericalDomainAnglesMax{ PiX2,PiDiv2 };

    float SphericalDomainRadiusMin = 1;
    float SphericalDomainRadiusMax = 2;

    uint SphericalDomainAnglesNumSamples = 16;
    uint SphericalDomainRadiusNumSamples = 8;

    float2 SphericalDomainAnglesDelta;
    float SphericalDomainRadiusDelta;
    
    bool bUnitSphere = false;

    bool bVisualizeMagnitudeWithColor = true;

    RVector3D<float3> SphericalInputGrid;


public:
    void VisualizeFunctionSphericalCoord(std::function<float(float3)> func, RGeometry& visualizer)
    {

        DrawUISphericalCoord();


        float2 angularLength = SphericalDomainAnglesMax - SphericalDomainAnglesMin;
        float RadiusLength = SphericalDomainRadiusMax - SphericalDomainRadiusMin;

        SphericalDomainAnglesDelta = angularLength / (SphericalDomainAnglesNumSamples - 1);
        SphericalDomainRadiusDelta = RadiusLength / (SphericalDomainRadiusNumSamples - 1);

        if (bUnitSphere)
        {
            SphericalDomainRadiusNumSamples = 1;
            SphericalDomainRadiusMin = SphericalDomainRadiusMax;
            SphericalDomainRadiusDelta = 1;
        }

        NumSamples = uint3{ SphericalDomainAnglesNumSamples, SphericalDomainAnglesNumSamples, SphericalDomainRadiusNumSamples };

        SamplingPeriodLength.x = SphericalDomainRadiusDelta;

        OutputGrid.Resize(NumSamples);
        InputGrid.Resize(NumSamples);
        SphericalInputGrid.Resize(NumSamples);

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    uint3 sampleCoord{ x,y,z };
                    float theta = SphericalDomainAnglesMin.x + float(x) * SphericalDomainAnglesDelta.x;
                    float phi = SphericalDomainAnglesMin.y + float(y) * SphericalDomainAnglesDelta.y;
                    float r = SphericalDomainRadiusMin + float(z) * SphericalDomainRadiusDelta;

                    OutputGrid[sampleCoord] = func(float3{ theta, phi, r });

                    InputGrid[sampleCoord] = MapSphericalToCartesian(theta, phi, r);
                    SphericalInputGrid[sampleCoord] = float3{theta, phi, r};
                }
            }
        }

        GenerateIntegralsSphericalCoords();


        Visualize(visualizer);
        //VisualizeUnitSphere()
    }

private:

    void DrawUISphericalCoord()
    {
        if (ImGui::Begin("RFunctionVisualizationGenerator"))
        {

            ImGui::Checkbox("UnitSphere", &bUnitSphere);


            ImGui::Checkbox("bVisualizeMagnitudeWithColor", &bVisualizeMagnitudeWithColor);


            uint2 logBase2Size;
            logBase2Size.x = RMath::Log2(SphericalDomainAnglesNumSamples);
            logBase2Size.y = RMath::Log2(SphericalDomainRadiusNumSamples);
            if (ImGui::SliderInt2("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
            {
                SphericalDomainAnglesNumSamples = (uint)pow(2u, logBase2Size.x);
                SphericalDomainRadiusNumSamples = (uint)pow(2u, logBase2Size.y);
            }

            ImGui::Text("SamplingPeriodLength: %f, %f", SamplingPeriodLength.x, SamplingPeriodLength.y);


            float2 thetaMinMax{ SphericalDomainAnglesMin.x, SphericalDomainAnglesMax.x };
            if (ImGui::SliderFloat2("thetaMinMax", &thetaMinMax.x, 0, PiX2))
            {
                SphericalDomainAnglesMin.x = thetaMinMax.x;
                SphericalDomainAnglesMax.x = thetaMinMax.y;
            }
            float2 phiMinMax{ SphericalDomainAnglesMin.y, SphericalDomainAnglesMax.y };
            if (ImGui::SliderFloat2("phiMinMax", &phiMinMax.x, -PiDiv2, PiDiv2))
            {
                SphericalDomainAnglesMin.y = phiMinMax.x;
                SphericalDomainAnglesMax.y = phiMinMax.y;
            }

            ImGui::DragFloat("RadiusMin", &SphericalDomainRadiusMin, 0.1);
            ImGui::DragFloat("RadiusMax", &SphericalDomainRadiusMax, 0.1, SphericalDomainRadiusMin + 0.1, 50);

            ImGui::Text("Integral (whole domain): %f", DefiniteIntegral);

        }
        ImGui::End();
    }


    void GenerateIntegralsSphericalCoords()
    {
        float Sum = 0;

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    uint3 coord{ x, y, z };
                    ///f(x,y,z) * dr * dTheta * dPhi * r*r * sin(phi)
                    const auto& curSphericalCoords = SphericalInputGrid[coord];
                    float r = curSphericalCoords.z;
                    float phi = curSphericalCoords.y;
                    float volume = OutputGrid[coord] * SphericalDomainRadiusDelta * SphericalDomainAnglesDelta.x * SphericalDomainAnglesDelta.y * pow(r, 2.f) * sin(phi) /*or cos(phi)*/;
                    volume = abs(volume);
                    Sum = Sum + volume;

                    if (SphericalDomainRadiusMin == 0.f && z == 0)
                    {
                        break;
                    }
                }
            }
        }

        DefiniteIntegral = Sum;

    }

    void VisualizeUnitSphere(float3 color, RGeometry& visualizer);



};












static inline float GaussianFunc(float2 uv, float2 scale)
{
    auto dotRes = VectorDot(uv, float2{ uv.x * scale.x, uv.y * scale.y });
    return pow(e, dotRes * (-0.5f));
}

static inline float GaussianFunc(float3 vec, float3 scale)
{
    auto dotRes = VectorDot(vec, float3{ vec.x * scale.x, vec.y * scale.y, vec.z * scale.z });
    return pow(e, dotRes * (-0.5f));
}


static inline float GaussianFunc(float x, float a)
{
    return pow(e, pow(x, 2.f) * a * (-0.5f));
}







static inline float WaveHeightAfterImpulse(float distanceFromImpulse, float time, float impulseStrength = 1.f)
{
    float x = distanceFromImpulse;
    return (pow(e, -pow(x - time + 2, 2.f)) / (x + 1)) * cos(x) * impulseStrength;
}




static inline float IntensityOverDistanceFunc(c_float3& particlePos, c_float3& sourcePos, float slope=1, float maxIntensity=1)
{
    float x = VectorGetLength(particlePos - sourcePos);
    float a = slope;
    float b = sqrt(1.f / maxIntensity);

    return 1.f / pow(a * x + b, 2.f);

}

static inline float IntensityOverDistanceFuncGaussian(c_float3& particlePos, c_float3& sourcePos, float slope = 1, float maxIntensity = 1)
{
    float x = VectorGetLength(particlePos - sourcePos);
    float a = slope;
    float s = maxIntensity;
    return s * pow(e, -pow(a * x, 2.f));

}











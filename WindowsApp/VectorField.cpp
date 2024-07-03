#include "VectorField.h"

#include "InterpolationAndCurves.h"

void RVectorFieldVisualizer::Visualize(RGeometry& visualizer)
{
    ///Calculate min/max magnitudes of all vectors
    float magnitudeMax = VectorGetLength(VectorField.GetSlice(0).GetElement(uint2{ 0,0 }));
    float magnitudeMin = magnitudeMax;

    float divMax = Divergence[uint3{ 0,0,0 }];
    float divMin = divMax;

    float curlMax = VectorGetLength(Curl[uint3{ 0,0,0 }]);
    float curlMin = curlMax;

    auto numSlices = b3D ? (NumSamples.z - 1) : 1;

    for (uint z = 0; z < numSlices; z++)
    {
        for (uint x = 0; x < NumSamples.x; x++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                uint3 coord{ x,y,z };

                float curMagnitude = VectorGetLength(VectorField[coord]);
                magnitudeMax = std::max(magnitudeMax, curMagnitude);
                magnitudeMin = std::min(magnitudeMin, curMagnitude);

                float curDivergence = Divergence[coord];
                divMax = std::max(divMax, curDivergence);
                divMin = std::min(divMin, curDivergence);

                float curCurl = VectorGetLength(Curl[coord]);
                curlMax = std::max(curlMax, curCurl);
                curlMin = std::min(curlMin, curCurl);
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

                float3 curInput = InputSet[coord];
                float3 curVec = VectorField[coord];
                float curMagnitude = VectorGetLength(curVec);
                float curDivergence = Divergence[coord];
                float curCurl = VectorGetLength(Curl[coord]);

                ///Set Color
                {
                    float t = 0;

                    auto compute_t = [](float curValue, float minValue, float maxValue)->float
                    {
                        ///Compute t only if there is a significant range
                        if (fabs(maxValue - minValue) >= 0.01)
                        {
                            return MapTo01Range(curValue, minValue, maxValue);
                        }
                        return 0.5;
                    };

                    if (CurColorVisMode == EColorVisualizeMode::Divergence)
                    {
                        t = compute_t(curDivergence, divMin, divMax);
                    }
                    else if(CurColorVisMode == EColorVisualizeMode::Magnitude)
                    {
                        t = compute_t(curMagnitude, magnitudeMin, magnitudeMax);
                    }
                    else if (CurColorVisMode == EColorVisualizeMode::Curl)
                    {
                        t = compute_t(curCurl, curlMin, curlMax);
                    }

                    ///Get color based on cur magnitude
                    static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
                    static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
                    static const float3 RedColor = GetColorRGBFloat(EColor::Red);

                    float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

                    visualizer.SetColorRaw(CurColor);
                }


                ///Render direction arrow
                // visualizer.AddArrowPolyboardDirection
                float3 dir = VectorNormalize(curVec);
                visualizer.AddArrowPolyboard(curInput - dir * SamplingPeriodLength.x * 0.5, curInput + dir * SamplingPeriodLength.x * 0.5);


            }
        }
    }

    ///Point On A Func
    float3 curInput = InputSet[PointOnAFunc];
    visualizer.SetColor(EColor::White);
    visualizer.AddCircleCameraFacing(curInput, 0.05);

}


#include "InterpolationEtc/InterpolationCommon.hlsli"
//===========================================================================
//
//                              SAMPLING INTERPOLATION
//
//===========================================================================

//Decompose TexCoord into int and frac parts
void TransformTexCoords(float2 uvGlobal, uint2 texSize, float scale, uint2 offset, out uint2 outTexCoord, out float2 outUVLocal)
{
    float2 texelCoordsFloat = uvGlobal * texSize / scale;
    ///Determine Spline Segment
    uint2 segmentIndex = int2(texelCoordsFloat.xy);
    float2 uvLocal = frac(texelCoordsFloat);

    //Apply offset
    segmentIndex = clamp(segmentIndex + offset, 0, texSize - 1);

    outTexCoord = segmentIndex;
    outUVLocal = uvLocal;
}

float3 SampleNearestNeighbor(Texture2D<float3> outputSet, float2 uv, uint2 texSize, float scale = 1.f, uint2 offset = 0)
{
    float2 texelCoordsFloat = uv * texSize / scale;
    uint2 segmentIndex = round(texelCoordsFloat.xy);
    segmentIndex = clamp(segmentIndex + offset, 0, texSize - 1);
    return outputSet[segmentIndex];
}

float3 SampleBilinear(Texture2D<float3> outputSet, float2 uv, uint2 texSize, float scale = 1.f, uint2 offset = 0)
{
    uint2 texelCoord;
    float2 uvLocal;
    TransformTexCoords(uv, texSize, scale, offset, texelCoord, uvLocal);

    float3 s00 = outputSet[texelCoord];
    float3 s10 = outputSet[texelCoord + uint2(1, 0)];
    float3 s01 = outputSet[texelCoord + uint2(0, 1)];
    float3 s11 = outputSet[texelCoord + uint2(1, 1)];
    return BilinearInterpolate(uvLocal, s00, s10, s01, s11);
}

float3 SampleBilinearSmooth(Texture2D<float3> outputSet, float2 uv, uint2 texSize, float scale = 1.f, uint2 offset = 0)
{
    uint2 texelCoord;
    float2 uvLocal;
    TransformTexCoords(uv, texSize, scale, offset, texelCoord, uvLocal);

    float3 s00 = outputSet[texelCoord];
    float3 s10 = outputSet[texelCoord + uint2(1, 0)];
    float3 s01 = outputSet[texelCoord + uint2(0, 1)];
    float3 s11 = outputSet[texelCoord + uint2(1, 1)];

    uvLocal = float2(ConvertParameterToSmoothstep(uvLocal.x), ConvertParameterToSmoothstep(uvLocal.y));

    return BilinearInterpolate(uvLocal, s00, s10, s01, s11);
}


float3 GetDerivativeCatmullRom(float3 pPrev, float3 pNext, float t = 0.5f)
{
    return (pNext - pPrev) * (1.f - t);
}
float3 CubicCatmullRomInterpolateSetX(Texture2D<float3> outputSet, uint2 segmentIndex, float u, float tension)
{
    float3 sPrev = outputSet[segmentIndex + uint2(-1, 0)];
    float3 s0 = outputSet[segmentIndex + uint2(0, 0)];
    float3 s1 = outputSet[segmentIndex + uint2(1, 0)];
    float3 sNext = outputSet[segmentIndex + uint2(2, 0)];

    float3 d0 = GetDerivativeCatmullRom(sPrev, s1, tension);
    float3 d1 = GetDerivativeCatmullRom(s0, sNext, tension);

    return GetPointOnHermiteCurveParametric(u, s0, d0, s1, d1);
}
float3 SampleBicubicCatmullRom(Texture2D<float3> outputSet, uint2 setSize, float2 uv, float tension, float scale = 1.f, int2 offset = 0)
{
    uint2 sampleCoord;
    float2 uvLocal;
    TransformTexCoords(uv, setSize, scale, offset, sampleCoord, uvLocal);

    float3 HorizontalLinePrev = CubicCatmullRomInterpolateSetX(outputSet, sampleCoord + uint2(0, -1), uvLocal.x, tension);
    float3 HorizontalLine0 = CubicCatmullRomInterpolateSetX(outputSet, sampleCoord + uint2(0, 0), uvLocal.x, tension);
    float3 HorizontalLine1 = CubicCatmullRomInterpolateSetX(outputSet, sampleCoord + uint2(0, 1), uvLocal.x, tension);
    float3 HorizontalLineNext = CubicCatmullRomInterpolateSetX(outputSet, sampleCoord + uint2(0, 2), uvLocal.x, tension);

    float3 d0 = GetDerivativeCatmullRom(HorizontalLinePrev, HorizontalLine1, tension);
    float3 d1 = GetDerivativeCatmullRom(HorizontalLine0, HorizontalLineNext, tension);

    return GetPointOnHermiteCurveParametric(uvLocal.y, HorizontalLine0, d0, HorizontalLine1, d1);
}


float3 CubicLagrangeInterpolateSetX(Texture2D<float3> outputSet, uint2 segmentIndex, float u)
{
    float3 sPrev = outputSet[segmentIndex + uint2(-1, 0)];
    float3 s0 = outputSet[segmentIndex + uint2(0, 0)];
    float3 s1 = outputSet[segmentIndex + uint2(1, 0)];
    float3 sNext = outputSet[segmentIndex + uint2(2, 0)];
    return CubicInterpolate(u, sPrev, s0, s1, sNext, -1.f, 0.f, 1.f, 2.f);
    //return InterpolateCubicSplineSegmentLagrange(u, sPrev, s0, s1, sNext);
}
float3 SampleBicubicLagrange(Texture2D<float3> outputSet, uint2 setSize, float2 uv, float scale = 1.f, int2 offset = 0)
{
    uint2 sampleCoord;
    float2 uvLocal;
    TransformTexCoords(uv, setSize, scale, offset, sampleCoord, uvLocal);

    float3 HorizontalLinePrev = CubicLagrangeInterpolateSetX(outputSet, sampleCoord + uint2(0, -1), uvLocal.x);
    float3 HorizontalLine0 = CubicLagrangeInterpolateSetX(outputSet, sampleCoord + uint2(0, 0), uvLocal.x);
    float3 HorizontalLine1 = CubicLagrangeInterpolateSetX(outputSet, sampleCoord + uint2(0, 1), uvLocal.x);
    float3 HorizontalLineNext = CubicLagrangeInterpolateSetX(outputSet, sampleCoord + uint2(0, 2), uvLocal.x);

    return CubicInterpolate(uvLocal.y, HorizontalLinePrev, HorizontalLine0, HorizontalLine1, HorizontalLineNext, -1.f, 0.f, 1.f, 2.f);
    //return InterpolateCubicSplineSegmentLagrange(uvLocal.y, HorizontalLinePrev, HorizontalLine0, HorizontalLine1, HorizontalLineNext);
}



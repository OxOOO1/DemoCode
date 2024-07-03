
#include "InterpolationEtc/InterpolationCommon.hlsli"
#include "InterpolationEtc/InterpolationCommonCPP.h"

/* 
	CBuffers
 */
ConstantBuffer<GConstants> GConstantsCB : register(b0);

/* 
	SRV
 */
StructuredBuffer<float3> ControlPointsBufferSRV : register(t0);
//StructuredBuffer<float3> IntermediateOutputBufferSRV : register(t1);

/* 
	UAV
 */
RWStructuredBuffer<float3> IntermediateOutputBufferUAV : register(u0);
RWStructuredBuffer<float3> OutputBufferUAV : register(u1);

void ComputeDerivativesForPoint(float3 pPrev, float3 pCurrent, float3 pNext, out float3 vIn, out float3 vOut)
{
	const float3 t = GConstantsCB.TCBTension;
	const float3 c = GConstantsCB.TCBContinuity;
	const float3 b = GConstantsCB.TCBBias;

	float3 d0 = (pCurrent - pPrev);
	float3 d1 = (pNext - pCurrent);

	///TCB extension of Catmull-Rom
	vIn = d0 * ((1 - t)*(1 + b)*(1 - c))*0.5f + d1 * ((1 - t)*(1 - b)*(1 + c))*0.5f;
	vOut = d0 * ((1 - t)*(1 + b)*(1 + c))*0.5f + d1 * ((1 - t)*(1 - b)*(1 - c))*0.5f;
}

float3 ComputeDerivativesForPointIn(float3 pPrev, float3 pCurrent, float3 pNext)
{
	float3 vIn;
	float3 dummy;
	ComputeDerivativesForPoint(pPrev, pCurrent, pNext, vIn, dummy);
	return vIn;
}
float3 ComputeDerivativesForPointOut(float3 pPrev, float3 pCurrent, float3 pNext)
{
	float3 dummy;
	float3 vOut;
	ComputeDerivativesForPoint(pPrev, pCurrent, pNext, dummy, vOut);
	return vOut;
}

float3 GetPointOnSplineSegment(float t, float3 A, float3 B, float3 C, float3 D)
{
#if INTERPOLATION_METHOD == INTERPOLATION_METHOD_CR
	return InterpolateCubicSplineSegmentCatmullRom(t, A, B, C, D);
#elif INTERPOLATION_METHOD == INTERPOLATION_METHOD_LAGRANGE
	return InterpolateCubicSplineSegmentLagrange(t, A, B, C, D);
#elif INTERPOLATION_METHOD == INTERPOLATION_METHOD_TCB

	//compute derivatives
	float3 firstKnotVOut = ComputeDerivativesForPointOut(A, B, C);
	float3 secondKnotVIn = ComputeDerivativesForPointOut(B, C, D);

	return GetPointOnHermiteCurveParametric(t, B, firstKnotVOut, C, secondKnotVIn);

#else //Linear
	return LinearInterpolate(t, B, C);//DEBUG
#endif//INTERPOLATION_METHOD
	
	
	
}

uint GetControlPointLinearIndex(uint2 coord)
{
	return coord.y * GConstantsCB.NumControlPoints.x + coord.x;
}
uint GetOutputPointLinearIndex(uint2 coord)
{
	return coord.y * GConstantsCB.NumOuputPoints.x + coord.x; 
}


void GetControlPointsH(uint2 segmentCoord, out float3 point0, out float3 point1, out float3 point2, out float3 point3)
{
	point1 = ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord)];
	point2 = ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord + uint2(1, 0))];
	//edges
	point0 = segmentCoord.x > 0 ? ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord - uint2(1, 0))] : point1;
	point3 = (segmentCoord.x + 2) < GConstantsCB.NumControlPoints.x ? ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord + uint2(2, 0))] : point2;
}

[numthreads( 64, 1, 1 )]
void InterpolateHorizontal( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	
	const uint stackIndex = SampleCoord.y;

	if(SampleCoord.x >= GConstantsCB.NumOuputPoints.x)
	{
		return;
	}

	//map outputSampleIndex to Spline tGlobal
	float tGlobal = MapToRange(float(SampleCoord.x), 0.f, float(GConstantsCB.NumOuputPoints.x - 1), 0.f, float(GConstantsCB.NumControlPoints.x - 1));
	uint segmentIndex = int(tGlobal);
	float tLocal = tGlobal - segmentIndex;

	//load control points
	float3 point0;
	float3 point1;//center
	float3 point2;//center
	float3 point3;
	GetControlPointsH(uint2(segmentIndex, stackIndex), point0, point1, point2, point3);

	//Get point on spline
	float3 pointOnSplineSegment = GetPointOnSplineSegment(tLocal, point0, point1, point2, point3);

	IntermediateOutputBufferUAV[GetOutputPointLinearIndex(SampleCoord)] = pointOnSplineSegment;

}

void GetControlPointsV(uint2 segmentCoord, out float3 point0, out float3 point1, out float3 point2, out float3 point3)
{
	point1 = ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord)];
	point2 = ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord + uint2(0, 1))];
	//edges
	point0 = segmentCoord.y > 0 ? ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord - uint2(0, 1))] : point1;
	point3 = (segmentCoord.y + 2) < GConstantsCB.NumControlPoints.y ? ControlPointsBufferSRV[GetControlPointLinearIndex(segmentCoord + uint2(0, 2))] : point2;
}

//TODO:We can just copy intermediate buffer into Output, to avoid redundant interpolation

[numthreads( 1, 64, 1 )]
void InterpolateVertical( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	const uint columnIndex = SampleCoord.x;

	if((SampleCoord.y >= GConstantsCB.NumOuputPoints.y))
	{
		return;
	}

	float tGlobal = MapToRange(float(SampleCoord.y), 0.f, float(GConstantsCB.NumOuputPoints.y - 1), 0.f, float(GConstantsCB.NumControlPoints.y - 1));
	uint segmentIndex = int(tGlobal);
	float tLocal = tGlobal - segmentIndex;

	//TODO: If tLocal is 0, just output control point

	//load control points
	float3 point0;
	float3 point1;//center
	float3 point2;//center
	float3 point3;
	GetControlPointsV(uint2(columnIndex, segmentIndex), point0, point1, point2, point3);

	//Get point on spline
	float3 pointOnSplineSegment = GetPointOnSplineSegment(tLocal, point0, point1, point2, point3);

	OutputBufferUAV[GetOutputPointLinearIndex(SampleCoord)] = pointOnSplineSegment;

}
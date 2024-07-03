#pragma once

#include "GenerateSimpleGeometry.h"

class RSimpleGeometryPass : public RPass
{
public:

	void InitResources() override;
	void InitPSOs();

	void SetCommonStates(RCommandListGraphics& CmdList);

	void Execute();

	void SubmitGeometry(RGeometry& into);

	RTexture NormalTexture;

	void RenderUI();


	void RenderGPUGeneratedMesh(RStructuredBuffer* pVertexBufferGPU, RStructuredBuffer* pNormalsBufferGPU, uint numVertices, RStructuredBuffer* pTexCoordsBufferGPU, RTexture* pColorTexture, RTexture* pNormalTexture);
	void RenderInstancedMesh(const RDynamicVector<float3>& vertexBuffer, TStructuredBuffer<float3>& instancesPositions, uint numInstances, TStructuredBuffer<float3>* pInstancesColors = nullptr);

	void RenderImGeometry(RCommandListGraphics& CmdList, RGeometry& GImGeometry);


	struct RSamplingTimer
	{
		float SamplingRate = 30;//Samples Per Sec
		bool bCanSample = true;
		float TimeFromLastSample = 0;
		void UpdateTime()
		{
			float SampleDurationMs = 1000.0 / SamplingRate;

			auto timeCurrent = Time::Get().GetTimeAfterLaunchMs();
			float HowMuchTimeHavePassed = timeCurrent - TimeFromLastSample;

			if (HowMuchTimeHavePassed > SampleDurationMs)
			{
				bCanSample = true;
				TimeFromLastSample = timeCurrent;
			}
		}
	};





};


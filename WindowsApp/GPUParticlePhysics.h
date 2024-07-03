#pragma once

#include "MathGeneral.h"

namespace GpuParticlesPhysicsPass
{
	struct OutputStruct
	{
		TStructuredBuffer<float3>* pPositionsBuffer;
		TStructuredBuffer<float3>* pColorsBuffer;
		uint NumInstances;
		float InstanceRadius;
		RStructuredBuffer* MeshBufferVerticesGPU;
		RStructuredBuffer* MeshBufferNormalsGPU = nullptr;
		RStructuredBuffer* MeshBufferTexCordsGPU = nullptr;
		RTexture* pColorTexture = nullptr;
		RTexture* pNormalTexture = nullptr;
		uint NumTriangles=0;
	};

	void Execute(OutputStruct& outResult);


}


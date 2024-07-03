#pragma once

#include "MathGeneral.h"

#include "shaders/MeshGenerator/MeshGeneratorCommonStructs.h"

namespace GPUMeshGeneration
{

	///In Compute Shader, generate mesh and write into buffers
	///Bind those buffers and raster

	class GeometryGeneratorRootSig : public RRootSignature
	{
	public:
		GeometryGeneratorRootSig()
		{
			AddConstantBuffer("GlobalInputs", 0);
			AddBufferUAV("outVertexBuffer", 0);
			AddBufferUAV("outColorsBuffer", 1);

			Finalize(L"GeometryGeneratorRootSig");
		}
	};

	class GeometryGenerationShader : public RShader
	{
	public:
		GeometryGenerationShader()
		{
			Create("MeshGenerator/MeshGenerator.hlsl", "GenerateMeshMain", L"cs_6_0");
		}
	};

	
	///Constants
	static int2 NumSamples{ 32,32 };
	static uint BufferSizeBytes;

	static float3 MeshColor{1,0,0};

	///CBuffer
	
	static TCBuffer<CFuncInputs> CBuffer;

	///GPU Resources
	static TStructuredBuffer<float3> GeneratedVertexBuffer;
	static TStructuredBuffer<float3> GeneratedColorsBuffer;

	void RenderUI()
	{
		if (ImGui::Begin("GPU Mesh Generator"))
		{
			int exp2 = RMath::Log2(NumSamples.x);
			if (ImGui::SliderInt("Num Samples Pow 2", &exp2, 3, 10))
			{
				NumSamples.x = NumSamples.y = ipow(2, exp2);
			}

			ImGui::ColorEdit3("Mesh Color", &MeshColor.x);

		}
		ImGui::End();
	}

	void GenerateMesh()
	{

		PROFILE_FUNC();

		///
		///	CBuffer
		{
			if (CBuffer.GetBufferSizeBytes() == 0)
			{
				CBuffer.Create(TEXT("GlobalInputs"));
			}

			///Upload Constants
			CFuncInputs CBufferInputs;
			CBufferInputs.GNumThreads = uint2(NumSamples.x, NumSamples.y);
			CBufferInputs.MeshColor = MeshColor;

			CBuffer.Update(CBufferInputs);
		}

		///
		///	GPU Resources
		///
		
		BufferSizeBytes = NumSamples.y * NumSamples.x * 6;

		if (GeneratedVertexBuffer.GetElementCount() < BufferSizeBytes)
		{
			GeneratedVertexBuffer.Create(TEXT("GPUGeneratedVertexBuffer"), BufferSizeBytes);
		}
		if (GeneratedColorsBuffer.GetElementCount() < BufferSizeBytes)
		{
			GeneratedColorsBuffer.Create(TEXT("GPUGeneratedColorsBuffer"), BufferSizeBytes);
		}
		

		///
		///	Shader & RootSig
		///
		static GeometryGeneratorRootSig RootSig;
		static GeometryGenerationShader ComputeShader;


		RCommandListCompute& CmdList = RCommandListCompute::BeginNew();
		{
			GPU_PROFILE_SCOPE(MeshGenerator, CmdList);

			CmdList.SetRootSignature(RootSig);

			CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("GlobalInputs"), CBuffer.GetGpuVirtualAddress());

			CmdList.TransitionResource(GeneratedVertexBuffer, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(GeneratedColorsBuffer, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

			CmdList.SetBufferAsUAV(RootSig.GetRootParamIndex("outVertexBuffer"), GeneratedVertexBuffer);
			CmdList.SetBufferAsUAV(RootSig.GetRootParamIndex("outColorsBuffer"), GeneratedColorsBuffer);

			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;

			CmdList.Dispatch(NumSamples.x / ThreadGroupSize, NumSamples.y / ThreadGroupSize, 1);
		}
		CmdList.ExecuteCmdListAndReleaseContext(false);

	}

	


}

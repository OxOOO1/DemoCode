#pragma once

#include "MathGeneral.h"

namespace GPUFunctions
{








	///Func Visualizer

	class FuncVisualizerRootSig : public RRootSignature
	{
	public:
		FuncVisualizerRootSig()
		{
			AddConstantBuffer("FuncInputs", 0);
			AddTextureUAV("OutputTex", 0);

			Finalize(L"RSimpleGeometryPassRootSignature");
		}
	};

	class RFuncVisualizerShader : public RShader
	{
	public:
		RFuncVisualizerShader()
		{
			Create("GPUFunctions.hlsl", "FuncVisualizer", L"cs_6_0");

		}
	};


	void FuncVisualizer(uint2 numSamples, float2 samplingPeriodLength, float2 inputStart)
	{
		///
		///	CBuffer
		///
		struct alignas(16) CFuncInputs
		{
			float2 SamplingPeriodLength;
			float2 InputStart;
		};

		static TCBuffer<CFuncInputs> CBuffer;

		if (CBuffer.GetBufferSizeBytes() == 0)
		{
			CBuffer.Create(TEXT("FuncInputs"));
		}

		///Upload Constants
		CFuncInputs CBufferInputs;
		CBufferInputs.InputStart = inputStart;
		CBufferInputs.SamplingPeriodLength = samplingPeriodLength;

		CBuffer.Update(CBufferInputs);

		///
		///	UAV Texture
		///
		static RTexture OutputTexture;

		if (OutputTexture.Size.x < numSamples.x || OutputTexture.Size.y < numSamples.y)
		{
			OutputTexture.Create(numSamples.x, numSamples.y, DXGI_FORMAT_R32_FLOAT, true);
		}
		
		///
		///	Shader & RootSig
		///
		static FuncVisualizerRootSig RootSig;
		static RFuncVisualizerShader ComputeShader;


		RCommandListCompute& CmdList = RCommandListCompute::BeginNew(L"FuncVisualizer");

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("FuncInputs"), CBuffer.GetGpuVirtualAddress());

		CmdList.TransitionResource(OutputTexture, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.SetDynamicDescriptor(RootSig.GetRootParamIndex("OutputTex"), 0, OutputTexture.GetUnorderedAccessView());

		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;

		CmdList.Dispatch(numSamples.x / ThreadGroupSize, numSamples.y / ThreadGroupSize, 1);

		CmdList.ExecuteCmdListAndReleaseContext(false);

	}


}
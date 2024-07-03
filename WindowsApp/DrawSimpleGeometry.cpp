#include "DrawSimpleGeometry.h"

#include "LLRenderer.h"
#include "Core/Renderer/Renderer.h"


#include "GPUGenerateGeometry.h"
#include "GPUParticlePhysics.h"

#include "Burner/BurnerMain.h"


bool GbFirstBoot = true;


void* __cdecl operator new[](size_t size, const char* name, int flags, unsigned debugFlags, const char* file, int line)
{
	return new uint8_t[size];
}

void* __cdecl operator new[](size_t size, size_t, size_t, char const*, int, unsigned int, char const*, int)
{
	return new uint8_t[size];
}


using namespace RCommandListHelper;

class RSimpleGeometryPassRootSignature : public RRootSignatureExt, public RSingleton<RSimpleGeometryPassRootSignature>
{
public:
	RSimpleGeometryPassRootSignature()
	{
		//Constants
		AddConstantBuffer("SceneUniformBufferCBuffer", 0);
		AddConstantBuffer("FuncInputs", 1);
		AddConstants("FlowControlConstants", 2, 1);

		//SRVs
		AddSRV("VertexBufferPosition", 0);
		AddSRV("ColorsBuffer", 1);
		AddSRV("NormalsBuffer", 2);
		AddSRV("InstancesPositionsBuffer", 3);
		AddSRV("TexCoordBuffer", 4);
		AddSRV("ColorTexture", 5);
		AddSRV("NormalTexture", 6);

		AddSRV("InstancesColorsBuffer", 7);


		//UAVs
		//AddUAV("AtomicOffsetBufferUAV", 4);

		//Samplers
		AddStaticSampler(0, RCommonResources::Get().SamplerLinearClampDesc);


		Finalize(L"RSimpleGeometryPassRootSignature");
	}
};



template <bool Normals, bool Instance, bool bSingleColorForDraw = false, bool bTexture = false>
class RSimpleGeometryPassVertexShader : public RShader
{
public:

	RSimpleGeometryPassVertexShader()
	{
		if constexpr (Normals)
		{
			PushBackDefine({ L"NORMALS", L"1" });
		}
		if constexpr (Instance)
		{
			PushBackDefine({ L"INSTANCED_DRAW", L"1" });
		}
		if constexpr (bSingleColorForDraw)
		{
			PushBackDefine({ L"SINGLE_COLOR_PER_DRAW", L"1" });
		}
		if constexpr (bTexture)
		{
			PushBackDefine({ L"USES_TEXTURE", L"1" });
		}
		CreateVS("DrawSimpleGeometry.hlsl", "MainVS", L"vs_6_0");
	}
};

template <bool Normals, bool bTexture = false>
class RSimpleGeometryPassPixelShader : public RShader
{
public:
	RSimpleGeometryPassPixelShader()
	{
		if constexpr (Normals)
		{
			PushBackDefine({ L"NORMALS", L"1" });
		}
		if constexpr (bTexture)
		{
			PushBackDefine({ L"USES_TEXTURE", L"1" });
		}
		CreatePS("DrawSimpleGeometry.hlsl", "MainPS", L"ps_6_0");
	}
};

RSimpleGeometryPassVertexShader<false, false> VertexShader;
RSimpleGeometryPassVertexShader<true, false> VertexShaderNormals;
RSimpleGeometryPassVertexShader<false, true> VertexShaderInstanced;
RSimpleGeometryPassVertexShader<true, true> VertexShaderNormalsInstanced;

RSimpleGeometryPassVertexShader<false, false, true> VertexShaderSingleColor;
RSimpleGeometryPassVertexShader<true, false, true> VertexShaderNormalsSingleColor;

RSimpleGeometryPassVertexShader<false, false, false, true> VertexShaderTexture;
RSimpleGeometryPassVertexShader<true, false, false, true> VertexShaderNormalsTexture;

RSimpleGeometryPassPixelShader<false> PixelShader;
RSimpleGeometryPassPixelShader<true> PixelShaderNormals;
RSimpleGeometryPassPixelShader<false, true> PixelShaderTexture;
RSimpleGeometryPassPixelShader<true, true> PixelShaderNormalsTexture;

void RSimpleGeometryPass::Execute()
{

	PROFILE_FUNC();

	auto& GImGeometry = RGeometry::GInstance;
	GImGeometry.BeginGeometrySubmission();


	SubmitGeometry(GImGeometry);

	//Render
	/// *************************************************** Render ******************************************************
	/// *************************************************** Render ******************************************************
	/// *************************************************** Render ******************************************************


	RenderUI();


#if 0
	///GPU Mesh Generator
	//GPUMeshGeneration::RenderUI();
	GPUMeshGeneration::GenerateMesh();

	//RenderGPUGeneratedMesh();
#endif
	
#if 1///GPU Particle Physics

	GpuParticlesPhysicsPass::OutputStruct ParticlesResultStruct;
	GpuParticlesPhysicsPass::Execute(ParticlesResultStruct);

	///Generate Circle Mesh
	static RDynamicVector<float3> circleMesh;
	GImGeometry.GenerateCircleCameraFacing(circleMesh, float3{}, ParticlesResultStruct.InstanceRadius);


	if (ParticlesResultStruct.NumTriangles > 0)
	{
		RenderGPUGeneratedMesh(ParticlesResultStruct.MeshBufferVerticesGPU, ParticlesResultStruct.MeshBufferNormalsGPU, ParticlesResultStruct.NumTriangles, ParticlesResultStruct.MeshBufferTexCordsGPU, ParticlesResultStruct.pColorTexture, ParticlesResultStruct.pNormalTexture);
	}
	else
	{

		//RenderInstancedMesh(circleMesh, *ParticlesResultStruct.pPositionsBuffer, ParticlesResultStruct.NumInstances, ParticlesResultStruct.pColorsBuffer);

	}


#endif

#if 0//BurnerParticles Particles
	
	BurnerMain();


#endif



#if 1 ///Immediate Geometry Render
	{
		RCommandListGraphics& CmdList = RCommandList::BeginNew().GetGraphicsContext();

		///Upload all Geometry to GPU
		RGeometry::GInstance.EndGeometrySubmission(CmdList);
		RenderImGeometry(CmdList, GImGeometry);

		//Execute
		CmdList.ExecuteCmdListAndReleaseContext();
	}
	
#endif



	GbFirstBoot = false;
}



void RSimpleGeometryPass::InitPSOs()
{

}

void RSimpleGeometryPass::InitResources()
{
	InitPSOs();

	NormalTexture.CreateFromBitmap(RBitmap::FromFile("brick_normal.png"), DXGI_FORMAT_R8G8B8A8_UNORM);

}

void RSimpleGeometryPass::SetCommonStates(RCommandListGraphics& CmdList)
{
	//Set global states 
	CmdList.SetViewport(RSceneRenderer::Get().GetViewport());
	CmdList.SetScissor(RSceneRenderer::Get().GetScissorRect());

	auto& RT = RCommonResources::Get().SceneColor;
	auto& DSV = RCommonResources::Get().SceneDepth;

	CmdList.SetRenderTarget(RT.GetRTV(), DSV.GetDSV());

	CmdList.SetRootSignature(RSimpleGeometryPassRootSignature::Get());

	CmdList.SetConstantBuffer(RSimpleGeometryPassRootSignature::Get().GetRootParamIndex("SceneUniformBufferCBuffer"), RScene::Get().UniformBufferParameters.GetGlobalUniformBuffer()->ConstantBuffer.GetResource()->GetGPUVirtualAddress());
}


void RenderFPS(bool animate = true)
{
	if (ImGui::Begin("FPS"))
	{
		static float values[30] = { 0 };
		static int values_offset = 0;
		static double refresh_time = 0.0;

		values[values_offset] = Time::Get().GetDeltaTimeMs();
		values_offset = (values_offset + 1) % IM_ARRAYSIZE(values);

		{
			float average = 0.0f;
			for (int n = 0; n < IM_ARRAYSIZE(values); n++)
			{
				average += values[n];
			}
			average /= (float)IM_ARRAYSIZE(values);
			char overlay[32];
			sprintf_s(overlay, "avg %f ms", average);
			ImGui::PlotLines("", values, IM_ARRAYSIZE(values), values_offset, overlay, 0.0f, 120.0f, ImVec2(0, 80.0f));
		}
	}
	ImGui::End();
}


void RSimpleGeometryPass::RenderUI()
{

	RenderAllGlobals();

	if (ImGui::Begin("MathGlobals"))
	{

	}
	ImGui::End();

	RenderFPS();

	RCPUProfileEvent::ShowProfileUI();

}


void RSimpleGeometryPass::RenderImGeometry(RCommandListGraphics& CmdList, RGeometry& GImGeometry)
{
	{
		GPU_PROFILE_SCOPE(RenderImGeometry, CmdList);

		SetCommonStates(CmdList);

		auto& RootSignature = RSimpleGeometryPassRootSignature::Get();

		///
		///	CBuffer
		///
		struct alignas(16) CFuncInputs
		{
			float3 GPreProjectionScale;
			float pad;
			float2 ViewSpaceClamp;
		};

		static TCBuffer<CFuncInputs> CBuffer;

		if (CBuffer.GetBufferSizeBytes() == 0)
		{
			CBuffer.Create(TEXT("FuncInputs"));
		}

		CBuffer.Update(CFuncInputs{});

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("FuncInputs"), CBuffer.GetGpuVirtualAddress());

		///Draw All Geometry

		auto BindPSOFromDesc = [&](RGeometry::RenderingResources& resources, RGeometry::PipelineStateDesc psoDesc)
		{
			RPSOFactory::CommonPipelineStates psoStates;

			if (psoDesc.twoSided)
			{
				psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;
			}
			else
			{
				psoStates.RasterDesc = RCommonResources::Get().RasterizerDefault;
			}
			if (psoDesc.depthEnabled)
			{
				psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;
			}
			else
			{
				psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;
			}

			///Shaders
			RShader* pVertexShader;
			RShader* pPixelShader;
			if (psoDesc.usesNormals)
			{
				pVertexShader = &VertexShaderNormals;
				pPixelShader = &PixelShaderNormals;
			}
			else if (psoDesc.usesTexture)
			{
				pVertexShader = &VertexShaderTexture;
				pPixelShader = &PixelShaderTexture;
			}
			else
			{
				pVertexShader = &VertexShader;
				pPixelShader = &PixelShader;
			}

			psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
			psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, *pVertexShader, pPixelShader, psoStates);

			///We can destroy created PSO since D3D12 PSO is stored in map
			CmdList.SetPipelineState(PSO);
		};

		auto RenderCurResources = [&](RGeometry::RenderingResources& resources, RGeometry::PipelineStateDesc psoDesc)
		{
			if (resources.VertexBuffer.size() > 0)
			{
				///Bind PSO from resource
				BindPSOFromDesc(resources, psoDesc);

				bool bHasNormalTexture = psoDesc.usesNormalTexture;
				CmdList.SetConstantArray(RootSignature.GetRootParamIndex("FlowControlConstants"), 1, &bHasNormalTexture);

				CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

				if (psoDesc.usesNormals)
				{
					SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalsBuffer"), resources.GpuBufferNormals, true);
				}

				if (psoDesc.usesTexture)
				{
					SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("TexCoordBuffer"), resources.GpuBufferTexCoords, true);
					SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), *resources.pColorTexture, true);

					if (psoDesc.usesNormalTexture)
					{
						SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), *resources.pNormalTexture, true);
					}
				}
				else
				{
					SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorsBuffer"), resources.GpuBufferColors, true);
				}
				

				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), resources.GpuBufferVertices, true);
				
				CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("FuncInputs"), CBuffer.GetGpuVirtualAddress());

				CmdList.Draw(resources.VertexBuffer.size());
			}
		};

		for (auto& [key, resource] : GImGeometry.PSOStateToResourcesMap)
		{
			RGeometry::PipelineStateDesc desc = RGeometry::GetDescFromHash(key);
			RenderCurResources(resource, desc);
		}


		//Finish Draw
		//CmdList.TransitionResource(RCommonResources::Get().SceneColor, D3D12_RESOURCE_STATE_PRESENT);
	}
}

void RSimpleGeometryPass::RenderGPUGeneratedMesh(RStructuredBuffer* pVertexBufferGPU, RStructuredBuffer* pNormalsBufferGPU, uint numVertices, RStructuredBuffer* pTexCoordsBufferGPU, RTexture* pColorTexture, RTexture* pNormalTexture)
{

	RCommandListGraphics& CmdList = RCommandList::BeginNew().GetGraphicsContext();
	{
		
		GPU_PROFILE_SCOPE(RenderGPUGeneratedMesh, CmdList);

		SetCommonStates(CmdList);

		auto& RootSignature = RSimpleGeometryPassRootSignature::Get();

		/*
		*	CBuffer
		*/
		struct alignas(16) CFuncInputs
		{
			float3 GPreProjectionScale;
			float pad;

			float2 ViewSpaceClamp;
			float2 pad2;

			float3 ColorPerDraw;
		};

		static TCBuffer<CFuncInputs> CBuffer2;

		if (CBuffer2.GetBufferSizeBytes() == 0)
		{
			CBuffer2.Create(TEXT("FuncInputs"));
		}

		CFuncInputs CBufferInputsCPU;
		static float3 DrawColor{ GenerateRandomColor() };
		CBufferInputsCPU.ColorPerDraw = DrawColor;
		CBuffer2.Update(CBufferInputsCPU);

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex(TEXT("FuncInputs")), CBuffer2.GetGpuVirtualAddress());

		bool bHasNormalTexture = pNormalTexture != nullptr;

		CmdList.SetConstantArray(RootSignature.GetRootParamIndex("FlowControlConstants"), 1, &bHasNormalTexture);

		///Draw Geometry

		///Bind PSO from resource
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;

		psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;
		//psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;

		///Shaders
		RShader* pVertexShader;
		RShader* pPixelShader;

		bool bHasNormals = pNormalsBufferGPU != nullptr;
		if (bHasNormals)
		{
			if (pColorTexture)
			{
				pVertexShader = &VertexShaderNormalsTexture;
				pPixelShader = &PixelShaderNormalsTexture;
			}
			else
			{
				pVertexShader = &VertexShaderNormalsSingleColor;
				pPixelShader = &PixelShaderNormals;
			}
			
		}
		else
		{
			pVertexShader = &VertexShaderSingleColor;
			pPixelShader = &PixelShader;
		}

		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

		auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, *pVertexShader, pPixelShader, psoStates);

		///We can destroy created PSO since D3D12 PSO is stored in map
		CmdList.SetPipelineState(PSO);

		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("FuncInputs"), CBuffer2.GetGpuVirtualAddress());

		if (bHasNormals)
		{
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalsBuffer"), *pNormalsBufferGPU, true);
		}
		if (pColorTexture)
		{
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("TexCoordBuffer"), *pTexCoordsBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), *pColorTexture, true);

			if (bHasNormalTexture)
			{
				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), *pNormalTexture, true);
			}
		}

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorsBuffer"), GPUMeshGeneration::GeneratedColorsBuffer, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), *pVertexBufferGPU, true);

		CmdList.Draw(numVertices);

	}
	//Execute
	CmdList.ExecuteCmdListAndReleaseContext();



}

void RSimpleGeometryPass::RenderInstancedMesh(const RDynamicVector<float3>& vertexBuffer, TStructuredBuffer<float3>& instancesPositions, uint numInstances, TStructuredBuffer<float3>* pInstancesColors)
{
	PROFILE_FUNC();

	RCommandListGraphics& CmdList = RCommandList::BeginNew().GetGraphicsContext();
	{
		GPU_PROFILE_SCOPE(RenderInstancedMesh, CmdList);

		static TStructuredBuffer<float3> VertexBufferGPU;

		if (VertexBufferGPU.GetElementCount() != vertexBuffer.size())
		{
			VertexBufferGPU.Create(TEXT("VertexBuffer"), vertexBuffer.size());
		}

		VertexBufferGPU.UploadElements(0, vertexBuffer.size(), vertexBuffer.data());

		SetCommonStates(CmdList);
		auto& RootSignature = RSimpleGeometryPassRootSignature::Get();

		///
		///	CBuffer
		///
		struct alignas(16) CFuncInputs
		{
			float3 GPreProjectionScale;
			float pad;
			float2 ViewSpaceClamp;
		};

		static TCBuffer<CFuncInputs> CBuffer2;

		if (CBuffer2.GetBufferSizeBytes() == 0)
		{
			CBuffer2.Create(TEXT("FuncInputs"));
		}

		CBuffer2.Update(CFuncInputs{});

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex(TEXT("FuncInputs")), CBuffer2.GetGpuVirtualAddress());

		///Draw Geometry

		///Bind PSO from resource
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;

		psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;
		//psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;

		///Shaders
		RShader* pVertexShader;
		RShader* pPixelShader;
		pVertexShader = &VertexShaderInstanced;
		pPixelShader = &PixelShader;


		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

		auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, *pVertexShader, pPixelShader, psoStates);

		///We can destroy created PSO since D3D12 PSO is stored in map
		CmdList.SetPipelineState(PSO);

		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), VertexBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorsBuffer"), GPUMeshGeneration::GeneratedColorsBuffer, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("InstancesPositionsBuffer"), instancesPositions, true);

		check(pInstancesColors != nullptr);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("InstancesColorsBuffer"), *pInstancesColors, true);
		
		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("FuncInputs"), CBuffer2.GetGpuVirtualAddress());

		CmdList.DrawInstanced(vertexBuffer.size(), numInstances);


		//Finish Draw
		//CmdList.TransitionResource(RCommonResources::Get().SceneColor, D3D12_RESOURCE_STATE_PRESENT);
	}
	//Execute
	CmdList.ExecuteCmdListAndReleaseContext();

}


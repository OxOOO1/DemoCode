
#include "GPUParticlePhysics.h"

#include "LLRenderer.h"

#include <algorithm>

#include "Intersections.h"
#include "LinearAlgebra.h"

#include "shaders/Physics/ParticlePhysicsCommonStructs.h"

#include "Renderer/Passes/ParallelSort.h"
#include "Renderer/Passes/Blur.h"
#include "Renderer/Passes/GPUSurface.h"

#include "InterpolationAndCurves.h"

#include "GenerateSimpleGeometry.h"
#include "Renderer/Passes/ParallelReduction.h"

#include "Core/Renderer/Renderer.h"

using namespace RCommandListHelper;




class CopyTextureRootSignatureGlobal : public RRootSignatureExt, public RSingleton<CopyTextureRootSignatureGlobal>
{
public:
	CopyTextureRootSignatureGlobal()
	{
		AddConstants("ConstArr", 0, 2);

		AddSRV("Source", 0);
		AddUAV("Destination", 0);

		AddStaticSampler(0, RCommonResources::Get().SamplerLinearClampDesc);

		Finalize(L"CopyTextureRootSignatureGlobal");
	}
};


class CopyTextureResourceShader : public RShader
{
public:
	CopyTextureResourceShader()
	{
		Create("GPUFunctions.hlsl", "CopyTextureResource", L"cs_6_0");
	}
};

void CopyTextureResourceHelper(RCommandListCompute& CmdList, RTexture& dest, RTexture& src)
{
	PROFILE_FUNC();
	GPU_PROFILE_SCOPE(CopyTextureHelper, CmdList);

	auto& RootSig = CopyTextureRootSignatureGlobal::Get();

	CmdList.SetRootSignature(RootSig);

	RCommandListHelper::SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("Source"), src, true);

	RCommandListHelper::SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("Destination"), dest, true);

	CmdList.SetConstantArray(RootSig.GetRootParamIndex("ConstArr"), 2, &dest.Size.x);

	static CopyTextureResourceShader ComputeShader;
	CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

	constexpr uint ThreadGroupSize = 8;
	uint2 numGroupsToDispatch{ RMath::DivideAndRoundUp(dest.Size.x, ThreadGroupSize), RMath::DivideAndRoundUp(dest.Size.y, ThreadGroupSize) };
	CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);

}


class ParticlePhysicsRootSig : public RRootSignatureExt, public RSingleton<ParticlePhysicsRootSig>
{
public:
	ParticlePhysicsRootSig()
	{
		///Constants
		AddConstantBuffer("CGConstantsCB", 0);
		AddConstants("GCollisionConstantsCB", 1, 3);
		AddConstantBuffer("CGSpatialControlConstantsCB", 2);
		AddConstants("GMeshGenerationConstantsCB", 3, 2);
		AddConstants("VectorProjectionPerPassConstantsCB", 4, 1);

		///SRVs
		AddSRV("PlanesBufferSRV", 0);
		AddSRV("ParticlesBufferPositionCurSRV", 1);
		AddSRV("CellIdOffsetAndCountBufferSRV", 2);
		AddSRV("ObjectListIndexBufferSRV", 3);
		AddSRV("CellListIndexBufferSRV", 4);
		AddSRV("PerObjectControlBitsBufferSRV", 5);
		AddSRV("ConstraintsBufferSRV", 6);
		AddSRV("ParticlesStateBufferSRV", 7);
		AddSRV("ConstraintsStatesBufferSRV", 8);

		AddSRV("ParticlesToRenderCoordListSRV", 9);

		AddSRV("ConstraintsStatesBufferBLUESRV", 10);
		AddSRV("ConstraintsStatesBufferGREENSRV", 11);
		AddSRV("ConstraintsStatesBufferYELLOWSRV", 12);
		AddSRV("ConstraintsStatesBufferREDSRV", 13);

		AddSRV("SphereCapsulesBufferSRV", 14);

		AddSRV("ParticlesBufferPositionAfterCreationSRV", 15);

		AddSRV("VelocityField_U_ComponentTextureSRV", 16);
		AddSRV("VelocityField_V_ComponentTextureSRV", 17);

		AddSRV("CurParticlesVelocityBufferSRV", 18);
		
		AddSRV("VelocityFieldStateTextureSRV", 19);

		AddSRV("DensityFieldTextureSRV", 20);

		AddSRV("PressureFieldTextureSRV", 21);

		AddSRV("VectorFieldDivergenceTextureSRV", 22);

		AddSRV("VelocityField_U_ComponentTextureNewSRV", 23);
		AddSRV("VelocityField_V_ComponentTextureNewSRV", 24);

		AddSRV("VectorFieldCurlTextureSRV", 25);

		AddSRV("DensityFieldTextureNewSRV", 26);

		AddSRV("DensitySourceColorTextureSRV", 27);

		AddSRV("InputTextureWithNegativeValuesSRV", 28);

		AddSRV("CurlGradientTextureSRV", 29);

		AddSRV("CurlGradientMinMaxTextureSRV", 30);

		AddSRV("DensityFieldMacCormackAdvectedTextureSRV", 31);
		AddSRV("DensityFieldMacCormackReversedTextureSRV", 32);

		AddSRV("VelocityField_U_ComponentWeightTextureSRV", 33);
		AddSRV("VelocityField_V_ComponentWeightTextureSRV", 34);

		AddSRV("AdvectedParticleDensityTextureSRV", 20);

		///UAVs
		AddUAV("ParticlesBufferPositionCurUAV", 0);
		AddUAV("ParticlesBufferPositionPrevUAV", 1);
		AddUAV("NumCellsObjectIntersectsUAV", 2);
		AddUAV("PerObjectControlBitsBufferUAV", 3);
		AddUAV("AtomicOffsetBufferUAV", 4);
		AddUAV("ObjectListIndexBufferUAV", 5);
		AddUAV("CellListIndexBufferUAV", 6);
		AddUAV("ParticlesStateBufferUAV", 7);

		AddUAV("IntersectedPointIndexUAV", 8);
		AddUAV("ClickedObjectViewSpaceZUAV", 9);

		AddUAV("ConstraintsStatesBufferUAV", 10);

		AddUAV("MeshBufferVerticesUAV", 11);
		AddUAV("MeshBufferNormalsUAV", 12);
		AddUAV("MeshBufferCounterUAV", 13);

		AddUAV("InterlockedAddBuffer", 14);

		AddUAV("CurParticlesVelocityBufferUAV", 15);

		AddUAV("VelocityField_U_ComponentTextureUAV", 16);
		AddUAV("VelocityField_V_ComponentTextureUAV", 17);

		AddUAV("VectorFieldDivergenceTextureUAV", 18);

		AddUAV("DensityFieldTextureUAV", 19);

		AddUAV("PressureFieldTextureUAV", 20);

		AddUAV("VectorFieldCurlTextureUAV", 21);

		AddUAV("MeshBufferTexCoordsUAV", 22);

		AddUAV("VelocityFieldStateTextureUAV", 23);

		AddUAV("VisualizerColorTextureUAV", 24);

		AddUAV("CurlGradientTextureUAV", 25);

		AddUAV("CurlGradientDirectionTextureUAV", 26);

		AddUAV("PressureGradientTextureUAV", 27);

		AddUAV("VelocityField_U_ComponentWeightTextureUAV", 28);
		AddUAV("VelocityField_V_ComponentWeightTextureUAV", 29);

		AddUAV("AdvectedParticleDensityTextureUAV", 30);

		AddUAV("AverageDensityBufferUAV", 31);

		AddUAV("ParticlesColorsBufferUAV", 20);

		AddUAV("CurParticlesExternalForceBufferUAV", 9);

		AddUAV("PrevParticlesVelocityBufferUAV", 22);


		//Samplers
		AddStaticSampler(0, RCommonResources::Get().SamplerLinearClampDesc);
		AddStaticSampler(1, RCommonResources::Get().SamplerPointClampDesc);

		Finalize(L"ParticlePhysicsRootSig");
	}
};


class GPUParticlePhysics
{
public:

	///***************************
	///		PHYSICS OBJECTS SoA
	///***************************
	TStructuredBuffer<float3> ParticlesBufferPositionCurGPU;
	TStructuredBuffer<float3> ParticlesBufferPositionPrevGPU;

	TStructuredBuffer<float3> CurParticlesVelocityBufferGPU;
	TStructuredBuffer<float3> PrevParticlesVelocityBufferGPU;

	TStructuredBuffer<float3> CurParticlesExternalForceBufferGPU;

	TStructuredBuffer<float3> ParticlesBufferPositionAfterCreationGPU;//particles might have a constraint to default pos

	TStructuredBuffer<uint> ParticlesStateBufferGPU;

	TStructuredBuffer<float3> ParticlesColorsBufferGPU;

	//For Particle Advection on CPU
	RDynamicVector<float3> ParticlesBufferPositionCurCPU;
	RDynamicVector<float3> CurParticlesVelocityBufferCPU;

	RDynamicVector<float3> ParticlesColorsBufferCPU;

	///******************************
	///		SIMULATION PARAMETERS
	///******************************


	///Spawning Objects
	static constexpr uint MaxNumParticles = 1000000;
	bool bResetSimulation = true;
	bool bResetParticles = true;

	bool bSimulationRunning = true;
	bool bSimulateOnce = false;

	bool bSpawnObjectOnce = false;

	int DeltaTimeScale = 4;

	float RandomImpulseScale = 50.f;

	///CBuffer Constants
	CGConstants CGConstantsCB;
	TCBuffer<CGConstants> CGConstantsCBufferGPU;

	uint NumCollisionResolveIterations = 1;
	uint NumGlobalCollisionResolveIterations = 1;
	uint NumConstraintResolveIterations = 1;

	bool bGenerateCloth = false;
	//bool bGenerateCloth = true;

	bool bEnableSelfCollisions = true;
	bool bEnableConstraints = false;
	bool bDebugVisualizeConstraints = false;
	bool bDebugVisualizeSimulationBounds = true;

	float GetCollisionGridLengthMax()
	{
		return GetCollisionGridCellSize() * 1024.f;///1024 is coming from "10 bits per Grid coord"
	}

	float GetCollisionGridCellSize()
	{
		return CGConstantsCB.PointRadius * CellSizeScale;
	}

	void RenderUI()
	{
		if (ImGui::Begin("GPU Particle Physics"))
		{

			bResetSimulation = ImGui::Button("Reset");
			bResetParticles = ImGui::Button("Reset Particles");
			if (Keyboard::IsKeyPressed(VK_SHIFT) && Keyboard::IsKeyPressed('R'))
			{
				bResetSimulation = true;
			}
			else if(Keyboard::IsKeyPressed('R'))
			{
				bResetParticles = true;
			}
			
			ImGui::Checkbox("DIFFUSE", &bNeablePaintDiffusion);
			
			ImGui::Checkbox("Simulation Running", &bSimulationRunning);
			ImGui::Checkbox("Simulation in 2D", &CGConstantsCB.bSimulateParticlesIn2D);
			ImGui::Checkbox("IntegrateVerlet", &bIntegrateVerlet);
			bSimulateOnce = ImGui::Button("Simulate Once");

			ImGui::NewLine();
			if (ImGui::DragFloat("Simulation Bounds Length", &CGConstantsCB.SimulationBoundsLength, 1., 1, GetCollisionGridLengthMax(), "%.1f", ImGuiSliderFlags_ClampOnInput))
			{
				if (bEnableSelfCollisions)
				{
					CGConstantsCB.SimulationBoundsLength = std::min(CGConstantsCB.SimulationBoundsLength, GetCollisionGridLengthMax());
				}
			}
			ImGui::Checkbox("Debug: Visualze Simulation Bounds", &bDebugVisualizeSimulationBounds);
			
			ImGui::NewLine();
			ImGui::SliderInt("Delta Time Scale (1 / 60 * x)", &DeltaTimeScale, 1, 8);

			ImGui::DragFloat("Particle Radius", &CGConstantsCB.PointRadius, 0.001f, 0.01f, 10.f, "%.3f", ImGuiSliderFlags_ClampOnInput);
			ImGui::DragFloat("Rendered Particle Radius", &CGConstantsCB.RasterizedParticleRadiusScale, 0.01f, 0.01f, 10.f, "%.2f", ImGuiSliderFlags_ClampOnInput);
			ImGui::DragFloat3("Default Point Pos ", &CGConstantsCB.DefaultParticlePos.x, 0.1f, 0, 250);
			
			ImGui::NewLine();
			ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Generation");
			
			ImGui::Checkbox("Spawn Objects", ((bool*)&CGConstantsCB.bSpawnObjectsMode));
			bSpawnObjectOnce = ImGui::Button("Spawn Object Once");
			ImGui::TextColored(ImVec4(1, 0.4, 0.2, 1), "Cur Num Particles: %u", CGConstantsCB.NumParticles);
			ImGui::DragInt2("Particle Spawn Area Size", (int*)&CGSpatialControlConstantsCB.ParticleSpawnAreaSize.x, 1, 1, 64, "%d", ImGuiSliderFlags_ClampOnInput);
			if (ImGui::Button("10x10 Spawn Area"))
			{
				CGSpatialControlConstantsCB.ParticleSpawnAreaSize.y = CGSpatialControlConstantsCB.ParticleSpawnAreaSize.x = 10;
			}
			ImGui::ColorEdit3("New Particle Color", &NewlySpawnedParticleColor.x);
			ImGui::ColorEdit3("New Particle Color 2", &NewlySpawnedParticleColor2.x);
			ImGui::ColorEdit3("New Particle Color 3", &NewlySpawnedParticleColor3.x);
			 
			 
			ImGui::NewLine();
			ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "External Forces");
			ImGui::DragFloat("Damping", &CGConstantsCB.Damping, 0.01f, 0.01f, 1.f);
			ImGui::DragFloat("Random Impulse Scale", &RandomImpulseScale, 1.f, 0.f, 250.f);
			ImGui::DragFloat3("ExternalForce", &CGConstantsCB.ExternalForce.x, 0.1f, -250, 250);
			
			if (ImGui::BeginTabBar("PhysTabBar"))
			{
				//Particles
				if (ImGui::BeginTabItem("Particles"))
				{
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Self Collision");
					ImGui::Checkbox("Enable Self Collisions", &bEnableSelfCollisions);
					ImGui::SliderInt("Num Collision Resolve Iterations", (int*)&NumCollisionResolveIterations, 1, 8);
					ImGui::SliderInt("Num Global Collision Resolve Iterations (Expensive)", (int*)&NumGlobalCollisionResolveIterations, 1, 4);

					ImGui::EndTabItem();
				}

				//Collisions
				if (ImGui::BeginTabItem("Collisions"))
				{

					ImGui::TextColored(ImVec4(1, 0.2, 0.2, 1), "External Collision");
					ImGui::Checkbox("External Collisions before Constraints", &bExternalCollisionBeforeConstraints);
					
					ImGui::DragFloat3("Capsule Pos", &SphereCapsulePos1.x, 0.1f, -250, 250);
					ImGui::DragFloat("Capsule Radius", &SphereCapsuleRadius, 0.1f, 0.1f, 100.f);
					ImGui::DragFloat("Capsule Height", &SphereCapsuleHeight, 0.1f, 0.1f, 100.f);

					ImGui::NewLine();
					ImGui::DragFloat("Collision Response Coef", &CGConstantsCB.CollisionResponseCoef, 0.01f, 0.01f, 2.f);

					float collisionGridCellSize = GetCollisionGridCellSize();
					ImGui::Text("Collision Grid Cell Length: %f", collisionGridCellSize);

					ImGui::EndTabItem();
				}


				//Control
				if (ImGui::BeginTabItem("Control"))
				{
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Spatial Control");
					ImGui::DragFloat("Intersection Sphere Scale", &CGSpatialControlConstantsCB.IntersectionSphereScale, 0.1f, 1.f, 10.f);
					ImGui::DragFloat("Scissor Affect Radius", &CGSpatialControlConstantsCB.ScissorAffectRadius, 0.1f, 0.1f, 10.f);
					ImGui::Checkbox("Drag State", (bool*)&CGSpatialControlConstantsCB.bInDragState);
					ImGui::Text("Cur Mouse Offset World: %f, %f, %f", CGSpatialControlConstantsCB.MouseOffsetWorld.x, CGSpatialControlConstantsCB.MouseOffsetWorld.y, CGSpatialControlConstantsCB.MouseOffsetWorld.z);

					ImGui::Text("Cur Mouse Pos NDC: %f, %f", CurMousePosNDC.x, CurMousePosNDC.y);
					ImGui::Text("Prev Mouse Pos NDC: %f, %f", PrevMousePosNDC.x, PrevMousePosNDC.y);
					ImGui::Text("Mouse Pos NDC Offset: %f, %f", CurMousePosNDCOffsetDir.x, CurMousePosNDCOffsetDir.y);

					ImGui::Text("Control World Pos Cur: %f, %f, %f", CursorWorldPosCur.x, CursorWorldPosCur.y, CursorWorldPosCur.z);
					ImGui::Text("Control World Pos Prev: %f, %f, %f", CursorWorldPosPrev.x, CursorWorldPosPrev.y, CursorWorldPosPrev.z);

					ImGui::Text("Clicked Object Pos World / View: %f / %f", IntersectedParticlePos.z, ClickedObjectViewSpaceZ);

					ImGui::Text("First Particle Velocity: %f, %f, %f", FirstParticleVelocity.x, FirstParticleVelocity.y, FirstParticleVelocity.z);

					ImGui::EndTabItem();
				}

				//Constraints & Knots
				if (ImGui::BeginTabItem("Constraints & Knots"))
				{
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Constraints & knots");
					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Generation");
					ImGui::Checkbox("Pin Generated Objects", &bClothGenerationPinGeneratedObjects);
					if (ImGui::Button("Generate Cloth"))
					{
						bGenerateCloth = true;
					}
					ImGui::Checkbox("Cloth Horizontal", &bGenerateClothHorizontal);

					ImGui::DragInt2("Num knots in Cloth", (int*)&ClothSizeCPU.x, 1, 2, 512);
					ImGui::NewLine();
					ImGui::Checkbox("Generate Triangle Mesh", &bGenerateTriangleMesh);
					ImGui::Checkbox("Generate Triangle Mesh Normals", &bGenerateTriangleMeshNormals);
					ImGui::Checkbox("Mesh Discard Broken Connections", &bMeshGeneratorDiscardBrokenConnections);
					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Constraint");
					ImGui::Checkbox("Enable Constraints", &bEnableConstraints);
					ImGui::Checkbox("Enable Soft Pin To Pos", &bEnableSoftPin);
					ImGui::DragFloat("Spring Rest Length", &CGConstantsCB.SpringRestLength, 0.01f, 0.1f, 10.f);
					ImGui::DragFloat("Spring Stgiffness", &CGConstantsCB.ConnectionConstraintStiffness, 0.01f, 0.01f, 10.f);
					ImGui::DragFloat("Spring Stgiffness To Default Pos", &CGConstantsCB.ConnectionConstraintToDefaultPosStiffness, 0.01f, 0.01f, 10.f);
					if (ImGui::Button("PointSize To SpringLength"))
					{
						CGConstantsCB.PointRadius = CGConstantsCB.SpringRestLength * 0.5f;
					}
					ImGui::SliderInt("Num Constraint Resolve Iterations", (int*)&NumConstraintResolveIterations, 1, 100);
					ImGui::Checkbox("Debug: Visualze Constraints", &bDebugVisualizeConstraints);
					

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Water Surface");
					ImGui::Checkbox("Simulate Water Surface", &bClothSimulateWaterSurface);
					ImGui::DragFloat("Pressure Controlled Offset Scale", &CGConstantsCB.PressureControlledOffsetScale, 0.01f, 0.1f, 10.f);
					ImGui::Checkbox("Simulate Waves", &bSimulateWavesSimple);
					if (ImGui::DragFloat("Wave Simulation Speed", &CGConstantsCB.WaveEquationWaveSpeed, 0.01f, 0.1f, 100.f))
					{
						//CGConstantsCB.WaveEquationWaveSpeed = std::min(CGConstantsCB.WaveEquationWaveSpeed, CGConstantsCB.DeltaTime * CGConstantsCB.SpringRestLength);
					}

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Subdivision Surface");
					ImGui::Checkbox("Enable Subdivision Mesh", &bSubdivisionMesh);
					ImGui::DragInt("Subdivision Scale", (int*)&SubdivisionScale, 1, 1, 1024, "%d", ImGuiSliderFlags_ClampOnInput);
					ImGui::Text("Interpolation Method:");
					ImGui::RadioButton("CatmullRom", (int*)&SurfaceSubdivisionSettings.InterpolationMethod, INTERPOLATION_METHOD_CR);
					ImGui::SameLine();
					ImGui::RadioButton("TCB", (int*)&SurfaceSubdivisionSettings.InterpolationMethod, INTERPOLATION_METHOD_TCB);
					ImGui::SameLine();
					ImGui::RadioButton("Lagrange", (int*)&SurfaceSubdivisionSettings.InterpolationMethod, INTERPOLATION_METHOD_LAGRANGE);
					if (SurfaceSubdivisionSettings.InterpolationMethod == INTERPOLATION_METHOD_TCB)
					{
						ImGui::DragFloat("TCBTension", &SurfaceSubdivisionSettings.tension, 0.01f, -5.f, 5.f);
						ImGui::DragFloat("TCBContinuity", &SurfaceSubdivisionSettings.continuity, 0.01f, -5.f, 5.f);
						ImGui::DragFloat("TCBBias", &SurfaceSubdivisionSettings.bias, 0.01f, -5.f, 5.f);
					}

					ImGui::EndTabItem();
				}

				//Velocity Field
				if (ImGui::BeginTabItem("Velocity Field"))
				{
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Velocity Field");

					ImGui::RadioButton("CPU", (int*)&bProcessVectorFieldOnCPU, 1);
					ImGui::SameLine();
					ImGui::RadioButton("GPU", (int*)&bProcessVectorFieldOnCPU, 0);

					ImGui::Checkbox("Generate VectorField On CPU", &bGenerateVectorFieldOnCPU);

					if (ImGui::DragFloat("Velocity Field Bounds Length", &CGConstantsCB.VelocityFieldBoundsLength, 1., 1, GetCollisionGridLengthMax(), "%.1f", ImGuiSliderFlags_ClampOnInput))
					{
						if(!CGConstantsCB.bUseParticleBasedAdvection)
						CGConstantsCB.VelocityFieldBoundsLength = std::min(CGConstantsCB.VelocityFieldBoundsLength, CGConstantsCB.SimulationBoundsLength);
					}
					if (ImGui::DragInt("Velocity Field Grid Size", (int*)&CGConstantsCB.VelocityFieldGridSize, 1, 4, 1024, "%d", ImGuiSliderFlags_ClampOnInput))
					{
						CGConstantsCB.VelocityFieldGridSize = RMath::AlignUp(CGConstantsCB.VelocityFieldGridSize, 2);
					}
					
					ImGui::TextColored(ImVec4(1,1,0.2,1.), "Velocity Field Grid Cell Length: %f", CGConstantsCB.VelocityFieldGridCellLength);

					ImGui::Text("Velocity Field Cur Magnitude Min: %f, Max: %f", VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Velocity Field Generation");
					ImGui::DragFloat("Generated Velocity Magnitude Scale", &CGConstantsCB.GeneratedVelocityMagnitudeScale, 0.01f, 0.01f, 100.f);
#if USE_3COLOR_DENSITY == 1
					ImGui::ColorEdit3("Initial Density Color", &CGConstantsCB.InitialDensityColor.x);
#else
					ImGui::DragFloat("Generated Density Amount", &CGConstantsCB.InitialDensityColor, 0.01f, 0.01f, 1);
#endif//USE_3COLOR_DENSITY

					if (bGenerateVectorFieldOnCPU)
					{
						ImGui::DragInt("Velocity Generation Pointer Radius", (int*)&VelocityGeneratorTriggerRadius, 1, 1, 32);
						ImGui::DragInt("Density Generation Pointer Radius", (int*)&DensityGeneratorTriggerRadius, 1, 1, 32);
					}
					else
					{
						ImGui::DragFloat("Velocity Generation Pointer Radius", &CGConstantsCB.VelocityGeneratorRadius, 0.01f, 0.01f, 10.f);
						ImGui::DragFloat("Density Generation Pointer Radius", &CGConstantsCB.DensityGeneratorRadius, 0.01f, 0.01f, 10.f);
					}

					


					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Velocity");

					ImGui::DragFloat("Velocity Field Damping", &VelocityFieldDamping, 0.01f, 0.01f, 1);
					ImGui::DragFloat("Velocity Field Damping Start Threshold", &VelocityFieldDampingStartThreshold, 0.1f, 0.f, 1000.f);
					ImGui::DragFloat2("Velocity Field Border Force", &VelocityFieldForceAtBorders.x, 0.1f, -250, 250);

					ImGui::Checkbox("Clear VectorField", &bClearVectorField);

					ImGui::Checkbox("Diffuse Velocity", &bDiffuseVelocity);
					ImGui::DragFloat("Velocity Diffuse Diff Scale", &CGConstantsCB.DiffusionFactor, 0.01f, 0.01f, 10000.f);
					
					ImGui::Checkbox("Enable Vector Field Projection", &bEnableVectorFieldProjection);
					ImGui::Checkbox("Vector Field Projection After Advection", &bEnableVectorFieldProjectionAfterAdvection);
					ImGui::Checkbox("Use Pressure Projection Method", &bProjectVectorFieldImplicit);
					ImGui::Checkbox("Use Overrelaxation", &CGConstantsCB.bSolversUseOverrelaxation);
					ImGui::DragFloat("Overrelaxation Constant", &CGConstantsCB.OverrelaxationConstant, 0.01f, 1.f, 5.f);
					ImGui::SliderInt("VecField Proj Num Iterations", (int*)&VectorFieldProjectionNumIterations, 1, 100);
					ImGui::NewLine();
					ImGui::Checkbox("Use Vorticity Confinement", &bEnableVorticityConfinement);
					ImGui::DragFloat("Vorticity Confinement Scale", &CGConstantsCB.VorticityConfinementFactor, 0.01f, 0.f, 1000.f);
					ImGui::SliderInt("Vorticity Num Iterations", (int*)&VorticityNumIterations, 1, 100);
					ImGui::Checkbox("Apply Vorticity Mask To Density", &bApplyVorticityMaskToDensityVisualizer);
					ImGui::NewLine();
					ImGui::Checkbox("Enable Advection", &bEnableAdvection);
					ImGui::Checkbox("Advection MacCormack Method", &bAdvectDensityMacCormack);
					ImGui::Checkbox("Enable Velocity Advection", &bEnableVelocityAdvection);
					ImGui::DragFloat("Velocity Advection Scale", &CGConstantsCB.VelocityAdvectionScaleFactor, 0.01f, 0.01f, 2.f);
					ImGui::Checkbox("Allow Velocity Writes to Solid Cells", &bAllowVelocityWritesToSolidCellsFacesDuringProjection);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Density");
					ImGui::Checkbox("Clear Density", &bClearDensityField);
					if (ImGui::Button("Fill Density with Color"))
					{
						bFillDensityFieldlWithColor = true;
					}
					ImGui::Checkbox("Apply Density To Borders CPU", &bGenerateDensityOnBorder);
					ImGui::Checkbox("Diffuse Density", &bDiffuseDensity);
					ImGui::DragFloat("Density Diffuse Diff Scale", &CGConstantsCB.DiffusionFactor, 0.01f, 0.01f, 10000.f);
					ImGui::SliderInt("Diffusion Num Iterations", (int*)&DiffusionNumIterations, 1, 100);

					ImGui::DragFloat("Density Diffusion Factor", &DensityAdvectionDiffusionFactor, 0.01f, 0.01f, 1);
					ImGui::DragFloat("Density Advection Scale", &DensityAdvectionScaleFactor, 0.01f, 0.f, 1);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Visualizer");
					ImGui::Checkbox("Visualize Velocity Field", &bDebugVisualizeVelocityField);
					ImGui::Checkbox("Visualize Velocity Field : Separate Axises", &bDebugVisualizeVelocityFieldSeparateAxises);
					ImGui::Checkbox("Visualize Velocity Magnitudes", &bDebugVisualizeVelocityFieldMagnitudes);
					ImGui::Checkbox("Visualize Velocity Field Cells", &bDebugVisualizeVelocityFieldCells);
					ImGui::Checkbox("Visualize Velocity Field Cells States", &bDebugVisualizeVelocityFieldCellsStates);
					ImGui::Checkbox("Visualize Divergence", &bDebugVisualizeDivergence);
					ImGui::Checkbox("Visualize Curl", &bDebugVisualizeCurl);
					ImGui::Checkbox("Visualize Curl Gradient", &bDebugVisualizeCurlGradient);
					ImGui::Checkbox("Visualize Pressure", &bDebugVisualizePressure);
					ImGui::Checkbox("Visualize Pressure Gradient", &bDebugVisualizePressureGradient);
					ImGui::Checkbox("Visualize Velocity Field Streamlines", &bDebugVisualizeVectorFieldStreamlines);
					if (bDebugVisualizeVectorFieldStreamlines)
					{
						ImGui::DragInt("Streamline Num Samples", (int*)&StreamLineMaxNumSamples, 1, 1, 2048);
						ImGui::DragInt("Streamline Stride", (int*)&StreamLineStride, 1, 1, 8);
						ImGui::DragFloat("Streamline Step Size", &StreamLineStepSize, 0.01f, 0.01f, 10.f);
					}
					
					ImGui::Checkbox("Visualize Density Field", &bVisualizeDensityField);
					

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Boundaries");

					ImGui::Checkbox("Enable Field Obstacle", &bEnableSphereObstacle);
					if (bEnableSphereObstacle)
					{
						ImGui::Checkbox("Obstacle as Density Source", &bObstacleAsDensitySource);
						ImGui::Checkbox("Obstacle Generates Velocity", &bObstacleGeneratesVelocity);
						ImGui::DragFloat2("Velocity Field Obstacle Pos", &VelocityFieldObstaclePos.x, 0.1f, 0, 250);
						ImGui::DragFloat("Velocity Field Obstacle Radius", &VelocityFieldObstacleRadius, 0.01f, 0.01f, 10.f);
					}

					ImGui::Checkbox("Copy Velocity to Borders", &bCopyVelocitiesToBorders);
					ImGui::Checkbox("Copy Reflected Velocity to Borders", &bCopyReflectedVelocitiesToBorders);
					ImGui::Checkbox("Remove Velocity From Solid Cells", &bClearVelocityFromSolidCell);

					ImGui::Text("Enabled Borders:");
					ImGui::Checkbox("Left Border", &bEnabledBorders[0]);
					ImGui::Checkbox("Up Border", &bEnabledBorders[1]);
					ImGui::Checkbox("Right Border", &bEnabledBorders[2]);
					ImGui::Checkbox("Down Border", &bEnabledBorders[3]);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Info At Point");
					const auto& velAtPoint = LoadFromVelocityField(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo);
					ImGui::Text("Velocity At Point: %f, %f", velAtPoint.x, velAtPoint.y);

					auto cellCenterWorldSpace = GetVelocityCellCenterWorldSpaceFromCoord(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo);
					ImGui::Text("World Pos At Point: %f, %f", cellCenterWorldSpace.x, cellCenterWorldSpace.y);

					ImGui::EndTabItem();
				}
				 

				//Particle based advection
				if (ImGui::BeginTabItem("Particle Based Velocity Advection"))
				{
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Particle Based Velocity Advection");

					ImGui::Checkbox("Use Particle Advection by Velocity Field", &CGConstantsCB.bUseParticleBasedAdvection);
					ImGui::Checkbox("Push Particles from dense regions", &bParticleAdvectionCompensateDrift);
					ImGui::DragFloat("Push Particles Strength Scale", &CGConstantsCB.AdvectedParticleDriftCompensationScale, 0.01f, 0.0f, 2.f);
					ImGui::Text("Average Particle Density: %f", AverageParticleDensity);
					ImGui::DragFloat("PIC-FLIP Ratio", &CGConstantsCB.ParticleAdvectionPICFLIPRatio, 0.01f, 0.0f, 1.f);

					ImGui::Checkbox("TransferParticles Use Compute Shader", &bTransferParticleDataWithComputeShader);
					ImGui::Checkbox("Advected Particles Mark all neighbor cells", &bAllAdvectedParticlesMarkDensityStateToCells);


					ImGui::DragFloat("Diffusion Bumpinnes", &CGSDFConstantsCB.ParticleBumpinnesPow, 0.01f, 0, 100.f);


					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "SDF Render");
					
					ImGui::Checkbox("Render Painting", &bRenderPainting);
					ImGui::DragFloat("SDF Threshold", &CGSDFConstantsCB.Threshold, 0.1f, 0.0f, 10000.f);
					ImGui::DragFloat("SDF Blend Radius", &CGSDFConstantsCB.BlendRadius, 0.01f, -1000.0f, 1000.f);
					ImGui::DragFloat("SDF Surface Depth Offset", &CGSDFConstantsCB.SurfaceDepthOffset, 0.01f, -10.f, 10.f);
					ImGui::DragFloat("SDF Ray Origin Depth Offset", &CGSDFConstantsCB.RayOriginDepthOffset, 0.01f, -100.f, 100.f);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Falling Particle");
					ImGui::Checkbox("Falling Particles", &bDisableParticlesIntersectingBounds);
					ImGui::DragFloat("Falling Particle Start Pos", &CGSDFConstantsCB.FalingParticleStartPosZ, 0.01f, 0.0f, 10.f);
					ImGui::DragFloat("Falling Particle Speed Towards Surface", &CGSDFConstantsCB.FallingParticleSpeedTowardsSurface, 0.01f, 0.0f, 250.f);
					 
					ImGui::Text("Shape Blend Method:");
					ImGui::RadioButton("Exp", (int*)&CGSDFConstantsCB.SDFShapeBlendMode, SDF_BLEND_EXP);
					ImGui::SameLine();
					ImGui::RadioButton("Root", (int*)&CGSDFConstantsCB.SDFShapeBlendMode, SDF_BLEND_ROOT);
					ImGui::SameLine();
					ImGui::RadioButton("Power", (int*)&CGSDFConstantsCB.SDFShapeBlendMode, SDF_BLEND_POWER);
					ImGui::SameLine();
					ImGui::RadioButton("Poly", (int*)&CGSDFConstantsCB.SDFShapeBlendMode, SDF_BLEND_POLY);
					ImGui::SameLine();
					ImGui::RadioButton("Cubic", (int*)&CGSDFConstantsCB.SDFShapeBlendMode, SDF_BLEND_CUBIC);

					ImGui::Text("Normals Blend Method:");
					ImGui::RadioButton("nExp", (int*)&CGSDFConstantsCB.SDFNormalsBlendMode, SDF_BLEND_EXP);
					ImGui::SameLine();
					ImGui::RadioButton("nRoot", (int*)&CGSDFConstantsCB.SDFNormalsBlendMode, SDF_BLEND_ROOT);
					ImGui::SameLine();
					ImGui::RadioButton("nPower", (int*)&CGSDFConstantsCB.SDFNormalsBlendMode, SDF_BLEND_POWER);
					ImGui::SameLine();
					ImGui::RadioButton("nPoly", (int*)&CGSDFConstantsCB.SDFNormalsBlendMode, SDF_BLEND_POLY);
					ImGui::SameLine();
					ImGui::RadioButton("nCubic", (int*)&CGSDFConstantsCB.SDFNormalsBlendMode, SDF_BLEND_CUBIC);


					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Normals");
					ImGui::DragFloat("SDF Normals Smoothness", &CGSDFConstantsCB.NormalsSmoothness, 0.01f, 0.0f, 1.f);
					ImGui::DragFloat("SDF Normals Blend Radius", &CGSDFConstantsCB.NormalsBlendRadius, 0.01f, -1000.0f, 1000.f);
					ImGui::DragFloat("Paint Wetness Contribution to Height", &CGSDFConstantsCB.WetnessContributionToHeight, 0.01f, 0, 1.f);
					ImGui::Checkbox("SDF Height Map Blur", &bEnableHeightMapBlur);
					ImGui::Checkbox("Height Map Blur Normals", &bHeightMapBlurNormals);
					ImGui::SliderInt("Height Map Blur Num Iterations", (int*)&HeightBlurPassNumIterations, 1, 10);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Surface Hit");
					ImGui::DragFloat("SurfaceHit Velocity Scale", &CGSDFConstantsCB.SurfaceHitVelocityScale, 0.01f, 0.0f, 1000.f);
					ImGui::DragFloat("SurfaceHit Velocity Max", &CGSDFConstantsCB.SurfaceHitVelocityMax, 0.1f, 0.0f, 1000.f);
					

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Velocity etc");
					ImGui::DragFloat("Velocity Damp On Dry Surface", &CGSDFConstantsCB.PaintDryAreasVelocityDampForceScale, 0.01f, 0.0f, 1000.f);
					ImGui::DragFloat("Dry Surface Damp Threshold", &CGSDFConstantsCB.DrySurfaceDampThreshold, 0.01f, 0.0f, 1.f);
					ImGui::DragFloat("Wetness Accumulation Amount", &CGSDFConstantsCB.WetnessAccumulationAmount, 0.001f, 0.001f, 100.f);
					ImGui::DragFloat("Wetness Drying Amount", &CGSDFConstantsCB.WetnessDryingOutAmount, 0.001f, 0.001f, 100.f);
					ImGui::DragFloat("Max Wetness Limit", &CGSDFConstantsCB.WetnessLimit, 0.1f, 1.f, 1000.f);

					
					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Color");
					ImGui::Checkbox("Clear Paint Texture", &bClearPaintTexture);
					ImGui::ColorEdit3("SDF Output Color", &CGSDFConstantsCB.OutputColor.x);
					ImGui::DragFloat("SDF Color Blend Radius", &CGSDFConstantsCB.ColorBlendRadius, 0.01f, -1000.0f, 1000.f);

					ImGui::NewLine();
					ImGui::DragFloat("SDF Color Fade Threshold", &CGSDFConstantsCB.ColorFadeThreshold, 0.01f, 0.0f, 1.f);
					ImGui::DragFloat("SDF Color Fade Scale Inv", &CGSDFConstantsCB.ColorFadeScaleInv, 0.01f, 0.0f, 1.f);
					ImGui::DragFloat("SDF Color-Particle Diffusion Scale", &CGSDFConstantsCB.ColorDiffusionScale, 0.01f, 0.0f, 1.f);
					ImGui::DragFloat("Paint Diffusion Factor", &CGSDFConstantsCB.PaintDiffusionFactor, 0.01f, 0.0f, 1000.f);
					ImGui::Checkbox("Blur Paint Texture", &bEnablePaintTextureBlur);
					ImGui::Checkbox("Only Blur Wetness", &bEnablePaintWetnessBlur);
					ImGui::SliderInt("Paint Blur Num Iterations", (int*)&PaintBlurPassNumIterations, 1, 10);

					ImGui::NewLine();
					ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Falling Particles");
					ImGui::DragFloat("SDF Falling Particles Blend Radius", &CGSDFConstantsCB.FallingParticleBlendRadius, 0.01f, -1000.0f, 1000.f);

					ImGui::EndTabItem();

				}


			}	
			ImGui::EndTabBar();
		}
		ImGui::End();
	}


	void AllocateResources()
	{
		///Main
		{
			ParticlesBufferPositionCurGPU.Create(TEXT("ParticlesBufferPositionCur"), MaxNumParticles);
			ParticlesBufferPositionAfterCreationGPU.Create(TEXT("ParticlesBufferPositionAfterCreation"), MaxNumParticles);
			ParticlesBufferPositionPrevGPU.Create(TEXT("ParticlesBufferPositionPrev"), MaxNumParticles);

			CurParticlesVelocityBufferGPU.Create(TEXT("CurParticlesVelocityBuffer"), MaxNumParticles);
			PrevParticlesVelocityBufferGPU.Create(TEXT("PrevParticlesVelocityBuffer"), MaxNumParticles);
			CurParticlesExternalForceBufferGPU.Create(TEXT("CurParticlesExternalForceBuffer"), MaxNumParticles);

			ParticlesStateBufferGPU.Create(TEXT("ParticlesStateBuffer"), MaxNumParticles);

			ParticlesColorsBufferGPU.Create(TEXT("ParticlesColorsBuffer"), MaxNumParticles);

			CGConstantsCBufferGPU.Create(TEXT("CGConstantsCB"));

		}


		///External Collisions
		{
			PlanesBufferGPU.Create(TEXT("Planes"), 6);
			SphereCapsulesBufferGPU.Create(TEXT("Sphere Capsules"), 2);
		}

		///Collisions
		{
			NumCellsObjectIntersectsGPU.Create(TEXT("NumCellsObjectIntersectsGPU"), MaxNumParticles);

			PerObjectControlBitsBufferGPU.Create(TEXT("PerObjectControlBitsBuffer"), MaxNumParticles);

			ObjectListIndexBufferGPU.Create(TEXT("ObjectListIndexBuffer"), MaxNumParticles * 8);
			CellListIndexBufferGPU.Create(TEXT("CellListIndexBuffer"), MaxNumParticles * 8);

			SortedObjectListIndexBufferGPU.Create(TEXT("SortedObjectListIndexBuffer"), MaxNumParticles * 8);
			SortedCellListIndexBufferGPU.Create(TEXT("SortedCellListIndexBuffer"), MaxNumParticles * 8);

			for (auto& buf : CellIdOffsetAndCountBufferGPUArr)
			{
				buf.Create(TEXT("CellIdOffsetAndCountBuffer"), MaxNumParticles * 8);
			}
			//CellIdOffsetAndCountBufferGPU.Create(TEXT("CellIdOffsetAndCountBuffer"), BURNER_MAX_PARTICLES * 8);

			CellCountBufferGPU.Create("CellCountBuffer", 1);

		}


		///Constraints
		{

		}

		///Spatial Control
		{
			CGSpatialControlConstantsCBGPU.Create(TEXT("CGSpatialControlConstantsCB"));

			IntersectedPointIndexGPU.Create(TEXT("IntersectedPointIndex"), 1);
			ClickedObjectViewSpaceZGPU.Create(TEXT("ClickedObjectViewSpaceZ"), 1);
		}

		///Mesh Generator
		{

		}

	}


	void Simulate()
	{

		if (bSimulationRunning || bSimulateOnce)
		{
			RCommandListCompute& CmdList = RCommandListCompute::BeginNew();

			{
				GPU_PROFILE_SCOPE(PhysicsSimulation, CmdList);

				if (GbFirstBoot)
				{
					///
					///	GPU Resources
					///
					AllocateResources();


					if (!bIntegrateVerlet)
					{
						CGConstantsCB.CollisionResponseCoef = 1.f;
					}

					bSpawnObjectOnce = true;


					#if PAINT_RENDER_PIXEL_ADVECT
					CGConstantsCB.bUseParticleBasedAdvection = false;
					bEnableDensity = true;
					bEnableAdvection = true;
					#endif

				}

				///	Update Main CBuffer
				{

					if (GbFirstBoot)
					{
						CGConstantsCB.SimulationBoundsLength = 20.f;
					}

					CGConstantsCB.DeltaTime = 1.f / (60.f * float(DeltaTimeScale));


					SpawnParticles(CmdList);

					CGConstantsCBufferGPU.Update(CGConstantsCB);
				}

				if (bResetSimulation || bResetParticles)
				{
					CGConstantsCB.NumParticles = 1;
					CmdList.ClearUAV(ParticlesStateBufferGPU);
				}

				///Generate Cloth
				if (bGenerateCloth)
				{
					GenerateKnotsConstraintsPhysObject(CmdList);
					bGenerateCloth = false;
				}

				if (!CGConstantsCB.bUseParticleBasedAdvection || GbFirstBoot)
				{

					CGConstantsCB.VelocityFieldGridCellLength = (CGConstantsCB.VelocityFieldBoundsLength) / float(CGConstantsCB.VelocityFieldGridSize);

#if PAINT_RENDER_PIXEL_ADVECT
					CGConstantsCB.SimulationBoundsLength = CGConstantsCB.VelocityFieldBoundsLength - CGConstantsCB.VelocityFieldGridCellLength * 2;
					CGConstantsCB.SimulationBoundsStartOffset = CGConstantsCB.VelocityFieldBoundsStartOffset + CGConstantsCB.VelocityFieldGridCellLength;
					//Fill GBuffer
					if(GbFirstBoot || bClearPaintTexture)
					RenderGBuffer(CmdList);
#endif

					//Generate Velocity Field
					ProcessVelocityField2D(CmdList);
				}

				if (CGConstantsCB.bUseParticleBasedAdvection)//Particles use velocity field bounds
				{
					CGConstantsCB.SimulationBoundsLength = CGConstantsCB.VelocityFieldBoundsLength - CGConstantsCB.VelocityFieldGridCellLength * 2;
					CGConstantsCB.SimulationBoundsStartOffset = CGConstantsCB.VelocityFieldBoundsStartOffset + CGConstantsCB.VelocityFieldGridCellLength;
				}

				if (bClothSimulateWaterSurface)
				{
					ApplyPressureFieldOffsetToSurfaceParticles(CmdList);
				}

				if (bSimulateWavesSimple)
				{
					SimulateWave(CmdList);
				}

				if (!CGConstantsCB.bUseParticleBasedAdvection)
				{
					//AdvectParticlesWithVelocityField(CmdList);
				}

				//if (1 && Keyboard::IsKeyPressed('V'))//Mouse controlled velocity
				if(1 && Mouse::Get().IsRButtonPressed())
				{
					ApplyPointerVelocityToParticles(CmdList);
				}
				
				
				IntegrateParticles(CmdList);

				if (bExternalCollisionBeforeConstraints)
				{
					ExternalCollisionPass(CmdList);
				}
				
				if (bEnableConstraints)
				{
					ResolveConstraints(CmdList);
				}

				if (bEnableSoftPin)
				{
					ResolveConstraintsToDefaultPos(CmdList);
				}

				for (int i = 0; i < NumGlobalCollisionResolveIterations; i++)
				{
					if (bEnableSelfCollisions)
					{
						PointPointCollisionPass(CmdList);
					}

					//if (!bExternalCollisionBeforeConstraints)
					{
						ExternalCollisionPass(CmdList);
					}
				}

				if (CGConstantsCB.bUseParticleBasedAdvection)
				{
#if PARTICLE_ADVECTION_DEBUG
					RGeometry::GInstance.SetColor(EColor::Yellow);
					RGeometry::GInstance.AddCircle(CGConstantsCB.DefaultParticlePos);
#endif//PARTICLE_ADVECTION_DEBUG

					if (bProcessVectorFieldOnCPU == 1)
					{
						ExtractParticlesDataFromGPU(CmdList);
						
						TransferDataFromParticlesToVelocityFieldCPU(CmdList);

						if (bParticleAdvectionCompensateDrift)
						{
							ComputeAdvectedParticleDensityCPU();
						}
					}
					else
					{
						TransferDataFromParticlesToVelocityField(CmdList);
					}

					//Projection
					ProcessVelocityField2D(CmdList);

					if (bProcessVectorFieldOnCPU == 1)
					{
						ComputeParticlesColorsCPU();

						TransferDataFromVelocityFieldToParticlesCPU();

						UploadParticlesDataToGPU(CmdList);
					}
					else
					{
						if (bAssignParticlesColorsBasedOnDensity)
						{
							ComputeParticlesColors(CmdList);
						}

						TransferDataFromVelocityFieldToParticles(CmdList);
					}
				}







				

				SpatialControl_Main(CmdList);

				MeshGeneratorMain(CmdList);

				if (bDebugVisualizeConstraints)
				{
					DebugVisualizeConstraints(CmdList);
				}
				if (bDebugVisualizeSimulationBounds)
				{
					DebugVisualizeSimulationBounds();
				}
				
				if (bRenderPainting)
				{
					RenderPainting(CmdList);
				}

				if (bDisableParticlesIntersectingBounds)
				{
					DisableParticlesIntersectingBounds(CmdList);
				}

				ReadFirstParticleVelocity(CmdList);

			}
			CmdList.ExecuteCmdListAndReleaseContext(false);
		}

	}

	float3 FirstParticleVelocity;
	void ReadFirstParticleVelocity(RCommandListCompute& CmdList)
	{
		RReadbackBuffer velReadbackBuffer;
		velReadbackBuffer.Create(TEXT("velReadbackBuffer"), 1, sizeof(float3));

		CmdList.CopyBufferRegion(velReadbackBuffer, 0, CurParticlesVelocityBufferGPU, 0, sizeof(float3));
		CmdList.ExecuteCmdList(true);

		velReadbackBuffer.ExtractContents(&FirstParticleVelocity, sizeof(float3));
	}




	float3 NewlySpawnedParticleColor = float3(1, 0.1, 0.2);
	float3 NewlySpawnedParticleColor2 = float3(0.2, 1, 1);
	float3 NewlySpawnedParticleColor3 = float3(0.5, 1, 0.2);


	float3 PaintEmitterPosition{ 2,2,0 };

	float2 PaintParticleSpawnPos{ 0,7.4 };
	float2 PaintParticleSpawnVel{ 0.7,0.09 };
	float PaintParticleSpawnVelScale = 276;
	float ParticleSpawnStepSizeScale = 1.f;
	float3 PointCloudMeshOrientation{ 0,0,0 };

	float2 BucketPaintSinCosFreqScale{ 0.1, 0.1 };
	float2 BucketPaintSinCosScale{ 1.f, 1.f };

	float VelocityControllerScale = 5.f;
	float VelocityControllerMax = 100.f;

	bool bEnableSingleClickSpawn = true;


	void AssignNewColors()
	{
		NewlySpawnedParticleColor = GenerateRandomColor();
		NewlySpawnedParticleColor2 = GenerateRandomColor();
		NewlySpawnedParticleColor3 = GenerateRandomColor();
	}


	void SpawnParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(SpawnParticles, CmdList);

		static bool bIsSpawnKeyPressedThisFrame{ false };
		bool bHasSpawnKeyBeenPressedPrevFrame = bIsSpawnKeyPressedThisFrame;
		bIsSpawnKeyPressedThisFrame = Keyboard::IsKeyPressed('T');

		bool bManualSpawnEnabled;

		if (bEnableSingleClickSpawn)
		{
			bManualSpawnEnabled = bHasSpawnKeyBeenPressedPrevFrame ? false : bIsSpawnKeyPressedThisFrame;
		}
		else
		{
			bManualSpawnEnabled = bIsSpawnKeyPressedThisFrame;
		}


		const bool bAutoSpawnEnabled = (bSpawnObjectOnce || CGConstantsCB.bSpawnObjectsMode) && CGConstantsCB.NumParticles != (MaxNumParticles - 1);


		//Generate New Colors On Demand
		if (Keyboard::IsKeyPressed('C'))
		{
			AssignNewColors();
		}

		if (1 || bManualSpawnEnabled)
		{
			//compute cur spawner pos
			auto& viewMat = RScene::Get().GetCamera().Matrices.View;
			float viewSpaceZ = float3{ float3SSE{ DirectX::XMVector3Transform(float3SSE{ CGConstantsCB.DefaultParticlePos }, viewMat) } }.z;

			///Calculate World Space Mouse Pos from ClickedObjectViewSpaceZ
			float3 spawnerWorldPosCur = RScene::Get().GetCamera().GetMouseCursorPositionWorld(viewSpaceZ);

			CGSpatialControlConstantsCB.ParticleSpawnerPosition = spawnerWorldPosCur.xy();

		}
		else if(bAutoSpawnEnabled)
		{
			CGSpatialControlConstantsCB.ParticleSpawnerPosition = CGConstantsCB.DefaultParticlePos.xy();

			auto randomDir = VectorNormalize(float3{ RandomNumberGenerator::Get().NextFloat(-1.f, 1.f), RandomNumberGenerator::Get().NextFloat(-1.f, 1.f), RandomNumberGenerator::Get().NextFloat(-1.f, 1.f) });
			CGConstantsCB.RandomImpulse = randomDir * RandomImpulseScale;
		}

		const float2 curSpawnerPos = CGSpatialControlConstantsCB.ParticleSpawnerPosition;
		const uint2 spawnAreaSize = CGSpatialControlConstantsCB.ParticleSpawnAreaSize;
		const float2 spawnAreaLength = float2{ spawnAreaSize.x * CGConstantsCB.PointRadius * 2, spawnAreaSize.y * CGConstantsCB.PointRadius * 2 };
		const float2 spawnAreaHalfLength = spawnAreaLength / 2.f;
		float2 spawnOriginPos = (curSpawnerPos - spawnAreaHalfLength) + float2{ CGConstantsCB.PointRadius, CGConstantsCB.PointRadius };

		if (0 && (bManualSpawnEnabled || bAutoSpawnEnabled))//visualize spawn area
		{
			RGeometry::GInstance.SetColor(EColor::Yellow);
			RGeometry::GInstance.AddCircleCameraFacingLineList(float3{ CGSpatialControlConstantsCB.ParticleSpawnerPosition, /*CGConstantsCB.DefaultParticlePos.z*/ 0 }, VectorGetLength(spawnAreaHalfLength));
			//RGeometry::GInstance.AddPlaneLineListXY(float3{ CGSpatialControlConstantsCB.ParticleSpawnerPosition, /*CGConstantsCB.DefaultParticlePos.z*/ 0}, spawnAreaLength);
		}

		//const float SpawnDepthPosition = CGConstantsCB.DefaultParticlePos.z;
		float SpawnDepthPosition = CGSDFConstantsCB.FalingParticleStartPosZ;
		float SpawnEndDepthPosition = std::min(CGConstantsCB.DefaultParticlePos.z - 0.5f, SpawnDepthPosition + CGConstantsCB.DefaultParticlePos.z * 0.25f);

		float spawnStepSize = CGConstantsCB.PointRadius * 2.f;
		spawnStepSize *= ParticleSpawnStepSizeScale;

		if (ImGui::Begin("PaintSpawn"))
		{
			ImGui::Checkbox("One Time Spawn", &bEnableSingleClickSpawn);
			ImGui::DragFloat2("Spawn Pos", &PaintParticleSpawnPos.x, 0.1, 0, 30);
			ImGui::DragFloat2("Spawn Vel", &PaintParticleSpawnVel.x, 0.01, 0, 1);
			ImGui::DragFloat("Spawn Vel Scale", &PaintParticleSpawnVelScale, 0.1f, 0.0f, 1000.f);
			ImGui::DragFloat("Vel Controller Scale", &VelocityControllerScale, 0.1f, 0.0f, 1000.f);
			ImGui::DragFloat("Vel Controller Max", &VelocityControllerMax, 0.1f, 0.0f, 1000.f);
			ImGui::DragFloat("Falling Particle Start Pos", &CGSDFConstantsCB.FalingParticleStartPosZ, 0.01f, 0.0f, 10.f);
			ImGui::DragFloat("Falling Particle Speed Towards Surface", &CGSDFConstantsCB.FallingParticleSpeedTowardsSurface, 0.1f, 0.0f, 100.f);
			ImGui::DragFloat("Spawn Step Size Scale", &ParticleSpawnStepSizeScale, 0.01, 0.1, 10);
			ImGui::SliderFloat3("Point Cloud Orientation", &PointCloudMeshOrientation.x, -PiX2, PiX2);
			ImGui::DragFloat2("BucketPaint SinCos Freq Scale", &BucketPaintSinCosFreqScale.x, 0.01, 0.01, 10);
			ImGui::DragFloat2("BucketPaint SinCos Scale", &BucketPaintSinCosScale.x, 0.01, 0.01, 10);

		}
		ImGui::End();

#if 0
		RGeometry::GInstance.SetColor(EColor::Yellow);
		RGeometry::GInstance.AddArrowPolyboardDirection(float3{ spawnOriginPos, 0 }, float3{ (CurMousePosNDCOffsetDir) * 10, 0.f });
#endif



		//bManualSpawnEnabled = true;//DEBUG
		
		if(bAutoSpawnEnabled || bManualSpawnEnabled)
		{


			float2 curSpawnPos = spawnOriginPos;






			static RDynamicVector<float3> pointCloudMesh;
			pointCloudMesh.clear();


			//Generate PointCloud Mesh
			/*{
				curSpawnPos = PaintParticleSpawnPos;
				float curDepth = CGSDFConstantsCB.FalingParticleStartPosZ;

				for (int x = 0; x < spawnAreaSize.x; x++)
				{
					for (int y = 0; y < spawnAreaSize.y; y++)
					{

						curDepth = CGSDFConstantsCB.FalingParticleStartPosZ + sin(curSpawnPos.x) + cos(curSpawnPos.y);

						float3 newParticlePos = float3{ curSpawnPos, curDepth };


						pointCloudMesh.push_back(newParticlePos);

						curSpawnPos.y += spawnStepSize;
					}

					curSpawnPos.y = PaintParticleSpawnPos.y;
					curSpawnPos.x += spawnStepSize;
				}

				RBoundingBox box;
				box.CreateFromVertices(pointCloudMesh.size(), (const char*)pointCloudMesh.data(), sizeof(float3));
				float3 boxPosition = box.GetPosition();

				for (auto& point : pointCloudMesh)
				{
					point = point - boxPosition;
				}

				auto mat = MatrixOrientationEuler(PointCloudMeshOrientation.x, PointCloudMeshOrientation.y, PointCloudMeshOrientation.z);
				for (auto& point : pointCloudMesh)
				{
					point = VectorMul(point, mat);
				}

				for (auto& point : pointCloudMesh)
				{
					point = point + boxPosition;
				}
			}*/

			




			//spawnOriginPos = PaintParticleSpawnPos;//DEBUG
			curSpawnPos = spawnOriginPos;

			//float3 curSpawnVel = float3{ PaintParticleSpawnVel, 0.f } *PaintParticleSpawnVelScale;

			float2 dir = (CurMousePosNDCOffsetDir);
			//float2 dir = (CursorWorldPosCur.xy() - CursorWorldPosPrev.xy());
			//float2 dir = CGSpatialControlConstantsCB.MouseOffsetWorld.xy();
			dir  = dir * VelocityControllerScale * 1000.f;
			dir.x = std::clamp(dir.x, -VelocityControllerMax, VelocityControllerMax);
			dir.y = std::clamp(dir.y, -VelocityControllerMax, VelocityControllerMax);
			float3 curSpawnVel = float3{ dir , 0.f };


			static RDynamicVector<float3> newParticlesPositions;
			newParticlesPositions.clear();

			static RDynamicVector<float3> newParticlesVelocities;
			newParticlesVelocities.clear();

			static RDynamicVector<float3> newParticlesColors;
			newParticlesColors.clear();

			for (int x = 0; x < spawnAreaSize.x; x++)
			{
				for (int y = 0; y < spawnAreaSize.y; y++)
				{
					float3 newParticlePos = float3{ curSpawnPos, SpawnDepthPosition };

					float a = x * spawnAreaSize.y;
					float t = (a + y) / float(spawnAreaSize.x * spawnAreaSize.y);

					newParticlePos.z = SpawnDepthPosition + ( (sin(float(x) * BucketPaintSinCosFreqScale.x) + 1.f) * BucketPaintSinCosScale.x + (cos(float(y) * BucketPaintSinCosFreqScale.y) + 1.f) * BucketPaintSinCosScale.y);

					newParticlesPositions.push_back(newParticlePos);

					newParticlesVelocities.push_back(curSpawnVel);

					//interpolate colors
					//float3 curColor = TLinearInterpolate(t, NewlySpawnedParticleColor, NewlySpawnedParticleColor2);
					//float3 curColor = TQuadraticInterpolateBezier(t, NewlySpawnedParticleColor, NewlySpawnedParticleColor2, NewlySpawnedParticleColor3);
					float3 curColor = TQuadraticInterpolate(t, NewlySpawnedParticleColor, NewlySpawnedParticleColor2, NewlySpawnedParticleColor3);
					newParticlesColors.push_back(curColor);
					//newParticlesColors.push_back(NewlySpawnedParticleColor);

					curSpawnPos.y += spawnStepSize;

				}

				curSpawnPos.y = spawnOriginPos.y;
				curSpawnPos.x += spawnStepSize;
			}

			const uint numNewParticles = newParticlesPositions.size();

			//CGConstantsCB.NumParticles = 0;//DEBUG

			//Upload new particles
			CmdList.UploadToBuffer(ParticlesBufferPositionCurGPU, CGConstantsCB.NumParticles * sizeof(float3), newParticlesPositions.data(), numNewParticles * sizeof(float3));
			CmdList.UploadToBuffer(ParticlesBufferPositionPrevGPU, CGConstantsCB.NumParticles * sizeof(float3), newParticlesPositions.data(), numNewParticles * sizeof(float3));
			CmdList.UploadToBuffer(CurParticlesVelocityBufferGPU, CGConstantsCB.NumParticles * sizeof(float3), newParticlesVelocities.data(), numNewParticles * sizeof(float3));
			/*CmdList.UploadToBuffer(ParticlesBufferPositionCurGPU, 0, pointCloudMesh.data(), numNewParticles * sizeof(float3));
			CmdList.UploadToBuffer(ParticlesBufferPositionPrevGPU,0, pointCloudMesh.data(), numNewParticles * sizeof(float3));*/

			CmdList.UploadToBuffer(ParticlesColorsBufferGPU, CGConstantsCB.NumParticles * sizeof(float3), newParticlesColors.data(), numNewParticles * sizeof(float3));
			//CmdList.UploadToBuffer(ParticlesColorsBufferGPU, 0, newParticlesColors.data(), numNewParticles * sizeof(float3));

			{
				static RDynamicVector<uint> newParticlesStates;
				newParticlesStates.clear();
				newParticlesStates.resize(numNewParticles);
				std::fill(newParticlesStates.begin(), newParticlesStates.end(), BITMASK_PARTICLE_STATE_FALLING | BITMASK_PARTICLE_STATE_NEW);
				CmdList.UploadToBuffer(ParticlesStateBufferGPU, CGConstantsCB.NumParticles * sizeof(uint), newParticlesStates.data(), numNewParticles * sizeof(uint));
				//CmdList.UploadToBuffer(ParticlesStateBufferGPU, 0 * sizeof(uint), newParticlesStates.data(), numNewParticles * sizeof(uint));
			}

			CGConstantsCB.NumParticles += numNewParticles;
			//CGConstantsCB.NumParticles = numNewParticles;//DEBUG
		}
	}




	template<bool bIntegrateVerlet = true>
	class IntegrateParticlesShader : public RShader
	{
	public:
		IntegrateParticlesShader()
		{
			if constexpr (bIntegrateVerlet)
			{
				PushBackDefine({ L"INTEGRATE_VERLET", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "SimulateParticles", L"cs_6_0");
		}
	};

	bool bIntegrateVerlet = false;

	void IntegrateParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ParticleIntegration, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionPrevUAV"), ParticlesBufferPositionPrevGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferUAV"), CurParticlesVelocityBufferGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PrevParticlesVelocityBufferUAV"), PrevParticlesVelocityBufferGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesExternalForceBufferUAV"), CurParticlesExternalForceBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PlanesBufferSRV"), PlanesBufferGPU, false);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

		///
		///	Shader
		///
		if (bIntegrateVerlet)
		{
			static IntegrateParticlesShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		else
		{
			static IntegrateParticlesShader<false> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		
		constexpr uint ThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize), 1, 1);

	}




	///*********************
	///		COLLISION
	///*********************
	
	///Colliders

	TStructuredBuffer<float4> PlanesBufferGPU;

	TStructuredBuffer<float4> SphereCapsulesBufferGPU;


	bool bExternalCollisionBeforeConstraints = false;

	//float3 SphereCapsulePos1{ 12,12,12 };
	float3 SphereCapsulePos1{ 0,10,0 };
	float3 SphereCapsulePos2;
	float SphereCapsuleRadius = 2.f;
	float SphereCapsuleHeight = 2.f;


	void GenerateBounds(RCommandListCompute& CmdList)
	{
		if (bEnableSelfCollisions)
		{
			checkf(CGConstantsCB.SimulationBoundsLength < GetCollisionGridLengthMax(), "Bounds bigger than acceleration grid size %f", GetCollisionGridLengthMax());
		}

		float initialBoundScale = CGConstantsCB.SimulationBoundsLength / 2.f;
		
		RArray<BoundingPlane, 6> BoundingPlanes;
		///6 planes
		///down
		BoundingPlanes[0] = BoundingPlane(float3{ 0, -initialBoundScale, 0 }, float3{ 0,1,0 });
		///up
		BoundingPlanes[1] = BoundingPlane(float3{ 0, initialBoundScale, 0 }, float3{ 0,-1,0 });
		///left
		BoundingPlanes[2] = BoundingPlane(float3{ -initialBoundScale, 0, 0 }, float3{ 1,0,0 });
		///right
		BoundingPlanes[3] = BoundingPlane(float3{ initialBoundScale, 0, 0 }, float3{ -1,0,0 });
		///forward
		BoundingPlanes[4] = BoundingPlane(float3{ 0, 0, initialBoundScale }, float3{ 0,0,-1 });
		///back
		BoundingPlanes[5] = BoundingPlane(float3{ 0, 0, -initialBoundScale }, float3{ 0,0,1 });

		std::array<float4, 6> PlanesArr;
		for (size_t i = 0; i < 6; i++)
		{
			//Translate the plane to the positive quadrant
			float3 offset = float3::MakeFloat3(initialBoundScale);
			BoundingPlanes[i].Translate(offset);

			//Additional Translate to avoid touching other quadrants
			BoundingPlanes[i].Translate(float3{ CGConstantsCB.SimulationBoundsStartOffset, CGConstantsCB.SimulationBoundsStartOffset, CGConstantsCB.SimulationBoundsStartOffset });

			PlanesArr[i] = float4{ BoundingPlanes[i].n, BoundingPlanes[i].d };
		}

		///PLanes
		CmdList.UploadToBuffer(PlanesBufferGPU, 0, PlanesArr.data(), 6 * sizeof(float4));


		///Capsules

		SphereCapsulePos2 = SphereCapsulePos1;
		SphereCapsulePos2.y += SphereCapsuleHeight;

		std::array<float4, 2> capsuleData;
		capsuleData[0] = float4{ SphereCapsulePos1, SphereCapsuleRadius };
		capsuleData[1] = float4{ SphereCapsulePos2, SphereCapsuleRadius };

		CmdList.UploadToBuffer(SphereCapsulesBufferGPU, 0, capsuleData.data(), 2 * sizeof(float4));


	}


	class ExternalCollisionShader : public RShader
	{
	public:
		ExternalCollisionShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ResolveExternalCollisions", L"cs_6_0");
		}
	};


	void ExternalCollisionPass(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ExternalCollisionPass, CmdList);


		GenerateBounds(CmdList);



		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PlanesBufferSRV"), PlanesBufferGPU, false);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SphereCapsulesBufferSRV"), SphereCapsulesBufferGPU, false);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferUAV"), CurParticlesVelocityBufferGPU, true);

		///
		///	Shader
		///
		static ExternalCollisionShader ComputeShader;

		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize), 1, 1);


		///Visualize Capsule
		RGeometry::GInstance.SetColor(EColor::Blue);
		RGeometry::GInstance.AddCapsule(SphereCapsulePos1, SphereCapsulePos2, SphereCapsuleRadius, SphereCapsuleRadius);

	}

	///************************
	///		SELF COLLISION
	///************************

	///Collision Pass Resources
	TStructuredBuffer<uint> NumCellsObjectIntersectsGPU;

	TStructuredBuffer<uint> PerObjectControlBitsBufferGPU;///TODO:Replace with uint16 or merge with other state bitmask

	TStructuredBuffer<uint> ObjectListIndexBufferGPU;
	TStructuredBuffer<uint> CellListIndexBufferGPU;

	TStructuredBuffer<uint> SortedObjectListIndexBufferGPU;
	TStructuredBuffer<uint> SortedCellListIndexBufferGPU;

	TStructuredBuffer<uint> CellCountBufferGPU;

	///CBuffer Constants
	GCollisionConstants GCollisionConstantsCB;
	GMeshGenerationConstants GMeshGenerationConstantsCB;

	//TStructuredBuffer<uint2> CellIdOffsetAndCountBufferGPU;
	//static std::vector<uint2> CellIdOffsetAndCountBuffer;
	std::array<RDynamicVector<uint2>, 8> CellIdOffsetAndCountBufferArr; ///Separate Cells based on local index and execute separately
	std::array<TStructuredBuffer<uint2>, 8>  CellIdOffsetAndCountBufferGPUArr;


	uint CellCount;

	///Collision Pass Shaders
	class GenerateIntersectionPairListShader : public RShader
	{
	public:
		GenerateIntersectionPairListShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "GenerateIntersectionPairList", L"cs_6_0");
		}
	};

	class ResolveCollisionsShader : public RShader
	{
	public:
		ResolveCollisionsShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ResolveCollisions", L"cs_6_0");
		}
	};


	void PointPointCollisionPass(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(PointPoint_Collision, CmdList);

		ParticlePhysicsRootSig& RootSig = ParticlePhysicsRootSig::Get();

		///Update Resources
		static const uint InitialValue = 0;
		CmdList.FillBuffer(CellCountBufferGPU, 0, InitialValue, 4);

		///Create {ObjectID, CellID} pairs List
		{
			GPU_PROFILE_SCOPE(Fill_Spatial_Grid, CmdList);

			CmdList.SetRootSignature(RootSig);

			CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ObjectListIndexBufferUAV"), ObjectListIndexBufferGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CellListIndexBufferUAV"), CellListIndexBufferGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AtomicOffsetBufferUAV"), CellCountBufferGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PerObjectControlBitsBufferUAV"), PerObjectControlBitsBufferGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

			///
			///	Shader
			///
			static GenerateIntersectionPairListShader ComputeShader;

			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
			constexpr uint ThreadGroupSize = 64;
			CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize), 1, 1);
		}

		///Sort by CellID
		{
			RReadbackBuffer TempBuffer;
			TempBuffer.Create(TEXT("CellCountReadBuffer"), 1, 4);

			CmdList.CopyResource(TempBuffer, CellCountBufferGPU);

			CmdList.ExecuteCmdList(true);

			TempBuffer.ExtractContents(&CellCount);

			if (CellCount == 0) //no points 
			{
				return;
			}


			///TODO: Use Indirect Buffer instead!


			{
				GPU_PROFILE_SCOPE_TEXT(CmdList, L"Sort Cells List (NumCells: %i)", CellCount);
				RadixSortGPU(CmdList, CellListIndexBufferGPU, SortedCellListIndexBufferGPU, ObjectListIndexBufferGPU, SortedObjectListIndexBufferGPU, CellCount);
			}
		}

		


		///TODO:Move the code Below to GPU! https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-32-broad-phase-collision-detection-cuda


		///Generate Offsets (CPU)
		{
			GPU_PROFILE_SCOPE(GenerateCollisionCellsPerLocalId, CmdList);
			///Get Sorted CellID buffer on CPU

			RReadbackBuffer CellIdBufferReadBack;
			CellIdBufferReadBack.Create(TEXT("CellIdReadBuffer"), CellCount, 4);

			CmdList.CopyBufferRegion(CellIdBufferReadBack, 0, CellListIndexBufferGPU, 0, CellCount * 4u);
			CmdList.ExecuteCmdList(true);

			static RDynamicVector<uint> CellIdBufferCPU;
			CellIdBufferCPU.resize(CellCount+1);
			CellIdBufferReadBack.ExtractContents(CellIdBufferCPU.data(), CellCount * 4u);
			CellIdBufferCPU.back() = -1;///Paste -1 here indicating the end of the cells

			///Count repeated CellIds and generate offsets and count

			for (auto& cpuBuf : CellIdOffsetAndCountBufferArr)
			{
				cpuBuf.clear();
			}
			//CellIdOffsetAndCountBuffer.resize(0);/// size == num unique cells

			uint curCount = 1;///Num Objects in this cell
			uint offset = 0; ///Cell Start index
			uint lastCellId = CellIdBufferCPU[0];
			for (int i = 1; i < CellCount+1; i++)///we need to reach "-1"
			{
				uint curCellId = CellIdBufferCPU[i];
				if (curCellId != lastCellId)///new cell or we reached the end
				{
					///New cell is started, write down offset and count for prev cell
					uint2 offsetAndCount;
					offsetAndCount.x = offset;
					offsetAndCount.y = curCount;
					uint cellLocalId = GetLocalCellIndex(GetCellCoord(lastCellId));
					if (curCount > 1)
					{
						CellIdOffsetAndCountBufferArr[cellLocalId].push_back(offsetAndCount);
						//CellIdOffsetAndCountBuffer.push_back(offsetAndCount);
					}
					///and update offset
					offset = i;
					curCount = 1;
				}
				else///cur cell repeats
				{
					curCount++;
				}
				lastCellId = curCellId;
			}

			///Update GPU buffers
			for (int i = 0; i < 8; i++)
			{
				auto& offsetCountBufferPerLocalId = CellIdOffsetAndCountBufferArr[i];
				if (!offsetCountBufferPerLocalId.empty())
				{
					uint numCellsWithThisLocalId = offsetCountBufferPerLocalId.size();
					//CellIdOffsetAndCountBufferGPUArr[i].UploadElements(0, numUniqueCells, bufCPU.data());
					CmdList.UploadToBuffer(CellIdOffsetAndCountBufferGPUArr[i], 0, offsetCountBufferPerLocalId.data(), numCellsWithThisLocalId * sizeof(uint2));
				}
			}
		}

		 
		///Generate Collision Resolve Work based on num unique cells
		for(int i = 0; i < NumCollisionResolveIterations; i++)
		{
			GPU_PROFILE_SCOPE(ResolveCollisions, CmdList);

			///Bind Const Resources
			{
				CmdList.SetRootSignature(RootSig);

				CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ObjectListIndexBufferSRV"), ObjectListIndexBufferGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PerObjectControlBitsBufferSRV"), PerObjectControlBitsBufferGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionPrevUAV"), ParticlesBufferPositionPrevGPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferUAV"), ParticlesColorsBufferGPU, true);

				static ResolveCollisionsShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
			}

			
			///8 passes per each local id
			for (int i = 0; i < 8; i++)
			{
				uint numUniqueCells = CellIdOffsetAndCountBufferArr[i].size();
				if (numUniqueCells > 0)
				{
					GPU_PROFILE_SCOPE_TEXT(CmdList, L"LocalCellGroup %i (Num Cells: %i)", i, numUniqueCells);

					SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

					///Update Constants
					GCollisionConstantsCB.CurCellLocalId = i;
					GCollisionConstantsCB.NumCells = numUniqueCells;

					CmdList.SetConstantArray(RootSig.GetRootParamIndex("GCollisionConstantsCB"), 2, &GCollisionConstantsCB);

					auto& CellIdOffsetAndCountBufferGPUCur = CellIdOffsetAndCountBufferGPUArr[i];
					SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CellIdOffsetAndCountBufferSRV"), CellIdOffsetAndCountBufferGPUCur, true);

					constexpr uint ThreadGroupSize = 64;
					CmdList.Dispatch(RMath::DivideAndRoundUp(numUniqueCells, ThreadGroupSize), 1, 1);

				}
			}

		}

	}

	///************************
	///		CONSTRAINTS
	///************************
	
	RDynamicVector<float3> KnotsBufferCPU;
	std::array<RDynamicVector<CConstraint>, CONSTRAINT_NUM_COLOR_TYPES> EdgesBufferArrCPU;

	RDynamicVector<uint> KnotsStatesBufferCPU;

	std::array<TStructuredBuffer<CConstraint>, CONSTRAINT_NUM_COLOR_TYPES> ConstraintsBufferArrGPU;
	std::array<TStructuredBuffer<uint>, CONSTRAINT_NUM_COLOR_TYPES> ConstraintsStatesBufferArrGPU;

	uint2 ClothSizeCPU{ 64,64 };

	bool bEnableSoftPin = false;
	bool bClothGenerationPinGeneratedObjects = false;
	bool bGenerateClothHorizontal = true;

	uint GetFlatIndexToClothArray(uint2 coord)
	{
		return coord.y * CGConstantsCB.ClothSize.x + coord.x;
	};

	///TODO: Move generation to GPU
	void GenerateKnotsConstraintsPhysObject(RCommandList& CmdList)
	{
		GPU_PROFILE_SCOPE(GenerateKnotsConstraintsPhysObject, CmdList);

		//Copy ClothSizeCPU to CBuffer version
		//Align ClothSize to multiple of 2 and add 1 so that we have the same amount of constraints per each color group
		static constexpr uint alignment = 2;
		CGConstantsCB.ClothSize.x = RMath::AlignUp(ClothSizeCPU.x, alignment) + 1;
		CGConstantsCB.ClothSize.y = RMath::AlignUp(ClothSizeCPU.y, alignment) + 1;

		if (bClothSimulateWaterSurface)
		{
			CGConstantsCB.DefaultParticlePos.y = CGConstantsCB.DefaultParticlePos.z = CGConstantsCB.DefaultParticlePos.x;
			CGConstantsCB.ClothSize.y = CGConstantsCB.ClothSize.x;
			bClothGenerationPinGeneratedObjects = true;
			bEnableSoftPin = true;
			CGConstantsCB.DefaultParticlePos.y = std::max(CGConstantsCB.DefaultParticlePos.y, WATER_SURFACE_SIMULATION_HEIGHT);
			CGConstantsCB.ConnectionConstraintToDefaultPosStiffness = 0.1f;
		}

		//Generate initial state on CPU first

		KnotsBufferCPU.resize(CGConstantsCB.ClothSize.x * CGConstantsCB.ClothSize.y);
		for (auto& edge : EdgesBufferArrCPU)
		{
			edge.clear();
		}

		//For pinned knots
		KnotsStatesBufferCPU.resize(CGConstantsCB.ClothSize.x * CGConstantsCB.ClothSize.y);
		std::fill(KnotsStatesBufferCPU.begin(), KnotsStatesBufferCPU.end(), 0);

		c_float3& startPos = CGConstantsCB.DefaultParticlePos;

		if (bEnableSelfCollisions)
		{
			check((startPos.x + CGConstantsCB.ClothSize.x * CGConstantsCB.SpringRestLength) < GetCollisionGridLengthMax(), "Cloth size should be less than BoundsMax: %f", GetCollisionGridLengthMax());
		}

		//Generate knots
		if (bGenerateClothHorizontal)//make it horizontal
		{
			for (uint y = 0; y < CGConstantsCB.ClothSize.y; y++)
			{
				for (uint x = 0; x < CGConstantsCB.ClothSize.x; x++)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					KnotsBufferCPU[i] = float3{ startPos.x + float(x) * CGConstantsCB.SpringRestLength, startPos.y, startPos.z + float(y) * CGConstantsCB.SpringRestLength };
				}
			}
		}
		else
		{
			for (uint y = 0; y < CGConstantsCB.ClothSize.y; y++)
			{
				for (uint x = 0; x < CGConstantsCB.ClothSize.x; x++)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					KnotsBufferCPU[i] = float3{ startPos.x + float(x) * CGConstantsCB.SpringRestLength, startPos.y + float(y) * CGConstantsCB.SpringRestLength, startPos.z };
				}
			}
		}
		
		

		if (bClothGenerationPinGeneratedObjects)
		{
			//Pin Objects
			if (1)//pin all borders
			{
				for (uint y = 0; y < CGConstantsCB.ClothSize.y; y++)
				{
					for (uint x = 0; x < CGConstantsCB.ClothSize.x; x++)
					{
						uint i = GetFlatIndexToClothArray(uint2{ x,y });
						if (y == (CGConstantsCB.ClothSize.y - 1) || x == (CGConstantsCB.ClothSize.x - 1) || (y == 0) || (x == 0))
						{
							KnotsStatesBufferCPU[i] = bEnableSoftPin ? BITMASK_PARTICLE_STATE_PINNED_SOFT : BITMASK_PARTICLE_STATE_PINNED;
						}
						else if(bClothSimulateWaterSurface)
						{
							KnotsStatesBufferCPU[i] = BITMASK_PARTICLE_STATE_PINNED_SOFT;
						}

					}
				}
			}
			else//pin only highest row
			{
				for (uint i = 0; i < CGConstantsCB.ClothSize.x; i++)
				{
					uint i1 = GetFlatIndexToClothArray(uint2{ i, CGConstantsCB.ClothSize.y - 1 });
					KnotsStatesBufferCPU[i1] = bEnableSoftPin ? BITMASK_PARTICLE_STATE_PINNED_SOFT : BITMASK_PARTICLE_STATE_PINNED;;
				}
			}
		}
		
		

		//Generate Blue edges (.x stride 2 offset 0)
		//Column major so that next constraint is above
		for (uint x = 0; x < CGConstantsCB.ClothSize.x; x += 2)
		{
			for (uint y = 0; y < CGConstantsCB.ClothSize.y; y++)
			{
				if ((x + 1) < CGConstantsCB.ClothSize.x)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					uint iRight = GetFlatIndexToClothArray(uint2{ x + 1, y });
					auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_BLUE];
					CConstraint newEdge;
					newEdge.IndexToKnot1 = i;
					newEdge.IndexToKnot2 = iRight;
					edgeArr.push_back(newEdge);
				}
			}
		}

		//Generate Green edges (.x stride 2 offset 1)
		//Column major so that next constraint is above
		for (uint x = 1; x < CGConstantsCB.ClothSize.x; x += 2)
		{
			for (uint y = 0; y < CGConstantsCB.ClothSize.y; y++)
			{
				if ((x + 1) < CGConstantsCB.ClothSize.x)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					uint iRight = GetFlatIndexToClothArray(uint2{ x + 1, y });
					auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_GREEN];
					CConstraint newEdge;
					newEdge.IndexToKnot1 = i;
					newEdge.IndexToKnot2 = iRight;
					edgeArr.push_back(newEdge);
				}
			}
		}

		//Generate Yellow edges (.y stride 2 offset 0)
		//Row major so that next constraint is to the right
		for (uint y = 0; y < CGConstantsCB.ClothSize.y; y += 2)
		{
			for (uint x = 0; x < CGConstantsCB.ClothSize.x; x++)
			{
				if ((y + 1) < CGConstantsCB.ClothSize.y)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					uint iUp = GetFlatIndexToClothArray(uint2{ x, y + 1 });
					auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_YELLOW];
					CConstraint newEdge;
					newEdge.IndexToKnot1 = i;
					newEdge.IndexToKnot2 = iUp;
					edgeArr.push_back(newEdge);
				}
			}
		}

		//Generate Red edges (.y stride 2 offset 1)
		//Row major so that next constraint is to the right
		for (uint y = 1; y < CGConstantsCB.ClothSize.y; y += 2)
		{
			for (uint x = 0; x < CGConstantsCB.ClothSize.x; x++)
			{
				if ((y + 1) < CGConstantsCB.ClothSize.y)
				{
					uint i = GetFlatIndexToClothArray(uint2{ x,y });
					uint iUp = GetFlatIndexToClothArray(uint2{ x, y + 1 });
					auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_RED];
					CConstraint newEdge;
					newEdge.IndexToKnot1 = i;
					newEdge.IndexToKnot2 = iUp;
					edgeArr.push_back(newEdge);
				}
			}
		}


		///Upload to GPU
		CmdList.UploadToBuffer(ParticlesBufferPositionCurGPU, 0, KnotsBufferCPU.data(), KnotsBufferCPU.size() * sizeof(float3));
		CmdList.UploadToBuffer(ParticlesBufferPositionPrevGPU, 0, KnotsBufferCPU.data(), KnotsBufferCPU.size() * sizeof(float3));

		CmdList.UploadToBuffer(ParticlesBufferPositionAfterCreationGPU, 0, KnotsBufferCPU.data(), KnotsBufferCPU.size() * sizeof(float3));

		CmdList.UploadToBuffer(ParticlesStateBufferGPU, 0, KnotsStatesBufferCPU.data(), KnotsStatesBufferCPU.size() * sizeof(uint));

		CGConstantsCB.NumParticles = KnotsBufferCPU.size();

		//Constraints
		for (size_t i = 0; i < CONSTRAINT_NUM_COLOR_TYPES; i++)
		{
			auto& edgeGroup = EdgesBufferArrCPU[i];
			auto& edgeGroupGPU = ConstraintsBufferArrGPU[i];
			auto& edgeGroupStatesGPU = ConstraintsStatesBufferArrGPU[i];

			if (!edgeGroup.empty())
			{
				if (edgeGroupGPU.GetElementCount() != edgeGroup.size())
				{
					edgeGroupGPU.Create("ConstraintsBuffer", edgeGroup.size());
					
					char bufName[256];
					sprintf_s(bufName, 256, "ConstraintsStatesBuffer %i ", i);

					edgeGroupStatesGPU.Create(bufName, edgeGroup.size());
				}

				CmdList.FillBuffer(edgeGroupStatesGPU, 0, 0, edgeGroup.size() * sizeof(uint));
				CmdList.UploadToBuffer(edgeGroupGPU, 0, edgeGroup.data(), edgeGroup.size() * sizeof(CConstraint));
			}
		}

		CGSpatialControlConstantsCB.IntersectionSphereScale = CGConstantsCB.SpringRestLength * 1.5f;

	}

	class ResolveConstraintsShader : public RShader
	{
	public:
		ResolveConstraintsShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ResolveConstraints", L"cs_6_0");
		}
	};

	class ResolveConstraintsToDefaultPosShader : public RShader
	{
	public:
		ResolveConstraintsToDefaultPosShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ResolveConstraintsToDefaultPos", L"cs_6_0");
		}
	};

	void ResolveConstraints(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ResolveConstraints, CmdList);

		///Bind Const Resources
		
		auto& RootSig = ParticlePhysicsRootSig::Get();
		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		static ResolveConstraintsShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		auto ResolveConstraintGroup = [&](int i, float curPassStringStiffness = 1.f)
		{
			//if (i == 3) continue;
			uint numConstraints = EdgesBufferArrCPU[i].size();
			if (numConstraints > 0)
			{
				GPU_PROFILE_SCOPE_TEXT(CmdList, L"Resolve Constraints Type %i: Num %i", i, numConstraints);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

				///Use Root Constants from Collision Pass
				GCollisionConstantsCB.CurCellLocalId = i;
				GCollisionConstantsCB.NumCells = numConstraints;
				GCollisionConstantsCB.PerPassSpringStiffnessScale = curPassStringStiffness;

				CmdList.SetConstantArray(RootSig.GetRootParamIndex("GCollisionConstantsCB"), 3, &GCollisionConstantsCB);

				auto& curConstraintBufferGPU = ConstraintsBufferArrGPU[i];
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsBufferSRV"), curConstraintBufferGPU, true);

				auto& curConstraintStatesBufferGPU = ConstraintsStatesBufferArrGPU[i];
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferSRV"), curConstraintStatesBufferGPU, true);

				constexpr uint ThreadGroupSize = 64;
				CmdList.Dispatch(RMath::DivideAndRoundUp(numConstraints, ThreadGroupSize), 1, 1);

			}
		};
		
		for (int j = 0; j < NumConstraintResolveIterations; j++)
		{
			///pass per each constraint group
			for (int i = 0; i < CONSTRAINT_NUM_COLOR_TYPES; i++)
			{
				ResolveConstraintGroup(i);
			}

			ResolveConstraintGroup(0, 0.5f);
			ResolveConstraintGroup(2, 0.5f);
		}

	}

	//Soft pin
	void ResolveConstraintsToDefaultPos(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ResolveConstraintsToDefaultPos, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();
		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		static ResolveConstraintsToDefaultPosShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionAfterCreationSRV"), ParticlesBufferPositionAfterCreationGPU, true);

		constexpr uint ThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.ClothSize.x * CGConstantsCB.ClothSize.y, ThreadGroupSize), 1, 1);
	}


	void DebugVisualizeConstraints(RCommandListCompute& CmdList)
	{
		///Export positions to CPU

		static RDynamicVector<float3> positionsBufferCPU;
		positionsBufferCPU.resize(CGConstantsCB.NumParticles);


		RReadbackBuffer positionsReadbackBuffer;
		positionsReadbackBuffer.Create(TEXT("PositionsReadbackBuffer"), CGConstantsCB.NumParticles, sizeof(float3));

		CmdList.CopyBufferRegion(positionsReadbackBuffer, 0, ParticlesBufferPositionCurGPU, 0, CGConstantsCB.NumParticles * sizeof(float3));
		CmdList.ExecuteCmdList(true);

		positionsReadbackBuffer.ExtractContents(positionsBufferCPU.data(), CGConstantsCB.NumParticles * sizeof(float3));

		///Visualize each constraint group

		///NB: Using Global here:
		auto& GImGeometry = RGeometry::GInstance;

		///Render Blue
		{
			auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_BLUE];
			GImGeometry.SetColor(EColor::Blue);
			for (auto& edge : edgeArr)
			{
				auto p1 = positionsBufferCPU[edge.IndexToKnot1];
				auto p2 = positionsBufferCPU[edge.IndexToKnot2];

				GImGeometry.AddLine(p1, p2);
			}
		}
		///Render Green
		{
			auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_GREEN];
			GImGeometry.SetColor(EColor::Green);
			for (auto& edge : edgeArr)
			{
				auto p1 = positionsBufferCPU[edge.IndexToKnot1];
				auto p2 = positionsBufferCPU[edge.IndexToKnot2];

				GImGeometry.AddLine(p1, p2);
			}
		}
		///Render Yellow]
		{
			auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_YELLOW];
			GImGeometry.SetColor(EColor::Yellow);
			for (auto& edge : edgeArr)
			{
				auto p1 = positionsBufferCPU[edge.IndexToKnot1];
				auto p2 = positionsBufferCPU[edge.IndexToKnot2];

				GImGeometry.AddLine(p1, p2);
			}
		}
		///Render Red
		{
			auto& edgeArr = EdgesBufferArrCPU[CONSTRAINT_COLOR_RED];
			GImGeometry.SetColor(EColor::Red);
			for (auto& edge : edgeArr)
			{
				auto p1 = positionsBufferCPU[edge.IndexToKnot1];
				auto p2 = positionsBufferCPU[edge.IndexToKnot2];

				GImGeometry.AddLine(p1, p2);
			}
		}


	}

	void DebugVisualizeSimulationBounds()
	{
		float initialBoundScale = CGConstantsCB.SimulationBoundsLength / 2.f;

		const float3 center = float3::MakeFloat3(CGConstantsCB.SimulationBoundsStartOffset + initialBoundScale);
		const float3 extents = float3::MakeFloat3(initialBoundScale);
		RGeometry::GInstance.SetColor(EColor::Red);
		RGeometry::GInstance.AddCubeLineList(center, extents);
	}

	///***************************
	///		SPATIAL CONTROL
	///***************************

	///CBuffer Constants
	SpatialControlConstants CGSpatialControlConstantsCB;
	TCBuffer<SpatialControlConstants> CGSpatialControlConstantsCBGPU;

	TStructuredBuffer<uint> IntersectedPointIndexGPU;
	TStructuredBuffer<float3> ClickedObjectViewSpaceZGPU;

	class UpdateIntersectedPointShader : public RShader
	{
	public:
		UpdateIntersectedPointShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "SpatialControl_UpdateIntersectedPoint", L"cs_6_0");
		}
	};

	uint IntersectedParticleIndex = UINT_MAX;
	float3 IntersectedParticlePos{};
	float ClickedObjectViewSpaceZ = 10.f;
	bool bIntersectedParticleHit = false;

	float3 CursorWorldPosPrev{};//TODO:Move it somewhere
	float3 CursorWorldPosCur{};

	float2 PrevMousePosNDC;
	float2 CurMousePosNDC;
	float2 CurMousePosNDCOffsetDir;

	void SpatialControl_Main(RCommandListCompute& CmdList)
	{
		CurMousePosNDC = RScene::Get().GetCamera().GetMouseCursorNDCCoords();
		CurMousePosNDCOffsetDir = CurMousePosNDC - PrevMousePosNDC;

		if (!CGSpatialControlConstantsCB.bInDragState)
		{
			SpatialControl_UpdateIntersectedPoint(CmdList);
		}
		
		///Debug Vis
		if (bIntersectedParticleHit)
		{
			RGeometry::GInstance.SetColor(EColor::RedOrange);
			RGeometry::GInstance.AddCircle(IntersectedParticlePos, CGSpatialControlConstantsCB.IntersectionSphereScale * CGConstantsCB.PointRadius);
		}

		///If not in Drag State, if Intersected Point is found, find ClickedObjectViewSpaceZ
		if (!CGSpatialControlConstantsCB.bInDragState)
		{
			if (bIntersectedParticleHit)
			{
				auto& viewMat = RScene::Get().GetCamera().Matrices.View;
				ClickedObjectViewSpaceZ = float3{ float3SSE{ DirectX::XMVector3Transform(float3SSE{ IntersectedParticlePos }, viewMat) } }.z;

				if (std::isnan(ClickedObjectViewSpaceZ))
				{
					ClickedObjectViewSpaceZ = float3{ float3SSE{ DirectX::XMVector3Transform(float3SSE{ CGConstantsCB.DefaultParticlePos }, viewMat) } }.z;
				}
			}
		}

		///Calculate World Space Mouse Pos from ClickedObjectViewSpaceZ
		CursorWorldPosCur = RScene::Get().GetCamera().GetMouseCursorPositionWorld(ClickedObjectViewSpaceZ);
		CGSpatialControlConstantsCB.ConstraintScissorPos = CursorWorldPosCur;
		CGSpatialControlConstantsCB.MouseOffsetWorld = CursorWorldPosCur - CursorWorldPosPrev;
		CursorWorldPosPrev = CursorWorldPosCur;

		///Register mouse click
		const bool bLClickThisFrame = Mouse::Get().IsLButtonPressed();
		const bool bRClickThisFrame = Mouse::Get().IsRButtonPressed();

		///Drag State Logic
		if (CGSpatialControlConstantsCB.bInDragState)
		{
			///Offset the object
			if (CGConstantsCB.bSimulateParticlesIn2D)
			{
				IntersectedParticlePos.x = IntersectedParticlePos.x + CGSpatialControlConstantsCB.MouseOffsetWorld.x;
				IntersectedParticlePos.y = IntersectedParticlePos.y + CGSpatialControlConstantsCB.MouseOffsetWorld.y;
			}
			else
			{
				IntersectedParticlePos = IntersectedParticlePos + CGSpatialControlConstantsCB.MouseOffsetWorld;
			}
			
			CmdList.UploadToBuffer(ParticlesBufferPositionCurGPU, IntersectedParticleIndex * sizeof(float3), &IntersectedParticlePos, sizeof(float3));
			
			///Disable Drag if Mouse was released
			if (!bLClickThisFrame)
			{
				CGSpatialControlConstantsCB.bInDragState = 0;
			}
		}
		else
		{
			///If Drag was disabled and we clicked -> start Dragging
			if (bLClickThisFrame && bIntersectedParticleHit)
			{
				CGSpatialControlConstantsCB.bInDragState = 1;
			}
		}

		if (bRClickThisFrame && (CGSpatialControlConstantsCB.bInDragState == 0))
		{
			///Cut constraints
			SpatialControl_CutConstraints(CmdList);
		}

		CGSpatialControlConstantsCB.PrevFrameIntersectedPointIndex = IntersectedParticleIndex;


		PrevMousePosNDC = CurMousePosNDC;
	}

	void SpatialControl_UpdateIntersectedPoint(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(UpdateIntersectedPoint, CmdList);

		///Update Constants
		{
			RScene::Get().GetCamera().GetMouseCursorRayWorldSpace(CGSpatialControlConstantsCB.RayOrigin, CGSpatialControlConstantsCB.RayEnd);

			CGSpatialControlConstantsCBGPU.Update(CGSpatialControlConstantsCB);
		}
		

		///Clear buffers
		//CmdList.ClearUAV(IntersectedPointIndexGPU);
		//CmdList.ClearUAV(ClickedObjectViewSpaceZGPU);
		CmdList.FillBuffer(IntersectedPointIndexGPU, 0, UINT_MAX, 4);
		CmdList.FillBuffer(ClickedObjectViewSpaceZGPU, 0, 0.f, 4 * 3);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSpatialControlConstantsCB"), CGSpatialControlConstantsCBGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferUAV"), ParticlesStateBufferGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("IntersectedPointIndexUAV"), IntersectedPointIndexGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ClickedObjectViewSpaceZUAV"), ClickedObjectViewSpaceZGPU, true);


		///
		///	Shader
		///
		static UpdateIntersectedPointShader ComputeShader;

		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize), 1, 1);


		///Readback

		RReadbackBuffer indexReadbackBuffer;
		RReadbackBuffer posReadbackBuffer;
		indexReadbackBuffer.Create(TEXT("indexReadbackBuffer"), 1, sizeof(uint));
		posReadbackBuffer.Create(TEXT("posReadbackBuffer"), 1, sizeof(float3));

		CmdList.CopyBufferRegion(indexReadbackBuffer, 0, IntersectedPointIndexGPU, 0, sizeof(uint));
		CmdList.CopyBufferRegion(posReadbackBuffer, 0, ClickedObjectViewSpaceZGPU, 0, sizeof(float3));
		CmdList.ExecuteCmdList(true);

		indexReadbackBuffer.ExtractContents(&IntersectedParticleIndex, sizeof(uint));
		posReadbackBuffer.ExtractContents(&IntersectedParticlePos, sizeof(float3));

		bIntersectedParticleHit = IntersectedParticleIndex != UINT_MAX;

	}

	class CutConstraintShader : public RShader
	{
	public:
		CutConstraintShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "SpatialControl_CutConstraints", L"cs_6_0");
		}
	};

	void SpatialControl_CutConstraints(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(CutConstraints, CmdList);

		///Update Constants
		{
			CGSpatialControlConstantsCB.ConstraintScissorPos = CursorWorldPosCur;
			CGSpatialControlConstantsCBGPU.Update(CGSpatialControlConstantsCB);
		}

		auto& RootSig = ParticlePhysicsRootSig::Get();
		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSpatialControlConstantsCB"), CGSpatialControlConstantsCBGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		static CutConstraintShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		auto ResolveConstraintGroup = [&](int i, float curPassStringStiffness = 1.f)
		{
			//if (i == 3) continue;
			uint numConstraints = EdgesBufferArrCPU[i].size();
			if (numConstraints > 0)
			{
				GPU_PROFILE_SCOPE_TEXT(CmdList, L"Resolve Constraints Type %i: Num %i", i, numConstraints);

				///Use Root Constants from Collision Pass
				GCollisionConstantsCB.CurCellLocalId = i;
				GCollisionConstantsCB.NumCells = numConstraints;
				GCollisionConstantsCB.PerPassSpringStiffnessScale = curPassStringStiffness;

				CmdList.SetConstantArray(RootSig.GetRootParamIndex("GCollisionConstantsCB"), 3, &GCollisionConstantsCB);

				auto& curConstraintBufferGPU = ConstraintsBufferArrGPU[i];
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsBufferSRV"), curConstraintBufferGPU, true);

				auto& curConstraintStatesBufferGPU = ConstraintsStatesBufferArrGPU[i];
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferUAV"), curConstraintStatesBufferGPU, true);

				constexpr uint ThreadGroupSize = 64;
				CmdList.Dispatch(RMath::DivideAndRoundUp(numConstraints, ThreadGroupSize), 1, 1);

			}
		};

		for (int j = 0; j < NumConstraintResolveIterations; j++)
		{
			//pass per each constraint group
			for (int i = 0; i < CONSTRAINT_NUM_COLOR_TYPES; i++)
			{
				ResolveConstraintGroup(i);
			}
		}
	}


	class ApplyPointerVelocityToParticlesShader : public RShader
	{
	public:
		ApplyPointerVelocityToParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ApplyPointerVelocityToParticles", L"cs_6_0");
		}
	};

	void ApplyPointerVelocityToParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ApplyPointerVelocityToParticles, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		CGSpatialControlConstantsCBGPU.Update(CGSpatialControlConstantsCB);
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSpatialControlConstantsCB"), CGSpatialControlConstantsCBGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferUAV"), CurParticlesVelocityBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		static ApplyPointerVelocityToParticlesShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint kThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, kThreadGroupSize), 1, 1);

	}

	///***************************
	///		MESH GENERATOR
	///***************************
	
	bool bGenerateTriangleMesh = false;
	bool bGenerateTriangleMeshNormals = false;
	
	uint MeshGeneratorCurNumVertices;

	GPUSurfaceGenerationSettings SurfaceSubdivisionSettings;
	
	template<bool bDiscardBrokenConstraints>
	class MeshGenerationShader : public RShader
	{
	public:
		MeshGenerationShader()
		{
			if constexpr (bDiscardBrokenConstraints)
			{
				PushBackDefine({ L"DISCARD_BROKEN_CONSTRAINTS", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "MeshGenerator", L"cs_6_0");
		}
	};

	TStructuredBuffer<float3> MeshBufferVerticesGPU;
	TStructuredBuffer<float3> MeshBufferNormalsGPU;
	TStructuredBuffer<float2> MeshBufferTexCoordsGPU;
	TStructuredBuffer<uint> TriangleCountBufferGPU;
	RTexture* pGeneratedMeshColorTexture;
	RTexture* pGeneratedMeshNormalTexture;

	uint2 OutputMeshSize;

	bool bMeshGeneratorDiscardBrokenConnections = false;

	bool bMeshGeneratorUseColorTexture = true;

	//Tesellation
	bool bSubdivisionMesh = false;//for non-destructible surfaces
	TStructuredBuffer<float3> InterpolatedParticlePositionsBufferGPU;
	uint SubdivisionScale = 16;

	void MeshGeneratorMain(RCommandListCompute& CmdList)
	{
		RStructuredBuffer* pInputBuffer = &ParticlesBufferPositionCurGPU;
		OutputMeshSize = CGConstantsCB.ClothSize;
		bool bDiscardBrokenConstraints = bMeshGeneratorDiscardBrokenConnections;

		//Additional topology by Interpolating 
		if (bSubdivisionMesh)
		{
			OutputMeshSize.x = CGConstantsCB.ClothSize.x * SubdivisionScale;
			OutputMeshSize.y = CGConstantsCB.ClothSize.y * SubdivisionScale;
			GPUSurfaceGeneration(CmdList, ParticlesBufferPositionCurGPU, CGConstantsCB.ClothSize, OutputMeshSize, InterpolatedParticlePositionsBufferGPU, SurfaceSubdivisionSettings);
			pInputBuffer = &InterpolatedParticlePositionsBufferGPU;
			bDiscardBrokenConstraints = false;
		}

		if (bGenerateTriangleMesh)
		{
			MeshGeneratorTriangulate(CmdList, *pInputBuffer, OutputMeshSize, bDiscardBrokenConstraints);

			if (bClothSimulateWaterSurface && bMeshGeneratorUseColorTexture)
			{
				if (bEnableVorticityConfinement && bApplyVorticityMaskToDensityVisualizer)
				{
					//pGeneratedMeshColorTexture = &DensityTextureCopyGPU;
					pGeneratedMeshColorTexture = &DensityTextureGPU;
					pGeneratedMeshNormalTexture = &CulGradientDirectionTexture;
				}
				else
				{
					pGeneratedMeshColorTexture = &DensityTextureGPU;
				}
				
				
			}
		}

	}

	void MeshGeneratorTriangulate(RCommandListCompute& CmdList, RStructuredBuffer& inputSamplesGPU, uint2 numInputSamples, bool bDiscardBrokenConstraints = true)
	{
		GPU_PROFILE_SCOPE(Generate_Triangle_Mesh, CmdList);

		//Compute min num of triangles required
		const uint MaxNumTriangles = (numInputSamples.x - 1) * (numInputSamples.y - 1) * 2;
		const uint MaxNumVertices = MaxNumTriangles * 3;

		if (MeshBufferVerticesGPU.GetElementCount() < MaxNumVertices)
		{
			MeshBufferVerticesGPU.Create(TEXT("GeneratedMeshBufferVertices"), MaxNumVertices);
			MeshBufferNormalsGPU.Create(TEXT("GeneratedMeshBufferNormals"), MaxNumVertices);
			MeshBufferTexCoordsGPU.Create(TEXT("GeneratedMeshBufferTexCoords"), MaxNumVertices);

			TriangleCountBufferGPU.Create(TEXT("TriangleCountBuffer"), 1);
		}

		//Clear Counter
		CmdList.FillBuffer(TriangleCountBufferGPU, 0, 0, 4);

		//Update Constants
		{
			CGSpatialControlConstantsCBGPU.Update(CGSpatialControlConstantsCB);
		}

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		GMeshGenerationConstantsCB.MeshSize = numInputSamples;

		CmdList.SetConstantArray(RootSig.GetRootParamIndex("GMeshGenerationConstantsCB"), 2, &GMeshGenerationConstantsCB);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), inputSamplesGPU, true);

		if (bDiscardBrokenConstraints)
		{
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferBLUESRV"), ConstraintsStatesBufferArrGPU[CONSTRAINT_COLOR_BLUE], true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferGREENSRV"), ConstraintsStatesBufferArrGPU[CONSTRAINT_COLOR_GREEN], true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferYELLOWSRV"), ConstraintsStatesBufferArrGPU[CONSTRAINT_COLOR_YELLOW], true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ConstraintsStatesBufferREDSRV"), ConstraintsStatesBufferArrGPU[CONSTRAINT_COLOR_RED], true);
		}
		
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("InterlockedAddBuffer"), TriangleCountBufferGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("MeshBufferVerticesUAV"), MeshBufferVerticesGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("MeshBufferNormalsUAV"), MeshBufferNormalsGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("MeshBufferTexCoordsUAV"), MeshBufferTexCoordsGPU, true);

		/*Shader*/
		if (bDiscardBrokenConstraints)
		{
			static MeshGenerationShader<true> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		else
		{
			static MeshGenerationShader<false> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		

		constexpr uint ThreadGroupSize = 8;
		CmdList.Dispatch(numInputSamples.x / ThreadGroupSize, numInputSamples.y / ThreadGroupSize, 1);


		//Readback
		RReadbackBuffer countReadbackBuffer;
		countReadbackBuffer.Create(TEXT("countReadbackBuffer"), 1, sizeof(uint));

		CmdList.CopyBufferRegion(countReadbackBuffer, 0, TriangleCountBufferGPU, 0, sizeof(uint));
		CmdList.ExecuteCmdList(true);

		countReadbackBuffer.ExtractContents(&MeshGeneratorCurNumVertices, sizeof(uint));
	}


	///*********************************
	///		VELOCITY FIELD GENERATOR
	///*********************************


	bool bVelocityField2DVerticalSliceSimulation = true;

	bool bClearVectorField = false;

	int bProcessVectorFieldOnCPU = 0;
	bool bGenerateVectorFieldOnCPU = false;
	bool bEnableVectorFieldProjection = true;
	bool bEnableVectorFieldProjectionAfterAdvection = false;
	bool bProjectVectorFieldImplicit = false;
	bool bOutputDivergenceAndPressureResults = true;//implicitly compute divergence and pressure anyway, regardless of projection method used
	bool bEnableAdvection = true;
	bool bEnableVelocityAdvection = true;
	uint VectorFieldProjectionNumIterations = 50;

	float2 VelocityFieldForceAtBorders = float2{ 0, 0 };

	bool bDebugVisualizeVelocityField = false;
	bool bDebugVisualizeVelocityFieldCells = false;
	bool bDebugVisualizeVelocityFieldCellsStates = false;
	bool bDebugVisualizeVelocityFieldSeparateAxises = false;
	bool bDebugVisualizeVelocityFieldMagnitudes = false;
	bool bDebugVisualizeVectorFieldStreamlines = false;
	bool bVisualizeDensityField = true;

	bool bDebugVisualizeDivergence = false;
	bool bDebugVisualizeCurl = false;
	bool bDebugVisualizeCurlGradient = false;
	bool bDebugVisualizePressure = false;
	bool bDebugVisualizePressureGradient = false;

	float VelocityFieldDamping = 0.99f;
	int VelocityGeneratorTriggerRadius = 5;

	int DensityGeneratorTriggerRadius = 2;

	 /*Resources*/
	//TStructuredBuffer<float2> VelocityFieldBufferGPU;
	//TStructuredBuffer<float2> VelocityFieldBufferCopyGPU;
	
	//Each Velocity component stored separately
	RTexture VelocityField_U_ComponentTextureGPU{ TEXT("VelocityField_U_ComponentTexture")};
	RTexture VelocityField_V_ComponentTextureGPU{ TEXT("VelocityField_V_ComponentTexture") };
	RTexture VelocityField_U_ComponentTexture2GPU{ TEXT("VelocityField_U_ComponentTexture2") };//double buffering
	RTexture VelocityField_V_ComponentTexture2GPU{ TEXT("VelocityField_V_ComponentTexture2") };//double buffering
	RTexture VelocityField_U_ComponentTextureCopyGPU{ TEXT("VelocityField_U_ComponentTextureCopy") };
	RTexture VelocityField_V_ComponentTextureCopyGPU{ TEXT("VelocityField_V_ComponentTextureCopy") };

	RTexture VelocityField_U_ComponentTextureRenderTargetGPU{ TEXT("VelocityField_U_ComponentTextureRenderTarget") };
	RTexture VelocityField_V_ComponentTextureRenderTargetGPU{ TEXT("VelocityField_V_ComponentTextureRenderTarget") };

	//For particle-based advection
	RTexture VelocityField_U_ComponentWeightTextureGPU{ TEXT("VelocityField_U_ComponentWeightTexture") };
	RTexture VelocityField_V_ComponentWeightTextureGPU{ TEXT("VelocityField_V_ComponentWeightTexture") };

	RTexture VelocityField_U_ComponentWeightTextureRenderTargetGPU{ TEXT("VelocityField_U_ComponentWeightTextureRenderTarget") };
	RTexture VelocityField_V_ComponentWeightTextureRenderTargetGPU{ TEXT("VelocityField_V_ComponentWeightTextureRenderTarget") };

	RTexture AdvectedParticleDensityTextureGPU{ TEXT("AdvectedParticleDensityTexture") };
	RTexture AdvectedParticleDensityTextureRenderTargetGPU{ TEXT("AdvectedParticleDensityTextureRenderTarget") };

	RDynamicVector<VELOCITY_FIELD_STATE_TYPE> VelocityFieldStateBufferCPU;
	RTexture VelocityFieldStateTextureGPU{ TEXT("VelocityFieldStateTexture") };

	RDynamicVector<float> VelocityField_U_ComponentBufferCPU;
	RDynamicVector<float> VelocityField_V_ComponentBufferCPU;

	RDynamicVector<float> VelocityField_U_ComponentBufferCopyCPU;
	RDynamicVector<float> VelocityField_V_ComponentBufferCopyCPU;

	RDynamicVector<float> DivergenceBufferCPU;
	RTexture DivergenceTextureGPU{ TEXT("DivergenceTexture") };

	RTexture CurlTextureGPU{ TEXT("CurlTexture") };
	RTexture CulGradientTexture{ TEXT("CurlGradientTexture") };
	RTexture CulGradientDirectionTexture{ TEXT("CurlGradientDirectionTexture") };
	RTexture CulGradientMinMaxTexture{ TEXT("CurlGradientMinMaxTexture") };

	RTexture VelocityFieldMagnitudesTextureGPU{ TEXT("MagnitudesTexture") };

	//Density
	RDynamicVector<DENSITY_TYPE> DensityBufferCPU;
	RTexture DensityTextureGPU{ TEXT("DensityFieldTexture") };
	RTexture DensityTexture2GPU{ TEXT("DensityFieldTexture2") };//double buffering
	RTexture DensityTexture3GPU{ TEXT("DensityFieldTexture3") };//triple buffering
	RTexture DensityTextureCopyGPU{ TEXT("DensityFieldTextureCopy") };

	RTexture PressureTextureGPU{ TEXT("PressureTexture") };
	RTexture PressureGradientTextureGPU{ TEXT("PressureTexture") };
	RTexture PressureTexture2GPU{ TEXT("PressureTexture2") };

	RTexture PressureTextureColoredGPU{ TEXT("PressureVisualizeTexture") };

	uint GetFlatIndexToVelocityArray2D(uint2 coord)
	{
		return coord.y * CGConstantsCB.VelocityFieldGridSize + coord.x;
	};

	bool bEnableSphereObstacle = true;
	bool bObstacleAsDensitySource = true;
	bool bObstacleGeneratesVelocity = true;

	float3 VelocityFieldObstaclePos{ 12,12, 0 };
	float2 VelocityFieldObstaclePosPrev{ 12,12 };
	float VelocityFieldObstacleRadius{1};

	float VelocityFieldDampingStartThreshold = 150.f;
	
	void VelocityFieldApplyDamping()
	{
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				auto curVelocity = LoadFromVelocityField(coord);
				if (VectorGetLengthSquared(curVelocity) > VelocityFieldDampingStartThreshold * VelocityFieldDampingStartThreshold)
				{
					StoreInVelocityField(coord, curVelocity * VelocityFieldDamping);
				}
			}
		}
	}

	bool bEnableDensity = true;
	bool bClearDensityField = false;
	bool bGenerateDensityOnBorder = false;

	void GenerateDensity()
	{
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				//generate density at near-border
				if ((y == 1 || y == 2) && (x > 0 && x < (CGConstantsCB.VelocityFieldGridSize - 1)))
				{
					DensityBufferCPU[i] = CGConstantsCB.InitialDensityColor;
				}
			}
		}
	}

	void ComputeCurVelocityFieldCellIdMousePointsTo()
	{
		//test for each cell if mouse ray intersects
		float3 rayOrigin;
		float3 rayEnd;
		float3 rayDirection;
		RScene::Get().GetCamera().GetMouseCursorRayWorldSpace(rayOrigin, rayEnd);
		rayDirection = VectorNormalize(rayEnd - rayOrigin);

		const float2 cellHalfSize = float2::MakeFloat2(CGConstantsCB.VelocityFieldGridCellLength * 0.5f);

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				
				float3 cellPosWorld{ GetVelocityCellCenterWorldSpaceFromCoord(coord), 0 };

				if (!bVelocityField2DVerticalSliceSimulation)
				{
					cellPosWorld.z = cellPosWorld.y;
					cellPosWorld.y = 0;
				}

				float tDummy;
				//if (RayAABBIntersection(rayOrigin, rayDirection, cellPosWorld, float3{ cellHalfSize, 0 }, tDummy))
				if (RaySphereIntersection(rayOrigin, rayDirection, cellPosWorld, cellHalfSize.x))
				{
					CGConstantsCB.CurVelocityFieldCellIdMousePointsTo = coord;
					return;
				}
			}
		}
	}

	void ProcessVelocityField2D(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(Process_Velocity_Field, CmdList);

		if (bClothSimulateWaterSurface)
		{
			bVelocityField2DVerticalSliceSimulation = false;
			CGConstantsCB.VelocityFieldBoundsLength = (CGConstantsCB.ClothSize.x - 1) * CGConstantsCB.SpringRestLength;
			CGConstantsCB.VelocityFieldBoundsStartOffset = CGConstantsCB.DefaultParticlePos.x;
		}

		auto curBufSize = VelocityField_U_ComponentBufferCPU.size();

		uint numVelocityGridCells = CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize;

		if (curBufSize < numVelocityGridCells)
		{
			VelocityField_U_ComponentBufferCPU.resize(numVelocityGridCells);
			VelocityField_V_ComponentBufferCPU.resize(numVelocityGridCells);

			VelocityField_U_ComponentBufferCopyCPU.resize(numVelocityGridCells);
			VelocityField_V_ComponentBufferCopyCPU.resize(numVelocityGridCells);

			VelocityFieldStateBufferCPU.resize(numVelocityGridCells);
			
			DensityBufferCPU.resize(numVelocityGridCells);
		}

		 /*GPU Resources*/
		if (VelocityField_U_ComponentTextureGPU.Size.x != CGConstantsCB.VelocityFieldGridSize)
		{
			/*VelocityFieldBufferGPU.Create("GeneratedVelocityFieldBuffer", numVelocityGridCells);
			VelocityFieldBufferCopyGPU.Create("GeneratedVelocityFieldBufferCopy", numVelocityGridCells);*/

			VelocityField_U_ComponentTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_V_ComponentTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_U_ComponentTexture2GPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_V_ComponentTexture2GPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_U_ComponentTextureCopyGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_V_ComponentTextureCopyGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);

			VelocityField_U_ComponentTextureRenderTargetGPU.CreateRenderTarget(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT);
			VelocityField_V_ComponentTextureRenderTargetGPU.CreateRenderTarget(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT);

			VelocityField_U_ComponentWeightTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			VelocityField_V_ComponentWeightTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);

			VelocityField_U_ComponentWeightTextureRenderTargetGPU.CreateRenderTarget(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT);
			VelocityField_V_ComponentWeightTextureRenderTargetGPU.CreateRenderTarget(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT);

			AverageDensityBufferGPU.Create(TEXT("AverageDensityBuffer"), 1);

			AdvectedParticleDensityTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			AdvectedParticleDensityTextureRenderTargetGPU.CreateRenderTarget(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT);

			VelocityFieldStateTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, std::is_same<VELOCITY_FIELD_STATE_TYPE, uint8>::value ? DXGI_FORMAT_R8_UINT : DXGI_FORMAT_R32_UINT, true);

			DivergenceTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			CurlTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			CulGradientTexture.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			CulGradientDirectionTexture.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32G32B32A32_FLOAT, true);
			
#if PAINT_RENDER_PIXEL_ADVECT
			DensityTextureGPU.Create(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R11G11B10_FLOAT, true);
			DensityTexture2GPU.Create(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R11G11B10_FLOAT, true);
			DensityTexture3GPU.Create(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R11G11B10_FLOAT, true);
			DensityTextureCopyGPU.Create(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R11G11B10_FLOAT, true);
#else
			DensityTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, USE_3COLOR_DENSITY ? DXGI_FORMAT_R32G32B32A32_FLOAT : DXGI_FORMAT_R32_FLOAT, true);
			DensityTexture2GPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, USE_3COLOR_DENSITY ? DXGI_FORMAT_R32G32B32A32_FLOAT : DXGI_FORMAT_R32_FLOAT, true);
			DensityTexture3GPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, USE_3COLOR_DENSITY ? DXGI_FORMAT_R32G32B32A32_FLOAT : DXGI_FORMAT_R32_FLOAT, true);
			DensityTextureCopyGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, USE_3COLOR_DENSITY ? DXGI_FORMAT_R32G32B32A32_FLOAT : DXGI_FORMAT_R32_FLOAT, true);
#endif
			
			PressureTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);
			PressureGradientTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32G32_FLOAT, true);
			PressureTexture2GPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R32_FLOAT, true);

			//VISUALIZERS
			VelocityFieldMagnitudesTextureGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R11G11B10_FLOAT, true);
			PressureTextureColoredGPU.Create(CGConstantsCB.VelocityFieldGridSize, CGConstantsCB.VelocityFieldGridSize, DXGI_FORMAT_R11G11B10_FLOAT, true);
		}

		//Compute cell id we point to
		ComputeCurVelocityFieldCellIdMousePointsTo();
		CGConstantsCB.CurPointerVelocity = CurMousePosNDCOffsetDir * CGConstantsCB.GeneratedVelocityMagnitudeScale * 100.f;

		//generate velocities on cpu
		if (!CGConstantsCB.bUseParticleBasedAdvection)
		{
			if (bGenerateVectorFieldOnCPU)
			{
				//Clean Velocity Field
				for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
				{
					for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
					{
						uint2 coord{ x,y };
						uint i = GetFlatIndexToVelocityArray2D(coord);

						if (bClearVectorField)
						{
							StoreInVelocityField(i, float2{ 0,0 });
						}

						if (bClearDensityField)
						{
#if USE_3COLOR_DENSITY == 1
							DensityBufferCPU[i] = float3{ 0,0,0 };
#else
							DensityBufferCPU[i] = 0;
#endif//USE_3COLOR_DENSITY

						}

					}
				}

				if (1 && Keyboard::IsKeyPressed('V'))//Mouse controlled velocity
				{
					//Write cur velocity

					int2 cellIdMin{ std::max(0, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.x) - VelocityGeneratorTriggerRadius), std::max(0, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.y) - VelocityGeneratorTriggerRadius) };
					int2 cellIdMax{ std::min(int(CGConstantsCB.VelocityFieldGridSize) - 1, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.x + VelocityGeneratorTriggerRadius)), std::min(int(CGConstantsCB.VelocityFieldGridSize) - 1, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.y + VelocityGeneratorTriggerRadius)) };
					int2 curCellId;

					for (curCellId.x = cellIdMin.x; curCellId.x <= cellIdMax.x; curCellId.x++)
					{
						for (curCellId.y = cellIdMin.y; curCellId.y <= cellIdMax.y; curCellId.y++)
						{
							auto v = LoadFromVelocityField(uint2{ uint(curCellId.x), uint(curCellId.y) });
							StoreInVelocityField(uint2{ uint(curCellId.x), uint(curCellId.y) }, v + CGConstantsCB.CurPointerVelocity);
						}
					}
				}

				if (bEnableDensity && Keyboard::IsKeyPressed('N'))
				{
					//generate density at a point

					int2 cellIdMin{ std::max(0, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.x) - DensityGeneratorTriggerRadius), std::max(0, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.y) - DensityGeneratorTriggerRadius) };
					int2 cellIdMax{ std::min(int(CGConstantsCB.VelocityFieldGridSize) - 1, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.x + DensityGeneratorTriggerRadius)), std::min(int(CGConstantsCB.VelocityFieldGridSize) - 1, int(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo.y + DensityGeneratorTriggerRadius)) };
					int2 curCellId;

					for (curCellId.x = cellIdMin.x; curCellId.x <= cellIdMax.x; curCellId.x++)
					{
						for (curCellId.y = cellIdMin.y; curCellId.y <= cellIdMax.y; curCellId.y++)
						{
							DensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ uint(curCellId.x), uint(curCellId.y) })] = CGConstantsCB.InitialDensityColor;
						}
					}
				}

				//Generate Bounds & Obstacles
				{
					VelocityFieldSetBoundaries();

					if (bEnableSphereObstacle || GbFirstBoot)
					{
						VelocityFieldSetObstacle();
					}
				}

				if ((VelocityFieldForceAtBorders.x != 0) || (VelocityFieldForceAtBorders.y != 0))
				{
					ApplyConstForceToVelocityFieldBorders();
				}

				ApplyConstExternalForceToVectorField();

				if (bClearVelocityFromSolidCell)
				{
					ClearVelocityFromSolidCells();
				}

				if (bGenerateDensityOnBorder)
				{
					GenerateDensity();
				}

			}
			else//GPU Generation
			{
				//Clear
				{
					if (bClearVectorField)
					{
						CmdList.TransitionResource(VelocityField_U_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
						CmdList.TransitionResource(VelocityField_V_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
						CmdList.ClearUAV(VelocityField_U_ComponentTextureGPU);
						CmdList.ClearUAV(VelocityField_V_ComponentTextureGPU);
					}

					if (bClearDensityField)
					{
						CmdList.TransitionResource(DensityTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
						CmdList.ClearUAV(DensityTextureGPU);
					}
				}


				if (bFillDensityFieldlWithColor)
				{
					FillDensityFieldGPU(CmdList);
					bFillDensityFieldlWithColor = false;
				}

				GenerateVelocityFieldBoundariesGPU(CmdList);

				GenerateVelocityFieldGPU(CmdList);

			}
		}
		

		//Upload CPU generated Data to GPU
		if (bGenerateVectorFieldOnCPU)
		{
			GPU_PROFILE_SCOPE(UploadFromCPUToGPU, CmdList);
			//Upload CPU Generated Velocity to Textures
			static RDynamicVector<float> intermediateBuffer1;
			CopyFromInvertedCPUBufferToGPUTexture(CmdList, VelocityField_U_ComponentBufferCPU, intermediateBuffer1, VelocityField_U_ComponentTextureGPU);
			static RDynamicVector<float> intermediateBuffer2;
			CopyFromInvertedCPUBufferToGPUTexture(CmdList, VelocityField_V_ComponentBufferCPU, intermediateBuffer2, VelocityField_V_ComponentTextureGPU);

			static RDynamicVector<VELOCITY_FIELD_STATE_TYPE> intermediateBuffer3;
			CopyFromInvertedCPUBufferToGPUTexture(CmdList, VelocityFieldStateBufferCPU, intermediateBuffer3, VelocityFieldStateTextureGPU);

			//Upload Density
			{
#if USE_3COLOR_DENSITY == 1
#else 
				static RDynamicVector<float> intermediateBuffer4;
				CopyFromInvertedCPUBufferToGPUTexture(CmdList, DensityBufferCPU, intermediateBuffer4, DensityTextureGPU);
#endif//USE_3COLOR_DENSITY
			}
		}

		//*************************************
		//			Process on CPU
		//*************************************
		if ((bProcessVectorFieldOnCPU == 1) && (CGConstantsCB.bUseParticleBasedAdvection || bGenerateVectorFieldOnCPU))
		{
			if (VelocityFieldDamping != 1.f)
			{
				VelocityFieldApplyDamping();
			}
			if (bDiffuseVelocity)
			{
				DiffuseVelocityCPU();
			}

			//1. Project
			if (bEnableVectorFieldProjection)
			{
				if (bProjectVectorFieldImplicit)
				{
					ProjectPressureCPU();
				}
				else
				{
					for (uint i = 0; i < VectorFieldProjectionNumIterations; i++)
					{
						//TODO: Exit when max & min divergence is close to 0
						ProjectVectorFieldCPU(true);
					}
				}
				
			}

			ComputeDivergenceCPU();
			
			if (bCopyVelocitiesToBorders || bCopyReflectedVelocitiesToBorders)
			{
				ExtrapolateVelocitiesToBordersCPU();
			}

			if (bClearVelocityFromSolidCell)
			{
				ClearVelocityFromSolidCells();
			}

			if (bDiffuseDensity)
			{
				DiffuseDensityCPU();
			}

			//2. Advect
			if (bEnableAdvection && !CGConstantsCB.bUseParticleBasedAdvection)
			{
				if (bEnableVelocityAdvection)
				{
					AdvectVelocityFieldCPU();

					if (bClearVelocityFromSolidCell)
					{
						ClearVelocityFromSolidCells();
					}
				}

				if (bEnableDensity)
				{
					AdvectDensityFieldCPU();
				}
			}

			if (bEnableVectorFieldProjection && bEnableVectorFieldProjectionAfterAdvection)
			{
				if (bProjectVectorFieldImplicit)
				{
					ProjectPressureCPU();
				}
				else
				{
					for (uint i = 0; i < VectorFieldProjectionNumIterations; i++)
					{
						ProjectVectorFieldCPU(false);
					}
				}

				if (bClearVelocityFromSolidCell)
				{
					ClearVelocityFromSolidCells();
				}

			}

			

		}
		
		

		//*************************************
		//			Process on GPU
		//*************************************
		if (bProcessVectorFieldOnCPU == 0)
		{

			if (bDiffuseVelocity)
			{
				DiffuseVelocity(CmdList);
			}

			if (bOutputDivergenceAndPressureResults)
			{
				//Compute Divergence
				ComputeDivergence(CmdList);

				//Compute Pressure
				ComputePressurePoisson(CmdList);
			}
			if (bEnableVectorFieldProjection)
			{
				if (bProjectVectorFieldImplicit)
				{
					ProjectVectorFieldImplicitMethod(CmdList);
				}
				else
				{
					GPU_PROFILE_SCOPE(ProjectVectorFieldSimple, CmdList);
					for (uint i = 0; i < VectorFieldProjectionNumIterations; i++)
					{
						GPU_PROFILE_SCOPE_TEXT(CmdList, L"Pass %i", i);
						ProjectVectorFieldSimple(CmdList);
					}
				}
			}
			
			if (bEnableAdvection && bEnableVelocityAdvection && !CGConstantsCB.bUseParticleBasedAdvection)
			{
				AdvectVelocity(CmdList);
			}

			Vorticity(CmdList);

			if (bEnableVectorFieldProjection && bEnableVectorFieldProjectionAfterAdvection)
			{
				if (bOutputDivergenceAndPressureResults)
				{
					//Compute Divergence
					ComputeDivergence(CmdList);

					//Compute Pressure
					ComputePressurePoisson(CmdList);
				}
				if (bProjectVectorFieldImplicit)
				{
					ProjectVectorFieldImplicitMethod(CmdList);
				}
				else
				{
					GPU_PROFILE_SCOPE(ProjectVectorFieldAfterAdvection, CmdList);
					for (uint i = 0; i < (VectorFieldProjectionNumIterations / 2); i++)
					{
						GPU_PROFILE_SCOPE_TEXT(CmdList, L"Pass %i", i);
						ProjectVectorFieldSimple(CmdList);
					}
				}
			}

			if (bEnableDensity && !CGConstantsCB.bUseParticleBasedAdvection)
			{
				if (bDiffuseDensity)
				{
					DiffuseDensity(CmdList);
				}

				if (bEnableAdvection)
				{
					AdvectDensity(CmdList);
				}

				if (bEnableVorticityConfinement && bApplyVorticityMaskToDensityVisualizer)
				{
					DensityApplyVorticityGradientMask(CmdList);
				}
			}

			

			
		}


		//Extract all GPU Data back to CPU
		if (bGenerateVectorFieldOnCPU && (bProcessVectorFieldOnCPU == 0))
		{
			//ExtractVelocityFieldDataFromGPU(CmdList);
		}


		/******************************
			Debug Visualization Stage
		*******************************/
		{
			
			//Calculate min/max magnitudes of all vectors
			ComputeVelocityFieldMagnitudesMinMax();

			//Debug Visualize
			if (bDebugVisualizeVelocityField)
			{
				DebugVisualizeVelocityField(CmdList, bDebugVisualizeVelocityFieldSeparateAxises);

				//DebugVisualizeAdvectionAtAPoint();
			}

			if (bDebugVisualizeVectorFieldStreamlines)
			{
				DebugVisualizeVectorFieldStreamLines();
			}

			if (bVisualizeDensityField)
			{
				VisualizeDensity();
			}
			else if (bDebugVisualizeDivergence)
			{
				VisualizeDivergence(CmdList);
			}
			else if (bDebugVisualizeCurl)
			{
				VisualizeCurl(CmdList);
			}
			else if (bDebugVisualizeCurlGradient)
			{
				//VisualizeVelocityFieldTexture(CulGradientTexture);
				VisualizeVelocityFieldTexture(CulGradientDirectionTexture);
			}
			else if (bDebugVisualizePressure)
			{
				VisualizePressure(CmdList);
			}
			else if (bDebugVisualizePressureGradient)
			{
				VisualizeVelocityFieldTexture(PressureGradientTextureGPU);
			}
			else if (bDebugVisualizeVelocityFieldMagnitudes)
			{
				VisualizeVelocityFieldMagnitudes();
			}

			//Visualize pointer
			if(0)
			{
				RGeometry::GInstance.SetColor(EColor::RedDark);
				const auto cellPos = GetVelocityCellCenterWorldSpaceFromCoord(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo);
				const float half = CGConstantsCB.VelocityFieldGridCellLength * 0.5f;
				
				if (bVelocityField2DVerticalSliceSimulation)
				{
					RGeometry::GInstance.AddPlane(float3{ cellPos, 0 }, float2{ half, half });
				}
				else
				{
					RGeometry::GInstance.AddPlaneXZ(float3{ cellPos.x, 0, cellPos.y }, float2{ half, half });
				}

			}
		}

	}

	float3 MapTo3ColorRange(float curValue, float valueMin, float valueMax, const float3& color1 = GetColorRGBFloat(EColor::Blue), const float3& color2 = GetColorRGBFloat(EColor::Yellow), const float3& color3 = GetColorRGBFloat(EColor::Red))
	{
		float interpolationParameter = 0;
		if (fabs(valueMax - valueMin) >= 0.01)
		{
			interpolationParameter = MapTo01Range(curValue, valueMin, valueMax);
		}

		///Get color based on cur magnitude
		float3 CurColor = TQuadraticInterpolateBezier(interpolationParameter, color1, color2, color3);
		//static std::vector<float3> outputSet{ color1, color2, color3 };
		//float3 CurColor = TAitkenInterpolate(outputSet, std::vector<float>{0.f,0.5f,1.f}, interpolationParameter);
		return CurColor;
	};

	template<typename T>
	void CopyFromInvertedCPUBufferToGPUTexture(RCommandList& CmdList, RDynamicVector<T>& sourceBuffer, RDynamicVector<T>& bufferIntermediate /*we cant kill this resource*/, RTexture& destTexture)
	{
		bufferIntermediate.resize(sourceBuffer.size());

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
				uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

				bufferIntermediate[destIndex] = sourceBuffer[sourceIndex];
			}
		}

		CmdList.UploadToTexture(destTexture, bufferIntermediate.data(), destTexture.Size.x, destTexture.Size.y, destTexture.BytesPerPixel);
		//RCommandList::UploadToTextureImmediate(destTexture, bufferIntermediate.data(), destTexture.Size.x, destTexture.Size.y, destTexture.BytesPerPixel);
	}


	void VisualizeVelocityFieldTexture(RTexture& textureToVisualize, RTexture* pNormalsTexture = nullptr)
	{
		float half = (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) * 0.5f;
		float center = CGConstantsCB.VelocityFieldBoundsStartOffset + half;

		if (bVelocityField2DVerticalSliceSimulation)
		{
			//XY
			RGeometry::GInstance.AddPlaneTextured(float3{ center, center, CGConstantsCB.VelocityFieldBoundsStartOffset }, float2{ half, half }, textureToVisualize, false, pNormalsTexture);
		}
		else
		{
			//XZ
			RGeometry::GInstance.AddPlaneTextured(float3{ center, CGConstantsCB.DefaultParticlePos.y, center }, float2{ half, half }, textureToVisualize, true, pNormalsTexture);
		}

		
	}

	void VisualizeDivergence(RCommandListCompute& CmdList)
	{
		if (bProcessVectorFieldOnCPU == 1)
		{
			static RDynamicVector<uint32> divergenceBufferColoredInverted;
			divergenceBufferColoredInverted.resize(DivergenceBufferCPU.size());

			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };
					uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
					uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

					float curDivergence = (DivergenceBufferCPU[sourceIndex]);

					float3 CurColor;
					if (curDivergence >= 0)
					{
						CurColor = float3{ curDivergence, 0, curDivergence };
					}
					else
					{
						float curDivAbs = std::fabs(curDivergence);
						CurColor = float3{ 0, curDivAbs, curDivAbs };
					}

					divergenceBufferColoredInverted[destIndex] = RColor(CurColor.x, CurColor.y, CurColor.z).ConvertToR11G11B10FAndReturn();
				}
			}

			RCommandList::UploadToTextureImmediate(DivergenceTextureGPU, divergenceBufferColoredInverted.data(), DivergenceTextureGPU.Size.x, DivergenceTextureGPU.Size.y, DivergenceTextureGPU.BytesPerPixel);

			VisualizeVelocityFieldTexture(DivergenceTextureGPU);
		}
		else
		{
			auto& visualizerTexture = VisualizerConvertTextureWithNegativeValuesToColor(CmdList, DivergenceTextureGPU);
			VisualizeVelocityFieldTexture(visualizerTexture);
		}
		

	}

	void VisualizeCurl(RCommandListCompute& CmdList)
	{
		auto& visualizerTexture = VisualizerConvertTextureWithNegativeValuesToColor(CmdList, CurlTextureGPU);
		VisualizeVelocityFieldTexture(visualizerTexture);
	}

	void VisualizePressure(RCommandListCompute& CmdList)
	{
		if (bProcessVectorFieldOnCPU == 1)
		{
			static RDynamicVector<uint32> pressureBufferInverted;
			pressureBufferInverted.resize(PressureBufferCPU.size());

			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };
					uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
					uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

					float curpressure = (PressureBufferCPU[sourceIndex]);

					float3 CurColor;
					if (curpressure >= 0)
					{
						CurColor = float3{ curpressure, 0, 0 };
					}
					else
					{
						float curPrAbs = std::fabs(curpressure);
						CurColor = float3{ 0, 0, curPrAbs };
					}

					pressureBufferInverted[destIndex] = RColor(CurColor.x, CurColor.y, CurColor.z).ConvertToR11G11B10FAndReturn();
				}
			}

			RCommandList::UploadToTextureImmediate(PressureTextureColoredGPU, pressureBufferInverted.data(), PressureTextureColoredGPU.Size.x, PressureTextureColoredGPU.Size.y, PressureTextureColoredGPU.BytesPerPixel);

			VisualizeVelocityFieldTexture(PressureTextureColoredGPU);
		}
		else
		{
			auto& visualizerTexture = VisualizerConvertTextureWithNegativeValuesToColor(CmdList, PressureTextureGPU);
			VisualizeVelocityFieldTexture(visualizerTexture);
		}
		

	}

	void VisualizeVelocityFieldMagnitudes()
	{
		static RDynamicVector<uint32> magnitudesBufferColoredInverted;
		magnitudesBufferColoredInverted.resize(DivergenceBufferCPU.size());

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
				uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

				float2 curVelocity = LoadFromVelocityField(coord);
				const float curMagnitude = VectorGetLength(curVelocity);
				float3 CurColor = MapTo3ColorRange(curMagnitude, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax);

				magnitudesBufferColoredInverted[destIndex] = RColor(CurColor.x, CurColor.y, CurColor.z).ConvertToR11G11B10FAndReturn();
			}
		}

		RCommandList::UploadToTextureImmediate(VelocityFieldMagnitudesTextureGPU, magnitudesBufferColoredInverted.data(), VelocityFieldMagnitudesTextureGPU.Size.x, VelocityFieldMagnitudesTextureGPU.Size.y, VelocityFieldMagnitudesTextureGPU.BytesPerPixel);

		VisualizeVelocityFieldTexture(VelocityFieldMagnitudesTextureGPU);
		
	}

	void VisualizeDensity()
	{
		//Upload Density field to texture
#if 0
		{
			//Compute min max density, convert interpolate value in each cell into uint32 RGB color
			static RDynamicVector<uint32> densityBufferColoredCPU;
			densityBufferColoredCPU.resize(DensityBufferCPU.size());
			float densityMax = VectorGetLength(DensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 0,0 })]);
			float densityMin = VelocityFieldCurMagnitudeMax;
			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };

					float curDensity = (DensityBufferCPU[GetFlatIndexToVelocityArray2D(coord)]);
					densityMax = std::max(densityMax, curDensity);
					densityMin = std::min(densityMin, curDensity);
				}
			}

			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };
					uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });

					//Invert density buffer vertically
					uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

					float curDensity = (DensityBufferCPU[sourceIndex]);
					
					//float3 CurColor = MapTo3ColorRange(curDensity, densityMin, densityMax);
					float3 CurColor = MapTo3ColorRange(curDensity, 0.f, 1.f, GetColorRGBFloat(EColor::BlueEgyptian), GetColorRGBFloat(EColor::YellowCanary), GetColorRGBFloat(EColor::Red));
					//float3 CurColor{ curDensity, curDensity, curDensity };

					densityBufferColoredCPU[destIndex] = RColor(CurColor.x, CurColor.y, CurColor.z).ConvertToR11G11B10FAndReturn();

				}
			}

			RCommandList::UploadToTextureImmediate(DensityTextureColoredGPU, densityBufferColoredCPU.data(), DensityTextureColoredGPU.Size.x, DensityTextureColoredGPU.Size.y, DensityTextureColoredGPU.BytesPerPixel);

		}
#endif

		if (CGConstantsCB.bUseParticleBasedAdvection)
		{
			VisualizeVelocityFieldTexture(AdvectedParticleDensityTextureGPU);
		}
		else if (bEnableVorticityConfinement && bApplyVorticityMaskToDensityVisualizer)
		{
			VisualizeVelocityFieldTexture(DensityTextureCopyGPU);
		}
		else
		{
			VisualizeVelocityFieldTexture(DensityTextureGPU);
		}

	}

	void ApplyConstExternalForceToVectorField()
	{
		auto& states = VelocityFieldStateBufferCPU;

		for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
		{
			for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
			{
				uint2 coord{ x,y };

				if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
					continue;

				float2 vCur = LoadFromVelocityField(coord);
				StoreInVelocityField(coord, vCur + CGConstantsCB.ExternalForce.xy() * CGConstantsCB.DeltaTime);
			}
		}
	}

	void ComputeDivergenceCPU(bool bMethod2 = false)
	{
		auto& states = VelocityFieldStateBufferCPU;

		const float samplingPeriodLength = CGConstantsCB.VelocityFieldGridCellLength;

		DivergenceBufferCPU.resize(CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize);

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				if (states[i] == 0)
				{
					DivergenceBufferCPU[i] = 0.f;
					continue;
				}


				float2 uvCur = LoadFromVelocityField(coord);
				float uRight = LoadFromVelocityField_U(coord + uint2{ 1, 0 });
				float vUp = LoadFromVelocityField_V(coord + uint2{ 0, 1 });

				float div;

				if (bMethod2)
				{
					div = -0.5f * samplingPeriodLength * (uRight - uvCur.x + vUp - uvCur.y); //divergence computation from Jos Stam paper
				}
				else
				{
					div = ((uRight - uvCur.x) + (vUp - uvCur.y)) / (2.f * samplingPeriodLength);//simple central diff
				}

				DivergenceBufferCPU[i] = div;

			}
		}
	}

	bool bAllowVelocityWritesToSolidCellsFacesDuringProjection = false;//experimental

	void ProjectVectorFieldCPU(bool bFirstPass)
	{

		auto& states = VelocityFieldStateBufferCPU;

		if (bFirstPass)
		{
			PressureBufferCPU.resize(CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize);
			std::fill(PressureBufferCPU.begin(), PressureBufferCPU.end(), 0);
		}

		for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize-1; y++)
		{
			for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize-1; x++)
			{
				uint2 coord{ x,y };

				if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
				{
					continue;
				}

				if (CGConstantsCB.bUseParticleBasedAdvection && states[GetFlatIndexToVelocityArray2D(coord)] != VELOCITY_FIELD_CELL_STATE_DENSITY)
				{
					continue;
				}

				uint sx0 = states[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })] != VELOCITY_FIELD_CELL_STATE_SOLID;
				uint sx1 = states[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })] != VELOCITY_FIELD_CELL_STATE_SOLID;
				uint sy0 = states[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })] != VELOCITY_FIELD_CELL_STATE_SOLID;
				uint sy1 = states[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })] != VELOCITY_FIELD_CELL_STATE_SOLID;

				if (bAllowVelocityWritesToSolidCellsFacesDuringProjection)
				{
					sx0 = 1;
					sx1 = 1;
					sy0 = 1;
					sy1 = 1;
				}

				auto s = sx0 + sx1 + sy0 + sy1;

				if (s == 0) continue;

				float2 uvCur = LoadFromVelocityField(coord);
				float uRight = LoadFromVelocityField_U(coord + uint2{ 1, 0 });
				float vUp = LoadFromVelocityField_V(coord + uint2{ 0, 1 });

				float div = uRight - uvCur.x + vUp - uvCur.y;

				div = div / float(s);
				
				PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord)] -= div * CGConstantsCB.VelocityFieldGridCellLength * CGConstantsCB.DeltaTime;

				#if 1 //OVERRELAXATION
				{
					div = div * 1.9f;
				}
				#endif//OVERRELAXATION

				if (CGConstantsCB.bUseParticleBasedAdvection && AverageParticleDensity > 0.0 && bParticleAdvectionCompensateDrift) 
				{
					//reduce divergence in dense regions
					float pressure = AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(coord)] - AverageParticleDensity;
					if (pressure > 0.0)
					{
						div = div - CGConstantsCB.AdvectedParticleDriftCompensationScale * pressure;
					}
				}

				uvCur.x += sx0 * div;
				uvCur.y += sy0 * div;
				StoreInVelocityField(coord, uvCur);

				uRight -= sx1 * div;
				StoreInVelocityField_U(coord + uint2{ 1, 0 }, uRight);

				vUp -= sy1 * div;
				StoreInVelocityField_V(coord + uint2{ 0, 1 }, vUp);
			}
		}

		if (bAllowVelocityWritesToSolidCellsFacesDuringProjection)
		{
			ClearVelocityFromSolidCells();
		}

	}

	RDynamicVector<float> PressureBufferCPU;

	//Projection for collocated grid
	void ProjectPressureCPU()
	{
		const float samplingPeriodLength = CGConstantsCB.VelocityFieldGridCellLength;

		//compute divergence at each point
		ComputeDivergenceCPU(true);

		ExtrapolateValueToBorders(DivergenceBufferCPU);

		//compute pressure at each point
		PressureBufferCPU.resize(CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize);
		std::fill(PressureBufferCPU.begin(), PressureBufferCPU.end(), 0.f);

		for (int k = 0; k < VectorFieldProjectionNumIterations; k++)
		{
			for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
			{
				for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
				{
					uint2 coord{ x,y };

					 PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord)] =
							(DivergenceBufferCPU[GetFlatIndexToVelocityArray2D(coord)]
							+ PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })]
							+ PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })]
							+ PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })]
							+ PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })]
							) / 4.f;
				}
			}

			ExtrapolateValueToBorders(PressureBufferCPU);
			
		}

		//remove pressure gradient (pointing towards highest pressure area) from velocity vector
		for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
		{
			for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
			{
				uint2 coord{ x,y };

				uint sx0, sy0;

				if (bAllowVelocityWritesToSolidCellsFacesDuringProjection)
				{
					sx0 = 1;
					sy0 = 1;
				}
				else
				{
					if (VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord)] != VELOCITY_FIELD_CELL_STATE_DENSITY)
					{
						continue;
					}

					sx0 = VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })];
					sy0 = VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })];
				}

				const float overrelaxationScale = 1.5f;

				auto curVelocity = LoadFromVelocityField(coord);
				curVelocity.x -= sx0 * overrelaxationScale * (PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord)] - PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 1,0 })]) / samplingPeriodLength;
				curVelocity.y -= sy0 * overrelaxationScale * (PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord)] - PressureBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 0,1 })]) / samplingPeriodLength;
				StoreInVelocityField(coord, curVelocity);
			}
		}

		if (bAllowVelocityWritesToSolidCellsFacesDuringProjection)
		{
			ClearVelocityFromSolidCells();
		}

	}


	/**********************
		   Boundaries
	**********************/

	bool4 bEnabledBorders{ 1,1,1,1 };//left, up, right, down 

	void VelocityFieldSetBoundaries()
	{
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				VelocityFieldStateBufferCPU[i] = 1;

				const bool borderLeft = bEnabledBorders[0] && coord.x == 0;
				const bool borderUp = bEnabledBorders[1] && coord.y == (CGConstantsCB.VelocityFieldGridSize - 1);
				const bool borderRight = bEnabledBorders[2] && coord.x == (CGConstantsCB.VelocityFieldGridSize - 1);
				const bool borderDown = bEnabledBorders[3] && coord.y == 0;

				bool isBorder = 0;
				isBorder |= borderLeft;
				isBorder |= borderRight;
				isBorder |= borderDown;
				isBorder |= borderUp;

				if (isBorder)
				{
					VelocityFieldStateBufferCPU[i] = 0;
					//DensityBufferCPU[i] = 0;
				}
			}
		}
	}

	uint IndexToObstacleControlPoint;

	void VelocityFieldSetObstacle()
	{
		const float h = CGConstantsCB.VelocityFieldGridCellLength;

		if (GbFirstBoot)
		{
			float center = CGConstantsCB.VelocityFieldBoundsStartOffset + CGConstantsCB.VelocityFieldBoundsLength / 2.f;
			VelocityFieldObstaclePos = float3{ center, center, 0 };
			VelocityFieldObstaclePosPrev = VelocityFieldObstaclePos.xy();

			IndexToObstacleControlPoint = RSpatialControlSet::Get().ControlPoints.size();
			auto& controlPoint = RSpatialControlSet::Get().AddPoint(&VelocityFieldObstaclePos, VelocityFieldObstacleRadius * 0.5f);

		}

		VelocityFieldObstaclePos.z = 0.f;

		bool bUpdateVelocity = false;

		float2 obstacleVelocity{ 0,0 };

		//update vel only when in drag mode
		if (RSpatialControlSet::Get().bInPositionDragState && (&RSpatialControlSet::Get().ControlPoints[IndexToObstacleControlPoint] == RSpatialControlSet::Get().pSelectedPoint))
		{
			obstacleVelocity = (VelocityFieldObstaclePos.xy() - VelocityFieldObstaclePosPrev) / CGConstantsCB.DeltaTime;
			VelocityFieldObstaclePosPrev = VelocityFieldObstaclePos.xy();
			bUpdateVelocity = true;
		}


		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				float2 curPos = GetVelocityCellCenterWorldSpaceFromCoord(coord);

				if (RectangleSphereIntersection(curPos, float2::MakeFloat2(h * 0.5f), VelocityFieldObstaclePos.xy(), VelocityFieldObstacleRadius))
				{
					VelocityFieldStateBufferCPU[i] = 0;

					if (bUpdateVelocity && bObstacleGeneratesVelocity)
					{
						//assign obstacle velocity to all vecs in staggered grid cell
						StoreInVelocityField(i, obstacleVelocity);

						auto iRight = GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 });
						if (iRight < VelocityField_U_ComponentBufferCPU.size())
						{
							StoreInVelocityField_U(iRight, obstacleVelocity.x);
						}
						auto iiUp = GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 });
						if (iiUp < VelocityField_V_ComponentBufferCPU.size())
						{
							StoreInVelocityField_V(iiUp, obstacleVelocity.y);
						}

					}
					else
					{
						//VelocityFieldBufferCPU[i] = float2{ 0,0 };
					}


					//Apply density to obstacle
					if (bObstacleAsDensitySource)
					{
						DensityBufferCPU[i] = CGConstantsCB.InitialDensityColor;
					}
					else
					{
#if USE_3COLOR_DENSITY == 1
						DensityBufferCPU[i] = float3{ 0,0,0 };
#else
						DensityBufferCPU[i] = 0;
#endif//USE_3COLOR_DENSITY
					}
				}
			}
		}

	}

	void ApplyConstForceToVelocityFieldBorders()
	{
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				const bool borderLeft = coord.x == 1;
				const bool borderRight = coord.x == (CGConstantsCB.VelocityFieldGridSize - 2);
				const bool borderDown = coord.y == 1;
				const bool borderUp = coord.y == (CGConstantsCB.VelocityFieldGridSize - 2);

				const bool bRequiredBorder = borderDown;

				if (bRequiredBorder)
				{
					StoreInVelocityField(i, VelocityFieldForceAtBorders);
				}
			}
		}
	}

	//Use when solid borders are disabled
	bool bCopyVelocitiesToBorders = false;
	bool bCopyReflectedVelocitiesToBorders = false;

	//Copy values to borders
	void ExtrapolateVelocitiesToBordersCPU()
	{
		uint n = CGConstantsCB.VelocityFieldGridSize;

#if 1 //only copy separate components
		for (uint x = 1; x < n; x++)
		{
			if (bCopyVelocitiesToBorders)
			{
				VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, 0 })] = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, 1 })];
				VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, n - 1 })] = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, n - 2 })];
			}
			
			if (bCopyReflectedVelocitiesToBorders)
			{
				VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, 0 })] = -VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, 1 })];
				VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, n - 1 })] = -VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ x, n - 2 })];
			}
		}
		for (uint y = 1; y < n; y++)
		{
			if (bCopyVelocitiesToBorders)
			{
				VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 0, y })] = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 1, y })];
				VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ n - 1, y })] = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ n - 2, y })];
			}
			
			if (bCopyReflectedVelocitiesToBorders)
			{
				VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 0, y })] = -VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 1, y })];
				VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ n - 1, y })] = -VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ n - 2, y })];
			}
		}

		StoreInVelocityField(uint2{ 0, 0 }, (LoadFromVelocityField(uint2{ 1, 0 }) + LoadFromVelocityField(uint2{ 0, 1 })) * 0.5f);
		StoreInVelocityField(uint2{ 0, n - 1 }, (LoadFromVelocityField(uint2{ 1, n - 1 }) + LoadFromVelocityField(uint2{ 0, n - 2 })) * 0.5f);
		StoreInVelocityField(uint2{ n - 1, 0 }, (LoadFromVelocityField(uint2{ n - 2, 0 }) + LoadFromVelocityField(uint2{ n - 1, 1 })) * 0.5f);
		StoreInVelocityField(uint2{ n - 1, n - 1 }, (LoadFromVelocityField(uint2{ n - 2, n - 1 }) + LoadFromVelocityField(uint2{ n - 1, n - 2 })) * 0.5f);
#else
		for (uint x = 0; x < n; x++)
		{
			velocities[GetFlatIndexToVelocityArray2D(uint2{ x, 0 })] = velocities[GetFlatIndexToVelocityArray2D(uint2{ x, 1 })];
			velocities[GetFlatIndexToVelocityArray2D(uint2{ x, n - 1 })] = velocities[GetFlatIndexToVelocityArray2D(uint2{ x, n - 2 })];
		}
		for (uint y = 0; y < n; y++)
		{
			velocities[GetFlatIndexToVelocityArray2D(uint2{ 0, y })] = velocities[GetFlatIndexToVelocityArray2D(uint2{ 1, y })];
			velocities[GetFlatIndexToVelocityArray2D(uint2{ n - 1, y })] = velocities[GetFlatIndexToVelocityArray2D(uint2{ n - 2, y })];
		}

		ClearVelocityFromSolidCells();

#endif
	}

	template <typename Buffer>
	void ExtrapolateValueToBorders(Buffer& buffer)
	{
		uint n = CGConstantsCB.VelocityFieldGridSize;
		for (uint i = 1; i < n-1; i++)
		{
			buffer[GetFlatIndexToVelocityArray2D(uint2{ i, 0 })] = buffer[GetFlatIndexToVelocityArray2D(uint2{ i, 1 })];
			buffer[GetFlatIndexToVelocityArray2D(uint2{ i, n - 1 })] = buffer[GetFlatIndexToVelocityArray2D(uint2{ i, n - 2 })];

			buffer[GetFlatIndexToVelocityArray2D(uint2{ 0, i })] = buffer[GetFlatIndexToVelocityArray2D(uint2{ 1, i })];
			buffer[GetFlatIndexToVelocityArray2D(uint2{ n - 1, i })] = buffer[GetFlatIndexToVelocityArray2D(uint2{ n - 2, i })];
		}

		buffer[GetFlatIndexToVelocityArray2D(uint2{ 0, 0 })] = (buffer[GetFlatIndexToVelocityArray2D(uint2{ 1, 0 })] + buffer[GetFlatIndexToVelocityArray2D(uint2{ 0, 1 })]) * 0.5f;
		buffer[GetFlatIndexToVelocityArray2D(uint2{ 0, n-1 })] = (buffer[GetFlatIndexToVelocityArray2D(uint2{ 1, n-1 })] + buffer[GetFlatIndexToVelocityArray2D(uint2{ 0, n-2 })]) * 0.5f;
		buffer[GetFlatIndexToVelocityArray2D(uint2{ n-1, 0})] = (buffer[GetFlatIndexToVelocityArray2D(uint2{ n-2, 0 })] + buffer[GetFlatIndexToVelocityArray2D(uint2{ n-1, 1})]) * 0.5f;
		buffer[GetFlatIndexToVelocityArray2D(uint2{ n-1, n-1})] = (buffer[GetFlatIndexToVelocityArray2D(uint2{ n-2, n-1 })] + buffer[GetFlatIndexToVelocityArray2D(uint2{ n-1, n-2})]) * 0.5f;

	}


	bool bClearVelocityFromSolidCell = true;

	void ClearVelocityFromSolidCells()
	{
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				float2 curVelocity = LoadFromVelocityField(i);

				//make sure cur velocity doesn't flow into solid

				//solid to the left
				if (x > 0)
				{
					if (VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })] == 0)
					{
						if (curVelocity.x < 0)
						{
							StoreInVelocityField_U(i, 0);
						}
					}
				}
				
				//solid down
				if (y > 0)
				{
					if (VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })] == 0)
					{
						if (curVelocity.y < 0)
						{
							StoreInVelocityField_V(i, 0);
						}
					}
				}
				
				if (VelocityFieldStateBufferCPU[i] == 0)
				{
					if (curVelocity.x > 0)
					{
						StoreInVelocityField_U(i, 0);
					}
					if (curVelocity.y > 0)
					{
						StoreInVelocityField_V(i, 0);
					}
				}
			}
		}
	}








	float2 GetVelocityFieldCellStartPosWorldSpaceFromCoord(uint2 coord)
	{
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		const float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };
		float2 curPos{ boundsStartPos.x + coord.x * h, boundsStartPos.y + coord.y * h };
		return curPos;
	}

	float2 GetVelocityCellCenterWorldSpaceFromCoord(uint2 coord)
	{
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
		curPos = curPos + float2{ h * 0.5f, h * 0.5f };
		return curPos;
	}

	float2 GetVelocityXComponentWorldPosFromCoord(uint2 coord)
	{
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
		curPos.y = curPos.y + h * 0.5f;
		return curPos;
	}
	float2 GetVelocityYComponentWorldPosFromCoord(uint2 coord)
	{
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
		curPos.x = curPos.x + h * 0.5f;
		return curPos;
	}

	void AdvectVelocityFieldCPU()
	{
		//for horizontal vector u
		//compute cur full velocity vector at world space position where u is stored
			//average neighbour v components to get cur u velocity .y
		//subtract this cur velocity from u world space to get u prev world space pos
			//use weighted average to get velocity for prev world space
		
		auto& states = VelocityFieldStateBufferCPU;

		static RDynamicVector<float> velocity_U_Copy;
		static RDynamicVector<float> velocity_V_Copy;
		velocity_U_Copy = VelocityField_U_ComponentBufferCPU;
		velocity_V_Copy = VelocityField_V_ComponentBufferCPU;

		float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };


		for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
		{
			for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
			{
				uint2 coord{ x,y };
				float2 curVelocity = LoadFromVelocityField(coord);

				if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
					continue;

				//u
				if (states[GetFlatIndexToVelocityArray2D(coord - uint2(1, 0))] != 0 && y < (CGConstantsCB.VelocityFieldGridSize - 1))
				{
					auto u = curVelocity.x;

					//get average v
					float vCur = curVelocity.y;
					float vUp = LoadFromVelocityField_V(coord + uint2(0, 1));
					float vLeft = LoadFromVelocityField_V(coord - uint2(1, 0));
					float vLeftUp = LoadFromVelocityField_V(coord - uint2(1, 0) + uint2(0, 1));
					float vAvrg = (vCur + vUp + vLeft + vLeftUp) * 0.25f;

					float2 velocityFull{ u, vAvrg };

					//get u world space pos
					float2 curPos = GetVelocityXComponentWorldPosFromCoord(coord);

					float2 prevPos = curPos - velocityFull * CGConstantsCB.VelocityAdvectionScaleFactor * CGConstantsCB.DeltaTime;

					float2 velocityNew = SampleVelocityField(prevPos);
					curVelocity.x = velocityNew.x;

				}

				//v
				if (states[GetFlatIndexToVelocityArray2D(coord - uint2(0, 1))] != 0 && x < (CGConstantsCB.VelocityFieldGridSize - 1))
				{
					auto& v = curVelocity.y;

					//get average u
					float uCur = curVelocity.x;
					float uDown =  LoadFromVelocityField_U(coord - uint2(0, 1));
					float uRight = LoadFromVelocityField_U(coord + uint2(1, 0));
					float uRightDown = LoadFromVelocityField_U(coord + uint2(1, 0) - uint2(0, 1));
					float uAvrg = (uCur + uDown + uRight + uRightDown) * 0.25f;

					float2 velocityFull{ uAvrg, v};

					//get u world space pos
					float2 curPos = GetVelocityYComponentWorldPosFromCoord(coord);

					float2 prevPos = curPos - velocityFull * CGConstantsCB.VelocityAdvectionScaleFactor * CGConstantsCB.DeltaTime;

					float2 velocityNew = SampleVelocityField(prevPos);
					curVelocity.y = velocityNew.y;

				}

				velocity_U_Copy[GetFlatIndexToVelocityArray2D(coord)] = curVelocity.x;
				velocity_V_Copy[GetFlatIndexToVelocityArray2D(coord)] = curVelocity.y;

			}
		}

		VelocityField_U_ComponentBufferCPU = velocity_U_Copy;
		VelocityField_V_ComponentBufferCPU = velocity_V_Copy;


	}

	bool bDiffuseVelocity = false;

	void DiffuseVelocityCPU()
	{
		auto& states = VelocityFieldStateBufferCPU;

		static RDynamicVector<float> velocityFieldNew_U;
		static RDynamicVector<float> velocityFieldNew_V;
		velocityFieldNew_U.resize(VelocityField_U_ComponentBufferCPU.size());
		std::fill(velocityFieldNew_U.begin(), velocityFieldNew_U.end(), 0.f);
		velocityFieldNew_V.resize(VelocityField_V_ComponentBufferCPU.size());
		std::fill(velocityFieldNew_V.begin(), velocityFieldNew_V.end(), 0.f);

		float a = CGConstantsCB.DeltaTime * (CGConstantsCB.DiffusionFactor * 10.f);

		for (int i = 0; i < DiffusionNumIterations; i++)
		{
			for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
			{
				for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
				{
					uint2 coord{ x,y };

					float2 curVelocity = LoadFromVelocityField(coord);

					if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
					{
						velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord)] = curVelocity.x;
						velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord)] = curVelocity.x;
						continue;
					}

					velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord)] =
						(curVelocity.x +
							(
								velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })]
								+ velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })]
								+ velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })]
								+ velocityFieldNew_U[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })]
								) * a
							) / (1 + 4 * a);

					velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord)] =
						(curVelocity.y +
							(
								velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })]
								+ velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })]
								+ velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })]
								+ velocityFieldNew_V[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })]
								) * a
							) / (1 + 4 * a);

				}
			}
		}

		VelocityField_U_ComponentBufferCPU = velocityFieldNew_U;
		VelocityField_V_ComponentBufferCPU = velocityFieldNew_V;
	}

	uint DiffusionNumIterations = 40;

	bool bDiffuseDensity = false;

	void DiffuseDensityCPU()
	{
		auto& states = VelocityFieldStateBufferCPU;
		auto& densityFieldCur = DensityBufferCPU;

		static RDynamicVector<DENSITY_TYPE> densityFieldNew;
		densityFieldNew.resize(densityFieldCur.size());

#if USE_3COLOR_DENSITY == 1
		std::fill(densityFieldNew.begin(), densityFieldNew.end(), float3{ 0,0,0 });
#else
		std::fill(densityFieldNew.begin(), densityFieldNew.end(), 0.f);
#endif//USE_3COLOR_DENSITY

		ExtrapolateValueToBorders(densityFieldCur);

		//float a = CGConstantsCB.DeltaTime * (DensityDiffuseDiffScale * 0.001f) * (CGConstantsCB.VelocityFieldGridSize - 1 ) * (CGConstantsCB.VelocityFieldGridSize - 1);
		float a = CGConstantsCB.DeltaTime * (CGConstantsCB.DiffusionFactor * 10.f);

		for (int i = 0; i < DiffusionNumIterations; i++)
		{
			for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
			{
				for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
				{
					uint2 coord{ x,y };
					/*if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
					{
						densityFieldNew[GetFlatIndexToVelocityArray2D(coord)] = densityFieldCur[GetFlatIndexToVelocityArray2D(coord)];
						continue;
					}*/
						
					densityFieldNew[GetFlatIndexToVelocityArray2D(coord)] =
						(densityFieldCur[GetFlatIndexToVelocityArray2D(coord)] +
								(
								densityFieldNew[GetFlatIndexToVelocityArray2D(coord - uint2{ 1, 0 })]
								+ densityFieldNew[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })]
								+ densityFieldNew[GetFlatIndexToVelocityArray2D(coord - uint2{ 0, 1 })]
								+ densityFieldNew[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })]
								) * a
						) / (1 + 4 * a);
				}
			}

			ExtrapolateValueToBorders(densityFieldNew);

		}

		densityFieldCur = densityFieldNew;

	}

	float DensityAdvectionDiffusionFactor = 1.0f;//damping
	float DensityAdvectionScaleFactor = 1.f;

	void AdvectDensityFieldCPU()
	{
		auto& states = VelocityFieldStateBufferCPU;

		static RDynamicVector<DENSITY_TYPE> densityFieldCopy;
		densityFieldCopy = DensityBufferCPU;

		float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };

		for (uint y = 1; y < CGConstantsCB.VelocityFieldGridSize - 1; y++)
		{
			for (uint x = 1; x < CGConstantsCB.VelocityFieldGridSize - 1; x++)
			{
				uint2 coord{ x,y };
				
				if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
					continue;

				//get cur velocity of cell center 
				float2 curVelocity;
				curVelocity.x = (VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] + VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2{ 1, 0 })]) * 0.5f;
				curVelocity.y = (VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] + VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2{ 0, 1 })]) * 0.5f;

				float2 curPosWorld = GetVelocityCellCenterWorldSpaceFromCoord(coord);

				float2 newPosWorld = curPosWorld - curVelocity * DensityAdvectionScaleFactor * CGConstantsCB.DeltaTime;

				//sample density at new location
				DENSITY_TYPE newdensity = SampleDensity(newPosWorld);

				densityFieldCopy[GetFlatIndexToVelocityArray2D(coord)] = newdensity * DensityAdvectionDiffusionFactor;

			}
		}

		DensityBufferCPU = densityFieldCopy;

	}

	void GetVelocityFieldCoordFromWorldPos(float2 worldPos, uint2& outCurCoord, uint2& outCurCoordNext, float2& outFrac)
	{
		const float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };
		const uint n = CGConstantsCB.VelocityFieldGridSize;
		const float h = CGConstantsCB.VelocityFieldGridCellLength;

		const float2 minLimitPos = boundsStartPos + float2{ h, h };
		const float2 maxLimitPos = boundsStartPos + float2{ h * n, h * n };
		worldPos.x = std::clamp(worldPos.x, minLimitPos.x, maxLimitPos.x);
		worldPos.y = std::clamp(worldPos.y, minLimitPos.y, maxLimitPos.y);

		float2 curCoordf = (worldPos - boundsStartPos) / float(CGConstantsCB.VelocityFieldGridCellLength);
		outCurCoord = uint2{ uint(std::floor(curCoordf.x)), uint(std::floor(curCoordf.y)) };

		outCurCoord.x = std::min(outCurCoord.x, n - 1);
		outCurCoord.y = std::min(outCurCoord.y, n - 1);

		outCurCoordNext = uint2{ std::min(outCurCoord.x + 1, n - 1) , std::min(outCurCoord.y + 1, n - 1) };

		outFrac = float2{ curCoordf.x - outCurCoord.x, curCoordf.y - outCurCoord.y };
	};


	DENSITY_TYPE SampleDensity(float2 worldPos)
	{
		const float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		const uint n = CGConstantsCB.VelocityFieldGridSize;

		//density values located at the centers of the grid

		float2 worldPosMinusHalf = worldPos - float2{ h * 0.5f, h * 0.5f };

		uint2 curCoord;
		uint2 curCoordNext;
		float2 frac;
		GetVelocityFieldCoordFromWorldPos(worldPosMinusHalf, curCoord, curCoordNext, frac);

		//sample
		DENSITY_TYPE dCur = DensityBufferCPU[GetFlatIndexToVelocityArray2D(curCoord)];
		DENSITY_TYPE dUp = DensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
		DENSITY_TYPE dRight = DensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
		DENSITY_TYPE dRightUp = DensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

		//interpolate
		//horizontal with frac.x
		DENSITY_TYPE HorizontalUp = dUp * (1.f - frac.x) + dRightUp * frac.x;
		DENSITY_TYPE HorizontalDown = dCur * (1.f - frac.x) + dRight * frac.x;

		//vertical
		DENSITY_TYPE res = HorizontalDown * (1.f - frac.y) + HorizontalUp * frac.y;

		return res;
	}

	float2 LoadFromVelocityField(uint2 coord)
	{
		auto& u = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
		auto& v = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
		return float2{ u,v };
	}

	float LoadFromVelocityField_U(uint2 coord)
	{
		return VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
	}
	float LoadFromVelocityField_V(uint2 coord)
	{
		return VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
	}

	float2 LoadFromVelocityField(uint flatIndex)
	{
		auto& u = VelocityField_U_ComponentBufferCPU[flatIndex];
		auto& v = VelocityField_V_ComponentBufferCPU[flatIndex];
		return float2{ u,v };
	}

	void StoreInVelocityField(uint2 coord, float2 value)
	{
		VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] = value.x;
		VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] = value.y;
	}

	void StoreInVelocityField_U(uint2 coord, float value)
	{
		VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] = value;
	}
	void StoreInVelocityField_U(uint flatIndex, float value)
	{
		VelocityField_U_ComponentBufferCPU[flatIndex] = value;
	}
	void StoreInVelocityField_V(uint2 coord, float value)
	{
		VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord)] = value;
	}
	void StoreInVelocityField_V(uint flatIndex, float value)
	{
		VelocityField_V_ComponentBufferCPU[flatIndex] = value;
	}
	void StoreInVelocityField(uint flatIndex, float2 value)
	{
		VelocityField_U_ComponentBufferCPU[flatIndex] = value.x;
		VelocityField_V_ComponentBufferCPU[flatIndex] = value.y;
	}

	float2 SampleVelocityField(float2 worldPos)
	{
		const float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };
		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		const uint n = CGConstantsCB.VelocityFieldGridSize;

		float2 velocityFinal;

		//get u
		{
			//offset half .y to convert into u space
			float2 worldPosMinusHalfX = worldPos - float2{ 0, h * 0.5f };

			uint2 curCoord;
			uint2 curCoordNext;
			float2 frac;
			GetVelocityFieldCoordFromWorldPos(worldPosMinusHalfX, curCoord, curCoordNext, frac);

			//sample
			float uCur = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(curCoord)];
			float uUp = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
			float uRight = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
			float uRightUp = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

			//interpolate
			//horizontal with frac.x
			float velocityHorizontalUp = uUp * (1.f - frac.x) + uRightUp * frac.x;
			float velocityHorizontalDown = uCur * (1.f - frac.x) + uRight * frac.x;

			//vertical
			velocityFinal.x = velocityHorizontalDown * (1.f - frac.y) + velocityHorizontalUp * frac.y;
			//velocityFinal.y = vCur * (1.f - frac.x) * (1.f - frac.y) + vRight * frac.x * (1.f - frac.y) + vUp * (1.f - frac.x) * frac.y + vRightUp * frac.x * frac.y;
		}

		//get v
		{
			//offset half .x to convert into v space
			float2 worldPosMinusHalfX = worldPos - float2{ h * 0.5f, 0};

			uint2 curCoord;
			uint2 curCoordNext;
			float2 frac;
			GetVelocityFieldCoordFromWorldPos(worldPosMinusHalfX, curCoord, curCoordNext, frac);

			//sample
			float vCur = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(curCoord)];
			float vUp = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
			float vRight = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
			float vRightUp = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

			//interpolate
			//horizontal with frac.x
			float velocityHorizontalUp = vUp * (1.f - frac.x) + vRightUp * frac.x;
			float velocityHorizontalDown = vCur * (1.f - frac.x) + vRight * frac.x;

			//vertical
			velocityFinal.y = velocityHorizontalDown * (1.f - frac.y) + velocityHorizontalUp * frac.y;
			//velocityFinal.y = vCur * (1.f - frac.x) * (1.f - frac.y) + vRight * frac.x * (1.f - frac.y) + vUp * (1.f - frac.x) * frac.y + vRightUp * frac.x * frac.y;
		}

		return velocityFinal;

	}


	uint2 AdvectionVisCoord{ 10,10 };
	bool bAdvectionVisX = true;

	void DebugVisualizeAdvectionAtAPoint()
	{
		const uint2 coord = AdvectionVisCoord;
		auto& states = VelocityFieldStateBufferCPU;

		const float h = CGConstantsCB.VelocityFieldGridCellLength;

		float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord) + float2{ bAdvectionVisX ? 0 : h * 0.5f, bAdvectionVisX ? h * 0.5f : 0};
		
		RGeometry::GInstance.SetColor(EColor::White);
		RGeometry::GInstance.AddCircleCameraFacing(float3{ curPos, 0 }, h * 0.1f);

		if (states[GetFlatIndexToVelocityArray2D(coord)] == 0)
		{
			return;
		}

		float2 curVelocityFull;

		float2 curVelocity = LoadFromVelocityField(coord);

		//advect either u or v
		if (bAdvectionVisX)//process x component
		{
			auto u = curVelocity.x;

			//get average v
			float vCur = curVelocity.y;
			float vUp = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2(0, 1))];
			float vLeft = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2(1, 0))];
			float vLeftUp = VelocityField_V_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2(1, 0) + uint2(0, 1))];
			float vAvrg = (vCur + vUp + vLeft + vLeftUp) * 0.25f;

			curVelocityFull = float2{ u, vAvrg };

			float2 velFromField = SampleVelocityField(curPos);

			float2 ss = curVelocityFull;

		}
		else
		{
			auto& v = curVelocity.y;

			//get average u
			float uCur = curVelocity.x;
			float uDown = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2(0, 1))];
			float uRight = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2(1, 0))];
			float uRightDown = VelocityField_U_ComponentBufferCPU[GetFlatIndexToVelocityArray2D(coord + uint2(1, 0) - uint2(0, 1))];
			float uAvrg = (uCur + uDown + uRight + uRightDown) * 0.25f;

			curVelocityFull = float2{ uAvrg, v };

			float2 velFromField = SampleVelocityField(curPos);

			float2 ss = curVelocityFull;

		}

		RGeometry::GInstance.AddArrowPolyboard(curPos, curPos + curVelocityFull * CGConstantsCB.DeltaTime);

		//step back

		float2 prevPos = curPos - curVelocityFull * CGConstantsCB.DeltaTime;
		
		RGeometry::GInstance.SetColor(EColor::RedOrange);
		RGeometry::GInstance.AddCircleCameraFacing(float3{ prevPos, 0 }, h * 0.1f);

		RGeometry::GInstance.SetColor(EColor::Red);
		RGeometry::GInstance.AddLine(float3{ curPos, 0 }, float3{ prevPos, 0 });


		//get new velocity based on pos
		float2 newVelocity = SampleVelocityField(prevPos);

		RGeometry::GInstance.AddArrowPolyboard(prevPos, prevPos + newVelocity * CGConstantsCB.DeltaTime);

	}

	uint StreamLineMaxNumSamples = 32;
	uint StreamLineStride = 2;
	float StreamLineStepSize{ 1.f / 60.f };

	void DebugVisualizeVectorFieldStreamLines()
	{
		RGeometry::GInstance.SetColor(EColor::White);

		const float h = CGConstantsCB.VelocityFieldGridCellLength;
		//uint y = 1;
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y += StreamLineStride)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x += StreamLineStride)
			{
				uint2 coord{ x,y };

				float2 startPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord) + float2{ h * 0.5f,h * 0.5f };

				float2 prevPos = startPos;

				const float2 minLimitPos = float2{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset } + float2{ h, h };
				const float2 maxLimitPos = float2{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset } + float2{ h * CGConstantsCB.VelocityFieldGridSize, h * CGConstantsCB.VelocityFieldGridSize };

				for (int i = 0; i < StreamLineMaxNumSamples; i++)
				{
					//get cur velocity
					float2 curVelocity = SampleVelocityField(prevPos);
					//advance
					float2 curPos = prevPos + curVelocity * StreamLineStepSize;

					if (curPos.x < minLimitPos.x || curPos.x > maxLimitPos.x || curPos.y < minLimitPos.y || curPos.y > maxLimitPos.y)
					{
						break;
					}

					//render line
					RGeometry::GInstance.AddLine(float3{ prevPos, 0 }, float3{ curPos, 0 });

					prevPos = curPos;
				}
			}
		}
		
		
	}

	float VelocityFieldCurMagnitudeMax = 0;
	float VelocityFieldCurMagnitudeMin = 0;

	void ComputeVelocityFieldMagnitudesMinMax(bool bSeparateAxises = false)
	{
		VelocityFieldCurMagnitudeMax = VectorGetLength(LoadFromVelocityField(uint2{ 0,0 }));
		VelocityFieldCurMagnitudeMin = VelocityFieldCurMagnitudeMax;
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };

				if (bSeparateAxises)
				{
					float2 curMagnitude = (LoadFromVelocityField(coord));
					curMagnitude.x = fabs(curMagnitude.x);
					curMagnitude.y = fabs(curMagnitude.y);
					VelocityFieldCurMagnitudeMax = std::max(VelocityFieldCurMagnitudeMax, curMagnitude.x);
					VelocityFieldCurMagnitudeMax = std::max(VelocityFieldCurMagnitudeMax, curMagnitude.y);
					VelocityFieldCurMagnitudeMin = std::min(VelocityFieldCurMagnitudeMin, curMagnitude.x);
					VelocityFieldCurMagnitudeMin = std::min(VelocityFieldCurMagnitudeMin, curMagnitude.y);
				}
				else
				{
					float curMagnitude = VectorGetLength(LoadFromVelocityField(coord));
					VelocityFieldCurMagnitudeMax = std::max(VelocityFieldCurMagnitudeMax, curMagnitude);
					VelocityFieldCurMagnitudeMin = std::min(VelocityFieldCurMagnitudeMin, curMagnitude);
				}

			}
		}

		VelocityFieldCurMagnitudeMin = 0.0;
		VelocityFieldCurMagnitudeMax = 500.0;
	}


	void ExtractVelocityFieldDataFromGPU(RCommandList& CmdList)
	{
		GPU_PROFILE_SCOPE(ExtractVelocityFieldDataFromGPUToCPU, CmdList);

		const auto numVelocityGridCells = CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize;
		RReadbackBuffer VelocityField_U_BufferReadBack;
		RReadbackBuffer VelocityField_V_BufferReadBack;

		VelocityField_U_BufferReadBack.Create(TEXT("VelocityField_U_BufferReadBack"), numVelocityGridCells, sizeof(float));
		VelocityField_V_BufferReadBack.Create(TEXT("VelocityField_V_BufferReadBack"), numVelocityGridCells, sizeof(float));

		CmdList.ReadbackTexture2D(VelocityField_U_BufferReadBack, VelocityField_U_ComponentTextureGPU);
		CmdList.ReadbackTexture2D(VelocityField_V_BufferReadBack, VelocityField_V_ComponentTextureGPU);

		CmdList.ExecuteCmdList(true);

		static RDynamicVector<float> VelocityField_U_ComponentBufferCPUInverted;
		VelocityField_U_ComponentBufferCPUInverted.resize(VelocityField_U_ComponentBufferCPU.size());
		static RDynamicVector<float> VelocityField_V_ComponentBufferCPUInverted;
		VelocityField_V_ComponentBufferCPUInverted.resize(VelocityField_V_ComponentBufferCPU.size());

		VelocityField_U_BufferReadBack.ExtractContents(VelocityField_U_ComponentBufferCPUInverted.data(), numVelocityGridCells * sizeof(float));
		VelocityField_V_BufferReadBack.ExtractContents(VelocityField_V_ComponentBufferCPUInverted.data(), numVelocityGridCells * sizeof(float));

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
				uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

				VelocityField_U_ComponentBufferCPU[destIndex] = VelocityField_U_ComponentBufferCPUInverted[sourceIndex];
				VelocityField_V_ComponentBufferCPU[destIndex] = VelocityField_V_ComponentBufferCPUInverted[sourceIndex];

			}
		}

		//extract field states
#if 1
		if (bDebugVisualizeVelocityFieldCellsStates)
		{
			RReadbackBuffer velocityStateBufferReadback;
			velocityStateBufferReadback.Create(TEXT("VelocityField_State_BufferReadBack"), numVelocityGridCells, sizeof(VELOCITY_FIELD_STATE_TYPE));
			CmdList.ReadbackTexture2D(velocityStateBufferReadback, VelocityFieldStateTextureGPU);
			CmdList.ExecuteCmdList(true);

			static RDynamicVector<VELOCITY_FIELD_STATE_TYPE> VelocityField_State_ComponentBufferCPUInverted;
			VelocityField_State_ComponentBufferCPUInverted.resize(VelocityFieldStateBufferCPU.size());
			velocityStateBufferReadback.ExtractContents(VelocityField_State_ComponentBufferCPUInverted.data(), numVelocityGridCells * sizeof(VELOCITY_FIELD_STATE_TYPE));
			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };
					uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
					uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });
					VelocityFieldStateBufferCPU[destIndex] = VelocityField_State_ComponentBufferCPUInverted[sourceIndex];
				}
			}
		}
#endif

		//density
#if USE_3COLOR_DENSITY == 1
		
#else

		RReadbackBuffer DensityBufferReadBack;
		DensityBufferReadBack.Create(TEXT("DensityBufferReadBack"), numVelocityGridCells, sizeof(float));

		CmdList.ReadbackTexture2D(DensityBufferReadBack, DensityTextureGPU);

		CmdList.ExecuteCmdList(true);

		static RDynamicVector<float> densityBufferCPUInverted;
		densityBufferCPUInverted.resize(DensityBufferCPU.size());

		DensityBufferReadBack.ExtractContents(densityBufferCPUInverted.data(), numVelocityGridCells * sizeof(float));

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint sourceIndex = GetFlatIndexToVelocityArray2D(uint2{ x, y });
				uint destIndex = GetFlatIndexToVelocityArray2D(uint2{ x, (CGConstantsCB.VelocityFieldGridSize - 1) - y });

				DensityBufferCPU[destIndex] = densityBufferCPUInverted[sourceIndex];
			}
		}
#endif
	}

	void DebugVisualizeVelocityField(RCommandListCompute& CmdList, bool bSeparateAxises = false)
	{

		if (bProcessVectorFieldOnCPU == 0)
		{
			//extract velocity
			ExtractVelocityFieldDataFromGPU(CmdList);
		}

		const uint visualizedVelocityFieldPrecision = 64;

		const uint stride = CGConstantsCB.VelocityFieldGridSize / visualizedVelocityFieldPrecision;

		float2 boundsStartPos{ CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset };

		if (bSeparateAxises)
		{
			ComputeVelocityFieldMagnitudesMinMax(bSeparateAxises);
		}

#if 1
		//Visualize states

		if (bDebugVisualizeVelocityFieldCellsStates)
		{
			uint gridSize = CGConstantsCB.VelocityFieldGridSize + (stride > 1 ? 1 : 0);

			for (uint y = 0; y < gridSize; y += stride)
			{
				for (uint x = 0; x < gridSize; x += stride)
				{
					if (stride > 1)
					{
						if (x == CGConstantsCB.VelocityFieldGridSize)
						{
							x -= 1;
						}
						if (y == CGConstantsCB.VelocityFieldGridSize)
						{
							y -= 1;
						}
					}
					uint2 coord{ x,y };

					auto curState = VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord)];

					if (curState == VELOCITY_FIELD_CELL_STATE_SOLID || curState == VELOCITY_FIELD_CELL_STATE_DENSITY)
					{
						if (curState == VELOCITY_FIELD_CELL_STATE_SOLID)
						{
							RGeometry::GInstance.SetColor(EColor::RedOrange);
						}
						else
						{
							RGeometry::GInstance.SetColor(EColor::BlueEgyptian);
						}

						//const float h = CGConstantsCB.VelocityFieldGridCellLength;
						const float h = (CGConstantsCB.VelocityFieldBoundsLength) / float(visualizedVelocityFieldPrecision);

						float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
						float2 xCurPos = GetVelocityXComponentWorldPosFromCoord(coord);
						float2 yCurPos = GetVelocityYComponentWorldPosFromCoord(coord);

						RGeometry::GInstance.AddPlane(float3{ curPos + float2{ h * 0.5f,h * 0.5f }, CGConstantsCB.SimulationBoundsStartOffset }, float2{ h * 0.5f, h * 0.5f });
					}
				}
			}
		}
#endif

#if 0 //visualize particle density
		float densityMax = AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(uint2{ 0,0 })];
		float densityMin = densityMax;
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };

				float curDensity = AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
				densityMax = std::max(densityMax, curDensity);
				densityMin = std::min(densityMin, curDensity);
			}
		}
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y += stride)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x += stride)
			{
				uint2 coord{ x,y };

				float curDensity = AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(coord)];
				float3 CurColor = MapTo3ColorRange(curDensity, densityMin, densityMax, GetColorRGBFloat(EColor::BlueEgyptian), GetColorRGBFloat(EColor::YellowCanary), GetColorRGBFloat(EColor::Red));

				RGeometry::GInstance.SetColorRaw(CurColor);

				const float h = (CGConstantsCB.VelocityFieldBoundsLength) / float(visualizedVelocityFieldPrecision);

				float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
				float2 xCurPos = GetVelocityXComponentWorldPosFromCoord(coord);
				float2 yCurPos = GetVelocityYComponentWorldPosFromCoord(coord);

				RGeometry::GInstance.AddPlane(float3{ curPos + float2{ h * 0.5f,h * 0.5f }, CGConstantsCB.SimulationBoundsStartOffset }, float2{ h * 0.5f, h * 0.5f });
			}
		}
#endif

		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y+= stride)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x+= stride)
			{
				uint2 coord{ x,y };

				//const float h = CGConstantsCB.VelocityFieldGridCellLength;
				const float h = (CGConstantsCB.VelocityFieldBoundsLength) / float(visualizedVelocityFieldPrecision);

				float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord);
				float2 xCurPos = GetVelocityXComponentWorldPosFromCoord(coord);
				float2 yCurPos = GetVelocityYComponentWorldPosFromCoord(coord);

				{
					const float2 curVelocity = LoadFromVelocityField(coord);
					const float2 curVelocityNormalized = VectorNormalize(curVelocity) * h;

					//Visualize Magnitude
					if (bSeparateAxises)
					{
						float2 curMagnitude = (curVelocity);

						curMagnitude.x = fabs(curMagnitude.x);
						curMagnitude.y = fabs(curMagnitude.y);

						float3 CurColorX = MapTo3ColorRange(curMagnitude.x, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax);
						float3 CurColorY = MapTo3ColorRange(curMagnitude.y, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax);

						float tX = MapToRange(curMagnitude.x, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax, 0.5f, 1.f);
						float tY = MapToRange(curMagnitude.y, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax, 0.5f, 1.f);

						float2 sign{ curVelocity.x < 0 ? -1.f : 1.f, curVelocity.y < 0 ? -1.f : 1.f };

						//x
						if (curMagnitude.x > 0.01)
						{
							RGeometry::GInstance.SetColorRaw(CurColorX);
							RGeometry::GInstance.AddArrowPolyboard(float3{ xCurPos, CGConstantsCB.SimulationBoundsStartOffset }, float3{ xCurPos + float2{ sign.x * h * tX, 0 }, CGConstantsCB.SimulationBoundsStartOffset });
						}
						
						//y
						if (curMagnitude.y > 0.01)
						{
							RGeometry::GInstance.SetColorRaw(CurColorY);
							RGeometry::GInstance.AddArrowPolyboard(float3{ yCurPos, CGConstantsCB.SimulationBoundsStartOffset }, float3{ yCurPos + float2{ 0, sign.y * h * tY }, CGConstantsCB.SimulationBoundsStartOffset });
						}
						
					}
					else
					{
						const float curMagnitude = VectorGetLength(curVelocity);
						float3 CurColor = MapTo3ColorRange(curMagnitude, VelocityFieldCurMagnitudeMin, VelocityFieldCurMagnitudeMax);

						RGeometry::GInstance.SetColorRaw(CurColor);

						RGeometry::GInstance.AddArrowPolyboard(float3{ curPos + float2{ h * 0.5f, 0.f},  CGConstantsCB.SimulationBoundsStartOffset }, float3{ curPos + float2{ h * 0.5f, 0.f } + curVelocityNormalized,  CGConstantsCB.SimulationBoundsStartOffset });
					}
					
					if (bDebugVisualizeVelocityFieldCells)
					{
						RGeometry::GInstance.SetColor(EColor::White);
						RGeometry::GInstance.AddPlaneLineList(
							float3{ curPos, 0 },
							float3{ curPos, 0 } + float3{ 0, h, 0 },
							float3{ curPos, 0 } + float3{ h, h, 0 },
							float3{ curPos, 0 } + float3{ h, 0, 0 }
						);
					}
					
				}

			}
		}

	}














	/**********************************
			 VELOCITY FIELD GPU
	**********************************/

	template <bool Velocity, bool Density>
	class GenerateVelocityFieldShader : public RShader
	{
	public:
		GenerateVelocityFieldShader()
		{
			if constexpr (Velocity)
			{
				PushBackDefine({ L"GENERATE_VELOCITY", L"1" });
			}
			if constexpr (Density)
			{
				PushBackDefine({ L"GENERATE_DENSITY", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "GenerateVelocityField", L"cs_6_0");
		}
	};


	void GenerateVelocityFieldGPU(RCommandListCompute& CmdList)
	{
		const bool bGenerateVelocity = 1 && Keyboard::IsKeyPressed('V');
		const bool bGenerateDensity = 1 && Keyboard::IsKeyPressed('N');

		if (bGenerateVelocity || bGenerateDensity)
		{
			GPU_PROFILE_SCOPE(GenerateVelocityFieldGPU, CmdList);

			auto& RootSig = ParticlePhysicsRootSig::Get();

			CmdList.SetRootSignature(RootSig);

			CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureGPU, true);

			if (bGenerateVelocity && bGenerateDensity)
			{
				static GenerateVelocityFieldShader<true,true> ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
			}
			else if (bGenerateVelocity)
			{
				static GenerateVelocityFieldShader<true, false> ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
			}
			else
			{
				static GenerateVelocityFieldShader<false, true> ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
			}
			
			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		}

		
	}

	template<bool Texture>
	class FillDensityFieldShader : public RShader
	{
	public:
		FillDensityFieldShader()
		{
			if constexpr (Texture)
			{
				PushBackDefine({ L"COPY_TEXTURE", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "FillDensityField", L"cs_6_0");
		}
	};

	bool bFillDensityFieldlWithColor = true;
	bool bFillDensityFieldWithTexture = true;

	void FillDensityFieldGPU(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(FillDensityFieldGPU, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureGPU, true);

#if USE_3COLOR_DENSITY == 1
		if (bFillDensityFieldWithTexture)
		{
			static RTexture colorSource;
			if (!colorSource.IsValid())
			{
				colorSource.CreateFromBitmap(RBitmap::FromFile(RTexturePool::filepath + std::string{ RTexturePool::Paintings[8] }), DXGI_FORMAT_R8G8B8A8_UNORM);
			}

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensitySourceColorTextureSRV"), colorSource, true);

			static FillDensityFieldShader<true> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		else
#endif//USE_3COLOR_DENSITY
		{
			static FillDensityFieldShader<false> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
	}


	template <bool ClearStates >
	class GenerateVelocityFieldBoundariesShader : public RShader
	{
	public:
		GenerateVelocityFieldBoundariesShader()
		{
			if constexpr (ClearStates)
			{
				PushBackDefine({ L"CLEAR_STATES", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "GenerateVelocityFieldBoundaries", L"cs_6_0");
		}
	};

	void GenerateVelocityFieldBoundariesGPU(RCommandListCompute& CmdList, bool bClearStates = true)
	{
		GPU_PROFILE_SCOPE(GenerateVelocityFieldBoundariesGPU, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureUAV"), VelocityFieldStateTextureGPU, true);

#if PAINT_RENDER_PIXEL_ADVECT
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("InputTextureWithNegativeValuesSRV"), GBufferMetal, true);
#endif

		if (bClearStates)
		{
			static GenerateVelocityFieldBoundariesShader<true> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		else
		{
			static GenerateVelocityFieldBoundariesShader<false> ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));
		}
		

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
	}






	class VectorFieldDivergenceShader : public RShader
	{
	public:
		VectorFieldDivergenceShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ProjectVelocityField", L"cs_6_0");
		}
	};

	void ProjectVectorFieldSimple(RCommandListCompute& CmdList)
	{
		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

		if (CGConstantsCB.bUseParticleBasedAdvection)
		{
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AverageDensityBufferUAV"), AverageDensityBufferGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureUAV"), AdvectedParticleDensityTextureGPU, true);
		}

		//Compute Divergence 

		auto encodeUint2 = [](uint2 val)
		{
			return uint(val.x | (val.y << 16));
		};

		//4 passes per each local id so that neighbors dont touch
		{
			//GPU_PROFILE_SCOPE(GenerateDivergence, CmdList);
			static VectorFieldDivergenceShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / 2 / ThreadGroupSize };

			VectorProjectionPerPassConstants PerPassConstants;

			auto ExecutePassFunc = [&](uint2 curLocalCoordOffset)
			{
				PerPassConstants.Offset = encodeUint2(curLocalCoordOffset);

				CmdList.SetConstantArray(RootSig.GetRootParamIndex("VectorProjectionPerPassConstantsCB"), 1, &PerPassConstants);

				CmdList.TransitionResource(VelocityField_U_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
				CmdList.TransitionResource(VelocityField_V_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			};

			ExecutePassFunc(uint2(0, 1));
			ExecutePassFunc(uint2(1, 1));
			ExecutePassFunc(uint2(0, 0));
			ExecutePassFunc(uint2(1, 0));

		}
	}

	class ComputeDivergenceShader : public RShader
	{
	public:
		ComputeDivergenceShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ComputeDivergence", L"cs_6_0");
		}
	};

	class ComputePressureShader : public RShader
	{
	public:
		ComputePressureShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ComputePressure", L"cs_6_0");
		}
	};

	class CopyPressureShader : public RShader
	{
	public:
		CopyPressureShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "CopyPressureToBorders", L"cs_6_0");
		}
	};

	class ProjectPressureFieldShader : public RShader
	{
	public:
		ProjectPressureFieldShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ProjectPressureField", L"cs_6_0");
		}
	};

	void ComputeDivergence(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ComputeDivergence, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VectorFieldDivergenceTextureUAV"), DivergenceTextureGPU, true);

		static ComputeDivergenceShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
	}

	void ComputePressurePoisson(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ComputePressurePoisson, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VectorFieldDivergenceTextureSRV"), DivergenceTextureGPU, true);

		CmdList.TransitionResource(PressureTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.TransitionResource(PressureTexture2GPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.ClearUAV(PressureTextureGPU);
		CmdList.ClearUAV(PressureTexture2GPU);

		RGpuResource* pCurSource = &PressureTextureGPU;
		RGpuResource* pCurDest = &PressureTexture2GPU;

		VectorFieldProjectionNumIterations = RMath::AlignUp(VectorFieldProjectionNumIterations, 2);

		//Jacobi solver
		for (int i = 0; i < VectorFieldProjectionNumIterations * 2; i++)
		{
			GPU_PROFILE_SCOPE_TEXT(CmdList, L"Pass %i", i);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PressureFieldTextureSRV"), *pCurSource, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PressureFieldTextureUAV"), *pCurDest, true);

			static ComputePressureShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };
			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);

			CmdList.TransitionResource(*pCurDest, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

			//Copy Pressure to Borders
			if(1)
			{
				static CopyPressureShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 64;
				uint numGroupsToDispatch = RMath::DivideAndRoundUp((CGConstantsCB.VelocityFieldGridSize - 1), ThreadGroupSize);
				CmdList.Dispatch(numGroupsToDispatch, 1, 1);
			}

			std::swap(pCurSource, pCurDest);
		}
	}

	void ProjectVectorFieldImplicitMethod(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ProjectVectorFieldImplicitMethod, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		if (!bOutputDivergenceAndPressureResults)
		{
			//Compute Divergence
			ComputeDivergence(CmdList);

			//Compute Pressure
			ComputePressurePoisson(CmdList);
		}

		//Project Velocity
		{
			GPU_PROFILE_SCOPE(ProjectVelocity, CmdList);

			static ProjectPressureFieldShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PressureFieldTextureSRV"), PressureTextureGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PressureGradientTextureUAV"), PressureGradientTextureGPU, true);

			if (CGConstantsCB.bUseParticleBasedAdvection)
			{
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AverageDensityBufferUAV"), AverageDensityBufferGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureUAV"), AdvectedParticleDensityTextureGPU, true);
			}

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		}

	}

	class VelocityAdvectionShader : public RShader
	{
	public:
		VelocityAdvectionShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "AdvectVelocity", L"cs_6_0");
		}
	};

	void AdvectVelocity(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(AdvectVelocity, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		//TODO: Swap instead of copying
		//Make a copy of Velocity Grid to use as SRV
		CmdList.CopyResource(VelocityField_U_ComponentTextureCopyGPU, VelocityField_U_ComponentTextureGPU);
		CmdList.CopyResource(VelocityField_V_ComponentTextureCopyGPU, VelocityField_V_ComponentTextureGPU);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureCopyGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureCopyGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

		static VelocityAdvectionShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
	}

	template <bool bBackwardAdvection = true>
	class DensityAdvectionShader : public RShader
	{
	public:
		DensityAdvectionShader()
		{
			if constexpr (bBackwardAdvection)
			{
				PushBackDefine({ L"BACKWARD_ADVECTION", L"1" });
			}
			Create("Physics/ParticlePhysics.hlsl", "AdvectDensity", L"cs_6_0");
		}
	}; 

	class MacCormackAdvectionShader : public RShader
	{
	public:
		MacCormackAdvectionShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "AdvectDensityMacCormackFinal", L"cs_6_0");
		}
	};

	bool bAdvectDensityMacCormack = false;

	void AdvectDensity(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(AdvectDensity, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());


		uint outputSize = CGConstantsCB.VelocityFieldGridSize;

#if PAINT_RENDER_PIXEL_ADVECT
		outputSize = PaintingBackgroundSizePixels.x;
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("InputTextureWithNegativeValuesSRV"), GBufferMetal, true);
#endif

		

		if (bAdvectDensityMacCormack)
		{
			GPU_PROFILE_SCOPE(AdvectDensityMacCormack, CmdList);

			//phi_n+1 -- DensityTextureCopyGPU
			//phi_n -- DensityTexture2GPU
			//density SRV -- DensityTexture3GPU

			//Make a copy of Density to use as SRV
#if PAINT_RENDER_PIXEL_ADVECT
			CmdList.CopyResource(DensityTexture3GPU, GBufferColor);
#else
			CmdList.CopyResource(DensityTexture3GPU, DensityTextureGPU);
#endif

			//1. advect as usual but into separate buffer
			{
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTexture3GPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureCopyGPU, true);

				static DensityAdvectionShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ outputSize / ThreadGroupSize };

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}

			//2. reverse advect forward into another separate buffer using prev buffer as SRV
			{
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTextureCopyGPU, true);
				
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTexture2GPU, true);

				static DensityAdvectionShader<false> ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ outputSize / ThreadGroupSize };

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}

			//3. combine results and clamp
			{
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTexture3GPU, true);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldMacCormackAdvectedTextureSRV"), DensityTextureCopyGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldMacCormackReversedTextureSRV"), DensityTexture2GPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureGPU, true);

				static MacCormackAdvectionShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ outputSize / ThreadGroupSize };

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}
		}
		else
		{
			GPU_PROFILE_SCOPE(AdvectDensity, CmdList);

			//Make a copy of Density to use as SRV
#if PAINT_RENDER_PIXEL_ADVECT
			CmdList.CopyResource(DensityTextureCopyGPU, GBufferColor);
#else
			CmdList.CopyResource(DensityTextureCopyGPU, DensityTextureGPU);
#endif

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTextureCopyGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureGPU, true);

			static DensityAdvectionShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ outputSize / ThreadGroupSize };

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		}


#if PAINT_RENDER_PIXEL_ADVECT
		CmdList.CopyResource(GBufferColor, DensityTextureGPU);
#endif
		
	}


	class DiffuseVelocityShader : public RShader
	{
	public:
		DiffuseVelocityShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "DiffuseVelocity", L"cs_6_0");
		}
	};

	void DiffuseVelocity(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(DiffuseVelocity__Apply_Viscocity, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		static DiffuseVelocityShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

		CmdList.TransitionResource(VelocityField_U_ComponentTexture2GPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.TransitionResource(VelocityField_V_ComponentTexture2GPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.TransitionResource(VelocityField_U_ComponentTextureCopyGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.TransitionResource(VelocityField_V_ComponentTextureCopyGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
		CmdList.ClearUAV(VelocityField_U_ComponentTexture2GPU);
		CmdList.ClearUAV(VelocityField_V_ComponentTexture2GPU);
		CmdList.ClearUAV(VelocityField_U_ComponentTextureCopyGPU);
		CmdList.ClearUAV(VelocityField_V_ComponentTextureCopyGPU);

		RGpuResource* pCurSourceU = &VelocityField_U_ComponentTexture2GPU;
		RGpuResource* pCurSourceV = &VelocityField_V_ComponentTexture2GPU;

		RGpuResource* pCurDestU = &VelocityField_U_ComponentTextureCopyGPU;
		RGpuResource* pCurDestV = &VelocityField_V_ComponentTextureCopyGPU;
		
		//Jacobi solver
		for (int i = 0; i < DiffusionNumIterations; i++)
		{
			GPU_PROFILE_SCOPE_TEXT(CmdList, L"Pass %i", i);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureNewSRV"), *pCurSourceU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureNewSRV"), *pCurSourceV, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), *pCurDestU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), *pCurDestV, true);

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);

			std::swap(pCurSourceU, pCurDestU);
			std::swap(pCurSourceV, pCurDestV);
		}

		CmdList.CopyResource(VelocityField_U_ComponentTextureGPU, *pCurSourceU);
		CmdList.CopyResource(VelocityField_V_ComponentTextureGPU, *pCurSourceV);

	}


	class DiffuseDensityShader : public RShader
	{
	public:
		DiffuseDensityShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "DiffuseDensity", L"cs_6_0");
		}
	};

	void DiffuseDensity(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(DiffuseDensity, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		static DiffuseDensityShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		CmdList.ClearUAV(DensityTexture2GPU);
		CmdList.ClearUAV(DensityTextureCopyGPU);

		RGpuResource* pCurSource = &DensityTexture2GPU;
		RGpuResource* pCurDest = &DensityTextureCopyGPU;

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTextureGPU, true);

		//Jacobi solver
		for (int i = 0; i < DiffusionNumIterations; i++)
		{
			GPU_PROFILE_SCOPE_TEXT(CmdList, L"Pass %i", i);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureNewSRV"), *pCurSource, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), *pCurDest, true);

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);

			std::swap(pCurSource, pCurDest);
		}

		CmdList.CopyResource(DensityTextureGPU, *pCurSource);

	}



	/*
		VORTICITY
	*/


	class ComputeCurlShader : public RShader
	{
	public:
		ComputeCurlShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ComputeCurl", L"cs_6_0");
		}
	};

	class ApplyVorticityConfinementForceShader : public RShader
	{
	public:
		ApplyVorticityConfinementForceShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ApplyVorticityConfinementForce", L"cs_6_0");
		}
	};

	bool bEnableVorticityConfinement = false;
	uint VorticityNumIterations = 1;

	void Vorticity(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(Vorticity, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		for (int i = 0; i < VorticityNumIterations; i++)
		{
			//ComputeVorticity
			{
				GPU_PROFILE_SCOPE(ComputeVorticity, CmdList);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VectorFieldCurlTextureUAV"), CurlTextureGPU, true);

				static ComputeCurlShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}

			//Apply Vorticity Confining Force
			if (bEnableVorticityConfinement)
			{
				GPU_PROFILE_SCOPE(VorticityConfinement, CmdList);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VectorFieldCurlTextureSRV"), CurlTextureGPU, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurlGradientTextureUAV"), CulGradientTexture, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurlGradientDirectionTextureUAV"), CulGradientDirectionTexture, true);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

				static ApplyVorticityConfinementForceShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}
		}

	}




	/*
	VELOCITY FIELD TO PARTICLES
	*/


	bool bParticleAdvectionCompensateDrift = true;

	class TransferDataFromParticlesToVelocityFieldShader : public RShader
	{
	public:
		TransferDataFromParticlesToVelocityFieldShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "TransferDataFromParticlesToVelocityField", L"cs_6_0");
		}
	};

	class InvScaleEachVelocityByItsWeightShader : public RShader
	{
	public:
		InvScaleEachVelocityByItsWeightShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "InvScaleEachVelocityByItsWeight", L"cs_6_0");
		}
	};

	class ComputeAverageFieldDensityShader : public RShader
	{
	public:
		ComputeAverageFieldDensityShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ComputeAverageFieldDensity", L"cs_6_0");
		}
	};

	TStructuredBuffer<float> AverageDensityBufferGPU;

	void TransferDataFromParticlesToVelocityFieldCS(RCommandListCompute& CmdList)
	{

		GPU_PROFILE_SCOPE(TransferDataFromParticlesToVelocityFieldCS, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentWeightTextureUAV"), VelocityField_U_ComponentWeightTextureGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentWeightTextureUAV"), VelocityField_V_ComponentWeightTextureGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureUAV"), VelocityFieldStateTextureGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureUAV"), AdvectedParticleDensityTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferSRV"), CurParticlesVelocityBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		auto encodeUint2 = [](uint2 val)
		{
			return uint(val.x | (val.y << 16));
		};

		static TransferDataFromParticlesToVelocityFieldShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / 2 / ThreadGroupSize };

		VectorProjectionPerPassConstants PerPassConstants;

		//4 passes per each local id so that neighbors dont collide
		auto ExecutePassFunc = [&](uint2 curLocalCoordOffset)
		{
			PerPassConstants.Offset = encodeUint2(curLocalCoordOffset);

			CmdList.SetConstantArray(RootSig.GetRootParamIndex("VectorProjectionPerPassConstantsCB"), 1, &PerPassConstants);

			CmdList.TransitionResource(VelocityField_U_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(VelocityField_V_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

			CmdList.TransitionResource(VelocityField_U_ComponentWeightTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(VelocityField_V_ComponentWeightTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

			CmdList.TransitionResource(VelocityFieldStateTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		};

		{
			GPU_PROFILE_SCOPE(AccumulateVelocityAndWeights, CmdList);
			ExecutePassFunc(uint2(0, 1));
			ExecutePassFunc(uint2(1, 1));
			ExecutePassFunc(uint2(0, 0));
			ExecutePassFunc(uint2(1, 0));
		}

	}


	template <uint32 OtputComponentIndex, bool bForceAllParticlesMarkStates = false>
	class TransferDataFromParticlesToVelocityFieldShaderVertexShader : public RShader
	{
	public:
		TransferDataFromParticlesToVelocityFieldShaderVertexShader()
		{
			if constexpr (OtputComponentIndex == 0)
			{
				PushBackDefine({ L"OUTPUT_U", L"1" });
			}
			else if constexpr (OtputComponentIndex == 1)
			{
				PushBackDefine({ L"OUTPUT_V", L"1" });
			}
			else if constexpr (OtputComponentIndex == 2)
			{
				PushBackDefine({ L"OUTPUT_C", L"1" });
			}
			if constexpr (bForceAllParticlesMarkStates)
			{
				PushBackDefine({ L"ALL_PARTICLES_MARK_STATE", L"1" });
			}
			static_assert(OtputComponentIndex < 3);
			CreateVS("Physics/ParticlePhysics.hlsl", "TransferDataFromParticlesToVelocityFieldVS", L"vs_6_0");
		}
	};

	template <bool InjectVelocity = true>
	class TransferDataFromParticlesToVelocityFieldShaderPixelShader : public RShader
	{
	public:
		TransferDataFromParticlesToVelocityFieldShaderPixelShader()
		{
			if constexpr (InjectVelocity)
			{
				PushBackDefine({ L"INJECT_VELOCITY", L"1" });
			}
			CreatePS("Physics/ParticlePhysics.hlsl", "TransferDataFromParticlesToVelocityFieldPS", L"ps_6_0");
		}
	};

	bool bAllAdvectedParticlesMarkDensityStateToCells = false;

	void TransferDataFromParticlesToVelocityFieldVSPS(RCommandListGraphics& CmdList)
	{

		GPU_PROFILE_SCOPE(TransferDataFromParticlesToVelocityField_VSPS, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		//define Viewport
		static CD3DX12_VIEWPORT Viewport;
		static CD3DX12_RECT ScissorRect;
		
		const uint2 viewportSize = VelocityField_U_ComponentTextureRenderTargetGPU.Size;
		
		Viewport = CD3DX12_VIEWPORT{ 0.0f, 0.0f, (FLOAT)viewportSize.x, (FLOAT)viewportSize.y };
		ScissorRect = CD3DX12_RECT{ 0, 0, (LONG)viewportSize.x, (LONG)viewportSize.y };

		//Set global states 
		{
			CmdList.SetViewport(Viewport);
			CmdList.SetScissor(ScissorRect);

			CmdList.SetRootSignature(RootSig);
		}

		//Resources
		{
			CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferSRV"), CurParticlesVelocityBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureUAV"), VelocityFieldStateTextureGPU, true);
		}

		{
			GPU_PROFILE_SCOPE(ClearRTs, CmdList);

			CmdList.TransitionResource(VelocityField_U_ComponentTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET);
			CmdList.TransitionResource(VelocityField_V_ComponentTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET);
			CmdList.TransitionResource(VelocityField_U_ComponentWeightTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET);
			CmdList.TransitionResource(VelocityField_V_ComponentWeightTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET);
			CmdList.TransitionResource(AdvectedParticleDensityTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET, true);

			CmdList.ClearRenderTarget(VelocityField_U_ComponentTextureRenderTargetGPU);
			CmdList.ClearRenderTarget(VelocityField_V_ComponentTextureRenderTargetGPU);

			CmdList.ClearRenderTarget(VelocityField_U_ComponentWeightTextureRenderTargetGPU);
			CmdList.ClearRenderTarget(VelocityField_V_ComponentWeightTextureRenderTargetGPU);

			CmdList.ClearRenderTarget(AdvectedParticleDensityTextureRenderTargetGPU);
		}

		RPSOFactory::CommonPipelineStates psoStates;

		//Pipeline definition
		{
			psoStates.RasterDesc = RCommonResources::Get().RasterizerDefault;

			psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;

			psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_POINT;

			psoStates.BlendDesc = RCommonResources::Get().BlendAdditive;

			DXGI_FORMAT rtFormats[2];
			rtFormats[0] = VelocityField_U_ComponentTextureRenderTargetGPU.Format;
			rtFormats[1] = VelocityField_U_ComponentTextureRenderTargetGPU.Format;

			psoStates.NumRenderTargets = 2;
			psoStates.pRenderTargetFormat = rtFormats;

			CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_POINTLIST);
		}

		//Inject U
		{
			D3D12_CPU_DESCRIPTOR_HANDLE rtvHandles[2];
			rtvHandles[0] = VelocityField_U_ComponentTextureRenderTargetGPU.GetRenderTargetView();
			rtvHandles[1] = VelocityField_U_ComponentWeightTextureRenderTargetGPU.GetRenderTargetView();

			CmdList.SetRenderTargets(2, rtvHandles);

			static TransferDataFromParticlesToVelocityFieldShaderVertexShader<0> VS;
			static TransferDataFromParticlesToVelocityFieldShaderPixelShader PS;
			auto PSO = RPSOFactory::GetGraphicsPSO(RootSig, VS, &PS, psoStates);
			CmdList.SetPipelineState(PSO);

			CmdList.DrawInstanced(1, CGConstantsCB.NumParticles * 4);
		}

		//Inject V
		{
			

			D3D12_CPU_DESCRIPTOR_HANDLE rtvHandles[2];
			rtvHandles[0] = VelocityField_V_ComponentTextureRenderTargetGPU.GetRenderTargetView();
			rtvHandles[1] = VelocityField_V_ComponentWeightTextureRenderTargetGPU.GetRenderTargetView();

			CmdList.SetRenderTargets(2, rtvHandles);

			static TransferDataFromParticlesToVelocityFieldShaderVertexShader<1> VS;
			static TransferDataFromParticlesToVelocityFieldShaderPixelShader PS;
			auto PSO = RPSOFactory::GetGraphicsPSO(RootSig, VS, &PS, psoStates);
			CmdList.SetPipelineState(PSO);

			CmdList.DrawInstanced(1, CGConstantsCB.NumParticles * 4);
		}
		
		//Inject density into cell center
		{

			CmdList.SetRenderTarget(AdvectedParticleDensityTextureRenderTargetGPU.GetRenderTargetView());

			if (bAllAdvectedParticlesMarkDensityStateToCells)
			{
				static TransferDataFromParticlesToVelocityFieldShaderVertexShader<2, true> VS;
				static TransferDataFromParticlesToVelocityFieldShaderPixelShader<false> PS;
				auto PSO = RPSOFactory::GetGraphicsPSO(RootSig, VS, &PS, psoStates);
				CmdList.SetPipelineState(PSO);
			}
			else
			{
				static TransferDataFromParticlesToVelocityFieldShaderVertexShader<2> VS;
				static TransferDataFromParticlesToVelocityFieldShaderPixelShader<false> PS;
				auto PSO = RPSOFactory::GetGraphicsPSO(RootSig, VS, &PS, psoStates);
				CmdList.SetPipelineState(PSO);
			}
			

			CmdList.DrawInstanced(1, CGConstantsCB.NumParticles * 4);
		}

		{
			GPU_PROFILE_SCOPE(CopyBackResults, CmdList);
			CmdList.CopyResource(VelocityField_U_ComponentTextureGPU, VelocityField_U_ComponentTextureRenderTargetGPU);
			CmdList.CopyResource(VelocityField_V_ComponentTextureGPU, VelocityField_V_ComponentTextureRenderTargetGPU);
			CmdList.CopyResource(VelocityField_U_ComponentWeightTextureGPU, VelocityField_U_ComponentWeightTextureRenderTargetGPU);
			CmdList.CopyResource(VelocityField_V_ComponentWeightTextureGPU, VelocityField_V_ComponentWeightTextureRenderTargetGPU);
			CmdList.CopyResource(AdvectedParticleDensityTextureGPU, AdvectedParticleDensityTextureRenderTargetGPU);
		}

	}

	bool bTransferParticleDataWithComputeShader = false;

	void TransferDataFromParticlesToVelocityField(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(TransferDataFromParticlesToVelocityField, CmdList);

		//Reset Boundaries
		GenerateVelocityFieldBoundariesGPU(CmdList);

		//Clear Velocity Field
		{
			GPU_PROFILE_SCOPE(ClearResources, CmdList);
			CmdList.TransitionResource(VelocityField_U_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(VelocityField_V_ComponentTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(VelocityField_U_ComponentWeightTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(VelocityField_V_ComponentWeightTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.TransitionResource(AdvectedParticleDensityTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
			CmdList.ClearUAV(VelocityField_U_ComponentTextureGPU);
			CmdList.ClearUAV(VelocityField_V_ComponentTextureGPU);
			CmdList.ClearUAV(VelocityField_U_ComponentWeightTextureGPU);
			CmdList.ClearUAV(VelocityField_V_ComponentWeightTextureGPU);
			CmdList.ClearUAV(AdvectedParticleDensityTextureGPU);
		}
		
		auto computeFinishedFence = CmdList.ExecuteCmdList();

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		//Compute Divergence
		{
			if (bTransferParticleDataWithComputeShader)
			{
				TransferDataFromParticlesToVelocityFieldCS(CmdList);
			}
			else
			{
				RCommandListGraphics& CmdListGraphics = RCommandList::BeginNew().GetGraphicsContext();
				SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdListGraphics.GetType()).QueueWaitForAnotherFence(computeFinishedFence);

				TransferDataFromParticlesToVelocityFieldVSPS(CmdListGraphics);
				
				auto graphicsFinishedFence = CmdListGraphics.ExecuteCmdListAndReleaseContext();
				SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdList.GetType()).QueueWaitForAnotherFence(graphicsFinishedFence);

				GenerateVelocityFieldBoundariesGPU(CmdList, false);
			}

			{
				GPU_PROFILE_SCOPE(DivideFinalVelocityByWeight, CmdList);
				//Divide final velocities by final weights
				static InvScaleEachVelocityByItsWeightShader ComputeShader2;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader2));
				CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureUAV"), VelocityField_U_ComponentTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureUAV"), VelocityField_V_ComponentTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentWeightTextureUAV"), VelocityField_U_ComponentWeightTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentWeightTextureUAV"), VelocityField_V_ComponentWeightTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentWeightTextureUAV"), VelocityField_U_ComponentWeightTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentWeightTextureUAV"), VelocityField_V_ComponentWeightTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);
				uint numGroupsToDispatch = CGConstantsCB.VelocityFieldGridSize / 8;
				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}

			{
				GPU_PROFILE_SCOPE(ComputeAverageDensity, CmdList);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AverageDensityBufferUAV"), AverageDensityBufferGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureUAV"), AdvectedParticleDensityTextureGPU, true);

				static ComputeAverageFieldDensityShader ComputeShader2;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader2));

				CmdList.Dispatch(1, 1, 1);
			}
			
			//Cache un-projected velocity field
			CmdList.CopyResource(VelocityField_U_ComponentTextureCopyGPU, VelocityField_U_ComponentTextureGPU);
			CmdList.CopyResource(VelocityField_V_ComponentTextureCopyGPU, VelocityField_V_ComponentTextureGPU);

		}
	}

	class TransferDataFromVelocityFieldToParticlesShader : public RShader
	{
	public:
		TransferDataFromVelocityFieldToParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "TransferDataFromVelocityFieldToParticles", L"cs_6_0");
		}
	};

	void TransferDataFromVelocityFieldToParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(TransferDataFromVelocityFieldToParticles, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferUAV"), CurParticlesVelocityBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureNewSRV"), VelocityField_U_ComponentTextureCopyGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureNewSRV"), VelocityField_V_ComponentTextureCopyGPU, true);
		
		static TransferDataFromVelocityFieldToParticlesShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 64;
		uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch, 1, 1);
	}


	class ComputeParticlesColorsShader : public RShader
	{
	public:
		ComputeParticlesColorsShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ComputeParticlesColors", L"cs_6_0");
		}
	};

	bool bAssignParticlesColorsBasedOnDensity = false;

	void ComputeParticlesColors(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ComputeParticlesColors, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferUAV"), ParticlesColorsBufferGPU, true);

		//SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureUAV"), AdvectedParticleDensityTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("AdvectedParticleDensityTextureSRV"), AdvectedParticleDensityTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		static ComputeParticlesColorsShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 64;
		uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch, 1, 1);
	}

	void ExtractParticlesDataFromGPU(RCommandList& CmdList)
	{
		GPU_PROFILE_SCOPE(ExtractParticlesDataFromGPU, CmdList);

		RReadbackBuffer particlesReadbackBuf;
		RReadbackBuffer particlesVelReadbackBuf;
		particlesReadbackBuf.Create(TEXT("particlesReadbackBuf"), CGConstantsCB.NumParticles, sizeof(float3));
		particlesVelReadbackBuf.Create(TEXT("particlesReadbackBuf"), CGConstantsCB.NumParticles, sizeof(float3));

		CmdList.CopyBufferRegion(particlesReadbackBuf, 0, ParticlesBufferPositionCurGPU, 0, CGConstantsCB.NumParticles * sizeof(float3));
		CmdList.CopyBufferRegion(particlesVelReadbackBuf, 0, CurParticlesVelocityBufferGPU, 0, CGConstantsCB.NumParticles * sizeof(float3));

		CmdList.ExecuteCmdList(true);

		ParticlesBufferPositionCurCPU.resize(CGConstantsCB.NumParticles);
		CurParticlesVelocityBufferCPU.resize(CGConstantsCB.NumParticles);

		particlesReadbackBuf.ExtractContents(ParticlesBufferPositionCurCPU.data());
		particlesVelReadbackBuf.ExtractContents(CurParticlesVelocityBufferCPU.data());
	}

	RDynamicVector<float> VelocityWeights_U_BufferCPU;
	RDynamicVector<float> VelocityWeights_V_BufferCPU;


	void TransferDataFromParticlesToVelocityFieldCPU(RCommandListCompute& CmdList)
	{
		
		if (VelocityWeights_U_BufferCPU.size() < VelocityField_U_ComponentBufferCPU.size())
		{
			VelocityWeights_U_BufferCPU.resize(VelocityField_U_ComponentBufferCPU.size());
			VelocityWeights_V_BufferCPU.resize(VelocityField_V_ComponentBufferCPU.size());
		}

		uint n = CGConstantsCB.VelocityFieldGridSize;
		float h = CGConstantsCB.VelocityFieldGridCellLength;
		float hHalf = h * 0.5f;

		//Clear Velocity Field
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);

				StoreInVelocityField(i, float2{ 0,0 });
			}
		}

		//Clear weights
		std::fill(VelocityWeights_U_BufferCPU.begin(), VelocityWeights_U_BufferCPU.end(), 0);
		std::fill(VelocityWeights_V_BufferCPU.begin(), VelocityWeights_V_BufferCPU.end(), 0);


		//Update Boundaries & States
		VelocityFieldSetBoundaries();
		for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
		{
			for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
			{
				uint2 coord{ x,y };
				uint i = GetFlatIndexToVelocityArray2D(coord);
				VelocityFieldStateBufferCPU[i] = VelocityFieldStateBufferCPU[i] == 0 ? VELOCITY_FIELD_CELL_STATE_SOLID : VELOCITY_FIELD_CELL_STATE_AIR;
			}
		}


		//Fill Density Cells with Particles
		for (int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			float3& curParticlePos = ParticlesBufferPositionCurCPU[i];

#if PARTICLE_ADVECTION_DEBUG
			{
				curParticlePos = CGConstantsCB.DefaultParticlePos;
			}
#endif//PARTICLE_ADVECTION_DEBUG

			uint2 curCoord;
			uint2 curCoordNext;
			float2 frac;
			GetVelocityFieldCoordFromWorldPos(curParticlePos.xy(), curCoord, curCoordNext, frac);

			uint flatIndex = GetFlatIndexToVelocityArray2D(curCoord);
			if (VelocityFieldStateBufferCPU[flatIndex] == VELOCITY_FIELD_CELL_STATE_AIR)
			{
				VelocityFieldStateBufferCPU[flatIndex] = VELOCITY_FIELD_CELL_STATE_DENSITY;
			}
		}


		//for each velocity vector dimension
		for (int c = 0; c < 2; c++)
		{
			float dx = c == 0 ? 0.0 : hHalf;
			float dy = c == 0 ? hHalf : 0.0;

			auto& f = c == 0 ? VelocityField_U_ComponentBufferCPU : VelocityField_V_ComponentBufferCPU;
			auto& d = c == 0 ? VelocityWeights_U_BufferCPU : VelocityWeights_V_BufferCPU;

			for (int i = 0; i < CGConstantsCB.NumParticles; i++)
			{
				float3& curParticlePos = ParticlesBufferPositionCurCPU[i];
				float3& curParticleVel = CurParticlesVelocityBufferCPU[i];

#if PARTICLE_ADVECTION_DEBUG
				float3 defaultVelocity = PARTICLE_ADVECTION_DEBUG_VELOCITY;
				{
					curParticlePos = CGConstantsCB.DefaultParticlePos;
					curParticleVel = defaultVelocity;
				}
#endif//PARTICLE_ADVECTION_DEBUG

				float curParticleVelComponent = c == 0 ? curParticleVel.x : curParticleVel.y;

				//convert into staggered grid space
				float2 curParticlePosMinusHalf = curParticlePos.xy() - float2{ dx, dy };

				uint2 curCoord;
				uint2 curCoordNext;
				float2 frac;
				GetVelocityFieldCoordFromWorldPos(curParticlePosMinusHalf, curCoord, curCoordNext, frac);

				float dCur = (1.f - frac.x) * (1.f - frac.y);
				float dUp = (1.f - frac.x) * frac.y;
				float dRight = frac.x * (1.f - frac.y);
				float dRightUp = frac.x * frac.y;

				//inject weighted particle velocity into the field

				auto& vCur = f[GetFlatIndexToVelocityArray2D(curCoord)];
				auto& vUp = f[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
				auto& vRight = f[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
				auto& vRightUp = f[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

				vCur += curParticleVelComponent * dCur;
				vUp += curParticleVelComponent * dUp;
				vRight += curParticleVelComponent * dRight;
				vRightUp += curParticleVelComponent * dRightUp;

				//accumulate weights
				d[GetFlatIndexToVelocityArray2D(curCoord)] += dCur;
				d[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })] += dUp;
				d[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })] += dRight;
				d[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })] += dRightUp;

			}


			//divide by accumulated weights
			for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
			{
				for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
				{
					uint2 coord{ x,y };
					uint i = GetFlatIndexToVelocityArray2D(coord);

					if (d[i] > 0.f)
					{
						f[i] /= d[i];
					}

					//restore solid cells velocities

					bool solid = VelocityFieldStateBufferCPU[i] == VELOCITY_FIELD_CELL_STATE_SOLID;

					//if cur cell solid or cell to the left is solid
					if (solid || x > 0 && VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2(1,0))] == VELOCITY_FIELD_CELL_STATE_SOLID)
					{
						VelocityField_U_ComponentBufferCPU[i] = 0.f;
					}
					//if cur cell solid or downward cell is solid
					if (solid || y > 0 && VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(coord - uint2(0, 1))] == VELOCITY_FIELD_CELL_STATE_SOLID)
					{
						VelocityField_V_ComponentBufferCPU[i] = 0.f;
					}
				}
			}

		}

		//Cache un-projected Velocity Field 
		VelocityField_U_ComponentBufferCopyCPU = VelocityField_U_ComponentBufferCPU;
		VelocityField_V_ComponentBufferCopyCPU = VelocityField_V_ComponentBufferCPU;

	}


	void ComputeParticlesColorsCPU()
	{

		auto curVecSize = ParticlesColorsBufferCPU.size();
		if (curVecSize < CGConstantsCB.NumParticles)
		{
			ParticlesColorsBufferCPU.resize(CGConstantsCB.NumParticles);
			const float3 defaultColor = float3{ 1.f,1.f,1.f };
			std::fill(ParticlesColorsBufferCPU.begin() + curVecSize, ParticlesColorsBufferCPU.end(), defaultColor);
		}


		for (int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			float s = 0.01f;

			float3& curParticleColor = ParticlesColorsBufferCPU[i];

			curParticleColor.x = RMath::Clamp(curParticleColor.x - s, 0.f, 1.f);
			curParticleColor.y = RMath::Clamp(curParticleColor.y - s, 0.f, 1.f);
			curParticleColor.z = RMath::Clamp(curParticleColor.z + s, 0.f, 1.f);

			float3& curParticlePos = ParticlesBufferPositionCurCPU[i];
			
			uint2 curCoord;
			uint2 curCoordNext;
			float2 frac;
			GetVelocityFieldCoordFromWorldPos(curParticlePos.xy(), curCoord, curCoordNext, frac);

			uint flatCellIndex = GetFlatIndexToVelocityArray2D(curCoord);

			float d0 = AverageParticleDensity;

			if (d0 > 0.0) 
			{
				float relDensity = AdvectedParticleDensityBufferCPU[flatCellIndex] / d0;
				if (relDensity < 0.7) 
				{
					const float sc = 0.8;
					curParticleColor.x = sc;
					curParticleColor.y = sc;
					curParticleColor.z = 1.0;
				}
			}
		}


	}


	void TransferDataFromVelocityFieldToParticlesCPU()
	{

		uint n = CGConstantsCB.VelocityFieldGridSize;
		float h = CGConstantsCB.VelocityFieldGridCellLength;
		float hHalf = h * 0.5f;

		//for each velocity vector dimension
		for (int c = 0; c < 2; c++)
		{
			float dx = c == 0 ? 0.0 : hHalf;
			float dy = c == 0 ? hHalf : 0.0;

			auto& f = c == 0 ? VelocityField_U_ComponentBufferCPU : VelocityField_V_ComponentBufferCPU;
			auto& prevF = c == 0 ? VelocityField_U_ComponentBufferCopyCPU : VelocityField_V_ComponentBufferCopyCPU;
			auto& d = c == 0 ? VelocityWeights_U_BufferCPU : VelocityWeights_V_BufferCPU;

			for (int i = 0; i < CGConstantsCB.NumParticles; i++)
			{
				float3& curParticlePos = ParticlesBufferPositionCurCPU[i];
				float3& curParticleVel = CurParticlesVelocityBufferCPU[i];

#if PARTICLE_ADVECTION_DEBUG
				{
					curParticlePos = CGConstantsCB.DefaultParticlePos;
				}
#endif//PARTICLE_ADVECTION_DEBUG

				float& curParticleVelComponent = c == 0 ? curParticleVel.x : curParticleVel.y;

				//convert into staggered grid space
				float2 curParticlePosMinusHalf = curParticlePos.xy() - float2{ dx, dy };

				uint2 curCoord;
				uint2 curCoordNext;
				float2 frac;
				GetVelocityFieldCoordFromWorldPos(curParticlePosMinusHalf, curCoord, curCoordNext, frac);

				float dCur = (1.f - frac.x) * (1.f - frac.y);
				float dUp = (1.f - frac.x) * frac.y;
				float dRight = frac.x * (1.f - frac.y);
				float dRightUp = frac.x * frac.y;

				//inject weighted particle velocity into the field

				auto& vCur = f[GetFlatIndexToVelocityArray2D(curCoord)];
				auto& vUp = f[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
				auto& vRight = f[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
				auto& vRightUp = f[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

				auto& vCurPrev = prevF[GetFlatIndexToVelocityArray2D(curCoord)];
				auto& vUpPrev = prevF[GetFlatIndexToVelocityArray2D(uint2{ curCoord.x, curCoordNext.y })];
				auto& vRightPrev = prevF[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoord.y })];
				auto& vRightUpPrev = prevF[GetFlatIndexToVelocityArray2D(uint2{ curCoordNext.x, curCoordNext.y })];

				auto coordUp = uint2{ curCoord.x, curCoordNext.y };
				auto coordRight = uint2{ curCoordNext.x, curCoord.y };
				auto coordRightUp = uint2{ curCoordNext.x, curCoordNext.y };

				auto GetCellState = [&](uint2 _coord)
				{
					return VelocityFieldStateBufferCPU[GetFlatIndexToVelocityArray2D(_coord)];
				};

				uint2 offset;
				offset.x = c == 0 ? 1 : 0;
				offset.y = c == 0 ? 0 : 1;
				int stateCur = GetCellState(curCoord) != VELOCITY_FIELD_CELL_STATE_AIR || GetCellState(curCoord - offset) != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
				int stateUp = GetCellState(coordUp) != VELOCITY_FIELD_CELL_STATE_AIR || GetCellState(coordUp - offset) != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
				int stateRight = GetCellState(coordRight) != VELOCITY_FIELD_CELL_STATE_AIR || GetCellState(coordRight - offset) != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
				int stateRightUp = GetCellState(coordRightUp) != VELOCITY_FIELD_CELL_STATE_AIR || GetCellState(coordRightUp - offset) != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;

				check((stateCur + stateUp + stateRight + stateRightUp) > 0, "at least 1 cell should contain density");

				float d = stateCur * dCur + stateUp * dUp + stateRight * dRight + stateRightUp * dRightUp;

				if (d > 0.0) 
				{

					float picV = (stateCur * dCur * vCur + stateUp * dUp * vUp + stateRight * dRight * vRight + stateRightUp * dRightUp * vRightUp) / d;

					float corr = (stateCur * dCur * (vCur - vCurPrev) + stateUp * dUp * (vUp - vUpPrev) + stateRight * dRight * (vRight - vRightPrev) + stateRightUp * dRightUp * (vRightUp - vRightUpPrev)) / d;
					float flipV = curParticleVelComponent + corr;

					curParticleVelComponent = (1.0f - CGConstantsCB.ParticleAdvectionPICFLIPRatio) * picV + CGConstantsCB.ParticleAdvectionPICFLIPRatio * flipV;
				}

			}
		}
	}



	RDynamicVector<float> AdvectedParticleDensityBufferCPU;

	float AverageParticleDensity = 0;

	void ComputeAdvectedParticleDensityCPU()
	{
		uint n = CGConstantsCB.VelocityFieldGridSize;
		float h = CGConstantsCB.VelocityFieldGridCellLength;
		float hHalf = h * 0.5f;

		if (AdvectedParticleDensityBufferCPU.empty())
		{
			AdvectedParticleDensityBufferCPU.resize(CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize);
		}

		std::fill(AdvectedParticleDensityBufferCPU.begin(), AdvectedParticleDensityBufferCPU.end(), 0.f);

		for (int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			float3& curParticlePos = ParticlesBufferPositionCurCPU[i];

			//offset before floor-ing
			float2 curParticlePosMinusHalf = curParticlePos.xy() - float2{ hHalf, hHalf };

			uint2 curCoord;
			uint2 curCoordNext;
			float2 frac;
			GetVelocityFieldCoordFromWorldPos(curParticlePosMinusHalf, curCoord, curCoordNext, frac);

			float dCur = (1.f - frac.x) * (1.f - frac.y);
			float dUp = (1.f - frac.x) * frac.y;
			float dRight = frac.x * (1.f - frac.y);
			float dRightUp = frac.x * frac.y;

			uint2 curCoordRight{ curCoordNext.x, curCoord.y };
			uint2 curCoordUp{ curCoord.x, curCoordNext.y };

			AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(curCoord)] += dCur;
			AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(curCoordRight)] += dRight;
			AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(curCoordUp)] += dUp;
			AdvectedParticleDensityBufferCPU[GetFlatIndexToVelocityArray2D(curCoordNext)] += dRightUp;
		}

		//if (AverageParticleDensity == 0.0) 
		{
			float sum = 0.0;
			float numFluidCells = 0;

			uint gridSize = CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridSize;

			for (int i = 0; i < gridSize; i++)
			{
				if (VelocityFieldStateBufferCPU[i] == VELOCITY_FIELD_CELL_STATE_DENSITY)
				{
					sum += AdvectedParticleDensityBufferCPU[i];
					numFluidCells++;
				}
			}

			if (numFluidCells > 0)
			{
				AverageParticleDensity = sum / numFluidCells;
			}
		}

	}

	void UploadParticlesDataToGPU(RCommandList& CmdList)
	{
		GPU_PROFILE_SCOPE(UploadParticlesDataToGPU, CmdList);

		//Upload Particle Velocity to GPU
		CmdList.UploadToBuffer(CurParticlesVelocityBufferGPU, 0, CurParticlesVelocityBufferCPU.data(), CGConstantsCB.NumParticles * sizeof(float3));

		CmdList.UploadToBuffer(ParticlesColorsBufferGPU, 0, ParticlesColorsBufferCPU.data(), CGConstantsCB.NumParticles * sizeof(float3));
	}





	class AdvectParticlesWithVelocityFieldShader : public RShader
	{
	public:
		AdvectParticlesWithVelocityFieldShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "AdvectParticlesWithVelocityField", L"cs_6_0");
		}
	};

	//Direct advection, not velocity transfer
	void AdvectParticlesWithVelocityField(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(AdvectParticlesWithVelocityField, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV"), VelocityField_U_ComponentTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV"), VelocityField_V_ComponentTextureGPU, true);

		///
		///	Shader
		///
		static AdvectParticlesWithVelocityFieldShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint kThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, kThreadGroupSize), 1, 1);
	}

	bool bClothSimulateWaterSurface = false;//combine cloth-constraint simulation with 2D velocity field fluid simulation
	bool bSimulateWavesSimple = false;

	class ApplyPressureFieldOffsetToClothParticlesShader : public RShader
	{
	public:
		ApplyPressureFieldOffsetToClothParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ApplyPressureFieldOffsetToClothParticles", L"cs_6_0");
		}
	};

	void ApplyPressureFieldOffsetToSurfaceParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ApplyPressureFieldOffsetToSurfaceParticles, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		//SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionAfterCreationGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("PressureFieldTextureSRV"), PressureTextureGPU, true);

		///
		///	Shader
		///
		static ApplyPressureFieldOffsetToClothParticlesShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint2 numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.ClothSize.x, ThreadGroupSize), RMath::DivideAndRoundUp(CGConstantsCB.ClothSize.y, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);
	}



	class WaveEquationSolverShader : public RShader
	{
	public:
		WaveEquationSolverShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "SolveWaveEquation", L"cs_6_0");
		}
	};

	void SimulateWave(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(SimulateWave, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

		///
		///	Shader
		///
		static WaveEquationSolverShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint2 numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.ClothSize.x, ThreadGroupSize), RMath::DivideAndRoundUp(CGConstantsCB.ClothSize.y, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);
	}


	class NegativeValuesConvertShader : public RShader
	{
	public:
		NegativeValuesConvertShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "VisualizerConvertTextureWithNegativeValuesToColor", L"cs_6_0");
		}
	};

	RTexture& VisualizerConvertTextureWithNegativeValuesToColor(RCommandListCompute& CmdList, RTexture& inputTexture)
	{
		static RTexture colorTexture;

		if (colorTexture.Size.x != inputTexture.Size.x || colorTexture.Size.y != inputTexture.Size.y)
		{
			colorTexture.Create(inputTexture.Size.x, inputTexture.Size.y, DXGI_FORMAT_R32G32_FLOAT, true);
		}

		GPU_PROFILE_SCOPE(ConvertTextureForVisualization, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("InputTextureWithNegativeValuesSRV"), inputTexture, true);
		
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("VisualizerColorTextureUAV"), colorTexture, true);


		///
		///	Shader
		///
		static NegativeValuesConvertShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint2 numGroupsToDispatch{ RMath::DivideAndRoundUp(inputTexture.Size.x, ThreadGroupSize), RMath::DivideAndRoundUp(inputTexture.Size.y, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);

		return colorTexture;
	}


	bool bApplyVorticityMaskToDensityVisualizer = true;

	class VorticityMaskShader : public RShader
	{
	public:
		VorticityMaskShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "DensityApplyVorticityGradientMask", L"cs_6_0");
		}
	};

	void DensityApplyVorticityGradientMask(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(DensityApplyVorticityGradientMask, CmdList);

		ParallelReductionPass::TextureMinMax(CmdList, CulGradientTexture, CulGradientMinMaxTexture);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		static VorticityMaskShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint numGroupsToDispatch{ CGConstantsCB.VelocityFieldGridSize / ThreadGroupSize };

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CurlGradientTextureSRV"), CulGradientTexture, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CurlGradientMinMaxTextureSRV"), CulGradientMinMaxTexture, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureSRV"), DensityTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV"), VelocityFieldStateTextureGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("DensityFieldTextureUAV"), DensityTextureCopyGPU, true);

		CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);

	}




	class DisableParticlesShader : public RShader
	{
	public:
		DisableParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "DisableParticles", L"cs_6_0");
		}
	};

	void DisableParticlesIntersectingBounds(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(DisableParticlesIntersectingBounds, CmdList);

		auto& RootSig = ParticlePhysicsRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferUAV"), ParticlesStateBufferGPU, true);

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferUAV"), ParticlesColorsBufferGPU, true);

		static DisableParticlesShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 64;
		CmdList.Dispatch(RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize), 1, 1);
	}


	//=============================================
	//					SDF Rendering
	//=============================================

	class SDFRenderingRootSig : public RRootSignatureExt, public RSingleton<SDFRenderingRootSig>
	{
	public:
		SDFRenderingRootSig()
		{
			///Constants
			AddConstantBuffer("CGConstantsCB", 0);
			AddConstantBuffer("CGSDFConstantsCB", 1);

			AddConstants("PaintDiffusionPerPassConstantsCB", 2, 2);

			AddConstantBuffer("CGPaintRenderConstantsCB", 10);

			///SRVs
			AddSRV("PlanesBufferSRV", 0);
			AddSRV("ParticlesBufferPositionCurSRV", 1);
			AddSRV("CellIdOffsetAndCountBufferSRV", 2);
			AddSRV("SDFTextureSRV", 3);
			AddSRV("ParticlesColorsBufferSRV", 4);
			AddSRV("SDFOutputTextureSRV", 5);
			AddSRV("SDFOutputTextureNewSRV", 6);
			AddSRV("WetnessParameterTextureSRV", 8);
			AddSRV("ParticlesStateBufferSRV", 7);
			AddSRV("CurParticlesVelocityBufferSRV", 18);

			AddSRV("GBufferColorTextureSRV", 10);
			AddSRV("GBufferNormalTextureSRV", 11);
			AddSRV("GBufferRoughnessTextureSRV", 12);
			AddSRV("GBufferWorldSpaceZSRV", 13);
			AddSRV("GBufferStencilSRV", 14);

			AddSRV("ParticleRenderVertexBufferSRV", 15);
			AddSRV("ParticleTexCoordsBufferSRV", 16);
			AddSRV("ParticleSourceTextureSRV", 17);

			AddSRV("AdvectedParticleDensityTextureSRV", 20);
			AddSRV("AverageDensityBufferSRV", 21);

			AddSRV("GBufferMetalTextureSRV", 22);

			///UAVs
			AddUAV("ParticlesBufferPositionCurUAV", 0);
			AddUAV("SDFOutputTextureUAV", 21);
			AddUAV("SDFTextureUAV", 1);
			AddUAV("SDFOutputNormalsTextureUAV", 2);
			AddUAV("PaintMaxWetnessTextureUAV", 3);
			AddUAV("WetnessParameterTextureUAV", 4);
			AddUAV("OriginalPaintingTextureUAV", 5);
			AddUAV("ParticleSourceTextureUAV", 6);

			AddUAV("ParticlesStateBufferUAV", 7);

			AddUAV("ParticlesColorsBufferUAV", 20);
			AddUAV("CurParticlesVelocityBufferUAV", 15);
			AddUAV("CurParticlesExternalForceBufferUAV", 9);



			//Samplers
			AddStaticSampler(0, RCommonResources::Get().SamplerLinearClampDesc);
			AddStaticSampler(1, RCommonResources::Get().SamplerPointClampDesc);

			Finalize(L"SDFRenderingRootSig");
		}
	};

	RTexture SDFOutputTextureGPU{ TEXT("SDFOutputTexture") };
	RTexture SDFOutputTextureRenderTargetGPU{ TEXT("SDFOutputTexture") };
	RTexture SDFOutputTexture2GPU{ TEXT("SDFOutputTexture2") };
	RTexture SDFOutputTextureCopyGPU{ TEXT("SDFOutputTextureCopy") };
	RTexture PaintMaxWetnessTextureGPU{ TEXT("PaintMaxWetnessTexture") };
	RTexture WetnessParameterTextureGPU{ TEXT("WetnessParameterTextureGPU") };
	RTexture WetnessParameterTexture2GPU{ TEXT("WetnessParameterTexture2GPU") };
	RTexture SDFOutputNormalsTextureGPU{ TEXT("SDFOutputNormalsTexture") };
	RTexture SDFOutputNormalsTexture2GPU{ TEXT("SDFOutputNormalsTexture") };
	RTexture SDFTextureGPU{ TEXT("SDFTexture") };
	RTexture SDFTextureRenderTargetGPU{ TEXT("SDFTexture") };
	RTexture SDFTexture2GPU{ TEXT("SDFTexture2") };

	RTexture SDFOutputFallingParticlesTextureGPU{ TEXT("SDFOutputFallingParticlesTexture") };
	RTexture SDFOutputFallingParticlesTextureRenderTargetGPU{ TEXT("SDFOutputFallingParticlesTextureRenderTarget") };
	RTexture SDFOutputFallingParticlesNormalsTextureGPU{ TEXT("SDFOutputFallingParticlesNormalsTexture") };

	RTexture PaintingSource;

	RDepthBuffer SDFOutputDepthStencilBufferGPU{ 1.f };

	///CBuffer Constants
	SDFRenderConstants CGSDFConstantsCB;
	TCBuffer<SDFRenderConstants> CGSDFConstantsCBufferGPU;


	bool bRenderPainting = false;
	bool bDisableParticlesIntersectingBounds = false;

	class GenerateSignedDistanceFieldShader : public RShader
	{
	public:
		GenerateSignedDistanceFieldShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "GenerateSignedDistanceField", L"cs_6_0");
		}
	};

	class RenderFallingParticlesShader : public RShader
	{
	public:
		RenderFallingParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "RenderFallingParticles", L"cs_6_0");
		}
	};

	class GenerateNormalsFromSDFShader : public RShader
	{
	public:
		GenerateNormalsFromSDFShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "GenerateNormalsFromSDFHeightMap", L"cs_6_0");
		}
	};

	class ModifyPaintParticlesShader : public RShader
	{
	public:
		ModifyPaintParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ModifyPaintParticles", L"cs_6_0");
		}
	};

	class ModifyPaintParticlesWithGBufferShader : public RShader
	{
	public:
		ModifyPaintParticlesWithGBufferShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ModifyPaintParticlesWithGBuffer", L"cs_6_0");
		}
	};

	class PaintDiffusionShader : public RShader
	{
	public:
		PaintDiffusionShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "PaintDiffusion", L"cs_6_0");
		}
	};


	template <bool bCopyFrom>
	class CopyWetnessParameterShader : public RShader
	{
	public:
		CopyWetnessParameterShader()
		{
			if constexpr (bCopyFrom)
			{
				PushBackDefine({ L"COPY_FROM", L"1" });
			}
			else
			{
				PushBackDefine({ L"COPY_FROM", L"0" });
			}
			Create("Physics/ParticlePhysics.hlsl", "CopyWetnessParameter", L"cs_6_0");
		}
	};

	bool bNeablePaintDiffusion = false;

	bool bClearPaintTexture = false;

	bool bEnableHeightMapBlur = true;
	bool bHeightMapBlurNormals = true;
	uint HeightBlurPassNumIterations = 8;

	bool bEnablePaintTextureBlur = false;
	bool bEnablePaintWetnessBlur = true;
	uint PaintBlurPassNumIterations = 1;

	void RenderGBuffer(RCommandListCompute& CmdList)
	{

		GPU_PROFILE_SCOPE(RenderGBuffer, CmdList);

		auto computeFinishedFence = CmdList.ExecuteCmdList();

		RCommandListGraphics& CmdListGraphics = RCommandList::BeginNew().GetGraphicsContext();
		SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdListGraphics.GetType()).QueueWaitForAnotherFence(computeFinishedFence);

		//Init Resources

		if (GbFirstBoot)
		{
			CGPaintRenderConstantsCBufferGPU.Create(TEXT("CGPaintRenderConstantsCBufferGPU"));
			PaintRenderAllocateResources();
			AllocateGBufferRenderingResources();
		}

		UploadSceneMeshes(CmdListGraphics);

#if PAINT_RENDER_ADVECT_SOURCE
		//if(GbFirstBoot)
		if (bClearPaintTexture || GbFirstBoot)
#endif//PAINT_RENDER_ADVECT_SOURCE
		{
			DrawPaintingBackground(CmdListGraphics);
		}

		RenderShadowMap(CmdListGraphics);


		auto graphicsFinishedFence = CmdListGraphics.ExecuteCmdListAndReleaseContext();
		SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdList.GetType()).QueueWaitForAnotherFence(graphicsFinishedFence);


		//Generate Shadow Visibility Textures per light
		GenerateShadowVisibilityTextures(CmdList);
	}

	void RenderPainting(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(RenderPainting, CmdList);
		
		if (SDFOutputTextureGPU.Size.x != CGSDFConstantsCB.SDFOutputSize.x)
		{
			SDFOutputTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);
			SDFOutputTextureRenderTargetGPU.CreateRenderTarget(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT);
			SDFOutputFallingParticlesTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);
			SDFOutputFallingParticlesTextureRenderTargetGPU.CreateRenderTarget(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, RColor(1.f,1.f,1.f,0.f));
			SDFOutputTexture2GPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);

#if PAINT_RENDER_ADVECT_SOURCE
			SDFOutputTextureCopyGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R11G11B10_FLOAT, true);
#else
			SDFOutputTextureCopyGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);
#endif

			SDFTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT, true);
			SDFTextureRenderTargetGPU.CreateRenderTarget(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT);
			SDFTexture2GPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT, true);
			SDFOutputNormalsTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);
			SDFOutputNormalsTexture2GPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);

			SDFOutputFallingParticlesNormalsTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32G32B32A32_FLOAT, true);

			CGSDFConstantsCBufferGPU.Create(TEXT("CGSDFConstantsCBufferGPU"));


			PaintMaxWetnessTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT, true);
			CmdList.ClearUAV(PaintMaxWetnessTextureGPU);

			WetnessParameterTextureGPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT, true);
			WetnessParameterTexture2GPU.Create(CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_R32_FLOAT, true);

			//SDFOutputDepthStencilBufferGPU.Create(TEXT("SDFOutputDepthStencil"), CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_D24_UNORM_S8_UINT);
			SDFOutputDepthStencilBufferGPU.Create(TEXT("SDFOutputDepthStencil"), CGSDFConstantsCB.SDFOutputSize.x, CGSDFConstantsCB.SDFOutputSize.y, DXGI_FORMAT_D32_FLOAT);

		}

#if !PAINT_RENDER_PIXEL_ADVECT
		//Fill GBuffer
		RenderGBuffer(CmdList);
#endif

		
		if (GbFirstBoot || bClearPaintTexture || bResetSimulation)
		{
#if PAINT_RENDER_ADVECT_SOURCE
			//CopyTextureResourceHelper(CmdList, SDFOutputTextureGPU, PaintingSource);
#else
			CmdList.ClearUAV(SDFOutputTextureGPU, float4{ 1,1,1,0 });
			//CmdList.ClearUAV(SDFOutputTextureGPU);
#endif//PAINT_RENDER_ADVECT_SOURCE
			
		}


	


		auto& RootSig = SDFRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		float3 camPos = RScene::Get().GetCamera().GetPosition();
		CGSDFConstantsCB.CameraPositionWorld = camPos.xy();
		 
		CGSDFConstantsCB.GBufferSize = PaintingBackgroundSizePixels;
		CGSDFConstantsCBufferGPU.Update(CGSDFConstantsCB);
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSDFConstantsCB"), CGSDFConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);
		
		//add velocity to the particles on wet surface
		//diffuse color for particles on wet surface
		{
			GPU_PROFILE_SCOPE(ModifyPaintParticles, CmdList);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFTextureUAV"), SDFTextureGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferUAV"), ParticlesColorsBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferSRV"), CurParticlesVelocityBufferGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesExternalForceBufferUAV"), CurParticlesExternalForceBufferGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

			static ModifyPaintParticlesShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 64;
			uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize) };
			CmdList.Dispatch(numGroupsToDispatch, 1, 1);
		}


		//Modify particles with GBuffer
		{
			GPU_PROFILE_SCOPE(ModifyPaintParticlesWithGBuffer, CmdList);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFTextureUAV"), SDFTextureGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferUAV"), ParticlesColorsBufferGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesExternalForceBufferUAV"), CurParticlesExternalForceBufferGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferUAV"), ParticlesStateBufferGPU, true);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurUAV"), ParticlesBufferPositionCurGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("CurParticlesVelocityBufferUAV"), CurParticlesVelocityBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferNormalTextureSRV"), GBufferNormal, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferWorldSpaceZSRV"), GBufferWorldSpaceZ, true);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferStencilSRV"), GBufferStencil, true);

			static ModifyPaintParticlesWithGBufferShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 64;
			uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGConstantsCB.NumParticles, ThreadGroupSize) };
			CmdList.Dispatch(numGroupsToDispatch, 1, 1);
		}
		


		AllocateGenerateParticleRasterResources(CmdList);

#if RENDER_PARTICLES_SDF == 0
		{
			auto computeFinishedFence = CmdList.ExecuteCmdList();

			RCommandListGraphics& CmdListGraphics = RCommandList::BeginNew().GetGraphicsContext();
			SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdListGraphics.GetType()).QueueWaitForAnotherFence(computeFinishedFence);

			RenderParticlesThroughRaster(CmdListGraphics);

			auto graphicsFinishedFence = CmdListGraphics.ExecuteCmdListAndReleaseContext();
			SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdList.GetType()).QueueWaitForAnotherFence(graphicsFinishedFence);

			ApplyWetnessToMarkedParticles(CmdList);

		}
#else
		{

			//Generate SDF
			{
				GPU_PROFILE_SCOPE(GenerateSDF, CmdList);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFTextureUAV"), SDFTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PaintMaxWetnessTextureUAV"), PaintMaxWetnessTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), SDFOutputTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferSRV"), ParticlesColorsBufferGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

#if PAINT_RENDER_ADVECT_SOURCE
				//CmdList.CopyResource(SDFOutputTextureCopyGPU, SDFOutputTextureGPU);
				CmdList.CopyResource(SDFOutputTextureCopyGPU, GBufferColor);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), GBufferColor, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferMetalTextureSRV"), GBufferMetal, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("OriginalPaintingTextureUAV"), SDFOutputTextureCopyGPU, true);
#endif//PAINT_RENDER_ADVECT_SOURCE

				static GenerateSignedDistanceFieldShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);

#if PAINT_RENDER_ADVECT_SOURCE
				CmdList.CopyResource(GBufferColor, SDFOutputTextureCopyGPU);
#endif//PAINT_RENDER_ADVECT_SOURCE
			}


			//Render Falling Particles
			{
				GPU_PROFILE_SCOPE(GenerateSDF_FallingParticles, CmdList);

				CmdList.SetRootSignature(RootSig);

				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), SDFOutputFallingParticlesTextureGPU, true);
				SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFTextureUAV"), SDFTextureGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesColorsBufferSRV"), ParticlesColorsBufferGPU, true);
				SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);

				static RenderFallingParticlesShader ComputeShader;
				//static GenerateSignedDistanceFieldShader ComputeShader;
				CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

				constexpr uint ThreadGroupSize = 8;
				uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
				CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
			}

		}
#endif //RENDER_PARTICLES_SDF



		RTexture* pHeightMapTexture = &SDFTextureGPU;

		//Blur heightmap
		if(bEnableHeightMapBlur && !bHeightMapBlurNormals)
		{
			GPU_PROFILE_SCOPE(HeightMapBlur, CmdList);
			RTexture* pSource = &SDFTextureGPU;
			RTexture* pDest = &SDFTexture2GPU;
			for (int i = 0; i < HeightBlurPassNumIterations; i++)
			{
				BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT, true);
				std::swap(pSource, pDest);
			}
			
			pHeightMapTexture = pSource;

		}

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		CGSDFConstantsCBufferGPU.Update(CGSDFConstantsCB);
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSDFConstantsCB"), CGSDFConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		//Compute Normals
		{
			GPU_PROFILE_SCOPE(Compute_Normals, CmdList);

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputNormalsTextureUAV"), SDFOutputNormalsTextureGPU, true);

			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFTextureSRV"), *pHeightMapTexture, true);
			
			static GenerateNormalsFromSDFShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		}

		RTexture* pNormalTexture = &SDFOutputNormalsTextureGPU;

		if (bEnableHeightMapBlur && bHeightMapBlurNormals)
		{
			GPU_PROFILE_SCOPE(NormalsTextureBlur, CmdList);
			RTexture* pSource = &SDFOutputNormalsTextureGPU;
			RTexture* pDest = &SDFOutputNormalsTexture2GPU;
			for (int i = 0; i < HeightBlurPassNumIterations; i++)
			{
				BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT3, true);
				std::swap(pSource, pDest);
			}

			pNormalTexture = pSource;

		}

		RTexture* pPaintTexture = &SDFOutputTextureGPU;

		if (bEnablePaintTextureBlur)
		{
			GPU_PROFILE_SCOPE(PaintTextureBlur, CmdList);

#if 1
			if (bEnablePaintWetnessBlur)
			{
				CmdList.CopyResource(SDFOutputTexture2GPU, SDFOutputTextureGPU);

				//Copy wetness parameter into single channel texture
				{
					CmdList.SetRootSignature(RootSig);
					SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("WetnessParameterTextureUAV"), WetnessParameterTextureGPU, true);
					SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureGPU, true);

					static CopyWetnessParameterShader<1> CopyFromComputeShader;
					CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, CopyFromComputeShader));

					constexpr uint ThreadGroupSize = 8;
					uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
					CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
				}
				//Blur
				RTexture* pSource = &WetnessParameterTextureGPU;
				RTexture* pDest = &WetnessParameterTexture2GPU;
				for (int i = 0; i < PaintBlurPassNumIterations; i++)
				{
					BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT, true);
					std::swap(pSource, pDest);
				}
				//Copy back
				{
					CmdList.SetRootSignature(RootSig);
					SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), SDFOutputTexture2GPU, true);
					SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("WetnessParameterTextureSRV"), *pSource, true);

					static CopyWetnessParameterShader<0> CopyToComputeShader;
					CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, CopyToComputeShader));

					constexpr uint ThreadGroupSize = 8;
					uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
					CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
				}

				pPaintTexture = &SDFOutputTexture2GPU;
			}
			else
			{
				//Copy wetness parameter into single channel texture
				{
					CmdList.SetRootSignature(RootSig);
					SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("WetnessParameterTextureUAV"), WetnessParameterTextureGPU, true);
					SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureGPU, true);

					static CopyWetnessParameterShader<1> CopyFromComputeShader;
					CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, CopyFromComputeShader));

					constexpr uint ThreadGroupSize = 8;
					uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
					CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
				}
				BlurPassGPU(CmdList, SDFOutputTextureGPU, SDFOutputTexture2GPU, BlurPassType::TYPE_FLOAT4, true);
				RTexture* pSource = &SDFOutputTexture2GPU;
				RTexture* pDest = &SDFOutputTextureCopyGPU;
				for (int i = 1; i < PaintBlurPassNumIterations; i++)
				{
					BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT4, true);
					std::swap(pSource, pDest);
				}
				

				pPaintTexture = pSource;

				//Copy back
				{
					CmdList.SetRootSignature(RootSig);
					SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), *pPaintTexture, true);
					SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("WetnessParameterTextureSRV"), WetnessParameterTextureGPU, true);

					static CopyWetnessParameterShader<0> CopyToComputeShader;
					CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, CopyToComputeShader));

					constexpr uint ThreadGroupSize = 8;
					uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize) };
					CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
				}
			}
#else
			RTexture* pSource = &SDFOutputFallingParticlesTextureGPU;
			RTexture* pDest = &SDFOutputTextureCopyGPU;
			for (int i = 0; i < PaintBlurPassNumIterations * 2; i++)
			{
				BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT3, true);
				std::swap(pSource, pDest);
			}

#endif

		}


		//Diffuse Paint
#if 1
		if( 1 && bNeablePaintDiffusion)
		{
			GPU_PROFILE_SCOPE(DiffusePaint, CmdList);

			static PaintDiffusionShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint kThreadGroupSize = 8;
			uint2 numGroupsToDispatch{ CGSDFConstantsCB.SDFOutputSize.x / 2 / kThreadGroupSize, CGSDFConstantsCB.SDFOutputSize.y / 2 / kThreadGroupSize };

			PaintDiffusionPerPassConstants PerPassConstants;

			CmdList.CopyResource(SDFOutputTextureCopyGPU, SDFOutputTextureGPU);
			SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureCopyGPU, true);
			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), SDFOutputTextureGPU, true);

			//4 passes per each local id so that neighbors dont collide
			auto ExecutePassFunc = [&](uint2 curLocalCoordOffset)
			{
				PerPassConstants.LocalOffset = (curLocalCoordOffset);

				CmdList.SetConstantArray(RootSig.GetRootParamIndex("PaintDiffusionPerPassConstantsCB"), 2, &PerPassConstants);

				CmdList.TransitionResource(SDFOutputTextureGPU, D3D12_RESOURCE_STATE_UNORDERED_ACCESS);

				CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);
			};

			ExecutePassFunc(uint2{ 0,1 });
			ExecutePassFunc(uint2{ 0,0 });
			ExecutePassFunc(uint2{ 1,0 });
			ExecutePassFunc(uint2{ 1,1 });


		}
#endif

		{
			GPU_PROFILE_SCOPE(VisualizeAdvectedPaintTexture, CmdList);
			//VisualizeVelocityFieldTexture(OriginalPaintingTextureGPU);
			//VisualizeVelocityFieldTexture(SDFOutputTextureGPU);
		}


		//Render
		{
			auto computeFinishedFence = CmdList.ExecuteCmdList();

			RCommandListGraphics& CmdListGraphics = RCommandList::BeginNew().GetGraphicsContext();
			SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdListGraphics.GetType()).QueueWaitForAnotherFence(computeFinishedFence);

			GeneratePaintSceneCameraViewProjection();

			UpdateLights(CmdListGraphics);

			PaintVisualizerPass(CmdListGraphics, pNormalTexture, pPaintTexture);

			auto graphicsFinishedFence = CmdListGraphics.ExecuteCmdListAndReleaseContext();
			SLLRenderer::GetCmdQueueManager().GetQueueBasedOnType(CmdList.GetType()).QueueWaitForAnotherFence(graphicsFinishedFence);
		}

	}






	//=============================================================================================================
	//												PARTICLES RASTERIZER
	//=============================================================================================================


	RTexture ParticleSphereSourceTexture; //R channel for alpha, G channel for Height

	TStructuredBuffer<float3> ParticleRenderVertexBufferGPU;
	TStructuredBuffer<float2> ParticleTexCoordsBufferGPU;

	class GenerateParticleSourceTextureShader : public RShader
	{
	public:
		GenerateParticleSourceTextureShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "GenerateParticleSourceTexture", L"cs_6_0");
		}
	};


	class ParticlesRenderVertexShader : public RShader
	{
	public:
		ParticlesRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "ParticleRenderVS", L"vs_6_0");
		}
	};

	class StickyParticlesRenderPixelShader : public RShader
	{
	public:
		StickyParticlesRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "StickyParticleRenderPS", L"ps_6_0");
		}
	};

	class FallingParticlesRenderPixelShader : public RShader
	{
	public:
		FallingParticlesRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "FallingParticleRenderPS", L"ps_6_0");
		}
	};

	class ApplyWetnessToMarkedParticlesShader : public RShader
	{
	public:
		ApplyWetnessToMarkedParticlesShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "ApplyWetnessToMarkedParticles", L"cs_6_0");
		}
	};

	void AllocateGenerateParticleRasterResources(RCommandListCompute& CmdList)
	{
		if (GbFirstBoot)
		{
			ParticleSphereSourceTexture.Create(CGSDFConstantsCB.ParticleSourceTextureSize, CGSDFConstantsCB.ParticleSourceTextureSize, DXGI_FORMAT_R32G32_FLOAT, true);

			//Sprite Geometry
			{
				RDynamicVector<float3> vertexDataCpu;
				RDynamicVector<float2> texCoordDataCpu;

				float3 center = float3{ 0,0,0 };
				float2 extent = float2{ 1, 1 };

				float3 p0{ center - float3{extent, 0} };
				float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
				float3 p2{ center + float3{extent, 0} };
				float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

				vertexDataCpu.push_back(p0);
				vertexDataCpu.push_back(p1);
				vertexDataCpu.push_back(p2);

				vertexDataCpu.push_back(p2);
				vertexDataCpu.push_back(p3);
				vertexDataCpu.push_back(p0);

				const float2 leftDown = float2{ 0.f, 1.f };
				const float2 leftUp = float2{ 0.f, 0.f };
				const float2 rightUp = float2{ 1.f, 0.f };
				const float2 rightDown = float2{ 1.f, 1.f };

				texCoordDataCpu.push_back(leftDown);
				texCoordDataCpu.push_back(leftUp);
				texCoordDataCpu.push_back(rightUp);

				texCoordDataCpu.push_back(rightUp);
				texCoordDataCpu.push_back(rightDown);
				texCoordDataCpu.push_back(leftDown);

				ParticleRenderVertexBufferGPU.Create(TEXT("PaintVertexBuffer"), 6, vertexDataCpu.data());
				ParticleTexCoordsBufferGPU.Create(TEXT("PaintTexCoordsBuffer"), 6, texCoordDataCpu.data());
			}
		}
		
		{
			GPU_PROFILE_SCOPE(GenerateParticleSourceTexture, CmdList);

			auto& RootSig = SDFRenderingRootSig::Get();

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ParticleSourceTextureUAV"), ParticleSphereSourceTexture, true);

			static GenerateParticleSourceTextureShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			uint numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.ParticleSourceTextureSize, ThreadGroupSize) };
			CmdList.Dispatch(numGroupsToDispatch, numGroupsToDispatch, 1);
		}


	}

	void ApplyWetnessToMarkedParticles(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(ApplyWetnessToMarkedParticles, CmdList);

		auto& RootSig = SDFRenderingRootSig::Get();

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());

		CGSDFConstantsCBufferGPU.Update(CGSDFConstantsCB);
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGSDFConstantsCB"), CGSDFConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFTextureUAV"), SDFTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFTextureSRV"), SDFTextureRenderTargetGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureUAV"), SDFOutputTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("SDFOutputTextureSRV"), SDFOutputTextureRenderTargetGPU, true);
		//SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("WetnessParameterTextureUAV"), WetnessParameterTextureGPU, true);
		SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("PaintMaxWetnessTextureUAV"), PaintMaxWetnessTextureGPU, true);

		static ApplyWetnessToMarkedParticlesShader ComputeShader;
		CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

		constexpr uint ThreadGroupSize = 8;
		uint2 numGroupsToDispatch{ RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.x, ThreadGroupSize), RMath::DivideAndRoundUp(CGSDFConstantsCB.SDFOutputSize.y, ThreadGroupSize) };
		CmdList.Dispatch(numGroupsToDispatch.x, numGroupsToDispatch.y, 1);
	}

	void RenderParticlesThroughRaster(RCommandListGraphics& CmdList)
	{
		GPU_PROFILE_SCOPE(RenderParticlesThroughRaster, CmdList);

		auto& RootSignature = SDFRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSignature);

		//define Viewport
		static CD3DX12_VIEWPORT Viewport;
		static CD3DX12_RECT ScissorRect;

		const uint2 viewportSize = SDFOutputTextureRenderTargetGPU.Size;

		Viewport = CD3DX12_VIEWPORT{ 0.0f, 0.0f, (FLOAT)viewportSize.x, (FLOAT)viewportSize.y };
		ScissorRect = CD3DX12_RECT{ 0, 0, (LONG)viewportSize.x, (LONG)viewportSize.y };

		//Set global states
		CmdList.SetViewport(Viewport);
		CmdList.SetScissor(ScissorRect);

		CmdList.TransitionResource(SDFOutputTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET, true);
		CmdList.TransitionResource(SDFTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET, true);
		CmdList.TransitionResource(SDFOutputFallingParticlesTextureRenderTargetGPU, D3D12_RESOURCE_STATE_RENDER_TARGET, true);
		CmdList.TransitionResource(DepthBuffer, D3D12_RESOURCE_STATE_DEPTH_WRITE, true);

		CmdList.ClearRenderTarget(SDFOutputFallingParticlesTextureRenderTargetGPU);
		CmdList.ClearRenderTarget(SDFTextureRenderTargetGPU);
		//if (bClearPaintTexture)
		{
			CmdList.ClearRenderTarget(SDFOutputTextureRenderTargetGPU);
		}

		CmdList.ClearDepthAndStencil(SDFOutputDepthStencilBufferGPU);

		//Set Descriptors

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticlesColorsBufferSRV"), ParticlesColorsBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticlesStateBufferSRV"), ParticlesStateBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU, true);

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());
		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("CGPaintRenderConstantsCB"), CGPaintRenderConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticleRenderVertexBufferSRV"), ParticleRenderVertexBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticleTexCoordsBufferSRV"), ParticleTexCoordsBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ParticleSourceTextureSRV"), ParticleSphereSourceTexture, true);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("AdvectedParticleDensityTextureSRV"), AdvectedParticleDensityTextureGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("AverageDensityBufferSRV"), AverageDensityBufferGPU, true);

		//PSO 
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;

		//DepthStencil
		{
			/*psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateStencilWrite;
			psoStates.DepthStencilDesc.StencilWriteMask = D3D12_DEFAULT_STENCIL_WRITE_MASK;
			psoStates.DepthStencilDesc.FrontFace.StencilPassOp = D3D12_STENCIL_OP_INCR;
			psoStates.DepthStencilDesc.FrontFace.StencilFailOp = D3D12_STENCIL_OP_INCR;
			psoStates.DepthStencilDesc.BackFace.StencilPassOp = D3D12_STENCIL_OP_INCR;
			psoStates.DepthStencilDesc.BackFace.StencilFailOp = D3D12_STENCIL_OP_INCR;*/

			//psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateWriteOnly;
			psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;
		}



		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		//psoStates.BlendDesc = RCommonResources::Get().BlendTraditional;
		//psoStates.BlendDesc = RCommonResources::Get().BlendDefault;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;


		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);


		static constexpr uint NumRenderTargets = 2;


		//Render Sticky Particles
		{
			D3D12_CPU_DESCRIPTOR_HANDLE rtvHandles[NumRenderTargets];
			rtvHandles[0] = SDFOutputTextureRenderTargetGPU.GetRenderTargetView();
			rtvHandles[1] = SDFTextureRenderTargetGPU.GetRenderTargetView();

			DXGI_FORMAT rtFormats[NumRenderTargets];
			rtFormats[0] = SDFOutputTextureRenderTargetGPU.Format;
			rtFormats[1] = SDFTextureRenderTargetGPU.Format;

			psoStates.NumRenderTargets = NumRenderTargets;
			psoStates.pRenderTargetFormat = rtFormats;

			CmdList.SetRenderTargets(NumRenderTargets, rtvHandles, SDFOutputDepthStencilBufferGPU.GetDSV());

			static ParticlesRenderVertexShader VS;
			static StickyParticlesRenderPixelShader PS;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

			CmdList.SetPipelineState(PSO);

			CmdList.DrawInstanced(6, CGConstantsCB.NumParticles);
		}

		//Render Falling Particles
		{
			D3D12_CPU_DESCRIPTOR_HANDLE rtvHandles[NumRenderTargets];
			rtvHandles[0] = SDFOutputFallingParticlesTextureRenderTargetGPU.GetRenderTargetView();
			rtvHandles[1] = SDFTextureRenderTargetGPU.GetRenderTargetView();

			DXGI_FORMAT rtFormats[NumRenderTargets];
			rtFormats[0] = SDFOutputTextureRenderTargetGPU.Format;
			rtFormats[1] = SDFTextureRenderTargetGPU.Format;
			//rtFormats[2] = SDFOutputThisFrameTextureRenderTargetGPU.Format;

			psoStates.NumRenderTargets = NumRenderTargets;
			psoStates.pRenderTargetFormat = rtFormats;

			CmdList.SetRenderTargets(NumRenderTargets, rtvHandles, SDFOutputDepthStencilBufferGPU.GetDSV());

			static ParticlesRenderVertexShader VS;
			static FallingParticlesRenderPixelShader PS;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

			CmdList.SetPipelineState(PSO);

			CmdList.DrawInstanced(6, CGConstantsCB.NumParticles);
		}

		//CmdList.CopyResource(SDFOutputTextureGPU, SDFOutputTextureRenderTargetGPU);
		CmdList.CopyResource(SDFOutputFallingParticlesTextureGPU, SDFOutputFallingParticlesTextureRenderTargetGPU);
		//CmdList.CopyResource(SDFTextureGPU, SDFTextureRenderTargetGPU);

	}
















	//=============================================================================================================
	//												PAINT SHADING PASS
	//=============================================================================================================


	class PaintRenderingRootSig : public RRootSignatureExt, public RSingleton<PaintRenderingRootSig>
	{
	public:
		PaintRenderingRootSig()
		{
			///Constants
			AddConstantBuffer("CGConstantsCB", 0);
			AddConstantBuffer("CGPaintRenderConstantsCB", 10);
			AddConstantBuffer("CGSDFConstantsCB", 1);

			AddConstantBuffer("SceneUniformBufferCBuffer", 3);

			AddConstants("UtilityConstants", 2, 1);

			///SRVs
			AddSRV("VertexBufferPosition", 0);
			AddSRV("TexCoordBuffer", 1);
			AddSRV("ColorTexture", 2);
			AddSRV("NormalTexture", 3);
			AddSRV("PaintingTexture", 4);
			AddSRV("LightsBuffer", 5);
			AddSRV("FallingParticlesTexture", 6);

			AddSRV("ShadowVisibilityTextureUAV", 7);


			AddSRV("IndexBuffer", 7);
			AddSRV("VertexBufferRest", 8);

			AddSRV("RoughnessTexture", 9);

			AddSRV("GBufferColorTexture", 10);
			AddSRV("GBufferNormalTexture", 11);
			AddSRV("GBufferRoughnessTexture", 12);

			AddSRV("GBufferWorldSpaceZSRV", 13);

			AddSRV("ShadowMapTextureSRV", 14);

			AddSRV("LightSpaceTransformMatricesPerLightBufferSRV", 15);
			//AddSRV("ShadowMapTextureSRV", 16);

			AddSRV("MetalnessTexture", 16);
			AddSRV("GBufferMetalTexture", 17);

			AddSRV("ShadowMapPerLightTexturesSRV0", 18);
			AddSRV("ShadowMapPerLightTexturesSRV1", 19);
			AddSRV("ShadowMapPerLightTexturesSRV2", 20);
			AddSRV("ShadowMapPerLightTexturesSRV3", 21);
			


			///UAVs


			//Samplers
			//AddStaticSampler(0, RCommonResources::Get().SamplerLinearClampDesc);
			AddStaticSampler(0, RCommonResources::Get().SamplerLinearRepeatDesc);
			AddStaticSampler(1, RCommonResources::Get().SamplerPointClampDesc);
			AddStaticSampler(2, RCommonResources::Get().SamplerLinearBorderBlackDesc);
			//AddStaticSampler(3, RCommonResources::Get().SamplerLinearRepeatDesc);

			Finalize(L"PaintRenderingRootSig");
		}
	};


	TStructuredBuffer<float3> PaintRenderVertexBufferGPU;
	TStructuredBuffer<float2> PaintTexCoordsBufferGPU;

	///CBuffer Constants
	PaintRenderConstants CGPaintRenderConstantsCB;
	TCBuffer<PaintRenderConstants> CGPaintRenderConstantsCBufferGPU;

	//Lights
	TStructuredBuffer<LightDescription> LightsBufferGPU;
	RDynamicVector<LightDescription> LightsBufferCPU;
	RDynamicVector<uint> LightIndexInSpatialControlSet;
	uint MaxNumLights = 100.f;

	LightDescription& AddLight(c_float3& pos, bool bSpot = false)
	{
		check(LightsBufferCPU.size() < MaxNumLights);
		uint lightIndex = RSpatialControlSet::Get().ControlPoints.size();
		auto& newLight = LightsBufferCPU.emplace_back();
		newLight.Position = pos;
		newLight.Type = bSpot ? LIGHT_TYPE_SPOT : LIGHT_TYPE_POINT;
		RSpatialControlSet::Get().AddPoint(&newLight.Position);
		LightIndexInSpatialControlSet.push_back(lightIndex);
		return newLight;
	}

	void UpdateLights(RCommandListGraphics& CmdList)
	{
		if (LightsBufferCPU.empty())
		{
			LightsBufferCPU.reserve(MaxNumLights);
			AddLight(float3{ -1.3,12.4,-20.3 });
			AddLight(float3{ 23.3,17.8,-10.5 });
			AddLight(float3{ 23.3,0.,-5.5 });

			auto& spotLight = AddLight(float3{ 11,22.4,-2.8 }, true);
			spotLight.Radius = 25;
			spotLight.SpotlightAngles = { 0.7f,0 }; //x = cos(Inner); y = cos(Outer);

			LightsBufferGPU.Create(TEXT("LightBuffer"), MaxNumLights);
		}

		//Update directions for each light
		for(auto& light: LightsBufferCPU)
		{
			//each light always points at the painting
			light.Direction = VectorNormalize(GetPaintSimulationCenterPos() - light.Position);
		}

		CmdList.UploadToBuffer(LightsBufferGPU, 0, LightsBufferCPU.data(), LightsBufferCPU.size() * sizeof(LightDescription));
		CGPaintRenderConstantsCB.NumLights = LightsBufferCPU.size();

	}

	void PaintRenderAllocateResources()
	{
		if (PaintRenderVertexBufferGPU.GetElementCount() == 0)
		{
			RDynamicVector<float3> vertexDataCpu;
			RDynamicVector<float2> texCoordDataCpu;
			
			float half = GetPaintSimulationLength() * 0.5f;

			float3 center = GetPaintSimulationCenterPos();
			float2 extent = float2{ half, half };

			float3 p0{ center - float3{extent, 0} };
			float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
			float3 p2{ center + float3{extent, 0} };
			float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

			vertexDataCpu.push_back(p0);
			vertexDataCpu.push_back(p1);
			vertexDataCpu.push_back(p2);

			vertexDataCpu.push_back(p2);
			vertexDataCpu.push_back(p3);
			vertexDataCpu.push_back(p0);

			const float2 leftDown = float2{ 0.f, 1.f };
			const float2 leftUp = float2{ 0.f, 0.f };
			const float2 rightUp = float2{ 1.f, 0.f };
			const float2 rightDown = float2{ 1.f, 1.f };

			texCoordDataCpu.push_back(leftDown);
			texCoordDataCpu.push_back(leftUp);
			texCoordDataCpu.push_back(rightUp);

			texCoordDataCpu.push_back(rightUp);
			texCoordDataCpu.push_back(rightDown);
			texCoordDataCpu.push_back(leftDown);

			PaintRenderVertexBufferGPU.Create(TEXT("PaintVertexBuffer"), 6, vertexDataCpu.data());
			PaintTexCoordsBufferGPU.Create(TEXT("PaintTexCoordsBuffer"), 6, texCoordDataCpu.data());

		}
	}

	void PaintVisualizerDrawUI()
	{

		if (ImGui::Begin("Painting Render"))
		{

			if (ImGui::BeginTabBar("PhysTabBar"))
			{

				//Scene
				if (ImGui::BeginTabItem("Scene"))
				{
					ImGui::DragFloat3("View Position", &PaintSceneCameraPosition.x);
					if (ImGui::DragFloat("View Scale", &PaintSceneCameraScale.x))
					{
						PaintSceneCameraScale.y = PaintSceneCameraScale.x;
					}

					bReUploadScene = false;

					bReUploadScene |= ImGui::DragFloat2("Frame Mesh Scale", &MeshScale.x, 0.01, 0.1, 10.f);
					bReUploadScene |= ImGui::DragFloat("Frame Mesh Orientation Offset", &MeshOrientationOffset);

					bReUploadScene |= ImGui::DragFloat3("Paint Source Pos Offset", &PaintSourcePositionOffset.x, 0.01, -10.f, 10.f);
					bReUploadScene |= ImGui::DragFloat2("Paint Source Scale", &PaintSourceScale.x, 0.01, 0.1, 10.f);

					ImGui::EndTabItem();
				}

				//Shading
				if (ImGui::BeginTabItem("Shading"))
				{

					ImGui::DragFloat("Paint Brightness Scale", &CGPaintRenderConstantsCB.PaintBrightnessScale, 0.01f, 0.01f, 10.f);

					ImGui::SliderFloat("Uniform Roughness", &CGPaintRenderConstantsCB.UniformRoughness, 0.0f, 1.f);
					ImGui::SliderFloat("Uniform Metalness", &CGPaintRenderConstantsCB.UniformMetalness, 0.0f, 1.f);
					ImGui::SliderFloat2("Roughness Min Max", &CGPaintRenderConstantsCB.RoughnessMinMax.x, 0.0f, 1.f);
					ImGui::SliderFloat2("Metalness Min Max", &CGPaintRenderConstantsCB.MetalnessMinMax.x, 0.0f, 1.f);

					ImGui::Checkbox("Shade Lambert", (bool*)&CGPaintRenderConstantsCB.bShadeLambert);
					ImGui::Checkbox("Shade Specular", (bool*)&CGPaintRenderConstantsCB.bShadeSpecular);
					ImGui::Checkbox("Shade Fresnel", (bool*)&CGPaintRenderConstantsCB.bShadeFresnel);
					ImGui::SliderInt("Specular Shading Method", (int*)&CGPaintRenderConstantsCB.SpecularShadingMethod, 0, 2);

					ImGui::DragFloat("SpecularPower", &CGPaintRenderConstantsCB.SpecularPower, 0.5f, 0.01f, 1000.f);
					ImGui::DragFloat("SpecularAmount", &CGPaintRenderConstantsCB.SpecularAmount, 0.01f, 0.0f, 10.f);

					ImGui::SliderFloat2("Global Light Direction Spherical", &GlobalLightDirectionSphericalCoords.x, -Pi, Pi);
					ImGui::DragFloat3("Global Light Pos", &GlobalLightPosition.x);
					ImGui::Checkbox("Global Light From Cam", &bAssignGlobalLightPosFromCamera);
					ImGui::DragFloat("ShadowMap Scale", &ShadowMapScale.x, 0.01, 1.f, 10.f);
					ImGui::DragFloat("ShadowMap Depth Bias", &CGPaintRenderConstantsCB.ShadowDepthBias, 0.001, 0.f, 1.f);
					ImGui::Checkbox("Slope Correct Bias", (bool*)&CGPaintRenderConstantsCB.bSlopeCorrectShadowBias);

					ImGui::DragFloat("Amgient Light Intensity", &CGPaintRenderConstantsCB.AmbientLight, 0.01f, 0.0f, 10.f);

					ImGui::Checkbox("Debug Visualize Normal Map", (bool*)&CGPaintRenderConstantsCB.bVisualizeNormals);

					ImGui::EndTabItem();
				}

				//Lights
				if (ImGui::BeginTabItem("Lights"))
				{
					{
						auto size = LightsBufferCPU.size();
						for (int i = 0; i < size; i++)
						{
							auto curLightIndexInSet = LightIndexInSpatialControlSet[i];

							//if cur light is selected
							if (&RSpatialControlSet::Get().ControlPoints[curLightIndexInSet] == RSpatialControlSet::Get().pSelectedPoint)
							{
								LightDescription& curLightDesc = LightsBufferCPU[i];

								ImGui::NewLine();
								ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Light");
								ImGui::DragFloat3("LightPosition", &curLightDesc.Position.x);
								ImGui::DragFloat("Light Radius", &curLightDesc.Radius, 0.01, 0.01, 1000);
								if (curLightDesc.Type == LIGHT_TYPE_SPOT)
								{
									float2 angles;
									angles.x = std::acosf(curLightDesc.SpotlightAngles.x);
									angles.y = std::acosf(curLightDesc.SpotlightAngles.y);
									if (ImGui::SliderAngle("Spotlight Inner Angle", &angles.x, 0, 90))
									{
										curLightDesc.SpotlightAngles.x = std::cos(angles.x);
									}
									if (ImGui::SliderAngle("Spotlight Outer Angle", &angles.y, 0, 90))
									{
										curLightDesc.SpotlightAngles.y = std::cos(angles.y);
									}
									if (ImGui::DragFloat3("Spotlight Direction", &curLightDesc.Direction.x, 0.1))
									{
										curLightDesc.Direction = VectorNormalize(curLightDesc.Direction);
									}
								}
								ImGui::DragFloat("Light Intensity", &curLightDesc.Intensity, 0.01, 0.0, 100);
								ImGui::ColorEdit3("Light Color", &curLightDesc.Color.x);

								//Shadow Map
								if (ImGui::Begin("ShadowMap"))
								{
									auto windowSize = ImGui::GetWindowSize();
									ImGui::Image((ImTextureID)DepthBufferDescHandleGPUArr[i].GPU.ptr, ImVec2((float)windowSize.x, (float)windowSize.y));
								}
								ImGui::End();

								if (ImGui::Begin("ShadowVisibilityMap"))
								{
									auto windowSize = ImGui::GetWindowSize();
									ImGui::Image((ImTextureID)ShadowVisibilityDescHandleGPUArr[i].GPU.ptr, ImVec2((float)windowSize.x, (float)windowSize.y));
								}
								ImGui::End();

							}
						}
					}
					ImGui::EndTabItem();
				}


			}
			ImGui::EndTabBar();

		}
		ImGui::End();

	}

	
	class PaintRenderVertexShader : public RShader
	{
	public:
		PaintRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "PaintRenderVS", L"vs_6_0");
		}
	};

	class PaintRenderPixelShader : public RShader
	{
	public:
		PaintRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "PaintRenderPS", L"ps_6_0");
		}
	};

	float2 GlobalLightDirectionSphericalCoords{ PiDiv2,-PiDiv2/2.f };

	bool bDebugVisualizeGlobalLightDirection = false;
	
	void PaintVisualizerPass(RCommandListGraphics& CmdList, RTexture* pNormalsTexture, RTexture* pPaintTexture)
	{
		GPU_PROFILE_SCOPE(PaintVisualizerPass, CmdList);

		PaintVisualizerDrawUI();

		auto& RootSignature = PaintRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSignature);


		//Set global states 
		CmdList.SetViewport(RSceneRenderer::Get().GetViewport());
		CmdList.SetScissor(RSceneRenderer::Get().GetScissorRect());

		auto& RT = RCommonResources::Get().SceneColor;
		auto& DSV = RCommonResources::Get().SceneDepth;

		CmdList.SetRenderTarget(RT.GetRTV(), DSV.GetDSV());

		//Set Constants

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("SceneUniformBufferCBuffer"), RScene::Get().UniformBufferParameters.GetGlobalUniformBuffer()->ConstantBuffer.GetResource()->GetGPUVirtualAddress());


		//Update values
		if (bAssignGlobalLightPosFromCamera)
		{
			CGPaintRenderConstantsCB.GlobalLightDirection = RScene::Get().GetCamera().GetDirection();
		}
		else
		{
			//CGPaintRenderConstantsCB.GlobalLightDirection = MapSphericalToCartesian(GlobalLightDirectionSphericalCoords.x, GlobalLightDirectionSphericalCoords.y, 1);
		}
		

		CGPaintRenderConstantsCBufferGPU.Update(CGPaintRenderConstantsCB);
		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("CGPaintRenderConstantsCB"), CGPaintRenderConstantsCBufferGPU.GetGpuVirtualAddress());

		//Set Resources
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), PaintRenderVertexBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("TexCoordBuffer"), PaintTexCoordsBufferGPU, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), *pPaintTexture, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), *pNormalsTexture, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("PaintingTexture"), PaintingSource, true);

		//SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ShadowMapTextureSRV"), ShadowMap, true);
		// 
		//CmdList.TransitionResource(GlobalLightDepth, D3D12_RESOURCE_STATE_GENERIC_READ);
		//CmdList.SetDynamicDescriptorInRange(RRootSignatureExt::GetSRVRootIndex(), RootSignature.GetDescriptorIndex("ShadowMapTextureSRV"), GlobalLightDepth.GetDepthSRV());

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("GBufferColorTexture"), GBufferColor, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("GBufferNormalTexture"), GBufferNormal, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("GBufferRoughnessTexture"), GBufferRoughness, true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("GBufferMetalTexture"), GBufferMetal, true);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("FallingParticlesTexture"), SDFOutputFallingParticlesTextureGPU, true);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("LightsBuffer"), LightsBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("GBufferWorldSpaceZSRV"), GBufferWorldSpaceZ, true);

		check((ShadowVisibilityPerLightTextureGPUArr.size() == LightsBufferCPU.size()) && (LightsBufferCPU.size() == 4));

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ShadowMapPerLightTexturesSRV0"), ShadowVisibilityPerLightTextureGPUArr[0], true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ShadowMapPerLightTexturesSRV1"), ShadowVisibilityPerLightTextureGPUArr[1], true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ShadowMapPerLightTexturesSRV2"), ShadowVisibilityPerLightTextureGPUArr[2], true);
		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ShadowMapPerLightTexturesSRV3"), ShadowVisibilityPerLightTextureGPUArr[3], true);
		
		
		//PSO 
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;
		//psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;
		psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;

		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

		///Shaders
		static PaintRenderVertexShader VS;
		static PaintRenderPixelShader PS;

		auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

		///We can destroy created PSO since D3D12 PSO is stored in map
		CmdList.SetPipelineState(PSO);

		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);


		CmdList.Draw(6);


		if (bDebugVisualizeGlobalLightDirection)
		{
			RGeometry::GInstance.SetColor(EColor::YellowCanary);

			float half = (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) * 0.5f;
			float centerf = CGConstantsCB.VelocityFieldBoundsStartOffset + half;

			float3 center = float3{ centerf, centerf, CGConstantsCB.VelocityFieldBoundsStartOffset };
			float2 extent = float2{ half, half };

			float3 start = center - CGPaintRenderConstantsCB.GlobalLightDirection * 5.f;

			RGeometry::GInstance.AddArrowPolyboard(start, center);
		}
	}















	//=============================================================================================================
	//												BACKGROUND RENDERING
	//=============================================================================================================


	RTexture GBufferColor{TEXT("GBufferColor")};
	RTexture GBufferNormal{ TEXT("GBufferNormal") }; //TODO: USE R32G32 format, we don't need .z ----------------!!!!!!!!!!!!!!!!!!!!!!!!!
	RTexture GBufferRoughness{ TEXT("GBufferRoughness") };
	RTexture GBufferMetal{ TEXT("GBufferMetal") };

	RTexture GBufferWorldSpaceZ{ TEXT("GBufferWorldSpaceZ") };
	RDepthBuffer DepthBuffer{ 1.f };

	RTexture GBufferStencil{ TEXT("GBufferStencil" )};

	float3 PaintSceneCameraPosition{ 10,10, 0 };
	float2 PaintSceneCameraScale{ 1,1 };

	void GeneratePaintSceneCameraViewProjection()
	{
		if (GbFirstBoot)
		{
			float half = (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) * 0.5f;
			float centerf = CGConstantsCB.VelocityFieldBoundsStartOffset + half;
			PaintSceneCameraPosition = float3{ centerf, centerf, 0 };
			PaintSceneCameraScale = float2::MakeFloat2(half * 2.f);
		}
		CGPaintRenderConstantsCB.PaintSceneViewMatrix = DirectX::XMMatrixLookToLH(float3SSE(PaintSceneCameraPosition), float3SSE{ 0,0,1 }, float3SSE{0,1,0});
		CGPaintRenderConstantsCB.PaintSceneProjectionMatrix = DirectX::XMMatrixOrthographicLH(PaintSceneCameraScale.x, PaintSceneCameraScale.y, 0.01, 20);
		CGPaintRenderConstantsCB.PaintSceneViewProjectionMatrix = DirectX::XMMatrixMultiply(CGPaintRenderConstantsCB.PaintSceneViewMatrix, CGPaintRenderConstantsCB.PaintSceneProjectionMatrix);
	}

	class BackgroundMeshRenderVertexShader : public RShader
	{
	public:
		BackgroundMeshRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "BackgroundMeshRenderVS", L"vs_6_0");
		}
	};

	class BackgroundWallRenderVertexShader : public RShader
	{
	public:
		BackgroundWallRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "BackgroundWallRenderVS", L"vs_6_0");
		}
	};

	class BackgroundFrameRenderPixelShader : public RShader
	{
	public:
		BackgroundFrameRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "BackgroundFrameRenderPS", L"ps_6_0");
		}
	};

	class BackgroundWallRenderPixelShader : public RShader
	{
	public:
		BackgroundWallRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "BackgroundWallRenderPS", L"ps_6_0");
		}
	};
	class PaintSourceRenderPixelShader : public RShader
	{
	public:
		PaintSourceRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "PaintSourceRenderPS", L"ps_6_0");
		}
	};


	class BackgroundDepthRenderVertexShader : public RShader
	{
	public:
		BackgroundDepthRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "DepthRenderVS", L"vs_6_0");
		}
	};
	class BackgroundDepthIndexedRenderVertexShader : public RShader
	{
	public:
		BackgroundDepthIndexedRenderVertexShader()
		{
			CreateVS("Physics/ParticlePhysics.hlsl", "DepthRenderIndexedVS", L"vs_6_0");
		}
	};

	class BackgroundDepthRenderPixelShader : public RShader
	{
	public:
		BackgroundDepthRenderPixelShader()
		{
			CreatePS("Physics/ParticlePhysics.hlsl", "DepthRenderPS", L"ps_6_0");
		}
	};


	//Render Mesh

	RStructuredBuffer MeshVertexBuffer;
	RStructuredBuffer MeshIndexBuffer;
	RStructuredBuffer MeshVertexBufferRest;

	RMeshDesc MeshDesc;
	RVertexBufferCPU FrameMeshVertexBufferCPU;

	RTexture MeshTextureAlbedo;
	RTexture MeshTextureRoughness;
	RTexture MeshTextureNormal;
	RTexture MeshTextureMetalness;

	bool bReUploadScene = true;

	float MeshOrientationOffset{ 0 };
	float2 MeshScale{ 1.19,1.83 };

	RTexture WallTextureAlbedo;
	RTexture WallTextureRoughness;
	RTexture WallTextureNormal;

	RTexture AdditionalTextureAlbedo;
	RTexture AdditionalTextureRoughness;
	RTexture AdditionalTextureNormal;
	RTexture AdditionalTextureMetalness;

	float2 PaintSourceScale{ 0.5,0.5 };
	float3 PaintSourcePositionOffset{ 0,0,-0.35f };
	float PaintSourceOrientationOffset{ 0 };

	float OriginalPaintSourceXYRatio{ 1 };

	TStructuredBuffer<float3> PaintSourceVertexBufferGPU;


	uint2 PaintingBackgroundSizePixels{ 1024,1024 };

	/*struct VertexInputRest
	{
		float3 Normal;
		float3 Tangent;
		float3 Bitangent;
		float2 TexCoord;
	};*/

	float3 GetPaintSimulationCenterPos()
	{
		float half = GetPaintSimulationLength() * 0.5f;
		float centerf = CGConstantsCB.VelocityFieldBoundsStartOffset + half;
		return float3{ centerf, centerf, CGConstantsCB.DefaultParticlePos.z };
	}

	float GetPaintSimulationLength()
	{
		return (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength);
	}


	void RescaleVertexBuffer(float3 scale, RVertexBufferCPU& vb)
	{
		for (int v = 0; v < vb.NumVertices; v++)
		{
			float3& vertex = vb.TGetVertexRef<float3>(v);

			float3SSE vertexXM = vertex;

			vertexXM = DirectX::XMVector3Transform(vertexXM, DirectX::XMMatrixScaling(scale.x, scale.y, scale.z));

			vertex = vertexXM;
		}
	}

	void RotateVertexBuffer(float3 pitchYawRoll, RVertexBufferCPU& vb)
	{
		for (int v = 0; v < vb.NumVertices; v++)
		{
			float3& vertex = vb.TGetVertexRef<float3>(v);

			float3SSE vertexXM = vertex;

			vertexXM = DirectX::XMVector3Transform(vertexXM, DirectX::XMMatrixRotationRollPitchYaw(pitchYawRoll.x, pitchYawRoll.y, pitchYawRoll.z));

			vertex = vertexXM;
		}
	}

	void TranslateVertexBuffer(float3 translateVec, RVertexBufferCPU& vb)
	{
		for (int v = 0; v < vb.NumVertices; v++)
		{
			float3& vertex = vb.TGetVertexRef<float3>(v);

			float3SSE vertexXM = vertex;

			vertexXM = DirectX::XMVector3Transform(vertexXM, DirectX::XMMatrixTranslation(translateVec.x, translateVec.y, translateVec.z));

			vertex = vertexXM;
		}
	}


	void UploadSceneMeshes(RCommandList& CmdList)
	{
		if (GbFirstBoot)
		{
			RDynamicVector<float3> vertexDataCpu;
			RDynamicVector<float2> texCoordDataCpu;

			float half = GetPaintSimulationLength() * 0.5f;

			float3 center = GetPaintSimulationCenterPos();
			float2 extent = float2{ half, half };

			center = center + PaintSourcePositionOffset;
			extent.x *= OriginalPaintSourceXYRatio;
			extent.x *= PaintSourceScale.x;
			extent.y *= PaintSourceScale.y;

			float3 p0{ center - float3{extent, 0} };
			float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
			float3 p2{ center + float3{extent, 0} };
			float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

			vertexDataCpu.push_back(p0);
			vertexDataCpu.push_back(p1);
			vertexDataCpu.push_back(p2);

			vertexDataCpu.push_back(p2);
			vertexDataCpu.push_back(p3);
			vertexDataCpu.push_back(p0);

			const float2 leftDown = float2{ 0.f, 1.f };
			const float2 leftUp = float2{ 0.f, 0.f };
			const float2 rightUp = float2{ 1.f, 0.f };
			const float2 rightDown = float2{ 1.f, 1.f };

			texCoordDataCpu.push_back(leftDown);
			texCoordDataCpu.push_back(leftUp);
			texCoordDataCpu.push_back(rightUp);

			texCoordDataCpu.push_back(rightUp);
			texCoordDataCpu.push_back(rightDown);
			texCoordDataCpu.push_back(leftDown);

			PaintSourceVertexBufferGPU.Create(TEXT("PaintVertexBuffer"), 6, vertexDataCpu.data());
			PaintTexCoordsBufferGPU.Create(TEXT("PaintTexCoordsBuffer"), 6, texCoordDataCpu.data());
		}

		//Update Vertices
		if(bReUploadScene)
		{
			RDynamicVector<float3> vertexDataCpu;

			float half = GetPaintSimulationLength() * 0.5f;

			float3 center = GetPaintSimulationCenterPos();
			float2 extent = float2{ half, half };

			center = center + PaintSourcePositionOffset;
			extent.x *= OriginalPaintSourceXYRatio;
			extent.x *= PaintSourceScale.x;
			extent.y *= PaintSourceScale.y;

			float3 p0{ center - float3{extent, 0} };
			float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
			float3 p2{ center + float3{extent, 0} };
			float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

			vertexDataCpu.push_back(p0);
			vertexDataCpu.push_back(p1);
			vertexDataCpu.push_back(p2);

			vertexDataCpu.push_back(p2);
			vertexDataCpu.push_back(p3);
			vertexDataCpu.push_back(p0);

			CmdList.UploadToBuffer(PaintSourceVertexBufferGPU, 0, vertexDataCpu.data(), vertexDataCpu.size() * sizeof(float3));
		}

		if (FrameMeshVertexBufferCPU.NumVertices == 0)
		{
			RInstancedStaticMesh FrameMesh;
			FrameMesh.Create("frame1");

			//get mesh cur center pos

			auto& aabb = RBoundingVolumesManagerS::GetBox(FrameMesh.Mesh.GlobalId);
			float3 curPos = aabb.GetPosition();

			FrameMesh.TranslateCPU(float3(-curPos.x, -curPos.y, -curPos.z));
			FrameMesh.RescaleCPU(0.1f);

			FrameMeshVertexBufferCPU = FrameMesh.GeometryBufferCPU.VertexBufferPositionsCPU;

			FrameMesh.RescaleCPU(float3{ MeshScale, 1.f });

			FrameMesh.RotateCPU(float3{ 0, 0, MeshOrientationOffset });

			//Translate to simulation center
			float3 center = GetPaintSimulationCenterPos();

			float frameZOffset = -0.25f;
			//frameZOffset *= 4.f;

			FrameMesh.TranslateCPU(float3(center.x, center.y, center.z + frameZOffset));

			MeshVertexBuffer = FrameMesh.GpuData.VertexBufferPositions;
			MeshIndexBuffer = FrameMesh.GpuData.IndexBuffer;
			MeshVertexBufferRest = FrameMesh.GpuData.VertexBufferRest;

			MeshDesc = FrameMesh.Mesh.Desc;
		}
		
		if (bReUploadScene)
		{
			auto tempVB = FrameMeshVertexBufferCPU;

			RescaleVertexBuffer(float3{ MeshScale, 1.f }, tempVB);

			RotateVertexBuffer(float3{ 0, 0, MeshOrientationOffset }, tempVB);

			//Translate to simulation center
			float3 center = GetPaintSimulationCenterPos();

			float frameZOffset = -0.25f;
			//frameZOffset *= 4.f;

			TranslateVertexBuffer(float3(center.x, center.y, center.z + frameZOffset), tempVB);

			CmdList.UploadToBuffer(MeshVertexBuffer, 0, tempVB.Data.data(), tempVB.GetNumVertices()* tempVB.GetVertexSize());
		}


	}

	void AllocateGBufferRenderingResources()
	{
		MeshTextureAlbedo.CreateFromBitmap(RBitmap::FromFile("assets//paint/paint1/Albedo.png"), DXGI_FORMAT_R8G8B8A8_UNORM);
		MeshTextureNormal.CreateFromBitmap(RBitmap::FromFile("assets/paint/paint1/Normal.png"), DXGI_FORMAT_R8G8B8A8_UNORM);
		MeshTextureRoughness.CreateFromBitmap(RBitmap::FromFile("assets/paint/paint1/RMA.png"), DXGI_FORMAT_R8G8B8A8_UNORM);

		WallTextureAlbedo.CreateFromBitmap(RBitmap::FromFile("assets/paint/wall/brown_wallpaper_21_73_diffuseSmall.jpg"), DXGI_FORMAT_R8G8B8A8_UNORM);
		WallTextureNormal.CreateFromBitmap(RBitmap::FromFile("assets/paint/wall/brown_wallpaper_21_73_normalSmall.jpg"), DXGI_FORMAT_R8G8B8A8_UNORM);
		WallTextureRoughness.CreateFromBitmap(RBitmap::FromFile("assets/paint/wall/brown_wallpaper_21_73_glossinessSmall.jpg"), DXGI_FORMAT_R8G8B8A8_UNORM);

		AdditionalTextureAlbedo.CreateFromBitmap(RBitmap::FromFile("assets//paint/additional/Gold_DIF.png"), DXGI_FORMAT_R8G8B8A8_UNORM);
		AdditionalTextureNormal.CreateFromBitmap(RBitmap::FromFile("assets/paint/additional/Gold_NRM.png"), DXGI_FORMAT_R8G8B8A8_UNORM);
		AdditionalTextureRoughness.CreateFromBitmap(RBitmap::FromFile("assets/paint/additional/Gold_RGH.png"), DXGI_FORMAT_R8G8B8A8_UNORM);
		AdditionalTextureMetalness.CreateFromBitmap(RBitmap::FromFile("assets/paint/additional/Gold_MTL.png"), DXGI_FORMAT_R8G8B8A8_UNORM);

		GBufferColor.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R11G11B10_FLOAT);
		GBufferNormal.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R32G32B32A32_FLOAT);
		GBufferRoughness.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R32_FLOAT);
		GBufferMetal.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R32_FLOAT);
		GBufferWorldSpaceZ.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R32_FLOAT);
		DepthBuffer.Create("PaintingDepthBuffer", PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_D32_FLOAT);

		GBufferStencil.CreateRenderTarget(PaintingBackgroundSizePixels.x, PaintingBackgroundSizePixels.y, DXGI_FORMAT_R8_UINT);


		PaintingSource.CreateFromBitmap(RBitmap::FromFile(RTexturePool::filepath +  std::string{ RTexturePool::Paintings[0] }), DXGI_FORMAT_R8G8B8A8_UNORM);
		//OriginalPaintSourceXYRatio = float(PaintingSource.Size.x) / float(PaintingSource.Size.y);

	}

	void DrawPaintingBackground(RCommandListGraphics& CmdList)
	{
		GPU_PROFILE_SCOPE(BackgroundRendering, CmdList);


		//CGPaintRenderConstantsCB.GlobalLightDirection = MapSphericalToCartesian(GlobalLightDirectionSphericalCoords.x, GlobalLightDirectionSphericalCoords.y, 1);

		auto& RootSignature = PaintRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSignature);

		//define Viewport
		static CD3DX12_VIEWPORT Viewport;
		static CD3DX12_RECT ScissorRect;

		const uint2 viewportSize = GBufferColor.Size;

		Viewport = CD3DX12_VIEWPORT{ 0.0f, 0.0f, (FLOAT)viewportSize.x, (FLOAT)viewportSize.y };
		ScissorRect = CD3DX12_RECT{ 0, 0, (LONG)viewportSize.x, (LONG)viewportSize.y };

		//Set global states 
		/*CmdList.SetViewport(RSceneRenderer::Get().GetViewport());
		CmdList.SetScissor(RSceneRenderer::Get().GetScissorRect());*/
		CmdList.SetViewport(Viewport);
		CmdList.SetScissor(ScissorRect);

		/*auto& RT = RCommonResources::Get().SceneColor;
		auto& DSV = RCommonResources::Get().SceneDepth;

		CmdList.SetRenderTarget(RT.GetRTV(), DSV.GetDSV());*/

		CmdList.TransitionResource(GBufferColor, D3D12_RESOURCE_STATE_RENDER_TARGET, true);
		CmdList.TransitionResource(DepthBuffer, D3D12_RESOURCE_STATE_DEPTH_WRITE, true);

		CmdList.ClearRenderTarget(GBufferColor);
		CmdList.ClearRenderTarget(GBufferMetal);
		CmdList.ClearRenderTarget(GBufferNormal);
		CmdList.ClearRenderTarget(GBufferRoughness);
		CmdList.ClearRenderTarget(GBufferWorldSpaceZ);
		CmdList.ClearRenderTarget(GBufferStencil);
		CmdList.ClearDepth(DepthBuffer);

		static constexpr uint NumRenderTargets = 6;

		D3D12_CPU_DESCRIPTOR_HANDLE rtvHandles[NumRenderTargets];
		rtvHandles[0] = GBufferColor.GetRenderTargetView();
		rtvHandles[1] = GBufferNormal.GetRenderTargetView();
		rtvHandles[2] = GBufferRoughness.GetRenderTargetView();
		rtvHandles[3] = GBufferWorldSpaceZ.GetRenderTargetView();
		rtvHandles[4] = GBufferStencil.GetRenderTargetView();
		rtvHandles[5] = GBufferMetal.GetRenderTargetView();

		CmdList.SetRenderTargets(NumRenderTargets, rtvHandles, DepthBuffer.GetDSV());

		//Set Constants

		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("SceneUniformBufferCBuffer"), RScene::Get().UniformBufferParameters.GetGlobalUniformBuffer()->ConstantBuffer.GetResource()->GetGPUVirtualAddress());


		//Update values

		CGPaintRenderConstantsCBufferGPU.Update(CGPaintRenderConstantsCB);
		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("CGPaintRenderConstantsCB"), CGPaintRenderConstantsCBufferGPU.GetGpuVirtualAddress());

		//SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("LightsBuffer"), LightsBufferGPU, true);

		//PSO 
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;
		//psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateDisabled;
		psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;

		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

		DXGI_FORMAT rtFormats[NumRenderTargets];
		rtFormats[0] = GBufferColor.Format;
		rtFormats[1] = GBufferNormal.Format;
		rtFormats[2] = GBufferRoughness.Format;
		rtFormats[3] = GBufferWorldSpaceZ.Format;
		rtFormats[4] = GBufferStencil.Format;
		rtFormats[5] = GBufferMetal.Format;

		psoStates.NumRenderTargets = NumRenderTargets;
		psoStates.pRenderTargetFormat = rtFormats;

		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);


		//Render Frame Mesh
		{
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("IndexBuffer"), MeshIndexBuffer, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), MeshVertexBuffer, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferRest"), MeshVertexBufferRest, true);

			//SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), MeshTextureAlbedo, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), MeshTextureNormal, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("RoughnessTexture"), MeshTextureRoughness, true);

			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), AdditionalTextureAlbedo, true);
			//SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), AdditionalTextureNormal, true);
			//SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("RoughnessTexture"), AdditionalTextureRoughness, true);

			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("MetalnessTexture"), AdditionalTextureMetalness, true);

			static BackgroundMeshRenderVertexShader VS;
			static BackgroundFrameRenderPixelShader PS;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

			CmdList.SetPipelineState(PSO);

			CmdList.Draw(MeshDesc.NumofIndices, 0);
		}
		 
		//Render Background wall
		{
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), PaintRenderVertexBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("TexCoordBuffer"), PaintTexCoordsBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), WallTextureAlbedo, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), WallTextureNormal, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("RoughnessTexture"), WallTextureRoughness, true);

			static BackgroundWallRenderVertexShader VS;
			static BackgroundWallRenderPixelShader PS;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

			CmdList.SetPipelineState(PSO);

			CmdList.Draw(6, 0);
		}

		//Render Painting
		{
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), PaintSourceVertexBufferGPU, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("TexCoordBuffer"), PaintTexCoordsBufferGPU, true);

			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("ColorTexture"), PaintingSource, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("NormalTexture"), WallTextureNormal, true);
			SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("RoughnessTexture"), WallTextureRoughness, true);

			static BackgroundWallRenderVertexShader VS;
			static PaintSourceRenderPixelShader PS;

			auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

			CmdList.SetPipelineState(PSO);

			CmdList.Draw(6, 0);
		}

	}


	

	/*
		ShadowMap Render
	*/

	//RDepthBuffer GlobalLightDepth{ 1.f };
	//RTexture ShadowMap{ TEXT("ShadowMap") };
	//RDescriptorHandle ShadowMapDescHandle;

	//uint2 ShadowMapResolution{ 1024,1024 };
	uint2 ShadowMapResolution{ 2048,2048 };
	float2 ShadowMapScale{ 1,1 };

	float3 GlobalLightPosition{ 0, 5, 0 };
	bool bAssignGlobalLightPosFromCamera = true;


	RDynamicVector<RDepthBuffer> DepthBufferPerLightGPUArr;
	RDynamicVector<RTexture> ShadowMapPerLightGPUArr;
	RDynamicVector<RTexture> ShadowVisibilityPerLightTextureGPUArr;
	RDynamicVector<RDescriptorHandle> DepthBufferDescHandleGPUArr;
	RDynamicVector<RDescriptorHandle> ShadowVisibilityDescHandleGPUArr;

	TStructuredBuffer<LightTransformStruct> LightSpaceTransformMatricesPerLightBufferGPU;
	RDynamicVector<LightTransformStruct> LightSpaceTransformMatricesPerLightBufferCPU;


	void UpdateShadowMapResources(RCommandListGraphics& CmdList)
	{
		//generate transform matrices for each light
		const auto numLights = LightsBufferCPU.size();
		LightSpaceTransformMatricesPerLightBufferCPU.resize(numLights);
		for (int i = 0; i < numLights; i++)
		{
			auto& curLight = LightsBufferCPU[i];

			LightTransformStruct lightMats;

			lightMats.View = DirectX::XMMatrixLookToLH(float3SSE(curLight.Position), float3SSE{ curLight.Direction }, float3SSE{ 0,1,0 });

			float nearPlane = 0.1f;
			float farPlane = 50.f;

#if 1
			matrix proj = DirectX::XMMatrixOrthographicLH(PaintSceneCameraScale.x * ShadowMapScale.x, PaintSceneCameraScale.y * ShadowMapScale.x, nearPlane, farPlane);

#else
			float FieldOfView{ 1.57079632679f };// Pi/2 angle doesn't scale XY 
			matrix proj = DirectX::XMMatrixPerspectiveFovLH(FieldOfView, 1, nearPlane, farPlane);
#endif
			lightMats.ViewProjection = DirectX::XMMatrixMultiply(lightMats.View, proj);
			LightSpaceTransformMatricesPerLightBufferCPU[i] = lightMats;
		}

		if (LightSpaceTransformMatricesPerLightBufferGPU.GetElementCount() < numLights)
		{
			LightSpaceTransformMatricesPerLightBufferGPU.Create("LightSpaceTransformMatricesPerLightBuffer", numLights, LightSpaceTransformMatricesPerLightBufferCPU.data());
		}

		CmdList.UploadToBuffer(LightSpaceTransformMatricesPerLightBufferGPU, 0, LightSpaceTransformMatricesPerLightBufferCPU.data(), LightSpaceTransformMatricesPerLightBufferCPU.size() * sizeof(LightTransformStruct));

		//shadow maps


		auto numSHadowMaps = DepthBufferPerLightGPUArr.size();

		if (numSHadowMaps < numLights)
		{
			DepthBufferPerLightGPUArr.resize(numLights);
			ShadowMapPerLightGPUArr.resize(numLights);
			DepthBufferDescHandleGPUArr.resize(numLights);
			ShadowVisibilityDescHandleGPUArr.resize(numLights);
			ShadowVisibilityPerLightTextureGPUArr.resize(numLights);
			for (int i = numSHadowMaps; i < numLights; i++)
			{
				auto& newDepthBuffer = DepthBufferPerLightGPUArr[i];
				newDepthBuffer.DepthClear = 1.f;
				newDepthBuffer.Create(TEXT("LightDepthBuffer"), ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_D32_FLOAT);

				auto& newShadowMap = ShadowMapPerLightGPUArr[i];
				newShadowMap.CreateRenderTarget(ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_R32_FLOAT);

				//create GPU visible handle for UI
				auto& newHandle = DepthBufferDescHandleGPUArr[i];
				newHandle = SLLRenderer::Get().GetDescHeapCBVUAVSRVShaderVisible().GetNextFreeHandle();
				SLLRenderer::GetDevice()->CreateShaderResourceView(newShadowMap.GetResource(), nullptr, newHandle.CPU);


				ShadowVisibilityPerLightTextureGPUArr[i].Create(ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_R32_FLOAT, true);
				//create GPU visible handle for UI
				auto& newHandle2 = ShadowVisibilityDescHandleGPUArr[i];
				newHandle2 = SLLRenderer::Get().GetDescHeapCBVUAVSRVShaderVisible().GetNextFreeHandle();
				SLLRenderer::GetDevice()->CreateShaderResourceView(ShadowVisibilityPerLightTextureGPUArr[i].GetResource(), nullptr, newHandle2.CPU);
			}
		}

	}


	void RenderShadowMap(RCommandListGraphics& CmdList)
	{

		GPU_PROFILE_SCOPE(RenderShadowMap, CmdList);

		//Generate Light Matrix
		if (GbFirstBoot)
		{

			//GlobalLightDepth.Create("PaintingDepthBuffer", ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_D32_FLOAT);

			//ShadowMap.CreateRenderTarget(ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_R32_FLOAT);


			////Allocate SRV Gpu-Visible
			//ShadowMapDescHandle = SLLRenderer::Get().GetDescHeapCBVUAVSRVShaderVisible().GetNextFreeHandle();
			//SLLRenderer::GetDevice()->CreateShaderResourceView(ShadowMap.GetResource(), nullptr, ShadowMapDescHandle.CPU);


		}

		UpdateLights(CmdList);
		UpdateShadowMapResources(CmdList);

		if (bAssignGlobalLightPosFromCamera)
		{
			CGPaintRenderConstantsCB.GlobalLightViewMatrix = RScene::Get().GetCamera().Matrices.View;
		}
		else
		{
			//CGPaintRenderConstantsCB.GlobalLightDirection = MapSphericalToCartesian(GlobalLightDirectionSphericalCoords.x, GlobalLightDirectionSphericalCoords.y, 1);
			//CGPaintRenderConstantsCB.GlobalLightViewMatrix = DirectX::XMMatrixLookToLH(float3SSE(GlobalLightPosition), float3SSE{ CGPaintRenderConstantsCB.GlobalLightDirection }, float3SSE{ 0,1,0 });
		}

		auto proj = DirectX::XMMatrixOrthographicLH(PaintSceneCameraScale.x * ShadowMapScale.x, PaintSceneCameraScale.y * ShadowMapScale.x, 0.01, 20);
		CGPaintRenderConstantsCB.GlobalLightViewProjectionMatrix = DirectX::XMMatrixMultiply(CGPaintRenderConstantsCB.GlobalLightViewMatrix, proj);

		auto& RootSignature = PaintRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSignature);

		//define Viewport
		static CD3DX12_VIEWPORT Viewport;
		static CD3DX12_RECT ScissorRect;

		const uint2 viewportSize = ShadowMapResolution;

		Viewport = CD3DX12_VIEWPORT{ 0.0f, 0.0f, (FLOAT)viewportSize.x, (FLOAT)viewportSize.y };
		ScissorRect = CD3DX12_RECT{ 0, 0, (LONG)viewportSize.x, (LONG)viewportSize.y };

		//Set global states 
		CmdList.SetViewport(Viewport);
		CmdList.SetScissor(ScissorRect);


		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("SceneUniformBufferCBuffer"), RScene::Get().UniformBufferParameters.GetGlobalUniformBuffer()->ConstantBuffer.GetResource()->GetGPUVirtualAddress());

		CGPaintRenderConstantsCBufferGPU.Update(CGPaintRenderConstantsCB);
		CmdList.SetConstantBuffer(RootSignature.GetRootParamIndex("CGPaintRenderConstantsCB"), CGPaintRenderConstantsCBufferGPU.GetGpuVirtualAddress());

		//PSO 
		RPSOFactory::CommonPipelineStates psoStates;

		psoStates.RasterDesc = RCommonResources::Get().RasterizerTwoSided;
		psoStates.DepthStencilDesc = RCommonResources::Get().DepthStateReadWrite;

		psoStates.BlendDesc = RCommonResources::Get().BlendDisable;
		psoStates.TopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;

		psoStates.pRenderTargetFormat = &ShadowMapPerLightGPUArr[0].Format;

		CmdList.SetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

		SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("LightSpaceTransformMatricesPerLightBufferSRV"), LightSpaceTransformMatricesPerLightBufferGPU, true);

		//Render the scene for each light

		const auto numLights = LightsBufferCPU.size();

		for (int i = 0; i < numLights; i++)
		{
			GPU_PROFILE_SCOPE_TEXT(CmdList, L"(Light %i)", i);
			auto& curDepthBuffer = DepthBufferPerLightGPUArr[i];
			auto& curShadowMapRT = ShadowMapPerLightGPUArr[i];

			CmdList.TransitionResource(curDepthBuffer, D3D12_RESOURCE_STATE_DEPTH_WRITE, true);

			CmdList.ClearRenderTarget(curShadowMapRT);
			CmdList.ClearDepth(curDepthBuffer);

			CmdList.SetRenderTarget(curShadowMapRT.GetRenderTargetView(), curDepthBuffer.GetDSV());

			CmdList.SetConstantArray(RootSignature.GetRootParamIndex("UtilityConstants"), 1, &i);

			//Render Frame Mesh
			{
				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("IndexBuffer"), MeshIndexBuffer, true);
				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), MeshVertexBuffer, true);
				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferRest"), MeshVertexBufferRest, true);

				static BackgroundDepthIndexedRenderVertexShader VS;
				static BackgroundDepthRenderPixelShader PS;

				auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

				CmdList.SetPipelineState(PSO);

				CmdList.Draw(MeshDesc.NumofIndices, 0);
			}

			//Render Background wall
			{
				SetDescriptorSRV(CmdList, RootSignature.GetDescriptorIndex("VertexBufferPosition"), PaintRenderVertexBufferGPU, true);

				static BackgroundDepthRenderVertexShader VS;
				static BackgroundDepthRenderPixelShader PS;

				auto PSO = RPSOFactory::GetGraphicsPSO(RootSignature, VS, &PS, psoStates);

				CmdList.SetPipelineState(PSO);


				CmdList.Draw(6, 0);
			}
		}
	}




	class GenerateShadowVisibilityShader : public RShader
	{
	public:
		GenerateShadowVisibilityShader()
		{
			Create("Physics/ParticlePhysics.hlsl", "GenerateShadowVisibility", L"cs_6_0");
		}
	};

	bool bBlurShadowVisibility = false;

	void GenerateShadowVisibilityTextures(RCommandListCompute& CmdList)
	{
		GPU_PROFILE_SCOPE(GenerateShadowVisibilityTextures, CmdList);

		auto& RootSig = PaintRenderingRootSig::Get();

		CmdList.SetRootSignature(RootSig);

		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGConstantsCB"), CGConstantsCBufferGPU.GetGpuVirtualAddress());
		CmdList.SetConstantBuffer(RootSig.GetRootParamIndex("CGPaintRenderConstantsCB"), CGPaintRenderConstantsCBufferGPU.GetGpuVirtualAddress());

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("LightSpaceTransformMatricesPerLightBufferSRV"), LightSpaceTransformMatricesPerLightBufferGPU, true);

		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferWorldSpaceZSRV"), GBufferWorldSpaceZ, true);
		SetDescriptorSRV(CmdList, RootSig.GetDescriptorIndex("GBufferNormalTexture"), GBufferNormal, true);


		const auto numLights = LightsBufferCPU.size();
		for (int i = 0; i < numLights; i++)
		{
			GPU_PROFILE_SCOPE_TEXT(CmdList, L"(Light %i)", i);

			CmdList.SetConstantArray(RootSig.GetRootParamIndex("UtilityConstants"), 1, &i);

			CmdList.TransitionResource(DepthBufferPerLightGPUArr[i], D3D12_RESOURCE_STATE_GENERIC_READ);
			CmdList.SetDynamicDescriptor(RRootSignatureExt::GetSRVRootIndex(), RootSig.GetDescriptorIndex("ShadowMapTextureSRV"), DepthBufferPerLightGPUArr[i].GetDepthSRV());

			SetDescriptorUAV(CmdList, RootSig.GetDescriptorIndex("ShadowVisibilityTextureUAV"), ShadowVisibilityPerLightTextureGPUArr[i], true);

			static GenerateShadowVisibilityShader ComputeShader;
			CmdList.SetPipelineState(RPSOFactory::GetComputePSO(RootSig, ComputeShader));

			constexpr uint ThreadGroupSize = 8;
			CmdList.Dispatch(RMath::DivideAndRoundUp(ShadowMapResolution.x, ThreadGroupSize), RMath::DivideAndRoundUp(ShadowMapResolution.y, ThreadGroupSize), 1);


		}

		if (bBlurShadowVisibility)
		{
			GPU_PROFILE_SCOPE(Blur, CmdList);

			for (int i = 0; i < numLights; i++)
			{
				static RTexture intermediateBlurTexture;
				if (intermediateBlurTexture.Size.x == 0)
				{
					intermediateBlurTexture.Create(ShadowMapResolution.x, ShadowMapResolution.y, DXGI_FORMAT_R32_FLOAT, true);
				}
				RTexture* pSource = &ShadowVisibilityPerLightTextureGPUArr[i];
				RTexture* pDest = &intermediateBlurTexture;
				for (int i = 0; i < PaintBlurPassNumIterations; i++)
				{
					BlurPassGPU(CmdList, *pSource, *pDest, BlurPassType::TYPE_FLOAT, true);
					std::swap(pSource, pDest);
				}
			}
		}


		
	}
























#if 0
	void AssignAllDescriptors()
	{
		auto& RootSig = ParticlePhysicsRootSig::Get();

		SetDescriptorSRV(RootSig.GetDescriptorIndex("PlanesBufferSRV"), PlanesBufferGPU);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ParticlesBufferPositionCurSRV"), ParticlesBufferPositionCurGPU);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("CellIdOffsetAndCountBufferSRV"), CellIdOffsetAndCountBufferGPUArr);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ObjectListIndexBufferSRV"), 3);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("CellListIndexBufferSRV", 4);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("PerObjectControlBitsBufferSRV", 5);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsBufferSRV", 6);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ParticlesStateBufferSRV", ParticlesStateBufferGPU);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsStatesBufferSRV", 8);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("ParticlesToRenderCoordListSRV", 9);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsStatesBufferBLUESRV", 10);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsStatesBufferGREENSRV", 11);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsStatesBufferYELLOWSRV", 12);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("ConstraintsStatesBufferREDSRV", 13);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("SphereCapsulesBufferSRV", 14);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("ParticlesBufferPositionAfterCreationSRV", 15);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureSRV", 16);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureSRV", 17);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("CurParticlesVelocityBufferSRV", 18);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("VelocityFieldStateTextureSRV", 19);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("DensityFieldTextureSRV", 20);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("PressureFieldTextureSRV", 21);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("VectorFieldDivergenceTextureSRV", 22);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("VelocityField_U_ComponentTextureNewSRV", 23);
		SetDescriptorSRV(RootSig.GetDescriptorIndex("VelocityField_V_ComponentTextureNewSRV", 24);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("VectorFieldCurlTextureSRV", 25);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("DensityFieldTextureNewSRV", 26);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("DensitySourceColorTextureSRV", 27);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("InputTextureWithNegativeValuesSRV", 28);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("CurlGradientTextureSRV", 29);

		SetDescriptorSRV(RootSig.GetDescriptorIndex("CurlGradientMinMaxTextureSRV", 30);
	}
#endif
};




void GpuParticlesPhysicsPass::Execute(GpuParticlesPhysicsPass::OutputStruct& outResult)
{
	PROFILE_FUNC();

	static GPUParticlePhysics GpuPhysics;

	GpuPhysics.RenderUI();
	GpuPhysics.Simulate();

	if (GpuPhysics.bSubdivisionMesh)
	{
		outResult.pPositionsBuffer = &GpuPhysics.InterpolatedParticlePositionsBufferGPU;
		outResult.NumInstances = GpuPhysics.OutputMeshSize.x * GpuPhysics.OutputMeshSize.y;
		outResult.InstanceRadius = GpuPhysics.CGConstantsCB.PointRadius;
	}
	else
	{
		outResult.pPositionsBuffer = &GpuPhysics.ParticlesBufferPositionCurGPU;
		outResult.NumInstances = GpuPhysics.CGConstantsCB.NumParticles;
		//static constexpr float krenderedParticleScale = 1.5f;
		static constexpr float krenderedParticleScale = 1.f;
		outResult.InstanceRadius = GpuPhysics.CGConstantsCB.PointRadius * krenderedParticleScale;

		outResult.pColorsBuffer = &GpuPhysics.ParticlesColorsBufferGPU;
	}

	

	outResult.MeshBufferVerticesGPU = &GpuPhysics.MeshBufferVerticesGPU;
	outResult.NumTriangles = GpuPhysics.MeshGeneratorCurNumVertices;

	if (GpuPhysics.bGenerateTriangleMeshNormals)
	{
		outResult.MeshBufferNormalsGPU = &GpuPhysics.MeshBufferNormalsGPU;
	}

	if (GpuPhysics.pGeneratedMeshColorTexture)
	{
		outResult.MeshBufferTexCordsGPU = &GpuPhysics.MeshBufferTexCoordsGPU;
		outResult.pColorTexture = GpuPhysics.pGeneratedMeshColorTexture;
	}

	if (GpuPhysics.pGeneratedMeshNormalTexture)
	{
		outResult.pNormalTexture = GpuPhysics.pGeneratedMeshNormalTexture;
	}

	GpuPhysics.MeshGeneratorCurNumVertices = 0;
}




#ifndef _GPU_PARTICLES_COMMON_STRUCTS_H
#define _GPU_PARTICLES_COMMON_STRUCTS_H

#define USE_3COLOR_DENSITY 1
#if USE_3COLOR_DENSITY == 1
#define DENSITY_TYPE float3
#else
#define DENSITY_TYPE float
#endif



#define RENDER_PARTICLES_SDF 0
#define PAINT_RENDER_ADVECT_SOURCE 0
#define PAINT_RENDER_PIXEL_ADVECT 0

struct
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
CGConstants
{
	float3 ExternalForce;
	uint NumParticles;

	float DeltaTime;
	float Damping;
	float PointRadius;
	float CollisionResponseCoef;

	float3 RandomImpulse;
	uint bSpawnObjectsMode;

	float3 DefaultParticlePos;
	float SpringRestLength;

	float ConnectionConstraintStiffness;
	uint2 ClothSize; //num Knots, not Springs
	float ConnectionConstraintToDefaultPosStiffness;

	float SimulationBoundsLength;
	float SimulationBoundsStartOffset;
	bool bSolversUseOverrelaxation;
	float OverrelaxationConstant;

	float VelocityFieldBoundsLength;
	float VelocityFieldBoundsStartOffset;
	uint VelocityFieldGridSize;
	float VelocityFieldGridCellLength;

	float DiffusionFactor;
	float VorticityConfinementFactor;
	float PressureControlledOffsetScale;
	float WaveEquationWaveSpeed;

	uint2 CurVelocityFieldCellIdMousePointsTo;
	float VelocityGeneratorRadius;
	float DensityGeneratorRadius;

	float2 CurPointerVelocity;
	float MaxParticleVelocityMagnitudeThreshold;
	bool bSimulateParticlesIn2D;

#if USE_3COLOR_DENSITY == 1
	float3 InitialDensityColor;
	float pad;
#else
	float InitialDensityColor;
	float3 pad;
#endif

	float VelocityAdvectionScaleFactor;
	float GeneratedVelocityMagnitudeScale;
	bool bUseParticleBasedAdvection;
	float AdvectedParticleDriftCompensationScale;

	float ParticleAdvectionPICFLIPRatio;
	bool bSimulateParticlesIn2DAllowNegativeZOffset;
	float RasterizedParticleRadiusScale;

#ifdef CPP_SHADER
	CGConstants()
	{
		//ExternalForce = float3{ 0,-100,0 };
		//ExternalForce = float3{ 0,-250,0 };
		ExternalForce = float3{ 0,0,0 };

		//NumParticles = 1;
		NumParticles = 0;
		DeltaTime = 1.f / 60.f;
		Damping = 0.99f;

		//PointRadius = 0.1f;
		PointRadius = 0.05f;
		//PointRadius = 1.0f;

		RasterizedParticleRadiusScale = 1.1f;

		RandomImpulse = float3{0,0,0};
		bSpawnObjectsMode = 0;
		CollisionResponseCoef = 0.75f;
		//CollisionResponseCoef = 1.f;

		//DefaultParticlePos = float3{ 1,1,1 };
		DefaultParticlePos = float3{ 3.5,5.5,4 };
		DefaultParticlePos = float3{ 3.5,5.5,8 };
		//DefaultParticlePos = float3{ 3.5,5.5,16 };

		#if RENDER_PARTICLES_SDF
		DefaultParticlePos = float3{ 3.5,5.5,4 };
		#endif

		SpringRestLength = PointRadius * 2.f;
		//SpringRestLength = 5.f;

		ConnectionConstraintStiffness = 1.f;
		ConnectionConstraintToDefaultPosStiffness = 0.5f;

		ClothSize = uint2(33,33);

		SimulationBoundsLength = 100.f;
		SimulationBoundsStartOffset = 1.f;

		VelocityFieldBoundsLength = 20.f;
		VelocityFieldBoundsStartOffset = 1.f;

		//VelocityFieldGridSize = 64;
		//VelocityFieldGridSize = 128;
		//VelocityFieldGridSize = 256;
		VelocityFieldGridSize = 512;

		#if PAINT_RENDER_PIXEL_ADVECT
		VelocityFieldGridSize = 512;
		#endif

		DiffusionFactor = 1.f;
		VorticityConfinementFactor = 1.f;

		PressureControlledOffsetScale = 1.f;

		WaveEquationWaveSpeed = 1.f;

		VelocityGeneratorRadius = 1.5f;
		DensityGeneratorRadius = 0.5f;

		#if USE_3COLOR_DENSITY == 1
		InitialDensityColor = float3(1,0.5,0.5);
		#else
		InitialDensityColor = 1;
		#endif

		VelocityAdvectionScaleFactor = 1.f;

		bSolversUseOverrelaxation = true;
		OverrelaxationConstant = 1.9f;

		bSimulateParticlesIn2D = true;
		bSimulateParticlesIn2DAllowNegativeZOffset = true;

		MaxParticleVelocityMagnitudeThreshold = 25.f;

		//GeneratedVelocityMagnitudeScale = 50.f;
		GeneratedVelocityMagnitudeScale = 25.f;

		bUseParticleBasedAdvection = false;

		//AdvectedParticleDriftCompensationScale = 0.25f;
		AdvectedParticleDriftCompensationScale = 0.f;

		//ParticleAdvectionPICFLIPRatio = 0.5f;
		ParticleAdvectionPICFLIPRatio = 0.25f;
		//ParticleAdvectionPICFLIPRatio = 0.0f;

	}
#endif//CPP_SHADER
};

struct GCollisionConstants
{
	uint CurCellLocalId;
	uint NumCells;
	float PerPassSpringStiffnessScale;
};

static const float CellSizeScale = 2.15f;
//static const float CellSizeScale = 2.f * 2.15f;

//From 0 to 7
uint GetLocalCellIndex(uint3 cellCoord)
{
    uint index = 0;
    index |= ((cellCoord.x & 1) << 0);
    index |= ((cellCoord.y & 1) << 1);
    index |= ((cellCoord.z & 1) << 2);
    return index;
}
uint GetLocalCellIndex(int3 cellCoord)
{
    uint index = 0;
    index |= ((cellCoord.x & 1) << 0);
    index |= ((cellCoord.y & 1) << 1);
    index |= ((cellCoord.z & 1) << 2);
    return index;
}

// Get the hash id for a cell given its (x, y, z) coordinates
#define XSHIFT 0
#define YSHIFT 10
#define ZSHIFT 20
#define CELL_MAX_INDEX 0x3FF
uint GetCellIndex(int3 cellId)
{
	return 0x3fffffffu & ((cellId.x << XSHIFT) | (cellId.y << YSHIFT) | (cellId.z << ZSHIFT));
}
#define BITMASK10 0x3FF
int3 GetCellCoord(uint cellIndex)
{
	int3 cellId;
	cellId.x = int((cellIndex >> XSHIFT) & BITMASK10);
	cellId.y = int((cellIndex >> YSHIFT) & BITMASK10);
	cellId.z = int((cellIndex >> ZSHIFT) & BITMASK10);
	return cellId;
}

//PARTICLE STATES
#define BITMASK_PARTICLE_STATE_PINNED (1u << 0)
#define BITMASK_PARTICLE_STATE_SELECTED (1u << 1)
#define BITMASK_PARTICLE_STATE_DISABLED (1u << 2)
#define BITMASK_PARTICLE_STATE_PINNED_SOFT (1u << 3)
#define BITMASK_PARTICLE_STATE_FALLING (1u << 4)
#define BITMASK_PARTICLE_STATE_NEW (1u << 5)
#define BITMASK_PARTICLE_STATE_STICKY (1u << 6)

///CONSTRAINT STATES
#define BITMASK_CONSTRAINT_STATE_BROKEN (1 << 0)

struct CConstraint
{
	uint IndexToKnot1;
	uint IndexToKnot2;
};

///Constraints divided into groups of 4
#define CONSTRAINT_COLOR_BLUE 0 //x0
#define CONSTRAINT_COLOR_GREEN 1 //x1
#define CONSTRAINT_COLOR_YELLOW 2 //y0
#define CONSTRAINT_COLOR_RED 3 //y1
#define CONSTRAINT_NUM_COLOR_TYPES 4

struct 
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
SpatialControlConstants 
{
	float3 RayOrigin;
	uint PrevFrameIntersectedPointIndex;

	float3 RayEnd;
	uint bInDragState;

	float3 MouseOffsetWorld;
	float IntersectionSphereScale;

	float3 ConstraintScissorPos;
	float ScissorAffectRadius;

	float2 ParticleSpawnerPosition;
	uint2 ParticleSpawnAreaSize;

#ifdef CPP_SHADER
	SpatialControlConstants()
	{
		IntersectionSphereScale = 1.5f;

		ScissorAffectRadius = 0.1f;

		ParticleSpawnAreaSize = uint2(1,1);
	}
#endif
};


//Mesh generator
struct GMeshGenerationConstants
{
	uint2 MeshSize;
};



 /* 
 	Velocity Field
 */

struct VectorProjectionPerPassConstants //TODO: Make those constants uniform for all passes that require sub-passes
{
	uint Offset;
};

#define WATER_SURFACE_SIMULATION_HEIGHT 5.f

#define VELOCITY_FIELD_STATE_TYPE uint32

#define VELOCITY_FIELD_CELL_STATE_SOLID 0
#define VELOCITY_FIELD_CELL_STATE_DENSITY 1
#define VELOCITY_FIELD_CELL_STATE_AIR 2

#define PARTICLE_ADVECTION_DEBUG 0
#define PARTICLE_ADVECTION_DEBUG_VELOCITY float3(1,1,0)

#define PARTICLE_ADVECTION_FFP_FOR_DATA_TRANSFER 1 //use VS PS to transfer particle data into velocity field


//===============================================
//
//				SDF Render
//
//===============================================

#define SDF_BLEND_EXP 0
#define SDF_BLEND_ROOT 1
#define SDF_BLEND_POWER 2
#define SDF_BLEND_POLY 3
#define SDF_BLEND_CUBIC 4




struct 
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
SDFRenderConstants 
{
	uint2 SDFOutputSize;
	float Threshold;
	float BlendRadius;

	float3 OutputColor;
	float NormalsSmoothness;

	float ColorBlendRadius;
	float ColorFadeThreshold;
	float ColorFadeScaleInv;
	float ColorDiffusionScale;

	float PaintDryAreasVelocityDampForceScale;
	float DrySurfaceDampThreshold;
	float SpecularAmount;
	float NormalsBlendRadius;

	float RayOriginDepthOffset;
	float SurfaceDepthOffset;
	uint SDFShapeBlendMode;
	uint SDFNormalsBlendMode;

	float WetnessAccumulationAmount;
	float WetnessDryingOutAmount;
	float WetnessLimit;
	float PaintDiffusionFactor;

	float WetnessContributionToHeight;
	float FallingParticleBlendRadius;
	uint2 GBufferSize;

	float FalingParticleStartPosZ;
	float FallingParticleSpeedTowardsSurface;
	float2 CameraPositionWorld;

	uint ParticleSourceTextureSize;
	float SurfaceHitVelocityScale;
	float SurfaceHitVelocityMax;
	float ParticleBumpinnesPow;

#ifdef CPP_SHADER
	SDFRenderConstants()
	{
		//SDFOutputSize = uint2(512,512);
		SDFOutputSize = uint2(1024,1024);

		Threshold = 0.001f;
		BlendRadius = 16.f;

		FallingParticleBlendRadius = 12.f;

		OutputColor = float3(0.2, 1., 0.2);

		NormalsSmoothness = 0.2f;

		ColorBlendRadius = 0.1f;

		ColorFadeThreshold = 0.3f;
		ColorFadeScaleInv = 0.995f;

		ColorDiffusionScale = 0.01f;

		PaintDryAreasVelocityDampForceScale = 50.f;

		DrySurfaceDampThreshold = 1.f;

		//NormalsBlendRadius = 4.f;
		NormalsBlendRadius = 9.f;

		//SurfaceDepthOffset = 0.2f;
		SurfaceDepthOffset = 0.f;

		SDFShapeBlendMode = SDF_BLEND_EXP;
		SDFNormalsBlendMode = SDF_BLEND_EXP;

		RayOriginDepthOffset = -0.5f;

		WetnessAccumulationAmount = 5.f;
		WetnessDryingOutAmount = 1.f;

		WetnessLimit = 4.f;

		PaintDiffusionFactor = 1.f;

		WetnessContributionToHeight = 0.2f;

		#if PAINT_RENDER_ADVECT_SOURCE
		WetnessContributionToHeight = 0.9f;
		NormalsSmoothness = 0.1f;
		#endif

		FalingParticleStartPosZ = 0.1f;
		//FallingParticleSpeedTowardsSurface = 100.f;
		FallingParticleSpeedTowardsSurface = 30.f;


		ParticleSourceTextureSize = 64;

		SurfaceHitVelocityScale = 5.f;
		SurfaceHitVelocityMax = 50.f;

		ParticleBumpinnesPow = 10.f;

	}
#endif
};

struct PaintDiffusionPerPassConstants
{
	uint2 LocalOffset;
};












//===============================================
//
//				Paint Visualizer
//
//===============================================



#define PBR_RENDER 1


struct 
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
PaintRenderConstants 
{
	uint bShadeLambert;
	uint bShadeSpecular;
	uint bShadeFresnel;
	uint SpecularShadingMethod;

	float SpecularPower;
	float SpecularAmount;
	float2 pad;

	float3 GlobalLightDirection;
	uint NumLights;

	matrix PaintSceneViewMatrix;
	matrix PaintSceneProjectionMatrix;
	matrix PaintSceneViewProjectionMatrix;

	matrix GlobalLightViewProjectionMatrix;
	matrix GlobalLightViewMatrix;

	float ShadowDepthBias;
	bool bSlopeCorrectShadowBias;
	float UniformRoughness;
	float UniformMetalness;

	float2 RoughnessMinMax;
	float2 MetalnessMinMax;

	float PaintBrightnessScale;
	float AmbientLight;
	bool bVisualizeNormals;

#ifdef CPP_SHADER
	PaintRenderConstants()
	{
		bShadeLambert = true;
		bShadeSpecular = true;
		bShadeFresnel = false;
		SpecularShadingMethod = 1;

		SpecularPower = 512.f;
		SpecularAmount = 1.f;

		NumLights = 2;

		ShadowDepthBias = 0.005f;
		bSlopeCorrectShadowBias = true;

		UniformRoughness = 0.1f;
		UniformMetalness = 0.f;

		RoughnessMinMax = float2{0, 1};
		MetalnessMinMax = float2{0, 1};

		PaintBrightnessScale = 1.1f;

		AmbientLight = 0.f;

		bVisualizeNormals = false;


	}
#endif
};

#define LIGHT_TYPE_POINT 0
#define LIGHT_TYPE_SPOT 1


struct LightTransformStruct
{
	matrix View;
	matrix ViewProjection;
};


struct 
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
LightDescription
{
	float3 Position;
	float Radius;

	float3 Color;
	uint Type;

	float3 Direction;
	float Intensity;

	float2 SpotlightAngles;
	float2 pad;

#ifdef CPP_SHADER
	LightDescription()
	{
		Position = float3(10, 10, 10);
		Radius = 40.f;

		Color = float3(1,1,1);

		Type = LIGHT_TYPE_POINT;

		Intensity = 0.75f;

		Direction = float3( 0, -1, 0 );
		SpotlightAngles = float2{ 0.7f,0 };
	}
#endif

};













#endif//_GPU_PARTICLES_COMMON_STRUCTS_H
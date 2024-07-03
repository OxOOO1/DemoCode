#ifndef _GPU_BURNER_COMMON_STRUCTS_H
#define _GPU_BURNER_COMMON_STRUCTS_H

#define BURNER_MAX_PARTICLES 250000u

#define PARTICLE_STATE_DEAD 0
#define PARTICLE_STATE_ALIVE 1


#define BLENDING_DISABLE 0
#define BLENDING_ADDITIVE 1
#define BLENDING_ADDITIVE_ALPHA 2
#define BLENDING_MODULATE 3
#define BLENDING_SUBTRACT 4
#define BLENDING_MAX 5
#define BLENDING_ALPHA 6
#define BLENDING_SUBTRACT_SOURCE 7
#define BLENDING_MIN 8
#define BLENDING_ALPHA_ADD_ALPHA 9

struct
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
CGParticleConstants
{
	float3 ExternalForce;
	uint NumParticlesToSpawn;

	float DeltaTime;
	float3 ParticleSpawnPosition;
	
	float3 ParticleSpawnVelocity;
	float CurTimeSec;

	float2 ParticleSize;
	float RandVelocityScale;
	float UpwardForceScale;

	float2 FlipBookUVScale;
	float NoiseVelocityScale;
	float ParticleLifetimeSec;

	uint2 NumSpawners2D;
	float2 DiminishedParticleSizeMin;

	float2 SpawnerFireStateMinMax;
	float ParticleSizeScalePow;
	float ParticleColorBlendPow;

	float ParticleBrightness;
	float ArtificialColorAmount;
	float ParticleSizeScaleAdditional;
	float BuoyancyForceScale;

	int ParticleAnimationNumLoops; 
	float AlphaFadeInPow;
	float AlphaFadeOutPow;
	int BlendingMode;

	int bPreMultiplyAlpha;
	float PreMultiplyAlphaScalePow;
	float RuntimeInterpolatorParameter;
	float ParticleSpawnInterval;

	float2 FadeInParticleSizeMin;
	float2 LowTempParticleSizeLimiter;

	float ParticleBrightness2Lerp;
	float AlphaScale;
	float2 SimulationDomainLengthWorldSpace;

	float3 ArtificialColor;
	float MinRandomSize;

	float RenderToTextureSizeScale;
	float RenderToTextureIntensity;
	float SpawnerFireStateMax;

	matrix SceneViewProjectionMatrix;

#ifdef CPP_SHADER
	CGParticleConstants()
	{
		ExternalForce = float3{ 0,100,0 };

		NumParticlesToSpawn = 1;
		DeltaTime = 1.f / 60.f;

		ParticleSpawnPosition = float3{0,0,0};
		ParticleSpawnVelocity = float3{0,1,0};
		
		ParticleSize = float2{1.f, 1.f};
		
		CurTimeSec = 0.f;

		RandVelocityScale = 2.5f;
		//RandVelocityScale = 0.f;

		UpwardForceScale = 10.f;
		//UpwardForceScale = 0.f;

		FlipBookUVScale = float2{16, 4};

		//ParticleLifetimeSec = 3.f;
		ParticleLifetimeSec = 1.f;

		NoiseVelocityScale = 40.f;

		NumSpawners2D = uint2{16,16};

		DiminishedParticleSizeMin = float2{0, 0.2};

		SpawnerFireStateMinMax = float2{0.5, 10.f};

		ParticleSizeScalePow = 8.f;
		ParticleColorBlendPow = 2.f;

		ParticleBrightness = 1.f;
		ParticleBrightness2Lerp = 1.f;

		ArtificialColorAmount = 0.5f;

		ParticleSizeScaleAdditional = 1.f;

		BuoyancyForceScale = 2.5f;

		AlphaFadeInPow = 4.f;
		AlphaFadeOutPow = 2.f;

		ParticleAnimationNumLoops = 3;

		BlendingMode = BLENDING_ALPHA;

		bPreMultiplyAlpha = 0;

		PreMultiplyAlphaScalePow = 4.f;

		RuntimeInterpolatorParameter = 0.f;

		ParticleSpawnInterval = 1.f;

		FadeInParticleSizeMin = float2{0,0};

		LowTempParticleSizeLimiter = float2{5, 1};

		AlphaScale = 1.f;

		ArtificialColor = float3(0.2, 0.2, 1.f);

		MinRandomSize = 1.5f;

		RenderToTextureSizeScale = 0.01;
		RenderToTextureIntensity = 0.25;

		SpawnerFireStateMax = 100.f;

	}
#endif//CPP_SHADER
};










struct
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
CGPlaneConstants
{
	float DeltaTime;
	float FireSpreadSpeed;
	float FuelConsumeSpeed;
	float FireDissipationSpeed;

	//Fire Applier
	uint2 FireApplierPositionPixelSpace;
	float FireApplierRadius;
	float FireApplierStrength;

	float NoiseTextureInterpolator;
	float NoiseAdvectedSpreadStrength;
	float MaxFireThreshold;
	float PlaneLengthWorldSpace;

	uint2 BurnerPlaneSizePixels;
	float VisualizerFireScale;
	float VisualizerFirePow;

	float Time;
	float VirtualLightRadius;
	float VirtualLightFlickerSpeed;
	float VirtualLightIntensity;

	float BurnedImageSharpness;
	float AshColorScale;

	matrix SceneViewProjectionMatrix;

#ifdef CPP_SHADER
	CGPlaneConstants()
	{
		DeltaTime = 1.f / 60.f;

		FireSpreadSpeed = 10.f;

		FuelConsumeSpeed = 0.75f;

		FireDissipationSpeed = 0.25f;



		FireApplierPositionPixelSpace = uint2{100,100};
		FireApplierRadius = 5.f;
		FireApplierStrength = 10.f;

		NoiseTextureInterpolator = 0.f;

		NoiseAdvectedSpreadStrength = 0.75f;

		MaxFireThreshold = 10.f;

		VisualizerFireScale = 0.75f;
		VisualizerFirePow = 0.9f;

		PlaneLengthWorldSpace = 20.f;

		VirtualLightRadius = 15.f;
		VirtualLightFlickerSpeed = 0.2f;
		VirtualLightIntensity = 0.25f;

		BurnedImageSharpness = 0.1f;

		AshColorScale = 1.f;

	}
#endif//CPP_SHADER
};




#define SMOKE_BLOOM_SEPARATE 1


struct
#ifdef CPP_SHADER
	alignas(16) 
#endif//CPP_SHADER
CGCombinerConstants
{
	float BloomCurvePow;
	float BloomThreshold;
	float BloomStrength;
	float PlaneBloomStrength;

	float2 SmokeBloomColorClampMinMax;
	float2 SmokeBloomAlphaClampMinMax;

	float NoiseTextureInterpolator;
	float ExposureFinal;
	float Time;
	float CharcoalIntensity;

	float BurnedImageIntensity;



#ifdef CPP_SHADER
	CGCombinerConstants()
	{
		BloomCurvePow = 2.f;
		BloomThreshold = 0.4f;

		BloomStrength = 2.f;

		PlaneBloomStrength = 0.05f;

		ExposureFinal = 1.f;

		SmokeBloomColorClampMinMax = float2{0.15, 1};
		SmokeBloomAlphaClampMinMax = float2{0,1};

		NoiseTextureInterpolator = 0.f;

		CharcoalIntensity = 5.f;
		BurnedImageIntensity = 2.f;

	}
#endif//CPP_SHADER
};










#endif//_GPU_BURNER_COMMON_STRUCTS_H
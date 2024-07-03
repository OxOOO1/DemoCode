
#include "Physics\BurnerCommonStructs.h"

#include "Physics\PhysicsCommon.hlsli"

#include "InterpolationEtc\InterpolationCommon.hlsli"
/*
*	Resources
*/

SamplerState LinearSampler : register(s0);
SamplerState PointSampler : register(s1);
SamplerState LinearWrapSampler : register(s2);

cbuffer PassInfo : register(b10)
{
    uint2 GTextureSize;
};


//CB

ConstantBuffer<CGParticleConstants> CGParticleConstantsCB : register(b0);

//SRV
StructuredBuffer<float3> ParticlesBufferPositionCurSRV : register(t0);
StructuredBuffer<float3> ParticlesBufferVelocityCurSRV : register(t1);
StructuredBuffer<uint> ParticleCountBufferSRV : register(t2);

StructuredBuffer<float3> ParticleVertexBufferSRV : register (t3);
StructuredBuffer<float2> ParticleTexCoordsBufferSRV : register (t4);

Texture1D<float2> RandomVelocityTextureSRV : register (t5);

Texture2D<float4> ParticleColorTextureSRV : register (t6);

StructuredBuffer<uint> DeadParticlesIndicesListBufferSRV : register (t7);
StructuredBuffer<uint> DeadParticleCountBufferSRV : register(t8);

StructuredBuffer<float> ParticleAgeBufferSRV : register(t9);

StructuredBuffer<uint> ParticleStateBufferSRV : register(t10);

Texture2D<float2> RandomVelocityTexture2DSRV : register (t11);

StructuredBuffer<float2> SpawnersPositionsBufferSRV : register(t12);

Texture2D<float> PlaneFireTextureSRV : register (t13);

Texture2D<float4> SpawnerStateTextureSRV : register (t14);

StructuredBuffer<float> ParticleTemperatureBufferSRV : register(t15);

Texture2D<float4> SpawnersNoiseTextureSRV : register (t16);

StructuredBuffer<float> ParticleSizeBufferSRV : register(t17);

//UAV
RWStructuredBuffer<float3> ParticlesBufferPositionCurUAV : register (u0);
RWStructuredBuffer<float3> ParticlesBufferVelocityCurUAV : register (u1);
RWStructuredBuffer<uint> ParticleCountBufferUAV : register (u2);

RWStructuredBuffer<uint> DeadParticlesIndicesListBufferUAV : register (u3);
RWStructuredBuffer<uint> DeadParticleCountBufferUAV : register (u4);

RWStructuredBuffer<float> ParticleAgeBufferUAV : register (u5);

RWStructuredBuffer<uint> ParticleStateBufferUAV : register (u6);

RWTexture2D<float4> SpawnerStateTextureUAV : register (u7);

RWStructuredBuffer<float> ParticleTemperatureBufferUAV : register (u8);

RWStructuredBuffer<float> ParticleSizeBufferUAV : register (u9);

RWStructuredBuffer<float> SpawnersAgeBufferUAV : register (u10);

[numthreads( 64, 1, 1 )]
void FillDeadParticlesIndicesList( 
	uint SampleCoord : SV_DispatchThreadID )
{
	if(SampleCoord < BURNER_MAX_PARTICLES)
	{
		if(SampleCoord == 0)
		{
			DeadParticleCountBufferUAV[0] = BURNER_MAX_PARTICLES;
		}
		DeadParticlesIndicesListBufferUAV[SampleCoord] = BURNER_MAX_PARTICLES - 1 - SampleCoord;

		ParticleAgeBufferUAV[SampleCoord] = 0.f;
		ParticleStateBufferUAV[SampleCoord] = PARTICLE_STATE_DEAD;
	}

}


uint2 GetSpawnerIndex2D(uint flatIndex)
{
	uint2 SpawnerIndex2D;
	SpawnerIndex2D.y = (flatIndex % (CGParticleConstantsCB.NumSpawners2D.y));
	SpawnerIndex2D.x = (flatIndex / (CGParticleConstantsCB.NumSpawners2D.x));
	return SpawnerIndex2D;
}

[numthreads( 8, 8, 1 )]
void FillSpawnStateTexture( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float2 uv = float2(SampleCoord.xy) / float2(CGParticleConstantsCB.NumSpawners2D.xy);

	//float curFire = PlaneFireTextureSRV.SampleLevel(LinearSampler, uv, 0);
	float curFire = PlaneFireTextureSRV.SampleLevel(PointSampler, uv, 0);

	float4 curSpawnerState = SpawnerStateTextureUAV[SampleCoord];

	curSpawnerState.r = curFire;

	SpawnerStateTextureUAV[SampleCoord] = curFire;

}


float GetSpawnersNoise(uint2 sampleCoord)
{	

	float2 sampleCoordf = float2(sampleCoord) + 0.5f;
	float2 uv = sampleCoordf / float2(CGParticleConstantsCB.NumSpawners2D);
	uv.y += CGParticleConstantsCB.CurTimeSec * 0.001;
	float4 noise = SpawnersNoiseTextureSRV.SampleLevel(LinearWrapSampler, uv, 0);

	//float4 noise = SpawnersNoiseTextureSRV[sampleCoord];

	float finalNoise;

	if(CGParticleConstantsCB.RuntimeInterpolatorParameter < 1.f)
	{
		finalNoise = lerp(noise.x, noise.y, CGParticleConstantsCB.RuntimeInterpolatorParameter);
	}
	else if(CGParticleConstantsCB.RuntimeInterpolatorParameter < 2.f)
	{
		finalNoise = lerp(noise.y, noise.z, frac(CGParticleConstantsCB.RuntimeInterpolatorParameter));
	}
	else if(CGParticleConstantsCB.RuntimeInterpolatorParameter < 3.f)
	{
		finalNoise = lerp(noise.z, noise.w, frac(CGParticleConstantsCB.RuntimeInterpolatorParameter));
	}
	else
	{
		finalNoise = lerp(noise.w, noise.x, frac(CGParticleConstantsCB.RuntimeInterpolatorParameter));
	}

	return finalNoise;
}

[numthreads( 64, 1, 1 )]
void SpawnParticles( 
	uint SampleCoord : SV_DispatchThreadID )
{
	//if(SampleCoord < CGParticleConstantsCB.NumParticlesToSpawn)
	const uint totalNumSpawners = CGParticleConstantsCB.NumSpawners2D.x * CGParticleConstantsCB.NumSpawners2D.y;
	if(SampleCoord < totalNumSpawners)
	{

		const float2 spawnPos = SpawnersPositionsBufferSRV[SampleCoord];

		#if 1//SPAWNER STATE
			//get 2d coord
			uint2 sampleCoord2D = GetSpawnerIndex2D(SampleCoord);
			sampleCoord2D.y = CGParticleConstantsCB.NumSpawners2D.y - 1 - sampleCoord2D.y;
			
			#if 1//FIRE FROM TEX
			float2 fireTextureUV = saturate(MapToRange(spawnPos.xy, 0, CGParticleConstantsCB.SimulationDomainLengthWorldSpace.x, 0, 1));
			fireTextureUV.y = 1.f - fireTextureUV.y;
			const float curFireState = PlaneFireTextureSRV.SampleLevel(LinearSampler, fireTextureUV, 0);
			#else
			float curFireState = SpawnerStateTextureSRV[sampleCoord2D];
			#endif
			if(curFireState < CGParticleConstantsCB.SpawnerFireStateMinMax.x || curFireState > CGParticleConstantsCB.SpawnerFireStateMax)
			{
				return;
			}
		#endif

		
		float spawnerNoiseValue = GetSpawnersNoise(sampleCoord2D);

		{
			float curSpawnerAge = SpawnersAgeBufferUAV[SampleCoord];
			const bool bForceSpawn = curSpawnerAge == 0.f;

			const float spawnAgeIncrement = CGParticleConstantsCB.DeltaTime + clamp(spawnerNoiseValue * 0.001, -CGParticleConstantsCB.DeltaTime * 0.5f, CGParticleConstantsCB.DeltaTime);

			if(bForceSpawn)//activate our spawner and let it pass
			{
				curSpawnerAge += spawnAgeIncrement;
				SpawnersAgeBufferUAV[SampleCoord] = curSpawnerAge;
			}
			else
			{
				if (curSpawnerAge < CGParticleConstantsCB.ParticleSpawnInterval)//we still aren't allowed to spawn
				{
					curSpawnerAge += spawnAgeIncrement;
					SpawnersAgeBufferUAV[SampleCoord] = curSpawnerAge;
					return;
				}
				else
				{
					//reset timer and allow to spawn
					//curSpawnerAge = 0.f;
					curSpawnerAge = spawnAgeIncrement;
					SpawnersAgeBufferUAV[SampleCoord] = curSpawnerAge;
				}
			}
		}
	


		#if 1//RANDOMLY DISABLE SPAWNERS
		if(spawnerNoiseValue <= 0.f)
		{
			return;
		}
		#endif

		//decrement dead particles list, get next available particle index
		uint indexToDeadList;
		InterlockedAdd(DeadParticleCountBufferUAV[0], -1, indexToDeadList);

		if(indexToDeadList < 1 || indexToDeadList > BURNER_MAX_PARTICLES)
		{
			//no particles cur available
			return;
		}

		indexToDeadList -= 1;
		uint outParticleIndex = DeadParticlesIndicesListBufferSRV[indexToDeadList];
	
		ParticleSizeBufferUAV[outParticleIndex] = max(CGParticleConstantsCB.MinRandomSize, spawnerNoiseValue);
		//ParticleSizeBufferUAV[outParticleIndex] = max(CGParticleConstantsCB.MinRandomSize, spawnerNoiseValue * spawnerNoiseValue);
		
		ParticleAgeBufferUAV[outParticleIndex] = 0.f;
		ParticleStateBufferUAV[outParticleIndex] = PARTICLE_STATE_ALIVE;

		ParticleTemperatureBufferUAV[outParticleIndex] = curFireState;

#if 0//DEFAULT VELOCITY
		ParticlesBufferVelocityCurUAV[outParticleIndex] = CGParticleConstantsCB.ParticleSpawnVelocity;
#else
		//Get Rand Velocity based on curTime and threadId
		float2 randVel = RandomVelocityTextureSRV.SampleLevel(LinearWrapSampler, CGParticleConstantsCB.CurTimeSec + 0.1f * SampleCoord, 0);
		ParticlesBufferVelocityCurUAV[outParticleIndex] = float3(randVel.xy * CGParticleConstantsCB.RandVelocityScale, 0);
#endif

#if 1//get spawn position
		ParticlesBufferPositionCurUAV[outParticleIndex] = float3(spawnPos.xy, CGParticleConstantsCB.ParticleSpawnPosition.z);
#else
		ParticlesBufferPositionCurUAV[outParticleIndex] = CGParticleConstantsCB.ParticleSpawnPosition;
#endif
	}

}

[numthreads( 64, 1, 1 )]
void UpdateParticles( 
	uint SampleCoord : SV_DispatchThreadID )
{
	//const uint particleCount = ParticleCountBufferSRV[0];

	uint curParticleState = ParticleStateBufferUAV[SampleCoord];

	if(curParticleState > PARTICLE_STATE_DEAD)
	{
		//check if particle is aged
		float curParticleAge = ParticleAgeBufferUAV[SampleCoord];
		if(curParticleAge < CGParticleConstantsCB.ParticleLifetimeSec)
		{
			//update alive particle

			//curParticleAge += CGParticleConstantsCB.DeltaTime;//from temp

			const float particleAgeNormalised = curParticleAge / CGParticleConstantsCB.ParticleLifetimeSec;

			float3 curPos = ParticlesBufferPositionCurUAV[SampleCoord];
			float3 curVel = ParticlesBufferVelocityCurUAV[SampleCoord];

			const float curParticleSize = ParticleSizeBufferSRV[SampleCoord];
			

			const float damping = 0.99f;
			curVel = curVel * damping;

		#if 1 //noise velocity
			//float2 randVel = RandomVelocityTexture2DSRV.SampleLevel(LinearWrapSampler, curPos.xy + CGParticleConstantsCB.CurTimeSec, 0);

			float2 randVelUV = curPos.xy / CGParticleConstantsCB.SimulationDomainLengthWorldSpace;
			randVelUV.y += CGParticleConstantsCB.CurTimeSec * 0.1;
			randVelUV *= 0.1;
			float2 randVel = RandomVelocityTexture2DSRV.SampleLevel(LinearWrapSampler, randVelUV, 0);
			/* 
			float2 randVel = RandomVelocityTexture2DSRV.SampleLevel(LinearWrapSampler, randVelUV, 0); */
			#if 1//PERLIN
			randVel = randVel * 2.f - 1.f;
			randVel = normalize(randVel);
			#endif

			float velocityAdditionalScale = 1.f;
			#if 1//APPLY MORE VELOCITY TO OLD PARTICLES
			velocityAdditionalScale *= clamp(particleAgeNormalised, 0.5, 1.f);
			#endif

			curVel.x += randVel.x * velocityAdditionalScale * CGParticleConstantsCB.NoiseVelocityScale * CGParticleConstantsCB.DeltaTime;
			curVel.y += randVel.y * velocityAdditionalScale * CGParticleConstantsCB.NoiseVelocityScale * CGParticleConstantsCB.DeltaTime;

		#endif

			float buoyancyForce = 0.f;

		#if 1
			//hot particles are faster
			//float curParticleTemp = ParticleTemperatureBufferUAV[SampleCoord];
			float curParticleTemp = ParticleTemperatureBufferUAV[SampleCoord];

			//sample fire texture for cur temp
			float2 fireTextureUV = MapToRange(curPos.xy, 0, CGParticleConstantsCB.SimulationDomainLengthWorldSpace.x, 0, 1);
			fireTextureUV.x = saturate(fireTextureUV.x);

			if(any(fireTextureUV > float2(1.0, 1.0)))
			{
				curParticleTemp *= 0.95f;
				
			}
			else
			{
				fireTextureUV.y = 1.f - fireTextureUV.y;
				const float particleTempFromTexture = PlaneFireTextureSRV.SampleLevel(LinearSampler, fireTextureUV, 0);

				curParticleTemp += ((particleTempFromTexture - curParticleTemp) * CGParticleConstantsCB.DeltaTime * 0.1f);

				if(curParticleTemp < 0.01)
				{
					//curParticleAge += 100.f;
				}

				
			}

			curParticleTemp = max(CGParticleConstantsCB.SpawnerFireStateMinMax.x, curParticleTemp);
			ParticleTemperatureBufferUAV[SampleCoord] = curParticleTemp;


			const float SlowParticleScale = 0.9f;
			const float FastParticleScale = 1.5f;

			const float particleTempNorm = MapToRange(curParticleTemp, 0, 10, SlowParticleScale, FastParticleScale);

			const float particleSizeNorm = max(1.f, MapToRange(CGParticleConstantsCB.ParticleSizeScaleAdditional * curParticleSize, 2, 10, 1, 2));

			curParticleAge += CGParticleConstantsCB.DeltaTime * particleTempNorm * particleSizeNorm;


			curParticleTemp = saturate(MapToRange(curParticleTemp, CGParticleConstantsCB.SpawnerFireStateMinMax.x, CGParticleConstantsCB.SpawnerFireStateMinMax.y, 0, 1));
			buoyancyForce = curParticleTemp * CGParticleConstantsCB.BuoyancyForceScale;
		#else
			curParticleAge += CGParticleConstantsCB.DeltaTime;
		#endif

			const float3 acceleration = float3(0, CGParticleConstantsCB.UpwardForceScale + buoyancyForce, 0);
			curVel = curVel + acceleration * CGParticleConstantsCB.DeltaTime;

			curPos = curPos + curVel * CGParticleConstantsCB.DeltaTime;

			ParticlesBufferVelocityCurUAV[SampleCoord] = curVel;
			ParticlesBufferPositionCurUAV[SampleCoord] = curPos;

			ParticleAgeBufferUAV[SampleCoord] = curParticleAge;
		}
		else
		{
		#if 1
			//mark particle as dead
			ParticleAgeBufferUAV[SampleCoord] = 0.f;
			ParticleStateBufferUAV[SampleCoord] = PARTICLE_STATE_DEAD;

			#if 1
			//add to the dead particles list
			uint indexToDeadList;
			InterlockedAdd(DeadParticleCountBufferUAV[0], 1, indexToDeadList);
			//indexToDeadList += 1;
			DeadParticlesIndicesListBufferUAV[indexToDeadList] = SampleCoord;
			#endif
		#endif
		}
	}
	

}





//Render Particles


struct ParticleRenderInterpolators
{
	float4 Position : SV_Position;
	float2 TexCoords : TEXCOORD;
	nointerpolation float ParticleAgeNormalised : AGE;
	nointerpolation float ParticleTempNormalised : TEMP;
	nointerpolation float FlipBookIndex : FLIPINDEX;
	nointerpolation float AlphaScale : SCALE;  
};

void ParticleRenderVS(
	in uint VertexIndex : SV_VertexID,
	in uint InstanceID : SV_InstanceID,
	out ParticleRenderInterpolators vsout
	)
{

	uint curParticleState = ParticleStateBufferSRV[InstanceID];

	if(curParticleState == PARTICLE_STATE_DEAD)
	{
		vsout.Position = float4(0,0,0,1);
		vsout.TexCoords = float2(0,0);
	} 
	else
	{

		/*
		*	Age etc
		*/
		const float curParticleAge = min(CGParticleConstantsCB.ParticleLifetimeSec, ParticleAgeBufferSRV[InstanceID]);


#if 1//LOOPS
		vsout.ParticleAgeNormalised = saturate(curParticleAge / CGParticleConstantsCB.ParticleLifetimeSec);

		float ParticleAnimationParameterNormalised;
		
		if(CGParticleConstantsCB.ParticleAnimationNumLoops > 1)
		{
			ParticleAnimationParameterNormalised = frac(curParticleAge / (CGParticleConstantsCB.ParticleLifetimeSec / float(CGParticleConstantsCB.ParticleAnimationNumLoops)));
		}
		else
		{
			ParticleAnimationParameterNormalised = vsout.ParticleAgeNormalised;
		}
		

		//freeze

#else
		vsout.ParticleAgeNormalised = curParticleAge / CGParticleConstantsCB.ParticleLifetimeSec;
		const float ParticleAnimationParameterNormalised = vsout.ParticleAgeNormalised;
#endif


		const float curParticleTemp = ParticleTemperatureBufferSRV[InstanceID];
		vsout.ParticleTempNormalised = saturate(MapToRange(curParticleTemp, CGParticleConstantsCB.SpawnerFireStateMinMax.x, CGParticleConstantsCB.SpawnerFireStateMinMax.y, 0, 1));

		const float TotalFlipFrames = CGParticleConstantsCB.FlipBookUVScale.x * CGParticleConstantsCB.FlipBookUVScale.y;
		float FrameIndexf = (ParticleAnimationParameterNormalised * TotalFlipFrames);
		vsout.FlipBookIndex = FrameIndexf;

		const float curParticleSize = ParticleSizeBufferSRV[InstanceID];

		/*
		*	Position
		*/
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[InstanceID];
		float3 positionWorld = ParticleVertexBufferSRV[VertexIndex];	
		positionWorld.x *= CGParticleConstantsCB.ParticleSize.x * CGParticleConstantsCB.ParticleSizeScaleAdditional;
		positionWorld.y *= CGParticleConstantsCB.ParticleSize.y * CGParticleConstantsCB.ParticleSizeScaleAdditional;
		
		#if 0
		positionWorld.x *= curParticleTemp < 1 ? min(curParticleSize, CGParticleConstantsCB.LowTempParticleSizeLimiter.x) : curParticleSize;
		positionWorld.y *= curParticleTemp < 1 ? min(curParticleSize, CGParticleConstantsCB.LowTempParticleSizeLimiter.y) : curParticleSize;
		#else
		positionWorld.xy *= curParticleSize;
		#endif


		vsout.AlphaScale = clamp(MapToRange(curParticlePos.y, 0, CGParticleConstantsCB.SimulationDomainLengthWorldSpace.y, 0, 1), 0.4, 1);

#if 1//FADE IN
		//float normAgePow = vsout.ParticleAgeNormalised * (1.0 - vsout.ParticleAgeNormalised) * (1.0 - vsout.ParticleAgeNormalised) * 6.743;
		//float normAgePow = pow(vsout.ParticleAgeNormalised, CGParticleConstantsCB.ParticleSizeScalePow);
		float normAgeFadeInPow = sqrt(1 - pow(1.f - vsout.ParticleAgeNormalised, CGParticleConstantsCB.ParticleSizeScalePow));
		positionWorld.x *= clamp(normAgeFadeInPow, CGParticleConstantsCB.FadeInParticleSizeMin.x, 1.f);
		positionWorld.y *= clamp(normAgeFadeInPow, CGParticleConstantsCB.FadeInParticleSizeMin.y, 1.f);
#endif

#if 1//FADE OUT
		float normAgePow = pow(vsout.ParticleAgeNormalised, CGParticleConstantsCB.ParticleSizeScalePow);
		positionWorld.x *= clamp(1.f - normAgePow, CGParticleConstantsCB.DiminishedParticleSizeMin.x, 1.f);
		positionWorld.y *= clamp(1.f - normAgePow, CGParticleConstantsCB.DiminishedParticleSizeMin.y, 1.f);
#endif

		//translate
		positionWorld += curParticlePos;	
		vsout.Position = mul(CGParticleConstantsCB.SceneViewProjectionMatrix, float4(positionWorld, 1.0f));

		/*
		*	TexCoord
		*/

		float2 uvLocal = ParticleTexCoordsBufferSRV[VertexIndex];

		#if 1//flip uv's
		if(InstanceID % 3 == 0)
		{
			//uvLocal = float2(1.f - uvLocal.x, 1.f - uvLocal.y);
			uvLocal.x = 1.f - uvLocal.x;
			//uvLocal.y = 1.f - uvLocal.y;
		}
		#endif

		vsout.TexCoords = uvLocal;


		

	}

	
}

float4 ParticleRenderPS
	(
	in ParticleRenderInterpolators psin
	) : SV_Target
{
	uint flipBookIndex1D = floor(psin.FlipBookIndex);

	uint2 FlipBookIndex2D;
	FlipBookIndex2D.x = (flipBookIndex1D % int(CGParticleConstantsCB.FlipBookUVScale.x));
	FlipBookIndex2D.y = (flipBookIndex1D / int(CGParticleConstantsCB.FlipBookUVScale.x));

	float2 frameSize = rcp(CGParticleConstantsCB.FlipBookUVScale);

	//Scale texCoords for flipbook
	float2 uv = psin.TexCoords * frameSize;

	//offset UV based on cur frame index
	//uv += frameSize * FlipBookIndex2D;
	uv.x += frameSize.x * FlipBookIndex2D.x;
	uv.y += frameSize.y * FlipBookIndex2D.y;

	float4 color = ParticleColorTextureSRV.SampleLevel(LinearSampler, uv, 0);

	//DEBUG TEMP
	//return float4(psin.ParticleTempNormalised, 0, 0, color.a);


	//Premultiply Alpha
	//color.rgb *= color.a;

	#if 1//FADE IN
	color.a *= sqrt(1 - pow(1.f - psin.ParticleAgeNormalised, CGParticleConstantsCB.AlphaFadeInPow));
	//color.a *= psin.ParticleAgeNormalised;

	#endif

	#if 1//FADE OUT
	color.a *= (1.f - pow(psin.ParticleAgeNormalised, CGParticleConstantsCB.AlphaFadeOutPow));
	#endif

	color.a *= CGParticleConstantsCB.AlphaScale;

	if(CGParticleConstantsCB.ParticleColorBlendPow == 0)
	{
		color.a *= psin.AlphaScale;
	}

#if 1 //artificial color
	if(CGParticleConstantsCB.ParticleColorBlendPow > 0)
	{
		const float3 artColor1 = CGParticleConstantsCB.ArtificialColor;
		//float3 artColor1 = float3(1, 0.5, 0.1f);
		float3 artColor2 = float3(1, 1, 1.f);
		//float3 artColor2 = float3(1, 0.5, 0.1f);

		float powFromTemp = clamp(1.f - psin.ParticleTempNormalised, 0.1, 1) * 2.f;
		//float powFinal = powFromTemp;
		float powFinal = CGParticleConstantsCB.ParticleColorBlendPow;

		//Blend based on Temperature
		//float3 finalColor = lerp(artColor1, artColor2, pow((psin.ParticleTempNormalised), CGParticleConstantsCB.ParticleColorBlendPow));
		
		color.rgb = lerp(saturate(artColor1 * color.a * 4.f), color.rgb, pow((psin.ParticleTempNormalised), CGParticleConstantsCB.ParticleColorBlendPow));
		//color.rgb = lerp(artColor1 * color.a, color.rgb, pow(saturate(psin.ParticleAgeNormalised + psin.ParticleTempNormalised), CGParticleConstantsCB.ParticleColorBlendPow));

		//for smoke
		//color.rgb = lerp(artColor1 * color.rgb * color.a * 4.f, color.rgb, saturate(psin.ParticleAgeNormalised + 1.f - pow(psin.ParticleTempNormalised, CGParticleConstantsCB.ParticleColorBlendPow)));

		//Blend based on Age
		//float3 finalColor = lerp(artColor1, artColor2, pow(psin.ParticleAgeNormalised, powFinal));

		//Blend based on Age & Temp
		//float3 finalColor = lerp(artColor1, artColor2, pow(saturate(psin.ParticleAgeNormalised + psin.ParticleTempNormalised), powFinal));

		//Mix with original color
		//color.rgb = lerp(color.rgb, finalColor * color.a, CGParticleConstantsCB.ArtificialColorAmount);
		//color.rgb = lerp(color.rgb, finalColor * color.rgb, CGParticleConstantsCB.ArtificialColorAmount);
		//color.rgb = lerp(color.rgb, color.rgb * finalColor, CGParticleConstantsCB.ArtificialColorAmount);
		//color.rgb = lerp(color.rgb, color.rgb * finalColor, min(CGParticleConstantsCB.ArtificialColorAmount, 1.f - psin.ParticleTempNormalised));
	}
	
#endif

	if(CGParticleConstantsCB.bPreMultiplyAlpha)
	{
		color.rgb *= sqrt(1.f - pow(1.f - color.a, CGParticleConstantsCB.PreMultiplyAlphaScalePow));
		//color.rgb *= color.a;
		//color.rgb *= color.a > 0.1 ? 1.f : 0.f;
	}

	if(CGParticleConstantsCB.ParticleBrightness2Lerp != CGParticleConstantsCB.ParticleBrightness)
	{
		color.rgb *= lerp(CGParticleConstantsCB.ParticleBrightness, CGParticleConstantsCB.ParticleBrightness2Lerp, psin.ParticleAgeNormalised);
	}
	else
	{
		color.rgb *= CGParticleConstantsCB.ParticleBrightness;
	}

	/* color.a = (saturate(color.r)) * color.a;
	color.rgb *= CGParticleConstantsCB.ParticleBrightness; */


	//color.rgb *= psin.ParticleTempNormalised;

	return color;


}










struct ParticleRenderToTextureInterpolators
{
	float4 Position : SV_Position;
	float2 TexCoords : TEXCOORD;
};

void ParticleRenderToTextureVS(
	in uint VertexIndex : SV_VertexID,
	in uint InstanceID : SV_InstanceID,
	out ParticleRenderToTextureInterpolators vsout
	)
{

	uint curParticleState = ParticleStateBufferSRV[InstanceID];

	if(curParticleState == PARTICLE_STATE_DEAD)
	{
		vsout.Position = float4(0,0,0,1);
		vsout.TexCoords = float2(0,0);
	} 
	else
	{	
		const float curParticleSize = ParticleSizeBufferSRV[InstanceID];		
		/*
		*	Position
		*/
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[InstanceID];
		float3 positionWorld = ParticleVertexBufferSRV[VertexIndex];	
		positionWorld.xy *= CGParticleConstantsCB.ParticleSizeScaleAdditional * curParticleSize * CGParticleConstantsCB.RenderToTextureSizeScale;
		positionWorld.y *= 0.5;
		//translate
		positionWorld += curParticlePos;
		vsout.Position = mul(CGParticleConstantsCB.SceneViewProjectionMatrix, float4(positionWorld, 1.0f));		
		/*
		*	TexCoord
		*/		
		float2 uvLocal = ParticleTexCoordsBufferSRV[VertexIndex];
		vsout.TexCoords = uvLocal;
	}
}

float ParticleRenderToTexturePS
	(
	in ParticleRenderToTextureInterpolators psin
	) : SV_Target
{

	float color = ParticleColorTextureSRV.SampleLevel(LinearSampler, psin.TexCoords, 0).r;

	color *= CGParticleConstantsCB.RenderToTextureIntensity;

	return color;

}








//==================================================================BURNER PLANE


/*
*	Resources
*/

//CB

ConstantBuffer<CGPlaneConstants> CGPlaneConstantsCB : register(b0);

//SRV
Texture2D<float> FireTextureSRV : register (t0);
Texture2D<float> FuelTextureSRV : register (t1);
Texture2D<float> CopyTextureSRV : register (t2);
Texture2D<float3> NoiseTextureSRV : register (t3);
Texture2D<float3> ImageTextureSRV : register (t4);
Texture2D<float3> AshesTextureSRV : register (t5);
Texture2D<float3> SpotlightTextureSRV : register (t6);
Texture2D<float> FireTextureDownsampledSRV : register (t7);

StructuredBuffer<float3> PlaneVertexBufferSRV : register (t4);
StructuredBuffer<float2> PlaneTexCoordsBufferSRV : register (t5);

//UAV
RWTexture2D<float> FireTextureUAV : register (u0);
RWTexture2D<float> FuelTextureUAV : register (u1);
RWTexture2D<float> CopyTextureUAV : register (u2);
RWTexture2D<float3> VisualizerTextureUAV : register (u3);

#define FIRE_CLAMP_TO_1 0


float GetNoise(uint2 sampleCoord)
{
	float3 noise = NoiseTextureSRV[sampleCoord];

	float finalNoise;

	if(CGPlaneConstantsCB.NoiseTextureInterpolator < 1.f)
	{
		finalNoise = lerp(noise.x, noise.y, CGPlaneConstantsCB.NoiseTextureInterpolator);
	}
	else if(CGPlaneConstantsCB.NoiseTextureInterpolator < 2.f)
	{
		finalNoise = lerp(noise.y, noise.z, frac(CGPlaneConstantsCB.NoiseTextureInterpolator));
	}
	else
	{
		finalNoise = lerp(noise.z, noise.x, frac(CGPlaneConstantsCB.NoiseTextureInterpolator));
	}

	return finalNoise;
}


[numthreads( 8, 8, 1 )]
void FireSpread( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	//get cur fire and fuel
	float curFire = FireTextureUAV[SampleCoord];
	float curFuel = FuelTextureUAV[SampleCoord];

	const float originalFireValue = curFire;

	if(curFire > CGPlaneConstantsCB.MaxFireThreshold)
	{
		//we are already burning out
		return;
	}

	if(curFuel < 0.01)
	{
		//we can't burn anymore
		return;
	}

	//read neighbor fire values

	float neighborFire;

	float weight = 1.f / 8.f;

	float accumulatedFire = 0.f;

	//right
	neighborFire = CopyTextureSRV[SampleCoord + uint2(1, 0)];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}

	if(SampleCoord.x > 0)
	{
		//left
		neighborFire = CopyTextureSRV[SampleCoord - uint2(1, 0)];
		if(neighborFire > originalFireValue)
		{
			accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
		}
	}

	if(SampleCoord.y > 0)
	{
		//up
		neighborFire = CopyTextureSRV[SampleCoord - uint2(0, 1)];
		if(neighborFire > originalFireValue)
		{
			accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
		}
	}

	//down
	neighborFire = CopyTextureSRV[SampleCoord + uint2(0, 1)];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}

#if 1//corners
	neighborFire = CopyTextureSRV[SampleCoord - uint2(1, 1)];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}

	neighborFire = CopyTextureSRV[SampleCoord + uint2(1, 1)];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}

	neighborFire = CopyTextureSRV[(SampleCoord + int2(-1, 1))];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}

	neighborFire = CopyTextureSRV[(SampleCoord + int2(1, -1))];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * weight * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}
#endif//corners


#if 1//NOISE BASED ADVECTION
	float2 noiseVec = NoiseTextureSRV[SampleCoord].xy;
	noiseVec = noiseVec * 2.f - 1.f;
	noiseVec = normalize(noiseVec) + 0.5f;
	int2 noiseVecOffset = floor(noiseVec);
	neighborFire = CopyTextureSRV[(SampleCoord + noiseVecOffset)];
	if(neighborFire > originalFireValue)
	{
		accumulatedFire += neighborFire * CGPlaneConstantsCB.NoiseAdvectedSpreadStrength * CGPlaneConstantsCB.FireSpreadSpeed * CGPlaneConstantsCB.DeltaTime;
	}
#endif

#if 1//NOISE
	float noise = GetNoise(SampleCoord);
	accumulatedFire *= noise;
#endif

	curFire += accumulatedFire;

	if(curFire > originalFireValue)
	{
		#if FIRE_CLAMP_TO_1
		curFire = saturate(curFire);
		#endif

		FireTextureUAV[SampleCoord] = curFire;
	}


}

[numthreads( 8, 8, 1 )]
void FireUpdate( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	//fire consumes fuel until no fuel left

	//get cur fire and fuel
	float curFire = FireTextureUAV[SampleCoord];
	float curFuel = FuelTextureUAV[SampleCoord];

	bool bUpdate = false;

	if(curFuel > 0)
	{
		if(curFire > 1.f)
		{
			//curFuel -= curFire * CGPlaneConstantsCB.FuelConsumeSpeed * CGPlaneConstantsCB.DeltaTime;
			curFuel -= clamp(curFire, 0.1, 0.5) * CGPlaneConstantsCB.FuelConsumeSpeed * CGPlaneConstantsCB.DeltaTime;

			//curFire += curFuel * CGPlaneConstantsCB.FuelConsumeSpeed * CGPlaneConstantsCB.DeltaTime;

			bUpdate = true;
		}
	}
	else
	{
		//fire dissipates when no fuel left
		//if(curFire > 0)
		if(curFire > 0.01)
		{
			curFire -= curFire * CGPlaneConstantsCB.FireDissipationSpeed * CGPlaneConstantsCB.DeltaTime;

			bUpdate = true;
		}
	}

	if(bUpdate)
	{
		FireTextureUAV[SampleCoord] = curFire;
		FuelTextureUAV[SampleCoord] = curFuel;
	}
	

}



[numthreads( 8, 8, 1 )]
void ApplyFire( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	float dist = length(float2(CGPlaneConstantsCB.FireApplierPositionPixelSpace) - float2(SampleCoord));

	if(dist < CGPlaneConstantsCB.FireApplierRadius)
	{
		//set the pixel on fire if we have fuel
		float curFire = FireTextureUAV[SampleCoord];
		float curFuel = FuelTextureUAV[SampleCoord];

		#if FIRE_CLAMP_TO_1
		if(curFire >= 1.f)
		{
			return;
		}
		#endif

		if(curFuel > 0)
		{
			curFire += CGPlaneConstantsCB.FireApplierStrength * CGPlaneConstantsCB.DeltaTime;
			FireTextureUAV[SampleCoord] = curFire;
		}

		//DEBUG Render Applier To Texture
		//CopyTextureUAV[SampleCoord] = 1.f;
	}

}







[numthreads( 8, 8, 1 )]
void DownsampleFireTexture( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float2 sampleCoordf = float2(SampleCoord) + 0.5;
	float2 uv = sampleCoordf / float2(GTextureSize);

	float4 source = FireTextureSRV.SampleLevel(LinearSampler, uv, 0);

	FireTextureUAV[SampleCoord] = source;
}




[numthreads( 8, 8, 1 )]
void BurnerPlaneVisualizer( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float curFire = FireTextureSRV[SampleCoord];
	float curFuel = FuelTextureSRV[SampleCoord];

	float2 uv = (float2(SampleCoord) + 0.5f) / float2(CGPlaneConstantsCB.BurnerPlaneSizePixels);
	float3 imageColor = ImageTextureSRV.SampleLevel(LinearSampler, uv, 0);
	float3 ashesColor = AshesTextureSRV.SampleLevel(LinearSampler, uv, 0) * CGPlaneConstantsCB.AshColorScale * float3(1,1,1.2);

	float3 paperColor = ashesColor;

	if(curFuel > 0)
	{
		//apply black spots to image

		float2 noises = NoiseTextureSRV.SampleLevel(LinearWrapSampler, float2(uv.x * 2.f, uv.y), 0).rg;
		float noise = max(noises.r, noises.g);
		noise = MapToRange(noise, 0.2, 0.8, 0, 1);

		float2 blackSpotUV = uv;
		blackSpotUV.y += 0.15f;
		blackSpotUV.y *= 0.9f;
		float blackSpot = CopyTextureSRV.SampleLevel(LinearSampler, blackSpotUV, 0);
		blackSpot *= max(0.25f, noise);
		blackSpot = saturate(blackSpot);
		imageColor *= (1.f - blackSpot);
		
		paperColor = lerp(ashesColor, imageColor, saturate(curFuel));

		//Lights
		float spotlight = SpotlightTextureSRV.SampleLevel(LinearSampler, uv, 0);
		paperColor *= spotlight;
	}
	else
	{
		

	}

	//Area Point Lights

	//TODO: Try replacing those with simple fake-light texture renders that we render on the scene texture

	//	!!--Or calculate Intensity/Flicker of all lights in a separate pre-pass--!!

	{
		float3 colorFinal = float3(0,0,0);

		//get cur pixel world pos
		float2 sampleCoordf = float2(SampleCoord) + 0.5;
		float2 uv = sampleCoordf / float2(GTextureSize);

		float2 posWorld = uv * CGPlaneConstantsCB.PlaneLengthWorldSpace;

		//for each light
		//get it's world pos
		//get it's intensity from fire tex
		const uint2 numLights = uint2(4,4);
		for(uint y = 0; y < numLights.y; y++)
		{
			for(uint x = 0; x < numLights.x; x++)
			{
				const float2 lightPosOffsetCell = float2(0.5,0.5);
				const uint2 curLightIndex = uint2(x,y);
				const float distBetweenLights = CGPlaneConstantsCB.PlaneLengthWorldSpace / float(numLights.x);
				const float2 curLightPosWorld = distBetweenLights * lightPosOffsetCell + distBetweenLights * curLightIndex;

				const float curLightIntesity = pow(FireTextureDownsampledSRV[curLightIndex], 2) * CGPlaneConstantsCB.VirtualLightIntensity * 0.1;

				const float distance = length(curLightPosWorld - posWorld);
				const float attenuation = saturate(1.f - pow(distance / CGPlaneConstantsCB.VirtualLightRadius, 2));

				float2 noiseuv = (float2(curLightIndex) + 0.5f) / float2(numLights);
				noiseuv.x += CGPlaneConstantsCB.Time * CGPlaneConstantsCB.VirtualLightFlickerSpeed * 0.1;
				//noiseuv.y += float(curLightIndex.y * curLightIndex.x) *0.01f;
				float noiseT = NoiseTextureSRV.SampleLevel(LinearWrapSampler, noiseuv, 0).r;

				noiseT = MapToRange(noiseT, 0.2, 0.8, 0, 1);

				float3 lightColor1 = float3(1,0.8,1) * 0.1;
				float3 lightColor2 = float3(1, 0.5, 0);

				const float3 lightColorFinal = lerp(lightColor1, lightColor2, noiseT); 

				colorFinal += (lightColorFinal * curLightIntesity * attenuation);

			}
		}

		paperColor += colorFinal;

	}

	


	curFire *= CGPlaneConstantsCB.VisualizerFireScale;
	curFire = pow(curFire, CGPlaneConstantsCB.VisualizerFirePow);

	float3 fireColor;
	if(curFire >= 1.f)
	{
		fireColor = float3(curFire, curFire * 0.2, curFire * 0.1);
	}
	else
	{
		fireColor = float3(curFire * 0.1f, curFire * 0.2, curFire);
	}


	float3 finalColor = max(paperColor, fireColor);
	
	#if 0//amplify fire transition edge
	{
		float h = 0.05 * 0.1;
		float n = FireTextureSRV.SampleLevel(LinearSampler, uv + float2(0, h), 0);
		float e = FireTextureSRV.SampleLevel(LinearSampler, uv + float2(h, 0), 0);
		float s = FireTextureSRV.SampleLevel(LinearSampler, uv + float2(0, -h), 0);
		float w = FireTextureSRV.SampleLevel(LinearSampler, uv + float2(-h, 0), 0);

    	float dy = (n - s)*.5;
    	float dx = (e - w)*.5;

    	float edge = sqrt(dx*dx + dy*dy);

		finalColor.rgb += edge * 0.1;
	}
	#endif


	VisualizerTextureUAV[SampleCoord] = finalColor;

}




//Render Plane

struct PlaneRenderInterpolators
{
	float4 Position : SV_Position;
	float2 TexCoords : TEXCOORD;
};
void PlaneRenderVS(
	in uint VertexIndex : SV_VertexID,
	out PlaneRenderInterpolators vsout
	)
{
	/*
	*	Position
	*/
	float3 positionWorld = PlaneVertexBufferSRV[VertexIndex];	
	vsout.Position = mul(CGPlaneConstantsCB.SceneViewProjectionMatrix, float4(positionWorld, 1.0f));
	/*
	*	TexCoord
	*/
	float2 uvLocal = PlaneTexCoordsBufferSRV[VertexIndex];
	vsout.TexCoords = uvLocal;
}
float4 PlaneRenderPS
	(
	in PlaneRenderInterpolators psin
	) : SV_Target
{

	return float4(NoiseTextureSRV.SampleLevel(LinearSampler, psin.TexCoords, 0).rgb, 1);

}

float BurnedImagePlaneRenderPS
	(
	in PlaneRenderInterpolators psin
	) : SV_Target
{
	//render image edges as embers

	if(FuelTextureSRV.SampleLevel(PointSampler, psin.TexCoords, 0) > 0.01)
	{
		return 0;
	}
	else
	{
		float h = CGPlaneConstantsCB.BurnedImageSharpness * 0.1;
		float3 n = NoiseTextureSRV.SampleLevel(LinearSampler, psin.TexCoords + float2(0, h), 0);
		float3 e = NoiseTextureSRV.SampleLevel(LinearSampler, psin.TexCoords + float2(h, 0), 0);
		float3 s = NoiseTextureSRV.SampleLevel(LinearSampler, psin.TexCoords + float2(0, -h), 0);
		float3 w = NoiseTextureSRV.SampleLevel(LinearSampler, psin.TexCoords + float2(-h, 0), 0);

    	float3 dy = (n - s)*.5;
    	float3 dx = (e - w)*.5;

    	float3 edge = sqrt(dx*dx + dy*dy);
		float luminance = dot( edge.rgb, float3( 0.299f, 0.587f, 0.114f ) );
    	return luminance * 5.0 * float3(1, 0.2, 0.1);
	}
}













//==================================================================BURNER PLANE


/*
*	Resources
*/

//CB

ConstantBuffer<CGCombinerConstants> CGCombinerConstantsCB : register(b0);

//SRV
Texture2D<float4> FlameTextureSRV : register (t0);
Texture2D<float4> SmokeTextureSRV : register (t1);
Texture2D<float4> PlaneTextureSRV : register (t2);
Texture2D<float4> BloomTextureSRV : register (t3);
Texture2D<float> CharcoalTextureSRV : register (t4);
Texture2D<float3> CharcoalNoiseTextureSRV : register (t5);
Texture2D<float3> BurnedImageTextureSRV : register (t6);
Texture2D<float> FlameNoiseTextureSRV : register (t7);
Texture2D<float> FlameNoise2TextureSRV : register (t8);
Texture2D<float4> LensDirtTextureSRV : register (t9);

//UAV
RWTexture2D<float4> CombineResultTextureUAV : register (u0);

[numthreads( 8, 8, 1 )]
void CombinerPassMain( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float2 sampleCoordf = float2(SampleCoord) + 0.5;
	const float2 uv = sampleCoordf / float2(GTextureSize);

	const float4 flame = FlameTextureSRV[SampleCoord];

	#if 1//Heat haze:
	float2 distortionUV = uv;
	distortionUV.y += CGCombinerConstantsCB.Time * 0.1;
	distortionUV.x *= 2.5;
	distortionUV.y *= 0.5;
	float3 distortionNoise = CharcoalNoiseTextureSRV.SampleLevel(LinearWrapSampler, distortionUV, 0);
	distortionNoise = (distortionNoise * 2.f) - 1.f;
	distortionUV = uv;
	float t = (CGCombinerConstantsCB.Time * 1.f) % 2.f;
	if(t < 1.f)
	{
		distortionNoise.r = lerp(distortionNoise.r, distortionNoise.g, t);
	}
	else
	{
		distortionNoise.r = lerp(distortionNoise.g, distortionNoise.r, t - 1.f);
	}
	float4 heat = FlameTextureSRV.SampleLevel(LinearWrapSampler, uv + float2(0, 0.08), 0);
	distortionNoise.x *= 0.0025;
	distortionNoise.y *= 0.005;
	distortionNoise *= 0.5f;
	distortionNoise *= (heat.a + heat.r);
	distortionUV.x += distortionNoise.r;
	distortionUV.y += distortionNoise.g;
	float4 plane = PlaneTextureSRV.SampleLevel(LinearSampler, distortionUV, 0);
	#else
	//float4 plane = PlaneTextureSRV[SampleCoord];
	#endif

	float4 smoke = SmokeTextureSRV[SampleCoord];
	float4 bloom = BloomTextureSRV[SampleCoord];

	float4 lensDirt = LensDirtTextureSRV.SampleLevel(LinearSampler, uv, 0);
	
	//bloom.rgb += ((bloom.rgb * lensDirt.rgb)*10.f); 

	//Add bloom to Plane
	plane.rgb += bloom * CGCombinerConstantsCB.BloomStrength * CGCombinerConstantsCB.PlaneBloomStrength;

#if SMOKE_BLOOM_SEPARATE
	smoke.rgb += bloom * CGCombinerConstantsCB.BloomStrength * saturate(1.f + clamp(smoke.r, CGCombinerConstantsCB.SmokeBloomColorClampMinMax.x, CGCombinerConstantsCB.SmokeBloomColorClampMinMax.y) - clamp(smoke.a, CGCombinerConstantsCB.SmokeBloomAlphaClampMinMax.x, CGCombinerConstantsCB.SmokeBloomAlphaClampMinMax.y)) * clamp(smoke.a, 0., 1);
#endif

	//Blend smoke onto Plane
	float4 final = float4(smoke.rgb * smoke.a + plane.rgb * clamp(1.f - smoke.a, 0.0, 1.f), 1.f);

#if SMOKE_BLOOM_SEPARATE == 0
	//Apply Bloom to Smoke and Plane
	final.rgb += bloom * CGCombinerConstantsCB.BloomStrength * saturate(1.f + clamp(smoke.r, CGCombinerConstantsCB.SmokeBloomColorClampMinMax.x, CGCombinerConstantsCB.SmokeBloomColorClampMinMax.y) - clamp(smoke.a, CGCombinerConstantsCB.SmokeBloomAlphaClampMinMax.x, CGCombinerConstantsCB.SmokeBloomAlphaClampMinMax.y)) * clamp(smoke.a, 0., 1);
#endif

	final.rgb += ((saturate(bloom.rgb) * lensDirt.rgb) * 0.75f);

	//Blend Flame
	//final = float4(flame.rgb * flame.a + final.rgb * (1.f - flame.a), 1.f);
	final.rgb = max(flame.rgb, final.rgb);

	

	//Tonemap Final
	final.rgb *= CGCombinerConstantsCB.ExposureFinal;
	
	CombineResultTextureUAV[SampleCoord] = final;


}

float GetCharcoalNoise(float2 uv)
{
	float3 noise = CharcoalNoiseTextureSRV.SampleLevel(LinearWrapSampler, uv, 0);

	float finalNoise;

	if(CGCombinerConstantsCB.NoiseTextureInterpolator < 1.f)
	{
		finalNoise = lerp(noise.x, noise.y, CGCombinerConstantsCB.NoiseTextureInterpolator);
	}
	else if(CGCombinerConstantsCB.NoiseTextureInterpolator < 2.f)
	{
		finalNoise = lerp(noise.y, noise.z, frac(CGCombinerConstantsCB.NoiseTextureInterpolator));
	}
	else
	{
		finalNoise = lerp(noise.z, noise.x, frac(CGCombinerConstantsCB.NoiseTextureInterpolator));
	}

	return finalNoise;
}

[numthreads( 8, 8, 1 )]
void CharcoalCombine( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float charcoalFlame = CharcoalTextureSRV[SampleCoord] * CGCombinerConstantsCB.CharcoalIntensity;
	float burnedImageColor = BurnedImageTextureSRV[SampleCoord] * CGCombinerConstantsCB.BurnedImageIntensity;

	float2 sampleCoordf = float2(SampleCoord) + 0.5;
	float2 uv = sampleCoordf / float2(GTextureSize);
	uv.y += CGCombinerConstantsCB.Time * 0.01f * 0.05f;
	uv.x += CGCombinerConstantsCB.Time * 0.0011f * 0.05f;
	float noise = GetCharcoalNoise(uv);
	noise = clamp(MapToRange(noise, 0.4, 0.6, 0, 1), 0.1, 1);

	charcoalFlame = max(charcoalFlame, burnedImageColor);

	charcoalFlame *= noise;

	float3 charcoalColor = (float3(1.f, 0.2, 0.1) * charcoalFlame);

	/* CombineResultTextureUAV[SampleCoord] = float4(charcoalColor.rgb, 1.f);
	return; */
	
	float4 plane = PlaneTextureSRV[SampleCoord];

	plane.rgb = max(plane.rgb, charcoalColor);

	CombineResultTextureUAV[SampleCoord] = plane;

}


[numthreads( 8, 8, 1 )]
void FlamePostProcess( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	
	float2 sampleCoordf = float2(SampleCoord) + 0.5;
	float2 flameNoiseUV = sampleCoordf / float2(GTextureSize);

	float2 flameSamplingUV = flameNoiseUV;
	flameSamplingUV.y -= CGCombinerConstantsCB.Time * 0.1;
	flameSamplingUV *= 0.1;
	float3 distortionNoise = CharcoalNoiseTextureSRV.SampleLevel(LinearWrapSampler, flameSamplingUV, 0);
	distortionNoise = (distortionNoise * 2.f) - 1.f;
	flameSamplingUV = flameNoiseUV;

	float t = (CGCombinerConstantsCB.Time * 1.f) % 2.f;
	if(t < 1.f)
	{
		distortionNoise.r = lerp(distortionNoise.r, distortionNoise.g, t);
	}
	else
	{
		distortionNoise.r = lerp(distortionNoise.g, distortionNoise.r, t - 1.f);
	}
	flameSamplingUV.x += distortionNoise.r * 0.0075;
	flameSamplingUV.y += distortionNoise.g * 0.001;

	flameNoiseUV = flameSamplingUV;

	//float4 flame = FlameTextureSRV[SampleCoord];
	float4 flame = FlameTextureSRV.SampleLevel(LinearSampler, flameSamplingUV, 0);

	float2 perlinNoiseUV = sampleCoordf / float2(GTextureSize);
	perlinNoiseUV.y -= CGCombinerConstantsCB.Time * 0.01;
	perlinNoiseUV.x += CGCombinerConstantsCB.Time * 0.05;
	perlinNoiseUV *= 0.1f;
	float perlinNoise = GetCharcoalNoise(perlinNoiseUV);
	perlinNoise = clamp(MapToRange(perlinNoise, 0.4, 0.6, 0, 1), 0.f, 1.f);
	

	//Pre-Translate Scale
	flameNoiseUV *= 0.9f;
	flameNoiseUV.x *= 5.f;

	//Translate
	const float flameSpeed = 0.4f; //TODO: USE VARYING SPEED [0.25,0.75]
	flameNoiseUV.y += CGCombinerConstantsCB.Time * flameSpeed;
	flameNoiseUV.x += CGCombinerConstantsCB.Time * 0.05f;

	//Post-Translate Scale
	flameNoiseUV *= 0.6f;

	float2 flameNoise2;
	flameNoise2.r = FlameNoiseTextureSRV.SampleLevel(LinearWrapSampler, flameNoiseUV, 0);
	flameNoise2.g = FlameNoise2TextureSRV.SampleLevel(LinearWrapSampler, flameNoiseUV, 0);

	
	float flameNoise;
	if(t < 1.f)
	{
		flameNoise = lerp(flameNoise2.r, flameNoise2.g, t);
	}
	else
	{
		flameNoise = lerp(flameNoise2.g, flameNoise2.r, t - 1.f);
	}

	//perlinNoise = min(0.5f, perlinNoise);
	//flameNoise = saturate(flameNoise - perlinNoise);
	flameNoise = 1.f - flameNoise;
	//flameNoise = saturate(flameNoise + perlinNoise);

	flame.rgb *= flameNoise;

	CombineResultTextureUAV[SampleCoord] = flame;


}

[numthreads( 8, 8, 1 )]
void BloomCombine( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float4 flame = FlameTextureSRV[SampleCoord];
	float4 plane = PlaneTextureSRV[SampleCoord];

	flame.rgb += 0.25;

	//tonemap fire plane
	plane.rgb = /* saturate */(log(plane.rgb + 1.f) /* * 0.75f */);
	// reinhard tone mapping
    //plane.rgb = plane.rgb / (plane.rgb + (1.0f));
    // gamma correction 
	//const float gamma = 2.2;
    //plane.rgb = pow(plane.rgb, (1.0 / gamma));

	float4 final = float4(0,0,0,1);

	//float3 fireColor = max(flame.rgb, plane.rgb);
	//float3 fireColor = flame.rgb + plane.rgb;
	float3 fireColor = flame.rgb + plane.rgb * (1.f - flame.a);

	// 1. Do a very simple mathematical average:
	//float luminance = dot( fireColor.rgb, float3( 0.33f, 0.33f, 0.33f ) );

	// 2. Perform a more accurately weighted average:
	//float luminance = dot( average.rgb, float3( 0.299f, 0.587f, 0.114f ) );

	// 3. Take the maximum value of the incoming, same as computing the
	// brightness/value for an HSV/HSB conversion:
	//float luminance = max( fireColor.r, max( fireColor.g, fireColor.b ) );

	// 4. Compute the luminance component as per the HSL colour space:
	//float luminance = 0.5f * ( max( average.r, max( average.g, average.b ) ) + min( average.r, min( average.g, average.b ) ) );

	// 5. Use the magnitude of the colour
	//float luminance = length( average.rgb );


	if(any(fireColor.rgb > CGCombinerConstantsCB.BloomThreshold))
	{
		final.rgb = fireColor;
		if(all(final.rgb < 1))
		{
			final.rgb = pow(final.rgb, CGCombinerConstantsCB.BloomCurvePow);
		}
	}

	CombineResultTextureUAV[SampleCoord] = final;

}

#define BLUR_DOWNSAMPLE 1

[numthreads( 8, 8, 1 )]
void Downsample( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float2 sampleCoordf = float2(SampleCoord) + 0.5;
	float2 uv = sampleCoordf / float2(GTextureSize);

#if BLUR_DOWNSAMPLE
	float2 srcTexelSize = 1.f / float2(GTextureSize << 1u);
	float4 leftup = FlameTextureSRV.SampleLevel(LinearSampler, float2(uv.x - srcTexelSize.x, uv.y - srcTexelSize.y), 0);
	float4 rightup = FlameTextureSRV.SampleLevel(LinearSampler, float2(uv.x + srcTexelSize.x, uv.y - srcTexelSize.y), 0);
	float4 leftdown = FlameTextureSRV.SampleLevel(LinearSampler, float2(uv.x - srcTexelSize.x, uv.y + srcTexelSize.y), 0);
	float4 rightdown = FlameTextureSRV.SampleLevel(LinearSampler, float2(uv.x + srcTexelSize.x, uv.y + srcTexelSize.y), 0);
	float4 source = (leftup + rightup + leftdown + rightdown) * 0.25f;
#else
	float4 source = FlameTextureSRV.SampleLevel(LinearSampler, uv, 0);
#endif
	 CombineResultTextureUAV[SampleCoord] = source;
}


#define BLUR_UPSAMPLE 1

[numthreads( 8, 8, 1 )]
void Upsample( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	float2 sampleCoordf = float2(SampleCoord) + 0.5;

	float2 uv = sampleCoordf / float2(GTextureSize);

	float4 source = FlameTextureSRV.SampleLevel(LinearSampler, uv, 0);

	CombineResultTextureUAV[SampleCoord] = source;
}
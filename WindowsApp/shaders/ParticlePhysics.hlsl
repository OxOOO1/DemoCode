
#include "Physics\ParticlePhysicsCommonStructs.h"

#include "Physics\PhysicsCommon.hlsli"

#include "InterpolationEtc\InterpolationCommon.hlsli"

//#include "mixbox.hlsl"

//===========================================================================
//
//                              PARTICLE PHYSICS
//
//===========================================================================

#define EPSILON 0.000001
#define UINT_MAX      0xffffffff
#define PI 3.141592654f

/*
*	Resources
*/

SamplerState LinearSampler : register(s0);
SamplerState PointSampler : register(s1);
SamplerState LinearSamplerBorder : register(s2);
//SamplerState LinearSamplerRepeat : register(s3);

ConstantBuffer<CGConstants> CGConstantsCB : register(b0);

ConstantBuffer<GCollisionConstants> GCollisionConstantsCB : register(b1);

ConstantBuffer<SpatialControlConstants> CGSpatialControlConstantsCB : register(b2);

ConstantBuffer<GMeshGenerationConstants> GMeshGenerationConstantsCB : register(b3);

ConstantBuffer<VectorProjectionPerPassConstants> VectorProjectionPerPassConstantsCB : register(b4);


//SRV
StructuredBuffer<float4> PlanesBufferSRV : register(t0);

StructuredBuffer<float3> ParticlesBufferPositionCurSRV : register(t1);

StructuredBuffer<uint2> CellIdOffsetAndCountBufferSRV : register(t2);

StructuredBuffer<uint> ObjectListIndexBufferSRV : register(t3);
StructuredBuffer<uint> CellListIndexBufferSRV : register(t4);
StructuredBuffer<uint> PerObjectControlBitsBufferSRV : register(t5);

StructuredBuffer<CConstraint> ConstraintsBufferSRV : register(t6);

StructuredBuffer<uint> ParticlesStateBufferSRV : register(t7);

StructuredBuffer<uint> ConstraintsStatesBufferSRV : register(t8);

StructuredBuffer<uint2> ParticlesToRenderCoordListSRV : register(t9);

StructuredBuffer<uint> ConstraintsStatesBufferBLUESRV : register(t10);
StructuredBuffer<uint> ConstraintsStatesBufferGREENSRV : register(t11);
StructuredBuffer<uint> ConstraintsStatesBufferYELLOWSRV : register(t12);
StructuredBuffer<uint> ConstraintsStatesBufferREDSRV : register(t13);

StructuredBuffer<float4> SphereCapsulesBufferSRV : register(t14);

StructuredBuffer<float3> ParticlesBufferPositionAfterCreationSRV : register(t15);

Texture2D<float> VelocityField_U_ComponentTextureSRV : register(t16);
Texture2D<float> VelocityField_V_ComponentTextureSRV : register(t17);

StructuredBuffer<float3> CurParticlesVelocityBufferSRV : register(t18);

Texture2D<uint> VelocityFieldStateTextureSRV : register(t19);

Texture2D<DENSITY_TYPE> DensityFieldTextureSRV : register(t20);

Texture2D<float> PressureFieldTextureSRV : register(t21);

Texture2D<float> VectorFieldDivergenceTextureSRV : register(t22);

Texture2D<float> VelocityField_U_ComponentTextureNewSRV : register(t23);
Texture2D<float> VelocityField_V_ComponentTextureNewSRV : register(t24);

Texture2D<float> VectorFieldCurlTextureSRV : register(t25);

Texture2D<DENSITY_TYPE> DensityFieldTextureNewSRV : register(t26);

Texture2D<float3> DensitySourceColorTextureSRV : register(t27);

Texture2D<float> InputTextureWithNegativeValuesSRV : register(t28);//for DEBUG Visualizer

Texture2D<float> CurlGradientTextureSRV : register(t29);

Texture2D<float2> CurlGradientMinMaxTextureSRV : register(t30);

Texture2D<DENSITY_TYPE> DensityFieldMacCormackAdvectedTextureSRV : register(t31);
Texture2D<DENSITY_TYPE> DensityFieldMacCormackReversedTextureSRV : register(t32);

Texture2D<float> VelocityField_U_ComponentWeightTextureSRV : register(t33);
Texture2D<float> VelocityField_V_ComponentWeightTextureSRV : register(t34);

Texture2D<float> AdvectedParticleDensityTextureSRV : register(t20);//TODO:Change Register ----!!!!

//UAV
RWStructuredBuffer<float3> ParticlesBufferPositionCurUAV : register (u0);
RWStructuredBuffer<float3> ParticlesBufferPositionPrevUAV : register (u1);

RWStructuredBuffer<uint> NumCellsObjectIntersectsUAV : register(u2);

//Assign each object a set of 11 control bits,
//3 bits hold the type of the cell where the centroid is,
//8 bits specify the types of the cells intersected
RWStructuredBuffer<uint> PerObjectControlBitsBufferUAV : register(u3);

globallycoherent RWStructuredBuffer<uint> AtomicOffsetBufferUAV : register(u4);

RWStructuredBuffer<uint> ObjectListIndexBufferUAV : register(u5);
RWStructuredBuffer<uint> CellListIndexBufferUAV : register(u6);

RWStructuredBuffer<uint> ParticlesStateBufferUAV : register(u7);

RWStructuredBuffer<uint> IntersectedPointIndexUAV : register(u8);
RWStructuredBuffer<float3> ClickedObjectViewSpaceZUAV : register(u9);

RWStructuredBuffer<uint> ConstraintsStatesBufferUAV : register(u10);

RWStructuredBuffer<float3> MeshBufferVerticesUAV : register(u11);
RWStructuredBuffer<float3> MeshBufferNormalsUAV : register(u12);

RWStructuredBuffer<uint2> ParticlesToRenderCoordListUAV : register(u13);

RWStructuredBuffer<uint> InterlockedAddBuffer : register(u14);

RWStructuredBuffer<float3> CurParticlesVelocityBufferUAV : register(u15);

RWTexture2D<float> VelocityField_U_ComponentTextureUAV : register(u16);
RWTexture2D<float> VelocityField_V_ComponentTextureUAV : register(u17);

RWTexture2D<float> VectorFieldDivergenceTextureUAV : register(u18);

RWTexture2D<DENSITY_TYPE> DensityFieldTextureUAV : register(u19);

RWTexture2D<float> PressureFieldTextureUAV : register(u20);

RWTexture2D<float> VectorFieldCurlTextureUAV : register(u21);

RWStructuredBuffer<float2> MeshBufferTexCoordsUAV : register(u22);

RWTexture2D<uint> VelocityFieldStateTextureUAV : register(u23);

RWTexture2D<float2> VisualizerColorTextureUAV : register(u24);//for DEBUG Visualizer

RWTexture2D<float> CurlGradientTextureUAV : register(u25);
RWTexture2D<float3> CurlGradientDirectionTextureUAV : register(u26);

RWTexture2D<float2> PressureGradientTextureUAV : register(u27);

RWTexture2D<float> VelocityField_U_ComponentWeightTextureUAV : register(u28);
RWTexture2D<float> VelocityField_V_ComponentWeightTextureUAV : register(u29);

RWTexture2D<float> AdvectedParticleDensityTextureUAV : register(u30);

RWStructuredBuffer<float> AverageDensityBufferUAV : register(u31);

RWStructuredBuffer<float3> ParticlesColorsBufferUAV : register(u20);//TODO:Change Register ----!!!!

RWStructuredBuffer<float3> CurParticlesExternalForceBufferUAV : register(u9);//TODO:Change Register ----!!!!

RWStructuredBuffer<float3> PrevParticlesVelocityBufferUAV : register(u22);

#define FORCE_POS 0 ///DEBUG
#define RESOLVE_NAN_FLOATS 0
#define RESOLVE_ZERO_DISTANCE 1 ///DISABLING IT MIGHT CRASH THE SIMULATION

uint GetFlatIndexToVelocityArray2D(uint2 coord)
{
	//invert v
	coord.y = (CGConstantsCB.VelocityFieldGridSize - 1) - coord.y;
	return coord.y * CGConstantsCB.VelocityFieldGridSize + coord.x;
};

float GetVelocityFieldBoundsStart()
{
	return CGConstantsCB.VelocityFieldBoundsStartOffset;
}
float GetVelocityFieldBoundsEnd()
{
	return CGConstantsCB.VelocityFieldBoundsStartOffset + CGConstantsCB.VelocityFieldBoundsLength;
}

float GetSimulationBoundsStart()
{
	return CGConstantsCB.SimulationBoundsStartOffset;
}
float GetSimulationBoundsEnd()
{
	return CGConstantsCB.SimulationBoundsStartOffset + CGConstantsCB.SimulationBoundsLength;
}

[numthreads( 64, 1, 1 )]
void SimulateParticles( 
	uint SampleCoord : SV_DispatchThreadID )
{

	if(SampleCoord < CGConstantsCB.NumParticles)
	{
		uint particleStateBitmask = ParticlesStateBufferSRV[SampleCoord];
		if((particleStateBitmask & BITMASK_PARTICLE_STATE_PINNED) != 0)
		{
			return;
		}

		float3 posCur = ParticlesBufferPositionCurUAV[SampleCoord];
		float3 posPrev = ParticlesBufferPositionPrevUAV[SampleCoord];

#if RESOLVE_NAN_FLOATS
		if(any(isnan(posCur)))
		{
			posCur = posPrev = CGConstantsCB.DefaultParticlePos;
		}
#endif//RESOLVE_NAN_FLOATS

		//Apply Forces
		float3 netForce = float3(0,0,0);

		//External
		netForce += CGConstantsCB.ExternalForce;

		//Per particle
		if((particleStateBitmask & BITMASK_PARTICLE_STATE_FALLING) == 0)
		{
			netForce += CurParticlesExternalForceBufferUAV[SampleCoord];
			CurParticlesExternalForceBufferUAV[SampleCoord] = float3(0,0,0);
		}
		
		float3 curVelocity;

		//Integrate
	#if INTEGRATE_VERLET
		
		IntegrateVerlet(posCur, posPrev, curVelocity, netForce, CGConstantsCB.Damping, CGConstantsCB.DeltaTime);
	#else//EULER
		/* if((particleStateBitmask & BITMASK_PARTICLE_STATE_FALLING) != 0)
		{
			IntegrateVerlet(posCur, posPrev, curVelocity, netForce, CGConstantsCB.Damping, CGConstantsCB.DeltaTime);
		}
		else */
		{
			curVelocity = CurParticlesVelocityBufferUAV[SampleCoord];
			IntegrateEuler(posCur, posPrev, curVelocity, netForce, CGConstantsCB.Damping, CGConstantsCB.DeltaTime);
			//IntegrateRK3(posCur, posPrev, curVelocity, netForce, CGConstantsCB.Damping, CGConstantsCB.DeltaTime);
		}
		
	#endif//INTEGRATE_VERLET
		
		if(CGConstantsCB.bSimulateParticlesIn2D)
		{
			if(CGConstantsCB.bSimulateParticlesIn2DAllowNegativeZOffset)
			{
				posCur.z = min(CGConstantsCB.DefaultParticlePos.z, posCur.z);
				posPrev.z = min(CGConstantsCB.DefaultParticlePos.z, posCur.z);
			}
			else
			{
				posCur.z = CGConstantsCB.DefaultParticlePos.z;
				posPrev.z = CGConstantsCB.DefaultParticlePos.z;
			}
			
		}

		//Cache prev-cur Velocity
		PrevParticlesVelocityBufferUAV[SampleCoord] = CurParticlesVelocityBufferUAV[SampleCoord];
		CurParticlesVelocityBufferUAV[SampleCoord] = curVelocity;

		//Write back
		ParticlesBufferPositionCurUAV[SampleCoord] = posCur;
		ParticlesBufferPositionPrevUAV[SampleCoord] = posPrev;

	}
	

}



//===========================================================================
//                              COLLISION
//===========================================================================

[numthreads( 64, 1, 1 )]
void ResolveExternalCollisions( 
	uint SampleCoord : SV_DispatchThreadID )
{

	if(SampleCoord < CGConstantsCB.NumParticles)
	{
		uint particleStateBitmask = ParticlesStateBufferSRV[SampleCoord];
		if((particleStateBitmask & BITMASK_PARTICLE_STATE_PINNED) != 0)
		{
			return;
		}

		float3 posCur = ParticlesBufferPositionCurUAV[SampleCoord];

		//Resolve Collisions etc

		const bool bDisableFloorPlane = true;

		//Planes
		const uint NumPlanes = 6;
		for(int i = 0; i < NumPlanes; i++)
		{
			const float4 planeVal = PlanesBufferSRV[i];
			BoundingPlane plane;
			plane.n = planeVal.xyz;
			plane.d = planeVal.w;
			
			if(bDisableFloorPlane)
			{
				if(abs(plane.n.y - 1.f) < 0.1)
				{
					continue;
				}
			}
			if(ResolvePointPlaneCollision(posCur, CGConstantsCB.PointRadius, plane, CGConstantsCB.CollisionResponseCoef))
			{
				float3 curV = CurParticlesVelocityBufferUAV[SampleCoord];
				curV = curV + plane.n * abs(dot(plane.n, curV));
				if(CGConstantsCB.bSimulateParticlesIn2D)
				{
					curV.z = 0.f;
				}
				CurParticlesVelocityBufferUAV[SampleCoord] = curV;
			}
		}


#if 1
		//Capsules

		//Load capsule
		const float4 capsulePos1AndRadius = SphereCapsulesBufferSRV[0];
		const float3 capsulePos2 = SphereCapsulesBufferSRV[1].xyz;
		BoundingCapsule capsule;
		capsule.Pos1 = capsulePos1AndRadius.xyz;
		capsule.Pos2 = capsulePos2;
		capsule.Radius = capsulePos1AndRadius.w;

		//Resolve
		ResolvePointCapsuleCollision(posCur, CGConstantsCB.PointRadius, capsule, CGConstantsCB.CollisionResponseCoef);
#endif

		if(CGConstantsCB.bSimulateParticlesIn2D)
		{
			if(CGConstantsCB.bSimulateParticlesIn2DAllowNegativeZOffset)
			{
				posCur.z = min(CGConstantsCB.DefaultParticlePos.z, posCur.z);
			}
			else
			{
				posCur.z = CGConstantsCB.DefaultParticlePos.z;
			}
		}

		//Write back
		ParticlesBufferPositionCurUAV[SampleCoord] = posCur;

	}
	

}





//===========================================================================
//                              SELF COLLISION
//===========================================================================

float GetCollisionGridCellSize()
{
	return CGConstantsCB.PointRadius * CellSizeScale;
}

[numthreads( 64, 1, 1 )]
void GenerateIntersectionPairList( 
	uint SampleCoord : SV_DispatchThreadID )
{
	if(SampleCoord < CGConstantsCB.NumParticles)
	{

		uint curParticleState = ParticlesStateBufferSRV[SampleCoord];

		if((curParticleState & BITMASK_PARTICLE_STATE_DISABLED) != 0)
		{
			//don't add this particle into collision simulation
			return;
		}

		uint ControlBits = 0;

		float3 posCur = ParticlesBufferPositionCurSRV[SampleCoord];

		const float cellSize = GetCollisionGridCellSize();

		///Home Cell
		int3 homeCell_Coord = floor((posCur) / cellSize);
		int homeCell_Id = GetCellIndex(homeCell_Coord);
		uint homeCellLocalId = GetLocalCellIndex(homeCell_Coord);
		//Write down home cell Id to control bits
		ControlBits |= (homeCellLocalId << 8);


		//Intersection Cells
		int3 cellIdMin = floor((posCur - CGConstantsCB.PointRadius) / cellSize);
		cellIdMin = max(int3(0,0,0), cellIdMin);
		int3 cellIdMax = floor((posCur + CGConstantsCB.PointRadius) / cellSize);
		cellIdMax = min(int3(CELL_MAX_INDEX,CELL_MAX_INDEX,CELL_MAX_INDEX), cellIdMax);
		int3 curCellId;

		//Compute Num Cells we output
		uint cellCounter = 0;

		for (curCellId.x = cellIdMin.x; curCellId.x <= cellIdMax.x; curCellId.x++)
		{
			for (curCellId.y = cellIdMin.y; curCellId.y <= cellIdMax.y; curCellId.y++)
			{
				for (curCellId.z = cellIdMin.z; curCellId.z <= cellIdMax.z; curCellId.z++)
				{
					cellCounter++;
				}
			}
		}


		uint writeOffset;
		InterlockedAdd(AtomicOffsetBufferUAV[0], cellCounter, writeOffset);

		int lastCellId = -1;

		///Output Cells
		cellCounter = 0;
		for (curCellId.x = cellIdMin.x; curCellId.x <= cellIdMax.x; curCellId.x++)
		{
			for (curCellId.y = cellIdMin.y; curCellId.y <= cellIdMax.y; curCellId.y++)
			{
				for (curCellId.z = cellIdMin.z; curCellId.z <= cellIdMax.z; curCellId.z++)
				{
					uint cell_id = GetCellIndex(curCellId);
					uint localCellId = GetLocalCellIndex(curCellId);
					ControlBits |= ( 1u << (localCellId));
					if(cell_id != lastCellId)
					{
						ObjectListIndexBufferUAV[writeOffset + cellCounter] = SampleCoord;
						CellListIndexBufferUAV[writeOffset + cellCounter] = cell_id;
					}

					lastCellId = cell_id;

					cellCounter++;
				}
			}
		}

		PerObjectControlBitsBufferUAV[SampleCoord] = ControlBits;

	}
}


/*
*	Compute Offsets
*/
#if 0
[numthreads( 64, 1, 1 )]
void ComputeOffsetsForCells( 
	uint SampleCoord : SV_DispatchThreadID )
{
	///Iterate Over sorted Cells list
	///Write down offset + num elements for each work item



}
#endif


#define COLLISION_SKIP_ENABLED 0
#define HOME_CELL_COLLISION_SKIP_ENABLED 0
#define UPDATE_POSITION_INSTANTLY 1//Gauss-Siedel

#define DIFFUSE_PARTICLES_COLORS_ON_COLLISION 1

[numthreads( 64, 1, 1 )]
void ResolveCollisions( 
	uint SampleCoord : SV_DispatchThreadID )
{

	if(SampleCoord < GCollisionConstantsCB.NumCells)
	{

		const uint2 offsetAndCount = CellIdOffsetAndCountBufferSRV[SampleCoord];
		const uint offset = offsetAndCount.x;
		const uint count = offsetAndCount.y;

		const uint GCurCellLocalId = GCollisionConstantsCB.CurCellLocalId;

		for(int i = 0; i < count-1; i++)
		{
			bool bCurObjectHadCollisions = false;
			const uint curObjectIndex = ObjectListIndexBufferSRV[offset + i];
			const uint curObjectControlBits = PerObjectControlBitsBufferSRV[curObjectIndex];
			float3 curObjectPos = ParticlesBufferPositionCurUAV[curObjectIndex];

			float3 curObjectColor = ParticlesColorsBufferUAV[curObjectIndex];

			float3 prevObjectPos = ParticlesBufferPositionPrevUAV[curObjectIndex];
			const uint curObjectHomeCellId = (curObjectControlBits >> 8);
			const bool bCurObjectInHomeCell = curObjectHomeCellId == GCurCellLocalId;
			const bool bCurObjectWasProcessedPreviously = curObjectHomeCellId < GCurCellLocalId;

			for(int j = i+1; j < count; j++)
			{
				const uint otherObjectIndex = ObjectListIndexBufferSRV[offset + j];
				uint otherObjectControlBits = PerObjectControlBitsBufferSRV[otherObjectIndex];

				const uint otherObjectHomeCellId = (otherObjectControlBits >> 8);
				const bool bOtherObjectInHomeCell = otherObjectHomeCellId == GCurCellLocalId;
				const bool bOtherObjectWasProcessedPreviously = curObjectHomeCellId < GCurCellLocalId;
#if HOME_CELL_COLLISION_SKIP_ENABLED
				if(bCurObjectInHomeCell || bOtherObjectInHomeCell)
#endif//HOME_CELL_COLLISION_SKIP_ENABLED
				{
					//Avoid collision resolve if at least one of the objects has home cell local id that is < than cur cell id
					//and other object has this cell id in ControlBits
					//(meaning this collision was resolved already in prev pass)
#if COLLISION_SKIP_ENABLED
					if(bCurObjectWasProcessedPreviously && (1u << curObjectHomeCellId) & (otherObjectControlBits & 0xF))
					{
						continue;
					}
					else if(bOtherObjectWasProcessedPreviously && (1u << otherObjectHomeCellId) & (curObjectControlBits & 0xF))
					{
						continue;
					}
					else
#endif//COLLISION_SKIP_ENABLED
					{
						///perform collision
						const float3 otherObjectPos = ParticlesBufferPositionCurUAV[otherObjectIndex];
						float3 collisionAxis = curObjectPos - otherObjectPos;
						float distanceSq =  dot(collisionAxis, collisionAxis);
						float radiusSum = CGConstantsCB.PointRadius * 2.f;
						if(distanceSq < radiusSum * radiusSum)
						{
#if RESOLVE_ZERO_DISTANCE
							if(distanceSq < EPSILON)
							{
								distanceSq = EPSILON;
								if(any(prevObjectPos != curObjectPos))
								{
									collisionAxis = normalize(prevObjectPos - curObjectPos);
								}
								else
								{
									collisionAxis = float3(0, 0, 1); //Just push them away along Z
								}
							}
							else
#endif//RESOLVE_ZERO_DISTANCE
							{
								distanceSq = sqrt(distanceSq);
								collisionAxis = collisionAxis / distanceSq;//normalize
							}
							float delta = (radiusSum - distanceSq) * 0.5f * CGConstantsCB.CollisionResponseCoef;
							collisionAxis*=delta;
							curObjectPos += collisionAxis;
							bCurObjectHadCollisions = true;
							
							uint otherObjectStateBitmask = ParticlesStateBufferSRV[otherObjectIndex];
							if((otherObjectStateBitmask & BITMASK_PARTICLE_STATE_PINNED) == 0)
							{
								float3 pFinal = otherObjectPos - collisionAxis;
								if(CGConstantsCB.bSimulateParticlesIn2D)
								{
									if(CGConstantsCB.bSimulateParticlesIn2DAllowNegativeZOffset)
									{
										pFinal.z = min(CGConstantsCB.DefaultParticlePos.z, pFinal.z);
									}
									else
									{
										pFinal.z = CGConstantsCB.DefaultParticlePos.z;
									}
									
								}
								ParticlesBufferPositionCurUAV[otherObjectIndex] = pFinal;
							}

#if DIFFUSE_PARTICLES_COLORS_ON_COLLISION
							const float ColorDiffusionCoeff = 0.1f;
							
							//diffuse colors
							float3 otherObjectColor = ParticlesColorsBufferUAV[otherObjectIndex];
							float3 color = (curObjectColor + otherObjectColor) * 0.5;
							curObjectColor = curObjectColor + (color - curObjectColor) * ColorDiffusionCoeff;
							ParticlesColorsBufferUAV[otherObjectIndex] = otherObjectColor + (color - otherObjectColor) * ColorDiffusionCoeff;
#endif//DIFFUSE_PARTICLES_COLORS_ON_COLLISION

						}
					}
				}
			}

			if(bCurObjectHadCollisions)
			{
				uint curObjectStateBitmask = ParticlesStateBufferSRV[curObjectIndex];
				if((curObjectStateBitmask & BITMASK_PARTICLE_STATE_PINNED) == 0)
				{
					if(CGConstantsCB.bSimulateParticlesIn2D)
					{
						if(CGConstantsCB.bSimulateParticlesIn2DAllowNegativeZOffset)
						{
							curObjectPos.z = min(CGConstantsCB.DefaultParticlePos.z, curObjectPos.z);
						}
						else
						{
							curObjectPos.z = CGConstantsCB.DefaultParticlePos.z;
						}
					}
					ParticlesBufferPositionCurUAV[curObjectIndex] = curObjectPos;
#if DIFFUSE_PARTICLES_COLORS_ON_COLLISION
					ParticlesColorsBufferUAV[curObjectIndex] = curObjectColor;
#endif//DIFFUSE_PARTICLES_COLORS_ON_COLLISION	
				}
			}

		}



	}


}




//===========================================================================
//                              CONSTRAINT
//===========================================================================

[numthreads( 64, 1, 1 )]
void ResolveConstraints( 
	uint SampleCoord : SV_DispatchThreadID )
{

	if(SampleCoord < GCollisionConstantsCB.NumCells /* NumConstraints */)
	{
		//return if broken
		if((ConstraintsStatesBufferSRV[SampleCoord] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			return;
		}

		//Get objects
		CConstraint curConstraint = ConstraintsBufferSRV[SampleCoord];

		float3 objectAPos = ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot1];
		float3 objectBPos = ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot2];

		float3 dir = objectBPos - objectAPos;

		float curDist = length(dir);

#if RESOLVE_ZERO_DISTANCE
		if(curDist < EPSILON)
		{
			curDist = EPSILON;
			dir = float3(0, 0, 1) * curDist; //Just push them away along Z
		}
#endif//RESOLVE_ZERO_DISTANCE

        float delta = curDist - CGConstantsCB.SpringRestLength;
		///TODO: Use different stiffness based on distance, e.g. if distance is < 50%ofRest, use hard Stifness, otherwise use soft
        dir = (dir / curDist) * delta * 0.5f * CGConstantsCB.ConnectionConstraintStiffness * GCollisionConstantsCB.PerPassSpringStiffnessScale;

		objectAPos += dir;
		objectBPos -= dir;

		if(CGConstantsCB.bSimulateParticlesIn2D)
		{
			if(CGConstantsCB.bSimulateParticlesIn2DAllowNegativeZOffset)
			{
				objectAPos.z = min(CGConstantsCB.DefaultParticlePos.z, objectAPos.z);
				objectBPos.z = min(CGConstantsCB.DefaultParticlePos.z, objectBPos.z);
			}
			else
			{
				objectAPos.z = CGConstantsCB.DefaultParticlePos.z;
				objectBPos.z = CGConstantsCB.DefaultParticlePos.z;
			}
		}

#if 1 //Pinned Objects
		uint objectAStateBitmask = ParticlesStateBufferSRV[curConstraint.IndexToKnot1];
		uint objectBStateBitmask = ParticlesStateBufferSRV[curConstraint.IndexToKnot2];
		if((objectAStateBitmask & BITMASK_PARTICLE_STATE_PINNED) == 0)
		{
			ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot1] = objectAPos;
		}

		if((objectBStateBitmask & BITMASK_PARTICLE_STATE_PINNED) == 0)
		{
			ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot2] = objectBPos;
		}
#else
		ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot1] = objectAPos;
		ParticlesBufferPositionCurUAV[curConstraint.IndexToKnot2] = objectBPos;
#endif
	}
}

//Soft Pin equivalent
[numthreads( 64, 1, 1 )]
void ResolveConstraintsToDefaultPos( 
	uint SampleCoord : SV_DispatchThreadID )
{

	if(SampleCoord < (CGConstantsCB.ClothSize.x * CGConstantsCB.ClothSize.y))
	{
		uint particleStateBitmask = ParticlesStateBufferSRV[SampleCoord];
		if((particleStateBitmask & BITMASK_PARTICLE_STATE_PINNED_SOFT) != 0)
		{
			float3 posCur = ParticlesBufferPositionCurUAV[SampleCoord];

			float3 posInitial = ParticlesBufferPositionAfterCreationSRV[SampleCoord];

			float3 dir = posInitial - posCur;

			float curDist = length(dir);

			float delta = curDist - CGConstantsCB.SpringRestLength;
		
        	dir = dir * CGConstantsCB.ConnectionConstraintToDefaultPosStiffness * CGConstantsCB.ConnectionConstraintStiffness;

			posCur += dir;

			//Write back
			ParticlesBufferPositionCurUAV[SampleCoord] = posCur;
		}

	}
	

}

//===========================================================================
//                              SPATIAL CONTROL
//===========================================================================

[numthreads( 64, 1, 1 )]
void SpatialControl_UpdateIntersectedPoint( 
	uint SampleCoord : SV_DispatchThreadID )
{
	if(SampleCoord < CGConstantsCB.NumParticles)
	{

		uint curParticleState = ParticlesStateBufferUAV[SampleCoord];
		//De-select all particles
		curParticleState = curParticleState & (~BITMASK_PARTICLE_STATE_SELECTED);

		bool bIntersectionInThisThread = false;
		float3 objectPos;

		if(IntersectedPointIndexUAV[0] == UINT_MAX)//intersection point still haven't been found
		{
			//Test Mouse Ray against object
			float3 objectPos = ParticlesBufferPositionCurUAV[SampleCoord];
			if(RaySphereIntersection(CGSpatialControlConstantsCB.RayOrigin, CGSpatialControlConstantsCB.RayEnd, objectPos, CGConstantsCB.PointRadius * CGSpatialControlConstantsCB.IntersectionSphereScale))
			{
				bIntersectionInThisThread = true;
				curParticleState = curParticleState | BITMASK_PARTICLE_STATE_SELECTED;
				//curParticleState = curParticleState | BITMASK_PARTICLE_STATE_PINNED;
				uint prevValue;
				InterlockedExchange(IntersectedPointIndexUAV[0], SampleCoord, prevValue);

				ClickedObjectViewSpaceZUAV[0] = objectPos;
			}
		}

		AllMemoryBarrierWithGroupSync();

		//In the end update the state
		ParticlesStateBufferUAV[SampleCoord] = curParticleState;

	}
}


[numthreads( 64, 1, 1 )]
void SpatialControl_CutConstraints( 
	uint SampleCoord : SV_DispatchThreadID )
{
	if(SampleCoord < GCollisionConstantsCB.NumCells /* NumConstraints */)
	{
		//return if already broken
		if((ConstraintsStatesBufferUAV[SampleCoord] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			return;
		}

		//get points
		//Get objects
		CConstraint curConstraint = ConstraintsBufferSRV[SampleCoord];

		float3 objectAPos = ParticlesBufferPositionCurSRV[curConstraint.IndexToKnot1];
		float3 objectBPos = ParticlesBufferPositionCurSRV[curConstraint.IndexToKnot2];

		if (LineSegmentSphereIntersection(objectAPos, objectBPos, CGSpatialControlConstantsCB.ConstraintScissorPos, CGSpatialControlConstantsCB.ScissorAffectRadius))
        {
            ConstraintsStatesBufferUAV[SampleCoord] = BITMASK_CONSTRAINT_STATE_BROKEN;
        }

	}
}

//===========================================================================
//                              MESH GENERATION
//===========================================================================

//Row major
uint GetFlatIndexToMeshArray(uint2 coord)
{
	return coord.y * GMeshGenerationConstantsCB.MeshSize.x + coord.x;
};
uint GetFlatIndexToClothArray(uint2 coord)
{
	return coord.y * CGConstantsCB.ClothSize.x + coord.x;
};



#if DISCARD_BROKEN_CONSTRAINTS
///column major
uint GetFirstConstraintStride()
{
	return CGConstantsCB.ClothSize.y;
}
//blue/yellow
uint GetFlatIndexToFirstConstraintArray(uint2 coord)
{
	coord.x = coord.x / 2;
	return coord.x * GetFirstConstraintStride() + coord.y;
};

///row major
uint GetSecondConstraintStride()
{
	return CGConstantsCB.ClothSize.x;
}
//green/red
uint GetFlatIndexToSecondConstraintArray(uint2 coord)
{
	coord.y = coord.y / 2;
	return coord.y * GetSecondConstraintStride() + coord.x;
};
#endif//DISCARD_BROKEN_CONSTRAINTS

#define COMPUTE_NORMALS_HQ 1 //wont work with ripped cloth

float3 GetVertexNormal(float3 p, uint2 coord)
{
	// Get the neighboring points in the u and v directions
	const float3 u1 = coord.x > 0 ? ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(coord - uint2(1,0))] : p;
	const float3 u2 = coord.x < (GMeshGenerationConstantsCB.MeshSize.x - 1) ? ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(coord + uint2(1,0))] : p;
	const float3 v1 = coord.y > 0 ? ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(coord - uint2(0,1))] : p;
	const float3 v2 = coord.y < (GMeshGenerationConstantsCB.MeshSize.y - 1) ? ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(coord + uint2(0,1))] : p;

    // Compute the central difference vectors in the u and v directions
    float3 du = u2 - u1;
    float3 dv = v2 - v1;
    //return normalize(cross(du, dv));
    return normalize(cross(dv, du));
}

[numthreads( 8, 8, 1 )]
void MeshGenerator( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	if((SampleCoord.x >= (GMeshGenerationConstantsCB.MeshSize.x - 1)) || (SampleCoord.y >= (GMeshGenerationConstantsCB.MeshSize.y - 1)))
	{
		return;
	}


	//For each quad we add 2 triangles
	//if any constraint is broken, don't add triangle to the list

	//triangle 1

	bool bShouldRenderTriangle = true;
	bool bShouldRenderSecondTriangle = true;
	
#if DISCARD_BROKEN_CONSTRAINTS
	//Get constraint color and index based on knot local pos
	uint2 knotLocalCoord = uint2(SampleCoord.x % 2, SampleCoord.y % 2);
	if(knotLocalCoord.x == 0)
	{
		//BLUE
		if((ConstraintsStatesBufferBLUESRV[GetFlatIndexToFirstConstraintArray(SampleCoord)] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderTriangle = false;
		}
		if((ConstraintsStatesBufferBLUESRV[GetFlatIndexToFirstConstraintArray(SampleCoord) + 1] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderSecondTriangle = false;
		}
	}
	else
	{
		//GREEN
		if((ConstraintsStatesBufferGREENSRV[GetFlatIndexToFirstConstraintArray(SampleCoord)] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderTriangle = false;
		}
		if((ConstraintsStatesBufferGREENSRV[GetFlatIndexToFirstConstraintArray(SampleCoord) + 1] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderSecondTriangle = false;
		}
	}
	if(knotLocalCoord.y == 0)
	{
		//YELLOW
		if((ConstraintsStatesBufferYELLOWSRV[GetFlatIndexToSecondConstraintArray(SampleCoord)] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderTriangle = false;
		}
		if((ConstraintsStatesBufferYELLOWSRV[GetFlatIndexToSecondConstraintArray(SampleCoord) + 1] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderSecondTriangle = false;
		}
	}
	else
	{
		//RED
		if((ConstraintsStatesBufferREDSRV[GetFlatIndexToSecondConstraintArray(SampleCoord)] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderTriangle = false;
		}
		if((ConstraintsStatesBufferREDSRV[GetFlatIndexToSecondConstraintArray(SampleCoord) + 1] & BITMASK_CONSTRAINT_STATE_BROKEN) != 0)
		{
			bShouldRenderSecondTriangle = false;
		}
	}
#endif//DISCARD_BROKEN_CONSTRAINTS

	//constraints valid, add triangle
	if(bShouldRenderTriangle || bShouldRenderSecondTriangle)
	{
		const float3 p = 	  	ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(SampleCoord + uint2(0,0))];
		const float3 pRight = 	ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(SampleCoord + uint2(1,0))];
		const float3 pUp = 	  	ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(SampleCoord + uint2(0,1))];
		const float3 pRightUp = ParticlesBufferPositionCurSRV[GetFlatIndexToMeshArray(SampleCoord + uint2(1,1))];

		//Compute normal for each vertex in the quad
#if COMPUTE_NORMALS_HQ
		float3 n = 			GetVertexNormal(p, SampleCoord + uint2(0,0));
		float3 nRight = 	GetVertexNormal(pRight, SampleCoord + uint2(1,0));
		float3 nUp = 		GetVertexNormal(pUp, SampleCoord + uint2(0,1));
		float3 nRightUp = 	GetVertexNormal(pRightUp, SampleCoord + uint2(1,1));
#else
		//clockwise cross
		float3 n = normalize(cross(pUp - p, pRight-p));
		float3 nRight = normalize(cross(p - pRight, pRightUp-pRight));
		float3 nUp = normalize(cross(pRightUp - pUp, p - pUp));
		float3 nRightUp = normalize(cross(pRight - pRightUp, pUp - pRightUp));
#endif//COMPUTE_NORMALS_HQ

		//Compute TexCoords
		const float2 meshSize = GMeshGenerationConstantsCB.MeshSize - uint2(1,1);
		float2 uv = float2(SampleCoord) / meshSize;
		float2 uvRight = float2(SampleCoord + uint2(1,0)) / meshSize;
		float2 uvUp = float2(SampleCoord + uint2(0,1)) / meshSize;
		float2 uvRightUp = float2(SampleCoord + uint2(1,1)) / meshSize;
		uv.y = 1.f - uv.y;
		uvRight.y = 1.f - uvRight.y;
		uvUp.y = 1.f - uvUp.y;
		uvRightUp.y = 1.f - uvRightUp.y;
		
		if(bShouldRenderTriangle)
		{
			//append left right up vertices
			uint prevIndex;
			InterlockedAdd(InterlockedAddBuffer[0], 3, prevIndex);

			MeshBufferVerticesUAV[prevIndex + 0] = p;
			MeshBufferVerticesUAV[prevIndex + 1] = pUp;
			MeshBufferVerticesUAV[prevIndex + 2] = pRight;

			MeshBufferNormalsUAV[prevIndex + 0] = n;
			MeshBufferNormalsUAV[prevIndex + 1] = nUp;
			MeshBufferNormalsUAV[prevIndex + 2] = nRight;

			MeshBufferTexCoordsUAV[prevIndex + 0] = uv;
			MeshBufferTexCoordsUAV[prevIndex + 1] = uvUp;
			MeshBufferTexCoordsUAV[prevIndex + 2] = uvRight;
		}

		if(bShouldRenderSecondTriangle)
		{
			uint prevIndex;
			InterlockedAdd(InterlockedAddBuffer[0], 3, prevIndex);

			MeshBufferVerticesUAV[prevIndex + 0] = pRightUp;
			MeshBufferVerticesUAV[prevIndex + 1] = pRight;
			MeshBufferVerticesUAV[prevIndex + 2] = pUp;

			MeshBufferNormalsUAV[prevIndex + 0] = nRightUp;
			MeshBufferNormalsUAV[prevIndex + 1] = nRight;
			MeshBufferNormalsUAV[prevIndex + 2] = nUp;

			MeshBufferTexCoordsUAV[prevIndex + 0] = uvRightUp;
			MeshBufferTexCoordsUAV[prevIndex + 1] = uvRight;
			MeshBufferTexCoordsUAV[prevIndex + 2] = uvUp;
		}
	}

}




//===========================================================================
//                              VELOCITY FIELD
//===========================================================================


float2 GetVelocityFieldCellStartPosWorldSpaceFromCoord(uint2 coord, bool bInvertY)
{
	if(bInvertY)
	{
		coord.y = (CGConstantsCB.VelocityFieldGridSize - 1) - coord.y;
	}
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset);
	return float2( boundsStartPos.x + coord.x * h, boundsStartPos.y + coord.y * h );
}
float2 GetVelocityCellCenterWorldSpaceFromCoord(uint2 coord, bool bInvertY)
{
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord, bInvertY);
	return curPos + float2(h * 0.5f, h * 0.5f );
}
float2 GetVelocityXComponentWorldPosFromCoord(uint2 coord, bool bInvertY)
{
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord, bInvertY);
	curPos.y += h * 0.5f;
	return curPos;
}
float2 GetVelocityYComponentWorldPosFromCoord(uint2 coord, bool bInvertY)
{
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	float2 curPos = GetVelocityFieldCellStartPosWorldSpaceFromCoord(coord, bInvertY);
	curPos.x += h * 0.5f;
	return curPos;
}

float2 LoadFromVelocityField(uint2 coord)
{
	float u = VelocityField_U_ComponentTextureSRV[coord];
	float v = VelocityField_V_ComponentTextureSRV[coord];
	return float2(u,v);
}

float LoadFromVelocityField_U(uint2 coord)
{
	return VelocityField_U_ComponentTextureSRV[coord];
}
float LoadFromVelocityField_V(uint2 coord)
{
	return VelocityField_V_ComponentTextureSRV[coord];
}

//return bool2 per each velocity component
bool2 bCanProcessCurrentFluidCellVelocity(uint2 SampleCoord)
{
	//velocities between air cells are undefined
	//we can process the cell if at least 1 of 2 cells has density
	bool2 result;
	result.x = VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_DENSITY || VelocityFieldStateTextureSRV[SampleCoord - uint2(1,0)] == VELOCITY_FIELD_CELL_STATE_DENSITY;
	result.y = VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_DENSITY || VelocityFieldStateTextureSRV[SampleCoord + uint2(0,1)] == VELOCITY_FIELD_CELL_STATE_DENSITY;
	return result;
}


[numthreads( 64, 1, 1 )]
void ApplyPointerVelocityToParticles( 
	uint SampleCoord : SV_DispatchThreadID )
{
	//const float2 pointerPos = GetVelocityCellCenterWorldSpaceFromCoord(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo, false);
	const float2 pointerPos = CGSpatialControlConstantsCB.ConstraintScissorPos.xy;
	const float poniterRadius = CGConstantsCB.VelocityGeneratorRadius;

	const float3 curParticlePos = ParticlesBufferPositionCurSRV[SampleCoord];
	const float particleRadius = CGConstantsCB.PointRadius;

	if(SphereSphereIntersection(
		float3(curParticlePos.x, curParticlePos.y, 0), particleRadius,
		float3(pointerPos.x, pointerPos.y, 0), poniterRadius)
		)
	{
		//apply velocity to particle
		float3 curParticleVel = CurParticlesVelocityBufferUAV[SampleCoord];
		//curParticleVel.xy = curParticleVel.xy + CGConstantsCB.CurPointerVelocity;
		curParticleVel.xy = curParticleVel.xy + CGSpatialControlConstantsCB.MouseOffsetWorld.xy * CGConstantsCB.GeneratedVelocityMagnitudeScale;
		CurParticlesVelocityBufferUAV[SampleCoord] = curParticleVel;
	}

}

[numthreads( 8, 8, 1 )]
void GenerateVelocityField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const float2 cellPos = GetVelocityCellCenterWorldSpaceFromCoord(SampleCoord, false);
	const float cellRadius = CGConstantsCB.VelocityFieldGridCellLength * 0.5f;//inner sphere

	const float2 pointerPos = GetVelocityCellCenterWorldSpaceFromCoord(CGConstantsCB.CurVelocityFieldCellIdMousePointsTo, false).xy;
	
	//apply if current sample is in pointer range

#if GENERATE_VELOCITY
	if(VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		if(SphereSphereIntersection(float3(cellPos.x, cellPos.y, 0), cellRadius,
		float3(pointerPos.x, pointerPos.y, 0), CGConstantsCB.VelocityGeneratorRadius))
		{
			const uint2 texcoord = uint2(SampleCoord.x, (CGConstantsCB.VelocityFieldGridSize-1) - SampleCoord.y);
			float2 curVelocity;
			curVelocity.x = VelocityField_U_ComponentTextureUAV[texcoord];
			curVelocity.y = VelocityField_V_ComponentTextureUAV[texcoord];
			curVelocity = curVelocity + CGConstantsCB.CurPointerVelocity;
			VelocityField_U_ComponentTextureUAV[texcoord] = curVelocity.x;
			VelocityField_V_ComponentTextureUAV[texcoord] = curVelocity.y;
		}
	}
#endif//GENERATE_VELOCITY

#if GENERATE_DENSITY
	if(SphereSphereIntersection(float3(cellPos.x, cellPos.y, 0), cellRadius,
	float3(pointerPos.x, pointerPos.y, 0), CGConstantsCB.DensityGeneratorRadius))
	{
		const uint2 texcoord = uint2(SampleCoord.x, (CGConstantsCB.VelocityFieldGridSize-1) - SampleCoord.y);
		DensityFieldTextureUAV[texcoord] = CGConstantsCB.InitialDensityColor;
	}
#endif//GENERATE_VELOCITY
}

[numthreads( 8, 8, 1 )]
void GenerateVelocityFieldBoundaries( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const bool borderLeft = /* bEnabledBorders[0] &&  */SampleCoord.x == 0;
	const bool borderUp = /* bEnabledBorders[1] &&  */SampleCoord.y == (CGConstantsCB.VelocityFieldGridSize - 1);
	const bool borderRight = /* bEnabledBorders[2] && */ SampleCoord.x == (CGConstantsCB.VelocityFieldGridSize - 1);
	const bool borderDown = /* bEnabledBorders[3] && */ SampleCoord.y == 0;
	bool isBorder = 0;
	isBorder |= borderLeft;
	isBorder |= borderRight;
	isBorder |= borderDown;
	//isBorder |= borderUp;
	if (isBorder)
	{
		VelocityFieldStateTextureUAV[SampleCoord] = VELOCITY_FIELD_CELL_STATE_SOLID;
	}
	else
	{
#if CLEAR_STATES
		VelocityFieldStateTextureUAV[SampleCoord] = CGConstantsCB.bUseParticleBasedAdvection ? VELOCITY_FIELD_CELL_STATE_AIR : VELOCITY_FIELD_CELL_STATE_DENSITY;
#endif//CLEAR_STATES

#if PAINT_RENDER_PIXEL_ADVECT
	float2 sampleCoordf = float2(SampleCoord) + float2(0.5,0.5);
	float2 fieldCoordUV = MapToRange(sampleCoordf, 0.f, CGConstantsCB.VelocityFieldGridSize - 1, 0.f, 1.f);
	float metal = InputTextureWithNegativeValuesSRV.SampleLevel(LinearSampler, fieldCoordUV, 0);
	if(metal > 0.f)
	{
		VelocityFieldStateTextureUAV[SampleCoord] = VELOCITY_FIELD_CELL_STATE_SOLID;
		return;
	}
#endif

	}
}

[numthreads( 8, 8, 1 )]
void FillDensityField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
#if COPY_TEXTURE
	float2 texCoord = (float2(SampleCoord) + float2(0.5,0.5)) / float(CGConstantsCB.VelocityFieldGridSize).xx;
	float3 color = DensitySourceColorTextureSRV.SampleLevel(LinearSampler, texCoord, 0);
	DensityFieldTextureUAV[SampleCoord] = color;
#else//solid color
	DensityFieldTextureUAV[SampleCoord] = CGConstantsCB.InitialDensityColor;
#endif//COPY_TEXTURE
}

#define SCALE_DERIVATIVES 0

#define USE_FORWARD_DIFFERENCE 1

#define SIMULATE_DOWNWARD_PUSH 1
#define SIMULATE_UPWARD_PUSH 1

uint2 DecodeUint(uint val)
{
	return uint2(val & 0xFFFF, (val >> 16) & 0xFFFF);
}

[numthreads( 8, 8, 1 )]
void ProjectVelocityField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	uint2 coordOffset = DecodeUint(VectorProjectionPerPassConstantsCB.Offset);
	const uint2 coordStride = uint2(2,2);
	SampleCoord = coordOffset + SampleCoord * coordStride;

	//invert .y
	SampleCoord.y = (CGConstantsCB.VelocityFieldGridSize - 1) - SampleCoord.y;

	uint curState = VelocityFieldStateTextureSRV[SampleCoord];

	float2 uvCur;
	uvCur.x = VelocityField_U_ComponentTextureUAV[SampleCoord];
	uvCur.y = VelocityField_V_ComponentTextureUAV[SampleCoord];

	float divergence = 0.f;

	if(curState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		uint leftState = VelocityFieldStateTextureSRV[SampleCoord - uint2(1,0)] != VELOCITY_FIELD_CELL_STATE_SOLID;
		uint rightState = VelocityFieldStateTextureSRV[SampleCoord + uint2(1,0)] != VELOCITY_FIELD_CELL_STATE_SOLID;
		uint downState = VelocityFieldStateTextureSRV[SampleCoord + uint2(0,1)] != VELOCITY_FIELD_CELL_STATE_SOLID;
		uint upState = VelocityFieldStateTextureSRV[SampleCoord - uint2(0,1)] != VELOCITY_FIELD_CELL_STATE_SOLID;

		//compute divisor
		uint divisor = leftState + rightState + upState + downState;
		if(divisor == 0) return;
	
		float uRight = VelocityField_U_ComponentTextureUAV[SampleCoord + uint2(1,0)];
		float vUp = VelocityField_V_ComponentTextureUAV[SampleCoord - uint2(0,1)];

		divergence = uRight - uvCur.x + vUp - uvCur.y;

		divergence = divergence / float(divisor);

		if(CGConstantsCB.bSolversUseOverrelaxation)
		{
			divergence = divergence * CGConstantsCB.OverrelaxationConstant;
		}

		if (CGConstantsCB.bUseParticleBasedAdvection) 
		{
			const float averageDensity = AverageDensityBufferUAV[0];
			if(averageDensity > 0.f)
			{
				//reduce divergence in dense regions
				float pressure = AdvectedParticleDensityTextureUAV[SampleCoord] - averageDensity;
				if (pressure > 0.0)
				{
					divergence = divergence - CGConstantsCB.AdvectedParticleDriftCompensationScale * pressure;
				}
			}
		}

		if(leftState)
		{
			uvCur.x += divergence;
			VelocityField_U_ComponentTextureUAV[SampleCoord] = uvCur.x;
		}
		if(downState)
		{
			uvCur.y += divergence;
			VelocityField_V_ComponentTextureUAV[SampleCoord] = uvCur.y;
		}
		if(rightState)
		{
			uRight -= divergence;
			VelocityField_U_ComponentTextureUAV[SampleCoord + uint2(1,0)] = uRight;
		}
		if(upState)
		{
			vUp -= divergence;
			VelocityField_V_ComponentTextureUAV[SampleCoord - uint2(0,1)] = vUp;
		}

	}

}

float LoadFromVelocityField_UInv(uint2 coord)
{
	coord.y = (CGConstantsCB.VelocityFieldGridSize - 1) - coord.y;
	return VelocityField_U_ComponentTextureSRV[coord];
}
float LoadFromVelocityField_VInv(uint2 coord)
{
	coord.y = (CGConstantsCB.VelocityFieldGridSize - 1) - coord.y;
	return VelocityField_V_ComponentTextureSRV[coord];
}


float2 SampleVelocityField(float2 worldPos)
{
	float texelSize = 1.f / (CGConstantsCB.VelocityFieldGridSize);

	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const uint n = CGConstantsCB.VelocityFieldGridSize;
	const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset);
	const float2 minLimitPos = boundsStartPos + float2( h, h );
	const float2 maxLimitPos = boundsStartPos + float2( h * n, h * n );
	worldPos = clamp(worldPos, minLimitPos, maxLimitPos);

	float2 res;
	//get u
	{
		// map world pos to 0-1
		float2 velocityGridCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
		velocityGridCoordUV.x += (texelSize * 0.5f);	
		//invert .y
		velocityGridCoordUV.y = 1.f - velocityGridCoordUV.y;	
		res.x = VelocityField_U_ComponentTextureSRV.SampleLevel(LinearSampler, velocityGridCoordUV, 0);
	}
	
	//get v
	{
		// map world pos to 0-1
		float2 velocityGridCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
		velocityGridCoordUV.y += (texelSize * 0.5f);	
		//invert .y
		velocityGridCoordUV.y = 1.f - velocityGridCoordUV.y;

		res.y = VelocityField_V_ComponentTextureSRV.SampleLevel(LinearSampler, velocityGridCoordUV, 0);
	}

	return res;

}


void GatherVelocityFieldMinMax(float2 worldPos, out float2 uMinMax, out float2 vMinMax)
{
	float texelSize = 1.f / (CGConstantsCB.VelocityFieldGridSize);

	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const uint n = CGConstantsCB.VelocityFieldGridSize;
	const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset);
	const float2 minLimitPos = boundsStartPos + float2( h, h );
	const float2 maxLimitPos = boundsStartPos + float2( h * n, h * n );
	worldPos = clamp(worldPos, minLimitPos, maxLimitPos);

	//get u
	{
		// map world pos to 0-1
		float2 velocityGridCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
		velocityGridCoordUV.x += (texelSize * 0.5f);	
		//invert .y
		velocityGridCoordUV.y = 1.f - velocityGridCoordUV.y;	
		float4 uAll = VelocityField_U_ComponentTextureSRV.Gather(LinearSampler, velocityGridCoordUV, 0);

		uMinMax.x = min(uAll.x, min(uAll.y, min(uAll.z, uAll.w)));
		uMinMax.y = max(uAll.x, max(uAll.y, max(uAll.z, uAll.w)));
	}
	
	//get v
	{
		// map world pos to 0-1
		float2 velocityGridCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
		velocityGridCoordUV.y += (texelSize * 0.5f);	
		//invert .y
		velocityGridCoordUV.y = 1.f - velocityGridCoordUV.y;

		float4 vAll = VelocityField_V_ComponentTextureSRV.Gather(LinearSampler, velocityGridCoordUV, 0);

		vMinMax.x = min(vAll.x, min(vAll.y, min(vAll.z, vAll.w)));
		vMinMax.y = max(vAll.x, max(vAll.y, max(vAll.z, vAll.w)));
	}

}


[numthreads( 8, 8, 1 )]
void ComputeDivergence( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const float samplingPeriodLength = CGConstantsCB.VelocityFieldGridCellLength;

	uint curState = VelocityFieldStateTextureSRV[SampleCoord];

	if(curState != VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		VectorFieldDivergenceTextureUAV[SampleCoord] = 0.f;
		return;
	}

	float2 uvCur = LoadFromVelocityField(SampleCoord);
	float uRight = LoadFromVelocityField_U(SampleCoord + uint2( 1, 0 ));
	float vUp = LoadFromVelocityField_V(SampleCoord - uint2( 0, 1 ));

#if 1
	float div = -0.5f * samplingPeriodLength * (uRight - uvCur.x + vUp - uvCur.y); //divergence computation from Jos Stam paper
#else
	float div = ((uRight - uvCur.x) + (vUp - uvCur.y)) / (-2.f * samplingPeriodLength);//simple central diff
#endif

	if (CGConstantsCB.bUseParticleBasedAdvection) 
	{
		const float averageDensity = AverageDensityBufferUAV[0];
		if(averageDensity > 0.f)
		{
			//reduce divergence in dense regions
			float pressure = AdvectedParticleDensityTextureUAV[SampleCoord] - averageDensity;
			if (pressure > 0.0)
			{
				div = div + CGConstantsCB.AdvectedParticleDriftCompensationScale * pressure;
			}
		}
	}

	VectorFieldDivergenceTextureUAV[SampleCoord] = div;

}


//TODO: Use prev frame pressure values as initial guess for iterator 
[numthreads( 8, 8, 1 )]
void ComputePressure( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	if(VelocityFieldStateTextureSRV[SampleCoord] != VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		PressureFieldTextureUAV[SampleCoord] = 0.f;
		return;
	}

	//load pressure srv
	const float p0 = PressureFieldTextureSRV[SampleCoord - uint2(1,0)];
	const float p1 = PressureFieldTextureSRV[SampleCoord + uint2(1,0)];
	const float p2 = PressureFieldTextureSRV[SampleCoord - uint2(0,1)];
	const float p3 = PressureFieldTextureSRV[SampleCoord + uint2(0,1)];

	//load divergence 
	const float d = VectorFieldDivergenceTextureSRV[SampleCoord];

	//write to pressure UAV
	PressureFieldTextureUAV[SampleCoord] = (d + p0 + p1 + p2 + p3) / 4.f;
}

[numthreads( 64, 1, 1 )]
void CopyPressureToBorders( 
	uint i : SV_DispatchThreadID )
{
	uint n = CGConstantsCB.VelocityFieldGridSize;

	if(i < (n-2))
	{
		PressureFieldTextureUAV[uint2(i+1,0)] = PressureFieldTextureUAV[uint2(i+1,1)];
		PressureFieldTextureUAV[(uint2( i+1, n - 1 ))] = PressureFieldTextureUAV[(uint2( i+1, n - 2 ))];
		PressureFieldTextureUAV[(uint2( 0, i+1 ))] = PressureFieldTextureUAV[(uint2( 1, i+1 ))];
		PressureFieldTextureUAV[(uint2( n - 1, i+1 ))] = PressureFieldTextureUAV[(uint2( n - 2, i+1 ))];
	}
	else if(i == (n-2))
	{
		PressureFieldTextureUAV[(uint2( 0, 0 ))] = PressureFieldTextureUAV[(uint2( 1, 1 ))];
		PressureFieldTextureUAV[(uint2( 0, n-1 ))] = PressureFieldTextureUAV[(uint2( 1, n-2 ))];
		PressureFieldTextureUAV[(uint2( n-1, 0))] = PressureFieldTextureUAV[(uint2( n-2, 1 ))];
		PressureFieldTextureUAV[(uint2( n-1, n-1))] = PressureFieldTextureUAV[(uint2( n-2, n-2 ))];
	}
}

[numthreads( 8, 8, 1 )]
void ProjectPressureField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	const float samplingPeriodLength = CGConstantsCB.VelocityFieldGridCellLength;

	if (VelocityFieldStateTextureSRV[SampleCoord] != VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		PressureGradientTextureUAV[SampleCoord] = float2(0,0);
		return;
	}

	float2 curPressureGradient;
	curPressureGradient.x = (PressureFieldTextureSRV[SampleCoord] - PressureFieldTextureSRV[SampleCoord - uint2(1,0)]) / samplingPeriodLength;
	curPressureGradient.y = (PressureFieldTextureSRV[SampleCoord] - PressureFieldTextureSRV[SampleCoord + uint2(0,1)]) / samplingPeriodLength;
	PressureGradientTextureUAV[SampleCoord] = curPressureGradient;

	uint leftState = VelocityFieldStateTextureSRV[SampleCoord - uint2(1,0)];
	uint downState = VelocityFieldStateTextureSRV[SampleCoord + uint2(0,1)];

	float relaxationConstant = 0.5f;
	if(CGConstantsCB.bSolversUseOverrelaxation)
	{
		relaxationConstant = CGConstantsCB.OverrelaxationConstant;
	}

	//staggered grid version:
	if(leftState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		float u = VelocityField_U_ComponentTextureUAV[SampleCoord];
		u -= relaxationConstant * curPressureGradient.x;
		VelocityField_U_ComponentTextureUAV[SampleCoord] = u;
	}

	if(downState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		float v = VelocityField_V_ComponentTextureUAV[SampleCoord];
		v -= relaxationConstant * curPressureGradient.y;
		VelocityField_V_ComponentTextureUAV[SampleCoord] = v;
	}

}

[numthreads( 8, 8, 1 )]
void AdvectVelocity( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	//borders can exit right away
	if(any(SampleCoord >= (CGConstantsCB.VelocityFieldGridSize - 1)) || any(SampleCoord == 0))
	{
		return;
	}

	const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset);

	uint curState = VelocityFieldStateTextureSRV[SampleCoord];

	if(curState == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		return;
	}

	const float2 curVelocity = LoadFromVelocityField(SampleCoord);

	//u
	const bool bLeftCellIsSolid = VelocityFieldStateTextureSRV[SampleCoord - uint2(1, 0)] == VELOCITY_FIELD_CELL_STATE_SOLID;
	if(!bLeftCellIsSolid && SampleCoord.y < (CGConstantsCB.VelocityFieldGridSize - 1))
	{
		const float u = curVelocity.x;

		//get average v
		const float vCur = curVelocity.y;
		const float vUp = LoadFromVelocityField_V(SampleCoord - uint2(0, 1));
		const float vLeft = LoadFromVelocityField_V(SampleCoord - uint2(1, 0));
		const float vLeftUp = LoadFromVelocityField_V(SampleCoord - uint2(1, 0) - uint2(0, 1));
		const float vAvrg = (vCur + vUp + vLeft + vLeftUp) * 0.25f;
		float2 velocityFull = float2(u, vAvrg);

		//get u world space pos
		const float2 curPos = GetVelocityXComponentWorldPosFromCoord(SampleCoord, true);

		const float2 prevPos = curPos - velocityFull * CGConstantsCB.DeltaTime * CGConstantsCB.VelocityAdvectionScaleFactor;

		float2 velocityNew = SampleVelocityField(prevPos);

		VelocityField_U_ComponentTextureUAV[SampleCoord] = velocityNew.x;

	}

	//v
	const bool bDownCellIsSolid = VelocityFieldStateTextureSRV[SampleCoord + uint2(0, 1)] == VELOCITY_FIELD_CELL_STATE_SOLID;
	if(!bDownCellIsSolid && SampleCoord.x < (CGConstantsCB.VelocityFieldGridSize - 1))
	{
		const float v = curVelocity.y;

		//get average u
		const float uCur = curVelocity.x;
		const float uDown =  LoadFromVelocityField_U(SampleCoord + uint2(0, 1));
		const float uRight = LoadFromVelocityField_U(SampleCoord + uint2(1, 0));
		const float uRightDown = LoadFromVelocityField_U(SampleCoord + uint2(1, 0) + uint2(0, 1));
		const float uAvrg = (uCur + uDown + uRight + uRightDown) * 0.25f;
		float2 velocityFull = float2(uAvrg, v);

		//get v world space pos
		const float2 curPos = GetVelocityYComponentWorldPosFromCoord(SampleCoord, true);
		
		const float2 prevPos = curPos - velocityFull * CGConstantsCB.DeltaTime * CGConstantsCB.VelocityAdvectionScaleFactor;

		float2 velocityNew = SampleVelocityField(prevPos);

		VelocityField_V_ComponentTextureUAV[SampleCoord] = velocityNew.y;
	}


}

DENSITY_TYPE SampleDensityField(float2 worldPos)
{
    const float h = CGConstantsCB.VelocityFieldGridCellLength;
    const uint n = CGConstantsCB.VelocityFieldGridSize;
    const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset, CGConstantsCB.VelocityFieldBoundsStartOffset);
    const float2 minLimitPos = boundsStartPos + float2(h, h);
    const float2 maxLimitPos = boundsStartPos + float2(h * n, h * n);

    worldPos = clamp(worldPos, minLimitPos, maxLimitPos);

    // Map world coords to UV coords
    float2 fieldCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
    // Invert .y
	#if PAINT_RENDER_PIXEL_ADVECT == 0
    fieldCoordUV.y = 1.f - fieldCoordUV.y;
	#endif
    return DensityFieldTextureSRV.SampleLevel(LinearSampler, fieldCoordUV, 0);
}

void GatherDensityFieldMinMax(float2 worldPos, out float2 redMinMax, out float2 greenMinMax, out float2 blueMinMax)
{
    const float h = CGConstantsCB.VelocityFieldGridCellLength;
    const uint n = CGConstantsCB.VelocityFieldGridSize;
    const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset, CGConstantsCB.VelocityFieldBoundsStartOffset);
    const float2 minLimitPos = boundsStartPos + float2(h, h);
    const float2 maxLimitPos = boundsStartPos + float2(h * n, h * n);

    worldPos = clamp(worldPos, minLimitPos, maxLimitPos);

    // Map world coords to UV coords
    float2 fieldCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
    // Invert .y
    fieldCoordUV.y = 1.f - fieldCoordUV.y;
    float4 red = DensityFieldTextureSRV.GatherRed(LinearSampler, fieldCoordUV, 0);
    float4 green = DensityFieldTextureSRV.GatherGreen(LinearSampler, fieldCoordUV, 0);
    float4 blue = DensityFieldTextureSRV.GatherBlue(LinearSampler, fieldCoordUV, 0);

	redMinMax.x = min(red.x, min(red.y, min(red.z, red.w)));
	redMinMax.y = max(red.x, max(red.y, max(red.z, red.w)));

	greenMinMax.x = min(green.x, min(green.y, min(green.z, green.w)));
	greenMinMax.y = max(green.x, max(green.y, max(green.z, green.w)));

	blueMinMax.x = min(blue.x, min(blue.y, min(blue.z, blue.w)));
	blueMinMax.y = max(blue.x, max(blue.y, max(blue.z, blue.w)));

}

[numthreads( 8, 8, 1 )]
void AdvectDensity( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	uint2 veloictySamplingCoord = SampleCoord;
	#if PAINT_RENDER_PIXEL_ADVECT
	veloictySamplingCoord = SampleCoord / 2;
	#endif

	#if PAINT_RENDER_PIXEL_ADVECT
	//get UV
	float2 sampleCoordf = float2(SampleCoord) + float2(0.5,0.5);
	const float densitySize = 1024.f;
	float2 fieldCoordUV = MapToRange(sampleCoordf, 0.f, densitySize, 0.f, 1.f);

	float metal = InputTextureWithNegativeValuesSRV.SampleLevel(LinearSampler, fieldCoordUV, 0);
	if(metal > 0.f)
	{
		DensityFieldTextureUAV[SampleCoord] = DensityFieldTextureSRV[SampleCoord];
		return;
	}
	#endif//PAINT_RENDER_PIXEL_ADVECT

	if(VelocityFieldStateTextureSRV[veloictySamplingCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		return;
	}

#if 1 && PAINT_RENDER_PIXEL_ADVECT
	float2 fieldCoordUVInv = float2(fieldCoordUV.x, 1.f - fieldCoordUV.y);
	float2 curVelocity;
	curVelocity.x = VelocityField_U_ComponentTextureSRV.SampleLevel(LinearSampler, fieldCoordUVInv, 0);
	curVelocity.y = VelocityField_V_ComponentTextureSRV.SampleLevel(LinearSampler, fieldCoordUVInv, 0);
#else
	//get cur velocity of cell center 
	float2 curVelocity;
	curVelocity.x = (VelocityField_U_ComponentTextureSRV[veloictySamplingCoord] + VelocityField_U_ComponentTextureSRV[veloictySamplingCoord + uint2(1, 0)]) * 0.5f;
	curVelocity.y = (VelocityField_V_ComponentTextureSRV[veloictySamplingCoord] + VelocityField_V_ComponentTextureSRV[veloictySamplingCoord - uint2(0, 1)]) * 0.5f;
#endif

	//get worldPos of the cell
#if PAINT_RENDER_PIXEL_ADVECT
	//map to particle simulation bounds
	float2 curPosWorld = MapToRange(fieldCoordUV, 0.f, 1.f, GetSimulationBoundsStart(), GetSimulationBoundsEnd());
#else
	float2 curPosWorld = GetVelocityCellCenterWorldSpaceFromCoord(veloictySamplingCoord, true);
#endif

#if BACKWARD_ADVECTION == 1
	float2 newPosWorld = curPosWorld - curVelocity * CGConstantsCB.DeltaTime;
#else
	float2 newPosWorld = curPosWorld + curVelocity * CGConstantsCB.DeltaTime;
#endif//BACKWARD_ADVECTION
	
	//sample density at new location
	DENSITY_TYPE advectedDensity = SampleDensityField(newPosWorld);

	DENSITY_TYPE newDensity = advectedDensity;

	DensityFieldTextureUAV[SampleCoord] = newDensity;
}

[numthreads( 8, 8, 1 )]
void AdvectDensityMacCormackFinal( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	if(VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		return;
	}

	DENSITY_TYPE curDensity = DensityFieldTextureSRV[SampleCoord];

	DENSITY_TYPE correctionDensity = DensityFieldMacCormackReversedTextureSRV[SampleCoord];

	DENSITY_TYPE newDensity = DensityFieldMacCormackAdvectedTextureSRV[SampleCoord];

	DENSITY_TYPE error = (curDensity - correctionDensity);
	newDensity = newDensity + 0.5 * error;

#if 1 //clamp
	float2 curVelocity;
	curVelocity.x = (VelocityField_U_ComponentTextureSRV[SampleCoord] + VelocityField_U_ComponentTextureSRV[SampleCoord + uint2(1, 0)]) * 0.5f;
	curVelocity.y = (VelocityField_V_ComponentTextureSRV[SampleCoord] + VelocityField_V_ComponentTextureSRV[SampleCoord - uint2(0, 1)]) * 0.5f;

	float2 curPosWorld = GetVelocityCellCenterWorldSpaceFromCoord(SampleCoord, true);
	float2 newPosWorld = curPosWorld - curVelocity * CGConstantsCB.DeltaTime;

	float2 redMinMax;
	float2 greenMinMax;
	float2 blueMinMax;
	GatherDensityFieldMinMax(newPosWorld, redMinMax, greenMinMax, blueMinMax);

	newDensity.r = clamp(newDensity.r, redMinMax.x, redMinMax.y);
	newDensity.g = clamp(newDensity.g, greenMinMax.x, greenMinMax.y);
	newDensity.b = clamp(newDensity.b, blueMinMax.x, blueMinMax.y);
#endif

	DensityFieldTextureUAV[SampleCoord] = newDensity;

}



[numthreads( 8, 8, 1 )]
void DiffuseVelocity( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	const float a = CGConstantsCB.DeltaTime * (CGConstantsCB.DiffusionFactor * 10.f);

	float2 curVelocity = LoadFromVelocityField(SampleCoord);

	if (VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		VelocityField_U_ComponentTextureUAV[SampleCoord] = curVelocity.x;
		VelocityField_V_ComponentTextureUAV[SampleCoord] = curVelocity.y;
		return;
	}
#if 1
	VelocityField_U_ComponentTextureUAV[SampleCoord] =
		(curVelocity.x +
			(
				VelocityField_U_ComponentTextureNewSRV[SampleCoord - uint2(1,0)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord + uint2(1,0)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord - uint2(0,1)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord + uint2(0,1)]
			) * a
		) / (1 + 4 * a);
	
	VelocityField_V_ComponentTextureUAV[SampleCoord] =
		(curVelocity.y +
			(
				VelocityField_V_ComponentTextureNewSRV[SampleCoord - uint2(1,0)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord + uint2(1,0)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord - uint2(0,1)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord + uint2(0,1)]
			) * a
		) / (1 + 4 * a);
#else

	VelocityField_U_ComponentTextureUAV[SampleCoord] =
		(curVelocity.x * 0.99 +
			(
				VelocityField_U_ComponentTextureNewSRV[SampleCoord - uint2(1,0)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord + uint2(1,0)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord - uint2(0,1)]
				+ VelocityField_U_ComponentTextureNewSRV[SampleCoord + uint2(0,1)]
			) 
		);
	
	VelocityField_V_ComponentTextureUAV[SampleCoord] =
		(curVelocity.y * 0.99 +
			(
				VelocityField_V_ComponentTextureNewSRV[SampleCoord - uint2(1,0)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord + uint2(1,0)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord - uint2(0,1)]
				+ VelocityField_V_ComponentTextureNewSRV[SampleCoord + uint2(0,1)]
			) 
		);

#endif

}

[numthreads( 8, 8, 1 )]
void DiffuseDensity( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	const float a = CGConstantsCB.DeltaTime * (CGConstantsCB.DiffusionFactor * 10.f);

	if (VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		DensityFieldTextureUAV[SampleCoord] = DensityFieldTextureSRV[SampleCoord];
		return;
	}

	DENSITY_TYPE curDensity = DensityFieldTextureSRV[SampleCoord];

	DensityFieldTextureUAV[SampleCoord] =
		(curDensity +
			(
				/* VelocityFieldStateTextureSRV[SampleCoord - uint2(1,0)] == 0 ? curDensity :  */DensityFieldTextureNewSRV[SampleCoord - uint2(1,0)]
				+ /* VelocityFieldStateTextureSRV[SampleCoord + uint2(1,0)] == 0 ? curDensity :  */DensityFieldTextureNewSRV[SampleCoord + uint2(1,0)]
				+ /* VelocityFieldStateTextureSRV[SampleCoord - uint2(0,1)] == 0 ? curDensity :  */DensityFieldTextureNewSRV[SampleCoord - uint2(0,1)]
				+ /* VelocityFieldStateTextureSRV[SampleCoord + uint2(0,1)] == 0 ? curDensity :  */DensityFieldTextureNewSRV[SampleCoord + uint2(0,1)]
			) * a
		) / (1 + 4 * a);

}



#define SCALE_CURL 0


[numthreads( 8, 8, 1 )]
void ComputeCurl( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const float samplingPeriodLength = CGConstantsCB.VelocityFieldGridCellLength;

	if(VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		VectorFieldCurlTextureUAV[SampleCoord] = 0.f;
		return;
	}

#if 1	//Curl 0 Near Borders
	bool bNearBorder = 0;

	bNearBorder |= VelocityFieldStateTextureSRV[SampleCoord + uint2( 1, 0 )] == 0;
	bNearBorder |= VelocityFieldStateTextureSRV[SampleCoord - uint2( 1, 0 )] == 0;
	bNearBorder |= VelocityFieldStateTextureSRV[SampleCoord + uint2( 0, 1 )] == 0;
	bNearBorder |= VelocityFieldStateTextureSRV[SampleCoord - uint2( 0, 1 )] == 0;

	if(bNearBorder)
	{
		VectorFieldCurlTextureUAV[SampleCoord] = 0.f;
		return;
	}
#endif

	float texelSize = 1.f / (CGConstantsCB.VelocityFieldGridSize);

	//map cur sample coord to tex coord
	float2 texCoordCenter = (float2(SampleCoord) + float(0.5f).xx) * texelSize.xx;

	//.v component from left right neighbours
	float vLeft = VelocityField_V_ComponentTextureSRV.SampleLevel(LinearSampler, texCoordCenter - float2(texelSize, 0.f), 0);
	float vRight = VelocityField_V_ComponentTextureSRV.SampleLevel(LinearSampler, texCoordCenter + float2(texelSize, 0.f), 0);

	//.u component from up down neighbours
	float uUp = VelocityField_U_ComponentTextureSRV.SampleLevel(LinearSampler, texCoordCenter - float2(0.f, texelSize), 0);
	float uDown = VelocityField_U_ComponentTextureSRV.SampleLevel(LinearSampler, texCoordCenter + float2(0.f, texelSize), 0);

	float VpartialX = vRight - vLeft;
	float UpartialY = uUp - uDown;

	float curl = VpartialX - UpartialY;

#if SCALE_CURL
	curl = curl / (2 * CGConstantsCB.VelocityFieldGridCellLength);
#endif

	VectorFieldCurlTextureUAV[SampleCoord] = (curl);

}

#define NORMALIZE_CURL_GRADIENT 0

//TODO: Write to Force Field first, then apply force in separate pass

[numthreads( 8, 8, 1 )]
void ApplyVorticityConfinementForce( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	
	if(VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		return;
	}

	const float h = CGConstantsCB.VelocityFieldGridCellLength;

	const float3 vorticityAxis = float3(0,0,1);
	
#if 0
	float curlRight = VectorFieldCurlTextureSRV[SampleCoord];
	float curlUp = VectorFieldCurlTextureSRV[SampleCoord];
	float curlLeft = VectorFieldCurlTextureSRV[SampleCoord - uint2(1,0)];
	float curlDown = VectorFieldCurlTextureSRV[SampleCoord + uint2(0,1)];
	
	float curlCur = (curlRight + curlLeft + curlUp + curlDown) * 0.25f;
#else
	float curlRight = VectorFieldCurlTextureSRV[SampleCoord + uint2(1,0)];
	float curlLeft = VectorFieldCurlTextureSRV[SampleCoord - uint2(1,0)];
	float curlDown = VectorFieldCurlTextureSRV[SampleCoord + uint2(0,1)];
	float curlUp = VectorFieldCurlTextureSRV[SampleCoord - uint2(0,1)];

	float curlCur = VectorFieldCurlTextureSRV[SampleCoord];
#endif


	//1. compute gradient of curl vector magnitude
	
	curlRight = abs(curlRight);
	curlUp = abs(curlUp);
	curlLeft = abs(curlLeft);
	curlDown = abs(curlDown);

	float2 gradient;
	gradient.x = (curlRight - curlLeft);
	gradient.y = (curlUp - curlDown);

	CurlGradientTextureUAV[SampleCoord] = length(gradient);
	//float3 curlGradientDirection = normalize(float3(0,0,1) + float3(gradient.xy, 0));
	float3 curlGradientDirection = normalize(float3(0,1,0) + float3(gradient.x, 0, gradient.y));
	CurlGradientDirectionTextureUAV[SampleCoord] = curlGradientDirection;

	//2. normalize the gradient
	if(any(gradient < EPSILON))
	{
		return;
	}

#if NORMALIZE_CURL_GRADIENT
	gradient = normalize(gradient);
#else
	gradient = log(gradient + 1);//produces more chaotic swirls
#endif

	//3. final force = OurConstant (>0) * CellLength * (Normalized Gradient (cross product) Curl Vector) 
#if 1
	float2 vorticityForce;
	float3 N = float3(gradient.x, gradient.y, 0);
	float3 w = vorticityAxis * curlCur;
	vorticityForce = (100 * CGConstantsCB.VorticityConfinementFactor * h * cross(N, w)).xy;
#else
	//cross product with (0,0,1)
	float2 vorticityForce;
	vorticityForce.x = gradient.y;
	vorticityForce.y = -gradient.x;

	vorticityForce *= curlCur * CGConstantsCB.VorticityConfinementFactor * h;
#endif
	//apply the force to velocities
	uint leftState = VelocityFieldStateTextureSRV[SampleCoord - uint2(1,0)];
	uint downState = VelocityFieldStateTextureSRV[SampleCoord + uint2(0,1)];

	if(leftState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		VelocityField_U_ComponentTextureUAV[SampleCoord] = VelocityField_U_ComponentTextureUAV[SampleCoord] + vorticityForce.x * CGConstantsCB.DeltaTime;
	}

	if(downState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		VelocityField_V_ComponentTextureUAV[SampleCoord] = VelocityField_V_ComponentTextureUAV[SampleCoord] + vorticityForce.y * CGConstantsCB.DeltaTime;
	}

}

[numthreads( 8, 8, 1 )]
void ApplyPressureFieldOffsetToClothParticles( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	
	if(any(SampleCoord >= CGConstantsCB.ClothSize))
	{
		return;
	}

	float2 sampleCoordf = SampleCoord;
	sampleCoordf += float2(0.5f,0.5f);

	//convert sample coord to UV
	float2 texCoord = sampleCoordf / CGConstantsCB.ClothSize;

	texCoord.y = 1.f - texCoord.y;

	float curPressure = PressureFieldTextureSRV.SampleLevel(LinearSampler, texCoord, 0);

	//offset particle .y based on cur pressure

	float offset = WATER_SURFACE_SIMULATION_HEIGHT + curPressure * CGConstantsCB.PressureControlledOffsetScale;

	float3 curPos = ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord)];

	curPos.y = offset;
	//curPos.y = 10.f;

	ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord)] = curPos;

}

[numthreads( 64, 1, 1 )]
void AdvectParticlesWithVelocityField( 
	uint SampleCoord : SV_DispatchThreadID )
{
	
	if((SampleCoord.x >= CGConstantsCB.NumParticles))
	{
		return;
	}

	float3 posCur = ParticlesBufferPositionCurUAV[SampleCoord.x];

	/* float2 velocityGridCoordf = MapToRange(posCur.xy, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, float(CGConstantsCB.VelocityFieldGridSize));
	uint2 velocityGridCoord = uint2(uint(velocityGridCoordf.x), uint(velocityGridCoordf.y));
	float2 curExternalVelocity2D =  VelocityFieldBufferSRV[GetFlatIndexToVelocityArray2D(velocityGridCoord)]; */

	float2 velocity = SampleVelocityField(posCur.xy);

	posCur.xy = posCur.xy + velocity * CGConstantsCB.DeltaTime/100.f;

	ParticlesBufferPositionCurUAV[SampleCoord.x] = posCur;

}


void GetVelocityFieldCoordFromWorldPos(float2 worldPos, out uint2 outCurCoord, out uint2 outCurCoordNext, out float2 outFrac)
{
	//most of the calculations here are spent on boundary checks
	const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset,CGConstantsCB.VelocityFieldBoundsStartOffset);
	const uint n = CGConstantsCB.VelocityFieldGridSize;
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const float2 minLimitPos = boundsStartPos + float2( h, h );
	const float2 maxLimitPos = boundsStartPos + float2( h * (n), h * (n) );
	worldPos.x = clamp(worldPos.x, minLimitPos.x, maxLimitPos.x);
	worldPos.y = clamp(worldPos.y, minLimitPos.y, maxLimitPos.y);
	float2 curCoordf = (worldPos - boundsStartPos) / float(CGConstantsCB.VelocityFieldGridCellLength);
	outCurCoord = uint2(uint(floor(curCoordf.x)), uint(floor(curCoordf.y)));
	outCurCoord.x = min(outCurCoord.x, n - 1);
	outCurCoord.y = min(outCurCoord.y, n - 1);
	outCurCoordNext = uint2( min(outCurCoord.x + 1, n - 1), min(outCurCoord.y + 1, n - 1) );
	outFrac = float2( curCoordf.x - outCurCoord.x, curCoordf.y - outCurCoord.y );
};




[numthreads( 8, 8, 1 )]
void TransferDataFromParticlesToVelocityField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	//accumulate scaled velocities, accumulate weights per velocity, accumulate weights per cell center 

	uint2 coordOffset = DecodeUint(VectorProjectionPerPassConstantsCB.Offset);
	const uint2 coordStride = uint2(2,2);
	SampleCoord = coordOffset + SampleCoord * coordStride;
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (CGConstantsCB.VelocityFieldGridSize - 1) - SampleCoord.y);

	if(VelocityFieldStateTextureUAV[SampleCoordInv] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		return;
	}

	const uint n = CGConstantsCB.VelocityFieldGridSize;
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const float hHalf = h * 0.5f;

	bool bCellContainsParticles = false;

	//U & V velocity coord is equal to cur cell coord
	const uint2 curCellCoordU = SampleCoord;
	const uint2 curCellCoordV = SampleCoord;

	//scalar-iterate over all particles
	for(int i = 0; i < CGConstantsCB.NumParticles; i++)
	{
#if PARTICLE_ADVECTION_DEBUG
		const float3 curParticlePos = CGConstantsCB.DefaultParticlePos;
		const float3 curParticleVel = PARTICLE_ADVECTION_DEBUG_VELOCITY;
#else
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
		const float3 curParticleVel = CurParticlesVelocityBufferSRV[i];
#endif//PARTICLE_ADVECTION_DEBUG


		//TODO: Create local scopes to minimize register pressure


		//compute particle's U & V cell's coord
		//if cur cell U or V coord alignes with a particle -> inject

		//U
		{
			uint2 curParticleCoordU;
			uint2 curParticleCoordUNext;
			float2 curParticlefracU;
			const float2 curParticlePosUCellSpace = curParticlePos.xy - float2(0, hHalf);
			GetVelocityFieldCoordFromWorldPos(curParticlePosUCellSpace, curParticleCoordU, curParticleCoordUNext, curParticlefracU);

			if(all(curParticleCoordU == curCellCoordU))
			{
				//inverse so that we are in inv texture space
				curParticleCoordU.y = (n-1) - curParticleCoordU.y;
				curParticleCoordUNext.y = (n-1) - curParticleCoordUNext.y;

				const float dCur = (1.f - curParticlefracU.x) * (1.f - curParticlefracU.y);
				const float dUp = (1.f - curParticlefracU.x) * curParticlefracU.y;
				const float dRight = curParticlefracU.x * (1.f - curParticlefracU.y);
				const float dRightUp = curParticlefracU.x * curParticlefracU.y;

				const uint2 curCoord = curParticleCoordU;
				const uint2 curCoordUp = uint2(curParticleCoordU.x, curParticleCoordUNext.y);
				const uint2 curCoordRight = uint2(curParticleCoordUNext.x, curParticleCoordU.y);
				const uint2 curCoordRightUp = uint2(curParticleCoordUNext.x, curParticleCoordUNext.y);

				//accumulate velocities
				VelocityField_U_ComponentTextureUAV[curCoord] = VelocityField_U_ComponentTextureUAV[curCoord] + curParticleVel.x * dCur;
				VelocityField_U_ComponentTextureUAV[curCoordUp] = VelocityField_U_ComponentTextureUAV[curCoordUp] + curParticleVel.x * dUp;
				VelocityField_U_ComponentTextureUAV[curCoordRight] = VelocityField_U_ComponentTextureUAV[curCoordRight] + curParticleVel.x * dRight;
				VelocityField_U_ComponentTextureUAV[curCoordRightUp] = VelocityField_U_ComponentTextureUAV[curCoordRightUp] + curParticleVel.x * dRightUp;

				//accumulate weights
				VelocityField_U_ComponentWeightTextureUAV[curCoord] = VelocityField_U_ComponentWeightTextureUAV[curCoord] + dCur;
				VelocityField_U_ComponentWeightTextureUAV[curCoordUp] = VelocityField_U_ComponentWeightTextureUAV[curCoordUp] + dUp;
				VelocityField_U_ComponentWeightTextureUAV[curCoordRight] = VelocityField_U_ComponentWeightTextureUAV[curCoordRight] + dRight;
				VelocityField_U_ComponentWeightTextureUAV[curCoordRightUp] = VelocityField_U_ComponentWeightTextureUAV[curCoordRightUp] + dRightUp;
			}
		}
		
		//V
		{
			uint2 curParticleCoordV;
			uint2 curParticleCoordVNext;
			float2 curParticlefracV;
			const float2 curParticlePosVCellSpace = curParticlePos.xy - float2(hHalf, 0);
			GetVelocityFieldCoordFromWorldPos(curParticlePosVCellSpace, curParticleCoordV, curParticleCoordVNext, curParticlefracV);

			if(all(curParticleCoordV == curCellCoordV))
			{
				float dCur = (1.f - curParticlefracV.x) * (1.f - curParticlefracV.y);
				float dUp = (1.f - curParticlefracV.x) * curParticlefracV.y;
				float dRight = curParticlefracV.x * (1.f - curParticlefracV.y);
				float dRightUp = curParticlefracV.x * curParticlefracV.y;

				//inverse so that we are in inv texture space
				curParticleCoordV.y = (n-1) - curParticleCoordV.y;
				curParticleCoordVNext.y = (n-1) - curParticleCoordVNext.y;

				const uint2 curCoord = curParticleCoordV;
				const uint2 curCoordUp = uint2(curParticleCoordV.x, curParticleCoordVNext.y);
				const uint2 curCoordRight = uint2(curParticleCoordVNext.x, curParticleCoordV.y);
				const uint2 curCoordRightUp = uint2(curParticleCoordVNext.x, curParticleCoordVNext.y);

				//accumulate velocities
				VelocityField_V_ComponentTextureUAV[curCoord] = VelocityField_V_ComponentTextureUAV[curCoord] + curParticleVel.y * dCur;
				VelocityField_V_ComponentTextureUAV[curCoordUp] = VelocityField_V_ComponentTextureUAV[curCoordUp] + curParticleVel.y * dUp;
				VelocityField_V_ComponentTextureUAV[curCoordRight] = VelocityField_V_ComponentTextureUAV[curCoordRight] + curParticleVel.y * dRight;
				VelocityField_V_ComponentTextureUAV[curCoordRightUp] = VelocityField_V_ComponentTextureUAV[curCoordRightUp] + curParticleVel.y * dRightUp;

				//accumulate weights
				VelocityField_V_ComponentWeightTextureUAV[curCoord] = VelocityField_V_ComponentWeightTextureUAV[curCoord] + dCur;
				VelocityField_V_ComponentWeightTextureUAV[curCoordUp] = VelocityField_V_ComponentWeightTextureUAV[curCoordUp] + dUp;
				VelocityField_V_ComponentWeightTextureUAV[curCoordRight] = VelocityField_V_ComponentWeightTextureUAV[curCoordRight] + dRight;
				VelocityField_V_ComponentWeightTextureUAV[curCoordRightUp] = VelocityField_V_ComponentWeightTextureUAV[curCoordRightUp] + dRightUp;
			}
		}
		
		
		//cell center
		{
			uint2 curParticleCoordCenter;
			uint2 curParticleCoordCenterNext;
			float2 curParticlefracCenter;
			const float2 curParticlePosCenterCellSpace = curParticlePos.xy - float2(hHalf, hHalf);
			GetVelocityFieldCoordFromWorldPos(curParticlePosCenterCellSpace, curParticleCoordCenter, curParticleCoordCenterNext, curParticlefracCenter);

			//inject density into cell center
			if(all(curParticleCoordCenter == SampleCoord))
			{
				float dCur = (1.f - curParticlefracCenter.x) * (1.f - curParticlefracCenter.y);
				float dUp = (1.f - curParticlefracCenter.x) * curParticlefracCenter.y;
				float dRight = curParticlefracCenter.x * (1.f - curParticlefracCenter.y);
				float dRightUp = curParticlefracCenter.x * curParticlefracCenter.y;

				//inverse so that we are in inv texture space
				curParticleCoordCenter.y = (n-1) - curParticleCoordCenter.y;
				curParticleCoordCenterNext.y = (n-1) - curParticleCoordCenterNext.y;

				const uint2 curCoord = curParticleCoordCenter;
				const uint2 curCoordUp = uint2(curParticleCoordCenter.x, curParticleCoordCenterNext.y);
				const uint2 curCoordRight = uint2(curParticleCoordCenterNext.x, curParticleCoordCenter.y);
				const uint2 curCoordRightUp = uint2(curParticleCoordCenterNext.x, curParticleCoordCenterNext.y);

				AdvectedParticleDensityTextureUAV[curCoord] = AdvectedParticleDensityTextureUAV[curCoord] + dCur;
				AdvectedParticleDensityTextureUAV[curCoordUp] = AdvectedParticleDensityTextureUAV[curCoordUp] + dUp;
				AdvectedParticleDensityTextureUAV[curCoordRight] = AdvectedParticleDensityTextureUAV[curCoordRight] + dRight;
				AdvectedParticleDensityTextureUAV[curCoordRightUp] = AdvectedParticleDensityTextureUAV[curCoordRightUp] + dRightUp;
			}
		}

		
		//mark nearest cell as density
		uint2 curParticleCoordNearest;
		uint2 dummy;
		float2 dumy;
		GetVelocityFieldCoordFromWorldPos(curParticlePos.xy, curParticleCoordNearest, dummy, dumy);
		if(all(curParticleCoordNearest == SampleCoord))
		{
			bCellContainsParticles = true;
		}

	}

	if(bCellContainsParticles)
	{
		VelocityFieldStateTextureUAV[SampleCoordInv] = VELOCITY_FIELD_CELL_STATE_DENSITY;
	}
	else
	{
		VelocityFieldStateTextureUAV[SampleCoordInv] = VELOCITY_FIELD_CELL_STATE_AIR;
	}

}


struct Interpolators
{
	float4 Position : SV_Position;
	nointerpolation float CurParticleVel : VELOCITY;
	nointerpolation float CurParticleWeight : WEIGHT;
	nointerpolation bool bIsDisabledParticle : PARTICLE_STATE;
#if INJECT_VELOCITY == 0
	nointerpolation uint CellState : STATE;
	nointerpolation uint2 CellCoord : CELLCOORD;
#endif//INJECT_VELOCITY == 0
	//nointerpolation float2 TexCoords : TEXCOORD;
};

void TransferDataFromParticlesToVelocityFieldVS(
	in uint InstanceID : SV_InstanceID,
	out Interpolators vsout
	)
{
	const uint n = CGConstantsCB.VelocityFieldGridSize;
	const float h = CGConstantsCB.VelocityFieldGridCellLength;
	const float hHalf = h * 0.5f;

	const uint curParticleIndex = InstanceID / 4;
	const uint curParticleLocalIndex = InstanceID % 4;

	//get cur particle pos & vel

#if PARTICLE_ADVECTION_DEBUG
		const float3 curParticlePos = CGConstantsCB.DefaultParticlePos;
		const float3 curParticleVel = PARTICLE_ADVECTION_DEBUG_VELOCITY;
#else
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[curParticleIndex];
		const float3 curParticleVel = CurParticlesVelocityBufferSRV[curParticleIndex];
#endif//PARTICLE_ADVECTION_DEBUG

	uint particleStateBitmask = ParticlesStateBufferSRV[InstanceID];
	vsout.bIsDisabledParticle = (particleStateBitmask & BITMASK_PARTICLE_STATE_DISABLED) != 0;

	float2 coordf = float(CGConstantsCB.VelocityFieldGridSize) / 2.f;

	//1.cur
	//2.up
	//3.right
	//4.rightUp
	float d[4];
	uint2 coords[4];

#if OUTPUT_U //0
	uint2 curParticleCoordU;
	uint2 curParticleCoordUNext;
	float2 curParticlefracU;
	const float2 curParticlePosUCellSpace = curParticlePos.xy - float2(0, hHalf);
	GetVelocityFieldCoordFromWorldPos(curParticlePosUCellSpace, curParticleCoordU, curParticleCoordUNext, curParticlefracU);

	d[0] = (1.f - curParticlefracU.x) * (1.f - curParticlefracU.y);
	d[1] = (1.f - curParticlefracU.x) * curParticlefracU.y;
	d[2] = curParticlefracU.x * (1.f - curParticlefracU.y);
	d[3] = curParticlefracU.x * curParticlefracU.y;

	coords[0] = curParticleCoordU;
	coords[1] = uint2(curParticleCoordU.x, curParticleCoordUNext.y);
	coords[2] = uint2(curParticleCoordUNext.x, curParticleCoordU.y);
	coords[3] = uint2(curParticleCoordUNext.x, curParticleCoordUNext.y);

	const float curD = d[curParticleLocalIndex];
	const uint2 curCoord = coords[curParticleLocalIndex];

	//accumulate velocities
	vsout.CurParticleVel = curParticleVel.x * curD;
	vsout.CurParticleWeight = curD;

	coordf = curCoord;
	
#elif OUTPUT_V //1

	uint2 curParticleCoordV;
	uint2 curParticleCoordVNext;
	float2 curParticlefracV;
	const float2 curParticlePosVCellSpace = curParticlePos.xy - float2(hHalf, 0);
	GetVelocityFieldCoordFromWorldPos(curParticlePosVCellSpace, curParticleCoordV, curParticleCoordVNext, curParticlefracV);

	d[0] = (1.f - curParticlefracV.x) * (1.f - curParticlefracV.y);
	d[1] = (1.f - curParticlefracV.x) * curParticlefracV.y;
	d[2] = curParticlefracV.x * (1.f - curParticlefracV.y);
	d[3] = curParticlefracV.x * curParticlefracV.y;

	coords[0] = curParticleCoordV;
	coords[1] = uint2(curParticleCoordV.x, curParticleCoordVNext.y);
	coords[2] = uint2(curParticleCoordVNext.x, curParticleCoordV.y);
	coords[3] = uint2(curParticleCoordVNext.x, curParticleCoordVNext.y);

	const float curD = d[curParticleLocalIndex];
	const uint2 curCoord = coords[curParticleLocalIndex];

	//accumulate velocities
	vsout.CurParticleVel = curParticleVel.y * curD;
	vsout.CurParticleWeight = curD;

	coordf = curCoord;

#else //OUTPUT_CENTER //2

	uint2 curParticleCoordCenter;
	uint2 curParticleCoordCenterNext;
	float2 curParticlefracCenter;
	const float2 curParticlePosCenterCellSpace = curParticlePos.xy - float2(hHalf, hHalf);
	GetVelocityFieldCoordFromWorldPos(curParticlePosCenterCellSpace, curParticleCoordCenter, curParticleCoordCenterNext, curParticlefracCenter);

	d[0] = (1.f - curParticlefracCenter.x) * (1.f - curParticlefracCenter.y);
	d[1] = (1.f - curParticlefracCenter.x) * curParticlefracCenter.y;
	d[2] = curParticlefracCenter.x * (1.f - curParticlefracCenter.y);
	d[3] = curParticlefracCenter.x * curParticlefracCenter.y;

	coords[0] = curParticleCoordCenter;
	coords[1] = uint2(curParticleCoordCenter.x, curParticleCoordCenterNext.y);
	coords[2] = uint2(curParticleCoordCenterNext.x, curParticleCoordCenter.y);
	coords[3] = uint2(curParticleCoordCenterNext.x, curParticleCoordCenterNext.y);

	const float curD = d[curParticleLocalIndex];
	uint2 curCoord = coords[curParticleLocalIndex];

	vsout.CurParticleWeight = curD;

	coordf = curCoord;

	#if INJECT_VELOCITY == 0
	#if ALL_PARTICLES_MARK_STATE == 0
		//mark nearest cell as density
		uint2 curParticleCoordNearest;
		uint2 dummy;
		float2 dumy;
		GetVelocityFieldCoordFromWorldPos(curParticlePos.xy, curParticleCoordNearest, dummy, dumy);
		vsout.CellState = curParticleLocalIndex == 0 ? VELOCITY_FIELD_CELL_STATE_DENSITY : VELOCITY_FIELD_CELL_STATE_AIR;
		curParticleCoordNearest.y = (n-1) - curParticleCoordNearest.y;
		vsout.CellCoord = curParticleCoordNearest;
	#else
		vsout.CellState = VELOCITY_FIELD_CELL_STATE_DENSITY;
		curCoord.y = (n-1) - curCoord.y;
		vsout.CellCoord = curCoord;
	#endif//ALL_PARTICLES_MARK_STATE
		
		
	#endif//INJECT_VELOCITY == 0

#endif//OUTPUT

	coordf += float2(0.5,0.5);
	//map to [-1, 1] NDC
	coordf = coordf / float(CGConstantsCB.VelocityFieldGridSize);
	coordf = coordf * 2.f - 1.f;

	vsout.Position = float4(coordf.x, coordf.y, 0, 1);

}


void TransferDataFromParticlesToVelocityFieldPS
	(
	in Interpolators psin,
#if INJECT_VELOCITY
	out float outVelocity : SV_Target0,
    out float outWeight : SV_Target1
#else
	out float outWeight : SV_Target0
#endif//INJECT_VELOCITY
	)
{
	if(psin.bIsDisabledParticle)
	{
		discard;
	}
#if INJECT_VELOCITY
	outVelocity = psin.CurParticleVel;
#else
	//we inject cell state
	if(psin.CellState == VELOCITY_FIELD_CELL_STATE_DENSITY)
	{
		VelocityFieldStateTextureUAV[psin.CellCoord] = psin.CellState;
	}
#endif//INJECT_VELOCITY
	outWeight = psin.CurParticleWeight;
}






[numthreads( 8, 8, 1 )]
void InvScaleEachVelocityByItsWeight( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	bool bsolid = VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID;

	//remove density from solid cells
	if(bsolid)
	{
		AdvectedParticleDensityTextureUAV[SampleCoord] = 0.f;
	}

	//Reset if cur cell solid or cell to the left is solid
	if (bsolid || (SampleCoord.x > 0 && VelocityFieldStateTextureUAV[(SampleCoord - uint2(1,0))] == VELOCITY_FIELD_CELL_STATE_SOLID))
	{
		VelocityField_U_ComponentTextureUAV[SampleCoord] = 0.f;
	}
	else
	{
		//float weightU = VelocityField_U_ComponentWeightTextureSRV[SampleCoord];
		float weightU = VelocityField_U_ComponentWeightTextureUAV[SampleCoord];
		if(weightU > 0.f)
		{
			VelocityField_U_ComponentTextureUAV[SampleCoord] = VelocityField_U_ComponentTextureUAV[SampleCoord] / weightU;
		}
	}
	//Reset if cur cell solid or downward cell is solid
	if (bsolid || (SampleCoord.y < (CGConstantsCB.VelocityFieldGridSize-1) && VelocityFieldStateTextureUAV[(SampleCoord + uint2(0, 1))] == VELOCITY_FIELD_CELL_STATE_SOLID))
	{
		VelocityField_V_ComponentTextureUAV[SampleCoord] = 0.f;
	}
	else
	{
		//float weightV = VelocityField_V_ComponentWeightTextureSRV[SampleCoord];
		float weightV = VelocityField_V_ComponentWeightTextureUAV[SampleCoord];
		if(weightV > 0.f)
		{
			VelocityField_V_ComponentTextureUAV[SampleCoord] = VelocityField_V_ComponentTextureUAV[SampleCoord] / weightV;
		}
	}
	
}


[numthreads( 1, 1, 1 )]
void ComputeAverageFieldDensity( 
	uint SampleCoord : SV_DispatchThreadID )
{
	//iterate over all field cells, gather densities and num density fields

	uint numDensityCells = 0;
	float densitySum = 0.f;

	for (uint y = 0; y < CGConstantsCB.VelocityFieldGridSize; y++)
	{
		for (uint x = 0; x < CGConstantsCB.VelocityFieldGridSize; x++)
		{
			uint2 coord = uint2(x, y);
			if (VelocityFieldStateTextureUAV[coord] == VELOCITY_FIELD_CELL_STATE_DENSITY)
			{
				densitySum += AdvectedParticleDensityTextureUAV[coord];
				numDensityCells++;
			}
		}
	}

	if (numDensityCells > 0)
	{
		AverageDensityBufferUAV[0] = densitySum / numDensityCells;
	}

}

#define SAMPLE_PARTICLE_CELL_DENSITY_PRECISE 0

float SampleAdvectedParticleDensity(float2 worldPos)
{
    const float h = CGConstantsCB.VelocityFieldGridCellLength;
    const uint n = CGConstantsCB.VelocityFieldGridSize;
    const float2 boundsStartPos = float2(CGConstantsCB.VelocityFieldBoundsStartOffset, CGConstantsCB.VelocityFieldBoundsStartOffset);
    const float2 minLimitPos = boundsStartPos + float2(h, h);
    const float2 maxLimitPos = boundsStartPos + float2(h * n, h * n);

    worldPos = clamp(worldPos, minLimitPos, maxLimitPos);

    // Map world coords to UV coords
    float2 fieldCoordUV = MapToRange(worldPos, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd(), 0.f, 1.f);
    // Invert .y
    fieldCoordUV.y = 1.f - fieldCoordUV.y;
    return AdvectedParticleDensityTextureSRV.SampleLevel(LinearSampler, fieldCoordUV, 0);
}

[numthreads( 64, 1, 1 )]
void ComputeParticlesColors( 
	uint SampleCoord : SV_DispatchThreadID )
{
	if(SampleCoord < CGConstantsCB.NumParticles)
	{
		float s = 0.01f;

		float3 curParticleColor = ParticlesColorsBufferUAV[SampleCoord];

		//curParticleColor.x = clamp(curParticleColor.x - s, 0.f, 1.f);
		//curParticleColor.y = clamp(curParticleColor.y - s, 0.f, 1.f);
		//curParticleColor.z = clamp(curParticleColor.z + s, 0.f, 1.f);

		curParticleColor.x = clamp(curParticleColor.x + s, 0.f, 1.f);
		curParticleColor.y = clamp(curParticleColor.y + s, 0.f, 1.f);
		curParticleColor.z = clamp(curParticleColor.z - s, 0.5f, 1.f);

		float3 curParticlePos = ParticlesBufferPositionCurSRV[SampleCoord];
		
#if !SAMPLE_PARTICLE_CELL_DENSITY_PRECISE
		uint2 nearestCellCoord;
		uint2 dummy;
		float2 dumy;
		GetVelocityFieldCoordFromWorldPos(curParticlePos.xy, nearestCellCoord, dummy, dumy);

		nearestCellCoord.y = (CGConstantsCB.VelocityFieldGridSize-1) - nearestCellCoord.y;
#endif//SAMPLE_PARTICLE_CELL_DENSITY_PRECISE
		float d0 = AverageDensityBufferUAV[0];

		if (d0 > 0.0) 
		{
		#if SAMPLE_PARTICLE_CELL_DENSITY_PRECISE
			float curParticleDensity = SampleAdvectedParticleDensity(curParticlePos.xy);
		#else
			float curParticleDensity = AdvectedParticleDensityTextureUAV[nearestCellCoord];
		#endif//SAMPLE_PARTICLE_CELL_DENSITY_PRECISE
			float relDensity = curParticleDensity / d0;
			if (relDensity < 0.7) 
			{
				const float sc = 0.8;
				/* curParticleColor.x = sc;
				curParticleColor.y = sc;
				curParticleColor.z = 1.0; */
				curParticleColor.x = 1.;
				curParticleColor.y = 1.;
				curParticleColor.z = 0.9;
			}
		}

		ParticlesColorsBufferUAV[SampleCoord] = curParticleColor;

	}
}

[numthreads( 64, 1, 1 )]
void TransferDataFromVelocityFieldToParticles( 
	uint SampleCoord : SV_DispatchThreadID )
{

	uint n = CGConstantsCB.VelocityFieldGridSize;
	float h = CGConstantsCB.VelocityFieldGridCellLength;
	float hHalf = h * 0.5f;

	if(SampleCoord < CGConstantsCB.NumParticles)
	{
		float3 particleCurPos = ParticlesBufferPositionCurSRV[SampleCoord];
		float3 particleCurVel = CurParticlesVelocityBufferUAV[SampleCoord];

		//compute U & V coords for cell particle is related to
		//U
		uint2 curParticleCoordU;
		uint2 curParticleCoordUNext;
		float2 curParticlefracU;
		const float2 curParticlePosUCellSpace = particleCurPos.xy - float2(0, hHalf);
		GetVelocityFieldCoordFromWorldPos(curParticlePosUCellSpace, curParticleCoordU, curParticleCoordUNext, curParticlefracU);

		//V
		uint2 curParticleCoordV;
		uint2 curParticleCoordVNext;
		float2 curParticlefracV;
		const float2 curParticlePosVCellSpace = particleCurPos.xy - float2(hHalf, 0);
		GetVelocityFieldCoordFromWorldPos(curParticlePosVCellSpace, curParticleCoordV, curParticleCoordVNext, curParticlefracV);

		//transfer U
		{
			curParticleCoordU.y = (n-1) - curParticleCoordU.y;
			curParticleCoordUNext.y = (n-1) - curParticleCoordUNext.y;

			float dCur = (1.f - curParticlefracU.x) * (1.f - curParticlefracU.y);
			float dUp = (1.f - curParticlefracU.x) * curParticlefracU.y;
			float dRight = curParticlefracU.x * (1.f - curParticlefracU.y);
			float dRightUp = curParticlefracU.x * curParticlefracU.y;

			const uint2 curCoord = curParticleCoordU;
			const uint2 curCoordUp = uint2(curParticleCoordU.x, curParticleCoordUNext.y);
			const uint2 curCoordRight = uint2(curParticleCoordUNext.x, curParticleCoordU.y);
			const uint2 curCoordRightUp = uint2(curParticleCoordUNext.x, curParticleCoordUNext.y);

			//get velocity
			const float vCur = VelocityField_U_ComponentTextureSRV[curCoord];
			const float vUp = VelocityField_U_ComponentTextureSRV[curCoordUp];
			const float vRight = VelocityField_U_ComponentTextureSRV[curCoordRight];
			const float vRightUp = VelocityField_U_ComponentTextureSRV[curCoordRightUp];

			//get unprojected velocity
			const float vCurPrev = VelocityField_U_ComponentTextureNewSRV[curCoord];
			const float vUpPrev = VelocityField_U_ComponentTextureNewSRV[curCoordUp];
			const float vRightPrev = VelocityField_U_ComponentTextureNewSRV[curCoordRight];
			const float vRightUpPrev = VelocityField_U_ComponentTextureNewSRV[curCoordRightUp];

			uint2 offset = uint2(1, 0);
			int stateCur = VelocityFieldStateTextureSRV[curCoord] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoord - offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateUp = VelocityFieldStateTextureSRV[curCoordUp] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordUp - offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateRight = VelocityFieldStateTextureSRV[curCoordRight] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordRight - offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateRightUp = VelocityFieldStateTextureSRV[curCoordRightUp] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordRightUp - offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;

			float d = stateCur * dCur + stateUp * dUp + stateRight * dRight + stateRightUp * dRightUp;

			if (d > 0.0) 
			{
				float picV = (stateCur * dCur * vCur + stateUp * dUp * vUp + stateRight * dRight * vRight + stateRightUp * dRightUp * vRightUp) / d;
				float corr = (stateCur * dCur * (vCur - vCurPrev) + stateUp * dUp * (vUp - vUpPrev) + stateRight * dRight * (vRight - vRightPrev) + stateRightUp * dRightUp * (vRightUp - vRightUpPrev)) / d;
				float flipV = particleCurVel.x + corr;
				particleCurVel.x = (1.0f - CGConstantsCB.ParticleAdvectionPICFLIPRatio) * picV + CGConstantsCB.ParticleAdvectionPICFLIPRatio * flipV;
			}
		}

		//transfer V
		{
			//inverse so that we are in inv texture space
			curParticleCoordV.y = (n-1) - curParticleCoordV.y;
			curParticleCoordVNext.y = (n-1) - curParticleCoordVNext.y;

			float dCur = (1.f - curParticlefracV.x) * (1.f - curParticlefracV.y);
			float dUp = (1.f - curParticlefracV.x) * curParticlefracV.y;
			float dRight = curParticlefracV.x * (1.f - curParticlefracV.y);
			float dRightUp = curParticlefracV.x * curParticlefracV.y;

			const uint2 curCoord = curParticleCoordV;
			const uint2 curCoordUp = uint2(curParticleCoordV.x, curParticleCoordVNext.y);
			const uint2 curCoordRight = uint2(curParticleCoordVNext.x, curParticleCoordV.y);
			const uint2 curCoordRightUp = uint2(curParticleCoordVNext.x, curParticleCoordVNext.y);

			//get velocity
			const float vCur = VelocityField_V_ComponentTextureSRV[curCoord];
			const float vUp = VelocityField_V_ComponentTextureSRV[curCoordUp];
			const float vRight = VelocityField_V_ComponentTextureSRV[curCoordRight];
			const float vRightUp = VelocityField_V_ComponentTextureSRV[curCoordRightUp];

			//get unprojected velocity
			const float vCurPrev = VelocityField_V_ComponentTextureNewSRV[curCoord];
			const float vUpPrev = VelocityField_V_ComponentTextureNewSRV[curCoordUp];
			const float vRightPrev = VelocityField_V_ComponentTextureNewSRV[curCoordRight];
			const float vRightUpPrev = VelocityField_V_ComponentTextureNewSRV[curCoordRightUp];

			uint2 offset = uint2(0, 1);
			int stateCur = VelocityFieldStateTextureSRV[curCoord] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoord + offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateUp = VelocityFieldStateTextureSRV[curCoordUp] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordUp + offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateRight = VelocityFieldStateTextureSRV[curCoordRight] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordRight + offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;
			int stateRightUp = VelocityFieldStateTextureSRV[curCoordRightUp] != VELOCITY_FIELD_CELL_STATE_AIR || VelocityFieldStateTextureSRV[curCoordRightUp + offset] != VELOCITY_FIELD_CELL_STATE_AIR ? 1 : 0;

			float d = stateCur * dCur + stateUp * dUp + stateRight * dRight + stateRightUp * dRightUp;

			if (d > 0.0) 
			{
				float picV = (stateCur * dCur * vCur + stateUp * dUp * vUp + stateRight * dRight * vRight + stateRightUp * dRightUp * vRightUp) / d;
				float corr = (stateCur * dCur * (vCur - vCurPrev) + stateUp * dUp * (vUp - vUpPrev) + stateRight * dRight * (vRight - vRightPrev) + stateRightUp * dRightUp * (vRightUp - vRightUpPrev)) / d;
				float flipV = particleCurVel.y + corr;
				particleCurVel.y = (1.0f - CGConstantsCB.ParticleAdvectionPICFLIPRatio) * picV + CGConstantsCB.ParticleAdvectionPICFLIPRatio * flipV;
			}
		}

		CurParticlesVelocityBufferUAV[SampleCoord] = particleCurVel;


	}
}

[numthreads( 8, 8, 1 )]
void DensityApplyVorticityGradientMask( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	if (VelocityFieldStateTextureSRV[SampleCoord] == VELOCITY_FIELD_CELL_STATE_SOLID)
	{
		DensityFieldTextureUAV[SampleCoord] = 0;
		return;
	}

	const float maskScale = 0.5f;

	const DENSITY_TYPE curDensity = DensityFieldTextureSRV[SampleCoord];

	float gradient = CurlGradientTextureSRV[SampleCoord];

#if 0 //map
	float2 gradientMinMax = CurlGradientMinMaxTextureSRV[uint2(0,0)];
	if(gradientMinMax.y != 0)
	{
		gradient = MapToRange(gradient, gradientMinMax.x, gradientMinMax.y, 0, 1);
	}
#endif

#if 1 //clamp
	const float gradientClampMin = 0.0f;
	const float gradientClampMax = 1.f; 
	gradient = clamp(gradient, gradientClampMin, gradientClampMax);
#endif


#if 0//filter
	const uint filterNumIterations = 2;

	for(int i = 0; i < filterNumIterations; i++)
	{
		gradient = 1 - sqrt(1 - gradient * gradient);//more blacks
		//gradient = sqrt(1 - pow(1 - gradient, 2));//more whites
	}
#endif
	
	//add
	//DENSITY_TYPE res = curDensity + gradient * maskScale;
	//scale
	//DENSITY_TYPE res = curDensity * gradient;
	//add scale
	DENSITY_TYPE res = curDensity + curDensity * gradient;

	DensityFieldTextureUAV[SampleCoord] = res;

}



[numthreads( 8, 8, 1 )]
void SolveWaveEquation( 
	uint2 SampleCoord : SV_DispatchThreadID )
{	

	uint particleStateBitmask = ParticlesStateBufferSRV[GetFlatIndexToClothArray(SampleCoord)];
	if((particleStateBitmask & BITMASK_PARTICLE_STATE_PINNED) != 0)
	{
		return;
	}

	//borders excluded
	if(any(SampleCoord >= (CGConstantsCB.ClothSize-2)))
	{
		return;
	}
	SampleCoord += uint2(1,1);
	
	float3 curPos = 	  	ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord + uint2(0,0))];

	float hRight = 	ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord + uint2(1,0))].y;
	float hLeft = ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord - uint2(1,0))].y;
	float hUp = 	  	ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord + uint2(0,1))].y;
	float hDown = ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord - uint2(0,1))].y;

	float k = pow(CGConstantsCB.WaveEquationWaveSpeed / CGConstantsCB.SpringRestLength, 2);

	float a = k * (hRight + hLeft + hUp + hDown - 4 * curPos.y);

	curPos.y += a * CGConstantsCB.DeltaTime;

	ParticlesBufferPositionCurUAV[GetFlatIndexToClothArray(SampleCoord + uint2(0,0))] = curPos;

}

[numthreads( 8, 8, 1 )]
void VisualizerConvertTextureWithNegativeValuesToColor( 
	uint2 SampleCoord : SV_DispatchThreadID )
{	

	float value = InputTextureWithNegativeValuesSRV[SampleCoord];

	float2 finalColor;

	if(value > 0)
	{
		finalColor.x = value;
	}
	else
	{
		finalColor.y = abs(value);
	}
	
	VisualizerColorTextureUAV[SampleCoord] = finalColor;

}




[numthreads( 64, 1, 1 )]
void DisableParticles(
	uint SampleIndex : SV_DispatchThreadID)
{
	if(SampleIndex < CGConstantsCB.NumParticles)
	{
		uint curParticleState = ParticlesStateBufferUAV[SampleIndex];

		if((curParticleState & BITMASK_PARTICLE_STATE_FALLING) != 0)
		{
			return;
		}

		//particles that intersect those planes are marked disabled

		//do distance check on .y axis

		const float3 curParticlePos = ParticlesBufferPositionCurSRV[SampleIndex];
		const float particleRadius = CGConstantsCB.PointRadius;

		/* const float killPlaneHeight1 = 0.75f * (GetSimulationBoundsEnd());
		const float killPlaneHeight2 = 0.25f * (GetSimulationBoundsEnd()); */

		const float killPlaneHeight1 = 0.75f * (GetSimulationBoundsEnd());
		const float killPlaneHeight2 = 0.287f * (GetSimulationBoundsEnd());

		if(abs(curParticlePos.y - killPlaneHeight1) < particleRadius)
		{
			curParticleState = curParticleState | BITMASK_PARTICLE_STATE_FALLING;
			ParticlesStateBufferUAV[SampleIndex] = curParticleState;

			//ParticlesColorsBufferUAV[SampleIndex] = float3(1,0,0);

		}
		else if(abs(curParticlePos.y - killPlaneHeight2) < particleRadius)
		{
			curParticleState = curParticleState | BITMASK_PARTICLE_STATE_FALLING;
			ParticlesStateBufferUAV[SampleIndex] = curParticleState;

			//ParticlesColorsBufferUAV[SampleIndex] = float3(1,0,0);

		}
	}
}



//===============================================
//
//				SDF Render
//
//===============================================

ConstantBuffer<SDFRenderConstants> CGSDFConstantsCB : register(b1);

ConstantBuffer<PaintDiffusionPerPassConstants> PaintDiffusionPerPassConstantsCB : register(b2);

ConstantBuffer<PaintRenderConstants> CGPaintRenderConstantsCB : register(b10);

//UAV
RWTexture2D<float4> SDFOutputTextureUAV : register(u21);
RWTexture2D<float> SDFTextureUAV : register(u1);
RWTexture2D<float4> SDFOutputNormalsTextureUAV : register(u2);
RWTexture2D<float> PaintMaxWetnessTextureUAV : register(u3);
RWTexture2D<float> WetnessParameterTextureUAV : register(u4);

RWTexture2D<float3> OriginalPaintingTextureUAV : register(u5);

RWTexture2D<float2> ParticleSourceTextureUAV : register(u6);

//SRV
Texture2D<float> SDFTextureSRV : register(t3);
StructuredBuffer<float3> ParticlesColorsBufferSRV : register(t4);
Texture2D<float4> SDFOutputTextureSRV : register(t5);
Texture2D<float4> SDFOutputTextureNewSRV : register(t6);
Texture2D<float> WetnessParameterTextureSRV : register(t8);

Texture2D<float3> GBufferColorTextureSRV : register (t10);
Texture2D<float4> GBufferNormalTextureSRV : register (t11);
Texture2D<float> GBufferRoughnessTextureSRV : register (t12);
Texture2D<float> GBufferWorldSpaceZSRV : register (t13);
Texture2D<uint> GBufferStencilSRV : register (t14);

StructuredBuffer<float3> ParticleRenderVertexBufferSRV : register (t15);
StructuredBuffer<float2> ParticleTexCoordsBufferSRV : register (t16);

Texture2D<float2> ParticleSourceTextureSRV : register (t17);

StructuredBuffer<float> AverageDensityBufferSRV : register (t21);

Texture2D<float> GBufferMetalTextureSRV : register (t22);



//Functions
float SphereSDF(float2 pointPos, float2 spherePos, float sphereradius)
{
	return length(pointPos - spherePos) - sphereradius;
}

float SphereSDF(float3 pointPos, float3 spherePos, float sphereradius)
{
	return length(pointPos - spherePos) - sphereradius;
}

float SurfaceSDF(float3 pointPos, float surfacePosZ)
{
	return surfacePosZ - pointPos.z;
}


float MaxSmooth(float a, float b, float k)
{
	return log(exp(k*a) + exp(k * b)) / k;
}
float MinSmooth(float a, float b, float k)
{
	return -MaxSmooth(-a, -b, k);
}

// exponential smooth min (k=32)
//produces the same result regardless of the order of the operations (min(min...))
/* float smin_e( float a, float b, float k )
{
    float res = exp2( -k*a ) + exp2( -k*b );
    return -log2( res )/k;
} */
float smin_e( float a, float b, float k )
{
    float res = exp( -k*a ) + exp( -k*b );
    return -log( res )/k;
}

float smin_e_colorT( in float a, in float b, in float k )
{
    return exp2( -k*b );
}

// root smooth min (k=0.01)
float smin_r( float a, float b, float k )
{
    float h = a-b;
    return 0.5*( (a+b) - sqrt(h*h+k) );
}

// power smooth min (k=8)
float smin_power( float a, float b, float k )
{
    a = pow( a, k ); 
	b = pow( b, k );
    return pow( (a*b)/(a+b), 1.0/k );
}



// polynomial smooth min 1 (k=0.1)
float smin_polynomial( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return lerp( b, a, h ) - k*h*(1.0-h);
}

// polynomial smooth min 2, faster (k=0.1)
float smin_polynomial2( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*k*(1.0/4.0);
}

// polynomial smooth min
float sminCubic( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}


float sminColor( float a, float b, float k )
{
    float h =  max( k-abs(a-b), 0.0 )/k;
    float m = h*h*0.5;
    return (a<b) ? m : 1.0-m;
}

float sminCubicColor( float a, float b, float k )
{
	float h =  max( k-abs(a-b), 0.0 )/k;
    float m = h*h*h*0.5;
    return (a<b) ? m : 1.0-m;
}

//The generalization to any power n being:
float2 sminN( float a, float b, float k, float n )
{
    float h =  max( k-abs(a-b), 0.0 )/k;
    float m = pow(h, n)*0.5;
    float s = m*k/n; 
    return (a<b) ? float2(a-s,m) : float2(b-s,1.0-m);
}

float GetSmoothMinSDF(float sdfMin, float sdfMin2, float blendRadius)
{
	if(CGSDFConstantsCB.SDFShapeBlendMode == SDF_BLEND_EXP)
	{
		return smin_e(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFShapeBlendMode == SDF_BLEND_ROOT)
	{
		return smin_r(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFShapeBlendMode == SDF_BLEND_POWER)
	{
		return smin_power(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFShapeBlendMode == SDF_BLEND_POLY)
	{
		return smin_polynomial2(sdfMin, sdfMin2, blendRadius);
	}
	else//SDF_BLEND_CUBIC
	{
		return sminCubic(sdfMin, sdfMin2, blendRadius);
	}
}

float GetSmoothMinNormals(float sdfMin, float sdfMin2, float blendRadius)
{
	if(CGSDFConstantsCB.SDFNormalsBlendMode == SDF_BLEND_EXP)
	{
		return smin_e(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFNormalsBlendMode == SDF_BLEND_ROOT)
	{
		return smin_r(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFNormalsBlendMode == SDF_BLEND_POWER)
	{
		return smin_power(sdfMin, sdfMin2, blendRadius);
	}
	else if(CGSDFConstantsCB.SDFNormalsBlendMode == SDF_BLEND_POLY)
	{
		return smin_polynomial2(sdfMin, sdfMin2, blendRadius);
	}
	else//SDF_BLEND_CUBIC
	{
		return sminCubic(sdfMin, sdfMin2, blendRadius);
	}
}

float GetColorBlendingWeight(float sdfMin, float sdfMin2)
{
	return sminColor(sdfMin, sdfMin2, CGSDFConstantsCB.ColorBlendRadius);
	return sminCubicColor(sdfMin, sdfMin2, CGSDFConstantsCB.ColorBlendRadius);
	return smin_e_colorT(sdfMin, sdfMin2, CGSDFConstantsCB.ColorBlendRadius);
}

#define ENABLE_BLENDING 1

float GetColorIntensity(float3 color)
{
    // Convert the color to grayscale by taking the average of its RGB components
    float grayscale = (color.r + color.g + color.b) / 3.0;

    // Normalize the grayscale intensity to the range [0, 1]
    float intensity = saturate(grayscale);

    return intensity;
}

#define SURFACE_PAINT_DIFFUSION 1

#define SDF_USE_RAYMARCH 1

#define SDF_NUM_PARTICLES_TO_BLEND 5

#define SDF_WRITE_INV_HEIGHT 1


[numthreads( 8, 8, 1 )]
void GenerateSignedDistanceField( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	//ortho raymarch until hit
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (CGSDFConstantsCB.SDFOutputSize.y - 1) - SampleCoord.y);
	float2 sampleCoordf = float2(SampleCoord) + float2(0.5,0.5);

	//get worldPos of the cell

	//map to 0-1
	float2 fieldCoordUV = MapToRange(sampleCoordf, 0.f, float(CGSDFConstantsCB.SDFOutputSize.x), 0.f, 1.f);
	//map to particle simulation bounds
	float2 cellCenterWorldPos = MapToRange(fieldCoordUV, 0.f, 1.f, GetSimulationBoundsStart(), GetSimulationBoundsEnd());

	const float3 rayOrigin = float3(cellCenterWorldPos.xy, CGConstantsCB.DefaultParticlePos.z - CGConstantsCB.PointRadius + CGSDFConstantsCB.RayOriginDepthOffset);
	const float3 rayDirection = float3(0,0,1);

	float sdfMin = 10000.f;
	float sdfMin2 = 10000.f;
	float sdfMinForNormals;

	int minParticleIndex = -1;
	int minParticleIndex2 = -1;

#if SDF_USE_RAYMARCH == 0

	//iterate over all particles 
	for(int i = 0; i < CGConstantsCB.NumParticles; i++)
	{
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
		const float particleRadius = CGConstantsCB.PointRadius;
		float sdfCur = SphereSDF(rayOrigin, curParticlePos, particleRadius);
#if ENABLE_BLENDING == 0
		sdfMin = min(sdfCur, sdfMin);
#else
		if (sdfCur < sdfMin)
		{
			sdfMin2 = sdfMin;
			sdfMin = sdfCur;

			minParticleIndex2 = minParticleIndex;
			minParticleIndex = i;
		}
		else if (sdfCur < sdfMin2)
		{
		    sdfMin2 = sdfCur;
			minParticleIndex2 = i;
		}
#endif//ENABLE_BLENDING
	}

#if ENABLE_BLENDING

	//get final blended value

	float colorLerpT;
	{
		colorLerpT = GetColorBlendingWeight(sdfMin, sdfMin2);

		sdfMinForNormals = GetSmoothMinNormals(sdfMin, sdfMin2, CGSDFConstantsCB.NormalsBlendRadius);

		sdfMin = GetSmoothMinSDF(sdfMin, sdfMin2, CGSDFConstantsCB.BlendRadius);
	}

#endif//ENABLE_BLENDING

#else //RAYMARCH

	float3 sdfMinColor = float3(1,1,1);

	float3 curPointOnRay = rayOrigin;

	//find 2 closest particles to iterate on
	{
		//put sample position from the up surface to to the middle
		curPointOnRay.z = rayOrigin.z + CGConstantsCB.PointRadius / 2.f;

		for(int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			if((ParticlesStateBufferSRV[i] & BITMASK_PARTICLE_STATE_FALLING) != 0)
			{
				continue;
			}
			/* if((ParticlesStateBufferSRV[i] & BITMASK_PARTICLE_STATE_FALLING) == 0)
			{
				continue;
			} */

			const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
			const float particleRadius = CGConstantsCB.PointRadius;
			float sdfCur = SphereSDF(curPointOnRay, curParticlePos, particleRadius);

			if (sdfCur < sdfMin)
			{
				sdfMin2 = sdfMin;
				sdfMin = sdfCur;

				minParticleIndex2 = minParticleIndex;
				minParticleIndex = i;
			}
			else if (sdfCur < sdfMin2)
			{
			    sdfMin2 = sdfCur;
				minParticleIndex2 = i;
			}
		}
	}

	if(isinf(sdfMin))
	{
		return;
	}

#if PAINT_RENDER_ADVECT_SOURCE == 0
	//blend colors
	{
		const float3 color1 = ParticlesColorsBufferSRV[minParticleIndex];
		sdfMinColor = color1;
		
		if(minParticleIndex2 > 0)//TODO: Don't run this branch if 2nd particle is too far away
		{
			float colorLerpT = GetColorBlendingWeight(sdfMin, sdfMin2);
			const float3 color2 = ParticlesColorsBufferSRV[minParticleIndex2];
			sdfMinColor = lerp(sdfMinColor, color2, colorLerpT); 
		}
	}
#endif	

	curPointOnRay = rayOrigin;
	const float surfacePositionZ = CGConstantsCB.DefaultParticlePos.z + CGSDFConstantsCB.SurfaceDepthOffset;

	float sdfDepthForNormals = curPointOnRay.z;

	const int kMaxIterations = 100;
	//const int kMaxIterations = 1;
	for(int i = 0; i < kMaxIterations; i++)
	{
		//blend particles and surface
		{
			const float sdfCur1 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex], CGConstantsCB.PointRadius);
			sdfMin = sdfCur1;
		}
		
		if(minParticleIndex2 >= 0)
		{
			const float sdfCur2 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex2], CGConstantsCB.PointRadius);
			//sdfMin = smin_e(sdfMin, sdfCur2, CGSDFConstantsCB.BlendRadius);
			sdfMin = GetSmoothMinSDF(sdfMin, sdfCur2, CGSDFConstantsCB.BlendRadius);
		}

		if(sdfMin < CGSDFConstantsCB.Threshold)
		{
			//we've hit a surface
			break;
		}

		curPointOnRay.z = curPointOnRay.z + sdfMin;
	}

	curPointOnRay = rayOrigin;
	for(int j = 0; j < kMaxIterations; j++)
	{
		{
			const float sdfCur1 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex], CGConstantsCB.PointRadius);
			sdfMinForNormals = sdfCur1;
		}

		if(minParticleIndex2 >= 0)
		{
			const float sdfCur2 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex2], CGConstantsCB.PointRadius);
			sdfMinForNormals = GetSmoothMinNormals(sdfMinForNormals, sdfCur2, CGSDFConstantsCB.NormalsBlendRadius);
		}

		{
			const float sdfToSurface = SurfaceSDF(curPointOnRay, surfacePositionZ);
			sdfMinForNormals = GetSmoothMinNormals(sdfMinForNormals, sdfToSurface, CGSDFConstantsCB.NormalsBlendRadius);
			//sdfMin = smin_e(sdfMin, sdfToSurface, CGSDFConstantsCB.BlendRadius);
		}

		if(sdfMinForNormals < CGSDFConstantsCB.Threshold)
		{
			//we've hit a surface
			break;
		}

		curPointOnRay.z = curPointOnRay.z + sdfMinForNormals;
	} 
	
	#endif//SDF_USE_RAYMARCH


	//output SDF to compute normals
	{
	#if SDF_USE_RAYMARCH
		const float4 curColor = SDFOutputTextureUAV[SampleCoordInv];
		#if SDF_WRITE_INV_HEIGHT
		SDFTextureUAV[SampleCoordInv] = (surfacePositionZ - curPointOnRay.z) + (curColor.a * CGSDFConstantsCB.WetnessContributionToHeight);
		#else
		SDFTextureUAV[SampleCoordInv] = curPointOnRay.z + (curColor.a * CGSDFConstantsCB.WetnessContributionToHeight);
		#endif
		//SDFTextureUAV[SampleCoordInv] = curPointOnRay.z;
	#else
		SDFTextureUAV[SampleCoordInv] = sdfMinForNormals;
	#endif//SDF_USE_RAYMARCH
	}


#if SDF_USE_RAYMARCH
	const float thres = CGSDFConstantsCB.Threshold;
#else
	const float thres = max(CGConstantsCB.PointRadius/2.f, CGSDFConstantsCB.Threshold);
#endif
	if(sdfMin <= thres)
	{
#if PAINT_RENDER_ADVECT_SOURCE == 0
		float3 finalColor;
		
	#if SDF_USE_RAYMARCH
		finalColor = sdfMinColor;
	#else
		const float3 color1 = ParticlesColorsBufferSRV[minParticleIndex];
		finalColor = color1;
		
		if(minParticleIndex2 > 0)
		{
			const float3 color2 = ParticlesColorsBufferSRV[minParticleIndex2];
			finalColor = lerp(finalColor, color2, colorLerpT); 
		}

	#endif//SDF_USE_RAYMARCH

#if SURFACE_PAINT_DIFFUSION
		//if cur color is intense, mix with cur surface
		const float4 curColor = SDFOutputTextureUAV[SampleCoordInv + uint2(0,1)];
		if(curColor.a > CGSDFConstantsCB.ColorFadeThreshold)
		{
			finalColor = lerp(curColor.rgb, finalColor, 0.5f);
			//finalColor = curColor;
		}
#endif//SURFACE_PAINT_DIFFUSION

		float curTileMaxWetness = PaintMaxWetnessTextureUAV[SampleCoordInv];
		float4 curTileColor = SDFOutputTextureUAV[SampleCoordInv];

		if(curTileColor.a < 1.f)
		{
			curTileColor.a = 1.f;
		}
		else if(curTileColor.a < CGSDFConstantsCB.WetnessLimit)
		{
			curTileColor.a += CGSDFConstantsCB.WetnessAccumulationAmount * CGConstantsCB.DeltaTime;
			if(curTileColor.a > curTileMaxWetness)
			{
				PaintMaxWetnessTextureUAV[SampleCoordInv] = curTileColor.a;
			}
		}
		
		SDFOutputTextureUAV[SampleCoordInv] = float4(finalColor.xyz, curTileColor.a);

#else//ADVECTION
		
		//based on cur particle velocity, sample paint color from prev particle pos


		//float metal = GBufferMetalTextureSRV[SampleCoordInv];
		float2 curFieldCoordUV = MapToRange(cellCenterWorldPos, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
		curFieldCoordUV.y = 1.f - curFieldCoordUV.y;
		float metal = GBufferMetalTextureSRV.SampleLevel(LinearSampler, curFieldCoordUV, 0);

		if(metal < EPSILON)
		{
			//get world pos of cur tile
			float2 curParticleVelocity = CurParticlesVelocityBufferSRV[minParticleIndex];

			float2 prevCellCenterWorldPos = cellCenterWorldPos - curParticleVelocity * CGConstantsCB.DeltaTime;
			//map to [0,1]
			float2 prevFieldCoordUV = MapToRange(prevCellCenterWorldPos, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
			prevFieldCoordUV.y = 1.f - prevFieldCoordUV.y;

			float2 prevCoordf = prevFieldCoordUV * CGSDFConstantsCB.SDFOutputSize;
			uint2 prevCoord = uint2(uint(floor(prevCoordf.x)), uint(floor(prevCoordf.y)));
			//float4 prevTileColor = SDFOutputTextureSRV[prevCoord];
			float4 prevTileColor = SDFOutputTextureSRV.SampleLevel(LinearSampler, prevFieldCoordUV, 0);

			//write prev tile color here
			float3 curTileColor = OriginalPaintingTextureUAV[SampleCoordInv];
			//curTileColor.rgb = prevTileColor;
			curTileColor.rgb = lerp(prevTileColor.rgb, curTileColor, 0.25f);

			#if 0
			float3 curParticleColor1 = ParticlesColorsBufferUAV[minParticleIndex];
			curTileColor.rgb = lerp(curParticleColor1.rgb, curTileColor, 0.75f);
			#else
			ParticlesColorsBufferUAV[minParticleIndex] = curTileColor;
			#endif

			OriginalPaintingTextureUAV[SampleCoordInv] = curTileColor;
			
			//OriginalPaintingTextureUAV[SampleCoordInv] = float3(1,0,1);
			//ParticlesColorsBufferUAV[minParticleIndex] = float3(1,0,1);
		}
		else
		{
			OriginalPaintingTextureUAV[SampleCoordInv] = ParticlesColorsBufferUAV[minParticleIndex];

			//OriginalPaintingTextureUAV[SampleCoordInv] = float3(0,1,0);
			//ParticlesColorsBufferUAV[minParticleIndex] = float3(0,1,0);
		}


		

		

		float4 curPaintColor = SDFOutputTextureUAV[SampleCoordInv];
		if(curPaintColor.a < 1.f)
		{
			curPaintColor.a = 1.f;
		}
		else if(curPaintColor.a < CGSDFConstantsCB.WetnessLimit)
		{
			curPaintColor.a += CGSDFConstantsCB.WetnessAccumulationAmount * CGConstantsCB.DeltaTime;
			float curTileMaxWetness = PaintMaxWetnessTextureUAV[SampleCoordInv];
			if(curPaintColor.a > curTileMaxWetness)
			{
				PaintMaxWetnessTextureUAV[SampleCoordInv] = curPaintColor.a;
			}
		}
		SDFOutputTextureUAV[SampleCoordInv] = curPaintColor;

		
		

#endif//ADVECTION

	}
	else
	{
		//SDFOutputTextureUAV[SampleCoordInv] = float4(0,0,1,1);

		//simulate color drying out
		float4 curSurfaceColor =  SDFOutputTextureUAV[SampleCoordInv];
		float curTileMaxWetness = PaintMaxWetnessTextureUAV[SampleCoordInv];
		if(curSurfaceColor.a > curTileMaxWetness * CGSDFConstantsCB.ColorFadeThreshold)
		{
		#if 0
			curSurfaceColor.rgb *= CGSDFConstantsCB.ColorFadeScaleInv;
			curSurfaceColor.a *= CGSDFConstantsCB.ColorFadeScaleInv;
		#else
			curSurfaceColor.a -= CGSDFConstantsCB.WetnessDryingOutAmount * CGConstantsCB.DeltaTime * 0.1f;
			//curSurfaceColor.rgb *= min(1.f, curSurfaceColor.a);
		#endif
			SDFOutputTextureUAV[SampleCoordInv] = curSurfaceColor;
		}
	}
}





#if 0
float extent = (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) * 0.25;

		float2 simulationCenter = GetPaintSimulationCenterPos();
		float simulationlengthHalf = GetPaintSimulationLength() * 0.5;

		if(
			cellCenterWorldPos.x > (CGConstantsCB.VelocityFieldBoundsStartOffset + extent) 
			&& cellCenterWorldPos.x < (CGConstantsCB.VelocityFieldBoundsStartOffset + (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) - extent)
			&& cellCenterWorldPos.y > (CGConstantsCB.VelocityFieldBoundsStartOffset + extent)  
			&& cellCenterWorldPos.y < (CGConstantsCB.VelocityFieldBoundsStartOffset + (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) - extent)
		)
		{
			//based on cur particle velocity, sample paint color from prev particle pos

			//get world pos of cur tile
			float2 curParticleVelocity = CurParticlesVelocityBufferSRV[minParticleIndex];

			float2 prevCellCenterWorldPos = cellCenterWorldPos - curParticleVelocity * CGConstantsCB.DeltaTime;
			//map to [0,1]
			float2 prevFieldCoordUV;
			//float2 prevFieldCoordUV = MapToRange(prevCellCenterWorldPos, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
			prevFieldCoordUV.x = MapToRange(prevCellCenterWorldPos.x, (simulationCenter.x - extent), simulationCenter.x + extent, 0.f, 1.f);
			prevFieldCoordUV.y = MapToRange(prevCellCenterWorldPos.y, (simulationCenter.y - extent), simulationCenter.y + extent, 0.f, 1.f);
			prevFieldCoordUV.y = 1.f - prevFieldCoordUV.y;

			//float2 prevCoordf = prevFieldCoordUV * CGSDFConstantsCB.SDFOutputSize;
			//uint2 prevCoord = uint2(uint(floor(prevCoordf.x)), uint(floor(prevCoordf.y)));
			//float4 prevTileColor = SDFOutputTextureSRV[prevCoord];

			float4 prevTileColor = SDFOutputTextureSRV.SampleLevel(LinearSampler, prevFieldCoordUV, 0);
			//float4 prevTileColor = SDFOutputTextureSRV[SampleCoordInv];

			//write prev tile color here
			float4 curTileColor = SDFOutputTextureForAdvectionUAV[SampleCoordInv];
			curTileColor.rgb = prevTileColor;
			//curTileColor.rgb = lerp(prevTileColor, curTileColor, 0.75f);


			
			float4 curTilWetnessColor = SDFOutputTextureUAV[SampleCoordInv];
			if(curTilWetnessColor.a < 1.f)
			{
				curTilWetnessColor.a = 1.f;
			}
			else if(curTilWetnessColor.a < CGSDFConstantsCB.WetnessLimit)
			{
				float curTileMaxWetness = PaintMaxWetnessTextureUAV[SampleCoordInv];
				curTilWetnessColor.a += CGSDFConstantsCB.WetnessAccumulationAmount * CGConstantsCB.DeltaTime;
				if(curTilWetnessColor.a > curTileMaxWetness)
				{
					PaintMaxWetnessTextureUAV[SampleCoordInv] = curTilWetnessColor.a;
				}
			}

			curTilWetnessColor.rgb = curTileColor.rgb;
			//SDFOutputTextureUAV[SampleCoordInv] = curTilWetnessColor;

			ParticlesColorsBufferUAV[minParticleIndex] = curTileColor;

			//SDFOutputTextureUAV[SampleCoordInv] = curTileColor;

		#if 1
			//float2 curTileUV = MapToRange(cellCenterWorldPos, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
			float2 curTileUV;
			curTileUV.x = MapToRange(cellCenterWorldPos.x, (simulationCenter.x - extent), simulationCenter.x + extent, 0.f, 1.f);
			curTileUV.y = MapToRange(cellCenterWorldPos.y, (simulationCenter.y - extent), simulationCenter.y + extent, 0.f, 1.f);
			//float2 curTileUV = MapToRange(cellCenterWorldPos, (CGConstantsCB.VelocityFieldBoundsStartOffset + extent) , (CGConstantsCB.VelocityFieldBoundsStartOffset + (CGConstantsCB.VelocityFieldGridSize * CGConstantsCB.VelocityFieldGridCellLength) - extent), 0.f, 1.f);
			curTileUV.y = 1.f - curTileUV.y; 

			//prevFieldCoordUV *= 0.5f;

			//curTileUV = MapToRange(curTileUV, 0.f, 1.f, 0.25f, 0.75f);	
			uint2 scaledSampleCoord = floor(curTileUV * float2(CGSDFConstantsCB.SDFOutputSize)); 	
			SDFOutputTextureForAdvectionUAV[scaledSampleCoord] = curTileColor;
		#else
			SDFOutputTextureForAdvectionUAV[SampleCoordInv] = curTileColor;
			//SDFOutputTextureForAdvectionUAV[SampleCoordInv] = float4(1,1,0, 1);
		#endif

#endif


[numthreads( 8, 8, 1 )]
void RenderFallingParticles( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	//ortho raymarch until hit
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (CGSDFConstantsCB.SDFOutputSize.y - 1) - SampleCoord.y);
	float2 sampleCoordf = float2(SampleCoord) + float2(0.5,0.5);

	//get worldPos of the cell

	//map to 0-1
	float2 fieldCoordUV = MapToRange(sampleCoordf, 0.f, float(CGSDFConstantsCB.SDFOutputSize.x), 0.f, 1.f);
	//map to particle simulation bounds
	float2 cellCenterWorldPos = MapToRange(fieldCoordUV, 0.f, 1.f, GetSimulationBoundsStart(), GetSimulationBoundsEnd());

	const float surfacePositionZ = CGConstantsCB.DefaultParticlePos.z;
	const float3 rayOrigin = float3(cellCenterWorldPos.xy, max(0.f, CGSDFConstantsCB.FalingParticleStartPosZ - CGConstantsCB.PointRadius));

	//const float3 rayDirection = float3(0,0,1);

	float sdfMin = 10000.f;
	float sdfMin2 = 10000.f;
	float sdfMinForNormals;

	int minParticleIndex = -1;
	int minParticleIndex2 = -1;

	float3 sdfMinColor = float3(1,1,1);

	float3 curPointOnRay = rayOrigin;

	#if 1
	
	//find 2 closest particles to iterate on
	{
		//put sample position from the up surface to to the middle
		curPointOnRay.z = rayOrigin.z + CGConstantsCB.PointRadius / 2.f;	
		for(int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			if((ParticlesStateBufferSRV[i] & BITMASK_PARTICLE_STATE_FALLING) == 0)
			{
				continue;
			}	
			const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
			const float particleRadius = CGConstantsCB.PointRadius;
			float sdfCur = SphereSDF(curPointOnRay, curParticlePos, particleRadius);	
			if (sdfCur < sdfMin)
			{
				sdfMin2 = sdfMin;
				sdfMin = sdfCur;	
				minParticleIndex2 = minParticleIndex;
				minParticleIndex = i;
			}
			else if (sdfCur < sdfMin2)
			{
			    sdfMin2 = sdfCur;
				minParticleIndex2 = i;
			}
		}
	}	

#if 1
	//first test
	const float kEarlyOutDistance = CGConstantsCB.DefaultParticlePos.z;
	if(sdfMin > kEarlyOutDistance)
	{
		SDFOutputTextureUAV[SampleCoordInv] = float4(0,0,0,0);
		//SDFOutputTextureUAV[SampleCoordInv] = float4(0,1,0,1);
		return;
	}
#endif//

	//blend colors
	{
		const float3 color1 = ParticlesColorsBufferSRV[minParticleIndex];
		sdfMinColor = color1;

		if(minParticleIndex2 > 0)//TODO: Don't run this branch if 2nd particle is too far away
		{
			float colorLerpT = GetColorBlendingWeight(sdfMin, sdfMin2);
			const float3 color2 = ParticlesColorsBufferSRV[minParticleIndex2];
			sdfMinColor = lerp(sdfMinColor, color2, colorLerpT); 
		}
	}

	curPointOnRay = rayOrigin;
	
	float minStepSize = 0.01f; 
	const int kMaxIterations = 16;
	for(int i = 0; i < kMaxIterations; i++)
	{
		sdfMin = 10000.f;

		for(int i = 1; i < CGConstantsCB.NumParticles; i++)
		{
			if((ParticlesStateBufferSRV[i] & BITMASK_PARTICLE_STATE_FALLING) == 0)
			{
				continue;
			}

			const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
			//float curParticleRadius = CGConstantsCB.PointRadius * pow(4.f - curParticlePos.z, 2);
			float curParticleRadius = CGConstantsCB.PointRadius * pow(4.f - curParticlePos.z, 2);

			const float sdfCur = SphereSDF(curPointOnRay, curParticlePos, curParticleRadius);
			sdfMin = GetSmoothMinSDF(sdfMin, sdfCur, CGSDFConstantsCB.FallingParticleBlendRadius);
		}

		//const float sdfToSurface = SurfaceSDF(curPointOnRay, surfacePositionZ);
		//sdfMin = GetSmoothMinSDF(sdfMin, sdfToSurface, CGSDFConstantsCB.FallingParticleBlendRadius);

		const bool bApproachedSurface = curPointOnRay.z >= surfacePositionZ;

		if(sdfMin < CGSDFConstantsCB.Threshold || bApproachedSurface)
		{
			//we've hit a surface
			break;
		}

		curPointOnRay.z = curPointOnRay.z + max(minStepSize, sdfMin);
		minStepSize += 0.02;

	}

	if(sdfMin < CGSDFConstantsCB.Threshold)
	{
		float curHeight = SDFTextureUAV[SampleCoordInv];
		#if SDF_WRITE_INV_HEIGHT
		SDFTextureUAV[SampleCoordInv] = curHeight + (surfacePositionZ - curPointOnRay.z);
		#else
		SDFTextureUAV[SampleCoordInv] = curHeight + curPointOnRay.z;
		#endif//SDF_WRITE_INV_HEIGHT
	}
	

	#else

	//find 2 closest particles to iterate on
	{
		//put sample position from the up surface to to the middle
		curPointOnRay.z = rayOrigin.z + CGConstantsCB.PointRadius;

		for(int i = 0; i < CGConstantsCB.NumParticles; i++)
		{
			if((ParticlesStateBufferSRV[i] & BITMASK_PARTICLE_STATE_FALLING) == 0)
			{
				continue;
			}

			const float3 curParticlePos = ParticlesBufferPositionCurSRV[i];
			const float particleRadius = CGConstantsCB.PointRadius;
			float sdfCur = SphereSDF(curPointOnRay, curParticlePos, particleRadius);

			if (sdfCur < sdfMin)
			{
				sdfMin2 = sdfMin;
				sdfMin = sdfCur;

				minParticleIndex2 = minParticleIndex;
				minParticleIndex = i;
			}
			else if (sdfCur < sdfMin2)
			{
			    sdfMin2 = sdfCur;
				minParticleIndex2 = i;
			}
		}
	}

	//blend colors
	{
		const float3 color1 = ParticlesColorsBufferSRV[minParticleIndex];
		sdfMinColor = color1;
		
		if(minParticleIndex2 > 0)//TODO: Don't run this branch if 2nd particle is too far away
		{
			float colorLerpT = GetColorBlendingWeight(sdfMin, sdfMin2);
			const float3 color2 = ParticlesColorsBufferSRV[minParticleIndex2];
			sdfMinColor = lerp(sdfMinColor, color2, colorLerpT); 
		}
	}	

	curPointOnRay = rayOrigin;
	//const float surfacePositionZ = CGConstantsCB.DefaultParticlePos.z + CGSDFConstantsCB.SurfaceDepthOffset;
	const float surfacePositionZ = CGConstantsCB.DefaultParticlePos.z + CGSDFConstantsCB.SurfaceDepthOffset + 0.4f;

	float sdfDepthForNormals = curPointOnRay.z;

	const int kMaxIterations = 100;
	for(int i = 0; i < kMaxIterations; i++)
	{
		//blend particles and surface
		{
			const float sdfCur1 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex], CGConstantsCB.PointRadius);
			sdfMin = sdfCur1;
		}
		
		if(minParticleIndex2 >= 0)
		{
			const float sdfCur2 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex2], CGConstantsCB.PointRadius);
			//sdfMin = smin_e(sdfMin, sdfCur2, CGSDFConstantsCB.BlendRadius);
			sdfMin = GetSmoothMinSDF(sdfMin, sdfCur2, CGSDFConstantsCB.FallingParticleBlendRadius);
		}

		if(sdfMin < CGSDFConstantsCB.Threshold)
		{
			//we've hit a surface
			break;
		}

		curPointOnRay.z = curPointOnRay.z + sdfMin;
	}

	curPointOnRay = rayOrigin;
	for(int j = 0; j < kMaxIterations; j++)
	{
		{
			const float sdfCur1 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex], CGConstantsCB.PointRadius);
			sdfMinForNormals = sdfCur1;
		}

		if(minParticleIndex2 >= 0)
		{
			const float sdfCur2 = SphereSDF(curPointOnRay, ParticlesBufferPositionCurSRV[minParticleIndex2], CGConstantsCB.PointRadius);
			sdfMinForNormals = GetSmoothMinNormals(sdfMinForNormals, sdfCur2, CGSDFConstantsCB.FallingParticleBlendRadius);
		}

		{
			const float sdfToSurface = SurfaceSDF(curPointOnRay, surfacePositionZ);
			sdfMinForNormals = GetSmoothMinNormals(sdfMinForNormals, sdfToSurface, CGSDFConstantsCB.FallingParticleBlendRadius);
			//sdfMin = smin_e(sdfMin, sdfToSurface, CGSDFConstantsCB.BlendRadius);
		}

		if(sdfMinForNormals < CGSDFConstantsCB.Threshold)
		{
			//we've hit a surface
			break;
		}

		curPointOnRay.z = curPointOnRay.z + sdfMinForNormals;
	} 
	
	float curHeight = SDFTextureUAV[SampleCoordInv];
	#if SDF_WRITE_INV_HEIGHT
	SDFTextureUAV[SampleCoordInv] = curHeight + (surfacePositionZ - curPointOnRay.z);
	#else
	SDFTextureUAV[SampleCoordInv] = curHeight + curPointOnRay.z;
	#endif//SDF_WRITE_INV_HEIGHT
	//
	
#endif

	if(sdfMin <= CGSDFConstantsCB.Threshold)
	{
		SDFOutputTextureUAV[SampleCoordInv] = float4(sdfMinColor.xyz, 1.f);
	}
	else
	{
		SDFOutputTextureUAV[SampleCoordInv] = float4(0,0,0,0);
	}
}





[numthreads( 8, 8, 1 )]
void GenerateParticleSourceTexture( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (CGSDFConstantsCB.ParticleSourceTextureSize - 1) - SampleCoord.y);
	
	float2 uv = float2(SampleCoord) / (CGSDFConstantsCB.ParticleSourceTextureSize-1).xx;
	float2 vec = MapToRange(uv, 0, 1, -1, 1);
	float l = length(vec);
	float res = max(0, 1 - l);

	/* float3 N;
	N.xy = uv*2.0-1.0;
	float r2 = dot(N.xy, N.xy);
	N.z = -sqrt(1.0 - r2); */

	//QuadraticThroughAGivenPoint()

	float resBumpy = sqrt(1 - pow(1 - res, 2));
	float resBumpy10 = sqrt(1 - pow(1 - res, 10));

	ParticleSourceTextureUAV[SampleCoordInv] = float2(sqrt(1 - pow(1 - res, CGSDFConstantsCB.ParticleBumpinnesPow)), resBumpy);


}


struct ParticleRenderInterpolators
{
	float4 Position : SV_Position;
	float3 PositionView : POSVIEW;
	nointerpolation float3 Color : COLOR;
	nointerpolation bool bIsFallingParticle : PARTICLE_STATE_FALLING;
	nointerpolation bool bIsDisabledParticle : PARTICLE_STATE_DISABLED;
	float2 TexCoords : TEXCOORD;
};



void ParticleRenderVS(
	in uint VertexIndex : SV_VertexID,
	in uint InstanceID : SV_InstanceID,
	out ParticleRenderInterpolators vsout
	)
{
	const float3 curParticlePos = ParticlesBufferPositionCurSRV[InstanceID];

	uint particleStateBitmask = ParticlesStateBufferSRV[InstanceID];
	vsout.bIsFallingParticle = (particleStateBitmask & BITMASK_PARTICLE_STATE_FALLING) != 0;

	vsout.bIsDisabledParticle = (particleStateBitmask & BITMASK_PARTICLE_STATE_DISABLED) != 0;

	const bool bIsNewlyGeneratedParticle = (particleStateBitmask & BITMASK_PARTICLE_STATE_NEW) != 0;

	/*
	*	Position
	*/
	float3 positionWorld = ParticleRenderVertexBufferSRV[VertexIndex];

	positionWorld.xy *= CGConstantsCB.PointRadius * CGConstantsCB.RasterizedParticleRadiusScale;

	#if 1//apply additional scale based on velocity
	if(!bIsNewlyGeneratedParticle)
	{
		const float3 curParticleVelocity = CurParticlesVelocityBufferSRV[InstanceID];
		float2 velNorm = normalize(abs(curParticleVelocity.xy));
		//positionWorld.x *= saturate(0.90 + velNorm.x);
		//positionWorld.y *= (1.f + velNorm.y * 1.5);
		//positionWorld.x *= (1.f + velNorm.x * 1.0);
		//positionWorld.y *= (1.f + velNorm.y * 2.0);
		//positionWorld.x *= (1.f + pow(velNorm.x,2));
		positionWorld.y *= (1.f + pow(velNorm.y,3) * 2.0f);

	}
	#endif

	#if 1//apply additional scale based on density
	if(vsout.bIsFallingParticle && !bIsNewlyGeneratedParticle)
	//if(!bIsNewlyGeneratedParticle)
	{
		const float curParticleDensity = SampleAdvectedParticleDensity(curParticlePos.xy);
		const float avrgDensity = AverageDensityBufferSRV[0];
		positionWorld.xy *= clamp(curParticleDensity, 0.5, 1.1);
		//positionWorld.xy *= clamp(curParticleDensity / avrgDensity, 0.5, 1.5);
		//float r = (curParticleDensity / avrgDensity);
		//vsout.Color = float3(r,r,r);
	}
	#endif
	vsout.Color = ParticlesColorsBufferSRV[InstanceID];
	
	//additionally scale XY based on distance
	if(vsout.bIsFallingParticle)
	{
		positionWorld.xy *= max(bIsNewlyGeneratedParticle ? 2.0f : 1.0f, CGConstantsCB.DefaultParticlePos.z - curParticlePos.z);
	}
	
	
	//translate
	positionWorld += curParticlePos;


	vsout.PositionView = mul(CGPaintRenderConstantsCB.PaintSceneViewMatrix, float4(positionWorld, 1.0f));
	//vsout.PositionWorld = positionWorld;
	vsout.Position = mul(CGPaintRenderConstantsCB.PaintSceneViewProjectionMatrix, float4(positionWorld, 1.0f));

	/*
	*	TexCoord
	*/
	vsout.TexCoords = ParticleTexCoordsBufferSRV[VertexIndex];
}


//throughraster
void StickyParticleRenderPS
	(
	in ParticleRenderInterpolators psin
	,out float outDepth : SV_Depth
	,out float4 outColor : SV_Target0
	,out float outHeightMap : SV_Target1
	)
{

	if(psin.bIsDisabledParticle)
	{
		discard;
	}

	if(psin.bIsFallingParticle)
	{
		discard;
	}

	float2 texture = ParticleSourceTextureSRV.SampleLevel(LinearSampler, psin.TexCoords, 0);

	float4 pixelPos;
	pixelPos.xyz = psin.PositionView;
	pixelPos.z = pixelPos.z - (texture.r) * CGConstantsCB.PointRadius;
	pixelPos = mul(CGPaintRenderConstantsCB.PaintSceneProjectionMatrix, float4(pixelPos.xyz, 1.0f));
	outDepth = pixelPos.z / pixelPos.w;

	if(texture.r > 0)
	{
		outColor = float4(psin.Color, texture.r);
		//outColor = float4(psin.Color, texture.g);
		//outColor = float4(psin.Color, 1.f);
		outHeightMap = texture.g;
	}
	else
	{
		discard;
	}
}

void FallingParticleRenderPS
	(
	in ParticleRenderInterpolators psin
	,out float outDepth : SV_Depth
	,out float4 outColor : SV_Target0
	,out float outHeightMap : SV_Target1
	)
{

	if(psin.bIsDisabledParticle)
	{
		discard;
	}

	if(!psin.bIsFallingParticle)
	{
		discard;
	}

	float2 texture = ParticleSourceTextureSRV.SampleLevel(LinearSampler, psin.TexCoords, 0);

	//convert view pos 
	float4 pixelPos;
	pixelPos.xyz = psin.PositionView;
	pixelPos.z = pixelPos.z - (texture.r) * CGConstantsCB.PointRadius;
	pixelPos = mul(CGPaintRenderConstantsCB.PaintSceneProjectionMatrix, float4(pixelPos.xyz, 1.0f));
	outDepth = pixelPos.z / pixelPos.w;

	if(texture.r > 0)
	{
		//outColor = float4(psin.Color, texture.r);
		outColor = float4(psin.Color, 1.f);
		outHeightMap = (CGConstantsCB.DefaultParticlePos.z - psin.PositionView.z) + texture.g;
	}
	else
	{
		discard;
	}
}

[numthreads( 8, 8, 1 )]
void ApplyWetnessToMarkedParticles( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (CGSDFConstantsCB.SDFOutputSize.y - 1) - SampleCoord.y);

	//Sample Color .a
	float4 curFramePaintTexture = SDFOutputTextureSRV[SampleCoordInv];
	float markedPixelValue = curFramePaintTexture.a;

	float4 prevFramePaintTexture = SDFOutputTextureUAV[SampleCoordInv];

	//float curWetness = WetnessParameterTextureUAV[SampleCoordInv];
	float curWetness = prevFramePaintTexture.a;

	float curTileMaxWetness = PaintMaxWetnessTextureUAV[SampleCoordInv];

	if(markedPixelValue > 0)
	{
		// if(curWetness < 1.f)
		// {
		// 	curWetness = 1.f;
		// }
		if(curWetness < markedPixelValue)
		{
			//curWetness = 1.f;
			curWetness = markedPixelValue;
		}
		else if(curWetness < CGSDFConstantsCB.WetnessLimit)
		{
			//curWetness += CGSDFConstantsCB.WetnessAccumulationAmount * CGConstantsCB.DeltaTime;
			curWetness += markedPixelValue * CGSDFConstantsCB.WetnessAccumulationAmount * CGConstantsCB.DeltaTime;
			if(curWetness > curTileMaxWetness)
			{
				PaintMaxWetnessTextureUAV[SampleCoordInv] = curWetness;
			}
		}
	}
	else
	{
		curFramePaintTexture.rgb = prevFramePaintTexture.rgb;

		//decrement
		if(curWetness > curTileMaxWetness * CGSDFConstantsCB.ColorFadeThreshold)
		{
			curWetness -= CGSDFConstantsCB.WetnessDryingOutAmount * CGConstantsCB.DeltaTime * 0.1f;
		}
	}

	//WetnessParameterTextureUAV[SampleCoordInv] = curWetness;

	//Mix paint
	#if 1
	{
		float4 prevFramePaintTexture = SDFOutputTextureUAV[SampleCoordInv];
		if(prevFramePaintTexture.a > CGSDFConstantsCB.ColorFadeThreshold)
		{
			float wetDiff = saturate(curFramePaintTexture.a - prevFramePaintTexture.a);
			wetDiff = 0.5f;
			curFramePaintTexture.rgb = lerp(prevFramePaintTexture.rgb, curFramePaintTexture.rgb, wetDiff);
			//finalColor = curColor;
		}
	}
	#endif

	SDFOutputTextureUAV[SampleCoordInv] = float4(curFramePaintTexture.rgb, curWetness);


	//add current wetness contribution to height map
	float curHeight = SDFTextureSRV[SampleCoordInv];
	#if 0
	if((curHeight - CGConstantsCB.PointRadius) <= curWetness * CGSDFConstantsCB.WetnessContributionToHeight)
	{
		SDFTextureUAV[SampleCoordInv] = curHeight + curWetness * CGSDFConstantsCB.WetnessContributionToHeight;
	}
	else
	{
		//falling particle
		SDFTextureUAV[SampleCoordInv] = curHeight;
	}
	#else
		SDFTextureUAV[SampleCoordInv] = curHeight + curWetness * CGSDFConstantsCB.WetnessContributionToHeight;
	#endif
	


}


[numthreads( 8, 8, 1 )]
void CopyWetnessParameter( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	#if COPY_FROM
	WetnessParameterTextureUAV[SampleCoord] = SDFOutputTextureSRV[SampleCoord].a;
	#else//COPY_INTO
	float4 c = SDFOutputTextureUAV[SampleCoord];
	c.a = WetnessParameterTextureSRV[SampleCoord];
	SDFOutputTextureUAV[SampleCoord] = c;
	#endif//COPY_FROM
}


#define PAINT_DIFFUSE_8_TAPS 0

[numthreads( 8, 8, 1 )]
void PaintDiffusion( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	uint2 coordOffset = (PaintDiffusionPerPassConstantsCB.LocalOffset);
	const uint2 coordStride = uint2(2,2);
	SampleCoord = coordOffset + SampleCoord * coordStride;

	//invert .y
	SampleCoord.y = (CGSDFConstantsCB.SDFOutputSize.y - 1) - SampleCoord.y;	

	float4 curColor = SDFOutputTextureSRV[SampleCoord];

	const float kThreshold = 0.1f;

	if(curColor.a < kThreshold)
	{
		return;
	}

	//Tile can only write to a neighbor who's Wetness is < than his

	//Sample nearby tiles, if any is dry -> propagate cur tile color to it
	float4 colorLeft = SDFOutputTextureSRV[SampleCoord - uint2(1,0)];
	float4 colorRight = SDFOutputTextureSRV[SampleCoord + uint2(1,0)];
	float4 colorUp = SDFOutputTextureSRV[SampleCoord - uint2(0,1)];
	float4 colorDown = SDFOutputTextureSRV[SampleCoord + uint2(0,1)];
#if PAINT_DIFFUSE_8_TAPS
	float4 colorLeftDown = SDFOutputTextureSRV[SampleCoord - uint2(1,0) + uint2(0,1)];
	float4 colorLeftUp = SDFOutputTextureSRV[SampleCoord - uint2(1,1)];
	float4 colorRightDown = SDFOutputTextureSRV[SampleCoord + uint2(1,1)];
	float4 colorRightUp = SDFOutputTextureSRV[SampleCoord + uint2(1,0) - uint2(0,1)];
#endif//PAINT_DIFFUSE_8_TAPS
	const float colorMixT = 0.5;

	curColor.a = saturate(curColor.a);

	if(colorLeft.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(1,0)] = lerp(colorLeft, curColor, colorMixT);
	}
	if(colorRight.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(1,0)] = lerp(colorRight, curColor, colorMixT);
	}
	if(colorUp.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(0,1)] = lerp(colorUp, curColor, colorMixT);
	}
	if(colorDown.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(0,1)] = lerp(colorDown, curColor, colorMixT);
	}

#if PAINT_DIFFUSE_8_TAPS
	if(colorLeftDown.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(1,0) + uint2(0,1)] = lerp(colorLeftDown, curColor, colorMixT);
	}
	if(colorLeftUp.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(1,1)] = lerp(colorLeftUp, curColor, colorMixT);
	}
	if(colorRightDown.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(1,1)] = lerp(colorRightDown, curColor, colorMixT);
	}
	if(colorRightUp.a < kThreshold)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(1,0) - uint2(0,1)] = lerp(colorRightUp, curColor, colorMixT);
	}
#endif//PAINT_DIFFUSE_8_TAPS

	/* const float kColorDampScale = 0.9999999f;

	if(colorLeft.a < 0.75)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(1,0)] = curColor * kColorDampScale;
	}
	if(colorRight.a < 0.75)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(1,0)] = curColor* kColorDampScale;
	}
	if(colorUp.a < 0.75)
	{
		SDFOutputTextureUAV[SampleCoord - uint2(0,1)] = curColor* kColorDampScale;
	}
	if(colorDown.a < 0.75)
	{
		SDFOutputTextureUAV[SampleCoord + uint2(0,1)] = curColor* kColorDampScale;
	} */

}

float4 SamplePaintTexture(float3 worldPos)
{
	float2 uvPos = MapToRange(worldPos.xy, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
	uvPos.y = 1.f - uvPos.y;

	/* float2 curCoordf = uvPos * CGSDFConstantsCB.SDFOutputSize;
	uint2 curCoord = uint2(uint(floor(curCoordf.x)), uint(floor(curCoordf.y))); */

	return SDFOutputTextureSRV.SampleLevel(LinearSampler, uvPos, 0);
}


#define PAINT_PARTICLES_DEBUG 0
#define PARTICLE_PAINT_DIFFUSION 1

[numthreads( 64, 1, 1 )]
void ModifyPaintParticles(
	uint SampleIndex : SV_DispatchThreadID)
{
	if(SampleIndex < CGConstantsCB.NumParticles)
	{
		
		uint curParticleState = ParticlesStateBufferSRV[SampleIndex];

		if((curParticleState & BITMASK_PARTICLE_STATE_DISABLED) != 0)
		{
			return;
		}
		if((curParticleState & BITMASK_PARTICLE_STATE_FALLING) != 0)
		{
			return;
		}


		float3 curParticleColor = ParticlesColorsBufferUAV[SampleIndex];
		const float3 curParticlePos = ParticlesBufferPositionCurSRV[SampleIndex];
		const float3 curParticleVelocity = CurParticlesVelocityBufferSRV[SampleIndex];
	
		//if particle is on the wet surface (and is surrounded by wet surfaces) -> boost it's velocity
		//boost depends on the wetness

		float4 paintColor = SamplePaintTexture(curParticlePos);

		//advect particle and check wetness in place where it's going

		const float kAdvectionCheckScale = 1.f;
		//const float kAdvectionCheckScale = 50.f;

		float velocityLength = length(curParticleVelocity);
		if(velocityLength < EPSILON)
		{
	#if PAINT_PARTICLES_DEBUG
			ParticlesColorsBufferUAV[SampleIndex] = float3(1,0,1);
	#endif//PAINT_PARTICLES_DEBUG
			return;
		}

		float3 advectedParticlePosition = curParticlePos + curParticleVelocity * CGConstantsCB.DeltaTime * kAdvectionCheckScale  + /* move outside the circle */ normalize(curParticleVelocity) * CGConstantsCB.PointRadius * 1.5;
		float4 nextPaintColor = SamplePaintTexture(advectedParticlePosition);

		//if surface is wet
		if(paintColor.a > CGSDFConstantsCB.DrySurfaceDampThreshold && nextPaintColor.a > CGSDFConstantsCB.DrySurfaceDampThreshold)
		{
	#if PAINT_PARTICLES_DEBUG
			ParticlesColorsBufferUAV[SampleIndex] = float3(1,0,0);
	#else

		#if PARTICLE_PAINT_DIFFUSION
			//mix curParticleColor with surface color
			//float3 colorAvrg = (curParticleColor + paintColor.rgb) * 0.5;
			float3 colorAvrg = (curParticleColor + paintColor.rgb + nextPaintColor.rgb) / 3.f;
			curParticleColor = curParticleColor + (colorAvrg - curParticleColor) * CGSDFConstantsCB.ColorDiffusionScale;
			ParticlesColorsBufferUAV[SampleIndex] = curParticleColor;

			//paintColor.rgb = paintColor.rgb + (colorAvrg - paintColor.rgb) * CGSDFConstantsCB.ColorDiffusionScale;
		#endif//PARTICLE_PAINT_DIFFUSION
	#endif//PAINT_PARTICLES_DEBUG	
		}
		else
		{
			//damp particles that are sliding down dry surface
	#if 1
			//TODO:Only apply upward force
			//CurParticlesExternalForceBufferUAV[SampleIndex] = (curParticleVelocity) * -CGSDFConstantsCB.PaintDryAreasVelocityDampForceScale;
			CurParticlesExternalForceBufferUAV[SampleIndex] = float3(curParticleVelocity.x * 0.25f, curParticleVelocity.y, 0) * -CGSDFConstantsCB.PaintDryAreasVelocityDampForceScale;
			//CurParticlesExternalForceBufferUAV[SampleIndex] = float3(curParticleVelocity.x, curParticleVelocity.y, 0) * -CGSDFConstantsCB.PaintDryAreasVelocityDampForceScale;
	#else //square method
			velocityLength = length(curParticleVelocity);
            float resistance = -CGSDFConstantsCB.PaintDryAreasVelocityDampForceScale * velocityLength * velocityLength;
            float3 resistanceWithDirection = (curParticleVelocity / velocityLength) * resistance;
            CurParticlesExternalForceBufferUAV[SampleIndex] = resistanceWithDirection;
	#endif

	#if PAINT_PARTICLES_DEBUG
			ParticlesColorsBufferUAV[SampleIndex] = float3(0,1,0);
	#endif//PAINT_PARTICLES_DEBUG
			
		}

	}

}

[numthreads( 64, 1, 1 )]
void ModifyPaintParticlesWithGBuffer(
	uint SampleIndex : SV_DispatchThreadID)
{
	if(SampleIndex < CGConstantsCB.NumParticles)
	{	
		
		uint curParticleState = ParticlesStateBufferUAV[SampleIndex];

		if((curParticleState & BITMASK_PARTICLE_STATE_DISABLED) != 0)
		{
			return;
		}

		float3 curParticlePos = ParticlesBufferPositionCurUAV[SampleIndex];

		//check if particle should be disabled
		//if(0)
		{
			if(curParticlePos.y <= GetSimulationBoundsStart() + EPSILON)
			{
				curParticleState = curParticleState | BITMASK_PARTICLE_STATE_DISABLED;
				ParticlesStateBufferUAV[SampleIndex] = curParticleState;
				return;
			}
		}

		float2 uvPos = MapToRange(curParticlePos.xy, GetSimulationBoundsStart(), GetSimulationBoundsEnd(), 0.f, 1.f);
		uvPos.y = 1.f - uvPos.y;

		//Sample Scene World Pos
		float curSceneZ = GBufferWorldSpaceZSRV.SampleLevel(LinearSampler, uvPos, 0);

		const float stickyParticleDetectionDistanceOffset = 0.75; //higher value more 
		const float fallingParticleDetectionDistanceOffset = EPSILON * 10; //offset surface forward and test intersection

		if((curParticleState & BITMASK_PARTICLE_STATE_FALLING) != 0)
		{
			//check if particle collides with a surface
			float particleOffset = CGConstantsCB.PointRadius;
			float surfaceOffset = stickyParticleDetectionDistanceOffset;

			const bool bIsNewParticle = (curParticleState & BITMASK_PARTICLE_STATE_NEW) != 0;

#if 0//DEBUG_FALLING_PARTICLES
			if(bIsNewParticle)
			{
				curParticlePos.z = CGSDFConstantsCB.FalingParticleStartPosZ;
				//curParticlePos.z = CGSDFConstantsCB.FalingParticleStartPosZ + float(SampleIndex) * 0.01;
				ParticlesBufferPositionCurUAV[SampleIndex] = curParticlePos;
				//return;
			}
#endif

			if(bIsNewParticle)
			{
				particleOffset = EPSILON * 10;
				surfaceOffset = 0.f;
			}

			//check for collisions
			if(((curParticlePos.z + particleOffset) > curSceneZ + surfaceOffset))
			{	
				//ParticlesColorsBufferUAV[SampleIndex] = float3(1,1,0);

				if(bIsNewParticle)
				{
					//if new particle landed on height surface -> mark it as sticky one
					uint2 coord = floor(float2(CGSDFConstantsCB.GBufferSize.xy) * uvPos);
					uint heightSurfaceBits = GBufferStencilSRV[coord];
					if(heightSurfaceBits > 0)
					{
						curParticleState = curParticleState | BITMASK_PARTICLE_STATE_STICKY;
					}


					//add velocity along normal that we've hit
					//if(1)
					{
						float3 curV = CurParticlesVelocityBufferUAV[SampleIndex];
						float curVStrength = length(curV);
						float3 curNormal = GBufferNormalTextureSRV.SampleLevel(PointSampler, uvPos, 0).xyz;
						if(abs(curNormal.z) > 0.99)
						{
							//curNormal.xy = (curNormal.xy) * 1.5f;
						}
						curNormal.z = 0.f;
						
						float clampValue = CGSDFConstantsCB.SurfaceHitVelocityMax;
						//curVStrangth = pow(curVStrangth, 2);
						curV = clamp(curNormal * CGSDFConstantsCB.SurfaceHitVelocityScale * curVStrength, -clampValue, clampValue);
						//curV = reflect(curV, curNormal) * 10.f * curVStrangth;
						CurParticlesVelocityBufferUAV[SampleIndex] = curV;
					}
					
				}

				curParticleState = curParticleState & (~BITMASK_PARTICLE_STATE_FALLING);
				curParticleState = curParticleState & (~BITMASK_PARTICLE_STATE_NEW);
				ParticlesStateBufferUAV[SampleIndex] = curParticleState;

				if(bIsNewParticle)
				{
					return;
				}
				//ParticlesColorsBufferUAV[SampleIndex] = float3(1,1,0);
			}
			else
			{
				
			#if 1//move falling particle towards the surface
				if((curParticleState & BITMASK_PARTICLE_STATE_NEW) != 0)
				{
					//curParticlePos.z = curParticlePos.z + CGConstantsCB.PointRadius * CGSDFConstantsCB.FallingParticleSpeedTowardsSurface * CGConstantsCB.DeltaTime;
					curParticlePos.z = curParticlePos.z + CGSDFConstantsCB.FallingParticleSpeedTowardsSurface * CGConstantsCB.DeltaTime;
					curParticlePos.z = min(curParticlePos.z, CGConstantsCB.DefaultParticlePos.z);

					//float t = curParticlePos.z / CGConstantsCB.DefaultParticlePos.z;
					//ParticlesColorsBufferUAV[SampleIndex] = lerp(float3(0,0,1), float3(1,0,0), t);
				}
				else
				{
					curParticlePos.z = curParticlePos.z + CGConstantsCB.PointRadius * 0.5 * CGConstantsCB.DeltaTime;

					//ParticlesColorsBufferUAV[SampleIndex] = float3(0,1,0);
				}
				ParticlesBufferPositionCurUAV[SampleIndex] = curParticlePos;
			#endif
				return;
			}
			
		}

#if RENDER_PARTICLES_SDF

#else
		/* if((curParticleState & BITMASK_PARTICLE_STATE_FALLING) == 0)
		{
			curParticleState = curParticleState | BITMASK_PARTICLE_STATE_DISABLED;
			ParticlesStateBufferUAV[SampleIndex] = curParticleState;
			return;
		} */

#endif

		const bool bIsStickyParticle = (curParticleState & BITMASK_PARTICLE_STATE_STICKY) != 0;

		if(bIsStickyParticle)
		{

			//when particle leaves the mesh territory-> mark it as falling one 
			uint2 coord = floor(float2(CGSDFConstantsCB.GBufferSize.xy) * uvPos);
			uint heightSurfaceBits = (GBufferStencilSRV[coord + uint2(0,2)] + GBufferStencilSRV[coord + uint2(0,3)]) * 0.5f;
			if(heightSurfaceBits < 1)
			{
				curParticleState = curParticleState | BITMASK_PARTICLE_STATE_FALLING;
				//curParticleState = curParticleState & (~BITMASK_PARTICLE_STATE_STICKY);
				ParticlesStateBufferUAV[SampleIndex] = curParticleState;

				//ParticlesColorsBufferUAV[SampleIndex] = float3(1,1,0);
				return;
			}
			
		}

#if 1//COLLISIONS
		if(!bIsStickyParticle && (curParticlePos.z - EPSILON) > curSceneZ)
		{
			//ParticlesColorsBufferUAV[SampleIndex] = float3(1,0,1);

			//TODO: Detect normal to push from based on particle pos on screen quadrant

			float3 curV = CurParticlesVelocityBufferUAV[SampleIndex];
			//push particle away, along it's velocity

			//currently disabled
			//curParticlePos = curParticlePos - normalize(curV) * CGConstantsCB.PointRadius;

			//push particle away from a surface
			curParticlePos.z = curParticlePos.z - CGConstantsCB.PointRadius * CGConstantsCB.DeltaTime * curSceneZ * 2.f;

			ParticlesBufferPositionCurUAV[SampleIndex] = curParticlePos;
			//curV.xy = float2(curV.x * -1.f, curV.y * -1.f);
			curV.xy = float2(curV.x * -0.5f, curV.y * -0.5f);
			CurParticlesVelocityBufferUAV[SampleIndex] = curV;
			return;
		}
		else
		{
			//ParticlesColorsBufferUAV[SampleIndex] = float3(0,1,0);
			if(bIsStickyParticle)
			{	
				//stick particles just follow the surface
				curParticlePos.z = curSceneZ;
			}
		}
#endif

#if 1 //advection along normal
		if((curParticleState & BITMASK_PARTICLE_STATE_FALLING) == 0)
		{
			const float normalAdvectionScale = 3.f /* 2.f */;
			//advect particle along surface normal
			float3 curV = CurParticlesVelocityBufferUAV[SampleIndex];
			float2 curNormal = GBufferNormalTextureSRV.SampleLevel(PointSampler, uvPos, 0).xy;
			//curV.y = min(0.f, curV.y + curNormal.y * 0.1);
			curV.y = curV.y + curNormal.y * 0.1;
			if(abs(curV.x) < abs(curNormal.x *  normalAdvectionScale * 5))
			{
				//curV.x = curNormal.x * normalAdvectionScale;
				curV.x += curNormal.x * normalAdvectionScale * 15.f;
			}
			CurParticlesVelocityBufferUAV[SampleIndex] = curV;
		}
#endif


#if 1
	//mark falling particles
	if((curParticlePos.z + CGConstantsCB.PointRadius) < (curSceneZ + fallingParticleDetectionDistanceOffset))
	{
		curParticleState = curParticleState | BITMASK_PARTICLE_STATE_FALLING;
		ParticlesStateBufferUAV[SampleIndex] = curParticleState;
	}

#endif

	}
}


#define OUTPUT_NORMALS 1

[numthreads( 8, 8, 1 )]
void GenerateNormalsFromSDFHeightMap( 
	uint2 SampleCoord : SV_DispatchThreadID )
{

	float curSDF = SDFTextureSRV[SampleCoord];

#if SDF_USE_RAYMARCH
	const float thres = CGSDFConstantsCB.Threshold;
#else
	const float thres = max(CGConstantsCB.PointRadius/2.f, CGSDFConstantsCB.Threshold);
#endif

#if SDF_USE_RAYMARCH == 0
	if(curSDF <= thres)
#endif//SDF_USE_RAYMARCH
	{
		//SDFOutputTextureUAV[SampleCoord] = float4(CGSDFConstantsCB.OutputColor.xyz, 1.f);
		float rightSDF = SDFTextureSRV[SampleCoord + uint2(1,0)];
		float upSDF = SDFTextureSRV[SampleCoord - uint2(0,1)];
		float leftSDF = SDFTextureSRV[SampleCoord - uint2(1,0)];
		float downSDF = SDFTextureSRV[SampleCoord + uint2(0,1)];

		float3 normal = float3(rightSDF - leftSDF, upSDF - downSDF, -CGSDFConstantsCB.NormalsSmoothness);
		normal = normalize(normal);
		SDFOutputNormalsTextureUAV[SampleCoord] = float4(normal.xyz, 1.f);
	}
#if SDF_USE_RAYMARCH == 0
	else
	{
		//SDFOutputTextureUAV[SampleCoord] = 0.f;
		//SDFOutputTextureUAV[SampleCoord] = SDFOutputTextureUAV[SampleCoord] * CGConstantsCB.Damping;
		SDFOutputNormalsTextureUAV[SampleCoord] = float4(0,0,-1.f, 1.f);
	}
#endif//SDF_USE_RAYMARCH

}







//===============================================
//
//				Painting Render
//
//===============================================

cbuffer SceneUniformBufferCBuffer : register(b3)
{
	matrix ViewMatrix;
	matrix ProjectionMatrix;
	matrix ViewProjectionMatrix;

	float3 DirectionWorldSpace;
	float pad;
	float3 DirectionViewSpace;
	float Intensity;
	float3 Color;
	float AmbientLightIntensity;

	float3 GCameraPosition;

	uint GNumLocalLights;
}

cbuffer UtilityConstants : register(b2)
{
    uint GPassIndex;
};

//ConstantBuffer<PaintRenderConstants> CGPaintRenderConstantsCB : register(b10);

//SRV
StructuredBuffer<float3> VertexBufferPosition : register (t0);
StructuredBuffer<float2> TexCoordBuffer : register (t1);
Texture2D<float4> ColorTexture : register (t2);
Texture2D<float3> NormalTexture : register (t3);
Texture2D<float3> PaintingTexture : register (t4);
StructuredBuffer<LightDescription> LightsBuffer : register(t5);
Texture2D<float4> FallingParticlesTexture : register (t6);

Texture2D<float3> RoughnessTexture : register (t9);

Texture2D<float3> GBufferColorTexture : register (t10);
Texture2D<float4> GBufferNormalTexture : register (t11);
Texture2D<float> GBufferRoughnessTexture : register (t12);


Texture2D<float> ShadowMapTextureSRV : register (t14);

StructuredBuffer<LightTransformStruct> LightSpaceTransformMatricesPerLightBufferSRV : register (t15);

Texture2D<float> MetalnessTexture : register (t16);
Texture2D<float> GBufferMetalTexture : register (t17);

Texture2D<float> ShadowMapPerLightTexturesSRV[] : register (t18);

//UAV

struct PaintRenderInterpolators
{
	float4 Position : SV_Position;
	float3 PositionWorld : POSWORLD;
	float3 PositionLightViewSpace : POSLIGHT;
	float4 PositionLightClipSpace : POSLIGHTCLIP;
	float2 TexCoords : TEXCOORD;
};

void PaintRenderVS(
	in uint VertexIndex : SV_VertexID,
	out PaintRenderInterpolators vsout
	)
{
	/*
	*	Position
	*/
	float3 positionWorld = VertexBufferPosition[VertexIndex];

	vsout.PositionWorld = positionWorld;
	vsout.Position = mul(float4(positionWorld, 1.0f), ViewProjectionMatrix);

	/*
	*	TexCoord
	*/
	vsout.TexCoords = TexCoordBuffer[VertexIndex];


	//Light
	vsout.PositionLightViewSpace = mul(CGPaintRenderConstantsCB.GlobalLightViewMatrix, float4(positionWorld, 1.0f));
	vsout.PositionLightClipSpace = mul(CGPaintRenderConstantsCB.GlobalLightViewProjectionMatrix, float4(positionWorld, 1.0f));


}

void Fresnel( inout float3 specular, inout float3 diffuse, float3 lightDir, float3 halfVec )
{
    float fresnel = pow(1.0 - saturate(dot(lightDir, halfVec)), 5.0);
    specular = lerp(specular, 1, fresnel);
    diffuse = lerp(diffuse, 0, fresnel);
}

float3 CalculatePointLight(LightDescription lightDesc, const float3 positionWorld, const float3 inColor, const float3 normal, const float gloss)
{
	//Dir
    float3 vToLight = lightDesc.Position - positionWorld;
    float vToLightLength = sqrt(dot(vToLight, vToLight));
    vToLight /= vToLightLength; //normalize

	float diffuse = 0.f;
	float specular = 0.f;

	//Fade
	float DistanceFalloff = saturate(1.f - pow(vToLightLength / lightDesc.Radius, 2));

	//Lambert
	if(CGPaintRenderConstantsCB.bShadeLambert)
	{
		diffuse = lightDesc.Color * lightDesc.Intensity * max(0.f, dot(vToLight, normal));
	}

	float3 halfVec = float3(0,0,0);

	if(CGPaintRenderConstantsCB.bShadeSpecular)
	{	
		if(CGPaintRenderConstantsCB.SpecularShadingMethod == 0)
		{
			//Phong specular calculation
    		float3 lightReflect = reflect(CGPaintRenderConstantsCB.GlobalLightDirection, (normal));
    		float3 viewDirection = normalize(GCameraPosition - positionWorld);
    		specular = pow(max(0, dot(lightReflect, viewDirection)), CGPaintRenderConstantsCB.SpecularPower);
		}
		else
		{
			// Blinn-Phong specular calculation
			float3 vToCam = normalize(GCameraPosition - positionWorld);
    		halfVec = normalize(vToLight + vToCam);
			
			const float specularPowerScaled = CGPaintRenderConstantsCB.SpecularPower * 8.f;
			specular = pow(max(0, dot(halfVec, normal)), specularPowerScaled);
			if(CGPaintRenderConstantsCB.SpecularShadingMethod == 2)
			{
				specular *= (specularPowerScaled + 2) / 8;
			}
		}
		specular *= CGPaintRenderConstantsCB.SpecularAmount * lightDesc.Intensity * clamp(gloss, 0., 1.);
	}

	if(CGPaintRenderConstantsCB.bShadeSpecular && CGPaintRenderConstantsCB.bShadeFresnel)
	{
		Fresnel(specular, diffuse, vToLight, halfVec);
	}

    return inColor * DistanceFalloff * (diffuse) + specular;

}

float3 FresnelSchlick(float cosTheta, float3 F0)
{
    return F0 + (1.0 - F0) * pow(saturate(1.0 - cosTheta), 5.0);
}

float3 FresnelSchlickWithRoughness(float cosTheta, float3 F0, float roughness)
{
    return F0 + (max(float(1.0 - roughness).rrr, F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float NormalDistributionGGX(float3 N, float3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
	
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
	
    return num / denom;
}

float GeometrySmith(float3 N, float3 V, float3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
	
    return ggx1 * ggx2;
}

float3 CalculatePointLightPBR(LightDescription lightDesc, const float3 positionWorld, const float3 albedo, const float3 normal, const float roughness, const float metalness)
{
	
	float3 n = normalize(normal);
	float3 v = normalize(GCameraPosition - positionWorld);
	float3 l = normalize(lightDesc.Position - positionWorld);
	float3 h = normalize(v + l);

	float distance    = length(l);
	/* 
    float attenuation = 1.0 / (distance * distance); */

	float attenuation = saturate(1.f - pow(distance / lightDesc.Radius, 2));

    float3 radiance     = lightDesc.Color * lightDesc.Intensity * attenuation;

	//calculate specular-diffuse ratio with Fresnel
	float3 F0 = float(0.04).rrr; //surface reflection at zero incidence
	F0 = lerp(F0, albedo, metalness);
	float3 F  = FresnelSchlick(max(dot(h, v), 0.0), F0);
	//float3 F  = FresnelSchlickWithRoughness(max(dot(h, v), 0.0), F0, roughness);

	//normal distribution function
	float NDF = NormalDistributionGGX(n, h, roughness);  

	//geometry function     
	float G = GeometrySmith(n, v, l, roughness);  

	//Cook-Torrance BRDF:
	float3 numerator    = NDF * G * F;
	float denominator = 4.0 * max(dot(n, v), 0.0) * max(dot(n, l), 0.0) + EPSILON;
	float3 specular     = numerator / denominator;  

	//return NDF; 
	//return G;//might be redundant, we don't have large angles
	//return F;

	//refracted diffuse vs reflected specular
	float3 kS = F;
	float3 kD = float(1.0).rrr - kS;
  
	//kD *= 1.0 - metalness;

    float NdotL = max(dot(n, l), 0.0);   
    return (kD * albedo / PI + specular) * radiance * NdotL; //final Radiance

}


float3 CalculateSpotLight(LightDescription lightDesc, const float3 positionWorld, const float3 inColor, const float3 normal, const float gloss)
{
	//Spotlight fade
	float3 vToLight = normalize(lightDesc.Position - positionWorld);
    //Spotlight fade
    float2 SpotlightAngles = lightDesc.SpotlightAngles;
    float SpotlightFalloff = dot(vToLight, -lightDesc.Direction);

     // x = 1.0f / (cos(coneInner) - cos(coneOuter)), y = cos(coneOuter)
    SpotlightAngles.x = 1.f / (SpotlightAngles.x - SpotlightAngles.y);
    SpotlightFalloff = saturate((SpotlightFalloff - SpotlightAngles.y) * SpotlightAngles.x);

    return SpotlightFalloff * CalculatePointLight(lightDesc, positionWorld, inColor, normal, gloss);
} 

float3 CalculateSpotLightPBR(LightDescription lightDesc, const float3 positionWorld, const float3 albedo, const float3 normal, const float roughness, const float metalness)
{
	//Spotlight fade
	float3 vToLight = normalize(lightDesc.Position - positionWorld);
    //Spotlight fade
    float2 SpotlightAngles = lightDesc.SpotlightAngles;
    float SpotlightFalloff = dot(vToLight, -lightDesc.Direction);

     // x = 1.0f / (cos(coneInner) - cos(coneOuter)), y = cos(coneOuter)
    SpotlightAngles.x = 1.f / (SpotlightAngles.x - SpotlightAngles.y);
    SpotlightFalloff = saturate((SpotlightFalloff - SpotlightAngles.y) * SpotlightAngles.x);

    return SpotlightFalloff * CalculatePointLightPBR(lightDesc, positionWorld, albedo, normal, roughness, metalness);
} 

float3 CalculateGlobalLight(const float3 positionWorld, const float3 color, const float3 normal, const float gloss)
{
	float3 DirToLight = -CGPaintRenderConstantsCB.GlobalLightDirection;

	float diffuse = 0.f;
	float specular = 0.f;

	if(CGPaintRenderConstantsCB.bShadeLambert)
	{
		//Back facing triangle:
		/* if (dot(normal, vsin.fragViewPos) >= 0.0f)
    	{
    	    normal = -normal;
    	} */
    	float dotProductRes = dot(DirToLight, (normal));
    	//dotProductRes = abs(dotProductRes);//For one-sided objects
		diffuse = max(0., dotProductRes);
	}

	if(CGPaintRenderConstantsCB.bShadeSpecular)
	{	
		if(CGPaintRenderConstantsCB.SpecularShadingMethod == 0)
		{
			// Blinn-Phong specular calculation
    		float3 lightReflect = reflect(CGPaintRenderConstantsCB.GlobalLightDirection, (normal));
    		float3 viewDirection = normalize(GCameraPosition - positionWorld);
    		specular = pow(saturate(dot(lightReflect, viewDirection)), CGPaintRenderConstantsCB.SpecularPower);
		}
		else
		{
			float3 vToCam = normalize(GCameraPosition - positionWorld);
    		float3 halfVec = normalize(DirToLight + vToCam);
			
			const float specularPowerScaled = CGPaintRenderConstantsCB.SpecularPower * 8.f;
			specular = pow(saturate(dot(halfVec, normal)), specularPowerScaled);
			if(CGPaintRenderConstantsCB.SpecularShadingMethod == 2)
			{
				specular *= (specularPowerScaled + 2) / 8;
			}
		}
		specular *= CGPaintRenderConstantsCB.SpecularAmount * clamp(gloss, 0., 1.);
	}

	return color * diffuse + specular;
}


float2 ScaleAndCenterTexture(float2 texCoord)
{
    // Center the texture coordinates
    float2 centeredTexCoord = texCoord - float2(0.5, 0.5);

    // Scale the texture coordinates to make it twice smaller
    float2 scaledTexCoord = centeredTexCoord * 2.f;

    // Translate the texture coordinates back to the original range [0, 1]
    float2 finalTexCoord = scaledTexCoord + float2(0.5, 0.5);

    return finalTexCoord;
}


float3 Unity_NormalBlend_Reoriented_float(float3 A, float3 B)
{
    float3 t = A.xyz + float3(0.0, 0.0, 1.0);
    float3 u = B.xyz * float3(-1.0, -1.0, 1.0);
    return (t / t.z) * dot(t, u) - u;
}

float3 BlendNormalsPD(float3 n1, float3 n2)
{
	float3 r = normalize(float3(n1.xy*n2.z + n2.xy*n1.z, n1.z*n2.z));
	return r;
	return r*0.5 + 0.5;
}

float3 BlendNormalsC(float3 n1, float3 n2)
{
	return normalize(float3(n1.xy*n2.z + n2.xy*n1.z, abs(n1.z * n2.z) * -1.f));
}

float3 BlendNormalsWhiteout(float3 n1, float3 n2)
{
	float2 pd = n1.xy/n1.z + n2.xy/n2.z; // Add the PDs
	float3 r = normalize(float3(n1.xy + n2.xy, n1.z*n2.z));
	return r;
	return r*0.5 + 0.5;
}

float3 BlendNormalsUDN(float3 n1, float3 n2)
{
	float2 pd = n1.xy/n1.z + n2.xy/n2.z; // Add the PDs
	float3 r = normalize(float3(n1.xy + n2.xy, n1.z));
	return r;
	return r*0.5 + 0.5;
}

float3 BlendNormalsLerp(float3 n1, float3 n2, float t)
{
	float2 pd = lerp(n1.xy/n1.z, n2.xy/n2.z, t);
	float3 r = normalize(float3(pd, 1));
	return r;
}

#define SHADE_GLOBAL_LIGHT 0
#define SHADE_LOCAL_LIGHT 1

float4 PaintRenderPS
	(
	in PaintRenderInterpolators psin
	) : SV_Target
{

	float3 Color = GBufferColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	float3 normal = GBufferNormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).xyz;

	float Roughness = GBufferRoughnessTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);

	#if PBR_RENDER == 0
	Roughness = (1.f - saturate(Roughness));
	#endif

#if PBR_RENDER
	float Metalness = saturate(GBufferMetalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0));
#endif

#if 1

	float4 PaintTexture = ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);


	#if PAINT_RENDER_ADVECT_SOURCE
	float originalPaintTextureA = PaintTexture.a;
	//PaintTexture.a = max(1.f, PaintTexture.a);
	
	//Metalness = 0;

	#endif//PAINT_RENDER_ADVECT_SOURCE

	PaintTexture.rgb *= CGPaintRenderConstantsCB.PaintBrightnessScale;

	//Gloss = Gloss + PaintTexture.a;
	if(PaintTexture.a > 0)
	{
	#if PBR_RENDER
		if(Metalness > 0.f)
		{
			#if PAINT_RENDER_ADVECT_SOURCE == 0
			if(PaintTexture.a < 1.f)
			{
				PaintTexture.rgb /= CGPaintRenderConstantsCB.PaintBrightnessScale;
			}
			else
			{
				Metalness = 0.1f;
			}
			#endif
			
		}
		else
		{
			//Roughness = 1.f - saturate(PaintTexture.a);

			#if PAINT_RENDER_ADVECT_SOURCE
			Roughness = max(0.05, 1.f - saturate(originalPaintTextureA));
			#else
			Roughness = max(0.05, 1.f - saturate(PaintTexture.a));
			#endif

			
			//Metalness = 0.f;
			//Metalness = CGPaintRenderConstantsCB.UniformMetalness;
		}
	#else
		Roughness = max(0.f, PaintTexture.a - 1.f);
	#endif//PBR_RENDER
	}
	
	const float alpha = saturate(PaintTexture.a);

#if PAINT_RENDER_ADVECT_SOURCE

#else
	Color = PaintTexture.rgb * alpha + Color * (1 - alpha);
#endif//PAINT_RENDER_ADVECT_SOURCE

	float3 paintNormal = NormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	paintNormal = normalize(paintNormal);
	paintNormal.xy *= -1.f;

	normal.z = abs(normal.z);
	paintNormal.z = abs(paintNormal.z);
	//normal = BlendNormalsPD(paintNormal, normal);
	normal = BlendNormalsLerp(paintNormal, normal, 1 - alpha);
	//normal = normalize(paintNormal.rgb * alpha + normal * (1 - alpha));
	normal.z *= -1.f;

	
#if PAINT_RENDER_ADVECT_SOURCE
	//normal = paintNormal;
	//normal.z *= -1.f;

#endif

#endif
	
#if 1
	//const float4 fallingPaintColor = FallingParticlesTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	float4 fallingPaintColor = FallingParticlesTexture.SampleLevel(PointSampler, psin.TexCoords, 0);
	/* if(fallingPaintColor.a > 0.)
	{
		Color = fallingPaintColor;
		Gloss = 1;
		normal = paintNormal;
		normal.z *= -1.f;
	} */

	const float fallingPaintAlpha = saturate(fallingPaintColor.a);
	const bool bIsFallingParticle = fallingPaintColor.a > 0;
	if(bIsFallingParticle)
	{
#if PAINT_RENDER_ADVECT_SOURCE
		Color = lerp(Color, fallingPaintColor.rgb, 0.75f);
#else
		fallingPaintColor.rgb *= CGPaintRenderConstantsCB.PaintBrightnessScale;
		Color = fallingPaintColor.rgb * fallingPaintAlpha + Color * (1 - fallingPaintAlpha);
#endif//PAINT_RENDER_ADVECT_SOURCE

	#if PBR_RENDER
		//Roughness = 1.f - saturate(fallingPaintColor.a);
		Roughness = CGPaintRenderConstantsCB.UniformRoughness;
		Metalness = CGPaintRenderConstantsCB.UniformMetalness;
	#else
		Roughness = fallingPaintColor.a;
	#endif//PBR_RENDER
		normal = paintNormal;
		normal.z *= -1.f;
	}

#endif

	float3 shadedColor = float3(0,0,0);

	//Direct Global Light
	#if SHADE_GLOBAL_LIGHT
	{
		float3 GlobalLight = CalculateGlobalLight(psin.PositionWorld, Color, normal, Gloss);
		shadedColor += GlobalLight;
	}
	#endif//SHADE_GLOBAL_LIGHT

	#if SHADE_LOCAL_LIGHT
	for(uint i = 0; i < CGPaintRenderConstantsCB.NumLights; i++)
	{
		float shadow = 1.f;

//#if PAINT_RENDER_ADVECT_SOURCE == 0
		if(!bIsFallingParticle)
		{
			shadow = ShadowMapPerLightTexturesSRV[i].SampleLevel(PointSampler, psin.TexCoords, 0);
			shadow = clamp(shadow, 0., 1.f);
		}
	//#endif

		#if 0 && PBR_RENDER
		Roughness = CGPaintRenderConstantsCB.UniformRoughness;
		Metalness = CGPaintRenderConstantsCB.UniformMetalness;
		#else
		Roughness = clamp(Roughness, CGPaintRenderConstantsCB.RoughnessMinMax.x, CGPaintRenderConstantsCB.RoughnessMinMax.y);
		Metalness = clamp(Metalness, CGPaintRenderConstantsCB.MetalnessMinMax.x, CGPaintRenderConstantsCB.MetalnessMinMax.y);
		#endif

		LightDescription curLightDesc = LightsBuffer[i];
		if(curLightDesc.Type == LIGHT_TYPE_POINT)
		{
		#if PBR_RENDER
			shadedColor += CalculatePointLightPBR(LightsBuffer[i], psin.PositionWorld, Color, normal, Roughness, Metalness) * shadow;
		#else
			shadedColor += CalculatePointLight(LightsBuffer[i], psin.PositionWorld, Color, normal, Roughness) * shadow;
		#endif//PBR_RENDER
		}
		else
		{
		#if PBR_RENDER
			shadedColor += CalculateSpotLightPBR(LightsBuffer[i], psin.PositionWorld, Color, normal, Roughness, Metalness) * shadow;
		#else
			shadedColor += CalculateSpotLight(LightsBuffer[i], psin.PositionWorld, Color, normal, Roughness) * shadow;
		#endif//PBR_RENDER
		}
	}
	#endif//SHADE_LOCAL_LIGHT

#if 1
	float3 ambient = Color * CGPaintRenderConstantsCB.AmbientLight /* * AO */;
    shadedColor += ambient;
#endif

	#if 1
	if(CGPaintRenderConstantsCB.bVisualizeNormals)
	{
		/* normal.xy = MapToRange(normal.xy, -1,1, 0.f,1.f);
		normal.z = MapToRange(normal.z, -1,1, 0.f,1.f); */
		shadedColor = normal;
	}
	#endif

	//shadedColor *= 0.25f;


	return float4(shadedColor.rgb,1);

}
















//===============================================
//
//				Background Render
//
//===============================================


struct GBufferRenderInterpolators
{
	float4 Position : SV_Position;
	float3 PositionWorld : POSWORLD;
	float2 TexCoords : TEXCOORD;
};

struct VertexInputRest
{
	float3 Normal;
	float3 Tangent;
	float3 Bitangent;
	float2 TexCoord;
};

//SRV
StructuredBuffer<uint> IndexBuffer : register (t7);
StructuredBuffer<VertexInputRest> VertexBufferRest : register (t8);

float3 DecodeNormalTextureRotated(float3 normTex)
{
	//swap
	float cachedY = normTex.y;
	normTex.y = normTex.x;
	normTex.x = cachedY;

	normTex.x = 1.0f - normTex.x;//Flip Y coord to align with DirectX UV Basis
	normTex.x = normTex.x * 2.0f - 1.0f;
    normTex.y = normTex.y * 2.0f - 1.0f;
    normTex.z = normTex.z * 2.0f - 1.0f;
	normTex.z *= -1.f;
    return normalize(normTex);
}

float3 DecodeNormalTexture(float3 normTex)
{
	normTex.y = 1.0f - normTex.y;//Flip Y coord to align with DirectX UV Basis
	normTex.x = normTex.x * 2.0f - 1.0f;
    normTex.y = normTex.y * 2.0f - 1.0f;
    normTex.z = normTex.z * 2.0f - 1.0f;
	normTex.z *= -1.f;
    return normalize(normTex);
}

void BackgroundMeshRenderVS(
	in uint VertexIndex : SV_VertexID,
	out GBufferRenderInterpolators vsout
	)
{
	uint VertexId = IndexBuffer[VertexIndex];
	/*
	*	Position
	*/
	float3 positionWorld = VertexBufferPosition[VertexId];
	vsout.PositionWorld = positionWorld;
	vsout.Position = mul(CGPaintRenderConstantsCB.PaintSceneViewProjectionMatrix, float4(positionWorld, 1.0f));
	/*
	*	Rest
	*/
	VertexInputRest inputRest = VertexBufferRest[VertexId];
	/*
	*	TexCoord
	*/
	vsout.TexCoords = inputRest.TexCoord;
}

void BackgroundWallRenderVS(
	in uint VertexIndex : SV_VertexID,
	out GBufferRenderInterpolators vsout
	)
{
	/*
	*	Position
	*/
	float3 positionWorld = VertexBufferPosition[VertexIndex];

	vsout.PositionWorld = positionWorld;
	vsout.Position = mul(CGPaintRenderConstantsCB.PaintSceneViewProjectionMatrix, float4(positionWorld, 1.0f));

	/*
	*	TexCoord
	*/
	vsout.TexCoords = TexCoordBuffer[VertexIndex];
}

void BackgroundFrameRenderPS
	(
	in GBufferRenderInterpolators psin,
	out float3 outColor : SV_Target0,
    out float3 outNormal : SV_Target1,
    out float outRoughness : SV_Target2,
	out float outZPosition : SV_Target3,
	out uint outStencil : SV_Target4
	,out float outMetal : SV_Target5
	) 
{

	//return ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	//return float4(DecodeNormalTexture(NormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).xyz), 1);

	float2 texCoordsScaled = psin.TexCoords * 2.f;

	float3 shadedColor = float3(0,0,0);

	float3 Color = ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).rgb;

	float3 normal = DecodeNormalTextureRotated(NormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).xyz);

	float Roughness = RoughnessTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).r;

	//float Metalness = MetalnessTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).r;
	float Metalness = RoughnessTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).g;

	//Metalness = max(0.25, Metalness);

	outColor = Color;
	outNormal = normal;
	outRoughness = Roughness;
	outZPosition = psin.PositionWorld.z;
	outStencil = 1;//mark frame pixels
	outMetal = Metalness;
}

void BackgroundWallRenderPS
	(
	in GBufferRenderInterpolators psin,
	out float3 outColor : SV_Target0,
    out float3 outNormal : SV_Target1,
    out float outRoughness : SV_Target2,
	out float outZPosition : SV_Target3
	)
{

	//return ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	//return float4(DecodeNormalTexture(NormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).xyz), 1);

	float3 Albedo = ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);
	//float3 Color = float3(0.5,0.5,0.5);
	//float3 Color = float3(0.01,0.35,0.15) * 0.75;
	//float3 Color = (Albedo + float3(0.01,0.35,0.15)) * 0.4f;
	//float3 Color = (Albedo + float3(0.13,0.67,0.41)) * 0.4f;
	float3 Color = (Albedo + float3(0.13,0.32,0.67)) * 0.4f;
	//float3 Color = (Albedo) * 0.9f;
	float3 normal = DecodeNormalTexture(NormalTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).xyz);
	float Roughness = RoughnessTexture.SampleLevel(LinearSampler, psin.TexCoords, 0).r;
	//const float Gloss = 1.f;

	outColor = Color;
	outNormal = normal;
	outRoughness = Roughness;
	outZPosition = psin.PositionWorld.z;
}

void PaintSourceRenderPS
	(
	in GBufferRenderInterpolators psin,
	out float3 outColor : SV_Target0,
    out float3 outNormal : SV_Target1,
    out float outRoughness : SV_Target2,
	out float outZPosition : SV_Target3
	)
{

	float3 Color = ColorTexture.SampleLevel(LinearSampler, psin.TexCoords, 0);

	outColor = Color;
}

/* 
*	 Depth Render
*/

struct DepthRenderInterpolants
{
	float4 Position : SV_Position;
	float3 PositionLightView : POSVIEW;
};

void DepthRenderVS(
	in uint VertexIndex : SV_VertexID,
	out DepthRenderInterpolants vsout
	)
{
	/*
	*	Position
	*/
	float3 positionWorld = VertexBufferPosition[VertexIndex];

	//get light matrix
	LightTransformStruct lightTransform = LightSpaceTransformMatricesPerLightBufferSRV[GPassIndex];
	vsout.PositionLightView = mul(lightTransform.View, float4(positionWorld, 1.0f));
	vsout.Position = mul(lightTransform.ViewProjection, float4(positionWorld, 1.0f));
}

void DepthRenderIndexedVS(
	in uint VertexIndex : SV_VertexID,
	out DepthRenderInterpolants vsout
	)
{
	uint VertexId = IndexBuffer[VertexIndex];
	/*
	*	Position
	*/
	float3 positionWorld = VertexBufferPosition[VertexId];

	//get light matrix
	LightTransformStruct lightTransform = LightSpaceTransformMatricesPerLightBufferSRV[GPassIndex];
	vsout.PositionLightView = mul(lightTransform.View, float4(positionWorld, 1.0f));
	vsout.Position = mul(lightTransform.ViewProjection, float4(positionWorld, 1.0f));
}

void DepthRenderPS
	(
	in DepthRenderInterpolants psin,
	out float outViewZ : SV_Target0
	) 
{

	outViewZ = psin.PositionLightView.z;

}



RWTexture2D<float> ShadowVisibilityTextureUAV : register(u7);

[numthreads( 8, 8, 1 )]
void GenerateShadowVisibility( 
	uint2 SampleCoord : SV_DispatchThreadID )
{
	
	const float shadowTextureSize = 2048;
	const uint2 SampleCoordInv = uint2(SampleCoord.x, (shadowTextureSize - 1) - SampleCoord.y);
	float2 uv = SampleCoord / (shadowTextureSize);
	float2 uvInv = float2(uv.x, 1.f - uv.y);

	float3 worldPos;

	worldPos.xy = MapToRange(uv.xy, 0.f, 1.f, GetVelocityFieldBoundsStart(), GetVelocityFieldBoundsEnd());

	//sample world space texture and map it to light space
	worldPos.z = GBufferWorldSpaceZSRV.SampleLevel(PointSampler, uvInv, 0).r;

	//map to light space
	//float4 lightViewPos = mul(CGPaintRenderConstantsCB.GlobalLightViewMatrix, float4(worldPos, 1.0f));
	//float curPixelDepthLightSpace = lightViewPos.z;

	LightTransformStruct lightTransform = LightSpaceTransformMatricesPerLightBufferSRV[GPassIndex];

	float4 lightClipPos = mul(lightTransform.ViewProjection, float4(worldPos, 1.0f));
	float curPixelDepthLightSpace = lightClipPos.z / lightClipPos.w;

	//NDC to TexCoord
	float2 projCoords = lightClipPos.xy / lightClipPos.w;
	projCoords.xy = projCoords.xy * 0.5 + 0.5;
	projCoords.y = 1.f - projCoords.y;
	const float kBias = CGPaintRenderConstantsCB.ShadowDepthBias;
	float bias = kBias;

	float3 normal = GBufferNormalTexture.SampleLevel(LinearSampler, uvInv, 0);

	if(CGPaintRenderConstantsCB.bSlopeCorrectShadowBias)
	{
		bias = max(kBias * 10.f * (1.0 - dot(normal, -CGPaintRenderConstantsCB.GlobalLightDirection)), kBias);  
	}

	float shadow = 1.f;

#if 0
	float lightDepth = ShadowMapTextureSRV.SampleLevel(PointSampler, projCoords.xy, 0);
	shadow = curPixelDepthLightSpace <= (lightDepth + bias) ? 1.f : 0.f;
#else //PCF
	shadow = 0.f;
	const float ShadowMapSize = 2048;
	float2 texelSize = 1.0 / ShadowMapSize;
	float counter = 0.f;
	int start = -4;
	int end = 4;
	for(int x = start; x <= end; ++x)
	{
	    for(int y = start; y <= end; ++y)
	    {
			float lightDepth = ShadowMapTextureSRV.SampleLevel(PointSampler, projCoords.xy + float2(x,y) * texelSize, 0);
			shadow += curPixelDepthLightSpace <= (lightDepth + bias) ? 1.f : 0.f;
			counter += 1.f;
	    }    
	}
	shadow /= counter;
#endif
	if(projCoords.x > 1.f || projCoords.x < 0.f || projCoords.y > 1.f || projCoords.y < 0.f)
	{
		shadow = 1.f;
	}
	shadow = clamp(shadow, 0.5f, 1.f);
	ShadowVisibilityTextureUAV[SampleCoordInv] = shadow;
}
import { init } from "@sentry/browser";
import { Vector2 } from "../core/types";
import { ParticleDesc } from "../scenes/payne/paintingMain/particleDesc";
import { scDeclareSceneUniformBuffer,  scGetBasePassInterpolatorDecl,  scGetBasePassUniformTexDeclPS, scGetBasePassUniformVarDeclPS } from "./shaderBasePass";
import { scGetCommonFuncsDeclarations, scGetGIFuncsDeclarations, scGetProceduralCubeVertsDecl, scGetSphericalFuncsDeclarations, scPBRDeclarations } from "./shaderCommon";
import { scGetSDFFuncsDeclarations, scGetSDFShadowDeclarations } from "./shaderSDF";
import { scGetBasePassLightShading, scGetDestrTileCommonFuncsDecl } from "./shaderDestrTiles";
import { scDeclareCurlNoiseFuncs } from "./shaderParticles";


//=============================================================================================================================
//
// 														PARTICLES
//
//=============================================================================================================================



//=============================================================================================================================
// 																											UPDATE
//=============================================================================================================================


function scParticleDefines(inDesc : ParticleDesc)
{
	return `#define MODEL_CUBE ` + (inDesc.bProceduralCube ? `1` : `0`)
	+
	` \n `
	+
	`#define ATLAS_TEXTURE ` + (inDesc.bAtlasTexture ? `1` : `0`)
	+
	` \n `
	+
	`#define MODEL_MESH ` + (inDesc.bUsesMesh ? `1` : `0`)
	+
	` \n `
	+
	`#define PATH_SMOKE ` + (inDesc.bSmokePath ? `1` : `0`)
	+
	` \n `
	+
	`#define CURL_NOISE ` + (inDesc.bUseCurlNoiseVelocity ? `1` : `0`)
	+
	` \n `
	
	;
}


export function GetDestParticleUpdateShaderVS(numTiles : Vector2, tileSize : Vector2, inDesc : ParticleDesc) {
    	return (
        /* glsl */ `#version 300 es
  
		precision highp float;
		precision mediump sampler2D;

		` +scParticleDefines(inDesc)+/* glsl */`

		` +scDeclareSceneUniformBuffer()+/* glsl */`
  
		layout(location = 0) in vec3 inPosition;
	  	layout(location = 1) in vec3 inVelocity;
	  	layout(location = 2) in float inAge;

		uniform float RandNumSeed;
		uniform ivec2 ParticleCurSpawnIndexRange; //Particles withing this range are force-respawned

		uniform sampler2D TileHealthStateTexture;
  
		out vec3 outPosition;
		out vec3 outVelocity;
		out float outAge;

		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetDestrTileCommonFuncsDecl(numTiles, tileSize)+ /* glsl */`
			
		vec3 GetSpawnPos()
		{
			//From input pos
			vec2 posNDC = SceneBuffer.UserInputPosNDC.xy;
			//vec3 spawnPos = TransformNDCToWorld(posNDC, -SceneBuffer.CameraPosition.z);
			vec3 spawnPos = SceneBuffer.IntersectionPosWS;
			spawnPos.z = min(-0.01, spawnPos.z);
			return spawnPos;
		}

	#if CURL_NOISE

		`+scDeclareCurlNoiseFuncs()+/* glsl */`

		vec3 GetCurlNoiseVelocity(vec3 pos, vec3 curVel, float age)
		{
			float curVelLength = length(curVel);
			vec3 noisePosition = pos;
			noisePosition *= 1.0;
			float noiseTime = SceneBuffer.Time * 0.00025000;
			vec4 xNoisePotentialDerivatives = vec4(0.0);
    		vec4 yNoisePotentialDerivatives = vec4(0.0);
    		vec4 zNoisePotentialDerivatives = vec4(0.0);
			float curl = 2.5; 

			float s = 1.0;
			if(age > 1.0)
			{
				s = mix(0.1, 0.01, min(1.0, (age - 1.0) * 0.25));
			}
			else
			{
				s = mix(0.5, 0.1, min(1.0, age));
			}

			float noiseVelScale = 10.0 * s;

			const int OCTAVES = 3;
			for (int i = 0; i < OCTAVES; ++i) {
				float scale = (1.0 / 2.0) * pow(2.0, float(i));
				float noiseScale = pow(curl, float(i));
				if (curl == 0.0 && i == 0) {
					noiseScale = 1.0;
				}
				xNoisePotentialDerivatives += simplexNoiseDerivatives(vec4(noisePosition * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				yNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(123.4, 129845.6, -1239.1)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				zNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(-9519.0, 9051.0, -123.0)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
			}
			vec3 noiseVelocity = vec3(
				zNoisePotentialDerivatives[1] - yNoisePotentialDerivatives[2],
				xNoisePotentialDerivatives[2] - zNoisePotentialDerivatives[0],
				yNoisePotentialDerivatives[0] - xNoisePotentialDerivatives[1]
				) * 0.07500000;

			noiseVelocity *= noiseVelScale * 0.5;

			if(curVelLength > 0.0)
			{
				float alignment = dot(normalize(noiseVelocity), (curVel / curVelLength));
				if(alignment < 0.5)
				{
					noiseVelocity *= 1.0 + curVelLength;
				}
				
			}

			return noiseVelocity;
		}
	#endif

		vec3 GetSpawnVelocity()
		{
			//Get Random Dir
			float thetaMin = 0.5;
			float thetaMax = max(thetaMin + 0.01, PIDIVTWO - 0.05);

			#if (ATLAS_TEXTURE && PATH_SMOKE)
			thetaMin = 0.0;
			thetaMax = 0.1 * 0.75;
			#endif

			float rnd = hashf(SceneBuffer.Time + float(gl_VertexID) * 23.731 + RandNumSeed * 1.12);
			float rnd2 = hashf(SceneBuffer.Time * 0.98 + float(gl_VertexID) * 7.12 + 1.1 + RandNumSeed * 17.3);
			float phi = 2.0 * PI * rnd;
			
			float theta = MapToRange(rnd2, 0.0, 1.0, thetaMin, thetaMax);
			//theta = 1.4;
			float cosTheta = cos(theta);
			float sinTheta = sin(theta);

			vec3 V;
	    	V.x = cos(phi) * sinTheta;
	    	V.y = sin(phi) * sinTheta;
	    	V.z = -cosTheta * float(`+inDesc.InitVelocityZScale+/* glsl */`);
			V *= 0.75;

			V.y += float(`+inDesc.InitVelHeightBoost+/* glsl */`) * (1.0 - max(0.0, V.y));

			V *= float(`+inDesc.InitVelMin+/* glsl */`) + hashf(SceneBuffer.Time + float(gl_VertexID) * 3.7082) * 3.0 * float(`+inDesc.RandVelScale+/* glsl */`);

		#if (PATH_SMOKE == 0)
			if(gl_VertexID % 23 == 0)
			{
				V *= 2.0;
			}
			else if(gl_VertexID % 17 == 0)
			{
				V.z *= 2.5;
			}
		#endif


			return V;
		}

		vec3 GetCurForce()
		{
			return vec3(0.0, -10.0 * float(`+inDesc.GravityForceScale+/* glsl */`), 0.0);
		}

		void main()
		{

			const vec2 kEmmiterTextureSpawnRange = vec2(float(` +
					inDesc.EmitterTextureSpawnRange.x +
					/* glsl */ `), float(` +
						inDesc.EmitterTextureSpawnRange.y +
					/* glsl */ `));


			float curAge = inAge;
			
			const float ParticleLifeTimeSec = float(`+inDesc.ParticleLifeTimeSec+/* glsl */`);
			int particleId = int(gl_VertexID);
			bool bCanSpawn = false;
			bCanSpawn = particleId >= ParticleCurSpawnIndexRange.x && particleId < ParticleCurSpawnIndexRange.y;

			vec3 spawnPos = GetSpawnPos();

			if(bCanSpawn)
			{
				ivec3 curTileId = GetTileIdFromPosWS(spawnPos);
				if(curTileId.z < 0)
				{
					bCanSpawn = false;
				}
				else
				{
					float health = texelFetch(TileHealthStateTexture, curTileId.xy, 0).r;
					if(health > kEmmiterTextureSpawnRange.x && health < kEmmiterTextureSpawnRange.y)
					{

					}
					else
					{
						bCanSpawn = false;
					}
				}
			}

			if(bCanSpawn)
			{
				//Respawn
				
				vec3 velocity = GetSpawnVelocity();
				outPosition = spawnPos + velocity * SceneBuffer.DeltaTime * 5.0 * float(`+inDesc.InitPosRandOffsetScale+/* glsl */`);
				outVelocity = velocity;
				outAge = 0.0;
				return;
			}
			else if(curAge <= ParticleLifeTimeSec)
			{

				float ageNorm = curAge / ParticleLifeTimeSec;

				//Update
				vec3 curPos = inPosition;
				vec3 curVel = inVelocity;

				vec3 netForce = vec3(0.0);
				netForce += GetCurForce();

			#if CURL_NOISE
				netForce += GetCurlNoiseVelocity(curPos, curVel, curAge);
			#endif


				curVel += netForce * SceneBuffer.DeltaTime;
				float dampingScale = 0.2 + float(`+inDesc.VelocityDampingAdd+/* glsl */`);
				curVel /= (1.0 + dampingScale * SceneBuffer.DeltaTime);
				curPos += curVel * SceneBuffer.DeltaTime;

				outPosition = curPos;
				outVelocity = curVel;
				outAge = curAge + SceneBuffer.DeltaTime;
			}
			else
			{
				outPosition = vec3(-100.0);
				outVelocity = vec3(0.0);
				outAge = curAge;
			}
		}`
    );
}


//=============================================================================================================================
// 																											RENDER
//=============================================================================================================================
export function scGetBasePassParticleInterpolatorDecl(bPixelShader : boolean)
{
	const inOrOut = bPixelShader ? `in` : `out`;
	return /* glsl */`
		`+inOrOut+/* glsl */` vec2 interpolatorTexCoords;
		`+inOrOut+/* glsl */` vec3 interpolatorWorldSpacePos;
		`+inOrOut+/* glsl */` vec3 interpolatorTBNNormal;
		flat `+inOrOut+/* glsl */` uint interpolatorFaceId;
	`;
}

export function GetShaderSourceBasePassParticleRenderVS(inDesc : ParticleDesc) {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;

		` +scParticleDefines(inDesc)+/* glsl */`

		layout(location = 0) in vec3 inPosition;
	  	layout(location = 1) in vec3 inVelocity;
	  	layout(location = 2) in float inAge;
	#if MODEL_MESH
		layout(location = 3) in vec3 inVBPosition;
		layout(location = 4) in vec2 inVBTexCoords;
		layout(location = 5) in vec3 inVBNormals;
		uniform vec3 PreTranslate;
	#endif//MODEL_MESH

		uniform float RandNumSeed;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		`+scGetBasePassParticleInterpolatorDecl(false)+ /* glsl */`

		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`
	
		void main()
		{
			vec3 pos = vec3(-100.0);
			vec2 texCoords = vec2(0.0);
			vec3 normal = vec3(0.0, 0.0, -1.0);
			

			const float ParticleLifeTimeSec = float(`+inDesc.ParticleLifeTimeSec+/* glsl */`);
			bool bParticleDead = inAge > ParticleLifeTimeSec;
			bParticleDead = false;

			//compute local space pos
			if(!bParticleDead)
			{
		#if (MODEL_MESH)
				pos = inVBPosition.xyz;
				pos -= PreTranslate;
				normal = inVBNormals;
				texCoords = inVBTexCoords;

				interpolatorFaceId = dot(normal, vec3(0.0, 0.0, -1.0)) > 0.9 ? CUBE_FACE_FRONT : CUBE_FACE_LEFT;
		#else
				uint vertId = uint(gl_VertexID);
				pos = CubePositions[vertId];
				normal = CubeNormals[vertId];
				texCoords = CubeTexCoords[vertId];
				
			#if (MODEL_CUBE)
				interpolatorFaceId = GetCurCubeFaceIndex(vertId);

				#if 1//Cube pos offset
					uint uniqueIndexId = VertexToIndexMap[vertId];
					if(uniqueIndexId == 0u || uniqueIndexId == 4u)
					{
						vec2 rnd2 = hash2(float(gl_InstanceID) * 9.671);

						vec2 offset;
						offset.x = MapToRange(rnd2.x, 0.0, 1.0, 0.1, 1.0);
						offset.y = MapToRange(rnd2.y, 0.0, 1.0, 0.0, 1.0);

						pos.x += offset.x;
						pos.y += offset.y;
					}
					else if(uniqueIndexId == 1u || uniqueIndexId == 5u)
					{
						vec2 rnd2 = hash2(float(gl_InstanceID) * 12.671);
						float offsetX = MapToRange(rnd2.x, 0.0, 1.0, 0.75, 1.85);
						float offsetY = MapToRange(rnd2.y, 0.0, 1.0, 0.0, 0.25);
						pos.x += offsetX;
						pos.y += offsetY;
					}
				#endif//Cube pos offset
			#endif//MODEL_CUBE
		#endif//MODEL_MESH

				//Get random texture from atlas
				#if (ATLAS_TEXTURE)
				{
					texCoords.x *= 1.0/3.0;
					float instanceIdf = float(gl_InstanceID % 3);
					texCoords.x += (1.0/3.0) * instanceIdf;
				}
				#endif

				vec3 size = vec3(0.05 * 0.8);
			#if (!MODEL_MESH)
				size.z *= 0.15;
			#else
				size.z *= 3.0;
			#endif
				size *= float(`+inDesc.SizeScale+/* glsl */`);
				//Random size per instance
			#if 1
				{
					float scale = MapToRange(hash(uint(gl_InstanceID)), 0.0, 1.0, float(`+inDesc.SizeRandMinMax.x+/* glsl */`), float(`+inDesc.SizeRandMinMax.y+/* glsl */`));
					size *= scale;

					scale = MapToRange(hash(uint(floor(float(gl_InstanceID) * 1.821))), 0.0, 1.0, float(`+inDesc.RandHeightSizeMinMax.x+/* glsl */`), float(`+inDesc.RandHeightSizeMinMax.y+/* glsl */`));
					size.y *= scale;
				}
			#endif
				
				pos.x *= size.x;
				pos.y *= size.y;
				pos.z *= size.z;

			#if (MODEL_CUBE)
				texCoords = pos.xy;
				texCoords.x *= 2.0;
			#endif

				//Get random instance orientation 2d
			#if 1
				{
					float angle = hash(uint(gl_InstanceID) * 3u + uint(RandNumSeed * 13.2)) * PI * 2.0;

					float rotationSpeed = MapToRange(hashf(float(gl_InstanceID) * 0.32), 0.0, 1.0, 0.0, 2.0 * float(`+inDesc.RotationSpeedScale+/* glsl */`));

					angle += float(gl_InstanceID % 2 == 0 ? -1.0 : 1.0) * inAge * 1.2 * rotationSpeed;

				#if (MODEL_CUBE || MODEL_MESH) //3D
				mat3 rotationMatrix = MatrixRotateAroundAxis(vec3(0.0, 0.0, 1.0), angle);
					float yawAngleSpeed = MapToRange(hashf(float(gl_InstanceID) * 7.12), 0.0, 1.0, 0.0, 0.5 * float(`+inDesc.RotationSpeedScale+/* glsl */`));
					rotationMatrix *= MatrixRotateAroundAxis(vec3(0.0, 1.0, 0.0), angle * yawAngleSpeed);
					pos = rotationMatrix * pos;
					normal = rotationMatrix * normal;
				#else //2D
					float cosTheta = cos(angle);
					float sinTheta = sin(angle);
					mat2 rotationMatrix = mat2(
						cosTheta, sinTheta,
						-sinTheta, cosTheta
					);
					pos.xy = rotationMatrix * pos.xy;
					normal.xy = rotationMatrix * normal.xy;
				#endif//(MODEL_CUBE || MODEL_MESH)
				}
			#endif

				//compute center pos based on instance
				{
					pos += inPosition;
				}

				interpolatorWorldSpacePos = pos;

				//Camera Space
				pos.xyz -= SceneBuffer.CameraPosition; 
				//Projection
				pos.xy *= SceneBuffer.CameraZoom;
				pos.x /= SceneBuffer.ScreenRatio;

				interpolatorTexCoords = texCoords;

				interpolatorTBNNormal = normal;
			}

			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

		}`;
}

export function GetShaderSourceBasePassParticlePS(inDesc : ParticleDesc) {
    return /* glsl */ `#version 300 es
		
		precision highp float;
		precision highp sampler2D;
		precision highp samplerCube;

		` +scParticleDefines(inDesc)+/* glsl */`
	
		layout(location = 0) out vec3 OutColor;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		` +scGetBasePassUniformVarDeclPS()+/* glsl */`

		uniform sampler2D AlbedoTexture;
		uniform sampler2D MaskTexture;

		`+scGetBasePassParticleInterpolatorDecl(true)+ /* glsl */`

		`+scGetCommonFuncsDeclarations()+ /* glsl */`
		`+scGetSphericalFuncsDeclarations()+ /* glsl */`
		`+scPBRDeclarations(true)+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`

		void main()
		{	
			vec3 light = vec3(0.0);

			PBRSettings s;

			s.Albedo = texture(AlbedoTexture, interpolatorTexCoords.xy).rgb;

			if(bEmitter > 0)
			{
				OutColor = s.Albedo * float(`+inDesc.Brightness+/* glsl */`); return;
			}

			s.Normal = normalize(interpolatorTBNNormal);

			s.Roughness = 0.65;

			s.Metalness = 0.0;

			s.PixelPos = interpolatorWorldSpacePos;
			s.vToCam = normalize(SceneBuffer.CameraPosition.xyz - interpolatorWorldSpacePos);
			

		#if (!MODEL_CUBE && !MODEL_MESH)
			float mask = texture(MaskTexture, interpolatorTexCoords.xy).r;
			if(mask < 0.1)
			{
				discard;
			}
			s.Roughness = 0.85;
		#else
			if (interpolatorFaceId != CUBE_FACE_FRONT)
			{
				s.Albedo = pow(s.Albedo, vec3(0.75));
				s.Roughness = 0.75;
			}
		#endif

			vec3 indirectDiffuse = SceneBuffer.LightColor;
			indirectDiffuse = vec3(0.5);
			light += vec3(indirectDiffuse * s.Albedo / PI * 1.0) * (max(-0.3, SceneBuffer.AmbientLight) + 0.3);

		#if (!MODEL_CUBE && !MODEL_MESH)

			s.Albedo *= vec3(1.0, 0.7, 0.5);

			s.LightPos = SceneBuffer.SpotlightPos2;
			vec3 lightTargetPos = s.LightPos + SceneBuffer.SpotlightTargetPos2;

			float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, normalize(lightTargetPos - s.LightPos), s.PixelPos, 1.0, 0.75);
			float d = length(s.LightPos - s.PixelPos);

			//float attenuation = pow(clamp(1.f - pow(d / s.LightRadius, 2.0), 0.f, 1.f), 2.0);
			float num = pow(clamp(1.0 - pow(d / 10.0, 4.0), 0.0, 1.0), 2.0);
			float denom =  (d * d + 1.0);
			float attenuation = num / denom;
			SpotlightFalloff = MapToRange(SpotlightFalloff, 0.0, 1.0, 0.25, 1.0);
			light = s.Albedo * SpotlightFalloff * SceneBuffer.LightIntensity1 * 0.5 * attenuation;

		#else

			s.Albedo *= 1.2;

			`+scGetBasePassLightShading()+ /* glsl */`

		#endif

			
			

			OutColor = light;

			
		}`;
	}




	
export function GetShaderSourceBasePassSmokeParticleRenderVS(inDesc : ParticleDesc) {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;

		` +scParticleDefines(inDesc)+/* glsl */`

		layout(location = 0) in vec3 inPosition;
	  	layout(location = 1) in vec3 inVelocity;
	  	layout(location = 2) in float inAge;
	#if MODEL_MESH
		layout(location = 3) in vec3 inVBPosition;
		layout(location = 4) in vec2 inVBTexCoords;
		layout(location = 5) in vec3 inVBNormals;
		uniform vec3 PreTranslate;
	#endif//MODEL_MESH

		uniform float RandNumSeed;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		`+scGetBasePassParticleInterpolatorDecl(false)+ /* glsl */`
		flat out float interpolatorAge;
		flat out float interpolatorFrameIndex;
		flat out int interpolatorInstanceId;

		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`
	
		void main()
		{

			interpolatorInstanceId = gl_InstanceID;

			vec3 pos = vec3(-100.0);
			vec2 texCoords = vec2(0.0);
			
			const float ParticleLifeTimeSec = float(`+inDesc.ParticleLifeTimeSec+/* glsl */`);
			bool bParticleDead = inAge > ParticleLifeTimeSec;
			bParticleDead = false;

			//compute local space pos
			if(!bParticleDead)
			{
				uint vertId = uint(gl_VertexID);
				pos = CubePositions[vertId];
				texCoords = CubeTexCoords[vertId];
				
				vec3 size = vec3(1.0);

				size *= float(`+inDesc.SizeScale+/* glsl */`);
				//Random size per instance
			#if 1
				{
					float scale = MapToRange(hash(uint(gl_InstanceID)), 0.0, 1.0, float(`+inDesc.SizeRandMinMax.x+/* glsl */`), float(`+inDesc.SizeRandMinMax.y+/* glsl */`));
					size *= scale;
				}
			#endif

			#if ATLAS_TEXTURE
				//fade in
				{

					const float minSize = 0.5;
					float ageT = minSize + min(0.5, pow(inAge, 0.5));
					size *= ageT;
				}
			#endif
				
				pos.x *= size.x;
				pos.y *= size.y;
				pos.z *= size.z;

				//Get random instance orientation 2d
			#if 1
				{
					float angle = hash(uint(gl_InstanceID) * 3u + uint(RandNumSeed * 13.2)) * PI * 2.0;

					float rotationSpeed = MapToRange(hashf(float(gl_InstanceID) * 0.32), 0.0, 1.0, 1.0, 2.0 * float(`+inDesc.RotationSpeedScale+/* glsl */`));

					#if ATLAS_TEXTURE
						rotationSpeed *= 0.1;
					#endif

					angle += float(gl_InstanceID % 2 == 0 ? -1.0 : 1.0) * inAge * 1.2 * rotationSpeed;

					float cosTheta = cos(angle);
					float sinTheta = sin(angle);
					mat2 rotationMatrix = mat2(
						cosTheta, sinTheta,
						-sinTheta, cosTheta
					);
					pos.xy = rotationMatrix * pos.xy;
				}
			#endif

				//compute center pos based on instance
				{
					pos += inPosition;
				}

				interpolatorWorldSpacePos = pos;

				//Camera Space
				pos.xyz -= SceneBuffer.CameraPosition; 
				//Projection
				pos.xy *= SceneBuffer.CameraZoom;
				pos.x /= SceneBuffer.ScreenRatio;

				#if ATLAS_TEXTURE
					if(gl_InstanceID % 3 == 0)
					{
						texCoords.x = 1.f - texCoords.x;
					}

				#endif

				const float ParticleLifeTimeSec = float(`+inDesc.ParticleLifeTimeSec+/* glsl */`);
				float ageNorm = min(0.999f, inAge / ParticleLifeTimeSec);
				interpolatorAge = ageNorm;

				float ageAnimation = min(0.999f, inAge / (ParticleLifeTimeSec * 0.5));
				ageAnimation = pow(ageAnimation, 0.2);
				float TotalFlipFrames = 25.0;
				interpolatorFrameIndex = (ageAnimation * TotalFlipFrames);

				interpolatorTexCoords = texCoords;
			}

			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

		}`;
}


export function GetShaderSourceBasePassSmokeParticlePS(inDesc : ParticleDesc) {
    return /* glsl */ `#version 300 es
		
		precision highp float;
		precision highp sampler2D;
		precision highp samplerCube;

		` +scParticleDefines(inDesc)+/* glsl */`
	
		layout(location = 0) out vec4 OutColor;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		` +scGetBasePassUniformVarDeclPS()+/* glsl */`

		uniform sampler2D AlbedoTexture;
		uniform sampler2D MaskTexture;

		`+scGetBasePassParticleInterpolatorDecl(true)+ /* glsl */`
		flat in float interpolatorAge;
		flat in float interpolatorFrameIndex;
		flat in int interpolatorInstanceId;

		`+scGetCommonFuncsDeclarations()+ /* glsl */`
		`+scGetSphericalFuncsDeclarations()+ /* glsl */`
		`+scPBRDeclarations(true)+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`

		void main()
		{	

			//vec4 color = texture(AlbedoTexture, interpolatorTexCoords.xy).rgba;

			#if ATLAS_TEXTURE
				const vec2 FlipbookSizeRC = vec2(5.0, 5.0);
				vec2 frameSize = 1.f / vec2(5.0, 5.0);
				vec2 uv = interpolatorTexCoords * frameSize;
			
				uint flipBookIndex1D = uint(floor(interpolatorFrameIndex));
				uvec2 FlipBookIndex2D;
				FlipBookIndex2D.x = (flipBookIndex1D % uint(5.0));
				FlipBookIndex2D.y = (flipBookIndex1D / uint(5.0));
			
				uv.x += (frameSize.x * float(FlipBookIndex2D.x));
				uv.y += (frameSize.y * float(FlipBookIndex2D.y));

				vec4 color = texture(AlbedoTexture, uv).rgba;

				#if 1//SMOOTH TRANSITION //TODO:COMPILE TIME CONDITIONAL
					float numFrames = float(FlipbookSizeRC.x * FlipbookSizeRC.y);
					if(ceil(interpolatorFrameIndex) < numFrames)
					{
						flipBookIndex1D = uint(ceil(interpolatorFrameIndex));
						FlipBookIndex2D.x = (flipBookIndex1D % uint(FlipbookSizeRC.x));
						FlipBookIndex2D.y = (flipBookIndex1D / uint(FlipbookSizeRC.x));
						uv = interpolatorTexCoords * frameSize;
						uv.x += (frameSize.x * float(FlipBookIndex2D.x));
						uv.y += (frameSize.y * float(FlipBookIndex2D.y));
						vec4 color2 = texture(AlbedoTexture, uv).rgba;
						color = mix(color, color2, fract(interpolatorFrameIndex));
					}
				#endif

				color.rgb *= color.a;

				//color.rgb *= vec3(1.0, 0.8, 0.5);
				color.rgb *= 0.5;

				color.a *= 0.01;
				//color.a *= mix(0.0, 0.1, hash(uint(interpolatorInstanceId)));

				color *= (1.0 - interpolatorAge);

				//color.a *= (0.25 + color.r * 0.75);

				#if 1
					float radialDistanceScale = length(interpolatorTexCoords - vec2(0.5, 0.5));
					radialDistanceScale = (1.f - clamp(radialDistanceScale, 0.f, 1.f));
					color.a *= radialDistanceScale * radialDistanceScale;
				#endif

			#else
				vec4 color = vec4(vec3(1.0, 0.8, 0.5), 0.1);

				color.rgb *= mix(0.1, 1.0, hash(uint(interpolatorInstanceId)));

				float thresRandParam = 0.5;
				float blurParam = 0.35;
				blurParam = min(thresRandParam - 0.05, blurParam);
				float cutThres = thresRandParam;
				float thres = cutThres - blurParam; // < 0.5 !
				float s = length(interpolatorTexCoords - vec2(0.5, 0.5));
				if(s > cutThres)
				{
					s = 1.0;
				}
				else if(s > thres)
				{
					s = MapToRange(s, thres, cutThres, 0.0, 1.0);
				}	
				else
				{
					s = 0.f;
				}
				s = (1.f - clamp(s, 0.f, 1.f));
				color *= s;
				color.a *= s;

			#endif

			OutColor = color;

			
		}`;
	}
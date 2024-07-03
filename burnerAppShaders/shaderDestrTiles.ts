
import { init } from "@sentry/browser";
import { Vector2 } from "../core/types";
import { ParticleDesc } from "../scenes/payne/paintingMain/particleDesc";
import { scDeclareSceneUniformBuffer,  scGetBasePassInterpolatorDecl,  scGetBasePassUniformTexDeclPS, scGetBasePassUniformVarDeclPS } from "./shaderBasePass";
import { scGetCommonFuncsDeclarations, scGetGIFuncsDeclarations, scGetProceduralCubeVertsDecl, scGetSphericalFuncsDeclarations, scPBRDeclarations } from "./shaderCommon";
import { scGetSDFFuncsDeclarations, scGetSDFShadowDeclarations } from "./shaderSDF";


export function scGetDestrTileCommonFuncsDecl(numTiles : Vector2, tileSize : Vector2)
{
	return /* glsl */`

	const float kTileHealthThres = 0.05;
	
	const ivec2 kNumTiles = ivec2(int(` +
		numTiles.x +
		/* glsl */ `), int(` +
		numTiles.y +
		/* glsl */ `));

	const vec2 kTileSizef = vec2(float(` +
			tileSize.x +
			/* glsl */ `), float(` +
			tileSize.y +
			/* glsl */ `));

		//.z [-1 : 1] for out of bounds check
		ivec3 GetTileIdFromPosWS(vec3 pos)
		{
			ivec3 res = ivec3(0);
			vec2 domainLength = vec2(float(kNumTiles.x) * kTileSizef.x, float(kNumTiles.y) * kTileSizef.y);
			vec2 posNorm;
			vec2 tileIdf;
			posNorm.x = MapToRange(pos.x, -domainLength.x * 0.5, domainLength.x * 0.5, 0.0, 1.0);
			tileIdf.x = posNorm.x * float(kNumTiles.x);
			
			float yOffset = 0.0;
			if (int(tileIdf.x) % 2 == 1)
			{
				yOffset += kTileSizef.y * 0.5;
			}

			posNorm.y = MapToRange(pos.y, -domainLength.y * 0.5 + yOffset, domainLength.y * 0.5 + yOffset, 0.0, 1.0);

			if(posNorm.x < 0.0 || posNorm.y < 0.0 || posNorm.x > 1.0 || posNorm.y > 1.0)
			{
				res.z = -1;
				return res;
			}
			else
			{
				
				tileIdf.y = posNorm.y * float(kNumTiles.y);

				ivec2 id = ivec2(tileIdf);
				res.x = id.x;
				res.y = id.y;
				res.z = 1;
				return res;
			}
		} 
	`;
}


//Collision Check
export function GetProcTileUpdatePS(numTiles : Vector2, tileSize : Vector2
) {
    	return (
        /* glsl */ `#version 300 es
  
		precision highp float;
		precision mediump sampler2D;

		` +scDeclareSceneUniformBuffer()+/* glsl */`
  
		layout(location = 0) out float outHealthState;

		uniform sampler2D TileHealthStateTexture;
  
		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetDestrTileCommonFuncsDecl(numTiles, tileSize)+ /* glsl */`
	  
		void main()
		{
			ivec2 instanceId2D = ivec2(gl_FragCoord.xy);
					
			vec3 posWS = SceneBuffer.IntersectionPosWS;

			ivec3 curTileId = GetTileIdFromPosWS(posWS);

			float res = texelFetch(TileHealthStateTexture, instanceId2D, 0).r;

			if(res > 0.0)
			{	
				//if(payload.x > 0.0)
				if(curTileId.x == instanceId2D.x && curTileId.y == instanceId2D.y)
				{
					const float bulletImpactScale = 0.33 * 0.75;// 0.33;
					res -= bulletImpactScale /* * SceneBuffer.DeltaTime */;
				}

				//Check neighbors
				#if 0
				if((instanceId2D.x > 0 && instanceId2D.x < kNumTiles.x - 1) &&
				(instanceId2D.y > 0 && instanceId2D.y < kNumTiles.y - 1))
				{	if(res > 0.0)
					{
						float resUp = texelFetch(TileHealthStateTexture, instanceId2D + ivec2(0, 1), 0).r;
						float resDown = texelFetch(TileHealthStateTexture, instanceId2D - ivec2(0, 1), 0).r;
						if((resUp + resDown) < 0.1)
						{
							res = 0.0;
						}
					}
				}
				#endif
			}

			outHealthState = res;
		  
  
		}`
    );
}

export function scGetDestrTileBasePassInterpolatorDecl(bPixelShader : boolean)
{
	const inOrOut = bPixelShader ? `in` : `out`;
	return /* glsl */`
		`+inOrOut+/* glsl */` vec2 interpolatorTexCoords;
		`+inOrOut+/* glsl */` vec2 interpolatorTexCoordsLocal;
		`+inOrOut+/* glsl */` vec3 interpolatorWorldSpacePos;
		`+inOrOut+/* glsl */` vec3 interpolatorTBNNormal;
		`+inOrOut+/* glsl */` vec3 interpolatorTBNTangent;
		flat `+inOrOut+/* glsl */` float interpolatorTileHealth;
	`;
}

export function GetShaderSourceBasePassProceduralTileRenderVS(numTiles : Vector2, tileSize : Vector2) {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		uniform sampler2D TileHealthStateTexture;

		`+scGetDestrTileBasePassInterpolatorDecl(false)+ /* glsl */`

		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`

		`+scGetDestrTileCommonFuncsDecl(numTiles, tileSize)+ /* glsl */`
	
		void main()
		{
			const ivec2 kNumTiles = ivec2(int(` +
				numTiles.x +
				/* glsl */ `), int(` +
				numTiles.y +
				/* glsl */ `));

			const vec2 kTileSizef = vec2(float(` +
				tileSize.x +
				/* glsl */ `), float(` +
				tileSize.y +
				/* glsl */ `));


			int instanceIdFlat = int(gl_InstanceID);
			ivec2 instanceId2D = ivec2(instanceIdFlat % kNumTiles.x, instanceIdFlat / kNumTiles.x);

			vec2 domainLength = vec2(float(kNumTiles.x) * kTileSizef.x, float(kNumTiles.y) * kTileSizef.y);
			vec2 startPos = -0.5 * domainLength;
			vec3 pos = vec3(-100.0);
			vec2 texCoords;
			vec3 normal;
			vec3 tangent;
			uint faceId = 0u;

			float health = texelFetch(TileHealthStateTexture, instanceId2D, 0).r;

			interpolatorTileHealth = health;

			//compute local space pos
			if(health > kTileHealthThres)
			{
				////Clockwise, starting from left down
				uint vertId = uint(gl_VertexID);
				pos = CubePositions[vertId];
				normal = CubeNormals[vertId];
				tangent = CubeTangents[vertId];
				texCoords = CubeTexCoords[vertId];

				interpolatorTexCoordsLocal = texCoords + 1.2 * float(gl_InstanceID);

				pos.z *= 0.1 * 0.1;

				const float scale = 0.5;
				pos.x *= kTileSizef.x * scale;
				pos.y *= kTileSizef.y * scale;

				faceId = GetCurCubeFaceIndex(vertId);
			}

			//compute center pos based on instance
			{
				vec2 halfSize = kTileSizef * 0.5;

				vec2 offset = vec2(0.0);
				if(instanceId2D.x % 2 == 1)
				{
					offset.y += halfSize.y;
				}

				vec2 centerPos = startPos + offset + halfSize + vec2(kTileSizef.x * float(instanceId2D.x), kTileSizef.y * float(instanceId2D.y));

				pos.xy += centerPos;
			}

			if(faceId == CUBE_FACE_BACK)
			{
				interpolatorTexCoords.x = MapToRange(pos.x, -domainLength.x * 0.5, domainLength.x * 0.5, 0.0, 1.0);
				interpolatorTexCoords.y = MapToRange(pos.y, -domainLength.y * 0.5, domainLength.y * 0.5, 0.0, 1.0);
				interpolatorTexCoords *= 2.0;
			}
			else
			{
				interpolatorTexCoords = texCoords * 0.1 * vec2(0.5, 1.0);
			}
			

			interpolatorWorldSpacePos = pos;

			//Camera Space
			pos.xyz -= SceneBuffer.CameraPosition; 
			//Projection
			pos.xy *= SceneBuffer.CameraZoom;
			pos.x /= SceneBuffer.ScreenRatio;
			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

			interpolatorTBNNormal = normal;
			interpolatorTBNTangent = tangent;

		}`;
}

export function scGetDestrTileBasePassUniformTexDeclPS()
{
	return /* glsl */`
		uniform sampler2D AlbedoTexture;
		uniform sampler2D NormalTexture;
		uniform sampler2D RoughnessTexture;
		uniform sampler2D MaskTexture;
		uniform sampler2D TileHealthStateTexture;
	`;
}


export function scGetBasePassLightShading()
{
	return (/* glsl */`
	if(bEmitter == 0)
		{
			s.LightColor = SceneBuffer.LightColor;
			s.LightRadius = SceneBuffer.PointLightRadius1;

			/* {
				s.LightPos = SceneBuffer.PointLightPos1;
				light += CalculateLightPBR(s);
			} */

			//Spot Light
			#if 1
			{
				s.LightPos = SceneBuffer.SpotlightPos2;
				vec3 lightTargetPos = s.LightPos + SceneBuffer.SpotlightTargetPos2;
				float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, normalize(lightTargetPos - s.LightPos), s.PixelPos, 1.0, 0.75);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				//lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				//SpotlightFalloff = 1.0;
				light += CalculateLightPBR(s) * SpotlightFalloff;

				{
					//emitted light
					s.LightPos.z = -0.85;
					s.LightPos.y += 0.5;
					s.LightColor *= 0.2;
					light += CalculateLightPBR(s);
				}
			}
			#endif

		}
	`);
}

export function scGetProcTileSpecificPSCode()
{
	return (/* glsl */`

	//damage texture
	{
		vec3 MaskColorTex = texture(MaskTexture, interpolatorTexCoordsLocal.xy * 0.5).rgb;

		if(interpolatorTileHealth < 1.0)
		{
			//cracks
			vec2 U = interpolatorTexCoordsLocal;
			vec2 V = U * vec2(1.0);
			const float CRACK_zebra_scale = 1.5;
			const float CRACK_zebra_amp = 1.7;
			 vec2 D = fbm22(CRACK_zebra_scale*U) / CRACK_zebra_scale / CRACK_zebra_amp;
			vec3  H = voronoiB( V + D );
			float d = H.x;                                     // distance to cracks
			  float r = voronoi(V).x;                            // distance to center
			const float thickness = 1.0;
			const float CRACK_slope = 2.0;
			float CRACK_width = mix(0.015, 0.0, interpolatorTileHealth * interpolatorTileHealth);
			const float CRACK_profile = 0.15;
			d = min(1.0, CRACK_slope * pow(max(0.0, d-CRACK_width), CRACK_profile));

			d *= thickness;
			d = clamp(d, 0.0, 1.0);
			float dCached = d;
			//s.Albedo = vec3(1.0);
			//s.Albedo *= d;
			d = 1.0 - d;
			d *= d * d;

			s.Albedo = mix(s.Albedo, MaskColorTex * vec3(1.0, 0.8, 0.5), d * 0.4);
			s.Roughness = mix(s.Roughness, 1.0, d);
		}
	}

	s.Albedo *= 1.2;


	#if 1 //SSR
	{
		vec3 sampleRayDir = normalize(reflect(-s.vToCam, s.Normal));
		float ssrScale = dot(sampleRayDir, vec3(0.0, 0.0, 1.0));	
		vec3 rayOrigin = interpolatorWorldSpacePos;
		float a = s.Roughness;
		float a2 = a * a;
		
		if(ssrScale > 0.0)
		{
			light += ssrScale * 0.5;

		}

	}
	#endif

	`);
}


export function GetShaderSourceBasePassDestrRenderPS(numTiles : Vector2, tileSize : Vector2, inSC_CUSTOM : string = ``) {
    return /* glsl */ `#version 300 es
		
		precision highp float;
		precision highp sampler2D;
		precision highp samplerCube;

		layout(location = 0) out vec3 OutColor;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		` +scGetBasePassUniformVarDeclPS()+/* glsl */`
		` +scGetDestrTileBasePassUniformTexDeclPS()+/* glsl */`

		`+scGetDestrTileBasePassInterpolatorDecl(true)+ /* glsl */`

		`+scGetCommonFuncsDeclarations()+ /* glsl */`
		`+scGetSphericalFuncsDeclarations()+ /* glsl */`
		`+scPBRDeclarations(true)+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`

		`+scGetDestrTileCommonFuncsDecl(numTiles, tileSize)+ /* glsl */`

		void main()
		{	

			vec3 light = vec3(0.0);

			PBRSettings s;

			//Albedo
			s.Albedo = texture(AlbedoTexture, interpolatorTexCoords.xy).rgb;
			

			//Normal
			s.Normal = texture(NormalTexture, interpolatorTexCoords.xy).rgb;
			s.Normal = DecodeNormalTexture(s.Normal, 0.8);
			{
				vec3 TBNNormal = normalize(interpolatorTBNNormal)/*  * -1.0 */;
				vec3 TBNTangent = normalize(interpolatorTBNTangent);
				vec3 BiTangent = cross(TBNNormal, TBNTangent);
				mat3 TangentToWorldMatrix = mat3(
				(TBNTangent),
				normalize(BiTangent),
				(TBNNormal)
				);
				s.Normal = normalize(TangentToWorldMatrix * s.Normal);
			#if 0//DEBUG NORMALS
				vec3 c = max(vec3(0.0), s.Normal) * 0.85;
				OutColor = c;
				return;
			#endif
			}

			//Roughness
			s.Roughness = max(0.0, texture(RoughnessTexture, interpolatorTexCoords.xy).r);
			s.Roughness = min(1.0, s.Roughness + 0.05);

			s.Metalness = Metalness;

			s.PixelPos = interpolatorWorldSpacePos;
			s.vToCam = normalize(SceneBuffer.CameraPosition.xyz - interpolatorWorldSpacePos);


			
			`+inSC_CUSTOM+ /* glsl */`



			vec3 indirectDiffuse = SceneBuffer.LightColor;
			indirectDiffuse = vec3(0.5);
			light += vec3(indirectDiffuse * s.Albedo / PI * 1.0) * (max(-0.3, SceneBuffer.AmbientLight) + 0.3);

			`+scGetBasePassLightShading()+ /* glsl */`


			//OutColor = s.Albedo; return;
			OutColor = light;

			
		}`;
	}



//=============================================================================================================================
// 																											DECAL
//=============================================================================================================================

export function GetDestDecalUpdateShaderVS() {
	return (
	/* glsl */ `#version 300 es

	precision highp float;
	precision mediump sampler2D;

	layout(location = 0) in vec3 inPosition;

	out vec3 outPosition;

	uniform ivec2 ParticleCurSpawnIndexRange; //Particles withing this range are force-respawned
	
	` +scDeclareSceneUniformBuffer()+/* glsl */`

	`+scGetCommonFuncsDeclarations()+ /* glsl */`
		
	vec3 GetSpawnPos()
	{
		//From input pos
		vec2 posNDC = SceneBuffer.UserInputPosNDC.xy;
		//vec3 spawnPos = TransformNDCToWorld(posNDC, -SceneBuffer.CameraPosition.z);
		vec3 spawnPos = SceneBuffer.IntersectionPosWS;
		//spawnPos.z = min(-0.01, spawnPos.z);
		return spawnPos;
	}

	void main()
	{
		int particleId = int(gl_VertexID);
		if(particleId >= ParticleCurSpawnIndexRange.x && particleId < ParticleCurSpawnIndexRange.y)
		{
			outPosition = GetSpawnPos();
		}
		else
		{
			outPosition = inPosition;
		}
	}`
	);
}


export function GetShaderSourceBasePassDecalRenderVS() {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;

		layout(location = 0) in vec3 inPosition;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		uniform sampler2D TileHealthStateTexture;

		`+scGetDestrTileBasePassInterpolatorDecl(false)+ /* glsl */`

		`+scGetCommonFuncsDeclarations()+ /* glsl */`

		`+scGetProceduralCubeVertsDecl()+ /* glsl */`
	
		void main()
		{
			
			vec3 pos = vec3(-100.0);
			vec2 texCoords;
			vec3 normal = vec3(0.0, 0.0, -1.0);
			vec3 tangent = vec3(1.0, 0.0, 0.0);

			//compute local space pos
			{
				////Clockwise, starting from left down
				uint vertId = uint(gl_VertexID);
				pos = CubePositions[vertId];
				//normal = CubeNormals[vertId];
				texCoords = CubeTexCoords[vertId];

				pos.z = -0.01;
				//pos.z = 0.0;

				const float scale = 0.05;
				pos.x *= scale;
				pos.y *= scale;

			}

			//Per Instance Scale and Orientation
			{
				float scale = MapToRange(hash(uint(gl_InstanceID)), 0.0, 1.0, 0.75, 1.25);
				pos.xy *= scale;


				float angle = hash(uint(gl_InstanceID) * 3u + uint(13.2)) * PI * 2.0;

				float cosTheta = cos(angle);
				float sinTheta = sin(angle);
				mat2 rotationMatrix = mat2(
					cosTheta, sinTheta,
					-sinTheta, cosTheta
				);
				pos.xy = rotationMatrix * pos.xy;
				tangent.xy = rotationMatrix * tangent.xy;

			}

			//compute center pos based on instance
			{
				pos.xy += inPosition.xy;
				//pos.z = inPosition.z;
			}

			interpolatorWorldSpacePos = pos;

			//Camera Space
			pos.xyz -= SceneBuffer.CameraPosition; 
			//Projection
			pos.xy *= SceneBuffer.CameraZoom;
			pos.x /= SceneBuffer.ScreenRatio;
			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

			interpolatorTBNNormal = normal;
			interpolatorTBNTangent = tangent;

			interpolatorTexCoords = texCoords;

		}`;
}


export function scGetDecalSpecificPSCode()
{
	return (/* glsl */`

	if(bHasMask > 0)
	{
		float mask = texture(MaskTexture, interpolatorTexCoords.xy).r;
		if(mask < 0.5)
		{
			discard;
		}
	}
	
	ivec3 curTileId = GetTileIdFromPosWS(interpolatorWorldSpacePos);
	if(curTileId.z < 0)
	{
		discard;
	}
	else
	{
		float health = texelFetch(TileHealthStateTexture, curTileId.xy, 0).r;
		if(health < kTileHealthThres)
		{
			discard;
		}
	}

	s.Albedo *= vec3(1.0, 0.7, 0.5);

	`);
}


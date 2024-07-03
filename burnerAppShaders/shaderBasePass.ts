import { GSettings } from "../core/settings";
import { Vector2 } from "../core/types";
import { SceneBase, SceneUniformBufferCPU } from "../scenes/scene";
import { scGetCommonFuncsDeclarations, scGetGIFuncsDeclarations, scGetSphericalFuncsDeclarations, scPBRDeclarations } from "./shaderCommon";
import { scGetSDFFuncsDeclarations, scGetSDFShadowDeclarations } from "./shaderSDF";

export function scDeclareSceneUniformBuffer()
{
	
	const decl =  /* glsl */`
	layout(std140) uniform SceneUniformBuffer {
		`+SceneBase.UniformBuffer.Params.GenerateParamsDeclarations()+/* glsl */`
	} SceneBuffer;

	`
	//console.log(decl);
	return decl;
	
	;
}

export function scGetBasePassInterpolatorDecl(bPixelShader : boolean)
{
	const inOrOut = bPixelShader ? `in` : `out`;
	return /* glsl */`
		`+inOrOut+/* glsl */` vec2 interpolatorTexCoords;
		`+inOrOut+/* glsl */` vec3 interpolatorWorldSpacePos;
		`+inOrOut+/* glsl */` vec3 interpolatorProjectorSpacePos;
		`+inOrOut+/* glsl */` vec3 interpolatorTBNNormal;
		`+inOrOut+/* glsl */` vec3 interpolatorTBNTangent;
		flat `+inOrOut+/* glsl */` float interpolatorTangentSign;
	`;
}

export function GetShaderSourceBasePassMainVS() {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;
	
		layout(location = 0) in vec3 VertexBuffer;
		layout(location = 1) in vec2 TexCoordsBuffer;
		layout(location = 2) in vec3 NormalsBuffer;
		layout(location = 3) in vec4 TangentsBuffer;

		` +scDeclareSceneUniformBuffer()+/* glsl */`


		uniform mat4 LocalToWorldMatrix;
		uniform vec3 PreTranslate;
		uniform mat4 WorldToViewMatrix;
		uniform float TextureScale;

		`+scGetBasePassInterpolatorDecl(false)+ /* glsl */`

		//out vec3 interpolatorColor;

		`+scGetCommonFuncsDeclarations()+ /* glsl */`
	
		void main()
		{
			vec3 pos = VertexBuffer.xyz;
			pos -= PreTranslate;

			pos.xyz = vec4(vec4(pos, 1.0) * LocalToWorldMatrix).xyz;

			mat4x4 rotMat = transpose(inverse(LocalToWorldMatrix));
			interpolatorTBNNormal = vec4(vec4(NormalsBuffer, 0.0) * rotMat).xyz;
			interpolatorTBNTangent = vec4(vec4(TangentsBuffer.xyz, 0.0) * rotMat).xyz;
			interpolatorTangentSign = TangentsBuffer.w;

			interpolatorWorldSpacePos = pos;

			//interpolatorProjectorSpacePos = pos - SceneBuffer.ProjectorPos + vec3(0.0, 0.0, 1.0);
			/* vec3 projPos = SceneBuffer.WorldToProjectorSpaceMat[3].xyz;
			interpolatorProjectorSpacePos = pos + projPos + vec3(0.0, 0.0, 1.0); */
			vec3 projSpace = (vec4(pos, 1.0) * SceneBuffer.WorldToProjectorSpaceMat).xyz;
			projSpace.xy *= SceneBuffer.ProjectorZoom;
			//projSpace.z /= SceneBuffer.ProjectorZoom;
			//projSpace.z += 1.0;
			projSpace.z += 0.01;
			interpolatorProjectorSpacePos = projSpace;
			

			//Camera Space
			pos.xyz = (vec4(pos, 1.0) * WorldToViewMatrix).xyz;
			//pos.xyz -= SceneBuffer.CameraPosition; 

			//Projection
			pos.xy *= SceneBuffer.CameraZoom;
			pos.x /= SceneBuffer.ScreenRatio;

			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

			interpolatorTexCoords = TexCoordsBuffer.xy * TextureScale;

			#if 0//Visualize triangles
			{
				float rnd = hash((uint(gl_VertexID) / 3u));
				vec3 hsv = vec3(rnd, 1.0, 1.0);
				vec3 c = hsv2rgb(hsv);
				interpolatorColor = c;
			}
			#endif

		}`;
}

function scGetBasePassLightShading()
{
	return (/* glsl */`
	if(bEmitter == 0)
		{
			s.LightColor = SceneBuffer.LightColor;
			s.LightRadius = SceneBuffer.PointLightRadius1;

			#ifdef BASEPASS_SHADOW
			vec3 shadowPixelPos = s.PixelPos;
			//shadowPixelPos = s.PixelPos + s.Normal * 1.01;
			#endif

			//Light1
		#ifndef LAMP_LIGHTS
			{
				s.LightPos = SceneBuffer.PointLightPos1;
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				light += CalculateLightPBR(s) * lightShadow;
			}
		#endif//LAMP_LIGHTS

			//Light2
			#if 1
			{
				//s.LightRadius = SceneBuffer.PointLightRadius1 * 1.5;
				s.LightPos = SceneBuffer.PointLightPos2;
				s.LightColor = mix(SceneBuffer.LightColor2, SceneBuffer.LightColor, 0.25);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				//lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				light += CalculateLightPBR(s) * 0.25 * lightShadow * AO;
			}
			#endif
		

			//Spot Light
			#if 1
			{
				s.LightPos = SceneBuffer.SpotlightPos2;
				vec3 lightTargetPos = SceneBuffer.SpotlightTargetPos2 - vec3(0.0, 0.5, 0.0);
				float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, normalize(lightTargetPos - s.LightPos), s.PixelPos, 1.0, 0.75);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				//SpotlightFalloff = 1.0;
				light += CalculateLightPBR(s) * 1.5 * lightShadow * SpotlightFalloff * AO;
			}
			#endif

			//Flashlight
			#ifndef LAMP_LIGHTS
			#if 1
			{
				s.LightColor = mix(vec3(1.0), SceneBuffer.LightColor2, 0.25);
				
				//s.LightPos = SceneBuffer.SpotlightTargetPos;
				vec3 lightTargetPos = SceneBuffer.SpotlightTargetPos;
				s.LightRadius = 30.0;
				#if 1
				s.LightPos = SceneBuffer.SpotlightPos;
				//s.LightPos = SceneBuffer.CameraPosition;
				//s.LightColor *= 2.0;
				#else
				const float spotOffset = 0.5;
				s.LightPos.xy = SceneBuffer.SpotlightTargetPos.xy * -spotOffset;
				s.LightPos.y += SceneBuffer.SpotlightPos.y;
				//s.LightPos.xy = -SceneBuffer.SpotlightTargetPos.xy;
				s.LightPos.z = SceneBuffer.SpotlightPos.z;
				#endif
				
				float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, normalize(lightTargetPos - s.LightPos), s.PixelPos, SceneBuffer.SpotlightOuterAngle, SceneBuffer.SpotlightInnerAngleOffset);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				if(SpotlightFalloff > 0.001)
				{
					lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
				}
			#endif//BASEPASS_SHADOW
				light += CalculateLightPBR(s) * 2.0 * lightShadow * SpotlightFalloff;
			}
			#endif
			#endif//LAMP_LIGHTS



			//Lamp Spotlights
			#ifdef LAMP_LIGHTS
			//Spot Light left
			{
				s.LightRadius = SceneBuffer.LampLightRadius;
				s.LightColor = SceneBuffer.LightColor * SceneBuffer.LampLightIntensity;
				s.LightPos = SceneBuffer.SpotlightPos3;
				float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, vec3(0.0, -1.0, 0.0), s.PixelPos, 1.25, 0.2);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				//lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				//SpotlightFalloff = 1.0;
				light += CalculateLightPBR(s) * 1.5 * lightShadow * SpotlightFalloff * AO;

				{
					//emitted light
					s.LightPos.z = -0.85;
					s.LightPos.y += 0.3;
					s.LightColor *= 0.6;
					light += CalculateLightPBR(s) * AO;
				}
			}
			//Spot Light Right
			{
				s.LightColor = SceneBuffer.LightColor;
				s.LightColor = mix(vec3(0.8), s.LightColor, 0.75) * SceneBuffer.LampLightIntensity;
				s.LightPos = SceneBuffer.SpotlightPos3;
				s.LightPos.x *= -1.0;
				s.LightPos.x += 0.85;
				float SpotlightFalloff = GetSpotlightFalloff(s.LightPos, vec3(0.0, -1.0, 0.0), s.PixelPos, 1.25, 0.2);
				float lightShadow = 1.0;
			#ifdef BASEPASS_SHADOW
				//lightShadow = SDFComputeShadow(shadowPixelPos, s.LightPos);
			#endif //BASEPASS_SHADOW
				//SpotlightFalloff = 1.0;
				/* float rndHash = hashf(mod(SceneBuffer.Time * 0.005, 10.0));
				float rnd = (sin(SceneBuffer.Time) + 1.0) * 0.5;
				const float flickThres = 0.06;
				if(!(rnd > (0.7 - flickThres) && rnd < 0.7))
				{
					rndHash = 1.0;
				}
				rnd = max(0.85, rnd);
				rndHash = max(0.7, rndHash);
				rnd *= rndHash; */
				float rnd = SceneBuffer.LampFlickerScale;
				light += CalculateLightPBR(s) * 1.5 * lightShadow * SpotlightFalloff * AO * rnd;
				{
					//emitted light
					s.LightPos.z = -0.85;
					s.LightPos.y += 0.3;
					s.LightColor *= 0.6;
					light += CalculateLightPBR(s) * AO * rnd;
				}
			}
			#endif

		}
	`);
}


export function GetShaderSourceBasePassAlphaRenderMainPS() {
	return /* glsl */ `#version 300 es

	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float OutColor;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`

	uniform sampler2D AlbedoTexture;

	in vec2 interpolatorTexCoords;
	in vec3 interpolatorWorldSpacePos;
	in vec3 interpolatorProjectorSpacePos;
	in vec3 interpolatorTBNNormal;
	in vec3 interpolatorTBNTangent;
	flat in float interpolatorTangentSign;

	void main()
	{

		vec4 c = texture(AlbedoTexture, interpolatorTexCoords.xy);
		c.a = MapToRange(c.a, 0.0, 1.0, 0.2, 1.0);
		c.r *= c.a;
		OutColor = c.r;
		//OutColor = vec4(light, 1.0);
	}`;
}

export function GetShaderSourceBasePassDensityRenderMainPS(numCells : Vector2) {
	return /* glsl */ `#version 300 es

	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 OutColor;
	layout(location = 1) out float OutAlpha;
	layout(location = 2) out float OutMask;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`

	uniform sampler2D AlbedoTexture;
	uniform sampler2D NormalTexture;
	uniform sampler2D PrevFrameSceneColor;
	uniform sampler2D PrevFrameSceneColorBlur;

	in vec2 interpolatorTexCoords;
	in vec3 interpolatorWorldSpacePos;
	in vec3 interpolatorProjectorSpacePos;
	in vec3 interpolatorTBNNormal;
	in vec3 interpolatorTBNTangent;
	flat in float interpolatorTangentSign;

	void main()
	{
	#if 1
		vec2 coords = vec2(gl_FragCoord.xy) / vec2(600.0);
		
		float scale = 1.0;
		{
			const float scaleMin = 0.75;
			const float scaleMax = 1.25;
			float t = mod(SceneBuffer.Time * 0.1, 2.f);
			if(t < 1.f)
			{
				scale = mix(scaleMin, scaleMax, t);
			}
			else
			{
				scale = mix(scaleMax, scaleMin, t - 1.f);
			}
		}
		coords *= 2.0;
		float a = SceneBuffer.Time * 0.01;
		coords.x = coords.x * cos(a) + coords.y * -sin(a);
		coords.y = coords.x * sin(a) + coords.y * cos(a);
		vec3 c = texture(AlbedoTexture, coords.xy).rgb;
		c *= 1.75;
		//c = Contrast(c, 0.2);

		float gradMask = 1.0;
		{
			coords = vec2(gl_FragCoord.xy) / vec2(300.0);
			coords.x += SceneBuffer.Time * 0.017;
			coords.y += SceneBuffer.Time * 0.017 * 0.3;
			gradMask = texture(NormalTexture, coords.xy).r;
		}
		vec3 final = vec3(0.0);
		//final += gradMask * vec3(1.0, 0.0, 0.0);
		{
			final.r = max(c.r, final.r);
			final.gb += c.gb;
		}
		//final *= vec3(1.0, 0.2, 0.1);
		OutColor = final;
		OutAlpha = 0.75;
		OutMask = 1.0;
	#else

		const vec2 kNumCells = vec2(float(` +
		numCells.x +
		/* glsl */ `), float(` +
		numCells.y +
		/* glsl */ `));

		vec2 coords = vec2(gl_FragCoord.xy) / kNumCells;

		vec3 c = texture(PrevFrameSceneColor, coords.xy).rgb;

		//OutColor = c;
		OutColor = vec3(1.0);
		OutAlpha = 1.0;
		OutMask = 1.0;

	#endif
	}`;
}

export function scGetBasePassUniformVarDeclPS()
{
	return /* glsl */`
		uniform float Brightness;
		uniform float Metalness;
		uniform int EUseLightprobe;//0-disabled, 1-default, 2-debug
		uniform int LightprobeIndex;
		uniform int LightprobesGridSize;
		uniform int bEmitter;
		uniform int bHasMask;
		uniform mat4 WorldToViewMatrix;
	`;
}

export function scGetBasePassUniformTexDeclPS()
{
	return /* glsl */`
		uniform sampler2D AlbedoTexture;
		uniform sampler2D NormalTexture;
		uniform sampler2D RoughnessTexture;
		uniform sampler2D MaskTexture;

		uniform sampler2D BRDF_LUT;

		uniform sampler2D ProjectorTexture;
		uniform sampler2D PrevFrameSceneColor;
		uniform sampler2D PrevFrameSceneColorBlur;
		uniform samplerCube CubemapTexture;
		uniform sampler2D SphericalTexture;
		uniform sampler2D SHCoefficientsTexture;
		uniform highp sampler3D SDFTexture;
		uniform highp sampler3D SDFTexture2;
		uniform sampler2D RopePositionsTextureSDF;
	`;
}





export function GetShaderSourceBasePassMainPS(sceneSDFFunc : string) {
    return /* glsl */ `#version 300 es
		
		precision highp float;
		precision highp sampler2D;
		precision highp samplerCube;
	
		layout(location = 0) out vec3 OutColor;
	
		` +scDeclareSceneUniformBuffer()+/* glsl */`

		` +scGetBasePassUniformVarDeclPS()+/* glsl */`
		` +scGetBasePassUniformTexDeclPS()+/* glsl */`

		`+scGetBasePassInterpolatorDecl(true)+ /* glsl */`

		#define BASEPASS_SHADOW
		
		`+scGetCommonFuncsDeclarations()+ /* glsl */`
		`+scGetSphericalFuncsDeclarations()+ /* glsl */`
		`+scPBRDeclarations()+ /* glsl */`
		`+scGetGIFuncsDeclarations()+ /* glsl */`
		`+scGetSDFFuncsDeclarations(sceneSDFFunc)+ /* glsl */`
		`+scGetSDFShadowDeclarations()+ /* glsl */`

		float GetCamDistFade()
		{
			float scale = 1.0;
			#if 1
			float d = length(SceneBuffer.CameraPosition.xyz - interpolatorWorldSpacePos);
			const float dThres = 2.5;
			const float shadow = 0.2;
			if(d < dThres)
			{
				float brightScale = clamp(MapToRange(d, 0.0, dThres, shadow, 1.0), shadow, 1.0);
				brightScale *= brightScale;
				scale *= brightScale;
			}
			#endif
			return scale;
		}

		void main()
		{
			vec3 light = vec3(0.0);
			PBRSettings s;
			
			/////////////
			//	ALBEDO
			////////////
			s.Albedo = texture(AlbedoTexture, interpolatorTexCoords.xy).rgb;

			if(bEmitter < 1)
			{
				s.Albedo += vec3(0.0, 0.0, 0.01);
				//s.Albedo = vec3(0.05, 0.16, 0.64);
				//s.Albedo = vec3(0.74, 0.05, 0.05);
				//s.Albedo = vec3(0.05, 0.64, 0.05);
				//s.Albedo = vec3(0.64, 0.64, 0.05);
			}
			//OutColor = s.Albedo; return;

			if(bEmitter > 0)
			{
				float scale = 1.0;
				#if 1 //EDGE FADE-OUT
					const float edgeShadow = 0.75;
					const highp float rectPow = 16.f;
					vec2 ndcPos = MapToRange(interpolatorTexCoords, 0.0, 1.0, -1.0, 1.0);
					highp float rectCircleLength = pow(abs(ndcPos.x), rectPow) + pow(abs(ndcPos.y), rectPow);
					const float thres = 0.1;
					if(rectCircleLength < thres)
					{
						scale = 1.0;
					}
					else
					{
						scale = MapToRange(rectCircleLength, thres, 2.0, 1.0, edgeShadow);
						//scale = smoothstep(0.0, 1.0, scale);
					}
					const float edgeThres = 0.99;

					if(abs(ndcPos.x) > edgeThres)
					{
						scale *= MapToRange(abs(ndcPos.x), edgeThres, 1.0, 1.0, edgeShadow);
					}
					if(abs(ndcPos.y) > edgeThres)
					{
						scale *= MapToRange(abs(ndcPos.y), edgeThres, 1.0, 1.0, edgeShadow);
					}
					scale *= scale;
					scale = min(1.0, scale + 0.1);
				#endif
				//Fade objects close to camera
				scale *= GetCamDistFade();

			#if 1//EMITTER IS SHADED
				light += s.Albedo * SceneBuffer.EmmissiveColor * Brightness * scale;
			#else
				OutColor = s.Albedo * SceneBuffer.EmmissiveColor * Brightness * scale;
				return;
			#endif
			}
			else
			{
				//s.Albedo = vec3(0.0, 0.0, 1.0);
				//s.Albedo *= vec3(0.5, 0.5, 1.0);
			}
		
			/////////////
			//	NORMAL
			////////////
			s.Normal = texture(NormalTexture, interpolatorTexCoords.xy).rgb;
			const float normalHarshness = 1.0;
			s.Normal = DecodeNormalTexture(s.Normal, normalHarshness);
			//s.Normal = vec3(0.0, 0.0, 1.0);
			{
				vec3 TBNNormal = normalize(interpolatorTBNNormal)/*  * -1.0 */;
				vec3 TBNTangent = normalize(interpolatorTBNTangent);
				vec3 BiTangent = cross(TBNNormal, TBNTangent) * interpolatorTangentSign;
				mat3 TangentToWorldMatrix = mat3(
				(TBNTangent),
				normalize(BiTangent),
				(TBNNormal)
				);
				s.Normal = normalize(TangentToWorldMatrix * s.Normal);
			#if 0//DEBUG NORMALS
				vec3 c = max(vec3(0.0), s.Normal) * 0.85;
				//vec3 c = max(vec3(0.0), TBNTangent) * 0.85;
				OutColor = c;
				return;
			#endif
			}
			
			////////////////
			//	ROUGHNESS
			///////////////
			const float roughnessMin = 0.05; 
			#if 1 //RGH from texture
			s.Roughness = max(0.0, texture(RoughnessTexture, interpolatorTexCoords.xy).r);
			s.Roughness = min(1.0, s.Roughness + roughnessMin);
			#else
			s.Roughness = SceneBuffer.UniformRoughness + roughnessMin;
			#endif
			//s.Roughness *= s.Roughness; //Filament rgh remap
			//OutColor = vec3(s.Roughness); return;

			//s.Albedo = vec3(s.Roughness);

			#if 1
			s.Metalness = Metalness;
			#else
			s.Metalness = SceneBuffer.UniformMetalness;
			#endif


			/////////////////////////
			//
			//	  	SHADING
			//
			////////////////////////
			

			s.PixelPos = interpolatorWorldSpacePos;
			s.vToCam = normalize(SceneBuffer.CameraPosition.xyz - interpolatorWorldSpacePos);


			/////////////////////////
			//	  INDIRECT LIGHT
			////////////////////////
			//Ambient Occlusion
			float AO = 1.0;
			#if 1
			if(bEmitter == 0)
			{
				float NDotV = max(dot(s.Normal, s.vToCam), 0.0);
				vec3 F0 = vec3(SceneBuffer.FresnelMin);
				F0 = mix(F0, s.Albedo, s.Metalness);
				vec3 F = FresnelSchlickRoughness(NDotV, F0, s.Roughness);
				vec3 kS = F;
				vec3 kD = 1.0 - kS;

				
				#if 1
				{
					//try using harsher normal
					AO *= CalcAOPrecise( interpolatorWorldSpacePos, s.Normal );
					float ran = hash( uint(gl_FragCoord.x) + 1920U* uint(gl_FragCoord.y) + (1920U*1080U) );  
					AO *= CalcAOSpread( interpolatorWorldSpacePos, s.Normal, ran );
					AO *= CalcAOSpreadHemisphere( interpolatorWorldSpacePos, s.Normal, 1., 0.5 );
					//AO *= AO;
					//OutColor = vec3(AO * 0.5); return;
				}
				#endif


				//Ambient Diffuse
				{
					#if 1
						vec3 samplePos = s.Normal;
						samplePos.y *= -1.0;
					
						//SH Basis Coefficients
						if(EUseLightprobe > 0)
						{
							if(EUseLightprobe == 1)
							{
								vec3 indirectDiffuse = GetIrradianceFromLightprobes(samplePos, interpolatorWorldSpacePos);
								indirectDiffuse = max(vec3(0.01), indirectDiffuse);
								indirectDiffuse *= vec3(0.8, 0.9, 1.0);
								//indirectDiffuse *= 10.0;
								//indirectDiffuse *= vec3(0.25, 0.5, 1.0) * 2.0;
								light += vec3(indirectDiffuse * s.Albedo / PI * kD) * AO * (1.0 + SceneBuffer.AmbientLight + 0.25);
								//light += vec3(vec3(0.5) * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight);
								//light += vec3(vec3(0.0, 0.0, 1.0) * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight);

							}
							else//TODO: Remove, used for GI Probes visualizer
							{
								light = SampleDebugLightprobe(samplePos);
								OutColor = light; return;
							}
						}
					#else //Regular diffuse
					{
						vec3 indirectDiffuse = SceneBuffer.LightColor;
						indirectDiffuse = vec3(0.5);
						light += vec3(indirectDiffuse * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight + 0.3);
					}

					#endif//Indirect diffuse
				}


				//Ambient Specular

				vec3 sampleRayDir = normalize(reflect(-s.vToCam, s.Normal));
				float ssrScale = dot(sampleRayDir, vec3(0.0, 0.0, 1.0));
				#if 0 //Reflection probe
				if(ssrScale < 0.0)
				{
					vec3 samplePos = sampleRayDir;
					samplePos.y *= -1.0;
					{
						float angle = SceneBuffer.EnvMapRotationAngle;
						mat3 rotationMatrix = mat3(
							cos(angle), 0.0, sin(angle),
							0.0,  1.0, 0.0,
						   -1.0*sin(angle), 0.0, cos(angle)
						);
						samplePos = rotationMatrix * samplePos;
					}

					//get relfected color
					int numMips = 10;
					float envMip = s.Roughness * float(numMips - 1);
					vec3 prefilteredColor = textureLod(CubemapTexture, samplePos, envMip).rgb;
					prefilteredColor.g *= 0.75;
					vec2 envBRDF = texture(BRDF_LUT, vec2(NDotV, (s.Roughness))).rg;
					vec3 fScale = vec3(F * envBRDF.x + envBRDF.y);
					vec3 indirectSpecular = prefilteredColor * fScale;
					light += indirectSpecular * AO;
					//OutColor = light; return;
				}
				#endif

				#if 1//SSR
				{
					
					vec3 rayOrigin = interpolatorWorldSpacePos;
					float a = s.Roughness;
					float a2 = a * a;
					
					if(ssrScale > 0.0)
					{
						//OutColor = vec3(ssrScale); return; 

				#if 0 //Importance sample rough surfaces
					//get low-discrepancy random value
					const uint SAMPLE_COUNT = 16u;
					float totalWeight = 0.0;   
    				vec3 prefilteredColor = vec3(0.0);  
					uint coordX = uint(gl_FragCoord.x);
					uint coordY = uint(gl_FragCoord.y);
					for(uint i = 0u; i < SAMPLE_COUNT; ++i)
					{
						vec2 Xi = Hammersley(i, SAMPLE_COUNT);
						//vec2 Xi = vec2(hashf(float(i + coordX * coordY)), hash(i));
						vec3 curLightDir = ImportanceSampleGGX(Xi, sampleRayDir, s.Roughness);
						curLightDir = normalize(curLightDir);
						rayOrigin = interpolatorWorldSpacePos + sampleRayDir * 0.01;
						float NdotL = max(dot(s.Normal, curLightDir), 0.0);
						SDFMarchResult sdfRes = MarchSDFScene(rayOrigin, curLightDir);
						if(sdfRes.bHit)
						{
							vec3 posNDC = (vec4(sdfRes.HitPosWS, 1.0) * WorldToViewMatrix).xyz;
							posNDC.xy *= SceneBuffer.CameraZoom;
							posNDC.x /= SceneBuffer.ScreenRatio;
							posNDC.xy /= (0.1f + posNDC.z);
							vec2 sceneUV = (posNDC.xy + 1.0) * 0.5;
							
							vec3 h = normalize(s.vToCam + curLightDir);
							float NdotH  = max(dot(s.Normal, h), 0.0);
							float HdotV  = max(dot(s.vToCam, h), 0.0);
							float NDF = NormalDistributionGGXFast(NdotH, a2);  
							float pdf = (NDF * NdotH / (4.0 * HdotV + 0.0001)); 

							vec3 prevSceneColor = textureLod(PrevFrameSceneColor, sceneUV.xy, 0.0).rgb;
							prevSceneColor = mix(prevSceneColor, textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb, s.Roughness * 0.75);
							//vec3 prevSceneColor = textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb;
							prefilteredColor += prevSceneColor * NdotL * pdf;
							totalWeight += NdotL * pdf;
						}
					}
					prefilteredColor = prefilteredColor / (totalWeight + 0.0001);
					vec2 envBRDF = texture(BRDF_LUT, vec2(NDotV, (s.Roughness))).rg;
					vec3 fresnelCustom = FresnelSchlickRoughness(NDotV, vec3(0.5), s.Roughness);
					vec3 fScale = vec3(fresnelCustom * envBRDF.x + envBRDF.y);
					if(!any(isnan(prefilteredColor)))
					{
						light += max(vec3(0.0), prefilteredColor) * fScale * AO;
					}

					//
					/////////////TODO DO SMTH ABOUT THIS, DISABLE SAMPLES OUTSIDE OF SCREEN
					//
					//light += clamp(prefilteredColor, vec3(0.0), vec3(10.0)) * fScale * AO;
					

				#else //Simple reflect-sample
					
					//bias the ray to avoid self-termination
					float rayOffset = 0.0;
					rayOffset += (1.0 - NDotV) * 0.25;
					rayOrigin += sampleRayDir * rayOffset;
					SDFMarchResult sdfRes = MarchSDFScene(rayOrigin, sampleRayDir);
					if(sdfRes.bHit)
					{
						float d = length(interpolatorWorldSpacePos.xyz - sdfRes.HitPosWS.xyz);
						/* light.r += 1.0 - clamp(d / 10.0, 0.0, 1.0); */
						//light.r = 0.3;
						//map hitPos to NDC
						vec3 posNDC = (vec4(sdfRes.HitPosWS, 1.0) * WorldToViewMatrix).xyz;
						posNDC.xy *= SceneBuffer.CameraZoom;
						posNDC.x /= SceneBuffer.ScreenRatio;
						posNDC.xy /= (0.1f + posNDC.z);
						vec2 sceneUV = (posNDC.xy + 1.0) * 0.5;
						d = MapToRange(d, 0.0, 20.0, 0.0, 1.0);
						d *= d;
						//d = 0.0;
						float reflectionLerpT = clamp(d + (s.Roughness * s.Roughness) - 0.1, 0.0, 1.0);
						//prevSceneMip = 0.0;
						vec3 prevSceneColor = textureLod(PrevFrameSceneColor, sceneUV.xy, 0.0).rgb;
						vec3 prevSceneColorBlur = textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb * 0.75;
						prevSceneColor = mix(prevSceneColor, prevSceneColorBlur, reflectionLerpT);
						/* vec2 envBRDF = texture(BRDF_LUT, vec2(NDotV, (s.Roughness))).rg;
						vec3 fScale = vec3(F * envBRDF.x + envBRDF.y); */
						vec3 fresnelCustom = FresnelSchlickRoughness(NDotV, vec3(0.5), s.Roughness);
						vec3 fScale = fresnelCustom;
						light += max(vec3(0.0), prevSceneColor * fScale * AO);
					}
					else
					{
						//light.g += sampleRayDir.y * 2.0;
						//light.r += min(1.0, abs(sdfRes.tRes - SDF_MAX_DISTANCE));
						//light.rb += vec2(0.6); //DEBUG
					}
				#endif //Importance Sample
				}
				}
				#endif //SSR

			}
			#endif//Indirect Light

			`+scGetBasePassLightShading()+ /* glsl */`
		
		#if 1//PROJECTOR
			{
				float divisor = interpolatorProjectorSpacePos.z;
				//divisor = 5.0;
				vec2 projNDCPos = vec2(interpolatorProjectorSpacePos.x / divisor, interpolatorProjectorSpacePos.y / divisor);
				projNDCPos.y *= SceneBuffer.ProjectorXYRatio;
				vec2 projUVPos = MapToRange(projNDCPos, -1.0, 1.0, 0.0, 1.0);
				vec3 projectorColor = vec3(SceneBuffer.ProjectorIntensity);
				//vec3 projectorColor = SceneBuffer.ProjectorColor * SceneBuffer.ProjectorIntensity;
				//s.Ambient = projectorColor * SceneBuffer.AmbientLight;
				//light += s.Albedo * s.Ambient;
				if(projUVPos.x < 1.0 && projUVPos.x > 0.0 && projUVPos.y > 0.0 && projUVPos.y < 1.0 && interpolatorProjectorSpacePos.z > 0.0)
				{

					const highp float rectPow = 32.f;
					highp float rectCircleLength = pow(abs(projNDCPos.x), rectPow) + pow(abs(projNDCPos.y), rectPow);
					float projMask = 1.0;
					if(rectCircleLength < 0.8)
					{
						projMask = 1.0;
					}
					else
					{
						const float maxTh = 1.0;
						float sm = MapToRange(min(maxTh, rectCircleLength), 0.8, maxTh, 1.0, 0.0);
						projMask = sm * sm;
						//projMask = 0.0;
					}

				#if 1 //EDGE FADE-OUT
					const float edgeThres = 0.99;

					float scale = 1.0;

					if(abs(projNDCPos.x) > edgeThres)
					{
						scale *= MapToRange(abs(projNDCPos.x), edgeThres, 1.0, 1.0, 0.0);
					}
					if(abs(projNDCPos.y) > edgeThres)
					{
						scale *= MapToRange(abs(projNDCPos.y), edgeThres, 1.0, 1.0, 0.0);
					}
					scale *= scale;
					projMask *= scale;
				#endif


					#if 0
					vec2 uvInv = vec2(projUVPos.x, 1.0 - projUVPos.y);
					#else
					vec2 uvInv = projUVPos;
					#endif
					vec3 projTexture = texture(ProjectorTexture, uvInv.xy).rgb;
					projectorColor *= vec3(0.00) + projTexture * 2.0;
					//projectorColor = Contrast(projectorColor, 1.5);
					projectorColor = pow(projectorColor, vec3(SceneBuffer.ProjectorIntensityCurve));
					projectorColor *= projMask;
					vec3 l = normalize(SceneBuffer.ProjectorPos - s.PixelPos);
					float NdotL = max(dot(s.Normal, l), 0.0);

					float shadow = 1.f; 
				#if 1//SHADOW
					{
						shadow = SDFComputeShadow(s.PixelPos, SceneBuffer.ProjectorPos);
					}
				#endif//SHADOW

					if(bEmitter > 0)
					{
						light += (s.Albedo * 0.75) * projectorColor / PI * 0.5 * (1.0 + 2.0 * (1.0 - min(1.0, Brightness)));
					}
					else
					{
						light += (/* s.Albedo * */ (1.0 - min(1.0, s.Metalness * 1.5))) * projectorColor / PI * NdotL * shadow * 0.5;
					}
					
					
				}

			}
		#endif//PROJECTOR
			
			//Fade objects close to camera
			light *= GetCamDistFade();

			OutColor = light;
			//OutColor = vec4(light, 1.0);
		}`;
}


export function GetShaderSourceIntegrateSpecularBRDFPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec2 outBRDF;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scPBRDeclarations()+ /* glsl */`

	in vec2 vsOutTexCoords;

	vec2 IntegrateBRDF(float NdotV, float roughness)
	{
	    vec3 V;
	    V.x = sqrt(1.0 - NdotV*NdotV);
	    V.y = 0.0;
	    V.z = NdotV;

	    float A = 0.0;
	    float B = 0.0;

	    vec3 N = vec3(0.0, 0.0, 1.0);

	    const uint SAMPLE_COUNT = 1024u;
	    for(uint i = 0u; i < SAMPLE_COUNT; ++i)
	    {
	        vec2 Xi = Hammersley(i, SAMPLE_COUNT);
	        vec3 H  = ImportanceSampleGGX(Xi, N, roughness);
	        vec3 L  = normalize(2.0 * dot(V, H) * H - V);

	        float NdotL = max(L.z, 0.0);
	        float NdotH = max(H.z, 0.0);
	        float VdotH = max(dot(V, H), 0.0);

	        if(NdotL > 0.0)
	        {
	            float G = GeometrySmith_IBL(N, V, L, roughness);
	            float G_Vis = (G * VdotH) / (NdotH * NdotV);
	            float Fc = pow(1.0 - VdotH, 5.0);

	            A += (1.0 - Fc) * G_Vis;
	            B += Fc * G_Vis;
	        }
	    }
	    A /= float(SAMPLE_COUNT);
	    B /= float(SAMPLE_COUNT);
	    return vec2(A, B);
	}

	void main()
	{
		vec2 uv = vec2(vsOutTexCoords.x, vsOutTexCoords.y);
		outBRDF = IntegrateBRDF(uv.x, uv.y);
	}`;
}

export function GetShaderSourceProceduralLineRenderVS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	
	layout(location = 0) in vec2 VertexBuffer;

	` +scDeclareSceneUniformBuffer()+/* glsl */`

	uniform vec3 LineStart;
	uniform vec3 LineEnd;
	uniform vec2 LineScale;

	out vec2 interpolatorTexCoords;
	out vec3 interpolatorNormal;
	out vec3 interpolatorTangent;
	out vec3 interpolatorWorldSpacePos;

	void main()
	{

		vec3 biTangent = LineEnd - LineStart;

		vec3 camPos = SceneBuffer.CameraPosition;

		#if 1
		vec3 dirToCamFromStart = normalize(camPos - LineStart);
		vec3 dirToCamFromEnd = normalize(camPos - LineEnd);
		#else
		
		#endif
		
		vec3 tangentStart = cross(biTangent, dirToCamFromStart);
		vec3 tangentEnd = cross(biTangent, dirToCamFromEnd);

		#if 1
		dirToCamFromStart = cross(tangentStart, biTangent);
		dirToCamFromEnd = cross(tangentEnd, biTangent);
		#endif

		//Fox perspective
		//tangentStart = tangentEnd;

		tangentStart = (normalize(tangentStart));
		tangentEnd = (normalize(tangentEnd));


		vec3 verts0 = LineStart - tangentStart * LineScale.x;
		vec3 verts1 = LineEnd - tangentEnd * LineScale.y;
		vec3 verts2 = LineEnd + tangentEnd * LineScale.y;
		vec3 verts3 = LineStart + tangentStart * LineScale.x;

		////Clockwise, starting from left down
		vec2 uv0 = vec2(0.0, 0.0);
		vec2 uv1 = vec2(0.0, 1.0);
		vec2 uv2 = vec2(1.0, 1.0);
		vec2 uv3 = vec2(1.0, 0.0);

		uint vertId = uint(gl_VertexID);
		vec3 pos = verts3;
		interpolatorTexCoords = uv3;
		interpolatorNormal = dirToCamFromStart;
		interpolatorTangent = tangentStart;
		if(vertId == 0u || vertId == 3u)
		{
			pos = verts0;
			interpolatorTexCoords = uv0;
			interpolatorNormal = dirToCamFromStart;
			interpolatorTangent = tangentStart;
		}
		else if(vertId == 1u)
		{
			pos = verts1;
			interpolatorTexCoords = uv1;
			interpolatorNormal = dirToCamFromEnd;
			interpolatorTangent = tangentEnd;
		}
		else if(vertId == 2u || vertId == 4u)
		{
			pos = verts2;
			interpolatorTexCoords = uv2;
			interpolatorNormal = dirToCamFromEnd;
			interpolatorTangent = tangentEnd;
		}

		interpolatorWorldSpacePos = pos;

		pos.xyz -= SceneBuffer.CameraPosition; 
		pos.xy *= SceneBuffer.CameraZoom;
		pos.x /= SceneBuffer.ScreenRatio;

		gl_Position = vec4(pos.xy,  pos.z / 20.0, (0.0001f + pos.z));
	}`;
}

export function GetShaderSourceProceduralRectRenderVS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	
	layout(location = 0) in vec2 VertexBuffer;

	` +scDeclareSceneUniformBuffer()+/* glsl */`

	uniform sampler2D RopePositionsTexture; 

	out vec2 interpolatorTexCoords;
	out vec3 interpolatorNormal;
	out vec3 interpolatorTangent;
	out vec3 interpolatorWorldSpacePos;

	void main()
	{

		int instanceId = int(gl_InstanceID);
		int offset = instanceId * 3;
		//get cur instance pos
		vec3 instancePos = vec3(
			texelFetch(RopePositionsTexture, ivec2(offset + 0, 0), 0).r,
			texelFetch(RopePositionsTexture, ivec2(offset + 1, 0), 0).r,
			texelFetch(RopePositionsTexture, ivec2(offset + 2, 0), 0).r
		);
		//vec3 instancePos = vec3(2.0, 0.1 * float(instanceId), -4.0);

		vec3 camPos = SceneBuffer.CameraPosition;
		float scale = SceneBuffer.SwitchRopeSphereRadius;
		
		vec3 biTangent = vec3(0.0, 1.0, 0.0);
		vec3 dirToCam = normalize(camPos - instancePos);
		vec3 tangent = normalize(cross(dirToCam, biTangent));
		biTangent = cross(dirToCam, tangent);

		vec3 verts0 = instancePos - tangent * scale - biTangent * scale;
		vec3 verts1 = instancePos - tangent * scale + biTangent * scale;
		vec3 verts2 = instancePos + tangent * scale + biTangent * scale;
		vec3 verts3 = instancePos + tangent * scale - biTangent * scale;

		////Clockwise, starting from left down
		vec2 uv0 = vec2(0.0, 0.0);
		vec2 uv1 = vec2(0.0, 1.0);
		vec2 uv2 = vec2(1.0, 1.0);
		vec2 uv3 = vec2(1.0, 0.0);

		uint vertId = uint(gl_VertexID);
		vec3 pos = verts3;
		interpolatorTexCoords = uv3;
		interpolatorNormal = dirToCam;
		interpolatorTangent = tangent;
		if(vertId == 0u || vertId == 3u)
		{
			pos = verts0;
			interpolatorTexCoords = uv0;
			interpolatorNormal = dirToCam;
			interpolatorTangent = tangent;
		}
		else if(vertId == 1u)
		{
			pos = verts1;
			interpolatorTexCoords = uv1;
			interpolatorNormal = dirToCam;
			interpolatorTangent = tangent;
		}
		else if(vertId == 2u || vertId == 4u)
		{
			pos = verts2;
			interpolatorTexCoords = uv2;
			interpolatorNormal = dirToCam;
			interpolatorTangent = tangent;
		}

		interpolatorWorldSpacePos = pos;

		pos.xyz -= SceneBuffer.CameraPosition; 
		pos.xy *= SceneBuffer.CameraZoom;
		pos.x /= SceneBuffer.ScreenRatio;

		gl_Position = vec4(pos.xy,  pos.z / 20.0, (0.0001f + pos.z));
	}`;
}

function scShapeRenderDefines( bCircle : boolean)
{
	const bcircle = `#define SHAPE_CIRCLE ` + (bCircle ? `1` : `0`);

	return bcircle;
}

export function GetShaderSourceProceduralShapeRenderPS(sceneSDFFunc : string, bCircle : boolean) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;
	precision highp samplerCube;

	`+scShapeRenderDefines(bCircle)+ /* glsl */`

	layout(location = 0) out vec3 OutColor;

	` +scDeclareSceneUniformBuffer()+/* glsl */`

	uniform float Brightness;
	uniform float Metalness;
	uniform int EUseLightprobe;//0-disabled, 1-default, 2-debug
	uniform int LightprobeIndex;
	uniform int LightprobesGridSize;
	uniform int bEmitter;
	uniform mat4 WorldToViewMatrix;

	uniform sampler2D AlbedoTexture;
	uniform sampler2D NormalTexture;
	uniform sampler2D RoughnessTexture;
	uniform sampler2D BRDF_LUT;

	uniform sampler2D ProjectorTexture;
	uniform sampler2D PrevFrameSceneColor;
	uniform sampler2D PrevFrameSceneColorBlur;
	uniform samplerCube CubemapTexture;
	uniform sampler2D SphericalTexture;
	uniform sampler2D SHCoefficientsTexture;
	uniform highp sampler3D SDFTexture;
	uniform highp sampler3D SDFTexture2;
	uniform sampler2D RopePositionsTextureSDF;
	
	in vec2 interpolatorTexCoords;
	in vec3 interpolatorNormal;
	in vec3 interpolatorTangent;
	in vec3 interpolatorWorldSpacePos;

	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`
	`+scPBRDeclarations()+ /* glsl */`
	`+scGetGIFuncsDeclarations()+ /* glsl */`
	`+scGetSDFFuncsDeclarations(sceneSDFFunc)+ /* glsl */`
	

	void main()
	{
		vec3 light = vec3(0.0);
		PBRSettings s;

		vec2 texSampleCoords = interpolatorTexCoords.xy * 0.25;

		s.Albedo = texture(AlbedoTexture, texSampleCoords.xy).rgb;
		//OutColor = vec3(0.8, 0.6, 0.1); return;
		
		
		{
		#if SHAPE_CIRCLE
			//hemiSphere
			{
				vec2 uv = interpolatorTexCoords.xy;
				uv = uv * 2.0 - 1.0;
				uv *= -1.0;
				if(length(uv) > 1.0)
				{
					discard;
				}
				s.Normal = normalize(vec3(uv, sqrt(1.0 - dot(uv, uv))));
			}
		#else
			{
				const highp float rectPow = 16.f;
				vec2 ndcPos = MapToRange(interpolatorTexCoords, 0.0, 1.0, -1.0, 1.0);
				highp float rectCircleLength = pow(abs(ndcPos.x), rectPow) + pow(abs(ndcPos.y), rectPow);
				const float thres = 0.9;
				if(rectCircleLength > thres)
				{
					discard;
				}
			}
			//in cylindrical coords, compute normal in tangent space
			{
				float angle = MapToRange(interpolatorTexCoords.x, 1.0, 0.0, 0.0, PI);
				s.Normal.x = -cos(angle);
				s.Normal.y = (1.0 - interpolatorTexCoords.y) * 0.2;
				#if 1//edge fade
					{
						vec2 ndcNorm = (interpolatorTexCoords - 0.5) * 2.0;
						const float edgeThres = 0.95;
						if(abs(ndcNorm.y) > edgeThres)
						{
							float tt = MapToRange(abs(ndcNorm.y), edgeThres, 1.0, 0.0, 1.0);
							tt *= tt;
							s.Normal.y = mix(s.Normal.y, sign(ndcNorm.y) * -0.75, tt);
						}
					}
				#endif
				//s.Normal.y = 0.0;
				s.Normal.z = sin(angle);
				s.Normal = normalize(s.Normal);
			}
		#endif//SHAPE_CIRCLE
			#if 1 //Detail normal
			{
				vec3 detailNormal = texture(NormalTexture, texSampleCoords.xy).rgb;
				detailNormal = DecodeNormalTexture(detailNormal, 0.1);
				//s.Normal = normalize(detailNormal);
				s.Normal = normalize(vec3(s.Normal.xy + detailNormal.xy, s.Normal.z));
			}
			#endif
			//s.Normal = vec3(0.0, 0.0, 1.0);
			vec3 TBNNormal = normalize(interpolatorNormal);
			vec3 TBNTangent = normalize(interpolatorTangent) * -1.0;
			vec3 BiTangent = cross(interpolatorTangent, interpolatorNormal);
			mat3 TangentToWorldMatrix = mat3(
				(TBNTangent),
				normalize(BiTangent),
				(TBNNormal)
				);
			s.Normal = normalize(TangentToWorldMatrix * s.Normal);
			#if 0//DEBUG NORMALS
				vec3 c = max(vec3(0.0), s.Normal) * 0.85;
				//vec3 c = max(vec3(0.0), TBNTangent) * 0.85;
				OutColor = c;
				return;
			#endif
		}

		//s.Normal = interpolatorNormal;

		const float roughnessMin = 0.05; 
		#if 1 //RGH from texture
		s.Roughness = max(0.0, texture(RoughnessTexture, interpolatorTexCoords.xy).r);
		s.Roughness = min(1.0, s.Roughness + roughnessMin);
		#else
		s.Roughness = SceneBuffer.UniformRoughness + roughnessMin;
		#endif

		s.Metalness = Metalness;

		s.PixelPos = interpolatorWorldSpacePos;
		s.vToCam = normalize(SceneBuffer.CameraPosition.xyz - interpolatorWorldSpacePos);

		s.LightColor = SceneBuffer.LightColor;

		//Ambient Occlusion
		float AO = 1.0;

		`+scGetBasePassLightShading()+ /* glsl */`

		float NDotV = max(dot(s.Normal, s.vToCam), 0.0);
		vec3 F0 = vec3(SceneBuffer.FresnelMin);
		F0 = mix(F0, s.Albedo, s.Metalness);
		vec3 F = FresnelSchlickRoughness(NDotV, F0, s.Roughness);
		vec3 kS = F;
		vec3 kD = 1.0 - kS;

		light += s.Albedo * 0.025 * kD;

		//Ambient Diffuse
		{
			#if 1
				vec3 samplePos = s.Normal;
				samplePos.y *= -1.0;
			
				//SH Basis Coefficients
				if(EUseLightprobe > 0)
				{
					if(EUseLightprobe == 1)
					{
						vec3 indirectDiffuse = GetIrradianceFromLightprobes(samplePos, interpolatorWorldSpacePos);
						indirectDiffuse *= vec3(0.8, 0.9, 1.0);
						//indirectDiffuse *= 10.0;
						//indirectDiffuse *= vec3(0.25, 0.5, 1.0) * 2.0;
						light += vec3(indirectDiffuse * s.Albedo / PI * kD) * AO * (1.0 + SceneBuffer.AmbientLight + 0.25);
						//light += vec3(vec3(0.5) * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight);
						light += vec3(vec3(0.0, 0.0, 1.0) * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight);

					}
					else//TODO: Remove, used for GI Probes visualizer
					{
						//light = SampleDebugLightprobe(samplePos);
					}
				}
			#else //Regular diffuse
			{
				vec3 indirectDiffuse = SceneBuffer.LightColor;
				indirectDiffuse = vec3(0.5);
				light += vec3(indirectDiffuse * s.Albedo / PI * kD) * AO * (SceneBuffer.AmbientLight + 0.3);
			}

			#endif//Indirect diffuse
		}

		#if 1//SSR
		{
			vec3 sampleRayDir = normalize(reflect(-s.vToCam, s.Normal));
			vec3 rayOrigin = interpolatorWorldSpacePos;
			float a = s.Roughness;
			float a2 = a * a;
			if(dot(sampleRayDir, vec3(0.0, 0.0, 1.0)) > 0.0)
			{

		#if 1 //Importance sample rough surfaces
			//get low-discrepancy random value
			const uint SAMPLE_COUNT = 16u;
			float totalWeight = 0.0;   
			vec3 prefilteredColor = vec3(0.0);  
			uint coordX = uint(gl_FragCoord.x);
			uint coordY = uint(gl_FragCoord.y);
			for(uint i = 0u; i < SAMPLE_COUNT; ++i)
			{
				vec2 Xi = Hammersley(i, SAMPLE_COUNT);
				//vec2 Xi = vec2(hashf(float(i + coordX * coordY)), hash(i));
				vec3 curLightDir = ImportanceSampleGGX(Xi, sampleRayDir, s.Roughness);
				curLightDir = normalize(curLightDir);
				rayOrigin = interpolatorWorldSpacePos + sampleRayDir * 0.01;
				float NdotL = max(dot(s.Normal, curLightDir), 0.0);
				SDFMarchResult sdfRes = MarchSDFScene(rayOrigin, curLightDir);
				if(sdfRes.bHit)
				{
					vec3 posNDC = (vec4(sdfRes.HitPosWS, 1.0) * WorldToViewMatrix).xyz;
					posNDC.xy *= SceneBuffer.CameraZoom;
					posNDC.x /= SceneBuffer.ScreenRatio;
					posNDC.xy /= (0.1f + posNDC.z);
					vec2 sceneUV = (posNDC.xy + 1.0) * 0.5;
					
					vec3 h = normalize(s.vToCam + curLightDir);
					float NdotH  = max(dot(s.Normal, h), 0.0);
					float HdotV  = max(dot(s.vToCam, h), 0.0);
					float NDF = NormalDistributionGGXFast(NdotH, a2);  
					float pdf = max(0.0, (NDF * NdotH / (4.0 * HdotV + 0.0001))); 
					
					vec3 prevSceneColor = textureLod(PrevFrameSceneColor, sceneUV.xy, 0.0).rgb;
					prevSceneColor = mix(prevSceneColor, textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb, s.Roughness * 0.75);
					//vec3 prevSceneColor = textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb;
					float wght = NdotL * pdf;
					prefilteredColor += prevSceneColor * wght;
					totalWeight += wght;
				}
			}
			prefilteredColor = prefilteredColor / (totalWeight + 0.0001);
			vec2 envBRDF = texture(BRDF_LUT, vec2(NDotV, (s.Roughness))).rg;
			vec3 fScale = vec3(F * envBRDF.x + envBRDF.y);
			if(!any(isnan(prefilteredColor)))
			{
				light += max(vec3(0.0), prefilteredColor) * fScale * AO * 0.5;
			}
			

			//
			/////////////TODO DO SMTH ABOUT THIS, DISABLE SAMPLES OUTSIDE OF SCREEN
			//
			//light += clamp(prefilteredColor, vec3(0.0), vec3(10.0)) * fScale * AO;
			

		#else //Simple reflect-sample
			
			//bias the ray to avoid self-termination
			float rayOffset = 0.0;
			rayOffset += (1.0 - NDotV) * 0.25;
			rayOrigin += sampleRayDir * rayOffset;
			SDFMarchResult sdfRes = MarchSDFScene(rayOrigin, sampleRayDir);
			if(sdfRes.bHit)
			{
				float d = length(interpolatorWorldSpacePos.xyz - sdfRes.HitPosWS.xyz);
				/* light.r += 1.0 - clamp(d / 10.0, 0.0, 1.0); */
				//light.r = 0.3;
				//map hitPos to NDC
				vec3 posNDC = (vec4(sdfRes.HitPosWS, 1.0) * WorldToViewMatrix).xyz;
				posNDC.xy *= SceneBuffer.CameraZoom;
				posNDC.x /= SceneBuffer.ScreenRatio;
				posNDC.xy /= (0.1f + posNDC.z);
				vec2 sceneUV = (posNDC.xy + 1.0) * 0.5;
				d = MapToRange(d, 0.0, 20.0, 0.0, 1.0);
				d *= d;
				//d = 0.0;
				float reflectionLerpT = clamp(d + (s.Roughness * s.Roughness) - 0.1, 0.0, 1.0);
				//prevSceneMip = 0.0;
				vec3 prevSceneColor = textureLod(PrevFrameSceneColor, sceneUV.xy, 0.0).rgb;
				vec3 prevSceneColorBlur = textureLod(PrevFrameSceneColorBlur, sceneUV.xy, 0.0).rgb * 0.75;
				prevSceneColor = mix(prevSceneColor, prevSceneColorBlur, reflectionLerpT);
				/* vec2 envBRDF = texture(BRDF_LUT, vec2(NDotV, (s.Roughness))).rg;
				vec3 fScale = vec3(F * envBRDF.x + envBRDF.y); */
				vec3 fresnelCustom = FresnelSchlickRoughness(NDotV, vec3(0.5), s.Roughness);
				vec3 fScale = fresnelCustom;
				light += max(vec3(0.0), prevSceneColor * fScale * AO);
			}
			else
			{
				//light.g += sampleRayDir.y * 2.0;
				//light.r += min(1.0, abs(sdfRes.tRes - SDF_MAX_DISTANCE));
				//light.rb += vec2(0.6); //DEBUG
			}
		#endif //Importance Sample
		}
		}
		#endif //SSR

		//outColor = vec3(0.7, 0.7, 0.2);

		//edge fade
		{
			const float fadeMax = 0.5;
			float f = 1.0;
			vec2 ndcNorm = abs(interpolatorTexCoords - 0.5) * 2.0;
		#if !SHAPE_CIRCLE
			const float edgeThres = 0.9;
			if((ndcNorm.x) > edgeThres)
			{
				f *= MapToRange((ndcNorm.x), edgeThres, 1.0, 1.0, fadeMax);
			}
		#else
			float sl = length(interpolatorTexCoords - 0.5) * 2.0;
			const float edgeThres = 0.9;
			if(sl > edgeThres)
			{
				f *= MapToRange(sl, edgeThres, 1.0, 1.0, fadeMax);
			}
		#endif
			light *= f;
		}

		OutColor = light;

	}`;
}
import { Vector2 } from "../core/types";
import { scDeclareSceneUniformBuffer } from "./shaderBasePass";
import { scGetCommonFuncsDeclarations, scPBRDeclarations } from "./shaderCommon";
import { scGetDestrTileCommonFuncsDecl } from "./shaderDestrTiles";


export function scGetSDFFuncsDeclarations(sceneSDFFunc : string)
{
	return (/* glsl */`

	#define SDF_MAX_STEPS 64
	#define SDF_EPSILON 0.05
	#define SDF_MAX_DISTANCE 20.0

	float sdSphere( vec3 pointPos, vec3 spherePos, float r )
	{
	  	return length(pointPos - spherePos) - r;
	}
	float sdPlane(vec3 p, vec3 n, float d )
	{
	  	return dot(p,n) + d;
	}

	float sdBox( vec3 curPos, vec3 p, vec3 b )
	{
	    vec3 d = abs(curPos - p) - b;
	    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
	}

	float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
	{
	  vec3 pa = p - a, ba = b - a;
	  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	  return length( pa - ba*h ) - r;
	}

	float sdCappedCone( vec3 p, vec3 a, vec3 b, float ra, float rb )
	{
	  float rba  = rb-ra;
	  float baba = dot(b-a,b-a);
	  float papa = dot(p-a,p-a);
	  float paba = dot(p-a,b-a)/baba;
	  float x = sqrt( papa - paba*paba*baba );
	  float cax = max(0.0,x-((paba<0.5)?ra:rb));
	  float cay = abs(paba-0.5)-0.5;
	  float k = rba*rba + baba;
	  float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
	  float cbx = x-ra - f*rba;
	  float cby = paba - f;
	  float s = (cbx<0.0 && cay<0.0) ? -1.0 : 1.0;
	  return s*sqrt( min(cax*cax + cay*cay*baba,
	                     cbx*cbx + cby*cby*baba) );
	}

	float sdVerticalCapsule( vec3 curPos, vec3 p, float h, float r )
	{
		//p.y -= h * 0.5;
		p = curPos - p;
	  	p.y -= clamp( p.y, 0.0, h );
	  	return length( p ) - r;
	}

	float sdCustomMesh(vec3 curPos, highp sampler3D texture, vec3 meshBoundsCenter, vec3 meshBoundsExtent, float meshScale, bool bRayMeshIntersection)
	{
		float curD = SDF_MAX_DISTANCE;
		vec3 sdfMeshBoxMin = meshBoundsCenter - meshBoundsExtent;
		vec3 sdfMeshBoxMax = meshBoundsCenter + meshBoundsExtent;
		//sample SDF texture
		if(bRayMeshIntersection)
		{
			vec3 uvw = MapToRange(curPos, sdfMeshBoxMin, sdfMeshBoxMax, vec3(0.0), vec3(1.0));
			float distToSDFMeshAABB = sdBox(curPos, meshBoundsCenter, meshBoundsExtent);
			distToSDFMeshAABB = max(0.0, distToSDFMeshAABB);
			uvw = clamp(uvw, vec3(0.0), vec3(1.0));
			//if(all(greaterThan(uvw, vec3(0.0))) && all(lessThan(uvw, vec3(1.0))))
			{
				float meshDist = textureLod(texture, uvw, 0.0).r;
				//meshDist = 0.0;
				meshDist *= meshScale;
				curD = min(curD, meshDist + distToSDFMeshAABB);
			}
		}

		return curD;

	}

	`+sceneSDFFunc+ /* glsl */`
	
	struct SDFMarchResult
	{
		bool bHit;
		vec3 HitPosWS;
		float tRes;
		int NumSteps; //Debug
	};

	SDFMarchResult MarchSDFScene(vec3 rayOrigin, vec3 rayDir)
	{
		bool bHit = false;
		vec3 hitPositionWS = rayOrigin;
		
		// March along the ray until we hit something or reach maximum distance
		float t = 0.01;
		
		bool bRayMeshIntersection = false;
		{
			vec4 rayMeshIntersectionRes = RayAABBIntersection(rayOrigin, rayDir, SceneBuffer.SDFMeshBoundsCenter, SceneBuffer.SDFMeshBoundsExtent);
			bRayMeshIntersection = rayMeshIntersectionRes.x > 0.0;
			if(bRayMeshIntersection)
			{
				t = rayMeshIntersectionRes.z;
			}
		}

	#ifdef SECOND_SDF_MESH
		bool bRayMeshIntersection2 = false;
		{
			vec4 rayMeshIntersectionRes = RayAABBIntersection(rayOrigin, rayDir, SceneBuffer.SDFMeshBoundsCenter2, SceneBuffer.SDFMeshBoundsExtent2);
			bRayMeshIntersection2 = rayMeshIntersectionRes.x > 0.0;
			if(bRayMeshIntersection2)
			{
				t = min(t, rayMeshIntersectionRes.z);
			}
		}
	#endif //SECOND_SDF_MESH
			
		int numSteps = 0;
		for (int i = 0; i < SDF_MAX_STEPS; i++) {
		#if 0
			float hitDist = 0.0001 * t;
		#else
			float hitDist = 0.0001 * t;
			hitDist += 0.001 * float(i);
		#endif
			hitPositionWS = rayOrigin + t * rayDir;

		#ifdef SECOND_SDF_MESH
			float dist = sceneSDF(hitPositionWS, bRayMeshIntersection, bRayMeshIntersection2);
		#else
			float dist = sceneSDF(hitPositionWS, bRayMeshIntersection, false);
		#endif

		#if 1 //negative dist is a hit
			if(dist < 0.0)
			{
				bHit = true;
				t += dist;
				hitPositionWS = rayOrigin + t * rayDir;
				break; 
			}
		#endif
			if (abs(dist) < abs(hitDist))
			{
				bHit = true;
				break; 
			}
			else if( t > (SDF_MAX_DISTANCE))
			{
				break;
			}
			t += dist; // Move along the ray by the distance to the nearest surface
			numSteps+=1;
		}

		SDFMarchResult res;
		res.bHit = bHit;
		res.HitPosWS = hitPositionWS;
		res.tRes = t;
		res.NumSteps = numSteps;

		return res;

	}

	// https://iquilezles.org/articles/nvscene2008/rwwtt.pdf
	float calcAO( in vec3 pos, in vec3 nor )
	{
		float occ = 0.0;
	    float sca = 1.0;
	    for( int i=0; i<7; i++ )
	    {
			//const float stepSize = 0.25;
			const float stepSize = 0.2;
	        float h = 0.01 + 0.12*float(i) * stepSize;
	        float d = sceneSDF( pos + h*nor, false, false );
	        occ += (h-d)*sca;
	        sca *= 0.85;
	        if( occ>0.35 ) break;
	    }
	    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );
	}

	float CalcAOPrecise( in vec3 pos, in vec3 nor )
	{
		float occ = 0.0;
	    float sca = 1.0;

		bool bRayMeshIntersection = RayAABBIntersection(pos, nor, SceneBuffer.SDFMeshBoundsCenter, SceneBuffer.SDFMeshBoundsExtent).x > 0.0;
	#ifdef SECOND_SDF_MESH
		bool bRayMeshIntersection2 = RayAABBIntersection(pos, nor, SceneBuffer.SDFMeshBoundsCenter2, SceneBuffer.SDFMeshBoundsExtent2).x > 0.0;
	#endif

	    for( int i=0; i<3; i++ )
	    {
			//const float stepSize = 0.25;
			float ran = hash( uint(gl_FragCoord.x) * uint(i) + 1920U* uint(gl_FragCoord.y) + (1920U*1080U) );  
			//ran = 1.0;
			float stepSize = 0.1 + ran * 0.01;
	        float h = stepSize + stepSize * float(i);
		#ifdef SECOND_SDF_MESH
			float d = sceneSDF( pos + h*nor, bRayMeshIntersection, bRayMeshIntersection2 );
	    #else
			float d = sceneSDF( pos + h*nor, bRayMeshIntersection, false );
		#endif//SECOND_SDF_MESH
			float diff = h-d;
	        occ += pow(diff, 1.0)*0.1*sca;
	        //sca *= 0.1;
	        if( occ>0.35 ) break;
	    }
	    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );
	}

	float CalcAOSpreadHemisphere( in vec3 p, in vec3 n, in float maxDist, in float falloff )
	{
		const int nbIte = 16;
		const float nbIteInv = 1./float(nbIte);
		const float rad = 1.-1. * nbIteInv; //Hemispherical factor (self occlusion correction)

		float ao = 0.0;

		for( int i=0; i<nbIte; i++ )
		{
			float l = hashf(float(i)) * maxDist;
			vec3 rd = normalize(n + randomHemisphereDir(n, l ) * rad /* * 0.0 */) * l; // mix direction with the normal for self occlusion problems!

			ao += (l - max(sceneSDF( p + rd, false, false ), 0.)) / maxDist * falloff;
		}

		return clamp( 1.-ao * nbIteInv, 0., 1.);
	}
	
	float CalcAOSpread( in vec3 pos, in vec3 nor, float ra )
	{
	    float occ = 0.0;
	    {
			int i = 1;
	        float h = 0.01 + 4.0*pow(float(i)/31.0,2.0);
	        vec2 an = hash2( ra + float(i)*13.1 )*vec2( 3.14159, 6.2831 );
	        vec3 dir = vec3( sin(an.x)*sin(an.y), sin(an.x)*cos(an.y), cos(an.x) );
			
			vec3 randomVec = vec3(0.0);
			float d = 1.0;

			float occScale = 0.25;
			
			float randDirScale = 0.1;
			vec3 tangentBasis = normalize(vec3(1.0, 0.0, 1.0) + dir * randDirScale);
			vec3 tangent   = normalize(tangentBasis - nor * dot(tangentBasis, nor));
			vec3 bitangent = cross(nor, tangent);
			mat3 TBN       = mat3(tangent, bitangent, nor);
			float stepSize = 0.5 + 0.75 * float(i) + ra * 0.02;
			
			randomVec = normalize(vec3(1.0, 0.0, 0.0) + vec3(dir.xy * 0.1, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += abs(d) * occScale;
			}

			randomVec = normalize(vec3(-1.0, 0.0, 0.0) + vec3(dir.xy * 0.1, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += abs(d) * occScale;
			}

			randomVec = normalize(vec3(0.0, 1.0, 0.0) + vec3(dir.xy * 0.1, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += abs(d) * occScale;
			}

			randomVec = normalize(vec3(0.0, -1.0, 0.0) + vec3(dir.xy * 0.1, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += abs(d) * occScale;
			}

			#if 1//ADDITIONAL CORNERS
			occScale *= 0.1;
			float rndVecScale = 0.1;
			randomVec = normalize(vec3(1.0, 1.0, 0.0) + vec3(dir.xy * rndVecScale, 0.0));
			randomVec = normalize(TBN * randomVec);
			stepSize *= 2.5;
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += (abs(d)) * occScale;
			}
			randomVec = normalize(vec3(-1.0, -1.0, 0.0) + vec3(dir.xy * rndVecScale, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += (abs(d)) * occScale;
			}
			randomVec = normalize(vec3(-1.0, 1.0, 0.0) + vec3(dir.xy * rndVecScale, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += (abs(d)) * occScale;
			}
			randomVec = normalize(vec3(1.0, -1.0, 0.0) + vec3(dir.xy * rndVecScale, 0.0));
			randomVec = normalize(TBN * randomVec);
			d = sceneSDF(pos + randomVec * stepSize, false, false);
			if(d < 0.0)
			{
				occ += (abs(d)) * occScale;
			}
			#endif


			occ *= occ * occ;
			//occ *= occ;


			
	        //occ += 10.0 * clamp( 5.0*sceneSDF( pos + h*dir ) * 0.1 /* / h */, -1.0, 1.0);
	    }

		occ = 1.0 - clamp( occ, 0.0, 1.0 );
	    return occ;

	}

	

	
	`);
}



export function scGetSDFShadowDeclarations()
{
	return (/* glsl */`

	/* float sceneSDFShadow(vec3 curPos, bool bRayMeshIntersection)
	{
		float curD = SDF_MAX_DISTANCE;
		
		//main painting pos
		{
			curD = min(curD, sdBox(curPos, SceneBuffer.OccluderPlaneCenter, SceneBuffer.OccluderPlaneExtent));
		}
		return curD;
	} */

	
	

	
	
	`);
}


export function GetShaderSourceSDFRenderPS(sceneSDFFunc : string) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	in vec2 vsOutTexCoords;

	uniform sampler2D ProjectorTexture;
	uniform sampler2D BlueNoiseTexture;
	uniform highp sampler3D SDFTexture;
	uniform highp sampler3D SDFTexture2;

	uniform sampler2D RopePositionsTextureSDF;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSDFFuncsDeclarations(sceneSDFFunc)+ /* glsl */`

	vec3 calcNormal( in vec3 pos, bool bRayMeshIntersection, bool bRayMeshIntersection2 )
	{
		// inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
	    vec3 n = vec3(0.0);
	    for( int i=0; i<4; i++ )
	    {
	        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
	        n += e*sceneSDF(pos+0.0005*e, bRayMeshIntersection, bRayMeshIntersection2);
	      //if( n.x+n.y+n.z>100.0 ) break;
	    }
	    return normalize(n);
	}

	float sdCone( in vec3 pointPos, in vec3 conePos, in vec2 c, float h )
	{
		vec3 p = pointPos - conePos;
    	vec2 q = h*vec2(c.x,-c.y)/c.y;
    	vec2 w = vec2( length(p.xz), p.y );
		
		vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    	vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    	float k = sign( q.y );
    	float d = min(dot( a, a ),dot(b, b));
    	float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
		return sqrt(d)*sign(s);
	}

	float bayer2(vec2 a){
		a = floor(a);
		return fract( dot(a, vec2(.5, a.y * .75)) );
	}
	#define bayer4(a)   (bayer2( .5*(a))*.25+bayer2(a))
	#define bayer8(a)   (bayer4( .5*(a))*.25+bayer2(a))
	#define bayer16(a)  (bayer8( .5*(a))*.25+bayer2(a))
	#define bayer32(a)  (bayer16(.5*(a))*.25+bayer2(a))
	#define bayer64(a)  (bayer32(.5*(a))*.25+bayer2(a))
	#define bayer128(a) (bayer64(.5*(a))*.25+bayer2(a))

	void main()
	{
		vec3 rayOrigin = SceneBuffer.CameraPosition;
		rayOrigin.z -= 0.1;

		vec3 color = vec3(0.0);
		

		//NDC to world
		vec2 posNDC = MapToRange(vsOutTexCoords, 0.0, 1.0, -1.0, 1.0);
		vec3 posWS = vec3(posNDC.x * SceneBuffer.ScreenRatio, posNDC.y, 0.0);
		posWS.xy /= SceneBuffer.CameraZoom;
		const float farPlane = 20.0;
		posWS.xy *= farPlane;
		posWS.xy += SceneBuffer.CameraPosition.xy;
		posWS.z = SceneBuffer.CameraPosition.z + farPlane;

		vec3 rayDir = normalize(posWS - rayOrigin);

		////////////////
		//	SDF
		///////////////
		bool bHit = false;
		vec3 hitPositionWS = rayOrigin;
		// March along the ray until we hit something or reach maximum distance
		float t = 0.01;
		int numSteps = 0;

		bool bRayMeshIntersection = false;
		{
			vec4 rayMeshIntersectionRes = RayAABBIntersection(rayOrigin, rayDir, SceneBuffer.SDFMeshBoundsCenter, SceneBuffer.SDFMeshBoundsExtent);
			bRayMeshIntersection = rayMeshIntersectionRes.x > 0.0;
			/* if(bRayMeshIntersection)
			{
				t = rayMeshIntersectionRes.z;
			} */
		}

		////////
		//Main
		////////
		SDFMarchResult payload = MarchSDFScene(rayOrigin, rayDir);
		bHit = payload.bHit;
		hitPositionWS = payload.HitPosWS;
		t = payload.tRes;
		numSteps = payload.NumSteps;

		float depth = clamp(t / 10.0, 0.0, 1.0);

		vec3 normal = vec3(0.0, 0.0, -1.0);
		normal = calcNormal(hitPositionWS, true, true);

		float occ = 1.0;
		float ran = hash( uint(gl_FragCoord.x) + 1920U* uint(gl_FragCoord.y) + (1920U*1080U) );  
		occ = CalcAOPrecise( hitPositionWS, normal );
		/* occ *= CalcAOSpread( hitPositionWS, normal, ran );
		occ *= CalcAOSpreadHemisphere( hitPositionWS, normal, 1., 0.5 ); */

		/* hitPositionWS.x = clamp(hitPositionWS.x, -3.5, 3.5); 
		hitPositionWS.y = clamp(hitPositionWS.y, -2.5, 2.5);  */
		
		
		if(bHit)
		{
			//color = vec3(1.0 - depth);
			//color = vec3(0.1);
			color = vec3(1.0, 0.5, 0.1);
			//color *= (1.0 - depth);
			//color = vec3(1.0);
			//color =  -normal;
			color *= occ;

			float NdotV = 1.0 - max(dot(normal, vec3(0.0, 0.0, -1.0)), 0.0);

			color *= NdotV;

			//Shade
			{
				//vec3 lightpos = vec3(0.0, 1.0, -5.0);
				vec3 lightpos = SceneBuffer.SpotlightPos2;
				vec3 l = normalize(lightpos - hitPositionWS);
				float distanceToCurLight = length(lightpos - hitPositionWS);
				float NdotL = max(dot(normal, l), 0.0);
				float ssThickness = 0.0;
				float sssRes = 0.0;
				#if 0
				//SSS
				{
					//compute thickness
					const float sstMin = 0.01;
					float sst = sstMin;
					vec3 ssrdir = vec3(0.0, 0.0, 1.0);
					//ssrdir = l;
					for( int i=0; i<32; i++ )
					{
						vec3 ssp = hitPositionWS + ssrdir * sst;

						float h = sceneSDF(ssp, true);

						if(h < 0.0)
						{
							ssThickness += 0.01;
						}

						sst += 0.01;
					}

					//color += (1.0 - min(1.0, ssThickness));
					sssRes = (1.0 - min(1.0, ssThickness));
					/* float wrap = sssRes;
					NdotL = max(0.0, (NdotL + wrap) / (1.0 + wrap)); */
				}
				#endif

				
				//NdotL = 1.0;
				//occ = 1.0;
				float lightShadow = SDFComputeShadow(hitPositionWS, lightpos);\
				float attenuation = clamp(1.f - pow(distanceToCurLight / 20.0, 2.0), 0.f, 1.f);
				color += vec3(1.0, 1.0, 0.8) * NdotL * attenuation * occ /* * lightShadow */;
				//color = vec3(1.0, 1.0, 0.8)  * attenuation  * lightShadow;
			}

			color *= 0.3;

			float stepsVis = float(numSteps) / float(SDF_MAX_STEPS);
			stepsVis *= stepsVis;
			//color = vec3(stepsVis);

		#if 0
			float edge = 0.0;
			{
				SDFMarchResult payloadT;
				const float offset = 0.001;
				const bool bMarchMesh = true;
				payloadT = MarchSDFScene(rayOrigin, rayDir + vec3(offset, 0.0, 0.0));
				edge += length(normal - calcNormal(payloadT.HitPosWS, bMarchMesh, bMarchMesh));
				payloadT = MarchSDFScene(rayOrigin, rayDir + vec3(0.0, offset, 0.0));
				edge += length(normal - calcNormal(payloadT.HitPosWS, bMarchMesh, bMarchMesh));
			}
			color = vec3(edge) * 0.5 * vec3(1.0, 0.75, 0.1);
		#endif

		}
		else
		{
			color = vec3(1.0, 0.0, 1.0) * 0.3;
		}

		color = max(vec3(0.0), normal) * 0.85; 

		
		outColor = color;
	

	}`;
}

export function GetShaderSourceSDFVolumeLightRenderPS(sceneSDFFunc : string) {
	return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	in vec2 vsOutTexCoords;

	uniform sampler2D ProjectorTexture;
	uniform sampler2D BlueNoiseTexture;
	uniform highp sampler3D SDFTexture;
	uniform highp sampler3D SDFTexture2;
	uniform sampler2D RopePositionsTextureSDF;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scPBRDeclarations()+ /* glsl */`
	`+scGetSDFFuncsDeclarations(sceneSDFFunc)+ /* glsl */`
	
	void main()
	{
		vec3 rayOrigin = SceneBuffer.CameraPosition;
		rayOrigin.z -= 0.1;

		vec3 color = vec3(0.0);

		//NDC to world
		vec2 posNDC = MapToRange(vsOutTexCoords, 0.0, 1.0, -1.0, 1.0);
		vec3 posWS = vec3(posNDC.x * SceneBuffer.ScreenRatio, posNDC.y, 0.0);
		posWS.xy /= SceneBuffer.CameraZoom;
		const float farPlane = 20.0;
		posWS.xy *= farPlane;
		posWS.xy += SceneBuffer.CameraPosition.xy;
		posWS.z = SceneBuffer.CameraPosition.z + farPlane;

		vec3 rayDir = normalize(posWS - rayOrigin);

		vec3 volumeCurPos = rayOrigin;
		const float marchLength = 10.0;
		const float marchStepSize = marchLength * 0.025;
		float tt = 0.0;
		vec3 scatter = vec3(0.0);
		float transmittance = 1.0;


		//start marching from the scene bbox start
		float sceneNearPlaneZ = SceneBuffer.SceneBoundsCenter.z + SceneBuffer.SceneBoundsExtent.z;
		if(rayOrigin.z < sceneNearPlaneZ)
		{
			float sdfNearPlane = abs(rayOrigin.z - sceneNearPlaneZ);
			float advT = sdfNearPlane / rayDir.z;
			//rayOrigin = rayOrigin + advT * rayDir;
			//rayOrigin = rayOrigin + sdfNearPlane * rayDir;
		}

		vec2 pixelPos = vec2(gl_FragCoord.xy);
		float jitter = textureLod(BlueNoiseTexture, pixelPos / 512.0 /* uvOffset + vsOutTexCoords * vec2(SceneBuffer.ScreenRatio, 1.0) * 1.0 */, 0.0).r;
		const float kGoldenRatioConjugate = 0.61803398875f;
        //jitter = fract(jitter + mod(SceneBuffer.Time * 100.0, 64.0) * kGoldenRatioConjugate);
		jitter = fract(jitter + hashf(SceneBuffer.Time) * 0.1);


		//float jitter = smoothNoise(rayDir);
		float rayStartOffset = jitter;
		for (float tt = rayStartOffset; tt < marchLength; ) {

			float od = marchStepSize;

			volumeCurPos = rayOrigin + tt * rayDir;

			#if 0 //Flashlight
			{
				vec3 lightColor = vec3(1.0, 1.0, 1.0);
				vec3 lightPos = SceneBuffer.SpotlightPos;
				//lightPos = SceneBuffer.CameraPosition;
				vec3 lightTargetPos = SceneBuffer.SpotlightTargetPos;
				vec3 p = volumeCurPos;
				float d = length(lightPos - p);
				vec3 l = (lightPos - p) / d;
				float SpotlightFalloff = GetSpotlightFalloff(lightPos, normalize(lightTargetPos - lightPos), volumeCurPos, SceneBuffer.SpotlightOuterAngle, SceneBuffer.SpotlightInnerAngleOffset);
				scatter += lightColor * od * SpotlightFalloff * 0.05;
			}
			#endif

			#if 0 //Spotlight
			{
				vec3 lightColor = vec3(1.0, 1.0, 1.0);
				vec3 lightPos = SceneBuffer.SpotlightPos2;
				vec3 lightTargetPos = SceneBuffer.SpotlightTargetPos2 - vec3(0.0, 0.5, 0.0);
				vec3 p = volumeCurPos;
				float d = length(lightPos - p);
				vec3 l = (lightPos - p) / d;
				float SpotlightFalloff = GetSpotlightFalloff(lightPos, normalize(lightTargetPos - lightPos), volumeCurPos, 1.0, 0.75);
				#if 0 //occlusion
				{
					float occlusion = 1.0;
					
					//cone
					{
						float minT = 0.1;
						float maxT = d;
						float t = minT;

						for( int i=0; i<32; i++ )
						{
							float h = SDF_MAX_DISTANCE;
							vec3 curPos = p + l*t;
							{
								const float handleRad = 0.1;
								h = min(h, sdCappedCone(curPos,
									 SceneBuffer.LightSwitchHandlePos + vec3(0.0, 0.35, 0.0),
									 SceneBuffer.LightSwitchHandlePos - vec3(0.0, 0.35, 0.0), 
									 handleRad * 0.7,
									handleRad));
							}

							t += h;

							if ((h)<(0.001*t))
							{
								occlusion *= 0.0;
								break; 
							}
							else if( t > maxT)
							{
								break;
							}
						}
					}

					#if 0//Rope
					{
						int numInstances = int(SceneBuffer.SwitchRopeNumSpheres);
						float radius = SceneBuffer.SwitchRopeSphereRadius * 2.0;
						for(int i = 0; i < numInstances; i++)
						{
							int offset = i * 3;
							vec3 instancePos = vec3(
								texelFetch(RopePositionsTextureSDF, ivec2(offset + 0, 0), 0).r,
								texelFetch(RopePositionsTextureSDF, ivec2(offset + 1, 0), 0).r,
								texelFetch(RopePositionsTextureSDF, ivec2(offset + 2, 0), 0).r
							);
							bool bIntersection = RaySphereIntersection(p, l, instancePos, radius);
							if(bIntersection)
							{
								occlusion *= 0.0;
								break; 
							}
						}
					}
					#endif
					SpotlightFalloff *= occlusion;
				}
				#endif
				SpotlightFalloff *= SpotlightFalloff;
				float att = GetLightAttenuation(d, SceneBuffer.PointLightRadius1 * 2.0);
				//att = 1.0;
				scatter += lightColor * od * SpotlightFalloff * att * 0.75;
			}
			#endif

			#if 1
			//projector
			{
				//get projector space pos
				vec3 posProjSpace = (vec4(volumeCurPos, 1.0) * SceneBuffer.WorldToProjectorSpaceMat).xyz;
				posProjSpace.xy *= SceneBuffer.ProjectorZoom;
				posProjSpace.z += 0.01;
				float divisor = posProjSpace.z;
				vec2 projNDCPos = vec2(posProjSpace.x / divisor, posProjSpace.y / divisor);
				projNDCPos.y *= SceneBuffer.ProjectorXYRatio;
				vec2 projUVPos = MapToRange(projNDCPos, -1.0, 1.0, 0.0, 1.0);
				vec3 projectorColor = SceneBuffer.ProjectorColor * SceneBuffer.ProjectorIntensity;
				
				if(projUVPos.x < 1.0 && projUVPos.x > 0.0 && projUVPos.y > 0.0 && projUVPos.y < 1.0 && posProjSpace.z > 0.0)
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
						projMask = 0.0;
					}

					#if 0
					vec2 uvInv = vec2(projUVPos.x, 1.0 - projUVPos.y);
					#else
					vec2 uvInv = projUVPos;
					#endif
					vec3 projTexture = texture(ProjectorTexture, uvInv.xy).rgb;
					projectorColor *= (vec3(0.01) + projTexture) * SceneBuffer.ProjectorIntensityVolume;
					projectorColor *= projMask;

					float ProjectorFalloff = SceneBuffer.ProjectorFalloffVolume;
					float d = posProjSpace.z;
					float num = pow(clamp(1.0 - pow(d /  ProjectorFalloff, 4.0), 0.0, 1.0), 2.0);
					float denom =  (d * d + 1.0);
					float attenuation = num / denom;

					scatter += projectorColor * od * attenuation /* * transmittance */;
				}
			}
			#endif

			#if 0
			//sdf mesh
			{
				vec3 sdfMeshBoxExtent = vec3(1.0, 1.0, 1.0);
				vec3 sdfMeshBoxCenter = vec3(0.0, 0.0, -4.0);
				vec3 sdfMeshBoxMin = sdfMeshBoxCenter - sdfMeshBoxExtent;
				vec3 sdfMeshBoxMax = sdfMeshBoxCenter + sdfMeshBoxExtent;
				{
					//sample SDF texture
					if(sdBox(volumeCurPos, sdfMeshBoxCenter, sdfMeshBoxExtent) < 0.1)
					{
						//we are inside mesh volume
						//map world pos to UVW space
						vec3 uvw = MapToRange(volumeCurPos, sdfMeshBoxMin, sdfMeshBoxMax, vec3(0.0), vec3(1.0));
						//uvw.z = 1.0 - uvw.z;
						//curD = min(curD, textureLod(SDFTexture, uvw, 0.0).r);
						float dToL = texture(SDFTexture, uvw).r;
						/* dToL = max(0.2, dToL);
						scatter += dToL * 0.3 * od; */
						if(dToL < 0.01)
						{
							scatter += 0.25 * od;
						}
						else
						{
							scatter += 0.025 * od;
						}
					}
				}
			}
			#endif

			//decay for distance to observer
			transmittance *= exp2(-od * 0.1); //exp decay

			//jitter = hashf(rayDir.x + rayDir.y + tt);
			//jitter = 0.0;
			//tt += marchStepSize + jitter * marchStepSize;
			tt += marchStepSize;

		}

		color = scatter;

		outColor = color;

	}
	
	`;
}


export function GetShaderSourceDestrSceneVolumeLightRenderPS(numTiles : Vector2, tileSize : Vector2, sceneSDFFunc : string) {
	return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	in vec2 vsOutTexCoords;

	uniform sampler2D ProjectorTexture;
	uniform sampler2D BlueNoiseTexture;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scPBRDeclarations()+ /* glsl */`
	`+scGetSDFFuncsDeclarations(sceneSDFFunc)+ /* glsl */`

	`+scGetDestrTileCommonFuncsDecl(numTiles, tileSize)+ /* glsl */`
	
	void main()
	{
		vec3 rayOrigin = SceneBuffer.CameraPosition;
		rayOrigin.z -= 0.1;

		vec3 color = vec3(0.0);

		//NDC to world
		vec2 posNDC = MapToRange(vsOutTexCoords, 0.0, 1.0, -1.0, 1.0);
		vec3 posWS = vec3(posNDC.x * SceneBuffer.ScreenRatio, posNDC.y, 0.0);
		posWS.xy /= SceneBuffer.CameraZoom;
		const float farPlane = 20.0;
		posWS.xy *= farPlane;
		posWS.xy += SceneBuffer.CameraPosition.xy;
		posWS.z = SceneBuffer.CameraPosition.z + farPlane;

		vec3 rayDir = normalize(posWS - rayOrigin);

		vec3 volumeCurPos = rayOrigin;
		const float marchLength = 5.0;
		const float marchStepSize = marchLength * 0.025 * 0.25;
		float tt = 0.0;
		vec3 scatter = vec3(0.0);
		float transmittance = 1.0;


		vec2 pixelPos = vec2(gl_FragCoord.xy);
		float jitter = textureLod(BlueNoiseTexture, pixelPos / 512.0 /* uvOffset + vsOutTexCoords * vec2(SceneBuffer.ScreenRatio, 1.0) * 1.0 */, 0.0).r;
		const float kGoldenRatioConjugate = 0.61803398875f;
        //jitter = fract(jitter + mod(SceneBuffer.Time * 100.0, 64.0) * kGoldenRatioConjugate);
		jitter = fract(jitter + hashf(SceneBuffer.Time) * 0.1);

		//jitter = 0.0;

		//float jitter = smoothNoise(rayDir);
		float rayStartOffset = jitter;
		for (float tt = rayStartOffset; tt < marchLength; ) {

			float od = marchStepSize;

			volumeCurPos = rayOrigin + tt * rayDir;

			#if 1 //Spotlight
			{
				vec3 lightColor = vec3(1.0, 1.0, 1.0);
				vec3 lightPos = SceneBuffer.SpotlightPos2;
				vec3 lightTargetPos = lightPos + SceneBuffer.SpotlightTargetPos2;
				vec3 p = volumeCurPos;
				float d = length(lightPos - p);
				vec3 l = (lightPos - p) / d;

				//float SpotlightFalloff = GetSpotlightFalloff(lightPos, normalize(lightTargetPos - lightPos), volumeCurPos, 1.0, 0.75);
				float SpotlightFalloff = GetSpotlightFalloff(lightPos, normalize(lightTargetPos - lightPos), volumeCurPos, 0.85, 0.15);
				
				//SpotlightFalloff *= SpotlightFalloff;

				float att = GetLightAttenuation(d, SceneBuffer.PointLightRadius1 * 2.0);
				att *= att;
				//att = 1.0;
				scatter += lightColor * od * SpotlightFalloff * att * 0.75;
			}
			#endif

			#if 1 //Direct Light
			if(volumeCurPos.z < 0.01)
			{
				vec3 lightDir = normalize(vec3(1.0, -1.0, -0.5));

				vec3 planeCenterWS = vec3(0.0, 0.0, 0.001 * 0.5);
				vec3 planeExtentWS = vec3(2.0, 2.0, 0.001);

				vec4 payload = RayAABBIntersection(volumeCurPos, -lightDir, planeCenterWS, planeExtentWS);
				bool bHit = payload.x > 0.0;

				if(bHit)
				{
					vec3 p = volumeCurPos - lightDir * payload.z;

					ivec3 curTileId = GetTileIdFromPosWS(p);

					if(curTileId.z > 0)
					{
						float health = texelFetch(ProjectorTexture, curTileId.xy, 0).r;
						if(health < kTileHealthThres)
						{
							scatter += vec3(1.0) * 0.01 * 2.0;
						}
					}
				}
				

			}
			#endif


			//decay for distance to observer
			transmittance *= exp2(-od * 0.1); //exp decay

			//jitter = hashf(rayDir.x + rayDir.y + tt);
			//jitter = 0.0;
			//tt += marchStepSize + jitter * marchStepSize;
			tt += marchStepSize;

		}

		color = scatter;

		outColor = color;

	}
	
	`;
}

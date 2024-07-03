import { scDeclareSceneUniformBuffer } from "./shaderBasePass";
import { scGetCommonFuncsDeclarations, scGetSphericalFuncsDeclarations, scPBRDeclarations } from "./shaderCommon";


export function GetShaderSourceComputeIrradianceRectPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec4 outIrradiance;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`

	uniform sampler2D SphericalTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		vec3 irradiance = vec3(0.0);

		float sampleDelta = 0.025;
		
		
		//get cur angle
		vec2 curSphCoord = UVToSpherical(vsOutTexCoords) + vec2(SceneBuffer.PhiOffset, SceneBuffer.ThetaOffset);
		float phiStart = curSphCoord.x - PI;
		float phiEnd = phiStart + TWOPI;
		float thetaStart = curSphCoord.y - PIDIVTWO;
		float thetaEnd = curSphCoord.y;

		#if 0
		vec2 sampleUV = SphericalToUVMod(vec2(curSphCoord.x + SceneBuffer.PhiOffset, curSphCoord.y + SceneBuffer.ThetaOffset));
		irradiance = texture(SphericalTexture, vec2(sampleUV.x, /* 1.0 -  */sampleUV.y)).rgb /* * cos(theta) * sin(theta) */;
		#endif

		#if 1
		float nrSamples = 0.0; 
		irradiance = vec3(0.0);
		sampleDelta = PI / 16.0;
		for(float phi = phiStart; phi < phiEnd; phi += sampleDelta)
		{
			float thetaLocalSpace = 0.0;
		    for(float theta = thetaStart; theta < thetaEnd; theta += sampleDelta)
		    {
		        vec2 sampleUV = SphericalToUVMod(vec2(phi, theta));
		        irradiance += texture(SphericalTexture, vec2(sampleUV.x, 1.0 - sampleUV.y)).rgb * cos(thetaLocalSpace) * sin(thetaLocalSpace);
		        nrSamples++;
				thetaLocalSpace += sampleDelta;
		    }
		}
		irradiance = PI * irradiance * (1.0 / float(nrSamples));
		#endif

		outIrradiance = vec4(irradiance.rgb, 1.0);

	}`;
}

export function GetShaderSourceComputeIrradianceCubemapPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec4 outIrradiance;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`

	uniform samplerCube CubemapTexture;
	uniform sampler2D SHCoefficientsTexture;

	uniform vec3 Normal;
	uniform vec3 Tangent;

	in vec2 vsOutTexCoords;

	#define SPHVEC 0.5773502588272095f

	const int n_bands = 3; // Assuming you have a predefined number of bands

	struct SHSample {
	    vec3 sph;
	    vec3 vec;
	    float coeff[n_bands*(n_bands+1)];
	};

	float rand(vec2 co) {
	    return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453);
	}

	float rand2( vec2 n ) { return fract(sin(vec2(n.x,n.y+1.0))*vec2(43758.5453123,22578.1459123)).x; }


	void SH_setup_spherical_samples(inout SHSample samples[1], int sqrt_n_samples) {
	    float oneoverN = 1.0 / float(sqrt_n_samples);
	    int i = 0;
	    for (int a = 0; a < sqrt_n_samples; a++) {
	        for (int b = 0; b < sqrt_n_samples; b++) {
	            float x = (float(a) + rand(vec2(float(a), float(b)))) * oneoverN;
	            float y = (float(b) + rand(vec2(float(a), float(b)))) * oneoverN;
	            float theta = 2.0 * acos(sqrt(1.0 - x));
	            float phi = 2.0 * PI * y;
	            samples[i].sph = vec3(theta, phi, 1.0);
	            vec3 vec = vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	            samples[i].vec = vec;
	            for (int l = 0; l < n_bands; ++l) {
	                for (int m = -l; m <= l; ++m) {
	                    int index = l*(l+1)+m;
	                    //samples[i].coeff[index] = SH(l, m, theta, phi); // Assuming SH is a function defined elsewhere
	                }
	            }
	            ++i;
	        }
	    }
	}

	void main()
	{

		//Get cur tangent dir from UV
		const vec3 nLeftDown = (vec3(-SPHVEC, -SPHVEC, SPHVEC));
		vec3 nRightDown = (vec3(SPHVEC, -SPHVEC, SPHVEC));
		vec3 nLeftUp = (vec3(-SPHVEC, SPHVEC, SPHVEC));
		vec3 nRightUp = (vec3(SPHVEC, SPHVEC, SPHVEC));
		vec3 dirUp = normalize(mix(nLeftUp, nRightUp, vsOutTexCoords.x));
		vec3 dirDown = normalize(mix(nLeftDown, nRightDown, vsOutTexCoords.x));
		vec3 dir = normalize(mix(dirUp, dirDown, vsOutTexCoords.y));
		vec3 Bitangent = cross(Normal, Tangent);
		mat3 TangentToWorldMatrix = mat3(
		(Tangent),
		normalize(Bitangent),
		(Normal)
		);

		
		dir = normalize(TangentToWorldMatrix * dir);
		
	#if 1 
	{
		vec3 up    = vec3(0.0, 1.0, 0.0);
		vec3 right = normalize(cross(up, dir));
		up         = normalize(cross(dir, right));

		//Transform to SH space

	#if 1//convolution kernels
		const float h0 = PI;
		const float h1 = 2.09439510239; //2pi / 3
		const float h2 = 0.785398163397; //pi / 4
		const float kernels[9] = float[9](h0, h1, h1, h1, h2, h2, h2, h2, h2);
	#else
		float h0 = SceneBuffer.SHCoef.x;
		float h1 = SceneBuffer.SHCoef.y; //2pi / 3
		float h2 = SceneBuffer.SHCoef.z; //pi / 4
		float kernels[9] = float[9](h0, h1, h1, h1, h2, h2, h2, h2, h2);
	#endif

		//SH Basis Coefficients
		float SHCoefArr[9] = float[9](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	#if 0 //MK Sample
		const int sqrt_n_samples = 64;
		float oneoverN = 1.0 / float(sqrt_n_samples);
	    for (int a = 0; a < sqrt_n_samples; a++) {
	        for (int b = 0; b < sqrt_n_samples; b++) {
	            float x = (float(a) + rand(vec2(float(a), float(b)))) * oneoverN;
	            float y = (float(b) + rand(vec2(float(b), float(a + b)))) * oneoverN;
	            float theta = 2.0 * acos(sqrt(1.0 - x));
	            float phi = 2.0 * PI * y;
	            //samples[i].sph = vec3(theta, phi, 1.0);
	            vec3 sampleVec = vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

				float funcOutput = dot(vec3(0.333), texture(CubemapTexture, sampleVec).rgb);

				float shBasisArr[9] = SHBasis2(sampleVec);

				for(int i = 0; i < 9; i++)
				{
					SHCoefArr[i] += funcOutput * shBasisArr[i] * kernels[i];
				}
	        }
	    }

		for(int i = 0; i < 9; i++)
		{
			SHCoefArr[i] *= (4.0 * PI) / float(sqrt_n_samples * sqrt_n_samples);
		}

	#else

		for(int i = 0; i < 9; i++)
		{
			SHCoefArr[i] = texelFetch(SHCoefficientsTexture, ivec2(i, 0), 0).r;
		}

	#endif
		

		float shBasisArr[9] = SHBasis2(dir);
		float res = 0.0;
		for(int i = 0; i < 9; i++)
		{
			res += shBasisArr[i] * SHCoefArr[i];
		}

		outIrradiance = vec4(vec3(res), 1.0);
		return;

	}
	#endif

	#if 0
		outIrradiance = texture(CubemapTexture, dir).rgba;
		return;
	#else
		vec3 irradiance = vec3(0.0);  

		vec3 up    = vec3(0.0, 1.0, 0.0);
		vec3 right = normalize(cross(up, dir));
		up         = normalize(cross(dir, right));


		const float phiRangeLength = TWOPI;
		const float thetaRangeLength = PIDIVTWO;
		float thetaEnd = thetaRangeLength * 0.75;


		float sampleDelta = 0.1;
		float phiSampleDeltaRad = phiRangeLength / 16.0;
	#if 1
		float nrSamples = 0.0; 
		
		for(float theta = 0.0; theta < thetaEnd; theta += sampleDelta)
		{
			for(float phi = 0.0; phi < phiRangeLength; phi += sampleDelta)
			{
		    	// spherical to cartesian (in tangent space)
		    	vec3 tangentSample = vec3(sin(theta) * cos(phi),  sin(theta) * sin(phi), cos(theta));
		    	// tangent space to world
		    	vec3 sampleVec = tangentSample.x * right + tangentSample.y * up + tangentSample.z * dir; 
				
		    	irradiance += texture(CubemapTexture, sampleVec).rgb * cos(theta) * sin(theta);
		    	nrSamples++;
			}
		}
		irradiance = PI * irradiance * (1.0 / float(nrSamples));
		outIrradiance = vec4(irradiance.rgb, 1.0);
	#else
				vec3 tangentSample = vec3(sin(0.0) * cos(0.0),  sin(0.0) * sin(0.0), cos(0.0));
		        // tangent space to world
		        vec3 sampleVec = tangentSample.x * right + tangentSample.y * up + tangentSample.z * dir; 
			
		        irradiance += texture(CubemapTexture, sampleVec).rgb;
				outIrradiance = vec4(irradiance.rgb, 1.0);
				return;
	#endif
		
	#endif

	}`;
}

export function GetShaderSourceConvertRectEnvToCubemapPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec4 outColor;

	#define SPHVEC 0.5773502588272095f

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`

	uniform sampler2D SphericalTexture;

	uniform vec3 Normal;
	uniform vec3 Tangent;

	in vec2 vsOutTexCoords;

	void main()
	{

		//Get cur tangent dir from UV
		const vec3 nLeftDown = (vec3(-SPHVEC, -SPHVEC, SPHVEC));
		vec3 nRightDown = (vec3(SPHVEC, -SPHVEC, SPHVEC));
		vec3 nLeftUp = (vec3(-SPHVEC, SPHVEC, SPHVEC));
		vec3 nRightUp = (vec3(SPHVEC, SPHVEC, SPHVEC));
		vec3 dirUp = normalize(mix(nLeftUp, nRightUp, vsOutTexCoords.x));
		vec3 dirDown = normalize(mix(nLeftDown, nRightDown, vsOutTexCoords.x));
		vec3 dir = normalize(mix(dirUp, dirDown, vsOutTexCoords.y));
		vec3 Bitangent = cross(Normal, Tangent);
		mat3 TangentToWorldMatrix = mat3(
		(Tangent),
		normalize(Bitangent),
		(Normal)
		);

		dir = normalize(TangentToWorldMatrix * dir);
		vec2 uv = DirToUV(dir);
		vec3 color = texture(SphericalTexture, uv).rgb;

		//Pseudo HDR
		const float scale = (0.25);
		const float powScale = 3.0;
		//color = pow(color + vec3(scale), vec3(powScale)) - vec3(pow(scale, powScale));
		color = pow(color, vec3(2.0)) * 2.0;

		outColor = vec4(color, 1.0);

	}`;
}

export function GetShaderSourceComputeReflectedIrradianceCubemapPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec4 outColor;

	#define SPHVEC 0.5773502588272095f

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`
	`+scPBRDeclarations()+ /* glsl */`


	uniform samplerCube CubemapTexture;

	uniform vec3 Normal;
	uniform vec3 Tangent;
	uniform float CurRoughness;
	uniform float SrcMip;

	in vec2 vsOutTexCoords;


	void main()
	{

		//Get cur tangent dir from UV
		const vec3 nLeftDown = (vec3(-SPHVEC, -SPHVEC, SPHVEC));
		vec3 nRightDown = (vec3(SPHVEC, -SPHVEC, SPHVEC));
		vec3 nLeftUp = (vec3(-SPHVEC, SPHVEC, SPHVEC));
		vec3 nRightUp = (vec3(SPHVEC, SPHVEC, SPHVEC));
		vec3 dirUp = normalize(mix(nLeftUp, nRightUp, vsOutTexCoords.x));
		vec3 dirDown = normalize(mix(nLeftDown, nRightDown, vsOutTexCoords.x));
		vec3 dir = normalize(mix(dirUp, dirDown, vsOutTexCoords.y));
		vec3 Bitangent = cross(Normal, Tangent);
		mat3 TangentToWorldMatrix = mat3(
		(Tangent),
		normalize(Bitangent),
		(Normal)
		);

		dir = normalize(TangentToWorldMatrix * dir);
		
		//Approximation has view and normal aligned
		vec3 N = dir;    
    	vec3 R = N;
    	vec3 V = R;
		//float CurRoughness = 0.04;
		float rgh = clamp(CurRoughness, 0.04, 1.0);

    	const uint SAMPLE_COUNT = 1024u;
    	float totalWeight = 0.0;   
    	vec3 prefilteredColor = vec3(0.0);     
    	for(uint i = 0u; i < SAMPLE_COUNT; ++i)
    	{
			//get low-discrepancy random value
    	    vec2 Xi = Hammersley(i, SAMPLE_COUNT);
			//get cur halfway vec, bias samples towards the specular lobe
			//for perfectly smooth surfaces every H = N
			//in reality the shape will also be affected by the viewing angle but which is V = N here for simplicity
    	    vec3 H  = ImportanceSampleGGX(Xi, N, rgh);
			//compute light dir based on halfway vec
			//for perfectly smooth surfaces we just sample alogn normal (N = L = R = V)
    	    vec3 L  = normalize(2.0 * dot(V, H) * H - V);
			//weight by cosine lobe
    	    float NdotL = max(dot(N, L), 0.0);
    	    if(NdotL > 0.0)
    	    {
    	        prefilteredColor += textureLod(CubemapTexture, L, 0.0).rgb * NdotL;
    	        //prefilteredColor += texture(CubemapTexture, L).rgb * NdotL;
    	        totalWeight      += NdotL;
    	    }
    	}
    	prefilteredColor = prefilteredColor / totalWeight;

    	outColor = vec4(prefilteredColor, 1.0);

	}`;
}


export function GetShaderSourceFullscreenLineRenderVS()
{
	return (/* glsl */ `#version 300 es

	precision highp float;

	void main()
	{
		vec3 pos = vec3(-1.0, 0.0, 1.0);
		if(gl_VertexID == 1)
		{
			pos.x = 1.0;
		}
		gl_Position = vec4(pos, 1.0);
	}`);
}


export function GetShaderSourceComputeSHCoefPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outCoefficient;

	` +scDeclareSceneUniformBuffer()+/* glsl */`
	`+scGetCommonFuncsDeclarations()+ /* glsl */`
	`+scGetSphericalFuncsDeclarations()+ /* glsl */`

	uniform samplerCube CubemapTexture;

	float rand(vec2 co) {
	    return fract(sin(dot(co.xy, vec2(12.9898, 78.233))) * 43758.5453);
	}

	float rand2( vec2 n ) { return fract(sin(vec2(n.x,n.y+1.0))*vec2(43758.5453123,22578.1459123)).x; }

	void main()
	{
		//Transform to SH space

	#if 1//convolution kernels
		const float h0 = PI;
		//const float h1 = 2.09439510239; //2pi / 3
		const float h1 = PIDIVTWO; //custom, looks better when l < 3
		const float h2 = 0.785398163397; //pi / 4
		const float kernels[9] = float[9](h0, h1, h1, h1, h2, h2, h2, h2, h2);
	#else
		float h0 = SceneBuffer.SHCoef.x;
		float h1 = SceneBuffer.SHCoef.y; //2pi / 3
		float h2 = SceneBuffer.SHCoef.z; //pi / 4
		float kernels[9] = float[9](h0, h1, h1, h1, h2, h2, h2, h2, h2);
	#endif

		//SH Basis Coefficients
		float SHCoefArrR[9] = float[9](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		float SHCoefArrG[9] = float[9](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		float SHCoefArrB[9] = float[9](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

		const int sqrt_n_samples = 64;
		float oneoverN = 1.0 / float(sqrt_n_samples);
	    for (int a = 0; a < sqrt_n_samples; a++) {
	        for (int b = 0; b < sqrt_n_samples; b++) {
	            float x = (float(a) + rand(vec2(float(a), float(b)))) * oneoverN;
	            float y = (float(b) + rand(vec2(float(b), float(a + b)))) * oneoverN;
	            float theta = 2.0 * acos(sqrt(1.0 - x));
	            float phi = 2.0 * PI * y;
	            //samples[i].sph = vec3(theta, phi, 1.0);
	            vec3 sampleVec = vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));

				vec3 radiance = texture(CubemapTexture, sampleVec).rgb;

				float shBasisArr[9] = SHBasis2(sampleVec);

				for(int i = 0; i < 9; i++)
				{
					SHCoefArrR[i] += radiance.r * shBasisArr[i] * kernels[i];
				}
				for(int i = 0; i < 9; i++)
				{
					SHCoefArrG[i] += radiance.g * shBasisArr[i] * kernels[i];
				}
				for(int i = 0; i < 9; i++)
				{
					SHCoefArrB[i] += radiance.b * shBasisArr[i] * kernels[i];
				}
	        }
	    }

		for(int i = 0; i < 9; i++)
		{
			SHCoefArrR[i] *= (4.0 * PI) / float(sqrt_n_samples * sqrt_n_samples);
			SHCoefArrG[i] *= (4.0 * PI) / float(sqrt_n_samples * sqrt_n_samples);
			SHCoefArrB[i] *= (4.0 * PI) / float(sqrt_n_samples * sqrt_n_samples);
		}

		//12 samples, 4 coeff per color
		//Pack by color
		float SHCoefArr[12];
		SHCoefArr[0] = SHCoefArrR[0];
		SHCoefArr[1] = SHCoefArrR[1];
		SHCoefArr[2] = SHCoefArrR[2];
		SHCoefArr[3] = SHCoefArrR[3];

		SHCoefArr[4] = SHCoefArrG[0];
		SHCoefArr[5] = SHCoefArrG[1];
		SHCoefArr[6] = SHCoefArrG[2];
		SHCoefArr[7] = SHCoefArrG[3];

		SHCoefArr[8] = SHCoefArrB[0];
		SHCoefArr[9] = SHCoefArrB[1];
		SHCoefArr[10] = SHCoefArrB[2];
		SHCoefArr[11] = SHCoefArrB[3];

		//Output sample based on cur coord
		int xCoord = int(gl_FragCoord.x);
		float outputCoef = SHCoefArr[xCoord];
		outCoefficient = outputCoef;
		//outCoefficient = float(xCoord);


	}`;
}
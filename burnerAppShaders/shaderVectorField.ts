import { Vector2 } from "../core/types";
import { scDeclareCurlNoiseFuncs } from "./shaderParticles";



export function GetShaderSourceComputeDivergencePS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outDiv;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		//if solid cell - quit
		/* if(sampleCoord.x == 0 || sampleCoord.x == (kNumCells.x - 1) || sampleCoord.y == 0 || sampleCoord.y == (kNumCells.y - 1))
		{
			discard;
		} */

	#if 1 //Pressure method
		float velXLeft = texelFetch(VecFieldXTexture, sampleCoord - ivec2(1, 0), 0).r;
		float velXRight = texelFetch(VecFieldXTexture, sampleCoord + ivec2(1, 0), 0).r;
		float velYDown = texelFetch(VecFieldYTexture, sampleCoord - ivec2(0, 1), 0).r;
		float velYUp = texelFetch(VecFieldYTexture, sampleCoord + ivec2(0, 1), 0).r;

		vec2 velCur;
		velCur.x = texelFetch(VecFieldXTexture, sampleCoord, 0).r;
		velCur.y = texelFetch(VecFieldYTexture, sampleCoord, 0).r;

		if(sampleCoord.x == 0)
		{
			velXLeft = -velCur.x;
		}
		if(sampleCoord.x == (kNumCells.x - 1))
		{
			velXRight = -velCur.x;
		}
		if(sampleCoord.y == 0)
		{
			velYDown = -velCur.y;
		}
		if(sampleCoord.y == (kNumCells.y - 1))
		{
			velYUp = -velCur.y;
		}

		float divergence = -0.5 * (velXRight - velXLeft + velYUp - velYDown);
	#else
		int leftState = int(sampleCoord.x != 1);
		int rightState = int(sampleCoord.x != (kNumCells.x - 2));
		int downState = int(sampleCoord.y != 1);
		int upState = int(sampleCoord.y != (kNumCells.y - 2));

		int divisor = leftState + rightState + upState + downState;

		float velXLeft = texelFetch(VecFieldXTexture, sampleCoord, 0).r;
		float velXRight = texelFetch(VecFieldXTexture, sampleCoord + ivec2(1, 0), 0).r;

		float velYDown = texelFetch(VecFieldYTexture, sampleCoord, 0).r;
		float velYUp = texelFetch(VecFieldYTexture, sampleCoord + ivec2(0, 1), 0).r; 

		float divergence = velXRight - velXLeft + velYUp - velYDown;

		divergence = divergence / float(divisor);

		const float OverrelaxationConstant = 1.9;
		divergence *= OverrelaxationConstant;
	#endif

		outDiv = divergence;

	}`;
}

export function GetShaderSourceComputeCurlPS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outCurl;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		#if 1
		float velXLeft = texelFetch(VecFieldYTexture, sampleCoord - ivec2(1, 0), 0).r;
		float velXRight = texelFetch(VecFieldYTexture, sampleCoord + ivec2(1, 0), 0).r;
		float velYDown = texelFetch(VecFieldXTexture, sampleCoord - ivec2(0, 1), 0).r;
		float velYUp = texelFetch(VecFieldXTexture, sampleCoord + ivec2(0, 1), 0).r;
		#else
		float velXLeft = texelFetch(VecFieldXTexture, sampleCoord - ivec2(1, 0), 0).r;
		float velXRight = texelFetch(VecFieldXTexture, sampleCoord + ivec2(1, 0), 0).r;
		float velYDown = texelFetch(VecFieldYTexture, sampleCoord - ivec2(0, 1), 0).r;
		float velYUp = texelFetch(VecFieldYTexture, sampleCoord + ivec2(0, 1), 0).r;
		#endif

		float curl = velXRight - velXLeft + velYDown - velYUp;

		outCurl = curl * 0.5;

	}`;
}

export function GetShaderSourceApplyVorticityPS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform float DeltaTime;
	uniform float VorticityScale;

	uniform sampler2D CurlTexture;
	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		float curlLeft = texelFetch(CurlTexture, sampleCoord - ivec2(1, 0), 0).r;
		float curlRight = texelFetch(CurlTexture, sampleCoord + ivec2(1, 0), 0).r;
		float curlDown = texelFetch(CurlTexture, sampleCoord - ivec2(0, 1), 0).r;
		float curlUp = texelFetch(CurlTexture, sampleCoord + ivec2(0, 1), 0).r;

		float curlCur = texelFetch(CurlTexture, sampleCoord, 0).r;

		vec2 grad = vec2(abs(curlUp) - abs(curlDown), abs(curlRight) - abs(curlLeft)) * 0.5;

	#if 1 //NORMALIZE_CURL_GRADIENT
		grad /= length(grad) + 0.0001;
	#else
		grad = log(grad + vec2(1.0));//produces more chaotic swirls
		curlCur *= 10.0;
	#endif

		grad *= curlCur * VorticityScale;

		grad.y *= -1.0;

		vec2 vorticityForce = grad;

		float curVelX = texelFetch(VecFieldXTexture, sampleCoord, 0).r;
		float curVelY = texelFetch(VecFieldYTexture, sampleCoord, 0).r;


		outX = curVelX + vorticityForce.x * DeltaTime;
		outY = curVelY + vorticityForce.y * DeltaTime;

	}`;
}

export function GetShaderSourceComputePressurePS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outPressure;

	uniform sampler2D PressureTexture;
	uniform sampler2D DivergenceTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		float p0 = texelFetch(PressureTexture, sampleCoord - ivec2(1, 0), 0).r;
		float p1 = texelFetch(PressureTexture, sampleCoord + ivec2(1, 0), 0).r;
		float p2 = texelFetch(PressureTexture, sampleCoord - ivec2(0, 1), 0).r;
		float p3 = texelFetch(PressureTexture, sampleCoord + ivec2(0, 1), 0).r;

		float d = texelFetch(DivergenceTexture, sampleCoord, 0).r;

		outPressure = (d + p0 + p1 + p2 + p3) * 0.25;
		//outPressure = (d);

	}`;
}

export function GetShaderSourceClearPressurePS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outPressure;

	uniform sampler2D PressureTexture;

	uniform float PressureClearValue;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		float p = texelFetch(PressureTexture, sampleCoord, 0).r;

		outPressure = p * PressureClearValue;

	}`;
}

export function GetShaderSourceProjectPressurePS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform sampler2D PressureTexture;
	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		float p0 = texelFetch(PressureTexture, sampleCoord - ivec2(1, 0), 0).r;
		float p1 = texelFetch(PressureTexture, sampleCoord + ivec2(1, 0), 0).r;
		float p2 = texelFetch(PressureTexture, sampleCoord - ivec2(0, 1), 0).r;
		float p3 = texelFetch(PressureTexture, sampleCoord + ivec2(0, 1), 0).r;
		vec2 curPresGrad = vec2(p1 - p0, p3 - p2);

		float relaxationConst = 1.0;
		#if 1
		{
			relaxationConst = 1.9;
		}
		#endif

		float curVelX = texelFetch(VecFieldXTexture, sampleCoord, 0).r;
		float curVelY = texelFetch(VecFieldYTexture, sampleCoord, 0).r;

		outX = curVelX - curPresGrad.x;
		outY = curVelY - curPresGrad.y;

	}`;
}


export function GetShaderSourceAdvectVelocityPS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	`+scDeclareCurlNoiseFuncs()+ /* glsl */`

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	uniform float DeltaTime;
	uniform float Decay;
	uniform float Time;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);
		
		float curVelX = textureLod(VecFieldXTexture, vsOutTexCoords, 0.0).r;
		float curVelY = textureLod(VecFieldYTexture, vsOutTexCoords, 0.0).r;
		
		//convert coord to WS, advect, convert back to texture space

		vec2 coord = vsOutTexCoords - vec2(curVelX, curVelY) * DeltaTime;

		float newVelX = textureLod(VecFieldXTexture, coord, 0.0).r;
		float newVelY = textureLod(VecFieldYTexture, coord, 0.0).r;

		//curl noise
		#if 0
		{
			vec2 inVelocity = vec2(curVelX, curVelY);
			vec3 noisePosition = vec3(1.0);
			noisePosition.xy = vec2(sampleCoord);
			noisePosition *= 0.1;
			float noiseTime = Time * 1.025000;
			//noiseTime = 1.0;
			noisePosition.z = Time * 0.1;

			vec4 xNoisePotentialDerivatives = vec4(0.0);
    		vec4 yNoisePotentialDerivatives = vec4(0.0);
    		vec4 zNoisePotentialDerivatives = vec4(0.0);
			/* float curl = 0.75 + min(1.0, length(inVelocity)); 
			float noiseVelScale = 1.0 + min(1.0, length(inVelocity)); */
			float curl = 0.5;
			float noiseVelScale = 0.5;
			noiseVelScale = 0.025 + min(1.0, length(inVelocity) * 2.0);
			//noiseVelScale *= 2.5;

			//LQ
			{
				xNoisePotentialDerivatives = vec4(0.0);
				yNoisePotentialDerivatives = vec4(0.0);
				zNoisePotentialDerivatives = vec4(0.0);
				noiseTime = Time * 0.25 * 0.5;
				noisePosition *= 0.02;
				noiseVelScale = 0.5 + min(1.0, length(inVelocity) * 1.0);
				noiseVelScale *= 0.75;
				curl = 0.9;
				int i = 0;
				float scale = (1.0 / 2.0) * pow(2.0, float(i));
				float noiseScale = 1.0;
				xNoisePotentialDerivatives += simplexNoiseDerivatives(vec4(noisePosition * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				yNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(123.4, 129845.6, -1239.1)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				zNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(-9519.0, 9051.0, -123.0)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				vec3 noiseVelocity = vec3(
					zNoisePotentialDerivatives[1] - yNoisePotentialDerivatives[2],
					xNoisePotentialDerivatives[2] - zNoisePotentialDerivatives[0],
					yNoisePotentialDerivatives[0] - xNoisePotentialDerivatives[1]
					) * 0.07500000 * noiseVelScale * 0.5;
				newVelX += noiseVelocity.x * DeltaTime;
				newVelY += noiseVelocity.y * DeltaTime;
			}
		}
		#endif

		float decay = 1.0 / (1.0 + Decay * DeltaTime);

		if(sampleCoord.x == 0)
		{
			newVelX *= 0.0;
		}
		if(sampleCoord.x == (kNumCells.x - 1))
		{
			newVelX *= 0.0;
		}
		if(sampleCoord.y == 0)
		{
			newVelY *= 0.0;
		}
		if(sampleCoord.y == (kNumCells.y - 1))
		{
			newVelY *= 0.0;
		}

		outX = newVelX * decay;
		outY = newVelY * decay;

	}`;
}

export function GetShaderSourceDiffuseVelocityPS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	uniform float DeltaTime;
	uniform float Decay;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		float a = DeltaTime * (Decay * 10.f);

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);
		
		float curVelX = texelFetch(VecFieldXTexture, sampleCoord, 0).r;
		float curVelY = texelFetch(VecFieldYTexture, sampleCoord, 0).r;

		float newVelX =
		(curVelX +
			(
				texelFetch(VecFieldXTexture, sampleCoord - ivec2(1, 0), 0).r
				+ texelFetch(VecFieldXTexture, sampleCoord + ivec2(1, 0), 0).r
				+ texelFetch(VecFieldXTexture, sampleCoord - ivec2(0, 1), 0).r
				+ texelFetch(VecFieldXTexture, sampleCoord + ivec2(0, 1), 0).r
			) * a
		) / (1.0 + 4.0 * a);

		float newVelY =
		(curVelY +
			(
				texelFetch(VecFieldYTexture, sampleCoord - ivec2(1, 0), 0).r
				+ texelFetch(VecFieldYTexture, sampleCoord + ivec2(1, 0), 0).r
				+ texelFetch(VecFieldYTexture, sampleCoord - ivec2(0, 1), 0).r
				+ texelFetch(VecFieldYTexture, sampleCoord + ivec2(0, 1), 0).r
			) * a
		) / (1.0 + 4.0 * a);

		outX = newVelX;
		outY = newVelY;

	}`;
}

function scAdvectionDefines(bReverseAdvection : boolean, bAlpha : boolean)
{
	const reverse = `#define ADVECTION_REVERSE ` + (bReverseAdvection ? `1` : `0`);
	const alpha = `#define ADVECTION_ALPHA ` + (bAlpha ? `1` : `0`);
	const dimtype = `#define DIM_TYPE ` + (bAlpha ? `float` : `vec3`); 

	return reverse + ` \n ` + alpha + ` \n ` + dimtype;
}



export function GetShaderSourceAdvectDensityPS(numCells : Vector2, bReverseAdvection : boolean, bAlpha : boolean) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	`+scAdvectionDefines(bReverseAdvection, bAlpha)+ /* glsl */`
	`+scDeclareCurlNoiseFuncs()+ /* glsl */`

	layout(location = 0) out DIM_TYPE outDensity;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;
	uniform sampler2D DensityTexture;

	uniform float DeltaTime;
	uniform float Time;
	uniform float Decay;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);
		
		float curVelX = textureLod(VecFieldXTexture, vsOutTexCoords, 0.0).r;
		float curVelY = textureLod(VecFieldYTexture, vsOutTexCoords, 0.0).r;

		//curl noise
		#if 0
		{
			vec2 inVelocity = vec2(curVelX, curVelY);
			vec3 noisePosition = vec3(1.0);
			noisePosition.xy = vec2(sampleCoord);
			noisePosition *= 0.1;
			float noiseTime = Time * 1.025000;
			//noiseTime = 1.0;
			noisePosition.z = Time * 0.1;

			vec4 xNoisePotentialDerivatives = vec4(0.0);
    		vec4 yNoisePotentialDerivatives = vec4(0.0);
    		vec4 zNoisePotentialDerivatives = vec4(0.0);
			/* float curl = 0.75 + min(1.0, length(inVelocity)); 
			float noiseVelScale = 1.0 + min(1.0, length(inVelocity)); */
			float curl = 0.5;
			float noiseVelScale = 0.5;
			noiseVelScale = 0.1 + min(1.0, length(inVelocity) * 2.0);
			//noiseVelScale *= 2.5;

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

			curVelX += noiseVelocity.x;
			curVelY += noiseVelocity.y;

			//LQ
			{
				xNoisePotentialDerivatives = vec4(0.0);
				yNoisePotentialDerivatives = vec4(0.0);
				zNoisePotentialDerivatives = vec4(0.0);
				noiseTime = Time * 0.25 * 0.5;
				//noisePosition *= 0.02;
				noisePosition *= 0.2;
				noiseVelScale = 0.2 + min(1.0, length(inVelocity) * 1.0);
				noiseVelScale *= 0.75;
				curl = 0.9;
				int i = 0;
				float scale = (1.0 / 2.0) * pow(2.0, float(i));
				float noiseScale = 1.0;
				xNoisePotentialDerivatives += simplexNoiseDerivatives(vec4(noisePosition * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				yNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(123.4, 129845.6, -1239.1)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				zNoisePotentialDerivatives += simplexNoiseDerivatives(vec4((noisePosition + vec3(-9519.0, 9051.0, -123.0)) * pow(2.0, float(i)), noiseTime)) * noiseScale * scale;
				noiseVelocity = vec3(
					zNoisePotentialDerivatives[1] - yNoisePotentialDerivatives[2],
					xNoisePotentialDerivatives[2] - zNoisePotentialDerivatives[0],
					yNoisePotentialDerivatives[0] - xNoisePotentialDerivatives[1]
					) * 0.07500000 * noiseVelScale * 0.5;
				curVelX += noiseVelocity.x;
				curVelY += noiseVelocity.y;
			}
		}
		#endif

		//vec4 curColor = textureLod(DensityTexture, vsOutTexCoords, 0.0).rgba;
		
		//convert coord to WS, advect, convert back to texture space

		#if ADVECTION_REVERSE
			vec2 coord = vsOutTexCoords + vec2(curVelX, curVelY) * DeltaTime;
		#else
			vec2 coord = vsOutTexCoords - vec2(curVelX, curVelY) * DeltaTime;
		#endif

		#if ADVECTION_ALPHA
		float newColor = textureLod(DensityTexture, coord, 0.0).r;
		#else
		vec3 newColor = textureLod(DensityTexture, coord, 0.0).rgb;
		#endif

		float dec = 1.0 / (1.0 + Decay * DeltaTime);

		const int borderOffset = 2;
		if(sampleCoord.x < borderOffset)
		{
			newColor *= 0.0;
		}
		if(sampleCoord.x >= (kNumCells.x - borderOffset))
		{
			newColor *= 0.0;
		}
		if(sampleCoord.y < borderOffset)
		{
			newColor *= 0.0;
		}
		if(sampleCoord.y >= (kNumCells.y - borderOffset))
		{
			newColor *= 0.0;
		}

		outDensity = newColor * dec;

	}`;
}

export function GetShaderSourceAdvectDensityMacKormackPS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outDensity;

	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;
	uniform sampler2D DensityTexture;
	uniform sampler2D DensityAdvectedTexture;
	uniform sampler2D DensityReversedTexture;

	uniform float DeltaTime;
	uniform float Decay;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		vec3 curDensity = texelFetch(DensityTexture, sampleCoord, 0).rgb;
		vec3 correctionDensity = texelFetch(DensityReversedTexture, sampleCoord, 0).rgb;
		vec3 newDensity = texelFetch(DensityAdvectedTexture, sampleCoord, 0).rgb;

		vec3 error = (curDensity - correctionDensity);
		newDensity = newDensity + 0.5 * error;

		float dec = 1.0 / (1.0 + Decay * DeltaTime);

		const int borderOffset = 2;
		if(sampleCoord.x < borderOffset)
		{
			newDensity *= 0.0;
		}
		if(sampleCoord.x >= (kNumCells.x - borderOffset))
		{
			newDensity *= 0.0;
		}
		if(sampleCoord.y < borderOffset)
		{
			newDensity *= 0.0;
		}
		if(sampleCoord.y >= (kNumCells.y - borderOffset))
		{
			newDensity *= 0.0;
		}

	#if 1//CLAMP
		float curVelX = textureLod(VecFieldXTexture, vsOutTexCoords, 0.0).r;
		float curVelY = textureLod(VecFieldYTexture, vsOutTexCoords, 0.0).r;

		vec2 coord = vsOutTexCoords - vec2(curVelX, curVelY) * DeltaTime;

		const float mO = 0.1;

		//red
		{
			vec4 comp = vec4(
				textureOffset(DensityTexture, coord, ivec2(0,0), 0.0).r,
				textureOffset(DensityTexture, coord, ivec2(1,0), 0.0).r,
				textureOffset(DensityTexture, coord, ivec2(0,1), 0.0).r,
				textureOffset(DensityTexture, coord, ivec2(1,1), 0.0).r
			);
			vec2 minMax;
			minMax.x = min(comp.x, min(comp.y, min(comp.z, comp.w)));
			minMax.y = max(comp.x, max(comp.y, max(comp.z, comp.w)));
			newDensity.r = clamp(newDensity.r, minMax.x - mO, minMax.y + mO);
		}
		//green
		{
			vec4 comp = vec4(
				textureOffset(DensityTexture, coord, ivec2(0,0), 0.0).g,
				textureOffset(DensityTexture, coord, ivec2(1,0), 0.0).g,
				textureOffset(DensityTexture, coord, ivec2(0,1), 0.0).g,
				textureOffset(DensityTexture, coord, ivec2(1,1), 0.0).g
			);
			vec2 minMax;
			minMax.x = min(comp.x, min(comp.y, min(comp.z, comp.w)));
			minMax.y = max(comp.x, max(comp.y, max(comp.z, comp.w)));
			newDensity.g = clamp(newDensity.g, minMax.x - mO, minMax.y + mO);
		}
		//blue
		{
			vec4 comp = vec4(
				textureOffset(DensityTexture, coord, ivec2(0,0), 0.0).b,
				textureOffset(DensityTexture, coord, ivec2(1,0), 0.0).b,
				textureOffset(DensityTexture, coord, ivec2(0,1), 0.0).b,
				textureOffset(DensityTexture, coord, ivec2(1,1), 0.0).b
			);
			vec2 minMax;
			minMax.x = min(comp.x, min(comp.y, min(comp.z, comp.w)));
			minMax.y = max(comp.x, max(comp.y, max(comp.z, comp.w)));
			newDensity.b = clamp(newDensity.b, minMax.x - mO, minMax.y + mO);
		}
	#endif
		
		outDensity = newDensity * dec;
		

	}`;
}




export function GetShaderSourceResolveDivergencePS(numCells : Vector2) {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform sampler2D DivergenceTexture;
	uniform sampler2D VecFieldXTexture;
	uniform sampler2D VecFieldYTexture;

	in vec2 vsOutTexCoords;

	void main()
	{
		const ivec2 kNumCells = ivec2(int(` +
        numCells.x +
        /* glsl */ `), int(` +
        numCells.y +
        /* glsl */ `));

		ivec2 sampleCoord = ivec2(gl_FragCoord.xy);

		//if solid cell - quit
		if(sampleCoord.x == 0 || sampleCoord.x == (kNumCells.x - 1) || sampleCoord.y == 0 || sampleCoord.y == (kNumCells.y - 1))
		{
			discard;
		}

		//for X

		float curVelX = texelFetch(VecFieldXTexture, sampleCoord, 0).r;

		//add right Div, negate Div from left
		float divRight = texelFetch(DivergenceTexture, sampleCoord, 0).r;
		float divLeft = texelFetch(DivergenceTexture, sampleCoord - ivec2(1, 0), 0).r;
		
		outX = curVelX + divRight - divLeft;

		//for Y

		float curVelY = texelFetch(VecFieldYTexture, sampleCoord, 0).r;

		//add up Div, negate Div from down
		float divUp = texelFetch(DivergenceTexture, sampleCoord, 0).r;
		float divDown = texelFetch(DivergenceTexture, sampleCoord - ivec2(0, 1), 0).r;
		
		outY = curVelY + divUp - divDown;
		

	}`;
}




export function GetShaderSourceApplyVelocityVS() {
    return (
        /* glsl */ `#version 300 es

	precision highp float;

	layout(location = 0) in vec2 VertexBuffer;

	uniform float ScreenRatio;

	uniform vec2 Position;
	uniform float Radius;

	out vec2 vsOutTexCoords;

	void main()
	{
		vec2 pos = VertexBuffer.xy;
		vec2 scale = vec2(1.0);
		scale *= Radius;
		scale.x /= ScreenRatio;
		
		//calculate offset
		vec2 posOffset = Position;
		
		pos = (pos.xy * scale.xy) + posOffset.xy;

		//pos.x /= ScreenRatio;
		
		gl_Position = vec4(pos.xy, 0.0, 1.0);
		vsOutTexCoords = (VertexBuffer.xy + 1.0) * 0.5; // Convert to [0, 1] range
	}`
    );
}

export function GetShaderSourceApplyVelocityPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform vec2 Velocity;

	in vec2 vsOutTexCoords;

	void main()
	{
		float s = 1.0;
		float d = length(vsOutTexCoords - vec2(0.5));
		s = clamp(1.0 - (d * 2.0), 0.0, 1.0);
		/* if(d > 0.5)
		{
			s = 0.0;
		} */
		outX = Velocity.x * s;
		outY = Velocity.y * s;
	}`;
}

export function GetShaderSourceInjectDensityPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	uniform vec3 ColorFilter;

	in vec2 vsOutTexCoords;

	void main()
	{
		vec2 flippedUVs = vec2(vsOutTexCoords.x, 1.f - vsOutTexCoords.y);

		vec3 color = ColorFilter;

		float s = 1.0;
		float d = length(vsOutTexCoords - vec2(0.5));
		s = clamp(1.0 - (d * 2.0), 0.0, 1.0);
		if(s < 0.25)
		{
			discard;
		}

		outColor = color;
	}`;
}

export function GetShaderSourceInjectImpulsePS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float outX;
	layout(location = 1) out float outY;

	uniform float ImpulseStrength;

	in vec2 vsOutTexCoords;

	void main()
	{
		float s = 1.0;
		vec2 dir = (vsOutTexCoords - vec2(0.5));
		float d = length(dir);
		s = clamp(1.0 - (d * 2.0), 0.0, 1.0);
		/* if(d > 0.5)
		{
			s = 0.0;
		} */
		outX = dir.x * s * ImpulseStrength;
		outY = dir.y * s * ImpulseStrength;
	}`;
}





export function GetShaderSourceInjectDensityColorVS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	
	layout(location = 0) in vec2 VertexBuffer;

	uniform float ViewportRatio;
	uniform vec2 PositionNDC;
	uniform vec2 Scale;

	out vec2 vsOutTexCoords;


	void main()
	{
		vec2 pos = vec3(VertexBuffer.xy, 0.0f).xy;
		pos.xy *= Scale;

		pos.x /= ViewportRatio;

		pos += PositionNDC;


		gl_Position = vec4(pos.xy, 1.0, 1.0);
		vsOutTexCoords = (VertexBuffer.xy + 1.0) * 0.5; // Convert to [0, 1] range
	}`;
}


export function GetShaderSourceInjectDensityTexturedPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	uniform sampler2D ColorTexture;
	uniform vec3 ColorFilter;

	in vec2 vsOutTexCoords;

	void main()
	{
		vec2 flippedUVs = vec2(vsOutTexCoords.x, 1.f - vsOutTexCoords.y);

		vec3 color = texture(ColorTexture, flippedUVs.xy).rgb;

		//color = pow(color, vec3(2.0)) * 5.0;
		color *= ColorFilter;
		color *= 1.5;

		outColor = color;
		//outColor = vec3(1.0, 0.0, 0.5);
	}`;
}

export function GetShaderSourceInjectDensitySceneColorPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec3 outColor;

	uniform sampler2D ColorTexture;
	uniform vec3 ColorFilter;
	uniform vec2 ViewportSize;

	in vec2 vsOutTexCoords;

	void main()
	{

		vec2 sampleCoordf = vec2(gl_FragCoord.xy) + vec2(0.5, 0.5);
		vec2 screenUV = sampleCoordf / ViewportSize;

		vec3 color = textureLod(ColorTexture, screenUV.xy, 0.0).rgb;

		// Calculate the dispersion effect
		vec3 dispersionColor;
		vec2 ndc = screenUV.xy * 2.0 - 1.0;
		ndc /= length(ndc) + 0.0001;
		const float DispOffsetScale = 1.0;
		vec2 dispOffset = ndc * DispOffsetScale * 0.01;

		dispersionColor.r = (textureLod(ColorTexture, screenUV.xy - dispOffset, 0.f).r);
		dispersionColor.g = (textureLod(ColorTexture, screenUV.xy, 0.f).g);
		dispersionColor.b = (textureLod(ColorTexture, screenUV.xy + dispOffset, 0.f).b);

		color = dispersionColor;

		

		//color = pow(color, vec3(2.0)) * 5.0;
		/* color *= ColorFilter;
		color *= 1.5; */

		outColor = color;
		outColor = vec3(1.0);
		//outColor = vec3(1.0, 0.0, 0.5);
	}`;
}


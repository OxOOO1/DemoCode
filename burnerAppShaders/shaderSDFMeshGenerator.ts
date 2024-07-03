import { scDeclareSceneUniformBuffer } from "./shaderBasePass";
import { scGetCommonFuncsDeclarations } from "./shaderCommon";

//Procedurally Render To Volume Texture 
export function GetShaderSourceVolumeTextureRenderPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float Output;

	uniform int LayerIndex;

	in vec2 vsOutTexCoords;

	`+scGetCommonFuncsDeclarations()+ /* glsl */`

	void main()
	{
		float meshSize = 1.0;
        vec3 meshExtent = vec3(meshSize, meshSize, meshSize);
        vec3 meshMin = vec3(-meshSize, -meshSize, -meshSize);
        float texPrec = 64.0;
        vec3 textureSize3D = vec3(texPrec, texPrec, 4);
        vec3 voxelSize = vec3(0, 0, 0);
        voxelSize.x = (meshExtent.x * 2.0) / textureSize3D.x;
        voxelSize.y = (meshExtent.y * 2.0) / textureSize3D.y;
        voxelSize.z = (meshExtent.z * 2.0) / textureSize3D.z;

		float sphereRadius = 0.5;

		ivec3 id = ivec3(int(gl_FragCoord.x), int(gl_FragCoord.y), LayerIndex);

		vec3 curVoxelPos;
		curVoxelPos.x = meshMin.x + voxelSize.x * float(id.x) + voxelSize.x * 0.5;
        curVoxelPos.y = meshMin.y + voxelSize.y * float(id.y) + voxelSize.y * 0.5;
        curVoxelPos.z = meshMin.z + voxelSize.z * float(id.z) + voxelSize.z * 0.5;

		float sdf = length(curVoxelPos) - sphereRadius;

		/* float val = 0.0;
		if(sdf < 0.1)
		{
			val = 1.0;
		} */

		Output = sdf;

	}`;
}



////////////////////////////Mesh Visualizer

export function GetShaderSourceMeshVisualizerVS() {
    return /* glsl */ `#version 300 es

		precision highp float;
		precision highp sampler2D;

		` +scDeclareSceneUniformBuffer()+/* glsl */`

		uniform sampler2D VertexBufferTexture; 

		out vec3 interpolatorColor;

		float hash( uint n ) 
		{
			n = (n << 13U) ^ n;
	    	n = n * (n * n * 15731U + 789221U) + 1376312589U;
	    	// floating point conversion from https://iquilezles.org/articles/sfrand
	    	return uintBitsToFloat( (n>>9U) | 0x3f800000U ) - 1.0;
		}

		vec3 hsv2rgb(vec3 hsv) {
			float h = hsv.x;
			float s = hsv.y;
			float v = hsv.z;
		
			float c = v * s;
			float x = c * (1.0 - abs(mod(h * 6.0, 2.0) - 1.0));
			float m = v - c;
		
			vec3 rgb;
		
			if (h < 1.0/6.0) {
				rgb = vec3(c, x, 0.0);
			} else if (h < 2.0/6.0) {
				rgb = vec3(x, c, 0.0);
			} else if (h < 3.0/6.0) {
				rgb = vec3(0.0, c, x);
			} else if (h < 4.0/6.0) {
				rgb = vec3(0.0, x, c);
			} else if (h < 5.0/6.0) {
				rgb = vec3(x, 0.0, c);
			} else {
				rgb = vec3(c, 0.0, x);
			}
		
			rgb = rgb + vec3(m);
		
			return rgb;
		}

		void main()
		{
			int vId = int(gl_VertexID) * 3;
			int offset = vId % 450;
			int row = int(floor(float(vId) / 450.0));

			vec3 pos = vec3(
				texelFetch(VertexBufferTexture, ivec2(offset + 0, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 1, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 2, row), 0).r
			);

			//Camera Space
			//pos.xyz = (vec4(pos, 1.0) * WorldToViewMatrix).xyz;
			pos.xyz -= SceneBuffer.CameraPosition; 

			//Projection
			pos.xy *= SceneBuffer.CameraZoom;
			pos.x /= SceneBuffer.ScreenRatio;

			gl_Position = vec4(pos.xy, pos.z / 20.0, (0.1f + pos.z));

			float rnd = hash((uint(gl_VertexID) / 3u));

			vec3 hsv = vec3(rnd, 1.0, 1.0);

			vec3 c = hsv2rgb(hsv);

			interpolatorColor = c;
		}`;
}

export function GetShaderSourceMeshVisualizerPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out vec4 outColor;

	in vec3 interpolatorColor;

	void main()
	{
		outColor = vec4(interpolatorColor, 1.0);
	}`;
}










export function GetShaderSourceMeshToSDFConverterPS() {
    return /* glsl */ `#version 300 es
	
	precision highp float;
	precision highp sampler2D;

	layout(location = 0) out float Output;

	uniform int LayerIndex;
	uniform int NumTriangles;
	uniform ivec3 SDFTextureSize;
	uniform vec3 SDFMeshExtent;

	uniform sampler2D VertexBufferTexture; 

	in vec2 vsOutTexCoords;

	`+scGetCommonFuncsDeclarations()+ /* glsl */`

	float dot2( in vec3 v ) { return dot(v,v); }

	float udTriangle( in vec3 v1, in vec3 v2, in vec3 v3, in vec3 p )
	{
	    // prepare data    
	    vec3 v21 = v2 - v1; vec3 p1 = p - v1;
	    vec3 v32 = v3 - v2; vec3 p2 = p - v2;
	    vec3 v13 = v1 - v3; vec3 p3 = p - v3;
	    vec3 nor = cross( v21, v13 );

		//float signEx = dot(p1, nor) < 0.0 ? 1.0 : -1.0;
		const float signEx = 1.0;
	
	    return ( signEx * sqrt( // inside/outside test    
	                 (sign(dot(cross(v21,nor),p1)) + 
	                  sign(dot(cross(v32,nor),p2)) + 
	                  sign(dot(cross(v13,nor),p3))<2.0) 
	                  ?
	                  // 3 edges    
	                  min( min( 
	                  dot2(v21*clamp(dot(v21,p1)/dot2(v21),0.0,1.0)-p1), 
	                  dot2(v32*clamp(dot(v32,p2)/dot2(v32),0.0,1.0)-p2) ), 
	                  dot2(v13*clamp(dot(v13,p3)/dot2(v13),0.0,1.0)-p3) )
	                  :
	                  // 1 face    
	                  dot(nor,p1)*dot(nor,p1)/dot2(nor) ));
	}

	struct Triangle {
		vec3 v0;
		vec3 v1;
		vec3 v2;
	};

	bool isTriangleFrontFacing(vec3 camPos, vec3 v1, vec3 v2, vec3 v3) {
	#if 0
		// Calculate the vectors for two edges of the triangle
		vec3 edge1 = v2 - v1;
		vec3 edge2 = v3 - v1;
	
		// Compute the normal of the triangle using the cross product
		vec3 normal = cross(edge1, edge2);
	
		// Compute the view vector from the camera to the triangle's first vertex
		vec3 viewVector = v1 - camPos;
	
		// Calculate the dot product of the normal and view vector
		float dotProduct = dot(normal, viewVector);
	
		// If the dot product is positive, the triangle is front facing; otherwise, it's back facing
		return dotProduct > 0.0;
	#else
		vec3 edge1 = v2 - v1;
		vec3 edge2 = v3 - v1;

		vec3 normal = cross(edge1, edge2);

		vec3 viewVector = camPos - v1;

		float dotProduct = dot(normal, viewVector);

		return (dotProduct > 0.0);
	#endif
	}

	void main()
	{
		/* const float meshSize = 1.0;
        const vec3 meshExtent = vec3(meshSize, meshSize, meshSize); */
		//vec3 meshExtent = vec3(1.0, 0.5, 1.0);
		vec3 meshExtent = SDFMeshExtent;
		vec3 meshMin = -meshExtent;
        //const vec3 meshMin = vec3(-meshSize, -meshSize, -meshSize);
        vec3 voxelSize = vec3(0, 0, 0);
        voxelSize.x = (meshExtent.x * 2.0) / float(SDFTextureSize.x);
        voxelSize.y = (meshExtent.y * 2.0) / float(SDFTextureSize.y);
        voxelSize.z = (meshExtent.z * 2.0) / float(SDFTextureSize.z);

		ivec3 id = ivec3(int(gl_FragCoord.x), int(gl_FragCoord.y), LayerIndex);

		vec3 curVoxelPos;
		curVoxelPos.x = meshMin.x + voxelSize.x * float(id.x) + voxelSize.x * 0.5;
        curVoxelPos.y = meshMin.y + voxelSize.y * float(id.y) + voxelSize.y * 0.5;
        curVoxelPos.z = meshMin.z + voxelSize.z * float(id.z) + voxelSize.z * 0.5;

		float minDist = 100.0; 

		int insideTriCount = 0;

		Triangle closestTriangle;

		for (int i = 0; i < NumTriangles; i++) {

			//int offset = i * 9;
			int vId = i * 9;
			int offset = vId % 450;
			int row = int(floor(float(vId) / 450.0));

			vec3 v0 = vec3(
				texelFetch(VertexBufferTexture, ivec2(offset + 0, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 1, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 2, row), 0).r
			);
			vec3 v1 = vec3(
				texelFetch(VertexBufferTexture, ivec2(offset + 3, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 4, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 5, row), 0).r
			);
			vec3 v2 = vec3(
				texelFetch(VertexBufferTexture, ivec2(offset + 6, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 7, row), 0).r,
				texelFetch(VertexBufferTexture, ivec2(offset + 8, row), 0).r
			);

			float dist = (udTriangle(v0, v1, v2, curVoxelPos));

			/* bool bFrontFace = isTriangleFrontFacing(curVoxelPos, v0, v1, v2);
			if(!bFrontFace)
			{
				insideTriCount += 1;
			} */

			if(dist < minDist)
			{
				minDist = dist;

				closestTriangle.v0 = v0;
				closestTriangle.v1 = v1;
				closestTriangle.v2 = v2;
			}
		}

		bool bFrontFace = isTriangleFrontFacing(curVoxelPos, closestTriangle.v0, closestTriangle.v1, closestTriangle.v2);

		#if 1
		if(!bFrontFace)
		{
			minDist *= -1.0;
		}
		else
		{
			//minDist -= 0.02;
		}
		#else
			//minDist -= 0.02;
			minDist -= 0.0125;
		#endif

		/* {
			float insideTriThres = 0.6;
			if(!bFrontFace && (insideTriCount > int(floor(float(NumTriangles) * insideTriThres))))
			{
				minDist *= -1.0;
			}
		} */


		/* {
			vec3 v0 = vec3(-0.5, -0.8, 0.2);
			vec3 v1 = vec3(0.0, 0.8, 0.0);
			vec3 v2 = vec3(0.9, -0.8, -0.2);

			Output = udTriangle(v0, v1, v2, curVoxelPos); return;
		} */

		Output = minDist;

	}`;
}


import { Vector3 } from "../core/types";

export function scGetVec3Decl(vec : Vector3)
{
	return /* glsl */`vec3(float(` +
	vec.x +
	/* glsl */ `), float(` +
	vec.y +
	/* glsl */ `),
	float(` +
	vec.z +
	/* glsl */ `)
	)`
}


export function scGetCommonFuncsDeclarations()
{
	return (/* glsl */`
	

	#define PI 3.141592654f
	#define TWOPI 6.283185308f
	#define PIDIVTWO 1.570796327f


	/* vec3 TransformNDCToWorld(vec2 posNDC, float viewPosZ)
	{
		vec3 posWS = vec3(posNDC.x * SceneBuffer.ScreenRatio, posNDC.y, 0.0);
		posWS.xy /= SceneBuffer.CameraZoom;

		const float nearPlane = 0.1;
		posWS.xy *= nearPlane + viewPosZ;
		posWS.xy += SceneBuffer.CameraPosition.xy;
		posWS.z = viewPosZ + SceneBuffer.CameraPosition.z + nearPlane;
		return posWS;
	} */

	mat3 MatrixRotateAroundAxis(vec3 axis, float angle) {
		vec3 n = /* normalize */(axis);
		float cosA = cos(angle);
		float sinA = sin(angle);
		float oneMinusCosA = 1.0 - cosA;
	
		mat3 mat = mat3(1.0); // Initialize as identity matrix
	
		mat[0][0] = (n.x * n.x) * oneMinusCosA + cosA;
		mat[0][1] = (n.x * n.y) * oneMinusCosA + n.z * sinA;
		mat[0][2] = (n.x * n.z) * oneMinusCosA - n.y * sinA;
	
		mat[1][0] = (n.x * n.y) * oneMinusCosA - n.z * sinA;
		mat[1][1] = (n.y * n.y) * oneMinusCosA + cosA;
		mat[1][2] = (n.y * n.z) * oneMinusCosA + n.x * sinA;
	
		mat[2][0] = (n.x * n.z) * oneMinusCosA + n.y * sinA;
		mat[2][1] = (n.y * n.z) * oneMinusCosA - n.x * sinA;
		mat[2][2] = (n.z * n.z) * oneMinusCosA + cosA;
	
		return mat;
	}


	//[0, 2]
	vec3 Contrast(vec3 color, float value) 
	{
		return vec3(0.5) + value * (color - vec3(0.5));
	}

	vec3 Saturation(vec3 color, float value)
	{
		vec3 weights_ = vec3(0.2125, 0.7154, 0.0721); // sums to 1
        float luminance_ = dot(color.rgb, weights_);
        return mix(vec3(luminance_), color, value);
	}
	
	float MapToRange(float t, float t0, float t1, float newt0, float newt1)
	{
		///Translate to origin, scale by ranges ratio, translate to new position
		return (t - t0) * ((newt1 - newt0) / (t1 - t0)) + newt0;
	}

	vec2 MapToRange(vec2 uv, float t0, float t1, float newt0, float newt1)
	{
		uv.x = MapToRange(uv.x, t0, t1, newt0, newt1);
		uv.y = MapToRange(uv.y, t0, t1, newt0, newt1);
		return uv;
	}

	vec3 MapToRange(vec3 uv, vec3 t0, vec3 t1, vec3 newt0, vec3 newt1)
	{
		vec3 uvw;
		uvw.x = MapToRange(uv.x, t0.x, t1.x, newt0.x, newt1.x);
		uvw.y = MapToRange(uv.y, t0.y, t1.y, newt0.y, newt1.y);
		uvw.z = MapToRange(uv.z, t0.z, t1.z, newt0.z, newt1.z);
		return uvw;
	}

	float hash( uint n ) 
	{
		n = (n << 13U) ^ n;
	    n = n * (n * n * 15731U + 789221U) + 1376312589U;
	    // floating point conversion from https://iquilezles.org/articles/sfrand
	    return uintBitsToFloat( (n>>9U) | 0x3f800000U ) - 1.0;
	}

	float hashf( float n ){return fract(sin(n)*43758.5453);}
	

	vec2 hash2( float n ) { return fract(sin(vec2(n,n+1.0))*vec2(43758.5453123,22578.1459123)); }

	#define Scale (vec3(0.8, 0.8, 0.8))
	#define K     (19.19)
	vec3 hash3(vec3 a) 
	{
    	a = fract(a * Scale);
    	a += dot(a, a.yxz + K);
    	return fract((a.xxy + a.yxx) * a.zyx);
	}

	// === Voronoi =====================================================
	// --- Base Voronoi. inspired by https://www.shadertoy.com/view/MslGD8

	#define hash22(p)  fract( 18.5453 * sin( p * mat2(127.1,311.7,269.5,183.3)) )
	#define VARIANT 0 
	const float ofs = 0.;
	#define disp(p) ( -ofs + (1.+2.*ofs) * hash22(p) )

	vec3 voronoi( vec2 u )  // returns len + id
	{
	    vec2 iu = floor(u), v;
		float m = 1e9,d;
	#if VARIANT
	    for( int k=0; k < 25; k++ ) {
	        vec2  p = iu + vec2(k%5-2,k/5-2),
	#else
	    for( int k=0; k < 9; k++ ) {
	        vec2  p = iu + vec2(k%3-1,k/3-1),
	#endif
	            o = disp(p),
	      	      r = p - u + o;
			d = dot(r,r);
	        if( d < m ) m = d, v = r;
	    }

	    return vec3( sqrt(m), v+u );
	}

	// --- Voronoi distance to borders. inspired by https://www.shadertoy.com/view/ldl3W8
	vec3 voronoiB( vec2 u )  // returns len + id
	{
	    vec2 iu = floor(u), C, P;
		float m = 1e9,d;
	#if VARIANT
	    for( int k=0; k < 25; k++ ) {
	        vec2  p = iu + vec2(k%5-2,k/5-2),
	#else
	    for( int k=0; k < 9; k++ ) {
	        vec2  p = iu + vec2(k%3-1,k/3-1),
	#endif
	              o = disp(p),
	      	      r = p - u + o;
			d = dot(r,r);
	        if( d < m ) m = d, C = p-iu, P = r;
	    }

	    m = 1e9;
	
	    for( int k=0; k < 25; k++ ) {
	        vec2 p = iu+C + vec2(k%5-2,k/5-2),
			     o = disp(p),
	             r = p-u + o;

	        if( dot(P-r,P-r)>1e-5 )
	        m = min( m, .5*dot( (P+r), normalize(r-P) ) );
	    }

	    return vec3( m, P+u );
	}


	// === pseudo Perlin noise =============================================
	#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
	int MOD = 1;  // type of Perlin noise
	
	// --- 2D
	#define hash21(p) fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453123)
	float noise2(vec2 p) {
	    vec2 i = floor(p);
	    vec2 f = fract(p); f = f*f*(3.-2.*f); // smoothstep

	    float v= mix( mix(hash21(i+vec2(0,0)),hash21(i+vec2(1,0)),f.x),
	                  mix(hash21(i+vec2(0,1)),hash21(i+vec2(1,1)),f.x), f.y);
		return   MOD==0 ? v
		       : MOD==1 ? 2.*v-1.
	           : MOD==2 ? abs(2.*v-1.)
	                    : 1.-abs(2.*v-1.);
	}

	float fbm2(vec2 p) {
	    float v = 0.,  a = .5;
	    mat2 R = rot(.37);

	    for (int i = 0; i < 9; i++, p*=2.,a/=2.) 
	        p *= R,
	        v += a * noise2(p);

	    return v;
	}
	#define noise22(p) vec2(noise2(p),noise2(p+17.7))
	vec2 fbm22(vec2 p) {
	    vec2 v = vec2(0);
	    float a = .5;
	    mat2 R = rot(.37);

	    for (int i = 0; i < 9; i++, p*=2.,a/=2.) 
	        p *= R,
	        v += a * noise22(p);

	    return v;
	}

	#if 0
	vec4 mod289(vec4 x){return x - floor(x * (1.0 / 289.0)) * 289.0;}
	vec4 perm(vec4 x){return mod289(((x * 34.0) + 1.0) * x);}

	float smoothNoise(vec3 p)
	{
    	vec3 a = floor(p);
    	vec3 d = p - a;
    	d = d * d * (3.0 - 2.0 * d);

    	vec4 b = a.xxyy + vec4(0.0, 1.0, 0.0, 1.0);
    	vec4 k1 = perm(b.xyxy);
    	vec4 k2 = perm(k1.xyxy + b.zzww);

    	vec4 c = k2 + a.zzzz;
    	vec4 k3 = perm(c);
    	vec4 k4 = perm(c + 1.0);

    	vec4 o1 = fract(k3 * (1.0 / 41.0));
    	vec4 o2 = fract(k4 * (1.0 / 41.0));

    	vec4 o3 = o2 * d.z + o1 * (1.0 - d.z);
    	vec2 o4 = o3.yw * d.x + o3.xz * (1.0 - d.x);

    	return o4.y * d.y + o4.x * (1.0 - d.y);
	}
	#endif

	float hashNoise (in vec3 x)
	{
		vec3 p = floor(x);
		vec3 f = fract(x);

		f = f*f*(3.0-2.0*f);

		float n = p.x + p.y*57.0 + 113.0*p.z;

		float res = mix(mix(mix( hashf(n+  0.0), hashf(n+  1.0),f.x),
							mix( hashf(n+ 57.0), hashf(n+ 58.0),f.x),f.y),
						mix(mix( hashf(n+113.0), hashf(n+114.0),f.x),
							mix( hashf(n+170.0), hashf(n+171.0),f.x),f.y),f.z);
		return res;
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

	//x:bool sign, y::distance, zw:tExit tEnter
	vec4 RayAABBIntersection(vec3 rayOrigin, vec3 rayDirection, vec3 aabbCenter, vec3 aabbExtent)
	{
		// Calculate the inverse direction of the ray
		vec3 invDirection;
		invDirection.x = 1.0f / rayDirection.x;
		invDirection.y = 1.0f / rayDirection.y;
		invDirection.z = 1.0f / rayDirection.z;
	
		// Calculate the minimum and maximum t values for each axis
		float tmin = (aabbCenter.x - aabbExtent.x - rayOrigin.x) * invDirection.x;
		float tmax = (aabbCenter.x + aabbExtent.x - rayOrigin.x) * invDirection.x;
		float tymin = (aabbCenter.y - aabbExtent.y - rayOrigin.y) * invDirection.y;
		float tymax = (aabbCenter.y + aabbExtent.y - rayOrigin.y) * invDirection.y;
		float tzmin = (aabbCenter.z - aabbExtent.z - rayOrigin.z) * invDirection.z;
		float tzmax = (aabbCenter.z + aabbExtent.z - rayOrigin.z) * invDirection.z;
	
		// Find the maximum and minimum t values for intersection
		float tEnter = max(max(min(tmin, tmax), min(tymin, tymax)), min(tzmin, tzmax));
		float tExit = min(min(max(tmin, tmax), max(tymin, tymax)), max(tzmin, tzmax));
		
		vec4 result;//x:bool, y::distance, xz:tExit tEnter
		// Check if the intersection is outside the valid range
		if ((tExit < 0.f) || (tEnter > tExit))
		{
			result.x = -1.f;
		}
		else
		{
			result.x = 1.f;
		}
	
		// Set the intersection distance
		result.y = (tEnter >= 0.f) ? tEnter : tExit;
		result.z = tEnter;
		result.w = tExit;
	
		return result;
	}

	bool RaySphereIntersection(vec3 rayOrigin, vec3 rayDir, vec3 spherePos, float sphereRadius)
	{
		vec3 m = rayOrigin - spherePos;
		float b = dot(m, rayDir);
		float c = dot(m, m) - sphereRadius * sphereRadius;
		if (c > 0.0f && b > 0.0f)
		{
			return false;
		}
		float discr = b * b - c;
		if (discr < 0.0f)
		{
			return false;
		}
		return true;
	}

	vec3 randomSphereDir(vec2 rnd)
	{
		float s = rnd.x * PI * 2.;
		float t = rnd.y * 2. - 1.;
		return vec3(sin(s), cos(s), t) / sqrt(1.0 + t * t);
	}

	vec3 randomHemisphereDir(vec3 dir, float i)
	{
		vec3 v = randomSphereDir( vec2(hashf(i+1.), hashf(i+2.)) );
		return v * sign(dot(v, dir));
	}
	
	`);
}

export function scGetSphericalFuncsDeclarations()
{
	return (/* glsl */`
	

	#define PI 3.141592654f
	#define TWOPI 6.283185308f
	#define PIDIVTWO 1.570796327f
	
	

	const vec2 invAtan = vec2(0.1591, 0.3183);
	const vec2 kAtan = vec2(TWOPI, PI);
	vec2 UVToSpherical(vec2 uv)
	{
		vec2 res = uv;
		res -= 0.5;
		res *= kAtan;
		return res;
	}

	vec2 SphericalToUV(vec2 phiTheta)
	{
		vec2 uv = phiTheta;
	    uv *= invAtan;//normalize
	    uv += 0.5;
	    return uv;
	}

	vec2 DirToUV(vec3 v)
	{
	    vec2 phiTheta = vec2(atan(v.z, v.x), asin(v.y));
	    return SphericalToUV(phiTheta);
	}

	//Default formula
	vec3 SphericalToCartesian(float phi, float theta)
	{
		vec3 res;
		res.x = cos(phi) * sin(theta); 
		res.y = sin(theta) * sin(phi);
		res.z = cos(theta);
		return res;
	}

	vec2 SphericalToUVMod(vec2 phiTheta)
	{
		vec3 dir = SphericalToCartesian(phiTheta.x, phiTheta.y);
		return DirToUV(dir);
		/* vec2 uv = phiTheta;
		uv.x = mod(uv.x + PI, 2.0 * PI) - PI;
		uv.y = asin(sin(uv.y));
	    uv *= invAtan;//normalize
	    uv += 0.5;
	    return uv; */
	}

	vec2 CartesianToSpherical(vec3 v)
	{
		return vec2(atan(v.z, v.x), asin(v.y));
	}

	float[4] SHBasis1(vec3 dir)
	{
		float x = dir.x;
		float y = dir.y;
		float z = dir.z;

		float z2 = z*z;

		float coefArr[4];

        /* m=0 */

        // l=0
        const float p_0_0 = (0.282094791773878140);
        coefArr[0] = p_0_0; // l=0,m=0
        // l=1
        float p_1_0 = (0.488602511902919920)*z;
        coefArr[2] = p_1_0; // l=1,m=0


        /* m=1 */

        float s1 = y;
        float c1 = x;

        // l=1
        const float p_1_1 = (-0.488602511902919920);
        coefArr[1] = p_1_1*s1; // l=1,m=-1
        coefArr[3] = p_1_1*c1; // l=1,m=+1


		return coefArr;
	}

	float[9] SHBasis2(vec3 dir)
	{
		float x = dir.x;
		float y = dir.y;
		float z = dir.z;

		float z2 = z*z;

		float coefArr[9];

        /* m=0 */

        // l=0
        const float p_0_0 = (0.282094791773878140);
        coefArr[0] = p_0_0; // l=0,m=0
        // l=1
        float p_1_0 = (0.488602511902919920)*z;
        coefArr[2] = p_1_0; // l=1,m=0
        // l=2
        float p_2_0 = (0.946174695757560080)*z2 + (-0.315391565252520050);
        coefArr[6] = p_2_0; // l=2,m=0


        /* m=1 */

        float s1 = y;
        float c1 = x;

        // l=1
        const float p_1_1 = (-0.488602511902919920);
        coefArr[1] = p_1_1*s1; // l=1,m=-1
        coefArr[3] = p_1_1*c1; // l=1,m=+1
        // l=2
        float p_2_1 = (-1.092548430592079200)*z;
        coefArr[5] = p_2_1*s1; // l=2,m=-1
        coefArr[7] = p_2_1*c1; // l=2,m=+1


        /* m=2 */

        float s2 = x*s1 + y*c1;
        float c2 = x*c1 - y*s1;

        // l=2
        const float p_2_2 = (0.546274215296039590);
        coefArr[4] = p_2_2*s2; // l=2,m=-2
        coefArr[8] = p_2_2*c2; // l=2,m=+2

		return coefArr;
	}
	

	

	
	
	`);
}

export function scPBRDeclarations(bHalfLambert = false)
{
	return (/* glsl */`
	

	float GetSpotlightFalloff(vec3 lightPos, vec3 lightDir, vec3 pixelPos, float angleOut, float angleIn)
	{
		vec3 vToL = normalize(lightPos - pixelPos);
		float curSpotlightAngle = dot(-vToL, lightDir);
		vec2 Angles = vec2(cos(angleOut - angleIn), cos(angleOut));
		vec2 SpotlightAngles;
		SpotlightAngles.x = 1.0f / (Angles.x - Angles.y);
		SpotlightAngles.y = Angles.y;
		float SpotlightFalloff = clamp((curSpotlightAngle - SpotlightAngles.y) * SpotlightAngles.x, 0.0, 1.0);
		return SpotlightFalloff;
	}


	float RadicalInverse_VdC(uint bits) 
	{
	    bits = (bits << 16u) | (bits >> 16u);
	    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	    return float(bits) * 2.3283064365386963e-10; // / 0x100000000
	}
	vec2 Hammersley(uint i, uint N)
	{
	    return vec2(float(i)/float(N), RadicalInverse_VdC(i));
	}  

	//Xi - uniform distribtuion rnd value
	//N - tangent space
	vec3 ImportanceSampleGGX(vec2 Xi, vec3 N, float roughness)
	{
	    float a = roughness*roughness;

	    float phi = 2.0 * PI * Xi.x;
		//tilt halfway vec based on cur roughness and rnd num
		//smooth surfaces map rnd[0,1] to [1, 0.995] - so small spread, only few samples needed
		//half rough surfaces map rnd[0,1] to [0.99, 0.197] - noticeable increase in spread
		//rough surfaces map rnd[0,1] to [0.99, 0.117] - less noticeable increse in spread compared to half roughness
	    float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (a*a - 1.0) * Xi.y)); //Higher value - more aligned with normal
	    float sinTheta = sqrt(1.0 - cosTheta*cosTheta);

	    // from spherical coordinates to cartesian coordinates
	    vec3 H;
	    H.x = cos(phi) * sinTheta;
	    H.y = sin(phi) * sinTheta;
	    H.z = cosTheta;

	    // from tangent-space vector to world-space sample vector
	    vec3 up        = abs(N.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
	    vec3 tangent   = normalize(cross(up, N));
	    vec3 bitangent = cross(N, tangent);

	    vec3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
	    return normalize(sampleVec);
	}  
	
	vec3 DecodeNormalTexture(vec3 normTex, float angleScale)
	{
		//normTex.y = 1.0f - normTex.y;
		normTex.x = normTex.x * 2.0f - 1.0f;
	    normTex.y = normTex.y * 2.0f - 1.0f;
	    normTex.z = normTex.z * 2.0f - 1.0f;
		normTex.z *= angleScale;
		//normTex.z *= -1.f;
	    return normalize(normTex);
	}

	

	//When F90 is 1
	vec3 F_Schlick(float cosTheta, vec3 F0) {
		float f = pow(1.0 - cosTheta, 5.0);
		return f + F0 * (1.0 - f);
	}

	vec3 FresnelSchlick(float cosTheta, vec3 F0)
	{
		float F90 = SceneBuffer.FresnelMax;
		return F0 + (F90 - F0) * pow(1.0 - cosTheta, 5.0);
	}

	//F90 depends on the roughness
	vec3 FresnelSchlickRoughness(float cosTheta, vec3 F0, float roughness)
	{
		float F90 = SceneBuffer.FresnelMax;
	    return F0 + (max(vec3(F90 - roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
	}   

	float NormalDistributionGGX(vec3 N, vec3 H, float roughness)
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

	float NormalDistributionGGXFast(float NdotH, float a2)
	{
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
		//float k = (r*r) / 2.0;

		float num   = NdotV;
		float denom = NdotV * (1.0 - k) + k;

		return num / denom;
	}

	float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
	{
		float NdotV = max(dot(N, V), 0.0);
		float NdotL = max(dot(N, L), 0.0);
		//Occlusion from view
		float ggx2  = GeometrySchlickGGX(NdotV, roughness);
		//Shadowing by neighbor
		float ggx1  = GeometrySchlickGGX(NdotL, roughness);

		return ggx1 * ggx2;
	}

	float GeometrySchlickGGX_IBL(float NdotV, float roughness)
	{
		float r = (roughness);
		float k = (r*r) / 2.0;

		float num   = NdotV;
		float denom = NdotV * (1.0 - k) + k;

		return num / denom;
	}

	float GeometrySmith_IBL(vec3 N, vec3 V, vec3 L, float roughness)
	{
		float NdotV = max(dot(N, V), 0.0);
		float NdotL = max(dot(N, L), 0.0);
		//Occlusion from view
		float ggx2  = GeometrySchlickGGX_IBL(NdotV, roughness);
		//Shadowing by neighbor
		float ggx1  = GeometrySchlickGGX_IBL(NdotL, roughness);

		return ggx1 * ggx2;
	}

	struct PBRSettings
	{
		highp vec3 Normal;
		vec3 LightPos;
		vec3 LightColor;
		float LightRadius;
		vec3 Radiance;
		vec3 vToCam;
		highp vec3 PixelPos;
		vec3 Albedo;
		vec3 Ambient;
		highp float Roughness;
		float Metalness;
	};

	float GetLightAttenuation(float d, float r)
	{
		float attenuation = 1.0;
	#if 0//Higher shoulder
		attenuation *= pow(clamp(1.f - pow(d / r, 2.0), 0.f, 1.f), 2.0);
	#else //UE4
		{
			float num = pow(clamp(1.0 - pow(d / r, 4.0), 0.0, 1.0), 2.0);
			float denom =  (d * d + 1.0);
			attenuation = num / denom;
		}
	#endif
		return attenuation;
	}

	vec3 CalculateLightPBR(PBRSettings s)
	{
		vec3 v = s.vToCam;
		float d = length(s.LightPos - s.PixelPos);
		vec3 l = (s.LightPos - s.PixelPos) / d;
		vec3 h = normalize(v + l);

		//Attenuation
		float attenuation = 1.0;
	#if 0//Higher shoulder
		attenuation *= pow(clamp(1.f - pow(d / s.LightRadius, 2.0), 0.f, 1.f), 2.0);
	#else //UE4
		{
			float num = pow(clamp(1.0 - pow(d / s.LightRadius, 4.0), 0.0, 1.0), 2.0);
			float denom =  (d * d + 1.0);
			attenuation = num / denom;
		}
	#endif
		
		s.Radiance = s.LightColor * SceneBuffer.LightIntensity1 * attenuation;

		//calculate specular-diffuse ratio with Fresnel
		//vec3 to support chromatic reflection for metals
		vec3 F0 = vec3(0.04); //surface reflection at zero incidence
		//vec3 F0 = vec3(SceneBuffer.FresnelMin);
		//For Metals Albedo texture represents base reflectivity (F0) 
		F0 = mix(F0, s.Albedo, s.Metalness);
		vec3 F  = F_Schlick(max(dot(h, v), 0.0), F0);

		//normal distribution function
		float NDF = NormalDistributionGGX(s.Normal, h, s.Roughness);  

		//geometry function     
		float G = GeometrySmith(s.Normal, v, l, s.Roughness);  

		//Cook-Torrance BRDF:
		vec3 numerator    = NDF * G * F;
		float denominator = 4.0 * max(dot(s.Normal, v), 0.0) * max(dot(s.Normal, l), 0.0) + 0.0001;
		vec3 specular     = numerator / denominator;  

		//refracted diffuse vs reflected specular
		vec3 kS = F;
		vec3 kD = vec3(1.0) - kS;

		
		float NdotL = max(dot(s.Normal, l), 0.0);
		
	`+(bHalfLambert ? 
		`#if 1//Half-Lambert
    	float NdotH = dot(s.Normal, h);
    	float halfLambert = 0.5 * (NdotL + NdotH);
		NdotL = halfLambert;
		const float softness = 0.5;
		NdotL = MapToRange(NdotL, 0.0, 1.0, softness, 1.0);
		#endif` 
	: ``)+/* glsl */`
	#if 0 //Softer NdotL
		const float softness = 0.5;
		NdotL = MapToRange(NdotL, 0.0, 1.0, softness, 1.0);
	#endif

		float specularScale = SceneBuffer.UniformSpecularScale;
		return (kD * s.Albedo / PI + specular * specularScale) * s.Radiance * NdotL;
	}


	//Interection
	
	#define EPSILON 0.000001

	struct Ray {
		vec3 origin;
		vec3 direction;
	};
	
	struct Triangle {
		vec3 v0;
		vec3 v1;
		vec3 v2;
		vec2 uv0;
		vec2 uv1;
		vec2 uv2;
	};
	
	struct TriangleIntersectionResult {
		bool hit;
		vec3 intersectionPoint;
		vec3 barycentricCoords;
	};
	
	TriangleIntersectionResult IntersectRayTriangle(Ray ray, Triangle triangle) {
		TriangleIntersectionResult result;
		result.hit = false;
	
		vec3 edge1 = triangle.v1 - triangle.v0;
		vec3 edge2 = triangle.v2 - triangle.v0;
		vec3 h = cross(ray.direction, edge2);
		float a = dot(edge1, h);
	
		if (a > -EPSILON && a < EPSILON) {
			return result;  // Ray is parallel to the triangle plane
		}
	
		float f = 1.0 / a;
		vec3 s = ray.origin - triangle.v0;
		result.barycentricCoords.x = f * dot(s, h);
	
		if (result.barycentricCoords.x < 0.0 || result.barycentricCoords.x > 1.0) {
			return result;
		}
	
		vec3 q = cross(s, edge1);
		result.barycentricCoords.y = f * dot(ray.direction, q);
	
		if (result.barycentricCoords.y < 0.0 || result.barycentricCoords.x + result.barycentricCoords.y > 1.0) {
			return result;
		}
	
		float t = f * dot(edge2, q);
	
		if (t > EPSILON) {
			result.hit = true;
			result.intersectionPoint = ray.origin + t * ray.direction;
		}
	
		return result;
	}

	vec2 BarycentricToUV(vec3 barycentricCoords, vec2 uv0, vec2 uv1, vec2 uv2) {
		float u = barycentricCoords.x;
		float v = barycentricCoords.y;
		float w = 1.0 - u - v;
	
		vec2 uv = u * uv0 + v * uv1 + w * uv2;
	
		return uv;
	}

#if 0//SHADOW
	float GetShadow()
	{
		//construct vertices
		vec3 planeCenterWS = SceneBuffer.OccluderPlaneCenter;
		vec3 planeExtentWS = SceneBuffer.OccluderPlaneExtent;
		vec3 vLeftDown = vec3(planeCenterWS.x - planeExtentWS.x, planeCenterWS.y - planeExtentWS.y, planeCenterWS.z);
		vec3 vLeftUp = vec3(planeCenterWS.x - planeExtentWS.x, planeCenterWS.y + planeExtentWS.y, planeCenterWS.z);
		vec3 vRightUp = vec3(planeCenterWS.x + planeExtentWS.x, planeCenterWS.y + planeExtentWS.y, planeCenterWS.z);
		vec3 vRightDown = vec3(planeCenterWS.x + planeExtentWS.x, planeCenterWS.y - planeExtentWS.y, planeCenterWS.z);
		/* vec3 vLeftDown = BurningPlaneVertex0;
		vec3 vLeftUp = BurningPlaneVertex1;
		vec3 vRightUp = BurningPlaneVertex2;
		vec3 vRightDown = BurningPlaneVertex3; */
		vec2 uvLeftDown = vec2(0.0, 0.0);
		vec2 uvLeftUp = vec2(0.0, 1.0);
		vec2 uvRightUp = vec2(1.0, 1.0);
		vec2 uvRightDown = vec2(1.0, 0.0);
		
		
		
		Ray rayToLight;
		rayToLight.origin = interpolatorWorldSpacePos;
		rayToLight.direction = l;
		
		//triangle right
		bool bHit = false;
		Triangle planeTriangle;
		planeTriangle.v0 = vRightUp;
		planeTriangle.v1 = vRightDown;
		planeTriangle.v2 = vLeftDown;
		
		planeTriangle.uv0 = uvLeftUp;
		planeTriangle.uv1 = uvRightUp;
		planeTriangle.uv2 = uvLeftDown;
		TriangleIntersectionResult triIntersectionRes = IntersectRayTriangle(rayToLight, planeTriangle);
		if(triIntersectionRes.hit)
		{
			bHit = true;
		}
		else
		{
			planeTriangle.v0 = vLeftDown;
			planeTriangle.v1 = vLeftUp;
			planeTriangle.v2 = vRightUp;
		
			planeTriangle.uv0 = uvRightDown;
			planeTriangle.uv1 = uvLeftDown;
			planeTriangle.uv2 = uvRightUp;
			triIntersectionRes = IntersectRayTriangle(rayToLight, planeTriangle);
			if(triIntersectionRes.hit)
			{
				bHit = true;
			}
		}

		if(bHit && s.PixelPos.z > -0.05)
		{
			vec2 hitTriangleUV = BarycentricToUV(triIntersectionRes.barycentricCoords, planeTriangle.uv0, planeTriangle.uv1, planeTriangle.uv2);
			/* float s = min(1.f, length(hitTriangleUV * 2.f - 1.f));
			s = pow(s, 4.f);
			s = max(0.f, s - hitTriangleUV.y);
			//shadow = s;
			shadow = mix(s, 0.0, hitTriangleUV.y); */

			float scale = 0.f;
			//penumbra size depends on the distance to occluder + dist to light
			float distanceToCurLight = length(SceneBuffer.ProjectorPos - s.PixelPos);
			float dToSource = distanceToCurLight / 10.f;
			dToSource = pow(clamp(dToSource, 0.0, 1.0), 4.f);
			float dToOccluder = length(triIntersectionRes.intersectionPoint - interpolatorWorldSpacePos);
			dToOccluder = clamp(dToOccluder * 0.5, 0.0, 1.0);
			//OutColor = vec4(dToOccluder, 0.0, 0.0, 1.0); return;
			const float penumbraScale = 0.25;
			const float dSourceFade = 0.25f; 
			const float dOccluderFade = 0.5f; 
			float horizonPenumbra = 1.f - (dToSource * dSourceFade + dToOccluder * dToOccluder * dOccluderFade) * penumbraScale ; //higher value sharper penumbra
			//horizonPenumbra = 1.0;
			vec2 rectSpacePos = hitTriangleUV * 2.f - 1.f;
			if(abs(rectSpacePos.x) > horizonPenumbra)
			{
				float m = MapToRange(abs(rectSpacePos.x), horizonPenumbra, 1.0, 0.0, 1.0);
				m = clamp(m, 0.0, 1.0);
				scale += m;
			}
		
			const float verPenumbraScale = 0.5f;
			float vertPenumbra = 1.f - (dToSource + dToOccluder) * verPenumbraScale; //higher value sharper penumbra
			//vertPenumbra = 1.0;
			vertPenumbra = horizonPenumbra;
			if(abs(rectSpacePos.y) > vertPenumbra)
			{
				float m = MapToRange(abs(rectSpacePos.y), vertPenumbra, 1.0, 0.0, 1.0);
				m = clamp(m, 0.0, 1.0);
				scale += m;
			}
		
			//SSAO
			//s *= min(1.f, dToOccluder * 10.f);

			shadow = clamp(scale, 0.0, 1.0);
			//shadow = clamp(shadow + (1.f - hitTriangleUV.y) * 0.1, 0.1, 1.0);

		}
	}
#endif//SHADOW

#if 0
	float GetBayerMask()
	{
		int bayerIndex = 0 % BAYER_LIMIT;
		const ivec2 projTextureSize = ivec2(1024, 512) / 2;
		ivec2 projCoord = ivec2(floor(float(projTextureSize.x) * projUVPos.x), floor(float(projTextureSize.y) * projUVPos.y));
		int bayerCoord = (projCoord.x + BAYER_LIMIT_H * projCoord.y) % BAYER_LIMIT;
		//if(bayerCoord == bayerFilter[bayerIndex])
		if(bayerCoord == 0)
		{
			projectorColor *= 2.0;
		}
		else
		{
			projectorColor *= 0.0f;
		}
	}
#endif


	
	`);
}


export function scGetGIFuncsDeclarations()
{
	return (/* glsl */`
	ivec3 GetSHGridIndex(vec3 posGridSpace, vec3 cellSize) {
		ivec3 index;
		index.x = int(floor(posGridSpace.x / cellSize.x));
		index.y = int(floor(posGridSpace.y / cellSize.y));
		index.z = int(floor(posGridSpace.z / cellSize.z));
		//index = clamp(index, ivec3(0), ivec3(2));
		return index;
	}

	int GetFlatProbeIndex(ivec3 cellId)
	{
		ivec3 numProbes = ivec3(LightprobesGridSize);
		return cellId.z * numProbes.x * numProbes.y + cellId.y * numProbes.x + cellId.x;
	}

	vec3 SampleLightprobe(float[4] shBasisArr, int probeIndex)
	{
		vec3 res;
		const int stride = 4;
		for(int i = 0; i < 4; i++)
		{
			res.r += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			res.g += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(stride + i, probeIndex), 0).r;
			res.b += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(stride*2 + i, probeIndex), 0).r;
		}
		return res;
	}

	vec3 SampleDebugLightprobe(vec3 curDir)
	{
		//Debug Visualizer

		//float shBasisArr[9] = SHBasis2(curDir);
		float shBasisArr[4] = SHBasis1(curDir);

		int probeIndex = LightprobeIndex;
		/* float SHCoefArr[9] = float[9](0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		for(int i = 0; i < 9; i++)
		{
			SHCoefArr[i] = texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
		}

		float res = 0.0;
		for(int i = 0; i < 9; i++)
		{
			res += shBasisArr[i] * SHCoefArr[i];
		} */

		return SampleLightprobe(shBasisArr, probeIndex);
	}

	vec3 GetIrradianceFromLightprobes(vec3 curDir, vec3 curWS)
	{
		//Compute SH Basis for cur normal
		//float shBasisArr[9] = SHBasis2(curDir);
		float shBasisArr[4] = SHBasis1(curDir);

		vec3 sceneExtent = SceneBuffer.SceneBoundsExtent;
		vec3 sceneCenterWS = SceneBuffer.SceneBoundsCenter;
		vec3 sceneLength = sceneExtent * 2.0;
		//const vec3 giGridCellSize = vec3(4.0, 4.0, 4.0);
		vec3 giGridCellSize = sceneLength / float(LightprobesGridSize - 1);

		//Map cur pos to positive scene box range 
		vec3 curPixelPosBoxNorm = curWS;
		curPixelPosBoxNorm = clamp(MapToRange(curPixelPosBoxNorm, sceneCenterWS - sceneExtent, sceneCenterWS + sceneExtent, vec3(0.0), sceneLength), vec3(0.0), sceneLength - vec3(0.05));
		ivec3 curCellCoord = GetSHGridIndex(curPixelPosBoxNorm, giGridCellSize);
		//Transform worldPos to uvw pos in cell space
		vec3 curPixelPosUVW = clamp(mod(curPixelPosBoxNorm, giGridCellSize) / giGridCellSize, vec3(0.0), vec3(1.0));
		curPixelPosUVW.x = smoothstep(0.0, 1.0, curPixelPosUVW.x);
		curPixelPosUVW.y = smoothstep(0.0, 1.0, curPixelPosUVW.y);
		curPixelPosUVW.z = smoothstep(0.0, 1.0, curPixelPosUVW.z);
		
		//Trilinear lerp between probes
		vec3 LineDownUpZ0 = vec3(0.0);
		vec3 LineDownUpZ1 = vec3(0.0);


		//Bilinear Lerp z=0
		{
			int z = 0;
			ivec3 probeCoord = ivec3(0);
			int probeIndex = 0;

			vec3 ProbeLeftDown = vec3(0.0);
			vec3 ProbeLeftUp = vec3(0.0);
			vec3 ProbeRightUp = vec3(0.0);
			vec3 ProbeRightDown = vec3(0.0);

			//Left Down
			probeCoord = curCellCoord + ivec3(0, 0, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeLeftDown = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeLeftDown += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Left Up
			probeCoord = curCellCoord + ivec3(0, 1, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeLeftUp = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeLeftUp += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Right Up
			probeCoord = curCellCoord + ivec3(1, 1, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeRightUp = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeRightUp += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Right Down
			probeCoord = curCellCoord + ivec3(1, 0, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeRightDown = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeRightDown += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */

			//Line Up
			vec3 LineDown = mix(ProbeLeftDown, ProbeRightDown, curPixelPosUVW.x);
			vec3 LineUp = mix(ProbeLeftUp, ProbeRightUp, curPixelPosUVW.x);
			LineDownUpZ0 = mix(LineDown, LineUp, curPixelPosUVW.y);
		}

		//Bilinear Lerp z=1
		{
			int z = 1;
			ivec3 probeCoord = ivec3(0);
			int probeIndex = 0;

			vec3 ProbeLeftDown = vec3(0.0);
			vec3 ProbeLeftUp = vec3(0.0);
			vec3 ProbeRightUp = vec3(0.0);
			vec3 ProbeRightDown = vec3(0.0);

			//Left Down
			probeCoord = curCellCoord + ivec3(0, 0, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeLeftDown = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeLeftDown += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Left Up
			probeCoord = curCellCoord + ivec3(0, 1, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeLeftUp = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeLeftUp += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Right Up
			probeCoord = curCellCoord + ivec3(1, 1, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeRightUp = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeRightUp += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */
			//Right Down
			probeCoord = curCellCoord + ivec3(1, 0, z);
			probeIndex = GetFlatProbeIndex(probeCoord);
			ProbeRightDown = SampleLightprobe(shBasisArr, probeIndex);
			/* for(int i = 0; i < 9; i++)
			{
				ProbeRightDown += shBasisArr[i] * texelFetch(SHCoefficientsTexture, ivec2(i, probeIndex), 0).r;
			} */

			//Line Up
			vec3 LineDown = mix(ProbeLeftDown, ProbeRightDown, curPixelPosUVW.x);
			vec3 LineUp = mix(ProbeLeftUp, ProbeRightUp, curPixelPosUVW.x);
			LineDownUpZ1 = mix(LineDown, LineUp, curPixelPosUVW.y);
		}

		//Z Lerp
		vec3 finalLerp = mix(LineDownUpZ0, LineDownUpZ1, curPixelPosUVW.z);

		return (finalLerp);
	}
	`);
}


export function scGetProceduralCubeVertsDecl()
{
	return /* glsl */`
	const vec3 CubePositions[36] = vec3[](
        // Front face
        vec3(-1.0, -1.0,  1.0), vec3(-1.0,  1.0,  1.0), vec3( 1.0,  1.0,  1.0),
        vec3(-1.0, -1.0,  1.0), vec3( 1.0,  1.0,  1.0), vec3( 1.0, -1.0,  1.0),
        // Back face
        vec3(-1.0, -1.0, -1.0), vec3( 1.0, -1.0, -1.0), vec3( 1.0,  1.0, -1.0),
        vec3(-1.0, -1.0, -1.0), vec3( 1.0,  1.0, -1.0), vec3(-1.0,  1.0, -1.0),
        // Left face
        vec3(-1.0, -1.0, -1.0), vec3(-1.0,  1.0, -1.0), vec3(-1.0,  1.0,  1.0),
        vec3(-1.0, -1.0, -1.0), vec3(-1.0,  1.0,  1.0), vec3(-1.0, -1.0,  1.0),
        // Right face
        vec3( 1.0, -1.0, -1.0), vec3( 1.0, -1.0,  1.0), vec3( 1.0,  1.0,  1.0),
        vec3( 1.0, -1.0, -1.0), vec3( 1.0,  1.0,  1.0), vec3( 1.0,  1.0, -1.0),
        // Top face
        vec3(-1.0,  1.0, -1.0), vec3( 1.0,  1.0, -1.0), vec3( 1.0,  1.0,  1.0),
        vec3(-1.0,  1.0, -1.0), vec3( 1.0,  1.0,  1.0), vec3(-1.0,  1.0,  1.0),
        // Bottom face
        vec3(-1.0, -1.0, -1.0), vec3(-1.0, -1.0,  1.0), vec3( 1.0, -1.0, -1.0),
        vec3( 1.0, -1.0, -1.0), vec3(-1.0, -1.0,  1.0), vec3( 1.0, -1.0,  1.0)
    );

	/* vec3(-1.0, -1.0,  1.0), // 0 */
    /* vec3(-1.0,  1.0,  1.0), // 1 */
    /* vec3( 1.0,  1.0,  1.0), // 2 */
    /* vec3( 1.0, -1.0,  1.0), // 3 */
    /* vec3(-1.0, -1.0, -1.0), // 4 */
    /* vec3(-1.0,  1.0, -1.0), // 5 */
    /* vec3( 1.0,  1.0, -1.0), // 6 */
    /* vec3( 1.0, -1.0, -1.0)  // 7 */

	// Mapping of 36 vertices to 8 unique vertices
	const uint VertexToIndexMap[36] = uint[](
    	0u, 1u, 2u,
		0u, 2u, 3u,
    	
    	4u, 7u, 6u,
		4u, 6u, 5u,
    	
    	4u, 5u, 1u,
		4u, 1u, 0u,
    	
    	7u, 3u, 2u,
		7u, 2u, 6u,
    	
    	5u, 6u, 2u,
		5u, 2u, 1u,
    	
    	4u, 0u, 7u,
		7u, 0u, 3u
	);
	
		const vec2 CubeTexCoords[36] = vec2[](
			// Front face
			vec2(0.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 1.0),
			vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(1.0, 0.0),
			// Back face
			vec2(0.0, 0.0), vec2(1.0, 0.0), vec2(1.0, 1.0),
			vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(0.0, 1.0),
			// Left face
			vec2(0.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 1.0),
			vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(1.0, 0.0),
			// Right face
			vec2(0.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 1.0),
			vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(1.0, 0.0),
			// Top face
			vec2(0.0, 0.0), vec2(1.0, 0.0), vec2(1.0, 1.0),
			vec2(0.0, 0.0), vec2(1.0, 1.0), vec2(0.0, 1.0),
			// Bottom face
			vec2(0.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 0.0),
			vec2(1.0, 0.0), vec2(0.0, 1.0), vec2(1.0, 1.0)
		);
	
		const vec3 CubeNormals[36] = vec3[](
			// Front face
			vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0),
			vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0),
			// Back face
			vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0),
			vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0),
			// Left face
			vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
			vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
			// Right face
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			// Top face
			vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0),
			vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 1.0, 0.0),
			// Bottom face
			vec3(0.0, -1.0, 0.0), vec3(0.0, -1.0, 0.0), vec3(0.0, -1.0, 0.0),
			vec3(0.0, -1.0, 0.0), vec3(0.0, -1.0, 0.0), vec3(0.0, -1.0, 0.0)
		);

		const vec3 CubeTangents[36] = vec3[](
			// Front face
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			// Back face
			vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
			vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
			// Left face
			vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0),
			vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0), vec3(0.0, 0.0, 1.0),
			// Right face
			vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0),
			vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0), vec3(0.0, 0.0, -1.0),
			// Top face
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			// Bottom face
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0),
			vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)
		);
		

		#define CUBE_FACE_FRONT 0u
		#define CUBE_FACE_BACK 1u
		#define CUBE_FACE_LEFT 2u
		#define CUBE_FACE_RIGHT 3u
		#define CUBE_FACE_TOP 4u
		#define CUBE_FACE_BOTTOM 5u

		const uint CubeFaces[6] = uint[](
			CUBE_FACE_FRONT,
			CUBE_FACE_BACK,
			CUBE_FACE_LEFT,
			CUBE_FACE_RIGHT,
			CUBE_FACE_TOP,
			CUBE_FACE_BOTTOM
		);

		uint GetCurCubeFaceIndex(uint vertId)
		{
			return CubeFaces[vertId / 6u];
		}

	`;
}
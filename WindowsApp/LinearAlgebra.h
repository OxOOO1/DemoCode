#pragma once

#include "MathGeneral.h"

inline float3 VectorMax(float val, c_float3& vec)
{
	return float3{ std::max(val, vec.x), std::max(val, vec.y), std::max(val, vec.z) };
}

inline float VectorDot(const float2& one, const float2& other)
{
	return
		one.x * other.x
		+ one.y * other.y
		;
}
inline float VectorDot(const float3& one, const float3& other)
{
	return
		one.x * other.x
		+ one.y * other.y
		+ one.z * other.z
		;
}
inline float VectorDot(const float4& one, const float4& other)
{
	return
		one.x * other.x
		+ one.y * other.y
		+ one.z * other.z
		+ one.w * other.w
		;
}
inline float3 VectorMul(const float3& one, float other) { return float3{ one.x * other, one.y * other, one.z * other }; }
inline float2 VectorMul(const float2& one, float other) { return float2{ one.x * other, one.y * other }; }
inline float3 VectorMul(const float3& vec, const Matrix4x4& mat)
{
	float3 result;
	result.x = vec.x * mat.m[0][0] + vec.y * mat.m[1][0] + vec.z * mat.m[2][0];
	result.y = vec.x * mat.m[0][1] + vec.y * mat.m[1][1] + vec.z * mat.m[2][1];
	result.z = vec.x * mat.m[0][2] + vec.y * mat.m[1][2] + vec.z * mat.m[2][2];
	return result;

}
inline float4 VectorMul(const float4& vec, const Matrix4x4& mat)
{
	float4 result;
	result.x = vec.x * mat.m[0][0] + vec.y * mat.m[1][0] + vec.z * mat.m[2][0] + vec.w * mat.m[3][0];
	result.y = vec.x * mat.m[0][1] + vec.y * mat.m[1][1] + vec.z * mat.m[2][1] + vec.w * mat.m[3][1];
	result.z = vec.x * mat.m[0][2] + vec.y * mat.m[1][2] + vec.z * mat.m[2][2] + vec.w * mat.m[3][2];
	result.w = vec.x * mat.m[0][3] + vec.y * mat.m[1][3] + vec.z * mat.m[2][3] + vec.w * mat.m[3][3];
	return result;

}
inline float3 VectorAdd(const float3& one, const float3& other) { return float3{ one.x + other.x, one.y + other.y, one.z + other.z }; }
inline float3 VectorDiv(const float3& one, float other) { return float3{ one.x / other, one.y / other, one.z / other }; }
inline float2 VectorDiv(const float2& one, float other) { return float2{ one.x / other, one.y / other }; }
inline float3 VectorNegate(const float3& vec)
{
	return float3{ -vec.x, -vec.y, -vec.z };
}
inline float3 VectorAbs(const float3& vec)
{
	return float3{ fabs(vec.x), fabs(vec.y), fabs(vec.z) };
}


inline float3 VectorCross(const float3& first, const float3& second)
{
	float3 result;
	result.x = first.y * second.z - second.y * first.z; //YZ for X
	result.y = first.z * second.x - second.z * first.x; //ZX for Y
	result.z = first.x * second.y - second.x * first.y; //XY for Z
	return result;
}


inline float VectorGetLength(c_float3& vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}
inline float VectorGetLength(const float2& vec)
{
	return sqrtf(vec.x * vec.x + vec.y * vec.y);
}
inline float VectorGetLength(float vec)
{
	return abs(vec);///sqrt(pow(vec,2))
}
static inline float VectorGetLengthSquared(c_float3& vec)
{
	return (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}
static inline float VectorGetLengthSquared(c_float2& vec)
{
	return (vec.x * vec.x + vec.y * vec.y);
}
inline float3 VectorNormalize(const float3& vec)
{
	//divide by length
	float length = VectorGetLength(vec);
	return VectorDiv(vec, length);
}
inline float2 VectorNormalize(const float2& vec)
{
	//divide by length
	float length = VectorGetLength(vec);
	return VectorDiv(vec, length);
}
inline float VectorNormalize(const float vec)
{
	return vec / fabs(vec);
}

inline Matrix4x4 MatrixMultiply(Matrix4x4 mat, float scalar)
{
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			mat.m[r][c] = mat.m[r][c] * scalar;
		}
	}
	return mat;
}

inline Matrix4x4 MatrixMultiply(const Matrix4x4& matA, const Matrix4x4& matB)
{
	Matrix4x4 matC;

	//matC._11 = matA._11 * matB._11 + matA._12 * matB._21 + matA._13 * matB._31 + matA._14 * matB._41;

	//For each element cij in the result, locate row i in A and column j in B
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matC.m[i][j] =
				matA.m[i][0] * matB.m[0][j]
				+ matA.m[i][1] * matB.m[1][j]
				+ matA.m[i][2] * matB.m[2][j]
				+ matA.m[i][3] * matB.m[3][j];
		}
	}

	return matC;

}

inline Matrix4x4 MatrixTranslate(const float3& translation)
{
	Matrix4x4 mat;


	mat._41 = translation.x; mat._42 = translation.y; mat._43 = translation.z; mat._44 = 1;

	return mat;

}
inline Matrix4x4 MatrixRotateAroundAxis(const float3& axis, float angle)
{
	float3 n = VectorNormalize(axis);
	float a = angle;

	Matrix4x4 mat;
	mat._44 = 1;

	mat._11 = (n.x * n.x) * (1 - cos(a)) + cos(a);			 mat._12 = (n.x * n.y) * (1 - cos(a)) + n.z * sin(a);  mat._13 = (n.x * n.z) * (1 - cos(a)) - n.y * sin(a);

	mat._21 = (n.x * n.y) * (1 - cos(a)) - n.z * sin(a);	 mat._22 = (n.y * n.y) * (1 - cos(a)) + cos(a);		   mat._23 = (n.y * n.z) * (1 - cos(a)) + n.x * sin(a);

	mat._31 = (n.x * n.z) * (1 - cos(a)) + n.y * sin(a);	 mat._32 = (n.y * n.z) * (1 - cos(a)) - n.x * sin(a);  mat._33 = (n.z * n.z) * (1 - cos(a)) + cos(a);


	return mat;
}


inline Matrix4x4 MatrixOrientationEuler(float pitch, float yaw, float roll)
{
	Matrix4x4 orientation;

	/*Step 1. Begin in the “identity” orientation—that is, with the object - space axes aligned with the upright axes.*/
	float3 dir = float3{ 0,0,1 };
	float3 up = float3{ 0, 1, 0 };
	float3 right = float3{ 1, 0, 0 };

	/*Step 2. Perform the heading rotation, rotating about the y - axis. Positive rotation rotates to the right(clockwise when viewed from above).*/

	orientation = MatrixRotateAroundAxis(up, yaw);

	/*Step 3. After heading has been applied, pitch measures the amount of rotation about the x - axis.
	This is the object - space x - axis, not the upright x - axis.
	Staying consistent with the left - hand rule, positive rotation rotates downward.
	In other words, pitch actually measures the angle of declination.*/

	right = VectorMul(right, orientation);

	orientation = MatrixMultiply(orientation, MatrixRotateAroundAxis(right, pitch));

	/*Step 4. After heading and pitch angles have been applied, bank measures the amount of rotation about the z - axis.
	Again, this is the objectspace z - axis, not the original upright - space z - axis.
	The left - hand rule dictates that positive bank rotates counterclockwise when viewed from the origin looking towards + z.*/


	dir = VectorMul(dir, orientation);

	orientation = MatrixMultiply(orientation, MatrixRotateAroundAxis(dir, roll));

	return orientation;
}



inline Matrix4x4 MatrixProjectPerspective(
	float FovAngleY,
	float AspectRatio,
	float NearZ,
	float FarZ
)
{
	/*float swapCopy = FarZ;
	FarZ = NearZ;
	NearZ = swapCopy;*/

	float    SinFov;
	float    CosFov;
	XMScalarSinCos(&SinFov, &CosFov, 0.5f * FovAngleY);

	float Height = CosFov / SinFov;
	float Width = Height / AspectRatio;
	float fRange = (FarZ) / (FarZ - NearZ);

	//z = z * fRange + (Z * -fRange * NearZ)

	Matrix4x4 M;
	M.m[0][0] = Width;
	M.m[0][1] = 0.0f;
	M.m[0][2] = 0.0f;
	M.m[0][3] = 0.0f;

	M.m[1][0] = 0.0f;
	M.m[1][1] = Height;
	M.m[1][2] = 0.0f;
	M.m[1][3] = 0.0f;

	M.m[2][0] = 0.0f;
	M.m[2][1] = 0.0f;
	M.m[2][2] = fRange;
	M.m[2][3] = 1.0f;//copy Z to W

	M.m[3][0] = 0.0f;
	M.m[3][1] = 0.0f;
	M.m[3][2] = -fRange * NearZ;
	M.m[3][3] = 0.0f;
	return M;

	//for New Z = Z * (FarZ) / (FarZ - NearZ) - (FarZ * NearZ) / (FarZ - NearZ)
}



inline float3 VectorReflect(const float3& d, const float3& n)
{
	auto normalizedN = VectorNormalize(n);
	float dot = VectorDot(d, normalizedN);
	float3 dParallel = VectorMul(normalizedN, dot);
	return d - dParallel - dParallel;

}

//Around Z means only XY values change
inline Matrix4x4 MatrixRotateAroundZ(float angle)
{
	Matrix4x4 mat;
	mat._44 = 1;

	//Z basis (row 3) remains unchanged
	mat._33 = 1;//Identity

	//row 1 is for X basis
	mat._11 = cos(angle); mat._12 = sin(angle);
	//row 2 for Y basis
	mat._21 = -sin(angle); mat._22 = cos(angle);


	return mat;
}

inline float2 VectorRotate(float2 vec, float angle)
{
	float2 res;
	res.x = vec.x * cos(angle) + vec.y * -sin(angle);
	res.y = vec.x * sin(angle) + vec.y * cos(angle);
	return res;
}

inline Matrix4x4 MatrixRotateAroundX(float angle)
{
	Matrix4x4 mat;
	mat._44 = 1;

	//X basis (row 1) remains unchanged
	mat._11 = 1;//Identity

	//row 2 is for Y basis
	mat._22 = cos(angle); mat._23 = sin(angle);
	//row 3 for Z basis
	mat._32 = -sin(angle); mat._33 = cos(angle);


	return mat;

}

inline Matrix4x4 MatrixRotateAroundY(float angle)
{
	Matrix4x4 mat;
	mat._44 = 1;

	//Y basis (row 2) remains unchanged
	mat._22 = 1;//Identity

	//row 1 is for X basis
	mat._11 = cos(angle); mat._13 = -sin(angle);
	//row 3 for Z basis
	mat._31 = sin(angle); mat._33 = cos(angle);


	return mat;

}

inline void GetLineFromDirectionVector(float3& outLineStart, float3& outLineEnd, float3 dir, float3 origin = float3{ 0,0,0 }, float scale = 1.f)
{
	outLineStart = origin - VectorMul(dir, 100 * scale);
	outLineEnd = origin + VectorMul(dir, 100 * scale);
}




inline Matrix4x4 MatrixScaleUniform(float scale)
{
	Matrix4x4 mat;
	mat._44 = 1;

	//	scale	0	0
	//{	0	scale	0 }
	//	0	0	scale

	mat._11 = scale;
	mat._22 = scale;
	mat._33 = scale;

	return mat;
}

inline Matrix4x4 MatrixScaleInDirection(const float3& direction, float k)
{
	float3 n = VectorNormalize(direction);

	Matrix4x4 mat;
	mat._44 = 1;

	mat._11 = 1 + (k - 1) * (n.x * n.x); mat._12 = (k - 1) * (n.x * n.y);	  mat._13 = (k - 1) * (n.x * n.z);
	mat._21 = (k - 1) * (n.x * n.y);	 mat._22 = 1 + (k - 1) * (n.y * n.y); mat._23 = (k - 1) * (n.y * n.z);
	mat._31 = (k - 1) * (n.x * n.z);	 mat._32 = (k - 1) * (n.y * n.z);	  mat._33 = 1 + (k - 1) * (n.z * n.z);


	return mat;
}



inline Matrix4x4 MatrixProjectOntoPlane(const float3& planeNormal)
{
	float3 n = VectorNormalize(planeNormal);

	Matrix4x4 mat;
	mat._44 = 1;

	mat._11 = 1 + (0 - 1) * (n.x * n.x); mat._12 = (0 - 1) * (n.x * n.y);	  mat._13 = (0 - 1) * (n.x * n.z);
	mat._21 = (0 - 1) * (n.x * n.y);	 mat._22 = 1 + (0 - 1) * (n.y * n.y); mat._23 = (0 - 1) * (n.y * n.z);
	mat._31 = (0 - 1) * (n.x * n.z);	 mat._32 = (0 - 1) * (n.y * n.z);	  mat._33 = 1 + (0 - 1) * (n.z * n.z);


	return mat;
}

///Map point from [Near, Far] Range to [0,1]
inline Matrix4x4 MatrixProjectionOrtho(float ViewWidth, float ViewHeight, float FarZ, float NearZ)
{
	float fRange = 1.0f / (FarZ - NearZ);

	Matrix4x4 M;
	M.m[0][0] = 2.0f / ViewWidth;
	M.m[0][1] = 0.0f;
	M.m[0][2] = 0.0f;
	M.m[0][3] = 0.0f;

	M.m[1][0] = 0.0f;
	M.m[1][1] = 2.0f / ViewHeight;
	M.m[1][2] = 0.0f;
	M.m[1][3] = 0.0f;

	M.m[2][0] = 0.0f;
	M.m[2][1] = 0.0f;
	M.m[2][2] = fRange;
	M.m[2][3] = 0.0f;

	M.m[3][0] = 0.0f;
	M.m[3][1] = 0.0f;
	M.m[3][2] = -fRange * NearZ;
	M.m[3][3] = 1.0f;
	return M;
}



inline Matrix4x4 MatrixTranspose(const Matrix4x4 m)
{
	Matrix4x4 ret;
	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			ret.m[c][r] = m.m[r][c];
		}
	}
	return ret;
}


inline Matrix4x4 MatrixLookAt(float3 position, float3 dir, float3 up, float3 right)
{
	float3 NegEyePosition = VectorNegate(position);

	right = VectorCross(up, dir);
	right = VectorNormalize(right);
	//right = VectorNegate(right);

	Matrix4x4 M;
	M.m[0][0] = right.x; M.m[0][1] = right.y; M.m[0][2] = right.z;
	M.m[1][0] = up.x;	 M.m[1][1] = up.y;	  M.m[1][2] = up.z;
	M.m[2][0] = dir.x;	 M.m[2][1] = dir.y;	  M.m[2][2] = dir.z;
	M.m[3][0] = 0;		 M.m[3][1] = 0;		  M.m[3][2] = 0;

	M.m[3][3] = 1;

	Matrix4x4 MatTranslation = MatrixTranslate(NegEyePosition);
	return MatrixMultiply(MatTranslation, MatrixTranspose(M));

}

inline Matrix4x4 MatrixBasisTransform(float3 position, float3 basisZ, float3 basisY, float3 basisX)
{
	basisX = VectorCross(basisY, basisZ);
	basisX = VectorNormalize(basisX);

	Matrix4x4 M;
	M.m[0][0] = basisX.x;	 M.m[0][1] = basisX.y;	  M.m[0][2] = basisX.z;
	M.m[1][0] = basisY.x;	 M.m[1][1] = basisY.y;	  M.m[1][2] = basisY.z;
	M.m[2][0] = basisZ.x;	 M.m[2][1] = basisZ.y;	  M.m[2][2] = basisZ.z;
	M.m[3][0] = 0;		     M.m[3][1] = 0;			  M.m[3][2] = 0;	      M.m[3][3] = 1;

	Matrix4x4 MatTranslation = MatrixTranslate(position);
	return MatrixMultiply(MatTranslation, M);
}

inline Matrix4x4 MatrixBasisTransformInversed(float3 position, float3 basisZ, float3 basisY, float3 basisX)
{
	float3 NegEyePosition = VectorNegate(position);

	basisX = VectorCross(basisY, basisZ);
	basisX = VectorNormalize(basisX);

	Matrix4x4 M;
	M.m[0][0] = basisX.x; M.m[0][1] = basisX.y; M.m[0][2] = basisX.z;
	M.m[1][0] = basisY.x;	 M.m[1][1] = basisY.y;	  M.m[1][2] = basisY.z;
	M.m[2][0] = basisZ.x;	 M.m[2][1] = basisZ.y;	  M.m[2][2] = basisZ.z;
	M.m[3][0] = 0;		 M.m[3][1] = 0;		  M.m[3][2] = 0;

	M.m[3][3] = 1;

	Matrix4x4 MatTranslation = MatrixTranslate(position);
	M = MatrixMultiply(MatTranslation, M);

	DirectX::XMFLOAT4X4 MDirect;

	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			MDirect.m[r][c] = M.m[r][c];
		}
	}

	DirectX::XMMATRIX mXM = DirectX::XMLoadFloat4x4(&MDirect);
	mXM = DirectX::XMMatrixInverse(nullptr, mXM);
	DirectX::XMStoreFloat4x4(&MDirect, mXM);

	for (int r = 0; r < 4; r++)
	{
		for (int c = 0; c < 4; c++)
		{
			M.m[r][c] = MDirect.m[r][c];
		}
	}

	return M;

}


inline float3 GetHalfVertex(float3 v1, float3 v2)
{
	float3 ret = VectorDiv(VectorAdd(v1, v2), 2);

	return VectorNormalize(ret);
}





static inline float MapToRange(float t, float t0, float t1, float newt0, float newt1)
{

	float scale = (newt1 - newt0) / (t1 - t0);
	float offset = -t0 * (newt1 - newt0) / (t1 - t0) + newt0;

	///Translate to origin, scale by ranges ratio, translate to new position
	return (t - t0) * ((newt1 - newt0) / (t1 - t0)) + newt0;

	//map to [0,1]
	t = (t - t0) / (t1 - t0);
	//Interpolate new range with normalised t
	return (1 - t) * newt0 + t * newt1;
}

static inline float2 MapToRange(c_float2& p, c_float2& p0, c_float2& p1, c_float2& newp0, c_float2& newp1)
{
	///Translate to origin, scale by ranges ratio, translate to new position
	///(p - p0) * (((newp1 - newp0) / (p1 - p0)) + newp0;
	float2 res;
	res.x = (p.x - p0.x) * ((newp1.x - newp0.x) / (p1.x - p0.x)) + newp0.x;
	res.y = (p.y - p0.y) * ((newp1.y - newp0.y) / (p1.y - p0.y)) + newp0.y;
	return res;
	///TODO:4x4 Matrix form
}

static inline float3 MapToRange(c_float3& p, c_float3& p0, c_float3& p1, c_float3& newp0, c_float3& newp1)
{
	///Translate to origin, scale by ranges ratio, translate to new position
	///(p - p0) * (((newp1 - newp0) / (p1 - p0)) + newp0;
	float3 res;
	res.x = (p.x - p0.x) * ((newp1.x - newp0.x) / (p1.x - p0.x)) + newp0.x;
	res.y = (p.y - p0.y) * ((newp1.y - newp0.y) / (p1.y - p0.y)) + newp0.y;
	res.z = (p.z - p0.z) * ((newp1.z - newp0.z) / (p1.z - p0.z)) + newp0.z;
	return res;
}

static inline float3 MapTo01Range(c_float3& p, c_float3& p0, c_float3& p1)
{
	///Translate to origin, scale by ranges ratio, translate to new position
	///(p - p0) * (((newp1 - newp0) / (p1 - p0)) + newp0;
	float3 res;
	res.x = (p.x - p0.x) * ((1) / (p1.x - p0.x));
	res.y = (p.y - p0.y) * ((1) / (p1.y - p0.y));
	res.z = (p.z - p0.z) * ((1) / (p1.z - p0.z));
	return res;
}

static inline float3 MapTo01RangeRescale(c_float3& p, c_float3& curSize)
{
	float3 res;
	res.x = p.x / curSize.x;
	res.y = p.y / curSize.y;
	res.z = p.z / curSize.z;
	return res;
}

static inline float MapTo01Range(float t, float t0, float t1)
{
	return (t - t0) / (t1 - t0);
}

static inline float MapTo01RangeClamped(float t, float t0, float t1)
{
	return std::clamp((t - t0) / (t1 - t0), 0.f, 1.f);
}

static inline float MapTo01RangeClampedSmooth(float t, float t0, float t1)
{
	float res = std::clamp((t - t0) / (t1 - t0), 0.f, 1.f);
	
}

//1/z mapping like the one done for DeviceZ
static inline float MapTo01RangeNonLinear(float t, float t0, float t1)
{
	//Z* (FarZ) / (FarZ - NearZ) - (FarZ * NearZ) / (FarZ - NearZ)
	check(t0 > 0);

	return (t * t1 / (t1 - t0) - (t1 * t0) / (t1 - t0)) / t;
}


///TODO: Map Mesh to Spherical Range->
///Take a VertexBuffer, construct an AABB out of it,
///translate VertexBuffer so that AABB center is at the origin,
///rescale-map the extents to [0,1], perform spherical transform, inverse the tranlation and scaling


static inline float3 MapCubeRangeToSphereRange(c_float3& vec, c_float3& aabbExtent)
{
	float3 res;
	float3 vec01 = MapTo01RangeRescale(vec, aabbExtent);///Rescale to [-1, 1]
	res.x = vec01.x * sqrt(1 - (pow(vec01.y, 2) / 2.f) - pow(vec01.z, 2) / 2.f + pow(vec01.y, 2) * pow(vec01.z, 2) / 3);
	res.y = vec01.y * sqrt(1 - (pow(vec01.z, 2) / 2.f) - pow(vec01.x, 2) / 2.f + pow(vec01.z, 2) * pow(vec01.x, 2) / 3);
	res.z = vec01.z * sqrt(1 - (pow(vec01.x, 2) / 2.f) - pow(vec01.y, 2) / 2.f + pow(vec01.x, 2) * pow(vec01.y, 2) / 3);
	res = float3{ aabbExtent.x * res.x, aabbExtent.y * res.y, aabbExtent.z * res.z };///Rescale back to original size
	return res;
}

static inline float2 MapPlaneRangeToCircle(float2 vec, float2 planeExtent)
{
	vec = float2{ vec.x / planeExtent.x, vec.y / planeExtent.y };///Rescale to [-1, 1]
	vec.x = vec.x * sqrt(1 - pow(vec.y, 2) / 2.f);
	vec.y = vec.y * sqrt(1 - pow(vec.x, 2) / 2.f);
	return float2{ vec.x * planeExtent.x, vec.y * planeExtent.y };///Rescale back to original size
}


static inline float3 ProjectCurvilinear(c_float3& vecViewSpace, c_float3& preProjectScale, float postProjectScale = 1.f)
{
	float3 vecScaled{ vecViewSpace.x * preProjectScale.x, vecViewSpace.y * preProjectScale.y, vecViewSpace.z * preProjectScale.z };
	///Divide by dist to origin, not distance along .z
	float distToOrigin = VectorGetLength(vecScaled);
	float3 vecProjected = float3{ vecViewSpace.x / distToOrigin, vecViewSpace.y / distToOrigin, vecViewSpace.z };
	/// Re-scale
	return vecProjected * postProjectScale;
}


static inline float2 MapPolarToCartesian(float angle, float r)
{
	float2 xy;
	xy.x = cos(angle) * r;
	xy.y = sin(angle) * r;
	return xy;
}

static inline float3 MapSphericalToCartesian(float theta, float phi, float r)
{
	float3 vec;
	vec.x = cos(phi) * cos(theta) * r;
	vec.y = sin(phi) * r;
	vec.z = cos(phi) * sin(theta) * r;
	return vec;
}


static inline void MapCartesianToSpherical(c_float3& dir, float& theta, float& phi)
{
	theta = std::atan2(dir.z, dir.x);
	phi = std::asin(-dir.y);
}

// A function to convert cartesian coordinates to barycentric coordinates
// p: the point in cartesian coordinates
// p0, p1, p2: the vertices of the triangle in cartesian coordinates
// returns a float3 containing the barycentric coordinates (s, t, 1 - s - t)
static inline float3 CartesianToBarycentric(c_float3& p, c_float3& p0, c_float3& p1, c_float3& p2)
{
	// Find the vectors u, v and w
	float3 u = p1 - p0;
	float3 v = p2 - p0;
	float3 w = p - p0;

	// Find the dot products uu, uv, vv, wu and wv
	float uu = VectorDot(u, u);
	float uv = VectorDot(u, v);
	float vv = VectorDot(v, v);
	float wu = VectorDot(w, u);
	float wv = VectorDot(w, v);

	// Calculate the barycentric coordinates s and t
	float s = (uv * wv - vv * wu) / (uv * uv - uu * vv);
	float t = (uv * wu - uu * wv) / (uv * uv - uu * vv);

	// Return the barycentric coordinates as a float3
	return float3(s, t, 1 - s - t);
}


static inline float3 GetArithmeticMean(const RDynamicVector<float3>& points)
{
	float3 sum{ 0,0,0 };
	for (auto& v : points)
	{
		sum = sum + v;
	}
	return VectorDiv(sum, points.size());
}

///Root Mean Square (RMS)
static inline float GetQuadraticMean(const RDynamicVector<float>& values)
{
	float sum{ 0 };
	for (float v : values)
	{
		sum += powf(v, 2);
	}
	float res = sum / (float)values.size();
	return sqrtf(res);
}

static inline float GetGeometricMean(const RDynamicVector<float>& values)
{
	float product{ 0 };
	for (float v : values)
	{
		product *= v;
	}
	return powf(product, 1.f / (float)values.size());
}

static inline float3 GetHarmonicMean(const RDynamicVector<float3>& points)
{
	float3 sum{ 0,0,0 };
	for (auto& v : points)
	{
		sum = sum + float3{ 1.f / v.x, 1.f / v.y, 1.f / v.z };
	}
	float n = points.size();
	return float3{ n / sum.x, n / sum.y, n / sum.z };
}





// A function to compute normal, tangent and bitangent using central difference method
static inline void ComputeNTBForGrid(const RVector2D<float3>& input, RVector2D<float3>& normal, RVector2D<float3>* tangent = nullptr, RVector2D<float3>* bitangent = nullptr)
{
	// Get the dimensions of the input
	uint2 size = input.Size;

	// Resize the output arrays to match the input
	normal.Resize(size);
	if (tangent)
	{
		tangent->Resize(size);
		if (bitangent)
		{
			bitangent->Resize(size);
		}
	}
	
	// Loop through each point in the input
	for (uint x = 0; x < size.x; x++)
	{
		for (uint y = 0; y < size.y; y++)
		{
			// Get the neighboring points in the u and v directions
			float3 u1 = (y > 0) ? input[uint2{ x, y - 1 }] : input[uint2{ x, y }];
			float3 u2 = (y < size.y - 1) ? input[uint2{ x, y + 1 }] : input[uint2{ x, y }];
			float3 v1 = (x > 0) ? input[uint2{ x - 1, y }] : input[uint2{ x, y }];
			float3 v2 = (x < size.x - 1) ? input[uint2{ x + 1, y }] : input[uint2{ x, y }];

			// Compute the central difference vectors in the u and v directions
			float3 du = u2 - u1;
			float3 dv = v2 - v1;

			// Compute the cross product of du and dv to get the normal vector
			float3 n = VectorCross(du, dv);

			// Normalize the normal vector
			n = VectorNormalize(n);

			// Store the normal, tangent and bitangent vectors in the output arrays
			normal[uint2{ x, y }] = n;
			if (tangent)
			{
				// Compute the tangent vector as the normalized du vector
				float3 t = VectorNormalize(du);
				tangent->GetElement(uint2{ x, y }) = t;

				if (bitangent)
				{
					// Compute the bitangent vector as the cross product of the normal and tangent vectors
					float3 b = VectorCross(n, t);
					bitangent->GetElement(uint2{ x, y }) = b;
				}
			}
			
		}
	}
}
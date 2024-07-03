#pragma once

#include "Utilities/Types.h"

#include <DirectXMath.h>
#include <d3d12.h>
#include <cmath>
#include <functional>
#include <type_traits>
#include "thirdParty/d3dx12.h"

#include "Resources/Texture/Texture.h"
#include "CommandContext.h"
#include "Scene/Scene.h"

#include "Resources/UniformBuffer.h"
#include "Resources/Shader.h"

#include "Core/Renderer/Passes/RenderPass.h"

#include "ShaderBindings.h"
#include "Sampler.h"

#include "Utilities/Math/Random.h"
#include "Utilities/Time.h"



extern bool GbFirstBoot;


///////////////////////////////////////
/////			Constants
///////////////////////////////////////
static constexpr float Pi = 3.141592654f;
static constexpr float PiX2 = 6.283185307f;
static constexpr float PiDiv2 = 1.570796327f;
static constexpr float PiDiv4 = 0.785398163f;
static constexpr float e = 2.718281828459;
constexpr float GoldenRatio = 1.61803398875;
constexpr float GoldenRatioX2 = 3.2360679775;


#define DEBUG_MEMORY (1 && _DEBUG)


///////////////////////////////////////
/////			Containers
///////////////////////////////////////

#define USE_EASTL 0

#if USE_EASTL
#include <EASTL/hash_map.h>
#include <EASTL/fixed_vector.h>
#include <EASTL/map.h>
#include <EASTL/vector.h>
#include <EASTL/string.h>
#include <EASTL/fixed_string.h>
#include <EASTL/string_hash_map.h>
#include <EASTL/fixed_list.h>
#include <EASTL/fixed_hash_map.h>
#endif

///Vector
#define USE_EASTL_VECTOR (1 && USE_EASTL)

template<typename V>
#if USE_EASTL_VECTOR
using RDynamicVector = eastl::vector<V>;
#else
using RDynamicVector = std::vector<V>;
#endif//USE_EASTL_VECTOR

template<typename V, size_t Size = 100>
#if USE_EASTL_VECTOR
using RStackVector = eastl::fixed_vector<V, Size, true>;
#else
using RStackVector = std::vector<V>;
#endif//USE_EASTL_VECTOR

template<typename V, size_t Size = 100>
#if USE_EASTL_VECTOR
using RStackString = eastl::fixed_string<V, Size, true>;
#else
using RStackString = std::vector<V>;
#endif//USE_EASTL_VECTOR

///Maps
#define USE_EASTL_MAP (1 && USE_EASTL)
///Unordered-map
template <typename TKey, typename TValue>
#if USE_EASTL_MAP
using RHashTable = eastl::hash_map<TKey, TValue>;
#else
using RHashTable = std::unordered_map<TKey, TValue>;
#endif//USE_EASTL_MAP

///Ordered-map
template <typename TKey, typename TValue>
#if USE_EASTL_MAP
using RHashTableOrdered = eastl::map<TKey, TValue>;
#else
using RHashTableOrdered = std::map<TKey, TValue>;
#endif//USE_EASTL_MAP


///EASTL Exclusive: fixed size containers:
#define USE_EASTL_ARRAY (1 && USE_EASTL)

///Contigious Linked-List
template <typename T, size_t Size = 100>
#if USE_EASTL_ARRAY
using RListArray = eastl::fixed_list<T, Size, true /*allowOverflow*/>;
#else
using RListArray = std::list<T>;
#endif//USE_EASTL

///Contigious Hash-table
template <typename TKey, typename TValue, size_t Size = 100>
#if USE_EASTL_ARRAY
using RHashTableArray = eastl::fixed_hash_map<TKey, TValue, Size>;
#else
using RHashTableArray = std::unordered_map<TKey, TValue>;
#endif//USE_EASTL

///Array
template <typename T, size_t Size>
using RArray = std::array<T, Size>;



enum EColor
{
	Black,
	White,
	Grey,
	Red,
	Green,
	Blue,
	Purple,
	Brown,
	Yellow,
	Orange,
	Pink,
	LightBlue,
	RedDark,
	BlueEgyptian,
	RedOrange,
	YellowCanary,
};

static inline constexpr float3 GetColorRGBFloat(EColor color)
{
	if (color == Red)
	{
		return float3{ 0.8f, 0.1f, 0.1f };
	}
	if (color == Green)
	{
		return float3{ 0.0f, 0.8f, 0.2f };
	}
	if (color == Blue)
	{
		return float3{ 0.1f, 0.1f, 0.8f };
	}
	if (color == Purple)
	{
		return float3{ 0.6f, 0.2f, 0.6f };
	}
	if (color == Brown)
	{
		return float3{ 0.8f, 0.4f, 0.1f };
	}
	if (color == Yellow)
	{
		return float3{ 1.f, 0.9f, 0.1f };
	}
	if (color == Orange)
	{
		return float3{ 0.9f, 0.4f, 0.0f };
	}
	if (color == Pink)
	{
		return float3{ 1.0f, 0.3f, 0.3f };
	}
	if (color == LightBlue)
	{
		return float3{ 0.04f, 0.65f, 0.99f };
	}
	if (color == Black)
	{
		return float3{ 0.05f, 0.05f, 0.05f };
	}
	if (color == White)
	{
		return float3{ 0.99f, 0.99f, 0.99f };
	}
	if (color == Grey)
	{
		return float3{ 0.5f, 0.5f, 0.5f };
	}
	if (color == RedDark)
	{
		return float3{ 0.541f, 0.0f, 0.0f };
	}
	if (color == BlueEgyptian)
	{
		return float3{ 0.06f, 0.2f, 0.65f };
	}
	if (color == RedOrange)
	{
		return float3{ 0.99f, 0.32f, 0.25f };
	}
	if (color == YellowCanary)
	{
		return float3{ 1.f, 1.f, 0.60f };
	}
}

enum ETopology
{
	Line,
	Triangle,
	Point,
	LineStrip,
	TriangleStrip
};


enum class EGlobalUIType
{
	Slider = 0,
	Drag,
	Text,
	Color
};

///Static Variable - can't be destroyed
///Each new type should be placed in SubmitUI() & RenderAllGlobals()
template <typename T>
class RStaticUIVar
{
public:

	static inline std::vector<RStaticUIVar<T>*> StaticVariableRefPool;

	///TODO: Pool should construct this variable, we can only receive const pointer to this var
	RStaticUIVar(const char* name, T value, EGlobalUIType uiType = EGlobalUIType::Drag, float min = 0, float max = 100, float precision = 0.1)
		:
		Name(name),
		Value(value),
		UIType(uiType),
		MinValue(min),
		MaxValue(max),
		Precision(precision)
	{
		StaticVariableRefPool.push_back(this);
		if constexpr (std::is_same_v<T, uint> || std::is_same_v<T, uint2> || std::is_same_v<T, uint3> || std::is_same_v<T, uint4>)
		{
			MinValue = 0;
		}
	}

	RStaticUIVar(const RStaticUIVar&) = delete;
	RStaticUIVar& operator=(const RStaticUIVar&) = delete;

	inline operator T() { return Value; }
	inline operator T& () { return Value; }
	T& operator=(T other)
	{
		Value = other;
		return Value;
	}


	void SubmitUI()
	{
		if constexpr (std::is_same_v<T, int> || std::is_same_v<T, uint>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderInt(Name.c_str(), (int*)&Value, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragInt(Name.c_str(), (int*)&Value, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %i" };
				ImGui::Text(s.data(), Value);
			}
			break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (/*std::is_same_v<T, int2> ||*/ std::is_same_v<T, uint2>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderInt2(Name.c_str(), (int*)&Value.x, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragInt2(Name.c_str(), (int*)&Value.x, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %i, %i" };
				ImGui::Text(s.data(), Value.x, Value.y);
			}
			break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (std::is_same_v<T, float>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderFloat(Name.c_str(), &Value, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragFloat(Name.c_str(), &Value, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %f" };
				ImGui::Text(s.data(), Value);
			}
			break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (std::is_same_v<T, float2>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderFloat2(Name.c_str(), &Value.x, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragFloat2(Name.c_str(), &Value.x, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %f, %f" };
				ImGui::Text(s.data(), Value.x, Value.y);
			}
			break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (std::is_same_v<T, float3>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderFloat3(Name.c_str(), &Value.x, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragFloat3(Name.c_str(), &Value.x, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %f, %f, %f" };
				ImGui::Text(s.data(), Value.x, Value.y, Value.z);
			}
			break;
			case EGlobalUIType::Color:
				ImGui::ColorEdit3(Name.c_str(), &Value.x);
				break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (std::is_same_v<T, float4>)
		{
			switch (UIType)
			{
			case EGlobalUIType::Slider:
				ImGui::SliderFloat4(Name.c_str(), &Value.x, MinValue, MaxValue);
				break;
			case EGlobalUIType::Drag:
				ImGui::DragFloat4(Name.c_str(), &Value.x, Precision, MinValue, MaxValue);
				break;
			case EGlobalUIType::Text:
			{
				const std::string s{ Name + ": %f, %f, %f, %f" };
				ImGui::Text(s.data(), Value.x, Value.y, Value.z, Value.w);
			}
			break;
			case EGlobalUIType::Color:
				ImGui::ColorEdit4(Name.c_str(), &Value.x);
				break;
			default:
				check(false);
				break;
			}
		}
		else if constexpr (std::is_same_v<T, bool>)
		{
			ImGui::Checkbox(Name.c_str(), &Value);
		}
		else
		{
			check(false, "incorrect type");
		}


	}

	EGlobalUIType UIType = EGlobalUIType::Slider;

	std::string Name;

	T Value;
	float MinValue;
	float MaxValue;
	float Precision;

};

///!!! --- Don't declare them in global scope in .h files --- !!!

typedef RStaticUIVar<float> RGlobalFloat;
typedef RStaticUIVar<float2> RGlobalFloat2;
typedef RStaticUIVar<float3> RGlobalFloat3;
typedef RStaticUIVar<float4> RGlobalFloat4;
#define GlobalFloat(name, ...) static RGlobalFloat name{#name, __VA_ARGS__}
#define GlobalFloat2(name, ...) static RGlobalFloat2 name{#name, __VA_ARGS__}
#define GlobalFloat3(name, ...) static RGlobalFloat3 name{#name, __VA_ARGS__}
#define GlobalFloat4(name, ...) static RGlobalFloat4 name{#name, __VA_ARGS__}

typedef RStaticUIVar<int> RGlobalInt;
#define GlobalInt(name, ...) static RGlobalInt name{#name, __VA_ARGS__}

typedef RStaticUIVar<uint> RGlobalUint;
typedef RStaticUIVar<uint2> RGlobalUint2;
#define GlobalUint(name, ...) static RGlobalUint name{#name, __VA_ARGS__}
#define GlobalUint2(name, ...) static RGlobalUint2 name{#name, __VA_ARGS__}

typedef RStaticUIVar<bool> RGlobalBool;
#define GlobalBool(name, ...) static RGlobalBool name{#name, __VA_ARGS__}

static inline void RenderAllGlobals()
{
	if (ImGui::Begin("MathGlobals"))
	{

		for (auto& global : RStaticUIVar<float>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<float2>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<float3>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<float4>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}

		for (auto& global : RStaticUIVar<int>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<uint>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<uint2>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
		for (auto& global : RStaticUIVar<bool>::StaticVariableRefPool)
		{
			global->SubmitUI();
		}
	}
	ImGui::End();
}






static inline unsigned int ipow(unsigned int base, unsigned int exponent) {
	unsigned int result = 1;
	while (exponent > 0) 
	{
		if (exponent & 1) 
		{  // check if exponent is odd
			result *= base;
		}
		base *= base;
		exponent >>= 1;  // divide exponent by 2
	}
	return result;
}




struct Matrix4x4
{
	Matrix4x4()
	{
		for (int r = 0; r < 4; r++)
		{
			for (int c = 0; c < 4; c++)
			{
				if (r == c)//Diagonal
				{
					m[r][c] = 1.f;
				}
				else
				{
					m[r][c] = 0.f;
				}
			}
		}
	}
	union
	{
		struct//Row Major
		{
			float _11, _12, _13, _14;
			float _21, _22, _23, _24;
			float _31, _32, _33, _34;
			float _41, _42, _43, _44;
		};
		float m[4][4];
	};
};



///Collects Pointers to All (Static!) Instances of T
template<typename T>
class RSingletonArray
{

public:

	static std::vector<T*>& GetAllInstances()
	{
		return AllInstances;
	}

protected:
	static inline std::vector<T*> AllInstances;
	int IndexToAllInstancesArray = -1;

	RSingletonArray(T* pInstance)
	{
		int n = AllInstances.size();
		for (int i = 0; i < n; i++)
		{
			if (AllInstances[i] == nullptr)
			{
				AllInstances[i] = pInstance;
				IndexToAllInstancesArray = i;
				break;
			}
		}

		if (IndexToAllInstancesArray == -1)
		{
			IndexToAllInstancesArray = AllInstances.size();
			AllInstances.push_back(pInstance);
		}
	}

	RSingletonArray& operator=(const RSingletonArray&) = delete;
	RSingletonArray(const RSingletonArray&) = delete;

	~RSingletonArray()
	{
		AllInstances[IndexToAllInstancesArray] = nullptr;
	}
};

/*
///No pointer invalidation
template<typename T, uint32 _size>
class RContigiousOwningPoolFiniteSize
{
	int32 Emplace(const T& other)
	{
		if (!NullIndexArray.empty())
		{
			auto index = NullIndexArray.back();
			Pool[index] = other;
			NullIndexArray.pop_back();
		}
		else
		{
			Pool[NextAvailableIndex] = other;
			NextAvailableIndex++;
			check(NextAvailableIndex < _size);
		}
	}

	int32 Erase(int32 index)
	{
		NullIndexArray.push_back(index);
	}


	void ForEach(std::function<void(const T&)> predicate)
	{
		for (int i = 0; i < NextAvailableIndex; i++)
		{
			bool bIsInvalid = false;
			for(int k = 0; k < NullIndexArray.size(); k++)
			{
				if (i == NullIndexArray[k])
				{
					bIsInvalid = true;
					break;
				}
			}
			if (!bIsInvalid)
			{
				predicate(Pool[i]);
			}
		}
	}


	std::vector<uint32> NullIndexArray;
	uint32 NextAvailableIndex = 0;//Real Size
	std::array<T, _size> Pool;//Capacity
};

template <typename T>
class RPointerRef
{
	RIndexedRef(const T& other)
	{
		Pointer = Pool.Emplace(other);
	}
	~RIndexedRef()
	{

	}

	T& Get()
	{

	}

private:
	T* Pointer;

	static inline RContigiousOwningPoolFiniteSize<T, 1000> Pool;
};


template<typename T>
class RContigiousOwningPoolInfiniteSize
{
	int32 Emplace(const T& other)
	{
		if (!NullIndexArray.empty())
		{
			auto index = NullIndexArray.back();
			Pool[index] = other;
			NullIndexArray.pop_back();
		}
		else
		{
			Pool.emplace_back(other);
		}
	}

	int32 Erase(int32 index)
	{
		Pool[index].~T();
		NullIndexArray.push_back(index);
	}


	void ForEach(std::function<void(const T&)> predicate)
	{
		for (int i = 0; i < NextAvailableIndex; i++)
		{
			bool bIsInvalid = false;
			for (int k = 0; k < NullIndexArray.size(); k++)
			{
				if (i == NullIndexArray[k])
				{
					bIsInvalid = true;
					break;
				}
			}
			if (!bIsInvalid)
			{
				predicate(Pool[i]);
			}
		}
	}


	std::vector<uint32> NullIndexArray;
	std::vector<T> Pool;//Capacity
};

template <typename T>
class RIndexedRef
{
	RIndexedRef(const T& other)
	{
		Index = Pool.Emplace(other);
	}
	~RIndexedRef()
	{

	}

	T& Get()
	{

	}

private:
	uint32 Index;

	static inline RContigiousOwningPoolInfiniteSize<T> Pool;
};
*/


static inline float frac(float value)
{
	value = abs(value);
	int i = int(value);
	return value - i;
}


static inline float3 GenerateRandomColor()
{
	return float3{ RandomNumberGenerator::Get().NextFloat(), RandomNumberGenerator::Get().NextFloat(), RandomNumberGenerator::Get().NextFloat() };
}


__forceinline float WrapValuef(const float X, const float Min, const float Max)
{
	if (X < Min)
		return Max - std::fmod((Min - X), (Max - Min));
	else if (X > Max)
		return Min + std::fmod((X - Min), (Max - Min));
}

/** Wraps X to be between Min and Max, inclusive */
template< class T >
__forceinline T WrapValue(const T X, const T Min, const T Max)
{
	static_assert(!std::is_same<T, float>::value, "Use WrapValuef");
	T Size = Max - Min;
	T EndVal = X;
	while (EndVal < Min)
	{
		EndVal += Size;
	}
	while (EndVal > Max)
	{
		EndVal -= Size;
	}
	return EndVal;
}

/** Snaps a value to the nearest grid multiple */
template< class T >
__forceinline T GridSnapValue(T Location, T Grid)
{
	return (Grid == T{}) ? Location : (floor((Location + (Grid / (T)2)) / Grid) * Grid);
}







template <typename V>
using RVector2DLinear = RDynamicVector<RDynamicVector<V>>;


#if 1
static inline uint32 SwizzleZOrder(uint2 coord)
{

#if DEBUG_MEMORY
	check(coord.x < UINT16_MAX && coord.y < UINT16_MAX, "Only 16 bits are supported");
#endif

	static constexpr uint32 MASKS[4] = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF };
	static constexpr uint32 SHIFTS[4] = { 1, 2, 4, 8 };

	uint32 x = coord.x;  // Interleave lower 16 bits of x and y, so the bits of x
	uint32 y = coord.y;  // are in the even positions and bits from y in the odd;

	x = (x | (x << SHIFTS[3])) & MASKS[3];
	x = (x | (x << SHIFTS[2])) & MASKS[2];
	x = (x | (x << SHIFTS[1])) & MASKS[1];
	x = (x | (x << SHIFTS[0])) & MASKS[0];

	y = (y | (y << SHIFTS[3])) & MASKS[3];
	y = (y | (y << SHIFTS[2])) & MASKS[2];
	y = (y | (y << SHIFTS[1])) & MASKS[1];
	y = (y | (y << SHIFTS[0])) & MASKS[0];

	const uint32 result = x | (y << 1);
	return result;
}
#else ///TODO
static inline uint32 SwizzleZOrder(uint2 coord)
{
#if DEBUG_MEMORY
	check(coord.x < UINT16_MAX&& coord.y < UINT16_MAX, "Only 16 bits are supported");
#endif

	__m128i x = _mm_cvtsi32_si128(coord.x);
	__m128i y = _mm_cvtsi32_si128(coord.y);

	__m128i masks = _mm_set_epi32(0x0F0F0F0F, 0x03030303, 0x000F000F, 0x000000FF);
	__m128i shifts = _mm_set_epi32(8, 4, 2, 1);

	// Interleave lower 16 bits of x and y using SIMD instructions
	x = _mm_shuffle_epi8(x, masks);
	y = _mm_shuffle_epi8(y, masks);

	x = _mm_srl_epi32(x, shifts);
	y = _mm_srl_epi32(y, shifts);

	x = _mm_or_si128(x, y);

	x = _mm_mul_epu32(x, _mm_set_epi32(0x00010001, 0x00010001, 0x00010001, 0x00010001));
	x = _mm_and_si128(x, _mm_set_epi32(0xFF, 0xFF00, 0xFF0000, 0xFF000000));

	return _mm_cvtsi128_si32(x);
}
#endif

template<typename V>
class RVector2D
{
public:

	uint2 Size;

	RVector2D() = default;

	RVector2D(uint2 size)
	{
		Resize(size);
	}

	void Resize(uint2 size)
	{
		Size = size;

		///Align size to pow 2
		uint linearSize = pow(RMath::AlignPowerOfTwo(std::max(size.x, size.y) + 1), 2);
		//uint lastIndex = SwizzleZOrder(Size);
		//check(lastIndex < linearSize);
		data.resize(linearSize);
	}
	void Resize(int2 size)
	{
		Resize(uint2{ uint(size.x), uint(size.y) });
	}

	V& GetElement(uint2 index)
	{
#if DEBUG_MEMORY
		check(index.x < Size.x && index.y < Size.y, "Out of bounds read");
#endif
		auto flatIndex = SwizzleZOrder(index);
		return data[flatIndex];
	}
	const V& GetElement(uint2 index) const
	{
#if DEBUG_MEMORY
		check(index.x < Size.x && index.y < Size.y, "Out of bounds read");
#endif
		auto flatIndex = SwizzleZOrder(index);
		return data[flatIndex];
	}

	V& operator [](uint2 index)
	{
		return GetElement(index);
	}
	V& operator [](int2 index)
	{
		return GetElement(uint2{ (uint)index.x, (uint)index.y });
	}

	const V& operator [](uint2 index) const
	{
		return GetElement(index);
	}
	const V& operator [](int2 index) const
	{
		return GetElement(uint2{ (uint)index.x, (uint)index.y });
	}

	const RDynamicVector<V>& GetLinearData() const
	{
		return data;
	}

	auto begin()
	{
		return data.begin();
	}
	auto end()
	{
		return data.end();
	}

private:
	RDynamicVector<V> data;
	
};

template <typename V>
class RVector3D
{
public:
	uint3 Size;

	RVector3D() = default;

	RVector3D(uint3 size)
	{
		Resize(size);
	}

	void Resize(uint3 size)
	{
		Size = size;
		data.resize(size.z);
		for (RVector2D<V>& slice : data)
		{
			slice.Resize(uint2{ size.x, size.y });
		}
	}

	V& GetElement(uint slice, uint2 index)
	{
#if DEBUG_MEMORY
		check(slice < Size.z, "Out of bounds read");
#endif
		return data[slice][index];
	}
	const V& GetElement(uint slice, uint2 index) const
	{
#if DEBUG_MEMORY
		check(slice < Size.z, "Out of bounds read");
#endif
		return data[slice][index];
	}

	V& operator [](uint3 index)
	{
#if DEBUG_MEMORY
		check(index.z < Size.z, "Out of bounds read");
#endif
		return GetElement(index.z, uint2{ index.x, index.y });
	}
	const V& operator [](uint3 index) const
	{
#if DEBUG_MEMORY
		check(index.z < Size.z, "Out of bounds read");
#endif
		return GetElement(index.z, uint2{ index.x, index.y });
	}

	const RVector2D<V>& GetSlice(uint z)
	{
		return data[z];
	}

private:
	RDynamicVector<RVector2D<V>> data;
};


static inline uint32_t SwizzleZOrder(uint3 coord)
{
	uint32_t x = coord.x;
	uint32_t y = coord.y;
	uint32_t z = coord.z;

	// Interleave the bits of x, y, and z to form the Morton code
	x = (x | (x << 16)) & 0x030000FF;
	x = (x | (x << 8)) & 0x0300F00F;
	x = (x | (x << 4)) & 0x030C30C3;
	x = (x | (x << 2)) & 0x09249249;

	y = (y | (y << 16)) & 0x030000FF;
	y = (y | (y << 8)) & 0x0300F00F;
	y = (y | (y << 4)) & 0x030C30C3;
	y = (y | (y << 2)) & 0x09249249;
	y <<= 1;

	z = (z | (z << 16)) & 0x030000FF;
	z = (z | (z << 8)) & 0x0300F00F;
	z = (z | (z << 4)) & 0x030C30C3;
	z = (z | (z << 2)) & 0x09249249;
	z <<= 2;

	return x | y | z;
}

template <typename V>
class RVector3DLinear
{
public:
	uint3 Size;

	RVector3DLinear() = default;

	RVector3DLinear(const uint3& size)
	{
		Resize(size);
	}

	void Resize(const uint3& size)
	{
		Size = size;
		uint linearSize = pow(RMath::AlignPowerOfTwo(std::max(size.x, std::max(size.y, size.z)) + 1), 3);
		data.resize(linearSize);
	}

	V& GetElement(const uint3& index)
	{
#if DEBUG_MEMORY
		check(index.x < Size.x && index.y < Size.y && index.z < Size.z, "Out of bounds read");
#endif
		const uint flatIndex = SwizzleZOrder(index);
		return data[flatIndex];
	}

	const V& GetElement(const uint3& index) const
	{
#if DEBUG_MEMORY
		check(index.x < Size.x&& index.y < Size.y&& index.z < Size.z, "Out of bounds read");
#endif
		const uint flatIndex = SwizzleZOrder(index);
		return data[flatIndex];
	}

	V& operator [](const uint3& index)
	{
#if DEBUG_MEMORY
		check(index.x < Size.x&& index.y < Size.y&& index.z < Size.z, "Out of bounds read");
#endif
		return GetElement(index);
	}
	const V& operator [](const uint3& index) const
	{
#if DEBUG_MEMORY
		check(index.x < Size.x&& index.y < Size.y&& index.z < Size.z, "Out of bounds read");
#endif
		return GetElement(index);
	}


	auto begin()
	{
		return data.begin();
	}
	auto end()
	{
		return data.end();
	}

private:
	RDynamicVector<V> data;
};

///Use in sequential computations performed on single thread to prevent constant heap allocations
///HowTo: auto& vectorRAII = GetTemporaryStorage<float()
template<typename V>
static std::vector<V>& GetTemporaryStorage(uint initialSize = 0)
{
	static std::vector<V> storage;
	if (initialSize > 0)
	{
		storage.resize(initialSize);
	}
	else
	{
		storage.clear();
	}
	return storage;
}





template<typename T>
class RSingleton
{
public:

	static T& Get()
	{
		static T rootsignature;
		return rootsignature;
	}

};
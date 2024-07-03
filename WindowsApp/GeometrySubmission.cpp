#include "DrawSimpleGeometry.h"

#include <execution>

#include "GenerateSimpleGeometry.h"
#include "InterpolationAndCurves.h"
#include "LinearAlgebra.h"
#include "Calculus.h"
#include "Intersections.h"
#include "MappingVisualizer.h"
#include "GPUFunctions.h"
#include "VectorField.h"

#include "Physics.h"


void RSimpleGeometryPass::SubmitGeometry(RGeometry& GImGeometry)
{

	PROFILE_FUNC();

	if (GbFirstBoot)
	{
		
	}

	{
		///Coordinate Space
		GImGeometry.AddCoordinateSpaceVisualization3D(float3(0, 0, 0), float3{ 50,50,50 });
	}

	///Update World Space GUI
	RSpatialControlSet::Get().OnUpdate();
	RSpatialControlSet::Get().GenerateVisualizers(GImGeometry);


	///************************************************ Surface ************************************************
#if 0
	static RSurface Surface;
	
	static const auto RandomColor = GenerateRandomColor();
	GImGeometry.SetColorRaw(RandomColor);
	Surface.Main(RandomColor, GImGeometry);
#endif 


	///************************************************ Vector Field ************************************************
	

#if 0
	static RVectorFieldVisualizer VectorFieldVisualizer;


	static float3 origin{ 1,1,0 };
	GlobalFloat(radius, 2.f);
	GlobalFloat(amplitude, 2.f);
	GlobalFloat(frequency, 2.f);
	GlobalFloat(phase, 2.f);

	if (bFirstBoot)
	{
		VectorFieldVisualizer.b3D = false;

		VectorFieldVisualizer.NumSamples = uint3{ 16, 16, 16 };
		VectorFieldVisualizer.DomainMin = float3{ -2,-2,-2 };
		VectorFieldVisualizer.DomainMax = float3{ 2,2,2 };

		auto& point = RSpatialControlSet::Get().ControlPoints.emplace_back();
		point.pPosition = &origin;

	}

	

	auto vectorFieldfunc = [&](float3 v)->float3
	{
		// Calculate the distance from v to the center of the pipe
		float r = sqrt(v.x * v.x + v.y * v.y);

		// If the distance is greater than the radius, return zero vector
		if (r > radius.Value)
		{
			return float3{ 0, 0, 0 };
		}

		// Calculate the velocity in the z-direction using Poiseuille's formula
		float vz = -amplitude.Value * (radius.Value * radius.Value - r * r) / (4 * frequency.Value);

		// Return a vector with only z-component
		return float3{ vz, 0, 0 };
	};

	

	VectorFieldVisualizer.VisualizeFunction(vectorFieldfunc, GImGeometry);



#endif

#if 0
	static RParametricFunctionVisualizer ParametricFunctionVisualizer;

	GlobalFloat(ParamCircleRadius, 1);

	auto paramFunc = [](float t)
	{
		return float3{ cos(t)* ParamCircleRadius.Value, 0 , sin(t) * ParamCircleRadius.Value };
		return float3{ pow(t, 3.f) /3.f, 0 , t};
	};

	if (GbFirstBoot)
	{
		ParametricFunctionVisualizer.DomainMinMax.x = 0;
		ParametricFunctionVisualizer.DomainMinMax.y = PiX2;
	}

	ParametricFunctionVisualizer.VisualizeFunction(paramFunc, GImGeometry);

	//ParametricFunctionVisualizer.LineIntegralOverVectorField(vectorFieldfunc);

#endif

#if 0
	static RParametricFunctionVisualizer2D SurfaceVisualizer;

	auto paramFunc = [](float2 uv)
	{
		return MapSphericalToCartesian(uv.x, uv.y, 1);
	};

	if (GbFirstBoot)
	{
		SurfaceVisualizer.DomainMin = float2{ 0, 0 };
		SurfaceVisualizer.DomainMax = float2{ PiX2, PiDiv2 };
	}

	SurfaceVisualizer.VisualizeFunction(paramFunc, GImGeometry);


	auto fieldFunc = [&](float3 pos)
	{
		return float3{ pos.y * pos.z, pos.x * pos.z, pos.x * pos.y };
	};

	SurfaceVisualizer.ComputeFlux(fieldFunc);


#endif



#if 0
	static RScalarFunctionVisualizer2D Visualizer;

	GlobalFloat2(GS, float2{ 1,1 });

	auto func = [&](float2 uv)
	{
		auto dotRes = VectorDot(uv, float2{ uv.x * GS.Value.x, uv.y * GS.Value.y });
		return pow(e, dotRes * (-0.5f));
		return pow(e, -uv.x * uv.x + uv.x * uv.y - 2.5f* uv.y* uv.y) / PiX2;
	};

	if (GbFirstBoot)
	{
		Visualizer.InputMin = float2{ -1,-1 };
	}

	Visualizer.Main(func, GImGeometry);
#endif

#if 0

	static RScalarFunctionVisualizer3D ScalarVisualizer;

	auto func = [](float3 vec)
	{
		return vec.y*vec.y;
	};

	if (GbFirstBoot)
	{
		ScalarVisualizer.DomainMin = float3{ 0,0,0 };
		ScalarVisualizer.DomainMax = float3{ 10,2,10 };
	}

	ScalarVisualizer.VisualizeFunction(func, GImGeometry);

#endif

#if 0

	static RMappingFunctionVisualizer MapVisualizer;

	auto mapFunc = [&](float3 vec)
	{
		return float3{ pow(e, vec.x),pow(e, vec.y), pow(e, vec.z) };
		float2 v0 = vec.xy();
		float t = vec.z;
		return float3{ 2 + v0.x * t, -t * t / 2.f + v0.y * t, 0 };
		return float3{ vec.x+(vec.y*vec.y), vec.y + (vec.x), 0 };
	};

	if (bFirstBoot)
	{
		MapVisualizer.bEnabledInputAxis.z = 0;
	}

	MapVisualizer.VisualizeFunction(mapFunc, GImGeometry);


#endif


#if 0

	static VelocityFieldGeneratorFunctor VelocityFieldGenerator;

	VelocityFieldGenerator.GenerateVelocityField2D();


#endif


	///************************************************ Physics ************************************************

	//float deltaT = Time::Get().GetDeltaTimeMs();
	float deltaT = Time::Get().GetDeltaTimeSec();


#if 0
	static RConnectedPointsSimulation Simulation;

	if (GbFirstBoot)
	{

		Simulation.ClothSize = uint2{ 32,32 };
		Simulation.GeneratePoints();


	}
	deltaT = 1.f / 60.f;
	Simulation.SimulateAndRender(GImGeometry, deltaT);

#endif


}

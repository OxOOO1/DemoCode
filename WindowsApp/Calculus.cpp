#include "Calculus.h"

#include "InterpolationAndCurves.h"

RStaticUIVar<float> GlobalSamplingPrecision("GlobalSamplingPrecision", 0.03, EGlobalUIType::Drag, 0.01, 1, 0.01);


void RScalarFunctionVisualizer2D::VisualizeOutput(float3 color, RGeometry& visualizer)
{
	///Merge Input-Output for Plot Visualizer
	const uint2& numSamples = InputSet.Size;
	static RVector2D<float3> visualizationSamples;
	visualizationSamples.Resize(numSamples);

	for (uint x = 0; x < numSamples.x; x++)
	{
		for (uint y = 0; y < numSamples.y; y++)
		{
			uint2 coord{ x, y };
			const auto& curInput = InputSet[coord];
			///Plot Graph
			visualizationSamples[coord] = float3{ curInput.x, OutputSet[coord], curInput.y };
		}

	}

	///Color Visualizer
	float outputMax = visualizationSamples.GetElement(uint2{ 0,0 }).y;
	float outputMin = outputMax;

	if (bVisualizeMagnitudeWithColor)
	{
		///Get highest and lowest magnitude
		for (uint u = 0; u < OutputSet.Size.x; u++)
		{
			for (uint v = 0; v < OutputSet.Size.y; v++)
			{
				uint2 coord{ u, v };

				const float value = visualizationSamples[coord].y;
				outputMax = std::max(outputMax, value);
				outputMin = std::min(outputMin, value);
			}
		}
	}

	for (uint u = 0; u < NumSamples.x-1; u++)
	{
		for (uint v = 0; v < NumSamples.y-1; v++)
		{
			uint2 coord{ u, v };

			c_float3& point = visualizationSamples.GetElement(coord);
			c_float3& pointRight = visualizationSamples.GetElement(coord + uint2{ 1, 0 });
			c_float3& pointUp = visualizationSamples.GetElement(coord + uint2{ 0, 1 });
			c_float3& pointRightUp = visualizationSamples.GetElement(coord + uint2{ 1, 1 });


			if (bVisualizeMagnitudeWithColor)
			{
				///Generate Color from Magnitude
				///Convert cur magnitude to color interpolation parameter
				auto curOutputMagnitude = point.y;
				
				float t = 0;

				if (fabs(outputMax - outputMin) >= 0.01)
				{
					t = MapTo01Range(curOutputMagnitude, outputMin, outputMax);
				}
				else
				{
					t = 0.5;
				}

				///Get color based on cur magnitude
				static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
				static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
				static const float3 RedColor = GetColorRGBFloat(EColor::Red);

				float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

				visualizer.SetColorRaw(CurColor);
			}
			else
			{
				visualizer.SetColorRaw(color);
			}


			const bool bWire = false;
			if (bWire)
			{
				visualizer.AddLine(point, pointRight);
				visualizer.AddLine(point, pointUp);
			}
			else
			{
				visualizer.AddQuadToBuffer(point, pointUp, pointRightUp, pointRight);
			}



			///Draw Gradients
			const bool bDrawGradient = false;
			if (bDrawGradient)
			{
				auto derivative = GradientSet[coord];
				float3 gradient{ derivative.x, 0, derivative.y };
				visualizer.SetColor(EColor::White);
				visualizer.AddArrowPolyboard(point, point + VectorNormalize(gradient));
			}
			
			
		}

	}

}

void RScalarFunctionVisualizer3D::VisualizeUnitSphere(float3 color, RGeometry& visualizer)
{
	RVector2D<float> outputSurface = OutputGrid.GetSlice(0);
	RVector2D<float3> sphereSurfaceCartesian = InputGrid.GetSlice(0);

	///Color Visualizer
	float outputMax = outputSurface[uint2{ 0,0 }];
	float outputMin = outputMax;

	if (bVisualizeMagnitudeWithColor)
	{
		///Get highest and lowest magnitude
		for (uint u = 0; u < outputSurface.Size.x; u++)
		{
			for (uint v = 0; v < outputSurface.Size.y; v++)
			{
				uint2 coord{ u, v };

				const float value = outputSurface[coord];
				outputMax = std::max(outputMax, value);
				outputMin = std::min(outputMin, value);
			}
		}
	}

	for (uint v = 0; v < sphereSurfaceCartesian.Size.x - 1; v++)
	{
		for (uint u = 0; u < sphereSurfaceCartesian.Size.y - 1; u++)
		{
			uint2 coord{ u, v };

			c_float3& point = sphereSurfaceCartesian.GetElement(coord);
			c_float3& pointRight = sphereSurfaceCartesian.GetElement(coord + uint2{ 1, 0 });
			c_float3& pointUp = sphereSurfaceCartesian.GetElement(coord + uint2{ 0, 1 });
			c_float3& pointRightUp = sphereSurfaceCartesian.GetElement(coord + uint2{ 1, 1 });

			auto curOutput = outputSurface[coord];

			if (bVisualizeMagnitudeWithColor)
			{
				///Generate Color from Magnitude
				///Convert cur magnitude to color interpolation parameter
				
				float t = MapTo01Range(curOutput, outputMin, outputMax);

				///Get color based on cur magnitude
				static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
				static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
				static const float3 RedColor = GetColorRGBFloat(EColor::Red);

				float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

				visualizer.SetColorRaw(CurColor);
			}
			else
			{
				visualizer.SetColorRaw(color);
			}


			const bool bWire = false;
			if (bWire)
			{
				visualizer.AddLine(point, pointRight);
				visualizer.AddLine(point, pointUp);
			}
			else
			{
				float outputRight = outputSurface.GetElement(coord + uint2{ 1, 0 });
				float outputUp = outputSurface.GetElement(coord + uint2{ 0, 1 });
				float outputRightUp = outputSurface.GetElement(coord + uint2{ 1, 1 });
				c_float3& modpoint = VectorNormalize(point) * curOutput * SphericalDomainRadiusMin;
				c_float3& modpointRight = VectorNormalize(pointRight) * outputRight * SphericalDomainRadiusMin;
				c_float3& modpointUp = VectorNormalize(pointUp) * outputUp * SphericalDomainRadiusMin;
				c_float3& modpointRightUp = VectorNormalize(pointRightUp) * outputRightUp * SphericalDomainRadiusMin;
				visualizer.AddQuadToBuffer(modpoint, modpointUp, modpointRightUp, modpointRight);
				//visualizer.AddQuadToBuffer(point, pointUp, pointRightUp, pointRight);
			}


		}

	}

}

///Derivative Direction & Magnitude Visualizer (with Color)

void RParametricFunctionVisualizer::VisualizeDerivative(const RDynamicVector<float3>& outputSet, const RDynamicVector<float3> derivativeSet, RGeometry& visualizer, bool bDrawLine)
{
	int numPoints = derivativeSet.size();

	float magnitudeMax = VectorGetLength(derivativeSet[0]);
	float magnitudeMin = magnitudeMax;
	if (bVisualizeMagnitudeWithColor)
	{
		///Get highest and lowest magnitude
		for (int i = 1; i < numPoints; i++)
		{
			magnitudeMax = std::max(magnitudeMax, VectorGetLength(derivativeSet[i]));
			magnitudeMin = std::min(magnitudeMin, VectorGetLength(derivativeSet[i]));
		}
	}

	RDynamicVector<float3> outputPlusDerivativeCurve;
	if (bDrawLine)
	{
		outputPlusDerivativeCurve.resize(numPoints);
	}

	for (int i = 0; i < numPoints; i += DerivativeVisualizationPrecision)
	{
		const float3 derivative = derivativeSet[i];

		///Get final Output and Derivative Direction
		float3 output;
		float3 derivativeDirection;
		bool bNormalizeDerivative = bVisualizeMagnitudeWithColor;
		output = outputSet[i];
		derivativeDirection = bNormalizeDerivative ? VectorNormalize(derivative) : VectorMul(derivative, VisualizedDerivativeMagnitudeScale);

		float derivativeMagnitude = VectorGetLength(derivative);

		if (bVisualizeMagnitudeWithColor)
		{
			///Generate Color from Magnitude
			///Convert cur magnitude to color interpolation parameter
			float t = MapTo01Range(derivativeMagnitude, magnitudeMin, magnitudeMax);

			///Get color based on cur magnitude
			static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
			static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
			static const float3 RedColor = GetColorRGBFloat(EColor::Red);

			float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

			visualizer.SetColorRaw(CurColor);
		}

		visualizer.AddArrowPolyboard(output, output + derivativeDirection);

		if (bDrawLine)
		{
			outputPlusDerivativeCurve[i] = output + derivativeDirection;
		}

	}

	if (bDrawLine)
	{
		visualizer.SetColor(EColor::White);
		visualizer.AddCurvePolyboard(outputPlusDerivativeCurve);
	}


}


void RScalarFunctionVisualizer2D::VisualizeGradientField(RGeometry& visualizer)
{
	///Calculate min/max magnitudes of all vectors
	GradientMax = VectorGetLength(GradientSet.GetElement(uint2{ 0,0 }));
	GradientMin = GradientMax;

	LaplacianMax = LaplacianSet[uint2{ 0,0 }];
	LaplacianMin = LaplacianMax;

	for (uint x = 0; x < NumSamples.x; x++)
	{
		for (uint y = 0; y < NumSamples.y; y++)
		{
			uint2 coord{ x,y };

			float curMagnitude = VectorGetLength(GradientSet[coord]);
			GradientMax = std::max(GradientMax, curMagnitude);
			GradientMin = std::min(GradientMin, curMagnitude);

			float curDivergence = LaplacianSet[coord];
			LaplacianMax = std::max(LaplacianMax, curDivergence);
			LaplacianMin = std::min(LaplacianMin, curDivergence);
		}
	}

	for (uint x = 0; x < NumSamples.x; x++)
	{
		for (uint y = 0; y < NumSamples.y; y++)
		{
			uint2 coord{ x,y };

			///Render Input Grid

			float2 curInput = InputSet[coord];
			float2 curVec = GradientSet[coord];
			float curMagnitude = VectorGetLength(curVec);
			float curDivergence = LaplacianSet[coord];

			///Set Color
			{
				float t = 0;

				auto compute_t = [](float curValue, float minValue, float maxValue)->float
				{
					///Compute t only if there is a significant range
					if (fabs(maxValue - minValue) >= 0.01)
					{
						return MapTo01Range(curValue, minValue, maxValue);
					}
					return 0.5;
				};

				if (bVisualizeGradientDivergence)
				{
					t = compute_t(curDivergence, LaplacianMin, LaplacianMax);
				}
				else
				{
					t = compute_t(curMagnitude, GradientMin, GradientMax);
				}
				

				///Get color based on cur magnitude
				static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
				static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
				static const float3 RedColor = GetColorRGBFloat(EColor::Red);

				float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

				visualizer.SetColorRaw(CurColor);
			}


			///Render direction arrow

			float3 curPos{ curInput.x, 0, curInput.y };
			float3 dir = VectorNormalize(float3{ curVec.x, 0, curVec.y });
			visualizer.AddArrowPolyboard(curPos - dir * SamplingPeriodLength.x * 0.5, curPos + dir * SamplingPeriodLength.x * 0.5);


		}
	}

}


void RScalarFunctionVisualizer3D::VisualizeGradientField(RGeometry& visualizer)
{
	///Calculate min/max magnitudes of all vectors
	float magnitudeMax = VectorGetLength(GradientSet.GetSlice(0).GetElement(uint2{ 0,0 }));
	float magnitudeMin = magnitudeMax;

	float divMax = LaplacianSet[uint3{ 0,0,0 }];
	float divMin = divMax;

	for (uint z = 0; z < NumSamples.z - 3; z++)
	{
		for (uint x = 0; x < NumSamples.x - 3; x++)
		{
			for (uint y = 0; y < NumSamples.y - 3; y++)
			{
				uint3 coord{ x,y,z };

				float curMagnitude = VectorGetLength(GradientSet[coord]);
				magnitudeMax = std::max(magnitudeMax, curMagnitude);
				magnitudeMin = std::min(magnitudeMin, curMagnitude);

				float curDivergence = LaplacianSet[coord];
				divMax = std::max(divMax, curDivergence);
				divMin = std::min(divMin, curDivergence);
			}
		}
	}

	for (uint z = 0; z < NumSamples.z - 1; z++)
	{
		for (uint x = 0; x < NumSamples.x - 1; x++)
		{
			for (uint y = 0; y < NumSamples.y - 1; y++)
			{
				uint3 coord{ x,y,z };

				///Render Input Grid

				float3 curInput = InputGrid[coord];
				float3 curVec = GradientSet[coord];
				float curMagnitude = VectorGetLength(curVec);
				float curDivergence = LaplacianSet[coord];

				///Set Color
				{
					float t = 0;

					auto compute_t = [](float curValue, float minValue, float maxValue)->float
					{
						///Compute t only if there is a significant range
						if (fabs(maxValue - minValue) >= 0.01)
						{
							return MapTo01Range(curValue, minValue, maxValue);
						}
						return 0.5;
					};

					if (bVisualizeGradientDivergence)
					{
						t = compute_t(curDivergence, divMin, divMax);
					}
					else
					{
						t = compute_t(curMagnitude, magnitudeMin, magnitudeMax);
					}

					///Get color based on cur magnitude
					static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
					static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
					static const float3 RedColor = GetColorRGBFloat(EColor::Red);

					float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);

					visualizer.SetColorRaw(CurColor);
				}


				///Render direction arrow
				float3 dir = VectorNormalize(curVec);
				visualizer.AddArrowPolyboard(curInput - dir * SamplingPeriodLength.x * 0.1, curInput + dir * SamplingPeriodLength.x * 0.1);


			}
		}
	}
}
#pragma once

#include "LinearAlgebra.h"
#include "GenerateSimpleGeometry.h"

///Map Local Space to Another Local Space=
/// 
///Map Range 1D
void MapToRangeVisualizer(float x, float xStart, float xEnd, float xStartNew, float xEndNew, RGeometry& visualizer)
{

	float2 Offset{ 0, 10 };

	float xLength = xEnd - xStart;
	float xLengthNew = xEndNew - xStartNew;
	float xStartForNew = Offset.x + xLength + 10;

	///Draw 2 Lines Representing Ranges
	{
		
		visualizer.SetColor(EColor::RedOrange);

		///Original
		visualizer.AddLine(float3(Offset.x, Offset.y, 0), float3(Offset.x + xLength, Offset.y, 0));

		///New
		visualizer.AddLine(float3{ xStartForNew, Offset.y, 0 }, float3(xStartForNew + xLengthNew, Offset.y, 0));
	}

	///Draw 2 Points 
	{
		///Cur

		float3 CurPos = float3{ Offset.x + x, Offset.y, 0 };

		visualizer.SetColor(EColor::Green);

		visualizer.AddCircleCameraFacing(CurPos, 0.5);

		///New
		float xNew = MapToRange(x, xStart, xEnd, xStartNew, xEndNew);

		float3 newPos = float3{ xStartForNew + xNew - xStartNew, Offset.y, 0 };

		visualizer.AddCircleCameraFacing(newPos, 0.5);
	}



}

///Map Range 2D
void MapToRangeVisualizer(c_float2& xy, c_float2& xyStart, c_float2& xyEnd, c_float2& xyStartNew, c_float2& xyEndNew, RGeometry& visualizer)
{

	float2 StartOffset{ 0, 10 };

	float2 xLength = xyEnd - xyStart;
	float2 xLengthNew = xyEndNew - xyStartNew;
	float xStartForNew = StartOffset.x + xLength.x + 10;

	///Draw 2 Planes Representing Ranges
	{
		visualizer.SetColor(EColor::RedOrange);

		///Original
		float2 extent = xLength / 2.f;
		visualizer.AddPlaneLineListXY(float3(StartOffset + extent, 0), extent);

		///New
		float2 extentNew = xLengthNew / 2.f;
		float2 center{ xStartForNew + extentNew.x, StartOffset.y + extentNew.y };

		visualizer.AddPlaneLineListXY(float3{ center, 0 }, extentNew);
	}

	///Draw 2 Points 
	{
		///Cur

		float3 CurPos = float3{ StartOffset + xy, 0 };
		visualizer.SetColor(EColor::Green);
		visualizer.AddCircleCameraFacing(CurPos, 0.5);

		///New
		float2 xNew = MapToRange(xy, xyStart, xyEnd, xyStartNew, xyEndNew);
		float3 newPos = float3{ StartOffset + float2{xStartForNew + xNew.x - xyStartNew.x, xNew.y - xyStartNew.y }, 0 };
		visualizer.AddCircleCameraFacing(newPos, 0.5);
	}

}

///Map Range 3D
void MapToRangeVisualizer(c_float3& xyz, c_float3& xyzStart, c_float3& xyzEnd, c_float3& xyzStartNew, c_float3& xyzEndNew, RGeometry& visualizer)
{

	float3 StartOffset{ 0, 10, 0 };

	float3 xLength = xyzEnd - xyzStart;
	float3 xLengthNew = xyzEndNew - xyzStartNew;
	float xStartForNew = StartOffset.x + xLength.x + 10;

	///Draw 2 Planes Representing Ranges
	{
		visualizer.SetColor(EColor::RedOrange);

		///Original
		float3 extent = xLength / 2.f;
		visualizer.AddCubeLineList(StartOffset + extent, extent);

		///New
		float3 extentNew = xLengthNew / 2.f;
		float3 center{ xStartForNew + extentNew.x, StartOffset.y + extentNew.y, extentNew.z };

		visualizer.AddCubeLineList(center, extentNew);
	}

	///Draw 2 Points 
	{
		///Cur

		float3 CurPos = StartOffset + xyz;
		visualizer.SetColor(EColor::Green);
		visualizer.AddCircleCameraFacing(CurPos, 0.5);

		///New
		float3 xNew = MapToRange(xyz, xyzStart, xyzEnd, xyzStartNew, xyzEndNew);

		float3 newPos = float3{ StartOffset.xy() + float2{xStartForNew + xNew.x - xyzStartNew.x, xNew.y - xyzStartNew.y }, xNew.z - xyzStartNew.z };
		visualizer.AddCircleCameraFacing(newPos, 0.5);
	}

}



class RMappingFunctionVisualizer
{
public:
	RVector3D<float3> InputGrid;
	RVector3D<float3> OutputGrid;

	///Derivatives
	struct PartialDerivatives
	{
		RVector3D<float3> PartialX;
		RVector3D<float3> PartialY;
		RVector3D<float3> PartialZ;
	};
	PartialDerivatives FirstDerivativeSet;
	struct PartialDerivativesMagnitudes
	{
		RVector3D<float> PartialX;
		RVector3D<float> PartialY;
		RVector3D<float> PartialZ;
	};
	PartialDerivativesMagnitudes FirstDerivativeMagnitudeSet;

	uint3 NumSamples;
	float SamplingPeriodLength{1};
	float DomainLength{10};

	float3 StartPosOffset{ 0,0,0 };

	uint3 bEnabledInputAxis{ 1,1,1 };

	///Visualizer Settings
	bool bDrawGridAsPoints = false;
	bool bDrawDefaultGrid = false;
	bool bVisualizeMapping = false;
	bool bVisualizeDerivativeMagnitudes = false;
	bool bVisualizePointOnAFunc = false;

	float VisualizationInterpolation{ 1 };

public:
	void VisualizeFunction(std::function<float3(float3)> mappingFunc, RGeometry& visualizer)
	{
		///UI
		DrawUI();

		///Sampling
		SampleFunction(mappingFunc);

		///Rendering
		Visualize(visualizer);

		if (bVisualizeMapping)
		{
			VisualizeMapping(visualizer);
		}
		if (bVisualizePointOnAFunc)
		{
			VisualizePointOnAFunc(visualizer);
		}

	}

private:
	void DrawUI()
	{
		if (ImGui::Begin("RMappingVisualization"))
		{
			ImGui::Checkbox("xAxis", (bool*)&bEnabledInputAxis.x);
			ImGui::SameLine();
			ImGui::Checkbox("yAxis", (bool*)&bEnabledInputAxis.y);
			ImGui::SameLine();
			ImGui::Checkbox("zAxis", (bool*)&bEnabledInputAxis.z);

			ImGui::Checkbox("AsPoints", &bDrawGridAsPoints);
			ImGui::Checkbox("DrawDefaultGrid", &bDrawDefaultGrid);
			ImGui::Checkbox("VisualizeMapping", &bVisualizeMapping);
			ImGui::Checkbox("VisualizeDerivativeMagnitudes", &bVisualizeDerivativeMagnitudes);
			ImGui::Checkbox("VisualizePointOnAFunc", &bVisualizePointOnAFunc);
			if (bVisualizePointOnAFunc)
			{
				ImGui::DragInt2("Coord", (int*)&PointCoord.x, 1, 0, NumSamples.x-1);
				PointCoord.y = std::min(PointCoord.y, NumSamples.y-1);
				float3 inputPos = InputGrid[PointCoord];
				float3 outputPos = OutputGrid[PointCoord];
				ImGui::Text("InputPos: %f, %f, %f", inputPos.x, inputPos.y, inputPos.z);
				ImGui::Text("OutputPos: %f, %f, %f", outputPos.x, outputPos.y, outputPos.z);
				ImGui::DragFloat3("PointOffset", &PointOffset.x, 0.1, -10.f, 10.f);
				ImGui::Text("PointOffsetAbsolutePosition: %f, %f, %f", PointStandardPos.x, PointStandardPos.y, PointStandardPos.z);
			}



			ImGui::NewLine();

			ImGui::SliderFloat("Interpolation", &VisualizationInterpolation, 0, 1);

			ImGui::DragFloat3("StartOffset", &StartPosOffset.x);

			ImGui::DragFloat("SamplingPeriodLength", &SamplingPeriodLength, 0.1, 0.1, 100);
			ImGui::DragFloat("DomainLength", &DomainLength, SamplingPeriodLength, SamplingPeriodLength);

			ImGui::Text("Num Samples: %u, %u, %u", NumSamples.x, NumSamples.y, NumSamples.z);

		}
		ImGui::End();
	}

	void SampleFunction(std::function<float3(float3)>& mappingFunc)
	{
		uint numSamplesPerAxis = (DomainLength / SamplingPeriodLength) + 1;
		NumSamples.x = bEnabledInputAxis.x > 0 ? numSamplesPerAxis : 1;
		NumSamples.y = bEnabledInputAxis.y > 0 ? numSamplesPerAxis : 1;
		NumSamples.z = bEnabledInputAxis.z > 0 ? numSamplesPerAxis : 1;

		NumSamples.x = std::max(NumSamples.x, 1u);
		NumSamples.y = std::max(NumSamples.y, 1u);
		NumSamples.z = std::max(NumSamples.z, 1u);

		float3 StartPos = StartPosOffset - float3{ DomainLength / 2.f, DomainLength / 2.f, DomainLength / 2.f };

		if (bEnabledInputAxis.z == 0)
		{
			StartPos.z = 0;
		}
		if (bEnabledInputAxis.y == 0)
		{
			StartPos.y = 0;
		}
		if (bEnabledInputAxis.x == 0)
		{
			StartPos.x = 0;
		}

		InputGrid.Resize(NumSamples);
		OutputGrid.Resize(NumSamples);

		for (uint z = 0; z < NumSamples.z; z++)
		{
			for (uint x = 0; x < NumSamples.x; x++)
			{
				for (uint y = 0; y < NumSamples.y; y++)
				{
					uint3 coord{ x,y,z };
					float3 inputPos = StartPos + float3{ (float)x, (float)y, (float)z } *SamplingPeriodLength;

					InputGrid[coord] = inputPos;

					float3 outputPos = mappingFunc(inputPos);

					OutputGrid[coord] = outputPos;

				}
			}
		}


		GenerateDerivatives();

	}


	void GenerateDerivatives()
	{
		FirstDerivativeSet.PartialX.Resize(OutputGrid.Size);
		FirstDerivativeSet.PartialY.Resize(OutputGrid.Size);
		FirstDerivativeSet.PartialZ.Resize(OutputGrid.Size);
		FirstDerivativeMagnitudeSet.PartialX.Resize(OutputGrid.Size);
		FirstDerivativeMagnitudeSet.PartialY.Resize(OutputGrid.Size);
		FirstDerivativeMagnitudeSet.PartialZ.Resize(OutputGrid.Size);
		check(OutputGrid.Size.x == NumSamples.x);

		uint numSlices = (bEnabledInputAxis.z > 0) ? (NumSamples.z-1) : 1;
		uint numStacks = (bEnabledInputAxis.y > 0) ? (NumSamples.y) : 1;
		uint numColumns = (bEnabledInputAxis.x > 0) ? (NumSamples.x) : 1;

		for (uint z = 0; z < numSlices; z++)
		{
			for (uint y = 0; y < numStacks; y++)
			{
				for (uint x = 0; x < numColumns; x++)
				{
					uint3 coord{ x,y,z };

					c_float3& vec = OutputGrid[coord];

					if ((bEnabledInputAxis.x > 0) && ((coord.x + 1) < NumSamples.x))
					{
						c_float3& vecNextX = OutputGrid[coord + uint3{ 1,0,0 }];
						c_float3 partialX = (vecNextX - vec) / SamplingPeriodLength;
						FirstDerivativeSet.PartialX[coord] = partialX;
						FirstDerivativeMagnitudeSet.PartialX[coord] = VectorGetLength(partialX);
					}
					
					if ((bEnabledInputAxis.y > 0) && ((coord.y + 1) < NumSamples.y))
					{
						c_float3& vecNextY = OutputGrid[coord + uint3{ 0,1,0 }];
						c_float3 partialY = (vecNextY - vec) / SamplingPeriodLength;
						FirstDerivativeSet.PartialY[coord] = partialY;
						FirstDerivativeMagnitudeSet.PartialY[coord] = VectorGetLength(partialY);
					}

					if ((bEnabledInputAxis.z > 0) && ((coord.z + 1) < NumSamples.z))
					{
						c_float3& vecNextZ = OutputGrid[coord + uint3{ 0,0,1 }];
						c_float3 partialZ = (vecNextZ - vec) / SamplingPeriodLength;
						FirstDerivativeSet.PartialZ[coord] = partialZ;
						FirstDerivativeMagnitudeSet.PartialZ[coord] = VectorGetLength(partialZ);
					}
				}
			}
		}


	}

	void VisualizeMapping(RGeometry& visualizer)
	{
		///Line between Input-Output Points

		visualizer.EnableDepthTest();

		for (uint z = 0; z < NumSamples.z; z++)
		{
			for (uint x = 0; x < NumSamples.x; x++)
			{
				for (uint y = 0; y < NumSamples.y; y++)
				{
					uint3 coord{ x,y,z };

					///Render Input Grid

					float3 pointCur = InputGrid[coord];

					float3 pointCurMapped = OutputGrid[coord];

					float3 pointCurInterpolated = TLinearInterpolate(VisualizationInterpolation, pointCur, pointCurMapped);

					visualizer.SetColor(EColor::Orange);

					visualizer.AddLine(pointCur, pointCurInterpolated);


				}
			}
		}


		visualizer.DisableDepthTest();
	}


	void Visualize(RGeometry& visualizer)
	{
		
		visualizer.EnableDepthTest();

		///Min-Max derivative magnitude
		float magnitudeMax = FirstDerivativeMagnitudeSet.PartialX.GetSlice(0).GetElement(uint2{ 0,0 });
		float magnitudeMin = magnitudeMax;
		for (uint z = 0; z < NumSamples.z; z++)
		{
			for (uint x = 0; x < NumSamples.x; x++)
			{
				for (uint y = 0; y < NumSamples.y; y++)
				{
					uint3 coord{ x,y,z };

					///For each sample: partialX, partialY, partialZ

					float curMagnitude = (FirstDerivativeMagnitudeSet.PartialX[coord]);
					magnitudeMax = std::max(magnitudeMax, curMagnitude);
					magnitudeMin = std::min(magnitudeMin, curMagnitude);

					if ((bEnabledInputAxis.y > 0))
					{
						curMagnitude = (FirstDerivativeMagnitudeSet.PartialY[coord]);
						magnitudeMax = std::max(magnitudeMax, curMagnitude);
						magnitudeMin = std::min(magnitudeMin, curMagnitude);
					}
					if ((bEnabledInputAxis.z > 0))
					{
						curMagnitude = (FirstDerivativeMagnitudeSet.PartialZ[coord]);
						magnitudeMax = std::max(magnitudeMax, curMagnitude);
						magnitudeMin = std::min(magnitudeMin, curMagnitude);
					}

				}
			}
		}

		for (uint z = 0; z < NumSamples.z; z++)
		{
			for (uint x = 0; x < NumSamples.x; x++)
			{
				for (uint y = 0; y < NumSamples.y; y++)
				{
					uint3 coord{ x,y,z };

					///Render Input Grid

					float3 pointCur = InputGrid[coord];

					float3 pointCurMapped = OutputGrid[coord];

					float3 pointCurInterpolated = TLinearInterpolate(VisualizationInterpolation, pointCur, pointCurMapped);

					if (bDrawGridAsPoints)
					{
						if (bDrawDefaultGrid)
						{
							visualizer.SetColor(EColor::Grey);
							visualizer.AddCircleCameraFacing(pointCur, 0.1);
						}
						
						visualizer.SetColor(EColor::LightBlue);
						visualizer.AddCircleCameraFacing(pointCurInterpolated, 0.1);

					}
					else
					{

						auto RenderLineForAxisNeighbor = [&](uint3 coordOffset, float3 color)
						{
							float3 curNeighbor = InputGrid[coord + coordOffset];

							if (bDrawDefaultGrid)
							{
								visualizer.SetColor(EColor::Grey);
								visualizer.AddLine(pointCur, curNeighbor);
							}
							if (bVisualizeDerivativeMagnitudes)
							{
								visualizer.SetColorRaw(color);
							}
							else
							{
								visualizer.SetColor(EColor::LightBlue);
							}
							float3 curNeighborMapped = OutputGrid[coord + coordOffset];
							float3 pointInterpolated = TLinearInterpolate(VisualizationInterpolation, curNeighbor, curNeighborMapped);
							visualizer.AddLine(pointCurInterpolated, pointInterpolated);
							//visualizer.AddLineCameraFacing(pointCurMapped, pointRightMapped);
						};

						///Get color based on cur magnitude
						static const float3 BlueColor = GetColorRGBFloat(EColor::Blue);
						static const float3 YellowColor = GetColorRGBFloat(EColor::Yellow);
						static const float3 RedColor = GetColorRGBFloat(EColor::Red);

						///x Axis
						if ((coord.x + 1) < NumSamples.x)
						{
							float partialXMagnitude = FirstDerivativeMagnitudeSet.PartialX[coord];
							float t = MapTo01Range(partialXMagnitude, magnitudeMin, magnitudeMax);
							float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);
							RenderLineForAxisNeighbor(uint3{ 1, 0, 0 }, CurColor);
						}
						///y Axis
						if ((coord.y + 1) < NumSamples.y)
						{
							float partialYMagnitude = FirstDerivativeMagnitudeSet.PartialY[coord];
							float t = MapTo01Range(partialYMagnitude, magnitudeMin, magnitudeMax);
							float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);
							RenderLineForAxisNeighbor(uint3{ 0, 1, 0 }, CurColor);
						}
						///z Axis
						if ((coord.z + 1) < NumSamples.z)
						{
							float partialZMagnitude = FirstDerivativeMagnitudeSet.PartialZ[coord];
							float t = MapTo01Range(partialZMagnitude, magnitudeMin, magnitudeMax);
							float3 CurColor = TQuadraticInterpolateBezier(t, BlueColor, YellowColor, RedColor);
							RenderLineForAxisNeighbor(uint3{ 0, 0, 1 }, CurColor);
						}

					}

					

				}
			}
		}
		

		visualizer.DisableDepthTest();

	}

	uint3 PointCoord{ 0,0,0 };
	float3 PointOffset{ 1,0,0 };
	float3 PointStandardPos{ 1,0,0 };

	void VisualizePointOnAFunc(RGeometry& visualizer)
	{
		///Linear approximation at a point
		///f(point) + f'(point)(pointoffset)

		///Jacobian
		float3 BasisX{ 1,0,0 };
		float3 BasisY{0,1,0};
		float3 BasisZ{0,0,1};

		float3 fPoint = OutputGrid[PointCoord];
		float3 input = InputGrid[PointCoord];

		///Visualize Input/Output
		{
			visualizer.SetColor(EColor::RedDark);
			visualizer.AddCircleCameraFacing(input, SamplingPeriodLength * 0.1);

			visualizer.SetColor(EColor::RedOrange);
			visualizer.AddCircleCameraFacing(fPoint, SamplingPeriodLength * 0.1);

			visualizer.SetColor(EColor::Green);
			visualizer.AddLine(input, fPoint);
		}

		if (bEnabledInputAxis.x > 0)
		{
			float3 partialX = FirstDerivativeSet.PartialX[PointCoord];
			auto start = fPoint;
			auto end = fPoint + partialX;
			visualizer.SetColor(EColor::Red);
			visualizer.AddArrowPolyboardDirection(start, VectorNormalize(end-start));

			BasisX = VectorNormalize(end - start);
		}
		
		if (bEnabledInputAxis.y > 0)
		{
			float3 partialY = FirstDerivativeSet.PartialY[PointCoord];
			auto start = fPoint;
			auto end = fPoint + partialY;
			visualizer.SetColor(EColor::YellowCanary);
			visualizer.AddArrowPolyboardDirection(start, VectorNormalize(end - start));

			BasisY = VectorNormalize(end - start);
		}
		
		///Visualize Offset
		//if (bEnabledInputAxis.x > 0 && bEnabledInputAxis.y > 0)
		{
			///Linear Combination of Basis
			PointStandardPos = fPoint + (BasisX * PointOffset.x + BasisY * PointOffset.y);

			visualizer.SetColor(EColor::Pink);
			visualizer.AddCircleCameraFacing(PointStandardPos, SamplingPeriodLength * 0.25);
		}


	}


private:

};



struct RCoordinateSpaceVisualizer
{

	static inline const float3 Right{ 1, 0, 0 };
	static inline const float3 Up{ 0, 1, 0 };
	static inline const float3 Forward{ 0, 0, 1 };
	static inline const float4 Origin{ 0, 0, 0, 1 };

	static void Visualize(const Matrix4x4& matrix, RGeometry& visualizer)
	{
		///Transform Reference Frame with a Matrix
		float3 RightMapped = VectorMul(Right, matrix);
		float3 UpMapped = VectorMul(Up, matrix);
		float3 ForwardMapped = VectorMul(Forward, matrix);
		float4 OriginMapped = VectorMul(Origin, matrix);

		visualizer.SetColor(EColor::Orange);
		visualizer.AddArrowPolyboardDirection(OriginMapped.xyz(), RightMapped);
		visualizer.SetColor(EColor::Purple);
		visualizer.AddArrowPolyboardDirection(OriginMapped.xyz(), UpMapped);
		visualizer.SetColor(EColor::LightBlue);
		visualizer.AddArrowPolyboardDirection(OriginMapped.xyz(), ForwardMapped);
	}

	//RCoordinateSpaceVisualizer* ReferenceFrame;///Parent


};
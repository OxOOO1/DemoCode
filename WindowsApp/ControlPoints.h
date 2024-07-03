#pragma once

#include "MathGeneral.h"
#include "GenerateSimpleGeometry.h"

///************************** World Space Spatial GUI *************************


///TODO: Make it purely Immediate Mode Interface that doesn't own anything
///TODO: Window for a set of control points


struct RControlSetWindow
{

	float3 Center;
	float3 Extent;

	static std::map<std::string, RControlSetWindow> Windows;
};


///Sptially controlled GUI Node
struct RControlPoint
{
public:
	RControlPoint(float3* position = nullptr, float inSize = 0.25f)
	{
		if (position)
		{
			pPosition = position;
		}
		Size = inSize;
		SetColor(EColor::Orange);
	}

	float3* pPosition{nullptr};

	std::function<void(void)> OnPointSelectedCallback;


	bool bHasDirection = false;
	float3 Direction{ 1, 0, 0 };
	float3 Color;
	float Size = 0.25f;

	void SetColor(EColor color)
	{
		Color = GetColorRGBFloat(color);
	}

	uint WindowId;
};

///World Space Interactable GUI allowing to spatially control such Data as Position/Velocity
struct RSpatialControlSet
{

	///TODO: Either force it as Singleton, or define boundaries of each ControlSet window, otherwise multiple Control Sets will overlap
	static RSpatialControlSet& Get()
	{
		static RSpatialControlSet instance;
		return instance;
	}

	RDynamicVector<RControlPoint> ControlPoints;

	RControlPoint& AddPoint(float3* pPosition, float size = 0.25f)
	{
		return ControlPoints.emplace_back(pPosition, size);
	}

	void OnUpdate();

	void GenerateVisualizers(RGeometry& visualizer)
	{
		///Update Colors
		if (pIntersectedPoint)
		{
			pIntersectedPoint->SetColor(EColor::Red);
		}
		if (pSelectedPoint)
		{
			pSelectedPoint->SetColor(EColor::Yellow);
			if (bInDirectionDragState)
			{
				pSelectedPoint->SetColor(EColor::Red);
			}
		}

		///Update Visualizers
		for (RControlPoint& point : ControlPoints)
		{
			///Draw all Visualizers
			
			visualizer.SetColorRaw(point.Color);
			visualizer.AddCircleCameraFacing(*point.pPosition, point.Size);

			if (point.bHasDirection)
			{
				visualizer.SetColor(EColor::Blue);
				visualizer.AddArrowPolyboardDirection(*point.pPosition, *point.pPosition + point.Direction);
			}

			point.SetColor(EColor::Orange);

		}

	}

	///TODO Remove
	/*void ExtractPositions(RDynamicVector<float3>& into)
	{
		check(!RSpatialControlSet::Get().ControlPoints.empty());

		into.clear();

		for (auto& point : RSpatialControlSet::Get().ControlPoints)
		{
			into.push_back(*point.pPosition);
		}

	}*/

	void VisualizeLinesBetweenPoints(RGeometry& visualizer)
	{
		visualizer.SetColor(EColor::White);
		float3& PointPrev = *ControlPoints[0].pPosition;
		for (int i = 1; i < ControlPoints.size(); i++)
		{
			float3& CurPoint = *ControlPoints[i].pPosition;
			visualizer.AddLine(PointPrev, CurPoint, 0.01f);
			PointPrev = CurPoint;
		}
	}

	/*void GenerateNewPoint(float3 pos)
	{
		RControlPoint& newPoint = ControlPoints.emplace_back();
		newPoint.Position = pos;
	}
	void GenerateNewPoint(float3 pos, float3 dir)
	{
		RControlPoint& newPoint = ControlPoints.emplace_back();
		newPoint.Position = pos;

		newPoint.bHasDirection = true;
		newPoint.Direction = dir;
	}*/
	void EraseAllPoints()
	{
		ControlPoints.clear();
	}

	bool IsPointCurrentllyInteractedWith(const RControlPoint& point)
	{
		return &point == pSelectedPoint || &point == pIntersectedPoint;
	}

	static void EraseControlPointFromList(std::list<RControlPoint>& controlPointsList, RControlPoint* pPointToDelete)
	{
		for (std::list<RControlPoint>::iterator it = controlPointsList.begin(); it != controlPointsList.end(); it++)
		{
			if (&(*it) == pPointToDelete)
			{
				controlPointsList.erase(it);
				return;
			}
		}
	}

	void RenderUI()
	{
		//if (ImGui::Button("AllocateNewPoint"))
		//{
		//	//GenerateNewPoint();
		//}
		//if (ImGui::Button("RemoveSelectedPoint"))
		//{
		//	EraseControlPointFromList(ControlPoints, pSelectedPoint);
		//}
		//ImGui::SliderFloat("PointsSize", &ControlPointsSize, 0, 1);

		bool bHasIntersection = pIntersectedPoint != nullptr;
		ImGui::Checkbox("Has Intersection", &bHasIntersection);
	}

	c_float3& GetControllerPosCur() const 
	{
		return CursorWorldPosCur;
	}

//private:

	RSpatialControlSet()
	{
	}

	float3 CursorWorldPosPrev{};//TODO:Move it somewhere
	float3 CursorWorldPosCur{};

	bool bInPositionDragState = false;
	bool bInDirectionDragState = false;
	RControlPoint* pIntersectedPoint = nullptr;
	RControlPoint* pSelectedPoint = nullptr;
	bool bLClickThisFrame = false;
	bool bRClickThisFrame = false;
	float ClickedObjectViewSpaceZ = 10.f;

	
	bool bCursorAsSphere = false;///In order to drag multiple objects
	float CursorSphereRadius = 1.f;

};
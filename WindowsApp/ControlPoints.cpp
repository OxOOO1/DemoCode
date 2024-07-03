#include "ControlPoints.h"

#include "Intersections.h"

void RSpatialControlSet::OnUpdate()
{

	///Construct Mouse Ray
	float3 RayOrigin;
	float3 RayEnd;
	float3 RayDir;
	RScene::Get().GetCamera().GetMouseCursorRayWorldSpace(RayOrigin, RayEnd);
	RayDir = VectorNormalize(RayEnd - RayOrigin);
	///Test Mouse Ray against objects

	if (!bInPositionDragState && !bInDirectionDragState)
	{
		pIntersectedPoint = nullptr;

		for (RControlPoint& point : ControlPoints)
		{
			if (RaySphereIntersection(RayOrigin, RayDir, *point.pPosition, point.Size))
			{
				pIntersectedPoint = &point;
			}
		}
	}

	///Update Cursor Depth
	if (!bInPositionDragState)
	{
		///Depth depends on intersected point depth
		if (pIntersectedPoint != nullptr)
		{
			c_float3& ClickedObjectPos = *pIntersectedPoint->pPosition;
			auto& viewMat = RScene::Get().GetCamera().Matrices.View;
			ClickedObjectViewSpaceZ = float3{ float3SSE{ DirectX::XMVector3Transform(float3SSE{ ClickedObjectPos }, viewMat) } }.z;
		}
		else///If we don't have an intersection, assign some Z value
		{
			//ClickedObjectViewSpaceZ = 100;//Or Far Plane
		}

	}

	///Calculate World Space Mouse Pos
	CursorWorldPosCur = RScene::Get().GetCamera().GetMouseCursorPositionWorld(ClickedObjectViewSpaceZ);
	float3 MouseOffsetWorld = CursorWorldPosCur - CursorWorldPosPrev;
	CursorWorldPosPrev = CursorWorldPosCur;

	///Register mouse click
	bLClickThisFrame = Mouse::Get().IsLButtonPressed();
	bRClickThisFrame = Mouse::Get().IsRButtonPressed();

	///Drag
	if (bInPositionDragState)
	{
		///Offset the object
		*pIntersectedPoint->pPosition = *pIntersectedPoint->pPosition + MouseOffsetWorld;

		if (pIntersectedPoint->OnPointSelectedCallback)
		{
			pIntersectedPoint->OnPointSelectedCallback();
		}

		///Disable Drag if Mouse was released
		if (!bLClickThisFrame)
		{
			bInPositionDragState = false;
		}
	}
	else
	{
		///If Drag was disabled and we clicked -> start Dragging
		if (bLClickThisFrame && pIntersectedPoint != nullptr)
		{
			bInPositionDragState = true;
			pSelectedPoint = pIntersectedPoint;
		}
	}

	///Direction
	/*if (bInDirectionDragState)
	{
		///Offset the object
		pIntersectedPoint->Direction = pIntersectedPoint->Direction + MouseOffsetWorld;

		///Disable Drag if Mouse was released
		if (!bRClickThisFrame)
		{
			bInDirectionDragState = false;
		}
	}
	else
	{
		///If Drag was disabled and we clicked -> start Dragging
		if (bRClickThisFrame && pIntersectedPoint != nullptr)
		{
			bInDirectionDragState = true;
			pSelectedPoint = pIntersectedPoint;
		}
	}*/


	///Generate new point from mouse pos
	/*if (pIntersectedPoint == nullptr)
	{
		if (bLClickThisFrame)
		{
			GenerateNewPoint(CurMousePosWorld);
		}
		else if (bRClickThisFrame)
		{
			GenerateNewPoint(CurMousePosWorld, float3{ 1, 0, 0 });
		}

	}*/


	

}

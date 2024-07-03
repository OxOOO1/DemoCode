#pragma once

#include "MathGeneral.h"


struct RTriangleSimple
{
	float3 V1;
	float3 V2;
	float3 V3;
};


namespace RGeometryGenerator
{

	struct GeometryRef
	{
		RDynamicVector<float3>* VertexBuffer;
		RDynamicVector<float3>* IndexBuffer;
		RDynamicVector<float3>* NormalsBuffer;
		RDynamicVector<float3>* TexCoordBuffer;
	};

}



struct RGeometry
{



	static RGeometry GInstance;


	///================= Description ====================

	struct PipelineStateDesc
	{
		PipelineStateDesc()
		{
			std::memset(this, 0, sizeof(PipelineStateDesc));
			twoSided = 1;
		}

		uint32 GetHash()
		{
			return (depthEnabled) | (twoSided << 1) | (usesNormals << 2) | (usesTexture << 3) | (usesNormalTexture << 4);
		}

		uint bitpad : 27;
		uint usesNormalTexture : 1;
		uint usesTexture : 1;
		uint usesNormals : 1;
		uint twoSided : 1;
		uint depthEnabled : 1;
		//uint shaderId : 16;

		enum
		{
			EDepthEnabled = (1 << 0),
			ETwoSided = (1 << 1),
			EUsesNormals = (1 << 2),
			EUsesTexture = (1 << 3),
			EUsesNormalTexture = (1 << 4)
		};

	};
	static_assert(sizeof(PipelineStateDesc) == sizeof(uint32));

	PipelineStateDesc CurPSOState;

	void EnableDepthTest()
	{
		CurPSOState.depthEnabled = 1;
	}
	void DisableDepthTest()
	{
		CurPSOState.depthEnabled = 0;
	}

	void EnableNormals()
	{
		CurPSOState.usesNormals = 1;
	}
	void DIsableNormals()
	{
		CurPSOState.usesNormals = 0;
	}

	void EnableTexture()
	{
		CurPSOState.usesTexture = 1;
	}
	void DisableTexture()
	{
		CurPSOState.usesTexture = 0;
	}

	void EnableNormalTexture()
	{
		CurPSOState.usesNormalTexture = 1;
	}
	void DisableNormalTexture()
	{
		CurPSOState.usesNormalTexture = 0;
	}

	static PipelineStateDesc GetDescFromHash(uint32 hash)
	{
		PipelineStateDesc desc;
		desc.depthEnabled = (hash & PipelineStateDesc::EDepthEnabled) != 0;
		desc.twoSided = (hash & PipelineStateDesc::ETwoSided) != 0;
		desc.usesNormals = (hash & PipelineStateDesc::EUsesNormals) != 0;
		desc.usesTexture = (hash & PipelineStateDesc::EUsesTexture) != 0;
		desc.usesNormalTexture = (hash & PipelineStateDesc::EUsesNormalTexture) != 0;
		return desc;
	}

	///================= Resources ====================

	struct RenderingResources
	{
		RDynamicVector<float3> VertexBuffer;
		TStructuredBuffer<float3> GpuBufferVertices;

		//TStructuredBuffer<uint> GpuBufferIndices;
		//RDynamicVector<uint> Indices;

		RDynamicVector<float3> Normals;
		TStructuredBuffer<float3> GpuBufferNormals;

		//RDynamicVector<float3> Tangents;
		//TStructuredBuffer<float3> GpuBufferTangents;

		RDynamicVector<float2> TexCoords;
		TStructuredBuffer<float2> GpuBufferTexCoords;

		/*RDynamicVector<float2> IndicesToTextures;
		TStructuredBuffer<float2> GpuBufferIndicesToTextures;*/

		RDynamicVector<float3> ColorsBuffer;
		TStructuredBuffer<float3> GpuBufferColors;

		RTexture* pColorTexture;
		RTexture* pNormalTexture;


		void Clear()
		{
			VertexBuffer.clear();
			ColorsBuffer.clear();
			Normals.clear();
			TexCoords.clear();
			pColorTexture = nullptr;
			pNormalTexture = nullptr;
		}

		void AddTriangle(const float3& p0, const float3& p1, const float3& p2)
		{
			VertexBuffer.push_back(p0);
			VertexBuffer.push_back(p1);
			VertexBuffer.push_back(p2);
		}

		void AddNormals(const float3& n0, const float3& n1, const float3& n2)
		{
			Normals.push_back(n0);
			Normals.push_back(n1);
			Normals.push_back(n2);
		}

		void AddTexCoords(const float2& t0, const float2& t1, const float2& t2)
		{
			TexCoords.push_back(t0);
			TexCoords.push_back(t1);
			TexCoords.push_back(t2);
		}

		void AddTexture(RTexture& texture)
		{
			check(&texture);
			pColorTexture = &texture;
		}

		void AddNormalTexture(RTexture& texture)
		{
			check(&texture);
			pNormalTexture = &texture;
		}

		void AddColor(const float3& color)
		{
			ColorsBuffer.push_back(color);
		}

		template<typename V, class Container>
		static inline void UploadCpuBufferToGpu(RCommandList& CmdList, const Container& cpuBuffer, TStructuredBuffer<V>& gpuBuffer
			, RStringView bufferName = TEXT("ImGeometryBuffer")
		)
		{
			if (cpuBuffer.size() > 0)//Is there anything to upload
			{
				if (cpuBuffer.size() > gpuBuffer.GetElementCount())
				{
					gpuBuffer.Create(bufferName, cpuBuffer.size());
				}
				check(gpuBuffer.GetElementCount() >= cpuBuffer.size());
				//gpuBuffer.UploadElements(0, cpuBuffer.size(), cpuBuffer.data());
				CmdList.UploadToBuffer(gpuBuffer, 0, cpuBuffer.data(), cpuBuffer.size() * sizeof(V));
			}
		}

		void UploadResourcesToGpu(RCommandList& CmdList)
		{
			UploadCpuBufferToGpu(CmdList, VertexBuffer, GpuBufferVertices, TEXT("VertexPositions"));
			UploadCpuBufferToGpu(CmdList, ColorsBuffer, GpuBufferColors, TEXT("VertexColors"));
			UploadCpuBufferToGpu(CmdList, Normals, GpuBufferNormals, TEXT("Normals"));
			UploadCpuBufferToGpu(CmdList, TexCoords, GpuBufferTexCoords, TEXT("TexCoords"));
		}

	};




	///================= Data Submission ====================

	void BeginGeometrySubmission()
	{
		for (auto& [hahs, resources] : PSOStateToResourcesMap)
		{
			resources.Clear();
		}

		EnableDepthTest();
	}

	///Front is clockwise
	void AddTriangleToBuffer(const float3& p0, const float3& p1, const float3& p2)
	{
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddTriangle(p0, p1, p2);
		curResources.AddColor(CurColor);
	}

	void AddQuadToBuffer(const float3& leftDown, const float3& leftUp, const float3& rightUp, const float3& rightDown)
	{
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddTriangle(leftDown, leftUp, rightUp);
		curResources.AddColor(CurColor);
		curResources.AddTriangle(rightUp, rightDown, leftDown);
		curResources.AddColor(CurColor);
	}

	void AddTriangleNormalsToBuffer(const float3& n0, const float3& n1, const float3& n2)
	{
		check(CurPSOState.usesNormals);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddNormals(n0, n1, n2);
	}

	void AddQuadNormalsToBuffer(const float3& leftDown, const float3& leftUp, const float3& rightUp, const float3& rightDown)
	{
		check(CurPSOState.usesNormals);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddNormals(leftDown, leftUp, rightUp);
		curResources.AddNormals(rightUp, rightDown, leftDown);
	}


	void AddTriangleTexCoordsToBuffer(const float2& t0, const float2& t1, const float2& t2)
	{
		check(CurPSOState.usesTexture);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddTexCoords(t0, t1, t2);
	}

	void AddQuadTexCoordsToBuffer(const float2& leftDown, const float2& leftUp, const float2& rightUp, const float2& rightDown)
	{
		check(CurPSOState.usesTexture);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddTexCoords(leftDown, leftUp, rightUp);
		curResources.AddTexCoords(rightUp, rightDown, leftDown);
	}

	void AddTexture(RTexture& texture)
	{
		check(CurPSOState.usesTexture);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddTexture(texture);
	}

	void AddNormalsTexture(RTexture& texture)
	{
		check(CurPSOState.usesNormalTexture);
		RenderingResources& curResources = GetResourceForCurrentState();
		curResources.AddNormalTexture(texture);
	}

	void EndGeometrySubmission(RCommandList& CmdList)
	{
		{
			GPU_PROFILE_SCOPE(UploadGeometryResourceToGPU, CmdList);
			for (auto& [hahs, resources] : PSOStateToResourcesMap)
			{
				resources.UploadResourcesToGpu(CmdList);
			}
		}

		//...
		
	}






	///================= Description Rest ====================


	float3 CurColor{ 1, 0, 1 };
	ETopology CurTopology{ ETopology::Triangle };

	void SetColor(EColor ecolor, float scale = 1.f)
	{
		CurColor = GetColorRGBFloat(ecolor) * scale;
	}
	void SetColorRaw(c_float3& colorRGB, float scale = 1.f)
	{
		CurColor = colorRGB * scale;
	}
	const float3 GetColor() const
	{
		return CurColor;
	}

	
	void SetTopology(ETopology topo)
	{
		CurTopology = topo;
	}

	bool bHighDrawPriority = false; //Render those first



	///================= PSO - Resource Map ====================

	RHashTable<uint32, RenderingResources> PSOStateToResourcesMap;

	RenderingResources& GetResourceFromMap(uint32 key)
	{
		RenderingResources& resource = PSOStateToResourcesMap[key];
		return resource;
	}
	RenderingResources& GetResourceForCurrentState()
	{
		uint32 hash = CurPSOState.GetHash();
		RenderingResources& resource = GetResourceFromMap(hash);
		return resource;
	}

	






	///================= Shapes Generator ====================

	void AddPlaneTextured(c_float3& center, c_float2& extent, RTexture& texture, bool bHorizontalXZPlane = false, RTexture* pNormalTexture = nullptr)
	{
		EnableTexture();
		if (pNormalTexture != nullptr)
		{
			EnableNormalTexture();
		}

		AddTexture(texture);
		if (pNormalTexture != nullptr)
		{
			AddNormalsTexture(*pNormalTexture);
		}

		if (!bHorizontalXZPlane)
		{
			AddPlane(center, extent);
		}
		else
		{
			AddPlaneXZ(center, extent);
		}

		AddQuadTexCoordsToBuffer(float2{ 0.f, 1.f }, float2{ 0.f, 0.f }, float2{ 1.f, 0.f }, float2{ 1.f, 1.f });

		DisableTexture();

		if (pNormalTexture != nullptr)
		{
			DisableNormalTexture();
		}
	};


	void AddPlane(c_float3& center, c_float2& extent)
	{
		float3 p0{ center - float3{extent, 0} };
		float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
		float3 p2{ center + float3{extent, 0} };
		float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

		AddQuadToBuffer(p0, p1, p2, p3);
	};

	void AddPlaneXZ(c_float3& center, c_float2& extent)
	{
		//extent .y is .z
		float3 p0{ center - float3{extent.x, 0, extent.y} };
		float3 p1{ center.x - extent.x, center.y, center.z + extent.y };
		float3 p2{ center + float3{extent.x, 0, extent.y} };
		float3 p3{ center.x + extent.x, center.y, center.z - extent.y };

		AddQuadToBuffer(p0, p1, p2, p3);
	};

	void AddPlaneBounding(class BoundingPlane& plane, float size = 1.f);

	void AddPlaneLineList(c_float3& p0, c_float3& p1, c_float3& p2, c_float3& p3, float lineThickness = 1.f)
	{
		AddLine(p0, p1, lineThickness);
		AddLine(p1, p2, lineThickness);
		AddLine(p2, p3, lineThickness);
		AddLine(p3, p0, lineThickness);
	}

	void AddPlaneLineListXY(c_float3& center, c_float2& extent, float lineThickness = 1.f)
	{
		float3 p0{ center - float3{extent, 0} };
		float3 p1{ center.x - extent.x, center.y + extent.y, center.z };
		float3 p2{ center + float3{extent, 0} };
		float3 p3{ center.x + extent.x, center.y - extent.y, center.z };

		AddLine(p0, p1, lineThickness);
		AddLine(p1, p2, lineThickness);
		AddLine(p2, p3, lineThickness);
		AddLine(p3, p0, lineThickness);
	}

	void AddPlaneLineListXZ(c_float3& center, c_float2& extent, float lineThickness = 1.f)
	{
		float3 p0{ center - float3{extent.x, 0, extent.y} };
		float3 p1{ center.x - extent.x, center.y, center.z + extent.y };
		float3 p2{ center + float3{extent.x, 0, extent.y} };
		float3 p3{ center.x + extent.x, center.y , center.z - extent.y };

		AddLine(p0, p1, lineThickness);
		AddLine(p1, p2, lineThickness);
		AddLine(p2, p3, lineThickness);
		AddLine(p3, p0, lineThickness);
	}

	void AddCubeLineList(c_float3& center, c_float3& extent, float lineThickness = 1.f);

	void AddCircle(c_float3& pos = float3{ 0,0,0 }, float radius = 0.05);
	void AddCircleCameraFacing(c_float3& pos = float3{ 0,0,0 }, float radius = 0.05);

	static void GenerateCircleCameraFacing(RDynamicVector<float3>& outVertexBuffer, c_float3& pos = float3{ 0,0,0 }, float radius = 0.05);

	void AddCircleCameraFacingLineList(c_float3& pos = float3{ 0,0,0 }, float radius = 0.05);

	void AddCapsule(c_float3& spherePos1, c_float3& spherePos2, float radius1, float radius2)
	{
		AddCircleCameraFacingLineList(spherePos1, radius1);
		AddCircleCameraFacingLineList(spherePos2, radius2);

		float3SSE SpherePos1View = spherePos1;
		float3SSE SpherePos2View = spherePos2;

		const auto& viewMat = RScene::Get().GetCamera().Matrices.View;
		SpherePos1View = SpherePos1View.TransformAndReturn(viewMat);
		SpherePos2View = SpherePos2View.TransformAndReturn(viewMat);

		auto dirVec = SpherePos1View - SpherePos2View;

		dirVec = DirectX::XMVectorSetZ(dirVec, 0.0f);
		dirVec.Normalize();

		///Rotate dirVec
		static const float3SSE zaxis{ 0,0,1 };
		auto rotMat = DirectX::XMMatrixRotationAxis(zaxis, PiDiv2);

		dirVec = dirVec.TransformAndReturn(rotMat);

		auto lineStart1 = SpherePos1View + dirVec * radius1;
		auto lineEnd1 = SpherePos2View + dirVec * radius2;

		auto lineStart2 = SpherePos1View - dirVec * radius1;
		auto lineEnd2 = SpherePos2View - dirVec * radius2;

		///Trasform back to world space
		auto viewMatInv = DirectX::XMMatrixInverse(nullptr, viewMat);
		lineStart1 = lineStart1.TransformAndReturn(viewMatInv);
		lineEnd1 = lineEnd1.TransformAndReturn(viewMatInv);
		lineStart2 = lineStart2.TransformAndReturn(viewMatInv);
		lineEnd2 = lineEnd2.TransformAndReturn(viewMatInv);

		AddLine(lineStart1, lineEnd1, 5.f);
		AddLine(lineStart2, lineEnd2, 5.f);
	}

	void AddLine(c_float3& startPos, c_float3& endPos, float widthScale = 1.f);
	void AddLineCameraFacingLong(c_float3& startPos, c_float3& endPos, float widthScale = 1.f, float lineSegmentLengthThreshold = 25.f);

	void AddLineScreenSpace(c_float3& startPosWorld, c_float3& endPosWorld, float ndcWidth = 0.1f);

	void AddTriangleCameraFacing(c_float3& baseCenterPos, c_float3& apexPos, float width = 0.01f);

	void AddLinePolyboardHermite(const float3& p0, const float3& p1, const float3& tangent0, const float3& tangent1, float r);
	void AddLinePolyboard(c_float3& pPrev, c_float3& p0, c_float3& p1, c_float3& pNext, float r);

	void AddCurvePolyboard(const RDynamicVector<float3>& samples, float lineThickness = 1, const RDynamicVector<float3>& colors = RDynamicVector<float3>{});
	void AddCurvePolyboard(const RDynamicVector<float3>& knots, const RDynamicVector<float3>& tangents, float lineThickness = 1)
	{
		const auto numKnots = knots.size();
		for (auto i = 0; i < numKnots - 1; i++)
		{
			AddLinePolyboardHermite(knots[i], knots[i + 1], tangents[i], tangents[i + 1], lineThickness);
		}
	}

	
	void AddArrowPolyboard(c_float3& startPos, c_float3& endPos);
	void AddArrowPolyboard(float2 startPos, float2 endPos)
	{
		AddArrowPolyboard(float3{ startPos, 0 }, float3{ endPos, 0 });
	}
	void AddArrowPolyboardDirection(c_float3& startPos, c_float3& direction, float sizeLimit = 0);


	void AddConvexHullPlane(const RDynamicVector<float3>& points)
	{
		///TODO: Graham Scan

		for (int i = 1, numPoints = points.size(); i < numPoints; i++)
		{
			c_float3& pPrev = points[i - 1];
			c_float3& pCur = points[i];

			AddLine(pPrev, pCur);
		}

		///Connect EndPoints
		AddLine(points.back(), points.front());

	}

	void AddAxisLines3D(c_float3& centerPosRelativeToWorld, float extent);

	void AddCoordinateSpaceVisualization3D(c_float3& centerPosRelativeToWorld, c_float3& extent);


	void AddGrid2D(c_float2& center, c_float2& size, c_float2& periodLength);


	///================= Etc ====================

	void TransformAndReturn(const Matrix4x4& TransformAndReturn);

	void TranslateVertices(const float3& newPos);
	inline void TranslateVertices(const float2& newPos)
	{
		TranslateVertices(float3{ newPos.x, newPos.y, 0 });
	}

	D3D_PRIMITIVE_TOPOLOGY GetTopology()
	{
		if (CurTopology == ETopology::Line)
		{
			return D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		}
		if (CurTopology == ETopology::Triangle)
		{
			return D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		}
		if (CurTopology == ETopology::LineStrip)
		{
			return D3D_PRIMITIVE_TOPOLOGY_LINESTRIP;
		}
		if (CurTopology == ETopology::TriangleStrip)
		{
			return D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP;
		}
		if (CurTopology == ETopology::Point)
		{
			return D3D_PRIMITIVE_TOPOLOGY_POINTLIST;
		}
	}

};

//Tetrahedron: n vertices, n faces, all sides have the same length
struct RTetrahedron : RGeometry
{
	int NumFaces = 16;
	bool bMirrorPolyhedron = false;

	float BaseRadius = 0.5f;
	float Height = 1.f;

	void DrawUI()
	{
		ImGui::Text("===============//Tetrahedr//================");

		ImGui::SliderFloat("Height:", &Height, 1.f, 10.f);
		ImGui::SliderFloat("Base Radius:", &BaseRadius, 0.1f, 10.f);
		ImGui::SliderInt("Num Faces:", &NumFaces, 3, 24);

	}

	void Generate()
	{
		auto& BaseVertices = GetTemporaryStorage<float3>();

		float3 ApexPos = float3{ 0, Height, 0 };
		float3 BaseCenter = bMirrorPolyhedron ? float3{ -ApexPos.x,-ApexPos.y,-ApexPos.z } : float3{ 0,0,0 };

		float3 XYZ;

		float FaceAngle = PiX2 / float(NumFaces);

		//Rotate around Y axis, sampling vertices on XZ plane
		for (float rad = 0; rad < PiX2 + FaceAngle; rad += FaceAngle)
		{
			XYZ.x = cos(rad) * BaseRadius;
			XYZ.z = sin(rad) * BaseRadius;

			BaseVertices.push_back(XYZ);
		}

		//Generate Triangles for sides
		for (int vertexId = 0; vertexId < NumFaces; vertexId++)
		{
			AddTriangleToBuffer(ApexPos, BaseVertices[vertexId], BaseVertices[vertexId + 1]);
		}

		//Generate triangles for Base
		//Generate Triangles for sides
		for (int vertexId = 0; vertexId < NumFaces; vertexId++)
		{
			AddTriangleToBuffer(BaseCenter, BaseVertices[vertexId], BaseVertices[vertexId + 1]);
		}
	}



};




//Line
struct RSphericalCoordsLine : RGeometry
{
	float RadialDistance = 1.f;
	float YAngle = 0;
	float XZAngle = 0;

	float3 CartesianPositionOnSphere;

	void DrawUI()
	{
		ImGui::SliderFloat("RadialDistance", &RadialDistance, 0.1, 100.0);
		//ImGui::SliderFloat("YAngle", &YAngle, -PiDiv2, PiDiv2);
		ImGui::SliderFloat("YAngle", &YAngle, 0, PiX2);
		ImGui::SliderFloat("XZAngle", &XZAngle, 0.f, PiX2);
		ImGui::Text("Position: %f, %f, %f", CartesianPositionOnSphere.x, CartesianPositionOnSphere.y, CartesianPositionOnSphere.z);
	}

	void GenerateSphericalCoordsLine()
	{
		CartesianPositionOnSphere.x = RadialDistance * cos(XZAngle) * cos(YAngle);
		CartesianPositionOnSphere.y = RadialDistance * sin(YAngle);
		CartesianPositionOnSphere.z = RadialDistance * sin(XZAngle) * cos(YAngle);

		AddLine(float3{ 0,0,0 }, CartesianPositionOnSphere);
	}

	void GeneratePolarCoordsLine()
	{
		CartesianPositionOnSphere.x = RadialDistance * cos(YAngle);
		CartesianPositionOnSphere.y = RadialDistance * sin(YAngle);
		CartesianPositionOnSphere.z = 0;

		AddLine(float3{ 0,0,0 }, CartesianPositionOnSphere);
	}

};



struct RViewFrustum : RGeometry
{
	using RGeometry::RGeometry;

	float FoV = Pi / 3.0;//Y
	float Ratio = 1920.0 / 1080.0;
	float NearPlane = 0.1f;
	float FarPlane = 5.f;


	float HalfWidth = 0;
	float HalfHeight = 0;

	void RenderUI()
	{
		ImGui::SliderFloat("FoV", &FoV, 0, Pi);
		ImGui::SliderFloat("NearPlane", &NearPlane, 0.1, 10.0);
		ImGui::SliderFloat("FarPlane", &FarPlane, NearPlane + 0.01, 100.0);
		ImGui::Text("HalfWidth: %f, HalfHeight: %f", HalfWidth, HalfHeight);
	}

	void GenerateForLineList();



};





//Cube
struct RCube : RGeometry
{
	using RGeometry::RGeometry;

	float Size = 1.0f;

	void GenerateForTriangleList();

	void GenerateForLineList();


};


//Sphere
struct RSphere : RGeometry
{
	using RGeometry::RGeometry;

	int SphereSamplingPrecision = 64;

	void DrawUI()
	{
		ImGui::SliderInt("SphereSamplingPrecision", &SphereSamplingPrecision, 1, 128);
	}

	void GenerateSphere();
};








//Icosahedron:12 vertices 20 triangles
struct RIcosahedron : RGeometry
{
	using RGeometry::RGeometry;

	TStructuredBuffer<float3> IcosahedronGeometryBuffer;

	RDynamicVector<RTriangleSimple> IcosahedronTriangles;

	int NumIcosahedronDivisions = 1;

	void DrawUI()
	{
		ImGui::SliderInt("SphereSamplingPrecision", &NumIcosahedronDivisions, 1, 128);
	}

	//Separate each icosahedron triangle into 4
	void TesselateIcosahedron();

	void GenerateIcosahedron();
};


float2 GetPointOnACircleParametric(float radius, float t);

float2 GetPointOnALineParametric(float2 p0, float2 p1, float t);

//static inline float3 GetPointOnTriangle(RTriangleSimple tri, float3 barycentrics)
//{
//	//b1v1 + b2v2 + b3v3
//	return (VectorMul(tri.V1, barycentrics.x) + VectorMul(tri.V2, barycentrics.y) + VectorMul(tri.V3, barycentrics.z));
//}






////////////////////////////////////////////////////////////////////////////////////////////////////


struct RRealtimeSignalGraph : RGeometry
{

	RRealtimeSignalGraph()
	{
		InitPlaybackBuffer();
	}

	RDynamicVector<float> Samples;

	//Every n ms
	float SamplingRate = 30;//Samples Per Sec
	float GraphLengthSecs = 10.f;//Secs

	uint NumSamples = 0;

	bool bCanSample = true;

	void InitPlaybackBuffer()
	{
		NumSamples = SamplingRate * GraphLengthSecs;//Num Samples Per Sec * Num Seconds
		Samples.resize(NumSamples);
	}

	void Sample(float value)
	{
		UpdateTime();

		if (bCanSample)
		{
			*(Samples.end() - 1) = value;

			GenerateGraph();

			Playback();

			bCanSample = false;
		}
	}

	float TimeFromLastSample = 0;
	void UpdateTime()
	{
		float SampleDurationMs = 1000.0 / SamplingRate;

		auto timeCurrent = Time::Get().GetTimeAfterLaunchMs();
		float HowMuchTimeHavePassed = timeCurrent - TimeFromLastSample;

		if (HowMuchTimeHavePassed > SampleDurationMs)
		{
			bCanSample = true;
			TimeFromLastSample = timeCurrent;
		}
	}

	void GenerateGraph()
	{

		float SamplingPrecision = 1.0 / SamplingRate;
		SamplingPrecision = -SamplingPrecision; //So that we start from 

		float x1, x2, y1, y2;

		for (int i = 0; i < NumSamples - 1; i++)
		{
			y1 = Samples[i];
			y2 = Samples[i + 1];

			float xStartIndex = NumSamples - i;
			x1 = xStartIndex * SamplingPrecision;
			x2 = (xStartIndex + 1) * SamplingPrecision;

			AddLine(float3{ x2, y1, 0 }, float3{ x1, y2, 0 });

		}

	}

	void Playback()
	{
		for (uint i = 0; i < NumSamples - 1; i++)
		{
			Samples[i] = Samples[i + 1];
		}
	}


	void DrawUI()
	{
		ImGui::Text("===============//Sampler//================");

		if (ImGui::SliderFloat("SamplingRate", &SamplingRate, 1, 120))
		{
			InitPlaybackBuffer();
		}

	}

};


struct RHistogram2D : RGeometry
{

	float SamplingPrecision = 1.0f;
	float EndPoint = 10.f;

	void GenerateFromFunc(std::function<float(float)> func)
	{
		float thickness = SamplingPrecision / 2.0f;

		for (float t = 0; t < EndPoint; t += SamplingPrecision)
		{
			float y = func(t);

			GenerateThickLine(float3{ t, 0,0 }, float3{ t, y, 0 }, thickness);
		}

	}


	void GenerateThickLine(c_float3& startPos, c_float3& endPos, float thickness)//,float3 RightDirVec
	{
		std::array<float3, 4> vertsT;

		//Clockwise, starting from left down
		vertsT[0] = float3{ startPos.x - thickness, startPos.y, 0 };
		vertsT[1] = float3{ endPos.x - thickness, endPos.y, 0 };
		vertsT[2] = float3{ endPos.x + thickness, endPos.y, 0 };
		vertsT[3] = float3{ startPos.x + thickness, startPos.y, 0 };

		AddQuadToBuffer(vertsT[0], vertsT[1], vertsT[2], vertsT[3]);
	}
};














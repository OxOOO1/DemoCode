#include "DrawSimpleGeometry.h"
#include "GenerateSimpleGeometry.h"

#include "LinearAlgebra.h"
//#include "InterpolationAndCurves.h"
#include "Intersections.h"



RGeometry RGeometry::GInstance;


void RenderVariadicPrecisionGrid3D(float3 center, float3 extent, RGeometry& visualizer)
{

	///Segment Length LOD
		///[Distance to camera, Grid segment length]
	//constexpr float DistanceScale = 2.f;
	static constexpr std::array<float2, 3> LODs =
	{
		//float2{5, 0.5},
		float2{15.f, 1.f},
		float2{30, 2.5},
		float2{50, 5},
		//float2{100, 10},
		//float2{1000, 50},
	};

	float3 gridLength = extent * 2;
	float3 startPos = center - extent;
	float3 endPos = center + extent;


	const auto& camPos = RScene::Get().GetCamera().GetPosition();

	///XY Plane Grid
	if (1 && 1)
	{
		float zOffset = endPos.z;


		float distToCameraXY = std::abs(GetPlanePointDistance(float3{ 0, 0, 1 }, -zOffset, camPos));

		float segmentLengthXY = LODs.back().y;

		for (int i = 0; i < LODs.size(); i++)
		{
			if (distToCameraXY < LODs[i].x)
			{
				segmentLengthXY = LODs[i].y;
				break;
			}
		}

		float3 gidLengthInPoints = (gridLength / segmentLengthXY) + float3{ 1,1,1 };

		visualizer.SetColorRaw(float3{ 0.1,0.1,0.1 });

		float yOffset = extent.y;
		startPos.y += yOffset;
		endPos.y += yOffset;

		for (int y = 0; y < gidLengthInPoints.y; y++)
		{
			float yPos = startPos.y + y * segmentLengthXY;
			float3 pointPrev = float3{ startPos.x, yPos, zOffset };
			float3 pointCur = float3{ endPos.x, yPos, zOffset };

			visualizer.AddLineCameraFacingLong(pointPrev, pointCur);
		}

		for (int x = 0; x < gidLengthInPoints.x; x++)
		{
			float xPos = startPos.x + x * segmentLengthXY;
			float3 pointPrev = float3{ xPos, startPos.y, zOffset };
			float3 pointCur = float3{ xPos, endPos.y, zOffset };

			visualizer.AddLineCameraFacingLong(pointPrev, pointCur);
		}
	}

	///XZ Plane Grid
	{
		float distToCameraXZ = std::abs(GetPlanePointDistance(float3{ 0, 1, 0 }, 0, camPos));

		float segmentLengthXZ = LODs.back().y;

		for (int i = 0; i < LODs.size(); i++)
		{
			if (distToCameraXZ < LODs[i].x)
			{
				segmentLengthXZ = LODs[i].y;
				break;
			}
		}

		float3 gidLengthInPoints = (gridLength / segmentLengthXZ) + float3{ 1,1,1 };

		///From Left To Right

		for (int z = 0; z < gidLengthInPoints.z; z++)
		{
			float zPos = startPos.z + z * segmentLengthXZ;
			float3 pointPrev = float3{ startPos.x, center.y, zPos };
			float3 pointCur = float3{ endPos.x, center.y, zPos };

			visualizer.AddLineCameraFacingLong(pointPrev, pointCur);
		}

		for (int x = 0; x < gidLengthInPoints.x; x++)
		{
			float xPos = startPos.x + x * segmentLengthXZ;
			float3 pointPrev = float3{ xPos, center.y, startPos.z };
			float3 pointCur = float3{ xPos, center.y, endPos.z };

			visualizer.AddLineCameraFacingLong(pointPrev, pointCur);
		}
	}

}


void RGeometry::AddPlaneBounding(BoundingPlane& plane, float size)
{
	float3 pointOnPlane = plane.GetPointOnPlaneClosestToOrigin();
	float3 n = plane.GetNormal();


	///Default unit plane
	RArray<float3, 4> vertices =
	{
		float3{-size, -size,0},
		float3{size, -size,0},
		float3{size, size,0},
		float3{-size, size,0},
	};


	Matrix4x4 mat;

	if(1)///Matrix from a basis
	{
		///Generate Basis Vectors
		///Build Rotation Matrix from Basis Vectors

		float3 right;
		float3 up;

		if (fabs(n.y) == 1.f)
		{
			float3SSE rightxm{ n };
			rightxm = DirectX::XMVector3Orthogonal(rightxm);

			right = rightxm;
			up = VectorNormalize(VectorCross(right, n));
		}
		else
		{
			right = VectorNormalize(VectorCross(n, float3{ 0,1.f,0 }));
			up = VectorNormalize(VectorCross(right, n));
		}

		

		mat = MatrixBasisTransform(float3{ 0,0,0 }, n, up, right);
	}
	else
	{
		///initial spherical coords pos
		vertices =
		{
			float3{0.f, -size,-size},
			float3{0.f, -size,size},
			float3{0.f, size,size},
			float3{0.f, size,-size},
		};
		///Build yaw pitch rotation matrix from normal direction

		float yaw = 0;
		float pitch = 0;
		MapCartesianToSpherical(n, yaw, pitch);

		mat = MatrixOrientationEuler(pitch, yaw, 0.f);
	}

	///Rotate default plane
	for (auto& v : vertices)
	{
		///Rotate
		v = VectorMul(v, mat);
		///Translate
		v = v + pointOnPlane;
	}

	SetColor(EColor::Purple);
	
	//AddQuadToBuffer(vertices[0], vertices[1], vertices[2], vertices[3]);
	AddPlaneLineList(vertices[0], vertices[1], vertices[2], vertices[3]);

	AddArrowPolyboard(pointOnPlane, pointOnPlane + n);


}

void RGeometry::AddGrid2D(c_float2& center, c_float2& size, c_float2& periodLength)
{
	

}

void RGeometry::AddAxisLines3D(c_float3& pos, float extent)
{
	///Add Lines
	//AxisLines.SetColorRaw(float3{ 0.5, 0.5, 0.5 });
	//AxisLines.SetColor(EColor::LightBlue);
	SetColor(EColor::White);
	float lineWidth = 1.f;
	AddLineCameraFacingLong(float3{pos.x - extent, pos.y, pos.z }, pos, lineWidth);
	AddLineCameraFacingLong(pos, float3{ pos.x + extent, pos.y, pos.z }, lineWidth);
	AddLineCameraFacingLong(float3{ pos.x, pos.y- extent, pos.z }, pos, lineWidth);
	AddLineCameraFacingLong(pos, float3{ pos.x, pos.y+ extent, pos.z }, lineWidth);
	AddLineCameraFacingLong(float3{ pos.x, pos.y, pos.z - extent }, pos, lineWidth);
	AddLineCameraFacingLong(pos, float3{ pos.x, pos.y, pos.z + extent }, lineWidth);

	///Add Arrows
	SetColor(EColor::Red);
	AddArrowPolyboardDirection(pos, float3{ 1, 0, 0 }, extent / 2);
	SetColor(EColor::Green);
	AddArrowPolyboardDirection(pos, float3{ 0, 1, 0 }, extent / 2);
	SetColor(EColor::Blue);
	AddArrowPolyboardDirection(pos, float3{ 0, 0, 1 }, extent / 2);

}

void RGeometry::AddCoordinateSpaceVisualization3D(c_float3& center, c_float3& extent)
{
	//RenderVariadicPrecisionGrid3D(center - float3{ 0,-1,0 }, extent, *this);
	AddAxisLines3D(center, extent.x);
}


float2 GetPointOnACircleParametric(float radius, float t)
{
	return float2{ cos(t * PiX2) * radius, sin(t * PiX2) * radius };
}

float2 GetPointOnALineParametric(float2 p0, float2 p1, float t)
{
	return p0 + VectorMul(p1 - p0, t);
}


GlobalFloat(CameraFacingObjectScale, 0.001f, EGlobalUIType::Drag, 0.001, 1, 0.001);
GlobalFloat(CameraFacingObjectScalePower, 1.1f, EGlobalUIType::Drag, -1, 5, 0.001);

static inline float GetScaleOfObjectFacingCamera(c_float3& vPos, float scaleAdditional)
{
#if 0
	return 0.1 * scaleAdditional;
#else
	c_float3& camPos = RScene::Get().GetCamera().GetPosition();
	float distanceToCam = VectorGetLength(vPos - camPos);
	float scale = pow(distanceToCam, CameraFacingObjectScalePower.Value) * CameraFacingObjectScale.Value;
	return scale * scaleAdditional;
#endif
}

GlobalFloat(DirectionArrowScale, 100, EGlobalUIType::Drag, 1, 500, 0.1);

void RGeometry::AddArrowPolyboardDirection(c_float3& startPos, c_float3& direction, float sizeLimit)
{
	///Compute EndPos based on distance from camera
	float scale = GetScaleOfObjectFacingCamera(startPos, DirectionArrowScale.Value);

	if (sizeLimit != 0 && scale > sizeLimit)
	{
		scale = sizeLimit;
	}

	float3 endPos = startPos + direction * scale;

	///Generate Board line with Board Triangle
	//Pointer triangle takes 20% of the whole line
	float t = 0.8f;
	float3 LineEndPos = startPos * (1.f - t) + endPos * t;

	float thickness = 2;

	AddLine(startPos, LineEndPos, thickness);

	AddTriangleCameraFacing(LineEndPos, endPos, thickness * GoldenRatioX2);
}

void RSphere::GenerateSphere()
{
	RDynamicVector<RDynamicVector<float3>> StacksSectors;
	StacksSectors.resize(SphereSamplingPrecision + 1);
	for (auto& sector : StacksSectors)
	{
		sector.resize(SphereSamplingPrecision + 1);
	}

	float SphereSamplingPrecisionFloat = SphereSamplingPrecision;

	///Divide sphere into Sectors (X) and Stacks (Y)
	///Sector angle theta
	///Stack angle phi

	///theta = sectorIndex / numSectors
	///phi = stackIndex / numStacks

	float3 XYZ;
	float3 XYZPrev;

	float sectorStep = PiX2 / SphereSamplingPrecisionFloat;///full circle
	float stackStep = Pi / SphereSamplingPrecisionFloat;///half circle

	float theta, phi;

	for (int stackIndex = 0; stackIndex <= SphereSamplingPrecision; stackIndex++)
	{
		phi = Pi / 2 - (stackIndex * stackStep);//starting from 90 degrees
		XYZ.y = sin(phi);

		///add sector count + 1 vertices per stack
		for (int sectorIndex = 0; sectorIndex < SphereSamplingPrecision; sectorIndex++)
		{
			theta = sectorIndex * sectorStep;

			XYZ.x = cos(phi) * cos(theta);
			XYZ.z = cos(phi) * sin(theta);

			StacksSectors[stackIndex][sectorIndex] = XYZ;


		}
	}


	//Generate Indices
	for (int stackIndex = 0; stackIndex < SphereSamplingPrecision; stackIndex++)
	{
		for (int sectorIndex = 0; sectorIndex < SphereSamplingPrecision; sectorIndex++)
		{
			//2 triangles per sector

			float3 leftUp = StacksSectors[stackIndex][sectorIndex];
			float3 leftDown = StacksSectors[stackIndex + 1][sectorIndex];
			float3 rightUp = StacksSectors[stackIndex][(sectorIndex + 1) % (SphereSamplingPrecision)];

			if (stackIndex != 0)
			{
				AddTriangleToBuffer(leftDown, leftUp, rightUp);
			}


			float3 rightDown = StacksSectors[stackIndex + 1][(sectorIndex + 1) % (SphereSamplingPrecision)];

			if (stackIndex != SphereSamplingPrecision - 1)
			{
				AddTriangleToBuffer(leftDown, rightUp, rightDown);
			}

			// vertex tex coord (s, t) range between [0, 1]
			/*s = (float)j / sectorCount;
			t = (float)i / stackCount;
			texCoords.push_back(s);
			texCoords.push_back(t);*/

		}
	}

}


//Separate each icosahedron triangle into 4

void RIcosahedron::TesselateIcosahedron()
{
	RDynamicVector<RTriangleSimple> IcosahedronTrianglesTesselated;

	auto NumOfTrianglesOriginal = IcosahedronTriangles.size();
	for (int t = 0; t < NumOfTrianglesOriginal; t++)
	{
		RTriangleSimple NewTriangle;
		auto& triangle = IcosahedronTriangles[t];

		// compute 3 new vertices by spliting half on each edge
		//         v1       
		//        / \       
		// newV1 *---* newV3
		//      / \ / \     
		//    v2---*---v3   
		//       newV2      
		float3 NewVertex1 = GetHalfVertex(triangle.V1, triangle.V2);
		float3 NewVertex2 = GetHalfVertex(triangle.V2, triangle.V3);
		float3 NewVertex3 = GetHalfVertex(triangle.V3, triangle.V1);

		// add 4 new triangles to vertex array
		IcosahedronTrianglesTesselated.push_back(RTriangleSimple{ triangle.V1, NewVertex1, NewVertex3 });
		IcosahedronTrianglesTesselated.push_back(RTriangleSimple{ NewVertex1, triangle.V2, NewVertex2 });
		IcosahedronTrianglesTesselated.push_back(RTriangleSimple{ NewVertex3, NewVertex2, triangle.V3 });
		IcosahedronTrianglesTesselated.push_back(RTriangleSimple{ NewVertex1, NewVertex2, NewVertex3 });
	}

	IcosahedronTriangles = std::move(IcosahedronTrianglesTesselated);


	RDynamicVector<float3> geometryBufferCPU;
	for (auto& tri : IcosahedronTriangles)
	{

		geometryBufferCPU.push_back(tri.V1);
		geometryBufferCPU.push_back(tri.V2);
		geometryBufferCPU.push_back(tri.V3);

	}
	IcosahedronGeometryBuffer.UploadElements(0, geometryBufferCPU.size(), geometryBufferCPU.data());

}

void RIcosahedron::GenerateIcosahedron()
{
	RDynamicVector<float3> geometryBufferCPU;

	float VertexLongitude{ float(atan(0.5)) }; //Angle From center to edge vertex
	float VertexLatitude{ PiX2 / 5 };//(5 angles per circle)
	float FirstRowStartAngle = -PiDiv2;
	float SecondRowOffset = VertexLatitude / 2.;
	float SecondRowStartAngle = -PiDiv2 - SecondRowOffset;
	float Radius = 1.f;

	float YPos = sin(VertexLongitude);

	float3 ApexVertexPos{ 0,1,0 };
	float3 BaseVertexPos{ ApexVertexPos.x, -ApexVertexPos.y, ApexVertexPos.z };

	RDynamicVector<float3> verticesFirstRow;
	RDynamicVector<float3> verticesSecondRow;

	//Generate Edge vertices
	for (int VertexId = 0; VertexId < 5; VertexId++)
	{
		float3 Vertex1;
		float3 Vertex2;

		Vertex1.x = cos(VertexLongitude) * cos(FirstRowStartAngle);
		Vertex1.y = YPos;//Near Apex
		Vertex1.z = cos(VertexLongitude) * sin(FirstRowStartAngle);

		Vertex2.x = cos(VertexLongitude) * cos(SecondRowStartAngle);
		Vertex2.y = -YPos;//Near Base
		Vertex2.z = cos(VertexLongitude) * sin(SecondRowStartAngle);

		verticesFirstRow.push_back(Vertex1);
		verticesSecondRow.push_back(Vertex2);

		FirstRowStartAngle += VertexLatitude;
		SecondRowStartAngle += VertexLatitude;

	}

	//Generate Triangles

	//5 triangles from Apex
	for (int i = 0; i < 5; i++)
	{
		RTriangleSimple tri;
		tri.V1 = ApexVertexPos;

		tri.V2 = verticesFirstRow[(i + 1) % 5];
		tri.V3 = verticesFirstRow[i];

		IcosahedronTriangles.push_back(tri);
	}

	//5 triangles from Upper vertices
	for (int i = 0; i < 5; i++)
	{
		RTriangleSimple tri;
		tri.V1 = verticesFirstRow[i];

		tri.V2 = verticesSecondRow[(i + 1) % 5];
		tri.V3 = verticesSecondRow[i];

		IcosahedronTriangles.push_back(tri);
	}

	//5 triangles from bottom vertices
	for (int i = 0; i < 5; i++)
	{
		RTriangleSimple tri;
		tri.V1 = verticesFirstRow[i];
		tri.V2 = verticesSecondRow[(i + 1) % 5];
		tri.V3 = verticesFirstRow[(i + 1) % 5];

		IcosahedronTriangles.push_back(tri);
	}

	//5 triangles from bottom vertex
	for (int i = 0; i < 5; i++)
	{
		RTriangleSimple tri;
		tri.V1 = BaseVertexPos;

		tri.V2 = verticesSecondRow[(i + 1) % 5];
		tri.V3 = verticesSecondRow[i];

		IcosahedronTriangles.push_back(tri);
	}

	for (auto& tri : IcosahedronTriangles)
	{

		geometryBufferCPU.push_back(tri.V1);
		geometryBufferCPU.push_back(tri.V2);
		geometryBufferCPU.push_back(tri.V3);



	}


	IcosahedronGeometryBuffer.UploadElements(0, geometryBufferCPU.size(), geometryBufferCPU.data());



}

static int CircleGeneratorPrecision = 16;

const RDynamicVector<float3>& GetGeneratedCircleDefault()
{
	float CircleSamplingPrecisionFLoat = CircleGeneratorPrecision;

	static RDynamicVector<float3> CachedCircleVertexBuffer;

	if (CachedCircleVertexBuffer.empty())
	{
		float samplingPeriodLength = PiX2 / float(CircleGeneratorPrecision);
		float rad = 0;
		for (int i = 0; i <= CircleGeneratorPrecision; i++)
		{
			float3 pos;
			pos.x = cos(rad);
			pos.y = sin(rad);
			pos.z = 0;

			CachedCircleVertexBuffer.emplace_back(pos);

			rad += samplingPeriodLength;
		}
	}

	return CachedCircleVertexBuffer;

}

void RGeometry::AddCircle(c_float3& offset, float radius)
{
	auto& circleBuffer = GetGeneratedCircleDefault();
	auto circleSize = circleBuffer.size();

	for (int i = 1; i < circleSize; i++)
	{
		AddTriangleToBuffer(offset, (circleBuffer[i-1]*radius) + offset, (circleBuffer[i] * radius) + offset);
	}

}

void RGeometry::AddCircleCameraFacing(c_float3& offset, float radius)
{
	auto& circleBuffer = GetGeneratedCircleDefault();
	auto circleSize = circleBuffer.size();

	float3 n = RScene::Get().GetCamera().GetDirection();
	float3 right = RScene::Get().GetCamera().GetRightVec();
	float3 up = VectorCross(n, right);

	auto mat = MatrixBasisTransform(float3{ 0,0,0 }, n, up, right);

	float3 posPrev = VectorMul(circleBuffer[0] * radius, mat) + offset;
	
	for (int i = 1; i < circleSize; i++)
	{
		///Rotate + Translate Default Circle
		float3 posCur = VectorMul(circleBuffer[i] * radius, mat) + offset;

		AddTriangleToBuffer(offset, posPrev , posCur);

		posPrev = posCur;
	}

}

void RGeometry::GenerateCircleCameraFacing(RDynamicVector<float3>& outVertexBuffer, c_float3& offset, float radius)
{
	outVertexBuffer.clear();

	auto& circleBuffer = GetGeneratedCircleDefault();
	auto circleSize = circleBuffer.size();

	float3 n = RScene::Get().GetCamera().GetDirection();
	float3 right = RScene::Get().GetCamera().GetRightVec();
	float3 up = VectorCross(n, right);

	auto mat = MatrixBasisTransform(float3{ 0,0,0 }, n, up, right);

	float3 posPrev = VectorMul(circleBuffer[0] * radius, mat) + offset;

	for (int i = 1; i < circleSize; i++)
	{
		///Rotate + Translate Default Circle
		float3 posCur = VectorMul(circleBuffer[i] * radius, mat) + offset;

		outVertexBuffer.emplace_back(offset);
		outVertexBuffer.emplace_back(posPrev);
		outVertexBuffer.emplace_back(posCur);

		posPrev = posCur;
	}
}

void RGeometry::AddCircleCameraFacingLineList(c_float3& pos, float radius)
{
	///Generate Circle Samples
	auto& circleBuffer = GetGeneratedCircleDefault();
	auto circleNumSamples = circleBuffer.size();

	float3 n = RScene::Get().GetCamera().GetDirection();
	float3 right = RScene::Get().GetCamera().GetRightVec();
	float3 up = VectorCross(n, right);

	auto mat = MatrixBasisTransform(float3{ 0,0,0 }, n, up, right);

	static RDynamicVector<float3> TransformedCircleSampes;
	TransformedCircleSampes.resize(circleNumSamples);

	for (size_t i = 0; i < circleNumSamples; i++)
	{
		TransformedCircleSampes[i] = VectorMul(circleBuffer[i] * radius, mat) + pos;
	}

	AddCurvePolyboard(TransformedCircleSampes);
}

void RGeometry::AddLine(c_float3& startPos, c_float3& endPos, float widthScale)
{
#if 0
	AddLineScreenSpace(startPos, endPos, 0.0025);

#else
	c_float3 dirVec = endPos - startPos;

	c_float3& camPos = RScene::Get().GetCamera().GetPosition();

	float3 dirToCamFromStart = VectorNormalize(startPos - camPos);
	float3 dirToCamFromEnd = VectorNormalize(endPos - camPos);

	float3 lineNormalStart = VectorCross(dirVec, dirToCamFromStart);
	float3 lineNormalEnd = VectorCross(dirVec, dirToCamFromEnd);

	float widthStart = GetScaleOfObjectFacingCamera(startPos, widthScale);
	float widthEnd = GetScaleOfObjectFacingCamera(endPos, widthScale);

	lineNormalStart = VectorMul(VectorNormalize(lineNormalStart), widthStart);
	lineNormalEnd = VectorMul(VectorNormalize(lineNormalEnd), widthEnd);

	std::array<float3, 4> vertsT;

	//Clockwise, starting from left down
	vertsT[0] = startPos - lineNormalStart;
	vertsT[1] = endPos - lineNormalEnd;
	vertsT[2] = endPos + lineNormalEnd;
	vertsT[3] = startPos + lineNormalStart;

	AddTriangleToBuffer(vertsT[0], vertsT[1], vertsT[2]);
	AddTriangleToBuffer(vertsT[0], vertsT[2], vertsT[3]);
#endif

}

void RGeometry::AddLineScreenSpace(c_float3& startPosWorld, c_float3& endPosWorld, float widthScaleNDC)
{


	GlobalFloat(NDCScale, 0.0025, EGlobalUIType::Drag, 0.001, 1, 0.001);

	widthScaleNDC = NDCScale.Value;


	///Convert Start-End Points to NDC Space

	float4SSE startPosNDCsse = float4{ startPosWorld, 1 };
	float4SSE endPosNDCsse = float4{ endPosWorld, 1 };

	auto viewProjMat = RScene::Get().GetCamera().Matrices.ViewProjection;

	startPosNDCsse = startPosNDCsse.TransformAndReturn(viewProjMat);
	endPosNDCsse = endPosNDCsse.TransformAndReturn(viewProjMat);

	float4 startPosNDC = startPosNDCsse;
	float4 endPosNDC = endPosNDCsse;

	startPosNDC = float4{ startPosNDC.x / startPosNDC.w, startPosNDC.y / startPosNDC.w, startPosNDC.z / startPosNDC.w, 1 };
	endPosNDC = float4{ endPosNDC.x / endPosNDC.w, endPosNDC.y / endPosNDC.w, endPosNDC.z / endPosNDC.w, 1 };

	float2 dirVecNDC = VectorNormalize(endPosNDC.xy() - startPosNDC.xy());

	float2 scaleDirection = float2{ -dirVecNDC.y, dirVecNDC.x } * widthScaleNDC;

	///Adapt to Screen Ratio
	scaleDirection.x = scaleDirection.x / RScene::Get().GetCamera().ScreenRatio;

	std::array<float3, 4> vertsT;

	//Clockwise, starting from left down
	vertsT[0] = startPosNDC.xyz() - float3{ scaleDirection, 0 };
	vertsT[1] = endPosNDC.xyz() - float3{ scaleDirection, 0 };
	vertsT[2] = endPosNDC.xyz() + float3{ scaleDirection, 0 };
	vertsT[3] = startPosNDC.xyz() + float3{ scaleDirection, 0 };

	///Transform back to WorldSpace
	auto& camera = RScene::Get().GetCamera();
	vertsT[0] = camera.GetObjectWorldSpacePosFromNDC(vertsT[0]);
	vertsT[1] = camera.GetObjectWorldSpacePosFromNDC(vertsT[1]);
	vertsT[2] = camera.GetObjectWorldSpacePosFromNDC(vertsT[2]);
	vertsT[3] = camera.GetObjectWorldSpacePosFromNDC(vertsT[3]);

	float3 startPosWorldReversed = camera.GetObjectWorldSpacePosFromNDC(startPosNDC.xyz());
	float3 endPosWorldReversed = camera.GetObjectWorldSpacePosFromNDC(endPosNDC.xyz());

	AddTriangleToBuffer(vertsT[0], vertsT[1], vertsT[2]);
	AddTriangleToBuffer(vertsT[0], vertsT[2], vertsT[3]);
}

void RGeometry::AddLineCameraFacingLong(c_float3& startPos, c_float3& endPos, float widthScale, float lineSegmentLengthThreshold)
{
	///We perform cam distance from endpoints, when distance between the endpoints themsleves is too big, the line isn't scaled properly
	float lineLength = VectorGetLength(endPos - startPos);
	float3 lineDirection = VectorNormalize(endPos - startPos);
	int numSegments = RMath::DivideAndRoundUp(lineLength, lineSegmentLengthThreshold);

	for (int i = 0; i < numSegments; i++)
	{
		float3 start = startPos + lineDirection * i * lineSegmentLengthThreshold;
		float3 end = startPos + lineDirection * (i + 1) * lineSegmentLengthThreshold;
		if (i == numSegments - 1)
		{
			end = endPos;
		}
		AddLine(start, end, widthScale);
	}

}

void RGeometry::AddTriangleCameraFacing(c_float3& baseCenterPos, c_float3& apexPos, float width)
{
	c_float3 dirVec = apexPos - baseCenterPos;

	c_float3 camPos = RScene::Get().GetCamera().GetPosition();
	c_float3 dirToCam = VectorNormalize(baseCenterPos - camPos);

	float3 lineNormal = VectorCross(dirVec, dirToCam);

	float widthStart = GetScaleOfObjectFacingCamera(baseCenterPos, width);

	lineNormal = VectorMul(VectorNormalize(lineNormal), widthStart);

	AddTriangleToBuffer(baseCenterPos + lineNormal, baseCenterPos - lineNormal, apexPos);

}

void RGeometry::AddLinePolyboardHermite(const float3& p0, const float3& p1, const float3& tangent0, const float3& tangent1, float r)
{
	//For each point
	float3 camPos = RScene::Get().GetCamera().GetPosition();

	auto getNormal = [&camPos, &r](const float3& p, const float3& v)
	{
		float3 dirToCam = VectorNormalize(camPos - p);
		float3 lineNormal = VectorCross(v, dirToCam);

		float scale = GetScaleOfObjectFacingCamera(p, 3.f * r);

		lineNormal = VectorMul((lineNormal), scale);
		return lineNormal;
	};

	float3 lineNormalp0 = getNormal(p0, tangent0);
	float3 lineNormalp1 = getNormal(p1, tangent1);

	std::array<float3, 4> vertsT;

	//Clockwise, starting from left down
	vertsT[0] = p0 - lineNormalp0;
	vertsT[1] = p1 - lineNormalp1;
	vertsT[2] = p1 + lineNormalp1;
	vertsT[3] = p0 + lineNormalp0;

	AddTriangleToBuffer(vertsT[0], vertsT[1], vertsT[2]);
	AddTriangleToBuffer(vertsT[0], vertsT[2], vertsT[3]);
}

void RGeometry::AddLinePolyboard(c_float3& pPrev, c_float3& p0, c_float3& p1, c_float3& pNext, float r)
{
	AddLinePolyboardHermite(p0, p1, VectorNormalize(p1 - pPrev), VectorNormalize(pNext - p0), r);
}

void RGeometry::AddCurvePolyboard(const RDynamicVector<float3>& samples, float lineThickness, const RDynamicVector<float3>& colors)
{
	check(samples.size() > 3);

	//For First Point
	if (!colors.empty())
	{
		SetColorRaw(colors[0]);
	}
	AddLinePolyboard(samples[0] + samples[0] - samples[1], samples[0], samples[1], samples[2], lineThickness);

	//For Inter Points
	for (int i = 1; i < samples.size() - 2; i++)
	{
		if (!colors.empty())
		{
			SetColorRaw(colors[i]);
		}
		AddLinePolyboard(samples[i - 1], samples[i], samples[i + 1], samples[i + 2], lineThickness);
	}

	//For EndPoint
	if (!colors.empty())
	{
		SetColorRaw(colors.back());
	}
	int n = samples.size();
	AddLinePolyboard(samples[n - 3], samples[n - 2], samples[n - 1], samples[n - 1] + samples[n - 1] - samples[n - 2], lineThickness);
}

void RGeometry::AddArrowPolyboard(c_float3& startPos, c_float3& endPos)
{
	///Generate Board line with Board Triangle
	//Pointer triangle takes 20% of the whole line
	float t = 0.8f;
	float3 LineEndPos = startPos * (1.f - t) + endPos * t;

	float thickness = 2.f;

	AddLine(startPos, LineEndPos, thickness);

	AddTriangleCameraFacing(LineEndPos, endPos, thickness * GoldenRatioX2);
}

void RGeometry::AddCubeLineList(c_float3& center, c_float3& extent, float lineThickness)
{
	RArray<float3, 4> topVertices;
	RArray<float3, 4> bottomVertices;

	float3 vTop;
	float3 vBottom;

	float3 TopCenter{ center + float3{ 0, extent.y, 0 } };
	float3 BottomCenter{ center - float3{ 0, extent.y, 0 } };

	//Starting from {-0.5, 0.5, -0.5} vertex
	//1
	vTop = float3{ TopCenter.x - extent.x, TopCenter.y, TopCenter.z - extent.z };
	vBottom = float3{ BottomCenter.x - extent.x, BottomCenter.y, BottomCenter.z - extent.z };
	topVertices[0] = (vTop);
	bottomVertices[0] = (vBottom);
	//2
	vTop = float3{ TopCenter.x + extent.x, TopCenter.y, TopCenter.z - extent.z };
	vBottom = float3{ BottomCenter.x + extent.x, BottomCenter.y, BottomCenter.z - extent.z };
	topVertices[1] = (vTop);
	bottomVertices[1] = (vBottom);
	//3
	vTop = float3{ TopCenter.x + extent.x, TopCenter.y,  TopCenter.z + extent.z };
	vBottom = float3{ BottomCenter.x + extent.x, BottomCenter.y, BottomCenter.z + extent.z };
	topVertices[2] = (vTop);
	bottomVertices[2] = (vBottom);
	//4
	vTop = float3{ TopCenter.x - extent.x, TopCenter.y,  TopCenter.z + extent.z };
	vBottom = float3{ BottomCenter.x - extent.x, BottomCenter.y, BottomCenter.z + extent.z };
	topVertices[3] = (vTop);
	bottomVertices[3] = (vBottom);


	//Create 12 Lines for 6 planes

	//Top plane
	{
		AddLine(topVertices[0], topVertices[1], lineThickness);
		AddLine(topVertices[0], topVertices[3], lineThickness);
		AddLine(topVertices[2], topVertices[1], lineThickness);
		AddLine(topVertices[2], topVertices[3], lineThickness);
	}

	//Sides
	/*for(int v = 0; v < 4; v++)
	{
	geometryBufferCPU.push_back(topVertices[v]);
	geometryBufferCPU.push_back(bottomVertices[v]);
	}*/

	for (int v = 0; v < 4; v++)
	{
		AddLine(topVertices[v], VectorDiv(topVertices[v] + bottomVertices[v], 2.0f), lineThickness);
	}

	for (int v = 0; v < 4; v++)
	{
		AddLine(VectorDiv(topVertices[v] + bottomVertices[v], 2.0f), bottomVertices[v], lineThickness);
	}

	//Bottom plane
	{
		AddLine(bottomVertices[0], bottomVertices[1], lineThickness);
		AddLine(bottomVertices[0], bottomVertices[3], lineThickness);
		AddLine(bottomVertices[2], bottomVertices[1], lineThickness);
		AddLine(bottomVertices[2], bottomVertices[3], lineThickness);
	}
}



void RGeometry::TransformAndReturn(const Matrix4x4& TransformAndReturn)
{
	/*for (auto& vec : VertexBuffer)
	{
		vec = VectorMul(float4{ vec, 1.f }, Transform).xyz();
	}

	for (auto& n : Normals)
	{
		n = VectorMul(n, Transform);
	}
	for (auto& t : Tangents)
	{
		t = VectorMul(t, Transform);
	}*/
}

void RGeometry::TranslateVertices(const float3& newPos)
{
	///*Matrix4x4 translationMat = MatrixTranslate(newPos);

	//for (auto& vec : Vertices)
	//{
	//	vec = VectorMul(float4{ vec }, translationMat);
	//}*/

	//for (auto& vec : VertexBuffer)
	//{
	//	vec = vec + newPos;
	//}


}

void RViewFrustum::GenerateForLineList()
{
	RDynamicVector<float3> verticesFarPlane;
	RDynamicVector<float3> verticesNearPlane;

	float3 Origin{ 0,0,0 };

	float tanY = tan(FoV / 2.0);
	float tanX = tanY * Ratio;

	float3 FrustumCornerVertex;

	FrustumCornerVertex.z = 1;
	FrustumCornerVertex.y = tanY;
	FrustumCornerVertex.x = tanX;

	HalfWidth = FrustumCornerVertex.x;
	HalfHeight = FrustumCornerVertex.y;

	float FoVY = 2 * atan(HalfHeight * NearPlane / NearPlane);
	float FoVX = 2 * atan(HalfWidth * NearPlane / NearPlane);

	//Far Plane Vertices
	float3 FarPlaneV = VectorMul(FrustumCornerVertex, FarPlane);
	//+X+Y
	verticesFarPlane.push_back(FarPlaneV);
	//-X+Y
	verticesFarPlane.push_back(float3{ -FarPlaneV.x, FarPlaneV.y, FarPlaneV.z });
	//-X-Y
	verticesFarPlane.push_back(float3{ -FarPlaneV.x, -FarPlaneV.y, FarPlaneV.z });
	//+X-Y
	verticesFarPlane.push_back(float3{ FarPlaneV.x, -FarPlaneV.y, FarPlaneV.z });

	//Near Plane Vertices
	float3 NearPlaneV = VectorMul(FrustumCornerVertex, NearPlane);
	//+X+Y
	verticesNearPlane.push_back(NearPlaneV);
	//-X+Y
	verticesNearPlane.push_back(float3{ -NearPlaneV.x, NearPlaneV.y, NearPlaneV.z });
	//-X-Y
	verticesNearPlane.push_back(float3{ -NearPlaneV.x, -NearPlaneV.y, NearPlaneV.z });
	//+X-Y
	verticesNearPlane.push_back(float3{ NearPlaneV.x, -NearPlaneV.y, NearPlaneV.z });


	//Line List
	//Top plane
	{
		AddLine(verticesFarPlane[0], verticesFarPlane[1]);
		AddLine(verticesFarPlane[0], verticesFarPlane[3]);
		AddLine(verticesFarPlane[2], verticesFarPlane[1]);
		AddLine(verticesFarPlane[2], verticesFarPlane[3]);
	}
	//Sides
	for (int v = 0; v < 4; v++)
	{
		AddLine(verticesFarPlane[v], verticesFarPlane[v]);
	}
	//Bottom plane
	{
		AddLine(verticesNearPlane[0], verticesNearPlane[1]);
		AddLine(verticesNearPlane[0], verticesNearPlane[3]);
		AddLine(verticesNearPlane[2], verticesNearPlane[1]);
		AddLine(verticesNearPlane[2], verticesNearPlane[3]);
	}
	//Origin
	for (int v = 0; v < 4; v++)
	{
		AddLine(Origin, verticesNearPlane[v]);
	}
}

void RCube::GenerateForTriangleList()
{

	RDynamicVector<float3> geometryBufferCPU;

	RDynamicVector<float3> topVertices;
	RDynamicVector<float3> bottomVertices;

	float half = 0.5f * Size;
	float3 vTop;
	float3 vBottom;

	//Starting from {-0.5, 0.5, -0.5} vertex
	//1
	vTop = float3{ -half, half, -half };
	vBottom = float3{ -half, -half, -half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//2
	vTop = float3{ half, half, -half };
	vBottom = float3{ half, -half, -half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//3
	vTop = float3{ half, half, half };
	vBottom = float3{ half, -half, half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//4
	vTop = float3{ -half, half, half };
	vBottom = float3{ -half, -half, half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);


	//Create 12 Triangles for 6 planes

	//Top plane
	{
		AddQuadToBuffer(topVertices[0], topVertices[1], topVertices[2], topVertices[3]);
	}

	//Sides
	//1 triangle per vertex
	for (int t = 0; t < 4; t++)
	{
		AddTriangleToBuffer(topVertices[t], bottomVertices[t], topVertices[(t + 1) % 4]);
	}

	for (int t = 0; t < 4; t++)
	{
		AddTriangleToBuffer(bottomVertices[t], bottomVertices[(t + 1) % 4], topVertices[(t + 1) % 4]);
	}

	//Bottom plane
	{
		AddQuadToBuffer(bottomVertices[0], bottomVertices[1], bottomVertices[2], bottomVertices[3]);

	}

}

void RCube::GenerateForLineList()
{

	RDynamicVector<float3> geometryBufferCPU;

	RDynamicVector<float3> topVertices;
	RDynamicVector<float3> bottomVertices;

	float half = 0.5f * Size;
	float3 vTop;
	float3 vBottom;

	//Starting from {-0.5, 0.5, -0.5} vertex
	//1
	vTop = float3{ -half, half, -half };
	vBottom = float3{ -half, -half, -half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//2
	vTop = float3{ half, half, -half };
	vBottom = float3{ half, -half, -half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//3
	vTop = float3{ half, half, half };
	vBottom = float3{ half, -half, half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);
	//4
	vTop = float3{ -half, half, half };
	vBottom = float3{ -half, -half, half };
	topVertices.push_back(vTop);
	bottomVertices.push_back(vBottom);


	//Create 12 Lines for 6 planes

	//Top plane
	{
		AddLine(topVertices[0], topVertices[1]);
		AddLine(topVertices[0], topVertices[3]);
		AddLine(topVertices[2], topVertices[1]);
		AddLine(topVertices[2], topVertices[3]);
	}
	//Sides

	/*for(int v = 0; v < 4; v++)
	{
	geometryBufferCPU.push_back(topVertices[v]);
	geometryBufferCPU.push_back(bottomVertices[v]);
	}*/

	for (int v = 0; v < 4; v++)
	{
		AddLine(topVertices[v], VectorDiv(topVertices[v] + bottomVertices[v], 2.0f));
	}

	for (int v = 0; v < 4; v++)
	{
		AddLine(VectorDiv(topVertices[v] + bottomVertices[v], 2.0f), bottomVertices[v]);
	}

	//Bottom plane
	{
		AddLine(bottomVertices[0], bottomVertices[1]);
		AddLine(bottomVertices[0], bottomVertices[3]);
		AddLine(bottomVertices[2], bottomVertices[1]);
		AddLine(bottomVertices[2], bottomVertices[3]);
	}

}

#pragma once

#include "LinearAlgebra.h"
#include "DrawSimpleGeometry.h"

using AABB = DirectX::BoundingBox;
using BoundingSphere = DirectX::BoundingSphere;


static inline float3 GetLinePointDistanceVector(const float3& LinePoint1, const float3& LinePoint2, const float3& Point)
{
	float3 PointVector = Point - LinePoint1;
	float3 LineVector = LinePoint2 - LinePoint1;

	float LengthSq = VectorGetLengthSquared(LineVector);

	float PointProjectionScale = VectorDot(PointVector, LineVector);
	PointProjectionScale = PointProjectionScale / LengthSq;

	float3 DistanceVector = VectorMul(LineVector, PointProjectionScale);
	DistanceVector = PointVector - DistanceVector;

	return (DistanceVector);
}

static inline float GetLinePointDistance(const float3& LinePoint1, const float3& LinePoint2, const float3& Point)
{
	return VectorGetLength(GetLinePointDistanceVector(LinePoint1, LinePoint2, Point));
}


static inline float GetClosestDistanceFromPointToLineSegment(c_float3& p, c_float3& p0, c_float3& p1, float3& closestPoint)
{
	float3 d = p1 - p0;
	// Find the parameter t that minimizes the distance from p to the edge
	float t = VectorDot(p - p0, d) / VectorDot(d, d);

	// Clamp t to [0, 1] to get the closest point on the edge
	t = std::max(0.0f, std::min(1.0f, t));

	// Find the distance from p to the closest point on the edge
	closestPoint = (p0 + d * t);
	return VectorGetLength(p - closestPoint);
}

static inline bool RaySphereIntersection(const float3& rayOrigin, const float3& rayDir, const float3& spherePos, float sphereRadius)
{
	float3 m = rayOrigin - spherePos;
	float b = VectorDot(m, rayDir);
	float c = VectorDot(m, m) - sphereRadius * sphereRadius;
	if (c > 0.0f && b > 0.0f)
		return false;
	float discr = b * b - c;
	if (discr < 0.0f)
		return false;
	return true;
}

static inline bool RayAABBIntersection(const float3& rayOrigin, const float3& rayDirection, const float3& aabbCenter, const float3& aabbExtent, float& t)
{
	// Calculate the inverse direction of the ray
	float3 invDirection;
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
	float tEnter = std::max(std::max(std::min(tmin, tmax), std::min(tymin, tymax)), std::min(tzmin, tzmax));
	float tExit = std::min(std::min(std::max(tmin, tmax), std::max(tymin, tymax)), std::max(tzmin, tzmax));

	// Check if the intersection is outside the valid range
	if (tExit < 0 || tEnter > tExit)
	{
		return false;
	}

	// Set the intersection distance
	t = (tEnter >= 0) ? tEnter : tExit;

	return true;
}

static inline bool LineSphereIntersection(const float3& rayOrigin, const float3& rayEnd, const float3& spherePos, float sphereRadius)
{
	return GetLinePointDistance(rayOrigin, rayEnd, spherePos) < sphereRadius;
}



static inline bool LineSegmentSphereIntersection(const float3& lineStart, const float3& lineEnd, const float3& spherePos, float sphereRadius)
{
	float3 dummy;
	return GetClosestDistanceFromPointToLineSegment(spherePos, lineStart, lineEnd, dummy) < sphereRadius;
}

static inline bool RayPlaneIntersection(float4 plane, float3 p1, float3 p2, float3& result, float& tResult)
{
	float3 RayOrigin = p1; //point on a line
	float d = plane.w;
	float3 n{ plane.x, plane.y, plane.z };

	float3 rayDir = VectorNormalize(p2 - RayOrigin);

	//p = S + tRay

	if (VectorDot(rayDir, n) == 0) //Line is parallel to plane
	{
		return false;
	}

	// Compute the t value for the directed line ray intersecting the plane
	float t = (d - VectorDot(n, RayOrigin)) / VectorDot(n, rayDir);

	result = RayOrigin + VectorMul(rayDir, t);
	tResult = t;

	if (t > 0)
		return true;
	else return false;

}

static inline bool RectangleSphereIntersection(const float2& rectCenter, const float2& extent, const float2& sphereCenter, float sphereRadius) {
	// Calculate the minimum and maximum values of the square's x and y coordinates
	float rectMinX = rectCenter.x - extent.x;
	float rectMaxX = rectCenter.x + extent.x;
	float rectMinY = rectCenter.y - extent.y;
	float rectMaxY = rectCenter.y + extent.y;

	// Find the closest point on the square to the center of the sphere
	float closestX = std::max(rectMinX, std::min(sphereCenter.x, rectMaxX));
	float closestY = std::max(rectMinY, std::min(sphereCenter.y, rectMaxY));

	// Calculate the distance between the closest point and the sphere's center
	float distanceX = closestX - sphereCenter.x;
	float distanceY = closestY - sphereCenter.y;
	float distanceSquared = distanceX * distanceX + distanceY * distanceY;

	// Compare the distance squared with the square of the sphere's radius
	float radiusSquared = sphereRadius * sphereRadius;
	return distanceSquared <= radiusSquared;
}

static inline float GetPlanePointDistance(float3 n, float d, float3 point)
{
	return VectorDot(n, point) + d;
}

//If absolute distance between closest point on a plane and sphere center is < radius --> intersection
static inline bool SpherePlaneIntersection(float4 plane, float3 spherePos, float sphereRadius)
{
	float d = plane.w;
	float3 n{ plane.x, plane.y, plane.z };

	//P* N(normalized)) + D < r

	float distance = VectorDot(n, spherePos) + d;

	return fabs(distance) <= sphereRadius;
}

static inline bool SphereSphereIntersection(float3 spherePos, float sphereRadius, float3 spherePos2, float sphereRadius2)
{
	// Distance  between centers.
	float d = VectorGetLength(spherePos - spherePos2);

	// Sum of the rad.
	float r = sphereRadius + sphereRadius2;

	return d <= r;
}

static inline float3 FindCentroid(std::vector<float3>& vertices)
{
	float3 sum{ 0,0,0 };
	for (auto& v : vertices)
	{
		sum = sum + v;
	}
	return VectorDiv(sum, vertices.size());
}

///N * P + D = 0
class BoundingPlane
{
public:

	float3 n;
	float d = 0.f;

	BoundingPlane() = default;

	BoundingPlane(c_float3& normal, float _d)
	{
		n = VectorNormalize(normal);
		d = _d;
	}

	BoundingPlane(c_float3& pointOnPlane, c_float3& normal)
	{
		Initialize(pointOnPlane, normal);
	}

	void Initialize(c_float3& pointOnPlane, c_float3& normal)
	{
		n = VectorNormalize(normal);
		d = -VectorDot(pointOnPlane, n);
	}

	BoundingPlane(c_float3& p0, c_float3& p1, c_float3& p2)
	{
		Initialize(p0, VectorCross(p1 - p0, p2 - p0));
	}

	// Returns the direction the plane is facing.  (Warning:  might not be normalized.)
	float3 GetNormal(void) { return n; }

	// Returns the point on the plane closest to the origin
	float3 GetPointOnPlaneClosestToOrigin(void) { return VectorNegate(n) * d; }

	float3 GetClosestPointOnAPlaneTo(c_float3& point) {
		// Calculate the signed distance from the given point to the plane
		// Subtract the normal vector scaled by the distance from the given point to get the closest point on the plane
		return point - n * GetSignedDistanceFromPoint(point);
	}

	///if dist == 0 --> the point lies in the plane
	///if dist > 0 --> the point lies on the positive side of the plane since point would be on the side in which the normal vector points
	///if dist < 0 the point lies on the negative side of the plane.
	float GetSignedDistanceFromPoint(c_float3& point)
	{
		return VectorDot(n, point) + d;
	}

	float GetSignedDistanceFromPointHomogeneous(float4 point)
	{
		float4 L{ n, d };
		return VectorDot(L, point);
	}

	void Translate(c_float3& offset)
	{
		///Get point on plane
		auto point = GetPointOnPlaneClosestToOrigin();

		///Tranlate the point
		point = point + offset;

		///Recreate the plane
		Initialize(point, n);
	}

	enum class EIntersection
	{
		Inside,
		Outside,
		Intersects
	};

	EIntersection IntersectsSphere(c_float3& center, float radius, float& absoluteDistance)
	{
		auto signedDist = GetSignedDistanceFromPoint(center);

		absoluteDistance = abs(signedDist);

		if (absoluteDistance <= radius)
		{
			return EIntersection::Intersects;
		}
		if (signedDist > 0)
		{
			return EIntersection::Inside;
		}
		return EIntersection::Outside;

	}


};

/// A function that checks if two lines intersect and returns true if they do
/// The point of intersection is passed by reference and is modified if there is an intersection
namespace RIntersectionCheck
{
	struct Edge
	{
		c_float3& p1, p2;
	};
}
static inline bool EdgeEdgeIntersection(RIntersectionCheck::Edge l1, RIntersectionCheck::Edge l2, float3& outIntersectionPoint1, float3& outIntersectionPoint2)
{
	/// Find the direction vectors of the lines
	float3 d1 = l1.p2 - l1.p1;
	float3 d2 = l2.p2 - l2.p1;

	/// Find the normal vector of the plane containing the projected lines
	float3 d1Xd2 = VectorCross(d1, d2);


	bool bLinesAreParallel = false;

	///TODO:Most likely this branch will never be executed and coincident lines will be skipped, imo this is more failsafe
	/// Check if the lines are parallel or coincident
	if (d1Xd2.x == 0 && d1Xd2.y == 0 && d1Xd2.z == 0)
	{
		/// Check if the lines are coincident by comparing one point of each line
		auto c = VectorCross(l1.p1 - l2.p1, d2);
		if (c.x == 0 && c.y == 0 && c.z == 0)
		{
			/// The lines are coincident, check if any point intersects edge
			float epsilon = 0.01;

			outIntersectionPoint1 = l1.p1;
			float curD = GetClosestDistanceFromPointToLineSegment(outIntersectionPoint1, l2.p1, l2.p2, outIntersectionPoint2);
			if (curD < epsilon && curD != 0.f)
			{
				return true;
			}
			outIntersectionPoint1 = l1.p2;
			curD = GetClosestDistanceFromPointToLineSegment(outIntersectionPoint1, l2.p1, l2.p2, outIntersectionPoint2);
			if (curD < epsilon && curD != 0.f)
			{
				return true;
			}
			outIntersectionPoint1 = l2.p1;
			curD = GetClosestDistanceFromPointToLineSegment(outIntersectionPoint1, l1.p1, l1.p2, outIntersectionPoint2);
			if (curD < epsilon && curD != 0.f)
			{
				return true;
			}
			outIntersectionPoint1 = l2.p2;
			curD = GetClosestDistanceFromPointToLineSegment(outIntersectionPoint1, l1.p1, l1.p2, outIntersectionPoint2);
			if (curD < epsilon && curD != 0.f)
			{
				return true;
			}
		}
	}

	c_float3& p1 = l1.p1;
	c_float3& p2 = l2.p1;


	float d1Xd2LengthSq = VectorGetLengthSquared(d1Xd2);
	float t1 = VectorDot(VectorCross(p2 - p1, d2), d1Xd2) / d1Xd2LengthSq;
	float t2 = VectorDot(VectorCross(p2 - p1, d1), d1Xd2) / d1Xd2LengthSq;

	outIntersectionPoint1 = p1 + d1 * t1;
	outIntersectionPoint2 = p2 + d2 * t2;

	//float2 minmaxT{ 0, 1 };
	float2 minmaxT{ 0.2, 0.8 };

	if ((t1 > minmaxT.x) && (t1 < minmaxT.y) && (t2 > minmaxT.x) && (t2 < minmaxT.y))///Don't include the very ends if the lines
	{
		return true;
	}

	return false;

}

static inline bool SphereTriangleIntersectionBarycentric(c_float3& spherePos, float sphereRadius, c_float3& p0, c_float3& p1, c_float3& p2, float3& closestPoint, float3& barycentric)
{
	// Find the normal vector of the triangle
	float3 n = VectorCross(p1 - p0, p2 - p0);
	n = VectorNormalize(n);

	// Find the distance from the sphere's center to the triangle's plane
	float d = VectorDot(spherePos - p0, n);

	// If the distance is larger than the sphere's radius, they do not intersect
	if (abs(d) > sphereRadius)
	{
		return false;
	}

	// Project the sphere's center onto the triangle's plane
	float3 p = spherePos - n * d;

	// Convert the projected point to barycentric coordinates using the function above
	float3 b = CartesianToBarycentric(p, p0, p1, p2);

	// Check if the barycentric coordinates are in [0, 1] and sum to 1
	if (b.x >= 0 && b.y >= 0 && b.z >= 0)
	{
		// The closest point on the triangle is the projected point
		closestPoint = p;
		barycentric = b;
		return true;
	}

	return false;
}


static inline bool SphereTriangleIntersection(c_float3& spherePos, float sphereRadius, c_float3& p0, c_float3& p1, c_float3& p2, float3& closestPoint)
{
	// Find the normal vector of the triangle
	float3 n = VectorCross(p1 - p0, p2 - p0);
	n = VectorNormalize(n);

	// Find the distance from the sphere's center to the triangle's plane
	float d = VectorDot(spherePos - p0, n);

	// If the distance is larger than the sphere's radius, they do not intersect
	if (abs(d) > sphereRadius)
	{
		return false;
	}

	// Project the sphere's center onto the triangle's plane
	float3 p = spherePos - n * d;

	// Convert the projected point to barycentric coordinates using the function above
	float3 b = CartesianToBarycentric(p, p0, p1, p2);

	// Check if the barycentric coordinates are in [0, 1] and sum to 1
	if (b.x >= 0 && b.y >= 0 && b.z >= 0)
	{
		// The closest point on the triangle is the projected point
		closestPoint = p;
		return true;
	}

	// Otherwise, find the closest point on each edge of the triangle and check if it is within the sphere's radius
	float3 edge[3] = { p1 - p0, p2 - p1, p0 - p2 };
	float3 point[3] = { p0, p1, p2 };

	// Initialize the minimum distance and index variables
	float minDist = sphereRadius;
	float3 potentiallyClosestPoint;
	bool bCollision = false;

	for (int i = 0; i < 3; i++)
	{
		// Find the parameter t that minimizes the distance from p to the edge
		float t = VectorDot(p - point[i], edge[i]) / VectorDot(edge[i], edge[i]);

		// Clamp t to [0, 1] to get the closest point on the edge
		t = std::max(0.0f, std::min(1.0f, t));

		// Find the distance from p to the closest point on the edge
		float3 pointCur = (point[i] + edge[i] * t);
		float dist = VectorGetLength(p - pointCur);

		if (dist < minDist)
		{
			minDist = dist;
			potentiallyClosestPoint = pointCur;
			bCollision = true;
		}

	}

	if (bCollision)
	{
		// Otherwise, they intersect and we write the closest point on the triangle to the ref parameter
		closestPoint = potentiallyClosestPoint;
		return true;
	}

	return false;
	
}


///Cut out anything outside of the Sphere defined by CamPos and ViewDistance
static inline bool TruncateLineByCameraSphereCulling(float3& p0, float3& p1, c_float3& camPos, float viewDistance = 100.f)
{
	float distToCam0 = VectorGetLength(p0 - camPos);
	float distToCam1 = VectorGetLength(p1 - camPos);

	///Use cases:

	if (distToCam0 <= viewDistance && distToCam1 <= viewDistance)
	{
		///Points within sphere
		return true;
	}

	else if (distToCam0 > viewDistance && distToCam1 > viewDistance)
	{
		///Line passes through sphere?
		if (RaySphereIntersection(p0, VectorNormalize(p1 - p0), camPos, viewDistance))
		{
			DirectX::BoundingSphere sphere;
			sphere.Center = camPos;
			sphere.Radius = viewDistance;

			///prepare rays
			float3SSE p0XM = p0;
			float3SSE p0Dir = VectorNormalize(p1 - p0);
			float3SSE p1XM = p1;
			float3SSE p1Dir = VectorNormalize(p0 - p1);
			float newDist = 0;
			///Test for p0
			check(sphere.Intersects(p0XM, p0Dir, newDist));
			p0 = p0 + float3(p0Dir) * newDist;
			///Test for p1
			check(sphere.Intersects(p1XM, p1Dir, newDist));
			p1 = p1 + float3(p1Dir) * newDist;
			return true;
		}
		else
		{
			///All points outside
			return false; ///discard the line
		}
	}
	///2. One Point outside
	else
	{
		DirectX::BoundingSphere sphere;
		sphere.Center = camPos;
		sphere.Radius = viewDistance;

		float newDist = 0;
		if (distToCam0 > viewDistance)
		{
			///point0 is outside
			float3SSE p0XM = p0;
			float3SSE p0Dir = VectorNormalize(p1 - p0);
			check(sphere.Intersects(p0XM, p0Dir, newDist));
			p0 = p0 + float3(p0Dir) * newDist;

		}
		else if (distToCam1 > viewDistance)
		{
			///point1 is outside
			float3SSE p1XM = p1;
			float3SSE p1Dir = VectorNormalize(p0 - p1);
			check(sphere.Intersects(p1XM, p1Dir, newDist));
			p1 = p1 + float3(p1Dir) * newDist;
		}
		else
		{
			check(false);
		}
		return true;
	}

	check(false);
	return false;

}






























static inline bool ClipTriangleAgainstPlane(float4 Plane, RTriangleSimple& triangle, std::vector<RTriangleSimple>& returnTriangles)
{
	float d = Plane.w;
	float3 n = float3{ Plane.x, Plane.y, Plane.z };
	float3 v1 = triangle.V1;
	float3 v2 = triangle.V2;
	float3 v3 = triangle.V3;

	//Detect what vertices are outside the plane
	//d = (P * N(normalized)) + D
	//d < 0 means outside the plane
	bool v1Out = (VectorDot(v1, n) - d) < 0;
	bool v2Out = (VectorDot(v2, n) - d) < 0;
	bool v3Out = (VectorDot(v3, n) - d) < 0;

	if (v1Out && v2Out && v3Out) //Triangle is outside and discarded
	{
		return false;
	}
	else if (!v1Out && !v2Out && !v3Out) //Triangle is fully inside and should be rendered as is
	{
		return false;
	}

	//What edges need to be tested
	//If 1 vertex is outside - edge should be cut
	//if 2 vertices are outside - edge is discarded
	//if 2 vertices are inside - edge is kept as is

	struct Edge
	{
		float3 rayOrigin;
		float3 rayDir;
	};

	std::vector<Edge> edgesToCut;

	Edge e1{ v1, VectorNormalize(v2 - v1) };
	Edge e12{ v2, VectorNormalize(v1 - v2) };

	Edge e2{ v2, VectorNormalize(v3 - v2) };
	Edge e22{ v3, VectorNormalize(v2 - v3) };

	Edge e3{ v3, VectorNormalize(v1 - v3) };
	Edge e32{ v1, VectorNormalize(v3 - v1) };


	std::vector<float3> resultVertices;


	if (!v1Out)
	{
		resultVertices.push_back(v1);
	}
	if (!v2Out)
	{
		resultVertices.push_back(v2);
	}
	if (!v3Out)
	{
		resultVertices.push_back(v3);
	}



	if (v1Out && !v2Out || !v1Out && v2Out)//cut e1
	{
		if (v1Out)
		{
			edgesToCut.push_back(e12);
		}
		else
		{
			edgesToCut.push_back(e1);
		}

	}
	if (v3Out && !v2Out || !v3Out && v2Out)//cut e2
	{
		if (v2Out)
		{
			edgesToCut.push_back(e22);
		}
		else
		{
			edgesToCut.push_back(e2);
		}

	}
	if (v1Out && !v3Out || !v1Out && v3Out)//cut e3
	{
		if (v3Out)
		{
			edgesToCut.push_back(e32);
		}
		else
		{
			edgesToCut.push_back(e3);
		}
	}
	check(edgesToCut.size() < 3);



	for (Edge& e : edgesToCut)
	{
		//t = (N dot EdgeStart) / N dot EdgeDir
		//new vertex = EdgeStart + abs(t) * EdgeDir

		float t = (d - VectorDot(n, e.rayOrigin)) / VectorDot(n, e.rayDir);
		float3 newVertex = e.rayOrigin + VectorMul(e.rayDir, abs(t));
		resultVertices.push_back(newVertex);
	}

	check(resultVertices.size() >= 3);

	const int numTriangles = resultVertices.size() - 2;

	struct FastPredicate {
		float3 Center;
		bool operator() (float3 i, float3 j)
		{
			//float angle1 = atan(double(i.y - Center.y) / double(i.x - Center.x));
			//float angle2 = atan(double(j.y - Center.y) / double(j.x - Center.x));
			float3 vec1 = VectorNormalize(i - Center);
			float3 vec2 = VectorNormalize(j - Center);
			/*float angle1 = atan2(vec1.y, vec1.x);
			float angle2 = atan2(vec2.y, vec2.x);*/
			float angle1 = atan2(vec1.z, vec1.x);
			float angle2 = atan2(vec2.z, vec2.x);
			//angle1 = angle1 + asin(-vec1.z);
			//angle2 = angle2 + asin(-vec2.z);
			return angle1 < angle2;
		}
	} mypredicate;


	//////////////////////////////////////////////////////////////////////////////////////////////////////

	{
		float3 CenterPos{ 0,0,0 };

		std::vector<float3> vertsToSort;
		vertsToSort.push_back(float3{ -1, 0, 0 });
		vertsToSort.push_back(float3{ 0, 1, 0 });
		vertsToSort.push_back(float3{ 0, -1, 0 });
		vertsToSort.push_back(float3{ 1, 0, 0 });

		CenterPos = FindCentroid(vertsToSort);
		mypredicate.Center = CenterPos;
		//mypredicate.ComparisonVector = VectorNormalize(vertsToSort[0] - CenterPos);

		std::sort(vertsToSort.begin(), vertsToSort.end(), mypredicate);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////


	//Sort Vertices
	float3 TriCenterPos = FindCentroid(resultVertices);
	mypredicate.Center = TriCenterPos;
	//mypredicate.ComparisonVector = VectorNormalize(resultVertices[0] - TriCenterPos);


	std::sort(resultVertices.begin(), resultVertices.end(), mypredicate);

	for (int v = 1; v <= numTriangles; v++)
	{
		RTriangleSimple tri;
		tri.V1 = resultVertices[0];
		tri.V2 = resultVertices[v];
		tri.V3 = resultVertices[v + 1];
		returnTriangles.push_back(tri);
	}


	return true;


}




///**************************************************************************************************************
///
///											ACCELERATION STRUCTURES
/// 
///**************************************************************************************************************


struct Particle
{
	int GlobalIndexToNode = -1;
	BoundingSphere Sphere;
};


///For dynamic objects, fast re-insertion, fast neighbor search
class RDynamicOctree
{
public:

	float3 RootCenter;
	float3 RootExtent;

	void Initialize(c_float3& _rootCenter, c_float3& _rootExtent)
	{


		///Deduce the size of the tree
		uint treeSize = GetDepthStartIndex(MaxDepth + 1);
		Nodes.resize(treeSize);

		Node& root = Nodes[0];
		root.LocalCoord = uint3{ 0,0,0 };
		root.Depth = 0;
		root.IndexToChildrenStart = 1;
		Nodes.push_back(root);

		for (int i = 1; i <= MaxDepth; i++)
		{
			InitializeDepthLevel(i);
		}
	}

	void InitializeDepthLevel(uint curDepth)
	{
		uint depthStart = GetDepthStartIndex(curDepth);
		uint numNodes = GetNumNodesFromDepth(curDepth);
		uint childrenDepthStart = GetDepthStartIndex(curDepth + 1);
		for (int i = 0; i < numNodes; i++)
		{
			///Initialize nodes
			auto curNodeIndex = depthStart + i;
			auto& curNode = Nodes[curNodeIndex];
			curNode.Depth = curDepth;
			curNode.LocalCoord = GetLocalCoordFromIndex(curDepth, i);

			GenerateBoundsForNode(curNode);

			if (curDepth != (MaxDepth))
			{
				curNode.IndexToChildrenStart = childrenDepthStart + i * 8;
				for (int j = 0; j < 8; j++)
				{
					auto& child = Nodes[curNode.IndexToChildrenStart + j];
					child.GlobalIndexToParent = curNodeIndex;
				}
			}
		}
	}


	static constexpr inline uint MaxDepth = 8;

	struct Node
	{
		uint3 LocalCoord;
		uint Depth;
		uint IndexToChildrenStart = 0;
		uint GlobalIndexToParent;
		std::vector<Particle*> Objects;
		AABB Bounds;
	};

	void Insert(Particle& object)
	{
		if (object.GlobalIndexToNode == -1)
		{
			InsertInner(0, object);
		}
		else
		{
			ReInsert(object);
		}
		
	}

	///TODO: Object should be present even in cells that only intersect it
	///And we need to check the neighbours of all cells that this object is present in
	/*void FindNeighbours(Particle& obj, RDynamicVector<Particle*>& outNodes)
	{
		outNodes.clear();

		auto& Node = Nodes[obj.GlobalIndexToNode];
		
		uint3 startCoord;
		startCoord.x = std::max(0, int(Node.LocalCoord.x) - 1);
		startCoord.y = std::max(0, int(Node.LocalCoord.y) - 1);
		startCoord.z = std::max(0, int(Node.LocalCoord.z) - 1);

		for (int x = 0; x < 2; x++)
		{
			for (int y = 0; y < 2; y++)
			{
				for (int z = 0; z < 2; z++)
				{
					uint3 coord{ startCoord.x + x, startCoord.y + y, startCoord.z + z };
					uint nodeId = GetLocalIndexFromCoord(Node.Depth, coord);
					auto& otherObjects = Nodes[nodeId].Objects;
					for (auto& otherObj : otherObjects)
					{
						outNodes.push_back(otherObj);
					}
				}
			}
		}
	}*/

	void ReInsert(Particle& obj)
	{
		auto& curNode = Nodes[obj.GlobalIndexToNode];
		///remove particle from node
		auto it = std::find(curNode.Objects.begin(), curNode.Objects.end(), &obj);
		check(it != curNode.Objects.end());
		curNode.Objects.erase(it);

		///TODO: Test object against neighboring boxes first
		
		///test object against all boxes from 1 depth level higher
		int curDepth = int(curNode.Depth);
		uint depthToCheck = std::max(0, curDepth - 1);
		uint NumNodes = GetNumNodesFromDepth(depthToCheck);
		uint startOffset = GetDepthStartIndex(depthToCheck);
		bool bEmplaced = false;
		for (int i = 0; i < NumNodes; i++)
		{
			if (InsertInner(startOffset + i, obj))
			{
				return;
			}
		}
		check(false);

	}

	bool InsertInner(uint curNodeIndex, Particle& obj)
	{
		Node& node = Nodes[curNodeIndex];
		auto containmentType = node.Bounds.Contains(obj.Sphere);
		
		bool bObjectEmplacedInChild = false;
		if (containmentType == DirectX::ContainmentType::CONTAINS)///node contains the object but what about children
		{
			for (int i = 0; i < 8; i++)
			{
				uint offset = node.IndexToChildrenStart;
				bObjectEmplacedInChild |= InsertInner(offset + i, obj);
				if (bObjectEmplacedInChild)
				{
					return true;
				}
			}
			if (!bObjectEmplacedInChild)
			{
				node.Objects.emplace_back(&obj);
				obj.GlobalIndexToNode = curNodeIndex;
				return true;
			}
		}
		else if (containmentType == DirectX::ContainmentType::INTERSECTS)
		{
			check(node.Depth != 0);
			return false;
		}
		else
		{
			return false;
		}
	}


	float3 GetExtentOfDepthLevel(uint depth)
	{
		return RootExtent / pow(0.5, depth);
	}

	float3 GetNodePositionFromLocalCoord(uint3 coord, c_float3& nodeExtents, int curGridLength)
	{
		float3 pos;
		pos.x = (int(coord.x) - curGridLength / 2) * nodeExtents.x * 2 + nodeExtents.x;
		pos.y = (int(coord.y) - curGridLength / 2) * nodeExtents.y * 2 + nodeExtents.y;
		pos.z = (int(coord.z) - curGridLength / 2) * nodeExtents.z * 2 + nodeExtents.z;
		pos = pos + RootCenter;
		return pos;
	}

	void GenerateBoundsForNode(Node& node)
	{
		node.Bounds.Extents = GetExtentOfDepthLevel(node.Depth);
		node.Bounds.Center = GetNodePositionFromLocalCoord(node.LocalCoord, *(const float3*)(&node.Bounds.Extents), GetSideLengthFromDepth(node.Depth));
	}

	uint GetSideLengthFromDepth(uint depth)
	{
		return ipow(2, depth);
	}
	uint GetNumNodesFromDepth(uint depth)
	{
		return ipow(8, depth);
	}

	uint GetLocalIndexFromCoord(uint depth, uint3 coord)
	{
		uint size = GetSideLengthFromDepth(depth);
		return (coord.y + size * coord.z) * size + coord.x;
	}
	uint3 GetLocalCoordFromIndex(uint depth, uint index)
	{
		uint size = GetSideLengthFromDepth(depth);
		uint3 coord;
		coord.z = index / (size * size);
		index -= coord.z * size * size;
		coord.y = index / size;
		coord.x = index % size;
		return coord;
	}
	
	uint GetDepthStartIndex(uint depth)
	{
		uint sum = (depth * (depth + 1)) / 2;
		return ipow(8, sum);
	}

	void GetGLobalNodeIndex(uint depth, uint localIndex)
	{
		uint offset = GetDepthStartIndex(depth);
		///First Depth Level Node
		uint id = offset + localIndex;
	}
	void GetGLobalNodeIndex(uint depth, uint3 coord)
	{
		uint offset = GetDepthStartIndex(depth);
		///First Depth Level Node
		uint id = offset + GetLocalIndexFromCoord(depth, coord);
	}

	uint3 GetChildCoord(uint3 p)
	{
		return uint3{ p.x * 2, p.y * 2, p.z * 2 };
	}
	uint3 GetParentCoord(uint3 c)
	{
		return uint3{ c.x / 2, c.y / 2, c.z / 2 };
	}

	RDynamicVector<Node> Nodes;
};


#if 0
class BVH
{

	struct Node
	{
		AABB bounds;
		Node* left;
		Node* right;
		std::vector<Object> objects;
	};

public:
	BVH(const std::vector<Object>& objects, int max_depth) :
		m_max_depth(max_depth), m_root(BuildTree(objects, 0)) {
	}

	~BVH() {
		DeleteTree(m_root);
	}

	void Traverse()
	{
		TraverseHelper(m_root);
	}

private:
	Node* BuildTree(const std::vector<Object>& objects, int depth)
	{
		Node* node = new Node();
		node->bounds = ComputeBoundingBox(objects);
		if (objects.size() == 1 || depth == m_max_depth)
		{
			node->objects = objects;
			DirectX::BoundingOrientedBox
		}
		else
		{
			// Find the axis with the largest variance
			int axis = GetMaxVarianceAxis(objects, node->bounds);
			// Sort the objects along that axis
			std::sort(objects.begin(), objects.end(),
				[axis](const Object& a, const Object& b) {
				return a.position[axis] < b.position[axis];
			});

			// Partition the objects into two groups
			auto mid = objects.size() / 2;
			auto left = std::vector<Object>(objects.begin(), objects.begin() + mid);
			auto right = std::vector<Object>(objects.begin() + mid, objects.end());

			// Recursively build the left and right subtrees
			node->left = BuildTree(left, depth + 1);
			node->right = BuildTree(right, depth + 1);
		}

		return node;
	}

	void DeleteTree(Node* node) {
		if (node) {
			DeleteTree(node->left);
			DeleteTree(node->right);
			delete node;
		}
	}

	void TraverseHelper(Node* node) {
		if (node) {
			// Process the objects in this node
			for (auto& obj : node->objects) {
				// Do something with obj
			}

			// Traverse the left and right subtrees
			TraverseHelper(node->left);
			TraverseHelper(node->right);
		}
	}

	int GetMaxVarianceAxis(const std::vector<Object>& objects, const AABB& bounds) {
		float3 extents = bounds.Extents;
		int axis = 0;
		float max_var = extents.x;
		if (extents.y > max_var) {
			axis = 1;
			max_var = extents.y;
		}
		if (extents.z > max_var) {
			axis = 2;
		}
		return axis;
	}

	AABB ComputeBoundingBox(const std::vector<Object>& objects) {
		AABB bbox;
		float3 min;
		float3 max;
		min = float3(INFINITY);
		max = float3(-INFINITY);
		for (auto& obj : objects)
		{
			min = Min(min, obj.position - obj.radius);
			max = Max(max, obj.position + obj.radius);
		}
		return bbox;
	}

	int m_max_depth;
	Node* m_root;
};
#endif


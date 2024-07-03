
void IntegrateVerlet(inout float3 curPos, inout float3 prevPos, out float3 outVelocity, float3 curAcceleration, float damping, float dt)
{
    float3 prevVelocity = (curPos - prevPos);
    prevPos = curPos;
    prevVelocity = prevVelocity * damping;
    curPos = curPos + prevVelocity + curAcceleration * dt * dt;
	outVelocity = (curPos - prevPos) / dt;
}

#define FORWARD_PBD_SCALE 10.f

void IntegrateEuler(inout float3 curPos, inout float3 prevPos, inout float3 curVelocity, float3 curAcceleration, float damping, float dt)
{
	//add a bit of PBD
	float3 prevVelocity = (curPos - prevPos);
	curVelocity += prevVelocity * FORWARD_PBD_SCALE;

	curVelocity = curVelocity * damping;

	curVelocity = curVelocity + curAcceleration * dt;
	curPos = curPos + curVelocity * dt;

	prevPos = curPos;//cache pos state before collisions
}

void IntegrateRK3(inout float3 curPos, inout float3 prevPos, inout float3 curVelocity, float3 curAcceleration, float damping, float dt)
{
	//add a bit of PBD
	float3 prevVelocity = (curPos - prevPos);
	curVelocity += prevVelocity * FORWARD_PBD_SCALE;

	curVelocity = curVelocity * damping;

	prevPos = curPos;

    float3 k1v = curAcceleration;
    float3 k2v = (curAcceleration + 0.5f * k1v * dt);
    float3 k3v = (curAcceleration + 0.75f * k2v * dt);

    curVelocity = curVelocity + (2.0f * k1v * dt + 3.0f * k2v * dt + 4.0f * k3v * dt) / 9.0f;

    float3 k1p = curVelocity;
    float3 k2p = (curVelocity + 0.5f * k1p * dt);
    float3 k3p = (curVelocity + 0.75f * k2p * dt);

    curPos = curPos + (2.0f * k1p * dt + 3.0f * k2p * dt + 4.0f * k3p * dt) / 9.0f;

	prevPos = curPos;//cache pos state before collisions
}

#define PLANE_INTERSECT_INSIDE 1
#define PLANE_INTERSECT_OUTSIDE 2
#define PLANE_INTERSECT_INTERSECT 3

struct BoundingPlane
{
    float3 n;
    float d;

    float GetSignedDistanceFromPoint(float3 pos)
    {
        return dot(n, pos) + d;
    }

    uint IntersectsSphere(float3 center, float radius, out float absoluteDistance)
	{
		float signedDist = GetSignedDistanceFromPoint(center);

		absoluteDistance = abs(signedDist);

		if (absoluteDistance <= radius)
		{
			return PLANE_INTERSECT_INTERSECT;
		}
		if (signedDist > 0)
		{
			return PLANE_INTERSECT_INSIDE;
		}
		return PLANE_INTERSECT_OUTSIDE;
	}
};

bool ResolvePointPlaneCollision(inout float3 curPosition, float radius, BoundingPlane plane, float CollisionResponseCoef)
{
    float curDist = 0;
    uint Intersectiontype = plane.IntersectsSphere(curPosition, radius, curDist);
    if (Intersectiontype != PLANE_INTERSECT_INSIDE)
    {
        float delta = radius + curDist;
        if (Intersectiontype == PLANE_INTERSECT_INTERSECT)
        {
            delta = radius - curDist;
        }
        curPosition += (plane.n * delta * CollisionResponseCoef);
		return true;
    }
	return false;
}

struct BoundingCapsule
{
	float3 Pos1;
	float3 Pos2;
	float Radius;
};

void ResolvePointCapsuleCollision(inout float3 curPosition, const float radius, const BoundingCapsule capsule, float CollisionResponseCoef) 
{
	const float3 capsuleStart = capsule.Pos1;
	const float3 capsuleEnd = capsule.Pos2;

    float3 diff = curPosition - capsuleStart;
    float3 capsuleDir = capsuleEnd - capsuleStart;
    float dirLengthSq = dot(capsuleDir,capsuleDir);
    float t = max(0.0f, min(1.0f, dot(diff,capsuleDir) / dirLengthSq));
    float3 closestPoint = { capsuleStart.x + t * capsuleDir.x, capsuleStart.y + t * capsuleDir.y, capsuleStart.z + t * capsuleDir.z };
	float3 dir = curPosition - closestPoint;
    float distSq = dot(dir,dir);
    if (distSq <= (radius + capsule.Radius) * (radius + capsule.Radius))
    {
        distSq = sqrt(distSq);
        dir = dir / distSq;
        float delta = ((radius + capsule.Radius) - distSq) * CollisionResponseCoef;
		curPosition += (dir * delta * CollisionResponseCoef);
    }
}

bool RaySphereIntersection(const float3 rayOrigin, const float3 pointOnRay, const float3 spherePos, float sphereRadius)
{
	float3 rayDir = normalize(pointOnRay - rayOrigin);
	float3 m = rayOrigin - spherePos;
	float b = dot(m, rayDir);
	float c = dot(m, m) - sphereRadius * sphereRadius;
	if (c > 0.0f && b > 0.0f)
	{
		return false;
	}
	float discr = b * b - c;
	if (discr < 0.0f)
	{
		return false;
	}
	return true;
}

float GetClosestDistanceFromPointToLineSegment(float3 p, float3 p0, float3 p1, out float3 closestPoint)
{
	float3 d = p1 - p0;
	// Find the parameter t that minimizes the distance from p to the edge
	float t = dot(p - p0, d) / dot(d, d);

	// Clamp t to [0, 1] to get the closest point on the edge
	t = max(0.0f, min(1.0f, t));

	// Find the distance from p to the closest point on the edge
	closestPoint = (p0 + d * t);
	return length(p - closestPoint);
}

bool LineSegmentSphereIntersection(float3 lineStart, float3 lineEnd, float3 spherePos, float sphereRadius)
{
	float3 dummy;
	return GetClosestDistanceFromPointToLineSegment(spherePos, lineStart, lineEnd, dummy) < sphereRadius;
}



bool SphereSphereIntersection(float3 spherePos, float sphereRadius, float3 spherePos2, float sphereRadius2)
{
	// Distance  between centers.
	float d = length(spherePos - spherePos2);

	// Sum of the rad.
	float r = sphereRadius + sphereRadius2;

	return d <= r;
}



///n-body gravity simulation
float3 GetGravityAxis(float3 dirVec)
{
	const float epsSquared = 1.f;//parameter
	return dirVec / pow(dot(dirVec, dirVec) + epsSquared, 3.f / 2.f);
}
				

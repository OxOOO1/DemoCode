#include "Physics.h"

#include <execution>

#include "InterpolationAndCurves.h"

float3 RConnectedPointsSimulation::ComputeCurWindForce()
{
    float t = Time::Get().GetTimeAfterLaunchSec();
    float main = sin(t * WindStrengthFrequency) + 1.f;
    float freq = main;
    freq += (sin(t * 2.3f * WindStrengthFrequency) + 1.f) * 0.5;
    freq += main * WindStrengthFrequency2Weight * 0.1f * (sin(t * 20.f) + 1.f) * 0.5f;
    ///Interpolate between angular endpoints
    const float2 pitchMinMax{ -PiDiv2, PiDiv2 };
    const float2 yawMinMax{ 0, PiX2 };
    float pitchAngleCur = TLinearInterpolate(sin(t / WindDirectionFrequency * 0.4f * 0.01f), pitchMinMax.x, pitchMinMax.y, -1.f, 1.f);
    float yawAngleCur = TLinearInterpolate(sin(t / WindDirectionFrequency * 0.01f), yawMinMax.x, yawMinMax.y, -1.f, 1.f);

    float3 windDirCur = MapSphericalToCartesian(yawAngleCur, pitchAngleCur, 1.f);

    WindForceCur = windDirCur * WindStrength * freq;
    return WindForceCur;

}

//void RConnectedPointsSimulation::SpatialHashGrid::ClearVectors()
//{
//    PROFILE_FUNC();
//#if GRID_USE_HASH_TABLE
//    CellToObjectsMap.clear();
//#else=
//    /*for (auto& vec : CellToObjectsMap)
//    {
//        vec.clear();
//    }*/
//    std::for_each(std::execution::par, CellToObjectsMap.begin(), CellToObjectsMap.end(), [](RDynamicVector<RPhysicsObject*>& vec) {vec.clear(); });
//#endif
//}

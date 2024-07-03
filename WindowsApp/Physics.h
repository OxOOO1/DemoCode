#pragma once
#include "Calculus.h"
#include "Intersections.h"
#include "ControlPoints.h"

static inline void UpdatePositionSmooth(float& CurPos, float DesiredPos,
    float AccelerationScale, float dtime
)
{
    float curDiff = (CurPos - DesiredPos);
    float CurVelocity = -AccelerationScale * curDiff;
    CurPos += CurVelocity * dtime;
}

static inline void UpdatePositionSmooth(float& CurPos, float& CurVelocity, float DesiredPos,
    float AccelerationScale, float DampingScale, float dtime
)
{
    float curDiff = (CurPos - DesiredPos);
    float acceleration = -AccelerationScale * curDiff - DampingScale * CurVelocity;

    CurVelocity += acceleration * dtime;
    CurPos += CurVelocity * dtime;
}





///Draw UI with pre-simulation settings
///Run simulation
///Draw UI only for running simulation settings

class RPhysicsSimulation
{
public:


    bool bSimulationRunning = false;
    bool bSimulateOnce = false;
    bool bSimulationStartedThisFrame = false;
    bool bResetSimulation = false;
    float SimulationSpeedScale = 1.0f;

    int NumSubSteps = 1;

    bool bSimulateIn2D = false;

    void DrawUI()
    {
        if (ImGui::Checkbox("SimulationRunning", &bSimulationRunning))
        {
            if (bSimulationRunning)
            {
                bSimulationStartedThisFrame = true;
            }
        }
        else
        {
            bSimulationStartedThisFrame = false;
        }
        bSimulateOnce = ImGui::Button("Simulate Once");
        ImGui::Checkbox("2D", &bSimulateIn2D);
        if (ImGui::Button("Reset Simulation"))
        {
            bResetSimulation = true;
        }
        ImGui::DragFloat("SimulationSpeed", &SimulationSpeedScale, 0.001);
        ImGui::SliderInt("Num SubSteps", (int*)&NumSubSteps, 1, 8);
    }

};

///Physics simulation using Position Based Dynamics with Verlet Integration
class RConnectedPointsSimulation : public RPhysicsSimulation
{
public:
    class RConnection;

    struct RPhysicsObject
    {

        RPhysicsObject() = default;

        explicit RPhysicsObject(float3 pos)
        {
            CurPosition = pos;
        }

        void SetPos(c_float3& pos)
        {
            CurPosition = pos;
        }

        void Reset()
        {
            CurAcceleration = float3{ 0,0,0 };
            CurVelocity = float3{ 0,0,0 };
            CurPosition = float3{ 0,0,0 };
            PrevPosition = float3{ 0,0,0 };
        }

        void Translate(c_float3& displacement)
        {
            if (!bIsPinned)
            {
                CurPosition = CurPosition + displacement;
            }
        }
        void TranslateNegate(c_float3& displacement)
        {
            if (!bIsPinned)
            {
                CurPosition = CurPosition - displacement;
            }
        }
        void AddAcceleration(c_float3& acceleration)
        {
            if (!bIsPinned)
            {
                CurAcceleration = CurAcceleration + acceleration;
            }
        }

        float3 CurPosition;
        float3 PrevPosition;
        float3 CurVelocity;
        float3 CurAcceleration;

        float3 SurfaceNormal{ 1,0,0 };

        bool bCurrentlyUnderCollision = false;

        bool bIsPinned = false;

        RDynamicVector<class RConnectedPointsSimulation::RConnection*> AllConnections;


        const RPhysicsObject* LastCollisionObject = nullptr;

        uint3 IndexInSpatialGrid;
    };

    static inline RDynamicVector<RPhysicsObject> PhysicsObjects;

    ///Control Points
    bool bPinObjects = true;

    ///External Force
    float3 ExternalForce{0, -100.f, 0};

    ///Acceleration Settings
    float AccelerationScale = 1.f;
    
    ///Damping
    bool bUniversalDamping = true;
    float VelocityDampingScale = 0.99f;

    ///Spring Settings
    static inline bool bSinglyLinkedEdges = false;///Use Links withing Points, this way force is single directional
    static inline bool bSpringDistanceSquared = false;
    static inline float ConnectionConstraintStiffness = 1.f;
    static inline float SpringRestLength = 0.5f;

    ///Spring Damping
    static inline bool bSpringDamping = true;
    static inline float SpringDampingFactor = 0.5f;

    ///Air Resistance
    bool bUseAirResistance = false;
    float AirResistanceScale = 1.f;

    void SimulateAndRender(RGeometry& visualizer, float dt)
    {
        PROFILE_FUNC();

        DrawUI();

        if (bEnableBounds)
        {
            GenerateBounds();///TODO:Execute only when needed
        }

        if (bResetSimulation)
        {
            ResetSimulation();
            bResetSimulation = false;
        }

        if (bSimulationRunning || bSimulateOnce)
        {
            PROFILE_SCOPE(Simulate);
            float subDt = dt / float(NumSubSteps);
            for (int i = 0; i < NumSubSteps; i++)
            {
                Simulate(subDt * SimulationSpeedScale);
            }
            
            bSimulateOnce = false;
        }

        if (bPinObjects)
        {
            PinObjects();
        }
        else
        {
            UnpinAllObjects();
        }

        Visualize(visualizer);
    }


///private:

    void DrawUI()
    {
        if (ImGui::Begin("Physics Simulation"))
        {
            RPhysicsSimulation::DrawUI();
            
            ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Integration Method:");
            ImGui::RadioButton("Euler", (int*)&IntegrationMethod, EIntegrationMethod::Euler);
            ImGui::SameLine();
            ImGui::RadioButton("Verlet", (int*)&IntegrationMethod, EIntegrationMethod::Verlet);

            ImGui::NewLine();
            ImGui::TextColored(ImVec4(1, 0.5, 0.5, 1), "Points Generation");
            if (ImGui::Button("Regenerate Points"))
            {
                GeneratePoints();
            }
            if (ImGui::Button("Spawn Point"))
            {
                SpawnPoint();
            }
            bool bGenModeChanged = false;
            bGenModeChanged |= ImGui::RadioButton("Rope", (int*)&GeneratorMode, EGeneratorMode::Rope);
            ImGui::SameLine();
            bGenModeChanged |= ImGui::RadioButton("Cloth", (int*)&GeneratorMode, EGeneratorMode::Cloth);
            ImGui::SameLine();
            bGenModeChanged |= ImGui::RadioButton("SoftBody", (int*)&GeneratorMode, EGeneratorMode::SoftBody);
            ImGui::SameLine();
            bGenModeChanged |= ImGui::RadioButton("Point", (int*)&GeneratorMode, EGeneratorMode::Point);
            if (GeneratorMode == EGeneratorMode::Rope)
            {
                bGenModeChanged |= ImGui::DragInt("Num Points in Rope", (int*)&NumRopePoints, 1, 2, 100);
            }
            if (GeneratorMode == EGeneratorMode::Cloth)
            {
                ImGui::Checkbox("Cloth: Normals", &bUseNormalsForCloth);
                ImGui::Checkbox("Cloth: Mesh", &bRenderClothMesh);
                bGenModeChanged |= ImGui::DragInt2("Num knots in Cloth", (int*)&ClothSize, 1, 2, 100);
            }
            if (bGenModeChanged)
            {
                GeneratePoints();
            }

            ImGui::NewLine();
            ImGui::TextColored(ImVec4(1, 0.2, 0.2, 1), "Spatial Control / Pin");
            ImGui::Checkbox("Pin Objects", &bPinObjects);

            ImGui::NewLine();

            ImGui::DragFloat("Acceleration Scale", &AccelerationScale, 0.001, 0.001, 1);

            ImGui::NewLine();

            if (ImGui::BeginTabBar("PhysTabBar"))
            {
                ImGui::NewLine();
                ///Colission&Bounds
                if (ImGui::BeginTabItem("Colission&Bounds"))
                {
                    ImGui::DragFloat("Point Size", &PointSize, 0.001, 0.1, 10.f);
                    if (ImGui::Button("PointSize To SpringLength"))
                    {
                        PointSize = SpringRestLength * 0.5f;
                    }
                    ImGui::Checkbox("Point Collision", &bPointPointCollision);
                    ImGui::Checkbox("Edge Collision", &bEdgeEdgeCollision);
                    ImGui::Checkbox("Triangle Collision", &bPointTriangleCollision);
                    ImGui::SliderInt("Num Collision Iterations", (int*)&NumCollisionIterations, 1, 8);
                    ImGui::DragFloat("Collision Response Coefficient", &CollisionResponseCoef, 0.001, 0.0, 10.f);
                    ImGui::Checkbox("Collisions Add Impulse", &bCollisionsAddImpulse);
                    if (bCollisionsAddImpulse)
                    {
                        ImGui::DragFloat("Collision Impulse Strength", &CollisionImpulseStrength, 0.001, 0.001, 10.f);
                    }
                    if (ImGui::Checkbox("Enable Bounds", &bEnableBounds))
                    {
                    }
                    if (bEnableBounds)
                    {
                        ImGui::Checkbox("Visualize Bounds", &bVisualizeBounds);
                        ImGui::DragFloat("Bounds Size", &BoundsSize, 1., 1, 50.f);
                        ImGui::DragFloat("Bounds Response Coefficient", &BoundsCollisionResponseCoef, 0.001, 0.0, 10.f);
                        ImGui::Checkbox("Rigid Sphere", &bEnableSphereCollider);
                        if (bEnableSphereCollider)
                        {
                            ImGui::DragFloat3("RSPhere Pos", &RigidSpherePos.x, 0.1);
                            ImGui::DragFloat("RSphere Radius", &RigidSphereRadius, 0.1, 0.1, 10.f);
                        }
                        
                        ImGui::NewLine();
                        if (bEnableCapsuleCollider)
                        {
                            ImGui::DragFloat3("Capsule Pos", &CapsulePos.x, 0.1);
                            ImGui::DragFloat3("Capsule Offset1", &CapsuleOffset1.x, 0.1);
                            ImGui::DragFloat3("Capsule Offset2", &CapsuleOffset2.x, 0.1);
                            ImGui::DragFloat("Capsule Radius", &CapsuleRadius, 0.1, 0.1, 10.f);
                        }
                        ImGui::Checkbox("Capsule Collider", &bEnableCapsuleCollider);
                        
                    }

                    ImGui::NewLine();
                    ImGui::TextColored(ImVec4(1, 0.5, 0.2, 1), "Acceleration Structure:");
                    ImGui::Checkbox("Use Acceleration Structure", &bUseAccelerationStructure);
                    ImGui::Checkbox("Store Object In Multiple Grids", &bHashGridCanWriteToMultipleGrids);
                    ImGui::DragFloat("Hash Grid Cell Size", &HashGrid.CellSize, 0.1, 0.1, 10);
                    if (ImGui::Button("Cell Size To SpringLength"))
                    {
                        HashGrid.CellSize = SpringRestLength;
                    }
                    ImGui::EndTabItem();
                }

                ///Constraints
                if (ImGui::BeginTabItem("Constraints"))
                {
                    ImGui::Checkbox("Use Constraints", &bEnableConnectionConstraints);
                    if (bEnableConnectionConstraints)
                    {
                        ImGui::SliderFloat("Tear Threshold", &SpringTearThreshold, 1.f, 25.f);
                        ImGui::DragFloat("Connection Rest Length", &SpringRestLength, 0.01, 0.0, 10);
                        ImGui::DragFloat("Connection Stiffness", &ConnectionConstraintStiffness, 0.001, -10, 10);
                        ImGui::Checkbox("Single Link", &bSinglyLinkedEdges);
                        ImGui::Checkbox("Use Spring Constraints", &bSpringConstraints);
                        if (bSpringConstraints)
                        {
                            ImGui::Checkbox("Use Distance Squared", &bSpringDistanceSquared);
                        }
                    }
                    ImGui::EndTabItem();
                }

                ///Damping
                if (ImGui::BeginTabItem("Damping"))
                {
                    ImGui::Checkbox("Enable universal Velocity Damping", &bUniversalDamping);
                    if (bUniversalDamping)
                    {
                        ImGui::DragFloat("Velocity Damping Scale", &VelocityDampingScale, 0.001, 0.0, 1);
                    }
                    if (bSpringConstraints)
                    {
                        ImGui::Checkbox("Enable Spring Damping", &bSpringDamping);
                        ImGui::DragFloat("Spring Damping Scale", &SpringDampingFactor, 0.001, -10, 10);
                    }
                    ImGui::Checkbox("Enable Drag / Air Resistance", &bUseAirResistance);
                    if (bUseAirResistance)
                    {
                        ImGui::DragFloat("Air Resistance Scale", &AirResistanceScale, 0.001, 0.001, 1);
                    }
                    ImGui::EndTabItem();
                }

                ///External Forces
                if (ImGui::BeginTabItem("External Forces"))
                {
                    bApplyRandomImpulse = ImGui::Button("Apply Random Impulses");
                    ImGui::DragFloat3("Additional Constant Force", &ExternalForce.x, 0.5, -250, 50);
                    ImGui::Checkbox("Wind", &bApplyWindForce);
                    if (bApplyWindForce)
                    {
                        ImGui::DragFloat("Wind Strength", &WindStrength, 0.1, 0.1, 50.f);
                        ImGui::DragFloat("Wind Strength Freq", &WindStrengthFrequency, 0.01, 0.01, 1.f);
                        ImGui::DragFloat("Wind Direction Freq", &WindDirectionFrequency, 0.01, 0.01, 1.f);
                    }
                    ImGui::EndTabItem();
                }


                // Add more tabs here if needed
                ImGui::EndTabBar();
            }

        }
        ImGui::End();

    }

    bool bVisualizeVelocity = false;
    float VisualizedVelocityArrowScale = 1.f;

    bool bRenderClothMesh = true;

    void Visualize(RGeometry& visualizer)
    {
        PROFILE_FUNC();

        for (auto& object : PhysicsObjects)
        {
            if (GeneratorMode != EGeneratorMode::Cloth || !bRenderClothMesh)
            {
                if (object.bCurrentlyUnderCollision)
                {
                    visualizer.SetColor(EColor::RedOrange);
                }
                else
                {
                    visualizer.SetColor(EColor::Orange);
                }

                visualizer.AddCircleCameraFacing(object.CurPosition, PointSize);
            }

            if (bVisualizeVelocity)
            {
                visualizer.SetColor(EColor::LightBlue);
                float length = VectorGetLength(object.CurVelocity);
                float3 v = length > 1 ? (object.CurVelocity / length) : object.CurVelocity;
                visualizer.AddArrowPolyboard(object.CurPosition, object.CurPosition + v * VisualizedVelocityArrowScale);
                visualizer.SetColor(EColor::Red);
                length = VectorGetLength(object.CurAcceleration);
                v = length > 1 ? (object.CurAcceleration / length) : object.CurAcceleration;
                visualizer.AddArrowPolyboard(object.CurPosition, object.CurPosition + v * VisualizedVelocityArrowScale);
            }
            
        }

        if (GeneratorMode == EGeneratorMode::Cloth && bRenderClothMesh)
        {
            static float3 clothColor = GenerateRandomColor();
            VisualizeClothMesh(clothColor, visualizer);
        }

        ///Connections
        if (bEnableConnectionConstraints)
        {
            if (GeneratorMode != EGeneratorMode::Cloth || !bRenderClothMesh)
            {
                for (auto& edge : ConnectionsBetweenPoints)
                {
                    if (!edge.bIsBroken)
                    {
                        if (edge.bIsTense)
                        {
                            visualizer.SetColor(EColor::RedOrange);
                        }
                        else if (edge.bIsIntersecting)
                        {
                            visualizer.SetColor(EColor::Red);
                        }
                        else
                        {
                            visualizer.SetColor(EColor::YellowCanary);
                        }

                        visualizer.AddLine(edge.pObjectA->CurPosition, edge.pObjectB->CurPosition);
                    }

                }
            }
        }
        
        if (bEnableBounds && bVisualizeBounds)
        {
            VisualizeBounds(visualizer);
        }

        if (bApplyWindForce)
        {
            float3 start = VectorNegate(WindForceCur * 0.5f);
            visualizer.AddArrowPolyboard(start, start + WindForceCur);
        }


    }


    void VisualizeClothMesh(c_float3& color, RGeometry& visualizer)
    {
        PROFILE_FUNC();

        visualizer.SetColorRaw(color);

        auto RenderTriangle = [&visualizer](RPhysicsObject& point, RPhysicsObject& pointRight, RPhysicsObject& pointUp, bool bUseNormals)
        {
            bool bCanRender = true;
            for (auto& connection : point.AllConnections)
            {
                ///if connection is related to his current neighbors --> check it
                if (connection->pObjectA == &pointRight || connection->pObjectB == &pointRight || connection->pObjectA == &pointUp || connection->pObjectB == &pointUp)
                {
                    if (connection->bIsBroken)
                    {
                        bCanRender = false;
                    }
                }
            }
            if (bCanRender)
            {
                visualizer.AddTriangleToBuffer(point.CurPosition, pointRight.CurPosition, pointUp.CurPosition);
                if (bUseNormals)
                {
                    visualizer.AddTriangleNormalsToBuffer(point.SurfaceNormal, pointRight.SurfaceNormal, pointUp.SurfaceNormal);
                }
            }
        };

        for (uint y = 0; y < ClothSize.y - 1; y++)
        {
            for (uint x = 0; x < ClothSize.x - 1; x++)
            {
                uint2 coord{ x,y };

                ///TODO: What if we render 4 triangles per point and sample 5 points
                auto& point = PhysicsObjects[GetFlatIndexToObjectArray(coord)];
                auto& pointRight = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 1, 0 })];
                auto& pointUp = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 0, 1 })];
                auto& pointRightUp = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 1, 1 })];


                RenderTriangle(point, pointRight, pointUp, bUseNormalsForCloth);
                RenderTriangle(pointRightUp, pointUp, pointRight, bUseNormalsForCloth);
                

                
                //visualizer.AddQuadToBuffer(point.CurPosition, pointUp.CurPosition, pointRightUp.CurPosition, pointRight.CurPosition);

                /*if (bUseNormalsForCloth)
                {
                    c_float3& normal = point.SurfaceNormal;
                    c_float3& normalRight = pointRight.SurfaceNormal;
                    c_float3& normalUp = pointUp.SurfaceNormal;
                    c_float3& normalightUp = pointRightUp.SurfaceNormal;

                    visualizer.AddQuadNormalsToBuffer(normal, normalUp, normalightUp, normalRight);
                }*/
            }
        }
    }


    ///Vector Field
    std::function<float3(float3)> VectorFieldFunc;
    float ExternalVectorFieldScale = 1.f;
    bool bVectorFieldAffectsAcceleration = true;


    ///********************************************************
    ///                     Simulation
    ///********************************************************


    float RandomImpulseScale = 10.f;
    bool bApplyRandomImpulse = false;
    
    void ApplyRandomImpulse(RPhysicsObject& object)
    {
        auto randomDir = VectorNormalize(float3{ RandomNumberGenerator::Get().NextFloat(-1.f, 1.f), RandomNumberGenerator::Get().NextFloat(-1.f, 1.f), bSimulateIn2D ? 0 : RandomNumberGenerator::Get().NextFloat(-1.f, 1.f) });
        object.AddAcceleration(randomDir * RandomImpulseScale);
    }

    void ResetSimulation()
    {
        bSimulationRunning = false;
        for (auto& object : PhysicsObjects)
        {
            object.Reset();
        }
    }

    void ApplyDamping(RPhysicsObject& object)
    {
        if (!object.bIsPinned)
        {
            if (bUniversalDamping)
            {
                object.CurVelocity = object.CurVelocity * VelocityDampingScale;
            }

            ///Apply Air Resistance
            if (bUseAirResistance)
            {
                float velocityLength = VectorGetLength(object.CurVelocity);
                float airResistance = -AirResistanceScale * velocityLength * velocityLength;
                float3 airResistanceWithDirection = (object.CurVelocity / velocityLength) * airResistance;
                object.AddAcceleration(airResistanceWithDirection);
            }
        }
    }

    bool bApplyWindForce = true;
    float WindStrength = 1.f;
    float3 WindDirectionInitial{ 1.f, 0.f, 0.f };
    float3 WindForceCur{};

    float WindStrengthFrequency = 0.5f;
    float WindStrengthFrequency2 = 1.f;
    float WindStrengthFrequency2Weight = 1.0f;
    float WindDirectionFrequency = 1.f;

    float3 ComputeCurWindForce();

    void ApplyExternalForces(RPhysicsObject& object)
    {
        ///Apply Force from Fields
        ///Const Force:
        object.AddAcceleration(ExternalForce);

        if (bApplyWindForce)
        {
            object.AddAcceleration(WindForceCur);
        }

        if (bApplyRandomImpulse)
        {
            ApplyRandomImpulse(object);
            bApplyRandomImpulse = false;
        }

        ///External Field
        if (VectorFieldFunc != nullptr)
        {
            float3 curFieldValue = VectorFieldFunc(object.CurPosition);

            if (bVectorFieldAffectsAcceleration)
            {
                object.AddAcceleration(curFieldValue * ExternalVectorFieldScale * 0.01);
            }
            else
            {
                ///TODO: Just shift position directly as in Hard Constraint
                ///object.CurVelocity = object.CurVelocity + curFieldValue * ExternalVectorFieldScale * 0.01;
            }
        }
    }

    void Simulate(float dt)
    {

        if (bApplyWindForce)
        {
            ComputeCurWindForce();
        }

        if (bSimulateScissors)
        {
            SimulateScissors();
        }

        for (auto& object : PhysicsObjects)
        {
            ///Remove all prev frame net forces
            object.CurAcceleration = float3{ 0,0,0 };
            object.bCurrentlyUnderCollision = false;

            ApplyDamping(object);

            ApplyExternalForces(object);
        }

        for (auto& object : PhysicsObjects)
        {
            IntegratePosVel(object, dt);
        }

        if (GeneratorMode == EGeneratorMode::Cloth && bPointTriangleCollision)
        {
            GenerateClothTrianglesForCollision();
        }

        if (bResolveCollisionsAfterIntegration)
        {
            ResolveAllColisionsAndConstraints();
        }
        
        if (bSimulateIn2D)
        {
            for (auto& object : PhysicsObjects)
            {
                object.CurPosition.z = 0.f;
                object.PrevPosition.z = 0.f;
            }
        }

        bool bNormalsRequiredThisFrame = bUseNormalsForCloth && (GeneratorMode == EGeneratorMode::Cloth);
        if (bNormalsRequiredThisFrame)
        {
            GenerateNormals();
        }
       
    }

    bool bSimulateScissors = true;
    float ScissorsAffectRadius = 0.5f;

    void SimulateScissors()
    {
        if (Mouse::Get().IsRButtonPressed())
        {
            float3 scissorPos = RSpatialControlSet::Get().GetControllerPosCur();
            
            for (auto& edge : ConnectionsBetweenPoints)
            {
                if (!edge.bIsBroken)
                {
                    if (LineSegmentSphereIntersection(edge.pObjectA->CurPosition, edge.pObjectB->CurPosition, scissorPos, ScissorsAffectRadius))
                    {
                        edge.bIsBroken = true;
                    }
                }
            }
        }
    }

    enum EIntegrationMethod : int
    {
        Euler,///TODO:Collisions & constraints should be resolved before integration
        Verlet,///Requires constant dt!
    } IntegrationMethod{ EIntegrationMethod::Verlet };


    inline void IntegratePosVel(RPhysicsObject& object, float dt)
    {
        if (!object.bIsPinned)
        {
            object.CurAcceleration = object.CurAcceleration * AccelerationScale;

            if ((IntegrationMethod != EIntegrationMethod::Euler) && !bSimulationStartedThisFrame)
            {
                float3 prevVelocity = (object.CurPosition - object.PrevPosition);
                object.PrevPosition = object.CurPosition;
                if (bUniversalDamping)
                {
                    prevVelocity = prevVelocity * VelocityDampingScale;
                }
                object.CurPosition = object.CurPosition + prevVelocity + object.CurAcceleration * dt * dt;
                object.CurVelocity = (object.CurPosition - object.PrevPosition) / dt;
            }
            else
            {
                object.PrevPosition = object.CurPosition;
                object.CurVelocity = object.CurVelocity + object.CurAcceleration * dt;
                object.CurPosition = object.CurPosition + object.CurVelocity * dt;
                ///object.CurPosition = object.CurPosition + object.CurVelocity * dt + object.CurAcceleration * dt * dt;
            }
        }
        else
        {
            object.PrevPosition = object.CurPosition;
        }
        
    }
    

    ///********************************************************
    ///                     Constraints
    ///********************************************************

    static inline bool bSpringConstraints = false;
    static inline bool bEnableConnectionConstraints = true;

    static inline float SpringTearThreshold = 25.f;

    bool bResolveCollisionsAfterIntegration = true;

    struct RConnection
    {

        RConnection() = default;

        RConnection(RPhysicsObject& a, RPhysicsObject& b)
            :pObjectA(&a), pObjectB(&b)
        {
        }

        ///After initialization of dynamic array
        void SetObjectsPointers()
        {
            pObjectA->AllConnections.push_back(this);
            pObjectB->AllConnections.push_back(this);
        }

        void AssignPoints(RPhysicsObject& a, RPhysicsObject& b)
        {
            pObjectA = &a;
            pObjectB = &b;
        }

        RPhysicsObject* pObjectA = nullptr;
        RPhysicsObject* pObjectB = nullptr;

        bool bIsBroken = false;
        bool bIsTense = false;

        bool bIsIntersecting = false;

        ///k(|B - A| - L)
        float GetSpringForce() const
        {
            float curDist = VectorGetLength(pObjectB->CurPosition - pObjectA->CurPosition);
            if (bSpringDistanceSquared)
            {
                curDist = 1.f / (pow(curDist, 2.f) + 1.f);
            }
            else
            {
                curDist = curDist - SpringRestLength;
            }
            return curDist * ConnectionConstraintStiffness * 1000.f;
        }

        float GetDampingForce() const
        {
            ///Determine how fast objects are moving towards or away from each other
            ///Whene objects move away from each other, damp force should act in direction facing each other

            float3 dirFromAtoB = VectorNormalize(pObjectB->CurPosition - pObjectA->CurPosition);
            float3 velBFromA = pObjectB->CurVelocity - pObjectA->CurVelocity;
            float dotAB = VectorDot(dirFromAtoB, velBFromA);
            return dotAB * SpringDampingFactor * 10.f;
        }

        void ResolveConnection()
        {
            if (bIsBroken)
            {
                return;
            }

            check(pObjectA != nullptr);
            check(pObjectB != nullptr);
            check(pObjectA != pObjectB);

            ///TODO:Optimise
            float distC = VectorGetLength(pObjectB->CurPosition - pObjectA->CurPosition);
            if (distC > SpringTearThreshold)
            {
                bIsBroken = true;
            }
            else if ((distC / SpringTearThreshold) > 0.75f )
            {
                bIsTense = true;
            }
            else
            {
                bIsTense = false;
            }



            if (bSpringConstraints)
            {
                float force = GetSpringForce();
                if (bSpringDamping)
                {
                    force += GetDampingForce();
                }
                if (!bSinglyLinkedEdges)///Whether only second endpoint is affected
                {
                    pObjectA->AddAcceleration(VectorNormalize(pObjectB->CurPosition - pObjectA->CurPosition) * force);
                }
                pObjectB->AddAcceleration(VectorNormalize(pObjectA->CurPosition - pObjectB->CurPosition) * force);
            }
            else///position based dynamics constraint
            {
                float3 dir = pObjectB->CurPosition - pObjectA->CurPosition;
                const bool bUseApproximation = false;
                if (bUseApproximation)
                {
                    //using sqrt Taylor approximation
                    const float restLengthSquared = SpringRestLength * SpringRestLength;
                    dir = dir * (restLengthSquared / (VectorDot(dir,dir) + restLengthSquared) - 0.5f);
                }
                else
                {
                    float curDist = VectorGetLength(dir);
                    float delta = curDist - SpringRestLength;
                    dir = (dir / curDist) * delta * 0.5f * ConnectionConstraintStiffness;///TODO: Use different stiffness based on distance, e.g. if distance is < 50%ofRest, use hard Stifness, otherwise use soft
                }

                if (!bSinglyLinkedEdges)///Whether only second endpoint is affected
                {
                    pObjectA->Translate(dir);
                }
                pObjectB->TranslateNegate(dir);
            }

        }
    };

    static inline RDynamicVector<RConnection> ConnectionsBetweenPoints;

    void ResolveAllConnectionConstraints()
    {
        PROFILE_FUNC();
        check(!ConnectionsBetweenPoints.empty());

        for (auto& edge : ConnectionsBetweenPoints)
        {
            edge.ResolveConnection();
        }
    }

    ///********************************************************
    ///                     Collisions
    ///********************************************************

    bool bCollisionsAddImpulse = false;

    bool bVisualizeBounds = true;
    bool bEnableBounds = true;
    bool bFloorBoundsOnly = true;
    float BoundsSize = 5.f;
    RArray<BoundingPlane, 6> BoundingPlanes;

    void GenerateBounds()
    {
        ///6 planes
        ///down
        BoundingPlanes[0] = BoundingPlane(float3{ 0, -BoundsSize, 0 }, float3{ 0,1,0 });
        ///up
        BoundingPlanes[1] = BoundingPlane(float3{ 0, BoundsSize, 0 }, float3{ 0,-1,0 });
        ///left
        BoundingPlanes[2] = BoundingPlane(float3{ -BoundsSize, 0, 0 }, float3{ 1,0,0 });
        ///right
        BoundingPlanes[3] = BoundingPlane(float3{ BoundsSize, 0, 0 }, float3{ -1,0,0 });
        ///forward
        BoundingPlanes[4] = BoundingPlane(float3{ 0, 0, BoundsSize }, float3{ 0,0,-1 });
        ///back
        BoundingPlanes[5] = BoundingPlane(float3{ 0, 0, -BoundsSize }, float3{ 0,0,1 });
    }

    void VisualizeBounds(RGeometry& visualizer)
    {
        for (int i = 0; i < 6; i++)
        {
            visualizer.AddPlaneBounding(BoundingPlanes[i], BoundsSize * 0.5f);
            if (bFloorBoundsOnly)
            {
                break;
            }
        }

        if (bEnableSphereCollider)
        {
            visualizer.AddCircleCameraFacing(RigidSpherePos, RigidSphereRadius);
        }

        if (bEnableCapsuleCollider)
        {
            visualizer.AddCapsule(CapsulePos + CapsuleOffset1, CapsulePos + CapsuleOffset2, CapsuleRadius, CapsuleRadius);
        }

    }

    int NumCollisionIterations = 1;
    float CollisionResponseCoef = 1.f;
    float BoundsCollisionResponseCoef = 0.75f;

    float CollisionImpulseStrength = 0.001f;

    bool bPointPointCollision = true;
    bool bEdgeEdgeCollision = false;


    bool bEnableSphereCollider = false;
    float3 RigidSpherePos;
    float RigidSphereRadius = 1.f;


    bool bEnableCapsuleCollider = true;
    float3 CapsulePos;
    float CapsuleRadius = 2.f;
    float3 CapsuleOffset1{0,2,0};
    float3 CapsuleOffset2;

    void ResolvePointPlaneCollision(RPhysicsObject& curObject, BoundingPlane& plane)
    {
        float curDist = 0;

        BoundingPlane::EIntersection Intersectiontype = plane.IntersectsSphere(curObject.CurPosition, PointSize, curDist);

        if (Intersectiontype != BoundingPlane::EIntersection::Inside)
        {
            curObject.bCurrentlyUnderCollision = true;

            float delta = PointSize + curDist;
            if (Intersectiontype == BoundingPlane::EIntersection::Intersects)
            {
                delta = PointSize - curDist;
            }
            float3 n = plane.GetNormal();
            curObject.Translate(n * delta * BoundsCollisionResponseCoef);

            if (bCollisionsAddImpulse)
            {
                curObject.AddAcceleration(n * delta * CollisionImpulseStrength);
            }
        }
    }

    void ResolvePointCapsuleCollision(RPhysicsObject& curObject, const float3& capsuleStart, const float3& capsuleEnd, float capsuleRadius) 
    {
        c_float3& sphereCenter = curObject.CurPosition;
        float sphereRadius = PointSize;

        float3 diff = sphereCenter - capsuleStart;
        float3 dir = capsuleEnd - capsuleStart;
        float dirLengthSq = VectorGetLengthSquared(dir);

        // If the capsule has zero length, treat it as a sphere
        /*if (dirLengthSq < 0.01f) 
        {
            float distSq = VectorGetLengthSquared(diff);
            return distSq <= (sphereRadius + capsuleRadius) * (sphereRadius + capsuleRadius);
        }*/

        float t = fmaxf(0.0f, fminf(1.0f, VectorDot(diff,dir) / dirLengthSq));
        float3 closestPoint = { capsuleStart.x + t * dir.x, capsuleStart.y + t * dir.y, capsuleStart.z + t * dir.z };
        float distSq = VectorGetLengthSquared(sphereCenter - closestPoint);

        if (distSq <= (sphereRadius + capsuleRadius) * (sphereRadius + capsuleRadius))
        {
            curObject.bCurrentlyUnderCollision = true;

            distSq = sqrtf(distSq);
            c_float3 dir = (sphereCenter - closestPoint) / distSq;
            float delta = ((sphereRadius + capsuleRadius) - distSq) * CollisionResponseCoef;

            curObject.Translate(dir * delta);
        }

    }

    static inline float PointSize = 0.1f;

    bool bPointTriangleCollision = false;

    void ResolveEdgeEdgeCollision()
    {
        for (auto& edge : ConnectionsBetweenPoints)
        {
            edge.bIsIntersecting = false;
        }

        ///edge-edge intersection
        auto numEdges = ConnectionsBetweenPoints.size();
        for (int i = 0; i < numEdges; i++)
        {
            auto& curEdge = ConnectionsBetweenPoints[i];
            RIntersectionCheck::Edge curLineSegment{ curEdge.pObjectA->CurPosition, curEdge.pObjectB->CurPosition };


            if (!curEdge.bIsBroken)
            {
                for (int j = i+1; j < numEdges; j++)
                {
                    auto& otherEdge = ConnectionsBetweenPoints[j];
                    RIntersectionCheck::Edge otherLineSegment{ otherEdge.pObjectA->CurPosition, otherEdge.pObjectB->CurPosition };

                    if (!otherEdge.bIsBroken)
                    {
                        float3 collisionPoint1;
                        float3 collisionPoint2;
                        if (EdgeEdgeIntersection(curLineSegment, otherLineSegment, collisionPoint1, collisionPoint2))
                        {
                            ///displace edges with collision points deltas
                            c_float3 collisionAxis = collisionPoint1 - collisionPoint2;
                            float distance = VectorGetLength(collisionAxis);

                            float radiusSum = PointSize * 0.1f;

                            if (distance < radiusSum)///push objects from each other
                            {

                                curEdge.bIsIntersecting = true;
                                otherEdge.bIsIntersecting = true;

                                c_float3 dir = collisionAxis / distance;
                                float delta = (radiusSum - distance) * 0.5f * CollisionResponseCoef;

                                curEdge.pObjectA->Translate(dir * delta);
                                curEdge.pObjectB->Translate(dir * delta);
                                otherEdge.pObjectA->TranslateNegate(dir * delta);
                                otherEdge.pObjectB->TranslateNegate(dir * delta);

                                //bSimulationRunning = false;

                            }
                        }
                    }
                }
            }

        }
    }

    struct Triangle
    {
        float3 p0, p1, p2;
        float3 n0, n1, n2;

        RPhysicsObject* o0;
        RPhysicsObject* o1;
        RPhysicsObject* o2;
    };

    RDynamicVector<Triangle> ClothTriangles;

    void GenerateClothTrianglesForCollision()
    {
        PROFILE_FUNC();

        ClothTriangles.clear();

        ///Only add triangles that don't have broken edges
        auto AddTriangle = [&](RPhysicsObject& point, RPhysicsObject& pointRight, RPhysicsObject& pointUp, bool bUseNormals)
        {
            bool bCanRender = true;
            for (auto& connection : point.AllConnections)
            {
                ///if connection is related to his current neighbors --> check it
                if (connection->pObjectA == &pointRight || connection->pObjectB == &pointRight || connection->pObjectA == &pointUp || connection->pObjectB == &pointUp)
                {
                    if (connection->bIsBroken)
                    {
                        bCanRender = false;
                    }
                }
            }
            if (bCanRender)
            {
                Triangle tr;
                tr.p0 = point.CurPosition;
                tr.p1 = pointRight.CurPosition;
                tr.p2 = pointUp.CurPosition;
                if (bUseNormals)
                {
                    tr.n0 = point.SurfaceNormal;
                    tr.n1 = pointRight.SurfaceNormal;
                    tr.n2 = pointUp.SurfaceNormal;
                }

                tr.o0 = &point;
                tr.o1 = &pointRight;
                tr.o2 = &pointUp;

                ClothTriangles.emplace_back(tr);
            }
        };

        for (uint y = 0; y < ClothSize.y - 1; y++)
        {
            for (uint x = 0; x < ClothSize.x - 1; x++)
            {
                uint2 coord{ x,y };

                ///TODO: What if we render 4 triangles per point and sample 5 points
                auto& point = PhysicsObjects[GetFlatIndexToObjectArray(coord)];
                auto& pointRight = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 1, 0 })];
                auto& pointUp = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 0, 1 })];
                auto& pointRightUp = PhysicsObjects[GetFlatIndexToObjectArray(coord + uint2{ 1, 1 })];

                AddTriangle(point, pointRight, pointUp, bUseNormalsForCloth);
                AddTriangle(pointRightUp, pointUp, pointRight, bUseNormalsForCloth);
            }
        }
    }


    void ResolvePointTriangleCollisions(RPhysicsObject& object, RConnectedPointsSimulation::Triangle& triangle)
    {
        bool triangleOfTheObject = false;
        if (triangle.o0 == &object || triangle.o1 == &object || triangle.o2 == &object)
        {
            triangleOfTheObject = true;
        }

        float radius = PointSize;

        if (!triangleOfTheObject)
        {
            float3 collisionPoint;
            float3 barycentric;
            if (SphereTriangleIntersectionBarycentric(object.CurPosition, radius, triangle.p0, triangle.p1, triangle.p2, collisionPoint, barycentric))
            {
                object.bCurrentlyUnderCollision = true;

                c_float3 collisionAxis = object.CurPosition - collisionPoint;
                float distance = VectorGetLength(collisionAxis);
                //c_float3 dir = collisionAxis / distance;
                c_float3 dir = VectorNegate(VectorNormalize(object.CurVelocity));
                float delta = (radius - distance) * CollisionResponseCoef;

                object.Translate(dir * delta);
            }
        }
    }

    void ResolvePointPointCollision(RPhysicsObject& curObject, RPhysicsObject& otherObject)
    {
        c_float3 collisionAxis = curObject.CurPosition - otherObject.CurPosition;
        float distanceSq = VectorGetLengthSquared(collisionAxis);

        float radiusSum = PointSize * 2.f;

        /*if (GeneratorMode == EGeneratorMode::Cloth)
        {
            radiusSum = SpringRestLength - 0.1;
        }*/

        if (distanceSq < radiusSum * radiusSum)///push objects from each other
        {
            curObject.bCurrentlyUnderCollision = true;
            otherObject.bCurrentlyUnderCollision = true;

            distanceSq = sqrtf(distanceSq);
            c_float3 dir = collisionAxis / distanceSq;
            float delta = (radiusSum - distanceSq) * 0.5f * CollisionResponseCoef;

            curObject.Translate(dir * delta);
            otherObject.TranslateNegate(dir * delta);
            if (bCollisionsAddImpulse)
            {
                curObject.AddAcceleration(dir * delta * CollisionImpulseStrength);
                otherObject.AddAcceleration(VectorNegate(dir * delta * CollisionImpulseStrength));
            }
        }
    }

    void ResolveAllColisionsAndConstraints()
    {
        PROFILE_FUNC();
        for (int i = 0; i < NumCollisionIterations; i++)
        {
            ///Resolve Connection Constraints
            if (bEnableConnectionConstraints)
            {
                ResolveAllConnectionConstraints();
            }

            if (bEnableConnectionConstraints && bEdgeEdgeCollision)
            {
                ResolveEdgeEdgeCollision();
            }

            if (bUseAccelerationStructure)
            {
                PROFILE_SCOPE(Hash_Grid_Update);
                HashGrid.CellSize = std::min(HashGrid.CellSize, PointSize * 2);
                HashGrid.ClearVectors(PhysicsObjects);
                for (auto& obj : PhysicsObjects)
                {
                    obj.LastCollisionObject = nullptr;

                    if (bHashGridCanWriteToMultipleGrids)
                    {
                        HashGrid.Add(obj, PointSize);
                    }
                    else
                    {
                        HashGrid.Add(obj);
                    }
                }
            }

            const auto numObjects = PhysicsObjects.size();
            for (auto i = 0; i < numObjects; i++)
            {
                auto& curObject = PhysicsObjects[i];

                ///Resolve collisions between points
                if(bPointPointCollision)
                {
                    if (bUseAccelerationStructure)
                    {
                        static RDynamicVector<RPhysicsObject*> pObjects;
                        pObjects.clear();
                        HashGrid.GetNearest(curObject, pObjects);
                        for (RPhysicsObject* pObj : pObjects)
                        {
                            if (pObj != &curObject)
                            {
                                RPhysicsObject& otherObject = *pObj;
                                ResolvePointPointCollision(curObject, otherObject);
                            }
                        }
                    }
                    else
                    {
                        for (auto j = i + 1; j < numObjects; j++)
                        {
                            auto& otherObject = PhysicsObjects[j];
                            ResolvePointPointCollision(curObject, otherObject);
                        }
                    }
                    
                }

                if (GeneratorMode == EGeneratorMode::Cloth && bPointTriangleCollision)
                {
                    ///test against all valid triangles
                    for (auto& triangle : ClothTriangles)
                    {
                        ResolvePointTriangleCollisions(curObject, triangle);
                    }
                }

                ///Resolve external collisions 
                if (bEnableBounds)
                {
                    ///Bounding planes
                    for (int i = 0; i < 6; i++)
                    {
                        ResolvePointPlaneCollision(curObject, BoundingPlanes[i]);
                        if (bFloorBoundsOnly)
                        {
                            break;
                        }
                    }

                    ///Bounding Sphere
                    if (bEnableSphereCollider)
                    {
                        c_float3 collisionAxis = curObject.CurPosition - RigidSpherePos;
                        float distance = VectorGetLength(collisionAxis);

                        float radiusSum = PointSize + RigidSphereRadius;
                        if (distance < radiusSum)///push objects from each other
                        {
                            curObject.bCurrentlyUnderCollision = true;

                            c_float3 dir = collisionAxis / distance;
                            float delta = (radiusSum - distance) * 0.5f * CollisionResponseCoef;

                            curObject.Translate(dir * delta);
                            if (bCollisionsAddImpulse)
                            {
                                curObject.AddAcceleration(dir * delta * CollisionImpulseStrength);
                            }

                        }
                    }

                    if (bEnableCapsuleCollider)
                    {
                        ResolvePointCapsuleCollision(curObject, CapsulePos + CapsuleOffset1, CapsulePos + CapsuleOffset2, CapsuleRadius);
                    }

                }

            }
        }

    }


    ///********************************************************
    ///                 Generator Settings
    ///********************************************************
    enum EGeneratorMode : int
    {
        Rope,
        Cloth,
        SoftBody,
        Point
    };
    EGeneratorMode GeneratorMode{ RConnectedPointsSimulation::EGeneratorMode::Cloth };

    uint NumRopePoints = 3;
    uint2 ClothSize{64,64};
    uint GetFlatIndexToObjectArray(uint2 coord)
    {
        return coord.y * ClothSize.x + coord.x;
    };

    void SpawnPoint()
    {
        PhysicsObjects.emplace_back();
    }

    bool bUseNormalsForCloth = false;

    void GenerateNormals()
    {
        for (uint y = 0; y < ClothSize.y; y++)
        {
            for (uint x = 0; x < ClothSize.x; x++)
            {
                // Get the neighboring points in the u and v directions
                float3 u1 = (y > 0) ? PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y - 1 })].CurPosition : PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y })].CurPosition;
                float3 u2 = (y < ClothSize.y - 1) ? PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y + 1 })].CurPosition : PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y })].CurPosition;
                float3 v1 = (x > 0) ? PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x - 1, y })].CurPosition : PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y })].CurPosition;
                float3 v2 = (x < ClothSize.x - 1) ? PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x + 1, y })].CurPosition : PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x, y })].CurPosition;

                // Compute the central difference vectors in the u and v directions
                float3 du = u2 - u1;
                float3 dv = v2 - v1;

                PhysicsObjects[GetFlatIndexToObjectArray(uint2{ x,y })].SurfaceNormal = VectorNormalize(VectorCross(du, dv));
            }
        }
    }

    void GeneratePoints()
    {
        ResetSimulation();

        if (GeneratorMode == EGeneratorMode::Rope)
        {
            PhysicsObjects.resize(NumRopePoints);
            ConnectionsBetweenPoints.resize(NumRopePoints-1);

            ///Each subsequent point connected to previous
            for (int i = 0; i < NumRopePoints; i++)
            {
                PhysicsObjects[i].SetPos(float3{float(2.f * i), 2.f, 0 });

                if (i < (NumRopePoints - 1))
                {
                    ConnectionsBetweenPoints[i].AssignPoints(PhysicsObjects[i], PhysicsObjects[i + 1]);
                }
            }
        }
        else if (GeneratorMode == EGeneratorMode::Cloth)
        {

            PhysicsObjects.resize(ClothSize.x * ClothSize.y);
            ConnectionsBetweenPoints.clear();

            float2 startPos{ SpringRestLength * float(ClothSize.x) / -2.f, SpringRestLength * float(ClothSize.y) / -2.f };

            ///Generate knots
            for (uint y = 0; y < ClothSize.y; y++)
            {
                for (uint x = 0; x < ClothSize.x; x++)
                {
                    uint i = GetFlatIndexToObjectArray(uint2{ x,y });
                    PhysicsObjects[i].SetPos(float3{ startPos.x + float(x) * SpringRestLength, startPos.y + float(y) * SpringRestLength, 0 });
                }
            }

            ///Generate edges
            for (uint y = 0; y < ClothSize.y; y++)
            {
                for (uint x = 0; x < ClothSize.x; x++)
                {
                    uint i = GetFlatIndexToObjectArray(uint2{ x,y });
                    if ((y+1) < ClothSize.y)
                    {
                        uint iUp = GetFlatIndexToObjectArray(uint2{ x, y + 1 });
                        ConnectionsBetweenPoints.emplace_back(PhysicsObjects[i], PhysicsObjects[iUp]);
                    }
                    if ((x + 1) < ClothSize.x)
                    {
                        uint iRight = GetFlatIndexToObjectArray(uint2{ x + 1, y });
                        ConnectionsBetweenPoints.emplace_back(PhysicsObjects[i], PhysicsObjects[iRight]);
                    }
                }
            }

            for (auto& edge : ConnectionsBetweenPoints)
            {
                edge.SetObjectsPointers();
            }

            ///Rotate around X
            auto mat = MatrixRotateAroundX(PiDiv2);
            for (auto& object : PhysicsObjects)
            {
                object.CurPosition = VectorMul(object.CurPosition, mat) + float3{ 0, 5.f, 0 };
            }

            GenerateNormals();

            PointSize = SpringRestLength * 0.5f;
            HashGrid.CellSize = SpringRestLength;

            HashGrid.InitMemory();

        }
        else if (GeneratorMode == EGeneratorMode::SoftBody)
        {
            PhysicsObjects.resize(4);
            ConnectionsBetweenPoints.resize(6);

            int i = 0;
            for (int y = 0; y < 2; y++)
            {
                for (int x = 0; x < 2; x++)
                {
                    PhysicsObjects[i].SetPos(float3{ 2.f + float(x), 2.f + float(y),0 });
                    i++;
                }
            }

            ///	2 - 3
            ///	| x |
            ///	0 - 1
            ConnectionsBetweenPoints[0].AssignPoints(PhysicsObjects[0], PhysicsObjects[2]);
            ConnectionsBetweenPoints[1].AssignPoints(PhysicsObjects[2], PhysicsObjects[3]);
            ConnectionsBetweenPoints[2].AssignPoints(PhysicsObjects[3], PhysicsObjects[1]);
            ConnectionsBetweenPoints[3].AssignPoints(PhysicsObjects[0], PhysicsObjects[1]);

            ConnectionsBetweenPoints[4].AssignPoints(PhysicsObjects[0], PhysicsObjects[3]);
            ConnectionsBetweenPoints[5].AssignPoints(PhysicsObjects[1], PhysicsObjects[2]);
        }
        else if (GeneratorMode == EGeneratorMode::Point)
        {
            PhysicsObjects.resize(1);
            ConnectionsBetweenPoints.resize(0);
        }
        else
        {
            check(false);
        }

        AttachSpatialControlPoints();

    }

    void AttachSpatialControlPoints()
    {
        check(!PhysicsObjects.empty());

        auto& controlPoints = RSpatialControlSet::Get().ControlPoints;

        auto numPoints = PhysicsObjects.size();

        controlPoints.resize(numPoints);

        for (int i = 0; i < numPoints; i++)
        {
            controlPoints[i].pPosition = &PhysicsObjects[i].CurPosition;

            RPhysicsObject& objectRef = PhysicsObjects[i];

            controlPoints[i].OnPointSelectedCallback = [&objectRef]()
            {
                objectRef.bIsPinned = true;
            };
        }

    }

    void UnpinAllObjects()
    {
        for (auto& object : PhysicsObjects)
        {
            object.bIsPinned = false;
        }
    }

    void PinObjects()
    {
        UnpinAllObjects();

if (GeneratorMode == EGeneratorMode::Cloth)
{
    for (uint i = 0; i < ClothSize.x; i++)
    {
#if 1
        uint i1 = GetFlatIndexToObjectArray(uint2{ i, ClothSize.y - 4 });
        PhysicsObjects[i1].bIsPinned = true;
#else
        uint i1 = GetFlatIndexToObjectArray(uint2{ i, ClothSize.y - 1 });
        PhysicsObjects[i1].bIsPinned = true;
        i1 = GetFlatIndexToObjectArray(uint2{ i, 4 });
        PhysicsObjects[i1].bIsPinned = true;
#endif

    }
}
else
{
    ///Assign Spatial Control to First Node
    auto& node = PhysicsObjects.front();
    node.bIsPinned = true;
}
    }


    ///********************************************************
    ///                 Acceleration Structures
    ///********************************************************

    bool bUseAccelerationStructure = true;
    bool bHashGridCanWriteToMultipleGrids = false;

#define GRID_USE_HASH_TABLE 1 ///Use HashTable when the Bounds of the Simulation are undefined, otherwise use 3D grid it's more scalable

    class SpatialHashGrid {
    public:

        SpatialHashGrid() = default;

#if GRID_USE_HASH_TABLE
        SpatialHashGrid(float cell_size) : CellSize(cell_size) {}
        void InitMemory()
        {
        }
#else
        SpatialHashGrid(float cell_size, float3 gridExtent) : CellSize(cell_size), GridBoxExtents(gridExtent) 
        {
            InitMemory();
        }

        void InitMemory()
        {
            ///Allocate Mem
            float3 gridSize = (GridBoxExtents * 2) / CellSize;
            uint3 dataSize{ uint(gridSize.x), uint(gridSize.y), uint(gridSize.z) };
            CellToObjectsMap.Resize(dataSize);
        }

        float3 GridBoxExtents;

#endif
        // Add a particle to the grid
        void Add(RPhysicsObject& p)
        {
            
#if GRID_USE_HASH_TABLE
            int cell_x = static_cast<int>(p.CurPosition.x / CellSize);
            int cell_y = static_cast<int>(p.CurPosition.y / CellSize);
            int cell_z = static_cast<int>(p.CurPosition.z / CellSize);
            int cell_id = GetCelllIndex(cell_x, cell_y, cell_z);
            CellToObjectsMap[cell_id].push_back(&p);
#else
            int cell_x = static_cast<int>((p.CurPosition.x + GridBoxExtents.x ) / CellSize);
            int cell_y = static_cast<int>((p.CurPosition.y + GridBoxExtents.y ) / CellSize);
            int cell_z = static_cast<int>((p.CurPosition.z + GridBoxExtents.z ) / CellSize);
            check((cell_x > 0) && (cell_y > 0) && (cell_z > 0));
            auto id = uint3{ uint(cell_x), uint(cell_y), uint(cell_z) };
            CellToObjectsMap[id].push_back(&p);
            p.IndexInSpatialGrid = id;
#endif
            
        }


        ///For objects of radius >0.5
        void Add(RPhysicsObject& p, float radius)
        {
#if GRID_USE_HASH_TABLE
            int cell_x_min = static_cast<int>((p.CurPosition.x - radius) / CellSize);
            int cell_x_max = static_cast<int>((p.CurPosition.x + radius) / CellSize);
            int cell_y_min = static_cast<int>((p.CurPosition.y - radius) / CellSize);
            int cell_y_max = static_cast<int>((p.CurPosition.y + radius) / CellSize);
            int cell_z_min = static_cast<int>((p.CurPosition.z - radius) / CellSize);
            int cell_z_max = static_cast<int>((p.CurPosition.z + radius) / CellSize);

            for (int cell_x = cell_x_min; cell_x <= cell_x_max; ++cell_x)
            {
                for (int cell_y = cell_y_min; cell_y <= cell_y_max; ++cell_y)
                {
                    for (int cell_z = cell_z_min; cell_z <= cell_z_max; ++cell_z)
                    {
                        int cell_id = GetCelllIndex(cell_x, cell_y, cell_z);
                        CellToObjectsMap[cell_id].push_back(&p);
                    }
                }
            }
#else

            float3 pos = p.CurPosition + GridBoxExtents;
            float3 gridSize = (GridBoxExtents * 2) / CellSize;
            int3 dataSize{ int(gridSize.x), int(gridSize.y), int(gridSize.z) };

            int cell_x_min = std::clamp(static_cast<int>((pos.x - radius) / CellSize), 0, dataSize.x);
            int cell_x_max = std::clamp(static_cast<int>((pos.x + radius) / CellSize), 0, dataSize.x);
            int cell_y_min = std::clamp(static_cast<int>((pos.y - radius) / CellSize), 0, dataSize.y);
            int cell_y_max = std::clamp(static_cast<int>((pos.y + radius) / CellSize), 0, dataSize.y);
            int cell_z_min = std::clamp(static_cast<int>((pos.z - radius) / CellSize), 0, dataSize.z);
            int cell_z_max = std::clamp(static_cast<int>((pos.z + radius) / CellSize), 0, dataSize.z);

            for (int cell_x = cell_x_min; cell_x <= cell_x_max; ++cell_x)
            {
                for (int cell_y = cell_y_min; cell_y <= cell_y_max; ++cell_y)
                {
                    for (int cell_z = cell_z_min; cell_z <= cell_z_max; ++cell_z)
                    {
                        CellToObjectsMap[uint3{ uint(cell_x), uint(cell_y), uint(cell_z) }].push_back(&p);
                    }
                }
            }
#endif
        }

        // Get all particles within a given radius of a point
        // The output particles are stored in the provided buffer
        // Returns the number of particles in the buffer
        size_t GetNearest(const RPhysicsObject& object, RDynamicVector<RPhysicsObject*>& out_particles) const
        {
            out_particles.clear();
            size_t num_particles = 0;

            const auto& point = object.CurPosition;
            float radius = PointSize;

#if GRID_USE_HASH_TABLE
            int min_cell_x = static_cast<int>((point.x - radius) / CellSize);
            int max_cell_x = static_cast<int>((point.x + radius) / CellSize);
            int min_cell_y = static_cast<int>((point.y - radius) / CellSize);
            int max_cell_y = static_cast<int>((point.y + radius) / CellSize);
            int min_cell_z = static_cast<int>((point.z - radius) / CellSize);
            int max_cell_z = static_cast<int>((point.z + radius) / CellSize);

            for (int cell_x = min_cell_x; cell_x <= max_cell_x; ++cell_x)
            {
                for (int cell_y = min_cell_y; cell_y <= max_cell_y; ++cell_y)
                {
                    for (int cell_z = min_cell_z; cell_z <= max_cell_z; ++cell_z)
                    {
                        int cell_id = GetCelllIndex(cell_x, cell_y, cell_z);
                        const auto& particles_in_cell = CellToObjectsMap.find(cell_id);
                        if (particles_in_cell != CellToObjectsMap.end())
                        {
                            for (RPhysicsObject* p : particles_in_cell->second)
                            {
                                if(p != &object && p->LastCollisionObject != &object && object.LastCollisionObject != p)
                                {
                                    out_particles.push_back(p);
                                    p->LastCollisionObject = &object;
                                    num_particles++;
                                }
                            }
                        }
                    }
                }
            }
#else
            float3 pos = point + GridBoxExtents;
            float3 gridSize = (GridBoxExtents * 2) / CellSize;
            int3 dataSize{ int(gridSize.x), int(gridSize.y), int(gridSize.z) };

            int3 centerId{ int(object.IndexInSpatialGrid.x), int(object.IndexInSpatialGrid.y), int(object.IndexInSpatialGrid.z) };

            for (int x = std::max(centerId.x - 1, 0); x < std::min(centerId.x + 2, dataSize.x); x++)
            {
                for (int y = std::max(centerId.y - 1, 0); y < std::min(centerId.y + 2, dataSize.y); y++)
                {
                    for (int z = std::max(centerId.z - 1, 0); z < std::min(centerId.z + 2, dataSize.z); z++)
                    {
                        uint3 coord{ x, y, z};
                        auto& vec = CellToObjectsMap[coord];
                        for (RPhysicsObject* p : vec)
                        {
                            if (p != &object && p->LastCollisionObject != &object && object.LastCollisionObject != p)
                            {
                                out_particles.push_back(p);
                                p->LastCollisionObject = &object;
                                num_particles++;
                            }
                        }
                    }
                }
            }
#endif
            return num_particles;
        }

        void ClearVectors(RDynamicVector<RPhysicsObject>& objects)
        {
                PROFILE_FUNC();
        #if GRID_USE_HASH_TABLE
                CellToObjectsMap.clear();
        #else
                for (auto& obj : objects)
                {
                    CellToObjectsMap[obj.IndexInSpatialGrid].clear();
                }
        #endif
        }

        float CellSize = 1.f;
#if GRID_USE_HASH_TABLE
        RHashTable<int, RStackVector<RPhysicsObject*, 5>> CellToObjectsMap;
#else
        RVector3DLinear<RDynamicVector<RPhysicsObject*>> CellToObjectsMap;
#endif
        // Get the hash id for a cell given its (x, y, z) coordinates
        int GetCelllIndex(int x, int y, int z) const
        {
            // Use Cantor pairing function to combine the three coordinates into a single integer
            return ((x + y) * (x + y + 1) / 2 + y) * (y + z + 1) + z;
        }

    };

#if GRID_USE_HASH_TABLE
    SpatialHashGrid HashGrid{ 1 };
#else
    SpatialHashGrid HashGrid{ 1, float3{20,20,20} };
#endif

};













///RHeatEquation Example
#if 0

static RScalarFunctionVisualizer2D Visualizer;

if (bFirstBoot)
{
    Visualizer.InputMin = float2{ -2.5,-2.5 };
    Visualizer.InputMax = float2{ 2.5,2.5 };
}

auto temperatureFunc = [&](float2 pos)
{
    return pos.x + pos.y;
    return pow(e, -pow(pos.x * pos.y, 2.f));
};

RHeatEquation::DrawUI();

//Visualizer.Main(temperatureFunc, GImGeometry);

Visualizer.DrawUI();

if (bFirstBoot || RHeatEquation::Reboot)
{
    std::function f = temperatureFunc;
    Visualizer.SampleFunction(f);

    RHeatEquation::Reboot = false;
}

RHeatEquation::Process(Visualizer, deltaT);

Visualizer.Visualize(GImGeometry);


#endif

class RHeatEquation
{
public:

    static inline bool Run = false;
    static inline bool RunOnce = false;
    static inline bool bReadOutOfBounds = false;
    static inline bool Reboot = false;

    static inline int numSteps = 1;

    static inline float aConst = 0.001f;

    static void Process(RScalarFunctionVisualizer2D& function, float dt)
    {

        check(function.Function != nullptr); ///Function should've been sampled

        DrawUI();

        for (int i = 0; i < numSteps; i++)
        {
            if (Run || RunOnce)
            {
                for (uint y = 0; y < function.NumSamples.y; y++)
                {
                    for (uint x = 0; x < function.NumSamples.x; x++)
                    {
                        uint2 coord{ x,y };

                        auto curOut = function.OutputSet[coord];

                        auto divCur = function.LaplacianSet[coord] * function.SamplingPeriodLength.x * function.SamplingPeriodLength.y;

                        auto velocity = aConst * divCur * dt;

                        function.OutputSet[coord] = curOut + velocity;
                        //function.OutputSet[coord] = curOut * 0.99;
                    }
                }

                RunOnce = false;

                function.GenerateDerivatives();

            }

        }
    }


    static void DrawUI()
    {
        if (ImGui::Begin("Heat Equation"))
        {
            RunOnce = ImGui::Button("RunOnce");
            ImGui::Checkbox("Run", &Run);
            ImGui::Checkbox("Reboot", &Reboot);
            ImGui::Checkbox("bReadOutOfBounds", &bReadOutOfBounds);

            ImGui::DragFloat("aConst", &aConst, 0.0001, -10, 10);

            ImGui::SliderInt("NumSteps", &numSteps, 1, 50);
        }
        ImGui::End();
    }



};
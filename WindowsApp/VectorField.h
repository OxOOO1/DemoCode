#pragma once

#include "LinearAlgebra.h"
#include "GenerateSimpleGeometry.h"

class RVectorFieldVisualizer
{
public:

	RVector3D<float3> InputSet;
	RVector3D<float3> VectorField;

    ///Derivatives
    struct PartialDerivatives
    {
        float3 PartialX;
        float3 PartialY;
        float3 PartialZ;
    };
    RVector3D<PartialDerivatives> FirstDerivativeSet;
    RVector3D<PartialDerivatives> SecondDerivativeSet;

    RVector3D<float> Divergence;
    RVector3D<float3> Curl;///Positive value - clockwise rotation

	float3 SamplingPeriodLength{};
	uint3 NumSamples{ 4, 4, 4 };
	float3 DomainMin{ -1,-1,-1 };
	float3 DomainMax{ 1,1,1 };
    float DomainScale = 1.f;
    float SliceOffset = 0.f;

    ///For specific Sample point
    uint3 PointOnAFunc;

    std::function<float3(float3)> Function;

    ///Visualizer Settings
    float3 FlowLineStartPos{ 0,0,0 };
    float3 CurFlowArrowPos{ 0,0,0 };
    bool bVisualizeFlowLine = false;
    bool bFlowLinePaused = false;
    float FlowLineSegmentLength = 0.01f;
    uint FlowLineNumSegments = 100;

    bool b3D = true;

    enum class EColorVisualizeMode : int32
    {
        Magnitude,
        Divergence,
        Curl
    };
    EColorVisualizeMode CurColorVisMode{ EColorVisualizeMode::Magnitude };


    void VisualizeFunction(std::function<float3(float3)> func, RGeometry& visualizer)
    {
        Function = func;

        DrawUI();

        SampleFunction(func);

        Visualize(visualizer);

        if (bVisualizeFlowLine)
        {
            VisualizeFlowLine(visualizer);
        }

        if (bSpanParticles)
        {
            SpanParticles();

            bSpanParticles = false;
        }

        if (bRemoveParticles)
        {
            Particles.clear();
            bRemoveParticles = false;
        }


        if (!Particles.empty())
        {
            VisualizeParticlesSimulation(visualizer);
        }

    }


    void VisualizeGradientField(RVector2D<float3> gradients, RGeometry& visualizer)
    {
        VectorField.Resize(uint3{ gradients.Size.x, gradients.Size.y, 1 });

        for (uint y = 0; y < NumSamples.y; y++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                uint2 coord{ x,y };
                const auto& v = gradients[coord];
                VectorField[uint3{x,y,0}] = float3{ v.x, v.y, 0 };
            }
        }

        Visualize(visualizer);

    }


private:

    void SampleFunction(std::function<float3(float3)> func)
    {
        auto DomainMinScaled = DomainMin * DomainScale;
        auto DomainMaxScaled = DomainMax * DomainScale;

        float3 DomainLength = DomainMaxScaled - DomainMinScaled;

        SamplingPeriodLength.x = DomainLength.x / (NumSamples.x - 1);
        SamplingPeriodLength.y = DomainLength.y / (NumSamples.y - 1);

        if (b3D)
        {
            SamplingPeriodLength.z = DomainLength.z / (NumSamples.z - 1);
        }
        else
        {
            SamplingPeriodLength.z = 0;
            DomainMinScaled.z = SliceOffset;
        }


        InputSet.Resize(NumSamples);
        VectorField.Resize(NumSamples);


        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    uint3 sampleCoord{ x,y,z };
                    float3 uvw;
                    uvw.x = DomainMinScaled.x + float(x) * SamplingPeriodLength.x;
                    uvw.y = DomainMinScaled.y + float(y) * SamplingPeriodLength.y;
                    uvw.z = DomainMinScaled.z + float(z) * SamplingPeriodLength.z;

                    VectorField[sampleCoord] = func(uvw);

                    InputSet[sampleCoord] = uvw;
                }
            }
        }


        GenerateDerivatives();

    }

    void GenerateDerivatives()
    {
        FirstDerivativeSet.Resize(VectorField.Size);

        Divergence.Resize(VectorField.Size);
        Curl.Resize(VectorField.Size);

        uint numSlices = b3D ? (NumSamples.z - 1) : 1;

        for (uint z = 0; z < numSlices; z++)
        {
            for (uint y = 0; y < NumSamples.y - 1; y++)
            {
                for (uint x = 0; x < NumSamples.x - 1; x++)
                {
                    uint3 coord{ x,y,z };

                    c_float3& vec = VectorField[coord];
                    c_float3& vecNextX = VectorField[coord + uint3{ 1,0,0 }];
                    c_float3& vecNextY = VectorField[coord + uint3{ 0,1,0 }];
                    
                    c_float3 partialX = (vecNextX - vec) / SamplingPeriodLength.x;
                    c_float3 partialY = (vecNextY - vec) / SamplingPeriodLength.y;

                    FirstDerivativeSet[coord].PartialX = partialX;
                    FirstDerivativeSet[coord].PartialY = partialY;
                    

                    if (b3D)
                    {
                        c_float3& vecNextZ = VectorField[coord + uint3{ 0,0,1 }];
                        c_float3 partialZ = (vecNextZ - vec) / SamplingPeriodLength.z;
                        FirstDerivativeSet[coord].PartialZ = partialZ;

                        Divergence[coord] = partialX.x + partialY.y + partialZ.z;
                        Curl[coord] = float3{partialY.z - partialZ.y, partialZ.x - partialX.z, partialX.y - partialY.x};
                    }
                    else
                    {
                        Divergence[coord] = partialX.x + partialY.y;
                        Curl[coord] = float3{0,0,1} * (partialX.y - partialY.x);
                    }
                }
            }
        }

        ///Copy last elements
        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint x = 0; x < NumSamples.x; x++)
            {
                uint y = NumSamples.y - 1;

                uint3 coord{ x,y,z };

                FirstDerivativeSet[coord].PartialX = FirstDerivativeSet[coord - uint3{ 0,1,0 }].PartialX;
                FirstDerivativeSet[coord].PartialY = FirstDerivativeSet[coord - uint3{ 0,1,0 }].PartialY;
                Divergence[coord] = Divergence[coord - uint3{ 0,1,0 }];
                Curl[coord] = Curl[coord - uint3{ 0,1,0 }];
            }
            for (uint y = 0; y < NumSamples.x; y++)
            {
                uint x = NumSamples.x - 1;

                uint3 coord{ x,y,z };

                FirstDerivativeSet[coord].PartialX = FirstDerivativeSet[coord - uint3{ 1,0,0 }].PartialX;
                FirstDerivativeSet[coord].PartialY = FirstDerivativeSet[coord - uint3{ 1,0,0 }].PartialY;
                Divergence[coord] = Divergence[coord - uint3{ 1,0,0 }];
                Curl[coord] = Curl[coord - uint3{ 1,0,0 }];
            }
        }
        
        

    }


    void DrawUI()
    {
        if (ImGui::Begin("Vector Field Visualizer"))
        {

            uint2 logBase2Size;
            logBase2Size.x = RMath::Log2(NumSamples.x);
            if (ImGui::SliderInt("NumSamplesPow2", (int*)&logBase2Size.x, 1, 11))
            {
                NumSamples.x = (uint)pow(2u, logBase2Size.x);
                NumSamples.y = NumSamples.x;
            }

            ImGui::Text("NumSamplesAll: %i", NumSamples.x * NumSamples.y * NumSamples.z);

            ImGui::Checkbox("3D", &b3D);
            if (!b3D)
            {
                NumSamples.z = 1;
                ImGui::DragFloat("SliceStart", &SliceOffset, 0.1, -100, 100);
            }
            else
            {
                NumSamples.z = NumSamples.x;
            }

            ImGui::DragFloat("DomainScale", &DomainScale, 0.1, 0.1, 2);
            if (ImGui::DragFloat("DomainMin", &DomainMin.x, 0.1, -100, 0))
            {
                DomainMin.y = DomainMin.z = DomainMin.x;
            }
            if (ImGui::DragFloat("DomainMax", &DomainMax.x, 0.1, 1, 100))
            {
                DomainMax.y = DomainMax.z = DomainMax.x;
            }


            ImGui::Text("SamplingPeriodLength: %f, %f", SamplingPeriodLength.x, SamplingPeriodLength.y);

            ImGui::TextColored(ImVec4(0.8, 0.8, 0.2, 1), "VisualizerSettings:");
            ImGui::RadioButton("Magnitude", (int*)&CurColorVisMode, (int)EColorVisualizeMode::Magnitude);
            ImGui::SameLine();
            ImGui::RadioButton("Divergence", (int*)&CurColorVisMode, (int)EColorVisualizeMode::Divergence);
            ImGui::SameLine();
            ImGui::RadioButton("Curl", (int*)&CurColorVisMode, (int)EColorVisualizeMode::Curl);

            ImGui::NewLine();
            ImGui::Checkbox("VisualizeFlowLine", &bVisualizeFlowLine);
            if (bVisualizeFlowLine)
            {
                if (ImGui::DragFloat3("FlowStart", &FlowLineStartPos.x, 0.1, DomainMin.x, DomainMax.x))
                {
                    CurFlowArrowPos = FlowLineStartPos;
                }
                ImGui::DragFloat("FlowLineSegmentLength", &FlowLineSegmentLength, 0.01, 0.01, 1);
                ImGui::DragInt("FlowLineNumSegments", (int*)&FlowLineNumSegments, 1, 1, 1000);
                ImGui::NewLine();

                if (ImGui::Button("RestartAnimation"))
                {
                    CurFlowArrowPos = FlowLineStartPos;
                }
                ImGui::Checkbox("PauseAnimation", &bFlowLinePaused);

                ImGui::NewLine();
            }


            ImGui::TextColored(ImVec4(1, 1, 0, 1), "PointOnAFunc");
            ImGui::SliderInt3("PointOnAFunc", (int*)&PointOnAFunc.x, 0, NumSamples.x - 1);
            if (VectorField.Size.x > 0)
            {
                float CurCurl = VectorGetLength(Curl[PointOnAFunc]);
                float CurDivergence = Divergence[PointOnAFunc];
                float CurMagnitude = VectorGetLength(VectorField[PointOnAFunc]);

                float3 partialX = (FirstDerivativeSet[PointOnAFunc].PartialX);
                ImGui::Text("PartialX: %f, %f, %f", partialX.x, partialX.y, partialX.z);

                float3 partialY = (FirstDerivativeSet[PointOnAFunc].PartialY);
                ImGui::Text("PartialY: %f, %f, %f", partialY.x, partialY.y, partialY.z);
                if (b3D)
                {
                    float3 partialZ = (FirstDerivativeSet[PointOnAFunc].PartialZ);
                    ImGui::Text("PartialZ: %f, %f, %f", partialZ.x, partialZ.y, partialZ.z);
                }


                ImGui::Text("Magnitude: %f, Divergence: %f, Curl: %f", CurMagnitude, CurDivergence, CurCurl);

                ImGui::NewLine();
            }

            ImGui::NewLine();
            ImGui::TextColored(ImVec4(1, 1, 0, 1), "ParticlesSimulation");
            if (ImGui::Button("Span Particles"))
            {
                bSpanParticles = true;
            }
            if (ImGui::Button("Destroy Particles"))
            {
                bRemoveParticles = true;
            }
            if (!Particles.empty())
            {
                ImGui::Checkbox("ParticleSimulationRunning", &bSimulateParticles);
                ImGui::Checkbox("SimulateIn2D", &bParticlesSimulation2D);
                ImGui::DragFloat("ParticleSimulationSpeed", &ParticleSimulationSpeed, 0.01, 0.01, 10);
                ImGui::SliderFloat("ParticleSizeScale", &ParticleSizeScale, 0.01, 2);
                if (ImGui::Button("Re-Import Particles"))
                {
                    bImportParticles = true;
                }
            }
        }
        ImGui::End();

    }


    void Visualize(RGeometry& visualizer);

    

    void VisualizeFlowLine(RGeometry& visualizer)
    {
        visualizer.SetColor(EColor::LightBlue);
        visualizer.AddCircleCameraFacing(FlowLineStartPos, VectorGetLength(SamplingPeriodLength) * 0.1);

        float3 curPos = FlowLineStartPos;

        visualizer.SetColor(EColor::Green);

        ///Generate flow line from a point
        for (int i = 0; i < FlowLineNumSegments; i++)
        {
            float3 prevPos = curPos;
            curPos = curPos + Function(curPos) * FlowLineSegmentLength;

            visualizer.AddLine(prevPos, curPos);

        }


        ///Animated arrow through a line
        static float3 positionCached;
        static float3 direction;
        
        if (!bFlowLinePaused)
        {
            positionCached = CurFlowArrowPos;

            auto dt = Time::Get().GetDeltaTimeSec();
            direction = Function(CurFlowArrowPos) * dt;
            CurFlowArrowPos = CurFlowArrowPos + direction;
        }
        

        visualizer.SetColor(EColor::LightBlue);
        visualizer.AddCircleCameraFacing(positionCached, 0.1);
        visualizer.SetColor(EColor::White);
        visualizer.AddArrowPolyboard(positionCached, positionCached + direction);
        
    }


public:

    RDynamicVector<float3> Particles;
    bool bSimulateParticles = false;
    float ParticleSizeScale = 0.5f;
    float ParticleSimulationSpeed = 1.f;
    bool bImportParticles = true;
    bool bSpanParticles = false;
    bool bRemoveParticles = false;
    bool bParticlesSimulation2D = true;

    void SpanParticles()
    {
        ///Populate domain with particles

        Particles.resize(NumSamples.x * NumSamples.y * NumSamples.z);

        uint i = 0;

        for (uint z = 0; z < NumSamples.z; z++)
        {
            for (uint y = 0; y < NumSamples.y; y++)
            {
                for (uint x = 0; x < NumSamples.x; x++)
                {
                    Particles[i] = InputSet[uint3{x,y,z}];
                    i++;
                }
            }
        }
    }

    void ImportParticlesIntoField(const RDynamicVector<float3>& particlesPositions)
    {
        if (bImportParticles)
        {
            Particles = particlesPositions;
        }
        
        bImportParticles = false;
    }

    void VisualizeParticlesSimulation(RGeometry& visualizer)
    {
        float dt = 1.f;

        if (bSimulateParticles)
        {
            dt = Time::Get().GetDeltaTimeSec();
        }
        

        {
            for (auto& particlePos : Particles)
            {
                auto curDirection = Function(particlePos) * dt * ParticleSimulationSpeed;
                if (bParticlesSimulation2D)
                {
                    curDirection.z = 0.f;
                }
                auto newParticlePos = particlePos + curDirection;

                visualizer.SetColor(EColor::White);
                visualizer.AddLine(particlePos, newParticlePos, ParticleSizeScale);

                if (bSimulateParticles)
                {
                    particlePos = newParticlePos;
                }

            }
        }


        
    }

private:

};

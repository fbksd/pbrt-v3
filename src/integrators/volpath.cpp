
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// integrators/volpath.cpp*
#include "integrators/volpath.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"
#include <fbksd/renderer/samples.h>

namespace pbrt {

STAT_INT_DISTRIBUTION("Integrator/Path length", pathLength);
STAT_COUNTER("Integrator/Volume interactions", volumeInteractions);
STAT_COUNTER("Integrator/Surface interactions", surfaceInteractions);

// VolPathIntegrator Method Definitions
void VolPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum VolPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth, SampleBuffer* sampleBuffer, Spectrum* diffuse) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f), beta(1.f);
    RayDifferential ray(r);
    bool specularBounce = false;
    bool hadNonSpecularBounce = false;
    int bounces;
    // Added after book publication: etaScale tracks the accumulated effect
    // of radiance scaling due to rays passing through refractive
    // boundaries (see the derivation on p. 527 of the third edition). We
    // track this value in order to remove it from beta when we apply
    // Russian roulette; this is worthwhile, since it lets us sometimes
    // avoid terminating refracted rays that are about to be refracted back
    // out of a medium and thus have their beta value increased.
    Float etaScale = 1;

    for (bounces = 0;; ++bounces) {
        // Intersect _ray_ with scene and store intersection in _isect_
        SurfaceInteraction isect;
        bool foundIntersection = scene.Intersect(ray, &isect);

        // Sample the participating medium, if present
        MediumInteraction mi;
        if (ray.medium) beta *= ray.medium->Sample(ray, sampler, arena, &mi);
        if (beta.IsBlack()) break;

        // Handle an interaction with a medium or a surface
        Spectrum directL;
        if (mi.IsValid()) {
            // Terminate path if ray escaped or _maxDepth_ was reached
            if (bounces >= maxDepth) break;

            if(bounces <= 1)
            {
                const Point3f &p = mi.p;
                const Normal3f &n = mi.n;
                sampleBuffer->set(WORLD_X, bounces, p.x);
                sampleBuffer->set(WORLD_Y, bounces, p.y);
                sampleBuffer->set(WORLD_Z, bounces, p.z);
                Normal3f nn = Faceforward(n, -ray.d);
                sampleBuffer->set(NORMAL_X, bounces, nn.x);
                sampleBuffer->set(NORMAL_Y, bounces, nn.y);
                sampleBuffer->set(NORMAL_Z, bounces, nn.z);

                Spectrum tex = beta;
                float rgb[3];
                tex.ToRGB(rgb);
                sampleBuffer->set(TEXTURE_COLOR_R, bounces, rgb[0]);
                sampleBuffer->set(TEXTURE_COLOR_G, bounces, rgb[1]);
                sampleBuffer->set(TEXTURE_COLOR_B, bounces, rgb[2]);
            }
            if(bounces == 0)
            {
                // copy the depth value because the benchmark needs it
                r.tMax = ray.tMax;
            }

            ++volumeInteractions;
            // Handle scattering at point in medium for volumetric path tracer
            const Distribution1D *lightDistrib =
                lightDistribution->Lookup(mi.p);
            Spectrum diffComp;
            L += beta * UniformSampleOneLight(mi, scene, arena, sampler, directL, bounces, sampleBuffer, diffComp,
                                              true, lightDistrib);

            Vector3f wo = -ray.d, wi;
            mi.phase->Sample_p(wo, &wi, sampler.Get2D());
            ray = mi.SpawnRay(wi);
            specularBounce = false;
        } else {
            ++surfaceInteractions;
            // Handle scattering at point on surface for volumetric path tracer

            // Possibly add emitted light at intersection
            if (bounces == 0 || specularBounce) {
                // Add emitted light at path vertex or from the environment
                if (foundIntersection)
                    L += beta * isect.Le(-ray.d);
                else
                    for (const auto &light : scene.infiniteLights)
                        L += beta * light->Le(ray);
            }

            // Terminate path if ray escaped or _maxDepth_ was reached
            if (!foundIntersection || bounces >= maxDepth) break;

            // Compute scattering functions and skip over medium boundaries
            isect.ComputeScatteringFunctions(ray, arena, true);
            if (!isect.bsdf) {
                ray = isect.SpawnRay(ray.d);
                bounces--;
                continue;
            }

            if(bounces <= 1)
            {
                const Point3f &p = isect.p;
                const Normal3f &n = isect.shading.n;
                sampleBuffer->set(WORLD_X, bounces, p.x);
                sampleBuffer->set(WORLD_Y, bounces, p.y);
                sampleBuffer->set(WORLD_Z, bounces, p.z);
                Normal3f nn = Faceforward(n, -ray.d);
                sampleBuffer->set(NORMAL_X, bounces, nn.x);
                sampleBuffer->set(NORMAL_Y, bounces, nn.y);
                sampleBuffer->set(NORMAL_Z, bounces, nn.z);

                Spectrum tex = isect.bsdf->getAlbedo();
                float rgb[3];
                tex.ToRGB(rgb);
                sampleBuffer->set(TEXTURE_COLOR_R, bounces, rgb[0]);
                sampleBuffer->set(TEXTURE_COLOR_G, bounces, rgb[1]);
                sampleBuffer->set(TEXTURE_COLOR_B, bounces, rgb[2]);
            }
            if(bounces == 0)
            {
                // copy the depth value because the benchmark needs it
                r.tMax = ray.tMax;
            }

            // Sample illumination from lights to find attenuated path
            // contribution
            const Distribution1D *lightDistrib =
                lightDistribution->Lookup(isect.p);
            Spectrum diffComp;
            L += beta * UniformSampleOneLight(isect, scene, arena, sampler, directL, bounces, sampleBuffer, diffComp,
                                              true, lightDistrib);

            if(bounces == 0)
            {
                float rgb[] = {0.f, 0.f, 0.f};
                directL.ToRGB(rgb);
                sampleBuffer->set(DIRECT_LIGHT_R, rgb[0]);
                sampleBuffer->set(DIRECT_LIGHT_G, rgb[1]);
                sampleBuffer->set(DIRECT_LIGHT_B, rgb[2]);
            }

            // Sample BSDF to get new path direction
            Vector3f wo = -ray.d, wi;
            Float pdf;
            BxDFType flags;
            Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                              BSDF_ALL, &flags);
            if (f.IsBlack() || pdf == 0.f) break;
            beta *= f * AbsDot(wi, isect.shading.n) / pdf;
            DCHECK(std::isinf(beta.y()) == false);
            specularBounce = (flags & BSDF_SPECULAR) != 0;
            if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
                Float eta = isect.bsdf->eta;
                // Update the term that tracks radiance scaling for refraction
                // depending on whether the ray is entering or leaving the
                // medium.
                etaScale *=
                    (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
            }
            ray = isect.SpawnRay(wi);

            if(!hadNonSpecularBounce && !specularBounce)
            {
                hadNonSpecularBounce = true;

                const Point3f &p = isect.p;
                const Normal3f &n = isect.shading.n;
                sampleBuffer->set(WORLD_X_NS, p.x);
                sampleBuffer->set(WORLD_Y_NS, p.y);
                sampleBuffer->set(WORLD_Z_NS, p.z);
                Normal3f nn = Faceforward(n, wo);
                sampleBuffer->set(NORMAL_X_NS, nn.x);
                sampleBuffer->set(NORMAL_Y_NS, nn.y);
                sampleBuffer->set(NORMAL_Z_NS, nn.z);
                float rgb[3];
                isect.bsdf->getAlbedo().ToRGB(rgb);
                sampleBuffer->set(TEXTURE_COLOR_R_NS, rgb[0]);
                sampleBuffer->set(TEXTURE_COLOR_G_NS, rgb[1]);
                sampleBuffer->set(TEXTURE_COLOR_B_NS, rgb[2]);
            }

            // Account for attenuated subsurface scattering, if applicable
            if (isect.bssrdf && (flags & BSDF_TRANSMISSION)) {
                // Importance sample the BSSRDF
                SurfaceInteraction pi;
                Spectrum S = isect.bssrdf->Sample_S(
                    scene, sampler.Get1D(), sampler.Get2D(), arena, &pi, &pdf);
                DCHECK(std::isinf(beta.y()) == false);
                if (S.IsBlack() || pdf == 0) break;
                beta *= S / pdf;

                // Account for the attenuated direct subsurface scattering
                // component
                Spectrum diffComp;
                L += beta *
                     UniformSampleOneLight(pi, scene, arena, sampler, directL, 0, nullptr, diffComp, true,
                                           lightDistribution->Lookup(pi.p));

                // Account for the indirect subsurface scattering component
                Spectrum f = pi.bsdf->Sample_f(pi.wo, &wi, sampler.Get2D(),
                                               &pdf, BSDF_ALL, &flags);
                if (f.IsBlack() || pdf == 0) break;
                beta *= f * AbsDot(wi, pi.shading.n) / pdf;
                DCHECK(std::isinf(beta.y()) == false);
                specularBounce = (flags & BSDF_SPECULAR) != 0;
                ray = pi.SpawnRay(wi);
            }
        }

        // Possibly terminate the path with Russian roulette
        // Factor out radiance scaling due to refraction in rrBeta.
        Spectrum rrBeta = beta * etaScale;
        if (rrBeta.MaxComponentValue() < rrThreshold && bounces > 3) {
            Float q = std::max((Float).05, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q) break;
            beta /= 1 - q;
            DCHECK(std::isinf(beta.y()) == false);
        }
    }
    ReportValue(pathLength, bounces);
    return L;
}

VolPathIntegrator *CreateVolPathIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");
    return new VolPathIntegrator(maxDepth, camera, sampler, pixelBounds,
                                 rrThreshold, lightStrategy);
}

}  // namespace pbrt

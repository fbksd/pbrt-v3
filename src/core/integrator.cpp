
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

// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "stats.h"
#include "samplers/random.h"
#include "samplers/zerotwosequence.h"


namespace pbrt {


class SparsePixelSampler
{
public:
    SparsePixelSampler(int beginx, int endx, int beginy, int endy, int n, bool isSPP)
    {
        m_beginX = beginx;
        m_endX = endx;
        m_beginY = beginy;
        m_endY = endy;
        m_width = m_endX - m_beginX;
        m_height = m_endY - m_beginY;
        m_xPos = m_beginX;
        m_yPos = m_beginY;
        m_sparseSampleIndex = 0;

        if(isSPP)
        {
            m_spp = n;
            m_sparseSampleCount = 0;
        }
        else
        {
            m_spp = n / (float)(m_width*m_height);
            m_sparseSampleCount = n % (m_width*m_height);
        }

        if(m_sparseSampleCount)
            m_stage = SPARSE;
        if(m_spp)
            m_stage = SPP;
    }

    virtual bool nextPixel(Point2i* pos)
    {
        int x = m_xPos, y = m_yPos;

        switch(m_stage)
        {
            case SPP:
                if(m_yPos == m_endY || m_xPos == m_endX)
                {
                    m_stage = SPARSE;
                    goto SPARSE_CASE;
                }

                if(++m_xPos == m_endX)
                {
                    m_xPos = m_beginX;
                    ++m_yPos;
                }
                break;
            case SPARSE:
                SPARSE_CASE:
                if(m_sparseSampleIndex++ == m_sparseSampleCount)
                    return false;
                x = m_beginX + m_random.UniformUInt32(m_width);
                y = m_beginY + m_random.UniformUInt32(m_height);
        }

        pos->x = x;
        pos->y = y;
        return true;
    }

protected:
    enum Stage
    {
        SPP,
        SPARSE
    };

    int m_beginX;
    int m_endX;
    int m_beginY;
    int m_endY;
    int m_width;
    int m_height;
    int m_xPos, m_yPos;
    int m_spp;
    int m_sparseSampleCount;
    int m_sparseSampleIndex;
    Stage m_stage;
    RNG m_random;
};

// Sampler used when the # of samples is not in SPP.
class SparseSampler: public Sampler
{
public:
    // New interface
    SparseSampler(int ns, int seed = 0):
        Sampler(ns)
    {
        m_numSamples = ns;
        if(RoundUpPow2(m_numSamples) == m_numSamples)
            m_sampler = new ZeroTwoSequenceSampler(m_numSamples);
        else
            m_sampler = new RandomSampler(m_numSamples, seed);
    }

    void StartPixel(const Point2i &p)
    {
        m_sampler->StartPixel(p);
    }

    virtual Float Get1D()
    {
        return m_sampler->Get1D();
    }

    virtual Point2f Get2D()
    {
        return m_sampler->Get2D();
    }

    virtual std::unique_ptr<Sampler> Clone(int seed)
    {
        return m_sampler->Clone(seed);
    }


    // Old interface
    /*
    SparseSampler(int xstart, int xend, int ystart, int yend, int ns, float sopen, float sclose):
        Sampler(xstart, xend, ystart, yend, ns, sopen, sclose)
    {
        m_numSamples = ns;
        if(isPerfectSquare(m_numSamples))
        {
            int ss = std::sqrt(m_numSamples);
            m_sampler = new StratifiedSampler(0, 1, 0, 1, ss, ss, true, sopen, sclose);
        }
        else
        {
            int xs = m_numSamples;
            m_sampler = new StratifiedSampler(0, 1, 0, 1, xs, 1, true, sopen, sclose);
        }
    }

    ~SparseSampler()
    { delete m_sampler; }

    int MaximumSampleCount() { return m_sampler->MaximumSampleCount(); }

    int GetMoreSamples(Sample *sample, RNG &rng)
    {
        int result = m_sampler->GetMoreSamples(sample, rng);
        for(size_t i = 0; i < result; ++i)
        {
            sample[i].imageX = Lerp(sample[i].imageX, xPixelStart, xPixelEnd);
            sample[i].imageY = Lerp(sample[i].imageY, yPixelStart, yPixelEnd);
        }

        return result;
    }

    int RoundSize(int sz) const { return sz; }

    Sampler *GetSubSampler(int num, int count)
    {
        int sps = getSubSamplerSize(num, count);
        if(sps)
            return new SparseSampler(xPixelStart, xPixelEnd, yPixelStart, yPixelEnd, sps, shutterOpen, shutterClose);
        else
            return nullptr;
    }

    int getSubSamplerSize(int num, int count) const
    {
        int sps = m_numSamples / count;
        int rest = m_numSamples % count;
        sps += num == 0 ? rest : 0;
        return sps;
    }
    */

private:
    Sampler* m_sampler;
    int m_numSamples; // this is the plain number of samples, not in spp.
};

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

// Integrator Method Definitions
Integrator::~Integrator() {}

// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
                                MemoryArena &arena, Sampler &sampler,
                                const std::vector<int> &nLightSamples,
                                bool handleMedia) {
    ProfilePhase p(Prof::DirectLighting);
    Spectrum L(0.f);
    Spectrum directL;
    for (size_t j = 0; j < scene.lights.size(); ++j) {
        // Accumulate contribution of _j_th light to _L_
        const std::shared_ptr<Light> &light = scene.lights[j];
        int nSamples = nLightSamples[j];
        const Point2f *uLightArray = sampler.Get2DArray(nSamples);
        const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
        if (!uLightArray || !uScatteringArray) {
            // Use a single sample for illumination from _light_
            Point2f uLight = sampler.Get2D();
            Point2f uScattering = sampler.Get2D();
            L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
                                arena, directL, handleMedia);
        } else {
            // Estimate direct lighting using sample arrays
            Spectrum Ld(0.f);
            for (int k = 0; k < nSamples; ++k)
                Ld += EstimateDirect(it, uScatteringArray[k], *light,
                                     uLightArray[k], scene, sampler, arena, directL,
                                     handleMedia);
            L += Ld / nSamples;
        }
    }
    return L;
}

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
                               MemoryArena &arena, Sampler &sampler, Spectrum& directL, int bounce,
                               SampleBuffer* sampleBuffer, bool handleMedia, const Distribution1D *lightDistrib) {
    ProfilePhase p(Prof::DirectLighting);
    // Randomly choose a single light to sample, _light_
    int nLights = int(scene.lights.size());
    if (nLights == 0) return Spectrum(0.f);
    int lightNum;
    Float lightPdf;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
        lightPdf = Float(1) / nLights;
    }
    const std::shared_ptr<Light> &light = scene.lights[lightNum];
    Point2f uLight = sampler.Get2D();
    Point2f uScattering = sampler.Get2D();
    if(bounce == 0 && sampleBuffer)
    {
        uLight.x = sampleBuffer->set(LIGHT_X, uLight.x);
        uLight.y = sampleBuffer->set(LIGHT_Y, uLight.y);
    }
    return EstimateDirect(it, uScattering, *light, uLight,
                          scene, sampler, arena, directL, handleMedia) / lightPdf;
}

Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, Spectrum& directL, bool handleMedia, bool specular) {
    BxDFType bsdfFlags =
        specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
    Spectrum Ld(0.f);
    // Sample light source with multiple importance sampling
    Vector3f wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
    VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
            << wi << ", pdf: " << lightPdf;
    if (lightPdf > 0 && !Li.IsBlack()) {
        // Compute BSDF or phase function's value for light sample
        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            // Evaluate BSDF for light sampling strategy
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
                AbsDot(wi, isect.shading.n);
            scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
            VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
        } else {
            // Evaluate phase function for light sampling strategy
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->p(mi.wo, wi);
            f = Spectrum(p);
            scatteringPdf = p;
            VLOG(2) << "  medium p: " << p;
        }
        if (!f.IsBlack()) {
            // Compute effect of visibility for light source sample
            if (handleMedia) {
                Li *= visibility.Tr(scene, sampler);
                VLOG(2) << "  after Tr, Li: " << Li;
            } else {
              if (!visibility.Unoccluded(scene)) {
                VLOG(2) << "  shadow ray blocked";
                Li = Spectrum(0.f);
              } else
                VLOG(2) << "  shadow ray unoccluded";
            }

            // Add light's contribution to reflected radiance
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags)) {
                    Ld += f * Li / lightPdf;
                    directL += Li / lightPdf;
                }
                else {
                    Float weight =
                        PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                    Ld += f * Li * weight / lightPdf;
                    directL += Li * weight / lightPdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!IsDeltaLight(light.flags)) {
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction()) {
            // Sample scattered direction for surface interactions
            BxDFType sampledType;
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
                                     bsdfFlags, &sampledType);
            f *= AbsDot(wi, isect.shading.n);
            sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
        } else {
            // Sample scattered direction for medium interactions
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
            f = Spectrum(p);
            scatteringPdf = p;
        }
        VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
            scatteringPdf;
        if (!f.IsBlack() && scatteringPdf > 0) {
            // Account for light contributions along sampled direction _wi_
            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);
                if (lightPdf == 0) return Ld;
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }

            // Find intersection and compute transmittance
            SurfaceInteraction lightIsect;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction =
                handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
                            : scene.Intersect(ray, &lightIsect);

            // Add light contribution from material sampling
            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIsect.primitive->GetAreaLight() == &light)
                    Li = lightIsect.Le(-wi);
            } else
                Li = light.Le(ray);
            if (!Li.IsBlack()) {
                Ld += f * Li * Tr * weight / scatteringPdf;
                directL += Li * Tr * weight / scatteringPdf;
            }
        }
    }
    return Ld;
}

std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
    const Scene &scene) {
    if (scene.lights.empty()) return nullptr;
    std::vector<Float> lightPower;
    for (const auto &light : scene.lights)
        lightPower.push_back(light->Power().y());
    return std::unique_ptr<Distribution1D>(
        new Distribution1D(&lightPower[0], lightPower.size()));
}

// SamplerIntegrator Method Definitions
/*
void SamplerIntegrator::Render(const Scene &scene) {
    Preprocess(scene, *sampler);
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);
    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            MemoryArena arena;

            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                {
                    ProfilePhase pp(Prof::StartPixel);
                    tileSampler->StartPixel(pixel);
                }

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;

                do {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // Generate camera ray for current sample
                    RayDifferential ray;
                    Float rayWeight =
                        camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(
                        1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                    ++nCameraRays;

                    // Evaluate radiance along camera ray
                    Spectrum L(0.f);
                    if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);

                    // Issue warning if unexpected radiance value returned
                    if (L.HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (L.y() < -1e-5) {
                        LOG(ERROR) << StringPrintf(
                            "Negative luminance value, %f, returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            L.y(), pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (std::isinf(L.y())) {
                          LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                        ray << " -> L = " << L;

                    // Add camera ray's contribution to image
                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
    camera->film->WriteImage();
}
*/

void SamplerIntegrator::Render(const Scene &scene)
{
    RenderingServer server;
    server.onSetParameters(
        [this](const SampleLayout& layout){
        m_layout = layout;
    });

    server.onGetSceneInfo([this, &scene](){
        SceneInfo info;
        this->getSceneInfo(scene, &info);
        return info;
    });

    server.onEvaluateSamples([this, &scene](int64_t spp, int64_t remainingCount){
        printf("Rendering ...\n");
        return this->evaluateSamples(scene, spp, remainingCount);
    });

    server.run(/*2227*/);
}

Spectrum SamplerIntegrator::SpecularReflect(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    // Compute specular reflection direction _wi_ and BSDF value
    Vector3f wo = isect.wo, wi;
    Float pdf;
    BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
    Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

    // Return contribution of specular reflection
    const Normal3f &ns = isect.shading.n;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential rd = isect.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = isect.p + isect.dpdx;
            rd.ryOrigin = isect.p + isect.dpdy;
            // Compute differential reflected directions
            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;
            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);
            rd.rxDirection =
                wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
            rd.ryDirection =
                wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
        }
        return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) /
               pdf;
    } else
        return Spectrum(0.f);
}

Spectrum SamplerIntegrator::SpecularTransmit(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    Vector3f wo = isect.wo, wi;
    Float pdf;
    const Point3f &p = isect.p;
    const BSDF &bsdf = *isect.bsdf;
    Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum L = Spectrum(0.f);
    Normal3f ns = isect.shading.n;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential rd = isect.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dpdx;
            rd.ryOrigin = p + isect.dpdy;

            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;

            // The BSDF stores the IOR of the interior of the object being
            // intersected.  Compute the relative IOR by first out by
            // assuming that the ray is entering the object.
            Float eta = 1 / bsdf.eta;
            if (Dot(wo, ns) < 0) {
                // If the ray isn't entering, then we need to invert the
                // relative IOR and negate the normal and its derivatives.
                eta = 1 / eta;
                ns = -ns;
                dndx = -dndx;
                dndy = -dndy;
            }

            /*
              Notes on the derivation:
              - pbrt computes the refracted ray as: \wi = -\eta \omega_o + [ \eta (\wo \cdot \N) - \cos \theta_t ] \N
                It flips the normal to lie in the same hemisphere as \wo, and then \eta is the relative IOR from
                \wo's medium to \wi's medium.
              - If we denote the term in brackets by \mu, then we have: \wi = -\eta \omega_o + \mu \N
              - Now let's take the partial derivative. (We'll use "d" for \partial in the following for brevity.)
                We get: -\eta d\omega_o / dx + \mu dN/dx + d\mu/dx N.
              - We have the values of all of these except for d\mu/dx (using bits from the derivation of specularly
                reflected ray deifferentials).
              - The first term of d\mu/dx is easy: \eta d(\wo \cdot N)/dx. We already have d(\wo \cdot N)/dx.
              - The second term takes a little more work. We have:
                 \cos \theta_i = \sqrt{1 - \eta^2 (1 - (\wo \cdot N)^2)}.
                 Starting from (\wo \cdot N)^2 and reading outward, we have \cos^2 \theta_o, then \sin^2 \theta_o,
                 then \sin^2 \theta_i (via Snell's law), then \cos^2 \theta_i and then \cos \theta_i.
              - Let's take the partial derivative of the sqrt expression. We get:
                1 / 2 * 1 / \cos \theta_i * d/dx (1 - \eta^2 (1 - (\wo \cdot N)^2)).
              - That partial derivatve is equal to:
                d/dx \eta^2 (\wo \cdot N)^2 = 2 \eta^2 (\wo \cdot N) d/dx (\wo \cdot N).
              - Plugging it in, we have d\mu/dx =
                \eta d(\wo \cdot N)/dx - (\eta^2 (\wo \cdot N) d/dx (\wo \cdot N))/(-\wi \cdot N).
             */
            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

            Float mu = eta * Dot(wo, ns) - AbsDot(wi, ns);
            Float dmudx =
                (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdx;
            Float dmudy =
                (eta - (eta * eta * Dot(wo, ns)) / AbsDot(wi, ns)) * dDNdy;

            rd.rxDirection =
                wi - eta * dwodx + Vector3f(mu * dndx + dmudx * ns);
            rd.ryDirection =
                wi - eta * dwody + Vector3f(mu * dndy + dmudy * ns);
        }
        L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
    }
    return L;
}


void SamplerIntegrator::getSceneInfo(const Scene &scene, SceneInfo *info)
{
    int width = camera->film->fullResolution.x;
    int height = camera->film->fullResolution.y;

    info->set<int>("width", width);
    info->set<int>("height", height);

    // NOTE: In PBRT, shutterOpen and shutterClose defaults are 0 and 1, respectively,
    // which makes them not good to decide if a scene has motion blur.
    // This info (`has_motion_blur`) is being set in the benchmark configuration file.
    info->set<float>("shutter_open", camera->shutterOpen);
    info->set<float>("shutter_close", camera->shutterClose);

    bool hasAreaLights = false;
    for(const auto& light : scene.lights)
    {
        if(light->flags & (int)LightFlags::Area)
        {
            hasAreaLights = true;
            break;
        }
    }
    info->set<bool>("has_area_lights", hasAreaLights);
}

bool SamplerIntegrator::evaluateSamples(const Scene& scene, int spp, int remaining)
{
    int64_t width = camera->film->fullResolution.x;
    int64_t height = camera->film->fullResolution.y;
    int64_t numPixels = width * height;
    int rest = remaining;
    float filmScale = camera->film->getScale();
    float maxSampleLuminance = camera->film->getMaxSampleLuminance();
    bool seekByPixel = !(m_layout.hasInput("IMAGE_X") || m_layout.hasInput("IMAGE_Y"));

    size_t offset = 0;
    if(spp)
    {
        std::unique_ptr<Sampler> sppSampler;
        sppSampler.reset(new RandomSampler(spp));
        if(seekByPixel)
            sppRender(scene, sppSampler.get(), width, height, filmScale, maxSampleLuminance);
        else
            offsetRender(scene, sppSampler.get(), width, height, filmScale, maxSampleLuminance);
        offset = spp * numPixels * m_layout.getSampleSize();
    }

    if(rest)
    {
        sparseRender(scene, rest, width, height, filmScale, maxSampleLuminance, offset);
    }

    LOG(INFO) << "Rendering finished";
}

void SamplerIntegrator::sppRender(const Scene& scene,
                                  Sampler* sppSampler,
                                  int width,
                                  int height,
                                  float filmScale,
                                  float maxSampleLuminance)
{
    Preprocess(scene, *sppSampler);
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    int spp = sppSampler->samplesPerPixel;
    Bounds2i sampleBounds(Point2i(0, 0), Point2i(width, height));
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);
    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            MemoryArena arena;

            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sppSampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            SamplesPipe pipe;
            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                {
                    ProfilePhase pp(Prof::StartPixel);
                    tileSampler->StartPixel(pixel);
                }

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;

                pipe.seek(pixel.x, pixel.y, spp, width);

                do {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // NOTE: I'm seeking by pixel only: not support for adaptive sampling
                    SampleBuffer sampleBuffer = pipe.getBuffer();
                    cameraSample.pFilm.x = sampleBuffer.set(IMAGE_X, cameraSample.pFilm.x);
                    cameraSample.pFilm.y = sampleBuffer.set(IMAGE_Y, cameraSample.pFilm.y);
                    cameraSample.pLens.x = sampleBuffer.set(LENS_U, cameraSample.pLens.x);
                    cameraSample.pLens.y = sampleBuffer.set(LENS_V, cameraSample.pLens.y);
                    cameraSample.time = sampleBuffer.set(TIME, cameraSample.time);

                    // Generate camera ray for current sample
                    RayDifferential ray;
                    Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                    ++nCameraRays;

                    // Evaluate radiance along camera ray
                    Spectrum L(0.f);
                    if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena, 0, &sampleBuffer);

                    // Issue warning if unexpected radiance value returned
                    if (L.HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (L.y() < -1e-5) {
                        LOG(ERROR) << StringPrintf(
                            "Negative luminance value, %f, returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            L.y(), pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (std::isinf(L.y())) {
                          LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                        ray << " -> L = " << L;

                    float lum = L.y();
                    if(lum > maxSampleLuminance)
                        L *= maxSampleLuminance / lum;
                    Float rgb[3] = {0.f, 0.f, 0.f};
                    L.ToRGB(rgb);
                    sampleBuffer.set(COLOR_R, rgb[0]*filmScale*rayWeight);
                    sampleBuffer.set(COLOR_G, rgb[1]*filmScale*rayWeight);
                    sampleBuffer.set(COLOR_B, rgb[2]*filmScale*rayWeight);
                    if(std::isfinite(ray.tMax))
                        sampleBuffer.set(DEPTH, ray.tMax);
                    else
                        sampleBuffer.set(DEPTH, 0.f);
                    pipe << sampleBuffer;

                    // Add camera ray's contribution to image
                    // filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
}

void SamplerIntegrator::offsetRender(const Scene& scene,
                                     Sampler* sppSampler,
                                     int width,
                                     int height,
                                     float filmScale,
                                     float maxSampleLuminance)
{
    Preprocess(scene, *sppSampler);
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    int spp = sppSampler->samplesPerPixel;
    Bounds2i sampleBounds(Point2i(0, 0), Point2i(width, height));
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);

    // Precompute the offset on the buffer for each tile
    std::vector<size_t> tileOffsets(nTiles.x * nTiles.y, 0);
    size_t currentTotalOffset = 0;
    for (int y = 0; y < nTiles.y; ++y)
    for (int x = 0; x < nTiles.x; ++x)
    {
        int x0 = sampleBounds.pMin.x + x * tileSize;
        int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
        int y0 = sampleBounds.pMin.y + y * tileSize;
        int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
        Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
        int index = y * nTiles.x + x;
        tileOffsets[index] = currentTotalOffset;
        currentTotalOffset += tileBounds.Area() * spp * m_layout.getSampleSize();
    }

    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            MemoryArena arena;

            // Get sampler instance for tile
            int seed = tile.y * nTiles.x + tile.x;
            std::unique_ptr<Sampler> tileSampler = sppSampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            SamplesPipe pipe;
            pipe.seek(tileOffsets[seed]);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                {
                    ProfilePhase pp(Prof::StartPixel);
                    tileSampler->StartPixel(pixel);
                }

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;


                do {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // NOTE: I'm seeking by pixel only: not support for adaptive sampling
                    SampleBuffer sampleBuffer = pipe.getBuffer();
                    cameraSample.pFilm.x = sampleBuffer.set(IMAGE_X, cameraSample.pFilm.x);
                    cameraSample.pFilm.y = sampleBuffer.set(IMAGE_Y, cameraSample.pFilm.y);
                    cameraSample.pLens.x = sampleBuffer.set(LENS_U, cameraSample.pLens.x);
                    cameraSample.pLens.y = sampleBuffer.set(LENS_V, cameraSample.pLens.y);
                    cameraSample.time = sampleBuffer.set(TIME, cameraSample.time);

                    // Generate camera ray for current sample
                    RayDifferential ray;
                    Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                    ++nCameraRays;

                    // Evaluate radiance along camera ray
                    Spectrum L(0.f);
                    if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena, 0, &sampleBuffer);

                    // Issue warning if unexpected radiance value returned
                    if (L.HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (L.y() < -1e-5) {
                        LOG(ERROR) << StringPrintf(
                            "Negative luminance value, %f, returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            L.y(), pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (std::isinf(L.y())) {
                          LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                        ray << " -> L = " << L;

                    float lum = L.y();
                    if(lum > maxSampleLuminance)
                        L *= maxSampleLuminance / lum;
                    Float rgb[3] = {0.f, 0.f, 0.f};
                    L.ToRGB(rgb);
                    sampleBuffer.set(COLOR_R, rgb[0]*filmScale*rayWeight);
                    sampleBuffer.set(COLOR_G, rgb[1]*filmScale*rayWeight);
                    sampleBuffer.set(COLOR_B, rgb[2]*filmScale*rayWeight);
                    if(std::isfinite(ray.tMax))
                        sampleBuffer.set(DEPTH, ray.tMax);
                    else
                        sampleBuffer.set(DEPTH, 0.f);
                    pipe << sampleBuffer;

                    // Add camera ray's contribution to image
                    // filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
}

void SamplerIntegrator::sparseRender(const Scene& scene,
                                     size_t rest,
                                     int width,
                                     int height,
                                     float filmScale,
                                     float maxSampleLuminance,
                                     size_t offset)
{
    SparsePixelSampler pixelSampler(0, width, 0, height, rest, false);
    RandomSampler sppSampler(1);

    SamplesPipe pipe;
    pipe.seek(offset);

    MemoryArena arena;
    Point2i pixel;
    while(pixelSampler.nextPixel(&pixel))
    {
        {
            ProfilePhase pp(Prof::StartPixel);
            sppSampler.StartPixel(pixel);
        }

        // Do this check after the StartPixel() call; this keeps
        // the usage of RNG values from (most) Samplers that use
        // RNGs consistent, which improves reproducability /
        // debugging.
        if (!InsideExclusive(pixel, pixelBounds))
            continue;

        do {
            // Initialize _CameraSample_ for current sample
            CameraSample cameraSample = sppSampler.GetCameraSample(pixel);

            // NOTE: I'm seeking by pixel only: not support for adaptive sampling
            SampleBuffer sampleBuffer = pipe.getBuffer();
            cameraSample.pFilm.x = sampleBuffer.set(IMAGE_X, cameraSample.pFilm.x);
            cameraSample.pFilm.y = sampleBuffer.set(IMAGE_Y, cameraSample.pFilm.y);
            cameraSample.pLens.x = sampleBuffer.set(LENS_U, cameraSample.pLens.x);
            cameraSample.pLens.y = sampleBuffer.set(LENS_V, cameraSample.pLens.y);
            cameraSample.time = sampleBuffer.set(TIME, cameraSample.time);

            // Generate camera ray for current sample
            RayDifferential ray;
            Float rayWeight = camera->GenerateRayDifferential(cameraSample, &ray);
            ray.ScaleDifferentials(1 / std::sqrt((Float)sppSampler.samplesPerPixel));
            ++nCameraRays;

            // Evaluate radiance along camera ray
            Spectrum L(0.f);
            if (rayWeight > 0) L = Li(ray, scene, sppSampler, arena, 0, &sampleBuffer);

            // Issue warning if unexpected radiance value returned
            if (L.HasNaNs()) {
                LOG(ERROR) << StringPrintf(
                                  "Not-a-number radiance value returned "
                                  "for pixel (%d, %d), sample %d. Setting to black.",
                                  pixel.x, pixel.y,
                                  (int)sppSampler.CurrentSampleNumber());
                L = Spectrum(0.f);
            } else if (L.y() < -1e-5) {
                LOG(ERROR) << StringPrintf(
                                  "Negative luminance value, %f, returned "
                                  "for pixel (%d, %d), sample %d. Setting to black.",
                                  L.y(), pixel.x, pixel.y,
                                  (int)sppSampler.CurrentSampleNumber());
                L = Spectrum(0.f);
            } else if (std::isinf(L.y())) {
                LOG(ERROR) << StringPrintf(
                                  "Infinite luminance value returned "
                                  "for pixel (%d, %d), sample %d. Setting to black.",
                                  pixel.x, pixel.y,
                                  (int)sppSampler.CurrentSampleNumber());
                L = Spectrum(0.f);
            }
            VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                       ray << " -> L = " << L;

            float lum = L.y();
            if(lum > maxSampleLuminance)
                L *= maxSampleLuminance / lum;
            Float rgb[3] = {0.f, 0.f, 0.f};
            L.ToRGB(rgb);
            sampleBuffer.set(COLOR_R, rgb[0]*filmScale*rayWeight);
            sampleBuffer.set(COLOR_G, rgb[1]*filmScale*rayWeight);
            sampleBuffer.set(COLOR_B, rgb[2]*filmScale*rayWeight);
            if(std::isfinite(ray.tMax))
                sampleBuffer.set(DEPTH, ray.tMax);
            else
                sampleBuffer.set(DEPTH, 0.f);
            pipe << sampleBuffer;

            // Add camera ray's contribution to image
            // filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

            // Free _MemoryArena_ memory from computing image sample
            // value
            arena.Reset();
        } while (sppSampler.StartNextSample());
    }
}

}  // namespace pbrt

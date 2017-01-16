/**
 * Copyright (c) 2011, 2012
 * Claudio Kopper <claudio.kopper@icecube.wisc.edu>
 * and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
 * OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 *
 * $Id: I3CLSimModuleHelper.cxx 136718 2015-08-21 20:23:49Z nega $
 *
 * @file I3CLSimModuleHelper.cxx
 * @version $Revision: 136718 $
 * @date $Date: 2015-08-21 22:23:49 +0200 (Fr, 21 Aug 2015) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include "clsim/I3CLSimModuleHelper.h"

#include "clsim/function/I3CLSimFunctionConstant.h"
#include "clsim/function/I3CLSimFunctionFromTable.h"
#include "clsim/function/I3CLSimFunctionDeltaPeak.h"
#include "clsim/random_value/I3CLSimRandomValueInterpolatedDistribution.h"
#include "clsim/random_value/I3CLSimRandomValueWlenCherenkovNoDispersion.h"
#include "clsim/random_value/I3CLSimRandomValueConstant.h"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/variant/get.hpp>

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3MCTreeUtils.h"


namespace I3CLSimModuleHelper {

    namespace {
        double CherenkovYieldDistribution(double wlen, I3CLSimMediumPropertiesConstPtr mediumProperties, double beta=1.)
        {
            I3CLSimFunctionConstPtr nPhaseDist =
            mediumProperties->GetPhaseRefractiveIndex(0); // this assumes the refractive index does not change between layers

            if (!nPhaseDist->HasNativeImplementation())
                log_fatal("The refractive index distribution needs a native implementation to be usable!");

            const double nPhase = nPhaseDist->GetValue(wlen);

            return (2.*M_PI/(137.*(wlen*wlen)))*(1. - 1./ ( std::pow(beta*nPhase,2.) ) ); // dN/dxdwlen
        }

        // the normalization will not be correct here
        double CherenkovYieldDistributionNoDispersion(double wlen)
        {
            return 1./(wlen*wlen); // dN/dxdwlen
        }

    };


    I3CLSimRandomValueConstPtr
    makeWavelengthGenerator(I3CLSimFunctionConstPtr unbiasedSpectrum,
                            I3CLSimFunctionConstPtr wavelengthGenerationBias,
                            I3CLSimMediumPropertiesConstPtr mediumProperties)
    {
        {
            // special handling for delta peaks
            I3CLSimFunctionDeltaPeakConstPtr deltaPeak =
            boost::dynamic_pointer_cast<const I3CLSimFunctionDeltaPeak>(unbiasedSpectrum);
            if (deltaPeak) {
                const double peakPosition = deltaPeak->GetPeakPosition();

                return I3CLSimRandomValueConstantConstPtr
                (new I3CLSimRandomValueConstant(peakPosition));
            }
        }

        // if we get here, it's not a delta peak

        double minWlen = unbiasedSpectrum->GetMinWlen();
        double maxWlen = unbiasedSpectrum->GetMaxWlen();

        // check if the spectrum is from a tabulated distribution (instead of
        // a parameterized one)
        I3CLSimFunctionFromTableConstPtr unbiasedSpectrumFromTable;
        unbiasedSpectrumFromTable =
        boost::dynamic_pointer_cast<const I3CLSimFunctionFromTable>(unbiasedSpectrum);

        if (!unbiasedSpectrumFromTable) {
            // do not clip wavelengths if they are from a tabulated distribution. In that case,
            // re-use the entire table binning
            if (mediumProperties->GetMinWavelength() > minWlen) minWlen=mediumProperties->GetMinWavelength();
            if (mediumProperties->GetMaxWavelength() < maxWlen) maxWlen=mediumProperties->GetMaxWavelength();
        }

        const double wlenRange = maxWlen-minWlen;
        if (wlenRange <= 0.) log_fatal("Internal error, wavelength range <= 0!");

        if (wavelengthGenerationBias->GetMinWlen() > minWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the spectrum wavelength range!");
        if (wavelengthGenerationBias->GetMaxWlen() < maxWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the spectrum wavelength range!");

        // Check if the spectrum values are from a tabulated distribution.
        // If yes, use the table binning, if not make up a binning.
        if (unbiasedSpectrumFromTable)
        {
            std::size_t wlenPoints = unbiasedSpectrumFromTable->GetNumEntries();
            const double firstWlen = unbiasedSpectrumFromTable->GetFirstWavelength();

            std::vector<double> spectrum(wlenPoints, NAN);
            std::vector<double> wavelengths(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = unbiasedSpectrumFromTable->GetEntryWavelength(i);
                const double entry = unbiasedSpectrumFromTable->GetEntryValue(i);
                const double bias = wavelengthGenerationBias->GetValue(wavelength);

                spectrum[i] = bias * entry;
                wavelengths[i] = wavelength;
            }

            if (unbiasedSpectrumFromTable->GetInEqualSpacingMode()) {
                const double wlenStep = unbiasedSpectrumFromTable->GetWavelengthStepping();

                return I3CLSimRandomValueInterpolatedDistributionConstPtr
                (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                                wlenStep,
                                                                spectrum));
            } else {
                // slightly less efficient if non-equally spaced
                return I3CLSimRandomValueInterpolatedDistributionConstPtr
                (new I3CLSimRandomValueInterpolatedDistribution(wavelengths,
                                                                spectrum));
            }
        }
        else
        {
            // use a pre-defined binning of 10ns (an arbitrary value..)

            std::size_t wlenPoints = static_cast<std::size_t>(wlenRange/(10.*I3Units::nanometer))+2;
            const double firstWlen = minWlen;
            const double wlenStep = wlenRange/static_cast<double>(wlenPoints-1);

            std::vector<double> spectrum(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = firstWlen + static_cast<double>(i)*wlenStep;
                const double entry = unbiasedSpectrum->GetValue(wavelength);
                const double bias = wavelengthGenerationBias->GetValue(wavelength);

                spectrum[i] = bias * entry;
            }

            return I3CLSimRandomValueInterpolatedDistributionConstPtr
            (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                            wlenStep,
                                                            spectrum));
        }
    }

    I3CLSimRandomValueConstPtr
    makeCherenkovWavelengthGenerator(I3CLSimFunctionConstPtr wavelengthGenerationBias,
                                     bool generateCherenkovPhotonsWithoutDispersion,
                                     I3CLSimMediumPropertiesConstPtr mediumProperties)
    {
        const double minWlen = mediumProperties->GetMinWavelength();
        const double maxWlen = mediumProperties->GetMaxWavelength();
        const double wlenRange = maxWlen-minWlen;
        if (wlenRange <= 0.) log_fatal("Internal error, wavelength range <= 0!");

        if (wavelengthGenerationBias->GetMinWlen() > minWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the medium property range!");
        if (wavelengthGenerationBias->GetMaxWlen() < maxWlen)
            log_fatal("wavelength generation bias has to have a wavelength range larger or equal to the medium property range!");

        bool noBias=false;
        //bool biasIsConstant=false;

        {
            I3CLSimFunctionConstantConstPtr wavelengthGenerationBiasConstant =
            boost::dynamic_pointer_cast<const I3CLSimFunctionConstant>(wavelengthGenerationBias);

            if (wavelengthGenerationBiasConstant)
            {
                //biasIsConstant=true;

                if ( std::abs(wavelengthGenerationBiasConstant->GetValue((minWlen+maxWlen)/2.)-1.) < 1e-10 )
                    noBias=true;
            }
        }

        I3CLSimFunctionFromTableConstPtr wavelengthGenerationBiasFromTable;
        wavelengthGenerationBiasFromTable =
        boost::dynamic_pointer_cast<const I3CLSimFunctionFromTable>(wavelengthGenerationBias);



        if ((!noBias) && (generateCherenkovPhotonsWithoutDispersion))
        {
            log_warn("**********");
            log_warn(" Using the \"GenerateCherenkovPhotonsWithoutDispersion\" option");
            log_warn(" with a biased photon spectrum generation does not yield a performance");
            log_warn(" increase. You might consider turning this option off to get a better");
            log_warn(" approximation of the Cherenkov spectrum.");
            log_warn("**********");
        }

        // Check if the bias values are from a tabulated distribution.
        // If yes, use the table binning, if not make up a binning.
        if (wavelengthGenerationBiasFromTable)
        {
            std::size_t wlenPoints = wavelengthGenerationBiasFromTable->GetNumEntries();
            const double firstWlen = wavelengthGenerationBiasFromTable->GetFirstWavelength();

            std::vector<double> spectrum(wlenPoints, NAN);
            std::vector<double> wavelengths(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = wavelengthGenerationBiasFromTable->GetEntryWavelength(i);
                const double bias = wavelengthGenerationBiasFromTable->GetEntryValue(i);

                if (generateCherenkovPhotonsWithoutDispersion)
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistributionNoDispersion(wavelength);
                }
                else
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistribution(wavelength, mediumProperties);
                }

                wavelengths[i] = wavelength;
            }

            if (wavelengthGenerationBiasFromTable->GetInEqualSpacingMode()) {
                const double wlenStep = wavelengthGenerationBiasFromTable->GetWavelengthStepping();

                return I3CLSimRandomValueInterpolatedDistributionConstPtr
                (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                                wlenStep,
                                                                spectrum));
            } else {
                // slightly less efficient if non-equally spaced
                return I3CLSimRandomValueInterpolatedDistributionConstPtr
                (new I3CLSimRandomValueInterpolatedDistribution(wavelengths,
                                                                spectrum));
            }
        }
        else if ((noBias) && (generateCherenkovPhotonsWithoutDispersion))
        {
            return I3CLSimRandomValueWlenCherenkovNoDispersionConstPtr
            (new I3CLSimRandomValueWlenCherenkovNoDispersion(minWlen, maxWlen));
        }
        else
        {
            std::size_t wlenPoints = static_cast<std::size_t>(wlenRange/(10.*I3Units::nanometer))+2;
            const double firstWlen = minWlen;
            const double wlenStep = wlenRange/static_cast<double>(wlenPoints-1);

            std::vector<double> spectrum(wlenPoints, NAN);
            for (std::size_t i=0;i<wlenPoints;++i)
            {
                const double wavelength = firstWlen + static_cast<double>(i)*wlenStep;
                const double bias = wavelengthGenerationBias->GetValue(wavelength);

                if (generateCherenkovPhotonsWithoutDispersion)
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistributionNoDispersion(wavelength);
                }
                else
                {
                    spectrum[i] =
                    bias * CherenkovYieldDistribution(wavelength, mediumProperties);
                }
            }

            return I3CLSimRandomValueInterpolatedDistributionConstPtr
            (new I3CLSimRandomValueInterpolatedDistribution(firstWlen,
                                                            wlenStep,
                                                            spectrum));
        }


    }


    // @param OpenCLInitOptions options
    //
    //     struct OpenCLInitOptions {
    //         const I3CLSimOpenCLDevice &device,
    //         I3RandomServicePtr rng,
    //         I3CLSimSimpleGeometryFromI3GeometryPtr geometry,
    //         I3CLSimMediumPropertiesConstPtr medium,
    //         I3CLSimFunctionConstPtr wavelengthGenerationBias,
    //         const std::vector<I3CLSimRandomValueConstPtr> &wavelengthGenerators,
    //         bool enableDoubleBuffering,
    //         bool doublePrecision,
    //         bool stopDetectedPhotons,
    //         bool saveAllPhotons,
    //         double saveAllPhotonsPrescale,
    //         double fixedNumberOfAbsorptionLengths,
    //         double pancakeFactor,
    //         uint32_t photonHistoryEntries,
    //         uint32_t limitWorkgroupSize
    //     }
    //
    //     The boost python bindings apparently do not support so many arguments
    //     passed directly into the method. Therefore, this options struct is used.
    //
    I3CLSimStepToPhotonConverterOpenCLPtr initializeOpenCL(OpenCLInitOptions options)
    {
        I3CLSimStepToPhotonConverterOpenCLPtr conv(new I3CLSimStepToPhotonConverterOpenCL(options.rng, options.device.GetUseNativeMath()));

        conv->SetDevice(options.device);

        conv->SetWlenGenerators(options.wavelengthGenerators);
        conv->SetWlenBias(options.wavelengthGenerationBias);

        conv->SetMediumProperties(options.medium);
        conv->SetGeometry(options.geometry);

        conv->SetEnableDoubleBuffering(options.enableDoubleBuffering);
        conv->SetDoublePrecision(options.doublePrecision);
        conv->SetStopDetectedPhotons(options.stopDetectedPhotons);
        conv->SetSaveAllPhotons(options.saveAllPhotons);
        conv->SetSaveAllPhotonsPrescale(options.saveAllPhotonsPrescale);

        conv->SetFixedNumberOfAbsorptionLengths(options.fixedNumberOfAbsorptionLengths);
        conv->SetDOMPancakeFactor(options.pancakeFactor);

        conv->SetPhotonHistoryEntries(options.photonHistoryEntries);

        conv->Compile();
        //log_trace("%s", conv.GetFullSource().c_str());

        std::size_t maxWorkgroupSize = conv->GetMaxWorkgroupSize();
        if (options.limitWorkgroupSize!=0) {
            maxWorkgroupSize = std::min(static_cast<std::size_t>(options.limitWorkgroupSize), maxWorkgroupSize);
        }

        conv->SetWorkgroupSize(maxWorkgroupSize);
        const std::size_t workgroupSize = conv->GetWorkgroupSize();

        // use approximately the given number of work items, convert to a multiple of the workgroup size
        std::size_t maxNumWorkitems = (static_cast<std::size_t>(options.device.GetApproximateNumberOfWorkItems())/workgroupSize)*workgroupSize;
        if (maxNumWorkitems==0) maxNumWorkitems=workgroupSize;

        conv->SetMaxNumWorkitems(maxNumWorkitems);

        log_info("maximum workgroup size is %zu", maxWorkgroupSize);
        log_info("configured workgroup size is %zu", workgroupSize);
        if (maxNumWorkitems != options.device.GetApproximateNumberOfWorkItems()) {
            log_info("maximum number of work items is %zu (user configured was %" PRIu32 ")", maxNumWorkitems, options.device.GetApproximateNumberOfWorkItems());
        } else {
            log_debug("maximum number of work items is %zu (user configured was %" PRIu32 ")", maxNumWorkitems, options.device.GetApproximateNumberOfWorkItems());
        }

        conv->Initialize();

        return conv;
    }

    I3CLSimLightSourceToStepConverterGeant4Ptr initializeGeant4(I3RandomServicePtr rng,
                                                             I3CLSimMediumPropertiesConstPtr medium,
                                                             I3CLSimFunctionConstPtr wavelengthGenerationBias,
                                                             uint64_t bunchSizeGranularity,
                                                             uint64_t maxBunchSize,
                                                             const I3CLSimLightSourceParameterizationSeries &parameterizationList,
                                                             const std::string &physicsListName,
                                                             double maxBetaChangePerStep,
                                                             uint32_t maxNumPhotonsPerStep,
                                                             bool multiprocessor)
    {
        I3CLSimLightSourceToStepConverterGeant4Ptr conv
        (
         new I3CLSimLightSourceToStepConverterGeant4
         (
          physicsListName,
          maxBetaChangePerStep,
          maxNumPhotonsPerStep
         )
        );

        conv->SetRandomService(rng);
        conv->SetWlenBias(wavelengthGenerationBias);
        conv->SetMediumProperties(medium);
        conv->SetMaxBunchSize(maxBunchSize);
        conv->SetBunchSizeGranularity(bunchSizeGranularity);

        conv->SetLightSourceParameterizationSeries(parameterizationList);

        conv->Initialize();

        return conv;
    }

}

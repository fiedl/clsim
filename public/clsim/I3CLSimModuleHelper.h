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
 * $Id: I3CLSimModuleHelper.h 108199 2013-07-12 21:33:08Z nwhitehorn $
 *
 * @file I3CLSimModuleHelper.h
 * @version $Revision: 108199 $
 * @date $Date: 2013-07-12 23:33:08 +0200 (Fr, 12 Jul 2013) $
 * @author Claudio Kopper
 */

#ifndef I3CLSIMMODULEHELPER_H_INCLUDED
#define I3CLSIMMODULEHELPER_H_INCLUDED


#include "phys-services/I3RandomService.h"

#include "clsim/random_value/I3CLSimRandomValue.h"
#include "clsim/function/I3CLSimFunction.h"

#include "clsim/I3CLSimMediumProperties.h"
#include "clsim/I3CLSimSimpleGeometryFromI3Geometry.h"

#include "clsim/I3CLSimStepToPhotonConverterOpenCL.h"
#include "clsim/I3CLSimLightSourceToStepConverterGeant4.h"

#include "clsim/I3CLSimLightSourceParameterization.h"

#include "clsim/I3CLSimOpenCLDevice.h"

#include <vector>
#include <string>

namespace I3CLSimModuleHelper {

    struct OpenCLInitOptions {
        const I3CLSimOpenCLDevice &device;
        I3RandomServicePtr rng;
        I3CLSimSimpleGeometryFromI3GeometryPtr geometry;
        I3CLSimMediumPropertiesConstPtr medium;
        I3CLSimFunctionConstPtr wavelengthGenerationBias;
        const std::vector<I3CLSimRandomValueConstPtr> &wavelengthGenerators;
        bool enableDoubleBuffering;
        bool doublePrecision;
        bool stopDetectedPhotons;
        bool saveAllPhotons;
        double saveAllPhotonsPrescale;
        double maxNumOutputPhotonsCorrectionFactor;
        bool simulateHoleIce;
        double holeIceScatteringLengthFactor;
        double holeIceAbsorptionLengthFactor;
        double fixedNumberOfAbsorptionLengths;
        double pancakeFactor;
        uint32_t photonHistoryEntries;
        uint32_t limitWorkgroupSize;
        I3Vector<I3Position> holeIceCylinderPositions;
        I3Vector<float> holeIceCylinderRadii;
        I3Vector<float> holeIceCylinderScatteringLengths;
        I3Vector<float> holeIceCylinderAbsorptionLengths;
    };

    I3CLSimStepToPhotonConverterOpenCLPtr
    initializeOpenCL(OpenCLInitOptions options);

    I3CLSimLightSourceToStepConverterGeant4Ptr
    initializeGeant4(I3RandomServicePtr rng,
                     I3CLSimMediumPropertiesConstPtr medium,
                     I3CLSimFunctionConstPtr wavelengthGenerationBias,
                     uint64_t bunchSizeGranularity,
                     uint64_t maxBunchSize,
                     const I3CLSimLightSourceParameterizationSeries &parameterizationList,
                     const std::string &physicsListName,
                     double maxBetaChangePerStep,
                     uint32_t maxNumPhotonsPerStep,
                     bool multiprocessor=false);

    I3CLSimRandomValueConstPtr
    makeCherenkovWavelengthGenerator(I3CLSimFunctionConstPtr wavelengthGenerationBias,
                                     bool generateCherenkovPhotonsWithoutDispersion,
                                     I3CLSimMediumPropertiesConstPtr mediumProperties);

    I3CLSimRandomValueConstPtr
    makeWavelengthGenerator(I3CLSimFunctionConstPtr unbiasedSpectrum,
                            I3CLSimFunctionConstPtr wavelengthGenerationBias,
                            I3CLSimMediumPropertiesConstPtr mediumProperties);

};


#endif //I3CLSIMMODULEHELPER_H_INCLUDED

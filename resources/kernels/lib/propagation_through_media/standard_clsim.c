// This is the standard clsim algorithm for propagation through different media.
//
// Source:
// https://github.com/fiedl/clsim/blob/icesim-v05-00-07/resources/kernels/propagation_kernel.c.cl
//

#ifndef STANDARD_CLSIM_C
#define STANDARD_CLSIM_C

inline void apply_propagation_through_different_media_with_standard_clsim(
    floating4_t photonPosAndTime, floating4_t photonDirAndWlen,
    floating_t *sca_step_left, floating_t *abs_lens_left,
    floating_t *distancePropagated, floating_t *distanceToAbsorption) {

  // this block is along the lines of the PPC kernel
  {
#ifdef getTiltZShift_IS_CONSTANT
#define effective_z (photonPosAndTime.z - getTiltZShift_IS_CONSTANT)
#else
    // apply ice tilt
    const floating_t effective_z =
        photonPosAndTime.z - getTiltZShift(photonPosAndTime);
    int currentPhotonLayer =
        min(max(findLayerForGivenZPos(effective_z), 0), MEDIUM_LAYERS - 1);
#endif

    const floating_t photon_dz = photonDirAndWlen.z;

    // add a correction factor to the number of absorption lengths *abs_lens_left
    // before the photon is absorbed. This factor will be taken out after this
    // propagation step. Usually the factor is 1 and thus has no effect, but it
    // is used in a direction-dependent way for our model of ice anisotropy.
    const floating_t abs_len_correction_factor =
        getDirectionalAbsLenCorrFactor(photonDirAndWlen);

    *abs_lens_left *= abs_len_correction_factor;

    // the "next" medium boundary (either top or bottom, depending on step
    // direction)
    floating_t mediumBoundary = (photon_dz < ZERO)
                                    ? (mediumLayerBoundary(currentPhotonLayer))
                                    : (mediumLayerBoundary(currentPhotonLayer) +
                                       (floating_t)MEDIUM_LAYER_THICKNESS);

    // track this thing to the next scattering point
#ifdef PRINTF_ENABLED
// dbg_printf("   - next scatter in %f scattering lengths\n", *sca_step_left);
#endif

    floating_t currentScaLen =
        getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
    floating_t currentAbsLen =
        getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);

    floating_t ais =
        (photon_dz * *sca_step_left -
         my_divide((mediumBoundary - effective_z), currentScaLen)) *
        (ONE / (floating_t)MEDIUM_LAYER_THICKNESS);
    floating_t aia =
        (photon_dz * *abs_lens_left -
         my_divide((mediumBoundary - effective_z), currentAbsLen)) *
        (ONE / (floating_t)MEDIUM_LAYER_THICKNESS);

#ifdef PRINTF_ENABLED
// dbg_printf("   - ais=%f, aia=%f, j_initial=%i\n", ais, aia,
// currentPhotonLayer);
#endif

    // propagate through layers
    int j = currentPhotonLayer;
    if (photon_dz < 0) {
      for (; (j > 0) && (ais < ZERO) && (aia < ZERO);
           mediumBoundary -= (floating_t)MEDIUM_LAYER_THICKNESS,
           currentScaLen = getScatteringLength(j, photonDirAndWlen.w),
           currentAbsLen = getAbsorptionLength(j, photonDirAndWlen.w),
           ais += my_recip(currentScaLen), aia += my_recip(currentAbsLen))
        --j;
    } else {
      for (; (j < MEDIUM_LAYERS - 1) && (ais > ZERO) && (aia > ZERO);
           mediumBoundary += (floating_t)MEDIUM_LAYER_THICKNESS,
           currentScaLen = getScatteringLength(j, photonDirAndWlen.w),
           currentAbsLen = getAbsorptionLength(j, photonDirAndWlen.w),
           ais -= my_recip(currentScaLen), aia -= my_recip(currentAbsLen))
        ++j;
    }

#ifdef PRINTF_ENABLED
// dbg_printf("   - j_final=%i\n", j);
#endif

    if ((currentPhotonLayer == j) || ((my_fabs(photon_dz)) < EPSILON)) {
      *distancePropagated = *sca_step_left * currentScaLen;
      *distanceToAbsorption = *abs_lens_left * currentAbsLen;
    } else {
      const floating_t recip_photon_dz = my_recip(photon_dz);
      *distancePropagated =
          (ais * ((floating_t)MEDIUM_LAYER_THICKNESS) * currentScaLen +
           mediumBoundary - effective_z) *
          recip_photon_dz;
      *distanceToAbsorption =
          (aia * ((floating_t)MEDIUM_LAYER_THICKNESS) * currentAbsLen +
           mediumBoundary - effective_z) *
          recip_photon_dz;
    }
#ifdef getTiltZShift_IS_CONSTANT
    currentPhotonLayer = j;
#endif

#ifdef PRINTF_ENABLED
// dbg_printf("   - *distancePropagated=%f\n", *distancePropagated);
#endif

    // get overburden for distance
    if (*distanceToAbsorption < *distancePropagated) {
      *distancePropagated = *distanceToAbsorption;
      *abs_lens_left = ZERO;
    } else {
      *abs_lens_left =
          my_divide(*distanceToAbsorption - *distancePropagated, currentAbsLen);
    }

    // hoist the correction factor back out of the absorption length
    *abs_lens_left = my_divide(*abs_lens_left, abs_len_correction_factor);
  }
}

#endif
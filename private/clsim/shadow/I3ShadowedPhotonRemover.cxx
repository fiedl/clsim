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
 * $Id: I3ShadowedPhotonRemover.cxx 108199 2013-07-12 21:33:08Z nwhitehorn $
 *
 * @file I3ShadowedPhotonRemover.cxx
 * @version $Revision: 108199 $
 * @date $Date: 2013-07-12 23:33:08 +0200 (Fr, 12 Jul 2013) $
 * @author Claudio Kopper
 */

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif 
#include <inttypes.h>

#include <limits>

#include "clsim/shadow/I3ShadowedPhotonRemover.h"

#include "dataclasses/I3Constants.h"

I3ShadowedPhotonRemover::I3ShadowedPhotonRemover() //(I3ExtraGeometry extraGeometry) 
//:
//extraGeometry_(extraGeometry)
{
    log_trace("%s", __PRETTY_FUNCTION__);

}

I3ShadowedPhotonRemover::~I3ShadowedPhotonRemover()
{
    log_trace("%s", __PRETTY_FUNCTION__);
}



bool I3ShadowedPhotonRemover::IsPhotonShadowed(const I3Photon &photon) const
{

    
    return false;
}

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
 * $Id: I3ShadowedPhotonRemover.cxx 136718 2015-08-21 20:23:49Z nega $
 *
 * @file I3ShadowedPhotonRemover.cxx
 * @version $Revision: 136718 $
 * @date $Date: 2015-08-21 22:23:49 +0200 (Fr, 21 Aug 2015) $
 * @author Claudio Kopper
 */

#include <sstream>

#include <clsim/shadow/I3ShadowedPhotonRemover.h>

#include <boost/preprocessor/seq.hpp>

using namespace boost::python;
namespace bp = boost::python;


void register_I3ShadowedPhotonRemover()
{
    {
        bp::scope I3ShadowedPhotonRemover_scope = 
        bp::class_<I3ShadowedPhotonRemover, boost::shared_ptr<I3ShadowedPhotonRemover>, boost::noncopyable>
        ("I3ShadowedPhotonRemover", 
         bp::init<>()
         //(
         // (
         //  bp::arg("random"),
         //  bp::arg("DOMRadius"),
         //  bp::arg("holeRadius"),
         //  bp::arg("mediumProperties"),
         //  bp::arg("holeIceAbsorptionLength"),
         //  bp::arg("holeIceScatteringLength")
         // )
         //)
        )

        .def("IsPhotonShadowed", &I3ShadowedPhotonRemover::IsPhotonShadowed)
        ;
    }
    
    bp::implicitly_convertible<boost::shared_ptr<I3ShadowedPhotonRemover>, boost::shared_ptr<const I3ShadowedPhotonRemover> >();
}

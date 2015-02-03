
#include "icetray/I3Units.h"
#include "clsim/tabulator/Axes.h"
#include "clsim/I3CLSimHelperToFloatString.h"
#include "opencl/I3CLSimHelperLoadProgramSource.h"

#include <boost/foreach.hpp>

namespace {

std::string 
loadKernel(const std::string& name, bool header=false)
{
    const std::string I3_SRC(getenv("I3_SRC"));
    const std::string kernelBaseDir = I3_SRC+"/clsim/resources/kernels/";
    const std::string ext = header ? ".h.cl" : ".c.cl";
    return I3CLSimHelper::LoadProgramSource(kernelBaseDir+name+ext);
}

}

namespace clsim {

namespace tabulator {

Axes::Axes(const std::vector<value_type> &axes) : axes_(axes), n_dim_(axes_.size()),
    shape_(n_dim_, 1), strides_(n_dim_, 1)
{
	int i = n_dim_-1;
	shape_[i] = axes_[i]->GetNBins();
	for (i--; i >= 0; i--)
		strides_[i] *= strides_[i+1]*shape_[i+1];
	n_bins_ = strides_[0]*shape_[0];
}

std::string
Axes::GetBinIndexFunction() const
{
	std::ostringstream ss;
	ss << "inline uint getBinIndex(coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return ";
	
	for (uint i = 0; i < axes_.size(); i++) {
		std::ostringstream var;
		var << "coords.s" << i;
		ss << strides_[i] << "*" << axes_[i]->GetIndexCode(var.str());
		if (i+1 < axes_.size())
			ss << "\n + ";
	}
	
	ss << ";\n}\n";
	
	return ss.str();
}

std::string
Axes::GenerateBinningCode() const
{
	return GetCoordinateFunction() + "\n"
	    + GetBoundsCheckFunction() + "\n"
	    + GetBinIndexFunction() + "\n";
}

double
Axes::GetBinVolume(size_t idx) const
{
	// unravel index
	size_t idxs[n_dim_];
	for (unsigned j=0; j < n_dim_; j++) 
		idxs[j] = idx/strides_[j] % shape_[j];

	return GetBinVolume(idxs);
}

std::string
SphericalAxes::GetCoordinateFunction() const
{
	return loadKernel("spherical_coordinates");
}

std::string
SphericalAxes::GetBoundsCheckFunction() const
{
	std::ostringstream ss;
	ss << "inline bool isOutOfBounds(const coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return (coords.s3 > "<<I3CLSimHelper::ToFloatString(this->at(3)->GetMax())<<");";
	ss << "\n}\n";

	return ss.str();
}

double
SphericalAxes::GetBinVolume(size_t *const idxs) const
{
	// NB: since we combine the bins at azimuth > 180 degrees with the
	// other half of the sphere, the true volume of an azimuthal bin is
	// twice its nominal value.
	return ((std::pow(at(0)->GetBinEdge(idxs[0]+1), 3) - std::pow(at(0)->GetBinEdge(idxs[0]), 3))/3.)
	    * 2*I3Units::degree*(at(1)->GetBinEdge(idxs[1]+1) - at(1)->GetBinEdge(idxs[1]))
	    * (at(2)->GetBinEdge(idxs[2]+1) - at(2)->GetBinEdge(idxs[2]));
}

std::string
CylindricalAxes::GetCoordinateFunction() const
{
	return loadKernel("cylindrical_coordinates");
}

std::string
CylindricalAxes::GetBoundsCheckFunction() const
{
	std::ostringstream ss;
	ss << "inline bool isOutOfBounds(const coordinate_t coords)";
	ss << "\n{\n";
	ss << "    return (coords.s3 > "<<I3CLSimHelper::ToFloatString(this->at(3)->GetMax())<<");";
	ss << "\n}\n";

	return ss.str();
}

double
CylindricalAxes::GetBinVolume(size_t *const idxs) const
{
	// NB: since we combine the bins at azimuth > pi with the
	// other half of the cylinder, the true volume of an azimuthal bin is
	// twice its nominal value.
	return ((std::pow(at(0)->GetBinEdge(idxs[0]+1), 2) - std::pow(at(0)->GetBinEdge(idxs[0]), 2))/2.)
	    * 2*(at(1)->GetBinEdge(idxs[1]+1) - at(1)->GetBinEdge(idxs[1]))
	    * (at(2)->GetBinEdge(idxs[2]+1) - at(2)->GetBinEdge(idxs[2]));
}

}

}
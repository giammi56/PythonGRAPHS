#include "CH_Electron.h"
#include "CH_Spectrometer.h"
#include "CH_Basics.h"

#include "Math.h"

namespace CH
{

	electron_class::electron_class()
	{
		this->raw.m = 1.0 / MASSAU; // all masses are in AMU!
		this->raw.q = -1.0;
	}
	electron_class::electron_class(particle_raw rawdata)
	{
		this->raw = rawdata;

		this->raw.m = 1.0 / MASSAU; // all masses are in AMU!
		this->raw.q = -1.0;
	}

	electron_class::~electron_class()
	{
	}

	void electron_class::set_raw(particle_raw rawdata) {
		this->raw.data = rawdata.data;
		this->raw.m = 1.0 / MASSAU;
		this->raw.q = -1.0;
		this->raw.method = rawdata.method;
		this->raw.phi = rawdata.phi;
	}


	void electron_class::cut_wiggle(spectrometer_class * spect, double width)
	{
		if(fabs(fmod(this->raw.data.tof,spect->Bfield_ns))<width)
			valid = false;
	}
}
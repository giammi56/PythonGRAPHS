#include "CH_Ion.h"
#include "CH_Basics.h"

#include "Math.h"

namespace CH
{

	ion_class::ion_class()
	{
		this->raw.m = 1.0; // all masses are in AMU, default ion is  proton.
		this->raw.q = 1.0;
	}

	ion_class::~ion_class()
	{
	}


	void ion_class::clone_from(ion_class * ion)
	{
		if (!ion) return;
		this->channel = ion->channel;
		this->cor = ion->cor;
		this->mom = ion->mom;
		this->raw = ion->raw;
		this->t_mean = ion->t_mean;
		this->t_width = ion->t_width;
		this->valid = ion->valid;
	}


}
#pragma warning(disable : 4996)
#include "CH_Ranges.h"

namespace CH 
{

	range_class::range_class() {
		this->range_map_from = new double[10,30];
		this->range_map_to = new double[10,30];
		prepare_tags();
		init_default_values();
	}

	range_class::~range_class() {
		delete this->range_map_from;
		delete this->range_map_to;
	}

	void range_class::prepare_tags() {
		// particle types
		this->particle_types.push_back("ELECTRON");	
		this->particle_types.push_back("ION");
		this->particle_types.push_back("PROJECTILE");
		this->particle_types.push_back("DIATOMIC");

		//single particle variables
		this->single_vars.push_back("TOF");
		this->single_vars.push_back("P");
		this->single_vars.push_back("PX");
		this->single_vars.push_back("PY");
		this->single_vars.push_back("PZ");
		this->single_vars.push_back("ENERGY");
		
		// diatomic molecule variables
		this->diatomic_vars.push_back("PREL");
		this->diatomic_vars.push_back("PRELX");
		this->diatomic_vars.push_back("PRELY");
		this->diatomic_vars.push_back("PRELZ");
		this->diatomic_vars.push_back("PCM");
		this->diatomic_vars.push_back("PCMX");
		this->diatomic_vars.push_back("PCMY");
		this->diatomic_vars.push_back("PCMZ");
		this->diatomic_vars.push_back("KER");
	
	}

	void range_class::init_default_values() {
		// electrons
		this->range_map_from[RA_ELECTRON, TOF] = 0.0;
		this->range_map_from[RA_ELECTRON, P] = 0.0;
		this->range_map_from[RA_ELECTRON, PX] = -2.0;
		this->range_map_from[RA_ELECTRON, PY] = -2.0;
		this->range_map_from[RA_ELECTRON, PZ] = -2.0;
		this->range_map_from[RA_ELECTRON, ENERGY] = 0.0;

		this->range_map_to[RA_ELECTRON, TOF] = 80.0;
		this->range_map_to[RA_ELECTRON, P] = 2.0;
		this->range_map_to[RA_ELECTRON, PX] = 2.0;
		this->range_map_to[RA_ELECTRON, PY] = 2.0;
		this->range_map_to[RA_ELECTRON, PZ] = 2.0;
		this->range_map_to[RA_ELECTRON, ENERGY] = 20.0;

		// ions
		this->range_map_from[RA_ION, TOF] = 0.0;
		this->range_map_from[RA_ION, P] = 0.0;
		this->range_map_from[RA_ION, PX] = -1.0;
		this->range_map_from[RA_ION, PY] = -1.0;
		this->range_map_from[RA_ION, PZ] = -1.0;
		this->range_map_from[RA_ION, ENERGY] = 0.0;

		this->range_map_to[RA_ION, TOF] = 8000.0;
		this->range_map_to[RA_ION, P] = 1.0;
		this->range_map_to[RA_ION, PX] = 1.0;
		this->range_map_to[RA_ION, PY] = 1.0;
		this->range_map_to[RA_ION, PZ] = 1.0;
		this->range_map_to[RA_ION, ENERGY] = 0.1;

		// projectile
		this->range_map_from[RA_PROJECTILE, TOF] = 0.0;
		this->range_map_from[RA_PROJECTILE, P] = 0.0;
		this->range_map_from[RA_PROJECTILE, PX] = -1.0;
		this->range_map_from[RA_PROJECTILE, PY] = -1.0;
		this->range_map_from[RA_PROJECTILE, PZ] = -1.0;
		this->range_map_from[RA_PROJECTILE, ENERGY] = 0.0;

		this->range_map_to[RA_PROJECTILE, TOF] = 8000.0;
		this->range_map_to[RA_PROJECTILE, P] = 1.0;
		this->range_map_to[RA_PROJECTILE, PX] = 1.0;
		this->range_map_to[RA_PROJECTILE, PY] = 1.0;
		this->range_map_to[RA_PROJECTILE, PZ] = 1.0;
		this->range_map_to[RA_PROJECTILE, ENERGY] = 0.1;

		// diatomic molecule
		this->range_map_from[RA_DIATOMIC, PREL] = 0.0;
		this->range_map_from[RA_DIATOMIC, PRELX] = -100.0;
		this->range_map_from[RA_DIATOMIC, PRELY] = -100.0;
		this->range_map_from[RA_DIATOMIC, PRELZ] = -100.0;
		this->range_map_from[RA_DIATOMIC, PCM] = 0.0;
		this->range_map_from[RA_DIATOMIC, PCMX] = -4.0;
		this->range_map_from[RA_DIATOMIC, PCMY] = -4.0;
		this->range_map_from[RA_DIATOMIC, PCMZ] = -4.0;
		this->range_map_from[RA_DIATOMIC, KINER] = 0.0;

		this->range_map_to[RA_DIATOMIC, PREL] = 100.0;
		this->range_map_to[RA_DIATOMIC, PRELX] = 100.0;
		this->range_map_to[RA_DIATOMIC, PRELY] = 100.0;
		this->range_map_to[RA_DIATOMIC, PRELZ] = 100.0;
		this->range_map_to[RA_DIATOMIC, PCM] = 4.0;
		this->range_map_to[RA_DIATOMIC, PCMX] = 4.0;
		this->range_map_to[RA_DIATOMIC, PCMY] = 4.0;
		this->range_map_to[RA_DIATOMIC, PCMZ] = 4.0;
		this->range_map_to[RA_DIATOMIC, KINER] = 20.0;
	
		array_size[RA_ELECTRON] = NUM_SINGLE_VARS;
		array_size[RA_ION] = NUM_SINGLE_VARS;
		array_size[RA_PROJECTILE] = NUM_SINGLE_VARS;
		array_size[RA_DIATOMIC] = NUM_DIATOMIC_VARS;

	}

	double range_class::get_from(int particle_type, int var) {
		// this is not supposed to happen, but prevents possible crashes...
		if(var >= array_size[particle_type])
			return 0.0;
		return this->range_map_from[particle_type, var];
	}

	double range_class::get_to(int particle_type, int var) {
		// this is not supposed to happen, but prevents possible crashes...
		if(var >= array_size[particle_type])
			return 0.0;
		return this->range_map_to[particle_type, var];
	}

	bool range_class::set_from(int particle_type, int var, double val) {
		if(var >= array_size[particle_type])
			return false;
		this->range_map_from[particle_type, var] = val;
		return true;
	}

	bool range_class::set_to(int particle_type, int var, double val) {
		if(var >= array_size[particle_type])
			return false;
		this->range_map_to[particle_type, var] = val;
		return true;
	}

}

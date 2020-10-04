#ifndef CH_POLYATOMIC_ALREADY_INCLUDED
	#define CH_POLYATOMIC_ALREADY_INCLUDED
	
#include "CH_Basics.h"
#include "CH_Ion.h"
#include "..\CH_Event.h"
#include "CH_coordinate_system.h"
	
#define USE_ION_MATRIX_RECURSIVE

namespace CH
{
class polyatomic_class{

private:
			mom_cor_param * mfac;

public:
//private:
		int nr_of_combinations;
		int max_nr_of_combinations;
		int* used_hits; //array of hits that are already matched with an ion
		double* value_to_sort; //array of the quantity that is used to decide which is the correct match
		int* index; //index to keep track during sorting the array "value_to_sort"
		double* temp_mass; //temporary array for swapping the masses
		double* temp_charge; //temporary array for swapping the charges
		ion_class ** temp_ion;

		void set_mom_fac(mom_cor_param * fac);

		bool check_validity(double var, double min, double max, bool not);
		void reset();
		void fill_ion_matrix(CH_event_struct* evti);
		bool sort_ion_matrix(int iteration, CH_event_struct *evti, spectrometer_class *spect, reaction_struct* reaction);
		bool sort_momenta();
		//void labframe_transformation(int direction, CH_vector a, CH_vector b);
		void trafo_to_molframe(Coordinate_System frame);

//public:
		double max_target_value;
		double ambiguity_parameter;
		bool valid;
		int channel;
		int number_of_ions;
		bool incomplete;
		//ion_class ** input_ion;
		ion_class ** ion;
		int** ion_matrix;
		int** matches;
		CH_vector * momenta_in_mol_frame;

		double KER();
		double total_mass;
		CH_vector mom_cm; //center-of-mass / sum momentum
		double momentum_magnitude_sum;

		polyatomic_class();
		~polyatomic_class();
		


		void set_ions(int number_of_ions, ion_class** ions);
		void set_channel(int channel);
		void set_channel(int channel, bool valid);
		void set_valid(bool valid);
	
		bool run_ion_matrix(reaction_struct *reaction, CH_event_struct* evti, spectrometer_class *spect);
		
 
		void process(spectrometer_class* spect);
		void calculate_sum_momenta();

		void invalidate(rtag_struct* rtag);

}; //end class Polyatomic

/*
public: 
  int ID;  //index for the channel
  int histindex;  //for the ID of the histograms
  char * channel_name; //for naming histograms and folders
  int nr_of_fragments; //number of fragments in the channel
  int * channel_species; //which ions from the array coulomb_explosion.ion_species are used in this channel?
  double * channel_mass; //masses of these ions
  double * channel_charge; //charges of these ions
  int max_nr_of_combinations; //how many combinations hit<->fragment should be tested in each channe
  double * psum; //sum momentum for each combination
  int * index; //array to keep track when sorting psum
  double pmax[3]; //maximum momentum (x,y,z) to write a combination into the psum array
  double psum_min[3]; //best sum momentum in this channel
  double pz_offset; //to correct psum_z if calibration was not precise
  double max_method; //which reconstruction method should still be accepted? (parameter 1349 )


private:
//variables for recursive algorithm
	int* used_hits; //during testing of a combination: which hits have already been used?
	int iteration; //how many instances of the recursive function have been called?+
	bool method_check; 
	

public: 
  ~Polyatomic(); 
  bool init(int id, double parameter[], double* mass, double* charge); //initializes the class with parameters from the config-file (called only once during the run; in LMF2root.cpp)
  bool reset();
  bool sort_ion_matrix(int nr_of_hits, int iteration, int& nr_of_combinations_in_channel, int** ion_matrix, int** matches, double *rx, double* ry, double* rtof, double* rmethod, double parameter[]); //recursive algorithm
  bool draw_IM_histos(Ueberstruct * Ueber); //draws some histograms for evaluating the recursive algorithm

}; //end class Polyatomic
*/
/*
class coulomb_explosion{
public:
	//********DATA MEMBERS*****************
  int nr_of_ion_species;
  int nr_of_fragmentation_channels;
  int max_nr_of_combinations; //how many combinations hit<->fragment should be tested in each channel?
  int max_nr_of_fragments; //number of fragments in the largest fragmentation channel
   int max_nr_of_hits; //maximum number of hits, later set to 16
  
	//arrays for characterization of all possible ions
	double * mass;
	double * charge;
	double * p0; //momentum that is used to calculate the time window
	double * tmin; //time window
	double * tmax;
  
	int ** ion_matrix;
	
	fragmentation_channel * break_up;//creates an instance of "fragmentation_channel" for every break-up channel that is specified in config.txt

	int ** matches;//2d-array; for each match, the order of the hits is stored in a 1d-array
	
	double * psum_all; //array for storing all the calculated sum momenta (concatenates all break_up[i].psum[] )	
	int * index_psum_all; //array for keeping track during sorting of the sum momenta (concatenates all break_up[i].index[] )
	
	double psum_max; //maximum sum momentum that the best channel is allowed to have
	double ambiguity_parameter_match; //in order to check if only the best channel has a low sum momentum; if not, then the event is discarded due to ambiguity
	double ambiguity_parameter_channel; //in order to check if only the best channel has a low sum momentum; if not, then the event is discarded due to ambiguity

	int final_channel; //channel that has the lowest sum momentum
	int nr_of_fragments_final;//number of fragments in that channel


	double* MASS;
	double* CHARGE;
	double* RX ;
	double* RY ;
	double* RTOF;

	double* px;
	double* py;
	double* pz;
	double* p2sum;
	double* psum;
	

	char hist_title[50];
	//********METHODS***********

  //coulomb_explosion(int nr_of_fragmentation_channels); 
  ~coulomb_explosion();
  bool init(double parameter[]); //initializes the class with parameters from the config-file (called only once during the run; in LMF2root.cpp)
  bool reset(); //resets certain values for a new event (called for every event, in analysis.cpp)
  bool calculate_time_windows(double parameter[]); //determines time windows for each ion species (called only once during the run; in LMF2root.cpp)
  bool fill_ion_matrix(int nr_of_hits, double* rtof); //decides which hit could be linked to which ion 
  bool find_best_match(Ueberstruct * Ueber); //sorts all possibles matches according to sum momentum and decides if the one with lowest sum momentum is accepted
  bool return_best_match(int& flag, double* rx, double* ry, double* rtof);//writes masses, charges and raw data into arrays for further analysis

  
    
};

*/

} //end namespace CH

#endif;
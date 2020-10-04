#ifndef IPA_ALREADY_INCLUDED
	#define IPA_ALREADY_INCLUDED


	#include "OS_Version.h"

	#include "../ColAHelL/src/ColAHelL.h"
	#include "../src/rootstuff.h"
	#include "../src/Histo.h"

//	#define IPA_MAX_HITS_REC 3
//	#define IPA_MAX_HITS_ELEC 2
//	#define IPA_MAX_HITS_PROJ 1

	#define IPAVERSION (0.26)
	using namespace CH;

	#define MAXTHREADS 8

	#define EDETPAGE 0
	#define RDETPAGE 1
	#define PDETPAGE 2
	#define ESPECPAGE 3
	#define RSPECPAGE 4
	#define EMOMPAGE 5
	#define RMOMPAGE 6
	#define PPCOPAGE 7
	#define USERPAGE 9

	struct ipa_parameter 
	{
		double value;
		double min;
		double max;
		double step;

		char *name;
	};

	struct bpoint
	{
		float x;
		float y;
	};

	class ipa_class
	{	
		private:
			int cores_in_use;

			int channel;
			int num_events;
			int page;

			int num_events_read;
			CH_event_struct *eventdata;
			CH_event_struct *ipa_evt;

			ipa_parameter *par;
			int num_par;

			char *keys;
			double precision;

			int NumPars[16];

			rootstuff *root;
			TCanvas * myCanvas[16];
			TH1D * H1D[MAXTHREADS][16];
			TH2D * H2D[MAXTHREADS][16];

//			bool interpolation;
			bool PPCO_enabled;
			bool PPCO_momcon;
			TCanvas * myCanvas2;
			TH2D * H_PIPICO;
			int pco_bins;
			double pco_max;
			double pco_min;
			double ppco_m1[4];
			double ppco_m2[4];
			double ppco_q1[4];
			double ppco_q2[4];

			// for autotuning
			TH2D *ctpeArray[MAXTHREADS][100];
			TObjArray *aSlices[100];
			TH1D *Fit[100];
			bool LUT_enabled;
			int LUT_electrons;
			double LUT_mom;
			int LUT_bins_mom;
			int LUT_bins_phi;
			int LUT_bins_ctheta;
			int LUT_zero;
			int LUT_mergeNbins;
			int LUT_slidingmerge;
			double LUT_momrange_min;
			double LUT_momrange_max;
			int LUT_Profile;
	
			// wiggle tool
			bool wiggletool_enabled;
			TH2D * H_wiggletool_x;
			TH2D *H_wiggletool_y;
			TH2D *H_wiggletool_r;
			double wiggletool_max;
			double wiggletool_min;
			double wiggletool_t_0_elec;
			double wiggletool_t_0_tof_calc;

			// fish tool
			bool fishtool_enabled;
			TH2D * H_fishtool_x;
			TH2D *H_fishtool_y;
			TH2D *H_fishtool_r;
			double fishtool_max;
			double fishtool_min;
			double fishtool_t_0_elec;
			double fishtool_ee1;
			double fishtool_ee2;
			double fishtool_ee3;
			TMarker *mrk[3][3][36];
			bpoint *knots[128];
			bpoint *firstControlPoints[128];
			bpoint *secondControlPoints[128];

			// ion tof tool
			bool iontoftool_enabled;
			TH2D * H_iontoftool_x;
			TH2D *H_iontoftool_y;
			double iontoftool_max;
			double iontoftool_min;
			double dx;
			double dy;
			double dt;
			double drot;

			processor_class* IPA_Proc[MAXTHREADS];
			tof_calc_class* IPA_Tof;

			//spectrometer_class *Spect;
			double len_e[10];
			double E_e[10];
			//GN
			double MeanTOFe[10];
			//GN_!
			double len_r[10];
			double E_r[10];

			//double cor_tab[36][30];
			
			//cor_param * e_fac;
			//cor_param * r_fac;
			//cor_param * p_fac;
			//mom_cor_param * e_mom_fac;
			//mom_cor_param * r_mom_fac;
			//mom_cor_param * p_mom_fac;
			//double e_rot;
			//double r_rot;
			//double p_rot;
			std::vector<double>    interpolation_adjustment_ceff;
			reaction_struct *cur_reaction;
			bool reaction_def_valid;
			
			bool cut_wiggle;
			double cut_wiggle_width;

			void std_calc(CH_event_struct *evti);
			void user_function(CH_event_struct *evti, int th_num);
			void WriteParToStruct();
			void combine_histos();
			void user_draw();
		
			// special IPAs
			void init_PIPICO();
			void draw_PIPICO();
			void draw_PIPICOline(int color, double MassP1_amu, double MassP2_amu, double ChargeP1_au, double ChargeP2_au);

			void init_extract_LUT();
			void fill_extract_LUT(CH_event_struct *evti, int th_num);
			void draw_extract_LUT();
			void write_correctiontable();

			void init_wiggletool();
			void draw_wiggletool();
			void draw_wiggletool_markers(bool r=false);

			void init_fishtool();
			void draw_fishtool();
			void draw_fishtool_markers(int type, double *U);
			double TTOT(double pz, double *U); 
			float getPt(float n1, float n2 , float perc);
			bpoint BezierQuad(bpoint *P0, bpoint *P1, bpoint *P2, float perc);
			bpoint BezierCube(bpoint *P0, bpoint *P1, bpoint *P2, bpoint *P3, float perc);
			void CalcControlPoints(int num_knots);
			void ipa_class::GetFirstControlPoints(double *rhs, int n, double *res);

			void init_iontoftool();
			void draw_iontoftool();
			void draw_iontoftool_markers();
			void draw_iontoftool_jet_markers(bool y);

			electron_class * e[MAXTHREADS][16];
			ion_class * r[MAXTHREADS][16];
			diatomic_class * mol[MAXTHREADS];
			polyatomic_class * big_mol[MAXTHREADS];

			particle_raw * e_raw[16];
			particle_raw * r_raw[16];
			particle_raw * p_raw[16];

			double *cor_tab[36][30];

			void init_user_histograms();

		public:

			ipa_class();
			~ipa_class();

			void init(int chan, int number_of_events, rootstuff *rt);
			void Setup(ColAHelL* MasterCH);
			void append_data(CH_event_struct *evnt);
			void FillPPCO(CH_event_struct *evnt);
			void Fill_wiggletool(CH_event_struct *evnt);
			void Fill_fishtool(CH_event_struct *evnt);
			void Fill_iontoftool(CH_event_struct *evnt);
			bool is_PPCO_tool_enabled();
			void enable_PPCO(int bins, double min, double max, bool momcon);
			void enable_LUT(double mom, int bins_mom, int bins_phi, int bins_ctheta, int electrons);
			void enable_wiggletool(int bins, double min, double max);
			void enable_fishtool(int bins, double min, double max);
			void enable_iontoftool(int bins, double min, double max);
			void add_par(double val, double minimum, double maximum, double steps, char *name, int n=-1);
			void set_fac_rot(int det, cor_param *fac);
			void set_mom_fac(int det, mom_cor_param *mom_fac, std::vector<double>  interpolation_adjustment_ceff = std::vector<double>(5,0) );
			void set_spec(spectrometer_class *Spect);
			void set_tof_calc_class(tof_calc_class *tof);
			void draw_screen();
			void init_interactive();
			int GetChannel();
			void main_loop();
			int events_in_memory();
			void set_reaction(reaction_struct *reac);

	};	// end ipa_class

#endif;
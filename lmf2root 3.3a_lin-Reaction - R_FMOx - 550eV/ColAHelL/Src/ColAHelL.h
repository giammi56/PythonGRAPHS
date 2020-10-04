#ifndef COLAHELL_ALREADY_INCLUDED
	#define COLAHELL_ALREADY_INCLUDED

	#pragma warning( disable : 4996 ) // disable strcpy-warning (doesn't work anymore?)
	#define _CRT_SECURE_NO_WARNINGS  // disable strcpy + sprintf warning
	#pragma warning( disable : 4018 ) // unsigned/signed-warning (it is not a really good idea to disable it)


	#include "CH_Basics.h"
	#include "CH_FUN_Lowlevel.h"
	#include "..\CH_Event.h"
	#include "CH_Reaction.h"
	#include "CH_Ranges.h"
	#include "CH_Vector.h"
	#include "CH_Coordinate_system.h"

	#include "CH_Spectrometer.h"
	#include "CH_Presorter.h"
	#include "CH_Tof.h"
	#include "CH_Electron.h"
	#include "CH_Ion.h"
	#include "CH_Diatomic.h"
	#include "CH_Polyatomic.h"
	#include "CH_Console.h"
	#include "CH_Hist.h"
	#include "CH_Processor.h"
	#include "CH_Histograms.h"

	#include <string>
	#include <vector>
	#include <fstream>

	#include "../RapidXML/rapidxml.hpp"
	using namespace rapidxml;
	using namespace std;
	using namespace CH;

	#define MAX_IONS 64
	#define MAX_ELECTRONS 64
	#define MAX_PROJECTILES 64

//	struct ColAHelL_cfg;
	
	struct Det_St;
	struct Momenta_St;
	struct Invalidate_if_St;
	struct Invalidate_ifnot_St;

	class ColAHelL {
	
		// those are needed e.g. by IPA, so make'em public..
		public:
			//ColAHelL_cfg *Cfg;

			spectrometer_class *Spect;
			tof_calc_class *CTof;
			processor_class * Proc;

		private:
			double CH_VERSION;
			xml_document<> doc;
			
			streambuf *psbuf, *backup;
			ofstream filestr;

			int Status;
			char CH_err[1024];

			presorter_class * presorters[64];
			int num_presorters;
			std::vector<CH_event_struct>  prs_evt_out;
			histo_handler * PrsHist;
			int maxNumPrsHists;

			histograms_class * HGrams;
			int maxNumCHHists;

			electron_class * e[64];
			ion_class * i[64];
			diatomic_class * mol;
			polyatomic_class * big_mol;
			reaction_struct * current_reaction;

			inline string ToUp(string Str) {for(int i=0;i < Str.size();i++){Str[i]=toupper(Str[i]);};return Str;};

			void CopyDetPara(xml_node<> *node, cor_param *e_fac);
			void CopyMomPara(xml_node<> *node, mom_cor_param *fac);

			void init_reac(reaction_struct *reac);
			void prepare_tag(int type, rtag_struct *temp_tag, xml_node<> *node, bool negate);
			//void prepare_tag(int type, rtag_struct *temp_tag, Invalidate_ifnot *Inv);
			bool OpenXML(char* FName);

			void SetupEventStruct();
			bool SetupTOFCalc();
			bool SetupSpectrometer();
			bool SetupDetectors();
			bool SetupMomentumParameters();
			bool SetupReactions();
			bool SetupPresorters();
			bool SetupHistogramming();
			bool SetupParticleClasses();
			void DisplayFortune();
			void DisplayConfucius();
			void Error(string s);

		public:

			ColAHelL();
			
			int InitAll();
			int InitAll(bool CrLogFile);
			int InitAll(const char* FName, bool CrLogFile);
			int ReplaceCHCfg(const char* FName);

			~ColAHelL();
			void SetVersion(double ver) {
				CH_VERSION = ver;
			}
			//current event ColAHelL is processing
			CH_event_struct * evt;
			void ResetEventStruct(CH_event_struct * evt, __int64 event_number);

			//calculate times-of-flight.
			void CalculateTOFs();
			
			// run all presorters on event and return a vector (array) of presorted events
			std::vector<CH_event_struct> PresortEvent();

			// check whether event belongs to presorter channel CHANNEL. 
			bool CheckEvent(CH_event_struct * evt, int channel);

			// calculate momenta etc.
			void ProcessEvent(CH_event_struct * evt, int ReactionDefNum);
			void ProcessEvent(int ReactionDefNum);
			
			// Fill std histograms. As an event can belong to several reactions we need to 
			// specify which reaction definition we want to use.
			void FillCHHists();
	
			bool IsReactionDefined(int reaction);
			string GetReactionName(int reaction);
			int EvtBelongsToReactions(int channel);

			cor_param* GetCorParam(int Det);
			mom_cor_param* GetMomCorParam(int Det);

			double GetPDetSize();
			double GetEDetSize();
			double GetIDetSize();
			int GetMaxPrsHists();
			histo_handler* GetPrsHists();
			int GetMaxUserHists();
			int GetMaxHists();
			histo_handler* GetHists(int Type);

			void SetMasterFolder(const char* folder);
	};

#endif
#pragma warning(disable : 4800)
#include "OS_Version.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtupleD.h"

#include <math.h>
#include "rootstuff.h"
#include "Histo.h"
#include "TF1.h"
//#include "TMinuit.h"

#include "functions.h"
#include "Ueberstruct.h"

//#include "polynomial_TOF_to_P.hpp"
//#include "resort64c.h"

#define AUTO __LINE__+1000 //returns the line number in the current file. Is used to make life easier for histogram numbering.

char ntuple_identifier_A[2048];
bool create_ntuple_identifier_A;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int analysis(__int64 eventcounter, double  parameter[], TTree * Data, Ueberstruct * Ueber)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	
	// do some renaming for less typing...
	Histo * Hist = Ueber->Hist;

	Ueber->start_new_root_file = false;

	// variable used when avoiding ColAHelL
	double reaction=-1., ehit=-1., rhit=-1., phit=-1., bunchmarker=-1., timestamp =-1., LMFStart=-1., LMFStop=-1.;
	double ex[64],ey[64],etof[64],etime[64],emethod[64];
	double rx[64],ry[64],rtof[64],rtime[64],rmethod[64];
	double px[64],py[64],ptof[64],ptime[64],pmethod[64];
	double scan_val[16];

	double pdet_size	= 5. + parameter[122];				// size of projectile detector (for displaying spectra properly)
	double rdet_size	= 5. + parameter[222];				// size of recoil detector (for displaying spectra properly)
	double edet_size	= 5. + parameter[322];				// size of electron detector (for displaying spectra properly)

	double NTupleData[250];
	bool WriteNTuple = false;
	char tmp[200];;

	int NRmax = (int)(parameter[61]+0.01);
	int NEmax = (int)(parameter[62]+0.01);
	int NPmax = (int)(parameter[60]+0.01);

	if (Ueber->error) {
		Ueber->stop_reading_files = true;
		return 0;
	}

	if(eventcounter == 0) {
		Ueber->EntriesInFile = 0;
		Ueber->eventswritten = 0;
		create_ntuple_identifier_A = true;
	}


	if(Ueber->EntriesInFile == 0) {

		int MaxRec = 0;
		int MaxElec = 0;
		int MaxProj = 0;

		char ntuple_identifier[500];
		bool Check = 0;
		// set branch addresses...

		Data->SetBranchAddress("reaction",&reaction);
		Data->SetBranchAddress("ehit",&ehit);
		Data->SetBranchAddress("rhit",&rhit);
		Data->SetBranchAddress("phit",&phit);
		Data->SetBranchAddress("bunchmarker",&bunchmarker);
		Data->SetBranchAddress("timestamp",&timestamp);
		Data->SetBranchAddress("LMFStart",&LMFStart);
		Data->SetBranchAddress("LMFStop",&LMFStop);


		// find out if/how many scan parameters are in Ntuple
		for(int i=15;i>=0;i--) {
			sprintf(ntuple_identifier,"scanval%i",i);
			if(Data->GetBranch(ntuple_identifier)) {
				Ueber->scan_in_data = true;
				Ueber->scan_num_vals = i+1;
				break;
			}
		}

		//check whether we have a scan value in the ntuple
		if(Ueber->scan_in_data) {
			for(int i=0;i<Ueber->scan_num_vals;i++) {
				sprintf(ntuple_identifier,"scanval%i",i);
				Data->SetBranchAddress(ntuple_identifier,&scan_val[i]);
			}
		}

		// find out how many recoils are in Ntuple
		do {
			sprintf(ntuple_identifier,"r%ix",++MaxRec);
		} while(Data->GetBranch(ntuple_identifier));
		MaxRec = (MaxRec<17 ? MaxRec : 16);
		Ueber->RecInNTuple = MaxRec;

		for(int i=0;i<MaxRec-1;i++) {
			sprintf(ntuple_identifier,"r%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rx[i]);
			sprintf(ntuple_identifier,"r%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ry[i]);
			sprintf(ntuple_identifier,"r%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rtime[i]);
			sprintf(ntuple_identifier,"r%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rtof[i]);
			sprintf(ntuple_identifier,"r%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rmethod[i]);
		}

		// find out how many electrons are in Ntuple
		do {
			sprintf(ntuple_identifier,"e%ix",++MaxElec);
		} while(Data->GetBranch(ntuple_identifier));
		MaxElec = (MaxElec<17 ? MaxElec : 16);
		Ueber->ElecInNTuple = MaxElec;

		for(int i=0;i<MaxElec-1;i++) {
			sprintf(ntuple_identifier,"e%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ex[i]);
			sprintf(ntuple_identifier,"e%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ey[i]);
			sprintf(ntuple_identifier,"e%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&etime[i]);
			sprintf(ntuple_identifier,"e%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&etof[i]);
			sprintf(ntuple_identifier,"e%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&emethod[i]);
		}

		// find out how many projectyles are in Ntuple
		do {
			sprintf(ntuple_identifier,"p%ix",++MaxProj);
		} while(Data->GetBranch(ntuple_identifier));
		MaxProj = (MaxProj<17 ? MaxProj : 16);
		Ueber->ProjInNTuple = MaxProj;
		
		for(int i=0;i<MaxProj-1;i++) {
			sprintf(ntuple_identifier,"p%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&px[i]);
			sprintf(ntuple_identifier,"p%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&py[i]);
			sprintf(ntuple_identifier,"p%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ptime[i]);
			sprintf(ntuple_identifier,"p%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ptof[i]);
			sprintf(ntuple_identifier,"p%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&pmethod[i]);
		}
	}

	Data->GetEntry(Ueber->EntriesInFile);

	if(Ueber->EntriesInFile < Data->GetEntries()-1) {
		++Ueber->EntriesInFile;
	} else {
		Ueber->EntriesInFile = 0;
	}

// ----------------------------------------------------------------------------------------------------------------------------
// If you do not want to use any ColAHelL things and work on the raw data, you may do that here.
//
// The raw data is stored in:
//
// electron: ehit, ex[64],ey[64],etof[64],etime[64],emethod[64];
// ion: rhit, rx[64],ry[64],rtof[64],rtime[64],rmethod[64];
// projectile: phit, px[64],py[64],ptof[64],ptime[64],pmethod[64];
// scanning_values: scan_val[16];
// ----------------------------------------------------------------------------------------------------------------------------

// A few very basic histograms
	Hist->fill(AUTO,"reaction_ntuple",reaction,1.,"Reaction channel flag (ntuple)",203,-1.25,100.25,"reaction_flag","ntuple");
	Hist->fill(AUTO,"rhit_ntuple",rhit,1.,"Recoil hits (ntuple)",43,-1.25,20.25,"rhit","ntuple");
	Hist->fill(AUTO,"phit_ntuple",phit,1.,"Projectile hits (ntuple)",43,-1.25,20.25,"phit","ntuple");
	Hist->fill(AUTO,"ehit_ntuple",ehit,1.,"Electron hits (ntuple)",43,-1.25,20.25,"ehit","ntuple");


// ----------------------------------------------------------------------------------------------------------------------------
// ColAHelL processing starts here. If you want to include some analysis and custom spectra,
// you can do that by changing the file "USER_Analysis.cpp" in the ColAHelL project. I do not 
// recommend to do that here... 

	if(Ueber->use_CH && Ueber->CH_status >-90) {
		CH_event_struct * CH_evt;

		double ds = -1.0;
		ds =  Ueber->CH->GetPDetSize();
		if(ds>0.0)
			pdet_size = ds + 5.0;
		ds =  Ueber->CH->GetIDetSize();
		if(ds>0.0)
			rdet_size = ds + 5.0;
		ds =  Ueber->CH->GetEDetSize();
		if(ds>0.0)
			edet_size = ds + 5.0;

		//copy data to ColAHelL
		Ueber->CH->ResetEvent(eventcounter);
		Ueber->CH->SetEvent_reaction(reaction);
		Ueber->CH->SetEvent_bm(bunchmarker);
		for(int i=0;i<ehit;i++)
			Ueber->CH->SetEvent_xyttof(EL,ex[i],ey[i],etime[i],etof[i],emethod[i]);
		for(int i=0;i<rhit;i++)
			Ueber->CH->SetEvent_xyttof(IO,rx[i],ry[i],rtime[i],rtof[i],rmethod[i]);
		for(int i=0;i<phit;i++)
			Ueber->CH->SetEvent_xyttof(PR,px[i],py[i],ptime[i],ptof[i],pmethod[i]);
		for(int i=0;i<Ueber->scan_num_vals;i++)
			Ueber->CH->SetEvent_scanvals(scan_val[i]);

		CH_evt = Ueber->CH->GetEvent();
// IPA needs food... Feeding is done here...
// -----------------------------------------------------------
		if(Ueber->ipa_enable) {
			Ueber->IPA->append_data(CH_evt);
			return 0;
		}

// Process the data and plot histograms. We might have several different
// reaction definitions belonging to the same "channel" (e.g. in case
// the "ANY" presorter was used), so we might need to loop over event 
// several times...

		// Let's see if we find a reaction definition for the current channel..
		// -----------------------------------------------------------	

		int EvtBelongsToReactions = Ueber->CH->EvtBelongsToReactions((int)CH_evt->reaction);
		// We need to loop, as there might be several reactions using the same channel number...
		// (e.g. when people use the ANY presorter...)
		for(int i=0;i<EvtBelongsToReactions;i++) {

			// ColAHelL kicks in here! ...
			// -----------------------------------------------------------
			// process particles...
			if(Ueber->CH_status>-40)
				Ueber->CH->ProcessEvent(i);
			// fill std. ColAHelL histograms
			if(Ueber->CH_status>-30)
				Ueber->CH->FillHistograms();

//			__int64 TotalNumEvents = (__int64)(Data->GetEntries()); //total number of events in the root file that is open

			if(WriteNTuple) {
				int NRec = (int)rhit;
				int NElec = (int)ehit;
				int NPro = (int)phit;
		
				if(NRec > NRmax)
					NRec = NRmax;
				if(NElec > NEmax)
					NElec = NEmax;
				if(NPro > NPmax)
					NPro = NPmax;

				NTupleData[0] = reaction;
				NTupleData[1] = double(NElec);
				NTupleData[2] = double(NRec);
				NTupleData[3] = double(NPro);
				NTupleData[4] = bunchmarker;
				NTupleData[5] = timestamp;
				NTupleData[6] = (double)Ueber->LMF_input->Starttime;
				NTupleData[7] = (double)Ueber->LMF_input->Stoptime; 

				int absind = 8;		

				if(Ueber->scan_in_data) {
					for(int i=0;i<Ueber->scan_num_vals;i++)
					{
						NTupleData[absind+i]=scan_val[i];
					}
					absind += Ueber->scan_num_vals;
				}

				for(int i=0;i<NRec;i++)
				{	// begin for: write all Recoils to NTuple
						NTupleData[absind+i*5] = CH_evt->r.x[i];
						NTupleData[absind+i*5+1] = CH_evt->r.y[i];
						NTupleData[absind+i*5+2] = CH_evt->r.time[i];
						NTupleData[absind+i*5+3] = CH_evt->r.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->r.method[i];
				}	// end for: write all Recoils to NTuple

				absind += NRmax*5;

				for(int i=0;i<NElec;i++)
				{	// begin for: write all Electrons to NTuple
						NTupleData[absind+i*5] = CH_evt->e.x[i];
						NTupleData[absind+i*5+1] = CH_evt->e.y[i];
						NTupleData[absind+i*5+2] = CH_evt->e.time[i];
						NTupleData[absind+i*5+3] = CH_evt->e.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->e.method[i];
				}	// end for: write all Electrons to NTuple

				absind += NEmax*5;

				for(int i=0;i<NPro;i++)
				{	// begin for: write all Projectiles to NTuple
						NTupleData[absind+i*5] = CH_evt->p.x[i];
						NTupleData[absind+i*5+1] = CH_evt->p.y[i];
						NTupleData[absind+i*5+2] = CH_evt->p.time[i];
						NTupleData[absind+i*5+3] = CH_evt->p.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->p.method[i];
				}	// end for: write all Projectiles to NTuple

			//////////////////////////////////////////////////////////////////////////
			//	write data to ntuple
			//////////////////////////////////////////////////////////////////////////

				if(create_ntuple_identifier_A) {
					create_ntuple_identifier_A = false;
			
					sprintf(ntuple_identifier_A,"reaction:ehit:rhit:phit:bunchmarker:timestamp:LMFStart:LMFStop");

					if(Ueber->scan_in_data) {

					for(int i=0;i<Ueber->scan_num_vals;i++)
					{
						sprintf(tmp,":scanval%d",i);
						strcat(ntuple_identifier_A,tmp);
					}
				}

				for(int i=1;i<=NRmax;i++)
				{
					sprintf(tmp,":r%dx:r%dy:r%dmcp:r%dtof:r%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
				}
				for(int i=1;i<=NEmax;i++)
				{
					sprintf(tmp,":e%dx:e%dy:e%dmcp:e%dtof:e%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
				}
				for(int i=1;i<=NPmax;i++)
				{
					sprintf(tmp,":p%dx:p%dy:p%dmcp:p%dtof:p%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
					}
				}

				Hist->NTupleD(0,"Data","BESSY2010",ntuple_identifier_A, 32000, NTupleData);
				++Ueber->eventswritten;

				if(parameter[57]>0.5) {
					unsigned __int64 max_events = (unsigned __int64)(parameter[56]+0.1);
					if(Ueber->eventswritten > (__int64)max_events && max_events > 0) {
						Ueber->start_new_root_file = true;
						Ueber->eventswritten = 0;
						Hist->Reset();
					}
				}
			}
		}
	} //end loop over event

	return 0;
}
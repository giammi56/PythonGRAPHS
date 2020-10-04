#include "OS_Version.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"

#include <math.h>
#include "rootstuff.h"
#include "Histo.h"
#include "TF1.h"
//#include "TMinuit.h"
#include "TNtupleD.h"

#include "functions.h"
#include "Ueberstruct.h"


//#include "ADC_analysis.h"
#include "../ADC/ADC_meta.h"

#include "console.h"

__int32 detector_stuff(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[]);

char ntuple_identifier[2048];
bool create_ntuple_identifier;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__int32 sort_and_write_NTuple(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[], double timestamp)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
#ifdef _DEBUG
	assert(Ueber);
#endif

	double NTupleData[250]; 
	double scan_val[16];

	bool writeNTuple = false;

	if (Ueber->error) {
		Ueber->stop_reading_files = true;
		return 0;
	}

	Ueber->write_this_event_into_sorted_LMF = false;
	Ueber->start_new_root_file = false;

	int NRec  = 0;
	int NElec = 0;
	int NPro  = 0;

	int NRmax = (int)(parameter[61]+0.01);
	int NEmax = (int)(parameter[62]+0.01);
	int NPmax = (int)(parameter[60]+0.01);

	if(eventcounter == 0) {
		Ueber->eventswritten = 0;
		create_ntuple_identifier = true;
	}

	Histo * Hist = Ueber->Hist;

	if (!detector_stuff(eventcounter,Ueber,parameter)) return 0;

	double * tdc_ns = Ueber->tdc_ns; // this line must be below   if (!detectorstuff()...
	__int32 * cnt = Ueber->cnt; // this line must be below   if (!detectorstuff()...

	double bunchspacing = 0.;
	int bunchmarker_channel = -1;
	bool useBM = false;

	double pdet_size	= 5. + parameter[122];				// size of projectile detector (for displaying spectra properly)
	double rdet_size	= 5. + parameter[222];				// size of recoil detector (for displaying spectra properly)
	double edet_size	= 5. + parameter[322];				// size of electron detector (for displaying spectra properly)

															// Scan data is encoded in TDC channel
	if(Ueber->scan_channel>-1 && Ueber->scan_channel<999) {
		Ueber->scan_in_data = true;
		for(int i=0;i<Ueber->scan_num_vals;i++)
			scan_val[i] = (double)tdc_ns[Ueber->scan_channel*NUM_IONS+i]/(double)Ueber->scan_factor;
	}
	// Scan data is encoded in LMF parameters
	if(Ueber->scan_channel == 999) {
		Ueber->scan_in_data = true;
		for(int i=0;i<Ueber->scan_num_vals;i++)
			scan_val[i] = Ueber->LMF_input->Parameter[Ueber->scan_par_nums[i]];
	}

	// Include your analysis of the TDC data here.
	//------------------------------------------------------------------------------------------------
	// Access the results from detector reconstruction using the detector struct.
	// The name of the struct is "proj", "rec" or "elec" 
	// depending on the detector. The struct is created as a part of "Ueber". 
	// Thus it is addressed using "Ueber->".
	//
	// For example the electron detector is accessed by "Ueber->elec".
	//
	// See Ueberstruct.h for the definition of the detector struct. However, the following 
	// are the important variables:
	//
	//		__int32		number_of_reconstructed_hits	- number of hits or this event
	//		__int32		method[NUM_HIT]					- reconstruction method that was used
	//		double	x[NUM_HIT]						- x position of impact
	//		double	y[NUM_HIT]						- y position of impact
	//		double	time[NUM_HIT]					- MCP time 
	//
	// The first hit is NUM_HIT=0.
	//
	// Another example: in order to write the x-position of the second hit on a 
	// recoil detector to xpos use:
	//		double xpos = Ueber->rec->x[1];
	//
	//------------------------------------------------------------------------------------------------



	// Feed things into ColAHelL if enabled...	
	// ("evt->"-struct is only available when using ColAHelL!)
	// -----------------------------------------------------------

	CH::CH_event_struct* evt = 0;
	if(Ueber->use_CH && Ueber->CH_status > -90) {

		useBM = Ueber->CH->GetUseBM();		
		if(useBM) {
			bunchspacing = Ueber->CH->GetBunchSpacing();
			bunchmarker_channel = Ueber->CH->GetBMChannel();
		}
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

		// Limit data to max. 64 hits 
		Ueber->elec->number_of_reconstructed_hits = (Ueber->elec->number_of_reconstructed_hits<65) ? Ueber->elec->number_of_reconstructed_hits : 64;
		Ueber->rec->number_of_reconstructed_hits = (Ueber->rec->number_of_reconstructed_hits<65) ? Ueber->rec->number_of_reconstructed_hits : 64;
		Ueber->proj->number_of_reconstructed_hits = (Ueber->proj->number_of_reconstructed_hits<65) ? Ueber->proj->number_of_reconstructed_hits : 64;

		// copy event data to ColAHelL
		Ueber->CH->ResetEvent(eventcounter);
		for(int i=0;i<Ueber->elec->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(EL, Ueber->elec->x[i], Ueber->elec->y[i], Ueber->elec->time[i], (double)Ueber->elec->method[i]); 	
		for(int i=0;i<Ueber->rec->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(IO, Ueber->rec->x[i], Ueber->rec->y[i], Ueber->rec->time[i], (double)Ueber->rec->method[i]); 
		for(int i=0;i<Ueber->proj->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(PR, Ueber->proj->x[i], Ueber->proj->y[i], Ueber->proj->time[i], (double)Ueber->proj->method[i]); 
		// get the bunchmarker time from the tdc array...
		Ueber->CH->SetBMTime(tdc_ns, cnt);
		// copy scan data to ColAHelL
		for(int i=0;i<Ueber->scan_num_vals;i++)
			Ueber->CH->SetEvent_scanvals(scan_val[i]);

		// calculate TOFs using ColAHelL definitions
		Ueber->CH->CalculateTOFs();
		evt = Ueber->CH->GetEvent();

		// IPA needs food, as well... Feeding is done here...
		// (IPA works only together with ColAHelL!)
		// -----------------------------------------------------------
		if(Ueber->ipa_enable && Ueber->IPA) { 
			if(Ueber->IPA->GetChannel() == -1) {
				Ueber->IPA->append_data(evt);
				return 0;
			} else {
				std::vector<CH::CH_event_struct> IPAPrsEvts = Ueber->CH->PresortEvent();
				for(int i=0;i<IPAPrsEvts.size();i++) {
					Ueber->IPA->append_data(&IPAPrsEvts[i]);
				}
			}
			return 0;
		}
	}

	//---------------------------------------------------------------------------------------------------------------------------
	// some typical spectra ...
	//---------------------------------------------------------------------------------------------------------------------------



	if(Ueber->scan_channel>-1 && Ueber->scan_channel<999) {
		Hist->fill1(108,"scan_hits",cnt[Ueber->scan_channel],1.,"number of hits scan channel",1000,-1.,25.,"cnt","raw");
		Hist->fill1(109,"scan value",scan_val[0],1.,"scan value",1000,0.,100.,"cnt","raw");
	}

	if(Ueber->use_CH) {
		if(useBM)
		{
			Hist->fill1(100,"bunchmarker1",tdc_ns[bunchmarker_channel*NUM_IONS+0],1.,"Time Bunchmarker #1",10000,-10000.,10000.,"bunchmarker [ns]","raw");
			Hist->fill1(101,"bunchmarker_hit",cnt[bunchmarker_channel],1.,"Bunchmarker hits",63,-1.25,30.25,"bunchmarker hits","raw");
		}

		Hist->fill1(105,"phit",evt->p.num_hits,1.,"projectile number of reconstructed hits",1000,-1.,25.,"evt->p.num_hits","raw");
		Hist->fill1(106,"ehit",evt->e.num_hits,1.,"electron number of reconstructed hits",1000,-1.,25.,"evt->e.num_hits","raw");
		Hist->fill1(107,"rhit",evt->r.num_hits,1.,"recoil number of reconstructed hits",1000,-1.,25.,"evt->r.num_hits","raw");

		if (evt->p.num_hits>0 )
		{	// begin if: one projectile
			Hist->fill2(110,"proj1xy",evt->p.x[0],evt->p.y[0],1.,"Projectile position #1",400,-1.*pdet_size,pdet_size,"proj1 x [mm]",400,-1.*pdet_size,pdet_size,"proj1 y [mm]","raw");
		}	// end if: one projectile


		if (evt->r.num_hits>0 )
		{	//  begin if: one recoil
			Hist->fill2(140,"rec1xy",evt->r.x[0],evt->r.y[0],1.,"Recoil position #1",400,-1.*rdet_size,rdet_size,"rec1 x [mm]",400,-1.*rdet_size,rdet_size,"rec1 y [mm]","raw");

			Hist->fill1(145,"rec1tof",evt->r.tof[0],1.,"Recoil TOF #1",10000,-10000.,60000.,"rec1 TOF [ns]","raw");
			Hist->fill1(1451,"rec1tof_short",evt->r.tof[0],1.,"Recoil TOF #1",10000,-2500.,7500.,"rec1 TOF [ns]","raw");
			Hist->fill2(146,"rec1wiggle",evt->r.tof[0],(sqrt(evt->r.x[0]*evt->r.x[0]+evt->r.y[0]*evt->r.y[0])),1.,"Recoil wiggle #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.,rdet_size,"rec1 radius [mm]","raw");
			Hist->fill2(147,"rec1xfish",evt->r.tof[0],evt->r.x[0],1.,"Recoil x-fish #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 x-pos [mm]","raw");
			Hist->fill2(1471,"rec1xfish_short",evt->r.tof[0],evt->r.x[0],1.,"Recoil x-fish #1",500,-2500.,7500.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 x-pos [mm]","raw");
			Hist->fill2(148,"rec1yfish",evt->r.tof[0],evt->r.y[0],1.,"Recoil y-fish #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 y-pos [mm]","raw");
			Hist->fill2(1481,"rec1yfish_short",evt->r.tof[0],evt->r.y[0],1.,"Recoil y-fish #1",500,-2500.,7500.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 y-pos [mm]","raw");
		}	//  end if: one recoil

		if (evt->r.num_hits>1 )
		{	//  begin if: two recoils
			Hist->fill2(150,"rec2xy",evt->r.x[1],evt->r.y[1],1.,"Recoil position #2",400,-1.*rdet_size,rdet_size,"rec2 x [mm]",400,-1.*rdet_size,rdet_size,"rec2 y [mm]","raw");

			Hist->fill1(155,"rec2tof",evt->r.tof[1],1.,"Recoil TOF #2",1000,-10000.,60000.,"rec2 TOF [ns]","raw");
			Hist->fill1(1551,"rec2tof_short",evt->r.tof[1],1.,"Recoil TOF #2",1000,-2500.,7500.,"rec2 TOF [ns]","raw");
			Hist->fill2(156,"rec2wiggle",evt->r.tof[1],(sqrt(evt->r.x[1]*evt->r.x[1]+evt->r.y[1]*evt->r.y[1])),1.,"Recoil wiggle #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.,rdet_size,"rec2 radius [mm]","raw");
			Hist->fill2(157,"rec2xfish",evt->r.tof[1],evt->r.x[1],1.,"Recoil x-fish #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 x-pos [mm]","raw");
			Hist->fill2(1571,"rec2xfish_short",evt->r.tof[1],evt->r.x[1],1.,"Recoil x-fish #2",500,-2500.,7500.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 x-pos [mm]","raw");
			Hist->fill2(158,"rec2yfish",evt->r.tof[1],evt->r.y[1],1.,"Recoil y-fish #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 y-pos [mm]","raw");
			Hist->fill2(1581,"rec2yfish_short",evt->r.tof[1],evt->r.y[1],1.,"Recoil y-fish #2",500,-2500.,7500.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 y-pos [mm]","raw");

			Hist->fill2(160,"recsumdiff",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff",1000,10000.,40000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1601,"recsumdiff_short",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff",1000,-5000.,15000.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(161,"recsumdiff2",2.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (2*r1 + r2)",800,0.,50000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1611,"recsumdiff2_short",2.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (2*r1 + r2)",800,-7500.,22500.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(162,"recsumdiff3",3.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (3*r1 + r2)",800,0.,60000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1621,"recsumdiff3_short",3.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (3*r1 + r2)",800,-10000.,30000.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");

			Hist->fill2(163,"pipico",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum",500,0.,50000.,"rec1 TOF [ns]",500,0.,50000.,"rec2 TOF [ns]","raw");
			Hist->fill2(1631,"pipico_short",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum",800,-2500.,7500.,"rec1 TOF [ns]",800,-2500.,7500.,"rec2 TOF [ns]","raw");

			double phi0 = atan2(evt->r.x[0],evt->r.y[0])*180./3.14152;
			double phi1 = atan2(evt->r.x[1],evt->r.y[1])*180./3.14152;

			if(fabs(deltaphi(phi0,phi1)-180.0)<10.0)
			{
				Hist->fill2(165,"pipico_momcon",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum (coarse) /w simple mom. check",500,0.,50000.,"rec1 TOF [ns]",500,0.,50000.,"rec2 TOF [ns]","raw");
				Hist->fill2(1651,"pipico_momcon_short",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum (fine) /w simple mom. check",800,-2500.,7500.,"rec1 TOF [ns]",800,-2500.,7500.,"rec2 TOF [ns]","raw");
			}
		}	//  end if: two recoils


		if ( evt->r.num_hits>2)
		{	// begin if: three recoils
			Hist->fill2(170,"rec3xy",evt->r.x[2],evt->r.y[2],1.,"Recoil position #3",400,-1.*rdet_size,rdet_size,"rec3 x [mm]",400,-1.*rdet_size,rdet_size,"rec3 y [mm]","raw");

			Hist->fill1(175,"rec3tof",evt->r.tof[2],1.,"Recoil TOF #3",1000,-10000.,60000.,"rec3 TOF [ns]","raw");
			Hist->fill1(1751,"rec3tof_short",evt->r.tof[2],1.,"Recoil TOF #3",1000,-2500.,7500.,"rec3 TOF [ns]","raw");
			Hist->fill2(176,"rec3wiggle",evt->r.tof[2],(sqrt(evt->r.x[2]*evt->r.x[2]+evt->r.y[2]*evt->r.y[2])),1.,"Recoil wiggle #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.,rdet_size,"rec3 radius [mm]","raw");
			Hist->fill2(177,"rec3xfish",evt->r.tof[2],evt->r.x[2],1.,"Recoil x-fish #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 x-pos [mm]","raw");
			Hist->fill2(1771,"rec3xfish_short",evt->r.tof[2],evt->r.x[2],1.,"Recoil x-fish #3",500,-2500.,7500.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 x-pos [mm]","raw");
			Hist->fill2(178,"rec3yfish",evt->r.tof[2],evt->r.y[2],1.,"Recoil y-fish #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 y-pos [mm]","raw");
			Hist->fill2(1781,"rec3yfish_short",evt->r.tof[2],evt->r.y[2],1.,"Recoil y-fish #3",500,-2500.,7500.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 y-pos [mm]","raw");
			Hist->fill2(179,"pipipico",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[2],1.,"PIPIPICO spectrum (coarse)",800,0.,80000.,"rec1+rec2 TOF [ns]",500,0.,50000.,"rec3 TOF [ns]","raw");
			Hist->fill2(1791,"pipipico_short",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[2],1.,"PIPIPICO spectrum (fine)",500,-5000.,15000.,"rec1+rec2 TOF [ns]",500,-2500.,7500.,"rec3 TOF [ns]","raw");
		}	// end if: three recoils


		if ( evt->e.num_hits>0 )
		{	// begin if: one electron
			Hist->fill2(180,"elec1xy",evt->e.x[0],evt->e.y[0],1.,"Electron position #1",400,-1.*edet_size,edet_size,"elec x [mm]",400,-1.*edet_size,edet_size,"elec y [mm]","raw");

			if(useBM && (bunchspacing > 0.0))
			{
				Hist->fill1(185,"elec1tof",evt->e.tof[0],1.,"Electron TOF #1",int(6*bunchspacing),-1.5*bunchspacing,1.5*bunchspacing,"elec1 TOF [ns]","raw");
				Hist->fill2(186,"elec1wiggle",evt->e.tof[0],(sqrt(evt->e.x[0]*evt->e.x[0]+evt->e.y[0]*evt->e.y[0])),1.,"Electron wiggle #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.,edet_size,"elec1 radius [mm]","raw");
				Hist->fill2(187,"elec1xfish",evt->e.tof[0],evt->e.x[0],1.,"Electron x-fish #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 x-pos [mm]","raw");
				Hist->fill2(188,"elec1yfish",evt->e.tof[0],evt->e.y[0],1.,"Electron y-fish #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 y-pos [mm]","raw");
			}
			else 
			{
				Hist->fill1(185,"elec1tof",evt->e.tof[0],1.,"Electron TOF #1",5000,-2500.,2500.,"elec1 TOF [ns]","raw");
				Hist->fill2(186,"elec1wiggle",evt->e.tof[0],(sqrt(evt->e.x[0]*evt->e.x[0]+evt->e.y[0]*evt->e.y[0])),1.,"Electron wiggle #1",500,0.,250.,"elec1 TOF [ns]",200,-1.,edet_size,"elec1 radius [mm]","raw");
				Hist->fill2(187,"elec1xfish",evt->e.tof[0],evt->e.x[0],1.,"Electron x-fish #1",500,0.,250.,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 x-pos [mm]","raw");
				Hist->fill2(188,"elec1yfish",evt->e.tof[0],evt->e.y[0],1.,"Electron y-fish #1",500,0.,250.,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 y-pos [mm]","raw");
			}

		}	// end if: one electron

		if ( evt->e.num_hits>1 )
		{	// begin if: two electrons
			Hist->fill2(190,"elec2xy",evt->e.x[1],evt->e.y[1],1.,"Electron position #2",400,-1.*edet_size,edet_size,"elec2 x [mm]",400,-1.*edet_size,edet_size,"elec2 y [mm]","raw");
			double dr = sqrt((evt->e.x[0] - evt->e.x[1])*(evt->e.x[0] - evt->e.x[1]) + (evt->e.y[0] - evt->e.y[1])*(evt->e.y[0] - evt->e.y[1]));
			//Hist->fill2(199, "dtof_vs_dr", fabs(evt->e.tof[0] - evt->e.tof[1]),dr,1.,"dToF vs dR",100,-50,50, )


			if(useBM && (bunchspacing > 0.0))
			{
				Hist->fill1(195,"elec2tof",evt->e.tof[1],1.,"Electron TOF #2",int(6*bunchspacing),-1.5*bunchspacing,1.5*bunchspacing,"elec2 TOF [ns]","raw");
				Hist->fill2(196,"elec2wiggle",evt->e.tof[1],(sqrt(evt->e.x[1]*evt->e.x[1]+evt->e.y[1]*evt->e.y[1])),1.,"Electron wiggle #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.,edet_size,"elec2 radius [mm]","raw");
				Hist->fill2(197,"elec2xfish",evt->e.tof[1],evt->e.x[1],1.,"Electron x-fish #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 x-pos [mm]","raw");
				Hist->fill2(198,"elec2yfish",evt->e.tof[1],evt->e.y[1],1.,"Electron y-fish #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 y-pos [mm]","raw");
			}
			else
			{
				Hist->fill1(195,"elec2tof",evt->e.tof[1],1.,"Electron TOF #2",2000,-1000.,1000.,"elec2 TOF [ns]","raw");
				Hist->fill2(196,"elec2wiggle",evt->e.tof[1],(sqrt(evt->e.x[1]*evt->e.x[1]+evt->e.y[1]*evt->e.y[1])),1.,"Electron wiggle #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.,edet_size,"elec2 radius [mm]","raw");
				Hist->fill2(197,"elec2xfish",evt->e.tof[1],evt->e.x[1],1.,"Electron x-fish #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 x-pos [mm]","raw");
				Hist->fill2(198,"elec2yfish",evt->e.tof[1],evt->e.y[1],1.,"Electron y-fish #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 y-pos [mm]","raw");
			}
		}	// end if: two electrons

		if (evt->e.num_hits == 2) {
			double dpos_mm = sqrt((evt->e.x[0] - evt->e.x[1])*(evt->e.x[0] - evt->e.x[1]) + (evt->e.y[0] - evt->e.y[1])*(evt->e.y[0] - evt->e.y[1]));
			double dtof_ns = fabs(evt->e.time[0] - evt->e.time[1]);

			Hist->fill1(80, "dpos", dpos_mm, 1., "dpos_mm", 2000, -50., 50., "dpos [mm]", "DPA");
			Hist->fill1(81, "dtof", dtof_ns, 1., "dtof_mm", 2000, 0., 1000., "dtof [ns]", "DPA");
			Hist->fill2(82, "dtof_vs_dpos", dtof_ns, fabs(dpos_mm), 1., "dtof_vs_dpos", 250, 0, 1000., "dtof [ns]", 200, 0, 50, "dpos [mm]", "DPA");
		}
	} // end if: "use CH"

	ADC_meta * adc = (ADC_meta*)(Ueber->adc_meta);

	if (adc->m_event->channels[0]) {
		if (adc->m_event->channels[0]->NbrPulses) {
			if (adc->m_event->channels[0]->Pulses[0].NbrPeaks) {
				//...
			}
		}
	}

	//---------------------------------------------------------------------------------------------------------------------------
	//	                            P R E S O R T E R                                      
	//---------------------------------------------------------------------------------------------------------------------------
	// Presorter definitions are part of ColAHelL!
	if(Ueber->CH && Ueber->CH_status > -50) {
		// Dublicate event in case it fits to several presorters... (default: true)
		bool dublicate = true;

		// presort event
		std::vector<CH::CH_event_struct> PrsEvts = Ueber->CH->PresortEvent();

		// write resulting data to file 
		for(unsigned int prsnum=0;prsnum<PrsEvts.size();prsnum++) {

			writeNTuple = true;

			NRec = PrsEvts[prsnum].r.num_hits;
			if(NRec > NRmax)
				NRec = NRmax;

			NElec = PrsEvts[prsnum].e.num_hits;
			if(NElec > NEmax)
				NElec = NEmax;

			NPro = PrsEvts[prsnum].p.num_hits;
			if(NPro > NPmax)
				NPro = NPmax;

			NTupleData[0] = double(PrsEvts[prsnum].reaction);
			NTupleData[1] = double(NElec);
			NTupleData[2] = double(NRec);
			NTupleData[3] = double(NPro);
			NTupleData[4] = PrsEvts[prsnum].bunchmarker;
			NTupleData[5] = timestamp;
			NTupleData[6] = (double)Ueber->LMF_input->Starttime;
			NTupleData[7] = (double)Ueber->LMF_input->Stoptime; 

			int absind = 8;

			if(Ueber->scan_channel>-1) {
				for(int i=0;i<Ueber->scan_num_vals;i++)
				{
					NTupleData[absind+i]=scan_val[i];
				}
				absind += Ueber->scan_num_vals;
			}

			for(int i=0;i<NRec;i++)
			{	
				NTupleData[absind+i*5]   = double(PrsEvts[prsnum].r.x[i]);
				NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].r.y[i]);
				NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].r.time[i]);
				NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].r.tof[i]);
				NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].r.method[i]);
			}	
			absind += NRmax*5;

			for(int i=0;i<NElec;i++)
			{	
				NTupleData[absind+i*5]   = double(PrsEvts[prsnum].e.x[i]);
				NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].e.y[i]);
				NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].e.time[i]);
				NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].e.tof[i]);
				NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].e.method[i]);
			}	
			absind += NEmax*5;

			for(int i=0;i<NPro;i++)
			{	
				NTupleData[absind+i*5]   = double(PrsEvts[prsnum].p.x[i]);
				NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].p.y[i]);
				NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].p.time[i]);
				NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].p.tof[i]);
				NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].p.method[i]);
			}	
			absind += NPmax*5;

			//////////////////////////////////////////////////////////////////////////
			//	write data to ntuple
			//////////////////////////////////////////////////////////////////////////

			char tmp[200];

			if(create_ntuple_identifier) {
				create_ntuple_identifier = false;

				sprintf(ntuple_identifier,"reaction:ehit:rhit:phit:bunchmarker:timestamp:LMFStart:LMFStop");

				if(Ueber->scan_channel>-1) {
					for(int i=0;i<Ueber->scan_num_vals;i++)
					{
						sprintf(tmp,":scanval%d",i);
						strcat(ntuple_identifier,tmp);
					}
				}

				for(int i=1;i<=NRmax;i++)
				{
					sprintf(tmp,":r%dx:r%dy:r%dmcp:r%dtof:r%dflag",i,i,i,i,i);
					strcat(ntuple_identifier,tmp);
				}
				for(int i=1;i<=NEmax;i++)
				{
					sprintf(tmp,":e%dx:e%dy:e%dmcp:e%dtof:e%dflag",i,i,i,i,i);
					strcat(ntuple_identifier,tmp);
				}
				for(int i=1;i<=NPmax;i++)
				{
					sprintf(tmp,":p%dx:p%dy:p%dmcp:p%dtof:p%dflag",i,i,i,i,i);
					strcat(ntuple_identifier,tmp);
				}
			}

			Hist->NTupleD(0,"Data","BESSY2012",ntuple_identifier, 32000, NTupleData);
			++Ueber->eventswritten;

			// plot some histogramms of all presorted events
			Hist->fill1(200,"reaction_ntuple",PrsEvts[prsnum].reaction,1.,"Reaction channel flag",203,-1.25,100.25,"reaction_flag","ntuple");

			Hist->fill1(205,"phit",PrsEvts[prsnum].p.num_hits,1.,"projectile number of reconstructed hits",1000,-1.,25.,"evt->p.num_hits","ntuple");
			Hist->fill1(206,"ehit",PrsEvts[prsnum].e.num_hits,1.,"electron number of reconstructed hits",1000,-1.,25.,"evt->e.num_hits","ntuple");
			Hist->fill1(207,"rhit",PrsEvts[prsnum].r.num_hits,1.,"recoil number of reconstructed hits",1000,-1.,25.,"evt->r.num_hits","ntuple");


			if (NElec>0)
			{	
				Hist->fill2(210,"elec1xy_ntuple",PrsEvts[prsnum].e.x[0],PrsEvts[prsnum].e.y[0],1.,"Electron position #1 (ntuple)",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[0] [mm]",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
				Hist->fill1(215,"elec1tof_ntuple",PrsEvts[prsnum].e.tof[0],1.,"Electron TOF #1 (ntuple)",1000,-500.0,500.0,"PrsEvts[prsnum].e.tof[0] [ns]","ntuple");
				Hist->fill2(216,"elec1wiggle_ntuple",PrsEvts[prsnum].e.tof[0],(sqrt(PrsEvts[prsnum].e.x[0]*PrsEvts[prsnum].e.x[0]+PrsEvts[prsnum].e.y[0]*PrsEvts[prsnum].e.y[0])),1.,"Electron wiggle #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.,edet_size,"elec_radius [mm]","ntuple");
				Hist->fill2(217,"elec1xfish_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.x[0],1.,"Electron x-fish #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[0] [mm]","ntuple");
				Hist->fill2(218,"elec1yfish_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.y[0],1.,"Electron y-fish #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
				if ( fabs(PrsEvts[prsnum].e.x[0])<10. ) Hist->fill2(219,"elec1yfishcut_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.y[0],1.,"Electron y-fish (x<10) #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",20,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
			}	

			if (NElec>1)
			{	
				Hist->fill2(220,"elec2xy_ntuple",PrsEvts[prsnum].e.x[1],PrsEvts[prsnum].e.y[1],1.,"Electron position #2 (ntuple)",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[1] [mm]",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[1] [mm]","ntuple");
				Hist->fill1(225,"elec2tof_ntuple",PrsEvts[prsnum].e.tof[1],1.,"Electron TOF #2 (ntuple)",1000,-500.,500.,"PrsEvts[prsnum].e.tof[1] [ns]","ntuple");
				Hist->fill2(226,"elec2wiggle_ntuple",PrsEvts[prsnum].e.tof[1],(sqrt(PrsEvts[prsnum].e.x[1]*PrsEvts[prsnum].e.x[1]+PrsEvts[prsnum].e.y[1]*PrsEvts[prsnum].e.y[1])),1.,"Electron wiggle #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.,edet_size,"elec2_radius [mm]","ntuple");
				Hist->fill2(227,"elec2xfish_ntuple",PrsEvts[prsnum].e.tof[1],PrsEvts[prsnum].e.x[1],1.,"Electron x-fish #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[1] [mm]","ntuple");
				Hist->fill2(228,"elec2yfish_ntuple",PrsEvts[prsnum].e.tof[1],PrsEvts[prsnum].e.y[1],1.,"Electron y-fish #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[1] [mm]","ntuple");
			}	

			if (NRec>0)
			{	
				Hist->fill2(250,"rec1xy_ntuple",PrsEvts[prsnum].r.x[0],PrsEvts[prsnum].r.y[0],1.,"Recoil position #1 (ntuple)",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[0] [mm]",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[0] [mm]","ntuple");
				Hist->fill1(255,"rec1tof_ntuple",PrsEvts[prsnum].r.tof[0],1.,"Recoil TOF #1 (ntuple)",30100,-100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]","ntuple");
				Hist->fill2(256,"rec1wiggle_ntuple",PrsEvts[prsnum].r.tof[0],(sqrt(PrsEvts[prsnum].r.x[0]*PrsEvts[prsnum].r.x[0]+PrsEvts[prsnum].r.y[0]*PrsEvts[prsnum].r.y[0])),1.,"Recoil wiggle #1 (ntuple)",500,10000.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.,rdet_size,"rec1_radius [mm]","ntuple");
				Hist->fill2(257,"rec1xfish_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.x[0],1.,"Recoil x-fish #1 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[0] [mm]","ntuple");
				Hist->fill2(258,"rec1yfish_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.y[0],1.,"Recoil y-fish #1 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[0] [mm]","ntuple");
				Hist->fill2(259,"rec1tof_flag_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].reaction,1.,"Recoil TOF Reaction channel #1 (ntuple)",100,5000.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",405,-1.125,100.125,"reaction_flag","ntuple");
			}	

			if (NRec>1)
			{	
				Hist->fill2(260,"rec2xy_ntuple",PrsEvts[prsnum].r.x[1],PrsEvts[prsnum].r.y[1],1.,"Recoil position #2 (ntuple)",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[1] [mm]",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[1] [mm]","ntuple");
				Hist->fill1(265,"rec2tof_ntuple",PrsEvts[prsnum].r.tof[1],1.,"Recoil TOF #2 (ntuple)",30100,-100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
				Hist->fill2(266,"rec2wiggle_ntuple",PrsEvts[prsnum].r.tof[1],(sqrt(PrsEvts[prsnum].r.x[1]*PrsEvts[prsnum].r.x[1]+PrsEvts[prsnum].r.y[1]*PrsEvts[prsnum].r.y[1])),1.,"Recoil wiggle #2 (ntuple)",500,10000.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.,rdet_size,"rec2_radius [mm]","ntuple");
				Hist->fill2(267,"rec2xfish_ntuple",PrsEvts[prsnum].r.tof[1],PrsEvts[prsnum].r.x[1],1.,"Recoil x-fish #2 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[1] [mm]","ntuple");
				Hist->fill2(268,"rec2yfish_ntuple",PrsEvts[prsnum].r.tof[1],PrsEvts[prsnum].r.y[1],1.,"Recoil y-fish #2 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[1] [mm]","ntuple");
				Hist->fill2(270,"pipico_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.tof[1],1.,"PIPICO spectrum (ntuple)",500,0.,50000.,"PrsEvts[prsnum].r.tof[0] [ns]",500,0.,50000.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
				Hist->fill2(271,"pipico_short_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.tof[1],1.,"PIPICO spectrum fine (ntuple)",800,-2500.,7500.,"PrsEvts[prsnum].r.tof[0] [ns]",800,-2500.,7500.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
			}		

			if (NPro>0)
			{	
				Hist->fill2(290,"proj1xy_ntuple",PrsEvts[prsnum].p.x[0],PrsEvts[prsnum].p.y[0],1.,"Projectile position #1 (ntuple)",400,-1.*pdet_size,pdet_size,"PrsEvts[prsnum].p.x[0] [mm]",400,-1.*pdet_size,pdet_size,"PrsEvts[prsnum].p.y[0] [mm]","ntuple");
			}	

			// ColAHelL momentum stuff kicks in here! 
			// (comment in if you don't mind processing time and need momenta etc. already at this point...)
			// -----------------------------------------------------------
			/*
			if(Ueber->use_CH && Ueber->CH_status > -40)  {
			int EvtBelongsToReactions = Ueber->CH->EvtBelongsToReactions((int)PrsEvts[prsnum].reaction);
			// We need to loop, as there might be several reactions using the same channel number...
			// (e.g. when people use the ANY presorter...)
			for(int i=0;i<EvtBelongsToReactions;i++) {
			// process particles...
			Ueber->CH->ProcessEvent(&PrsEvts[prsnum],i);
			// fill std. ColAHelL histograms
			if(Ueber->CH_status >-30)
			Ueber->CH->FillHistograms();
			}
			}
			*/
			if(dublicate == false)
				break;	
		}

		if(parameter[57]>0.5) {
			unsigned __int32 max_events = (unsigned __int32)(parameter[56]+0.1);
			if(Ueber->eventswritten > max_events && max_events > 0) {
				Ueber->start_new_root_file = true;
				Ueber->eventswritten = 0;
			}
		}
	} // end if: ist ColAHelL enabled/defined?
	return 1;
}

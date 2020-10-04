#include "OS_Version.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TProfile.h"
#include "TCutG.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TMarker.h"
#include "ipa.h"
#include "./src/console.h"
#include <omp.h>

__int32 my_kbhit();

ipa_class::ipa_class(){}

ipa_class::~ipa_class(){
	if(!reaction_def_valid)
	{
		delete this->cur_reaction->e_fac;
		delete this->cur_reaction->e_mom_fac;
		delete this->cur_reaction->r_fac;
		delete this->cur_reaction->r_mom_fac;
		delete this->cur_reaction->p_fac;
		delete this->cur_reaction->p_mom_fac;
		delete this->cur_reaction;
	}

	if(this->fishtool_enabled) {
		for(int j=0;j<3;j++) { 
			for(int i=0;i<36;i++) {
				if(mrk[0][j][i])
					delete mrk[0][j][i];
				if(mrk[1][j][i])
					delete mrk[1][j][i];
				if(mrk[2][j][i])
					delete mrk[2][j][i];
			}
		}
		for(int i=0;i<128;i++) {
			if(knots[i])
				delete knots[i];
			if(firstControlPoints[i])
				delete firstControlPoints[i]; 
			if(secondControlPoints[i])
				delete secondControlPoints[i];
		}
	}
}

void ipa_class::init(int chan, int number_of_events, rootstuff *rt)
{

#ifdef _DEBUG
//	printf("IPA: USING MEGA SLOW DEBUGGING MODE !!!\n");
//	_CrtSetDbgFlag(_CrtSetDbgFlag( _CRTDBG_REPORT_FLAG ) | _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_CRT_DF | _CRTDBG_CHECK_ALWAYS_DF);
#endif

	// Thread handling...
	SYSTEM_INFO sysinfo;
	GetSystemInfo( &sysinfo );
	cores_in_use = sysinfo.dwNumberOfProcessors;

	if(cores_in_use>MAXTHREADS)
		cores_in_use = MAXTHREADS;

	omp_set_num_threads(cores_in_use);
	this->precision = 1.0;

	this->channel = chan;
	this->num_events  = number_of_events;
	this->num_events_read = 0;
	this->page = EDETPAGE;

	this->NumPars[EDETPAGE]= 8;
	this->NumPars[RDETPAGE]= 8;
	this->NumPars[PDETPAGE]= 8;
	//GN
	this->NumPars[ESPECPAGE]= 10;
	//GN_!
	this->NumPars[RSPECPAGE]= 10;
	this->NumPars[EMOMPAGE]= 7;
	this->NumPars[RMOMPAGE]= 7+5; //+5 for interpolation spect
	this->NumPars[PPCOPAGE]= 12;
	this->NumPars[USERPAGE]= 0;


	this->PPCO_enabled = false;
	this->PPCO_momcon = false;
//	this->interpolation = false;
	
	this->LUT_enabled = false;

	this->wiggletool_enabled = false;

	this->fishtool_enabled = false;

	this->iontoftool_enabled = false;

	this->ipa_evt = new CH_event_struct();
	for(int i = 0; i<16;i++) {  //for polyatomic
		this->ipa_evt->e.method.push_back(0.0);
		this->ipa_evt->e.x.push_back(0.0);
		this->ipa_evt->e.y.push_back(0.0);
		this->ipa_evt->e.tof.push_back(0.0);
		this->ipa_evt->e.time.push_back(0.0);
		this->ipa_evt->r.method.push_back(0.0);
		this->ipa_evt->r.x.push_back(0.0);
		this->ipa_evt->r.y.push_back(0.0);
		this->ipa_evt->r.tof.push_back(0.0);
		this->ipa_evt->r.time.push_back(0.0);
		this->ipa_evt->p.method.push_back(0.0);
		this->ipa_evt->p.x.push_back(0.0);
		this->ipa_evt->p.y.push_back(0.0);
		this->ipa_evt->p.tof.push_back(0.0);
		this->ipa_evt->p.time.push_back(0.0);
	}

	this->eventdata = new CH_event_struct[number_of_events];
	
	this->par = new ipa_parameter[256];
	this->num_par = 0;
	this->keys = "1q2w3e4r5t7u8i9o0psxdcfvgbhnjmzh";

	this->root = rt;

	for(int i = 0; i<16;i++) {
		for(int j=0;j<MAXTHREADS;j++) {
			this->e[j][i] = new electron_class();
			this->r[j][i] = new ion_class();
		}
		this->e_raw[i] = new particle_raw();
		this->r_raw[i] = new particle_raw();
		this->p_raw[i] = new particle_raw();
	}

	for(int j=0;j<MAXTHREADS;j++) {
		this->mol[j] = new diatomic_class();
		this->big_mol[j] = new polyatomic_class();
		this->IPA_Proc[j] = new processor_class();
		this->IPA_Proc[j]->Spect = new spectrometer_class();
	}

	this->cur_reaction = new reaction_struct();
	this->reaction_def_valid = false;
			
	this->cut_wiggle = false;
	this->cut_wiggle_width = 0.0;

	for(int i=0;i<16;i++) {
		for(int j=0;j<MAXTHREADS;j++) {
			H1D[j][i] = 0;
			H2D[j][i] = 0;
		}
	}
	
}

void ipa_class::Setup(ColAHelL* MasterCH) {

	this->set_fac_rot(0, MasterCH->GetCorParam(0));
	this->set_mom_fac(0, MasterCH->GetMomCorParam(0));
	this->set_fac_rot(1, MasterCH->GetCorParam(1));
	this->set_mom_fac(1, MasterCH->GetMomCorParam(1), std::vector<double>(5, 0));
/*	Sorry Josh, we need to fix this... :O(
	if (Ueber->Proc->Spect->ion_side->interpolation_approximation) {
		this->set_mom_fac(1, Ueber->Proc->r_mom_fac, Ueber->Proc->reaction_list[(int)this->_channel]->interpolation_adjustment_variables[0]);
	}
	else {
		this->set_mom_fac(1, Ueber->Proc->r_mom_fac, std::vector<double>(5, 0));
	}
*/
	this->set_fac_rot(2, MasterCH->GetCorParam(2));
	this->set_mom_fac(2, MasterCH->GetMomCorParam(2));

	this->set_spec(MasterCH->Spect);
	this->set_tof_calc_class(MasterCH->CTof);

	if(MasterCH->IsReactionDefined(this->channel)) {
		reaction_struct* Reac = MasterCH->Proc->reaction_list[MasterCH->Proc->find_reaction(this->channel)];
		this->set_reaction(Reac);
		
		this->set_fac_rot(0, Reac->e_fac);
		this->set_mom_fac(0, Reac->e_mom_fac);
		this->set_fac_rot(1, Reac->r_fac);
		this->set_mom_fac(1, Reac->r_mom_fac);
/*
		if (Ueber->Proc->Spect->ion_side->interpolation_approximation) {
				this->set_mom_fac(1, Ueber->Proc->reaction_list[i]->r_mom_fac, Ueber->Proc->reaction_list[(int)this->_channel]->interpolation_adjustment_variables[0]);
			}
			else {
				this->set_mom_fac(1, Ueber->Proc->reaction_list[i]->r_mom_fac, std::vector<double>(5, 0));
			}
*/		
		this->set_fac_rot(2, Reac->p_fac);
		this->set_mom_fac(2, Reac->p_mom_fac);
	} else {
		this->cur_reaction->e_fac = new cor_param();
		this->cur_reaction->e_mom_fac = new mom_cor_param();
		this->cur_reaction->r_fac = new cor_param();
		this->cur_reaction->r_mom_fac = new mom_cor_param();
		this->cur_reaction->p_fac = new cor_param();
		this->cur_reaction->p_mom_fac = new mom_cor_param();
	}
	WriteParToStruct();
}

int ipa_class::GetChannel() {
	return this->channel;
}

bool ipa_class::is_PPCO_tool_enabled(){
	return this->PPCO_enabled;
}

void ipa_class::set_fac_rot(int det, cor_param *fac) 
{
	if(det==0) {
		add_par(fac->dx,-80.0,80.0,0.1,"Shift pos x by",0);
		add_par(fac->dy,-180.0,180.0,0.1,"Shift pos y by",1);
		add_par(fac->dt,-1800.0,1800.0,0.1,"Shift t by",2);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch pos x by",3);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch pos y by",4);
		add_par(fac->overall_stretch,0.2,10.0,0.01,"Stretch pos by",5);
		add_par(fac->rot_ang,-180.0,180.0,1.0,"Rotate pos by",6);
		if(fac->raw_order)
			add_par(1.0,0.0,1.0,1.0,"Shift/Str then Rot",7);
		else
			add_par(0.0,0.0,1.0,1.0,"Shift/Str then Rot",7);
	}
	if(det==1) {
		add_par(fac->dx,-80.0,80.0,0.1,"Shift pos x by",10);
		add_par(fac->dy,-80.0,80.0,0.1,"Shift pos y by",11);
		add_par(fac->dt,-1800.0,1800.0,1.0,"Shift t by",12);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch pos x by",13);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch pos y by",14);
		add_par(fac->overall_stretch,0.2,10.0,0.01,"Stretch pos by",15);
		add_par(fac->rot_ang,-180.0,180.0,1.0,"Rotate pos by",16);
		add_par(fac->raw_order,0.0,1.0,1.0,"Shift/Str then Rot",17);
	}

	if(det==2) {
		add_par(fac->dx,-80.0,80.0,0.1,"Shift pos x by",20);
		add_par(fac->dy,-80.0,80.0,0.1,"Shift pos y by",21);
		add_par(fac->dt,-1800.0,1800.0,0.1,"Shift t by",22);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch pos x by",23);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch pos y by",24);
		add_par(fac->overall_stretch,0.2,10.0,0.01,"Stretch pos by",25);
		add_par(fac->rot_ang,-180.0,180.0,1.0,"Rotate pos by",26);
		add_par(fac->raw_order,0.0,1.0,1.0,"Shift/Str then Rot",27);
	}
}

void ipa_class::set_mom_fac(int det, mom_cor_param *fac,  std::vector<double>  interpolation_adjustment_ceff )
{
	if(det==0) {
		add_par(fac->dx,-80.0,80.0,0.1,"Shift x by",50);
		add_par(fac->dy,-180.0,180.0,0.1,"Shift y by",51);
		add_par(fac->dz,-1800.0,1800.0,0.1,"Shift z by",52);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch x by",53);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch y by",54);
		add_par(fac->z_stretch,0.2,2.0,0.01,"Stretch z by",55);
		add_par(fac->overall_stretch,0.2,10.0,0.01,"Stretch ALL by",56);
	}
	if(det==1) {
		add_par(fac->dx,-80.0,80.0,0.1,"Shift x by",60);
		add_par(fac->dy,-180.0,180.0,0.1,"Shift y by",61);
		add_par(fac->dz,-1800.0,1800.0,0.1,"Shift z by",62);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch x by",63);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch y by",64);
		add_par(fac->z_stretch,0.2,2.0,0.01,"Stretch z by",65);
		add_par(fac->overall_stretch,0.2,10.0,0.01,"Stretch ALL by",66);
		//
		add_par(interpolation_adjustment_ceff[0], -1., 1.0, 0.001, "interpol. ceff x      ", 67);
		add_par(interpolation_adjustment_ceff[1], -1., 1.0, 0.001, "interpol. ceff x^2    ", 68);
		add_par(interpolation_adjustment_ceff[2], -1., 1.0, 0.001, "interpol. ceff y      ", 69);
		add_par(interpolation_adjustment_ceff[3], -1., 1.0, 0.001, "interpol. ceff y^2    ", 70);
		add_par(interpolation_adjustment_ceff[4], -1., 1.0, 0.001, "interpol. ceff r      ", 71);
	}

	if(det==2) {
	/*	add_par(fac->dx,-80.0,80.0,0.1,"Shift pos x by",20);
		add_par(fac->dy,-80.0,80.0,0.1,"Shift pos y by",21);
		add_par(fac->dz,-1800.0,1800.0,0.1,"Shift t by",22);
		add_par(fac->x_stretch,0.2,2.0,0.01,"Stretch pos x by",23);
		add_par(fac->y_stretch,0.2,2.0,0.01,"Stretch pos y by",24);
		add_par(fac->overall_stretch,0.2,2.0,0.01,"Stretch pos by",25);*/
	}


}

void ipa_class::set_spec(spectrometer_class *Spect) {
	add_par(Spect->electron_side->number_of_regions,0,3.0,1.0,"Number of regions",30);
	add_par(Spect->electron_side->Efields[0],-1000,1000.0,0.01,"E(0) [V/cm]",31);
	add_par(Spect->electron_side->Efields[1],-1000,1000.0,0.01,"E(1) [V/cm]",32);
	add_par(Spect->electron_side->Efields[2],-1000,1000.0,0.01,"E(2) [V/cm]",33);
	add_par(Spect->electron_side->lengths[0],0,1000.0,0.1,"Length(0) [mm]",34);
	add_par(Spect->electron_side->lengths[1],0,1000.0,0.1,"Length(1) [mm]",35);
	add_par(Spect->electron_side->lengths[2],0,1000.0,0.1,"Length(2) [mm]",36);
	add_par(Spect->Bfield_ns,0,1000.0,0.1,"B-Field [ns]",37);
	add_par(Spect->Bfield_clockwise,0,1.0,1.,"B-Field is clk-wise",38);
	//GN
	add_par(Spect->MeanTOFe,-1000,1000.,0.01,"LIN: MeanTOF [ns]",39);
	//GN_!

	add_par(Spect->ion_side->number_of_regions,0,3.0,1.0,"Number of regions",40);
	add_par(Spect->ion_side->linear_approximation,0,1.0,1.0,"Linear approximation",41);
	add_par(Spect->ion_side->Efields[0],-1000,1000.0,0.01,"E(0) [V/cm]",42);
	add_par(Spect->ion_side->Efields[1],-1000,1000.0,0.01,"E(1) [V/cm]",43);
	add_par(Spect->ion_side->Efields[2],-1000,1000.0,0.01,"E(2) [V/cm]",44);
	add_par(Spect->ion_side->lengths[0],0,1000.0,0.1,"Length(0) [mm]",45);
	add_par(Spect->ion_side->lengths[1],0,1000.0,0.1,"Length(1) [mm]",46);
	add_par(Spect->ion_side->lengths[2],0,1000.0,0.1,"Length(2) [mm]",47);
	add_par(Spect->VJet,0.0,5000.0,10.0,"Jet velocity [m/s]",48);
	add_par(Spect->AngJet,-180.0,360.0,5.0,"Jet direction [deg]",49);
	
	//if (Spect->ion_side->interpolation_approximation == 1) {
	//	this->interpolation = true;
	//}
	
}

void ipa_class::set_tof_calc_class(tof_calc_class *tof) {
	this->IPA_Tof = tof;
}

void ipa_class::set_reaction(reaction_struct *reac){
	this->reaction_def_valid = true;
	this->cur_reaction = reac;
}

void ipa_class::append_data(CH_event_struct *evnt)
{
	if(this->wiggletool_enabled) { 
		Fill_wiggletool(evnt);
		return;
	}
	if(this->fishtool_enabled) { 
		Fill_fishtool(evnt);
		return;
	}
	if(this->iontoftool_enabled) { 
		Fill_iontoftool(evnt);
		return;
	}
	if(this->PPCO_enabled) {
		FillPPCO(evnt);
		return;
	}
	if(this->num_events_read < this->num_events) {
		if(this->channel == -1 || this->channel == evnt->reaction) {
			ipa_evt->bunchmarker = evnt->bunchmarker;
			ipa_evt->reaction = evnt->reaction;
			
			ipa_evt->bunchmarker = evnt->bunchmarker;
			ipa_evt->e.num_hits = (evnt->e.num_hits<17) ? evnt->e.num_hits : 16; //polyatomic, before: 3
			ipa_evt->r.num_hits = (evnt->r.num_hits<17) ? evnt->r.num_hits : 16;  
			ipa_evt->p.num_hits = (evnt->p.num_hits<17) ? evnt->p.num_hits : 16;
			for(int i=0;i<ipa_evt->e.num_hits;i++) {
				ipa_evt->e.method[i] = evnt->e.method[i];
				ipa_evt->e.x[i] = evnt->e.x[i];
				ipa_evt->e.y[i] = evnt->e.y[i];
				ipa_evt->e.tof[i] = evnt->e.tof[i];
				ipa_evt->e.time[i] = evnt->e.time[i];
			}
			for(int i=0;i<ipa_evt->r.num_hits;i++) {
				ipa_evt->r.method[i] = evnt->r.method[i];
				ipa_evt->r.x[i] = evnt->r.x[i];
				ipa_evt->r.y[i] = evnt->r.y[i];
				ipa_evt->r.tof[i] = evnt->r.tof[i];
				ipa_evt->r.time[i] = evnt->r.time[i];
			}
			for(int i=0;i<ipa_evt->p.num_hits;i++) {
				ipa_evt->p.method[i] = evnt->p.method[i];
				ipa_evt->p.x[i] = evnt->p.x[i];
				ipa_evt->p.y[i] = evnt->p.y[i];
				ipa_evt->p.tof[i] = evnt->p.tof[i];
				ipa_evt->p.time[i] = evnt->p.time[i];
			}
			
			std::copy(ipa_evt, ipa_evt + 1, &eventdata[this->num_events_read]);
			num_events_read++;
		}
	}
}

void ipa_class::add_par(double val, double minimum, double maximum, double steps, char *name, int n)
{
	if(n==-1) {
		n = USERPAGE*10 + this->num_par++;
		this->NumPars[USERPAGE]++;
	}
	this->par[n].value = val;
	this->par[n].min = minimum;
	this->par[n].max = maximum;
	this->par[n].step = steps;
	this->par[n].name = name;
}

void ipa_class::main_loop()
{
	draw_screen();
	init_interactive();
	
	if(this->PPCO_enabled) {
		init_PIPICO();
		WriteParToStruct();
		draw_PIPICO();
	}

	if(this->LUT_enabled) {
		init_extract_LUT();
		WriteParToStruct();
	}

	if(this->wiggletool_enabled) {
		init_wiggletool();
		WriteParToStruct();
	}

	if(this->fishtool_enabled) {
		init_fishtool();
		WriteParToStruct();
	}

	if(this->iontoftool_enabled) {
		init_iontoftool();
		WriteParToStruct();
	}

	while(true)
	{
		__int32 c = my_kbhit();
		
		if (c) {
			while (my_kbhit());
			switch (c) {
				case 0x3b:
					this->page = EDETPAGE;
					break;
				case 0x3c:
					this->page = RDETPAGE;
					break;
				case 0x3d:
					this->page = PDETPAGE;
					break;
				case 0x3e:
					this->page = ESPECPAGE;
					break;
				case 0x3f:
					this->page = RSPECPAGE;
					break;
				case 0x40:
					this->page = EMOMPAGE;
					break;
				case 0x41:
					this->page = RMOMPAGE;
					break;
				case 0x42:
					this->page = USERPAGE;
					break;
				case 0x43:
					if(this->PPCO_enabled || this->LUT_enabled || this->wiggletool_enabled || this->fishtool_enabled || this->iontoftool_enabled)
						this->page = PPCOPAGE;
					break;
				default:
					break;
			}
			
			for(int i=0;i<this->NumPars[this->page];i++) {
				if(c == this->keys[i*2]) {
					if(par[this->page*10+i].value<par[this->page*10+i].max)
						par[this->page*10+i].value += par[this->page*10+i].step*this->precision;
				}
				if(c == this->keys[i*2+1]) {
					if(par[this->page*10+i].value>par[this->page*10+i].min)
						par[this->page*10+i].value -= par[this->page*10+i].step*this->precision;
				}
			}
			
			if(c=='Q')
				break;
			
			if(c==' ') {
				WriteParToStruct();
				
				if(this->PPCO_enabled)
					draw_PIPICO();
				if(this->wiggletool_enabled) 
					draw_wiggletool();
				if(this->fishtool_enabled) 
					draw_fishtool();
				if(this->iontoftool_enabled) 
					draw_iontoftool();

				// only loop over data if we don't use the PPCO or wiggle tool or ... 
				if(!this->PPCO_enabled && !this->wiggletool_enabled && !this->fishtool_enabled && !this->iontoftool_enabled) {
					for(int i=0;i<16;i++) {	
						for(int j=0;j<MAXTHREADS;j++) {
							if(H1D[j][i]!=0)
								H1D[j][i]->Reset();
							if(H2D[j][i]!=0)
								H2D[j][i]->Reset();
						}
					}
					int th_num;
					int idx;
					#pragma omp parallel default(shared) private(idx, th_num)
					{
						bool valid_ions[MAXTHREADS];
						
						//loop here
						#pragma omp for schedule(dynamic) /*default(shared) private(idx, th_num) */
						for(idx=0;idx<this->num_events_read;idx++) {

							// Do processing.. ColAHelL can only work on ions in case the reaction is properly defined...
							th_num = omp_get_thread_num();
							valid_ions[th_num] = false;
							if(reaction_def_valid) {
								for(int i=0;i<eventdata[idx].e.num_hits;i++)
									e[th_num][i]->valid = false;

								valid_ions[th_num] = IPA_Proc[th_num]->process_ions(cur_reaction, &eventdata[idx], r[th_num], mol[th_num], big_mol[th_num]);
							}
							if(valid_ions || !reaction_def_valid) 
								bool dummy = IPA_Proc[th_num]->process_electrons(cur_reaction, &eventdata[idx], e[th_num]);
							
							if(this->LUT_enabled)
								fill_extract_LUT(&eventdata[idx],th_num);
								
							user_function(&eventdata[idx],th_num);
						}
					}//end of paralell region

					combine_histos();
					user_draw();
				}
			}
			
			if(c=='+') 
				precision*=10.0;
			
			if(c=='-') 
				precision/=10.0;
			
			// Lookuptable tool special keys
			if(this->LUT_enabled) {
				if(c=='C')
					write_correctiontable();

				if(c=='V') {
					WriteParToStruct();
					draw_extract_LUT();			
				}
			}

			draw_screen();
			
		}
		gSystem->Sleep(20);
		gSystem->ProcessEvents();
	}
	delete this->ipa_evt;
//	delete[] this->eventdata;
}

int ipa_class::events_in_memory() {
	return this->num_events_read;
}

void ipa_class::WriteParToStruct()
{
	
// correction factors
	for(int j=0;j<MAXTHREADS;j++) {
		IPA_Proc[j]->e_fac->dx = par[0].value;
		IPA_Proc[j]->e_fac->dy = par[1].value;
		IPA_Proc[j]->e_fac->dt = par[2].value;
		IPA_Proc[j]->e_fac->x_stretch = par[3].value;
		IPA_Proc[j]->e_fac->y_stretch = par[4].value;
		IPA_Proc[j]->e_fac->overall_stretch = par[5].value;
		IPA_Proc[j]->e_fac->rot_ang = par[6].value;
		if(par[7].value>0.5)
			IPA_Proc[j]->e_fac->raw_order = true;
		else
			IPA_Proc[j]->e_fac->raw_order = false;

		cur_reaction->e_fac->dx = par[0].value;
		cur_reaction->e_fac->dy = par[1].value;
		cur_reaction->e_fac->dt = par[2].value;
		cur_reaction->e_fac->x_stretch = par[3].value;
		cur_reaction->e_fac->y_stretch = par[4].value;
		cur_reaction->e_fac->overall_stretch = par[5].value;
		cur_reaction->e_fac->rot_ang = par[6].value;
		cur_reaction->e_fac->raw_order = IPA_Proc[j]->e_fac->raw_order;

		IPA_Proc[j]->r_fac->dx = par[10].value;
		IPA_Proc[j]->r_fac->dy = par[11].value;
		IPA_Proc[j]->r_fac->dt = par[12].value;
		IPA_Proc[j]->r_fac->x_stretch = par[13].value;
		IPA_Proc[j]->r_fac->y_stretch = par[14].value;
		IPA_Proc[j]->r_fac->overall_stretch = par[15].value;
		IPA_Proc[j]->r_fac->rot_ang = par[16].value;
		if(par[17].value>0.5)
			IPA_Proc[j]->r_fac->raw_order = true;
		else
			IPA_Proc[j]->r_fac->raw_order = false;

		if(reaction_def_valid) {
			cur_reaction->r_fac->dx = par[10].value;
			cur_reaction->r_fac->dy = par[11].value;
			cur_reaction->r_fac->dt = par[12].value;
			cur_reaction->r_fac->x_stretch = par[13].value;
			cur_reaction->r_fac->y_stretch = par[14].value;
			cur_reaction->r_fac->overall_stretch = par[15].value;
			cur_reaction->r_fac->rot_ang = par[16].value;
			cur_reaction->r_fac->raw_order = IPA_Proc[j]->r_fac->raw_order;
		}

		IPA_Proc[j]->p_fac->dx = par[20].value;
		IPA_Proc[j]->p_fac->dy = par[21].value;
		IPA_Proc[j]->p_fac->dt = par[22].value;
		IPA_Proc[j]->p_fac->x_stretch = par[23].value;
		IPA_Proc[j]->p_fac->y_stretch = par[24].value;
		IPA_Proc[j]->p_fac->overall_stretch = par[25].value;
		IPA_Proc[j]->p_fac->rot_ang = par[26].value;
		if(par[27].value>0.5)
			IPA_Proc[j]->p_fac->raw_order = true;
		else
			IPA_Proc[j]->p_fac->raw_order = false;

		cur_reaction->p_fac->dx = par[20].value;
		cur_reaction->p_fac->dy = par[21].value;
		cur_reaction->p_fac->dt = par[22].value;
		cur_reaction->p_fac->x_stretch = par[23].value;
		cur_reaction->p_fac->y_stretch = par[24].value;
		cur_reaction->p_fac->overall_stretch = par[25].value;
		cur_reaction->p_fac->rot_ang = par[26].value;
		cur_reaction->p_fac->raw_order = IPA_Proc[j]->p_fac->raw_order;

	// spectrometer	
	
		E_e[0] = par[31].value;
		E_e[1] = par[32].value;
		E_e[2] = par[33].value;
		len_e[0] = par[34].value;
		len_e[1] = par[35].value;
		len_e[2] = par[36].value;
		//GN
		MeanTOFe[0] = par[39].value;
		//GN_!
		E_r[0] = par[42].value;
		E_r[1] = par[43].value;
		E_r[2] = par[44].value;
		len_r[0] = par[45].value;
		len_r[1] = par[46].value;
		len_r[2] = par[47].value;

		IPA_Proc[j]->Spect->set_VJet(par[48].value, par[49].value);
		//GN
		if (par[34].value == 0 && par[35].value == 0 && par[36].value == 0) {
			IPA_Proc[j]->Spect->set_electron_arm_lin(E_e[0]);
			IPA_Proc[j]->Spect->set_mean_tof_e(MeanTOFe[0]);
		}
		else
			IPA_Proc[j]->Spect->set_electron_arm((unsigned short)(par[30].value + 0.1),len_e,E_e);
		//GN_!
		IPA_Proc[j]->Spect->set_Bfield(par[37].value,(bool)par[38].value,true);
		IPA_Proc[j]->Spect->set_ion_arm((unsigned short)(par[40].value + 0.1),len_r,E_r);
		if(par[41].value > 0.5)
			IPA_Proc[j]->Spect->set_ion_arm_lin(E_r[0]);
//		if(this->interpolation)
//			IPA_Proc[j]->Spect->set_ion_arm_interpolation();

	// mom correction factors
		IPA_Proc[j]->e_mom_fac->dx = par[50].value;
		IPA_Proc[j]->e_mom_fac->dy = par[51].value;
		IPA_Proc[j]->e_mom_fac->dz = par[52].value;
		IPA_Proc[j]->e_mom_fac->x_stretch = par[53].value;
		IPA_Proc[j]->e_mom_fac->y_stretch = par[54].value;
		IPA_Proc[j]->e_mom_fac->z_stretch = par[55].value;
		IPA_Proc[j]->e_mom_fac->overall_stretch = par[56].value;
		if(reaction_def_valid){
			cur_reaction->e_mom_fac->dx = par[50].value;
			cur_reaction->e_mom_fac->dy = par[51].value;
			cur_reaction->e_mom_fac->dz = par[52].value;
			cur_reaction->e_mom_fac->x_stretch = par[53].value;
			cur_reaction->e_mom_fac->y_stretch = par[54].value;
			cur_reaction->e_mom_fac->z_stretch = par[55].value;
			cur_reaction->e_mom_fac->overall_stretch = par[56].value;
		}

		IPA_Proc[j]->r_mom_fac->dx = par[60].value;
		IPA_Proc[j]->r_mom_fac->dy = par[61].value;
		IPA_Proc[j]->r_mom_fac->dz = par[62].value;
		IPA_Proc[j]->r_mom_fac->x_stretch = par[63].value;
		IPA_Proc[j]->r_mom_fac->y_stretch = par[64].value;
		IPA_Proc[j]->r_mom_fac->z_stretch = par[65].value;
		IPA_Proc[j]->r_mom_fac->overall_stretch = par[66].value;
		if(reaction_def_valid) {
			cur_reaction->r_mom_fac->dx = par[60].value;
			cur_reaction->r_mom_fac->dy = par[61].value;
			cur_reaction->r_mom_fac->dz = par[62].value;
			cur_reaction->r_mom_fac->x_stretch = par[63].value;
			cur_reaction->r_mom_fac->y_stretch = par[64].value;
			cur_reaction->r_mom_fac->z_stretch = par[65].value;
			cur_reaction->r_mom_fac->overall_stretch = par[66].value;
		}
	}

	if(this->PPCO_enabled) {
		this->ppco_m1[0] = par[70].value;
		this->ppco_m2[0] = par[71].value;
		this->ppco_q1[0] = par[72].value;
		this->ppco_q2[0] = par[73].value;

		this->ppco_m1[1] = par[74].value;
		this->ppco_m2[1] = par[75].value;
		this->ppco_q1[1] = par[76].value;
		this->ppco_q2[1] = par[77].value;

		this->ppco_m1[2] = par[78].value;
		this->ppco_m2[2] = par[79].value;
		this->ppco_q1[2] = par[80].value;
		this->ppco_q2[2] = par[81].value;

	/*	this->ppco_m1[3] = par[82].value;
		this->ppco_m2[3] = par[83].value;
		this->ppco_q1[3] = par[84].value;
		this->ppco_q2[3] = par[85].value;*/
	}
	if(this->LUT_enabled) {
		this->LUT_zero = (int)par[70].value;
		this->LUT_mergeNbins = (int)par[71].value;
		this->LUT_slidingmerge = (int)par[72].value;
		this->LUT_momrange_min = par[73].value;
		this->LUT_momrange_max = par[74].value;
		this->LUT_mom = par[75].value;	
		this->LUT_Profile =(int)par[76].value;
	}
}

void ipa_class::draw_screen()
{
	cls();
	
	int from = 10*this->page;
	int to = 10*this->page + this->NumPars[this->page];
	
//	if ( !(this->interpolation) && this->page== RMOMPAGE)//display extra paramenters if interpolation is used
//		to -= 5;

	Green(true);
	printf("\nInteractive Parameter Adjustment (IPA) Version %.2f\n",IPAVERSION);
	printf("========================================================================\n");
	White(false);
	printf("Events: %d",this->num_events_read);
	printf(", Thrds: %d", this->cores_in_use);
	printf(", prec: %.2f", this->precision);

	White(true);
	switch (this->page) {
	case EDETPAGE:
		printf("		Electron detector parameters\n\n");
		break;
	case RDETPAGE:
		printf("		Recoil detector parameters\n\n");
		break;
	case PDETPAGE:
		printf("		Projectile detector parameters\n\n");
		break;
	case ESPECPAGE:
		printf("		Spectrometer parameters electron arm\n\n");
		break;
	case RSPECPAGE:
		printf("		Spectrometer parameters ion arm\n\n");
		break;
	case EMOMPAGE:
		printf("		Electron momentum parameters\n\n");
		break;
	case RMOMPAGE:
		printf("		Recoil momentum parameters\n\n");
		break;
	case PPCOPAGE:
		if(this->PPCO_enabled)
			printf("		M and Q for PIPICO drawing tool\n\n");
		if(this->LUT_enabled)
			printf("		Lookuptable tool parameters\n\n");
		if(this->wiggletool_enabled)
			printf("		Wiggle tool parameters\n\n");
		if(this->fishtool_enabled)
			printf("		Fish tool parameters\n\n");
		if(this->iontoftool_enabled)
			printf("		Ion tof tool parameters\n\n");
		break;
	case USERPAGE:
		printf("		User defined parameters\n\n");
		to = 10*this->page + this->num_par;
		break;
	default:
		break;
	}
	White(false);


	printf("Adjustable parameters are:\n");
	printf("-------------------------\n");

	int k=0;
	for(int i=from;i<to;i++)
	{
		char inc[2] = " ";
		inc[0] = keys[2*(k)];
		char dec[2] = " "; 
		dec[0] = keys[2*(k++)+1];

		printf("Name: ");
		White(true);
		printf("%20.20s: ",par[i].name);
		printf("%.4f ",par[i].value);
		White(false);
		printf("(press '%.1s' to inc , '%.1s' to dec)\n",&inc,&dec);
	}

	White(true);
	printf("\n==================== Press SPACE to update canvas! ===================== \n");
	Green(true);
	printf("Press F1..F3 for electron, ion, or projectile detector page.\n");
	printf("Press F4,F5 for spectrometer definitions.\n");
	printf("Press F6,F7 for momentum tweaking.\n");
	printf("...or F8 for user parameters.");
	if(this->PPCO_enabled)
		printf("...or, finally, F9 for PIPICO tool.\n");
	if(this->LUT_enabled)
		printf("...or, finally, F9 for LUT tool.\n");
	if(this->wiggletool_enabled)
		printf("...or, finally, F9 for wiggle tool.\n");
	if(this->fishtool_enabled)
		printf("...or, finally, F9 for fish tool.\n");
	if(this->iontoftool_enabled)
		printf("...or, finally, F9 for ion tof tool.\n");	
	Red(true);
	printf("..(or press Q (i.e. 'SHIFT + q') to quit...)\n");
	White(false);

}

void ipa_class::std_calc(CH_event_struct *evti) {

	// Do processing.. ColAHelL can only work on ions in case the reaction is properly defined...
/*	bool valid_ions = false;
	if(reaction_def_valid) {
		for(int i=0;i<evti->e.num_hits;i++)
			e[i]->valid = false;

		//this still only works for one ion!!!
		//maybe that is good enough
/*		cur_reaction->interpolation_adjustment_variables[0][0]=this->par[67].value;
		cur_reaction->interpolation_adjustment_variables[0][1]=this->par[68].value;
		cur_reaction->interpolation_adjustment_variables[0][2]=this->par[69].value;
		cur_reaction->interpolation_adjustment_variables[0][3]=this->par[70].value;
		cur_reaction->interpolation_adjustment_variables[0][4]=this->par[71].value;
*/
/*		valid_ions = IPA_Proc->process_ions( cur_reaction, evti, r, mol, big_mol);
	}
	if(valid_ions || !reaction_def_valid) 
		bool dummy = IPA_Proc->process_electrons( cur_reaction, evti, e);
*/				
}

void ipa_class::init_interactive() {
	
	gStyle->SetOptStat("eou");
	gStyle->SetNdivisions(508,"XY");
	gStyle->SetOptFit(0000);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetCanvasColor(0);

/*	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetTitleColor(0);
	gStyle->SetStatColor(0);
	
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameFillStyle(1001);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameBorderMode(0);
	*/
	gStyle->SetStatColor(0);
	gStyle->SetStatTextColor(1);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatFont(62);
	gStyle->SetStatFontSize(0.03f);
	gStyle->SetStatStyle(0);
	gStyle->SetStatFormat("6.4g");
	gStyle->SetStatX(0.9f);
	gStyle->SetStatY(0.9f);
	gStyle->SetStatW(0.5f);
	gStyle->SetStatH(0.3f);

	gStyle->SetPadGridX(kTRUE);
	gStyle->SetPadGridY(kTRUE);

	root->color();

	if(!PPCO_enabled && !wiggletool_enabled && !fishtool_enabled && !iontoftool_enabled)
		init_user_histograms();

}

void ipa_class::combine_histos() {
	for(int i=0;i<16;i++) {
		if(H1D[0][i]) {
			for(int th=1;th<MAXTHREADS;th++) {
				H1D[0][i]->Add(H1D[th][i]);
			}
		}
		if(H2D[0][i]) {
			for(int th=1;th<MAXTHREADS;th++) {
				H2D[0][i]->Add(H2D[th][i]);
			}
		}
	}
	//merge histos from multicore support..
	if(LUT_enabled) {
		for(int i=0;i<this->LUT_bins_phi;i++) {
			for(int th=1;th<MAXTHREADS;th++) {
				ctpeArray[0][i]->Add(ctpeArray[th][i]);
			}
		}
	}
}
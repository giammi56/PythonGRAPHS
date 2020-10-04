#include "CH_Histograms.h"

// Standard histograms for reactions of "ION" and "ION_ELEC" type
namespace CH
{
	
	// there is nothing to plot in case we did not have a reaction definition
	void histograms_class::plot_ions(reaction_struct *cur_reaction, ion_class ** r, int rhit, electron_class ** e, int ehit, double *scan_val) 
	{
		// limit to 16 hits.. 
		if(rhit>16)
			rhit=16;
		
		histo_handler *Hist = this->Hist_ions;
		// histogram index stuff
		int hmax = HISTS_PER_CHANNEL; // overall number of ColAHelL histograms per reaction/channel
		int hoff = hmax*16*(int)cur_reaction->reac_num;

		char reaction_dir[360];
		//char x_axis_title[80];
		//char y_axis_title[80];

		char rootdir[360];
		char ndir[360];
		char tdir[360];

		if(this->use_master_folder) {
			strcpy(reaction_dir,this->CH_master_folder);
			strcat(reaction_dir,cur_reaction->name);
		} else
			strcpy(reaction_dir,cur_reaction->name);			

		strcpy(rootdir,reaction_dir);
		strcat(rootdir,"/ions");

		for(int i=0;i<(int)rhit;i++) {
			strcpy(ndir,rootdir);
			sprintf(tdir,"/hit_%i",i);
			strcat(ndir,tdir);
			strcpy(tdir,ndir);
						
			if(r[i]->valid) {												

				Hist->fill(hoff+hmax*i+0,"ion energy",r[i]->energy(),1.,"ion energy",50,0.0,0.250,"ion energy [eV]",strcat(ndir,"/energy"));
				Hist->fill(hoff+hmax*i+1,"phi vs ion energy",r[i]->raw.phi,r[i]->energy(),1.,"ion energy vs. phi on detector",72,-180.0,180.0,"phi [deg]",50,0.0,25.0,"ion energy [eV]",ndir); 
				Hist->fill(hoff+hmax*i+2,"ctheta vs ion energy",r[i]->mom.Cos_Theta(),r[i]->energy(),1.,"ion energy vs. cos(theta_z)",72,-1.0,1.0,"ctheta",50,0.0,25.0,"ion energy [eV]",ndir); 
						
				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+3,"p_x vs p_y",r[i]->mom.x,r[i]->mom.y,1.,"p_x vs. p_y",50,-10.5,10.5,"p_x [a.u.]",50,-10.5,10.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
				Hist->fill(hoff+hmax*i+4,"p_x vs p_z",r[i]->mom.x,r[i]->mom.z,1.,"p_x vs. p_z",50,-10.5,10.5,"p_x [a.u.]",50,-10.5,10.5,"p_z [a.u.]",ndir);
				Hist->fill(hoff+hmax*i+5,"p_y vs p_z",r[i]->mom.y,r[i]->mom.z,1.,"p_y vs. p_z",50,-10.5,10.5,"p_y [a.u.]",50,-10.5,10.5,"p_z [a.u.]",ndir);
				Hist->fill(hoff+hmax*i+6,"p_mag",r[i]->mom.Mag(),1.,"|p|",50,0.0,0.25,"|p| [a.u.]",ndir);

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+7,"phi",r[i]->mom.Phi_deg(),1.,"phi in the labframe",72,-180.0,180.0,"phi [deg]",strcat(ndir,"/angles_labframe"));
				Hist->fill(hoff+hmax*i+10,"ctheta",r[i]->mom.Cos_Theta(),1.,"cos(theta) in the labframe",72,-1.0,1.0,"cos(theta)",ndir);

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+13,"position",r[i]->raw.data.x,r[i]->raw.data.y,1.,"position",400,-1.*rdet_size,rdet_size,"x [mm]",400,-1.*rdet_size,rdet_size,"y [mm]",strcat(ndir,"/raw"));
				Hist->fill(hoff+hmax*i+14,"tof",r[i]->raw.data.tof,1.,"time-of-flight",4000,-100.,10000.,"tof [ns]",ndir);
				Hist->fill(hoff+hmax*i+15,"wiggle",r[i]->raw.data.tof,(sqrt(r[i]->raw.data.x*r[i]->raw.data.x+r[i]->raw.data.y*r[i]->raw.data.y)),1.,"wiggles",4000,-100.,10000.,"tof [ns]",200,-1.,rdet_size,"r [mm]",ndir);
				Hist->fill(hoff+hmax*i+16,"fish x",r[i]->raw.data.tof,r[i]->raw.data.x,1.,"x-fish",4000,-100.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"x [mm]",ndir);
				Hist->fill(hoff+hmax*i+17,"fish y",r[i]->raw.data.tof,r[i]->raw.data.y,1.,"y-fish",4000,-100.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"y [mm]",ndir);
				if ( fabs(r[i]->raw.data.y)<10. ) 
					Hist->fill(hoff+hmax*i+18,"fish filet x",r[i]->raw.data.tof,r[i]->raw.data.x,1.,"x-fish filet",2000,-100.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"x [mm]",ndir);
				if ( fabs(r[i]->raw.data.x)<10. ) 
					Hist->fill(hoff+hmax*i+19,"fish filet y",r[i]->raw.data.tof,r[i]->raw.data.y,1.,"y-fish filet",2000,-100.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"y [mm]",ndir);

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+20,"pr_x vs pe_x",r[i]->mom.x,e[i]->mom.x,1.,"p_x vs. p_x",50,-2.5,2.5,"pr_x [a.u.]",50,-2.5,2.5,"pe_x [a.u.]",strcat(ndir,"/ion coincidence"));
				Hist->fill(hoff+hmax*i+21,"pr_y vs pe_y",r[i]->mom.y,e[i]->mom.y,1.,"p_y vs. p_y",50,-2.5,2.5,"pr_y [a.u.]",50,-2.5,2.5,"pe_y [a.u.]",ndir);
				Hist->fill(hoff+hmax*i+22,"pr_z vs pe_z",r[i]->mom.z,e[i]->mom.z,1.,"p_z vs. p_z",50,-2.5,2.5,"pr_z [a.u.]",50,-2.5,2.5,"pe_z [a.u.]",ndir);

			}
		}
	}
}
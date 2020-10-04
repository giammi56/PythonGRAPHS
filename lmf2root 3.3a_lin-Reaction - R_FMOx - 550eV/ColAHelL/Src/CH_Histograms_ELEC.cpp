#include "CH_Histograms.h"

namespace CH
{
	void histograms_class::plot_electrons(reaction_struct *cur_reaction, electron_class ** e, int ehit, double *scan_val) 
	{
		histo_handler *Hist = this->Hist_e;

		// limit to 16 hits.. 
		if(ehit>15)
			ehit=15;

		// histogram index stuff
		int hmax = HISTS_PER_CHANNEL; // overall number of ColAHelL histograms per reaction/channel
		int hoff = hmax*15*(int)cur_reaction->reac_num;
		int ID = cur_reaction->ID;

		char reaction_dir[360];
		//		char x_axis_title[80];
		//		char y_axis_title[80];
		char rootdir[360];
		char ndir[360];
		char tdir[360];

		//SET HERE THE VALUE!!!
		
		double emin_el = 9.;
		double emax_el = 14.;
		//

		if(this->use_master_folder) {
			strcpy(reaction_dir,this->CH_master_folder);
			strcat(reaction_dir,cur_reaction->name);
		} else
			strcpy(reaction_dir,cur_reaction->name);		

		strcpy(rootdir,reaction_dir);
		strcat(rootdir,"/electrons");

		for(int i=0;i<(int)ehit;i++) {

			strcpy(ndir,rootdir);
			sprintf(tdir,"/hit_%i",i);
			strcat(ndir,tdir);
			strcpy(tdir,ndir);

			if(e[i]->valid) {
				// fill standard histograms
				// Ar series
				if (cur_reaction->channel == 248) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 12.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 1000, 20.0, 6.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 4.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 4.0,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 252) {				
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 1000, 0.0, 12.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 8.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 8.0,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 254) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 12.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 10.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 10.0,"electron energy [eV]", ndir);
				} 
				else if (cur_reaction->channel == 256) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 14.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 12.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 12.0,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 260) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 14.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 14.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 14.0,"electron energy [eV]", ndir);
				}
				// Molecular brake-ups series
				else if (cur_reaction->channel == 9) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 4.0, 18.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 4.0, 18.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 4.0, 18.0,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 10) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 20.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 14.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 14.0,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 14.0,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 11) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 100.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 25.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy_small", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+11,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy_small", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+12,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 12) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 100.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 25.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy_small", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+11,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy_small", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+12,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 13) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 100.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 25.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy_small", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+11,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy_small", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+12,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
				}
				else if (cur_reaction->channel == 14) {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 800, 0.0, 100.0,"electron energy [eV]",strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 0.0, 25.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy_small", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+11,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy_small", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el,"electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+12,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, emin_el, emax_el*1.5,"electron energy [eV]", ndir);

				}
				else {
					Hist->fill(hoff+hmax*i+9,"electron energy", e[i]->energy(), 1., "electron energy", 360, 0.0, 20.0, "electron energy [eV]", strcat(ndir, "/energy"));
					Hist->fill(hoff+hmax*i+0,"electron energy_small", e[i]->energy(), 1., "electron energy", 360, 6.0, 14.0, "electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+1,"phi vs electron energy", e[i]->raw.phi, e[i]->energy(), 1., "electron energy vs. phi on detector", 140, -180.0, 180.0, "phi [deg]", 400, 0.0, 20.0, "electron energy [eV]", ndir);
					Hist->fill(hoff+hmax*i+2,"ctheta vs electron energy", e[i]->mom.Cos_Theta(), e[i]->energy(), 1., "electron energy vs. cos(theta_z)", 140, -1.0, 1.0, "ctheta", 400, 0.0, 20.0, "electron energy [eV]", ndir);
				}

				strcpy(ndir,tdir); 
				// Ar series
				if (cur_reaction->channel == 248) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-.1,.1,"p_x [a.u.]",220,-.1,.1,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-.1,.1,"p_x [a.u.]",220,-.1,.1,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-.1,.1,"p_y [a.u.]",220,-.1,.1,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 252) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.,1.,"p_x [a.u.]",220,-1.,1.,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.,1.,"p_x [a.u.]",220,-1.,1.,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.,1.,"p_y [a.u.]",220,-1.,1.,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,1.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 254) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.1,1.1,"p_x [a.u.]",220,-1.1,1.1,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.1,1.1,"p_x [a.u.]",220,-1.1,1.1,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.1,1.1,"p_y [a.u.]",220,-1.1,1.1,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,.8,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 256) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				} 
				else if (cur_reaction->channel == 260) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				// Molecular brake-ups series
				else if (cur_reaction->channel == 9) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-2.5,2.5,"p_x [a.u.]",220,-2.5,2.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-2.5,2.5,"p_x [a.u.]",220,-2.5,2.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-2.5,2.5,"p_y [a.u.]",220,-2.5,2.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 10) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 11) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 12) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 13) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 14) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-2.5,2.5,"p_x [a.u.]",220,-2.5,2.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-2.5,2.5,"p_x [a.u.]",220,-2.5,2.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-2.5,2.5,"p_y [a.u.]",220,-2.5,2.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else if (cur_reaction->channel == 14 && ID == 16) {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-1.5,1.5,"p_x [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-1.5,1.5,"p_y [a.u.]",220,-1.5,1.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}
				else {
					Hist->fill(hoff+hmax*i+3,"p_x vs p_y",e[i]->mom.x,e[i]->mom.y,1.,"p_x vs. p_y",220,-0.5,0.5,"p_x [a.u.]",220,-0.5,0.5,"p_y [a.u.]",strcat(ndir,"/momenta"));
					Hist->fill(hoff+hmax*i+4,"p_x vs p_z",e[i]->mom.x,e[i]->mom.z,1.,"p_x vs. p_z",220,-0.5,0.5,"p_x [a.u.]",220,-0.5,0.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+5,"p_y vs p_z",e[i]->mom.y,e[i]->mom.z,1.,"p_y vs. p_z",220,-0.5,0.5,"p_y [a.u.]",220,-0.5,0.5,"p_z [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+6,"p_mag",e[i]->mom.Mag(),1.,"|p|",50,0.0,2.5,"|p| [a.u.]",ndir);
				}

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+7,"phi",e[i]->mom.Phi_deg(),1.,"phi in the labframe",72,-180.0,180.0,"phi [deg]",strcat(ndir,"/angles_labframe"));
				Hist->fill(hoff+hmax*i+8,"ctheta",e[i]->mom.Cos_Theta(),1.,"cos(theta) in the labframe",72,-1.0,1.0,"cos(theta)",ndir);
				// labframe frame ADs for beta etc.
				if(cur_reaction->LF_cond->type > -1) {

					if((e[i]->energy() > cur_reaction->LF_cond->emin) && (e[i]->energy() < cur_reaction->LF_cond->emax)) {
						CH_vector photon_dir = CH_vector(1.0,0.0,0);
						CH_vector pol_dir = CH_vector(0,1.0,0);
						CH_vector perp_dir = CH_vector(0,0,1.0);

						Coordinate_System lf = Coordinate_System(photon_dir,pol_dir);
						CH_vector e_LF;
						e_LF = lf.project_vector(e[i]->mom);

						strcat(ndir,"/photo-electron");

						switch(cur_reaction->LF_cond->type) { 						
						case 0: // linear light 
							Hist->fill(hoff+hmax*i+9,"ctheta to pol (dipole)",cos(e[i]->mom.Angle(pol_dir)),1.,"ctheta /w resp. to light propagation",72,-1.0,1.0,"phi [deg]",ndir);	
							// We need to restrict to cases where we are in the plane photon_dir x pol_dir
							if(fabs(e[i]->mom.Angle_deg(perp_dir)-90.0)<20.0) {
								Hist->fill(hoff+hmax*i+10,"ctheta to light (non-dipole)",cos(e[i]->mom.Angle(photon_dir)),1.,"ctheta /w resp. to light propagation",72,-1.0,1.0,"phi [deg]",ndir);	
							}
							break;					
						case 1: // circular light
								// We can integrate over all planes which include the photon_dir
								//Hist->fill(hoff+hmax*i+11,"ctheta to light",cos(e[i]->mom.Angle(photon_dir)),1.,"cos(theta) /w resp. to light propagation",36,-1.0,1.0,"cos(theta)",ndir);
							Hist->fill(hoff + hmax*i + 11, "ctheta to light vs ee", cos(e[i]->mom.Angle(photon_dir)), e[i]->energy(), 1., "cos(theta) /w resp. to light propagation vs electron energy [eV]", 36, -1.0, 1.0, "cos(theta)", 100, 0.0, 15.0, "electron energy [eV]", ndir);
							break;
						}
						Hist->fill(hoff+hmax*i+12,"LFAD3D",e_LF.Phi_deg(),e_LF.Cos_Theta(),1.,"phi vs. cos(theta) in the photon frame",72,-180.0,180.0,"phi [deg]",72,-1.0,1.0,"cos(theta)",ndir);
					}
				}

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*i+13,"p_x vs tof (valid)",e[i]->raw.data.tof,e[i]->mom.x,1.,"B-check x",250,-25.,150.,"tof [ns]",200,-1.5,10.5,"p_x [a.u.]",strcat(ndir,"/raw/B-orientation"));
				Hist->fill(hoff+hmax*i+14,"p_y vs tof (valid)",e[i]->raw.data.tof,e[i]->mom.y,1.,"B-check y",250,-25.,150.,"tof [ns]",200,-1.5,10.5,"p_y [a.u.]",ndir);
			}

			strcpy(ndir,tdir);
			//Hist->fill(hoff+hmax*i+15,"position",e[i]->raw.data.x,e[i]->raw.data.y,1.,"position",400,-1.*edet_size*1.5,edet_size*1.5,"x [mm]",400,-1.*edet_size*1.5,edet_size*1.5,"y [mm]",strcat(ndir,"/raw"));
			Hist->fill(hoff+hmax*i+15,"position",e[i]->raw.data.x,e[i]->raw.data.y,1.,"position",400,-180,180,"x [mm]",400,-180,180,"y [mm]",strcat(ndir,"/raw"));
			Hist->fill(hoff+hmax*i+16,"tof",e[i]->raw.data.tof,1.,"time-of-flight",250,-25.,150.,"tof [ns]",ndir);
			Hist->fill(hoff+hmax*i+25,"tof_red",e[i]->raw.data.tof,1.,"time-of-flight",250,32.,44.,"tof [ns]",ndir);
			
			//GN
			//if (cur_reaction->channel == 14) {
			Hist->fill(hoff + hmax*i + 17, "wiggle", e[i]->raw.data.tof, (sqrt(e[i]->raw.data.x*e[i]->raw.data.x + e[i]->raw.data.y*e[i]->raw.data.y)), 1., "wiggles", 250, -25., 150., "tof [ns]", 200, -1., 90, "r [mm]", ndir);
			Hist->fill(hoff + hmax*i + 18, "fish x", e[i]->raw.data.tof, e[i]->raw.data.x, 1., "x-fish", 250, 0., 50., "tof [ns]", 200, -90, 90, "x [mm]", ndir);
			Hist->fill(hoff + hmax*i + 19, "fish y", e[i]->raw.data.tof, e[i]->raw.data.y, 1., "y-fish", 250, 0., 50., "tof [ns]", 200, -90, 90, "y [mm]", ndir);
			if (fabs(e[i]->raw.data.y) < 10.)
				Hist->fill(hoff + hmax*i + 20, "fish filet x", e[i]->raw.data.tof, e[i]->raw.data.x, 1., "x-fish filet", 250, 0., 50., "tof [ns]", 200, -90, 90, "x [mm]", ndir);
			if (fabs(e[i]->raw.data.x) < 10.)
				Hist->fill(hoff + hmax*i + 21, "fish filet y", e[i]->raw.data.tof, e[i]->raw.data.y, 1., "y-fish filet", 250, 0., 50., "tof [ns]", 200, -90, 90, "y [mm]", ndir);
			//}
			//else {
			//	Hist->fill(hoff + hmax*i + 17, "wiggle", e[i]->raw.data.tof, (sqrt(e[i]->raw.data.x*e[i]->raw.data.x + e[i]->raw.data.y*e[i]->raw.data.y)), 1., "wiggles", 250, -25., 150., "tof [ns]", 200, -1., edet_size, "r [mm]", ndir);
			//	Hist->fill(hoff + hmax*i + 18, "fish x", e[i]->raw.data.tof, e[i]->raw.data.x, 1., "x-fish", 250, 0., 50., "tof [ns]", 200, -1.*edet_size, edet_size, "x [mm]", ndir);
			//	Hist->fill(hoff + hmax*i + 19, "fish y", e[i]->raw.data.tof, e[i]->raw.data.y, 1., "y-fish", 250, 0., 50., "tof [ns]", 200, -1.*edet_size, edet_size, "y [mm]", ndir);
			//	if (fabs(e[i]->raw.data.y) < 10.)
			//		Hist->fill(hoff + hmax*i + 20, "fish filet x", e[i]->raw.data.tof, e[i]->raw.data.x, 1., "x-fish filet", 250, 0., 50., "tof [ns]", 200, -1.*edet_size, edet_size, "x [mm]", ndir);
			//	if (fabs(e[i]->raw.data.x) < 10.)
			//		Hist->fill(hoff + hmax*i + 21, "fish filet y", e[i]->raw.data.tof, e[i]->raw.data.y, 1., "y-fish filet", 250, 0., 50., "tof [ns]", 200, -1.*edet_size, edet_size, "y [mm]", ndir);
			//}
			//GN_!
			strcpy(ndir,tdir);
			Hist->fill(hoff+hmax*i+22,"p_x vs tof (all)",e[i]->raw.data.tof,e[i]->mom.x,1.,"B-check x",250,-25.,150.,"tof [ns]",200,-1.5,1.5,"p_x [a.u.]",strcat(ndir,"/raw/B-orientation"));
			Hist->fill(hoff+hmax*i+23,"p_y vs tof (all)",e[i]->raw.data.tof,e[i]->mom.y,1.,"B-check y",250,-25.,150.,"tof [ns]",200,-1.5,1.5,"p_y [a.u.]",ndir);			
		}

		Hist->fill(hoff+hmax+24,"ehit",ehit,1.,"number of hits electron",43,-1.25,20.25,"number of hits",rootdir);

		// electron multihit histograms 
		int e1 = 0;
		int e2 = 0;
		bool twoelec = false;
		CH_vector pesum;
		if(cur_reaction->two_elec) {
			strcpy(ndir,rootdir);
			sprintf(tdir,"/coincidence");
			strcat(ndir,tdir);
			strcpy(tdir,ndir);

			// look for first electron
			for(int i=0;i<ehit;i++) {
				if(e[i]->valid) {
					if(e[i]->energy()>cur_reaction->e_master_min && e[i]->energy()<cur_reaction->e_master_max) {
						e1 = i;
						twoelec = true;
						break;
					}
				}
			}
			// look for second electron
			if(twoelec) {
				twoelec = false;
				for(int i=0;i<ehit;i++) {
					if(e[i]->valid && i!=e1) {
						e2 = i;
						twoelec = true;
						break;
					}
				}
			}

			if(twoelec) {
				Hist->fill(hoff+hmax*15+25,"electron energy 1 vs 2",e[e1]->energy(),e[e2]->energy(),1.0,"electron energy 1st electron vs. 2nd electron",100,0.0,20.0,"energy electron 1 [eV]",100,0.0,20.0,"energy electron 2 [eV]",ndir);
				Hist->fill(hoff+hmax*15+26,"electron energy 1 vs 2 sym",e[e1]->energy(),e[e2]->energy(),0.5,"electron energy 1st electron vs. 2nd electron",100,0.0,20.0,"energy electron 1 [eV]",100,0.0,20.0,"energy electron 2 [eV]",ndir);
				Hist->fill(hoff+hmax*15+26,"electron energy 1 vs 2 sym",e[e2]->energy(),e[e1]->energy(),0.5,"electron energy 1st electron vs. 2nd electron",100,0.0,20.0,"energy electron 1 [eV]",100,0.0,20.0,"energy electron 2 [eV]",ndir);

				Hist->fill(hoff+hmax*15+27,"energy sharing",e[e1]->energy()/(e[e1]->energy() + e[e2]->energy()),0.5,"electron energy sharing",50,0.0,1.0,"energy sharing",ndir);
				Hist->fill(hoff+hmax*15+27,"energy sharing",e[e2]->energy()/(e[e1]->energy() + e[e2]->energy()),0.5,"electron energy sharing",50,0.0,1.0,"energy sharing",ndir);	

				if(cur_reaction->e_master_min<0.0) {
					Hist->fill(hoff+hmax*15+28,"energy sharing vs ctheta",cos(e[e1]->mom.Angle(e[e2]->mom)),e[e1]->energy()/(e[e1]->energy() + e[e2]->energy()),0.5,"electron energy sharing vs. angle between electrons",36,-1.0,1.0,"cos(theta_e1_e2)",50,0.0,1.0,"energy sharing",ndir);
					Hist->fill(hoff+hmax*15+28,"energy sharing vs ctheta",cos(e[e2]->mom.Angle(e[e1]->mom)),e[e2]->energy()/(e[e1]->energy() + e[e2]->energy()),0.5,"electron energy sharing vs. angle between electrons",36,-1.0,1.0,"cos(theta_e1_e2)",50,0.0,1.0,"energy sharing",ndir);	
				} else {
					Hist->fill(hoff+hmax*15+29,"energy sharing vs ctheta",cos(e[e2]->mom.Angle(e[e1]->mom)),e[e2]->energy()/(e[e1]->energy() + e[e2]->energy()),0.5,"electron energy sharing vs. angle between electrons",36,-1.0,1.0,"cos(theta_e1_e2)",50,0.0,1.0,"energy sharing",ndir);	
				}

				Hist->fill(hoff+hmax*15+30,"electron sum energy",e[e1]->energy()+e[e2]->energy(),1.0,"electron sum energy",100,0.0,25.0,"energy sum electron [eV]",ndir);

				strcpy(ndir,tdir);
				pesum = e[e1]->mom + e[e2]->mom;
				Hist->fill(hoff+hmax*15+31,"p_sum_x vs p_sum_y",pesum.x,pesum.y,1.,"p_sum_x vs. p_sum_y",90,-1.5,1.5,"p_sum_x [a.u.]",90,-1.5,1.5,"p_sum_y [a.u.]",strcat(ndir,"/electron sum momenta"));
				Hist->fill(hoff+hmax*15+32,"p_sum_x vs p_sum_z",pesum.x,pesum.z,1.,"p_sum_x vs. p_sum_z",90,-1.5,1.5,"p_sum_x [a.u.]",90,-1.5,1.5,"p_sum_z [a.u.]",ndir);
				Hist->fill(hoff+hmax*15+33,"p_sum_y vs p_sum_z",pesum.y,pesum.z,1.,"p_sum_y vs. p_sum_z",90,-1.5,1.5,"p_sum_y [a.u.]",90,-1.5,1.5,"p_sum_z [a.u.]",ndir);				
				Hist->fill(hoff+hmax*15+34,"p_sum_mag",pesum.Mag(),1.,"|p_sum|",50,0.0,2.5,"|p_sum| [a.u.]",ndir);

				strcpy(ndir,tdir);
				// Create Newton/elec frame plots 
				Coordinate_System eframe = Coordinate_System(e[e1]->mom, e[e2]->mom);
				CH_vector eFrame_mom;
				eFrame_mom = eframe.project_vector(e[2]->mom);
				Hist->fill(hoff+hmax*15+35,"electron frame norm. all",eFrame_mom.z/e[e1]->mom.Mag(),eFrame_mom.x/e[e1]->mom.Mag(),1.,"p_e2 relative to p_e1 (normalized to e1)",80,-2.5,1.5,"pz",60,0.0,3.0,"py",strcat(ndir,"/relative momenta"));
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())<0.25)
					Hist->fill(hoff+hmax*15+36,"electron frame norm. slow",eFrame_mom.z/e[e1]->mom.Mag(),eFrame_mom.x/e[e1]->mom.Mag(),1.,"p_e2 relative to p_e1 (norm. to e1), e1 slow",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())>0.25 && e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())<0.75)
					Hist->fill(hoff+hmax*15+37,"electron frame norm. mid",eFrame_mom.z/e[e1]->mom.Mag(),eFrame_mom.x/e[e1]->mom.Mag(),1.,"p_e2 relative to p_e1 (norm. to e1), ee1=ee2",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())>0.75)
					Hist->fill(hoff+hmax*15+38,"electron frame norm. fast",eFrame_mom.z/e[e1]->mom.Mag(),eFrame_mom.x/e[e1]->mom.Mag(),1.,"p_e2 relative to p_e1 (norm. to e1), e1 fast",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);

				Hist->fill(hoff+hmax*15+39,"electron frame all",eFrame_mom.z,eFrame_mom.x,1.,"p_e2 relative to p_e1",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())<0.25)
					Hist->fill(hoff+hmax*15+40,"electron frame slow",eFrame_mom.z,eFrame_mom.x,1.,"p_e2 relative to p_e1, e1 slow",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())>0.25 && e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())<0.75)
					Hist->fill(hoff+hmax*15+41,"electron frame mid",eFrame_mom.z,eFrame_mom.x,1.,"p_e2 relative to p_e1, ee1=ee2",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);
				if(e[e1]->energy()/(e[e1]->energy() + e[e2]->energy())>0.75)
					Hist->fill(hoff+hmax*15+42,"electron frame fast",eFrame_mom.z,eFrame_mom.x,1.,"p_e2 relative to p_e1, e1 fast",80,-2.5,1.5,"pz",60,0.0,3.0,"py",ndir);

				strcpy(ndir,tdir);
				Hist->fill(hoff+hmax*15+43,"dt dx",e[e2]->raw.data.time-e[e1]->raw.data.time,e[e2]->raw.data.x-e[e1]->raw.data.x,1.,"Multihit check dt vs. dx",88,-2.0,20,"dt [ns]",120,-30.0,30.0,"dx [mm]",strcat(ndir,"/multihit checks"));
				Hist->fill(hoff+hmax*15+44,"dt dy",e[e2]->raw.data.time-e[e1]->raw.data.time,e[e2]->raw.data.y-e[e1]->raw.data.y,1.,"Multihit check dt vs. dy",88,-2.0,20,"dt [ns]",120,-30.0,30.0,"dy [mm]",ndir);
				double dR=sqrt(pow((e[e2]->raw.data.x - e[e1]->raw.data.x),2)-pow((e[e2]->raw.data.y - e[e1]->raw.data.y),2)); 
				Hist->fill(hoff+hmax*15+45,"dt dR",e[e2]->raw.data.time-e[e1]->raw.data.time,dR,1.,"Multihit check dt vs. dR",88,-2.0,20,"dt [ns]",60,0.0,30.0,"dR [mm]",ndir);

				// two electron labframe frame ADs for beta etc.
				if (cur_reaction->LF_cond->type > -1) {
					strcpy(ndir, tdir);
					double esum = e[e1]->energy() + e[e2]->energy();
					CH_vector DiElec = e[e1]->mom + e[e2]->mom;
					if ((esum > cur_reaction->LF_cond->emin) && (esum < cur_reaction->LF_cond->emax)) {
						CH_vector photon_dir = CH_vector(1.0, 0.0, 0);
						CH_vector pol_dir = CH_vector(0, 1.0, 0);
						CH_vector perp_dir = CH_vector(0, 0, 1.0);

						Coordinate_System lf = Coordinate_System(photon_dir, pol_dir);
						CH_vector e_LF;
						e_LF = lf.project_vector(DiElec);

						strcat(ndir, "/angles_labframe");

						switch (cur_reaction->LF_cond->type) {
						case 0: // linear light 
							Hist->fill(hoff + hmax*15 + 46, "phi to pol (dipole)", DiElec.Angle_deg(pol_dir), 1., "phi /w resp. to light propagation", 72, -180.0, 180.0, "phi [deg]", ndir);
							// We need to restrict to cases where we are in the plane photon_dir x pol_dir
							if (fabs(DiElec.Angle_deg(perp_dir))<20.0) {
								Hist->fill(hoff + hmax*15 + 47, "phi to light (non-dipole)", DiElec.Angle_deg(photon_dir), 1., "phi /w resp. to light propagation", 72, -180.0, 180.0, "phi [deg]", ndir);
							}
							break;
						case 1: // circular light
								// We can integrate over all planes which include the photon_dir
							Hist->fill(hoff + hmax*15 + 48, "ctheta to light", cos(DiElec.Angle(photon_dir)), 1., "cos(theta) /w resp. to light propagation", 24, -1.0, 1.0, "cos(theta)", ndir);
							break;
						}
						Hist->fill(hoff + hmax*15 + 49, "LFAD3D", e_LF.Phi_deg(), e_LF.Cos_Theta(), 1., "phi vs. cos(theta) in the photon frame", 36, -180.0, 180.0, "phi [deg]", 24, -1.0, 1.0, "cos(theta)", ndir);
					}
				}
			}
		}

		if(cur_reaction->photon_scan) {
			strcpy(ndir,reaction_dir);
			sprintf(tdir,"/photon energy scan/electrons");
			strcat(ndir,tdir);
			strcpy(tdir,ndir);

			int phbins = int((cur_reaction->ph_max - cur_reaction->ph_min)/cur_reaction->ph_step);

			// We will plot only the first valid electron here, so find it... 
			int i;
			for(i=0;i<(int)ehit;i++) {
				if(e[i]->valid)
					break;
			}

			Hist->fill(hoff+hmax*15+50,"hv vs p",scan_val[cur_reaction->ph_scan_channel],e[i]->mom.Mag(),1.,"Photon energy vs. |p|",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",50,0.0,2.5,"|p| [a.u.]",ndir);
			Hist->fill(hoff+hmax*15+51,"hv vs energy",scan_val[cur_reaction->ph_scan_channel],e[i]->energy(),1.,"Photon energy vs. electron energy",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",50,0.0,20.0,"electron energy [eV]",ndir);
			if(twoelec) {
				Hist->fill(hoff+hmax*15+52,"hv vs psum",scan_val[cur_reaction->ph_scan_channel],pesum.Mag(),1.,"Photon energy vs. |p|",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",50,0.0,2.5,"|psum| [a.u.]",ndir);
				Hist->fill(hoff+hmax*15+53,"hv vs electron sum energy",scan_val[cur_reaction->ph_scan_channel],e[e1]->energy()+e[e2]->energy(),1.,"Photon energy vs. electron sum energy",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",50,0.0,20.0,"electron sum energy [eV]",ndir);			
			}
		}
	}
}
#include "CH_Histograms.h"

// Standard histograms for reactions of "DIATOMIC", "DIATOMIC_DISS", "DIATOMIC_ELEC" and "DIATOMIC_DISS" type
namespace CH
{
	// there is nothing to plot in case we did not have a reaction definition
	void histograms_class::plot_diatomic(reaction_struct *cur_reaction, diatomic_class * mol, electron_class ** e, int ehit, double *scan_val)
	{
		histo_handler *Hist = this->Hist_ions;
		// histogram index stuff
		int hmax = HISTS_PER_CHANNEL; // overall number of ColAHelL histograms per reaction/channel
		int hoff = hmax * 15 * (int)cur_reaction->reac_num;
		int CH = cur_reaction->channel;

		// limit to 16 hits.. 
		if (ehit > 15)
			ehit = 15;

		char reaction_dir[360];
		//		char x_axis_title[80];
		//		char y_axis_title[80];

		char rootdir[360];
		char ndir[360];
		char alldir[360];
		char tdir[360];
		char talldir[360];

		if (this->use_master_folder) {
			strcpy(reaction_dir, this->CH_master_folder);
			strcat(reaction_dir, cur_reaction->name);
		}
		else
			strcpy(reaction_dir, cur_reaction->name);

		strcpy(rootdir, reaction_dir);

		// fill standard ion histograms.
		ion_class *r[2];

		r[0] = mol->ion[0];
		r[1] = mol->ion[1];

		for (int i = 0; i < 2; i++) {
			strcpy(ndir, rootdir);
			sprintf(tdir, "/ion_%i", i);
			strcat(ndir, tdir);
			strcpy(tdir, ndir);

			Hist->fill(hoff + hmax*i + 0, "ion energy", r[i]->energy(), 1., "ion energy", 50, 0.0, 25.0, "ion energy [eV]", strcat(ndir, "/energy"));
			Hist->fill(hoff + hmax*i + 1, "phi vs ion energy", r[i]->raw.phi, r[i]->energy(), 1., "ion energy vs. phi on detector", 148, -180.0, 180.0, "phi [deg]", 150, 0.0, 10.0, "ion energy [eV]", ndir);
			Hist->fill(hoff + hmax*i + 2, "ctheta vs ion energy", r[i]->mom.Cos_Theta(), r[i]->energy(), 1., "ion energy vs. cos(theta_z)", 148, -1.0, 1.0, "ctheta", 150, 0.0, 10.0, "ion energy [eV]", ndir);

			strcpy(ndir, tdir);
			Hist->fill(hoff + hmax*i + 3, "p_x vs p_y", r[i]->mom.x, r[i]->mom.y, 1., "p_x vs. p_y", 200, -200., 200., "p_x [a.u.]", 200, -200., 200., "p_y [a.u.]", strcat(ndir, "/momenta"));
			Hist->fill(hoff + hmax*i + 4, "p_x vs p_z", r[i]->mom.x, r[i]->mom.z, 1., "p_x vs. p_z", 200, -200., 200., "p_x [a.u.]", 200, -200., 200., "p_z [a.u.]", ndir);
			Hist->fill(hoff + hmax*i + 5, "p_y vs p_z", r[i]->mom.y, r[i]->mom.z, 1., "p_y vs. p_z", 200, -200., 200., "p_y [a.u.]", 200, -200., 200., "p_z [a.u.]", ndir);
			Hist->fill(hoff + hmax*i + 6, "p_mag", r[i]->mom.Mag(), 1., "|p|", 60, 20., 200., "|p| [a.u.]", ndir);

			strcpy(ndir, tdir);
			Hist->fill(hoff + hmax*i + 7, "phi", r[i]->mom.Phi_deg(), 1., "phi in the labframe", 72, -180.0, 180.0, "phi [deg]", strcat(ndir, "/angles_labframe"));
			Hist->fill(hoff + hmax*i + 8, "ctheta", r[i]->mom.Cos_Theta(), 1., "cos(theta) in the labframe", 72, -1.0, 1.0, "cos(theta)", ndir);

			strcpy(ndir, tdir);
			Hist->fill(hoff + hmax*i + 9, "position", r[i]->raw.data.x, r[i]->raw.data.y, 1., "position", 400, -1.*rdet_size, rdet_size, "x [mm]", 400, -1.*rdet_size, rdet_size, "y [mm]", strcat(ndir, "/raw"));
			Hist->fill(hoff + hmax*i + 10, "tof", r[i]->raw.data.tof, 1., "time-of-flight", 4000, -1000., 20000., "tof [ns]", ndir);
			Hist->fill(hoff + hmax*i + 11, "wiggle", r[i]->raw.data.tof, (sqrt(r[i]->raw.data.x*r[i]->raw.data.x + r[i]->raw.data.y*r[i]->raw.data.y)), 1., "wiggles", 500, -25., 10000., "tof [ns]", 200, -1., rdet_size, "r [mm]", ndir);
			Hist->fill(hoff + hmax*i + 12, "fish x", r[i]->raw.data.tof, r[i]->raw.data.x, 1., "x-fish", 500, -25., 10000., "tof [ns]", 200, -1.*rdet_size, rdet_size, "x [mm]", ndir);
			Hist->fill(hoff + hmax*i + 13, "fish y", r[i]->raw.data.tof, r[i]->raw.data.y, 1., "y-fish", 500, -25., 10000., "tof [ns]", 200, -1.*rdet_size, rdet_size, "y [mm]", ndir);
			if (fabs(r[i]->raw.data.y) < 10.)
				Hist->fill(hoff + hmax*i + 14, "fish filet x", r[i]->raw.data.tof, r[i]->raw.data.x, 1., "x-fish filet", 500, -25., 10000., "tof [ns]", 200, -1.*rdet_size, rdet_size, "x [mm]", ndir);
			if (fabs(r[i]->raw.data.x) < 10.)
				Hist->fill(hoff + hmax*i + 15, "fish filet y", r[i]->raw.data.tof, r[i]->raw.data.y, 1., "y-fish filet", 500, -25., 10000., "tof [ns]", 200, -1.*rdet_size, rdet_size, "y [mm]", ndir);

		}

		// fill standard diatomic histograms
		strcpy(ndir, rootdir);
		sprintf(tdir, "/mol_ion");
		strcat(ndir, tdir);
		strcpy(tdir, ndir);

		Hist->fill(hoff + hmax + 16, "p_relx vs p_rely", mol->mom_rel.x, mol->mom_rel.y, 1., "relative momentum x vs y", 200, -200., 200., "p_relx [a.u.]", 200, -200., 200., "p_rely [a.u.]", strcat(ndir, "/momenta"));
		Hist->fill(hoff + hmax + 17, "p_relx vs p_relz", mol->mom_rel.x, mol->mom_rel.z, 1., "relative momentum x vs z", 200, -200., 200., "p_relx [a.u.]", 200, -200., 200., "p_relz [a.u.]", ndir);
		Hist->fill(hoff + hmax + 18, "p_rely vs p_relz", mol->mom_rel.y, mol->mom_rel.z, 1., "relative momentum y vs z", 200, -200., 200., "p_rely [a.u.]", 200, -200., 200., "p_relz [a.u.]", ndir);
		if (CH == 9 || CH == 11) {
			Hist->fill(hoff + hmax + 19, "p_sumx vs p_sumy", mol->mom_cm.x, mol->mom_cm.y, 1., "center of mass momentum x vs y", 100, -120., 120., "p_cmx [a.u.]", 100, -120., 120., "p_cmy [a.u.]", ndir);
			Hist->fill(hoff + hmax + 20, "p_sumx vs p_sumz", mol->mom_cm.x, mol->mom_cm.z, 1., "center of mass momentum x vs z", 100, -120., 120., "p_cmx [a.u.]", 100, -120., 120., "p_cmz [a.u.]", ndir);
			Hist->fill(hoff + hmax + 21, "p_sumy vs p_sumz", mol->mom_cm.y, mol->mom_cm.z, 1., "center of mass momentum y vs z", 100, -120., 120., "p_cmy [a.u.]", 100, -120., 120., "p_cmz [a.u.]", ndir);
		}
		else {
			Hist->fill(hoff + hmax + 19, "p_sumx vs p_sumy", mol->mom_cm.x, mol->mom_cm.y, 1., "center of mass momentum x vs y", 100, -30., 30., "p_cmx [a.u.]", 100, -30., 30., "p_cmy [a.u.]", ndir);
			Hist->fill(hoff + hmax + 20, "p_sumx vs p_sumz", mol->mom_cm.x, mol->mom_cm.z, 1., "center of mass momentum x vs z", 100, -30., 30., "p_cmx [a.u.]", 100, -30., 30., "p_cmz [a.u.]", ndir);
			Hist->fill(hoff + hmax + 21, "p_sumy vs p_sumz", mol->mom_cm.y, mol->mom_cm.z, 1., "center of mass momentum y vs z", 100, -30., 30., "p_cmy [a.u.]", 100, -30., 30., "p_cmz [a.u.]", ndir);
		}
		Hist->fill(hoff+hmax+22,"ion1_x vs ion2_x",mol->ion[0]->mom.x,mol->ion[1]->mom.x,1.,"momentum ion1 vs ion2 x-direction",150,-200.,200.,"p_1x [a.u.]",150,-200.,200.,"p_2x [a.u.]",ndir);
		Hist->fill(hoff+hmax+23,"ion1_y vs ion2_y",mol->ion[0]->mom.y,mol->ion[1]->mom.y,1.,"momentum ion1 vs ion2 y-direction",150,-200.,200.,"p_1y [a.u.]",150,-200.,200.,"p_2y [a.u.]",ndir);
		Hist->fill(hoff+hmax+24,"ion1_z vs ion2_z",mol->ion[0]->mom.z,mol->ion[1]->mom.z,1.,"momentum ion1 vs ion2 z-direction",150,-200.,200.,"p_1z [a.u.]",150,-200.,200.,"p_2z [a.u.]",ndir);

		Hist->fill(hoff+hmax+25,"p_relx vs p_sumz",mol->mom_rel.x,mol->mom_cm.z,1.,"rel. mom. x vs cm. mom. z",200,-200.,200.,"p_relx [a.u.]",100,-20.,20.,"p_cmz [a.u.]",ndir);
		Hist->fill(hoff+hmax+26,"p_rely vs p_sumz",mol->mom_rel.y,mol->mom_cm.z,1.,"rel. mom. y vs cm. mom. z",200,-200.,200.,"p_rely [a.u.]",100,-20.,20.,"p_cmz [a.u.]",ndir);

		Hist->fill(hoff+hmax+27,"p_sum vs KER",mol->mom_cm.Mag(),mol->KER(),1.,"center of mass momentum vs. KER",150,0.0,80.0,"p_cm [a.u.]",150,0.0,10.0,"KER [eV]",ndir); 

		strcpy(ndir,tdir);
		Hist->fill(hoff+hmax+28,"phi vs KER",mol->mom_rel.Phi_deg(),mol->KER(),1.,"phi rel. mom. vs. KER",148,-180.0,180.0,"phi [deg]",150,0.0,10.0,"KER [eV]",ndir); 
		Hist->fill(hoff+hmax+29,"ctheta vs KER",mol->mom_rel.Cos_Theta(),mol->KER(),1.,"cos(rel. mom.) vs. KER",148,-1.0,1.0,"ctheta",150,0.0,10.0,"KER [eV]",ndir); 
		Hist->fill(hoff+hmax+30,"KER",mol->KER(),1.,"kinetic energy release",400,0.,10.,"KER [eV]",ndir);

		// fill standard diatomic coincidence histograms
		if(ehit>0) {
			strcpy(ndir,rootdir);
			sprintf(tdir,"/mol_coincidence");
			strcat(ndir,tdir);
			strcpy(tdir,ndir);

			strcpy(alldir,rootdir);
			strcat(alldir,"/mol_coincidence/all_hits");
			strcpy(talldir,alldir);

			double e_sum_en = 0.0;
			int valid_electrons = 0;

			for(int i=0;i<(int)ehit;i++) {

				strcpy(ndir,rootdir);
				sprintf(tdir,"/mol_coincidence/elec_hit_%i",i);
				strcat(ndir,tdir);
				strcpy(tdir,ndir);

				if(e[i]->valid) {
					e_sum_en += e[i]->energy();
					valid_electrons++;

					Hist->fill(hoff+hmax*i+31,"KER vs electron energy",mol->KER(),e[i]->energy(),1.,"kinetic energy release vs electron energy",400,0.,20.,"KER [eV]",100,0.0,25.0,"electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*15+32,"KER vs electron energy",mol->KER(),e[i]->energy(),1.,"kinetic energy release vs electron energy",400,0.,20.,"KER [eV]",100,0.0,25.0,"electron energy [eV]",alldir);

					Hist->fill(hoff+hmax*i+33,"KER vs KER+electron energy",mol->KER(),mol->KER()+e[i]->energy(),1.,"kinetic energy release vs sum energy",400,0.,20.,"KER [eV]",100,0.0,35.0,"KER + electron energy [eV]",ndir);
					Hist->fill(hoff+hmax*15+34,"KER vs KER+electron energy",mol->KER(),mol->KER()+e[i]->energy(),1.,"kinetic energy release vs sum energy",400,0.,20.,"KER [eV]",100,0.0,35.0,"KER + electron energy [eV]",alldir);

					Hist->fill(hoff+hmax*i+35,"p_sumx vs p_ex",mol->mom_cm.x,e[i]->mom.x,1.,"psum_x vs. p_ex",50,-50.,50.,"psum_x [a.u.]",50,-1.5,1.5,"pe_x [a.u.]",strcat(ndir,"/mol_psum"));
					Hist->fill(hoff+hmax*i+36,"p_sumy vs p_ey",mol->mom_cm.y,e[i]->mom.y,1.,"psum_y vs. p_ey",50,-50.,50.,"psum_y [a.u.]",50,-1.5,1.5,"pe_y [a.u.]",ndir);
					Hist->fill(hoff+hmax*i+37,"p_sumz vs p_ez",mol->mom_cm.z,e[i]->mom.z,1.,"psum_z vs. p_ez",50,-50.,50.,"psum_z [a.u.]",50,-1.5,1.5,"pe_z [a.u.]",ndir);

					Hist->fill(hoff+hmax*15+38,"p_sumx vs p_ex",mol->mom_cm.x,e[i]->mom.x,1.,"psum_x vs. p_ex",50,-50.,50.,"psum_x [a.u.]",50,-1.5,1.5,"pe_x [a.u.]",strcat(alldir,"/mol_psum"));
					Hist->fill(hoff+hmax*15+39,"p_sumy vs p_ey",mol->mom_cm.y,e[i]->mom.y,1.,"psum_y vs. p_ey",50,-50.,50.,"psum_y [a.u.]",50,-1.5,1.5,"pe_y [a.u.]",alldir);
					Hist->fill(hoff+hmax*15+40,"p_sumz vs p_ez",mol->mom_cm.z,e[i]->mom.z,1.,"psum_z vs. p_ez",50,-50.,50.,"psum_z [a.u.]",50,-1.5,1.5,"pe_z [a.u.]",alldir);

					strcpy(ndir,tdir);
					strcpy(alldir,talldir);

					// molecular frame
					if(cur_reaction->MF_cond->type > -1) {

						CH_vector photon_dir = CH_vector(1.0,0.0,0);
						CH_vector pol_dir = CH_vector(0,1.0,0);
						Coordinate_System mf = Coordinate_System(mol->mom_rel, pol_dir);
						Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
						CH_vector e_MF;

						strcat(ndir,"/mol_frame");
						strcpy(tdir,ndir);
						strcat(alldir,"/mol_frame");

						switch(cur_reaction->MF_cond->type) { 
						case 0:// linear light
							e_MF = mf.project_vector(e[i]->mom);

							if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
								Hist->fill(hoff+hmax*i+41,"MFPAD3D all orientations",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);
								Hist->fill(hoff+hmax*15+42,"MFPAD3D all orientations",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",alldir);
							}
							if(mol->mom_rel.Angle_deg(pol_dir) < 20) {
								if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									Hist->fill(hoff+hmax*i+43,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(ndir,"/parallel to polarization"));
									Hist->fill(hoff+hmax*i+44,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",ndir);

									Hist->fill(hoff+hmax*15+45,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(alldir,"/parallel to polarization"));
									Hist->fill(hoff+hmax*15+46,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",alldir);
								}
								if((e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									Hist->fill(hoff+hmax*i+47,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",ndir);

									Hist->fill(hoff+hmax*15+48,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",alldir);
								}
							}
							strcpy(ndir,tdir);
							strcpy(alldir,talldir);
							if(fabs(mol->mom_rel.Angle_deg(pol_dir)-90) < 20) {
								if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									Hist->fill(hoff+hmax*i+49,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(ndir,"/perpendicular to polarization"));

									Hist->fill(hoff+hmax*15+50,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(alldir,"/perpendicular to polarization"));			
									if(fabs(e_MF.Phi_deg())<15 || e_MF.Phi_deg()>185 || e_MF.Phi_deg() < -165) {
										Hist->fill(hoff+hmax*i+51,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",ndir);

										Hist->fill(hoff+hmax*15+52,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",alldir);
									}
								}
								if((e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									if(fabs(e_MF.Phi_deg())<15 || e_MF.Phi_deg()>185 || e_MF.Phi_deg() < -165) {
										Hist->fill(hoff+hmax*i+53,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",ndir);

										Hist->fill(hoff+hmax*15+54,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",alldir);
									}
								}
							}
							break;
						case 1:// circular light				
							e_MF = mf_circ.project_vector(e[i]->mom);

							if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
								Hist->fill(hoff+hmax*i+55,"MFPAD3D all orientations",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);
								Hist->fill(hoff+hmax*15+56,"MFPAD3D all orientations",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",alldir);
							}
							if((fabs(mol->mom_rel.Angle_deg(photon_dir)-90.0) < 15)) {
								if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									Hist->fill(hoff+hmax*i+57,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",48,-180.,180.,"Phi",36,-1.0,1.0,"cos(theta)",strcat(ndir,"/in polarization plane"));

									Hist->fill(hoff+hmax*15+58,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",48,-180.,180.,"Phi",36,-1.0,1.0,"cos(theta)",strcat(alldir,"/in polarization plane"));
									if(fabs(e[i]->mom.Angle_deg(photon_dir)-90.0) < 15) {
										Hist->fill(hoff+hmax*i+59,"MFPAD",e_MF.Phi_deg(),1.,"molecular frame angular distribution",72,-180.,180.,"Phi",ndir);						

										Hist->fill(hoff+hmax*15+60,"MFPAD",e_MF.Phi_deg(),1.,"molecular frame angular distribution",72,-180.,180.,"Phi",alldir);						
									}
								}
								if((e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									if(fabs(e[i]->mom.Angle_deg(photon_dir)-90.0) < 15) {
										Hist->fill(hoff+hmax*i+61,"KER vs. MFPAD",mol->KER(),e_MF.Phi_deg(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-180.,180.,"Phi",ndir);

										Hist->fill(hoff+hmax*15+62,"KER vs. MFPAD",mol->KER(),e_MF.Phi_deg(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-180.,180.,"Phi",alldir);																	
									}
								}
							}
							strcpy(ndir,tdir);
							strcpy(alldir,talldir);
							if((mol->mom_rel.Angle_deg(photon_dir) < 20)) {
								if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									Hist->fill(hoff+hmax*i+63,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(ndir,"/perpendicular to polarization plane"));

									Hist->fill(hoff+hmax*15+64,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",strcat(alldir,"/perpendicular to polarization plane"));
									if(fabs(e_MF.Phi_deg())<15 || e_MF.Phi_deg()>185 || e_MF.Phi_deg() < -165) {
										Hist->fill(hoff+hmax*i+65,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",ndir);

										Hist->fill(hoff+hmax*15+66,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",alldir);
									}
								}
								if((e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
									if(fabs(e_MF.Phi_deg())<15 || e_MF.Phi_deg()>185 || e_MF.Phi_deg() < -165) {
										Hist->fill(hoff+hmax*i+67,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",ndir);

										Hist->fill(hoff+hmax*15+68,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",alldir);
									}
								}
							}
							break;
						case 2: // Auger electron
							e_MF = mf.project_vector(e[i]->mom);
							if((mol->KER() > cur_reaction->MF_cond->rmin) && (mol->KER() < cur_reaction->MF_cond->rmax) && (e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
								Hist->fill(hoff+hmax*i+69,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);
								Hist->fill(hoff+hmax*i+70,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",ndir);

								Hist->fill(hoff+hmax*15+71,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",alldir);
								Hist->fill(hoff+hmax*15+72,"MFPAD",e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-1.,1.,"cos(theta)",alldir);
							}
							if((e[i]->energy() > cur_reaction->MF_cond->emin) && (e[i]->energy() < cur_reaction->MF_cond->emax)) {
								Hist->fill(hoff+hmax*i+73,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",ndir);

								Hist->fill(hoff+hmax*15+74,"KER vs. MFPAD",mol->KER(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution vs. KER",150,0.,15.,"KER[eV]",36,-1.,1.,"cos(theta)",alldir);
							}	
							break;
						}
					}
				}
			}
			strcpy(ndir,rootdir); sprintf(tdir,"/mol_coincidence");	strcat(ndir,tdir); strcpy(tdir,ndir); 
			strcpy(alldir,rootdir);	strcat(alldir,"/mol_coincidence/all_hits");	strcpy(talldir,alldir);

			if (valid_electrons > 1) {
				Hist->fill(hoff + hmax * 15 + 75, "KER vs electron sum energy", mol->KER(), e_sum_en, 1., "kinetic energy release vs electron sum energy", 200, 0., 20., "KER [eV]", 200, 0.0, 50.0, "esum energy [eV]", alldir);
				Hist->fill(hoff + hmax * 15 + 76, "KER vs total kin. energy", mol->KER(), e_sum_en + mol->KER(), 1., "kinetic energy release vs total kinetic energy", 200, 0., 20., "KER [eV]", 200, 0.0, 50.0, "total energy [eV]", alldir);
				Hist->fill(hoff + hmax * 15 + 77, "Total kin. energy", e_sum_en + mol->KER(), 1., "Total kinetic energy", 200, 0.0, 50.0, "total energy [eV]", alldir);
			}
		}

		if(cur_reaction->photon_scan) {
			int phbins = int((cur_reaction->ph_max - cur_reaction->ph_min)/cur_reaction->ph_step);

			strcpy(ndir,reaction_dir);
			sprintf(tdir,"/photon energy scan");
			strcat(ndir,tdir);
			Hist->fill(hoff+hmax*15+78,"hv",scan_val[cur_reaction->ph_scan_channel],1.,"Photon energy",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",ndir);

			strcpy(ndir,reaction_dir);
			sprintf(tdir,"/photon energy scan/ions");
			strcat(ndir,tdir);
			strcpy(tdir,ndir);
			Hist->fill(hoff+hmax*15+79,"hv vs KER",scan_val[cur_reaction->ph_scan_channel],mol->KER(),1.,"Photon energy vs. KER",phbins, cur_reaction->ph_min, cur_reaction->ph_max,"hv [eV]",80,0.0,20.0,"KER [eV]",ndir);
		}
	}
}
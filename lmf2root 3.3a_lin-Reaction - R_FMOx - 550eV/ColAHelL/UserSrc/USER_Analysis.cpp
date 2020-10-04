#include "..\Src\CH_Histograms.h"
#include "assert.h"
#include "Math.h"
#include <vector>
#include ".\Src\functions.h"
//TEST and MFPDA cannot be 1 both at the same time becuase of AUTO__LINE limitations
#define MFPAD 1
#define TEST 0

namespace CH
{
#define AUTO __LINE__ //returns the line number in the current file. Can be used to make life easier for histogram numbering.

	void histograms_class::UserAnalysis(reaction_struct *cur_reaction, electron_class ** e, ion_class ** r, diatomic_class * mol, polyatomic_class * big_mol, CH_event_struct *evt) {

		int ID = cur_reaction->ID;
		int CH = cur_reaction->channel;

		CH_vector photon_dir_till = CH_vector(1.0,0,0); ///light dir
		CH_vector photon_dir; photon_dir.x = 1.; photon_dir.y = 0.; photon_dir.z = 0.;
		double d_photon_dir[3]; d_photon_dir[0]=1.; d_photon_dir[1]=0.; d_photon_dir[2]=0.;
		double d_el[3];   d_el[0]=e[0]->mom.x;     d_el[1]=e[0]->mom.y;     d_el[2]=e[0]->mom.z;
		double d_sum[3]; d_sum[0]=-mol->mom_cm.x; d_sum[1]=-mol->mom_cm.y; d_sum[2]=-mol->mom_cm.z;

		double d_frag1[3];  d_frag1[0]=mol->ion[0]->mom.x; d_frag1[1]=mol->ion[0]->mom.y; d_frag1[2]=mol->ion[0]->mom.z;
		CH_vector frag1;       frag1.x=mol->ion[0]->mom.x;    frag1.y=mol->ion[0]->mom.y;    frag1.z=mol->ion[0]->mom.z;
		double d_frag2[3];  d_frag2[0]=mol->ion[1]->mom.x; d_frag2[1]=mol->ion[1]->mom.y; d_frag2[2]=mol->ion[1]->mom.z;
		CH_vector frag2;       frag2.x=mol->ion[1]->mom.x;    frag2.y=mol->ion[1]->mom.y;    frag2.z=mol->ion[1]->mom.z;

		CH_vector mol_axsis; mol_axsis.x = mol->ion[0]->mom.x - mol->ion[1]->mom.x; mol_axsis.y = mol->ion[0]->mom.y - mol->ion[1]->mom.y; mol_axsis.z = mol->ion[0]->mom.z - mol->ion[1]->mom.z;

		CH_vector pol_dir_till = CH_vector(0, 1.0, 0);
		CH_vector pol_dir; pol_dir.x = 0.; pol_dir.y = 1.; pol_dir.z = 0.;

		CH_vector e_new;
		//Correction to the e mom z due to momentum conservation (in a.u. )
		e_new.z = copysign(e[0]->mom.z,e[0]->mom.z) * sqrt(11.66*2/27.211-e[0]->mom.x*e[0]->mom.x-e[0]->mom.y*e[0]->mom.y);

		//From standard diatomic
		Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
		CH_vector e_MF;
		e_MF = mf_circ.project_vector(e[0]->mom);

		// char reaction_dir[800];
		// char rootdir[800];
		char ndir[800];
		char tdir[800];
		char sdir[800];
		char title[200];
		char title_el[200];
		int k = 1; //used to indetify the coordiante system from Kilian

				// Number of possible user histograms.
		MaxUserHists = 10000;

		histo_handler *Hist = this->Hist_User;
		int ehit = evt->e.num_hits;
		int rhit = evt->r.num_hits;

		double relphi;
		double phi_p, phi_m;

		//VALUES FOR MFPADS
		//NOTE: the light binnign is ALWAYS (6,12) cos(theta),phi
		int nbin_phi = 24;
		int nbin_cos = 12;

		//SET HERE THE VALUE!!
		double emin_el = 9.;
		double emax_el = 14.;

		//to limt the momentum of the neutral
		double sum_min = 10.;
		double sum_max = 55.;
		//
		double phi_yz = atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI;

		// example: process reaction with ID 3, which is a diatomic

		if (ID == 424) {
			strcpy(ndir, "check_orientation/He_24.63eV");
			Hist->fill(AUTO, "orientation_check_x", r[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -1.5, 1.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_y", r[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -1.5, 1.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_z", r[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -1.5, 1.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_x", r[0]->mom.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -0.5, 0.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_y", r[0]->mom.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -0.5, 0.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_z", r[0]->mom.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -0.5, 0.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
		}
		if (ID == 427) {
			strcpy(ndir, "check_orientation/He_27.0eV");
			Hist->fill(AUTO, "orientation_check_x", r[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -1.5, 1.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_y", r[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -1.5, 1.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_z", r[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -1.5, 1.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_x", r[0]->mom.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -0.5, 0.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_y", r[0]->mom.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -0.5, 0.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_z", r[0]->mom.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -0.5, 0.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
		}
		if (ID == 430) {
			strcpy(ndir, "check_orientation/He_30.0eV");
			Hist->fill(AUTO, "orientation_check_x", r[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -1.5, 1.0, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_y", r[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -1.5, 1.0, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_z", r[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -1.5, 1.0, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_x", r[0]->mom.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -0.5, 0.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_y", r[0]->mom.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -0.5, 0.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_z", r[0]->mom.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -0.5, 0.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
		}
		if (ID == 432) {
			strcpy(ndir, "check_orientation/He_32.0eV");
			Hist->fill(AUTO, "orientation_check_x", r[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -1.5, 1.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_y", r[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -1.5, 1.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "orientation_check_z", r[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -1.5, 1.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_x", r[0]->mom.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -0.5, 0.5, "p_rel x [a.u]", 350, -1., 1., "p_el x [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_y", r[0]->mom.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -0.5, 0.5, "p_rel y [a.u]", 350, -1., 1., "p_el y [a.u.]", ndir);
			Hist->fill(AUTO, "seccond_orientation_check_z", r[0]->mom.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -0.5, 0.5, "p_rel z [a.u]", 350, -1., 1., "p_el z [a.u.]", ndir);
		}

		// problem in energy sharing plot: there is one line with many counts		
		// special invalidates
		/*if (fabs(mol->mom_cm.x) < 7. && fabs(mol->mom_cm.y) < 7. && mol->mom_cm.z > 3.5 && mol->mom_cm.z < 16.)
		mol->valid = false;*/

		// parent ion: test for flipping
		strcpy(ndir, "angular_distr_el/Parent_r0_valid");
		if (CH == 10 && ID == 5 && e[0]->valid && r[0]->valid) {  //IN
			Hist->fill(AUTO, "recoil_xy", r[0]->raw.data.x, r[0]->raw.data.y, 1., "Recoil position", 400, -70., 70., "rx [mm]", 400, -70., 70., "ry [mm]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_px vs p_z", e[0]->mom.x, e[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -1.5, 1.5, "p_x [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			if (fabs(e[0]->mom.Mag())>0.4 && e[0]->energy()<30. && e[0]->energy()>10.){
				Hist->fill(AUTO, "check_momentum_py vs p_z_mod", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "check_momentum_px vs p_z_mod", e[0]->mom.x, e[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -1.5, 1.5, "p_x [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			}
			Hist->fill(AUTO, "check_cthetaEE_full", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 40, -1., 1., "cos(theta)", ndir);
			if (fabs(e[0]->mom.Mag())>0.4 && e[0]->energy()<30. && e[0]->energy()>10.) {
				Hist->fill(AUTO, "cos(theta)_e[0]_mod_x", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 40, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "cos(theta)_e[0]_mod_y", (e[0]->mom.y) / (e[0]->mom.Mag()), 1., "cos(theta) = py_e0/||p|| [adm]", 40, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "cos(theta)_e[0]_mod_z", (e[0]->mom.z) / (e[0]->mom.Mag()), 1., "cos(theta) = pz_e0/||p|| [adm]", 40, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "cos(theta)_e[0]_mod", e[0]->mom.Cos_Theta(), 1., "cos(theta) = pz_e0/||p|| [adm]", 40, -1., 1., "cos(theta)", ndir);
			}
		} //OUT

		  ////////////////////////////////////////////////////  CH 14  ////////////////////////////////////////////////////

		strcpy(ndir, "angular_distr_el/CH14/ID_ALL_mol_e0_valid/EN_gate");
		if (CH == 14 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(sdir,ndir);
			//Definiton from Kilian coordinate system
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			//Classical definition from Till
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mod definition from Till
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 4500., 7500., "tof1 [ns]", 350, 4500., 7000., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF_36", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 36, -1., 1., "cos(beta)", 36, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_72", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 72, -1., 1., "cos(beta)", 72, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 36, -180., 180., "phi", 36, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 36, -1., 1., "cos(beta)", 36, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_36", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 36, -1., 1., "cos(beta)", 36, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_72", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 72, -1., 1., "cos(beta)", 72, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) { // polarization plane
					Hist->fill(AUTO, "cos(theta)_e[0]_pol_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_pol_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
				if (mol->mom_rel.Angle_deg(pol_dir) < 20) {// propagation plane
					Hist->fill(AUTO, "cos(theta)_e[0]_prop_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_prop_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
				strcpy(ndir,sdir);
				strcat(ndir,"/b1_fitting_LF_redPHI");
				for (int j = 0; j < 16; j++) {
					if (fabs(j*0.125-0.875-mol_axsis.x/mol_axsis.Mag())<0.0625) {  //fabs(j*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
						sprintf_s(title,100,"cos(theta)_e[0]_cosbeta_%4.2f",j*0.125-1);
						Hist->fill(AUTO+2000+j, title, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					}
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/CDAD");
			//Hist->fill(AUTO, "cos(phi)_e[0]", (e[0]->mom.z) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 36, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "cos(phi)_e[0]_alldir", cos(atan2(e[0]->mom.y,e[0]->mom.z)), 1., "cos(phi) = cos(atan2(py,pz)) [adm]", 36, -1., 1., "cos(phi)", ndir);
			Hist->fill(AUTO, "phi_e[0]_alldir", e[0]->mom.Phi()* 180/PI, 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) {// polarization plane
				Hist->fill(AUTO, "cos(phi)_e[0]_pol",  cos(atan2(e[0]->mom.y,e[0]->mom.z)), 1., "cos(phi) = cos(atan2(py,pz)) [adm]", 36, -1., 1., "cos(phi)", ndir);
				Hist->fill(AUTO, "phi_e[0]_pol", e[0]->mom.Phi()*180/PI, 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_MF_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_multinew_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_MF_mix");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH14/ID1_mol_e0_valid/EN_gate");
		if (CH == 14 && ID == 1 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			//Definiton from Kilian coordinate system
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			//Classical definition from Till
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mod definition from Till
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 4500., 7500., "tof1 [ns]", 350, 4500., 7000., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag()) * 180 / PI , acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag()) * 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) { // polarization plane
					Hist->fill(AUTO, "cos(theta)_e[0]_pol_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_pol_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
				if (mol->mom_rel.Angle_deg(pol_dir) < 20) {// propagation plane
					Hist->fill(AUTO, "cos(theta)_e[0]_prop_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_pro_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
			}
			strcpy(sdir,ndir);
			strcat(ndir,"/MFPAD_e_MF_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_multinew_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_MF_mix");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH14/ID2_mol_e0_valid/EN_gate");
		if (CH == 14 && ID == 2 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			//Definiton from Kilian coordinate system
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			//Classical definition from Till
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mod definition from Till
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 4500., 7500., "tof1 [ns]", 350, 4500., 7000., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) { // polarization plane
					Hist->fill(AUTO, "cos(theta)_e[0]_pol_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_pol_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
				if (mol->mom_rel.Angle_deg(pol_dir) < 20) {// propagation plane
					Hist->fill(AUTO, "cos(theta)_e[0]_prop_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_prop_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
			}
			strcpy(sdir,ndir);
			strcat(ndir,"/MFPAD_e_MF_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_multinew_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_MF_mix");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH14/ID3_mol_e0_valid/EN_gate");
		if (CH == 14 && ID == 3 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			//Definiton from Kilian coordinate system
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			//Classical definition from Till
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mod definition from Till
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 4500., 7500., "tof1 [ns]", 350, 4500., 7000., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE", e[0]->mom.Cos_Theta(), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) { // polarization plane
					Hist->fill(AUTO, "cos(theta)_e[0]_pol_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_pol_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
				if (mol->mom_rel.Angle_deg(pol_dir) < 20) {// propagation plane
					Hist->fill(AUTO, "cos(theta)_e[0]_prop_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
					Hist->fill(AUTO, "theta_e[0]_prop_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				}
			}
			strcpy(sdir,ndir);
			strcat(ndir,"/MFPAD_e_MF_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_multinew_coord");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_30", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
			strcpy(ndir,sdir);
			strcat(ndir,"/MFPAD_e_MF_mix");
			Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
				Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
				if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
					Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
				}
			}
		} // from line 120!

		  ////////////////////////////////////////////////////  CH 11  ////////////////////////////////////////////////////

		strcpy(ndir, "angular_distr_el/CH11/ID_ALL_mol_e0_valid/EN_gate");
		if (CH == 11 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//mf_circ: Classical Till ref system
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mf_circ_mod: mod Till ref system
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			//definition of the third neutral particle
			CH_vector neutral_mom_ALL; double d_sum_ALL[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_ALL.x = -mol->mom_cm.x; neutral_mom_ALL.y = -mol->mom_cm.y; neutral_mom_ALL.z = -mol->mom_cm.z;
				d_sum_ALL[0]=-mol->mom_cm.x; d_sum_ALL[1]=-mol->mom_cm.y; d_sum_ALL[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_ALL);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//multinew_mf: Kilian ref system
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_ALL);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_ALL);
			//NFrame0: the same multinew_mf BUT has the modified cos(theta) and phi fucntions
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_ALL);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3800, 4600., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			strcpy(ndir,tdir);
			strcat(ndir,"/CDAD");
			//Hist->fill(AUTO, "cos(phi)_e[0]", (e[0]->mom.z) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 36, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "cos(phi)_e[0]_alldir", cos(e[0]->mom.Phi()), 1., "cos(phi) = cos(atan2(py,px)) [adm]", 36, -1., 1., "cos(phi)", ndir);
			Hist->fill(AUTO, "phi_e[0]_alldir", e[0]->mom.Phi(), 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) {// polarization plane
				Hist->fill(AUTO, "cos(phi)_e[0]_pol", cos(e[0]->mom.Phi()), 1., "cos(phi) = cos(atan2(py,px)) [adm]", 36, -1., 1., "cos(phi)", ndir);
				Hist->fill(AUTO, "phi_e[0]_pol", e[0]->mom.Phi(), 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<0.6  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.050 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.6 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() < 40) {
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_multinew_mf ref sysyem
					//USE OF A VECTOR JOINTED TO THE MF 
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_MF ref system
					//USE OF A VECTOR JOINTED TO THE MF
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_MF ref system with flipped coordinates
					//USE OF A VECTOR JOINTED TO THE MF
					strcpy(ndir,tdir);
					strcat(ndir,"/MFPAD_e_MF_mix");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
						Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
							Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					//std angles (just a try)
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_std");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2000 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2000 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2300 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2300 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2600 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2600 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//NORMAL PLANE INSTEAD OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_CH_VECTOR");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2900 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2900 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}//NORMAL PLANE INSTEAD OF LIGHT_FLIPPED ACCORDING TO KILIAN
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_CH_VECTOR_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3200 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3200 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		  //molecular break-up
		strcpy(ndir, "angular_distr_el/CH11/ID1_mol_e0_valid/EN_gate");
		if (CH == 11 && ID == 1 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH1; double d_sum_CH1[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH1.x = -mol->mom_cm.x; neutral_mom_CH1.y = -mol->mom_cm.y; neutral_mom_CH1.z = -mol->mom_cm.z;
				d_sum_CH1[0]=-mol->mom_cm.x; d_sum_CH1[1]=-mol->mom_cm.y; d_sum_CH1[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH1);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH1);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH1);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH1);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH1.y, neutral_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3800, 4600., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<0.6  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.050 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.6 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH1.y, neutral_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH1.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod() < 40) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH1dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH1.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH1dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3500 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3500 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3800 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3800 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH11/ID2_mol_e0_valid/EN_gate");
		if (CH == 11 && ID == 2 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH2; double d_sum_CH2[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH2.x = -mol->mom_cm.x; neutral_mom_CH2.y = -mol->mom_cm.y; neutral_mom_CH2.z = -mol->mom_cm.z;
				d_sum_CH2[0]=-mol->mom_cm.x; d_sum_CH2[1]=-mol->mom_cm.y; d_sum_CH2[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH2);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH2);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH2);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH2);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH2.y, neutral_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3800, 4600., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<0.6  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.050 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.6 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH2.y, neutral_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH2.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod() < 40) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH2dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH2.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH2dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4100 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4100 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4400 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4400 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH11/ID3_mol_e0_valid/EN_gate");
		if (CH == 11 && ID == 3 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH3; double d_sum_CH3[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH3.x = -mol->mom_cm.x; neutral_mom_CH3.y = -mol->mom_cm.y; neutral_mom_CH3.z = -mol->mom_cm.z;
				d_sum_CH3[0]=-mol->mom_cm.x; d_sum_CH3[1]=-mol->mom_cm.y; d_sum_CH3[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH3);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH3);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH3);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH3);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH3.y, neutral_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3800, 4600., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<0.6  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.050 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.6 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH3.y, neutral_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH3.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod() < 40) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH3dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH3.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH3dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4700 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4700 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(5000 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(5000 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!


		  ////////////////////////////////////////////////////  CH 9  ////////////////////////////////////////////////////
		strcpy(ndir, "angular_distr_el/CH9/ID_ALL_mol_e0_valid/EN_gate");
		if (CH == 9 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//mf_circ: Classical Till ref system
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//mf_circ_mod: mod Till ref system
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel,photon_dir,k);
			CH_vector e_MF_mod=mf_circ_mod.project_vector(e[0]->mom);
			//definition of the third neutral particle
			CH_vector neutral_mom_ALL; double d_sum_ALL[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_ALL.x = -mol->mom_cm.x; neutral_mom_ALL.y = -mol->mom_cm.y; neutral_mom_ALL.z = -mol->mom_cm.z;
				d_sum_ALL[0]=-mol->mom_cm.x; d_sum_ALL[1]=-mol->mom_cm.y; d_sum_ALL[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_ALL);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//multinew_mf: Kilian ref system
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_ALL);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_ALL);
			//NFrame0: the same multinew_mf BUT has the modified cos(theta) and phi fucntions
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_ALL);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 2700, 3300., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			strcpy(ndir,tdir);
			strcat(ndir,"/CDAD");
			//Hist->fill(AUTO, "cos(phi)_e[0]", (e[0]->mom.z) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 36, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "cos(phi)_e[0]_alldir", cos(e[0]->mom.Phi()), 1., "cos(phi) = cos(atan2(py,px)) [adm]", 36, -1., 1., "cos(phi)", ndir);
			Hist->fill(AUTO, "phi_e[0]_alldir", e[0]->mom.Phi(), 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 20) {// polarization plane
				Hist->fill(AUTO, "cos(phi)_e[0]_pol", cos(e[0]->mom.Phi()), 1., "cos(phi) = cos(atan2(py,px)) [adm]", 36, -1., 1., "cos(phi)", ndir);
				Hist->fill(AUTO, "phi_e[0]_pol", e[0]->mom.Phi(), 1., "phi = atan2(py,px) [deg]", 36, -180., 180., "phi", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<1.00  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.10 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.00  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.7 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() < -5) {
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_multinew_mf ref sysyem
					//USE OF A VECTOR JOINTED TO THE MF 
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_MF ref system
					//USE OF A VECTOR JOINTED TO THE MF
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_MF ref system with flipped coordinates
					//USE OF A VECTOR JOINTED TO THE MF
					strcpy(ndir,tdir);
					strcat(ndir,"/MFPAD_e_MF_mix");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 15) {
						Hist->fill(AUTO, "MFPAD3D_polplane_usr_15", e_MF.Phi_deg_mod(), e_MF.Cos_Theta_mod(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 15) {
							Hist->fill(AUTO, "MFPAD_polplane_usr_15", e_MF.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					//std angles (just a try)
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_std");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2000 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2000 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2300 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2300 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2600 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2600 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//NORMAL PLANE INSTEAD OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_CH_VECTOR");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(2900 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(2900 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}//NORMAL PLANE INSTEAD OF LIGHT_FLIPPED ACCORDING TO KILIAN
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_CH_VECTOR_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3200 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3200 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH9/ID1_mol_e0_valid/EN_gate");
		if (CH == 9 && ID == 1 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH1; double d_sum_CH1[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH1.x = -mol->mom_cm.x; neutral_mom_CH1.y = -mol->mom_cm.y; neutral_mom_CH1.z = -mol->mom_cm.z;
				d_sum_CH1[0]=-mol->mom_cm.x; d_sum_CH1[1]=-mol->mom_cm.y; d_sum_CH1[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH1);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH1);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH1);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH1);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH1.y, neutral_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (CH1)", 350, 2700, 3300., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<1.00  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.10 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.00  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.7 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH1.y, neutral_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH1.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod() < -5) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH1dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH1.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH1.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH1dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3500 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3500 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(3800 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(3800 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH9/ID2_mol_e0_valid/EN_gate");
		if (CH == 9 && ID == 2 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH2; double d_sum_CH2[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH2.x = -mol->mom_cm.x; neutral_mom_CH2.y = -mol->mom_cm.y; neutral_mom_CH2.z = -mol->mom_cm.z;
				d_sum_CH2[0]=-mol->mom_cm.x; d_sum_CH2[1]=-mol->mom_cm.y; d_sum_CH2[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH2);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH2);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH2);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH2);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH2.y, neutral_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (CH2)", 350, 2700, 3300., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<1.00  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.10 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.00  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.7 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH2.y, neutral_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH2.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod() < -5) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH2dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH2.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH2.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH2dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4100 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4100 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4400 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4400 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		strcpy(ndir, "angular_distr_el/CH9/ID3_mol_e0_valid/EN_gate");
		if (CH == 9 && ID == 3 && mol->valid && e[0]->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el)	 {
			strcpy(tdir, ndir);
			//Classical Till ref
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel,photon_dir);
			CH_vector e_MF=mf_circ.project_vector(e[0]->mom);
			//Third neutral particle
			CH_vector neutral_mom_CH3; double d_sum_CH3[3]; double Nummertosort=-2;
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_mom_CH3.x = -mol->mom_cm.x; neutral_mom_CH3.y = -mol->mom_cm.y; neutral_mom_CH3.z = -mol->mom_cm.z;
				d_sum_CH3[0]=-mol->mom_cm.x; d_sum_CH3[1]=-mol->mom_cm.y; d_sum_CH3[2]=-mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_CH3);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_CH3);
			Coordinate_System multinew_mf = Coordinate_System(vcross,mol->ion[1]->mom,k);
			CH_vector e_multinew_mf=multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf=multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_CH3);
			//Newton plots: Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
			//NOTE: mol->ion[0]->mom instead of cross01
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom,k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_CH3);
			Hist->fill(AUTO,"Newton plot",NFrame0_mom[1].x/NFrame0_mom[0].Mag(),-NFrame0_mom[1].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",NFrame0_mom[2].x/NFrame0_mom[0].Mag(),-NFrame0_mom[2].y/NFrame0_mom[0].Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_CH3.y, neutral_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (CH3)", 350, 2700, 3300., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag())* 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag())* 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y,e[0]->mom.z)*180/PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir,tdir);
			if (NFrame0_mom[1].x/NFrame0_mom[0].Mag()>-1.8  &&  NFrame0_mom[1].x/NFrame0_mom[0].Mag()<-0.7 &&
				-NFrame0_mom[1].y/NFrame0_mom[0].Mag()<1.00  && -NFrame0_mom[1].y/NFrame0_mom[0].Mag()>0.10 &&
				NFrame0_mom[2].x/NFrame0_mom[0].Mag()>-0.3  &&  NFrame0_mom[2].x/NFrame0_mom[0].Mag()<0.5  && 
				-NFrame0_mom[2].y/NFrame0_mom[0].Mag()<0.00  && -NFrame0_mom[2].y/NFrame0_mom[0].Mag()>-0.7 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_CH3.y, neutral_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_MF_gated",NFrame0_mom[1].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_MF_gated",NFrame0_mom[2].x/mol->ion[0]->mom.Mag(),-NFrame0_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH3.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);
#if MFPAD == 1
				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod() < -5) {
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_CH3dir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_CH3.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_CH3.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					strcpy(ndir, tdir);
					strcat(ndir,"/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_CH3dir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					//SINGLE MFPAD WITH MOLECULAR ORIENTATION
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(4700 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(4700 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
					//ACCORDING TO KILIAN, HERE +COS(Theta) and MIRRORED_PHI FOR THE CHANGE OF ellicity OF LIGHT
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_multinew_flipped");
					for (int i = 0; i < 12; i++) {
						for (int j = 0; j < 6; j++) { //fabs(i*bin-(period/2-bin/2)-x)<bin/2  ; period=bin*i
							if (fabs(i*30-165.+light_multinew_mf.Phi_deg_mod_mirr())<15 && abs(j*0.3333333-0.8333335+light_multinew_mf.Cos_Theta_mod())<0.166666) {
								sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								sprintf_s(title_el, 100, "cos(theta)_e[0]_costheta_%4.2f_phi_%d", j*0.333 - 1, i * 30 - 180);
								//cout << "i value = " << i << ", j value = " << j << ".\n\n";
								Nummertosort=i*6+j;
								Hist->fill3(AUTO,"MFPAD_Mathematica",e_multinew_mf.Phi_deg_mod(),e_multinew_mf.Cos_Theta_mod(),Nummertosort,1.,"phiLaser",nbin_phi,-180.,180.,"PhiELE",nbin_cos,-1.,1.,"ThetELE", 72,-0.5,71-5,"Nummertosort", ndir);
								Hist->fill(5000 + i * 10 + j, title, e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", nbin_phi, -180., 180., "Phi", nbin_cos, -1., 1., "cos(theta)", ndir);
								Hist->fill(5000 + i * 10 + j + 130, title_el, e_multinew_mf.Cos_Theta_mod(), 1., "cos(theta) = px_e0/||p|| [adm]", nbin_cos+nbin_cos/2, -1., 1., "cos(theta)", ndir);
							}
						}
					}
				}
#endif
			}
		} // from line 120!

		  //test section
#if TEST == 1
		if (mol->valid && e[0]->energy() > emin_el && e[0]->energy() < emax_el) {
			strcpy(ndir, "angular_distr_el/ALL_CH/test");
			strcpy(tdir, ndir);
			//mf_circ: Classical Till ref system
			Coordinate_System mf_circ = Coordinate_System(mol->mom_rel, photon_dir);
			CH_vector e_MF = mf_circ.project_vector(e[0]->mom);
			//mf_circ_mod: mod Till ref system
			Coordinate_System mf_circ_mod = Coordinate_System(mol->mom_rel, photon_dir, k);
			CH_vector e_MF_mod = mf_circ_mod.project_vector(e[0]->mom);
			//definition of the third neutral particle
			CH_vector neutral_mom_ALL; double d_sum_ALL[3]; double Nummertosort = -2;
			if (mol->mom_cm.Mag() > sum_min && mol->mom_cm.Mag() < sum_max) {
				neutral_mom_ALL.x = -mol->mom_cm.x; neutral_mom_ALL.y = -mol->mom_cm.y; neutral_mom_ALL.z = -mol->mom_cm.z;
				d_sum_ALL[0] = -mol->mom_cm.x; d_sum_ALL[1] = -mol->mom_cm.y; d_sum_ALL[2] = -mol->mom_cm.z;
			}
			CH_vector mom_relsum = mol->mom_rel.operator-(neutral_mom_ALL);
			CH_vector normal_plane = mol->mom_rel.Cross(mom_relsum);
			CH_vector vcross = mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//multinew_mf: Kilian ref system
			//Coordinate_System multinew_mf = Coordinate_System(vcross,neutral_mom_ALL);
			Coordinate_System multinew_mf = Coordinate_System(vcross, mol->ion[1]->mom, k);
			CH_vector e_multinew_mf = multinew_mf.project_vector(e[0]->mom);
			CH_vector light_multinew_mf = multinew_mf.project_vector(photon_dir);
			CH_vector mf_mom[3];
			mf_mom[0] = multinew_mf.project_vector(mol->ion[0]->mom);
			mf_mom[1] = multinew_mf.project_vector(mol->ion[1]->mom);
			mf_mom[2] = multinew_mf.project_vector(neutral_mom_ALL);
			//NFrame0: the same multinew_mf BUT has the modified cos(theta) and phi fucntions
			Coordinate_System NFrame0 = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom, k); //NOTE: mol->ion[0]->mom instead of CROSS
			CH_vector NFrame0_mom[3];
			NFrame0_mom[0] = NFrame0.project_vector(mol->ion[0]->mom);
			NFrame0_mom[1] = NFrame0.project_vector(mol->ion[1]->mom);
			NFrame0_mom[2] = NFrame0.project_vector(neutral_mom_ALL);
			Hist->fill(AUTO, "Newton plot", NFrame0_mom[1].x / NFrame0_mom[0].Mag(), -NFrame0_mom[1].y / NFrame0_mom[0].Mag(), 1., "Newton plot", 80, -2.5, 2.5, "px", 120, -3.0, 3.0, "py", ndir);
			Hist->fill(AUTO - 1, "Newton plot", NFrame0_mom[2].x / NFrame0_mom[0].Mag(), -NFrame0_mom[2].y / NFrame0_mom[0].Mag(), 1., "Newton plot", 80, -2.5, 2.5, "px", 120, -3.0, 3.0, "py", ndir);
			Hist->fill(AUTO, "Normal_plane_phi", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "Normal_plane_costheta", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "light_mf.Phi", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
			Hist->fill(AUTO, "light_mf.theta", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);

			//custom
			Hist->fill(AUTO, "check_neutral_momentum_py vs p_z", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "check_momentum_py vs p_z", e[0]->mom.y, e[0]->mom.z, 1., "p_y vs. p_z (no pol)", 100, -1.5, 1.5, "p_y [a.u.]", 100, -1.5, 1.5, "p_z [a.u.]", ndir);
			Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 2700, 3300., "tof1 [ns]", 350, 6300., 6900., "tof2 [ns]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_x", (e[0]->mom.x / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons x dir", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "check_cthetaEE_z", (e[0]->mom.z / e[0]->mom.Mag()), e[0]->energy(), 1., "theta angular distribution electrons", 20, -1., 1., "cos(theta)", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			Hist->fill(AUTO, "cos(theta)_e[0]", (e[0]->mom.x / e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
			Hist->fill(AUTO, "PECD_LF", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_36", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 36, -1., 1., "cos(beta)", 36, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_72", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 72, -1., 1., "cos(beta)", 72, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_ang", acos(mol_axsis.x / mol_axsis.Mag()) * 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
			Hist->fill(AUTO, "PECD_LF_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
			Hist->fill(AUTO, "PECD_LF_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			if (phi_yz > -140 && phi_yz < -40 || phi_yz > 40 && phi_yz < 140) {
				Hist->fill(AUTO, "cos(theta)_e[0]_redPHI", (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "cos(theta) = px_e0/||p|| [adm]", 15, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "theta_e[0]_redPHI", acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "theta = acos(px_e0/||p||) [adm]", 15, 0., 180., "cos(theta)", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_ang", acos(mol_axsis.x / mol_axsis.Mag()) * 180 / PI, acos(e[0]->mom.x / e[0]->mom.Mag()) * 180 / PI, 1., "PECD map theta e[0] vs mol orientation", 24, 0., 180., "beta", 24, 0., 180., "theta", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_phi", atan2(e[0]->mom.y, e[0]->mom.z) * 180 / PI, (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -180., 180., "phi", 16, -1., 1., "cos(theta)e[0]", ndir);
				Hist->fill(AUTO, "PECD_LF_redPHI_en", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), e[0]->energy(), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
			}
			
			//CUSTOM GATE on the Newtonplot
			// x and y are different becasue of the coordinate_system definition (according to Kilian one)
			strcpy(ndir, tdir);
			if (NFrame0_mom[1].x / NFrame0_mom[0].Mag() > -1.4  &&  NFrame0_mom[1].x / NFrame0_mom[0].Mag() < -0.6 &&
				-NFrame0_mom[1].y / NFrame0_mom[0].Mag() < 0.5  && -NFrame0_mom[1].y / NFrame0_mom[0].Mag() > -0.1 &&
				NFrame0_mom[2].x / NFrame0_mom[0].Mag() > -0.4  &&  NFrame0_mom[2].x / NFrame0_mom[0].Mag() < 0.6  &&
				-NFrame0_mom[2].y / NFrame0_mom[0].Mag() < 0.1  && -NFrame0_mom[2].y / NFrame0_mom[0].Mag() > -0.5 &&
				mol->mom_cm.Mag() < sum_max) {
				strcat(ndir, "/mol_orient_GATED");
				//custom
				Hist->fill(AUTO, "check_neutral_momentum_py vs p_z_gated", neutral_mom_ALL.y, neutral_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "Newton plot_MF_gated", NFrame0_mom[1].x / mol->ion[0]->mom.Mag(), -NFrame0_mom[1].y / mol->ion[0]->mom.Mag(), 1., "Newton plot", 80, -2.5, 2.5, "px", 120, -3.0, 3.0, "py", ndir);
				Hist->fill(AUTO - 1, "Newton plot_MF_gated", NFrame0_mom[2].x / mol->ion[0]->mom.Mag(), -NFrame0_mom[2].y / mol->ion[0]->mom.Mag(), 1., "Newton plot", 80, -2.5, 2.5, "px", 120, -3.0, 3.0, "py", ndir);
				Hist->fill(AUTO, "Normal_plane_phi_gated", normal_plane.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "Normal_plane_costheta_gated", normal_plane.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				Hist->fill(AUTO, "light_mf.Phi_gated", light_multinew_mf.Phi_deg_mod(), 1., "test", 100, -180., 180., "Phi", ndir);
				Hist->fill(AUTO, "light_mf.theta_gated", light_multinew_mf.Cos_Theta_mod(), 1., "test", 100, -1., 1., "cos(theta)", ndir);
				//Phi_deg_mod: NOTE the Newtonplot PLANE changes upond the change of coordinates: Phi is the angle btw -y and z, cos(theta) z/norm
				Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_mod", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_mod", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_mod", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
				Hist->fill(AUTO, "check_ion[0]_xz", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "x mom", 200, -200., 200., "z mom", ndir);
				Hist->fill(AUTO, "check_ion[1]_yz", mol->ion[1]->mom.y, mol->ion[1]->mom.z, 1., "pos. in lab frame", 200, -200., 200., "y mom", 200, -200., 200., "z mom", ndir);

				if (multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() > -125 && multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod() < -95 &&
					multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() > -80 && multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod() < -5) {
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_multinew_mf ref sysyem
					//USE OF A VECTOR JOINTED TO THE MF 
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPADs_e_multinew_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), e_multinew_mf.Cos_Theta_mod(), 1., "molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e_multinew_mf.Angle_deg(light_multinew_mf)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_multinew_mf.Phi_deg_mod(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					Hist->fill(AUTO, "check_pos_ion[0].multinew_mf_gated", multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[0]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[1].multinew_mf_gated", multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Phi_deg_mod(), multinew_mf.project_vector(mol->ion[1]->mom.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					Hist->fill(AUTO, "check_pos_ion[2].multinew_mf_gated", multinew_mf.project_vector(neutral_mom_ALL.Norm()).Phi_deg_mod(), multinew_mf.project_vector(neutral_mom_ALL.Norm()).Cos_Theta_mod(), 1., "pos. in lab frame", 36, -180., 180, "phi", 18, -1., 1, "cos(theta)", ndir);
					//INTEGRATED MFPAD ALONG ALL MOLECULAR DIRECTIONS e_MF ref system
					//USE OF A VECTOR JOINTED TO THE MF
					strcpy(ndir, tdir);
					strcat(ndir, "/MFPAD_e_MF_check");
					Hist->fill(AUTO, "MFPAD3D_ALLdir_usr", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
					if (fabs(mol->mom_rel.Angle_deg(photon_dir) - 90.0) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_molrel_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(e[0]->mom.Angle_deg(photon_dir) - 90.0) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_molrel_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
					if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
						Hist->fill(AUTO, "MFPAD3D_polplane_vnormal_30_gated", e_MF.Phi_deg(), e_MF.Cos_Theta(), 1., "molecular frame angular distribution", 48, -180., 180., "Phi", 36, -1.0, 1.0, "cos(theta)", ndir);
						if (fabs(normal_plane.Angle_deg(photon_dir)) < 30) {
							Hist->fill(AUTO, "MFPAD_polplane_vnormal_30_gated", e_MF.Phi_deg(), 1., "molecular frame angular distribution", 72, -180., 180., "Phi", ndir);
						}
					}
				}
			}
		}

		if (CH == 14 && ID == 1 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH14/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/ID=1");
			if (mol->valid && rhit > 1 && mol->ion[0]->raw.data.tof < mol->ion[1]->raw.data.tof) {
				phi_p = mol->ion[0]->raw.phi;
				phi_m = mol->ion[1]->raw.phi;
				relphi = phi_p - phi_m;
				if (phi_p > phi_m) {
					if (relphi > 180.) {
						relphi = relphi - 360.;
					}
				}
				else {
					if (phi_p <= phi_m) {
						if (relphi < -180.) {
							relphi = relphi + 360.;
						}
					}
				}
				Hist->fill(AUTO, "tof_ion1 vs. event_number", (int)(evt->event_number), mol->ion[0]->raw.data.tof, 1., "tof_ion1 vs #events", 200, 0., 6e+5, "#events [adm]", 200, 4900., 5500., "tof [ns]", ndir);
				Hist->fill(AUTO, "tof_ion2 vs. event_number", (int)(evt->event_number), mol->ion[1]->raw.data.tof, 1., "tof_ion2 vs #events",200, 0., 6e+5, "#events [adm]", 200, 6300, 7000., "tof [ns]", ndir);
				Hist->fill(AUTO, "tof_el vs. event_number", e[0]->raw.data.tof, (int)(evt->event_number), 1., "tof_el vs #events", 50, 34., 44., "tof [ns]", 200, 0., 6e+5, "#events [adm]", ndir);

				Hist->fill(AUTO, "orientation_check_x", mol->mom_cm.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -200., 200., "p_rel x [a.u]", 350, -1.5, 1.5, "p_el x [a.u.]", ndir);
				Hist->fill(AUTO, "orientation_check_y", mol->mom_cm.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -200., 200., "p_rel y [a.u]", 350, -1.5, 1.5, "p_el y [a.u.]", ndir);
				Hist->fill(AUTO, "orientation_check_z", mol->mom_cm.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -200., 200., "p_rel z [a.u]", 350, -1.5, 1.5, "p_el z [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_x", mol->mom_cm.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -200., 200., "p_rel x [a.u]", 350, -1.5, 1.5, "p_el x [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_y", mol->mom_cm.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -200., 200., "p_rel y [a.u]", 350, -1.5, 1.5, "p_el y [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_z", mol->mom_cm.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -200., 200., "p_rel z [a.u]", 350, -1.5, 1.5, "p_el z [a.u.]", ndir);

				// Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3000., 13000., "tof1 [ns]", 350, 3000., 13000., "tof2 [ns]", ndir);
				if ((fabs(relphi - 180.0) < 10.0)) {
					Hist->fill(AUTO, "pipico_momcon_short", evt->r.tof[0], evt->r.tof[1], 1., "PIPICO spectrum (fine) /w simple mom. check", 350, 4500, 7800., "rec1 TOF [ns]", 350, 4500., 7800., "rec2 TOF [ns]", ndir);
					Hist->fill(AUTO, "PIPICO_phi", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (phi condition)", 350, 4500., 7800., "tof1 [ns]", 350, 4500., 7800., "tof2 [ns]", ndir);
					Hist->fill(AUTO, "sum_diff_phi", mol->ion[0]->raw.data.tof + mol->ion[1]->raw.data.tof, mol->ion[1]->raw.data.tof - mol->ion[0]->raw.data.tof, 1., "SUM-DIFF (phi condition)", 350, 11650, 12000., "tof 1+2 [ns]", 350, 1200., 1600., "tof 2-1 [ns]", ndir);
					//if (e[0]->energy() < 10.) {
					//	Hist->fill(AUTO, "PIPICO_phi_e10", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (phi condition)", 350, 3000., 4800., "tof1 [ns]", 350, 4000., 5500., "tof2 [ns]", ndir);
					//	Hist->fill(AUTO, "sum_diff_phi_e10", mol->ion[0]->raw.data.tof + mol->ion[1]->raw.data.tof, mol->ion[1]->raw.data.tof - mol->ion[0]->raw.data.tof, 1., "SUM-DIFF (phi condition)", 350, 11650, 12000., "tof 1+2 [ns]", 350, 1200., 1600., "tof 2-1 [ns]", ndir);
					//}
				}
				Hist->fill(AUTO, "KER vs electron sum energy", mol->KER(), e[0]->energy(), 1, "kinetic energy release vs |p_sum_el|", 100, 0., 15., "KER [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "ion1 energy vs electron 1 energy", mol->ion[0]->energy(), e[0]->energy(), 1., "energy ion 1 vs electron 1", 100, 0., 10., "ion 1 en [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "ion2 energy vs electron 1 energy", mol->ion[1]->energy(), e[0]->energy(), 1., "energy ion 2 vs electron 1", 100, 0., 10., "ion 2 en [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "mol_ang_distrib", mol_axsis.x / mol_axsis.Mag(), 1., "angle respect to photon dir [deg]", 32, -1., 1., "beta", ndir);

				if (e[0]->raw.data.y> 0.) {
					Hist->fill(AUTO, "PECD_LF_y_upDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.y< 0.) {
					Hist->fill(AUTO, "PECD_LF_y_downDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.x> 0.) {
					Hist->fill(AUTO, "PECD_LF_x_upDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.x< 0.) {
					Hist->fill(AUTO, "PECD_LF_x_downDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
			}
		}

		if (CH == 14 && ID == 2 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH14/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/ID=2");
			if (mol->valid && rhit > 1 && mol->ion[0]->raw.data.tof < mol->ion[1]->raw.data.tof) {
				phi_p = mol->ion[0]->raw.phi;
				phi_m = mol->ion[1]->raw.phi;
				relphi = phi_p - phi_m;
				if (phi_p > phi_m) {
					if (relphi > 180.) {
						relphi = relphi - 360.;
					}
				}
				else {
					if (phi_p <= phi_m) {
						if (relphi < -180.) {
							relphi = relphi + 360.;
						}
					}
				}
				Hist->fill(AUTO, "tof_ion1 vs. event_number", (int)(evt->event_number), mol->ion[0]->raw.data.tof, 1., "tof_ion1 vs #events", 200, 0., 6e+5, "#events [adm]", 200, 4900., 5500., "tof [ns]", ndir);
				Hist->fill(AUTO, "tof_ion2 vs. event_number", (int)(evt->event_number), mol->ion[1]->raw.data.tof, 1., "tof_ion2 vs #events",200, 0., 6e+5, "#events [adm]", 200, 6300, 7000., "tof [ns]", ndir);
				Hist->fill(AUTO, "tof_el vs. event_number", e[0]->raw.data.tof, (int)(evt->event_number), 1., "tof_el vs #events", 50, 34., 44., "tof [ns]", 200, 0., 6e+5, "#events [adm]", ndir);

				Hist->fill(AUTO, "orientation_check_x", mol->mom_cm.x, e[0]->mom.x, 1., "p_rel x vs pel x", 350, -200., 200., "p_rel x [a.u]", 350, -1.5, 1.5, "p_el x [a.u.]", ndir);
				Hist->fill(AUTO, "orientation_check_y", mol->mom_cm.y, e[0]->mom.y, 1., "p_rel y vs pel y", 350, -200., 200., "p_rel y [a.u]", 350, -1.5, 1.5, "p_el y [a.u.]", ndir);
				Hist->fill(AUTO, "orientation_check_z", mol->mom_cm.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -200., 200., "p_rel z [a.u]", 350, -1.5, 1.5, "p_el z [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_x", mol->mom_cm.x + e[0]->mom.x, e[0]->mom.x, 1., "p_rel x vs pel z", 350, -200., 200., "p_rel x [a.u]", 350, -1.5, 1.5, "p_el x [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_y", mol->mom_cm.y + e[0]->mom.y, e[0]->mom.y, 1., "p_rel y vs pel z", 350, -200., 200., "p_rel y [a.u]", 350, -1.5, 1.5, "p_el y [a.u.]", ndir);
				Hist->fill(AUTO, "seccond_orientation_check_z", mol->mom_cm.z + e[0]->mom.z, e[0]->mom.z, 1., "p_rel z vs pel z", 350, -200., 200., "p_rel z [a.u]", 350, -1.5, 1.5, "p_el z [a.u.]", ndir);

				// Hist->fill(AUTO, "PIPICO", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (all)", 350, 3000., 13000., "tof1 [ns]", 350, 3000., 13000., "tof2 [ns]", ndir);
				if ((fabs(relphi - 180.0) < 10.0)) {
					Hist->fill(AUTO, "pipico_momcon_short", evt->r.tof[0], evt->r.tof[1], 1., "PIPICO spectrum (fine) /w simple mom. check", 350, 4500, 7800., "rec1 TOF [ns]", 350, 4500., 7800., "rec2 TOF [ns]", ndir);
					Hist->fill(AUTO, "PIPICO_phi", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (phi condition)", 350, 4500., 7800., "tof1 [ns]", 350, 4500., 7800., "tof2 [ns]", ndir);
					Hist->fill(AUTO, "sum_diff_phi", mol->ion[0]->raw.data.tof + mol->ion[1]->raw.data.tof, mol->ion[1]->raw.data.tof - mol->ion[0]->raw.data.tof, 1., "SUM-DIFF (phi condition)", 350, 11650, 12000., "tof 1+2 [ns]", 350, 1200., 1600., "tof 2-1 [ns]", ndir);
					//if (e[0]->energy() < 10.) {
					//	Hist->fill(AUTO, "PIPICO_phi_e10", mol->ion[0]->raw.data.tof, mol->ion[1]->raw.data.tof, 1., "PIPICO (phi condition)", 350, 3000., 4800., "tof1 [ns]", 350, 4000., 5500., "tof2 [ns]", ndir);
					//	Hist->fill(AUTO, "sum_diff_phi_e10", mol->ion[0]->raw.data.tof + mol->ion[1]->raw.data.tof, mol->ion[1]->raw.data.tof - mol->ion[0]->raw.data.tof, 1., "SUM-DIFF (phi condition)", 350, 11650, 12000., "tof 1+2 [ns]", 350, 1200., 1600., "tof 2-1 [ns]", ndir);
					//}
				}
				Hist->fill(AUTO, "KER vs electron sum energy", mol->KER(), e[0]->energy(), 1, "kinetic energy release vs |p_sum_el|", 100, 0., 15., "KER [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "ion1 energy vs electron 1 energy", mol->ion[0]->energy(), e[0]->energy(), 1., "energy ion 1 vs electron 1", 100, 0., 10., "ion 1 en [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "ion2 energy vs electron 1 energy", mol->ion[1]->energy(), e[0]->energy(), 1., "energy ion 2 vs electron 1", 100, 0., 10., "ion 2 en [eV]", 100, emin_el, emax_el, "electron energy [eV]", ndir);
				Hist->fill(AUTO, "mol_ang_distrib", mol_axsis.x / mol_axsis.Mag(), 1., "angle respect to photon dir [deg]", 32, -1., 1., "beta", ndir);

				if (e[0]->raw.data.y> 0.) {
					Hist->fill(AUTO, "PECD_LF_y_upDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.y< 0.) {
					Hist->fill(AUTO, "PECD_LF_y_downDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.x> 0.) {
					Hist->fill(AUTO, "PECD_LF_x_upDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
				if (e[0]->raw.data.x< 0.) {
					Hist->fill(AUTO, "PECD_LF_x_downDETECT", mol_axsis.x / mol_axsis.Mag(), (e[0]->mom.x) / (e[0]->mom.Mag()), 1., "PECD map cos(theta)e[0] vs mol orientation", 16, -1., 1., "cos(beta)", 16, -1., 1., "cos(theta)e[0]", ndir);
				}
			}
		}

		if (CH == 9 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH9/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/ALL+1");
			CH_vector neutral_test_mom_ALL;

			double d_sum_test_ALL[3];
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_test_mom_ALL.x = -mol->mom_cm.x; neutral_test_mom_ALL.y = -mol->mom_cm.y; neutral_test_mom_ALL.z = -mol->mom_cm.z;
				d_sum_test_ALL[0]=-mol->mom_cm.x; d_sum_test_ALL[1]=-mol->mom_cm.y; d_sum_test_ALL[2]=-mol->mom_cm.z;
			}
			CH_vector cross01; cross01=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			double d_cross01[3]; d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;
			labframe_transformation(1, d_cross01, d_frag1, d_photon_dir);  d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_sum_test_ALL);d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_frag2);       d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_el);          d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;

			double phiphot  = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2]);
			double phifrag1 = 180/PI*atan2(d_frag1[1],d_frag1[2]);
			double phifrag2 = 180/PI*atan2(d_frag2[1],d_frag2[2]);
			double phifrag3 = 180/PI*atan2(d_sum_test_ALL[1],d_sum_test_ALL[2]);
			double phiel =    180/PI*atan2(d_el[1],d_el[2]);
			double mirrorphiphot = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2])+180;
			if(mirrorphiphot>180){mirrorphiphot=mirrorphiphot-360;}

			double thetaphot = -(d_photon_dir[0]/(sqrt(d_photon_dir[0]*d_photon_dir[0]+d_photon_dir[1]*d_photon_dir[1]+d_photon_dir[2]*d_photon_dir[2])));
			double thetafrag1= -(d_frag1[0]/(sqrt(d_frag1[0]*d_frag1[0]+d_frag1[1]*d_frag1[1]+d_frag1[2]*d_frag1[2])));
			double thetafrag2= -(d_frag2[0]/(sqrt(d_frag2[0]*d_frag2[0]+d_frag2[1]*d_frag2[1]+d_frag2[2]*d_frag2[2])));
			double thetafrag3= -(d_sum_test_ALL[0]/(sqrt(d_sum_test_ALL[0]*d_sum_test_ALL[0]+d_sum_test_ALL[1]*d_sum_test_ALL[1]+d_sum_test_ALL[2]*d_sum_test_ALL[2])));
			double thetael =   -(d_el[0]/(sqrt(d_el[0]*d_el[0]+d_el[1]*d_el[1]+d_el[2]*d_el[2])));

			Hist->fill(AUTO,"MFAPD P_all",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			Hist->fill(AUTO,"MFAPD P_all_-theta",phiel ,-thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			strcpy(ndir, tdir);
			strcat(ndir, "/ALL+1/MFPADs");

			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333 - 0.833333 - thetaphot) < 0.1666666 && fabs(j * 30 - 165 - phiphot) < 15) {
						sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", k*0.333 - 1, j * 30 - 180);
						Hist->fill(5300+j*10+k, title, phiel, thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}
			strcpy(ndir, tdir);
			strcat(ndir, "/ALL+1/MFPADs_flipped");
			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333-0.833333+thetaphot)<0.1666666 && fabs(j*30-165-mirrorphiphot)<15) {
						sprintf_s(title, 100, "MFPAD3D_engate_-costheta_%4.2f_mirrphi_%d",k*0.333-1, j*30-180);
						Hist->fill(5800+j*10+k,title,phiel,thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}

			strcpy(ndir, tdir);
			strcat(ndir, "/ALL+1");
			//Newton plots
			CH_vector cross02; cross02=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_ALL);
			Coordinate_System multitest_mf = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom);
			CH_vector e_multitest_mf=multitest_mf.project_vector(e[0]->mom);
			CH_vector light_multitest_mf=multitest_mf.project_vector(photon_dir);
			CH_vector multitest_mf_mom[3];
			multitest_mf_mom[0] = multitest_mf.project_vector(mol->ion[0]->mom);
			multitest_mf_mom[1] = multitest_mf.project_vector(mol->ion[1]->mom);
			multitest_mf_mom[2] = multitest_mf.project_vector(neutral_test_mom_ALL);

			Hist->fill(AUTO,"Newton plot",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);

			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-1.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.7 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00 &&
				multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()>-0.3 &&  multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()<0.5  && 
				-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()<0.00 && -multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()>-0.7 &&
				mol->mom_cm.Mag() < sum_max) {
				Hist->fill(AUTO,"MFAPD P_all_gated_1",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_1",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_1",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_1", neutral_test_mom_ALL.y, neutral_test_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-0.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.5 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00) {
				Hist->fill(AUTO,"MFAPD P_all_gated_2",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_2",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_2",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_2", neutral_test_mom_ALL.y, neutral_test_mom_ALL.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
		}

		if (CH == 9 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH9/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/ALL-1");
			CH_vector neutral_test_mom_ALL;

			double d_sum_test_ALL[3];
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_test_mom_ALL.x = -mol->mom_cm.x; neutral_test_mom_ALL.y = -mol->mom_cm.y; neutral_test_mom_ALL.z = -mol->mom_cm.z;
				d_sum_test_ALL[0]=-mol->mom_cm.x; d_sum_test_ALL[1]=-mol->mom_cm.y; d_sum_test_ALL[2]=-mol->mom_cm.z;
			}

			CH_vector cross01; cross01=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_ALL);
			Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_ALL);
			CH_vector e_multitest_mf=multitest_mf.project_vector(e[0]->mom);
			CH_vector light_multitest_mf=multitest_mf.project_vector(photon_dir);
			double d_cross01[3]; d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;
			labframe_transformation(-1, d_cross01, d_frag1, d_photon_dir);  d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(-1, d_cross01, d_frag1, d_sum_test_ALL);d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(-1, d_cross01, d_frag1, d_frag2);       d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(-1, d_cross01, d_frag1, d_el);          d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;

			double phiphot  = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2]);
			double phifrag1 = 180/PI*atan2(d_frag1[1],d_frag1[2]);
			double phifrag2 = 180/PI*atan2(d_frag2[1],d_frag2[2]);
			double phifrag3 = 180/PI*atan2(d_sum_test_ALL[1],d_sum_test_ALL[2]);
			double phiel =    180/PI*atan2(d_el[1],d_el[2]);
			double mirrorphiphot = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2])+180;
			if(mirrorphiphot>180){mirrorphiphot=mirrorphiphot-360;}

			double thetaphot = -(d_photon_dir[0]/(sqrt(d_photon_dir[0]*d_photon_dir[0]+d_photon_dir[1]*d_photon_dir[1]+d_photon_dir[2]*d_photon_dir[2])));
			double thetafrag1= -(d_frag1[0]/(sqrt(d_frag1[0]*d_frag1[0]+d_frag1[1]*d_frag1[1]+d_frag1[2]*d_frag1[2])));
			double thetafrag2= -(d_frag2[0]/(sqrt(d_frag2[0]*d_frag2[0]+d_frag2[1]*d_frag2[1]+d_frag2[2]*d_frag2[2])));
			double thetafrag3= -(d_sum_test_ALL[0]/(sqrt(d_sum_test_ALL[0]*d_sum_test_ALL[0]+d_sum_test_ALL[1]*d_sum_test_ALL[1]+d_sum_test_ALL[2]*d_sum_test_ALL[2])));
			double thetael =   -(d_el[0]/(sqrt(d_el[0]*d_el[0]+d_el[1]*d_el[1]+d_el[2]*d_el[2])));

			Hist->fill(AUTO,"MFAPD P_all",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			Hist->fill(AUTO,"MFAPD P_all_-theta",phiel ,-thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
		}

		if (CH == 9 && ID == 1 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH9/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/CH1+1");
			CH_vector neutral_test_mom_CH1;

			double d_sum_test_CH1[3];
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_test_mom_CH1.x = -mol->mom_cm.x; neutral_test_mom_CH1.y = -mol->mom_cm.y; neutral_test_mom_CH1.z = -mol->mom_cm.z;
				d_sum_test_CH1[0]=-mol->mom_cm.x; d_sum_test_CH1[1]=-mol->mom_cm.y; d_sum_test_CH1[2]=-mol->mom_cm.z;
			}
			CH_vector cross01; cross01=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			double d_cross01[3]; d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;
			labframe_transformation(1, d_cross01, d_frag1, d_photon_dir);  d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_sum_test_CH1);d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_frag2);       d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_el);          d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;

			double phiphot  = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2]);
			double phifrag1 = 180/PI*atan2(d_frag1[1],d_frag1[2]);
			double phifrag2 = 180/PI*atan2(d_frag2[1],d_frag2[2]);
			double phifrag3 = 180/PI*atan2(d_sum_test_CH1[1],d_sum_test_CH1[2]);
			double phiel =    180/PI*atan2(d_el[1],d_el[2]);
			double mirrorphiphot = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2])+180;
			if(mirrorphiphot>180){mirrorphiphot=mirrorphiphot-360;}

			double thetaphot = -(d_photon_dir[0]/(sqrt(d_photon_dir[0]*d_photon_dir[0]+d_photon_dir[1]*d_photon_dir[1]+d_photon_dir[2]*d_photon_dir[2])));
			double thetafrag1= -(d_frag1[0]/(sqrt(d_frag1[0]*d_frag1[0]+d_frag1[1]*d_frag1[1]+d_frag1[2]*d_frag1[2])));
			double thetafrag2= -(d_frag2[0]/(sqrt(d_frag2[0]*d_frag2[0]+d_frag2[1]*d_frag2[1]+d_frag2[2]*d_frag2[2])));
			double thetafrag3= -(d_sum_test_CH1[0]/(sqrt(d_sum_test_CH1[0]*d_sum_test_CH1[0]+d_sum_test_CH1[1]*d_sum_test_CH1[1]+d_sum_test_CH1[2]*d_sum_test_CH1[2])));
			double thetael =   -(d_el[0]/(sqrt(d_el[0]*d_el[0]+d_el[1]*d_el[1]+d_el[2]*d_el[2])));

			Hist->fill(AUTO,"MFAPD P_all",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			Hist->fill(AUTO,"MFAPD P_all_-theta",phiel ,-thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			strcpy(ndir, tdir);
			strcat(ndir, "/CH1+1/MFPADs");

			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333 - 0.833333 - thetaphot) < 0.1666666 && fabs(j * 30 - 165 - phiphot) < 15) {
						sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", k*0.333 - 1, j * 30 - 180);
						Hist->fill(6300+j*10+k, title, phiel, thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}
			strcpy(ndir, tdir);
			strcat(ndir, "/CH1+1/MFPADs_flipped");
			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333-0.833333+thetaphot)<0.1666666 && fabs(j*30-165-mirrorphiphot)<15) {
						sprintf_s(title, 100, "MFPAD3D_engate_-costheta_%4.2f_mirrphi_%d",k*0.333-1, j*30-180);
						Hist->fill(6800+j*10+k,title,phiel,thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}

			strcpy(ndir, tdir);
			strcat(ndir, "/CH1+1");
			//Newton plots
			CH_vector cross02; cross02=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_CH1);
			Coordinate_System multitest_mf = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom);
			CH_vector e_multitest_mf=multitest_mf.project_vector(e[0]->mom);
			CH_vector light_multitest_mf=multitest_mf.project_vector(photon_dir);
			CH_vector multitest_mf_mom[3];
			multitest_mf_mom[0] = multitest_mf.project_vector(mol->ion[0]->mom);
			multitest_mf_mom[1] = multitest_mf.project_vector(mol->ion[1]->mom);
			multitest_mf_mom[2] = multitest_mf.project_vector(neutral_test_mom_CH1);

			Hist->fill(AUTO,"Newton plot",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);

			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-1.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.7 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00 &&
				multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()>-0.3 &&  multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()<0.5  && 
				-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()<0.00 && -multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()>-0.7 &&
				//multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() > 40 && multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() < 120 &&
				//multitest_mf.project_vector(neutral_test_mom_CH1.Norm()).Phi_deg() > -155 && multitest_mf.project_vector(neutral_test_mom_CH1.Norm()).Phi_deg() < -90 &&
				mol->mom_cm.Mag() < sum_max) {
				Hist->fill(AUTO,"MFAPD P_all_gated_1",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_1",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_1",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_1", neutral_test_mom_CH1.y, neutral_test_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-0.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.5 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00) {
				Hist->fill(AUTO,"MFAPD P_all_gated_2",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_2",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_2",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_2", neutral_test_mom_CH1.y, neutral_test_mom_CH1.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
		}

		if (CH == 9 && ID == 2 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH9/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/CH2+1");
			CH_vector neutral_test_mom_CH2;

			double d_sum_test_CH2[3];
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_test_mom_CH2.x = -mol->mom_cm.x; neutral_test_mom_CH2.y = -mol->mom_cm.y; neutral_test_mom_CH2.z = -mol->mom_cm.z;
				d_sum_test_CH2[0]=-mol->mom_cm.x; d_sum_test_CH2[1]=-mol->mom_cm.y; d_sum_test_CH2[2]=-mol->mom_cm.z;
			}
			CH_vector cross01; cross01=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			double d_cross01[3]; d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;
			labframe_transformation(1, d_cross01, d_frag1, d_photon_dir);  d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_sum_test_CH2);d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_frag2);       d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_el);          d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;

			double phiphot  = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2]);
			double phifrag1 = 180/PI*atan2(d_frag1[1],d_frag1[2]);
			double phifrag2 = 180/PI*atan2(d_frag2[1],d_frag2[2]);
			double phifrag3 = 180/PI*atan2(d_sum_test_CH2[1],d_sum_test_CH2[2]);
			double phiel =    180/PI*atan2(d_el[1],d_el[2]);
			double mirrorphiphot = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2])+180;
			if(mirrorphiphot>180){mirrorphiphot=mirrorphiphot-360;}

			double thetaphot = -(d_photon_dir[0]/(sqrt(d_photon_dir[0]*d_photon_dir[0]+d_photon_dir[1]*d_photon_dir[1]+d_photon_dir[2]*d_photon_dir[2])));
			double thetafrag1= -(d_frag1[0]/(sqrt(d_frag1[0]*d_frag1[0]+d_frag1[1]*d_frag1[1]+d_frag1[2]*d_frag1[2])));
			double thetafrag2= -(d_frag2[0]/(sqrt(d_frag2[0]*d_frag2[0]+d_frag2[1]*d_frag2[1]+d_frag2[2]*d_frag2[2])));
			double thetafrag3= -(d_sum_test_CH2[0]/(sqrt(d_sum_test_CH2[0]*d_sum_test_CH2[0]+d_sum_test_CH2[1]*d_sum_test_CH2[1]+d_sum_test_CH2[2]*d_sum_test_CH2[2])));
			double thetael =   -(d_el[0]/(sqrt(d_el[0]*d_el[0]+d_el[1]*d_el[1]+d_el[2]*d_el[2])));

			Hist->fill(AUTO,"MFAPD P_all",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			Hist->fill(AUTO,"MFAPD P_all_-theta",phiel ,-thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			strcpy(ndir, tdir);
			strcat(ndir, "/CH2+1/MFPADs");

			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333 - 0.833333 - thetaphot) < 0.1666666 && fabs(j * 30 - 165 - phiphot) < 15) {
						sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", k*0.333 - 1, j * 30 - 180);
						Hist->fill(7300+j*10+k, title, phiel, thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}
			strcpy(ndir, tdir);
			strcat(ndir, "/CH2+1/MFPADs_flipped");
			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333-0.833333+thetaphot)<0.1666666 && fabs(j*30-165-mirrorphiphot)<15) {
						sprintf_s(title, 100, "MFPAD3D_engate_-costheta_%4.2f_mirrphi_%d",k*0.333-1, j*30-180);
						Hist->fill(7800+j*10+k,title,phiel,thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}

			strcpy(ndir, tdir);
			strcat(ndir, "/CH2+1");
			//Newton plots
			CH_vector cross02; cross02=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_CH2);
			Coordinate_System multitest_mf = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom);
			CH_vector e_multitest_mf=multitest_mf.project_vector(e[0]->mom);
			CH_vector light_multitest_mf=multitest_mf.project_vector(photon_dir);
			CH_vector multitest_mf_mom[3];
			multitest_mf_mom[0] = multitest_mf.project_vector(mol->ion[0]->mom);
			multitest_mf_mom[1] = multitest_mf.project_vector(mol->ion[1]->mom);
			multitest_mf_mom[2] = multitest_mf.project_vector(neutral_test_mom_CH2);

			Hist->fill(AUTO,"Newton plot",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);

			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-1.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.7 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00 &&
				multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()>-0.3 &&  multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()<0.5  && 
				-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()<0.00 && -multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()>-0.7 &&
				//multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() > 40 && multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() < 120 &&
				//multitest_mf.project_vector(neutral_test_mom_CH2.Norm()).Phi_deg() > -155 && multitest_mf.project_vector(neutral_test_mom_ALL.Norm()).Phi_deg() < -90 &&
				mol->mom_cm.Mag() < sum_max) {
				Hist->fill(AUTO,"MFAPD P_all_gated_1",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_1",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_1",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_1", neutral_test_mom_CH2.y, neutral_test_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-0.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.5 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00) {
				Hist->fill(AUTO,"MFAPD P_all_gated_2",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_2",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_2",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_2", neutral_test_mom_CH2.y, neutral_test_mom_CH2.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
		}

		if (CH == 9 && ID == 3 && mol->valid && e[0]->energy()> emin_el && e[0]->energy()< emax_el) {
			strcpy(ndir, "angular_distr_el/CH9/test");
			strcpy(tdir, ndir);
			strcat(ndir, "/CH3+1");
			CH_vector neutral_test_mom_CH3;

			double d_sum_test_CH3[3];
			if (mol->mom_cm.Mag()>sum_min && mol->mom_cm.Mag()<sum_max) {
				neutral_test_mom_CH3.x = -mol->mom_cm.x; neutral_test_mom_CH3.y = -mol->mom_cm.y; neutral_test_mom_CH3.z = -mol->mom_cm.z;
				d_sum_test_CH3[0]=-mol->mom_cm.x; d_sum_test_CH3[1]=-mol->mom_cm.y; d_sum_test_CH3[2]=-mol->mom_cm.z;
			}
			CH_vector cross01; cross01=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			double d_cross01[3]; d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;
			labframe_transformation(1, d_cross01, d_frag1, d_photon_dir);  d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_sum_test_CH3);d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_frag2);       d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z; d_frag1[0]=mol->ion[1]->mom.x; d_frag1[1]=mol->ion[1]->mom.y; d_frag1[2]=mol->ion[1]->mom.z;
			labframe_transformation(1, d_cross01, d_frag1, d_el);          d_cross01[0]=cross01.x; d_cross01[1]=cross01.y; d_cross01[2]=cross01.z;

			double phiphot  = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2]);
			double phifrag1 = 180/PI*atan2(d_frag1[1],d_frag1[2]);
			double phifrag2 = 180/PI*atan2(d_frag2[1],d_frag2[2]);
			double phifrag3 = 180/PI*atan2(d_sum_test_CH3[1],d_sum_test_CH3[2]);
			double phiel =    180/PI*atan2(d_el[1],d_el[2]);
			double mirrorphiphot = 180/PI*atan2(d_photon_dir[1],d_photon_dir[2])+180;
			if(mirrorphiphot>180){mirrorphiphot=mirrorphiphot-360;}

			double thetaphot = -(d_photon_dir[0]/(sqrt(d_photon_dir[0]*d_photon_dir[0]+d_photon_dir[1]*d_photon_dir[1]+d_photon_dir[2]*d_photon_dir[2])));
			double thetafrag1= -(d_frag1[0]/(sqrt(d_frag1[0]*d_frag1[0]+d_frag1[1]*d_frag1[1]+d_frag1[2]*d_frag1[2])));
			double thetafrag2= -(d_frag2[0]/(sqrt(d_frag2[0]*d_frag2[0]+d_frag2[1]*d_frag2[1]+d_frag2[2]*d_frag2[2])));
			double thetafrag3= -(d_sum_test_CH3[0]/(sqrt(d_sum_test_CH3[0]*d_sum_test_CH3[0]+d_sum_test_CH3[1]*d_sum_test_CH3[1]+d_sum_test_CH3[2]*d_sum_test_CH3[2])));
			double thetael =   -(d_el[0]/(sqrt(d_el[0]*d_el[0]+d_el[1]*d_el[1]+d_el[2]*d_el[2])));

			Hist->fill(AUTO,"MFAPD P_all",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			Hist->fill(AUTO,"MFAPD P_all_-theta",phiel ,-thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
			strcpy(ndir, tdir);
			strcat(ndir, "/CH3+1/MFPADs");

			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333 - 0.833333 - thetaphot) < 0.1666666 && fabs(j * 30 - 165 - phiphot) < 15) {
						sprintf_s(title, 100, "MFPAD3D_engate_costheta_%4.2f_phi_%d", k*0.333 - 1, j * 30 - 180);
						Hist->fill(8300+j*10+k, title, phiel, thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}
			strcpy(ndir, tdir);
			strcat(ndir, "/CH3+1/MFPADs_flipped");
			for (int k = 0; k < 6; k++) {
				for (int j = 0; j < 12; j++) {
					if (fabs(k*0.33333-0.833333+thetaphot)<0.1666666 && fabs(j*30-165-mirrorphiphot)<15) {
						sprintf_s(title, 100, "MFPAD3D_engate_-costheta_%4.2f_mirrphi_%d",k*0.333-1, j*30-180);
						Hist->fill(8800+j*10+k,title,phiel,thetael, 1., "molecular frame angular distribution", 12, -180., 180., "Phi", 6, -1., 1., "cos(theta)", ndir);
					}
				}
			}

			strcpy(ndir, tdir);
			strcat(ndir, "/CH3+1");
			//Newton plots
			CH_vector cross02; cross02=mol->ion[0]->mom.Cross(mol->ion[1]->mom);
			//Coordinate_System multitest_mf = Coordinate_System(cross01,neutral_test_mom_CH3);
			Coordinate_System multitest_mf = Coordinate_System(mol->ion[0]->mom, mol->ion[1]->mom);
			CH_vector e_multitest_mf=multitest_mf.project_vector(e[0]->mom);
			CH_vector light_multitest_mf=multitest_mf.project_vector(photon_dir);
			CH_vector multitest_mf_mom[3];
			multitest_mf_mom[0] = multitest_mf.project_vector(mol->ion[0]->mom);
			multitest_mf_mom[1] = multitest_mf.project_vector(mol->ion[1]->mom);
			multitest_mf_mom[2] = multitest_mf.project_vector(neutral_test_mom_CH3);

			Hist->fill(AUTO,"Newton plot",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
			Hist->fill(AUTO-1,"Newton plot",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);

			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-1.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.7 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00 &&
				multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()>-0.3 &&  multitest_mf_mom[2].x/mol->ion[0]->mom.Mag()<0.5  && 
				-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()<0.00 && -multitest_mf_mom[2].y/mol->ion[0]->mom.Mag()>-0.7 &&
				//multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() > 40 && multitest_mf.project_vector(mol->ion[0]->mom.Norm()).Phi_deg() < 120 &&
				//multitest_mf.project_vector(neutral_test_mom_CH3.Norm()).Phi_deg() > -155 && multitest_mf.project_vector(neutral_test_mom_CH3.Norm()).Phi_deg() < -90 &&
				mol->mom_cm.Mag() < sum_max) {
				Hist->fill(AUTO,"MFAPD P_all_gated_1",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_1", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_1", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_1",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_1",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_1", neutral_test_mom_CH3.y, neutral_test_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
			if (multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()>-0.8 &&  multitest_mf_mom[1].x/mol->ion[0]->mom.Mag()<-0.5 &&
				-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()<1.00 && -multitest_mf_mom[1].y/mol->ion[0]->mom.Mag()>0.00) {
				Hist->fill(AUTO,"MFAPD P_all_gated_2",phiel ,thetael ,1.,"molecular frame angular distribution", 36, -180., 180., "Phi", 18, -1.0, 1.0, "cos(theta)", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs x_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[0] y vs z_gated_2", mol->ion[0]->mom.x, mol->ion[0]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs x_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.y, 1., "p_x vs. p_y (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_y [a.u.]", ndir);
				Hist->fill(AUTO, "mom ion[1] y vs z_gated_2", mol->ion[1]->mom.x, mol->ion[1]->mom.z, 1., "p_x vs. p_z (no pol)", 100, -200., 200, "p_x [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
				Hist->fill(AUTO,"Newton plot_2",multitest_mf_mom[1].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[1].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO-1,"Newton plot_2",multitest_mf_mom[2].x/mol->ion[0]->mom.Mag(),-multitest_mf_mom[2].y/mol->ion[0]->mom.Mag(),1.,"Newton plot",80,-2.5,2.5,"px",120,-3.0,3.0,"py",ndir);
				Hist->fill(AUTO, "check_neutral_momentum_py vs pz_2", neutral_test_mom_CH3.y, neutral_test_mom_CH3.z, 1., "p_y vs. p_z (no pol)", 100, -200., 200, "p_y [a.u.]", 100, -200., 200, "p_z [a.u.]", ndir);
			}
		}
#endif
	} //void
} //NAMESPACE_CH		
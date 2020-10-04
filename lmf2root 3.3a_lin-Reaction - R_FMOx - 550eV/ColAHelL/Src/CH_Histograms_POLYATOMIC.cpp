#include "CH_Histograms.h"

// Standard histograms for reactions of "POLYATOMIC" and "POLYATOMIC_ELEC" type
namespace CH
{
	// there is nothing to plot in case we did not have a reaction definition
	void histograms_class::plot_polyatomic(reaction_struct *cur_reaction, polyatomic_class * big_mol, electron_class ** e, int ehit, double *scan_val)
	{
		if(!big_mol->valid)
			return;

		histo_handler *Hist = this->Hist_ions;
		// histogram index stuff
		int hmax = HISTS_PER_CHANNEL; // overall number of ColAHelL histograms per reaction/channel
		int hoff = hmax*16*(int)cur_reaction->reac_num;

		// limit to 16 hits.. 
		if(ehit>16)
			ehit=16;

		char reaction_dir[360];
		char x_axis_title[80];
		char y_axis_title[80];

		char rootdir[360];
		char ndir[360];
		char alldir[360];
		char tdir[360];

		if(this->use_master_folder) {
			strcpy(reaction_dir,this->CH_master_folder);
			strcat(reaction_dir,cur_reaction->name);
		} else
			strcpy(reaction_dir,cur_reaction->name);		
		
		strcpy(rootdir,reaction_dir);

		double rtof_sum1 = 0.;
		double rtof_sum2 = 0.;

		for(int i=0; i<big_mol->number_of_ions;i++) {
			strcpy(ndir,rootdir);
			sprintf(tdir,"/ion%i_m%i_q%i", i+1, int(big_mol->ion[i]->raw.m+0.001), int(big_mol->ion[i]->raw.q + 0.001));
			strcat(ndir,tdir);
			strcpy(tdir,ndir);

			Hist->fill(hoff+hmax*i+0,"ion energy",big_mol->ion[i]->energy(),1.,"ion energy",50,0.0,25.0,"ion energy [eV]",strcat(ndir,"/energy"));
			Hist->fill(hoff+hmax*i+1,"phi vs ion energy",big_mol->ion[i]->raw.phi,big_mol->ion[i]->energy(),1.,"ion energy vs. phi on detector",72,-180.0,180.0,"phi [deg]",50,0.0,25.0,"ion energy [eV]",ndir); 
			Hist->fill(hoff+hmax*i+2,"ctheta vs ion energy",big_mol->ion[i]->mom.Cos_Theta(),big_mol->ion[i]->energy(),1.,"ion energy vs. cos(theta_z)",72,-1.0,1.0,"ctheta",50,0.0,25.0,"ion energy [eV]",ndir); 
							
			strcpy(ndir,tdir);
			Hist->fill(hoff+hmax*i+3,"p_x vs p_y",big_mol->ion[i]->mom.x,big_mol->ion[i]->mom.y,1.,"p_x vs. p_y",300,-300.,300.,"p_x [a.u.]",300,-300.,300.,"p_y [a.u.]",strcat(ndir,"/momenta"));
			Hist->fill(hoff+hmax*i+4,"p_x vs p_z",big_mol->ion[i]->mom.x,big_mol->ion[i]->mom.z,1.,"p_x vs. p_z",300,-300.,300.,"p_x [a.u.]",300,-300.,300.,"p_z [a.u.]",ndir);
			Hist->fill(hoff+hmax*i+5,"p_y vs p_z",big_mol->ion[i]->mom.y,big_mol->ion[i]->mom.z,1.,"p_y vs. p_z",300,-300.,300.,"p_y [a.u.]",300,-300.,300.,"p_z [a.u.]",ndir);
			Hist->fill(hoff+hmax*i+6,"p_mag",big_mol->ion[i]->mom.Mag(),1.,"|p|",400,0.0,400,"|p| [a.u.]",ndir);

			strcpy(ndir,tdir);
			Hist->fill(hoff+hmax*i+7,"phi",big_mol->ion[i]->mom.Phi_deg(),1.,"phi_x in the labframe",72,-180.0,180.0,"phi_x [deg]",strcat(ndir,"/angles_labframe"));
			Hist->fill(hoff+hmax*i+10,"ctheta",big_mol->ion[i]->mom.Cos_Theta(),1.,"cos(theta_x) in the labframe",72,-1.0,1.0,"cos(theta_x)",ndir);

			strcpy(ndir,tdir);
			Hist->fill(hoff+hmax*i+13,"position",big_mol->ion[i]->raw.data.x,big_mol->ion[i]->raw.data.y,1.,"position",400,-1.*rdet_size,rdet_size,"x [mm]",400,-1.*rdet_size,rdet_size,"y [mm]",strcat(ndir,"/raw"));
			Hist->fill(hoff+hmax*i+14,"tof",big_mol->ion[i]->raw.data.tof,1.,"time-of-flight",4000,-1000.,20000.,"tof [ns]",ndir);
			Hist->fill(hoff+hmax*i+15,"wiggle",big_mol->ion[i]->raw.data.tof,(sqrt(big_mol->ion[i]->raw.data.x*big_mol->ion[i]->raw.data.x+big_mol->ion[i]->raw.data.y*big_mol->ion[i]->raw.data.y)),1.,"wiggles",500,-25.,10000.,"tof [ns]",200,-1.,rdet_size,"r [mm]",ndir);
			Hist->fill(hoff+hmax*i+16,"fish x",big_mol->ion[i]->raw.data.tof,big_mol->ion[i]->raw.data.x,1.,"x-fish",500,-25.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"x [mm]",ndir);
			Hist->fill(hoff+hmax*i+17,"fish y",big_mol->ion[i]->raw.data.tof,big_mol->ion[i]->raw.data.y,1.,"y-fish",500,-25.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"y [mm]",ndir);
			if ( fabs(big_mol->ion[i]->raw.data.y)<10. ) 
				Hist->fill(hoff+hmax*i+18,"fish filet x",big_mol->ion[i]->raw.data.tof,big_mol->ion[i]->raw.data.x,1.,"x-fish filet",500,-25.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"x [mm]",ndir);
			if ( fabs(big_mol->ion[i]->raw.data.x)<10. ) 
				Hist->fill(hoff+hmax*i+19,"fish filet y",big_mol->ion[i]->raw.data.tof,big_mol->ion[i]->raw.data.y,1.,"y-fish filet",500,-25.,10000.,"tof [ns]",200,-1.*rdet_size,rdet_size,"y [mm]",ndir);
				
			if(i<(big_mol->number_of_ions-0.01)/2)
				rtof_sum1 += big_mol->ion[i]->raw.data.tof;
			else 
				rtof_sum2 += big_mol->ion[i]->raw.data.tof;
		}

		// fill standard polyatomic coincidence histograms
		strcpy(ndir,rootdir);
		sprintf(tdir,"/polyatomic");
		strcat(ndir,tdir);
		strcpy(tdir,ndir);

		if (cur_reaction->use_ion_matrix == 1)
		{
			for (int i=0; i<16; i++){
				for (int j=0; j<16; j++) //goes through each ion species
				{
					if (big_mol->ion_matrix[i][j] == 1) //if hit tof is within the bounds of the ion species
						Hist->fill2(hoff+hmax+20,"PIM",i,j,1.,"Ion matrix",16,0,16,"Hit number",16,0,16,"ion species",ndir);
				}
			}

			Hist->fill(hoff+hmax+21,"pzsum_[0]",big_mol->value_to_sort[0],1.,"sum_momentum",400,0.,100.,"momentum (a.u.)",ndir);
			Hist->fill(hoff+hmax+22,"pzsum_[1]",big_mol->value_to_sort[1],1.,"sum_momentum",400,0.,100.,"momentum (a.u.)",ndir);
		}

		Hist->fill(hoff+hmax+23,"KER",big_mol->KER(),1.,"kinetic energy release",400,0.,100.,"KER [eV]",ndir);

		Hist->fill(hoff+hmax+24,"p_sumx vs p_sumy large",big_mol->mom_cm.x,big_mol->mom_cm.y,1.,"center of mass momentum x vs y",300,-300.,300.,"p_cmx [a.u.]",300,-300.,300.,"p_cmx [a.u.]",ndir);
		Hist->fill(hoff+hmax+25,"p_sumx vs p_sumz large",big_mol->mom_cm.x,big_mol->mom_cm.z,1.,"center of mass momentum x vs z",300,-300.,300.,"p_cmx [a.u.]",300,-300.,300.,"p_cmz [a.u.]",ndir);
		Hist->fill(hoff+hmax+26,"p_sumy vs p_sumz large",big_mol->mom_cm.y,big_mol->mom_cm.z,1.,"center of mass momentum y vs z",300,-300.,300.,"p_cmy [a.u.]",300,-300.,300.,"p_cmz [a.u.]",ndir);
		Hist->fill(hoff+hmax+27,"p_sumx vs p_sumy",big_mol->mom_cm.x,big_mol->mom_cm.y,1.,"center of mass momentum x vs y",100,-20.,20.,"p_cmx [a.u.]",100,-20.,20.,"p_cmx [a.u.]",ndir);
		Hist->fill(hoff+hmax+28,"p_sumx vs p_sumz",big_mol->mom_cm.x,big_mol->mom_cm.z,1.,"center of mass momentum x vs z",100,-20.,20.,"p_cmx [a.u.]",100,-20.,20.,"p_cmz [a.u.]",ndir);
		Hist->fill(hoff+hmax+29,"p_sumy vs p_sumz",big_mol->mom_cm.y,big_mol->mom_cm.z,1.,"center of mass momentum y vs z",100,-20.,20.,"p_cmy [a.u.]",100,-20.,20.,"p_cmz [a.u.]",ndir);
					
		Hist->fill(hoff+hmax+30,"p_sum_magnitude",big_mol->mom_cm.Mag(),1.,"p_sum_magnitude",400,0.,400.,"momentum [a.u.]",ndir);
		if (big_mol->number_of_ions > 2)
		{
			// Create Dalitz plot with the ions provided by the Dalitz array
			sprintf(x_axis_title,"(|p|_%i - |p|_%i)/sqrt(3)",int(big_mol->ion[cur_reaction->Dalitz_array[1]]->raw.m+0.01),int(big_mol->ion[cur_reaction->Dalitz_array[2]]->raw.m+0.01));
			sprintf(y_axis_title,"|p|_%i - 1/3",int(big_mol->ion[cur_reaction->Dalitz_array[0]]->raw.m+0.01));
			Hist->fill(hoff+hmax+32,"Dalitz plot",(big_mol->ion[cur_reaction->Dalitz_array[1]]->mom.Mag() - big_mol->ion[cur_reaction->Dalitz_array[2]]->mom.Mag())/(big_mol->momentum_magnitude_sum*sqrt(3.)),big_mol->ion[cur_reaction->Dalitz_array[0]]->mom.Mag()/big_mol->momentum_magnitude_sum - 1./3.,1.,"Dalitz plot for momenta",100,-0.5,0.5,x_axis_title,100,-0.5,0.5,y_axis_title,ndir);
			
			// Create Newton plot with the ions provided by the Newton array
			Coordinate_System Nframe = Coordinate_System(big_mol->ion[cur_reaction->Newton_array[0]]->mom, big_mol->ion[cur_reaction->Newton_array[1]]->mom);
			CH_vector NFrame_mom[3];
			for (int i = 0; i<3; i++)
			{
				NFrame_mom[i] = Nframe.project_vector(big_mol->ion[cur_reaction->Newton_array[i]]->mom);
			}
			Hist->fill(hoff+hmax+33,"Newton plot",NFrame_mom[1].z/big_mol->ion[cur_reaction->Newton_array[0]]->mom.Mag(),NFrame_mom[1].x/big_mol->ion[cur_reaction->Newton_array[0]]->mom.Mag(),1.,"Newton plot",80,-2.5,1.5,"pz",120,-3.0,3.0,"py",ndir);
			Hist->fill(hoff+hmax+33,"Newton plot",NFrame_mom[2].z/big_mol->ion[cur_reaction->Newton_array[0]]->mom.Mag(),NFrame_mom[2].x/big_mol->ion[cur_reaction->Newton_array[0]]->mom.Mag(),1.,"Newton plot",80,-2.5,1.5,"pz",120,-3.0,3.0,"py",ndir);

			Hist->fill(hoff+hmax+34,"Ion Molframe",NFrame_mom[1].z,NFrame_mom[1].x,1.,"Ions Molframe",80,-240.0,240,"pz",80,-240.0,240.0,"py",ndir);
			Hist->fill(hoff+hmax+34,"Ion Molframe",NFrame_mom[2].z,NFrame_mom[2].x,1.,"Ions Molframe",80,-240.0,240,"pz",80,-240.0,240.0,"py",ndir);

			Hist->fill(hoff+hmax+35,"Ion Molframe EN",NFrame_mom[1].z,NFrame_mom[1].x,NFrame_mom[0].Mag(),1.,"Ions Molframe",50,-80.0,80,"pz",50,-80.0,80.0,"py",10,0.0,200.0,"pmag",ndir);
			Hist->fill(hoff+hmax+35,"Ion Molframe EN",NFrame_mom[2].z,NFrame_mom[2].x,NFrame_mom[0].Mag(),1.,"Ions Molframe",50,-80.0,80,"pz",50,-80.0,80.0,"py",10,0.0,200.0,"pmag",ndir);

			double BondAngle = big_mol->ion[cur_reaction->Newton_array[0]]->mom.Angle_deg(big_mol->ion[cur_reaction->Newton_array[1]]->mom);
			Hist->fill(hoff+hmax+36,"Bond angle vs KER",BondAngle, big_mol->KER(),1.,"",72,-180.0,180.0,"theta",200,0.0,50.0,"KER [eV]",ndir);
			Hist->fill(hoff+hmax+37,"cos(Bond angle) vs KER",cos(BondAngle/180.*PI), big_mol->KER(),1.0,"",36,-1.0,1.0,"cos(theta)",200,0.0,50.0,"KER [eV]",ndir);

			Coordinate_System mf = Coordinate_System(big_mol->ion[0]->mom.Norm() + big_mol->ion[1]->mom.Norm(), big_mol->ion[0]->mom.Norm() - big_mol->ion[1]->mom.Norm());			
			CH_vector i1_MF;
			i1_MF = mf.project_vector(big_mol->ion[0]->mom);
			CH_vector i2_MF;
			i2_MF = mf.project_vector(big_mol->ion[1]->mom);
			Hist->fill(hoff+hmax+38,"Ion 3D-Molframe",i1_MF.Phi_deg(),i1_MF.Cos_Theta(),1.,"Molecular frame angular distribution (Ion Position)",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);
			Hist->fill(hoff+hmax+38,"Ion 3D-Molframe",i2_MF.Phi_deg(),i2_MF.Cos_Theta(),1.,"Molecular frame angular distribution (Ion Position)",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);

			char title[50];
			sprintf_s(title,50,"%s%li%s","PI",big_mol->number_of_ions,"CO");
			sprintf_s(x_axis_title,80,"%s%li%s","Sum 1...",int((big_mol->number_of_ions-0.01)/2)+1," ions");
			sprintf_s(y_axis_title,80,"%s%li%s%li%s","Sum ",int((big_mol->number_of_ions-0.01)/2)+2,"...",big_mol->number_of_ions," ions");
			Hist->fill2(hoff+hmax+43,title,rtof_sum1,rtof_sum2,1.,title,500,00.,20000.,x_axis_title,500,0.,20000.,y_axis_title,ndir);
		}

		// ion/electron coincidences
		strcpy(alldir,rootdir);
		strcat(alldir,"/polyatomic_coincidence/all_hits");

		for(int i=0;i<(int)ehit;i++) {	
			strcpy(ndir,rootdir);
			sprintf(tdir,"/polyatomic_coincidence/elec_hit_%i",i);
			strcat(ndir,tdir);
//			strcpy(tdir,ndir);
		
			if(e[i]->valid) {
				Hist->fill(hoff+hmax*i+40,"KER vs electron energy",big_mol->KER(),e[i]->energy(),1.,"kinetic energy release vs electron energy",60,0.,60.,"KER [eV]",60,0.0,20.0,"electron energy [eV]",ndir);

				Coordinate_System mf = Coordinate_System(big_mol->ion[0]->mom, big_mol->ion[1]->mom);
				CH_vector e_MF;
				e_MF = mf.project_vector(e[i]->mom);
				Hist->fill(hoff+hmax*i+41,"MFPAD3D",e_MF.Phi_deg(),e_MF.Cos_Theta(),1.,"molecular frame angular distribution",36,-180.,180.,"Phi",24,-1.0,1.0,"cos(theta)",ndir);
			}

			// Merge all electrons in one plot
			if(e[i]->valid) {
				Hist->fill(hoff+hmax+42,"KER vs electron energy",big_mol->KER(),e[i]->energy(),1.,"kinetic energy release vs electron energy",60,0.,60.,"KER [eV]",60,0.0,20.0,"electron energy [eV]",alldir);
			}
		}

	}
}
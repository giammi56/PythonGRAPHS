#include "CH_Polyatomic.h"
#include "CH_Basics.h"
#include "CH_Ion.h"
#include "CH_FUN_Lowlevel.h"




namespace CH{

	polyatomic_class::polyatomic_class(){

	//PRIVATE
	this->nr_of_combinations = 0;
	this->max_nr_of_combinations = 30;

	
	//used hits
	this->used_hits = new int[16];
	for (int i = 0; i< 16; i++)
		this->used_hits[i] = -1;
	//value to sort
	this->value_to_sort = new double[max_nr_of_combinations];
	for (int i = 0; i< max_nr_of_combinations; i++)
		this->value_to_sort[i] = 9999.;
	//index
	this->index = new int[this->max_nr_of_combinations];
	for (int i = 0; i< this->max_nr_of_combinations; i++)
		this->index[i] = i;
	

	//temp_mass
	this->temp_mass = new double[16];
	for (int i = 0; i< 16; i++)
		temp_mass[i] = 0.;
	//temp_charge
	this->temp_charge = new double[16];
	for (int i = 0; i< 16; i++)
		temp_charge[i] = 0.;
	//temp_ion
	this->temp_ion = new ion_class*[16];
	
	memset(temp_ion,0,16*sizeof(ion_class*));

	//matches
	this->matches = new int*[this->max_nr_of_combinations]; //max number of combinations
	for (int i = 0; i< this->max_nr_of_combinations; i++)
	{
		matches[i] = new int[16];
		for (int j = 0; j < 16; j++)
			matches[i][j] = -1;
	}



	//PUBLIC
	this->ion = new ion_class*[16];
	this->ion_matrix = new int*[16];
	//ion_matrix
	for (int i = 0; i< 16; i++)
	{
		this->ion[i] = new ion_class();
		ion_matrix[i] = new int[16];
		for (int j = 0; j <16; j++)
			ion_matrix[i][j] = 0;
	}
	
	this->momenta_in_mol_frame = new CH_vector[16];

	this->mfac = new mom_cor_param();

	this->mfac->dx = 0.0;
	this->mfac->dy = 0.0;
	this->mfac->dz = 0.0;

	this->mfac->x_stretch = 1.0;
	this->mfac->y_stretch = 1.0;
	this->mfac->z_stretch = 1.0;
	this->mfac->overall_stretch = 1.0;

	this->mfac->mir_x = false;
	this->mfac->mir_y = false;
	this->mfac->mir_z = false;
}



 polyatomic_class::~polyatomic_class(){

		if (used_hits){ delete[] used_hits; used_hits =0;}
		if (value_to_sort){ delete[] value_to_sort; value_to_sort =0;}
		if (index){ delete[] index; index =0;}
		if (temp_mass){ delete[] temp_mass; temp_mass =0;}
		if (temp_charge){ delete[] temp_charge; temp_charge =0;}
		if (temp_ion) {
				
			for (__int32 i=0;i<16;i++) {
				if (temp_ion[i]) {
					delete temp_ion[i];
					temp_ion[i] = 0;
				}
			}
			delete[] temp_ion; temp_ion = 0;
		}

	for (int i = 0; i < this->max_nr_of_combinations; i++)
	{
		if (matches[i]) {delete[] matches[i]; matches[i] = 0;}	
	}	
	if (matches){ delete[] matches; matches =0;}


	for (int i = 0; i < 16; i++)
	{
		if (ion[i]) {delete ion[i]; ion[i] = 0;}
		if (ion_matrix[i]) {delete[] ion_matrix[i]; ion_matrix[i] = 0;}
	}

	if (ion) {delete[] ion; ion = 0;}
	if (ion_matrix) {delete[] ion_matrix; ion_matrix = 0;}
	if (momenta_in_mol_frame) {delete[] momenta_in_mol_frame; momenta_in_mol_frame = 0;}

	if(mfac) {delete mfac; mfac=0;}

//	if (cm) {delete cm; cm = 0;}
}



 void polyatomic_class::set_ions(int number_of_ions, ion_class** ions){

	this->number_of_ions = number_of_ions;
	
	for (int i = 0; i < number_of_ions; i++)
	{
		ion[i]->clone_from(ions[i]);
	}
}


 void polyatomic_class::set_channel(int channel) {
		this->channel = channel;
	}

void polyatomic_class::set_channel(int channel, bool valid) {
		this->channel = channel;
		this->valid = valid;
	}

void polyatomic_class::set_valid(bool valid) {
		this->valid = valid;
	}

void polyatomic_class::set_mom_fac(mom_cor_param * fac) {
	this->mfac->dx = fac->dx;
	this->mfac->dy = fac->dy;
	this->mfac->dz = fac->dz;

	this->mfac->x_stretch = fac->x_stretch;
	this->mfac->y_stretch = fac->y_stretch;
	this->mfac->z_stretch = fac->z_stretch;
	this->mfac->overall_stretch = fac->overall_stretch;

	this->mfac->mir_x = fac->mir_x;
	this->mfac->mir_y = fac->mir_y;
	this->mfac->mir_z = fac->mir_z;
}

bool polyatomic_class::run_ion_matrix(reaction_struct* reaction, CH_event_struct *evti, spectrometer_class *spect){
	
	this->max_target_value = reaction->max_target_value;
	this->ambiguity_parameter = reaction->ambiguity_parameter;

	bool match_found = false;
	
	this->fill_ion_matrix(evti);
	this->sort_ion_matrix(0,evti,spect,reaction);
	match_found = this->sort_momenta();
	//at this point this->ion[n]->phy->mom are not correct --> recalculate momenta with best match
	if (match_found)
		this->process(spect);

	return match_found;
}




void polyatomic_class::reset() { //need to reset all arrays of polyatomic before going to next event
	this->incomplete = false;
	this->nr_of_combinations = 0;
	//ion matrix
	for (int i=0; i<16; i++) //Goes through each hit
	{
		this->ion[i]->raw.data.x = 99999.;
		this->ion[i]->raw.data.y = 99999.;
		this->ion[i]->raw.data.tof = 99999.;

		for (int j=0; j<16; j++) //goes through each ion species
			this->ion_matrix[i][j]= 0; 
	}
	//used_hits
	for (int i = 0; i< 16; i++)
		this->used_hits[i] = -1;
	//matches
	for (int i = 0; i< this->max_nr_of_combinations; i++)
	{
		for (int j = 0; j < 16; j++)
			this->matches[i][j] = -1;
	}
	//value_to_sort
	for (int i = 0; i< max_nr_of_combinations; i++)
		this->value_to_sort[i] = 99999.;
	//index
	for (int i = 0; i< max_nr_of_combinations; i++)
		this->index[i] = i;
	//temp mass
	for (int i = 0; i< 16; i++)
		this->temp_mass[i] = 0.;
	//temp_charge
	for (int i = 0; i< 16; i++)
		this->temp_charge[i] = 0.;
	/*
	for (int i = 0; i< 16; i++)
		this->molecular_frame[i] = {1.,1.,1.};
*/

}

void polyatomic_class::fill_ion_matrix(CH_event_struct *evti)
{


	for (int i=0; i < evti->r.num_hits; i++) //Goes through each hit
	{
		for (int j=0; j < this->number_of_ions; j++) //goes through each ion species
		{
			if ( fabs(evti->r.tof[i] - this->ion[j]->t_mean)  < this->ion[j]->t_width) //if hit tof is within the bounds of the ion species
			{ this->ion_matrix[i][j]= 1; //positive mark in PIM
			//Ueber->Hist->fill2(3,"PIM",i,j,1.,"Ion matrix",16,0,16,"Hit number",16,0,16,"ion species");
			}
			else
				this->ion_matrix[i][j]= 0;
		}
	}


}


#ifdef USE_ION_MATRIX_RECURSIVE   // switch between Martin's and Kilian's code...

bool polyatomic_class::sort_ion_matrix(int iteration, CH_event_struct *evti, spectrometer_class *spect, reaction_struct* reaction){
 
  bool Tag = false; //Just to keep track if a hit is usable
	int i = 0;
	bool method_check = true;
	
	//for (int i=0; i<nr_of_hits; i++) //loop through each hit, and check if it matches the ion we're looking for (the recursion goes through the ion species, the loop goes through the hits)
	while ((i < evti->r.num_hits) && (this->nr_of_combinations < this->max_nr_of_combinations))
	{
	if (this->ion_matrix[i][iteration] == 1)  //if hit i matches the iteration-th ion in channel (funktioniert mit [i][iteration], solange jeder Kanal eigene ion_matrix hat
		{
			Tag = true; //succesful hit

			for (int j=0; j< iteration; j++){ // go through previous hits, see if the same hit used twice
				if (i == this->used_hits[j])
				{ Tag = false; //unsuccessful hit
				}
			}  //end for (int j=0; j< iteration; j++)
			if (Tag) { //if no hit reused
				this->used_hits[iteration] = i; //records the succesful match: create a list of the used hits
				if (iteration < (this->number_of_ions - 1) )//check if we need more hits to complete the channel
				{ 
					int iter = iteration + 1; //on to the next iteration (and thus ion species) and call the recursion
					this->sort_ion_matrix(iter, evti, spect, reaction);
 
				} //#1, where algorithm arrives if a recursive call has been finished -->(drops out of all instances?)-->goes to next step of the while-loop
				
				else //if we don't need more hits, i.e. we're complete
				{
					
					//"merit function"					
					double pxsum = 0.;
					double pysum = 0.;
					double pzsum = 0.;
					
				
					for (int k = 0; k< this->number_of_ions; k++)
					{
						
						this->ion[k]->process(spect,this->ion[used_hits[k]]->raw.data.x,this->ion[used_hits[k]]->raw.data.y,this->ion[used_hits[k]]->raw.data.tof); //, reaction->interpolation_points[k], reaction->interpolation_adjustment_variables[k]);
						pxsum += this->ion[k]->mom.x;
						pysum += this->ion[k]->mom.y;
						pzsum += this->ion[k]->mom.z;
						
					}
					for (int m=0; m < this->number_of_ions; m++)
							{
							//records the found combination to matches 
							this->matches[this->nr_of_combinations][m] = this->used_hits[m];
							}
						//folgende zwei Zeilen nicht unbedingt nötig(?)
						//if ( fabs(this->mom_cm.z) < psum_min[2]) //is pz of the current match better than the previous one?
						//	{ psum_min[0] = this->mom_cm.x; psum_min[1] = this->mom_cm.y; psum_min[2] = this->mom_cm.z;}  //store the momenta with the best pz momentum conservation for later histograms
						//value_to_sort[this->nr_of_combinations] = fabs(pzsum);//fabs(rec_p[2]-pz_offset); //corrected for pz_offset in case of bad calibration
						value_to_sort[this->nr_of_combinations] = fabs(pxsum) + fabs(pysum) + fabs(pzsum); 
						this->nr_of_combinations+=1; 
						
										
				}//end else; leaves the current instance of sort_PIM and goes to #1 in the superordinate instance
			}// end: if (Tag)
	}//end: if (PIM[i][channel[iteration]] == 1)
	i++;
		}//end loop
	
	return true;

}

#else

bool polyatomic_class::sort_ion_matrix(int iteration, CH_event_struct *evti, spectrometer_class *spect, reaction_struct* reaction){
		
	int Allemöglichkeiten=1;
	for (int z=0; z< number_of_ions; z++)
	{
		Allemöglichkeiten*=evti->r.num_hits;
	}
	for (int j=0; j< Allemöglichkeiten; j++){     //idee: alle Kombinationen durchzähen, die es gibt und als zahlen
		//in der Basis b= number of hits darstellen und dabei alle kombinationen haben.
		int Zahlenhelfer1=j;									//Dann muss man nur noch die Doofen weg werfen 
		int Zahlenhelfer2=j;
		int Matrix_test=0;
		bool Tag= false;

		for (int i=0; i< number_of_ions; i++){					//

			this->used_hits[number_of_ions-1-i]= Zahlenhelfer1%evti->r.num_hits;
			Zahlenhelfer2 = Zahlenhelfer2/evti->r.num_hits; //integer durch einander teilen ist teilen und abrunden.
			Zahlenhelfer1=Zahlenhelfer2;
		}

		//		 printf ("%s \n", "Test1");
		for (int i=0; i< number_of_ions; i++)
		{																			//Test, ob nur Matrixeinträge benutzt werden, die eine 1 haben.
			Matrix_test=Matrix_test+this->ion_matrix[this->used_hits[i]][i];
			//		printf (" %d ", this->used_hits[i]);									//
		}																			//Sonst wird die Kombination verworfen
		if(Matrix_test != number_of_ions)											//
		{continue;}

		for (int i=0; i< number_of_ions; i++){										//Test, ob ein Hit mehrfach verwendet wird.
			for (int j=i+1; j< number_of_ions; j++){								//
				if (this->used_hits[i]==this->used_hits[j])							//
				{																	//
					Tag= true;
				}
			}
		}
		if (Tag==true)
		{continue;}

		double pxsum = 0.;
		double pysum = 0.;
		double pzsum = 0.;


		for (int k = 0; k< this->number_of_ions; k++)
		{
			this->ion[k]->process(spect,this->ion[used_hits[k]]->raw.data.x,this->ion[used_hits[k]]->raw.data.y,this->ion[used_hits[k]]->raw.data.tof); //, reaction->interpolation_points[k], reaction->interpolation_adjustment_variables[k]);
			pxsum += this->ion[k]->mom.x;
			pysum += this->ion[k]->mom.y;
			pzsum += this->ion[k]->mom.z;
		}

		// printf ("\n %s Number of combinations %d", "Kombinationen", this->nr_of_combinations );
		for (int m=0; m < this->number_of_ions; m++)
		{
			//records the found combination to matches 
			this->matches[this->nr_of_combinations][m] = this->used_hits[m];
									//printf ("  %d ", this->used_hits[m]);

		}
								// printf (" \n");
		this->value_to_sort[this->nr_of_combinations] = fabs(pxsum) + fabs(pysum) + fabs(pzsum); 

						//	printf ("    %f  ", value_to_sort[this->nr_of_combinations]);
		this->nr_of_combinations+=1; 
		if (nr_of_combinations>29){break;}//damit die begrenzte Mathes Box nicht überläuft..
	}

	return true;
}

#endif

bool polyatomic_class::sort_momenta(){

	double temp1; int temp2;
	//sort the matches 
	
	for (int i=0; i < this->max_nr_of_combinations; i++)
	{
		for (int j=i+1; j < this->max_nr_of_combinations; j++)
		{
			if (this->value_to_sort[j] < this->value_to_sort[i]) 
			{	temp1 = this->value_to_sort[i]; this->value_to_sort[i] = this->value_to_sort[j]; this->value_to_sort[j]=temp1;
				temp2 = this->index[i]; this->index[i]=this->index[j]; this->index[j]=temp2;	
			}
		}
	}
	
	//checking validity of solution 	
	for (int j = 0; j < this->max_nr_of_combinations; j++)
	{	
	if ( (this->matches[this->index[j]][0] != -1) && ( fabs(this->value_to_sort[j]) > fabs(this->value_to_sort[0]) )) //find an entry that differs from the first one (necessary in the case there are identical fragments in the list that lead to identical values)
	{
		if (  ( fabs(this->value_to_sort[0]) < this->max_target_value) && (fabs(this->value_to_sort[j]) > this->ambiguity_parameter*this->value_to_sort[0]))//if at least one match was found and the criteria from the config-file are fulfilled
		{
			//swapping ions (analog to diatomic class)
			//first assign masses and charges of to temporary array 
			for (int i = 0; i < this->number_of_ions; i++)
			{
				this->temp_mass[i] = this->ion[this->matches[this->index[0]][i]]->raw.m; 
				this->temp_charge[i] = this->ion[this->matches[this->index[0]][i]]->raw.q;
			}
			for (int i = 0; i< this->number_of_ions; i++)
			{
				this->ion[i]->raw.m = this->temp_mass[i]; 
				this->ion[i]->raw.q = this->temp_charge[i]; 			
			}

	
			for (int i = 0; i< this->number_of_ions; i++)
				temp_ion[this->matches[this->index[0]][i]] = this->ion[i]; 
	
			for (int i = 0; i< this->number_of_ions; i++)
				this->ion[i] = temp_ion[i]; 
	
			return true;
		}
		else
			return false;
		}
	//continue for-loop over j if both entries are the same
	}
	return false;
}




void polyatomic_class::process(spectrometer_class *spect) {
		
		for (int i = 0; i < number_of_ions; i++)
		{
		this->ion[i]->process(spect);
		//this->ion[i]->shift_stretch_mom(this->mfac);
		
		}

		calculate_sum_momenta();
		 
		// Create neutral ion from cm momentum
		if(this->incomplete) {
			this->ion[this->number_of_ions]->raw.data.mcp = 0.0;
			this->ion[this->number_of_ions]->raw.data.time = 0.0;
			this->ion[this->number_of_ions]->raw.data.tof = 0.0;
			this->ion[this->number_of_ions]->raw.data.x = 0.0;
			this->ion[this->number_of_ions]->raw.data.y = 0.0;
			this->ion[this->number_of_ions]->raw.phi = 0.0;
			this->ion[this->number_of_ions]->raw.method = -1.0;

			this->ion[this->number_of_ions]->mom.x = -this->mom_cm.x;
			this->ion[this->number_of_ions]->mom.y = -this->mom_cm.y;
			this->ion[this->number_of_ions]->mom.z = -this->mom_cm.z;

			this->mom_cm.x = 0.0;
			this->mom_cm.y = 0.0;
			this->mom_cm.z = 0.0;

			this->total_mass += this->ion[this->number_of_ions]->raw.m;
			this->momentum_magnitude_sum += this->ion[this->number_of_ions]->mom.Mag();
			this->number_of_ions++;
		}
	}



void polyatomic_class::calculate_sum_momenta(){
		
	this->mom_cm.x = 0.;
	this->mom_cm.y = 0.;
	this->mom_cm.z = 0.;

	this->total_mass = 0.;
	this->momentum_magnitude_sum = 0.;

	for (int i = 0; i < this->number_of_ions; i++)
	{
		this->mom_cm.x += this->ion[i]->mom.x;
		this->mom_cm.y += this->ion[i]->mom.y;
		this->mom_cm.z += this->ion[i]->mom.z;
		
		this->total_mass += this->ion[i]->raw.m;
		this->momentum_magnitude_sum += this->ion[i]->mom.Mag();
	}
		
}
double polyatomic_class::KER() {
	double res = 0.;
	for (int i = 0; i < this->number_of_ions; i++)
	{
		res += this->ion[i]->mom.Mag() * this->ion[i]->mom.Mag() / (2.0 * this->ion[i]->raw.m * MASSAU) * EVAU;
	}
return res;
}

bool polyatomic_class::check_validity(double var, double min, double max, bool not)
	{
		if(not) {
			if(var < min || var >= max)
				return false;
			else 
				return this->valid;
		} else {
			if(var > min && var <= max)
				return false;
			else 
				return this->valid;
		}
	}

void polyatomic_class::invalidate(rtag_struct* rtag) 
	{
		switch (rtag->tag_type) {
		
		case RTAG_PSUMX:
			this->valid = check_validity(this->mom_cm.x,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUMY:
			this->valid = check_validity(this->mom_cm.y,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUMZ:
			this->valid = check_validity(this->mom_cm.z,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUM:
			this->valid = check_validity(this->mom_cm.Mag(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_KER:
			this->valid = check_validity(this->KER(),rtag->min,rtag->max,rtag->negate);
			break;

		default:
			break;
		}

	}


void polyatomic_class::trafo_to_molframe(Coordinate_System frame){
	
	for (int i = 0; i<this->number_of_ions; i++)
	{
		momenta_in_mol_frame[i] = frame.project_vector(this->ion[i]->mom);

	}



}

/*
void polyatomic_class::labframe_transformation(int direction, CH_vector a, CH_vector b) //adapt Achim Czasch's labframe transformation
{


	CH_vector norm_vector_x;
	CH_vector norm_vector_y;
	CH_vector norm_vector_z;

	

	double transformation_matrix[3][3];
	
	norm_vector_x = a.Norm();

	norm_vector_z = a.Cross(b); //-b[1]*a[2]+b[2]*a[1];
	norm_vector_z = norm_vector_z.Norm();
	

	norm_vector_y = norm_vector_z.Cross(a); //-a[1]*norm_vector_z[2]+a[2]*norm_vector_z[1];
	norm_vector_y = norm_vector_y.Norm();
	

	transformation_matrix[0][0] = norm_vector_x.x;
	transformation_matrix[0][1] = norm_vector_x.y;
	transformation_matrix[0][2] = norm_vector_x.z;

	transformation_matrix[1][0] = norm_vector_y.x;
	transformation_matrix[1][1] = norm_vector_y.y;
	transformation_matrix[1][2] = norm_vector_y.z;

	transformation_matrix[2][0] = norm_vector_z.x;
	transformation_matrix[2][1] = norm_vector_z.y;
	transformation_matrix[2][2] = norm_vector_z.z;



//         If direction =  +1: This transforms into the new system
//         If direction =  -1: Reverse transformation

	
	

	if (direction==-1) {

			for (int i=0; i<this->number_of_ions; i++)
			{
			this->molecular_frame[i].x=this->ion[i]->phy->mom.x * transformation_matrix[0][0] + this->ion[i]->phy->mom.y * transformation_matrix[1][0] + this->ion[i]->phy->mom.z * transformation_matrix[2][0];
			this->molecular_frame[i].y=this->ion[i]->phy->mom.x * transformation_matrix[0][1] + this->ion[i]->phy->mom.y * transformation_matrix[1][1] + this->ion[i]->phy->mom.z * transformation_matrix[2][1];
			this->molecular_frame[i].z=this->ion[i]->phy->mom.x * transformation_matrix[0][2] + this->ion[i]->phy->mom.y * transformation_matrix[1][2] + this->ion[i]->phy->mom.z * transformation_matrix[2][2];
			this->molecular_frame[i] = this->molecular_frame[i]/a.Mag();
			
			}


	} else {
		if (direction==+1) {
				

			for (int i=0; i<this->number_of_ions; i++)
			{
			this->molecular_frame[i].x=this->ion[i]->phy->mom.x * transformation_matrix[0][0] + this->ion[i]->phy->mom.y * transformation_matrix[0][1] + this->ion[i]->phy->mom.z * transformation_matrix[0][2];
			this->molecular_frame[i].y=this->ion[i]->phy->mom.x * transformation_matrix[1][0] + this->ion[i]->phy->mom.y * transformation_matrix[1][1] + this->ion[i]->phy->mom.z * transformation_matrix[1][2];
			this->molecular_frame[i].z=this->ion[i]->phy->mom.x * transformation_matrix[2][0] + this->ion[i]->phy->mom.y * transformation_matrix[2][1] + this->ion[i]->phy->mom.z * transformation_matrix[2][2];
			this->molecular_frame[i] = this->molecular_frame[i]/a.Mag();
			
			}
				
		} else {

			for (int i=0; i<this->number_of_ions; i++)
			{
			this->molecular_frame[i].x = 0;	this->molecular_frame[i].y = 0;	this->molecular_frame[i].z = 0;
			}
				
		}
	}

	



} // end polyatomic_class::labframe_transformation(int a, int b)

*/
} //end namespace CH
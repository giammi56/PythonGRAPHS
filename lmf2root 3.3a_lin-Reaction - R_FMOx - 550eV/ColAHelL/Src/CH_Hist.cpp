//#include "stdafx.h"
#include "CH_Hist.h"
#include <string>
#include <vector>
#include "CH_console.h"
#include <map>
//#include <thread>
//#include <mutex>


namespace CH
{
using namespace std;



////////////////////////////////////////////  axis  ////////////////////////////////////////////
axis::axis(const char * Name, const char * Title,int N_bins, double Min, double Max){
	//std::lock_guard<std::mutex> guard(mutex); //auto lock thread
	n_bins = N_bins;
	max = Max; 
	min = Min;
	strcpy_s(title , 256, Title);
	strcpy_s(name , 256, Name);
	//overflow =0;
	//underflow =0;

	return;
 }
axis::~axis(){
	//delete(title);
	//delete(name);
}

inline int axis::get_bin_address(double x){
	//std::lock_guard<std::mutex> guard(mutex); //auto lock thread
	if( _finite(x) != 0.0 ){ 
		if(x>=max){
			return -1;
		}
		if(x<min){
			return -10;
		}
		int bin = int (n_bins*(x-min)/(max-min) );
		if(bin>=n_bins)
			return -1;
		return bin;
	}else{
		return -1; 
	}
}
////////////////////////////////////////////  end  ////////////////////////////////////////////


////////////////////////////////////////////  H1d  ////////////////////////////////////////////

H1d::H1d(){};

H1d::H1d(const char * Name, const char * Title,int N_bins_, double Min, double Max, const char * X_axis_label, const char * Dir){
	//cout << "passed in " << Dir << endl;	
	strcpy_s(title , 256, Title);
	strcpy_s(name , 256, Name);
	strcpy_s(dir , 256, Dir);

	N_bins = N_bins_;

	//printf( "dir=%s, Dir=%s \n", dir, Dir);
	Xaxis = new axis("Xaxis", X_axis_label, N_bins, Min, Max );

	bins = new double[N_bins];
	memset(bins, 0, N_bins * sizeof(double));

//	vector<double> v1DVector(N_bins, 0.0);
//	bins=v1DVector;
	
	Entries=0;
	overflow_x			=0;
	underflow_x			=0;
	
}

H1d::~H1d(){
	if (bins) {
		delete[] bins;
		bins = 0;
	}
	delete(Xaxis);
	//delete(title);
	//delete(name);
	//delete(dir);
}

//void H1d::fill(double x){
//	int bin_id = Xaxis->get_bin_address(x);
//	++Entries;
//	if (bin_id !=-1) ++bins[bin_id];
//}

void H1d::weighted_fill(double x, double weight){
	int bin_id = Xaxis->get_bin_address(x);

/*	if (bin_id >= N_bins) return;
	if (bin_id < 0) return;
*/
	++Entries;
	if(bin_id ==-1){
		overflow_x++;
	}else if(bin_id ==-10){
		underflow_x++;
	}else{
		bins[bin_id]+=weight;
	}
}

void H1d::print_bin_contents(){
	//std::lock_guard<std::mutex> guard(mutex); //auto lock thread
	printf("\n");
	//printf("\nnumber of bins=%i\n",bins.size());
	double bin_width = (Xaxis->max-Xaxis->min)/Xaxis->n_bins;
	//printf("bin widh=%f\n", bin_width);
	for (int i=0; i<(Xaxis->n_bins); i++)
	{
		if (bins[i]!=0)	Red(true);
			printf("%2.f=<(bin%i)<%2.f:%2.f   ",(double)i*bin_width+Xaxis->min,  i  ,((double)i+1.)*bin_width+Xaxis->min,bins[i]);
		if (bins[i]!=0)	White(false);
	
		if((i+1)%10==0 ) printf("\n");
		
	}
	
	printf("\n");
	return;
}

//bool H1d::match(string NAME, string TITLE, int N_BINS, double MIN, double MAX, string X_LABEL, string DIR){
//	return NAME==this->name && TITLE==this->title && N_BINS==this->Xaxis->n_bins && MIN==this->Xaxis->min && MAX==this->Xaxis->max && X_LABEL==this->Xaxis->title && DIR==this->dir; 
//}

void H1d::print_info(){
	cout<<name<< ", "<< title<< ", "<< Xaxis->n_bins << ", "<< Xaxis->min<< ", "<< Xaxis->max<< ", "<< Xaxis->title<< ", "<< dir <<endl;

}
////////////////////////////////////////////  H1d end  ////////////////////////////////////////////

////////////////////////////////////////////  H2d  ////////////////////////////////////////////

H2d::H2d(){};

H2d::H2d(const char * Name, const char * Title,int X_N_bins_, double X_Min, double X_Max, const char * X_axis_label,int Y_N_bins_, double Y_Min, double Y_Max, const char * Y_axis_label, const char * Dir){
		
	strcpy_s(title , 256, Title);
	strcpy_s(name , 256, Name);
	strcpy_s(dir , 256, Dir);

	X_N_bins = X_N_bins_;
	Y_N_bins = Y_N_bins_;

	Xaxis = new axis("Xaxis", X_axis_label, X_N_bins, X_Min, X_Max );
	Yaxis = new axis("Yaxis", Y_axis_label, Y_N_bins, Y_Min, Y_Max );

	Entries=0;
	
	overflow_x			=0;
	overflow_y			=0;
	overflow_xy			=0;
	underflow_x			=0;
	underflow_y			=0;
	underflow_xy		=0;
	flow_under_x_over_y	=0;
	flow_over_x_under_y =0;

	bins = new double*[X_N_bins];
	for (int i=0; i<X_N_bins; i++) {
			bins[i] = new double[Y_N_bins];
			memset(bins[i], 0, Y_N_bins * sizeof(double));
	}

	//initialise all bins to 0
	//vector<vector<double>> v2DVector(X_N_bins, vector<double>(Y_N_bins,0.0));
	//bins=v2DVector;

}

H2d::~H2d(){
	if (bins) {
		for (int i=0; i<X_N_bins; i++) {
			if (bins[i]) {
				delete[] bins[i];
				bins[i] = 0;
			}
		}
		delete[] bins;
		bins = 0;
	}
	delete(Xaxis);
	delete(Yaxis);
	//delete(title);
	//delete(name);
	//delete(dir);
}

//void H2d::fill(double x, double y){
//	int bin_id_x = Xaxis->get_bin_address(x);
//	int bin_id_y = Yaxis->get_bin_address(y);
//	++Entries;
//	// add one to the bin if it isn't overflow are underflow
//	if (bin_id_x !=-1 && bin_id_y !=-1) ++bins[bin_id_x][bin_id_y];
//	return;
//}

void H2d::weighted_fill(double x, double y, double weight){
	int bin_id_x = Xaxis->get_bin_address(x);
	int bin_id_y = Yaxis->get_bin_address(y);

/*	if (bin_id_x >= X_N_bins) return;
	if (bin_id_y >= Y_N_bins) return;

	if (bin_id_x < 0) return;
	if (bin_id_y < 0) return;
*/
	++Entries;

	if(bin_id_x >= 0 && bin_id_y >= 0){
		bins[bin_id_x][bin_id_y]+=weight;
	}
	else if(bin_id_x ==-1 && bin_id_y >= 0){
		overflow_x++;
	}
	else if(bin_id_x ==-10 && bin_id_y >= 0){
		underflow_x++;
	}
	else if(bin_id_y ==-1 && bin_id_x >= 0){
		overflow_y++;
	}
	else if(bin_id_y ==-10 && bin_id_x >= 0){
		underflow_y++;
	}
	
	else if(bin_id_x ==-1 && bin_id_y == -1){
		overflow_xy++;
	}
	else if(bin_id_x ==-10 && bin_id_y == -10){
		underflow_xy++;
	}

	else if(bin_id_x ==-10 && bin_id_y == -1){
		flow_under_x_over_y++;
	}
	else if(bin_id_x ==-1 && bin_id_y == -10){
		flow_over_x_under_y++;
	}

	return;
}
//void H2d::print_bin_contents(){
//	//std::lock_guard<std::mutex> guard(mutex); //auto lock thread
//	printf("\n");
//	//printf("\nnumber of bins=%i\n",bins.size());
//	double bin_width = (Xaxis->max-Xaxis->min)/Xaxis->n_bins;
//	//printf("bin widh=%f\n", bin_width);
//	for (int i=0; i<((signed int)bins.size()); i++)
//	{
//		if (bins[i]!=0)	Red(true);
//			printf("%2.f=<(bin%i)<%2.f:%i   ",(double)i*bin_width+Xaxis->min,i,((double)i+1.)*bin_width+Xaxis->min,bins[i]);
//		if (bins[i]!=0)	White(false);
//	
//		if((i+1)%10==0 ) printf("\n");
//		
//	}
//	
//	printf("\n");
//	return;
//}

//bool H2d::match(string NAME, string TITLE, int X_N_BINS, double X_MIN, double X_MAX, string X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, string Y_LABEL,  string DIR){
//	return NAME==this->name && TITLE==this->title && X_N_BINS==this->Xaxis->n_bins && X_MIN==this->Xaxis->min && X_MAX==this->Xaxis->max && X_LABEL==this->Xaxis->title && Y_N_BINS==this->Yaxis->n_bins && Y_MIN==this->Yaxis->min && Y_MAX==this->Yaxis->max && Y_LABEL==this->Yaxis->title && DIR==this->dir; 
//}

void H2d::print_info(){
	cout<<name<< ", "<< title<< ", "<< Xaxis->n_bins << ", "<< Xaxis->min<< ", "<< Xaxis->max<< ", "<< Xaxis->title<< ", ";
	cout<< Yaxis->n_bins << ", "<< Yaxis->min<< ", "<< Yaxis->max<< ", "<< Yaxis->title<< ", "<< dir <<endl;
}
////////////////////////////////////////////  H2d end  ////////////////////////////////////////////




////////////////////////////////////////////  H3d  ////////////////////////////////////////////

H3d::H3d(){};

H3d::H3d(const char * Name, const char * Title,int X_N_bins_, double X_Min, double X_Max, const char * X_axis_label,int Y_N_bins_, double Y_Min, double Y_Max, const char * Y_axis_label, int Z_N_bins_, double Z_Min, double Z_Max, const char * Z_axis_label, const char * Dir){
		
	strcpy_s(title , 256, Title);
	strcpy_s(name , 256, Name);
	strcpy_s(dir , 256, Dir);

	X_N_bins = X_N_bins_;
	Y_N_bins = Y_N_bins_;
	Z_N_bins = Z_N_bins_;

	Xaxis = new axis("Xaxis", X_axis_label, X_N_bins, X_Min, X_Max );
	Yaxis = new axis("Yaxis", Y_axis_label, Y_N_bins, Y_Min, Y_Max );
	Zaxis = new axis("Zaxis", Z_axis_label, Z_N_bins, Z_Min, Z_Max );

	Entries=0;


	vector<vector<vector<double>>> v3DVector(X_N_bins, vector<vector<double>>(Y_N_bins, vector<double>(Z_N_bins, 0.0 )));
	bins=v3DVector;
}

H3d::~H3d(){
	bins.clear();
	bins.shrink_to_fit();
	delete(Xaxis);
	delete(Yaxis);
	delete(Zaxis);
}

//void H3d::fill(double x, double y, double z){
//	int bin_id_x = Xaxis->get_bin_address(x);
//	int bin_id_y = Yaxis->get_bin_address(y);
//	int bin_id_z = Zaxis->get_bin_address(z);
//	++Entries;
//	// add one to the bin if it isn't overflow are underflow
//	if (bin_id_x !=-1 && bin_id_y !=-1 && bin_id_z !=-1) 
//		bins[bin_id_x][bin_id_y][bin_id_z] = ++bins[bin_id_x][bin_id_y][bin_id_z];
//
//	return;
//}
void H3d::weighted_fill(double x, double y, double z, double weight){
	int bin_id_x = Xaxis->get_bin_address(x);
	int bin_id_y = Yaxis->get_bin_address(y);
	int bin_id_z = Zaxis->get_bin_address(z);

	if (bin_id_x >= X_N_bins) return;
	if (bin_id_y >= Y_N_bins) return;
	if (bin_id_z >= Z_N_bins) return;

	if (bin_id_x < 0) return;
	if (bin_id_y < 0) return;
	if (bin_id_z < 0) return;

	++Entries;

	if(bin_id_x >= 0 && bin_id_y >= 0 && bin_id_z >= 0){
		bins[bin_id_x][bin_id_y][bin_id_z]=bins[bin_id_x][bin_id_y][bin_id_z]+weight;
	}

	/////// over and underflow are not shown in root so I didn't bother to finish this part //////////
	//else if(bin_id_x ==-1 && bin_id_y >= 0 && bin_id_z >= 0){
	//	overflow_x++;
	//}
	//else if(bin_id_y ==-1 && bin_id_x >= 0 && bin_id_z >= 0){
	//	overflow_y++;
	//}
	//else if(bin_id_y >= 0 && bin_id_x >= 0 && bin_id_z ==-1){
	//	overflow_z++;
	//}
	//
	//else if(bin_id_x ==-1 && bin_id_y == -1 && bin_id_z >= 0){
	//	overflow_xy++;
	//}
	//else if(bin_id_x ==-1 && bin_id_y >= 0 && bin_id_z == -1){
	//	overflow_xz++;
	//}	
	//else if(bin_id_x >= 0 && bin_id_y == -1 && bin_id_z == -1){
	//	overflow_yz++;
	//}	
	//else if(bin_id_x == -1 && bin_id_y == -1 && bin_id_z == -1){
	//	overflow_xyz++;
	//}	
	//
	//
	//
	//
	//else if(bin_id_x ==-10 && bin_id_y >= 0 && bin_id_z >= 0){
	//	underflow_x++;
	//}
	//else if(bin_id_y ==-10 && bin_id_x >= 0 && bin_id_z >= 0){
	//	underflow_y++;
	//}
	//	else if(bin_id_x ==-10 && bin_id_y == -10 && bin_id_z >= 0){
	//	underflow_xy++;
	//}
	//else if(bin_id_x ==-10 && bin_id_y >= 0 && bin_id_z == -10){
	//	underflow_xz++;
	//}		
	//else if(bin_id_x >= 0 && bin_id_y == -10 && bin_id_z == -10){
	//	underflow_yz++;
	//}	
	//else if(bin_id_x == -10 && bin_id_y == -10 && bin_id_z == -10){
	//	underflow_xyz++;
	//}
	//
	//
	//
	//else if(bin_id_x ==-10 && bin_id_y == -1 ){
	//	flow_over_y_under_x++;
	//}
	//else if(bin_id_x ==-1 && bin_id_y == -10 ){
	//	flow_under_x_over_y++;
	//}




	return;
}

void H3d::print_info(){
	cout<<name<< ", "<< title<< ", "<< Xaxis->n_bins << ", "<< Xaxis->min<< ", "<< Xaxis->max<< ", "<< Xaxis->title<< ", ";
	cout<< Yaxis->n_bins << ", "<< Yaxis->min<< ", "<< Yaxis->max<< ", "<< Yaxis->title<< ", "<< dir <<endl;
}
////////////////////////////////////////////  H2d end  ////////////////////////////////////////////


////////////////////////////////////////////  histo_container  ////////////////////////////////////////////

histo_handler::histo_handler(){
	hist_check_counter=0;

	for(int i=0; i<MAX_NUM_HIST; ++i){
		H1d_vector.push_back(0);
		H2d_vector.push_back(0);
		H3d_vector.push_back(0);
	}
}

histo_handler::~histo_handler(){
	for(int i=0; i<MAX_NUM_HIST; ++i){
		if(H1d_vector[i]) delete(H1d_vector[i]);
		if(H2d_vector[i]) delete(H2d_vector[i]);
		if(H3d_vector[i]) delete(H3d_vector[i]);
	}
}


inline bool compare(const char* s1,const char* s2){
	for(int i=0; i<50; i++){
		if ( *(s1+i) != *(s2+i) )	{
		//	cout << "return false"<< endl;
			return false;
		}
		if ( *(s1+i) == '\0') {
		//	cout << "return true"<< endl;	
			return true;
		}
		//cout << *(s1+i) << ", "<< *(s2+i) << "i=" << i << endl;
	}
}

//void histo_handler::fill1(int id, const char * NAME, double x, const char * TITLE, int N_BINS, double MIN, double MAX, const char * X_LABEL, const char * DIR){
//	if( id > MAX_NUM_HIST) {
//		int max_num_hist = MAX_NUM_HIST;
//		printf("The current max number of histograms are %d.  Please change the max (in Simple_Hist.cpp -> histo_handler::histo_handler()).\n", max_num_hist);
//		return;
//	}
//
//	//char* temp;
//	////strcpy(temp, DIR);
//	////printf("id=%d", id);
//	////
//	//unsigned int hash = 0;
//	//int c;
//
//	////while (c = *DIR++)
//	////	hash += c;
//	//id = id + hash;
//	//printf(", DIR=%s, hash=%d, new id=%d \n", DIR, hash, id);
//
//
//	if( H1d_vector[id] == 0){
//		H1d * hist = new H1d(NAME, TITLE, N_BINS, MIN, MAX, X_LABEL, DIR);
//		H1d_vector[id]=hist;
//		H1d_vector[id]->fill(x);
//	}else{
//		//checking to see if the histogram matches is slow so we will only do this in the very begining 
//		if ( strcmp( H1d_vector[id]->get_name(), NAME)!=0 || H1d_vector[id]->get_X_n_bins() != N_BINS || H1d_vector[id]->get_X_max() != MAX || H1d_vector[id]->get_X_min() != MIN)// !compare( H1d_vector[id]->get_name(),NAME) strcmp( H1d_vector[id]->get_name(), NAME)!=0 )
//		{
//			Red(true);
//			cout<<endl;
//			cout<<"ERROR: Histogram is:" << NAME<< ", "<< TITLE<< ", "<< N_BINS<< ", "<< MIN<< ", "<< MAX<< ", "<< X_LABEL<< ", "<< DIR<<endl;
//			cout<<"And it should be   :";
//			H1d_vector[id]->print_info();
//			cout<<endl;
//			White(false);
//		}
//		else {
//			H1d_vector[id]->fill(x);
//			hist_check_counter++;
//		}
//	}
//	return;
//
//}



void histo_handler::fill1(int id, const char * NAME, double x,  double weight, const char * TITLE, int N_BINS, double MIN, double MAX, const char * X_LABEL, const char * DIR){
	if( id > MAX_NUM_HIST) {
		int max_num_hist = MAX_NUM_HIST;
		printf("The current max number of histograms are %d.  Please change the max (in Simple_Hist.cpp -> histo_handler::histo_handler()).\n", max_num_hist);
		return;
	}

	//printf("id=%d", id);
	//cout << "dir=" << DIR ;
	//printf(", DIR=%s , ", DIR);
	//cout << "dir=" << DIR << endl;
	//char temp[256];
	//strcpy(temp,DIR);
	////sprintf(temp, "%s\0",DIR);
	//unsigned int hash = 0;
	//int c=0;
	//int b=0;
	//for(int i=1;i<256;i++){	
	//	//cout<<temp[i];
	//	c=temp[i];
	//	hash = hash/33 + c;
	//	if(temp[i+1] == '\0')
	//		i=1000;
	//	b=temp[i-1];
	//}
	////cout << endl;
	//id = id + hash;
	//printf(", DIR=%s, hash=%d, new id=%d \n", DIR, hash, id);
	//cout << "dir=" << DIR << endl;
	//cout << "temp=" << temp << endl;
	if( H1d_vector[id] == 0){
		H1d * hist = new H1d(NAME, TITLE, N_BINS, MIN, MAX, X_LABEL, DIR);
		H1d_vector[id]=hist;
		H1d_vector[id]->weighted_fill(x,weight);
	}else{
		//	printf(", DIR=%s \n", H1d_vector[id]->get_dir());
		//cout<<"Histogram is "<< id << NAME << ", "<< DIR<< " and " << H1d_vector[id]->get_name() << ", " << H1d_vector[id]->get_dir() << endl;
	
		//checking to see if the histogram matches is slow so we will only do this in the very begining 
		if ( strcmp( H1d_vector[id]->get_dir(), DIR)!=0 || strcmp( H1d_vector[id]->get_name(), NAME)!=0 || H1d_vector[id]->get_X_n_bins() != N_BINS || H1d_vector[id]->get_X_max() != MAX || H1d_vector[id]->get_X_min() != MIN)// !compare( H1d_vector[id]->get_name(),NAME) strcmp( H1d_vector[id]->get_name(), NAME)!=0 )
		{
			Red(true);
			cout<<endl;
			cout<<"ERROR: Histogram is:" << NAME<< ", "<< TITLE<< ", "<< N_BINS<< ", "<< MIN<< ", "<< MAX<< ", "<< X_LABEL<< ", "<< DIR<<endl;
			cout<<"And it should be   :";
			H1d_vector[id]->print_info();
			cout<<endl;
			White(false);
		}
		else {
			H1d_vector[id]->weighted_fill(x,weight);
			hist_check_counter++;
		}
	}
	return;

}

void histo_handler::fill1(int id, double x, double weight){
	H1d_vector[id]->weighted_fill(x,weight);
}


void histo_handler::combine_hist(H1d * hist2, int pos){
	
	if(hist2 == 0) return;

	if(H1d_vector[pos] == 0){ //if it dosen't exit just add it
		H1d_vector[pos]=hist2;
	}
	else{
		for(int i=0; i < hist2->Xaxis->n_bins; ++i){
			H1d_vector[pos]->bins[i] = H1d_vector[pos]->bins[i] + hist2->bins[i];
		}
		H1d_vector[pos]->set_X_overflow(  H1d_vector[pos]->get_X_overflow()  +	hist2->get_X_overflow() );
		H1d_vector[pos]->set_X_underflow( H1d_vector[pos]->get_X_underflow() +	hist2->get_X_underflow() );
		
		H1d_vector[pos]->set_Entries(	H1d_vector[pos]->get_Entries() + hist2->get_Entries());
	}
}


	
//void histo_handler::fill2(int id, const char * NAME, double x, double y, const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, const char * DIR){
//
//	if( id > MAX_NUM_HIST) {
//		int max_num_hist = MAX_NUM_HIST;
//		printf("The current max number of histograms are %d.  Please change the MAX_NUM_HIST (in Simple_Hist.cpp -> histo_handler::histo_handler()).\n", max_num_hist);
//		return;
//	}
//
//	//printf("id=%d", id);
//	//cout << DIR ;
//
//	//unsigned int hash = 0;
//	//int c;
//
//	////while (c = *DIR++)
//	////	hash += c;
//	//id = id + hash;
//	//printf(", DIR=%s, hash=%d, new id=%d \n", DIR, hash, id);
//
//	if( H2d_vector[id] == 0){
//		H2d * hist = new H2d(NAME, TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, DIR);
//		H2d_vector[id]=hist;
//		H2d_vector[id]->fill(x,y);
//	}else{
//		//checking to see if the histogram matches is slow so we will only do this in the very begining 
//		if ( strcmp( H2d_vector[id]->get_name(), NAME)!=0         || H2d_vector[id]->get_X_n_bins() != X_N_BINS || H2d_vector[id]->get_X_max() != X_MAX || H2d_vector[id]->get_X_min() != X_MIN      || H2d_vector[id]->get_Y_n_bins() != Y_N_BINS || H2d_vector[id]->get_Y_max() != Y_MAX || H2d_vector[id]->get_Y_min() != Y_MIN)// !compare( H2d_vector[id]->get_name(),NAME) strcmp( H2d_vector[id]->get_name(), NAME)!=0 )
//		{
//			Red(true);
//			cout<<endl;
//			cout<<"ERROR: Histogram is:"<<NAME<< ", "<< TITLE<< ", "<< X_N_BINS<< ", "<< X_MIN<< ", "<< X_MAX<< ", "<< X_LABEL<< ", "<< DIR<<endl;
//			cout<<"And it should be   :";
//			H2d_vector[id]->print_info();
//			cout<<endl;
//			White(false);
//		}
//		else {
//			H2d_vector[id]->fill(x,y);
//			hist_check_counter++;
//		}
//	}
//	return;
//
//}


void histo_handler::fill2(int id, const char * NAME, double x, double y, double weight,const  char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, const char * DIR){

	if( id > MAX_NUM_HIST) {
		int max_num_hist = MAX_NUM_HIST;
		printf("The current max number of histograms are %d.  Please change the max (in Simple_Hist.cpp -> histo_handler::histo_handler()).\n", max_num_hist);
		return;
	}

	//printf("id=%d", id);
	
	//unsigned int hash = 0;
	//int c;

	////while (c = *DIR++)
	////	hash += c;
	//id = id + hash;
	
	//printf(", DIR=%s, hash=%d, new id=%d \n", DIR, hash, id);


	if( H2d_vector[id] == 0){
		H2d * hist = new H2d(NAME, TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, DIR);
		H2d_vector[id]=hist;
		H2d_vector[id]->weighted_fill(x,y,weight);
	}else{
		//checking to see if the histogram matches is slow so we will only do this in the very begining 
		if ( strcmp( H2d_vector[id]->get_name(), NAME)!=0         || H2d_vector[id]->get_X_n_bins() != X_N_BINS || H2d_vector[id]->get_X_max() != X_MAX || H2d_vector[id]->get_X_min() != X_MIN      || H2d_vector[id]->get_Y_n_bins() != Y_N_BINS || H2d_vector[id]->get_Y_max() != Y_MAX || H2d_vector[id]->get_Y_min() != Y_MIN)// !compare( H2d_vector[id]->get_name(),NAME) strcmp( H2d_vector[id]->get_name(), NAME)!=0 )
		{
			Red(true);
			cout<<endl;
			cout<<"ERROR: Histogram is:"<<NAME<< ", "<< TITLE<< ", "<< X_N_BINS<< ", "<< X_MIN<< ", "<< X_MAX<< ", "<< X_LABEL<< ", "<< DIR<<endl;
			cout<<"And it should be   :";
			H2d_vector[id]->print_info();
			cout<<endl;
			White(false);
		}
		else {
			H2d_vector[id]->weighted_fill(x,y,weight);
			hist_check_counter++;
		}
	}
	return;

}




void histo_handler::combine_hist(H2d * hist2, int pos){
	if(hist2 == 0) return;
	
	if(H2d_vector[pos] == 0){ //if it dosen't exit just add it
		H2d_vector[pos]=hist2;
	}
	else{
		for(int i=0; i < (int)hist2->get_X_n_bins(); ++i){
			for(int j=0; j < (int)hist2->get_Y_n_bins(); ++j){
				H2d_vector[pos]->bins[i][j] = H2d_vector[pos]->bins[i][j] + hist2->bins[i][j];
			}
		}

		H2d_vector[pos]->set_X_overflow(  H2d_vector[pos]->get_X_overflow()  +	hist2->get_X_overflow() );
		H2d_vector[pos]->set_X_underflow( H2d_vector[pos]->get_X_underflow() +	hist2->get_X_underflow() );
		
		H2d_vector[pos]->set_Y_overflow(  H2d_vector[pos]->get_Y_overflow()  +	hist2->get_Y_overflow() );
		H2d_vector[pos]->set_Y_underflow( H2d_vector[pos]->get_Y_underflow() +	hist2->get_Y_underflow() );

		std::cout << "fix over and underflows" << endl;
		
		H2d_vector[pos]->set_Entries(	H2d_vector[pos]->get_Entries() + hist2->get_Entries());

	}
}	

void histo_handler::fill3(int id, const char * NAME, double x, double y, double z, double weight,const  char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, int Z_N_BINS, double Z_MIN, double Z_MAX, const char * Z_LABEL, const char * DIR){

	if( id > MAX_NUM_HIST) {
		int max_num_hist = MAX_NUM_HIST;
		printf("The current max number of histograms are %d.  Please change the max (in Simple_Hist.cpp -> histo_handler::histo_handler()).\n", max_num_hist);
		return;
	}


	if( H3d_vector[id] == 0){
		H3d * hist = new H3d(NAME, TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, Z_N_BINS, Z_MIN, Z_MAX, Z_LABEL, DIR);
		H3d_vector[id]=hist;
		H3d_vector[id]->weighted_fill(x,y,z, weight);
		
		//printf("Hist (%s) is %d bytes\n", NAME, sizeof(*hist) );

	}else{
		if ( strcmp( H3d_vector[id]->get_name(), NAME)!=0   || H3d_vector[id]->get_X_n_bins() != X_N_BINS || H3d_vector[id]->get_X_max() != X_MAX || H3d_vector[id]->get_X_min() != X_MIN      || H3d_vector[id]->get_Y_n_bins() != Y_N_BINS || H3d_vector[id]->get_Y_max() != Y_MAX || H3d_vector[id]->get_Y_min() != Y_MIN)// !compare( H3d_vector[id]->get_name(),NAME) strcmp( H3d_vector[id]->get_name(), NAME)!=0 )
		{
			Red(true);
			cout<<endl;
			cout<<"ERROR: Histogram is:"<<NAME<< ", "<< TITLE<< ", "<< X_N_BINS<< ", "<< X_MIN<< ", "<< X_MAX<< ", "<< X_LABEL<< ", "<< DIR<<endl;
			cout<<"And it should be   :";
			H3d_vector[id]->print_info();
			cout<<endl;
			White(false);
		}
		else {
			H3d_vector[id]->weighted_fill(x,y,z,weight);
			hist_check_counter++;
		}
	}
	//return;
}				


void histo_handler::combine_hist(H3d * hist2, int pos){
	if(hist2 == 0) return;
	
	if(H3d_vector[pos] == 0){ //if it dosen't exit just add it
		H3d_vector[pos]=hist2;
	}
	else{
		for(int i=0; i < (int)hist2->get_X_n_bins(); ++i){
			for(int j=0; j < (int)hist2->get_Y_n_bins(); ++j){
				for(int k=0; k < (int)hist2->get_Z_n_bins(); ++k){
					H3d_vector[pos]->bins[i][j][k] = H3d_vector[pos]->bins[i][j][k] + hist2->bins[i][j][k];
				}
			}
		}
//		H3d_vector[pos]->set_X_overflow(  H3d_vector[pos]->get_X_overflow()  +	hist2->get_X_overflow() );
//		H3d_vector[pos]->set_X_underflow( H3d_vector[pos]->get_X_underflow() +	hist2->get_X_underflow() );
//		
//		H3d_vector[pos]->set_Y_overflow(  H3d_vector[pos]->get_Y_overflow()  +	hist2->get_Y_overflow() );
//		H3d_vector[pos]->set_Y_underflow( H3d_vector[pos]->get_Y_underflow() +	hist2->get_Y_underflow() );
//		
//		H3d_vector[pos]->set_Z_overflow(  H3d_vector[pos]->get_Z_overflow()  +	hist2->get_Z_overflow() );
//		H3d_vector[pos]->set_Z_underflow( H3d_vector[pos]->get_Z_underflow() +	hist2->get_Z_underflow() );
			
		H3d_vector[pos]->set_Entries(	H3d_vector[pos]->get_Entries() + hist2->get_Entries());

		}
}
}
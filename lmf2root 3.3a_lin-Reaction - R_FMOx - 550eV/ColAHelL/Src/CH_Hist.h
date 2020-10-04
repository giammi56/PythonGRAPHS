#ifndef SIMPLE_HIST_ALREADY_INCLUDED
	#define SIMPLE_HIST_ALREADY_INCLUDED
#include <string>
#include <vector>
#include <unordered_map> 

//#include <thread>
//#include <mutex>
//#include <functional> // for std::hash<int>

namespace CH
{
	using namespace std;


	#define MAX_NUM_HIST 100000


	class axis
	{
	//	std::mutex mutex;


	public:
		int n_bins;

		double max, min;
		char title[256];
		char name[256];
		//string axis_label;

		axis(const char * Name, const char * Title,int N_bins, double Min, double Max);
		~axis();

		int get_bin_address(double x);

		//vector * get_bins();
		//void print_bin_contents();

		//inline __int32 get_overflow(){return overflow;}
		//inline __int32 get_underflow(){return underflow;}
		//inline void set_overflow(__int32 value){ overflow = value;}
		//inline void set_underflow(__int32 value){ underflow = value;}

	};



	/////////////////////////////////////////////////////
	//1d histogram of doubles
	class H1d
	{

	
		char title[256];
		char name[256];
		char dir[256];

		int Entries;

		int overflow_x,  underflow_x;

		int N_bins;

	public:
		H1d();
		H1d(const char * Name, const char * Title,int N_bins, double Min, double Max, const char * X_axis_label, const char * Dir);
		~H1d();

		//char * key;

		//bool operator == (const H1d& rhs) const {
		//	return this->name == rhs.name && rhs.title==this->title && rhs.Xaxis->n_bins==this->Xaxis->n_bins && rhs.Xaxis->min==this->Xaxis->min && rhs.Xaxis->max==this->Xaxis->max && rhs.Xaxis->title==this->Xaxis->title && rhs.dir==this->dir;
		//} 
		axis * Xaxis;
		double * bins;
		
		//vector<double> bins;
		void print_bin_contents();
	//	void fill(double x);
		void weighted_fill(double x, double weight);
	//	bool match(string NAME, string TITLE, int N_BINS, double MIN, double MAX, string X_LABEL, string DIR);
		void print_info();

		inline const char * get_name(){return name;}
		inline const char * get_dir(){return dir;}	
		inline const char * get_title(){return title;}	
		inline const char * get_X_title(){return Xaxis->title;}	

		inline double get_X_max(){return Xaxis->max;}	
		inline double get_X_min(){return Xaxis->min;}	
		inline int	  get_X_n_bins(){return Xaxis->n_bins;}

		inline __int32 get_X_overflow(){return		overflow_x;}
		inline __int32 get_X_underflow(){return		underflow_x;}
		inline void set_X_overflow(__int32 value){	overflow_x = value;}
		inline void set_X_underflow(__int32 value){ underflow_x = value;}

		inline void set_Entries(int value){ Entries = value;}
		inline int get_Entries(){ return Entries ;}

	};


	/////////////////////////////////////////////////////
	//2d histogram of doubles
	class H2d 
	{
		axis * Xaxis;
		axis * Yaxis;
		char title[256];
		char name[256];
		char dir[256];

		int	X_N_bins, Y_N_bins;
		int Entries;
		int overflow_x , overflow_y, overflow_xy;
		int underflow_x, underflow_y, underflow_xy;
		int flow_under_x_over_y, flow_over_x_under_y; 


	public:
		H2d();
		H2d(const char * Name, const char * Title,int X_N_bins, double X_Min, double X_Max, const char * X_axis_label,int Y_N_bins, double Y_Min, double Y_Max, const char * Y_axis_label, const char * Dir);
		~H2d();

		//char * key;
		//vector<vector<double>> bins;
	
		double ** bins;

		//void print_bin_contents();
	//	void fill(double x, double y);
		void weighted_fill(double x, double y, double weight);
		//bool match(string NAME, string TITLE, int X_N_BINS, double X_MIN, double X_MAX, string X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, string Y_LABEL,  string DIR);
		void print_info();


		inline const char * get_name(){return				name;}
		inline const char * get_dir(){return				dir;}	
		inline const char * get_title(){return			title;}	

		inline const char * get_X_title(){return			Xaxis->title;}	


		inline double get_X_max(){return			Xaxis->max;}	
		inline double get_X_min(){return			Xaxis->min;}	
		inline int	  get_X_n_bins(){return			Xaxis->n_bins;}

		inline __int32 get_X_overflow(){return		overflow_x;}
		inline __int32 get_X_underflow(){return		underflow_x;}
		inline void set_X_overflow(__int32 value){	overflow_x = value;}
		inline void set_X_underflow(__int32 value){ underflow_x = value;}

		inline const char * get_Y_title(){return			Yaxis->title;}	

		inline double get_Y_max(){return			Yaxis->max;}	
		inline double get_Y_min(){return			Yaxis->min;}	
		inline int	  get_Y_n_bins(){return			Yaxis->n_bins;}

		inline __int32 get_Y_overflow(){return		overflow_y;}
		inline __int32 get_Y_underflow(){return		underflow_y;}
		inline void set_Y_overflow(__int32 value){	overflow_y = value;}
		inline void set_Y_underflow(__int32 value){ underflow_y = value;}

		inline __int32 get_overflow_xy(){return		overflow_xy;}
		inline __int32 get_underflow_xy(){return	underflow_xy;}

		inline void set_overflow_xy	(int value)	{overflow_xy  = value;}
		inline void set_underflow_xy(int value)	{underflow_xy = value;}

		inline __int32 get_flow_under_x_over_y(){return		flow_under_x_over_y;}
		inline __int32 get_flow_over_x_under_y(){return		flow_over_x_under_y;}
	
		inline void set_flow_under_x_over_y(int value){flow_under_x_over_y = value;}
		inline void set_flow_over_x_under_y(int value){flow_over_x_under_y = value;}
	

		inline void set_Entries(int value){ Entries = value;}
		inline int get_Entries(){ return Entries ;}

	};


	/////////////////////////////////////////////////////
	//3d histogram of doubles
	class H3d 
	{
		axis * Xaxis;
		axis * Yaxis;
		axis * Zaxis;
	
		char title[256];
		char name[256];
		char dir[256];

		int Entries;

		int X_N_bins, Y_N_bins, Z_N_bins;

		//int overflow_x , overflow_y, overflow_z, overflow_xy, overflow_xz, overflow_yz, overflow_xyz;
		//int underflow_x, underflow_y, underflow_z, underflow_xy, underflow_xz, underflow_yz, underflow_xyz;
		//int flow_under_x_over_y, flow_over_y_under_x; 

	public:
		H3d();
		H3d(const char * Name, const char * Title,int X_N_bins, double X_Min, double X_Max, const char * X_axis_label,int Y_N_bins, double Y_Min, double Y_Max, const char * Y_axis_label,int Z_N_bins, double Z_Min, double Z_Max, const char * Z_axis_label, const char * Dir);
		~H3d();

		vector<vector<vector<double>>> bins;
	
		//void fill(double x, double y, double z);
		void weighted_fill(double x, double y, double z, double weight);
		//bool match(string NAME, string TITLE, int X_N_BINS, double X_MIN, double X_MAX, string X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, string Y_LABEL,  string DIR);
		void print_info();


		inline const char * get_name(){return				name;}
		inline const char * get_dir(){return				dir;}	
		inline const char * get_title(){return			title;}	

		inline const char * get_X_title(){return			Xaxis->title;}	


		inline double get_X_max(){return			Xaxis->max;}	
		inline double get_X_min(){return			Xaxis->min;}	
		inline int	  get_X_n_bins(){return			Xaxis->n_bins;}

		//inline __int32 get_X_overflow(){return		Xaxis->overflow;}
		//inline __int32 get_X_underflow(){return		Xaxis->underflow;}
		//inline void set_X_overflow(__int32 value){	Xaxis->overflow = value;}
		//inline void set_X_underflow(__int32 value){ Xaxis->underflow = value;}

		inline const char * get_Y_title(){return			Yaxis->title;}	

		inline double get_Y_max(){return			Yaxis->max;}	
		inline double get_Y_min(){return			Yaxis->min;}	
		inline int	  get_Y_n_bins(){return			Yaxis->n_bins;}

		//inline __int32 get_Y_overflow(){return		Yaxis->overflow;}
		//inline __int32 get_Y_underflow(){return		Yaxis->underflow;}
		//inline void set_Y_overflow(__int32 value){	Yaxis->overflow = value;}
		//inline void set_Y_underflow(__int32 value){ Yaxis->underflow = value;}

		inline const char * get_Z_title(){return				Zaxis->title;}	
		inline double		get_Z_max(){return					Zaxis->max;}	
		inline double		get_Z_min(){return					Zaxis->min;}	
		inline int			get_Z_n_bins(){return				Zaxis->n_bins;}
		//inline __int32		get_Z_overflow(){return				Zaxis->overflow;}
		//inline __int32		get_Z_underflow(){return			Zaxis->underflow;}
		//inline void			set_Z_overflow(__int32 value){		Zaxis->overflow = value;}
		//inline void			set_Z_underflow(__int32 value){		Zaxis->underflow = value;}

		inline void set_Entries(int value){ Entries = value;}
		inline int get_Entries(){ return Entries ;}
	};



	///////////////////////////////////////////////////////
	//1d histogram of doubles
	class histo_handler
	{
		string name;
//		int size_H1d;
		int hist_check_counter;
	public:
		histo_handler();
		~histo_handler();



		vector<H1d*>	H1d_vector;
		void fill1(int id, double x, double weight);
		inline void fill1(int id, const char * NAME, double x, const char * TITLE, int N_BINS, double MIN, double MAX, const char * X_LABEL, const char * DIR){
			fill1( id, NAME,  x, 1.,  TITLE,  N_BINS,  MIN,  MAX,  X_LABEL,  DIR);}
		void fill1(int id, const char * NAME, double x, double weight, const char * TITLE, int N_BINS, double MIN, double MAX, const char * X_LABEL, const char * DIR);
		// good old lmf2root style
		void fill(int id, const char * NAME, double x, double weight, const char * TITLE, int N_BINS, double MIN, double MAX, const char * X_LABEL, const char * DIR) {
			fill1(id, NAME, x, weight, TITLE, N_BINS, MIN, MAX, X_LABEL, DIR);};

		void combine_hist(H1d * hist2, int pos);



		vector<H2d*>	H2d_vector;
		void fill2(int id, const char * NAME, double x, double y,					const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, const char * DIR){
			fill2( id, NAME, x,  y, 1.,	 TITLE,  X_N_BINS,  X_MIN, X_MAX,  X_LABEL,  Y_N_BINS,  Y_MIN,  Y_MAX,  Y_LABEL, DIR);}
		void fill2(int id, const char * NAME, double x, double y, double weight,	const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, const char * DIR);
		// good old lmf2root style
		void fill(int id, const char * NAME, double x, double y, double weight,	const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, const char * DIR) {
			fill2(id, NAME, x, y, weight, TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, DIR);};
		
		void combine_hist(H2d * hist2, int pos);


		vector<H3d*>	H3d_vector;
		inline void fill3(int id, const char * NAME, double x, double y, double z, 	const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, int Z_N_BINS, double Z_MIN, double Z_MAX, const char * Z_LABEL, const char * DIR){
			fill3( id, NAME,  x,  y, z, 1.,	TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, Z_N_BINS, Z_MIN, Z_MAX, Z_LABEL, DIR);}
		void fill3(int id, const char * NAME, double x, double y, double z, double weight,	const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, int Z_N_BINS, double Z_MIN, double Z_MAX, const char * Z_LABEL, const char * DIR);
		// good old lmf2root style	
		void fill(int id, const char * NAME, double x, double y, double z, double weight,	const char * TITLE, int X_N_BINS, double X_MIN, double X_MAX, const char * X_LABEL, int Y_N_BINS, double Y_MIN, double Y_MAX, const char * Y_LABEL, int Z_N_BINS, double Z_MIN, double Z_MAX, const char * Z_LABEL, const char * DIR) {
			fill3(id, NAME, x, y, z, weight, TITLE, X_N_BINS, X_MIN, X_MAX, X_LABEL, Y_N_BINS, Y_MIN, Y_MAX, Y_LABEL, Z_N_BINS, Z_MIN, Z_MAX, Z_LABEL, DIR);};

		void combine_hist(H3d * hist3, int pos);
	};
}
#endif
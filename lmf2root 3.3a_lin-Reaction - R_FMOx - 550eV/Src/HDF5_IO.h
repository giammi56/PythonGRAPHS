#ifndef HFD5_is_included_s786dgfsd67
	#define HFD5_is_included_s786dgfsd67

	// test for Visual Studio version above 2010:
	#ifdef _MSC_VER
		#if _MSC_VER > 1600
			#define USE_HDF5_STUFF
		#endif
	#endif

	#include "LMF_IO.h"
	
	#ifdef USE_HDF5_STUFF
	#include "..\HDF5\include\hdf5.h"

	#pragma comment(lib,"libszip.lib")
	#pragma comment(lib,"libhdf5.lib")

	//#pragma comment(lib,"libzlib.lib")
	//#pragma comment(lib,"libhdf5_cpp.lib")
	//#pragma comment(lib,"libhdf5_hl.lib")
	//#pragma comment(lib,"libhdf5_hl_cpp.lib")
	//#pragma comment(lib,"libhdf5_tools.lib")
	// libszip.lib;libzlib.lib;libhdf5.lib;libhdf5_cpp.lib;libhdf5_hl.lib;libhdf5_hl_cpp.lib;libhdf5_tools.lib;
	#endif

	#define Num_HDF5_ADC_Channels 16

	class HDF5_IO
	{
	public:
		HDF5_IO(int trigger_channel, double group_range_start_ns, double group_range_end_ns, __int64 time_slice_width_ns, int ADC_binsize_ps, double full_scale_mV, bool use_grouping);
		~HDF5_IO();

		void initialize();
		void CloseInputHDF5();
		void reset();
		bool OpenInputHDF5(const char* file_name);
		bool ReadNextHDF5packet_grouping(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words);
		bool ReadNextHDF5packet_slicing(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words);
		void seek_to_begin();

		bool	eof;
		__int32 TriggerChannel;
		__int32 error_flag;
		bool	input_file_is_open;
		__int32 number_of_channels;
		__int32	max_number_of_hits;
		__int64 int64_Numberofevents;
		__int64 int64_number_of_read_events;

		double bits_per_mVolt;
		double	sample_period_ps;
		double full_scale_mV;
		bool do_grouping;
		double	group_range_start_ns;
		double	group_range_end_ns;
		__int64	group_range_start_bins;
		__int64	group_range_end_bins;
		bool	trigger_exists;
		__int64 last_trigger_train_ID;

		ndigo_packet * local_ndigo_packet;

		__int64 current_train_ID;
		__int64 first_train_ID;
		__int64 current_packet_index[Num_HDF5_ADC_Channels];
		__int64 group_slice_size_in_bins;
		__int64 trigger_bin_position;
		__int64 TrainsInFile;
		bool init_was_successful;

		#ifdef USE_HDF5_STUFF
		hid_t       HDF5_file, memtype;
		herr_t      status;
		hsize_t     dims[2] = { 2,0 };
			__int32 *ptr_32;
		int ndims;

			
		// ADC channels
		// ----------------------------------------------------
		int MaxNumberOfPackets[Num_HDF5_ADC_Channels];
		// ADC packet data
		hid_t space_ADC[Num_HDF5_ADC_Channels];
		hid_t ds_ADC[Num_HDF5_ADC_Channels];
		hvl_t *ADCdata[Num_HDF5_ADC_Channels];
		// ADC timestamps
		hid_t space_ADC_TS[Num_HDF5_ADC_Channels];
		hid_t ds_ADC_TS[Num_HDF5_ADC_Channels];
		int *ADC_TS[Num_HDF5_ADC_Channels];
		// ADC packet TrainIDs
		hid_t space_ADC_TID[Num_HDF5_ADC_Channels];
		hid_t ds_ADC_TID[Num_HDF5_ADC_Channels];
		int *ADC_TID[Num_HDF5_ADC_Channels];

		// Photon properties
		// ----------------------------------------------------
		// photons/TrainID
		hid_t ds_ph_TID;
		int ph_TID[500];
		// photons/intensityTD
		hid_t space_ph_iTD;
		hid_t ds_ph_iTD;
		float **ph_iTD;
		#else
			void * HDF5_file;
			void * memtype;
		#endif

	};
#endif
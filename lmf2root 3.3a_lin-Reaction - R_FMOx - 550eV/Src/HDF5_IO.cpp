#include "HDF5_IO.h"

#ifdef HFD5_is_included_s786dgfsd67

/////////////////////////////////////////////////////////////////
HDF5_IO::HDF5_IO(int trigger_channel, double group_range_start_ns_, double group_range_end_ns_, __int64 time_slice_width_ns, int ADC_binsize_ps, double full_scale_mV_, bool use_grouping)
/////////////////////////////////////////////////////////////////
{
	group_range_start_ns = group_range_start_ns_;
	group_range_end_ns = group_range_end_ns_;
	if (ADC_binsize_ps > 0.) group_range_start_bins = (1000*(int)group_range_start_ns_)/ADC_binsize_ps;
	if (ADC_binsize_ps > 0.) group_range_end_bins   = (1000*(int)group_range_end_ns_)  /ADC_binsize_ps;
	full_scale_mV = full_scale_mV_;
	do_grouping = use_grouping;
	TriggerChannel = trigger_channel;
	group_slice_size_in_bins = 0;
	sample_period_ps = ADC_binsize_ps;
	if (ADC_binsize_ps > 0.) group_slice_size_in_bins = (1000 * time_slice_width_ns) / ADC_binsize_ps;
	local_ndigo_packet = 0;
	input_file_is_open = false;
	HDF5_file = 0;
	memtype = 0;
	error_flag = 0;
	init_was_successful = false;
	eof = false;
	initialize();
}



/////////////////////////////////////////////////////////////////
HDF5_IO::~HDF5_IO()
/////////////////////////////////////////////////////////////////
{
	if (local_ndigo_packet) {
		delete local_ndigo_packet;
		local_ndigo_packet = 0;
	}
}


/////////////////////////////////////////////////////////////////
// This function is called once during the construction of the instance
void HDF5_IO::initialize()
/////////////////////////////////////////////////////////////////
{
	init_was_successful = true;
	#ifdef USE_HDF5_STUFF
	if (!local_ndigo_packet) local_ndigo_packet = new ndigo_packet();
		if (do_grouping && TriggerChannel < 0) init_was_successful = false;
		if (do_grouping && (group_range_start_ns > group_range_end_ns))  init_was_successful = false;
		if (sample_period_ps < 1) init_was_successful = false;
		if ((!do_grouping) && (group_slice_size_in_bins < 1)) init_was_successful = false;
		if (!init_was_successful) error_flag = 1;
	#endif
}

/////////////////////////////////////////////////////////////////
void HDF5_IO::reset()
/////////////////////////////////////////////////////////////////
{
	error_flag = 0;
	eof = false;
	number_of_channels = 0;
	max_number_of_hits = 0;
	int64_Numberofevents = 0;
	int64_number_of_read_events = 0;
	if (input_file_is_open) CloseInputHDF5();
	TrainsInFile = 0;
	trigger_exists = false;
	trigger_bin_position = -1000000000000000000LL;
	last_trigger_train_ID = -1;
}

/////////////////////////////////////////////////////////////////
bool HDF5_IO::OpenInputHDF5(const char* file_name)
/////////////////////////////////////////////////////////////////
{
	#ifdef USE_HDF5_STUFF
		if (!init_was_successful) return false;
	reset();

	// Create variable length memory datatype:
	memtype = H5Tvlen_create(H5T_NATIVE_INT);
	if (!memtype) {error_flag = 1; return false;}

	// Open file
	HDF5_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (!HDF5_file) {
		if (memtype) {H5Tclose(memtype); memtype = 0;}
		error_flag = 1;
		return false;
	}

	// Open ADC channels datasets:	
	for (int i = 0; i < Num_HDF5_ADC_Channels; i++) {
			ds_ADC[i] = 0;
			ds_ADC_TS[i] = 0;
			ds_ADC_TID[i] = 0;
		char name[50];

			sprintf_s(name, "digitizer/channel_%d_%c/samples", i / 4, 65 + (i % 4));
		ds_ADC[i]		= H5Dopen(HDF5_file, name	, H5P_DEFAULT);
			if (ds_ADC[i] > 0) {
				sprintf_s(name, "digitizer/channel_%d_%c/timestamp", i / 4, 65 + (i % 4));
		ds_ADC_TS[i]	= H5Dopen(HDF5_file, name	, H5P_DEFAULT);
				sprintf_s(name, "digitizer/channel_%d_%c/TrainID", i / 4, 65 + (i % 4));
		ds_ADC_TID[i]	= H5Dopen(HDF5_file, name	, H5P_DEFAULT);
			} else ds_ADC[i] = 0;
	}

	// full file
	/*
	char tbuf[1000];
	for (int i = 0; i < Num_HDF5_ADC_Channels;i++) {
		sprintf_s(tbuf, "digitizer/channel_%d_%c/samples", i / 4, 65 + (i % 4));
		ds_ADC[0] = H5Dopen(file, tbuf, H5P_DEFAULT);
		sprintf_s(tbuf, "digitizer/channel_%d_%c/timestamp", i / 4, 65 + (i % 4));
		ds_ADC_TS[0] = H5Dopen(file, tbuf, H5P_DEFAULT);
		sprintf_s(tbuf, "digitizer/channel_%d_%c/TrainID", i / 4, 65 + (i % 4));
		ds_ADC_TID[0] = H5Dopen(file, tbuf, H5P_DEFAULT);
	}
	*/

	// Allocate memory for ADC data
	for (int i = 0; i < Num_HDF5_ADC_Channels; i++) {
			space_ADC[i] = 0;
			space_ADC_TS[i] = 0;
			space_ADC_TID[i] = 0;
			MaxNumberOfPackets[i] = 0;
			ADCdata[i] = 0;
			ADC_TS[i]  = 0;
			ADC_TID[i] = 0;

			if (ds_ADC[i] > 0) {
		space_ADC[i] = H5Dget_space(ds_ADC[i]);
		ndims = H5Sget_simple_extent_dims(space_ADC[i], dims, NULL);
				ADCdata[i] = new hvl_t[size_t(dims[0])];
				// ADCdata[i] = (hvl_t *)malloc(size_t(dims[0]) * sizeof(hvl_t));
		MaxNumberOfPackets[i] = (int)dims[0];
			}


			if (ds_ADC_TS[i] > 0) {
		space_ADC_TS[i] = H5Dget_space(ds_ADC_TS[i]);
		ndims = H5Sget_simple_extent_dims(space_ADC_TS[i], dims, NULL);
				ADC_TS[i] = new int[size_t(dims[0])];
				// ADC_TS[i] = (int *)malloc(size_t(dims[0]) * sizeof(int));
			}


			if (ds_ADC_TID[i] > 0) {
		space_ADC_TID[i] = H5Dget_space(ds_ADC_TID[i]);
		ndims = H5Sget_simple_extent_dims(space_ADC_TID[i], dims, NULL);
				ADC_TID[i] = new int[size_t(dims[0])];
				// ADC_TID[i] = (int *)malloc(size_t(dims[0]) * sizeof(int));
			}
	}

	// Open photon datasets
	ds_ph_TID = H5Dopen2(HDF5_file, "photons/TrainID", H5P_DEFAULT);
	ds_ph_iTD = H5Dopen2(HDF5_file, "photons/intensityTD", H5P_DEFAULT);
	// Allocate mem. for photon intensity array
	space_ph_iTD = H5Dget_space(ds_ph_iTD);
	ndims = H5Sget_simple_extent_dims(space_ph_iTD, dims, NULL);

		ph_iTD    = new float * [size_t(dims[0])];
		ph_iTD[0] = new float[size_t(dims[0] * dims[1])];

		//ph_iTD = (float **)malloc(size_t(dims[0] * sizeof(float *)));
		//ph_iTD[0] = (float *)malloc(size_t(dims[0] * dims[1] * sizeof(float)));
	for (int i = 1; i < dims[0]; i++) ph_iTD[i] = ph_iTD[0] + i * dims[1];

	TrainsInFile = (int)dims[0];
		current_train_ID = -1;
	first_train_ID = ph_TID[0];
		trigger_bin_position = -1000000000000LL;

		memset(current_packet_index, 0, sizeof(__int64)*Num_HDF5_ADC_Channels);

	// ADC channels
	for (int i = 0; i < Num_HDF5_ADC_Channels; i++) {
			if (ds_ADC[i]	  > 0) status = H5Dread(ds_ADC[i],		memtype,		 H5S_ALL, H5S_ALL, H5P_DEFAULT, ADCdata[i]);
			if (ds_ADC_TS[i]  > 0) status = H5Dread(ds_ADC_TS[i],	 H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ADC_TS[i]);
			if (ds_ADC_TID[i] > 0) status = H5Dread(ds_ADC_TID[i],	 H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ADC_TID[i]);
	}

	// Photon properties
	status = H5Dread(ds_ph_TID, H5T_NATIVE_INT,   H5S_ALL, H5S_ALL, H5P_DEFAULT, ph_TID);
	status = H5Dread(ds_ph_iTD, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ph_iTD[0]);

	number_of_channels = Num_HDF5_ADC_Channels;

		// TODO: check the following 3 variables with values:
		max_number_of_hits = 25;
		int64_Numberofevents = 0;
		bits_per_mVolt = 65536/500; // 16 bits for 500 mV full scale?

	input_file_is_open = true;
	#endif
	return true;
}



/////////////////////////////////////////////////////////////////
void HDF5_IO::CloseInputHDF5()
/////////////////////////////////////////////////////////////////
{
	#ifdef USE_HDF5_STUFF
	if (!input_file_is_open) return;

	// Close and release resources.  
	for (int i = 0; i < Num_HDF5_ADC_Channels; i++) {
		status = H5Dvlen_reclaim(memtype, space_ADC[i], H5P_DEFAULT, ADCdata[i]);
			if (ADCdata[i]) delete[] ADCdata[i];
			if (ADC_TS[i])  delete[] ADC_TS[i];
			if (ADC_TID[i]) delete[] ADC_TID[i];

			if (ds_ADC[i])			status = H5Dclose(ds_ADC[i]);
			if (space_ADC[i])		status = H5Sclose(space_ADC[i]);
			if (ds_ADC_TS[i])		status = H5Dclose(ds_ADC_TS[i]);
			if (space_ADC_TS[i])	status = H5Sclose(space_ADC_TS[i]);
			if (ds_ADC_TID[i])		status = H5Dclose(ds_ADC_TID[i]);
			if (space_ADC_TID[i])	status = H5Sclose(space_ADC_TID[i]);
	}

		if (ph_iTD[0]) delete[] ph_iTD[0];
		if (ph_iTD)    delete[] ph_iTD;
	status = H5Sclose(space_ph_iTD);

		if (memtype) status = H5Tclose(memtype);
		if (HDF5_file) status = H5Fclose(HDF5_file);

	memtype = 0;
	HDF5_file = 0;

	input_file_is_open = false;
	#endif
}



/////////////////////////////////////////////////////////////////
bool HDF5_IO::ReadNextHDF5packet_slicing(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words)
/////////////////////////////////////////////////////////////////
{
	#ifndef USE_HDF5_STUFF
		return true;
	#else
		if (!init_was_successful) return false;
	if (error_flag) return false;
		if (eof) return false;
	if (!input_file_is_open) {
		error_flag = 1;
		return false;
	}
	bEnd_of_group_detected = false;
		static bool group_is_not_empty = false;

		if (do_grouping) {
			error_flag = 1;
			return false;
		}


	L100:	

		if (trigger_bin_position < 0 || current_train_ID < 0) {
			current_train_ID = 9220000000000000000LL;
			for (int ch = 0; ch < Num_HDF5_ADC_Channels; ch++) {
				if (current_packet_index[ch] < MaxNumberOfPackets[ch]) {
					if (current_train_ID > ADC_TID[ch][current_packet_index[ch]]) current_train_ID = ADC_TID[ch][current_packet_index[ch]];
				}
			}
			trigger_bin_position = 0;
		}

	// TODO: handle Photon intensity here.
		// ph_TID[current_train_ID];

	for (int ch = 0; ch < Num_HDF5_ADC_Channels; ch++) {
		if (current_packet_index[ch] >= MaxNumberOfPackets[ch]) continue;

			if (ADC_TID[ch][current_packet_index[ch]] < current_train_ID) {
				current_packet_index[ch]++;
				ch--;
				continue;
			}
			if (ADC_TID[ch][current_packet_index[ch]] > current_train_ID) continue;
			
			__int64 TS = ADC_TS[ch][current_packet_index[ch]];
		if (TS - trigger_bin_position < 0) {
			current_packet_index[ch]++;
			ch--;
			continue;
		}
		if (TS - trigger_bin_position > group_slice_size_in_bins) continue;

			group_is_not_empty	= true;
			packet->channel		= ch%4;
			packet->card		= ch/4;
		packet->flags = 0;
		packet->type = CRONO_PACKET_TYPE_16_BIT_SIGNED;
				packet->timestamp = TS - trigger_bin_position;

			packet->length = ADCdata[ch][current_packet_index[ch]].len;

			ptr_32 = (__int32*)ADCdata[ch][current_packet_index[ch]].p;
		
			int len = packet->length; // number of 16 bit_words to copy
		if (len > buffersize_in_16bit_words) len = buffersize_in_16bit_words;
			for (__int32 i = 0; i < len; i++) {
				__int32 val = ptr_32[i];
				if (val >  32767) val =  32767;
				if (val < -32767) val = -32767;
				i16buffer[i] = __int16(val);
			}
			current_packet_index[ch]++;
		return true;
	}

		if (!group_is_not_empty) {
			bool at_least_one_channel_still_has_data = false;
			for (int ch = 0; ch < Num_HDF5_ADC_Channels; ch++) {
				if (current_packet_index[ch] < MaxNumberOfPackets[ch]) {
					at_least_one_channel_still_has_data = true;
					break;
				}
			}
			if (!at_least_one_channel_still_has_data) {
				eof = true;
				group_is_not_empty = false;
				return false;
			}

			bool at_least_one_channel_still_has_data_in_current_train = false;
			for (int ch = 0; ch < Num_HDF5_ADC_Channels; ch++) {
				if (ADC_TID[ch][current_packet_index[ch]] == current_train_ID) {
					at_least_one_channel_still_has_data_in_current_train = true;
					break;
				}
			}

			if (at_least_one_channel_still_has_data_in_current_train) {
		trigger_bin_position += group_slice_size_in_bins;
			} else {
				trigger_bin_position = -1000000000000000000LL;
	}
			goto L100;
		}

		bEnd_of_group_detected = true;
		group_is_not_empty = false;

				int64_number_of_read_events++;
		trigger_bin_position += group_slice_size_in_bins;
			
		return true;
	#endif
			}





/////////////////////////////////////////////////////////////////
bool HDF5_IO::ReadNextHDF5packet_grouping(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words)
/////////////////////////////////////////////////////////////////
{
	#ifndef USE_HDF5_STUFF
			return true;
	#else
		if (!init_was_successful) return false;
		if (error_flag) return false;
		if (eof) return false;
		if (!input_file_is_open) {
			error_flag = 1;
			return false;
		}
		bEnd_of_group_detected = false;

		if ((!do_grouping) || TriggerChannel < 0) {
			error_flag = 1;
		return false;
}




		int ch = 0;
		__int64 TS = 0;

		while (!trigger_exists) {
			ch = TriggerChannel;
			if (current_packet_index[ch] >= MaxNumberOfPackets[ch]) {
				eof = true;
				return false;
			}
			TS = ADC_TS[ch][current_packet_index[ch]];
			current_train_ID = ADC_TID[ch][current_packet_index[ch]];
			if (current_train_ID == last_trigger_train_ID) {
				if (TS - trigger_bin_position < group_range_end_bins) {
					current_packet_index[ch]++;
					continue;
				}
			}
			last_trigger_train_ID = current_train_ID;
			trigger_exists = true;
			trigger_bin_position = TS;
			int64_number_of_read_events++;

			// TODO: handle Photon intensity here.
			//ph_TID[current_train_ID];
			
			packet->channel		= ch%4;
			packet->card		= ch/4;
			packet->flags		= 0;
			packet->type		= CRONO_PACKET_TYPE_16_BIT_SIGNED;
			packet->timestamp	= TS - trigger_bin_position;
			packet->length = ADCdata[ch][current_packet_index[ch]].len;
			ptr_32 = (__int32*)ADCdata[ch][current_packet_index[ch]].p;
			int len = packet->length; // number of 16 bit_words to copy
			if (len > buffersize_in_16bit_words) len = buffersize_in_16bit_words;
			for (__int32 i = 0; i < len; i++) {
				__int32 val = ptr_32[i];
				if (val >  32767) val =  32767;
				if (val < -32767) val = -32767;
				i16buffer[i] = __int16(val);
			}

			current_packet_index[ch]++;
			return true;
		}

		for (ch = 0; ch < Num_HDF5_ADC_Channels; ch++) {
			if (current_packet_index[ch] >= MaxNumberOfPackets[ch]) continue;

			if (ADC_TID[ch][current_packet_index[ch]] > current_train_ID) continue;
			if (ADC_TID[ch][current_packet_index[ch]] < current_train_ID) {
				current_packet_index[ch]++;
				ch--;
				continue;
			}

			TS = ADC_TS[ch][current_packet_index[ch]];
			if (TS - trigger_bin_position < group_range_start_bins) {
				current_packet_index[ch]++;
				ch--;
				continue;
			}
			if (TS - trigger_bin_position > group_range_end_bins) continue;

			packet->channel		= ch%4;
			packet->card		= ch/4;
			packet->flags		= 0;
			packet->type		= CRONO_PACKET_TYPE_16_BIT_SIGNED;
			packet->timestamp	= TS - trigger_bin_position;
			packet->length = ADCdata[ch][current_packet_index[ch]].len;
			ptr_32 = (__int32*)ADCdata[ch][current_packet_index[ch]].p;
			int len = packet->length; // number of 16 bit_words to copy
			if (len > buffersize_in_16bit_words) len = buffersize_in_16bit_words;
			for (__int32 i = 0; i < len; i++) {
				__int32 val = ptr_32[i];
				if (val >  32767) val =  32767;
				if (val < -32767) val = -32767;
				i16buffer[i] = __int16(val);
		}

			current_packet_index[ch]++;
			return true;
	}
	
		bEnd_of_group_detected = true;
		trigger_exists = false;
		return true;

	#endif
}



void HDF5_IO::seek_to_begin()
{
	eof = false;
	current_train_ID = -1;
	last_trigger_train_ID = -1;
	trigger_exists = false;
	trigger_bin_position = -1000000000000000000LL;
	int64_number_of_read_events = 0;
	memset(current_packet_index, 0, sizeof(__int64)*Num_HDF5_ADC_Channels);
}

#endif
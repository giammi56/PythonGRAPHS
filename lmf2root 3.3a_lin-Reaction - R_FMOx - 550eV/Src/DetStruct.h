#pragma once

#ifndef DETSTRUCT_IS_ALREADY_INCLUDED
	#define DETSTRUCT_IS_ALREADY_INCLUDED

	#include "OS_Version.h"
	#include "resort64c.h"

	struct Det_struct{
		bool	use_this_detector;
		double	offset_timesum_U;
		double  offset_timesum_V;
		double  offset_timesum_W;
		double	hex_offset_w;
		double	center_X;
		double	center_Y;
		__int32	auto_calibration;
		bool	use_reconstruction;
		bool	use_sum_tracker_layer_u;
		bool	use_sum_tracker_layer_v;
		bool	use_sum_tracker_layer_w;
		sort_class * sorter;

		// outputs
		__int32	number_of_reconstructed_hits;
		__int32	method[NUM_IONS];
		double	x[NUM_IONS];
		double	y[NUM_IONS];
		double	time[NUM_IONS];
		double	tof[NUM_IONS];
		double	phi[NUM_IONS];
		peak_tracker_class * sumu_watchdog;
		peak_tracker_class * sumv_watchdog;
		peak_tracker_class * sumw_watchdog;

		void * Ueberstruct_pointer;
	};

#endif
#pragma once

#include "OS_Version.h"
#ifndef dont_use_MFC

#ifndef PARALLEL_SORT_CLASS_IS_DEFINED
#define PARALLEL_SORT_CLASS_IS_DEFINED


#include "resort64c.h"


#ifndef NO_RESORT64C_DLL_IMPORT
#define RESORT64C_API
#endif



class parallel_sort_class;


class event_class{
public:
	event_class(__int32 _number_of_channels, __int32 _number_of_hits);
	~event_class();

	volatile bool	please_terminate_this_thread;
	volatile bool	thread_is_running;

	volatile void *  Sort_Cthread;
	volatile __int32 ThreadID;

	double * tdc_ns;
	__int32 * cnt;

	volatile __int32 number_of_channels;
	volatile __int32 number_of_hits;

	volatile bool	sorter_proj_use_sorting;
	volatile bool	sorter_rec_use_sorting;
	volatile bool	sorter_elec_use_sorting;

	volatile parallel_sort_class * mother_parallel_sort_class;

	sort_class * sorter_proj;
	sort_class * sorter_rec;
	sort_class * sorter_elec;

	void * Event_do_sort_data;
	void * Event_do_read_sorted_data;
};





class RESORT64C_API parallel_sort_class {
public:
	parallel_sort_class(__int32 number_of_threads, __int32 number_of_channels, __int32 number_of_hits, sort_class * sorter1, bool use_sorting1_algorithm, sort_class * sorter2, bool use_sorting2_algorithm, sort_class * sorter3,bool use_sorting3_algorithm);
	~parallel_sort_class();

	void	start_sort_threads();
	void	terminate_sort_threads();
	void	wait_for_thread_with_ID(__int32 ID);
	void	do_sort_with_thread_ID(__int32 ID);
	static unsigned __int32	sort_thread(void* pointer_to_event_class);

	event_class ** events;

	volatile __int32		number_of_channels;
	volatile __int32		number_of_hits;
	volatile __int32		number_of_threads;
	//double * temp_tdc;
	//__int32 *	temp_cnt;

	volatile bool			sorter_proj_use_sorting;
	volatile bool			sorter_rec_use_sorting;
	volatile bool			sorter_elec_use_sorting;

	__int32					filled_threads;
	__int32					next_thread_ID_to_be_sorted;

private:
	void			reset_tdc_pointers(sort_class * sorter);
	void			init(__int32 number_of_threads, __int32 number_of_channels, __int32 number_of_hits, sort_class * sorter1, sort_class * sorter2, sort_class * sorter3);
};

#endif
#endif

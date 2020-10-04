#include "OS_Version.h"

#ifndef dont_use_MFC

	#include "parallel_sort_class.h"
	#include "afxmt.h"
	#include "afxwin.h"
	
	#ifdef _DEBUG
	#include <assert.h>
	#endif
	
	
	#define _Cu1 (0)
	#define _Cu2 (1)
	#define _Cv1 (2)
	#define _Cv2 (3)
	#define _Cw1 (4)
	#define _Cw2 (5)
	#define _Cmcp (6)
	
	
	
	
	
	/////////////////////////////////////////////////////////
	event_class::event_class(__int32 _number_of_channels, __int32 _number_of_hits)
	/////////////////////////////////////////////////////////
	{
		CEvent *ce1 = new CEvent(false,true,0,0);
		Event_do_sort_data = (void*)ce1;
		ce1->ResetEvent();

		CEvent *ce2 = new CEvent(false,true,0,0);
		Event_do_read_sorted_data = (void*)ce2;
		ce2->ResetEvent();

		thread_is_running = false;
		please_terminate_this_thread = false;
		sorter_elec = sorter_proj = sorter_rec = 0;
		Sort_Cthread = 0;
		number_of_channels	= _number_of_channels;
		number_of_hits		= _number_of_hits;
		tdc_ns = new double[number_of_channels*number_of_hits];
		memset(tdc_ns, 0, number_of_channels*number_of_hits * sizeof(double));
		cnt = new int[number_of_channels];
		memset(cnt, 0, number_of_channels * sizeof(__int32));
	}

	
	/////////////////////////////////////////////////////////
	event_class::~event_class()
	/////////////////////////////////////////////////////////
	{
		if (tdc_ns) {delete[] tdc_ns; tdc_ns = 0;}
		if (cnt) {delete[] cnt; cnt = 0;}
		if (this->Event_do_sort_data) {
			CEvent *ce = (CEvent*)this->Event_do_sort_data;
			if (ce) delete ce;
			ce = 0;
			this->Event_do_sort_data = 0;
		}
		if (this->Event_do_read_sorted_data) {
			CEvent *ce = (CEvent*)this->Event_do_read_sorted_data;
			if (ce) delete ce;
			ce = 0;
			this->Event_do_read_sorted_data = 0;
		}
	}
	
	
	
	
	
	
	/////////////////////////////////////////////////////////
	parallel_sort_class::parallel_sort_class(__int32 _number_of_threads, __int32 _number_of_channels, __int32 _number_of_hits, sort_class * _sorter_proj, bool _sorter_proj_use_sorting, sort_class * _sorter_rec, bool _sorter_rec_use_sorting, sort_class * _sorter_elec, bool _sorter_elec_use_sorting)
	/////////////////////////////////////////////////////////
	{
		events = 0;
		sorter_proj_use_sorting = _sorter_proj_use_sorting;
		sorter_rec_use_sorting = _sorter_rec_use_sorting;
		sorter_elec_use_sorting = _sorter_elec_use_sorting;
		init(_number_of_threads, _number_of_channels, _number_of_hits, _sorter_proj, _sorter_rec, _sorter_elec);
	}
	
	
	
	
	
	/////////////////////////////////////////////////////////
	parallel_sort_class::~parallel_sort_class()
	/////////////////////////////////////////////////////////
	{
		this->terminate_sort_threads();
	
		bool thread_still_running = false;
		if (events) {
			for (__int32 i=0;i<number_of_threads;++i) {
				if (!events[i]) continue;
				if (events[i]->thread_is_running) thread_still_running = true;
			}
		}
	
		if (thread_still_running) this->terminate_sort_threads();
	
		if (events) {
			for (__int32 i=0;i<number_of_threads;++i) {
				if (!events[i]) continue;
		
				if (events[i]->sorter_proj) {delete events[i]->sorter_proj; events[i]->sorter_proj = 0;}
				if (events[i]->sorter_rec)  {delete events[i]->sorter_rec;  events[i]->sorter_rec = 0; }
				if (events[i]->sorter_elec) {delete events[i]->sorter_elec; events[i]->sorter_elec = 0;}
	
				if (events[i]->Sort_Cthread) {
					CWinThread * t = (CWinThread *)events[i]->Sort_Cthread;
					if (t) delete t;
					t=0;
					events[i]->Sort_Cthread = 0;
				}
				if (events[i]) delete events[i];
				events[i] = 0;
			}
			delete[] events;
			events = 0;
		}
	}
	
	
	
	
	
		
	
	
	/////////////////////////////////////////////////////////
	void parallel_sort_class::init(__int32 _number_of_threads, __int32 _number_of_channels, __int32 _number_of_hits, sort_class * _sorter_proj, sort_class * _sorter_rec, sort_class * _sorter_elec)
	/////////////////////////////////////////////////////////
	{
		number_of_threads = _number_of_threads;
		events = new event_class*[number_of_threads];
		memset(events,0,number_of_threads*sizeof(event_class*));

		number_of_channels = _number_of_channels;
		number_of_hits = _number_of_hits;
		filled_threads = 0;
		next_thread_ID_to_be_sorted = 0;
	
		for (__int32 i = 0;i < number_of_threads; ++i) {
			events[i] = new event_class(number_of_channels, number_of_hits);
			events[i]->ThreadID = i;
			events[i]->mother_parallel_sort_class = this;
	
			if (_sorter_proj) {
				events[i]->sorter_proj = new sort_class();
				_sorter_proj->clone(events[i]->sorter_proj); 
				events[i]->sorter_proj->tdc_pointer = events[i]->tdc_ns; 
				events[i]->sorter_proj->count = events[i]->cnt; 
				reset_tdc_pointers(events[i]->sorter_proj);
			}
			
			if (_sorter_rec) {
				events[i]->sorter_rec = new sort_class();
				_sorter_rec->clone(events[i]->sorter_rec); 
				events[i]->sorter_rec->tdc_pointer = events[i]->tdc_ns;
				events[i]->sorter_rec->count = events[i]->cnt; 
				reset_tdc_pointers(events[i]->sorter_rec);
			}

			if (_sorter_elec) {
				events[i]->sorter_elec = new sort_class();
				_sorter_elec->clone(events[i]->sorter_elec);
				events[i]->sorter_elec->tdc_pointer = events[i]->tdc_ns;
				events[i]->sorter_elec->count = events[i]->cnt;
				reset_tdc_pointers(events[i]->sorter_elec);
			}
	
			events[i]->sorter_proj_use_sorting = sorter_proj_use_sorting;
			events[i]->sorter_rec_use_sorting = sorter_rec_use_sorting;
			events[i]->sorter_elec_use_sorting = sorter_elec_use_sorting;
		}
	}
	
	
	
	
	
	
	/////////////////////////////////////////////////////////
	void parallel_sort_class::start_sort_threads()
	/////////////////////////////////////////////////////////
	{
		for (__int32 i=0;i<number_of_threads;++i) {
			if (!events[i]->thread_is_running) {
				static CWinThread * this_thread;
				this_thread = AfxBeginThread(&sort_thread	,(LPVOID)events[i],THREAD_PRIORITY_NORMAL,0,0,NULL);
	//			this_thread = AfxBeginThread(&sort_thread	,(LPVOID)events[i],THREAD_PRIORITY_IDLE,0,0,NULL);
				events[i]->Sort_Cthread = (void*)this_thread;
				Sleep(10);
			}
		}
	}
	
	
	
	
	
	
	/////////////////////////////////////////////////////////
	void parallel_sort_class::terminate_sort_threads()
	/////////////////////////////////////////////////////////
	{
		bool terminating_threads = true;
		if (!events) return;
		while (terminating_threads) {
			terminating_threads = false;
			for (__int32 i=0;i<number_of_threads;++i) {
				if (!events[i]) continue;
				events[i]->please_terminate_this_thread = true;
				if (events[i]->thread_is_running) terminating_threads = true;
			}
			Sleep(5);
		}
	}
	
	
	
		

	
	
	/////////////////////////////////////////////////////////
	void parallel_sort_class::reset_tdc_pointers(sort_class * sorter)
	/////////////////////////////////////////////////////////
	{
		if (sorter) {
			if (sorter->initialization_successful) {
				for (__int32 i=0;i<7;i++) sorter->tdc[i] = (double*)0;

				if (sorter->Cu1 >= 0 && sorter->use_u1) sorter->tdc[_Cu1] = sorter->tdc_pointer + sorter->Cu1 * sorter->tdc_array_row_length;
				if (sorter->Cu2 >= 0 && sorter->use_u2) sorter->tdc[_Cu2] = sorter->tdc_pointer + sorter->Cu2 * sorter->tdc_array_row_length;
				if (sorter->Cv1 >= 0 && sorter->use_v1) sorter->tdc[_Cv1] = sorter->tdc_pointer + sorter->Cv1 * sorter->tdc_array_row_length;
				if (sorter->Cv2 >= 0 && sorter->use_v2) sorter->tdc[_Cv2] = sorter->tdc_pointer + sorter->Cv2 * sorter->tdc_array_row_length;
	
				if (sorter->use_HEX) {
					if (sorter->Cw1 >= 0 && sorter->use_w1) sorter->tdc[_Cw1] = sorter->tdc_pointer + sorter->Cw1 * sorter->tdc_array_row_length;
					if (sorter->Cw2 >= 0 && sorter->use_w2) sorter->tdc[_Cw2] = sorter->tdc_pointer + sorter->Cw2 * sorter->tdc_array_row_length;
				} 
	
				if (sorter->Cmcp >= 0 && sorter->use_MCP) {
					sorter->tdc[_Cmcp] = sorter->tdc_pointer + sorter->Cmcp * sorter->tdc_array_row_length;
				}
				sorter->tdc_Cu1 = sorter->tdc[_Cu1];
				sorter->tdc_Cu2 = sorter->tdc[_Cu2];
				sorter->tdc_Cv1 = sorter->tdc[_Cv1];
				sorter->tdc_Cv2 = sorter->tdc[_Cv2];
				sorter->tdc_Cw1 = sorter->tdc[_Cw1];
				sorter->tdc_Cw2 = sorter->tdc[_Cw2];
				sorter->tdc_Cmcp = sorter->tdc[_Cmcp];
			}
		}
	}
	
	





	/////////////////////////////////////////////////////////
	void parallel_sort_class::do_sort_with_thread_ID(__int32 ID)
	/////////////////////////////////////////////////////////
	{
		CEvent * ce =  (CEvent*)events[ID]->Event_do_sort_data;
		ce->SetEvent();
		filled_threads++;
		if (filled_threads > number_of_threads) filled_threads = number_of_threads;
	}
	



	/////////////////////////////////////////////////////////
	void parallel_sort_class::wait_for_thread_with_ID(__int32 ID)
	/////////////////////////////////////////////////////////
	{
		while (!events[ID]->please_terminate_this_thread) {
			CEvent * ce = (CEvent*)events[ID]->Event_do_read_sorted_data;
			if (::WaitForSingleObject(ce->m_hObject, 51) == WAIT_TIMEOUT) continue;
			ce->ResetEvent();
			break;
		}
	}
	


	
	
	
	
		
	
	/////////////////////////////////////////////////////////
	UINT parallel_sort_class::sort_thread(LPVOID voidpointer)
	/////////////////////////////////////////////////////////
	{
		event_class * my_event = (event_class*)voidpointer;
		my_event->thread_is_running = true;
	
		bool xxx = true;

		while(!my_event->please_terminate_this_thread) {
			if (xxx) {
				xxx = false;
				//TRACE("in thread %i: waiting to start sorting\n", my_event->ThreadID);
			}
			while(!my_event->please_terminate_this_thread) {
				CEvent * ce = (CEvent*)my_event->Event_do_sort_data;
				if (::WaitForSingleObject(ce->m_hObject, 51) == WAIT_TIMEOUT) continue;
				ce->ResetEvent();
				break;
			}
			if (my_event->please_terminate_this_thread) break;
	
			xxx = false;
			//TRACE("in thread %i: sorting started\n",my_event->ThreadID);
		
			if (my_event->sorter_proj) {
				if (my_event->sorter_proj_use_sorting) my_event->sorter_proj->sort(); else my_event->sorter_proj->run_without_sorting();
			}
	
			if (my_event->sorter_rec) {
				if (my_event->sorter_rec_use_sorting) my_event->sorter_rec->sort(); else my_event->sorter_rec->run_without_sorting();
			}

			if (my_event->sorter_elec) {
				if (my_event->sorter_elec_use_sorting) my_event->sorter_elec->sort(); else my_event->sorter_elec->run_without_sorting();
			}
	
			//TRACE("in thread %i: sort finished\n", my_event->ThreadID);
			CEvent * ce = (CEvent*)my_event->Event_do_read_sorted_data; ce->SetEvent();
		}
	
		my_event->thread_is_running = false;
		AfxEndThread(1,FALSE);
		return TRUE;
	}

#endif
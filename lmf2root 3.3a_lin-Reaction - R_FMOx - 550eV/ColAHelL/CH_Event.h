#ifndef CH_EVENT_ALREADY_INCLUDED
	#define CH_EVENT_ALREADY_INCLUDED

	namespace CH
	{
		struct CH_det_struct{
			int num_hits;
			std::vector<double> method;
			std::vector<double> x;
			std::vector<double> y;
			std::vector<double> time;
			std::vector<double> tof;
		};

		struct CH_event_struct{
		public:	
			__int64 event_number;
			double reaction;
			double bunchmarker;
			CH_det_struct r;
			CH_det_struct e;
			CH_det_struct p;
			std::vector<double> scanval;
		};

		struct CH_presorter_result_struct{
		public:	
			int channel;
			bool *valid_array_e;
			bool *valid_array_r;
			bool *valid_array_p;
		};
	}
#endif;
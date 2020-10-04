
#include "OS_Version.h"

int lmf2root_main(__int32 argc, char* argv[]);

#ifndef dont_use_mfc
	bool init_MFC();
#endif

/////////////////////////////////////////////////////////////////////////////
int main(__int32 argc, char* argv[], char* envp[])
/////////////////////////////////////////////////////////////////////////////
{
	#ifndef dont_use_MFC
		if (!init_MFC()) return 1;
	#endif
	lmf2root_main(argc, argv);
}

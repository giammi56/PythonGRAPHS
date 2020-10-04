#pragma once






#ifndef OS_VERSION_DEFINED_FOR_LMF2ROOT
	#define OS_VERSION_DEFINED_FOR_LMF2ROOT

#ifdef _DEBUG
	#define _CRTDBG_MAP_ALLOC
#endif


#define NO_WARN_MBCS_MFC_DEPRECATION

#ifndef WINVER                          // Specifies that the minimum required platform is Windows Vista.
#define WINVER 0x0600           // Change this to the appropriate value to target other versions of Windows.
#endif

#ifndef _WIN32_WINNT            // Specifies that the minimum required platform is Windows Vista.
#define _WIN32_WINNT 0x0600     // Change this to the appropriate value to target other versions of Windows.
#endif

#ifndef _WIN32_WINDOWS          // Specifies that the minimum required platform is Windows 98.
#define _WIN32_WINDOWS 0x0410 // Change this to the appropriate value to target Windows Me or later.
#endif

#ifndef _WIN32_IE                       // Specifies that the minimum required platform is Internet Explorer 7.0.
#define _WIN32_IE 0x0700        // Change this to the appropriate value to target other versions of IE.
#endif



#ifndef WIN32
	#define LINUX
#endif

	#pragma warning(disable : 4996)
	#pragma warning(disable : 4800)

	//ADC_analysis definitions:
	#define NUM_CHANNELS (100)
	#define NUM_IONS (100)
	#define MEAN_PULSE_LENGTH (100)

	#ifdef LINUX
		#define dont_use_MFC
		#ifndef __int32_IS_DEFINED
			#define __int32_IS_DEFINED
		#define __int32 int
		#define __int16 short
		#define __int64 long long
		#define __int8 char
		#endif
	#endif

#endif



#ifndef LINUX
#ifdef dont_use_MFC
#error: Do not use Release_std or Debug_std without a good reason. Ask Achim. Then remove this line.
#endif
#endif


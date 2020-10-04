project (lmf-to-root)

SET ( BUILD_TYPE Release)
SET ( BUILD_TYPE Debug)

SET ( lmf2root_SOURCES main_console_app.cpp LMF_IO.cpp mySleep.cpp rootstuff.cpp LMF2root.cpp analysis.cpp console.cpp presorter.cpp detector_stuff.cpp functions.cpp sort_and_write_NTuple.cpp Histo.cpp config_file_parser.cpp help.cpp )

FOREACH(cxxfile ${lmf2root_SOURCES})
	STRING(REPLACE ".cpp" ".h" hfile ${cxxfile})
	SET ( lmf2root_HEADERS ${lmf2root_HEADERS} ${hfile})
ENDFOREACH(cxxfile)
#MESSAGE("[i] Headers: ${lmf2root_HEADERS}")

#MESSAGE("[i] ROOT set: $ENV{ROOTSYS}")
#MESSAGE("[i] AliROOT set: $ENV{ALICE_ROOT}")


EXECUTE_PROCESS ( COMMAND root-config --cflags WORKING_DIRECTORY /tmp OUTPUT_VARIABLE ROOTCFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE )
MESSAGE("[i] ROOTCFLAGS set ${ROOTCFLAGS} ")
EXECUTE_PROCESS ( COMMAND root-config --incdir WORKING_DIRECTORY /tmp OUTPUT_VARIABLE ROOTINCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE )
MESSAGE("[i] ROOTINCLUDES is ${ROOTINCLUDES} ")

MESSAGE("[i] ROOTLIBS set ${ROOTLIBS} ")
EXECUTE_PROCESS ( COMMAND root-config --libdir WORKING_DIRECTORY /tmp OUTPUT_VARIABLE ROOTLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE )
MESSAGE("[i] ROOTLIBDIR set ${ROOTLIBDIR} ")
#MESSAGE("[i] Sources are : ${lmf2root_SOURCES}")


#Make sure the compiler can find include files from our Hello library.
# include_directories (${HELLO_SOURCE_DIR}/Hello)
#INCLUDE_DIRECTORIES ( ${lmf2root_SOURCE_DIR}/src )
INCLUDE_DIRECTORIES (${ROOTINCLUDES})

SET(CMAKE_CXX_FLAGS "${ROOTCFLAGS}")
SET(CMAKE_CPP_FLAGS "${ROOTCFLAGS}")
SET(CMAKE_LD_FLAGS "${ROOTLIBS}" )

# Make sure the linker can find the Hello library once it is built.
#link_directories (${lmf2root_BINARY_DIR} ${ROOTLIBS} ${ALICE_ROOT}/lib/tgt_${ALICE_TARGET} )

LINK_DIRECTORIES (${lmf2root_BINARY_DIR} ${ROOTLIBDIR} )
LINK_LIBRARIES ( /home/czasch/Linux/lmf2root/linux_lib/libResort64c.a )

ADD_LIBRARY ( lmf2root SHARED ${lmf2root_SOURCES} )
#ADD_LIBRARY ( lmf2root SHARED ${CINTDICTS} )
#ADD_LIBRARY ( lmf2root SHARED ${lmf2root_SOURCES} )

#SET_TARGET_PROPERTIES( lmf2root PROPERTIES LINK_FLAGS "${ROOTLIBS}" SO_VERSION 0.0.1)
#TARGET_LINK_LIBRARIES( lmf2root Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore freetype Gui)

SET (ROOTLIBLINK Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore freetype Gui VMC EG Geom Thread dl)

TARGET_LINK_LIBRARIES( lmf2root ${ROOTLIBLINK} )

#SET(optionalSources)
#SET(optionalLibs)

#IF(ROOT_FOUND)
#   SET(optionalSources ${optionalSources} ROOTview.cpp)
#   INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIR} )
#   SET(optionalLibs ${optionalLibs} ${ROOT_LIBRARIES} )
#ENDIF(ROOT_FOUND)

#--------------- EXECS
add_executable ( lmf2root ${lmf2root_SOURCES})
# Link the executable
target_link_libraries( lmf2root ${ROOTLIBLINK}  )


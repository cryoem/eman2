FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

# MESSAGE("original pgm list: ${e2programs}")

FILE(GLOB exclusion_list "e2plotFSC.py" "e2refine_evenodd.py" "e2refine_evenodd.py" "e2refmerge.py" "e2resolution.py" "e2montecarlorecon.py" "e2markbadparticles.py" "e22.py" "e2boxer21.py" "e2classesbyref.py" "e2ali2img.py" "e2avg_ffts.py" "e2fileinfo.py"  "e2flick.py" "e2findsubunit.py" "e2helical_recons.py" "e2modeleval.py" "e2preferences.py" "e2proc3d_huge.py" "e2scp.py" "e2tomoallvall.py" "e2tomoaverage.py" "e2tomohunter.py" "e2remoted.py" "e2tomosim.py")
LIST(REMOVE_ITEM e2programs ${exclusion_list})
#######MESSAGE(SEND_ERROR "${e2programs}")    #for debug purpose ######

# MESSAGE("cut back pgm list: ${e2programs}")

# SET(ED  $ENV{EMANDIR})
# SET(E2D $ENV{EMAN2DIR})
# SET(HOM $ENV{HOME})

SET(DBG $ENV{DEBUG})

# MESSAGE("envs: '${ED}' '${E2D}' '${HOM}' '${DBG}'")

FOREACH(f ${e2programs})
	IF(DBG MATCHES "y")
		MESSAGE("Copied file ${f} to bin.")
	ENDIF(DBG MATCHES "y")
ENDFOREACH(f)

INSTALL(PROGRAMS ${e2programs}
  DESTINATION    bin
)

#INSTALL(PROGRAMS    e2boxer.py 
#                    e2version.py 
#                    e2ctf.py 
#                    e2history.py
#                    #....
#        DESTINATION bin
#)

IF(WIN32)
	INSTALL(FILES e2boxer.py
		DESTINATION bin
		RENAME	e2boxer.pyw
		)
		
	INSTALL(FILES e2display.py
		DESTINATION bin
		RENAME	e2display.pyw
		)
ENDIF(WIN32)

# INSTALL(PROGRAMS e2ctf_auto.py
#       DESTINATION bin
#       )

# MESSAGE("Installed e2ctf_auto.py in bin manually (failed earlier for some reason.")

string(TIMESTAMP EMAN_TIMESTAMP "%Y-%m-%d %H:%M")
set(EMAN_GITHASH $ENV{GIT_FULL_HASH})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/e2version.py.in
				${CMAKE_INSTALL_PREFIX}/bin/e2version.py
				)

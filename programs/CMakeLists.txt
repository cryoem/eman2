FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

FILE(GLOB exclusion_list "e2plotFSC.py" "e2refine_evenodd.py" "e2refine_evenodd.py" "e2refmerge.py" "e2resolution.py" "e2montecarlorecon.py" "e2markbadparticles.py" "e22.py" "e2boxer21.py" "e2classesbyref.py" "e2ali2img.py" "e2avg_ffts.py" "e2fileinfo.py"  "e2flick.py" "e2findsubunit.py" "e2helical_recons.py" "e2modeleval.py" "e2preferences.py" "e2proc3d_huge.py" "e2scp.py" "e2tomoallvall.py" "e2tomoaverage.py" "e2tomohunter.py" "e2remoted.py" "e2tomosim.py")
LIST(REMOVE_ITEM e2programs ${exclusion_list})

SET(DBG $ENV{DEBUG})

FOREACH(f ${e2programs})
	IF(DBG MATCHES "y")
		MESSAGE("Copied file ${f} to bin.")
	ENDIF(DBG MATCHES "y")
ENDFOREACH(f)

INSTALL(PROGRAMS ${e2programs}
  DESTINATION    bin
)

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

string(TIMESTAMP EMAN_TIMESTAMP "%Y-%m-%d %H:%M")
string(SUBSTRING $ENV{GIT_FULL_HASH} 0 7 EMAN_GITHASH)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/e2version.py.in
				${CMAKE_INSTALL_PREFIX}/bin/e2version.py
				)

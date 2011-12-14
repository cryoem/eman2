FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
FILE(GLOB exclusion_list "e2ali2img.py" "e2avg_ffts.py" "e2fileinfo.py"  "e2flick.py" "e2findsubunit.py" "e2foldfitter.py" "e2helical_recons.py" "e2modeleval.py" "e2preferences.py" "e2proc3d_huge.py" "e2ptcltrace.py" "e2scp.py" "e2tomoallvall.py" "e2tomoaverage.py" "e2tomohunter.py" "e2remoted.py")
LIST(REMOVE_ITEM e2programs ${exclusion_list})
#######MESSAGE(SEND_ERROR "${e2programs}")    #for debug purpose ######
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
		
	INSTALL(FILES e2workflow.py
		DESTINATION bin
		RENAME	e2workflow.pyw
		)
ENDIF(WIN32)

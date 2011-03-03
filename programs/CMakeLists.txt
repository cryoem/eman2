FILE(GLOB e2programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
#FILE(GLOB exclusion_list "file1.py" "file2.py" "file3.py")
#LIST(REMOVE_ITEM e2programs ${exclusion_list})
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
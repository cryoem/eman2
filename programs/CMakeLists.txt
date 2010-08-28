LINK_LIBRARIES(EM2)
INSTALL_PROGRAMS(/bin "^.+\\.py$")
INSTALL_PROGRAMS(/bin golem)

IF(WIN32)
	INSTALL(FILES e2boxer.py
		DESTINATION bin
		RENAME	e2boxer.pyw
		)
		
	INSTALL(FILES e2display.py
		DESTINATION bin
		RENAME	e2display.pyw
		)
		
	INSTALL(FILES e2desktop.py
		DESTINATION bin
		RENAME	e2desktop.pyw
		)
		
	INSTALL(FILES e2workflow.py
		DESTINATION bin
		RENAME	e2workflow.pyw
		)
ENDIF(WIN32)

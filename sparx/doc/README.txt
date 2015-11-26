==========================================================================================
TITLE  : README - How to generate sxgui.py from Wiki documents using wikiparser.py
AUTHOR : Toshio Moriya
DATE   : 2015/11/25
==========================================================================================
------------------------------------------------------------------------------------------
Contents
------------------------------------------------------------------------------------------
1. Syntax of Wiki document and parsing rule 
2. How to generate sxgui.py


------------------------------------------------------------------------------------------
1. Parsing rule and syntax of Wiki document 
------------------------------------------------------------------------------------------
1-1. Parsing Rules

a. Empty lines are always ignored.

b. Only information in the target sections will be used by the parser (wikiparser.py).
   The target sections are shown below. Each section must start with "="



1-2. Syntax of Wiki document

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
= Name =

${SXSCRIPT} - ${SHORT_INFO}

= Usage =

${USAGE_IN_COMMAND_LINE}


=== Typical usage ===

# If there is keyword "mpirun" in this section, the parser assumes that this command supports MPI.


== Input ==
# The format of this section can be one of the followings. 
# (1) for command token where user is required enter the value:

    ${COMMDNA_TOKEN_KEY_BASE}:: ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default required ${COMMDNA_TOKEN_TYPE})

# (2) for command token where default value is provided by the sx script

    ${COMMDNA_TOKEN_KEY_BASE}:: ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default ${COMMDNA_TOKEN_DEFAULT})

# If a line does not follow these format, the parser assumes that it is a comment line 
# Put 4 spaces at the head of each line even though this is not requirement.
# Do not include output related command tokens here. Instead, put them in "== Output =="  
# The command tokens appears here must exist in "= Usage =" section.

# Line that has "*" and word "optional" indicates the remaining options are advanced


== Output ==
# Exactly same format as in "== Input =="

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

${SXSCRIPT}
   Name of the sx script file without ".py" extension.
   e.g. "sxviper"

   Used for text of the button in main window.

${SHORT_INFO} 
   A definition (description) for GUI button. 
   e.g. "Validated ''ab initio'' 3D structure determination"

   Used for a help info of the button and label at the top of main tab.
   
${USAGE_IN_COMMAND_LINE}
   The usage in command line.
   e.g "sxviper.py stack  directory  --ir=inner_radius --radius=outer_radius"
   
   Used to extract command tokens (i.e. arguments and options of the command). 
   If token includes no prefix, it is assigned as an argument.
   If token includes prefix of "--" or "-", the parser assume this is an option.
   The order of command tokens here is also kept in the GUI.
   The command tokens appears here must exist either in "== Input ==" or "== Output ==" section.

${COMMDNA_TOKEN_KEY_BASE}
	Key base or name of a command token (i.e. arguments and options) without any prefix "--" or "-".
	e.g. "radius"

	GUI uses this to generate the command line.
	
	If key base contains the keyword defined in construct_token_list_from_wiki function of wikiparser.py, 
	the associated special data type (listed below) will be assigned to this command token.
	
	
	"output"      : Line edit box and output info button
	                GUI also checks the existence of output directory/file before execution of the sx*.py
	                GUI abort the execution if the directory/file exists already
	"directory"   : Line edit box and open directory button
	"image"       : Line edit box and open file buttons for .hdf and .bdb 
	"parameters"  : Line edit box and open file button for all file types 
	"pdb"         : Line edit box and open file button for .pdb 
	"function"    : Two line edit boxes (function name & file path of the container script)
	                and open file button for .py
	
${COMMDNA_TOKEN_LABEL}
	User-friendly name of a command token.  
	e.g. "radius of the particle" (instead of "radius")

	Used as the label of this command token in GUI.

${COMMDNA_TOKEN_HELP}
	Help information about a command token.
	e.g. "has to be less than < int(nx/2)-1" (for "radius")

	Used as the help information in GUI.

${COMMDNA_TOKEN_DEFAULT}
	Default value of a command token.
	
	Used as the default value of this command token in GUI.
	
	GUI also uses default value to find its data type (either int, float, string, bool), 
	please make sure to use appropriate expression for the value.
	For example, use "1.0" instead of "1" for float type.

	In case of string type, it can be multiple words.
	e.g. "current directory"
	
	Please do not use "required" for the 1st word, since it is a part of syntax.

	If the value in GUI is not changed from the default value, GUI excludes 
	this command token from the command line upon its generation.
	
${COMMDNA_TOKEN_TYPE}
	Data type of a command token. Valid types are int, float, string, and bool.
	
	If the key base contains a keyword explained above, the type provided here will be 
	override by a special data type in the Wiki document parser.


------------------------------------------------------------------------------------------
2. How to generate sxgui.py
------------------------------------------------------------------------------------------
(1) Go to the directory contains the Wiki document parser script.

$ cd ~/EMAN2/src/eman2/sparx/templates/wikiparser.py

(2) Execute the Wiki document parser script

$ ./wikiparser.py

The script automatically updates sxgui.py in sparx/bin directory by copying 
sxgui_template.py to the sxgui.py, and inserting the extracted information into it.
Please check the print out of the Wiki document parser.  
It will tell you about incorrect formats of some Wiki documents. 
If this happens, please edit the Wiki documents.

(3) Copy the generated sxgui.py to your EMAN2/bin directory.

$ cp ~/EMAN2/src/eman2/sparx/bin/sxgui.py ~/EMAN2/bin

(4) Run sxgui.py in your project directory

$ sxgui.py &

That's it!

==========================================================================================
END OF README
==========================================================================================

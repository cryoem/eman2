==========================================================================================
TITLE  : README - How to generate sp_gui.py from Wiki documents using wikiparser.py
AUTHOR : Toshio Moriya
DATE   : 2018/03/02
==========================================================================================

------------------------------------------------------------------------------------------
Contents
------------------------------------------------------------------------------------------
1. Syntax of Wiki document and parsing rule.
2. How to generate sp_gui.py.

------------------------------------------------------------------------------------------
1. Basic parsing rules and syntax of Wiki document 
------------------------------------------------------------------------------------------
1-1. Basic parsing rules

a. Following lines are always ignored.
   - Empty line.
   - Line contains only a DokuWiki linebreak "\\". 
   - Line contains only a DokuWiki control macro "~~NOTOC~~"

b. Only information in the target sections will be used by the parser (wikiparser.py).
   The target sections are shown below. Each section must have follow the format of "===== ${SECTION_NAME} ====="


1-2. Syntax of Wiki document

As the basis, the Wiki document must follow the syntax of DokuWiki engine.

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
===== ${SXSCRIPT} =====
${SXSCRIPT_LABEL}: ${SXSCRIPT_SHORT_INFO}

===== Usage =====

${USAGE_IN_COMMAND_LINE}

# The command tokens appears here must exist also in "===== Input =====" section.
# However, the order of command tokens here is not necessary to be the same and not used in GUI.

===== Typical usage =====
# If there is the keyword "mpirun" in this section, the parser assumes that this command supports MPI.

===== Input =====
=== Main Parameters ===
# The command tokens listed in this subsection will be added to "Main" tab in GUI.

# The information in this section are used to extract command tokens (i.e. arguments and options of the command). 
#   If token includes no prefix (empty string), it is assigned as an argument.
#   If token includes prefix of "%%--%%" or "-", the parser assume this is an option.
#   The order of command tokens here is also kept in the GUI.
#   The command tokens appears here must exist either in "===== Usage =====" section.

# The format of this section can be one of the followings. 
# (1) for command token where user is required enter the value:

  ; ${COMMDNA_TOKEN_KEY_PREFIX}${COMMDNA_TOKEN_KEY_BASE} : ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default required ${COMMDNA_TOKEN_TYPE})

# (2) for command token where default value is provided by the sp_ script

  ; ${COMMDNA_TOKEN_KEY_PREFIX}${COMMDNA_TOKEN_KEY_BASE} : ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default ${COMMDNA_TOKEN_DEFAULT})

# (3) for command token of bool type (i.e. a flag option) where default value is provided by the sp_ script but the question is reversed in GUI.

  ; ${COMMDNA_TOKEN_KEY_PREFIX}${COMMDNA_TOKEN_KEY_BASE} : ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default ${COMMDNA_TOKEN_DEFAULT} question reversed in GUI)

# (4) for command token of bool type (i.e. a flag option) where default value is provided by the sp_ script but the default value is reversed in GUI.

  ; ${COMMDNA_TOKEN_KEY_PREFIX}${COMMDNA_TOKEN_KEY_BASE} : ${COMMDNA_TOKEN_LABEL}: ${COMMDNA_TOKEN_HELP} (default ${COMMDNA_TOKEN_DEFAULT} value reversed in GUI)

# If a line does not follow these formats, the parser assumes that it is a comment line. 
# Make sure to put 2 spaces before ";" and 1 spaces after ";" at the head of each line.
# The command tokens appears here must exist in "===== Usage =====" section.

=== Advanced Parameters ===
# The command tokens listed in this subsection will be added to "Advanced" tab in GUI.

# Use exactly same format as in "=== Main Parameters ===" subsection

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

${SXSCRIPT}
	Name of the sp_ script file without ".py" extension.
	e.g. "sp_viper"

	GUI uses this to generate the command line.
	Used for title in popup window tabs for this sp_ script.

${SXSCRIPT_LABEL}
	User friendly name of the sp_ script.
	e.g. "sp_viper"

	Used for text of the button in main window.

${SXSCRIPT_SHORT_INFO} 
	A definition (description) for GUI button. 
	e.g. "Validated //ab initio// 3D structure determination"

	Used for a help info of the button and label at the top of main tab.

${USAGE_IN_COMMAND_LINE}
	The usage in command line.
	e.g "sp_viper.py  stack  directory  --ir=INNER_RADIUS  --radius=OUTER_RADIUS"

	Used to check the consistency between command tokens (i.e. arguments and options of the command) 
	in "===== Input =====" and "===== Usage =====" sections.
	
	The command tokens appears here must exist also in "===== Input =====" section.
	However, the order of command tokens here is not necessary to be the same and not used in GUI.

	NOTE: Unlike in "===== Input =====" section, use "--" here since "%%--%%" is not required by DokuWiki engine.

${COMMDNA_TOKEN_KEY_PREFIX}
	This can be either, no prefix (empty string), "%%--%%", or "-".
	
	GUI uses this to generate the command line.
	- If token has no prefix (empty string), it is assigned as an argument.
	- If token has prefix of "%%--%%" or "-", the parser assume this is an option.
	
	NOTE: "%%" in "%%--%%" is necessary to suppress automatic character conversion of DokuWiki.

${COMMDNA_TOKEN_KEY_BASE}
	Key base or name of a command token (i.e. arguments and options) without any prefix "--" or "-".
	e.g. "radius"

	GUI uses this to generate the command line.
	
	If key base contains the keyword defined in construct_keyword_dict() function of wikiparser.py, 
	the associated special data type will be assigned to this command token.

${COMMDNA_TOKEN_LABEL}
	User-friendly name of a command token.  
	e.g. "radius of the particle" (instead of "radius")

	Used as the label of this command token in GUI.

${COMMDNA_TOKEN_HELP}
	Help information about a command token.
	e.g. "has to be less than int(nx/2)-1" (for "radius")

	Used as the help information in GUI.

${COMMDNA_TOKEN_DEFAULT}
	Default value of a command token.

	Used as the default value of this command token in GUI.

	GUI also uses default value to find its data type (either int, float, string, bool), 
	please make sure to use appropriate expression for the value.
	For example, use "1.0" instead of "1" for float type.

	In case of using string type to express a list of int or float type, 
	please use the following formula for wikiparser.py to deduce the type from the default value:
	
	1     ==>  int
	1.0   ==>  float
	‘1.0' ==>  string. please use single quotations (') instead of double (")
	Aaa   ==>  string

	In case of string type, it can be multiple words.
	e.g. "current directory"

	Please do not use 
	- "required" fir the 1st word, 
	- "question reversed in GUI" substring anywhere
	- "value reversed in GUI" substring anywhere
	since these are keyword in Wiki document parser syntax.

	If the value in GUI is not changed from the default value, 
	GUI excludes this command token from the command line upon its generation.

${COMMDNA_TOKEN_TYPE}
	Data type of a command token. Valid types are int, float, string, and bool.
	
	If the key base contains a keyword explained above, the type provided here will be 
	override by a special data type in the Wiki document parser.

------------------------------------------------------------------------------------------
2. How to generate sp_gui.py
------------------------------------------------------------------------------------------
(1) Go to the directory contains the Wiki document parser script.

$ cd ~/EMAN2/src/eman2/sphire/templates/wikiparser.py

(2) Execute the Wiki document parser script

$ ./wikiparser.py

The script generates sp_gui_auto.py in the same directory by copying sp_gui_template.py to the sp_gui_auto.py, 
and inserting the extracted information into sp_gui_auto.py.

Please check the print out of the Wiki document parser.  
It will tell you about incorrect formats of some Wiki documents. 
If this happens, please edit the Wiki documents.

(3) Copy the generated sp_gui_auto.py to the EMAN2/bin directory of your installation as sp_gui.py.

$ cp ~/EMAN2/src/eman2/sphire/bin/sp_gui_auto.py ~/EMAN2/bin/sp_gui.py

(4) Run sp_gui.py in your project directory

$ sp_gui.py &


That's it!

==========================================================================================
END OF README
==========================================================================================

# Tutorial


## How to create the SPHIRE tutorial


1. Create a word docx document out of the latest .gfm version stored in the versions folder:  
 **make latest\_to\_docx**  
 This will create a new **tutorial\_${time\_stemp}.docx**.  

2. Edit the file following the rules listed below.

3. Rename the edited file to **tutorial.docx**.

4. Check if you did not do any mistakes in the syntax:  
 **make check\_docx**  
 This will create a new file **tutorial\_check.gfm**, that can be used to check the occuring errors.  
 It also can also be used to create a diff with the .gfm file used in step 1.

5. Fix occuring errors and repeat step 4.  
 The errors:  
 ('|||URLEND|||>', 'not found')  
 ('|||URL|||<', 'not found')  
 ('>|||URLEND|||', 'not found')  
 Can be ignored.

6. Build the PDF  
 **make clean**  
 **make**..
 This will create the **build** folder and after the program finished, please check the **build/tutorial.log** file for more advanced errors.

7. Fix the last errors and after no errors occure (*Latexmk: All targets (../build/tutorial.pdf) are up-to-date* in log.txt), distribute the PDF


## How to create this tutorial


1. Edit README.md in **git flavoured markdown** syntax:  
- **Bold** - \*\*Bold\*\*
- *italic* - \*italic\*
- ***Bolditalic*** - \*\*\*italic\*\*\*
- \_ need to be replaced by \\\_

2. Build the pdf  
 **make readme**


## Installation


1. Download and install pandoc:  
 **https://pandoc.org/**

2. Download and install TeX:  
 Linux: **https://www.tug.org/texlive/**  
 Mac: **https://www.tug.org/mactex/**

3. Run a python 2.7 environment while compiling the tutorial


## How to edit the tutorial


### Format

- There are 4 format options available:
    - Tags
    - ***Bolditalic***
    - **Bold**
    - *italic*
- Tags - For postprocessing |||TAG|||; some tags require |||TAG||| text |||TAGEND|||.
- ***Bolditalic*** - Used for button names and inline terminal commands.
- **Bold** - Used for file names, directory names, GUI labels and to highlight other texts.
- *italic* - Used for file extensions.
- All other formattings should be avoided.


### Tags

Listed are just the ones for editing existing chapters/sections.  
In case you want to see all possibilities checkout **scripts/gfm\_to\_latex\_tags.txt**.


#### Inline

- |||NEWLINE|||  
Put this one between paragraphs.  
- |||SI|||number unit unit|||SIEND|||  
Put this around numbers that have a unit, e.g. |||SI||1.14 1 per pixel per angstrom^2|||SIEND|||.  
- |||UNIT|||unit unit|||UNITEND|||  
Put this around a unit , e.g. |||UNIT||1 per pixel per angstrom^2|||UNITEND|||.  
- |||URL|||url||URLEND|||  
Put this around numbers that have a unit, e.g. |||URL||https://....|||URLEND|||.  
- |||WIKI|||  
Put this to insert a link to the SPHIRE wiki  
- |||REDNUM||| text |||REDNUMEND|||  
Color text red  
- |||PURPLENUM||| text |||PURPLENUMEND|||  
Color text purple  

#### Environments

- |||ISSUE|||{Issue header} text |||ISSUEEND|||  
Start an issue section  
- |||ADVANCED|||{Advanced header} text |||ADVANCEDEND|||  
Start an Advanced section  
- |||EQUATION|||Equation in TeX syntax|||EQUATIONEND|||  
Put this for nice formatted equations.  
- |||TIP||| text |||TIPEND|||  
Tip section  
- |||NOTE||| text |||NOTEEND|||  
Note section  
- |||IMPORTANT||| text |||IMPORTANTEND|||  
Important section  
- |||TERMINAL||| text |||TERMINALEND|||  
Terminal input section  
- |||STDOUT||| text |||STDOUTEND|||  
Standard output/log file content section  
- |||CENTER||| text |||CENTEREND|||  
Centered text

#### Itemization

- |||ITEMIZEDASH||| text |||ITEMIZEDASHEND|||  
List with dashes as labelitem  
- |||ITEMIZEARROW||| text |||ITEMIZEARROWEND|||  
List with arrows as labelitem  
- |||ITEMIZEGUI||| text |||ITEMIZEGUIEND|||  
List with dots as labelitem  
- |||ENUMERATE||| text ||ENUMERATEEND|||  
List with numbers as labelitem  

#### Inside Itemization

- |||ITEM|||  
Put this before items in an |||ITEMIZEXXX||| |||ITEM||| text |||ITEMIZEXXXEND||| and |||ENUMERATE||| |||ITEM||| text |||ENUMERATEEND||| environment.  
- |||ITEMTIP||| text |||ITEMTIPEND|||  
Tip section in an |||ITEMIZEXXX||| |||ITEM||| text |||ITEMIZEXXXEND||| and |||ENUMERATE||| |||ITEM||| text |||ENUMERATEEND||| environment.  
- |||ITEMNOTE||| text |||ITEMNOTEEND|||  
Note section in an |||ITEMIZEXXX||| |||ITEM||| text |||ITEMIZEXXXEND||| and |||ENUMERATE||| |||ITEM||| text |||ENUMERATEEND||| environment.  
- |||ITEMIMPORTANT||| text |||ITEMIMPORTANTEND|||  
Important section in an |||ITEMIZEXXX||| |||ITEM||| text |||ITEMIZEXXXEND||| and |||ENUMERATE||| |||ITEM||| text |||ENUMERATEEND||| environment.  

#### Figures

- |||FIGUREGUI|||figure\_name|||FIGUREGUIEND|||  
Large figure with a width of the page.  
- |||FIGUREMEDIUM|||figure\_name|||FIGUREMEDIUMEND|||  
Medium figure with a width of 70% of the page.  
- |||FIGURESMALL|||figure\_name|||FIGURESMALLEND|||  
Small figure with a width of 40% of the page.  
- |||FIGUREADV|||figure\_name|||FIGUREADVEND|||  
Figure used in a advanced or issue environment.  
- |||FIGURESUB|||figure\_name\_1|||FIGURESUBMIDDLE|||figure\_name\_2|||FIGURESUBEND|||  
Two figures next to each other.  
- |||FIGUREMINI|||figure\_name\_1|||FIGUREMINIMIDDLE||| text |||FIGUREMINIEND|||  
Figure and text on the same level.  
- |||FIGURECAP||| Caption  
Start for an figure caption, e.g.  
|||FIGURESUB|||figure\_name\_1|||FIGURECAP||| CAP1 |||FIGURESUBMIDDLE|||figure\_name\_2|||FIGURECAP|||CAP2|||FIGURESUBEND|||  

### Acronyms
- |||ANOVA|||  
analysis of variance  
- |||BDB|||  
Berkley Data Bank  
- |||CRYO-EM|||  
cryo-electron microscopy  
- |||CS|||  
Spherical Aberration  
- |||CTF|||  
Contrast Transfer Function  
- |||ELM|||  
Electron Microscopy  
- |||FSC|||  
Fourier shell correlation  
- |||GUI|||  
graphical user interface  
- |||HDF|||  
Hierarchical Data Format  
- |||ISAC|||  
iterative stable alignment and clustering  
- |||ML|||  
Maximum Likelihood  
- |||MPI|||  
Message Passing Interface  
- |||MTF|||  
modulation transfer function  
- |||PDB|||  
protein data bank  
- |||SD|||  
standard deviation  
- |||SNR|||  
Signal to noise ratio  
- |||SPANA|||  
single particle analysis  
- |||SPARX|||  
single particle analysis for resolution extension  
- |||SPHIRE|||  
SParx for HIgh-REsolution electron microscopy  
- |||XFEG|||  
X-Field Emission Gun  
- |||SGE|||  
Sun Grid Engine  


### Basic Latex 

#### Inside equation

- \\frac{numerator}{denominator}
- \\text{Plain text}
- \\times
- \\approx

### Outside equation

- \\parencite{citationname}  
Name of the citation located in latex/Tex\_global/lit.bib.
- \\\\  
Force line break at this point.

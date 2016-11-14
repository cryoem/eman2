#-----------------------------------------------
def process(x):
    try:
        istream=open(x,"r");
    except IOError:
        print "Cannot find file "+x
        return

    ostream=open(x+".tex","w");

    line=istream.readline();    #take out the preamble
    while((line.find("\\begin{document}")==-1)and(line!="")):
        ostream.write("%"+line);
        line=istream.readline();

    if (line==""):  #no \begin{doc} found, so reset files
        try:            
            istream.close();
            ostream.close();
            istream=open(x,"r");
            ostream=open(x+".tex","w");
            line=istream.readline();
        except IOError:
            pass

    remove=["\\title", "\maketitle", "\documentclass", "\\begin{document}",
            "\end{document}", "\\tableofcontents"]; 
        
    while(line!=""):   #remove tags we don't want
        ok=1;
        for word in remove:
            if (line.find(word)!=-1):
                ok=0;
        if (ok==1):
            ostream.write(line);
        else:
            ostream.write("%"+line);
        line=istream.readline();

    try:
        istream.close();
        ostream.close();
        print "Input File: "+x+"\t\tOutput file: "+x+".tex"
    except IOError:
        print "Had trouble closing the file "+x+" or "+x+".tex"            

#------------------------------------------------


import sys

if (len(sys.argv)<2):
    sys.stderr.write("Need a source LaTeX file");
    exit(1);

sys.argv.pop(0);
for x in sys.argv:
    process(x)

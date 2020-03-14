import os
import EMAN2


files = [
    "fonts/DejaVuSerif.ttf", 
    "images/feather.png", 
    "images/EMAN2Icon.png", 
    "images/SirEMAN2.png",
    'lib/pmconfig/icons.json', 
    'lib/pmconfig/spr.json', 
    'lib/pmconfig/tomo.json', 
    'lib/pmconfig/tomosegpanel.json', 
]

eman2dir = EMAN2.e2getinstalldir()
print("EMAN2DIR: %s" % eman2dir)

for f in files:
    filepath = os.path.abspath(os.path.join(eman2dir, f))
    print("testing for existence of file via envar EMAN2DIR: %s" % filepath)
    
    assert os.path.isfile(filepath) == 1

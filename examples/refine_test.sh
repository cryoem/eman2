#!/usr/bin/env bash

e2rawdata.py orig_micrographs/*hdf --edgenorm --xraypixel --ctfest --apix=1.275 --voltage=300.0 --cs=2.7 --ac=10.0 --threads=32 --defocusmin=0.6 --defocusmax=4.0

e2import.py orig_boxes/*box --import_boxes --box_type=boxes

e2boxer.py --allmicrographs --boxsize=288 --ptclsize=192 --write_ptcls --apix=-1.0

e2ctf_auto.py --hires --apix=1.275 --voltage=300.0 --cs=2.7 --invartype=auto --defocusmin=0.6 --defocusmax=4.0 --constbfactor=80.0 --ac=10.0 --threads=32 --minqual=0

e2refine2d_bispec.py --input=sets/all__ctf_flip_fullres.lst --ncls=96 --alignsort --normproj --nbasisfp=8 --parallel=thread:32 --center=xform.center --classkeep=0.8 --classiter=4 --classalign=rotate_translate_tree:flip=1 --classaligncmp=ccc --classraligncmp=ccc --classaverager=ctf.weight.autofilt --classcmp=ccc --classnormproc=normalize.edgemean

e2initialmodel.py --input=r2db_01/classes_sort_00.hdf --iter=12 --tries=30 --shrink=3 --sym=d2 --automaskexpand=-1 --parallel=thread:32

e2refine_easy.py --input=sets/all__ctf_flip_lp5.lst --model=initial_models/model_00_01.hdf --targetres=8.0 --speed=5 --sym=d2 --iter=4 --mass=400.0 --invar --apix=0.0 --classkeep=0.9 --m3dkeep=0.8 --parallel=thread:32 --threads=32 --automaskexpand=-1 --ampcorrect=auto

e2refine_easy.py --input=sets/all__ctf_flip_fullres.lst --model=refine_02/threed_04.hdf --targetres=4.0 --speed=5 --sym=d2 --iter=3 --mass=400.0 --apix=0.0 --classkeep=0.9 --m3dkeep=0.8 --parallel=thread:32 --threads=32 --automaskexpand=-1 --ampcorrect=auto


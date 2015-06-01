

from EMAN2 import *
import numpy as np
import shutil
import os
import pandas as pd
from pandas import DataFrame
import subprocess as sp
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
import scipy
import scipy.cluster.hierarchy as sch

class MicrographSampler:

    def __init__(self,path,boxsize,nboxes=100):
        fname = os.path.basename(path)
        name,ext = fname.split('.')
        ndir = self.mkpth(os.path.abspath(path)[:-4])
        npath = ndir+'/'+fname
        shutil.copyfile(path,npath)
        self.path = npath
        self.proc_file = self.path[:-4] + '_proc.hdf'
        self.samp_file = self.path[:-4] + '_samp.hdf'
        self.img = EMData(self.path)
        self.hdr = self.img.get_attr_dict()
        self.boxsize = int(boxsize)
        self.boxes = []
        self.nboxes = 0
        self._generated = False
        self._processed = False
        self.gen_random_boxes(nboxes)

    def __repr__(self):
        return "Micrograph(%s,boxsize=%i)"%(self.path,self.boxsize)

    def get_box_at(self,x,y,asarray=False):
        region = Region(x-self.boxsize/2,y-self.boxsize/2,self.boxsize,self.boxsize)
        if asarray: return region.numpy()
        else: return region

    def add_box_at(self,x,y,meta={}):
        region = self.get_box_at(x,y)
        self.add_box(region,meta)

    def add_box(self,region,meta={}):
        if region not in self.boxes:
            self.boxes.append({'region':region,'meta':meta})

    def gen_random_boxes(self,nsamps):
        self.nboxes = self.nboxes + nsamps
        if not self._generated:
            self._random = True
            self._generated = True
        else: print("Appending %i regions for a total of %i"%(nsamps,self.nboxes))
        xlist = np.random.random_integers(int(self.boxsize/2),int(self.hdr['nx']-self.boxsize/2),nsamps)
        ylist = np.random.random_integers(int(self.boxsize/2),int(self.hdr['ny']-self.boxsize/2),nsamps)
        for x,y in zip(xlist,ylist):
            self.add_box_at(x,y)

    def gen_exhaustive_boxes(self,fx=0,fy=0,lx=-1,ly=-1,sx=1,sy=1):
        if not self._generated:
            self._exhaustive = True
            self._generated = True
            if fx == 0: fx = int(self.boxsize/2)
            if fy == 0: fy = int(self.boxsize/2)
            if lx == -1: lx = int(self.hdr['nx']-self.boxsize/2)
            if ly == -1: ly = int(self.hdr['ny']-self.boxsize/2)
            for y in xrange(fy,ly,sy):
                for x in xrange(fx,lx,sx):
                    self.add_box_at(x,y)
            self.nboxes = len(self.boxes)
        else: print("Exhaustive regions have already been generated.")

    def display(self):
        display(self.img)

    def show_box_at(self,x,y):
        r = self.get_box_at(x,y)
        display(self.img.get_clip(r))

    def show_boxes(self,infile=None):
        if not infile: infile = self.proc_file
        sp.call(['e2display.py',infile])

    def get_box(self,box_number):
        return self.boxes[box_number]

    def get_box_coords(self,box_number):
        return self.boxes[box_number]['region'].get_origin()[:2]

    def process_boxes(self,processor,params=None,infile=None,outfile=None):
        if not self._generated:
            print("You must first generate regions randomly or exhaustively.")
            return
        if infile == None:
            if self._processed: infile = self.proc_file
            else:
                infile = self.samp_file
                self._processed = True
        # perform additional processing in place
        if not outfile: outfile = self.proc_file
        elif infile == self.proc_file: outfile = infile
        for imgnum in xrange(EMUtil.get_image_count_c(infile)):
            self.process_box(infile,imgnum,processor,params,outfile)

    @staticmethod
    def process_box(infile,imgnum,processor,params=None,outfile=None):
        if not outfile: outfile = infile # process in place
        box = EMData(infile,imgnum)
        if params: box.process_inplace(processor,params)
        else: box.process_inplace(processor)
        box.write_image_c(outfile,imgnum)

    def write_stack(self,outfile=None,imgperfile=100000):
        if not self._generated: print("You must first generate regions randomly or exhaustively.")
        if self._random:
            if not outfile: outfile = self.samp_file
            for i,b in enumerate(self.boxes):
                box = self.img.get_clip(b['region'])
                box.write_image_c(outfile,i)
        elif self._exhaustive:
            path = self.path[:-4]
            os.makedirs(path)
            for i,b in enumerate(self.boxes):
                box = self.img.get_clip(b['region'])
                outname = path + "/r%s.hdf"%(i/imgperfile)
                box.write_image_c(outname,i)

    def get_feature_matrix(self,fname,asdf=False):
        mat = []
        for i in xrange(self.nboxes):
            img = EMData(fname,i)
            mat.append(img.numpy().flatten())
        if asdf:
            try:
                import pandas as pd
                return pd.DataFrame(np.vstack(mat))
            except: return np.vstack(mat)
        else: return np.vstack(mat)

    @staticmethod
    def mkpth(path, stem=''):
        containing = os.path.dirname(os.path.realpath(path))
        contents = os.listdir(containing)
        basename = os.path.basename(path)
        if '_' not in basename: basename += '_00'
        while basename in contents:
            components=basename.split('_')
            if components[-1].isdigit(): components[-1] = str(int(components[-1])+1).zfill(2)
            else: components.append('00')
            basename = '_'.join(components)
        if basename not in contents:
            path = containing + '/' + basename
            os.makedirs(path)
        return path


mgs = MicrographSampler(path='./micrographs/dh3962.hdf',
                        boxsize=96,
                        nboxes=10000)

# add labeled background boxes
mgs.add_box_at(2606,4506,meta={'type':"background",'intensity':"darker"})
mgs.add_box_at(1724,5078,meta={'type':"background",'intensity':"lighter"})
mgs.add_box_at(669,503,meta={'type':"background",'intensity':"light"})
mgs.add_box_at(447,5688,meta={'type':"background",'intensity':"dark"})

# add labeled particle boxes
mgs.add_box_at(1938,3916,meta={'type':"particle",'view':"top"})
mgs.add_box_at(2112,4417,meta={'type':"particle",'view':"side"})
mgs.add_box_at(3882,3099,meta={'type':"particle",'view':"side"})
mgs.add_box_at(2644,2445,meta={'type':"particle",'view':"side"})

# write samples to disk
mgs.write_stack()

# process boxes. function takes care of file handling.
mgs.process_boxes('normalize.edgemean')
mgs.process_boxes('math.meanshrink',{'n':2})


#mgs.show_boxes(mgs.proc_file)
#mgs.get_box_coords(2)


# testing PCA...nothing to see here :\
A = mgs.get_feature_matrix(mgs.proc_file)
A = np.asmatrix(A.T) * np.asmatrix(A)
U, S, V = np.linalg.svd(A)
eigvals = S**2 / np.cumsum(S)[-1]
sing_vals = np.arange(A.shape[1]) + 1

ncomps=10
fig = plt.figure(figsize=(8,5))
plt.plot(sing_vals[:ncomps], eigvals[:ncomps], 'ro-', linewidth=2)
plt.title('Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Eigenvalue')
plt.show()

D = mgs.get_feature_matrix(mgs.proc_file,asdf=True)
data = DataFrame()
data['mu'] = D.mean(axis=1)
data['sigma'] = D.std(axis=1)
data['kurtosis'] = D.kurtosis(axis=1)
data['kurtosis'] = D.kurtosis(axis=1)
data['min'] = D.min(axis=1)
data['max'] = D.max(axis=1)
data['median'] = D.median(axis=1)
data['sum'] = D.sum(axis=1)
data['skewness'] = D.skew(axis=1)
#data['mode'] = D.astype(int).mode(axis=1)
data['mad'] = D.mad(axis=1)
data['var'] = D.var(axis=1)
data['sem'] = D.sem(axis=1)
data['kurt'] = D.kurt(axis=1)
data['quantile_25'] = D.quantile(0.25,axis=1)
data['quantile_50'] = D.quantile(0.50,axis=1)
data['quantile_75'] = D.quantile(0.75,axis=1)

#data = data/data.max(axis=0) #normalize wrt max
df = data.as_matrix()

# Generate random features and distance matrix.
dist = np.zeros([df.shape[0],df.shape[0]])
for i in range(df.shape[0]):
    for j in range(df.shape[0]):
        dist[i,j] = np.sum(np.abs(df[i]-df[j]))

# Compute and plot first dendrogram.
fig = plt.figure(figsize=(20,20))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
Y = sch.linkage(dist, method='centroid')
Z1 = sch.dendrogram(Y, orientation='right')
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
Y = sch.linkage(dist, method='single')
Z2 = sch.dendrogram(Y)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
dist = dist[idx1,:]
dist = dist[:,idx2]
im = axmatrix.matshow(dist, aspect='auto', origin='lower', cmap=pylab.cm.Greys)
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
plt.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('./micrographs/dendrogram.png')


#!/usr/bin/env python
# Author: Lan Dang, 02/21/2024 (dlan@bcm.edu)

import sys
import weakref
import re
from PyQt5 import QtGui, QtWidgets, QtCore, QtOpenGL
from PyQt5.QtCore import Qt
from EMAN2 import *
from EMAN2_utils import interp_points, base_name
from eman2_gui.emapplication import EMApp

try:
	from cryoet_data_portal import Client, Tomogram, Run
	print("Importing cryoet_data_portal")

except:
	print("cryoet_data_portal library required to download data from czi database. Install the library in your conda environment using cmd: \npip install -U cryoet-data-portal")
	sys.exit(0)

def main():
	usage="""a popup app to download tomogram from czi cryoet data portal and import into EMAN2 as correct format for e2tomo_annotate.py
	[prog]
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--no_gui",action="store_true", help="Run the tomogram from GUI", default=False)
	parser.add_argument("--inquire_only",action="store_true", help="Inquire the dataset only", default=False)
	parser.add_argument("--import_to_eman2",action="store_true", help="Import to EMAN2 with format compatible with e2tomo_annotate", default=False)
	parser.add_argument("--tomogram",action="store_true", help="Download tomogram", default=True)
	parser.add_argument("--annotation",action="store_true", help="Download annotation", default=False)
	parser.add_argument("--imod",action="store_true", help="Data from old imod", default=False)
	parser.add_argument("--set_indices",type=str, help="Indices of tomogram to download. Default download all tomograms in the dataset", default="all")
	parser.add_argument("--set_id",type=int, help="Data set ID to download", default=None)

	(options, args) = parser.parse_args()

	em_app = EMApp()
	czi_loader = CZIDataLoader(em_app, options)
	if not options.no_gui:
		czi_loader.show()
		x=em_app.exec_()
		sys.exit(0)

	if options.inquire_only:
		czi_loader.inquire_dataset()
	if options.tomogram or options.annotation:
		czi_loader.download_dataset()
	if options.import_to_eman2:
		czi_loader.import_data_to_eman()



class CZIDataLoader(QtWidgets.QWidget):
	def __init__(self,application,options):
		super().__init__()

		self.setWindowTitle("CZI CryoET Data Portal")
		self.setMinimumSize(400, 250)
		self.app = weakref.ref(application)
		self.options = options
		self.dataset_id = self.options.set_id
		self.set_indices = self.options.set_indices
		self.download_tomo = self.options.tomogram
		self.download_anno = self.options.annotation
		self.imod = self.options.imod

		if not self.options.no_gui:
			self.init_gui()
			return

	def init_gui(self):
		self.dataset_id_le= QtWidgets.QLineEdit()
		if self.dataset_id:
			self.dataset_id_le.setText(str(self.dataset_id))
		self.inquire_bt= QtWidgets.QPushButton("Inquire")
		self.download_tomo_cb = QtWidgets.QCheckBox("Tomogram")
		self.download_anno_cb = QtWidgets.QCheckBox("Annotation")
		self.download_tomo_cb.setChecked(self.options.tomogram)
		self.download_anno_cb.setChecked(self.options.annotation)
		self.data_download_bt= QtWidgets.QPushButton("Download")
		self.download_tomo_num= QtWidgets.QLineEdit()
		self.download_tomo_num.setText("all")
		self.imod_data_cb = QtWidgets.QCheckBox("IMOD")
		self.imod_data_cb.setChecked(self.options.imod)
		if not os.path.exists("./CZI_data"):
			os.mkdir("./CZI_data")

		gbl = QtWidgets.QGridLayout(self)
		gbl.addWidget(QtWidgets.QLabel("Dataset ID"),0,0,1,1)
		gbl.addWidget(self.dataset_id_le,0,1,1,1)
		gbl.addWidget(self.inquire_bt,1,0,1,1)
		gbl.addWidget(self.data_download_bt,2,0,1,1)
		gbl.addWidget(self.download_tomo_num,2,1,1,1)

		gbl.addWidget(self.download_tomo_cb,3,1,1,1)
		gbl.addWidget(self.download_anno_cb,3,2,1,1)
		gbl.addWidget(self.imod_data_cb,3,0,1,1)

		self.import_to_eman2_bt = QtWidgets.QPushButton("Import")
		gbl.addWidget(self.import_to_eman2_bt,4,0,1,1)
		self.import_to_eman2_bt.clicked[bool].connect(self.import_data_to_eman)


		self.inquire_bt.clicked[bool].connect(self.inquire_dataset)
		self.data_download_bt.clicked[bool].connect(self.download_dataset)
		self.download_tomo_cb.stateChanged[int].connect(self.download_tomo_cb_changed)
		self.download_anno_cb.stateChanged[int].connect(self.download_anno_cb_changed)
		self.imod_data_cb.stateChanged[int].connect(self.imod_cb_changed)


	def create_data_folder(self,set_id):
		dataset_fname ="./CZI_data/{}".format(set_id)
		if not os.path.exists(dataset_fname):
			os.mkdir(dataset_fname)
		self.tomo_fname = "./CZI_data/{}/tomos".format(set_id)
		if not os.path.exists(self.tomo_fname):
			os.mkdir(self.tomo_fname)
		if self.download_anno:
			self.seg_fname = "./CZI_data/{}/segs".format(set_id)
			if not os.path.exists(self.seg_fname):
				os.mkdir(self.seg_fname)
			return

	def download_tomo_cb_changed(self,int):
		self.download_tomo=int
	def download_anno_cb_changed(self,int):
		self.download_anno=int
	def imod_cb_changed(self,int):
		self.imod=int

	def download_dataset(self):
		self.inquire_dataset()
		self.create_data_folder(self.get_dataset_id())
		try:
			self.set_indices = self.download_tomo_num.text()
		except:
			print("No set indices provide. Download all tomogram in the dataset")
			return
		indices = []
		tomo_num= self.set_indices
		if tomo_num =="all":
			print("Download all tomograms in the current dataset")
			indices = list(range(len(self.tomos)))
		elif re.search("[0-9]+(-[0-9]+){0,1}(,[0-9]+(-[0-9]+)*)*",tomo_num).group(0) == tomo_num:
			for s in tomo_num.split(","):
				l_small  =sorted([int(index) for index in s.split("-")])
				if len(l_small) == 1:
					indices.append(int(l_small[0]))
					continue
				else:
					for i in range(l_small[0],l_small[1]+1):
						indices.append(i)
			print("Download tomograms at index", sorted(list(set(indices))))
		else:
			print(re.search("[0-9]+(-[0-9]+){0,1}(,[0-9]+(-[0-9]+)*)*",tomo_num).group(0))
			print("Please provide a valid indices range to download tomograms. \nEx: 1-3,5-7 to download tomogram 1,2,3,5,6,7; or \n 'all' to download the whole dataset")
			return

		for index in sorted(list(set(indices))):
			if index not in range(len(self.tomos)):
				print("Index", index, "is out of range. Skip")
				continue
			tomo=self.tomos[index]
			js_fname = os.path.join(self.tomo_fname,tomo.name[:-4]+".json")

			#if self.download_tomo_cb.isChecked():
			if self.download_tomo:
				with open(js_fname, "w") as f:
					json.dump(tomo.to_dict(), f)
				print("Downloading Tomogram",tomo.name)
				tomo.download_mrcfile(dest_path=self.tomo_fname)
			# if self.download_anno_cb.isChecked():
			if self.download_anno:
				seg_dest = os.path.join(self.seg_fname,tomo.name)
				os.mkdir(seg_dest)
				print("Downloading Annotation(s) for tomogram",tomo.name)
				tomo.download_all_annotations(dest_path=seg_dest,shape="SegmentationMask",format="mrc")
		print("Finish downloading data")


	def inquire_dataset(self):
		try:
			if not self.dataset_id:
				try:
					self.dataset_id = int(self.dataset_id_le.text())
				except:
					print("Please provide a valid set_id")
					return
			client = Client()
			self.tomos = Tomogram.find(
			client,
			[Tomogram.tomogram_voxel_spacing.run.dataset.id==self.dataset_id],
			)
		except Exception as e:
			print("Invalid dataset id or",e,". Abort.")
			return
		print("Dataset ID {} includes {} tomograms".format(str(self.dataset_id),str(len(self.tomos))))

	def get_dataset_id(self):
		try:
			return self.dataset_id_le.text()
		except:
			return self.dataset_id

	def import_data_to_eman(self):
		self.create_data_folder(self.get_dataset_id())
		eman2_tomo_fname ="./{}_eman".format(self.get_dataset_id())
		if not os.path.exists(eman2_tomo_fname):
			os.mkdir(eman2_tomo_fname)
		else:
			print("Data folder already exist. Return")
			self.close()
			return
		for tomo_f in os.listdir(self.tomo_fname):
			if tomo_f.endswith("mrc"):
				ori_f = os.path.join(self.tomo_fname,tomo_f)
				print("Importing", tomo_f, "from", ori_f)
				eman2_f = os.path.join(eman2_tomo_fname,tomo_f[0:-4]+".hdf")
				if self.imod:
					imod_import = ":8 --process math.fixmode:byte_utos=1"
				else:
					imod_import = ""

				print("e2proc3d.py {} {}{} --process=normalize.maxmin".format(ori_f,eman2_f,imod_import))
				os.system("e2proc3d.py {} {}{} --process=normalize.maxmin".format(ori_f,eman2_f,imod_import))

				ori_seg_fold = os.path.join("./CZI_data/{}/segs".format(self.dataset_id),tomo_f[0:-4])
				if os.path.exists(ori_seg_fold) and len(os.listdir(ori_seg_fold)) != 0:
					try:
						iter = 1
						for seg_f in os.listdir(ori_seg_fold):
							if seg_f.endswith("mrc"):
								ori_seg_f = os.path.join(ori_seg_fold,seg_f)
								print("Importing segmentations",str(iter),"for",tomo_f[0:-4])
								eman2_seg_f = os.path.join("./segs/","{}_{}_seg_from_czi.hdf".format(base_name(eman2_f),str(iter)))
								os.system("e2proc3d.py {} {}".format(ori_seg_f,eman2_seg_f))
								iter += 1
					except Exception as e:
						print(e)
						continue
		self.close()



if __name__ == '__main__':
	main()

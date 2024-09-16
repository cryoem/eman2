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
from eman2_gui.valslider import StringBox

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
		self.setWindowTitle("CZII CryoET Data Portal")
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

	def populate_table(self,data_l):
		self.data_table.setRowCount(0)
		for i,data_name in enumerate(data_l):
			self.data_table.insertRow(i)
			self.data_table.setItem(i,0,QtWidgets.QTableWidgetItem(data_name))
		print("Finish populating data table")

	def init_gui(self):
		self.dataset_id_le= QtWidgets.QLineEdit()
		if self.dataset_id:
			self.dataset_id_le.setText(str(self.dataset_id))
		self.inquire_bt= QtWidgets.QPushButton("Inquire")
		self.data_table = QtWidgets.QTableWidget(0,1)
		self.data_table.setHorizontalHeaderLabels(["Tomogram Name"])
		self.data_table.setColumnWidth(0, 240)
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
		gbl.addWidget(self.dataset_id_le,0,1,1,2)
		gbl.addWidget(self.inquire_bt,1,0,1,1)
		gbl.addWidget(self.data_table,1,1,5,2)
		gbl.addWidget(self.data_download_bt,2,0,1,1)
		gbl.addWidget(self.imod_data_cb,3,0,1,1)
		gbl.addWidget(self.download_tomo_cb,4,0,1,1)
		gbl.addWidget(self.download_anno_cb,5,0,1,1)

		self.import_to_eman2_bt = QtWidgets.QPushButton("Import")
		self.binary_label_checkbox = QtWidgets.QCheckBox("Binary")
		self.multiclass_label_checkbox = QtWidgets.QCheckBox("Multiclass")
		self.binary_label_checkbox.setChecked(True)
		self.multiclass_label_checkbox.setChecked(True)

		gbl.addWidget(self.import_to_eman2_bt,9,0,1,1)
		gbl.addWidget(self.binary_label_checkbox,9,1,1,1)
		gbl.addWidget(self.multiclass_label_checkbox,9,2,1,1)
		self.import_to_eman2_bt.clicked[bool].connect(self.import_data_to_eman)

		self.annotate_eman2_bt = QtWidgets.QPushButton("Segmentation")
		self.region_sz_sb = StringBox(label="Region Sz",value="500",showenable=-1)
		self.zthick_sb = StringBox(label="Zthick",value="-1",showenable=-1)
		self.enable_undo_checkbox = QtWidgets.QCheckBox("Enable Undo")
		self.enable_undo_checkbox.setChecked(False)

		gbl.addWidget(self.annotate_eman2_bt,10,0,1,1)
		gbl.addWidget(self.region_sz_sb,11,0,1,1)
		gbl.addWidget(self.zthick_sb,11,1,1,1)
		gbl.addWidget(self.enable_undo_checkbox,11,2,1,1)

		self.inquire_bt.clicked[bool].connect(self.inquire_dataset)
		self.data_download_bt.clicked[bool].connect(self.download_dataset)
		self.download_tomo_cb.stateChanged[int].connect(self.download_tomo_cb_changed)
		self.download_anno_cb.stateChanged[int].connect(self.download_anno_cb_changed)
		self.imod_data_cb.stateChanged[int].connect(self.imod_cb_changed)
		self.annotate_eman2_bt.clicked[bool].connect(self.launch_e2tomo_annotate)

	def show_question_box(self,msg):
		msg = QMessageBox()
		msg.setIcon(QMessageBox.Question)
		msg.setText(msg)
		msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
		return(msg)

	def overwrite_folder(self,fname):
		msg = "Data folder "+fname+" already exists.\nDo you want to overwrite it?"
		result = QtWidgets.QMessageBox.question(self,"Overwrite",msg,QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
		if result == QtWidgets.QMessageBox.Yes:
			os.system(f"rm -rf {fname}")
			print("%s folder has been removed successfully" %fname)
			os.mkdir(fname)
			return True
		else:
			print(fname+" is not overwritten. Continue.")
			return False

	def download_tomo_cb_changed(self,int):
		self.download_tomo=int
	def download_anno_cb_changed(self,int):
		self.download_anno=int
	def imod_cb_changed(self,int):
		self.imod=int

	def download_dataset(self):
		#self.inquire_dataset()
		if not self.create_data_folder(self.get_dataset_id()):
			print("Data already downloaded. Please continue with the pipeline.")
			return
		indices = [item.row() for item in self.data_table.selectedItems()]
		for index in indices:
			tomo=self.tomos[index]
			js_fname = os.path.join(self.tomo_fname,tomo.name[:-4]+".json")
			if self.download_tomo:
				with open(js_fname, "w") as f:
					json.dump(tomo.to_dict(), f)
				print("Downloading Tomogram",tomo.name)
				tomo.download_mrcfile(dest_path=self.tomo_fname)
			if self.download_anno:
				seg_dest = os.path.join(self.seg_fname,tomo.name)
				os.mkdir(seg_dest)
				try:
					print("Downloading Annotation(s) for tomogram",tomo.name)
					tomo.download_all_annotations(dest_path=seg_dest,shape="SegmentationMask",format="mrc")
				except Exception as e:
					print("Error download annotation(s) for tomogram {} due to {}".format(tomo.name,e))

		print("Finish downloading data")

	def inquire_dataset(self):
		try:
			self.dataset_id = int(self.dataset_id_le.text())
		except:
			print("Please provide a valid set_id")
			return
		try:
			client = Client()
			tomos = Tomogram.find(
			client,
			[Tomogram.tomogram_voxel_spacing.run.dataset.id==self.dataset_id],
			)
			self.tomos = [tomo for tomo in tomos]
			self.tomos.sort(key=lambda x:x.name)

		except Exception as e:
			print("Invalid dataset id or",e,". Abort.")
			return
		print("Dataset ID {} includes {} tomograms".format(str(self.dataset_id),str(len(self.tomos))))
		tomo_l = [tomo.name for tomo in self.tomos]
		self.populate_table(data_l=tomo_l)


	def get_dataset_id(self):
		try:
			return self.dataset_id_le.text()
		except:
			return self.dataset_id


	def create_data_folder(self,set_id, set_params_only=False):
		dataset_fname ="./CZI_data/{}".format(set_id)
		if not os.path.exists(dataset_fname):
			os.mkdir(dataset_fname)

		self.tomo_fname = "./CZI_data/{}/tomos".format(set_id)
		self.seg_fname = "./CZI_data/{}/segs".format(set_id)
		if set_params_only:
			return
		if self.download_tomo:
			if not os.path.exists(dataset_fname):
				os.mkdir(dataset_fname)
			if not os.path.exists(self.tomo_fname):
				os.mkdir(self.tomo_fname)
			else:
				if not self.overwrite_folder(self.tomo_fname):
					return False
		if self.download_anno:
			if not os.path.exists(dataset_fname):
				os.mkdir(dataset_fname)
			if not os.path.exists(self.seg_fname):
				os.mkdir(self.seg_fname)
			else:
				if not self.overwrite_folder(self.seg_fname):
					return False
		return True

	def import_data_to_eman(self):
		self.create_data_folder(self.get_dataset_id(),set_params_only=True)
		self.eman2_tomo_fname ="./{}_eman".format(self.get_dataset_id())
		if not os.path.exists(self.eman2_tomo_fname):
			os.mkdir(self.eman2_tomo_fname)
		else:
			if not self.overwrite_folder(self.eman2_tomo_fname):
				return
			else:
				pass
		for tomo_f in os.listdir(self.tomo_fname):
			if self.imod:
				imod_import = ":8 --process math.fixmode:byte_utos=1"
			else:
				imod_import = ""
			if tomo_f.endswith("mrc"):
				ori_f = os.path.join(self.tomo_fname,tomo_f)
				print("Importing", tomo_f, "from", ori_f)
				eman2_f = os.path.join(self.eman2_tomo_fname,tomo_f[0:-4]+".hdf")
				os.system("e2proc3d.py {} {}{} --process=normalize.maxmin".format(ori_f,eman2_f,imod_import))
				print(f"Finish importing tomogram {tomo_f} to EMAN2 project.")
				if (self.binary_label_checkbox.isChecked() or self.multiclass_label_checkbox.isChecked()):
					ori_seg_fold = os.path.join("./CZI_data/{}/segs".format(self.get_dataset_id()),tomo_f[0:-4])
					if os.path.exists(ori_seg_fold) and len(os.listdir(ori_seg_fold)) != 0:
						annf_l = []
						try:
							for seg_f in os.listdir(ori_seg_fold):
								if seg_f.endswith("mrc"):
									ori_seg_f = os.path.join(ori_seg_fold,seg_f)
									if (self.binary_label_checkbox.isChecked()):
										sname  = "_".join(base_name(seg_f).split("_")[:-1])
										print("Importing segmentations",sname,"for",tomo_f[0:-4],"as a binary mask")
										eman2_seg_f = os.path.join("./segs/","{}_{}_czi__seg.hdf".format(base_name(eman2_f),sname))
										json_file = eman2_seg_f[0:-4]+".json"
										os.system("e2proc3d.py {} {} --compressbits=8".format(ori_seg_f,eman2_seg_f))
										#os.system("e2proc3d.py {} {}".format(ori_seg_f,eman2_seg_f))

										ser_text =  json.dumps(["1",sname,"-1"], default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))
										json_str = {ser_text:None}
										js=js_open_dict(json_file)
										js['tree_dict'] = json_str
									if (self.multiclass_label_checkbox.isChecked()):
										annf_l.append(ori_seg_f)
							if len(annf_l) > 0:
								eman2_seg_f = os.path.join("./segs/","{}_{}_from_czi__seg.hdf".format(base_name(eman2_f),"000-multi"))
								print("Generating the multicolor annotation")
								ann_out, json_str = self.write_multiclass_annotate(annf_l)
								ann_out.write_compressed(eman2_seg_f,0,8)
								json_file = eman2_seg_f[0:-4]+".json"
								js=js_open_dict(json_file)
								js['tree_dict'] = json_str
								print("Done importing multiclass segmentations for",tomo_f[0:-4],"as",eman2_seg_f)
						except Exception as e:
							print(e)
							continue
				else:
					continue
		print("Finish importing data to EMAN2 project.")

	def write_multiclass_annotate(self,annf_l):
		ann_out = EMData(annf_l[0])
		json_dict = {}
		for i,annf in enumerate(annf_l):
			if i>0:
				ann = EMData(annf)
				bg_ann = 1-ann.process("threshold.binary",{"value":0.1})
				ann_out = (i+1)*ann + ann_out*bg_ann
			sname  = "_".join(base_name(os.path.basename(annf)).split("_")[:-1])
			text = [str(i+1),sname,"-1"]
			ser_text =  json.dumps(text, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))
			json_dict[ser_text] =  None
		return ann_out,json_dict

	def launch_e2tomo_annotate(self):
		self.eman2_tomo_fname ="./{}_eman".format(self.get_dataset_id())
		zthick=self.zthick_sb.getValue()
		reg_sz=self.region_sz_sb.getValue()
		tomo_fname = self.eman2_tomo_fname
		if self.enable_undo_checkbox.isChecked():
			notmp = ""
		else:
			notmp = "--no_tmp"
		os.system(f"e2tomo_annotate.py  --zthick={zthick}  {notmp} --region_sz={reg_sz} --folder={tomo_fname} &")
		self.close()


if __name__ == '__main__':
	main()

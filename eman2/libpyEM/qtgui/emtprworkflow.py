#!/usr/bin/env python
#
# Author: David Woolford 11/10/08 (woolford@bcm.edu)
# Copyright (c) 2000-2008 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#


from emsprworkflow import *
from emform import *
from emsave import EMFileTypeValidator
from emapplication import error, EMErrorMessageDisplay
from EMAN2db import db_open_dict
	
class EMBaseTomoChooseFilteredPtclsTask(WorkFlowTask):
	"""Choose the data""" 
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Choose Filtered Images To Display"
		self.preferred_size = (480,300)
		
	def get_params(self):
		ptcl_opts = EMPartSetOptions(tpr_ptcls_dict)
		self.particles_map, self.particles_name_map, choices, self.name_map = ptcl_opts.get_particle_options()
		
		#if as_string:
		#params.append(ParamDef(name="particle_set_choice",vartype="string",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short=title,property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		db = db_open_dict(self.form_db_name)
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		if len(choices) > 0:
			if len(choices) == 1:
				vartype = "choice"
			else:
				vartype = "string"
				
			params.append(ParamDef(name="tomo_filt_choice",vartype=vartype,desc_long="Choose from the filtered tomogram particles",desc_short="Choose data",property=None,defaultunits=db.get("tomo_filt_choice",dfl=""),choices=choices))				
		else:
			params.append(ParamDef(name="blurb2",vartype="text",desc_short="",desc_long="",property=None,defaultunits="There are no particles in the project. Go back to earlier stages and box/import particles",choices=None))
		return params

	def on_form_ok(self,params):
		raise NotImplementedException("Inheriting classes must implement this function")

class EMTomoChooseFilteredPtclsTask(EMBaseTomoChooseFilteredPtclsTask):
	"""Choose the particle set you wish to filter. The available sets inlcude the raw particles, and any filtered sets you have previously generated.""" 
	def __init__(self):
		EMBaseTomoChooseFilteredPtclsTask.__init__(self)
		self.form_db_name ="bdb:tomo.choose.filtered"

	def on_form_ok(self,params):
		if not params.has_key("tomo_filt_choice") or params["tomo_filt_choice"] == None:
			error("Please choose some data")
			return
		choice = params["tomo_filt_choice"]
		
		task = EMTomoGenericReportTask(self.particles_map[self.particles_name_map[choice]])
		self.emit(QtCore.SIGNAL("replace_task"),task,"Filter Tomo Particles")
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)



class E2TomoFilterParticlesTask(WorkFlowTask):	
	"""This task is for Fourier filtering and/or rotating your data. If you choose to perform both of these operations, the Fourier filtering is performed before the rotation."""
	
	preprocessor_cache = None
	def __init__(self,ptcls_list=[],name_map={}):
		WorkFlowTask.__init__(self)
		self.window_title = "Filter Tomogram Particles"
		self.output_formats = ["bdb","hdf"] # disable img from the workflow because in EMAN2 we want to store more metadata in the header
		self.form_db_name = "bdb:emform.tomo.filter_particles"
		self.project_dict = tpr_ptcls_dict
		self.ptcls_list = ptcls_list
		self.name_map = name_map
		
	def get_params(self):
		params = []
		
		params = []
		
		self.table_tool = EMTomoPtclReportTool(self.project_dict,self.window_title)
		table = self.table_tool.get_particle_table_with_ptcls(self.ptcls_list)
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(table)
		self.add_filt_params(params)
		db = db_open_dict(self.form_db_name)
		pname =  ParamDef("name",vartype="string",desc_short="Filtered set name",desc_long="The processed sets will be referred to by this name",property=None,defaultunits=db.get("name",dfl="filt"),choices=[])
		params.append(pname)
		return params
	
	def add_filt_params(self,params):
		'''
		'''
		db = db_open_dict(self.form_db_name)
		az = ParamDef(name="az",vartype="float",desc_short="Az rotation",desc_long="Rotate your model about the z axis",property=None,defaultunits=db.get("az",0.0),choices=None)
		alt = ParamDef(name="alt",vartype="float",desc_short="Alt rotation",desc_long="Rotate your model about the x axis",property=None,defaultunits=db.get("alt",0.0),choices=None)
		phi = ParamDef(name="phi",vartype="float",desc_short="Phi rotation",desc_long="Rotate your model about the z' axis",property=None,defaultunits=db.get("phi",0.0),choices=None)
		ppostproc =  ParamDef("filter",vartype="string",desc_short="Filter",desc_long="A post processor applied to the reconstructed model",property=None,defaultunits=db.get("filter",dfl="None"),choices=self.get_postprocess_filt_options())
		ppostprocargs =  ParamDef(name="filterargs",vartype="string",desc_short="params",desc_long="Parameters for the post processor see \"e2help.py processors\"",property=None,defaultunits=db.get("filterargs",dfl=""),choices=[])	
		
		params.append([ppostproc,ppostprocargs])
		params.append([az,alt,phi])
		
		
		
	def get_postprocess_filt_options(self):
		if E2TomoFilterParticlesTask.preprocessor_cache == None:
			a = dump_processors_list()
			l = ["None"]
			for key in a.keys():
				if len(key) > 5 and key[:6] == "filter":
					vals = key.split(".")
					if len(vals) > 1:
						if vals[1] in ["lowpass","highpass"]:
							l.append(key)
							
			E2TomoFilterParticlesTask.preprocessor_cache = l
			 
		return E2TomoFilterParticlesTask.preprocessor_cache
	
	def check_params(self,params):
		error_message = []
		
		if params["filter"] != "None":
			if len(params["filterargs"]) == 0:
				error_message.append("If you supply a post processor for make3d, you have to supply one of the cutoff_abs, cutoff_freq, or cutoff_pixels parameters")
			else:
				s = params["filter"] + ":"
				s += params["filterargs"]
				p = parsemodopt(s)
				if p[0] == None:
					error_message.append("Error can't interpret the make3d post processor string (%s)" %(s))
				else:
					try:
						Processors.get(p[0], p[1])
						params["filter_processor"] = s
					except:
						error_message.append("Error, can't interpret parameters for the make3d post processor (%s)" %(s))
						vals = dump_processors_list()
						values = vals[p[0]]
						s = "The parameters for the %s processor are:"  %p[0]
						
						for i in xrange(1,len(values),3):
							s += " " + values[i] +","
						s = s[:-1] # get rid of the last column
						error_message.append(s)
						
		error_message.extend(self.check_name_param(params))		
		return error_message
	
	def check_name_param(self,params):
		error_message = []
		
		if not params.has_key("name"):
			error_message.append("You must supply a name for your filtered set")
		else:
			if params["name"] in self.get_previous_filtered_set_names():
				error_message.append("There is a previously filtered set with the name %s. Please choose another name" %params["name"])
		
		return error_message
	
	def on_form_ok(self,params):	
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		error_message = self.check_params(params)
		if len(error_message) >0: 
			self.show_error_message(error_message)
			return
		
		if params["alt"] != 0 or params["az"] != 0 or params["phi"] != 0:
			params["rotate"] = "%.2f,%.2f,%.2f" %(params["az"],params["alt"],params["phi"])

		if params.has_key("rotate") or params.has_key("filter_processor"):
			success,cmd = self.process_all(params)
			if not success:
				error("Command failed:"+cmd)
				return
		else:
			error("You have to supply a filter or a non zero rotation for any filtering to occur")
			return
		
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
	
		self.write_db_entries(params)	
	
	def get_previous_filtered_set_names(self):
		data_dict = EMProjectDataDict(tpr_ptcls_dict)
		db_map = data_dict.get_data_dict()
		
#		project_db = db_open_dict("bdb:project")
#		name = tpr_ptcls_dict
#		if not project_db.has_key(name):
#			return []
#		
#		db_map = project_db.get(name)
		previous_sets = []
		
		for root_name,dict in db_map.items():
			for filt,name in dict.items():
				if 	previous_sets.count(filt) == 0:
					previous_sets.append(filt)
					
		return previous_sets
	
	def convert_to_root(self,name):
		if self.name_map.has_key(name): return self.name_map[name]
		else:return name
	
	def output_names(self,params):
		'''
		Only works if you the params dictionary has the filenames and name keys 
		'''
		return ["bdb:tomo_particles#"+base_name(self.convert_to_root(f))+"_"+params["name"] for f in params["filenames"]]
		
	def process_all(self,params):
		
		outnames = self.output_names(params)
		
		progress = QtGui.QProgressDialog("Processing files...", "Abort", 0, len(params["filenames"]),None)
		progress.show()
	
		i = 0
		for i,name in enumerate(params["filenames"]):
			cmd = "e2proc3d.py"
 			cmd += " "+name
 			cmd += " "+outnames[i]
 			if params.has_key("filter_processor"):
 				cmd += " --process="+params["filter_processor"]
 			if params.has_key("rotate"):
 				cmd += " --rot="+params["rotate"]
 			success = (os.system(cmd) in (0,12))
 			if not success:
 				progress.close()
 				return False,cmd
 			else:
 				progress.setValue(i+1)
 				get_application().processEvents()
		
		progress.close()
		
		self.save_to_filt_ptcls_map(params, outnames)
		return True,"Success"
	
	def save_to_filt_ptcls_map(self,params,outnames):
		data_dict = EMProjectDataDict(tpr_ptcls_dict)
		db_map = data_dict.get_data_dict()
		
#		project_db = db_open_dict("")
#		project_list = tpr_ptcls_dict
#		db_map = project_db.get(project_list,dfl={})
		
		for i,name in enumerate(params["filenames"]):
			real_name = self.convert_to_root(name)
			if db_map.has_key(real_name):
				d = db_map[real_name]
				d[params["name"]] = outnames[i]
				db_map[real_name] = d
			else:
				d = {}
				d[params["name"]] = outnames[i]
				db_map[real_name] = d
		data_dict.update(db_map)



class EMTomoChooseFilteredPtclsForFiltTask(EMBaseTomoChooseFilteredPtclsTask):
	"""Choose the data you wish to filter""" 
	def __init__(self,task_type=E2TomoFilterParticlesTask):
		EMBaseTomoChooseFilteredPtclsTask.__init__(self)
		self.form_db_name ="bdb:tomo.choose.filtered.forfilt"
		self.task_type = task_type

	def on_form_ok(self,params):
		if not params.has_key("tomo_filt_choice") or params["tomo_filt_choice"] == None:
			error("Please choose some data")
			return
		choice = params["tomo_filt_choice"]
		
		task = self.task_type(self.particles_map[self.particles_name_map[choice]],self.name_map)
		self.emit(QtCore.SIGNAL("replace_task"),task,"Filter Tomo Particles")
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)

class EMTomoBootStapChoosePtclsTask(EMBaseTomoChooseFilteredPtclsTask):
	"""Choose the particle set you wish to use to generate the bootstrapped probe. The sets available will include the raw particles and any filtered sets you have generated.""" 

	def __init__(self):
		EMBaseTomoChooseFilteredPtclsTask.__init__(self)
		self.form_db_name ="bdb:tomo.choose.forbootstraptomo"

	def on_form_ok(self,params):
		if not params.has_key("tomo_filt_choice") or params["tomo_filt_choice"] == None:
			error("Please choose some data")
			return
		choice = params["tomo_filt_choice"]
		task = EMTomoBootstrapTask(self.particles_map[self.particles_name_map[choice]],self.name_map)
		self.emit(QtCore.SIGNAL("replace_task"),task,"Filter Tomo Particles")
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)

class EMTomoAlignParams:
	'''
	A class that get parameters commonly used by alignment based programs
	'''
class EMSubTomoDataReportTask(EMRawDataReportTask):
	'''Gets subtomogram stacks. ''' 
	def __init__(self):
		EMRawDataReportTask.__init__(self)
		self.project_list = tpr_subtomo_stacks
		
	def get_raw_data_table(self):
		'''
		Gets an EMTomographicFileTable - this is type of class that the emform knows how to handle 
		'''
		data_dict = EMProjectDataDict(self.project_list)
		project_data = data_dict.get_data_dict()
		project_names = project_data.keys()
		self.project_data_at_init = project_data # so if the user hits cancel this can be reset

		from emform import EMTomographicFileTable,EMFileTable
		table = EMTomographicFileTable(project_names,desc_short="Sub Tomograms",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(self.project_list)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
	
		#p.append(pdims) # don't think this is really necessary
		return table,len(project_names)

class EMRefDataReportTask(EMRawDataReportTask):
	'''Gets subtomogram references. ''' 
	def __init__(self):
		EMRawDataReportTask.__init__(self)
		self.project_list = tpr_subtomo_ref
		
	def get_raw_data_table(self):
		'''
		Gets an EMTomographicFileTable - this is type of class that the emform knows how to handle 
		'''
		data_dict = EMProjectDataDict(self.project_list)
		project_data = data_dict.get_data_dict()
		project_names = project_data.keys()
		self.project_data_at_init = project_data # so if the user hits cancel this can be reset

		from emform import EMTomographicFileTable,EMFileTable
		table = EMTomographicFileTable(project_names,name="refnames", desc_short="Sub Tomogram References",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(self.project_list)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
	
		#p.append(pdims) # don't think this is really necessary
		return table,len(project_names)	

class EMTomoBootstrapTask(WorkFlowTask):
	'''Use this task for running e2spt_classaverage'''
	
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.tomo_boxer_module = None
		self.form_db_name = "bdb:emform.tomo.classavg3d"
		self.window_title = "Launch e2spt_classaverage"
		self.report_task = None
		
	def get_tomo_hunter_basic_table(self):
		'''
		'''
		self.report_task = EMSubTomoDataReportTask()
		table,n = self.report_task.get_raw_data_table()# now p is a EMParamTable with rows for as many files as there in the project
		from emform import EMFileTable,int_lt
		return table, n
		
	def get_tomo_hunter_ref_table(self):
		'''
		'''
		self.report_task = EMRefDataReportTask()
		table,n = self.report_task.get_raw_data_table()# now p is a EMParamTable with rows for as many files as there in the project
		from emform import EMFileTable,int_lt
		return table, n
		
	def run_form(self):
		self.form = EMTableFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)	
		
	def get_params(self):
		table_params = []
		params = []
		db = db_open_dict(self.form_db_name)
		
		p,n = self.get_tomo_hunter_basic_table() # note n is unused, it's a refactoring residual		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of tomohunter",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		psavesteps = ParamDef(name="savesteps",vartype="boolean",desc_short="Savesteps",desc_long="Save the steps",property=None,defaultunits=db.get("savesteps",dfl=True),choices=None)
		psaveali = ParamDef(name="saveali",vartype="boolean",desc_short="Saveali",desc_long="Save the alignments",property=None,defaultunits=db.get("saveali",dfl=True),choices=None)
		params.append([psavesteps, psaveali])
		piter = ParamDef(name="number of iterations",vartype="int",desc_short="Number of iterations", desc_long="Number of iterations",property=None,defaultunits=db.get("number of iterations",dfl=5),choices=None )
		pncoarse = ParamDef(name="coarse number",vartype="int",desc_short="Coarse Number", desc_long="Coarse number",property=None,defaultunits=db.get("coarse number",dfl=6),choices=None )
		params.append([piter, pncoarse])
		pshrink = ParamDef(name="Percentage to shrink",vartype="int",desc_short="Shrink", desc_long="Percentage to shrink",property=None,defaultunits=db.get("Percentage to shrink",dfl=2),choices=None )
		pshrinkrefine = ParamDef(name="Percentage to shrink, refinement",vartype="int",desc_short="Shrink refine", desc_long="Percentage to shrink for refienment",property=None,defaultunits=db.get("Percentage to shrink, refinement",dfl=2),choices=None )
		params.append([pshrink, pshrinkrefine])
		
		proc_data = dump_processors_list()
		masks = {}
		for key in proc_data.keys():
			if len(key) >= 5 and key[:5] == "mask.":
				masks[key] = proc_data[key]
		masks["None"] = ["Choose this to stop masking from occuring"]
		pmask = ParamDef("mask",vartype="string",desc_short="Mask",desc_long="The mask to apply to the subtomos",property=None,defaultunits=db.get("mask",dfl="None"),choices=masks)
		pmaskparams = ParamDef("maskparams",vartype="string",desc_short="Params",desc_long="Parameters for the mask",property=None,defaultunits=db.get("maskparams",dfl=""))
		params.append([pmask, pmaskparams])
		
		filters = {}
		for key in proc_data.keys():
			if len(key) >= 7 and key[:7] == "filter.":
				filters[key] = proc_data[key]
		filters["None"] = ["Choose this to stop filtering from occuring"]
		pfilter = ParamDef("filter",vartype="string",desc_short="Filter",desc_long="The Filter to apply to the subtomos",property=None,defaultunits=db.get("filter",dfl="None"),choices=filters)
		pfilterparams = ParamDef("filterparams",vartype="string",desc_short="Params",desc_long="Parameters for the filter",property=None,defaultunits=db.get("filterparams",dfl=""))
		params.append([pfilter, pfilterparams])

		ali_data = dump_aligners_list()
		caligners = {}
		for key in ali_data.keys():
			if len(key) >= 19 and key[:19] == "rotate_translate_3d":
				caligners[key] = ali_data[key]
		pali = ParamDef("aligner3D",vartype="string",desc_short="Aligner3D",desc_long="The 3D course aligner",property=None,defaultunits=db.get("aligner3D",dfl="rotate_translate_3d"),choices=caligners)
		paliparams = ParamDef("ali3dparams",vartype="string",desc_short="Params",desc_long="Parameters for the 3D aligner",property=None,defaultunits=db.get("ali3dparams",dfl="search=10:delta=15:dphi=15:verbose=1"))
		params.append([pali, paliparams])
		
		craligners = {}
		for key in ali_data.keys():
			if len(key) >= 9 and key[:9] == "refine_3d":
				craligners[key] = ali_data[key]
		prali = ParamDef("raligner3D",vartype="string",desc_short="RAligner3D",desc_long="The 3D refine aligner",property=None,defaultunits=db.get("raligner3D",dfl="refine_3d_grid"),choices=craligners)
		praliparams = ParamDef("rali3dparams",vartype="string",desc_short="Params",desc_long="Parameters for the 3D refine aligner",property=None,defaultunits=db.get("rali3dparams",dfl="verbose=1"))
		params.append([prali, praliparams])
		
		ppostfilter = ParamDef("postfilter",vartype="string",desc_short="PostFilter",desc_long="The Filter to apply to the average",property=None,defaultunits=db.get("postfilter",dfl="None"),choices=filters)
		ppostfilterparams = ParamDef("postfilterparams",vartype="string",desc_short="Params",desc_long="Parameters for the postfilter",property=None,defaultunits=db.get("postfilterparams",dfl=""))
		params.append([ppostfilter, ppostfilterparams])
		
		pparallel = ParamDef("parallel",vartype="string",desc_short="Parallel",desc_long="Parallalization parameters",property=None,defaultunits=db.get("parallel",dfl=""))
		params.append(pparallel)
			
		#pylong = ParamDef(name="yshort",vartype="boolean",desc_short="yshort",desc_long="Use Z axis as normal",property=None,defaultunits=1,choices=None)
		#pinmem = ParamDef(name="inmemory",vartype="boolean",desc_short="inmemory",desc_long="Load the tomo into memory",property=None,defaultunits=1,choices=None)
		#papix = ParamDef(name="apix",vartype="float",desc_short=u"\u212B per pixel", desc_long="Angstroms per pixel",property=None,defaultunits=1.0,choices=None )
		#params.append([pylong, pinmem, papix])
		#params.append(papix)
		table_params.append(["Main",params])
		
		advanced_params = []
		r,rn = self.get_tomo_hunter_ref_table() # note rn is unused, it's a refactoring residual
		advanced_params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use tomo refs",desc_long="",property=None,defaultunits="Use this for selecting references. If no reference is chosen, then reference free alignment will be executed",choices=None))
		advanced_params.append(r)
		
		table_params.append(["References",advanced_params])
#		db = db_open_dict(self.form_db_name)
#		params.append(ParamDef(name="interface_boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=db.get("interface_boxsize",dfl=128),choices=[]))
#		#db_close_dict(self.form_db_name)
		return table_params
		
	def on_form_ok(self,params):
		
		if not params.has_key("filenames"):
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return
		if len(params["refnames"]) > 1:
			EMErrorMessageDisplay.run(["Only one referecne can be used"])
			return
		
		e23dcalist = "e2spt_classaverage.py"
		e23dcalist += " --input="+params['filenames'][0]
		if(params['filenames'][0][-4:-3] == "."):
			e23dcalist += " --output="+params['filenames'][0][:-4]+"_3DAVG."+params['filenames'][0][-3:] # output file hack
		else:
			e23dcalist += " --output="+params['filenames'][0]+"_3DAVG"
		
		spacer = ""
		try:
			e23dcalist += " --ref="+params["refnames"][0]
		except:
			pass
		if params["savesteps"]:
			e23dcalist += " --savesteps"
		if params["saveali"]:
			e23dcalist += " --saveali"
		if params["number of iterations"]:
			e23dcalist += " --iter="+str(params["number of iterations"])
		if params["Percentage to shrink"]:
			e23dcalist += " --shrink="+str(params["Percentage to shrink"])
		if params["Percentage to shrink, refinement"]:
			e23dcalist += " --shrinkrefine="+str(params["Percentage to shrink, refinement"])
		if params["coarse number"]:
			e23dcalist += " --npeakstorefine="+str(params["coarse number"])
		if params["mask"] != "None":
			spacer=""
			if params["maskparams"]: spacer=":"
			e23dcalist += " --mask="+params["mask"]+spacer+params["maskparams"]
		if params["filter"] != "None":
			spacer=""
			if params["filterparams"]: spacer=":"  
			e23dcalist += " --preprocess="+params["filter"]+spacer+params["filterparams"]
		spacer=""
		if params["ali3dparams"]: spacer=":" 
		e23dcalist += " --align="+params["aligner3D"]+spacer+params["ali3dparams"]
		spacer=""
		if params["rali3dparams"]: spacer=":" 
		e23dcalist += " --ralign="+params["raligner3D"]+spacer+params["rali3dparams"]
		if params["postfilter"] != "None":
			if params["postfilterparams"]: spacer=":" 
			e23dcalist += " --postprocess="+params["postfilter"]+spacer+params["postfilterparams"]
		if params["parallel"]:
			e23dcalist += " --parallel="+params["parallel"]
		print e23dcalist
		
		child = subprocess.Popen(e23dcalist, shell=True)
		
		self.form.close()
		self.form = None
		self.write_db_entries(params)

class EMTomoRawDataReportTask(EMRawDataReportTask):
	'''The tools under this tab are highly experimental. ''' 
	def __init__(self):
		EMRawDataReportTask.__init__(self)
		self.project_list = tpr_raw_data_dict
		
	def get_raw_data_table(self):
		'''
		Gets an EMTomographicFileTable - this is type of class that the emform knows how to handle 
		'''
		data_dict = EMProjectDataDict(self.project_list)
		project_data = data_dict.get_data_dict()
		project_names = project_data.keys()
		self.project_data_at_init = project_data # so if the user hits cancel this can be reset

		from emform import EMTomographicFileTable,EMFileTable
		table = EMTomographicFileTable(project_names,desc_short="Tomograms",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(self.project_list)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
	
		#p.append(pdims) # don't think this is really necessary
		return table,len(project_names)
		
class E2TomoBoxerGuiTask(WorkFlowTask):
	"""Select the file you want to process and hit okay, this will launch e2spt_boxer. The yshort option sets the Z axis normal to the screen, and inmemory load the tomo into memory for fast access"""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.tomo_boxer_module = None
		self.form_db_name = "bdb:emform.tomo.boxer"
		self.window_title = "Launch e2spt_boxer"
		self.report_task = None
		
	def get_tomo_boxer_basic_table(self):
		'''
		'''
		
		self.report_task = EMTomoRawDataReportTask()
		table,n = self.report_task.get_raw_data_table()# now p is a EMParamTable with rows for as many files as there in the project
		from emform import EMFileTable,int_lt
		table.insert_column_data(0,EMFileTable.EMColumnData("Stored Boxes",E2TomoBoxerGuiTask.get_tomo_boxes_in_database,"Boxes currently stored in the EMAN2 database",int_lt))
		
		return table, n

	def get_tomo_boxes_in_database(name):
		print "checking for boxes, but this aspect of things is not working yet...."+base_name(name)+" "+name
		#from e2spt_boxer import tomo_db_name
		#if db_check_dict(tomo_db_name):
			#tomo_db = db_open_dict(tomo_db_name)
			#image_dict = tomo_db.get(base_name(name),dfl={})
			#if image_dict.has_key("coords"):
				#return str(len(image_dict["coords"]))
		
		return "0"
	
	get_tomo_boxes_in_database = staticmethod(get_tomo_boxes_in_database)
	
	def get_params(self):
		params = []
		db = db_open_dict(self.form_db_name)
		
		p,n = self.get_tomo_boxer_basic_table() # note n is unused, it's a refactoring residual		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of e2spt_boxer",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		pylong = ParamDef(name="yshort",vartype="boolean",desc_short="yshort",desc_long="Use Z axis as normal",property=None,defaultunits=db.get("yshort",dfl=True),choices=None)
		pinmem = ParamDef(name="inmemory",vartype="boolean",desc_short="inmemory",desc_long="Load the tomo into memory",property=None,defaultunits=db.get("inmemory",dfl=True),choices=None)
		papix = ParamDef(name="apix",vartype="float",desc_short=u"\u212B per pixel", desc_long="Angstroms per pixel",property=None,defaultunits=db.get("apix",dfl=1.0),choices=None )
		params.append([pylong, pinmem, papix])
#		db = db_open_dict(self.form_db_name)
#		params.append(ParamDef(name="interface_boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=db.get("interface_boxsize",dfl=128),choices=[]))
#		#db_close_dict(self.form_db_name)
		return params
	
	def on_form_ok(self,params):
		
		if not params.has_key("filenames"):
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return

		self.write_db_entries(params)
		
		e2tblist = "e2spt_boxer.py"
		e2tblist += " "+params['filenames'][0]
		if params["yshort"]:
			e2tblist += " --yshort"
		if params["inmemory"]:
			e2tblist += " --inmemory"
		e2tblist += " --apix="+str(params["apix"])
		
		child = subprocess.Popen(e2tblist, shell=True)
		
		self.form.close()
		self.form = None
		self.write_db_entries(params)
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.tomo_boxer_module == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass
	
	def on_boxer_closed(self): 
		if self.tomo_boxer_module != None:
			self.tomo_boxer_module = None
			self.emit(QtCore.SIGNAL("gui_exit"))
	
	def on_boxer_idle(self):
		'''
		Presently this means boxer did stuff but never opened any guis, so it's safe just to emit the signal
		'''
		self.tomo_boxer_module = None
		self.emit(QtCore.SIGNAL("gui_exit"))
		
	def on_form_cancel(self):
		if self.report_task:
			self.report_task.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))

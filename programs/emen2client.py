#!/usr/bin/env python
import time
import os
import sys
import subprocess
import optparse
import getpass
import glob
import collections
import operator
import fnmatch
import struct
import copy
import shutil
import urllib
import urllib2
import urlparse
import httplib
import xmlrpclib
import gzip
import zlib
import math
import StringIO
import socket
import tempfile


VERSION = 20100928
USER_AGENT = "emen2client/%s"%VERSION


try:
	import json
except ImportError:
	try:
		import simplejson as json
	except ImportError:
		print """Note: No JSON module found. This is not necessary, but is recommended. Please upgrade to Python 2.6 or install simplejson, using "sudo easy_install simplejson" """
		json = None


try:
	import EMAN2
	import hashlib
except ImportError:
	# print "EMAN2 support not available"
	EMAN2 = None



##################################
# Compression schemes
##################################

# Compression is handled on the fly now during upload / download with a gzip pipe

# Based on http://code.activestate.com/recipes/551784/
class GzipPipe(StringIO.StringIO):
	"""This class implements a compression pipe suitable for asynchronous 
	process."""

	# Size of the internal buffer
	CHUNCK_SIZE = 1024*1024
	COMPRESSLEVEL = 3

	def __init__(self, filename=None) :
		"""Streaming compression using gzip.GzipFile
		
		@param filename source file
		
		"""
		
		self.filename = filename
		self.source = open(self.filename, "rb")
		self.source_eof = False
		self.buffer = ""
		self.name = self.filename + ".gz"

		self.pos = 0
		
		StringIO.StringIO.__init__(self)
		#super(GzipPipe, self).__init__(self)
		self.zipfile = gzip.GzipFile(filename=str(os.path.basename(filename)), mode='wb', compresslevel=self.COMPRESSLEVEL, fileobj=self)
	
	
	def write(self, data):
		"""Writes data to internal buffer. Do not call from outside."""
		self.buffer += data
		

	def read(self, size=-1) :
		"""Calling read() on a zip pipe will suck data from the source stream.
		@param	size Maximum size to read - Read whole compressed file if not specified.
		@return Compressed data
		"""

		# Feed the zipped buffer by writing source data to the zip stream
		while ((len(self.buffer) < size) or (size == -1)) and not self.source_eof:
		   
			# Get a chunk of source data
			chunk = self.source.read(GzipPipe.CHUNCK_SIZE)
			self.pos += len(chunk)
			
			# Feed the source zip file (that fills the compressed buffer)
			self.zipfile.write(chunk)
			
			# End of source file ?
			if (len(chunk) < GzipPipe.CHUNCK_SIZE) :
				self.source_eof = True
				self.zipfile.flush()
				self.zipfile.close()
				self.source.close()
				break


		# We have enough data in the buffer (or source file is EOF): Give it to the output
		if size == 0:
			result = ""
		if size >= 1:
			result = self.buffer[0:size]
			self.buffer = self.buffer[size:]
		else: # size < 0 : All requested
			result = self.buffer
			self.buffer = ""

		return result



##################################
# Exceptions
##################################
	
class OverwriteError(Exception):
	"""File exists, and flags are set to prevent overwriting"""
		
	

##################################
# Options
##################################

class DBRPCOptions(optparse.OptionParser):
	"""Default options to be used by EMEN2 clients"""

	def __init__(self, *args, **kwargs):
		#super(DBRPCOptions, self).__init__(self, *args, **kwargs)
		optparse.OptionParser.__init__(self, *args, **kwargs)

		self.add_option("--username","-U", type="string", help="Username")
		self.add_option("--password","-P", type="string", help="Password (Note: specifying passwords in shell commands is not secure)")

		self.add_option("--host","-H",type="string",help="Host endpoint URI",default="http://ncmidb.bcm.edu")
		self.add_option("--ctxid","-C",type="string",help="Valid Context ID")

	#def load_config(self):
	#	"""Load from config file"""
	#	pass



##################################
# Main exported DB interface classes
##################################


# Todo: switch to JSON-RPC which allows for positional/named arguments


class DBJSONRPCProxy(object):
	def __init__(self, host=None, username=None, password=None, ctxid=None, version=None):
		pass
		
		
		

class DBXMLRPCProxy(xmlrpclib.ServerProxy):
	"""Provide a more feature-rich DB XMLRPC interface than the bare python xmlrpclib"""

	def __init__(self, host=None, username=None, password=None, ctxid=None, version=None):
		#super(DBXMLRPCProxy, self).__init__(host, allow_none=True)
		xmlrpclib.ServerProxy.__init__(self, host, allow_none=True)

		self.__username = username
		self.__password = password

		self._ServerProxy__allow_none = True
		self._host = host
		self._handler = self._ServerProxy__handler

		self._setContext(ctxid)


	def __wrap_request(self, method, params):
		try:
			ret = self._ServerProxy__request(method, params)
		except xmlrpclib.Error, err:
			#raise xmlrpclib.Fault(0, "Error: %s"%(err.headers.get('X-Error')))
			raise ValueError, str(err.headers.get('X-Error'))

		except Exception, e:
			raise
			
		return ret


	def __getattr__(self, name):
		# print "-> ", name
		return xmlrpclib._Method(self.__wrap_request, name)
		
		
	def _login(self, username=None, password=None):
		ctxid = self.__wrap_request("_login", (username or self.__username, password or self.__password))
		self.__password = None
		return self._setContext(ctxid)


	def _setContext(self, ctxid):
		# ian: todo: make handler a dynamically generated property? maybe not..
		self._ctxid = ctxid
		if self._ctxid != None:
			self._ServerProxy__handler = "%s?ctxid=%s"%(self._handler, self._ctxid)
		return self._ctxid


	def _getContext():
		return self._ctxid




##################################
# File Transport
##################################

		
class FileTransport(object):
	"""Transfer files from/to DB"""

	def __init__(self, db=None, bdo=None, bdoid=None, param="file_binary", recid=None, newrecord=None, compress=True, filename=None, filedata=None, report=None, pos=None):
		"""Needs a BDO or BDOID for download, or a filename and recid | newrecord for upload
		
		@keyparam db DBProxy
		@keyparam bdo BDO object
		@keyparam bdoid BDO Id
		@keyparam filename Filename to upload
		@keyparam filedata Upload a string buffer instead of file reference. In this case, filename will be used for the name of the file.
		@keyparam report Upload/download status report hook/callback
		@keyparam compress Compress file
		
		"""
		
		self.db = db
		self.bdo = bdo
		self.bdoid = bdoid
					
		self.filename = filename
		self.filedata = filedata
		
		self.newrecord = newrecord
		self.param = param
		self.recid = recid
		
		self.compress = compress
		
		self._report = report
		self.pos = pos
		
		
	def _retryhandler(self, method, *args, **kwargs):
		"""Automatic retry wrapper
		@param method self._upload or self._download
		"""
		
		count = 0
		ret = None
		while True:
			count += 1
			try:
				ret = method(*args, **kwargs)
				break
			
			except KeyboardInterrupt, e:
				print "Keyboard Interrupt"
				raise
			
			# socket.error is recoverable, others are not recoverable!	
			except socket.error, e:
				self.report("Network error, trying again in 10 seconds (%s attempts): %s"%(count, e), i=1)
				time.sleep(10)

			except Exception, e:
				self.report("Unrecoverable error, aborting: %s"%(e), i=1)
				raise
				
		return ret		


	def _downloadretry(self, *args, **kwargs):
		return self._retryhandler(self._download, *args, **kwargs)


	def _uploadretry(self, *args, **kwargs):
		return self._retryhandler(self._upload, *args, **kwargs)


	def report(self, msg=None, progress=None, pos=None, i=0):
		try:
			self._report(msg=msg, progress=progress, f=self.filename or self.bdo.get("name"), pos=pos, i=i+1)
		except:
			pass
			

	def sidecar_read(self, filename=None):
		filename = filename or self.filename
		try:
			return json.load(file(filename+".json","r"))
		except:
			pass


	def sidecar_write(self, filename=None):
		filename = filename or self.filename
		try: 
			json.dump(self.bdo, file(filename+".json", "w"), indent=True)
			self.report("Wrote sidecar: %s.json"%filename, i=1)
		except:
			pass


	def setdb(self, db):
		self.db = db
		
		
	def setreport(self, report):
		self._report = report
		
		
	def setrecmap(self, recmap):
		self.recid = recmap.get(self.recid, self.recid)
			
		if self.newrecord:
			self.newrecord["children"] = [recmap.get(i,i) for i in self.newrecord.get("children",[])]
			self.newrecord["parents"] = [recmap.get(i,i) for i in self.newrecord.get("parents",[])]
		
			

	########################
	# Download methods
	########################

	def download(self, overwrite=False, rename=False, retry=False, sidecar=False):
		"""Download BDO to local disk
		
		@keyparam overwrite Overwrite an existing file
		@keyparam rename If file exists, rename file to avoid overwriting
		@keyparam retry	If there is a network error, try to download again

		"""
		
		self.report(pos=self.pos)
		
		# self._download_checkbdo()

		dl = self._download
		if retry:
			dl = self._downloadretry

		try:
			# Check that we can download the file, and then get it
			self._download_checkfile(overwrite=overwrite, rename=rename)
			self.report("Downloading %s -> %s..."%(self.bdo.get("filename"), self.filename), i=1)
			dl()
			if sidecar:
				self.sidecar_write()

		except OverwriteError, e:
			# Or report that we can't download it
			self.report("Skipping file: %s"%e, i=1)
			return
		
		
		return self.filename
		
			
			

	def _download_checkfile(self, overwrite=False, rename=False):
		"""Check if we can write the file, and if to use compression or not. See download() for params"""

		filename = os.path.basename(self.bdo.get("filename")) or self.bdo.get("name")		

		fsplit = os.path.splitext(filename)
		if self.compress and fsplit[-1] == ".gz":
			compress = True
			filename = fsplit[0]
		else:
			compress = False

		# print "reading sidecar file"
		# print self.sidecar_read(filename)

		# this may raise OverwriteError Exception, which should be caught by caller
		if os.path.exists(filename) and not overwrite:
			if rename:
				filename = "duplicate.%s:%s"%(self.bdo.get("recid"), filename)
			else:
				raise OverwriteError, "File exists: %s"%filename

		self.filename = filename
		self.compress = compress



	def _download(self):
		"""Actual download."""
		
		# Setup connection
		ctxid = urllib.urlencode({"ctxid":self.db._ctxid})

		# ian: changed .netloc to [1] because it caused failure on python 2.4
		host = urlparse.urlparse(self.db._host)[1] # .netloc
		path = "/download/%s?%s"%(self.bdo.get("name"), ctxid)		

		# Connect
		http = httplib.HTTPConnection(host)		
		http.request("GET", path)
		resp = http.getresponse()
		clength = float(resp.getheader('content-length'))

		# Download and pipe through gzip
		stdout = None
		if self.compress:
			stdout = open(self.filename, "wb")
			gz = subprocess.Popen(["gzip","-d"], stdout=stdout, stdin=subprocess.PIPE)
			fio = gz.stdin
		else:
			fio = open(self.filename,"wb")
						
		chunksize = 1024*1024
		outlen = 0
		while outlen < clength:
			chunk = resp.read(chunksize)
			outlen += len(chunk)
			fio.write(chunk)
			self.report(progress=(outlen/clength))


		fio.close()
		if stdout: stdout.close()
		http.close()
		
		return self.filename
		
	


	########################
	# Upload methods
	########################

	def upload(self, retry=False, sidecar=False): #, recid=None, newrecord=None, param='file_binary_image', compress=True, retry=True
		"""Handler provides records and files to put in DB

		@keyparam retry Enable automatic upload retry
		
		"""
		
		ul = self._upload
		if retry:
			ul = self._uploadretry

		self.report(pos=self.pos)

		s = self.sidecar_read()
		if s != None:
			self.report("File already exists in database. Check record %s."%s.get("recid"), i=1)			
			self.bdo = s
			return

		ret = ul()

		if sidecar:
			self.sidecar_write()
		
	
	
	def _upload(self, mimetype=None):
		"""Implements Transfer-Encoding:Chunked over a PUT request, with optional gzip pipe streaming. See upload() for arguments."""
		
		# Upload a string buffer as a file
		if self.filedata:
			filesize = len(self.filedata)
			fio = StringIO.StringIO(self.filedata)
			filename = self.filename
			
			
		# Or upload a file from disk	
		else:
			# Get file size for progress bar
			filesize = float(os.stat(self.filename).st_size)

			# Are we compressing this file?
			if self.compress:
				fio = GzipPipe(filename=self.filename)
			else:
				fio = open(self.filename, "rb")

			filename = os.path.basename(fio.name)


		if filesize == 0:
			self.report("Warning: Uploading empty file!")
			filesize = 1


		# Prepare connection info
		# recid = recid or self.newrecord.get('recid')
		qs = {
			"param": self.param,
			"ctxid": self.db._ctxid,
			"filename": filename
		}
		

		# If rec is set, with json support, use one-step newrecord/putbinary
		if self.newrecord:						
			if json:
				qs["newrecord"] = json.dumps(self.newrecord)
				self.recid = None
			else:
				self.report("Warning: JSON support is strongly recommended for upload; using workaround..")
				self.newrecord = self.db.putrecord(self.newrecord)
				self.recid = self.newrecord.get("recid")				
				self.report("... got recid %s"%self.recid)
				

		if self.newrecord == None and self.recid == None:
			raise ValueError, "No new record or target recid specified for upload"
			
		# ._put
		qs = urllib.urlencode(qs)
		host = self.db._host.split("://")[-1]
		path = "/upload/%s?%s"%(self.recid or "", qs)
						
		# Set headers			
		headers = {"User-Agent":USER_AGENT}
		headers["Transfer-Encoding"] = "chunked"
		if self.compress:
			mimetype = "application/x-gzip"
		if mimetype:
			headers["Content-Type"] = mimetype
			

		# Open connection
		http = httplib.HTTPConnection(host)
		http.request("PUT", path, headers=headers)			
			
		t = time.time()	
			
		# Upload in chunks
		chunksize = 1024*1024
		
		while True:
			try:
				chunk = fio.read(chunksize)
				cs = len(chunk)
				http.send("%X\r\n"%(cs))
				http.send(chunk)
				http.send("\r\n")			
				self.report(progress=(fio.tell()/filesize))

			except socket.error:
				if fio: fio.close()
				raise

			if not chunk:
				if fio: fio.close()
				break	

		# Responses..
		resp = http.getresponse()
		status, reason, response = resp.status, resp.reason, resp.read()
		http.close()
				
		if status not in range(200, 300):
			raise httplib.HTTPException, "Bad response code: %s: %s"%(status, reason)

		kbsec = (filesize / (time.time() - t))/1024

		if json:
			response = json.loads(response)
			self.bdo = response
			self.report("Done. Uploaded %s to record %s @ %0.2f kb/sec"%(self.filename, self.bdo.get("recid"), kbsec), i=1)
		else:
			self.report("Done. Uploaded %s @ %0.2f kb/sec"%(self.filename, kbsec), i=1)

		return response
		







##################################
# Base FileHandler
##################################




class FileHandler(object):
	
	request_params = []
	
	def __init__(self, db=None, recid=None, filename=None, filedata=None, rectype=None, applyparams=None, pos=None, options=None, **kwargs):
		"""Provide special support for different filetypes
		
		@keyparam db DBProxy
		@keyparam recid Record ID for download or upload
		@keyparam filename Filename for upload
		@keyparam rectype Record type for default-setting record creation
		@keyparam applyparams Parameters to overlay on new records (see: CCDHandler)
		@keyparam pos Display this FileHandler's position in the overall queue, tuple (position, total)
		@keyparam options Options. See UploadController and DownloadController OptionParsers for details.
		
		"""
		
		# @keyparam recurse For default behavior download, level of recursion from recid to look for files
		# @keyparam overwrite Overwrite files when downloading (default=False)
		# @keyparam rename Rename downloaded files to avoid existing files (default=False). If overwrite and rename are false, and a file exists, OverwriteError will be raised
		# @keyparam compress Enable compression for upload/download
		# @keyparam retry In the event of network error, keep retrying the upload/download action
		# @keyparam sidecar Write a sidecar of the BDO after upload or download
		# @keyparam force Force upload even if sidecar is found

		
		self.db = db
		
		# target recid for upload or download
		self.recid = recid
		
		# filename for upload..
		self.filename = filename
		self.filedata = filedata
		
		self.bdos = []

		# new records to apply to uploads
		self.rectype = rectype
		self.applyparams = applyparams or {}

		# default and passed options
		self.options = {
			"recurse":50,
			"overwrite": False,
			"rename": False,
			"compress": True,
			"retry": True,
			"sidecar": False
		}
		if options:
			self.options.update(options) 
		if kwargs:
			self.options.update(kwargs)


		# internal attributes
		self.tmpfiles = []
		self.checked = False
		self.quiet = False
		self.pos = pos
		self.progress = 0.0
		self.listeners = []
		
		self.init()
		

	def init(self):
		pass


	def newrecordinit(self, newrecord):
		"""Initialize and validate a new record"""
		
		tmprec = self.db.newrecord(newrecord.get('rectype'), self.recid)
		tmprec.update(newrecord)
		self.db.validaterecord(tmprec)
		return tmprec


	######################
	# Upload
	######################				
		
	def upload(self, retry=None, sidecar=None):
		"""Start the upload process for this handler"""
		
		if retry == None:
			retry = self.options.get("retry")

		if sidecar == None:
			sidecar = self.options.get("sidecar")

		ret = {}
		recmap = {}
		
		
		self.report(pos=self.pos, f=self.filename)
		self.report("Preparing for upload", i=1)

		self.prepare_upload()
		newrecords, files = self.get_upload_items()

		self.report("Checking and committing records", i=1)


		# Check records, call newrecord to properly init, validate, set permissions, etc.
		# If a particular record doesn't validate, will throw exception
		newrecords = [self.newrecordinit(i) for i in newrecords]
		for f in files:
			if f.newrecord:
				f.newrecord = self.newrecordinit(newrecord=f.newrecord)								
								
		# ian: killed a very ugly block of code.. this is much simpler
				
		# Commit any non-single-item files
		if newrecords:
			records = self.db.putrecord(newrecords)	
			recmap = dict(zip([i['recid'] for i in newrecords], [i['recid'] for i in records]))

		
		# Upload other files/newrecords
		filecount = len(files)
		for count, f in enumerate(files):
			f.setdb(self.db)
			f.setreport(self.report)
			f.setrecmap(recmap)

			if filecount > 1:
				f.pos = (count, filecount)

			f.upload(retry=retry, sidecar=sidecar)
			ret[f.filename] = f.bdo


		self.upload_postprocess(ret)
		self.cleanup()	
		return ret
		

	def get_upload_items(self):
		"""Returns records and files to upload"""
		f = FileTransport(recid=self.recid, filename=self.filename)
		return [], [f]


	def prepare_upload(self):
		"""Do any required prep work before files can be uploaded"""
		self.check_db()

	def upload_postprocess(self, ret):
		pass


	######################
	# Download
	######################

	def download(self, overwrite=None, rename=None, compress=None, retry=None, sidecar=None):
		"""Start the download process for this handler"""
		
		if overwrite == None:
			overwrite = self.options.get("overwrite")

		if rename == None:
			rename = self.options.get("rename")

		if compress == None:
			compress = self.options.get("compress")

		if retry == None:
			retry = self.options.get("retry")

		if sidecar == None:
			sidecar = self.options.get("sidecar")
		
		
		self.report(pos=self.pos, f=self.recid)
		self.report("Checking for items to download", i=1)
		
		recs, bdos = self.get_download_items()

		bdolen = len(bdos)

		if bdos:
			self.report("Found %s items in %s records"%(bdolen, len(recs)), i=1)
		else:
			self.report("Found no downloadable items in %s records"%len(recs), i=1)
			return			
			
		ret = {}
		
		for count, bdo in enumerate(bdos):
			u = FileTransport(db=self.db, bdo=bdo, report=self.report, compress=self.options.get("compress"), pos=(count, bdolen))
			u.download(overwrite=overwrite, rename=rename, retry=retry, sidecar=sidecar)
			ret[u.filename] = u.bdo
						
		self.download_postprocess(recs, ret)
				
		return ret


	def download_postprocess(self, recs, ret):
		pass
		
	

	def get_download_items(self):
		"""Find child recods of self.recid and get their BDOs"""
		
		#recs = self.db.getchildren(self.recid, "record", self.options.get("recurse"), None, True, True)		
		recs = self.db.getchildren(self.recid, self.options.get("recurse"))
		recs.append(self.recid)
		bdos = self.db.getbinary(recs)
		return recs, bdos

	

	########################
	# Status updates
	########################

	def addlistener(self, callback):
		self.listeners.append(callback)


	def report(self, msg=None, progress=None, f=None, pos=None, i=0):

		for listener in self.listeners:
			listener(msg=msg, progress=progress, f=f)
		
		if self.quiet:
			return
	
		# Normally suppress this...
		if progress != None and msg == None:
			return	

		if not msg:		
			if pos:
				print "\n%s%s of %s: %s"%("\t"*i, pos[0]+1, pos[1], f)
			else:
				print "\n%s%s"%("\t"*i, f)
		
		else:
			print "%s%s"%("\t"*i, msg)
											
													
	########################
	# Utility methods
	########################
		

	def _shortsleep(self):
		time.sleep(10)


	def check_db(self):
		pass


	########################
	# Overrideable methods for upload control
	########################


	def cleanup(self):
		for i in self.tmpfiles:
			try:
				self.report("Cleanup temp file: %s"%i, i=1)
				os.unlink(i)
			except Exception, e:
				self.report("Warning: Could not remove %s: %s"%(i, e), i=1)






#################################
# Filetype-specific handlers
##################################


class NewRecordFileHandler(FileHandler):

	def get_upload_items(self):
		newrecord = {}
		newrecord["recid"] = -1
		newrecord["rectype"] = self.rectype
		newrecord["parents"] = [self.recid]
		return [], [FileTransport(newrecord=newrecord, filename=self.filename, param="file_binary_image")]


class VolumeHandler(NewRecordFileHandler):

	def get_upload_items(self):
		newrecord = {}
		newrecord["recid"] = -1
		newrecord["title_structure"] = self.filename
		newrecord["rectype"] = self.rectype
		newrecord["parents"] = [self.recid]
		return [], [FileTransport(newrecord=newrecord, filename=self.filename, param="file_binary_image")]



class BoxHandler(FileHandler):
	def get_upload_items(self):
		newrecord = {}
		newrecord["recid"] = -1
		newrecord["rectype"] = "box"
		newrecord["parents"] = [self.recid]
		newrecord.update(self.applyparams)
		return [], [FileTransport(newrecord=newrecord, filename=self.filename, filedata=self.filedata, compress=False, param="box")]



class ParticleHandler(FileHandler):
	def get_upload_items(self):
		newrecord = {}
		newrecord["recid"] = -1
		newrecord["rectype"] = "particles"
		newrecord["parents"] = [self.recid]
		return [], [FileTransport(newrecord=newrecord, filename=self.filename, compress=False, param="ptcl")]



class GridImagingFileHandler(FileHandler):
	
	def check_db(self):		
		gi = self.db.getrecord(self.recid)

		if gi.get("rectype") != "grid_imaging":
			raise Exception, "This action may only be used with grid_imaging sessions!"

		microscopy = self.db.getparents(self.recid, 1, ["microscopy"])
		if not microscopy:
			self.report("\tWARNING! No microscopy record present for grid_imaging session!")
			#raise Exception, "No microscopy record present for grid_imaging session!"



class CCDHandler(GridImagingFileHandler):

	request_params = [
		'tem_magnification_set',
		'ctf_defocus_set',
		'angstroms_per_pixel',
		'time_exposure_tem',
		'tem_dose_rate',
		'current_screen',
		'ice_type',
		'assess_ice_thick',
		'assess_ice_comments',
		'assess_image_quality',
		'ccd_id',
		'status_energy_filter',
		'binning'
		]


	def get_upload_items(self):
		newrecord = {}
		newrecord["recid"] = -1
		newrecord["rectype"] = "ccd"
		newrecord["id_ccd_frame"] = os.path.basename(self.filename)
		newrecord["assess_image_quality"] = 5
		newrecord.update(self.applyparams)
		newrecord["parents"] = [self.recid]
		# read metadata
		return [], [FileTransport(newrecord=newrecord, filename=self.filename, param="file_binary_image")]




class ScanHandler(GridImagingFileHandler):
	pass



class StackHandler(GridImagingFileHandler):

	# See notes at bottom of file
	header_labels = [
		# X,Y,Z size
		['stack_size_nx', 'i', 0],
		['stack_size_ny', 'i', 0],
		['stack_size_nz', 'i', 0],
		## Mode: mode, 0 = unsigned byte, 1 = short int, 2 = float, 3 = short*2, 4 = float*2
		['stack_data_mode', 'i', 0],
		## Starting point of sub image: nxstart, nystart, nzstart
		['stack_start_nx', 'i', 0],
		['stack_start_ny', 'i', 0],
		['stack_start_nz', 'i', 0],
		## Grid Size: mx, my, mz
		['stack_size_mx', 'i', 0],
		['stack_size_my', 'i', 0],
		['stack_size_mz', 'i', 0],
		## Cell size: xlen, ylen, zlen... pixel spacing = xlen/mx
		['stack_size_xlen', 'i', 0],
		['stack_size_ylen', 'i', 0],
		['stack_size_zlen', 'i', 0],
		## Cell angles: alpha, beta, gamma
		['stack_angle_alpha', 'f', 0.0],
		['stack_angle_beta', 'f', 0.0],
		['stack_angle_gamma', 'f', 0.0],
		## Map column. Ignored by IMOD: mapc, mapr, maps. (1=x, 2=y, 3=z)
		['stack_map_mapc', 'i', 0],
		['stack_map_mapr', 'i', 0],
		['stack_map_maps', 'i', 0],
		## Minimum pixel value, amin
		['stack_pixel_min', 'f', 0.0],
		## Maximum pixel value, amax
		['stack_pixel_max', 'f', 0.0],
		## Mean pixel value, amean
		['stack_pixel_mean', 'f', 0.0],
		## Image type, ispg
		['stack_data_ispg', 'h', 0],
		## Space group number, nsymbt
		['stack_data_nsymbt', 'h', 0],
		## Number of bytes in extended header, next
		['stack_data_extheadersize', 'i', 0],
		## Creator ID, creatid
		['stack_data_creatid', 'h', 0],
		# 30 bytes of unused data: See explanation above
		['stack_unused', '30s', ''],
		## Number of bytes per section (SerialEM format), nint
		['stack_data_bytespersection', 'h', 0],
		## Flags for which types of short data (SerialEM format), nreal. (See above documentation.)
		['stack_data_extheaderflags', 'h', 0],
		# 28 bytes of unused data: See explanation above
		['stack_unused2', '28s', ''],
		['stack_data_idtype', 'h', 0],
		['stack_data_lens', 'h', 0],
		['stack_data_nd1', 'h', 0],
		['stack_data_nd2', 'h', 0],
		['stack_data_vd1', 'h', 0],
		['stack_data_vd2', 'h', 0],
		['stack_data_tiltangles_orig', 'fff', [0.0, 0.0, 0.0]],
		['stack_data_tiltangles_current', 'fff', [0.0, 0.0, 0.0]],
		['stack_data_xorg', 'f', 0.0],
		['stack_data_yorg', 'f', 0.0],
		['stack_data_zorg', 'f', 0.0],
		['stack_data_cmap', '4s', ''],
		# big/little endian flag
		['stack_data_stamp', '4s', ''],
		['stack_data_rms', 'f', 0.0],
		['stack_data_nlabl', 'i', 0],
		# spliced seperately: see after
		['stack_data_labels', '80s'*10, '']
	]


	extheader_flags = {
		1: {
			'pack': 'h',
			'load': lambda x:[x[0] / 100.0],
			'dump': lambda x:[x[0] * 100],
			'dest': ['specimen_tilt'],
			},

		2: {
			'pack': 'hhh',
			'load': lambda x:[x],
			'dump': lambda x:[x],
			'dest': ['stack_data_montage']
			},

		4: {
			'pack': 'hh',
			'load': lambda x:[x[0] / 25.0 , x[1] / 25.0],
			'dump': lambda x:[x[0] * 25   , x[1] * 25],
			'dest': ['stack_stagepos_x', 'stack_stagepos_y']
			},

		8: {
			'pack': 'h',
			'load': lambda x:[x[0] / 10.0],
			'dump': lambda x:[x[0] * 10],
			'dest': ['tem_magnification_set']
			},

		16: {
			'pack': 'h',
			'load': lambda x:[x[0] / 25000.0],
			'dump': lambda x:[x[0] * 25000],
			'dest': ['stack_intensity']
			},

		32: {
			'pack': 'hh',
			'load': lambda x:[sign(x[0])*(math.fabs(x[0])*256+(math.fabs(x[1])%256))*2**(sign(x[1])*(int(math.fabs(x[1]))/256))],
			'dump': lambda x:[0,0],
			'dest': ['stack_dose']
			}
	}



	def download_postprocess(self, recids, ret):

		rec = self.db.getrecord(self.recid)
		print "Rebuilding .st file. It is generally a good idea to run this download in a new directory."
		
		workfile = rec.get('stack_filename')	
		if os.access(workfile, os.F_OK):
			raise ValueError, "File %s already exists -- not going to overwrite. Exiting."
		
		f = open(workfile, "wb")
		
		print "Making header"
		a = self._make_header(rec)
		f.write(a)

		print "Making extheader"
		slices = self.db.getrecord(sorted(self.db.getchildren(self.recid, 1, "stackimage"), reverse=True))


		hl = rec.get('stack_data_extheaderflags')
		hs = rec.get('stack_data_extheadersize')
		b = self._make_extheader(slices, stack_data_extheaderflags=hl)
		padsize = hs - len(b)
		f.write(b)
		print "\t...size was %s, padding by %s to %s / %s"%(len(b), padsize, len(b)+padsize, hs)
		f.write("\0"*padsize)

		print "Adding images to stack"
		
		reti = {}
		for k,v in ret.items():
			reti[v.get("recid")] = k
			
		for rec in slices:
			fn = reti.get(rec.get("recid"))
			print fn
			fnopen = open(fn, "rb")
			fnopen.seek(1024)
			data = fnopen.read()
			f.write(data)
			print "\t...wrote %s bytes"%len(data)
			fnopen.close()
			os.unlink(fn)

		f.close()
		
		

	def init(self):
		pass


	def prepare_upload(self):
		try:
			self.header = self.get_header()
			self.slices = self.get_extheader()
			self.update_maxangles()
			self.header['stack_filename'] = self.filename
		except Exception, inst:
			raise ValueError, "This does not appear to be a SerielEM stack: %s"%inst


		# Extract headers for backup, unstack file and compress, add to tmpfiles to be removed after upload
		self.check_db()
		
		self.tmpdir = tempfile.mkdtemp()
		
		self.raw_headers = self.backup_headers()
		self.tmpfiles.extend(self.raw_headers)

		self.unstacked_files = self.unstack()
		self.tmpfiles.extend(self.unstacked_files)


	def get_upload_items(self):
		"""Get basic records and files to put in database. These will be filled out more fully with db.newrecord, etc."""

		stackrec = copy.copy(self.header)
		stackrec["recid"] = -100
		stackrec["parents"] = [self.recid]
		# print self.raw_headers
		files = [FileTransport(recid=-100, filename=i, compress=False) for i in self.raw_headers]

		for offset, (s, s_filename) in enumerate(zip(self.slices, self.unstacked_files)):
			s["recid"] = -(offset+101)
			s["parents"] = [-100]
			files.append(FileTransport(newrecord=s, filename=s_filename, param="file_binary_image"))

		return [stackrec], files


	########
	# Stack specific methods, mostly header related

	def update_maxangles(self):
		if self.slices[0].has_key("specimen_tilt"):
			tilts = [i['specimen_tilt'] for i in self.slices]

		self.header["stack_maxangle"] = max(tilts)
		self.header["stack_minangle"] = min(tilts)


	def get_header(self):
		"""Read an IMOD/SerialEM/MRC stack header"""
		f = open(self.filename,"rb")
		h = f.read(1024)
		f.close()

		return self._handle_header(h)


	def get_extheader(self):
		if not self.header:
			self.header = self.get_header()

		if not self.header["stack_data_extheadersize"]:
			return []

		f = open(self.filename, "rb")
		h = f.read(1024)
		eh = f.read(self.header["stack_data_extheadersize"])
		f.close()

		return self._handle_extheader(eh)


	def backup_headers(self):
		fin = open(self.filename, 'rb')
		raw_header = fin.read(1024)
		raw_extheader = fin.read(self.header["stack_data_extheadersize"])
		fin.close()

		ret = []

		header_outfile = "%s/%s.header"%(self.tmpdir, os.path.basename(self.filename))

		fout = open(header_outfile, 'wb')
		fout.write(raw_header)
		fout.close()
		ret.append(header_outfile)

		if raw_extheader:
			extheader_outfile = "%s/%s.header"%(self.tmpdir, os.path.basename(self.filename))
			fout = open(extheader_outfile, 'wb')
			fout.write(raw_extheader)
			fout.close()
			ret.append(extheader_outfile)

		return ret


	def unstack(self, slices=None):
		"""Unstack slices into new files.
		@keyparam slices Single slice or list of slices to unstack; default is all.
		@return Filenames of unstacked files
		"""

		raw_header = self._make_unstacked_header()
		header_offset, size = self._get_offset_and_size()
		fstack = open(self.filename, "rb")
				
		ret = []

		if slices and not hasattr(slices, "__iter__"):
			slices = [slices]
		if not slices:
			slices = range(len(self.slices))

		for n in slices:
			read_offset = header_offset + (size*n)
			
			fout_name = "%s/%s.%s.mrc"%(self.tmpdir, os.path.basename(self.filename), n)

			self.report("Unstacking frame %s into %s"%(n, fout_name), i=1)
			
			# check_filenames = [fout_name]
			# check_filenames.append(fout_name+".gz")
			# if not filter(None, map(os.path.exists, check_filenames)):

			#ret.append(fout_name)
			#continue

			fout = open(fout_name, "wb")
			fout.write(raw_header)
			fstack.seek(read_offset)

			blen = 1024*1024*8
			for buffer_size in [blen]*(size/blen)+[size%blen]:
				copy_buffer = fstack.read(buffer_size)
				fout.write(copy_buffer)

			fout.close()

			ret.append(fout_name)

		fstack.close()
		return ret


	def stack(self, outfile=None):
		pass


	def _get_offset_and_size(self):
		# 4    int     mode;      Types of pixel in image.  Values used by IMOD:
		# 	 0 = unsigned bytes,
		# 	 1 = signed short integers (16 bits),
		# 	 2 = float,
		# 	 3 = short * 2, (used for complex data)
		# 	 4 = float * 2, (used for complex data)
		#                          6 = unsigned 16-bit integers (non-standard)
		# 	16 = unsigned char * 3 (for rgb data, non-standard)

		header_offset = 1024 + self.header["stack_data_extheadersize"]
		depth = 2
		dmode = self.header['stack_data_mode']
		if dmode in [0]:
			depth = 1
		elif dmode in [2, 3]:
			depth = 4
		elif dmode in [4]:
			depth = 8
		elif dmode in [16]:
			depth = 3

		slicesize = self.header["stack_size_nx"] * self.header["stack_size_ny"] * depth

		return header_offset, slicesize


	def _handle_header(self, h):
		"""Extract data from header string (1024 bytes) and process"""

		d={}
		offset = 0

		for dest, format, default in self.header_labels:
			size = struct.calcsize(format)
			value = struct.unpack(format, h[offset:offset+size])
			if len(value) == 1: d[dest] = value[0]
			else: d[dest] = value
			offset += size

		# Process data labels
		n = d["stack_data_nlabl"]
		d["stack_data_labels"] = [i.strip() for i in d["stack_data_labels"][:n]]+[""]*(10-n)

		# Include rectype
		d["rectype"] = "stack"
		d["recid"] = -1

		# This guy gives us trouble
		d['stack_data_stamp'] = d['stack_data_stamp'].partition('\x00')[0]

		# Delete unused items
		try:
			del d["stack_unused"]
			del d["stack_unused2"]
		except:
			pass

		return d


	def _make_unstacked_header(self):

		h = copy.deepcopy(self.header)
		for i in ['stack_size_nz', 'stack_size_zlen']:
			h[i] = 1
		for i in ['stack_data_extheadersize', 'stack_data_bytespersection', 'stack_data_extheaderflags']:
			h[i] = 0

		rawheader = self._make_header(h)
		return rawheader


	def _make_header(self, header):
		rawheader = []

		header["stack_unused"] = ""
		header["stack_unused2"] = ""

		for dest, format, default in self.header_labels:
			value = header.get(dest, default)
			size = struct.calcsize(format)
			if hasattr(value,"__iter__"):
				value = struct.pack(format, *value)
			else:
				value = struct.pack(format, value)
			# print dest, format, header.get(dest), value
			rawheader.append(value)

		return "".join(rawheader)


	def _make_extheader(self, slices, stack_data_extheaderflags=0):
		eh = []

		# ian: todo: calculate keys based on info in DB instead of extheaderflags
		keys = self._extheader_getkeys(stack_data_extheaderflags)

		for s in slices:
			for key in keys:
				_k = self.extheader_flags.get(key)

				value = [s[i] for i in _k['dest']]
				value = _k['dump'](value)
				if hasattr(value, "__iter__"):
					value = struct.pack(_k['pack'], *value)
				else:
					value = struct.pack(_k['pack'], value)

				# print key, _k, [s[i] for i in _k['dest']]
				eh.append(value)

		return "".join(eh)


	def _handle_extheader(self, eh):
		"""Process extended header"""

		d = self.header
		ed = []
		keys = self._extheader_getkeys(d["stack_data_extheaderflags"])

		offset = 0
		for i in range(0, d["stack_size_nz"]):
			sslice = {}
			for key in keys:
				_k = self.extheader_flags.get(key)
				size = struct.calcsize(_k['pack'])
				# print "Consuming %s bytes (%s:%s) for %s"%(size, i+offset, i+offset+size, _k['dest'])
				value = struct.unpack(_k['pack'], eh[offset: offset+size])
				value = _k['load'](value)

				for _i in zip(_k['dest'], value):
					sslice[_i[0]] = _i[1]

				offset += size

			sslice["rectype"] = "stackimage"

			ed.append(sslice)

		return ed


	def _extheader_getkeys(self, flags):
		keys = []
		for i in sorted(self.extheader_flags.keys()):
			if flags & i:
				keys.append(i)
		return keys





##################################
# Application controllers
##################################


class EMEN2ClientController(object):
	def __init__(self, args=None, options=None, version=None):
		"""EMEN2 Client. DBRPCProxy available as 'db' attr after login."""
		self.version = version
		self.options = options
		self.args = None
		self.parser = None

		self.db = None

		self.setparser()
		if args != None:
			self.parse_args(args)


	def setparser(self):
		self.parser = DBRPCOptions()
		self.setparser_add_option()


	def setparser_add_option(self):
		pass


	def parse_args(self, inargs):
		self.options, self.args = self.parser.parse_args(inargs)
		self.check_args()


	def check_args(self):
		pass


	def run(self):
		pass


	def login(self):
		self.options.username = self.options.username or raw_input("Username: ")
		self.options.password = self.options.password or getpass.getpass("Password: ")

		self.db = DBXMLRPCProxy(host=self.options.host)
		self.checkversion()
		try:
			self.ctxid = self.db._login(username=self.options.username, password=self.options.password)
		except Exception, e:
			print "Unable to login: %s"%(e)
			sys.exit(0)
	

	def checkversion(self):
		VERSION_CLIENT = self.db.checkversion("emen2client")

		if self.version < VERSION_CLIENT:
			print """
Note: emen2client version %s is available; installed version is %s.

EMAN2 now includes emen2client.py, so you may download a new nightly build of EMAN2,
	http://ncmi.bcm.edu/ncmi/

...or visit the EMEN2 wiki and download emen2client.py directly:
	http://blake.grid.bcm.edu/emanwiki/EMEN2/emen2client.py

"""%(VERSION_CLIENT, self.version)

		# else:
		# 	print "emen2client version %s is up to date"%(self.version)




class UploadController(EMEN2ClientController):
	def setparser_add_option(self):
		self.parser.add_option("--noninteractive","-q", action="store_true", help="Do not prompt for parameter values")
 		self.parser.add_option("--sidecar", "-s", action="store_true", help="Write sidecar file after upload", default=True)
 		self.parser.add_option("--force", "-f", action="store_true", help="Force re-upload even if a sidecar is found")
		self.parser.add_option("--metafile", action="store_true", dest="metafile", help="Attempt to read JAMES/JADAS metadata files (default)", default=True)
		self.parser.add_option("--no-metafile", action="store_false", dest="metafile", help="Ignore metadata files")

		usage = """%prog upload [options] <record type> <recid> <files to upload>

Record type can be any valid database protocol.
Some record types have special, application-specific handlers, e.g.:

	ccd			CCD Frames
	scan		Scanned micrographs
	stack		Tomograms

Other values, e.g. "volume", will create child records of that type, with 1 file per child record.

Alternatively, you can use "none" for record type and the files will be attached directly to the specified record ID.
		
		"""

		self.parser.set_usage(usage)



	def _get_param_values(self, params):

		ret = {}
		if self.options.noninteractive:
			return ret

		pds = self.db.getparamdef(params)
		pds = dict([[i.get('name'), i] for i in pds])


		for param in params:
			pd = pds.get(param)
			ret[param] = self._get_param_value_console(pd)
		return ret


	#########################
	# Ugly console input code
	#########################

	def _get_param_value_console(self, param):
		name = param.get('name')

		print "\n----- %s -----"%name
		print "\tDescription: %s"%param.get('desc_long')
		if param.get('defaultunits'):
			print "\tUnits: %s"%param.get('defaultunits')


		punits = param.get('defaultunits') or ''
		lines = None
		selections = []

		try:
			data = {"ctxid": self.db._ctxid, "limit":10, "showchoices":0}		#,
			f = urllib2.urlopen('%s/find/value/%s/'%(self.db._host, name), urllib.urlencode(data))
			lines = f.read()
			f.close()
		except:
			pass

		count = 0
		if param.get('choices'):
			print "\n\tSuggested values:"
			for choice in param.get('choices', []):
				print "\t\t%s) %s"%(count, choice)
				selections.append(choice)
				count += 1

		if param.get('vartype') != "choice" and lines:
			if param.get('choices'):
				print "\n\tOther common values:"
			else:
				print "\n\tCommon values:"

			lines = lines.split("\n")
			for line in lines:
				print "\t\t%s) %s"%(count, line)
				selections.append(line)
				count += 1

		print "\n\t\t%s) None or N/A"%count
		selections.append('')
		count += 1

		if param.get('vartype') != "choice":
			print "\t\t%s) Enter a different not listed above"%count


		while True:
			inp = raw_input("\n\tSelection (0-%s): "%count)
			try:
				inp = int(inp)

				if inp == count:
					inp = raw_input("\tValue: ")
				else:
					inp = selections[inp]
			except:
				continue

			break

		return inp


	def upload(self):
		if self.action not in map_rectype_handler:
			try:
				rd = self.db.getrecorddef(self.action)
			except:
				print "Invalid record type: %s"%self.action
				return


		print ""
		print "%s Files to upload:"%len(self.filenames)
		for i in self.filenames:
			print "\t",i

		handlerclass = gethandler(rectype=self.action) # or NewRecordFileHandler
		applyparams = self._get_param_values(handlerclass.request_params)

		# There isn't a dict(options)...
		options = {}
		for k,v in self.options.__dict__.items():
			options[k]=v


		dbts = []
		filecount = len(self.filenames)
		for count, filename in enumerate(self.filenames):
			dbt = handlerclass(db=self.db, rectype=self.action, filename=filename, recid=self.recid, applyparams=applyparams, pos=(count,filecount), options=options)
			dbt.upload()



	def check_args(self):
		self.action = self.args[0]

		try:
			self.recid = int(self.args[1])
		except:
			raise Exception, "Record ID required as argument after upload action"


		self.filenames = reduce(operator.concat, map(glob.glob, self.args[2:]))

		if not self.filenames:
			raise Exception, "No files specified"


	def run(self):
		self.login()
		self.upload()






class DownloadController(EMEN2ClientController):

	def setparser_add_option(self):
		self.parser.add_option("--recurse", type="int",help="Recursion level",default=3)
		self.parser.add_option("--overwrite", "-o", action="store_true", help="Overwrite existing files (default is to skip)", default=False)
		self.parser.add_option("--rename", "-r", action="store_true", help="If a file already exists, save with format 'duplicate.recid:original_filename.dm3'", default=False)
		self.parser.add_option("--sidecar", "-s", action="store_true", help="Include sidecar file with EMEN2 metadata in JSON format")
		self.parser.add_option("--gzip", action="store_true", dest="compress", help="Decompress gzip'd files. Requires gzip in path (default)", default=True)
		self.parser.add_option("--no-gzip", action="store_false", dest="compress", help="Do not decompress gzip'd files.")
		
		usage = """%prog download [options] <recid> [filename-pattern]"""

		self.parser.set_usage(usage)
		
		# self.parser.add_option("--convert","-c",action="store_true",help="Convert files to mrc",default=0)
		# self.parser.add_option("--invert","-i",action="store_true",help="Invert contrast",default=0)


	def check_args(self):

		if self.options.sidecar and not json:
			print "Warning: sidecar support requires JSON"
			self.options.sidecar = False

		if len(self.args) < 1:
			raise  optparse.OptionError, "Record ID Required"

		try:
			self.recid = int(self.args[0])
			self.patterns = self.args[1:] or []
		except Exception, inst:
			raise  optparse.OptionError, "Record ID Required"


	def run(self):
		self.check_args()
		self.login()

		options = {}
		for k,v in self.options.__dict__.items():
			options[k] = v

		rec = self.db.getrecord(self.recid)
		handler = gethandler(rectype=rec.get('rectype'))

		dbt = handler(db=self.db, recid=self.recid, options=options, pos=(0,1))
		dbt.download()

		return







class SyncController(EMEN2ClientController):
	
	def setparser_add_option(self):

		self.parser.add_option("--check", action="store_true", help="Do not upload anything; just check file mappings")		
		self.parser.add_option("--ctf", action="store_true", help="Upload CTF Parameters")
		self.parser.add_option("--boxes", action="store_true", help="Upload Boxes")
		self.parser.add_option("--eman1", action="store_true", help="Look in EMAN1-style files instead of EMAN2 project database")		
		self.parser.add_option("--ptcls", action="store_true", help="Upload Particle Sets")
		self.parser.add_option("--confirm", action="store_true", help="Request confirmation of mappings before proceeding")		
		self.parser.add_option("--clip_filename", type="string", help="Remove a substr from source filenames to aid matching, e.g. _filt_ptcls")		
		self.parser.add_option("--match", type="string", help="Restrict particle sets to this substring, e.g. for quick testing")
		#self.parser.add_option("--check_boxsize", type="int", help="Check if a micrograph has been shrunk; if box_size < check_boxsize, zoom by box_size / check_boxsize")
		self.parser.add_option("--shrink_factor", type="float", help="Specify a shrink factor (e.g. 0.25 if the boxed micrograph was reduced by a factor of 4)", default=1.0)


		usage = """%prog sync [options] <project record ID>
		
This program will upload CTF parameters, boxes, and particle sets into an EMEN2 database. Because filenames are not always globally unique, you must specify a base project that will be searched for files.

If run with "--check", it will only test to see if it can make the correct file mappings between the local EMAN2 project and the remote database. Nothing will be written. This is a good way to test to see if the correct mappings can be made before you attempt to commit any changes.
		
		"""

		self.parser.set_usage(usage)
		

	def check_args(self):
		try:
			self.project = int(self.args[0])
		except:
			raise Exception, "Project record ID required"			
			
		if EMAN2 == None:
			raise Exception, "EMAN2 support is required"	
			
		if self.options.clip_filename:
			print "Debug: Clipping '%s' from filenames"%self.options.clip_filename
			test = "filename_ok_filt.mrc"
			test_clip = test.replace(self.options.clip_filename or '', "")
			print "\t%s -> %s"%(test, test_clip)
			
			
			
	def run(self):
		
		self.sources = set()
		self.source_bdo = {}
		self.source_ctf = {}
		self.source_boxes = {}
		self.source_box_size = {}
		self.source_quality = {}
		
		self.tmpdirs = []
		self.projrecs = []
				
		# Check arguments and auth
		self.check_args()	
		self.login()

		# Actually read EMAN1 or EMAN2 CTF/boxes
		if self.options.eman1:
			self.readEMAN1()
		else:
			self.readEMAN2()

		# Check project and get remote records
		self.getremoteinfo()

		# Make mappings
		self.findbinary()

		# Upload 
		if self.options.ctf:
			self.uploadctf()

		if self.options.boxes:
			self.uploadboxes()

		if self.options.ptcls:
			self.uploadptcls()

		self.cleanup()


	def cleanup(self):
		if not self.tmpdirs:
			return
			
		print "You will want to remove the following tmp directories. They are not automatically removed to prevent accidental file loss in the event of a bug or other error."
		for k in self.tmpdirs:
			print k



	def getremoteinfo(self):
		# Since filenames are not guaranteed unique, restrict to a project...
		print "\nChecking project %s on %s"%(self.project, self.db._host)
		
		self.projrecs = self.db.getchildren(self.project, -1, ["ccd","scan"])
		print "\tFound %s ccds/scans in remote project"%len(self.projrecs)

		self.bdos = self.db.getbinary(self.projrecs)
		print "\tFound %s binaries"%len(self.bdos)
		

		
	def readEMAN1(self):
		if self.options.ctf:
			self.readEMAN1ctf()
		if self.options.boxes:
			self.readEMAN1boxes()


	def readEMAN1ctf(self):
		try:
			f = open("ctfparm.txt", "r")
			r = [i.strip() for i in f.readlines()]
			f.close()
		except:
			print "No ctfparm.txt present; could not load EMAN1 CTF parameters"
			return
		
		indexes = {
			"ctf_defocus_measured": 0, # multiply by -1
			"ctf_bfactor": 3,
			"tem_voltage": 10,
			"aberration_spherical": 11,
			"angstroms_per_pixel": 12,
			"ctf_astig_angle": 2 ,
			"ctf_astig_defocus_diff": 1,
			"ctf_ampcont": 5, # multiply by 100
		}

		for i in r:
			try:
				source, params = i.split("\t")
				params = params.split(",")
			
				# it seems that if there are 14 params, astig params are inserted as 1,2
				shift = 0
				if len(params) == 14:
					shift = 2
			
				ctf = EMAN2.EMAN2Ctf()
				ctf.defocus = float(params[0]) * -1
				ctf.bfactor = float(params[1+shift])
				ctf.ampcont = float(params[3+shift]) * 100
			
				# print "defocus %s, bfactor %s, ampcont %s"%(ctf.defocus, ctf.bfactor, ctf.ampcont)
				self.source_ctf[source] = ctf
			except:
				print "Unable to parse CTF parameters, skipping"
				

	def readEMAN1boxes(self):
		boxes = glob.glob("*.box")
		
		if not boxes:
			print "No EMAN1 .box files found"
			return
		
		print "Found EMAN1 boxes: %s"%boxes
		
		for box in boxes:
			try:
				source = os.path.splitext(box)[0]
				box_size, coords = self.readEMAN1box(box)
			except:
				print "Could not load data from box %s, skipping"%box
				continue

			self.source_box_size[source] = box_size
			self.source_boxes[source] = coords
				

	def readEMAN1box(self, box):
		f = open(box, "r")
		r = [i.split() for i in f.readlines()]
		f.close()
	
		coords = []
		for b in r:
			box_size = int(b[2])
			xc = (int(b[0]) + (box_size/2)) / self.options.shrink_factor
			yc = (int(b[1]) + (box_size/2)) / self.options.shrink_factor
			coords.append([int(xc),int(yc)])

		box_size = int(box_size / self.options.shrink_factor)

		return box_size, coords
			
			
	def readEMAN2(self):

		###############
		# Find all the raw images, particle sets, and CTF parameters in the local DB
		###############

		print "\nOpening EMAN2 local project"

		projdb = EMAN2.db_open_dict("bdb:project")
		ptclsets = projdb.get("global.spr_ptcls_dict", {})
		e2ctfparms = EMAN2.db_open_dict("bdb:e2ctf.parms")
		total = len(ptclsets)
		
		# Read the EMAN2 managed particle sets
		for count, (k,v) in enumerate(ptclsets.items()):
			ref = v.get('Phase flipped') or v.get('Original Data')
	
			if not ref:
				print "No particles found for %s, skipping"%k
				continue
		
			print "%s of %s: Getting info for particle set %s from %s"%(count+1, total, k, ref)
		
			d = EMAN2.db_open_dict(ref)

			coords = []
			maxrec = d['maxrec']
			if maxrec == None:
				print "No particles in %s, skipping"%(k)
				d.close()				
				continue				
				
				
			# Get info from first particle in stack	
			ptcl = d[0]
			
			try: ctf = ptcl.get_attr("ctf")
			except: ctf = None
			
			try: box_size = int(ptcl.get_attr("nx") / self.options.shrink_factor)
			except: box_size = None

			# Mangle source name if necessary
			try:
				source = os.path.basename(ptcl.get_attr('ptcl_source_image'))
			except:
				source = k
				source = source.split("_ptcl")[0]
				print "\tUnable to get source image %s, using %s for the filename search"%(k,k)
				

			# Try to read boxes from particle headers
			if self.options.boxes:
				for i in range(maxrec+1):
					dptcl = d[i]
					try:
						x, y = dptcl.get_attr('ptcl_source_coord')
						x /= self.options.shrink_factor
						y /= self.options.shrink_factor
						coords.append([int(x), int(y)])
					except:
						coords.append(None)


				if None in coords:
					print "\tSome particles for %s did not specify coordinates"%k
					coords = []
					# self.options.boxes = False
				
				print "Got box_size %s and coords %s"%(box_size, coords)
					
					

			# Get alternative CTF and quality from e2ctfit settings
			ctf2, quality = self.readEMAN2e2ctfparms(e2ctfparms, k)
			if not ctf and ctf2:
				print "\tUsing CTF parameters from bdb:e2ctf.parms#%s"%k
				ctf = ctf2	
			
			if ctf:		
				self.source_ctf[source] = ctf
				print "Got CTF: defocus %s, B-factor %s"%(ctf.defocus, ctf.bfactor)
			else:
				print "\tNo CTF for %s"%k				
			
			if box_size and coords:
				self.source_box_size[source] = box_size
				self.source_boxes[source] = coords

			if quality:
				self.source_quality[source] = quality
				
			d.close()	
	
	
		# If we can't find any EMAN2 managed particle sets, at least check e2ctf.parms for any CTF params	
		if not ptclsets:
			print "No EMAN2 managed particle sets found; checking e2ctf.parms for CTF instead"
			# self.options.boxes = False
			
			for k,v in e2ctfparms.items():
				ctf, quality = self.readEMAN2e2ctfparms(e2ctfparms, k)
				source = k.split("_ptcl")[0]
				
				if ctf:
					self.source_ctf[source] = ctf
				if quality:
					self.source_quality[source] = quality

		print "\n%s Files in local project: "%len(self.source_ctf), self.source_ctf.keys()

		projdb.close()
		e2ctfparms.close()
		
		

	def readEMAN2e2ctfparms(self, e2ctfparms, item):
		ctf2 = None
		quality = None
		
		ctf_str = e2ctfparms.get(item)
		if ctf_str:
			try:
				ctf2 = EMAN2.EMAN2Ctf()
				ctf2.from_string(ctf_str[0])
				quality = ctf_str[-1]
			except:
				print "\tProblem reading CTF from bdb:e2ctf.parms for %s"%item
	
		return ctf2, quality


	def __filematch(self, i, j):
		if i in j: return True


	def findbinary(self):
		
		###############
		# For each file in the local project, search for a matching BDO and record ID in the database
		###############

		self.sources = set(self.source_ctf.keys() + self.source_boxes.keys())
		total = len(self.sources)

		bdosbyfilename = dict([[i["filename"], i] for i in self.bdos])
		filenames = bdosbyfilename.keys()

		ambiguous = {}
		nomatches = {}

		# ian: new style: I fixed getbinary, so just do one big BDO lookup instead of many db.findbinary
		
		print "remote filenames?:"
		print filenames

		for count, source in enumerate(self.sources):
		 	print "\n%s of %s: Looking up file %s in project %s"%(count+1, total, source, self.project)

			q = source.replace(self.options.clip_filename or '', "")
			gzipstrip = ".gz" in source

			# ian: replace this with functools.partial
			# matches = map(bdosbyfilename.get, filter(lambda x:q in x.split("."), filenames))
			matches = []
			for i in filenames:
				i2 = i
				if gzipstrip:
					i2 = i.replace('.gz','')
				#if q in i2:
				if self.__filematch(q, i2):
					matches.append(bdosbyfilename.get(i))
			
			
			if len(matches) > 1:
				print "\tAmbiguous match for %s: %s"%(source, [match["filename"] for match in matches])
				ambiguous[source] = matches
				continue	
		
		
			if len(matches) == 0:
				print "\tNo matches for %s"%source
				nomatches[source] = []
				continue
		
		
			match = matches[0]
			print "\tFound %s in record %s, matching filename is %s"%(match["name"], match["recid"], match["filename"])
			
			self.source_bdo[source] = match
		
						
		print "\n\nSuccessful Local-Remote Mappings:"
		for k,v in self.source_bdo.items():
			print "\t%s -> %s (recid %s)"%(k, v["filename"], v["recid"])
			
			
		if ambiguous:
			print "\n\nAmbiguous matches:"
			for k,v in ambiguous.items():
				print "\t%s -> %s"%(k, ",".join(v))
		
		if nomatches:	
			print "\n\nNo matches:"
			for k in nomatches:
				print "\t%s"%k
			
			
		if self.options.confirm:
			answer = raw_input("\nContinue with these mappings Y/N? ")
			if not answer.lower() in ["y","yes"]:
				print "Mappings rejected; SyncController exiting"
				# ian: do this in a better way, but not a harsh sys.exit
				self.options.ctf = False
				self.options.boxes = False


	def uploadctf(self):

		###############
		# Now that we have found all source images, CTF parameters, and database references, prepare to save in DB.
		###############

		recs = [i["recid"] for i in self.source_bdo.values()]
		recs = self.db.getrecord(recs)
		recs = dict([(i["recid"],i) for i in recs])
		putrecs = []
		
		print "\nGot %s records for CTF parameter upload"%len(recs)

		for source, match in self.source_bdo.items():
	
			ctf = self.source_ctf.get(source)
			quality = self.source_quality.get(source)

			if not ctf:
				print "No CTF found for %s"%source
				continue
				
			rec = recs[match["recid"]]

			if rec.get("ctf_defocus_measured") == ctf.defocus and rec.get("ctf_bfactor") == ctf.bfactor and rec.get("ctf_ampcont") == ctf.ampcont and rec.get("assess_image_quality") == quality:
				print "%s already has current CTF parameters"%source
				continue

			rec["ctf_defocus_measured"] = ctf.defocus
			rec["ctf_bfactor"] = ctf.bfactor
			rec["ctf_ampcont"] = ctf.ampcont

			if quality != None:
				rec["assess_image_quality"] = quality

			putrecs.append(rec)
	
		###############
		# Store / upload
		###############

		print "\nCommitting %s updated records with changed CTF..."%(len(putrecs))

		self.db.putrecord(putrecs)



	def uploadboxes(self):
		
		print "\nPreparing to upload boxes..."
		newboxes = []
		
		for source, boxes in self.source_boxes.items():
			
			bdo = self.source_bdo.get(source)
			box_size = self.source_box_size[source]
			
			if not bdo:
				print "\tNo bdo for source %s..."%source
				continue
			
			recid = bdo.get("recid")	
			
			# Check remote site for existing boxes
			remoteboxes = self.db.getchildren(recid, -1, ["box"])

			if len(remoteboxes) == 1:
				print "\tUpdating existing box record"
				newbox = self.db.getrecord(remoteboxes.pop())				

			else:
				if len(remoteboxes) > 1:
					print "\tNote: more than one box record already specified!"
				
				print "\tCreating new box record"
				newbox = self.db.newrecord("box", recid)
			
			print "\t%s / %s has %s boxes with box size %s"%(source, recid, len(boxes), box_size)

			newbox["box_coords"] = boxes
			newbox["box_size"] = box_size
			
			newboxes.append(newbox)
			
			# print "-----"
			# print newbox
			# print "-----"

		print "\tCommitting %s box records"%(len(newboxes))
		newrecs = self.db.putrecord(newboxes)
		# print "\t... %s"%", ".join([i.get('recid') for i in newrecs])
			
		
		
	def uploadptcls(self):
		
		print "\nUploading particles..."
		
		print "Not implemented yet"
		
		for source in self.source_ctf.keys():
			pass
		





##################################
# Mappings
##################################

map_filename_rectype = {
	".st":"stack",
	".dm3":"ccd"
}

map_rectype_handler = {
	"ccd": CCDHandler,
	"stack": StackHandler,
	"scan": ScanHandler,
	"volume": VolumeHandler,
	"none": FileHandler
}


def gethandler(filename=None, rectype=None):
	if not rectype:
		ext = os.path.splitext(filename)[-1]
		rectype = map_filename_rectype.get(ext)	
	return map_rectype_handler.get(rectype, FileHandler)




##################################
# Exported utility methods
##################################

def quickdb(host='http://ncmidb.bcm.edu', username=None, password=None):
	username = username or raw_input("Username: ")
	password = password or getpass.getpass("Password: ")
	db = DBXMLRPCProxy(host=host)
	db._login(username=username, password=password)
	return db


def sign(a):
	if a>0: return 1
	return -1




##################################
# Main()
##################################


def print_help():
	print """%s <action>
	
Actions available: upload, download

For detailed help: %s <action> --help
	
	"""%(sys.argv[0],sys.argv[0])




def main():
	parser = DBRPCOptions(version=VERSION)

	controllers = {
		"download":DownloadController,
		"upload":UploadController,
		"sync":SyncController
	}

	try:
		action = sys.argv[1]
	except:
		action = "help"

	if action in ["-h","--help","help"]:
		return print_help()

	try:
		if len(sys.argv) == 2:
			sys.argv.append("-h")
		controller = controllers.get(action)(args=sys.argv[2:], version=VERSION)
		
	except Exception, inst:
		print "Error: %s"%inst
		sys.exit(0)

	controller.run()



if __name__ == "__main__":
	main()



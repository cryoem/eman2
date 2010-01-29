#!/usr/bin/env python

#Author: Ross Coleman

from EMAN2 import test_image, get_image_directory, Transform, Region, \
    EMANVERSION, EMData, E2init, E2end
from EMAN2db import db_open_dict, db_check_dict, db_close_dict
from PyQt4 import QtGui, QtCore
from emapplication import EMStandAloneApplication, get_application
from emimage2d import EMImage2DModule
from emselector import EMSelectorModule
from emshape import EMShape, EMShapeDict
from math import *
from optparse import OptionParser
import weakref
import sys
import os
#from emimagemx import EMImageMXModule


"""
This program is used to box segments from a frame and extract particles from the segments.
A frame is the the EM micrograph image that contains the biological macromolecules that will be boxed.
A segment is a rectangular region that shows a helical region of a biological macromolecule. 
Different segments from the same frame will generally have different lengths but the same width.
A particle is a square or rectangle from inside a segment. Particles are chosen so that they overlap each other
within a segment. Usually, all particles from a frame will have the same dimensions.
"""


E2HELIXBOXER_DB = "bdb:e2helixboxercache"
def counterGen():
    """
    Calling this function will create a counter object.
    Ex: "counter = counterGen()"
    Then calling "counter.next()" will return the next number in {1, 2, 3, ...}
    """
    i = 0
    while True:
        i += 1
        yield i

def get_particles_from_segment( segment, px_overlap, px_length = None, px_width = None ):
    """
    Gets the square/rectangular "particles" inside a rectangular segment.
    @segment: EMData object that contains the rectangular region boxed by the user
    @px_overlap: pixels of overlap of the rectangular particles taken from the segment
    @px_length: length of the particles in pixels, defaults to the width of the segment
    @px_width: width of the particles in pixels, defaults to px_length
    @return: a list of EMData particles
    """
    if not px_width:
        px_width = segment.get_xsize()
    if not px_length:
        px_length = px_width
    assert px_length > px_overlap, "The overlap must be smaller than the particle length"
    
    px_step = px_length - px_overlap
    
    (xsize, ysize) = (px_width, px_length)
    (x0,y0) = (0,0)
    particles = []
    while y0 + ysize <= segment.get_ysize():
        particles.append( segment.get_clip(Region(x0, y0, xsize, ysize)) )
        y0 += px_step
    
    return particles

def get_segment_from_coords(frame, x1, y1, x2, y2, width):
    """
    Gets the rectangular segment of the image specified by the coordinates.
    @frame: the EMData object which holds the micrograph from which segments and particles are chosen
    @x1: x-coordinate in pixels of the first end-point along the long axis of symmetry of the rectangle
    @y1: y-coordinate in pixels of the first end-point along the long axis of symmetry of the rectangle
    @x2: x-coordinate in pixels of the second end-point along the long axis of symmetry of the rectangle
    @y2: y-coordinate in pixels of the second end-point along the long axis of symmetry of the rectangle
    @width: the width in pixels of the rectangle
    @return: the rectangular EMData segment specified by the coordinates and width 
    """
    l_vect = (x2-x1, y2-y1)
    length = sqrt(l_vect[0]**2+l_vect[1]**2)
    assert length != 0
    l_uvect = (l_vect[0]/length, l_vect[1]/length)
    centroid = ( (x1+x2)/2.0, (y1+y2)/2.0 )
    
    #Rotate so that the length is parallel to the y-axis
    #Angle between l_uvect and y-axis: l_uvect (dot) j_hat = cos (rot_angle)
    rot_angle = 180/pi*acos( l_uvect[1] )
    #Whether we rotate clockwise or counterclockwise depends on the sign of l_uvect[0] (the x-component)
    if l_uvect[0] < 0:
        rot_angle *= -1 #rotate the box clockwise
    
    tr = Transform()
    tr.set_trans(centroid)
    tr.set_rotation({"type":"2d", "alpha":rot_angle})
    segment_dimensions = ( int(round(width)), int(round(length)), 1 )
    segment = frame.get_rotated_clip( tr, segment_dimensions )
    return segment

def load_coords(coords_filepath):
    """
    load coordinates from a *.box or *.hbox file
    @return a list of tuples [(x0, x1, y1, y2, width), ...] 
    """
    (path, filename) = os.path.split(coords_filepath)
    (basename, extension) = os.path.splitext(filename)
    if extension == ".box":
        data = []
        datum = [None]*5
        for line in open(coords_filepath):
            line = line.split("\t")
            for i in range(len(line)):
                line[i] = int(line[i])
            if line[4] == -1:
                w = line[2]
                r = w / 2.0
                datum[0] = line[0] + r
                datum[1] = line[1] + r
                datum[4] = w
            elif line[4] == -2:
                assert line[2] == w
                r = w / 2.0
                datum[2] = line[0] + r
                datum[3] = line[1] + r
                w = None
                r = None
                data.append(tuple(datum))
                datum = [None]*5
    elif extension == ".hbox":
        data = []
        for line in open(coords_filepath):
            line = line.split("\t")
            for i in range(len(line)):
                line[i] = int(line[i])
            data.append(tuple(line))
    else:
        pass
    
    return data
def save_coords(coords_list, frame_name, output_type = "box", output_dir=""):
    """
    Saves coordinates and widths of the boxed segments to a file.
    @coords_list: a list of tuples (x1, y1, x2, y2, width), with each tuple corresponding to a segment
    @frame_name: the file name of the frame without the file extension
    @output_type: the file format used for saving coordinates. Can be "box" or "hbox".
    "*.box": the EMAN1 file format -- r is half the width (w) of the boxes
        x1-r    y1-r    w    w    -1
        x2-r    y2-r    w    w    -2
    "*.hbox": the EMAN2 file format
        x1    y1    x2    y2    w
    @output_dir: the file will be saved in {output_dir}/box.xyz/ where xyz is the width of the boxes     
    """
    if output_type == "box":
        output_filepath = os.path.join(output_dir, "box." + str(coords_list[0][4]))
        if not os.access(output_filepath, os.F_OK):
            os.mkdir(output_filepath)
        output_filepath = os.path.join(output_filepath, frame_name + ".box")
        out_file = open(output_filepath, "w")
        
        for coords in coords_list:
            (x1, y1) = (coords[0], coords[1])
            (x2, y2) = (coords[2], coords[3])
            width = coords[4]
            r = width / 2.0
                        
            #For some reason, EMAN1 subtracts half the box width from each coordinate
            #EMAN1 uses <cstdio> fprintf() and "%1.0f", which rounds half-integers away from zero
            #the string format operator works the same in Python as it does in C for decimal floats
            out_file.write( "%1.0f\t%1.0f\t%1.0f\t%1.0f\t-1\n" % (x1 - r, y1 - r, width, width) )
            out_file.write( "%1.0f\t%1.0f\t%1.0f\t%1.0f\t-2\n" % (x2 - r, y2 - r, width, width) )
        out_file.close()
        
    elif output_type == "hbox":
        output_filepath = os.path.join(output_dir, "box." + str(coords_list[0][4]))
        if not os.access(output_filepath, os.F_OK):
            os.mkdir(output_filepath)
        output_filepath = os.path.join(output_filepath, frame_name + ".hbox")
        out_file = open(output_filepath, "w")
        for coords in coords_list:
            out_file.write( "%i\t%i\t%i\t%i\t%i\n" % (coords[0], coords[1], coords[2], coords[3], coords[4]) )
        out_file.close()
        
    else:
        pass

def db_save_coords(frame_filepath, output_type = "box", output_dir=""):
    frame_filename = os.path.basename(frame_filepath)
    frame_name = os.path.splitext( frame_filename )[0]
    db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
    box_coords_list = db[frame_filename]
    save_coords(box_coords_list, frame_name, output_type, output_dir)
    
def db_get_segments_dict(frame_filepath):
    frame = EMData(frame_filepath)
    frame_filename = os.path.basename(frame_filepath)
    frame_name = os.path.splitext( frame_filename )[0]
    db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
    box_coords_list = db[frame_filename]
    if not box_coords_list:
        return {}
    segments_dict = {}
    for coords in box_coords_list:
        segment = get_segment_from_coords(frame, *coords)
        segments_dict[tuple(coords)] = segment
    return segments_dict

def save_segment(segment_emdata, frame_name, segment_num, output_type = "hdf", output_dir=""):
    """
    Saves a boxed segment to an image file.
    @segment_emdata: the EMData object that holds the image data for the segment
    @frame_name: the file name of the frame without the file extension
    @segment_num: the number that identifies this segment among those boxed from this frame
    @output_type: the image file type to write out
    @output_dir: the file will be saved in {output_dir}/seg.xyz/ where xyz is the width of the segment 
    """
    output_fname = "%s.%i.seg.%s" % (frame_name, segment_num, output_type)
    output_filepath = os.path.join(output_dir, "seg." + str(segment_emdata.get_xsize()))
    if not os.access(output_filepath, os.F_OK):
        os.mkdir(output_filepath)
    output_filepath = os.path.join(output_filepath, output_fname)
    if os.access(output_filepath, os.F_OK):
        os.remove(output_filepath)
    segment_emdata.write_image( output_filepath )

def db_save_segments(frame_filepath, output_type = "hdf", output_dir=""):
    frame_filename = os.path.basename(frame_filepath)
    frame_name = os.path.splitext( frame_filename )[0]
    segments_dict = db_get_segments(frame_filepath)
    i = 1
    for coords in segments_dict:
        segment = segments_dict[coords]
        save_segment(segment, frame_name, i, output_type, output_dir)
        i+=1

def save_particles(segment_emdata, segment_num, frame_name, px_overlap, px_length = None, px_width = None, output_type = "hdf", output_dir=""):
    """
    saves the particles in a segment to a stack file
    @segment_emdata: the EMData object containing the image data for the segment
    @segment_num: the number used to identify this segment among all those boxed in the frame
    @frame_name: the file name of the frame without the file extension
    @px_overlap, @px_length, @px_width: see get_particles_from_segment() function
    @output_type: the output image file type
    @output_dir: the output file will be saved in {output_dir}/ptcl.xyz where xyz is the width of each of the particles
    """
    particles = get_particles_from_segment(segment_emdata, px_overlap, px_length, px_width)
    output_filename = "%s.%i.ptcl.%s" % (frame_name, segment_num, output_type)
    output_filepath = os.path.join(output_dir, "ptcl." + str(px_width))
    if not os.access(output_filepath, os.F_OK):
        os.mkdir(output_filepath)
    output_filepath = os.path.join(output_filepath, output_filename)
    if os.access(output_filepath, os.F_OK):
        os.remove(output_filepath)
    for i in range(len(particles)):
        ptcl = particles[i]
        ptcl.write_image(output_filepath, i) #appending to the image stack


def db_save_particles(frame_filepath, px_overlap, px_length = None, px_width = None, output_type = "hdf", output_dir=""):
    pass

class EMHelixBoxerWidget(QtGui.QWidget):
    def __init__(self, frame_filepath, app):
        QtGui.QWidget.__init__(self)
        self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
        self.setWindowTitle("e2helixboxer")
        
        if not frame_filepath:
            self.frame_filepath = "test_image"
            img = test_image()
        else:
            self.frame_filepath = os.path.relpath(frame_filepath)
            img = EMData(frame_filepath)
        
        self.main_image = EMImage2DModule(img, application=app)
        self.main_image.set_file_name( os.path.split(self.frame_filepath)[1] ) # TODO: determine if this should use the entire file path
        self.main_image.shapes = EMShapeDict()
        self.segment_viewer = None #Will be an EMImage2DModule instance
        self.edit_mode = None #Values are in {None, "new", "move", "2nd_point", "1st_point", "delete"}
        self.current_boxkey = None
        self.initial_helix_box_data_tuple = None
        self.click_loc = None #Will be (x,y) tuple
        self.segments_dict = db_get_segments_dict(frame_filepath) #Will be like {(x1,y1,x2,y2,width): emdata}
        self.color = (0, 0, 1)
        self.selected_color = (0, 1, 0)
        self.counter = counterGen()
        self.coords_file_extension_dict = {"EMAN1":"box", "EMAN2": "hbox"}
        self.image_file_extension_dict = {"MRC":"mrc", "Spider":"spi", "Imagic": "img", "HDF5": "hdf"}

        if self.get_db_item("boxes") == None:
            self.set_db_item("boxes", [])
        else:
            boxList = self.get_db_item("boxes")
            for box_coords in boxList:
                key = self.generate_emshape_key()
                emshape_list = ["rectline"]
                emshape_list.extend(list(self.color))
                emshape_list.extend(list(box_coords))
                emshape_list.append(2)
                emshape = EMShape( emshape_list )
                self.main_image.add_shape(key, emshape)
            self.main_image.updateGL()

        self.__create_ui()
        
        qual = self.get_image_quality()
        if qual:
            self.img_quality_combobox.setCurrentIndex( qual )
        else:
            self.img_quality_combobox.setCurrentIndex( 2 )
        
        self.main_image.optimally_resize()
        
        self.coords_ftype_combobox.addItems( sorted(self.coords_file_extension_dict.keys()) )
        self.segs_ftype_combobox.addItems( sorted(self.image_file_extension_dict.keys()) )
        self.ptcls_ftype_combobox.addItems( sorted(self.image_file_extension_dict.keys()) )
        width = 100
        if self.segments_dict:
            first_coords = self.segments_dict.keys()[0]
            width = first_coords[4]
        self.box_width_spinbox.setValue(width)
        self.ptcls_width_spinbox.setValue( width )
        self.ptcls_length_spinbox.setValue( width )
        self.ptcls_overlap_spinbox.setValue( int(0.9*width) )
        
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedown"), self.mouse_down)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedrag"), self.mouse_drag)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mouseup"), self.mouse_up)
        self.connect(self.load_boxes_button, QtCore.SIGNAL("clicked()"), self.load_boxes )
        self.connect(self.output_dir_pushbutton, QtCore.SIGNAL("clicked()"), self.choose_dir )
        self.connect(self.write_output_button, QtCore.SIGNAL("clicked()"), self.write_ouput )
        self.connect(self.box_width_spinbox, QtCore.SIGNAL("valueChanged(int)"), self.width_changed)
        self.connect(self.done_but, QtCore.SIGNAL("clicked()"), self.exit_app )
        self.connect( self.img_quality_combobox, QtCore.SIGNAL("currentIndexChanged(int)"), self.set_image_quality )
        get_application().show_specific(self.main_image)        

        
    def __create_ui(self):
        self.box_width_label = QtGui.QLabel(self.tr("Box &Width:"))
        self.box_width_spinbox = QtGui.QSpinBox()
        self.box_width_spinbox.setMaximum(1000)
        self.box_width_label.setBuddy(self.box_width_spinbox)
        
        self.img_quality_label = QtGui.QLabel(self.tr("Image &Quality:"))
        self.img_quality_combobox = QtGui.QComboBox()
        qualities = [str(i) for i in range(5)]
        self.img_quality_combobox.addItems(qualities)
        self.img_quality_combobox.setCurrentIndex(2)
        self.img_quality_label.setBuddy(self.img_quality_combobox)
        #self.img_quality_combobox.setEnabled(False)

        self.load_boxes_button = QtGui.QPushButton(self.tr("&Load Boxes"))
        
        self.coords_groupbox = QtGui.QGroupBox(self.tr("Write &Coordinates:"))
        self.coords_groupbox.setCheckable(True)
        coords_ftype_label = QtGui.QLabel(self.tr("&File Type:"))
        self.coords_ftype_combobox = QtGui.QComboBox()
        coords_ftype_label.setBuddy(self.coords_ftype_combobox)
        
        self.segs_groupbox = QtGui.QGroupBox(self.tr("Write &Segments:"))
        self.segs_groupbox.setCheckable(True)
        segs_ftype_label = QtGui.QLabel(self.tr("File &Type:"))
        self.segs_ftype_combobox = QtGui.QComboBox()
        segs_ftype_label.setBuddy(self.segs_ftype_combobox)
        
        self.ptcls_groupbox = QtGui.QGroupBox(self.tr("Write &Particles:"))
        self.ptcls_groupbox.setCheckable(True)
        ptcls_ftype_label = QtGui.QLabel(self.tr("File T&ype:"))
        self.ptcls_ftype_combobox = QtGui.QComboBox()
        ptcls_ftype_label.setBuddy(self.ptcls_ftype_combobox)
        ptcls_overlap_label = QtGui.QLabel(self.tr("&Overlap:"))
        self.ptcls_overlap_spinbox = QtGui.QSpinBox()
        self.ptcls_overlap_spinbox.setMaximum(1000)
        ptcls_overlap_label.setBuddy(self.ptcls_overlap_spinbox)
        ptcls_width_label = QtGui.QLabel(self.tr("W&idth:"))
        self.ptcls_width_spinbox = QtGui.QSpinBox()
        self.ptcls_width_spinbox.setMaximum(1000)
        ptcls_width_label.setBuddy(self.ptcls_width_spinbox)
        ptcls_length_label = QtGui.QLabel(self.tr("&Length:"))
        self.ptcls_length_spinbox = QtGui.QSpinBox()
        self.ptcls_length_spinbox.setMaximum(1000)
        ptcls_length_label.setBuddy(self.ptcls_length_spinbox)
        
        output_dir_label = QtGui.QLabel(self.tr("Output &Directory"))
        self.output_dir_line_edit = QtGui.QLineEdit()
        output_dir_label.setBuddy(self.output_dir_line_edit)
        self.output_dir_pushbutton = QtGui.QPushButton(self.tr("&Browse"))
        self.write_output_button = QtGui.QPushButton(self.tr("W&rite Output"))

        self.done_but=QtGui.QPushButton(self.tr("&Done"))
        self.status_bar = QtGui.QStatusBar()
        self.status_bar.showMessage("Ready",10000)
        
        
        
        widthLayout = QtGui.QHBoxLayout()
        widthLayout.addWidget(self.box_width_label)
        widthLayout.addWidget(self.box_width_spinbox)
        
        qualityLayout = QtGui.QHBoxLayout()
        qualityLayout.addWidget(self.img_quality_label)
        qualityLayout.addWidget(self.img_quality_combobox)
        
        coords_ftype_layout = QtGui.QHBoxLayout()
        coords_ftype_layout.addWidget(coords_ftype_label)
        coords_ftype_layout.addWidget(self.coords_ftype_combobox)
        self.coords_groupbox.setLayout(coords_ftype_layout)
        
        segs_ftype_layout = QtGui.QHBoxLayout()
        segs_ftype_layout.addWidget(segs_ftype_label)
        segs_ftype_layout.addWidget(self.segs_ftype_combobox)
        self.segs_groupbox.setLayout(segs_ftype_layout)
        
        ptcls_ftype_layout = QtGui.QHBoxLayout()
        ptcls_ftype_layout.addWidget(ptcls_ftype_label)
        ptcls_ftype_layout.addWidget(self.ptcls_ftype_combobox)
        
        ptcls_overlap_layout = QtGui.QHBoxLayout()
        ptcls_overlap_layout.addWidget(ptcls_overlap_label)
        ptcls_overlap_layout.addWidget(self.ptcls_overlap_spinbox)
        
        ptcls_width_layout = QtGui.QHBoxLayout()
        ptcls_width_layout.addWidget(ptcls_width_label)
        ptcls_width_layout.addWidget(self.ptcls_width_spinbox)
        
        ptcls_length_layout = QtGui.QHBoxLayout()
        ptcls_length_layout.addWidget(ptcls_length_label)
        ptcls_length_layout.addWidget(self.ptcls_length_spinbox)
        
        ptcls_opts_layout = QtGui.QVBoxLayout()
        ptcls_opts_layout.addLayout(ptcls_ftype_layout)
        ptcls_opts_layout.addLayout(ptcls_overlap_layout)
        ptcls_opts_layout.addLayout(ptcls_width_layout)
        ptcls_opts_layout.addLayout(ptcls_length_layout)
        self.ptcls_groupbox.setLayout(ptcls_opts_layout)
        
        directory_layout = QtGui.QHBoxLayout()
        directory_layout.addWidget(output_dir_label)
        directory_layout.addWidget(self.output_dir_line_edit)
        
        button_layout = QtGui.QHBoxLayout()
        button_layout.addWidget(self.output_dir_pushbutton)
        button_layout.addWidget(self.write_output_button)
        
        layout = QtGui.QVBoxLayout()

        
        self.vbl = QtGui.QVBoxLayout(self)
        self.vbl.setMargin(0)
        self.vbl.setSpacing(6)
        self.vbl.setObjectName("vbl")
        self.vbl.addLayout(widthLayout)
        self.vbl.addLayout(qualityLayout)
        self.vbl.addWidget(self.load_boxes_button)
        self.vbl.addWidget(self.coords_groupbox)
        self.vbl.addWidget(self.segs_groupbox)
        self.vbl.addWidget(self.ptcls_groupbox)
        self.vbl.addLayout(directory_layout)
        self.vbl.addLayout(button_layout)
        self.vbl.addWidget(self.done_but)
        self.vbl.addWidget(self.status_bar)

    def choose_dir(self):
        """
        launches a file-browser dialog to select the base directory in which to save output coordinates, segments, and particles 
        """
        #selector = EMSelectorModule(save_as_mode=False)
        #selector.widget.save_as_line_edit.setEnabled(False)
        #path = selector.exec_()
        path = QtGui.QFileDialog.getExistingDirectory(self)
        self.output_dir_line_edit.setText(path)
    def color_boxes(self):
        emshapes_dict = self.main_image.get_shapes()
        for key in emshapes_dict:
            shape = emshapes_dict.get(key).shape
            for i in range(3):
                shape[i+1] = self.color[i]
        current_shape = emshapes_dict.get(self.current_boxkey)
        if current_shape:
            for i in range(3):
                current_shape.shape[i+1] = self.selected_color[i]
        self.main_image.shapechange=1
        self.main_image.updateGL()
    def display_segment(self, segment_emdata):
        """
        launches or updates an EMImage2DModule to display segment_emdata
        @segment_emdata: an EMData object that stores the image data for a segment
        """
        self.color_boxes()
        if not self.segment_viewer:
            self.segment_viewer = EMImage2DModule(application=get_application())
            self.segment_viewer.desktop_hint = "rotor" # this is to make it work in the desktop
            self.segment_viewer.setWindowTitle("Current Boxed Segment")
            self.segment_viewer.get_qt_widget().resize(200,800)
        self.segment_viewer.set_data(segment_emdata)

        get_application().show_specific(self.segment_viewer)
        scale = 100.0 / self.get_width()
        self.segment_viewer.set_scale(scale)
        if self.segment_viewer.inspector:
            self.segment_viewer.inspector.set_scale(scale)
        self.segment_viewer.updateGL()
    def exit_app(self):
        """
        quits the program
        """
        app = get_application()
        app.quit()
    def generate_emshape_key(self):
        """
        creates a unique key for a new "rectline" EMShape, which is used for boxing a segment
        @return: a string that is the key for the new "rectline" EMShape
        """
        i = self.counter.next()
        return "rectline%i" % i
    def get_width(self):
        """
        returns the current width for the segments
        """
        return self.box_width_spinbox.value()
    def load_boxes(self):
        """
        load boxes from a file selected in a file browser dialog
        """
        path = QtGui.QFileDialog.getOpenFileName(self, self.tr("Open Box Coordinates File"), "", self.tr("Boxes (*.box *.hbox)"))
        path = str(path)
        coords_list = load_coords(path)
        
        keep_boxes_msgbox = QtGui.QMessageBox()
        keep_boxes_msgbox.setText(self.tr("Keep current boxes?"))
        keep_boxes_msgbox.setInformativeText(self.tr("Do you want to keep your current boxes?"))
        keep_boxes_msgbox.setStandardButtons(QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        keep_boxes_msgbox.setDefaultButton(QtGui.QMessageBox.Yes)
        keep_current_boxes = keep_boxes_msgbox.exec_()

        if keep_current_boxes == QtGui.QMessageBox.No:
            self.main_image.shapes = EMShapeDict()
            self.set_db_item("boxes", [])
            self.segments_dict = {}
            if self.segment_viewer:
                self.display_segment(EMData(10,10))
        
        for coords in coords_list:
            emshape = EMShape(["rectline", self.color[0], self.color[1], self.color[2], coords[0], coords[1], coords[2], coords[3], coords[4], 2])
            key = self.generate_emshape_key()
            self.main_image.add_shape(key, emshape)
            segment = get_segment_from_coords(self.main_image.get_data(), *coords)
            self.segments_dict[coords] = segment
            self.add_box_to_db(coords)

        self.main_image.updateGL()
        
    def width_changed(self, width):
        """
        updates the widths of the boxed segments when the user changes the width to use for segments
        """
        if width < 1:
            return
        #other widget updates
        self.ptcls_length_spinbox.setValue(width)
        self.ptcls_overlap_spinbox.setValue( int(0.9*width) )
        self.ptcls_width_spinbox.setValue( width )
        
        #resize current boxes
        #TODO: this is similar to part of self.mouse_up ==> make both methods call a function with common code
        shapes = self.main_image.get_shapes() #an EMShapeDict of EMShapes
        for box_key in shapes.keys():
            old_emshape = shapes.get(box_key)
            old_coords = old_emshape.getShape()[4:9]
            new_coords = (old_coords[0], old_coords[1], old_coords[2], old_coords[3], width)
            segment = get_segment_from_coords( self.main_image.get_data(), *new_coords )
                        
            self.remove_box_from_db(old_coords)
            self.add_box_to_db(new_coords)
            self.segments_dict.pop(tuple(old_coords))
            self.segments_dict[new_coords] = segment
                        
            new_emshape = EMShape( ["rectline", self.color[0], self.color[1], self.color[2], new_coords[0], new_coords[1], new_coords[2], new_coords[3], new_coords[4], 2] )
            shapes[box_key] = new_emshape
            
        self.main_image.shapechange=1
        self.main_image.updateGL()
        
        if self.segment_viewer:
            self.display_segment(EMData(10,10))
        
    def write_ouput(self):
        """
        writes the coordinates for the segments, the image data for the segments, and the image data 
        for the particles to files if each of those options are checked
        """
        frame_filename = os.path.basename(self.frame_filepath)
        frame_name = os.path.splitext( frame_filename )[0]
        output_dir = str(self.output_dir_line_edit.text())
        if self.coords_groupbox.isChecked():
            coords_out_type = unicode( self.coords_ftype_combobox.currentText() )
            coords_out_type = self.coords_file_extension_dict[coords_out_type]
            save_coords(self.segments_dict.keys(), frame_name, coords_out_type, output_dir)
        if self.ptcls_groupbox.isChecked():
            i = 1
            for coords_key in self.segments_dict:
                px_overlap = self.ptcls_overlap_spinbox.value()
                px_length = self.ptcls_length_spinbox.value()
                px_width = self.ptcls_width_spinbox.value()
                output_type = self.image_file_extension_dict[unicode(self.ptcls_ftype_combobox.currentText())]
                
                seg = self.segments_dict[coords_key]
                save_particles(seg, i, frame_name, px_overlap, px_length, px_width, output_type, output_dir)
                i += 1
            pass
        if self.segs_groupbox.isChecked():
            seg_file_extension = self.image_file_extension_dict[unicode(self.segs_ftype_combobox.currentText())]
            print seg_file_extension
            i = 1
            for coords_key in self.segments_dict:
                print coords_key
                seg = self.segments_dict[coords_key]
                save_segment(seg, frame_name, i, seg_file_extension, output_dir)
                i += 1
        
    def get_db_item(self, key):
        """
        gets the value stored in the e2helixboxer database for the specified key and the current frame 
        """
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        val = db[self.frame_filepath]
        db_close_dict(db_name)
        return val
    def remove_db_item(self, key):
        """
        removes the key and its value from the e2helixboxer database for the current frame
        """
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        db.pop(key)
    def set_db_item(self, key, value):
        """
        sets the value stored in the e2helixboxer database for the specified key and the current frame 
        """
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        db[self.frame_filepath] = value
        db_close_dict(db_name)
    def get_image_quality(self):
        return self.get_db_item("quality")
    def set_image_quality(self, quality):
        self.set_db_item("quality", quality)
    def add_box_to_db(self, box_coords):
        """
        adds the coordinates for a segment to the e2helixboxer database for the current frame
        """
        assert len(box_coords) == 5, "box_coords must have 5 items"
        db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
        boxList = db[self.frame_filepath] #Get a copy of the db in memory
        boxList.append(tuple(box_coords))
        db[self.frame_filepath] = boxList #Needed to save changes to disk
    def remove_box_from_db(self, box_coords):
        """
        removes the coordinates for a segment in the e2helixboxer database for the current frame
        """
        assert len(box_coords) == 5, "box_coords must have 5 items"
        db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
        boxList = db[self.frame_filepath] #Get a copy of the db in memory
        boxList.remove(tuple(box_coords))
        db[self.frame_filepath] = boxList #Needed to save changes to disk

    def mouse_down(self, event, click_loc):
        """
        If the shift key is pressed and the click is inside a box, delete it.
        Otherwise, either create a new box or edit an existing one depending on click location.
        Imagine drawing two (infinite) lines through the long sides of each box.
        If the click is not between two of the lines for a box, we will create a new box.
        Then the behavior depends on distance from the shorter axis of symmetry--in other
        words, how far up or down the length of the box. Clicking in the middle 3/4 of the box
        (3/8 L from the shorter axis of symmetry) will result in moving the entire box.
        Clicking on a point betwen 3/8 L and 5/8 L from the shorter axis of symmetry
        results in moving that end of the box while keeping the midpoint of the other end fixed.
        
        @event: the mouse click event that causes a box to be added, removed, or modified
        @click_loc: the coordinates in Angstroms of the mouse click on the image
        """

        self.click_loc = click_loc
        box_key = None
        
        if self.main_image.get_shapes(): #helix boxes already exist
            box_key = self.main_image.get_shapes().closest_collision(click_loc[0], click_loc[1], fuzzy=True)
            if event.modifiers()&QtCore.Qt.ShiftModifier:
                if not box_key:
                    self.edit_mode = None #Nothing to delete
                else:
                    box_key = self.main_image.get_shapes().closest_collision(click_loc[0], click_loc[1], fuzzy=False)
                    self.edit_mode = "delete"
            else:
                if not box_key:
                    self.edit_mode = "new"
                else:
                    control_points = self.main_image.get_shapes().get(box_key).control_pts()
                    closest_pt_ix = 0
                    point = control_points[0]
                    min_squared_dist = (click_loc[0] - point[0])**2 + (click_loc[1] - point[1])**2
                    for i in (1,2):
                        point = control_points[i]
                        dist_squared = (click_loc[0] - point[0])**2 + (click_loc[1] - point[1])**2
                        if dist_squared < min_squared_dist:
                            min_squared_dist = dist_squared
                            closest_pt_ix = i
                    if closest_pt_ix == 0: #first endpoint
                        self.edit_mode = "1st_point"
                    elif closest_pt_ix == 1: #second endpoint
                        self.edit_mode = "2nd_point"
                    elif closest_pt_ix == 2: #midpoint
                        self.edit_mode = "move"
                    else:
                        self.edit_mode = "error"
            
        else: #no boxes exist
            if event.modifiers()&QtCore.Qt.ShiftModifier: #nothing to delete
                self.edit_mode = None
            else:
                self.edit_mode = "new" #create new box
       
                
        
        if self.edit_mode == "new" or not self.edit_mode:
            self.current_boxkey = None
            self.initial_helix_box_data_tuple = None
        elif self.edit_mode == "delete":
            box_coords = self.main_image.get_shapes().get(box_key).getShape()[4:9]
            self.remove_box_from_db(box_coords)
            self.main_image.del_shape(box_key)
            self.main_image.updateGL()
            self.current_boxkey = None
        else:
            self.current_boxkey = box_key
            self.initial_helix_box_data_tuple = tuple( self.main_image.get_shapes().get(box_key).getShape()[4:9] )

    def mouse_drag(self, event, cursor_loc):
        """
        Boxes are deleted in mouse_down, and the decision of how to edit is made there.
        However, new boxes are made and existing boxes are edited here.
        @event: the mouse click event that causes a box to be added, removed, or modified
        @click_loc: the coordinates in Angstroms of the mouse click on the image
        """
        
        if self.click_loc and self.edit_mode: #self.click_loc and self.edit_mode are set in mouse_down
            if self.edit_mode == "new":
                if self.click_loc[0] != cursor_loc[0] or self.click_loc[1] != cursor_loc[1]: #Don't make a zero-sized box
                    self.current_boxkey = self.generate_emshape_key()                
                    emshape_tuple = ( "rectline",self.color[0], self.color[1], self.color[2], 
                                           self.click_loc[0], self.click_loc[1], cursor_loc[0], cursor_loc[1], self.get_width(), 2 )
                    
                    emshape_box = EMShape(emshape_tuple)
                    self.main_image.add_shape(self.current_boxkey, emshape_box)
                    self.main_image.updateGL()
                    self.initial_helix_box_data_tuple = emshape_tuple[4:9]
                    self.edit_mode = "2nd_point"
                
            elif self.edit_mode == "delete":
                pass
            else:
                first = self.initial_helix_box_data_tuple[:2]
                second = self.initial_helix_box_data_tuple[2:4]
                width = self.initial_helix_box_data_tuple[4]
                move = (cursor_loc[0] - self.click_loc[0], cursor_loc[1]-self.click_loc[1])

                if self.edit_mode == "move":
                    first = (move[0]+first[0], move[1]+first[1])
                    second = (move[0]+second[0], move[1]+second[1])
                elif self.edit_mode == '1st_point': #move first point
                    first = (move[0]+first[0], move[1]+first[1])
                elif self.edit_mode == "2nd_point":
                    second = (move[0]+second[0], move[1]+second[1])
                
                box = self.main_image.get_shapes().get(self.current_boxkey)
                box.getShape()[4] = first[0]
                box.getShape()[5] = first[1]
                box.getShape()[6] = second[0]
                box.getShape()[7] = second[1]
                self.main_image.shapechange=1
                self.main_image.updateGL()

    def mouse_up(self, event, cursor_loc):
        """
        Once the mouse button comes back up, creating a new box, or editing
        an existing box is complete, so we need only clear variables relevant
        to creating or editing boxes, and get the image data from the boxed area.
        @event: the mouse click event that causes a box to be added, removed, or modified
        @click_loc: the coordinates in Angstroms of the mouse click on the image
        """

        if self.current_boxkey and self.edit_mode != "delete": 
            if self.initial_helix_box_data_tuple in self.get_db_item("boxes"):
                self.remove_box_from_db(self.initial_helix_box_data_tuple)
            box = self.main_image.get_shapes().get(self.current_boxkey)
            box_coords = box.getShape()[4:9]
            self.add_box_to_db(box_coords)
            segment = get_segment_from_coords( self.main_image.get_data(), *box_coords )
            data_tuple = tuple(box.getShape()[4:9])
            self.segments_dict[data_tuple] = segment
            
            self.display_segment(segment)
        
        self.click_loc = None
        self.edit_mode = None
        self.current_boxkey = None #We are done editing the box
        self.initial_helix_box_data_tuple = None

        
def main():
    progname = os.path.basename(sys.argv[0])
    usage = """%prog [options] <image>....

used to box alpha helices

For example:

e2helixboxer.py ????.mrc --boxwidth=256
"""
    parser = OptionParser(usage=usage,version=EMANVERSION)
    parser.add_option("--boxwidth","-B",type="int",help="Box width in pixels",default=128)
    (options, args) = parser.parse_args()
    if len(args) > 0:
        filename= args[0]
        (path, basename) = os.path.split(filename)
    else:
        filename = None
        (path, basename) = (None, None)
    logid=E2init(sys.argv)
    
    app = EMStandAloneApplication()
    helixboxer = EMHelixBoxerWidget(filename, app)
    helixboxer.show()
    app.execute()
    E2end(logid)

if __name__ == '__main__':
    main()

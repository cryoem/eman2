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
    if not px_width:
        px_width = segment.get_xsize()
    if not px_length:
        px_length = px_width
    assert px_length > px_overlap, "The overlap must be smaller than the particle length"
    
    px_step = px_length - px_overlap
    
    (xsize, ysize) = (px_width, px_length)
    (x0,y0) = (0,0)
    particles = []
    while y0 + ysize <= segment.get_xsize():
        particles.append(segment.get_clip(x0, y0, xsize, ysize))
        y0 += px_step
    
    return particles

def get_segment_from_coords(frame, x1, y1, x2, y2, width):
    l_vect = (x2-x1, y2-y1)
    length = sqrt(l_vect[0]**2+l_vect[1]**2)
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

def save_coords(coords_list, frame_name, output_type = "box", output_dir=""):
    if output_type == "box":
        path = os.path.join(output_dir, frame_name + ".box")
        out_file = open(path, "w")
        for coords in coords_list:
            out_file.write( "%i\t%i\t%i\t%i\t%i\n" % (coords[0], coords[1], coords[4], coords[4], -1) )
            out_file.write( "%i\t%i\t%i\t%i\t%i\n" % (coords[2], coords[3], coords[4], coords[4], -2) )
        out_file.close()
    elif output_type == "hbox":
        path = os.path.join(output_dir, frame_name + ".hbox")
        out_file = open(path, "w")
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
    segments_dict = {}
    for coords in box_coords_list:
        segment = get_segment_from_coords(frame, *coords)
        segments_dict[tuple(coords)] = segment
    return segments_dict

def save_segment(segment_emdata, frame_name, segment_num, output_type = "hdf", output_dir=""):
    output_fname = "%s_seg%i.%s" % (frame_name, segment_num, output_type)
    output_filepath = os.path.join(output_dir, output_fname)
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

def save_particles(segment, frame_name, px_overlap, px_length = None, px_width = None, output_dir=""):
    particles = get_particles_from_segment(segment, px_overlap, px_length, px_width)
    #TODO: save as a stack!

def db_save_particles(frame_filepath, px_overlap, px_length = None, px_width = None, output_type = "hdf", output_dir=""):
    pass

class EMHelixBoxerWidget(QtGui.QWidget):
    def __init__(self, frame_filepath, app):
        QtGui.QWidget.__init__(self)
        self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
        self.setWindowTitle("e2helixboxer")
        
        if not frame_filepath:
            basename = "test_image"
            img = test_image()
        else:
            (path, basename) = os.path.split(frame_filepath)
            img = EMData(frame_filepath)
        
        self.filename = basename
        self.main_image = EMImage2DModule(img, application=app)
        self.main_image.set_file_name(basename) # TODO: determine if this should use the entire file path
        self.main_image.shapes = EMShapeDict()
        self.segment_viewer = None #Will be an EMImage2DModule instance
        self.edit_mode = None #Values are in {None, "new", "move", "2nd_point", "1st_point", "delete"}
        self.current_boxkey = None
        self.initial_helix_box_data_tuple = None
        self.click_loc = None #Will be (x,y) tuple
        self.segments_dict = db_get_segments_dict(self.filename) #Will be like {(x1,y1,x2,y2,width): emdata}
        self.color = (1, 1, 1)
        self.counter = counterGen()
        self.coords_file_extension_dict = {"EMAN1":"box", "EMAN2": "hbox"}
        self.image_file_extension_dict = {"MRC":"mrc", "Spider":"spi", "Imagic": "img", "HDF5": "hdf"}

        if self.get_db_item("boxes") == None:
            self.set_db_item("boxes", [])
        else:
            boxList = self.get_db_item("boxes")
            for boxCoords in boxList:
                key = self.generate_emshape_key()
                emshape_list = ["rectline"]
                emshape_list.extend(list(self.color))
                emshape_list.extend(list(boxCoords))
                emshape_list.append(2)
                emshape = EMShape( emshape_list )
                self.main_image.add_shape(key, emshape)
            self.main_image.updateGL()

        self.__create_ui()
        
        self.coords_ftype_combobox.addItems( sorted(self.coords_file_extension_dict.keys()) )
        self.segs_ftype_combobox.addItems( sorted(self.image_file_extension_dict.keys()) )
        self.ptcls_ftype_combobox.addItems( sorted(self.image_file_extension_dict.keys()) )
        width = 100
        self.box_width_spinbox.setValue(width)
        self.ptcls_width_spinbox.setValue( width )
        self.ptcls_length_spinbox.setValue( width )
        self.ptcls_overlap_spinbox.setValue( int(0.9*width) )
        
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedown"), self.mouse_down)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedrag"), self.mouse_drag)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mouseup"), self.mouse_up)
        self.connect(self.output_dir_pushbutton, QtCore.SIGNAL("clicked()"), self.choose_dir )
        self.connect(self.write_output_button, QtCore.SIGNAL("clicked()"), self.write_ouput )
        self.connect(self.box_width_spinbox, QtCore.SIGNAL("valueChanged(int)"), self.width_changed)
        self.connect(self.done_but, QtCore.SIGNAL("clicked()"), self.exit_app )

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
        self.img_quality_combobox.setEnabled(False)
        
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
        
        self.ptcls_groupbox.setChecked(False)
        self.ptcls_groupbox.setEnabled(False)
        
        
        
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
        self.vbl.addWidget(self.coords_groupbox)
        self.vbl.addWidget(self.segs_groupbox)
        self.vbl.addWidget(self.ptcls_groupbox)
        self.vbl.addLayout(directory_layout)
        self.vbl.addLayout(button_layout)
        self.vbl.addWidget(self.done_but)
        self.vbl.addWidget(self.status_bar)

    def box_coords_string(self):
        print self.main_image.get_shapes()
        coords_list = [ shape.getShape()[4:9] for shape in self.main_image.get_shapes().values() ]
        lines = [ "\t".join( [str(i) for i in coords] ) for coords in coords_list ]
        ret = "\n".join(lines)
        return ret
    def choose_dir(self):
        #selector = EMSelectorModule(save_as_mode=False)
        #selector.widget.save_as_line_edit.setEnabled(False)
        #path = selector.exec_()
        
        path = QtGui.QFileDialog.getExistingDirectory(self)
        
        self.output_dir_line_edit.setText(path)
    def exit_app(self):
        app = get_application()
        app.quit()
    def generate_emshape_key(self):
        i = self.counter.next()
        return "rectline%i" % i
    def get_width(self):
        return self.box_width_spinbox.value()
    def width_changed(self, width):
        self.ptcls_length_spinbox.setValue(width)
        self.ptcls_overlap_spinbox.setValue( int(0.9*width) )
        self.ptcls_width_spinbox.setValue( width )
    def write_coords(self):
        em_selector_module = EMSelectorModule()
        file_path = em_selector_module.exec_()
        print file_path
        print self.box_coords_string()
    def write_particles(self):
        em_selector_module = EMSelectorModule()
        file_path = em_selector_module.exec_()
        print file_path
        pass
    def write_ouput(self):
        frame_filename = os.path.basename(self.filename)
        frame_name = os.path.splitext( frame_filename )[0]
        if self.coords_groupbox.isChecked():
            coords_out_type = unicode( self.coords_ftype_combobox.currentText() )
            coords_out_type = self.coords_file_extension_dict[coords_out_type]
            save_coords(self.segments_dict.keys(), frame_name, coords_out_type)
        if self.ptcls_groupbox.isChecked():
            #save_particles(self.filename)
            pass
        if self.segs_groupbox.isChecked():
            seg_file_extension = self.image_file_extension_dict[unicode(self.segs_ftype_combobox.currentText())]
            print seg_file_extension
            i = 1
            for coords_key in self.segments_dict:
                print coords_key
                seg = self.segments_dict[coords_key]
                save_segment(seg, frame_name, i, seg_file_extension)
                i += 1
        
    def get_db_item(self, key):
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        val = db[self.filename]
        db_close_dict(db_name)
        return val
    def remove_db_item(self, key):
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        db.pop(key)
    def set_db_item(self, key, value):
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        db[self.filename] = value
        db_close_dict(db_name)
    def add_box_to_db(self, boxCoords):
        assert len(boxCoords) == 5, "boxCoords must have 5 items"
        db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
        boxList = db[self.filename] #Get a copy of the db in memory
        boxList.append(tuple(boxCoords))
        db[self.filename] = boxList #Needed to save changes to disk
    def remove_box_from_db(self, boxCoords):
        assert len(boxCoords) == 5, "boxCoords must have 5 items"
        db = db_open_dict(E2HELIXBOXER_DB + "#boxes")
        boxList = db[self.filename] #Get a copy of the db in memory
        boxList.remove(tuple(boxCoords))
        db[self.filename] = boxList #Needed to save changes to disk

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
            boxCoords = self.main_image.get_shapes().get(box_key).getShape()[4:9]
            self.remove_box_from_db(boxCoords)
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
        """
        
        if self.click_loc and self.edit_mode: #self.click_loc and self.edit_mode are set in mouse_down
            if self.edit_mode == "new":
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
        """

        if self.current_boxkey and self.edit_mode != "delete":
            if self.initial_helix_box_data_tuple in self.get_db_item("boxes"):
                self.remove_box_from_db(self.initial_helix_box_data_tuple)
            box = self.main_image.get_shapes().get(self.current_boxkey)
            boxCoords = box.getShape()[4:9]
            self.add_box_to_db(boxCoords)
            segment = get_segment_from_coords( self.main_image.get_data(), *boxCoords )
            data_tuple = tuple(box.getShape()[4:9])
            self.segments_dict[data_tuple] = segment
            
            if not self.segment_viewer:
                self.segment_viewer = EMImage2DModule(application=get_application())
                self.segment_viewer.desktop_hint = "rotor" # this is to make it work in the desktop
                self.segment_viewer.setWindowTitle("Current Boxed Segment")
                self.segment_viewer.get_qt_widget().resize(200,800)
            self.segment_viewer.set_data(segment)
            #w = self.segment_viewer.width()
            #self.segment_viewer.set_origin(w/2,0)
            self.segment_viewer.updateGL()
            get_application().show_specific(self.segment_viewer)
            scale = 100 / self.get_width()
            self.segment_viewer.set_scale(scale)
            if self.segment_viewer.inspector:
                self.segment_viewer.inspector.set_scale(scale)
            self.segment_viewer.updateGL()

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

#!/usr/bin/env python

#Author: Ross Coleman

from EMAN2 import test_image, get_image_directory, Transform, Region, EMANVERSION, EMData, E2init, E2end
from EMAN2db import db_open_dict, db_check_dict, db_close_dict
from emapplication import EMStandAloneApplication, get_application
from emimage2d import EMImage2DModule
#from emimagemx import EMImageMXModule
from emshape import EMShape, EMShapeDict
from optparse import OptionParser

from math import *
import weakref, sys, os
from PyQt4 import QtGui, QtCore


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
        
class EMHelixBoxerWidget(QtGui.QWidget):
    def __init__(self, filename, app):
        QtGui.QWidget.__init__(self)
        self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
        self.setWindowTitle("e2helixboxer")
        
        if not filename:
            filename = "test_image"
            image = test_image()
        img = EMData(filename)
        self.filename = filename
        self.main_image = EMImage2DModule(img, application=app)
        self.main_image.set_file_name(filename) # TODO: set the filename
        self.main_image.shapes = EMShapeDict()
        self.particle_viewer = None #Will be an EMImage2DModule instance
        #self.box_list = HelixBoxList()
        self.edit_mode = None #Values are in {None, "new", "move", "2nd_point", "1st_point", "delete"}
        self.current_boxkey = None
        self.initial_helix_box_data_tuple = None
        self.click_loc = None #Will be (x,y) tuple
        self.particles_dict = {} #Will be like {(x1,y1,x2,y2,width): emdata}
        self.color = (1, 1, 1)
        self.counter = counterGen()
        
        self.__create_ui()
        
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedown"), self.mouse_down)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mousedrag"), self.mouse_drag)
        QtCore.QObject.connect( self.main_image.emitter(), QtCore.SIGNAL("mouseup"), self.mouse_up)
        self.connect( self.gen_output_but, QtCore.SIGNAL("clicked()"), self.box_coords_string )

        get_application().show_specific(self.main_image)
        
        
    def __create_ui(self):
        self.boxWidthLabel = QtGui.QLabel(self.tr("Box &Width"))
        self.boxWidthSpinBox = QtGui.QSpinBox()
        self.boxWidthSpinBox.setMaximum(1000)
        self.boxWidthSpinBox.setValue(100)
        self.boxWidthLabel.setBuddy(self.boxWidthSpinBox)
        
        self.imgQualityLabel = QtGui.QLabel(self.tr("Image &Quality"))
        self.imgQualityComboBox = QtGui.QComboBox()
        qualities = [str(i) for i in range(5)]
        self.imgQualityComboBox.addItems(qualities)
        self.imgQualityComboBox.setCurrentIndex(2)
        self.imgQualityLabel.setBuddy(self.imgQualityComboBox)
             
        self.gen_output_but=QtGui.QPushButton(self.tr("&Write Output"))
        self.done_but=QtGui.QPushButton(self.tr("&Done"))
        
        self.status_bar = QtGui.QStatusBar()
        self.status_bar.showMessage("Ready",10000)
        
        self.imgQualityComboBox.setEnabled(False)
        self.done_but.setEnabled(False)
        #self.gen_output_but.setEnabled(False)
        
        widthLayout = QtGui.QHBoxLayout()
        widthLayout.addWidget(self.boxWidthLabel)
        widthLayout.addWidget(self.boxWidthSpinBox)
        
        qualityLayout = QtGui.QHBoxLayout()
        qualityLayout.addWidget(self.imgQualityLabel)
        qualityLayout.addWidget(self.imgQualityComboBox)
        
        self.vbl = QtGui.QVBoxLayout(self)
        self.vbl.setMargin(0)
        self.vbl.setSpacing(6)
        self.vbl.setObjectName("vbl")
        self.vbl.addLayout(widthLayout)
        self.vbl.addLayout(qualityLayout)        
        self.vbl.addWidget(self.gen_output_but)
        self.vbl.addWidget(self.done_but)
        self.vbl.addWidget(self.status_bar)
        
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
            self.main_image.del_shape(box_key)
            self.main_image.updateGL()
            self.current_boxkey = None
        else:
            self.current_boxkey = box_key
            self.initial_helix_box_data_tuple = list( self.main_image.get_shapes().get(box_key).getShape()[4:] )

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
                self.initial_helix_box_data_tuple = emshape_tuple[4:]
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
            box = self.main_image.get_shapes().get(self.current_boxkey)
            
            #Now we'll get the area of the image that was boxed, and display it
            control_pts = box.control_pts()
            box_centroid = control_pts[2]
            p0 = control_pts[0]
            p1 = control_pts[1]
            l_vect = (p1[0]-p0[0], p1[1]-p0[1])
            box_length = sqrt(l_vect[0]**2+l_vect[1]**2)
            box_width = box.getShape()[8]
            l_uvect = (l_vect[0]/box_length, l_vect[1]/box_length)
            #Rotate so that the length is parallel to the y-axis
            #Angle between l_uvect and y-axis: l_uvect (dot) j_hat = cos (rot_angle)
            rot_angle = 180/pi*acos( l_uvect[1] )
            #Whether we rotate clockwise or counterclockwise depends on the sign of l_uvect[0] (the x-component)
            if l_uvect[0] < 0:
                rot_angle *= -1 #We want to rotate the box clockwise, so the angle is negative

            em_image = self.main_image.get_data()
            tr = Transform()
            tr.set_trans(box_centroid)
            tr.set_rotation({"type":"2d", "alpha":rot_angle})
            particle_dimensions = ( int(round(box_width)), int(round(box_length)), 1 )
            particle = em_image.get_rotated_clip( tr, particle_dimensions )
            data_tuple = tuple(box.getShape()[4:9])
            self.particles_dict[data_tuple] = particle
            self.set_db_item("boxes", self.particles_dict.keys() )
            
            if not self.particle_viewer:
                self.particle_viewer = EMImage2DModule(application=get_application())
                self.particle_viewer.desktop_hint = "rotor" # this is to make it work in the desktop
                self.particle_viewer.setWindowTitle("Current Boxed Particle")
                self.particle_viewer.get_qt_widget().resize(200,800)
            self.particle_viewer.set_data(particle)
            #w = self.particle_viewer.width()
            #self.particle_viewer.set_origin(w/2,0)
            self.particle_viewer.updateGL()
            get_application().show_specific(self.particle_viewer)
            scale = 100 / self.get_width()
            self.particle_viewer.set_scale(scale)
            if self.particle_viewer.inspector:
                self.particle_viewer.inspector.set_scale(scale)
            self.particle_viewer.updateGL()

        self.click_loc = None
        self.edit_mode = None
        self.current_boxkey = None #We are done editing the box
        self.initial_helix_box_data_tuple = None
    def generate_emshape_key(self):
        i = self.counter.next()
        return "rectline%i" % i
    def get_width(self):
        return self.boxWidthSpinBox.value()
    def get_db_item(self, key):
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        val = db[self.filename]
        db_close_dict(db_name)
        return val
    def set_db_item(self, key, value):
        db_name = E2HELIXBOXER_DB + "#" + key
        db = db_open_dict(db_name)
        db[self.filename] = value
        db_close_dict(db_name)
    def box_coords_string(self):
        lines = [ "\t".join( [str(i) for i in coords_size_tuple] ) for coords_size_tuple in self.particles_dict.keys() ]
        ret = "\n".join(lines)
        print ret
        return ret

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
    else:
        filename = None
    logid=E2init(sys.argv)
    db = db_open_dict(E2HELIXBOXER_DB)
    
    app = EMStandAloneApplication()
    helixboxer = EMHelixBoxerWidget(filename, app)
    helixboxer.show()
    app.execute()
    E2end(logid)

if __name__ == '__main__':
    main()

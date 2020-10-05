try:
  from pyStarDB import sp_pystardb as star
  STAR_AVAILABLE = True
except ImportError:
  STAR_AVAILABLE = False

import EMAN2
import EMAN2db
import EMAN2_cppwrap
from EMAN2 import Transform
from EMAN2 import EMUtil
from EMAN2 import EMAN2Ctf

import time
import unittest
import numpy  as np


inputfile = '/home/adnan/Desktop/star_file_project/particles_optics.star'
inputfile = '/home/adnan/PycharmProjects/newrelion/Refine3D/job056/run_data.star'
outputfile = '/home/adnan/Desktop/star_file_project/newtest.star'
bdbfile = 'bdb:/home/adnan/PycharmProjects/DoseWeighting/Relion_Stackv5/MotionCorr/job049/MOVIES_RELION/Particles/sphire_relion_stack'

star_file = star.StarFile(inputfile)
try:
  star_file['particles']
except:
  star_file['']


class Test_EMDB_functions(unittest.TestCase):

  def test_get_all_attributes(self):
    attributes = EMAN2db.db_get_all_attributes(inputfile, "ctf")
    self.assertEqual(str(attributes[0].to_dict()), str({'ampcont': 10.0, 'apix': 0.8849999904632568, 'background': [], 'bfactor': 0.0, 'cs': 1.399999976158142, 'defocus': 1.0729365348815918, 'dfang': 149.37908935546875, 'dfdiff': -0.026136621832847595, 'dsbg': 0.0, 'snr': [], 'voltage': 200.0} ))

    attributes = EMAN2db.db_get_all_attributes(inputfile, 'ptcl_source_coord')
    self.assertEqual(attributes[0], [2937.0, 3437.0])

    attributes = EMAN2db.db_get_all_attributes(inputfile, 'ptcl_source_image')
    self.assertEqual(attributes[0], 'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc')


  def test_get_image_count(self):
    count_info  = EMAN2db.db_get_image_count(inputfile)
    self.assertEqual(count_info , 5090)

  def test_read_image(self):
    img = EMAN2.EMData()
    img.read_image(inputfile, 0)

    self.assertEqual(str(img.get_attr('ctf').to_dict()), str(
      {'ampcont': 10.0, 'apix': 0.8849999904632568, 'background': [], 'bfactor': 0.0, 'cs': 1.399999976158142,
       'defocus': 1.0729365348815918, 'dfang': 149.37908935546875, 'dfdiff': -0.026136621832847595, 'dsbg': 0.0,
       'snr': [], 'voltage': 200.0}))
    self.assertEqual(img.get_attr('ptcl_source_coord'), [2937.0, 3437.0])
    self.assertEqual(img.get_attr('ptcl_source_image'),
                     'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc')

  def test_read_images(self):
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, [0, 4, 5, 7])

    self.assertEqual(len(pp), 4)
    self.assertEqual(str(pp[0].get_attr('ctf').to_dict()), str(
      {'ampcont': 10.0, 'apix': 0.8849999904632568, 'background': [], 'bfactor': 0.0, 'cs': 1.399999976158142,
       'defocus': 1.0729365348815918, 'dfang': 149.37908935546875, 'dfdiff': -0.026136621832847595, 'dsbg': 0.0,
       'snr': [], 'voltage': 200.0}))
    self.assertEqual(pp[0].get_attr('ptcl_source_coord'), [2937.0, 3437.0])
    self.assertEqual(pp[0].get_attr('ptcl_source_image'),
                     'MotionCorr/job049/MOVIES_RELION/20170629_00021_frameImage.mrc')


  def test_write_image(self):
    img = EMAN2.EMData()
    dd = [0, 2, 4, 5, 7]
    pp = img.read_images(inputfile, dd)

    for i, emdata in enumerate(pp):
         emdata.write_image(outputfile, dd[i])

    img1 = EMAN2.EMData()
    dd = [0, 1, 2, 3, 4]   # this is index of the row not the particle number
    pp_new = img1.read_images(outputfile, dd)

    for i in range(len(pp_new)):
      self.assertTrue(np.allclose(pp[i].get_2dview() , pp_new[i].get_2dview(), rtol =0.5))


  def test_read_header(self):
    dd = np.arange(5000).tolist()
    inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    start_time1 = time.time()
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)
    print("Star file read time--- %s seconds ---" % (time.time() - start_time1))

    start_time2 = time.time()
    img1 = EMAN2.EMData()
    pp1 = img1.read_images(bdbfile, dd, True)
    print("BDB file time--- %s seconds ---" % (time.time() - start_time2))




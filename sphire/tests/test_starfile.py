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
import os

import time
import unittest
import numpy  as np

# path = os.path.join(os.path.dirname(__file__), "../resources_test/TcdA1-0187_frames_ptcls.star.star")

inputfile = os.path.join(os.path.dirname(__file__), "../resources_tests/03_PARTICLES/mpi_proc_007/TcdA1-0187_frames_ptcls.star")
second_file = os.path.join(os.path.dirname(__file__), "../resources_tests/03_PARTICLES/mpi_proc_007/TcdA1-0188_frames_ptcls.star")
outputfile = os.path.join(os.path.dirname(__file__), "../resources_tests/output_test.star")
bdbfile = 'bdb:'+os.path.join(os.path.dirname(__file__), "../resources_tests/03_PARTICLES_BDB/mpi_proc_007/TcdA1-0187_frames_ptcls")

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
    self.assertEqual(count_info , 22)

  def test_read_image(self):
    img = EMAN2.EMData()
    img.read_image(inputfile, 0, True)
    print(img.get_attr('ptcl_source_coord'))

    self.assertEqual(img.get_attr('ptcl_source_coord'), [241,3868])
    self.assertEqual(img.get_attr('ptcl_source_image'),
                     'CorrectedSums/corrsum_dw/TcdA1-0187_frames.mrc')

  def test_read_images(self):
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, [0, 4, 5, 7])

    self.assertEqual(len(pp), 4)
    self.assertEqual(str(pp[0].get_attr('ctf').to_dict()), str(
      {'ampcont': 10.0, 'apix': 1.1399999856948853, 'background': [], 'bfactor': 0.0, 'cs': 0.0, 'defocus': 1.3080999851226807, 'dfang': 330.41998291015625, 'dfdiff': 0.034550998359918594, 'dsbg': 0.0, 'snr': [], 'voltage': 300.0}))
    self.assertEqual(pp[0].get_attr('ptcl_source_coord'), [241, 3868])
    self.assertEqual(pp[0].get_attr('ptcl_source_image'),
                     'CorrectedSums/corrsum_dw/TcdA1-0187_frames.mrc')


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
    dd = np.arange(5).tolist()
    # inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    start_time1 = time.time()
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)
    print("Star file read time--- %s seconds ---" % (time.time() - start_time1))

    # start_time2 = time.time()
    # img1 = EMAN2.EMData()
    # pp1 = img1.read_images(bdbfile, dd, True)
    # print("BDB file time--- %s seconds ---" % (time.time() - start_time2))



  def test_writing_test_again(self):
    dd = np.arange(50).tolist()
    # inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    start_time1 = time.time()
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)

    for i in range(len(pp)):
      pp[i].write_image(outputfile, i)


  def test_writing_images(self):
    print("You are awesome")
    dd = np.arange(50).tolist()
    # inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)

    # start_time1 = time.time()
    # for i in range(len(pp)):
    #   pp[i].write_image('/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/saved_in_loop.star', i)
    # print("Time used for saving data using write_image", time.time() - start_time1)

    start_time2 = time.time()
    # outputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/saved_in_onego.star'
    EMAN2db.write_images(pp, outputfile, np.arange(50).tolist() )

    print("Time used for saving data using write_images", time.time() - start_time2)
    print("Execution done")


  def test_bdb_star(self):

    dd = np.arange(5).tolist()
    img1 = EMAN2.EMData()
    bdb_stack = img1.read_images(bdbfile, dd, False)

    img2 = EMAN2.EMData()
    pp = img2.read_images(inputfile, dd, False)

    found_lost_keys = ['is_complex_ri']
    bdb_dont_have = ['data_path']

    for ind in range(len(pp)):
      for key, val in bdb_stack[ind].get_attr_dict().items():
        print(key)
        if key == "changecount" :
          print("Change count value")
          print("Value in BDB ", val)
          print("Value in star",  pp[ind].get_attr_dict()[key])

        if key == 'ctf':
          dict_bdb = bdb_stack[ind]['ctf'].to_dict()
          dict_star = pp[ind]['ctf'].to_dict()
      self.assertTrue(np.array_equal(pp[ind].get_2dview() , bdb_stack[ind].get_2dview() ))

  def test_bdb_star_reading(self):

    dd = np.arange(5).tolist()
    img2 = EMAN2.EMData()
    # pp1 = img2.read_images(star_file_loc)
    pp2 = img2.read_images(inputfile, dd, False)
    pp4 = img2.read_images(inputfile, dd)
    img1 = EMAN2.EMData()
    img2 = EMAN2.EMData()
    img1.read_image(inputfile, 0)
    img2.read_image(inputfile)
    self.assertTrue(np.array_equal(pp2[0].get_2dview(), pp4[0].get_2dview()))

    self.assertTrue(np.array_equal(img1.get_2dview(), img2.get_2dview()))


  def test_star_two_file_writing(self):

    img2 = EMAN2.EMData()
    img1 = EMAN2.EMData()
    img1.read_image(inputfile, 0)
    img2.read_image(second_file, 0)
    img1.write_image(outputfile, 0)
    img2.write_image(outputfile, 1)
    img2.write_image(outputfile, 2)

    pp = [img1, img2]
    EMAN2.EMData.write_images(pp, outputfile, [0,1])















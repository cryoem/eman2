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
    img.read_image(inputfile, 0, True)

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



  def test_writing_test_again(self):
    dd = np.arange(5000).tolist()
    inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    start_time1 = time.time()
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)

    for i in range(len(pp)):
      pp[i].write_image('/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/newstarfile.star', i)


  def test_writing_images(self):
    print("You are awesome")
    dd = np.arange(5000).tolist()
    inputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/test_starfile.star'
    img = EMAN2.EMData()
    pp = img.read_images(inputfile, dd, True)

    # start_time1 = time.time()
    # for i in range(len(pp)):
    #   pp[i].write_image('/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/saved_in_loop.star', i)
    #
    # print("Time used for saving data using write_image", time.time() - start_time1)


    start_time2 = time.time()
    outputfile = '/home/adnan/PycharmProjects/starfile_conv/sphire/tests/resources/saved_in_onego.star'
    EMAN2db.write_images(pp, outputfile, np.arange(5000).tolist() )

    print("Time used for saving data using write_images", time.time() - start_time2)
    print("Execution done")


    # for a , b in enumerate(range(25,55)):
    #   print(a,b)



  def test_bdb_star(self):
    bdb_file_loc = 'bdb:/home/adnan/PycharmProjects/Starfile_test_demo/04_ISAC_BDB/stack_ali2d'

    star_file_loc = '/home/adnan/PycharmProjects/Starfile_test_demo/04_ISAC/stack_ali2d.star'

    # star_file_loc = 'bdb:/home/adnan/PycharmProjects/Starfile_test_demo/04_ISAC/stack_ali2d'

    # bdb_file_loc = 'bdb:/home/adnan/PycharmProjects/Starfile_test_demo/03_PARTICLES_BDB/data'
    #
    # star_file_loc = '/home/adnan/PycharmProjects/Starfile_test_demo/03_PARTICLES_STAR/data.star'
    #

    dd = np.arange(50).tolist()
    img1 = EMAN2.EMData()
    bdb_stack = img1.read_images(bdb_file_loc, dd, False)

    img2 = EMAN2.EMData()
    pp = img2.read_images(star_file_loc, dd, False)

    # print(bdb_stack[0].get_attr('ctf'))
    # print(pp[0].get_attr('ctf'))

    found_lost_keys = ['is_complex_ri'
                       ]
    bdb_dont_have = ['data_path']


    import inspect
    def f1():
     f2()

    def f2():
      print('caller name', inspect.stack())

    f1()


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

          for ctfkey , ctfvalue in  dict_bdb.items():
            print(ctfkey , ctfvalue)
            print(ctfkey, dict_star[ctfkey])

        try:
          if val == pp[ind].get_attr_dict()[key]:
            print('They are same')
          else:
            print(val, pp[ind].get_attr_dict()[key])
        except KeyError :
          print('Not present' , val)
        print('')
      # print(bdb_stack[ind].get_2dview()[50:150, 50:150])
      self.assertTrue(np.array_equal(pp[ind].get_2dview() , bdb_stack[ind].get_2dview() ))

    # for ind , micro in enumerate(pp):
    #   for keys in bdb_stack[ind].get_attr_dict():

        # if keys in found_lost_keys:
        #   pass
        # elif keys in star.StarFile(star_file_loc).ignored_keys:
        #   continue
        #
        # elif keys in bdb_dont_have:
        #   pass
        # else:
        #   if keys == 'ctf':
        #     # print(pp[ind][keys].to_dict())
        #     # print(bdb_stack[ind][keys].to_dict())
        #     pass
        #
        #   else:
        #     # print(keys)
        #     # print(pp[ind][keys])
        #     # print(bdb_stack[ind][keys])
        #     self.assertEqual(pp[ind][keys] , bdb_stack[ind][keys])

  def test_bdb_star_writing(self):

    bdb_file_loc = 'bdb:/home/adnan/PycharmProjects/Starfile_test_demo/03_PARTICLES_NEWBDB/data'
    star_file_loc = '/home/adnan/PycharmProjects/Starfile_test_demo/03_PARTICLES_STAR/data.star'

    dd = np.arange(50).tolist()
    img2 = EMAN2.EMData()
    pp = img2.read_images(star_file_loc, dd, False)

    # for i in range(len(pp)):
    EMAN2db.write_images(pp, bdb_file_loc,range(len(pp)))






try:
  from starfile_db import sp_starfile_db as star
  STAR_AVAILABLE = True
except ImportError:
  STAR_AVAILABLE = False

import EMAN2
import EMAN2db
import EMAN2_cppwrap


from EMAN2 import EMUtil

from EMAN2 import EMAN2Ctf

# elif fsp.endswith('.star'):
#   if not STAR_AVAILABLE:
#     print('READING STAR FILES NOT AVAILABLE')
#     return
#   else:


inputfile = '/home/adnan/Desktop/star_file_project/particles_optics.star'

inputfile = '/home/adnan/PycharmProjects/newrelion/Refine3D/job056/run_data.star'


outputfile = '/home/adnan/Desktop/star_file_project/newtest.star'
star_file = star.StarFile(inputfile)


try:
  star_file.imported_content['particles']
except:
  star_file.imported_content['']

# img = EMAN2.EMData()
# img.read_image(inputfile, 0 )



pp = EMUtil.get_all_attributes(inputfile, 'ctf')

# a = EMAN2.EMData.read_images(inputfile)
#
#
#
# a[0].write_image(outputfile)


print(STAR_AVAILABLE)
print(star_file)


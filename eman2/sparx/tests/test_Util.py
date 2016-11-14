
import unittest
from optparse import OptionParser

IS_TEST_EXCEPTION = False

# ====================================================================================================================
class Test_Util_diff_between_matrix_of_3D_parameters_angles(unittest.TestCase):
	"""this is unit test of diff_between_matrix_of_3D_parameters_angles function from Util"""
	
	def internal_calculate_avg_diff(self, projs, matrix_rot):
		from utilities import angle_between_projections_directions
		n = len(projs[0])
		sc = len(projs)
		matrix_diff = [0]*sc
		for i in xrange(sc):
			matrix_diff[i] = [0]*i
		for i in xrange(sc):
			for j in xrange(i):
				diff = []
				for k in xrange(n):
					
					diff.append( angle_between_projections_directions( self.mult_transform(projs[i][k],matrix_rot[i][j]), projs[j][k] ) )
				matrix_diff[i][j] = diff
		avg_diff_per_image = []
		for k in xrange(n):
			avg = 0.0
			for i in xrange(sc):
				for j in xrange(i):
					avg += matrix_diff[i][j][k]
			avg /= (sc * (sc-1) / 2)
			avg_diff_per_image.append(avg)
		return avg_diff_per_image

	def internal_calculate_avg_diff_Util(self, projs, matrix_rot):
		from global_def import Util
		sc = len(projs)
		trans_matrix = []
		for i in xrange(sc):
			for j in xrange(i):
				trans_matrix.extend(matrix_rot[i][j][0:3])
		trans_projs = []
		for iConf in xrange(sc):
			for i in xrange(len(projs[0])):
				ttt = projs[iConf][i][0:5]
				while len(ttt) < 5:
					ttt.append(0.0)
				trans_projs.extend(ttt)
		avg_diff_per_image = Util.diff_between_matrix_of_3D_parameters_angles(trans_projs, trans_matrix)
		return avg_diff_per_image

	def wrap_rotation_between_anglesets(self, ang1, ang2):
		from utilities import rotation_between_anglesets
		
		phi, theta, psi = rotation_between_anglesets(ang1, ang2)
		return [phi, theta, psi]
	
	def mult_transform(self, v1, v2):
		from EMAN2 import Transform
		T1 = Transform({"type":"spider","phi":v1[0],"theta":v1[1],"psi":v1[2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
		T2 = Transform({"type":"spider","phi":v2[0],"theta":v2[1],"psi":v2[2],"tx":0.0,"ty":0.0,"tz":0.0,"mirror":0,"scale":1.0})
		T = T1*T2
		return [ T.get_params("spider")["phi"], T.get_params("spider")["theta"], T.get_params("spider")["psi"], T.get_params("spider")["tx"], T.get_params("spider")["ty"]  ]
	
	
	def calculate_matrix_rot(self, projs):
		sc = len(projs)
		matrix_rot  = [0]*sc
		for i in xrange(sc):
			matrix_rot [i] = [0]*i
		for i in xrange(sc):
			for j in xrange(i):
				matrix_rot[i][j] = self.wrap_rotation_between_anglesets(projs[i], projs[j])
		return matrix_rot


	def test_two_configurations(self):
		"""test_two_configurations........................"""
		params = [
				[ [10.0, 20.0, 30.0], [240.0, 40.0, 23.0], [150.0, 7.0, 80.0], [34.0, 70.0, 144.0] ] ,
				[ [150.0, 70.0, 35.0], [230.0, 45.0, 25.0], [160.0, 9.0, 75.0], [32.0, 68.0, 147.0] ]
				]
		
		matrix_rot = self.calculate_matrix_rot(params)
		results = self.internal_calculate_avg_diff_Util(params, matrix_rot)
		control = self.internal_calculate_avg_diff(params, matrix_rot)
		self.assertEqual(len(results), len(control))
		for i in xrange(len(results)):
			self.assertAlmostEquals( results[i], control[i], delta=1.0 )

	def test_four_configurations(self):
		"""test_four_configurations........................"""
		params = [
				[ [ 10.0, 20.0, 30.0], [240.0, 40.0, 23.0], [150.0, 7.0, 80.0], [34.0, 70.0, 144.0] ],
				[ [120.0, 70.0, 35.0], [230.0, 45.0, 25.0], [160.0, 9.0, 75.0], [32.0, 68.0, 147.0] ],
				[ [ 10.0, 26.0, 32.0], [ 40.0, 70.0, 23.0], [ 50.0, 7.0, 80.0], [14.0,  7.0, 144.0] ],
				[ [ 50.0, 75.0, 35.0], [130.0, 25.0, 25.0], [ 60.0, 0.0, 75.0], [39.0,  8.0, 147.0] ]
				]
		
		matrix_rot = self.calculate_matrix_rot(params)
		results = self.internal_calculate_avg_diff_Util(params, matrix_rot)
		control = self.internal_calculate_avg_diff(params, matrix_rot)
		self.assertEqual(len(results), len(control))
		for i in xrange(len(results)):
			self.assertAlmostEquals( results[i], control[i], delta=1.0 )


def test_main():
	from EMAN2 import Log
	p = OptionParser()
	p.add_option('--t', action='store_true', help='test exception', default=False )
	global IS_TEST_EXCEPTION
	opt, args = p.parse_args()
	if opt.t:
		IS_TEST_EXCEPTION = True
	Log.logger().set_level(-1)  #perfect solution for quenching the Log error information, thank Liwei
	suite = unittest.TestLoader().loadTestsFromTestCase(Test_Util_diff_between_matrix_of_3D_parameters_angles)
	unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
	test_main()

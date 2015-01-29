

def results_from_tests():

	def comparison_between_results_from_old_code_and_new_code():
		pass#
		pass#

		new code location:
			/home/hvoicu/Analysis/test_rviper
		def Comments_from_Pawel_about_the_two_sets_of_results():
			pass#
			pass#

			"""
			Conclusions from your runs.

			The result of the two runs, the first main iteration, are very similar, this is prior to elimination.  After elimination they get progressively less similar.

			It can mean one of two things:

			- too many images are eliminated, meaning percentile is too high

			- the three structures from one run of the program are more similar to each other than to three structure from the other run f the program.
			These two runs are totally independent.  It may mean that there is too mach cross-talk between run000-2.

			We will get more class averages and retest the overall strategy.
			However, please write it down somewhere so we have a starting point for discussions.

			Pawel.

			"""

def current_action_items():

	def Fix_main_iteration_starts_with_zero():
		pass#
		pass#
		To avoid checking the special condition of the first iteration, create main000 and put all relevant files there

	def Rules_for_importing_the_data():
		pass#
		pass#

		#  1.0  Looking at test600 structure it is unclear to me what the data structure is.
		There is the original g20.hdf, which then turns into bdb in test600, and a copy of it existis in master.
		The rule should be: if the input is hdf, you copy it into master as a bdb.  If the input is bdb, you
		create a copy as a virtual stack in master.  So the data part is only duplicated if input is hdf.

	def Find_a_better_fitness_function():
		pass#
		pass#

		L2 norm of reference volumes

	def Create_log_file_for_r_viper():
		pass#
		pass#

	def Create_documentation_on_web_page():
		pass#
		pass#


def finished_action_items():

	def Histogram_and_outliers_calculate_in_parallel():
		pass#
		pass#

		calculation of histogram takes time, it should be done in parallel with finding outliers

		Done, remains to be tested

	def Reported_Indices_refer_to_the_original_index_in_the_stack():
		pass#
		pass#
		Reporting of the images that are excluded from the stack must be done using the original index
		Solution: attach a flag to the original vstack that keeps track of changes

	def Is_rotation_done_properly():
		pass#
		pass#

	def Volume_reconstruction():
		pass#
		pass#
		use code from grandref

	def Make_projections_faster():
		pass#
		pass#

		def Write_a_transfer_function_for_refrings_that_concatenates_np_arrays():
			pass#
			pass#
			Done already: test_ring_transfer_function.py
			Needs to be inserted in the code

		def email_001():
			pass#
			pass#

			"""
			Hi,

			comments;

			- the program is waaaaaaaaay too slow.  We will not get anything this way.  It is too slow even for routine use, but it is
			impossible to test it if it takes so long to run.

			- please move to stampede immediately.  I do not recall what is the number of CPUs one can take there,
			but do not use larger than the number permitted by the development queue, even though you obviously must
			use normal queue.  In my observations requesting more doe not help; to the contrary, the waiting time
			increases.  In practice it means 256 or something like that.  You can quickly run couple of simple tests
			in the development queue to see how it performs.

			- while the times reported in the log file are not entirely reliable, it is quite obvious that generation of reprojections
			is the bottleneck.  The time is constant and is ~48s, while the matching (the other slow step) is initially very fast (0.05),
			later increases (see below).

			- In light of the above, we have to find a way to make projections faster, which in practice means figuring how to broadcast
			them safely.  I am afraid it would call for writing dedicated code and possibly making modifications in the wrapper.
			The idea would be to put generated projections in a longer consecutive array.  Here are my suggestions:

			1.  Ideally, one would put all refrings (reprojections converted to rings kept in one 1D array) in a very long 1D array
			and broadcast it to all nodes.  One would also one to skip headers, which I think is done now.  Still, for unknown reason
			the mpi wrapper converts everything to numpy array, so there is yet another level of copying.  I do not know why it is necessary.
			At the bottom, wrapper uses a C MPI function, maybe this is why the conversion is necessary.  However, I hope you’ll
			agree this is all awkward, as at the origin the data is a 32-float.  Something to look into, but not right away.

			The refrings are supposed to to have header information, it can be however added there after the transfer, as I believe is done in current code.

			The main problem is that there is a buffer limit in MPI, i.e., there are only that many bytes that MPI can transfer in one shot.
			The limitation might be system dependent.  I fought with it in the context of volume transfer (function reduce_EMData_to_root).
			However, maybe it is only for reduce-type function.  Please check.

			In short, please try to write a function that will:
			1. put all 1D refrings on a given node in a one long 1D array
			2. broadcasts number of refrings in it
			3. broadcasts the long array
			4 unpacks it.

			This would be of course memory-inefficient.  One should do putting on one long 1D array DURING generation of refrings, but thi
			can be done later.

			I believe it would solve the problem and significantly accelerates the code.
			"""

	def debug_restart():
		pass#
		pass#
		program picks up from where it stopped
		test branching




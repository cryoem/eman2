# import os
# from datetime import datetime
#
# """ the 'test_user_functions.py' is not tested because, at the 30/07/2019, there is nothing to test"""
# FILE_LIST_NOSE_TEST = [
#     "test_morphology.py",
#     "test_projection.py",
#     "test_reconstruction.py",
#     "test_utilities.py",
#     "test_alignment.py",
#     "test_filter.py",
#     "test_fundamentals.py",
#     "test_statistics.py",
#     "test_pixel_error.py",
#     "test_multi_shc.py",
# ]
# CURRENT_PATH = os.path.dirname(os.path.abspath(__file__))
#
# # in order to be able to run it from pycharm too, we have to specify the whole path of the commands. change them with your own path or set to FALSE and run it from the commandLine
# PYCHARM_RUN = True
# PYTEST_PATH = "/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/" #"/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/"
# MPI_PATH = "/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/" #"/home/adnan/applications/sphire/v1.1/envs/conda_fresh/bin/"
#
#
# OUTPUTFILE = "ALL_TEST_LOG_FILE"
#
#
# # it will save the result on the log file. NB: nosetest logs on stderr (see: https://stackoverflow.com/questions/9519717/how-do-i-redirect-the-output-of-nosetests-to-a-textfile)
# def main():
#     date = str(datetime.now()).split(".")[0].replace(" ", "_").replace(":", "")
#     outputfile = OUTPUTFILE + "_" + date + ".txt"
#     pytest_path = "" if PYCHARM_RUN is False else PYTEST_PATH
#     mpi_path = "" if PYCHARM_RUN is False else MPI_PATH
#     print("running 'test_systemExit.py' ...")
#     os.system(
#         pytest_path
#         + "pytest -v "
#         + os.path.join(CURRENT_PATH, "test_systemExit.py")
#         + " >"
#         + os.path.join(CURRENT_PATH, outputfile)
#     )
#
#     for f in FILE_LIST_NOSE_TEST:
#         print("running '" + f + "' ...")
#         os.system(
#             mpi_path
#             + "mpirun -np 1 nosetests -v "
#             + os.path.join(CURRENT_PATH, f)
#             + " 2>>"
#             + os.path.join(CURRENT_PATH, outputfile)
#         )
#
#
# if __name__ == "__main__":
#     main()

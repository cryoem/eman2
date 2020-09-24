from setuptools import setup, find_packages
import re
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name='sphire',
    version=find_version("sphire", "__init__.py"),
    python_requires='>3.4.0',
    packages=find_packages(),
    url='---',
    license='---',
    author='SPHIRE Development Team',
    author_email='sphire-devel@mpi-dortmund.mpg.de',
    description='The SPHIRE cryo-EM package',
)
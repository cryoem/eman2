from setuptools import setup, find_packages
import re
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))
os.chdir(here)

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

def create_entrypoints():
    import glob
    import os
    all_files = glob.glob("sphire/bin/*.py")

    entries = []
    for path in all_files:
        filename = os.path.basename(path)
        with open(path) as file:
            #https://regex101.com/r/b8NlEJ/2
            result = re.findall("^def\s+(_?main_?)\(\)",file.read(),re.M)
        if len(result) == 1:
            main_name = result[0]
            entry = filename + " = sphire.bin." + os.path.splitext(filename)[0]+":"+main_name
            entries.append(entry)
        elif filename in {"__init__.py"}:
            pass
        else:
            assert False,filename + " Main method not found."


    return entries

def get_longdescription():
    from os import path
    this_directory = path.abspath(path.dirname(__file__))
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
    return long_description
setup(
    name='sphire',
    version=find_version("sphire", "__init__.py"),
    python_requires='>3.4.0',
    packages=find_packages(),
    url='https://github.com/cryoem/eman2',
    license='GPLv3',
    author='SPHIRE Development Team',
    author_email='sphire-devel@mpi-dortmund.mpg.de',
    description='Python binaries and libraries for the SPHIRE cryo-EM package',
    long_description=get_longdescription(),
    long_description_content_type="text/markdown",
    entry_points={
        'console_scripts': create_entrypoints()},
)

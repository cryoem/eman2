Home: http://blake.bcm.edu/emanwiki

License: GPL2, BSD 3-Clause

Summary: A scientific image processing software suite with a focus on CryoEM and CryoET.



Current build status
====================

<a href="https://circleci.com/gh/cryoem/eman2">
<img alt="Linux" src="https://img.shields.io/circleci/project/github/cryoem/eman2/master.svg?label=Linux">
</a>

<a href="https://travis-ci.org/cryoem/eman2">
<img alt="macOS" src="https://img.shields.io/travis/cryoem/eman2/master.svg?label=macOS">
</a>

<a href="https://ci.appveyor.com/project/cryoem/eman2/branch/master">
<img alt="windows" src="https://img.shields.io/appveyor/ci/cryoem/eman2/master.svg?label=Windows">
</a>

Installation
====================
EMAN2 is now built within Anaconda, trying to compile using system dependencies (without Anaconda) is not
supported. Building from source with Anaconda takes just a few simple steps and is normally fairly foolproof.
Instructions at <a href="http://eman2.org">http://eman2.org</a>

Coding Style
====================
1) EMAN2 follows the GNU coding style with minor differences. We use
   GNU indent to make the proper indentation.
2) The naming styles are:
   1) All source code files use lower cases.
   2) All classes and data types use uppercase in the first letter.
   3) All functions use lower cases with '_' as the word separator.
   4) All parameter names for modular class use lower cases with '_' 
   as the word separator.

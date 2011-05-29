from distutils.core import setup
from distutils.extension import Extension
import os.path
import sys, glob, re, shutil

def main():

    include_dirs = re.findall( "-I\s*(\S+)", "-I/usr/include -I/home/andreas/test_install/include -I/usr/include/mysql  -g -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -fno-strict-aliasing -fwrapv -fPIC   -DUNIV_LINUX" ) +\
        [".","../pairsdblib", ".." ]

    libraries=["boost_python-mt", "pairsdblib"]
    library_dirs=[ "../pairsdblib/.libs" ]

    files = glob.glob( "modules/*.cpp" )

    setup(name="pairsdblib",    
          version='1.0',
          description="pairsdb utility functions in python",
          author='Andreas Heger',
          author_email='andreas.heger@gmail.com',
          ext_modules=[
                       Extension("pairsdblib",
                             files,
                             library_dirs=library_dirs,
                             libraries=libraries,
                             include_dirs=include_dirs,
                             depends=[]),
                             ],
          data_files = [ ('include/pairsdblib', ('exposed_decl.pypp.txt',)) ],
          )

if __name__ == "__main__":
    sys.exit(main())

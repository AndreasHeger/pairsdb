
USAGE = """python setup.py [options] command [command [...]]

The following commands are available:

* build: build and compile python extension to pairsdblib
* test: run some tests
* install: install the extension
"""

import re, sys, os, optparse, subprocess

from types import *

try:
    from pyplusplus import module_builder, messages, function_transformers
    from pyplusplus.module_builder import call_policies
    from pyplusplus.decl_wrappers import \
    return_value_policy, manage_new_object, copy_const_reference, reference_existing_object, return_self, return_arg

    from pygccxml import declarations
    GLOBAL_HAS_PYPLUSPLUS = True
except ImportError:
    GLOBAL_HAS_PYPLUSPLUS = False

import distutils.sysconfig
import os.path
import shutil

def findBoost( options ):
    """find location of boost."""
    if not options.boost_dir:
        if "BOOST_ROOT" in os.environ:
            options.boost_dir = os.environ['BOOST_ROOT']
            return 

        for x in ("/usr/local/boost", "/opt/boost" ):
            if os.path.exists( x ):
               options.boost_dir = x 
               return

        raise "could not find BOOST. Please specify location of BOOST as option or set BOOST_ROOT environment variable."

def checkRequisites( options ):
    """check if boost is present."""
    
    nerrors = 0
    if not os.path.exists(options.boost_dir):
        nerrors += 1
        print "could not find boost directory %s" % options.boost_dir
                    
    if not os.path.exists(options.pairsdblib_lib_dir ):                
        nerrors += 1
        print "could not find library directory %s" % options.pairsdblib_lib_dir
    elif not os.path.exists( options.pairsdblib_lib_dir + "/libpairsdblib.so" ):
        nerrors += 1
        print "could not find pairsdblib shared library %s/libpairsdblib.so" % options.pairsdblib_lib_dir
        
    return nerrors

def exportFunctions( mb ):
    """export utility functions."""

    ## include all free functions
    mb.namespace("pairsdblib").free_functions().include()
   
def exportHandles( mb ):
    """include handle classes. 
    
    These are shared_ptr<> typedefs.     
    """
    handles_to_export = []
    
    for handle in handles_to_export:
        
        # for c in mb.classes( lambda x: handle[1:]in x.name ):
        #     print "class=", c.name 
        # for c in mb.decls( lambda x: handle[1:] in x.name):
        #    print "decl=", c.name
            
        pointer = "boost::shared_ptr<pairsdblib::%s>" % handle[1:]
        try:
            d = mb.decl( pointer )
        except RuntimeError:
            print "could not find handle class %s, searching with %s" % (handle, pointer)
            continue
        
        # print "exporting %s" % handle
        
        d.include()
        # d.rename( handle )
        # d.alias = handle

def buildModule( include_paths, dest, options) :
    """build module using py++."""
    
    if not GLOBAL_HAS_PYPLUSPLUS:
        raise "can not build the interface, as py++ is not installed."
    
    if options.force:
        if os.path.exists( "cache" ):
            os.remove( "cache" )
            
    # add the include paths explicitely. This is a patch, as 
    # gccxml 0.9 was not able to pickup the directories when
    # using include_paths
    cflags = '--gccxml-cxxflags "%s"' % " ".join( [ "-I %s" % x for x in include_paths ] )
            
    # creating an instance of class that will help you to expose your declarations
    mb = module_builder.module_builder_t( [r"includes.h"]
                                          , gccxml_path=r""
                                          , cflags = options.gccxml_options + " " + cflags
                                          , cache="cache"
                                          , start_with_declarations=( "pairsdblib","py_details" )
                                          , working_directory=r"."
                                          , include_paths=include_paths
                                          , define_symbols=[]
                                          , )
    
    ## exclude py_details namespace, only used to instantiate template classes
    mb.namespace( 'py_details' ).exclude()

    if options.verbose:
        print "# declarations before building interface."
        mb.print_declarations()
    
    exportFunctions( mb )
        
    ## Every declaration will be exposed at its own line
    mb.classes().always_expose_using_scope = True

    my_exception = mb.class_( 'PairsdblibException' )
    my_exception.translate_exception_to_string( 'PyExc_RuntimeError', 'exc.what()')
    
    ######################################################################
    #Well, don't you want to see what is going on?
    if options.verbose:
        print "# declarations after building interface."        
        mb.print_declarations()
    
    # creating code creator. After this step you should not modify/customize declarations.
    mb.build_code_creator( module_name='pairsdblib' )
    mb.code_creator.add_include( "iostream" )
    mb.code_creator.add_include( "cstdio" )
    mb.split_module( "modules" )

    # check for alignlib installation
    alignlib_decl = "exposed_decl.pypp.txt"
    for p in sys.path:
        decls = os.path.join( p, "alignlib/exposed_decl.pypp.txt")
        if os.path.exists( decls ):
            mb.register_module_dependency( decls )
            break
    else:
        raise IOError( "could not find file %s with alignlib declarations" % alignlib_decl )
  
    # writing code to file.
    mb.write_module( dest )


if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "-f", "--force", dest="force", action="store_true",
                      help="force complete rebuilt.. [%default]")
    
    parser.add_option( "--boost-dir", dest="boost_dir", type="string",
                       help="location of boost [%default]." )
    
    parser.add_option( "--verbose", dest="verbose", action="store_true",
                       help="output details [%default]." )
    
    parser.add_option( "--gccxml-options", dest="gccxml_options", type="string",
                        help="flags to be passed to gccxml [%default]." )

    parser.add_option( "--cflags", dest="cflags", type="string",
                       help="flags to be passed on [%default]." )

    parser.set_defaults( extension_name = "pairsdblib",
                         force = False, 
                         src_dir = "../pairsdblib",
                         boost_dir = None,
                         compiler = None,
                         cflags = "-I/home/andreas/test/include -I/usr/include/mysql  -g -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -fno-strict-aliasing -fwrapv -fPIC",
                         gccxml_options = "",
                         pairsdblib_lib_dir = "../pairsdblib/.libs",
                         pairsdblib_include_dir = "../pairsdblib",
                         build_dir = ".",
                         verbose = False,
                         )
    
    (options, args) = parser.parse_args()

    if len(args) < 1:
        print USAGE
        raise "please supply a command"

    ## find boost
    findBoost( options )

    commands = map( lambda x: x.lower(), args)
    
    for command in commands:
        if command not in ("build", "test", "install", "generate-interface", "compile-interface" ):
            print USAGE
            raise "unknown command %s" % command
        
    for command in commands:

        if command in ("build", "generate-interface" ):
    
            nerrors = checkRequisites( options )
            
            if nerrors:
                print "found %i errors - aborting build." % (nerrors)
            
            ## installation directory of pairsdblib
            src_dir=os.path.abspath( options.src_dir )
            
            module_name = "%s.cpp" % options.extension_name
            
            if options.force or not os.path.exists( module_name):
                print "building module %s" % module_name     
                include_paths = [src_dir]
                dd = re.findall( "-I\s*(\S+)", options.cflags )
                include_paths.extend( dd )
                # patch for alignlib layout
                for d in dd:
                    p = os.path.join( d, "alignlib" )
                    if os.path.exists( p ):
                        include_paths.append( p )
                buildModule( include_paths = include_paths, dest = module_name, options = options )
                
            if command == "generate-interface": break                
        
            if command == "install":
            
                python_lib = distutils.sysconfig.get_python_lib()
                python_lib_data = python_lib + "/pairsdblib"
                results = []
                for root, dirs, files in os.walk('.'):
        
                    if "pairsdblib.so" in files:
                        results.append( os.path.join(root, "pairsdblib.so") )
                        
                if len(results) == 0:
                    print "could not find pairsdblib.so"
                    print "please run setup.py build first."
                    sys.exit(0)
                
                if len(results) > 1:
                    print "found more than one pairsdblib.so in %s." % str(results)
                    print "confused and thus not installing."
            
                print "installing %s in %s" % (results[0], python_lib )
                try:
                    shutil.copy( results[0], python_lib )
                    if not os.path.exists( python_lib_data ):  
                        os.mkdir( python_lib_data ) 
                    shutil.copy( "exposed_decl.pypp.txt", python_lib_data )
                except IOError, msg:
                    print "installation failed: %s" % str(msg)
                    
        elif command == "test":
            pass

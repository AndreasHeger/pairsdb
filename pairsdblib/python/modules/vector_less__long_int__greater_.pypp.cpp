// This file has been generated by Py++.

#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "includes.h"
#include "iostream"
#include "cstdio"
#include "vector_less__long_int__greater_.pypp.hpp"

namespace bp = boost::python;

void register_vector_less__long_int__greater__class(){

    bp::class_< std::vector< long int > >("vector_less__long_int__greater_")    
        .def( bp::vector_indexing_suite< ::std::vector< long int >, true >() );

}
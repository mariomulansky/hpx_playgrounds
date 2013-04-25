// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_OPERATIONS_GLOBAL_HPP
#define DATAFLOW_SHARED_OPERATIONS_GLOBAL_HPP

#include <vector>
#include <memory>

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/runtime/actions/component_action.hpp>

#include "serialization.hpp"

typedef std::vector<double> dvec;

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vector;

shared_vector operation3_global( shared_vector x1 , const shared_vector &x2 , 
                                 const shared_vector &x3 , 
                                 const double a1 , const double a2 )
{
    for( size_t i=0 ; i<x1->size() ; ++i )
        (*x1)[i] = a1*(*x2)[i] + a2*(*x3)[i];
    return x1;
}

HPX_PLAIN_ACTION( operation3_global , operation3_global_action );

#endif

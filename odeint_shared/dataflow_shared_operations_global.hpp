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

template< typename Operation >
shared_vector operation3( shared_vector x1 , const shared_vector &x2 , 
                          const shared_vector &x3 , Operation op )
{
    for( size_t i=0 ; i<x1->size() ; ++i )
        op( (*x1)[i] , (*x2)[i] , (*x3)[i] );
    return x1;
}


template< typename Operation >
struct operation3_action
    : hpx::actions::make_action<
    shared_vector (*)( shared_vector , 
                       const shared_vector& , 
                       const shared_vector& , 
                       Operation op ) , 
    &operation3<Operation>, 
    operation3_action<Operation> >
{};


HPX_REGISTER_PLAIN_ACTION_TEMPLATE((template < typename Operation>),(operation3_action< Operation>))


//HPX_PLAIN_ACTION( operation3_global , operation3_global_action );

#endif

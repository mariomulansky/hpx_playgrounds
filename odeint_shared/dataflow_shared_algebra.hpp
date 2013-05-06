// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_ALGEBRA_HPP
#define DATAFLOW_SHARED_ALGEBRA_HPP

#include "dataflow_shared_operations.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::lcos::future;
using hpx::find_here;

struct dataflow_shared_algebra
{

    //template< class S1 , class S2 , class S3 , class Op >
    // for now just a single state  type
    template< typename S , typename Op >
    void for_each3( S &s1 , const S &s2 , const S &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        for( size_t i=0 ; i<N ; ++i )
        {
            typedef std::vector< double > dvec;
            typedef std::shared_ptr< dvec > shared_vec;
            typedef hpx_odeint_actions::operation3_action<typename S::value_type::result_type,Op> action;
            s1[i] = dataflow< action >( find_here() ,
                                        s1[i] , s2[i] , s3[i] , 
                                        op );
        }
    }
};

#endif

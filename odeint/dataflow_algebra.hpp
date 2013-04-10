// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_ALGEBRA_HPP
#define DATAFLOW_ALGEBRA_HPP

#include "dataflow_operations.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

struct dataflow_algebra
{

    template< class S1 , class S2 , class S3 , class Op >
    static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        for( size_t i=0 ; i<N ; ++i )
        {
            s1[i] = dataflow< hpx_scale_sum2_action >( find_here() , s2[i] , s3[i] , 
                                                       op.m_alpha1 , op.m_alpha2 );
        }
    }

};


#endif

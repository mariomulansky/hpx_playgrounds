// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_ALGEBRA_HPP
#define DATAFLOW_ALGEBRA_HPP

#include "dataflow_operations.hpp"

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::find_here;

// struct dataflow_algebra
// {

//     template< class S1 , class S2 , class S3 , class Op >
//     static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
//     {
//         const size_t N = boost::size( s1 );
//         for( size_t i=0 ; i<N ; ++i )
//         {
//             s1[i] = dataflow< hpx_scale_sum2_action >( find_here() , s2[i] , s3[i] , 
//                                                        op.m_alpha1 , op.m_alpha2 );
//         }
//     }

// };

struct dataflow_sub_algebra
{
    template< class S1 , class S2 , class S3 , class Op >
    static void for_each3( S1 &s1 , const S2 &s2 , const S3 &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        for( size_t i=0 ; i<N ; ++i )
        {
            hpx::id_type target_locality = hpx::find_here();
            hpx::id_type o = hpx::components::new_< dataflow_actions::operation3 >( target_locality ).get();

            s1[i] = dataflow< typename dataflow_actions::operation3::operation3_action< std::vector<double> , Op> >( o , s1[i] , s2[i] , s3[i] , op );
        }
    }
};
/*
struct dataflow_ref_algebra
{
    template< class S1 , class S2 , class S3 , class Op >
    static void for_each3( S1 &s1 , S2 &s2 , S3 &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        for( size_t i=0 ; i<N ; ++i )
        {
            s1[i] = dataflow< hpx_scale_sum2_ref_action >( find_here() , s2[i] , s3[i] ,
                                                           op.m_alpha1 , op.m_alpha2 );
        }
    }
};
*/
#endif

// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_ALGEBRA_HPP
#define DATAFLOW_SHARED_ALGEBRA_HPP

#include <hpx/lcos/local/dataflow.hpp>
#include <hpx/util/unwrap.hpp>

using hpx::lcos::local::dataflow;
using hpx::lcos::future;
using hpx::util::unwrap;

struct local_dataflow_algebra
{

    //template< class S1 , class S2 , class S3 , class Op >
    // for now just a single state  type
    template< typename S , typename Op >
    void for_each3( S &s1 , const S &s2 , const S &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        for( size_t i=0 ; i<N ; ++i )
        {
            s1[i] = dataflow( hpx::launch::sync , unwrap(op) , s1[i] , s2[i] , s3[i] );
        }
    }
};

#endif

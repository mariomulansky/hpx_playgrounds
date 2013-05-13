// Copyright 2013 Mario Mulansky
#ifndef DATAFLOW_SHARED_ALGEBRA_HPP
#define DATAFLOW_SHARED_ALGEBRA_HPP

#include <hpx/include/iostreams.hpp>

using hpx::lcos::dataflow;
using hpx::lcos::dataflow_base;
using hpx::lcos::future;
using hpx::find_here;

namespace hpx_odeint_actions {

    template< typename S , typename Operation >
    S operation2d_3( S x1 , const S &x2 , const S &x3 , 
                     Operation op )
    {
        //hpx::cout << (boost::format("sizes: %d , %d , %d\n") % (x1->size()) % (x2->size()) % (x3->size()) ) << hpx::flush;
        for( size_t i=0 ; i<x1->size() ; ++i )
            for( size_t j=0 ; j<(*x1)[i].size() ; ++j )
                op( (*x1)[i][j] , (*x2)[i][j] , (*x3)[i][j] );
        return x1;
    }


    template< typename S , typename Operation >
    struct operation2d_3_action
        : hpx::actions::make_action<
        S (*)( S , 
               const S& , 
               const S& , 
               Operation op ) , 
        &operation2d_3<S,Operation>, 
        operation2d_3_action<S,Operation> >
    //        boost::mpl::true_ >
    {};

}

HPX_REGISTER_PLAIN_ACTION_TEMPLATE(
    (template< typename S , typename Operation >),
    (hpx_odeint_actions::operation2d_3_action< S , Operation >))


struct dataflow_shared_algebra_2d
{

    //template< class S1 , class S2 , class S3 , class Op >
    // for now just a single state  type
    template< typename S , typename Op >
    void for_each3( S &s1 , const S &s2 , const S &s3 , Op op )
    {
        const size_t N = boost::size( s1 );
        // std::cout << "dataflow sizes: " << boost::size( s1 ) << " , ";
        // std::cout << boost::size( s2 ) << " , " << boost::size( s3 ) << std::endl;
        // std::cout << s1[0].get_future().get()->size() << " , " << s2[0].get_future().get()->size() << " , " << s3[0].get_future().get()->size() << std::endl;
        // std::cout << std::flush;
        for( size_t i=0 ; i<N ; ++i )
        {
            typedef std::vector< double > dvec;
            typedef std::shared_ptr< dvec > shared_vec;
            typedef hpx_odeint_actions::operation2d_3_action< typename S::value_type::result_type,Op> action;
            s1[i] = dataflow< action >( find_here() ,
                                        s1[i] , s2[i] , s3[i] , 
                                        op );
        }
    }
};

#endif

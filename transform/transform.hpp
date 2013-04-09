//  Copyright (c) 2013 Mario Mulansky
//
// based on the for_each implementation in hpx/for_each.hpp by Thomas Heller
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#include <iostream>
#include <vector>
#include <algorithm>

#include <hpx/config.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/async.hpp>
#include <hpx/util/move.hpp>
#include <hpx/for_each.hpp>

#include <boost/range.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/ref.hpp>


namespace hpx { 

    namespace detail {

        template< typename F , typename Value >
        struct binder2 {

            binder2( F f , Value &x1 , const Value &x2 )
                : m_f( f ) , m_x1( x1 ) , m_x2( x2 )
            { }

            void operator()( void )
            {
                m_f( m_x1 , m_x2 );
            }

            F m_f;
            Value &m_x1;
            const Value &m_x2;
        };


        template< typename F , typename InOutIterator , typename InIterator >
        struct partial_transform {

            partial_transform( InOutIterator inout_begin , InOutIterator inout_end , InIterator in_begin , F f )
                : m_inout_begin( inout_begin ) , m_inout_end( inout_end ) , 
                  m_in_begin( in_begin ) , m_f( f )
            { }

            void operator()( void )
            {
                //std::clog << "running for " << m_inout_end-m_inout_begin << " elements" << std::endl;
                while( m_inout_begin != m_inout_end )
                {
                    m_f( *m_inout_begin++ , *m_in_begin++ );
                }
            }

            InOutIterator m_inout_begin , m_inout_end;
            InIterator m_in_begin;
            F m_f;
        };
    }


    template< typename RangeInOut , typename RangeIn , typename F >
    std::vector< lcos::future<void> > 
    transform_inplace( RangeInOut &range_inout , const RangeIn &range_in , F f , 
               std::size_t granularity = 0 )
    {

        typedef typename boost::range_value<RangeInOut>::type value_type;
        typedef void result_type;
        typedef lcos::future<result_type> future_type;
        typedef std::vector<future_type> futures_type;
        typedef typename boost::range_iterator<RangeInOut>::type inout_iterator_type;
        typedef typename boost::range_iterator<const RangeIn>::type in_iterator_type;

        futures_type futures(boost::size(range_inout));

        // defualt value granularity = size / threadcount, that 
        if( granularity == 0 )
        {
            granularity = boost::size(range_inout) / hpx::get_os_thread_count();
        }

        if(futures.size() <= granularity)
        {
            inout_iterator_type inout_begin( boost::begin( range_inout ) );
            in_iterator_type in_begin( boost::begin( range_in ) );
            std::size_t i = 0;
            while( inout_begin != boost::end(range_inout) )
            {
                
                //futures[i] = hpx::async(HPX_STD_BIND(hpx::util::protect(f), boost::ref(*inout_begin) , *in_begin));
                futures[i] = hpx::async( detail::binder2< hpx::util::detail::protected_bind<F> , value_type >( hpx::util::protect(f) , *inout_begin , *in_begin ) );

                // we have to use our own binder here to ensure const correctness of the arguments
                // what does protect do? do we need that?
                //futures[i] = hpx::async( detail::binder2< F , value_type >( f , *inout_begin , *in_begin ) );
                inout_begin++;
                in_begin++;
                ++i;
            }
            return futures;
        }

        inout_iterator_type inout_begin = boost::begin(range_inout);
        inout_iterator_type inout_mid = inout_begin + boost::size(range_inout)/2;
        inout_iterator_type inout_end = boost::end(range_inout);

        in_iterator_type in_begin = boost::begin(range_in);
        in_iterator_type in_mid = in_begin + boost::size(range_in)/2;
        in_iterator_type in_end = boost::end(range_in);

        futures_type (*transform_impl)( boost::iterator_range<inout_iterator_type> &, 
                                        const boost::iterator_range<in_iterator_type> & , 
                                        F, const std::size_t ) = transform_inplace;

        lcos::future<futures_type>
            left_future(
                hpx::async(
                    HPX_STD_BIND(
                        transform_impl
                      , boost::make_iterator_range(inout_begin, inout_mid)
                      , boost::make_iterator_range(in_begin, in_mid)
                      , f
                      , granularity
                    )
                )
            );

        lcos::future<futures_type>
            right_future(
                hpx::async(
                    HPX_STD_BIND(
                        transform_impl
                      , boost::make_iterator_range(inout_mid, inout_end)
                      , boost::make_iterator_range(in_mid, in_end)
                      , f 
                      , granularity
                    )
                )
            );

        typedef typename futures_type::iterator futures_iterator_type;

        futures_iterator_type futures_begin = futures.begin();
        futures_iterator_type futures_mid = futures.begin() + futures.size()/2;
        
        std::vector<hpx::lcos::future<void> > v;
        v.reserve(2);
        v.push_back(left_future.then(detail::for_each_copy_future<futures_iterator_type>(futures_begin)));
        v.push_back(right_future.then(detail::for_each_copy_future<futures_iterator_type>(futures_mid)));
        hpx::lcos::wait(v);

        return futures;
    }


    /* flat transform algorithm
     * divides the range into parts of size granularity and runs the 
     * transformation over each part in parallel.
     * granularity = 0 (default) let's the algorithm guess a chunksize 
     * based on size of the range and the thread count.
     */
    template< typename RangeInOut , typename RangeIn , typename F >
    std::vector< lcos::future<void> > 
    transform_inplace_flat( RangeInOut &range_inout , const RangeIn &range_in , F f , 
                    std::size_t granularity = 0 )
    {

        typedef typename boost::range_value<RangeInOut>::type value_type;
        typedef void result_type;
        typedef lcos::future<result_type> future_type;
        typedef std::vector<future_type> futures_type;
        typedef typename boost::range_iterator<RangeInOut>::type inout_iterator_type;
        typedef typename boost::range_iterator<const RangeIn>::type in_iterator_type;

        // default value granularity = size / threadcount
        if( granularity == 0 )
        {
            granularity = boost::size(range_inout) / hpx::get_os_thread_count();
        }

        // number of chunks the data is divided into
        std::size_t chunks = boost::size(range_inout)/granularity;
        if( boost::size(range_inout) % granularity != 0 )
        {
            // we have some remainder
            chunks++;
        }

        futures_type futures( chunks );

        inout_iterator_type inout_begin( boost::begin( range_inout ) );
        in_iterator_type in_begin( boost::begin( range_in ) );
        std::size_t n = 0;
        while( n<futures.size()-1 )
        {
            futures[n] = hpx::async( detail::partial_transform< F , inout_iterator_type , in_iterator_type >( inout_begin , inout_begin+granularity , in_begin , f ) );
            inout_begin += granularity;
            in_begin += granularity;
            n++;
        }
        futures[n] = hpx::async( detail::partial_transform< F , inout_iterator_type , in_iterator_type >( inout_begin , boost::end( range_inout ) , in_begin , f ) );

        return futures;
    }

} // namespace hpx

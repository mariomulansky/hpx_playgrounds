#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include <memory>

typedef std::vector< double > dvec;
typedef std::vector< dvec > dvecvec;
typedef std::shared_ptr< dvecvec > shared_vec;

struct initialize_zero
{
    const size_t m_N1;
    const size_t m_N2;

    initialize_zero( const size_t N1 , const size_t N2 )
        : m_N1( N1 ) , m_N2( N2 )
    { }

    shared_vec operator()( shared_vec v )
    {
        v->resize( m_N1 );
        for( size_t n=0 ; n<m_N1 ; ++n )
        {
            (*v)[n].resize( m_N2 );
            std::fill( (*v)[n].begin() , (*v)[n].end() , 0.0 );
        }
        return v;
    }
};


struct initialize_copy
{
    const dvecvec &m_data;
    const size_t m_index;
    const size_t m_len;

    initialize_copy( const dvecvec data , const size_t index , const size_t len )
        : m_data( data ) , m_index( index ) , m_len( len )
    { }

    shared_vec operator()( shared_vec v )
    {
        v->resize( m_len );
        for( size_t n=0 ; n<m_len ; ++n )
        {
            (*v)[n].resize( m_data[m_index+n].size() );
            std::copy( m_data[m_index+n].begin() , m_data[m_index+n].end() , (*v)[n].begin() );
        }
        return v;
    }
};

#endif

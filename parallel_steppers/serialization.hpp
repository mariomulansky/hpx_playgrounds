// Copyright 2013 Mario Mulansky

#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <array>
#include <memory>

#include <boost/serialization/serialization.hpp>

// add serialization to shared_ptr
namespace boost { namespace serialization {

template<class Archive , typename T , size_t N>
void serialize( Archive & ar, std::shared_ptr< std::array<T,N> > &v , 
                const unsigned int version )
{
    // empty for now
}

template<class Archive , typename T , size_t N>
void serialize( Archive & ar , std::array<T,N> &v , const unsigned int version )
{
    // empty for now
}

// template<class Archive>
// void serialize(Archive & ar, typename dvecvec::iterator iter , 
//                const unsigned int version)
// {
//     // empty for now
// }

} }

#endif

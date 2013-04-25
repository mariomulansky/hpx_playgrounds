// Copyright 2013 Mario Mulansky

#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <vector>
#include <memory>

#include <boost/serialization/serialization.hpp>

typedef std::vector< double > dvec;
typedef std::shared_ptr< dvec > shared_vec;


// add serialization to shared_ptr
namespace boost { namespace serialization {
template<class Archive>
void serialize(Archive & ar, shared_vec &v , const unsigned int version)
{
    // empty for now
}
} }

#endif

#include <boost/python.hpp>
#include "DmpAlgBgoRawTrack.h"

BOOST_PYTHON_MODULE(libDmpAlgBgoRawTrack){
  using namespace boost::python;

  class_<DmpAlgBgoRawTrack,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoRawTrack",init<>());
}

#ifndef TYPES_H
#define TYPES_H
#include <memory>
namespace jmat {
    typedef unsigned char uchar;
    typedef unsigned int  uint;
    typedef std::shared_ptr<uchar> shr_ptr;
    enum    DTYP {UCHAR=0 , INT=1 , FLOAT=2, DOUBLE=3};
}
#endif // TYPES_H

#ifndef TYPES_H
#define TYPES_H
#include <memory>
#include <cstdint>

namespace jmat {    
    typedef unsigned char  uchar;
    typedef uint8_t        uint8;
    typedef int8_t          int8;
    typedef uint32_t      uint32;
    typedef int32_t        int32;
    typedef uint64_t      uint64;
    typedef int64_t        int64;

    typedef std::shared_ptr<uchar> shr_ptr;
    enum    DTYP {UCHAR=0 , INT=1 , FLOAT=2, DOUBLE=3};
}
#endif // TYPES_H

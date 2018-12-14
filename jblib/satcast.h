#ifndef SATCAST_H
#define SATCAST_H
//
// bring similar code from opencv
//
#include "types.h"

namespace jmat {
    // template function prototypes
    template<typename T> static inline T sat_cast(int32  A){  return (T)A;   }
    template<typename T> static inline T sat_cast(float  A){  return (T)A;   }
    template<typename T> static inline T sat_cast(double A){  return (T)A;   }

    // template specializers
    template<> inline uchar sat_cast<uchar>(int32  A){ return ( (A >> 8) & 0x80000000) ? 0 : (A >> 8) ? UINT8_MAX : (uchar)A;   }
    template<> inline uchar sat_cast<uchar>(float  A){ return ( A < 0 ) ? 0 : (A > UINT8_MAX) ? UINT8_MAX : (uchar)A; }
    template<> inline uchar sat_cast<uchar>(double A){ return ( A < 0 ) ? 0 : (A > UINT8_MAX) ? UINT8_MAX : (uchar)A; }

    template<> inline int32 sat_cast<int32 >(float  A){ return ( A < INT32_MIN ) ? INT32_MIN : (A > INT32_MAX) ? INT32_MAX : (int32)A; }
    template<> inline int32 sat_cast<int32 >(double A){ return ( A < INT32_MIN ) ? INT32_MIN : (A > INT32_MAX) ? INT32_MAX : (int32)A; }
}

#endif // SATCAST_H

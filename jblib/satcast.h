/*
 * Copyright (C) 2020. Jong B. Choi
 * License : MIT License
 * contact : neoflame99@naver.com
 */

#ifndef SATCAST_H
#define SATCAST_H
//
// bring similar code from opencv
//
#include "types.h"

namespace jmat {
// template function prototypes
template<typename T> inline T sat_cast(int32  A){ return static_cast<T>(A); }
template<typename T> inline T sat_cast(float  A){ return static_cast<T>(A); }
template<typename T> inline T sat_cast(double A){ return static_cast<T>(A); }
template<typename T> inline T sat_cast(cmplx  A){ return static_cast<T>(A); }

// template specializers
template<> inline uchar sat_cast<uchar>(int32  A){ return ( (A >> 8) & 0x80000000) ? 0 : (A >> 8) ? UINT8_MAX : static_cast<uchar>(A); }
template<> inline uchar sat_cast<uchar>(float  A){ return ( A < 0 ) ? 0 : (A > UINT8_MAX) ? UINT8_MAX : static_cast<uchar>(A); }
template<> inline uchar sat_cast<uchar>(double A){ return ( A < 0 ) ? 0 : (A > UINT8_MAX) ? UINT8_MAX : static_cast<uchar>(A); }
template<> inline uchar sat_cast<uchar>(cmplx  A){ return ( A.re < 0 ) ? 0 : (A.re > UINT8_MAX) ? UINT8_MAX : static_cast<uchar>(A.re); }

template<> inline int32 sat_cast<int32>(float  A){ return ( A < INT32_MIN ) ? INT32_MIN : (A > INT32_MAX) ? INT32_MAX : static_cast<int32>(A); }
template<> inline int32 sat_cast<int32>(double A){ return ( A < INT32_MIN ) ? INT32_MIN : (A > INT32_MAX) ? INT32_MAX : static_cast<int32>(A); }
template<> inline int32 sat_cast<int32>(cmplx  A){ return static_cast<int32>(A.re); }

template<> inline float  sat_cast<float >(cmplx  A){ return static_cast<float >(A.re); }
template<> inline double sat_cast<double>(cmplx  A){ return A.re; }

}

#endif // SATCAST_H

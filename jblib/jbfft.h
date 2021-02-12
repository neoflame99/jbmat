#ifndef FFT_H
#define FFT_H

#include "jbMat.h"
#include "jbmath.h"
#include "satcast.h"
#ifdef _WIN_
    #define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <vector>

#define FFT_EXP_TABLE
#define MALLOC_F

namespace jmat{

    // fft
    void fft_dit2(_complex* dat, int32 len, bool backward=false);
    void fft_czt( _complex *dat, int32 len, bool inverse=false);
    void fft(_complex* dat, int32 len);
    void ifft(_complex* dat, int32 len);
    void fft2d(_complex* dat, int32 r_len, int32 c_len);
    void ifft2d(_complex* dat, int32 r_len, int32 c_len);
    inline void bitrev_permute(_complex* dat, int32 len);

    void fft_dif4(_complex *dat, int32 len, bool backward=false);
    void permute_radix4(_complex *a, int32 len);
    inline int32 digit4_rev(int x, int32 ldn, int32 radbit, int32 andBit);
    void fft_compositN(_complex *dat, int32 len, std::vector<int32>&fac, bool backward=false);
    void factorizeN(int32 N, std::vector<int32>& fac );

}
#endif // FFT_H

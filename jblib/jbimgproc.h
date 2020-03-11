#ifndef JBIMGPROC_H
#define JBIMGPROC_H
#include "jbMat.h"
#include "jbmath.h"
#include <math.h>
#include <vector>

#define FFT_EXP_TABLE
#define MALLOC_F

namespace jmat {
namespace imgproc {

    // color space conversion functions
    Mat rgb2ycc(const Mat& rgbIm, const int32 sel_eq = 0);
    Mat ycc2rgb(const Mat& yccIm, const int32 sel_eq = 0);
    Mat rgb2gray(const Mat& rgbIm, const int32 HowToGray = 0);
    Mat rgb2xyz(const Mat& rgbIm);
    Mat xyz2rgb(const Mat& xyzIm);
    Mat rgb2Yxy(const Mat& rgbIm);
    Mat Yxy2rgb(const Mat& YxyIm);
    template <typename _T> inline Mat _rgb2ycc(const Mat& rgbIm, const int32 sel_eq );
    template <typename _T> inline Mat _ycc2rgb(const Mat& yccIm, const int32 sel_eq );
    template <typename _T> inline Mat _rgb2gray(const Mat& rgbIm, const int32 HowToGray);
    template <typename _T> inline Mat _conv_rgb2xyz(const Mat& rgbIm);
    template <typename _T> inline Mat _conv_xyz2rgb(const Mat& xyzIm);
    template <typename _T> inline Mat _conv_rgb2Yxy(const Mat& rgbIm);
    template <typename _T> inline Mat _conv_Yxy2rgb(const Mat& YxyIm);

    // histogram functions
    Mat histoPmf(const Mat& src,const int32 bins, const double step);
    Mat histoCmf(const Mat& src,const int32 bins, const double step);
    Mat clip_HistoPmf(const Mat& src, const int32 clipVal,const int32 bins, const int32 step);
    Mat clip_HistoCmf(const Mat& src, const int32 clipVal,const int32 bins, const int32 step);
    Mat clip_HistoEqual(const Mat& src, const Mat& histCmf, const int32 step);
    template <typename _T> inline Mat _histoPmf(const Mat& src, const int32 bins, const double step);
    template <typename _T> inline Mat _histoCmf(const Mat& src, const int32 bins, const double step);
    template <typename _T> inline Mat _clip_HistoPmf(const Mat& src,const int32 clipVal, const int32 bins, const int32 step);
    template <typename _T> inline Mat _clip_HistoCmf(const Mat& src,const int32 clipVal, const int32 bins, const int32 step);
    template <typename _T> inline Mat _clip_HistoEqual(const Mat& src, const Mat& histCmf, const int32 step);

    // gamma function
    Mat gamma(const Mat& src, const double gmval);
    template <typename _T> inline Mat _gamma(const Mat& src, const double gmval);

    // tone mapping functions
    double inline nakaSigmoid (const double X, const double X0, const double Xmax );
    Mat gaussMaskGen (const double sigma, const double factor = 6, const uint32 ch=1);
    Mat boxMaskGen (const uint32 sz, const uint32 ch=1);
    Mat inline localMeanMat ( const Mat& src, const Mat& mask);
    Mat nakaSigTonemap( Mat& src, Mat& localmean, const double globalmean, const double Imax);
    Mat nakaSig3MeanTonemap( Mat& src, Mat& s_localmean, Mat& l_localmean, const double globalmean, const double Imax);
    Mat logRetinexTonemap( Mat& src, Mat& surround);
    template <typename _T> inline Mat _nakaSigTm(const Mat& Im, const Mat& localmean, const double globalmean, const double Imax);
    template <typename _T> inline Mat _nakaSig3MeanTm(const Mat& Im, const Mat& smLocalMean, const Mat& lgLocalMean, const double globalMean, const double Imax);
    template <typename _T> inline Mat _logRetinexTm(const Mat& Im, const Mat& localmean);

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

    // image resizing
    inline float cubic1d(float ma1, float a0, float a1, float a2, float t);
    Mat bicubicIntp(const Mat& src,const int32 sw, const int32 sh);
    Mat bilinearIntp(const Mat& src,const int32 sw, const int32 sh);
    Mat nearestIntp(const Mat& src,const int32 sw, const int32 sh);
    Mat decimate(const Mat& src,const int32 sw, const int32 sh);
    Mat copy_padding(const Mat& src, int32 pad_size=2);
    template <typename _T> inline Mat _bicubicIntp(const Mat& src, const int32 sw, const int32 sh);
    template <typename _T> inline Mat _nearestIntp(const Mat& src, const int32 sw, const int32 sh);
    template <typename _T> inline Mat _decimate(const Mat& src, const int32 sw, const int32 sh);

    // Image pyramid
    template <typename _T> inline Mat _gaussPyramid(const Mat& src, const int32 level);
    template <typename _T> inline Mat _laplPyramid( Mat& src, const int32 level);
    Mat gaussPyramid(const Mat& src, const int32 level=4);
    Mat laplPyramid( Mat& src, const int32 level=4);
}


namespace  imgproc{

    template <typename _T> Mat _rgb2ycc(const Mat& rgbIm, const int32 sel_eq ){
    /*
     *   if sel_eq = 0 (BT 601)
     *  Y = [  0.299    ,  0.587    ,  0.114    ]   [ r ]
     * Pb = [ -0.168736 , -0.331264 ,  0.5      ] * [ g ]
     * Pr = [  0.5      , -0.418688 , -0.081312 ]   [ b ]
     *
     *   if sel_eq = 1 (BT 709)
     *  Y = [  0.2126  ,  0.7152  ,  0.0722  ]   [ r ]
     * Pb = [ -0.11457 , -0.38543 ,	 0.5     ] * [ g ]
     * Pr = [  0.5     , -0.45415 ,	-0.04585 ]   [ b ]
     *
     */

        uint32 row = rgbIm.getRow();
        uint32 col = rgbIm.getCol();
        uint32 imsize = row * col;
        uint32 chsize = rgbIm.getChannel();

        Mat  A(rgbIm.getDatType(), row, col, chsize);
        _T* srcDat_pt = rgbIm.getDataPtr<_T>();
        _T* tarDat_pt = A.getDataPtr<_T>();

        uint32 ch_offset1 = imsize;
        uint32 ch_offset2 = imsize << 1;

        uint32 x,x2,x3;
        if( sel_eq == 0){
            for(x= 0 ; x < imsize; x++ ){
                x2 = x+ch_offset1;
                x3 = x+ch_offset2;
                tarDat_pt[x ] = bt601_r2y[0][0] * srcDat_pt[x  ] + bt601_r2y[0][1] * srcDat_pt[x2] + bt601_r2y[0][2] * srcDat_pt[x3];
                tarDat_pt[x2] = bt601_r2y[1][0] * srcDat_pt[x  ] + bt601_r2y[1][1] * srcDat_pt[x2] + bt601_r2y[1][2] * srcDat_pt[x3];
                tarDat_pt[x3] = bt601_r2y[2][0] * srcDat_pt[x  ] + bt601_r2y[2][1] * srcDat_pt[x2] + bt601_r2y[2][2] * srcDat_pt[x3];
            }
        }else if( sel_eq == 1){
            for(x= 0 ; x < imsize; x++ ){
                x2 = x+ch_offset1;
                x3 = x+ch_offset2;
                tarDat_pt[x ] = bt709_r2y[0][0] * srcDat_pt[x  ] + bt709_r2y[0][1] * srcDat_pt[x2] + bt709_r2y[0][2] * srcDat_pt[x3];
                tarDat_pt[x2] = bt709_r2y[1][0] * srcDat_pt[x  ] + bt709_r2y[1][1] * srcDat_pt[x2] + bt709_r2y[1][2] * srcDat_pt[x3];
                tarDat_pt[x3] = bt709_r2y[2][0] * srcDat_pt[x  ] + bt709_r2y[2][1] * srcDat_pt[x2] + bt709_r2y[2][2] * srcDat_pt[x3];
            }
        }
        return A;
    }
    template <typename _T> Mat _ycc2rgb(const Mat& yccIm, const int32 sel_eq ){
        /*
         *   if sel_eq = 0 (BT 601)
         *
         *  r = [ 1.0000 , -0.0000  ,  1.4020 ]   [ Y  ]
         *  g = [ 1.0000 , -0.3441  , -0.7141 ] * [ Pb ]
         *  b = [ 1.0000 ,  1.7720  ,  0.0000 ]   [ Pr ]
         *
         *   if sel_eq = 1 (BT 709)
         *  r = [ 1.0000 ,  0.0000  ,   1.5748  ]   [ Y  ]
         *  g = [ 1.0000 , -0.1873  ,  -0.4681  ] * [ Pb ]
         *  b = [ 1.0000 ,  1.8556  ,  -0.0000  ]   [ Pr ]
         *
         */

         uint32 row = yccIm.getRow();
         uint32 col = yccIm.getCol();
         uint32 imsize = row * col;
         uint32 chsize = yccIm.getChannel();

         Mat A(yccIm.getDatType(), row, col, chsize);
         _T *srcDat_pt = yccIm.getDataPtr<_T>();
         _T *tarDat_pt = A.getDataPtr<_T>();

         uint32 ch_offset1 = imsize ;
         uint32 ch_offset2 = imsize << 1;
         uint32 x,x2,x3;
         _T tmp1, tmp2, tmp3;
         if( sel_eq == 0){
             for(x= 0 ; x < imsize; x++){
                 x2 = x+ch_offset1;
                 x3 = x+ch_offset2;
                 tmp1 = bt601_y2r[0][0] * srcDat_pt[x  ] + bt601_y2r[0][1] * srcDat_pt[x2] + bt601_y2r[0][2] * srcDat_pt[x3];
                 tmp2 = bt601_y2r[1][0] * srcDat_pt[x  ] + bt601_y2r[1][1] * srcDat_pt[x2] + bt601_y2r[1][2] * srcDat_pt[x3];
                 tmp3 = bt601_y2r[2][0] * srcDat_pt[x  ] + bt601_y2r[2][1] * srcDat_pt[x2] + bt601_y2r[2][2] * srcDat_pt[x3];

                 tarDat_pt[x ] = (tmp1 < 0 ) ? 0 : tmp1;
                 tarDat_pt[x2] = (tmp2 < 0 ) ? 0 : tmp2;
                 tarDat_pt[x3] = (tmp3 < 0 ) ? 0 : tmp3;
             }
         }else if( sel_eq == 1){
             for(x= 0 ; x < imsize; x++){
                 x2 = x+ch_offset1;
                 x3 = x+ch_offset2;
                 tmp1 = bt709_y2r[0][0] * srcDat_pt[x  ] + bt709_y2r[0][1] * srcDat_pt[x2] + bt709_y2r[0][2] * srcDat_pt[x3];
                 tmp2 = bt709_y2r[1][0] * srcDat_pt[x  ] + bt709_y2r[1][1] * srcDat_pt[x2] + bt709_y2r[1][2] * srcDat_pt[x3];
                 tmp3 = bt709_y2r[2][0] * srcDat_pt[x  ] + bt709_y2r[2][1] * srcDat_pt[x2] + bt709_y2r[2][2] * srcDat_pt[x3];

                 tarDat_pt[x ] = (tmp1 < 0 ) ? 0 : tmp1;
                 tarDat_pt[x2] = (tmp2 < 0 ) ? 0 : tmp2;
                 tarDat_pt[x3] = (tmp3 < 0 ) ? 0 : tmp3;

             }
         }
         return A;
    }


    template <typename _T> Mat _rgb2gray(const Mat& rgbIm, const int32 HowToGray){
    /*
     *  if HowToGray = 0 (BT 601)
     *  Y =  0.299 * r   +  0.587 * g +  0.114 * b
     *
     *  if HowToGray = 1 (BT 709)
     *  Y =  0.2126 * r +  0.7152 * g +  0.0722 * b
     *
     *  if HowToGray = 2 (3 equal-weight)
     *  Y = 0.3333 * r + 0.3334 * g + 0.3333 * b
     */

        uint32 row = rgbIm.getRow();
        uint32 col = rgbIm.getCol();
        uint32 imsize = row * col;                

        Mat A(rgbIm.getDatType(), row, col, 1);
        _T *srcDat_pt = rgbIm.getDataPtr<_T>();
        _T *tarDat_pt = A.getDataPtr<_T>();

        uint32 ch_offset1 = imsize;
        uint32 ch_offset2 = imsize << 1;
        uint32 x,k;

        if( HowToGray==0){
            for(x= 0, k=0 ; x < imsize; x++, k++)
                tarDat_pt[k  ] =  bt601_r2y[0][0] * srcDat_pt[x  ] + bt601_r2y[0][1] * srcDat_pt[x+ch_offset1] + bt601_r2y[0][2] * srcDat_pt[x+ch_offset2];
        }else if( HowToGray==1){
            for(x= 0, k=0 ; x < imsize; x++, k++)
                tarDat_pt[k  ] =  bt709_r2y[0][0] * srcDat_pt[x  ] + bt709_r2y[0][1] * srcDat_pt[x+ch_offset1] + bt709_r2y[0][2] * srcDat_pt[x+ch_offset2];
        }else if( HowToGray==2){
            for(x= 0, k=0 ; x < imsize; x++, k++)
                tarDat_pt[k  ] =  0.3333 * srcDat_pt[x  ] + 0.3334 * srcDat_pt[x+ch_offset1] + 0.3333 * srcDat_pt[x+ch_offset2];
        }

        return A;
    }


    template <typename _T> inline Mat _conv_rgb2xyz(const Mat& rgbIm){
        uint32 row    = rgbIm.getRow();
        uint32 col    = rgbIm.getCol();
        uint32 chsize = rgbIm.getChannel();
        if( chsize != 3) return Mat();

        uint32 rcsize = rgbIm.getRowColSize();
        uint32 ch2    = rcsize;
        uint32 ch3    = rcsize << 1;

        Mat A(rgbIm.getDatType(), row, col, chsize);
        _T* dat64f = A.getDataPtr<_T>();
        _T x, y, z;
        int32 i2, i3;
        for(uint32 i=0; i < rcsize ; i++ ){
            i2 = i+ch2;
            i3 = i+ch3;
            x = rgb2xyz_bt709[0][0] * dat64f[i] + rgb2xyz_bt709[0][1] *dat64f[i2] + rgb2xyz_bt709[0][2] *dat64f[i3];
            y = rgb2xyz_bt709[1][0] * dat64f[i] + rgb2xyz_bt709[1][1] *dat64f[i2] + rgb2xyz_bt709[1][2] *dat64f[i3];
            z = rgb2xyz_bt709[2][0] * dat64f[i] + rgb2xyz_bt709[2][1] *dat64f[i2] + rgb2xyz_bt709[2][2] *dat64f[i3];
            dat64f[i ]= x;
            dat64f[i2]= y;
            dat64f[i3]= z;
        }

        return A;
    }
    template <typename _T> inline Mat _conv_xyz2rgb(const Mat& xyzIm){
        uint32 row    = xyzIm.getRow();
        uint32 col    = xyzIm.getCol();
        uint32 chsize = xyzIm.getChannel();
        if( chsize != 3) return Mat();

        uint32 rcsize = row * col;
        uint32 ch2    = rcsize;
        uint32 ch3    = rcsize << 1;

        Mat A(xyzIm.getDatType(), row, col, chsize);
        _T* dat64f = A.getDataPtr<_T>();
        _T r, g, b;
        int32 i2, i3;
        for(size_t i=0; i < rcsize; i++){
            i2 = i+ch2;
            i3 = i+ch3;
            r = xyz2rgb_bt709[0][0] * dat64f[i] + xyz2rgb_bt709[0][1] *dat64f[i2] + xyz2rgb_bt709[0][2] *dat64f[i3];
            g = xyz2rgb_bt709[1][0] * dat64f[i] + xyz2rgb_bt709[1][1] *dat64f[i2] + xyz2rgb_bt709[1][2] *dat64f[i3];
            b = xyz2rgb_bt709[2][0] * dat64f[i] + xyz2rgb_bt709[2][1] *dat64f[i2] + xyz2rgb_bt709[2][2] *dat64f[i3];
            dat64f[i ]= r;
            dat64f[i2]= g;
            dat64f[i3]= b;
        }

        return A;
    }

    template <typename _T> inline Mat _conv_rgb2Yxy(const Mat& rgbIm){
        uint32 row    = rgbIm.getRow();
        uint32 col    = rgbIm.getCol();
        uint32 chsize = rgbIm.getChannel();
        if( chsize != 3) return Mat();

        uint32 rcsize = rgbIm.getRowColSize();
        uint32 ch2    = rcsize;
        uint32 ch3    = rcsize << 1;

        Mat A(rgbIm.getDatType(), row, col, chsize);
        _T* dat64f = A.getDataPtr<_T>();

        _T X, Y, Z, W, x, y;
        //_T r,g,b;
        uint32 i2, i3;
        for(uint32 i=0; i < rcsize; i++){
            //r = dat64f[i    ];
            //g = dat64f[i+ch2];
            //b = dat64f[i+ch3];
            i2 = i+ch2;
            i3 = i+ch3;
            X = rgb2xyz_bt709[0][0] * dat64f[i] + rgb2xyz_bt709[0][1] *dat64f[i2] + rgb2xyz_bt709[0][2] *dat64f[i3];
            Y = rgb2xyz_bt709[1][0] * dat64f[i] + rgb2xyz_bt709[1][1] *dat64f[i2] + rgb2xyz_bt709[1][2] *dat64f[i3];
            Z = rgb2xyz_bt709[2][0] * dat64f[i] + rgb2xyz_bt709[2][1] *dat64f[i2] + rgb2xyz_bt709[2][2] *dat64f[i3];
            W = X + Y + Z;
            if( W <= 0.0) {
                x = 0.0;
                y = 0.0;
            }else{
                x = X/W;
                y = Y/W;
            }

            dat64f[i ]= Y;
            dat64f[i2]= x;
            dat64f[i3]= y;
        }

        return A;
    }

    template <typename _T> inline Mat _conv_Yxy2rgb(const Mat& YxyIm){
        uint32 row    = YxyIm.getRow();
        uint32 col    = YxyIm.getCol();
        uint32 chsize = YxyIm.getChannel();
        if( chsize != 3) return Mat();

        uint32 rcsize = YxyIm.getRowColSize();
        uint32 ch2    = rcsize;
        uint32 ch3    = rcsize << 1;

        Mat A(YxyIm.getDatType(), row, col, chsize);
        _T* dat64f = A.getDataPtr<_T>();

        _T  r, g, b;
        _T  X,Y,Z,x,y,W;
        uint32 i2, i3;
        for(uint32 i=0; i < rcsize; i++){
            i2 = i+ch2;
            i3 = i+ch3;
            Y = dat64f[i ];
            x = dat64f[i2];
            y = dat64f[i3];
            W = Y/y;
            X = x * W;
            Z = W-Y-X;
            r = xyz2rgb_bt709[0][0] *X + xyz2rgb_bt709[0][1] *Y + xyz2rgb_bt709[0][2] *Z;
            g = xyz2rgb_bt709[1][0] *X + xyz2rgb_bt709[1][1] *Y + xyz2rgb_bt709[1][2] *Z;
            b = xyz2rgb_bt709[2][0] *X + xyz2rgb_bt709[2][1] *Y + xyz2rgb_bt709[2][2] *Z;
            dat64f[i ]= r;
            dat64f[i2]= g;
            dat64f[i3]= b;
        }

        return A;
    }


    template <typename _T> inline Mat _histoPmf(const Mat& src, const int32 bins , const double step){

        Mat A = Mat::zeros(1,bins,1,DTYP::DOUBLE);
        double *tarDat_pt = A.getDataPtr<double>();
        _T *srcDat_pt = src.getDataPtr<_T>();

        int32 len = src.getLength();
        int32 k, d ;

        for(k=0; k < len ; k++){
            d = int32(srcDat_pt[k]/step);
            d = (d < 0) ? 0 : (d >= bins) ? bins-1: d;  // bin index range : 0 ~ bins-1
            tarDat_pt[d]++;
        }
        return A;
    }

    template <typename _T> inline Mat _histoCmf(const Mat& src, const int32 bins, const double step){

        Mat cmf = _histoPmf<_T>(src, bins, step);

        double *srcDat_pt = cmf.getDataPtr<double>();

        for(int32 k=1; k < bins ; k++)
            srcDat_pt[k] += srcDat_pt[k-1];

        return cmf;
    }

    template <typename _T> inline Mat _clip_HistoPmf(const Mat& src,const int32 clipVal, const int32 bins, const int32 step){

        Mat pmf = _histoPmf<_T>(src, bins, step );
        double *srcDat_pt = pmf.getDataPtr<double>();

        // clipping
        int32 sum_clipped =0;
        int32 binval;
        int32 k;
        for( k=0; k < bins ; k++){
            binval = int32(srcDat_pt[k]);
            if( binval > clipVal){
                sum_clipped += binval - clipVal;
                srcDat_pt[k] = clipVal;
            }
        }
        sum_clipped /= bins;
        // distributing the clipped sum
        for( k=0; k < bins ; k++){
            srcDat_pt[k] += sum_clipped;
        }

        return pmf;
    }

    template <typename _T> inline Mat _clip_HistoCmf(const Mat& src,const int32 clipVal, const int32 bins, const int32 step){

        Mat cmf = _clip_HistoPmf<_T>(src, clipVal, bins, step);
        double *srcDat_pt = cmf.getDataPtr<double>();

        // making cumiltive data
        for(int32 k=1; k < bins ; k++)
            srcDat_pt[k] += srcDat_pt[k-1];

        return cmf;
    }

    template <typename _T> inline Mat _clip_HistoEqual(const Mat& src, const Mat& histCmf, const int32 step){

        Mat A = src.copy();

        _T *srcDat_pt = src.getDataPtr<_T>();
        _T *tarDat_pt = A.getDataPtr<_T>();
        double *mapDat_pt = histCmf.getDataPtr<double>();

        int32 d0, d1, d2 , d3;
        int32 i0, i1;
        double mp1, mp2, mv;
        int32 bins = histCmf.getLength();
        int32 halfstep = step >> 1;
        int32 lowlmt = halfstep;
        int32 upplmt = (bins-1)*step + halfstep;
        for(uint32 i=0; i < src.getLength(); i++){
            d0 = int32(srcDat_pt[i]);
            d1 = d0 / step;
            d3 = d0 - d1*step;
            if( d0 < lowlmt ){
                tarDat_pt[i] = _T(mapDat_pt[0]);
            }else if( d0 >= upplmt ){
                tarDat_pt[i] = _T(mapDat_pt[bins-1]);
            }else {
                if( d3 < halfstep ){
                    i0 = d1-1;
                    i1 = d1;
                    d2 = halfstep + d3;
                }else{
                    i0 = d1;
                    i1 = d1+1;
                    d2 = d3 - halfstep;
                }
                mp1 = mapDat_pt[i0];
                mp2 = mapDat_pt[i1];
                mv  = (mp2-mp1)/d2;
                tarDat_pt[i] = _T(mp1 + mv);
            }
        }
        return A;
    }

template <typename _T> inline Mat _nakaSigTm(const Mat& Im, const Mat& localmean, const double globalmean, const double Imax){

    Mat Imt(Im.getDatType(), Im.getRow(), Im.getCol(), Im.getChannel());

    _T *ImDat_pt  = Im.getDataPtr<_T>();
    _T *ImtDat_pt = Imt.getDataPtr<_T>();
    _T *lmnDat_pt = localmean.getDataPtr<_T>();

    uint32 rc = Im.getRowColSize();
    uint32 ch = Im.getChannel();
    uint32 i, ich;
    for ( ich =0 ; ich < ch ; ++ich){
        for( i =0 ; i < rc; ++i )
            ImtDat_pt[i] = nakaSigmoid( ImDat_pt[i], lmnDat_pt[i] + globalmean, Imax);
    }
    return Imt;
}

template <typename _T> inline Mat _nakaSig3MeanTm(const Mat& Im, const Mat& smLocalMean, const Mat& lgLocalMean, const double globalmean, const double Imax){

    Mat Imt = Im.copy();

    _T *ImDat_pt  = Im.getDataPtr<_T>();
    _T *ImtDat_pt = Imt.getDataPtr<_T>();
    _T *slmnDat_pt = smLocalMean.getDataPtr<_T>();
    _T *llmnDat_pt = lgLocalMean.getDataPtr<_T>();

    uint32 rc = Im.getRowColSize();
    uint32 ch = Im.getChannel();
    uint32 i, ich;
    for ( ich =0 ; ich < ch ; ++ich){
        for( i =0 ; i < rc; ++i )
            ImtDat_pt[i] = nakaSigmoid( ImDat_pt[i], slmnDat_pt[i] + llmnDat_pt[i] + globalmean, Imax);
    }
    return Imt;
}

Mat inline localMeanMat ( const Mat& src, const Mat& mask){
    return  conv2d(src,mask,"symm","same");
}

double inline nakaSigmoid (const double X, const double X0, const double Xmax ){
    return (Xmax + X0 ) * X / (X + X0 ); // maximum return value is 0
}

template <typename _T> inline Mat _gamma(const Mat& src, const double gmval){
    Mat A = src.copy();
    double maxA = A.max().max().at<double>(0);
    A /= maxA;

    _T* adat_ptr = A.getDataPtr<_T>();
    uint32 len = A.getLength();

    for(uint32 i=0; i < len; ++i){
        adat_ptr[i] = pow( adat_ptr[i], gmval );
    }

    A *= maxA;
    return A;
}

template <typename _T> inline Mat _logRetinexTm(const Mat& Im, const Mat& surround){
    Mat Imt(Im.getDatType(),Im.getRow(), Im.getCol(), Im.getChannel()) ;
    _T *ImtDat_pt = Imt.getDataPtr<_T>();
    _T *ImDat_pt = Im.getDataPtr<_T>();
    _T *SrDat_pt = surround.getDataPtr<_T>();

    uint32 rc = Im.getRowColSize();
    uint32 ch = Im.getChannel();
    uint32 i, ich;
    double mn = 0;
    for ( ich =0 ; ich < ch ; ++ich){
        for( i =0 ; i < rc; ++i )
            mn += SrDat_pt[i];
    }
    mn /= rc*ch;

    for ( ich =0 ; ich < ch ; ++ich){
        for( i =0 ; i < rc; ++i )
            ImtDat_pt[i] = log(ImDat_pt[i])+log(ImDat_pt[i]/SrDat_pt[i]);//log(ImDat_pt[i]) - log(SrDat_pt[i]);
    }
    return Imt;
}

template <typename _T> inline Mat _bicubicIntp(const Mat& src, const int32 sw, const int32 sh){
    int32 mr = src.getRow();
    int32 mc = src.getCol();
    int32 mch = src.getChannel();
    DTYP  mdt = src.getDatType();
    int32 dr = mr*sh;
    int32 dc = mc*sw;
    int32 padw = 2;
    //Mat for boundary padding
    Mat A =	copy_padding(src,padw);
    int32 ar = A.getRow();
    int32 ac = A.getCol();
    int32 nc = ac*sw;
    int32 nr = ar*sh;
    Mat B = Mat::zeros(nr,nc,mch,mdt);

    // bicubic interpolation
    int32 y,x, oy, ox;
    _T *srcDat_p= A.getDataPtr<_T>();
    _T *desDat_p= B.getDataPtr<_T>();
    float xt,yt;
    float a_1[3],a0[3],a1[3],a2[3];  // the range 3 is for color channels
    float b_1[3],b0[3],b1[3],b2[3];
    int32 och_offset1 = ar*ac;
    int32 desCh_offset1 = nr*nc;
    int32 stp,k;
    int32 orow,pr_orow, n1_orow, n2_orow;
    int32 yRowByteStep;
    int32 xOff = padw*sw, yOff = padw*sh;
    for(y=yOff, yRowByteStep=nc*yOff ; y < nr-yOff; ++y, yRowByteStep+=nc){
        oy   = y/sh;
        orow = oy * ac;
        pr_orow = orow-ac;
        n1_orow = orow+ac;
        n2_orow = n1_orow+ac;
        for(x=xOff ; x < nc-xOff; ++x){
            ox = x/sw;
            xt = (float)(x-ox*sw)/sw;
            // f(x,-1)
            for( k=0, stp=pr_orow+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b_1[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }
            // f(x,0)
            for( k=0, stp=orow+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b0[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }
            // f(x,1)
            for( k=0, stp=n1_orow+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b1[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }
            // f(x,2)
            for( k=0, stp=n2_orow+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b2[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }

            yt =(float)(y-oy*sh)/sh;
            // f(x,y)
            for( k=0, stp=yRowByteStep+x; k < mch;++k, stp+=desCh_offset1 ){
                desDat_p[stp] =(_T) cubic1d(b_1[k],b0[k],b1[k],b2[k],yt);
            }
        }
    }

    Mat C(mdt,dr,dc,mch);
    matRect srcR(yOff,xOff,dr+yOff-1,dc+xOff-1);
    matRect tarR(0,0,dr-1,dc-1);
    Mat::sliceCopyMat(B,srcR,C,tarR);
    return C;
}

template <typename _T> inline Mat _bilinearIntp(const Mat& src, const int32 sw, const int32 sh){
    int32 mr = src.getRow();
    int32 mc = src.getCol();
    int32 mch = src.getChannel();
    DTYP  mdt = src.getDatType();
    int32 dr = mr*sh;
    int32 dc = mc*sw;
    int32 padw = 2;
    //Mat for boundary padding
    Mat A =	copy_padding(src,padw);
    int32 ar = A.getRow();
    int32 ac = A.getCol();
    int32 nc = ac*sw;
    int32 nr = ar*sh;
    Mat B = Mat::zeros(nr,nc,mch,mdt);

    // bilinear interpolation
    int32 y,x, oy, ox;
    _T *srcDat_p= A.getDataPtr<_T>();
    _T *desDat_p= B.getDataPtr<_T>();
    int32 pxt, pyt;
    double xt, yt;
    double a0[3],a1[3]; // 3 is for the color channel
    int32 och_offset = ar*ac;
    int32 desCh_offset = nr*nc;
    int32 stp,dstp,k;
    int32 orow, n1_orow;
    int32 yRowByteStep;
    int32 xOff = padw*sw, yOff = padw*sh;
    for(y=yOff, yRowByteStep=nc*yOff ; y < nr-yOff; ++y, yRowByteStep+=nc){
        oy   = y/sh;
        pyt  = y-oy*sh;
        orow = oy * ac;
        n1_orow = orow+ac;
        for(x=xOff ; x < nc-xOff; ++x){
            ox = x/sw;
            pxt = x-ox*sw;
            if(pxt==0 && pyt==0){
                for( k=0, dstp=yRowByteStep+x, stp=orow+ox; k < mch; ++k, stp+= och_offset, dstp+=desCh_offset)
                    desDat_p[dstp] = srcDat_p[stp];
            }else{
                xt = (double)pxt/sw;
                // f(x,0)
                for( k=0, stp=orow+ox; k < mch; ++k, stp+= och_offset){
                    a0[k] = ((double)srcDat_p[stp+1] -(double)srcDat_p[stp])*xt+(double)srcDat_p[stp];
                }
                // f(x,1)
                for( k=0, stp=n1_orow+ox; k < mch; ++k, stp+= och_offset){
                    a1[k] = ((double)srcDat_p[stp+1] -(double)srcDat_p[stp])*xt+(double)srcDat_p[stp];
                }

                yt =(double)pyt/sh;
                // f(x,y)
                for( k=0, dstp=yRowByteStep+x; k < mch;++k, dstp+=desCh_offset ){
                    desDat_p[dstp] = (_T)((a1[k]-a0[k])*yt+a0[k]);
                }
            }
        }
    }

    Mat C(mdt,dr,dc,mch);
    matRect srcR(yOff,xOff,dr+yOff-1,dc+xOff-1);
    matRect tarR(0,0,dr-1,dc-1);
    Mat::sliceCopyMat(B,srcR,C,tarR);
    return C;
}
template <typename _T> inline Mat _nearestIntp(const Mat& src, const int32 sw, const int32 sh){
    int32 mr = src.getRow();
    int32 mc = src.getCol();
    int32 mch = src.getChannel();
    DTYP  mdt = src.getDatType();
    int32 nr = mr*sh;
    int32 nc = mc*sw;
    Mat A = Mat::zeros(nr,nc,mch,mdt);

    // nearest interpolation
    int32 y,x, oy, ox;
    _T *srcDat_p= src.getDataPtr<_T>();
    _T *desDat_p= A.getDataPtr<_T>();
    int32 och_offset1   = mr*mc;
    int32 desCh_offset1 = nr*nc;
    int32 stp,k, orow, ostp;
    int32 yRowByteStep;
    for(y=0, yRowByteStep=0 ; y < nr; ++y, yRowByteStep+=nc){
        oy   = y/sh;
        orow = oy * mc;
        for(x=0 ; x < nc; ++x){
            ox = x/sw;
            ostp = orow;
            for(k=0, stp=yRowByteStep+x; k < mch; ++k, stp+=desCh_offset1){
                desDat_p[stp] =srcDat_p[ostp+ox];
                ostp += och_offset1;
            }
        }
    }

    return A;
}

template <typename _T> inline Mat _decimate(const Mat& src, const int32 sw, const int32 sh){
    int32 mr = src.getRow();
    int32 mc = src.getCol();
    int32 mch = src.getChannel();
    DTYP  mdt = src.getDatType();
    int32 nr = mr/sh;
    int32 nc = mc/sw;
    Mat A = Mat::zeros(nr,nc,mch,mdt);

    // nearest interpolation
    int32 y,x, oy, ox;
    _T *srcDat_p= src.getDataPtr<_T>();
    _T *desDat_p= A.getDataPtr<_T>();
    int32 och_offset   = src.getRowColSize(); // mr*mc
    int32 desCh_offset = A.getRowColSize(); //nr*nc;
    int32 stp,k, orow, ostp;
    int32 yRowByteStep;
    for(y=0, yRowByteStep=0 ; y < nr; ++y, yRowByteStep+=nc){
        oy   = y*sh;
        orow = oy * mc;
        for(x=0 ; x < nc; ++x){
            ox = x*sw;
            ostp = orow;
            for(k=0, stp=yRowByteStep+x; k < mch; ++k, stp+=desCh_offset){
                desDat_p[stp] =srcDat_p[ostp+ox];
                ostp += och_offset;
            }
        }
    }
    return A;
}

template <typename _T> inline Mat _gaussPyramid(const Mat& src, const int32 level){
    int32 r=src.getRow();
    int32 c=src.getCol();
    int32 ch=src.getChannel();

    Mat A=Mat::zeros(r,c+c/2,ch,src.getDatType());
    _T *srcDat_p = src.getDataPtr<_T>();
    _T *desDat_p = A.getDataPtr<_T>();

    int32 k;
    int32 sr, sc, lr, lc;
    Mat flt = boxMaskGen(3,ch);
    Mat B = src;
    Mat C ;
    matRect srcR(0,0,r-1,c-1);
    matRect tarR(0,0,r-1,c-1);
    Mat::sliceCopyMat(src,srcR,A,tarR);
    lr = r; lc = c;
    sr = 0;
    sc = 0;
    for(k=0; k< level; ++k){
        C = jmat::conv2d(B,flt,"symm");
        B = decimate(C,2,2);
        if( k & 0x00000001 ) sr +=lr;
        else                 sc +=lc;
        lr = B.getRow();
        lc = B.getCol();
        srcR.set(0,0,lr-1,lc-1);
        tarR.set(sr,sc,sr+lr-1,sc+lc-1);
        Mat::sliceCopyMat(B,srcR,A,tarR);
    }
    return A;
}

template <typename _T> inline Mat _laplPyramid(Mat& src, const int32 level){
    int32 r=src.getRow();
    int32 c=src.getCol();
    int32 ch=src.getChannel();

    Mat A=Mat::zeros(r,c+c/2,ch,src.getDatType());
    _T *srcDat_p = src.getDataPtr<_T>();
    _T *desDat_p = A.getDataPtr<_T>();

    int32 k;
    int32 sr, sc, lr, lc;
    Mat flt = boxMaskGen(3,ch);
    Mat B = src;
    Mat C, D, E, F;
    matRect srcR(0,0,r-1,c-1);
    matRect tarR(0,0,r-1,c-1);
    sr = 0;
    sc = 0;
    for(k=0; k< level; ++k){
        C = jmat::conv2d(B,flt,"symm");
        D = decimate(C,2,2);
        E = bilinearIntp(D,2,2);
        F = B - E + 128;
        lr = F.getRow();
        lc = F.getCol();
        srcR.set(0,0,lr-1,lc-1);
        tarR.set(sr,sc,sr+lr-1,sc+lc-1);
        Mat::sliceCopyMat(F,srcR,A,tarR);
        if( k & 0x00000001 ) sr +=lr;
        else                 sc +=lc;
        B = D.copy();
    }

    lr = D.getRow();
    lc = D.getCol();
    srcR.set(0,0,lr-1,lc-1);
    tarR.set(sr,sc,sr+lr-1,sc+lc-1);
    Mat::sliceCopyMat(D,srcR,A,tarR);
    return A;
}
} // end of imgproc namespace
} // end of jmat namespace
#endif // JBIMGPROC_H

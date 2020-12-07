/*
 * Copyright (C) 2020. Jong B. Choi
 * License : MIT License
 * contact : neoflame99@naver.com
 */

#ifndef JBIMGPROC_H
#define JBIMGPROC_H
#include "jbMat.h"
#include "jbmath.h"
#ifdef _WIN_
    #define _USE_MATH_DEFINES
#endif
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
    float cubic1d(float ma1, float a0, float a1, float a2, float t);
    Mat bicubicIntp(const Mat& src,const int32 s);
    Mat copy_padding(const Mat& src, int32 pad_size=2);
    template <typename _T> inline Mat _bicubicIntp(const Mat& src, int32 s);

    template <typename T> class bgr_g{
        public: T B, G, R;
    };
    template <typename T> class yuv_g{
        public: T Y, U, V;
    };
    template <typename T> class xyz_g{
        public: T X, Y, Z;
    };
    template <typename T> class Yxy_g{
        public: T Y, x, y;
    };
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
        uint32 chsize = rgbIm.getChannel();
        assert(chsize==3);

        Mat  A(rgbIm.getDatType(), row, col, chsize);
        bgr_g<_T>* bgr = (bgr_g<_T> *)rgbIm.getDataPtr();
        yuv_g<_T>* yuv = (yuv_g<_T> *)A.getDataPtr();

        uint32 x;
        if( sel_eq == 0){
            for(x=0 ; x < rgbIm.getRowColSize(); ++x){
                yuv[x].Y = (299*bgr[x].R + 587*bgr[x].G + 114*bgr[x].B)/1000; // Y
                yuv[x].U = (-168736*bgr[x].R -331264*bgr[x].G + 500000*bgr[x].B)/1000000; // U
                yuv[x].V = (500000*bgr[x].R -418688*bgr[x].G - 81312*bgr[x].B)/1000000; // V
            }
        }else if( sel_eq == 1){
            for(x=0 ; x < rgbIm.getRowColSize(); ++x){
                yuv[x].Y = (2126*bgr[x].R + 7152*bgr[x].G + 722*bgr[x].B)/10000; // Y
                yuv[x].U = (-11457*bgr[x].R -38543*bgr[x].G + 50000*bgr[x].B)/100000; // U
                yuv[x].V = (50000*bgr[x].R -45415*bgr[x].G -4585*bgr[x].B)/100000; // V
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
         uint32 chsize = yccIm.getChannel();
         assert(chsize==3);

         Mat A(yccIm.getDatType(), row, col, chsize);
         yuv_g<_T> *yuv = (yuv_g<_T> *)yccIm.getDataPtr();
         bgr_g<_T> *bgr = (bgr_g<_T> *)A.getDataPtr();
         uint32 x ;
         if( sel_eq == 0){
             for(x= 0 ; x < yccIm.getRowColSize(); ++x ){
                 bgr[x].R = (10000*yuv[x].Y + 14020*yuv[x].V)/10000;
                 bgr[x].G = (10000*yuv[x].Y - 3441*yuv[x].U -7141*yuv[x].V)/10000;
                 bgr[x].B = (10000*yuv[x].Y + 17720*yuv[x].U)/10000;

                 bgr[x].B = (bgr[x].B < 0 ) ? 0 : bgr[x].B;
                 bgr[x].G = (bgr[x].G < 0 ) ? 0 : bgr[x].G;
                 bgr[x].R = (bgr[x].R < 0 ) ? 0 : bgr[x].R;
             }
         }else if( sel_eq == 1){
             for(x= 0 ; x < yccIm.getRowColSize(); ++x){
                 bgr[x].R = (10000*yuv[x].Y + 15748*yuv[x].V)/10000;
                 bgr[x].G = (10000*yuv[x].Y - 1873*yuv[x].U - 4681*yuv[x].V)/10000;
                 bgr[x].B = (10000*yuv[x].Y + 18556*yuv[x].U)/10000;

                 bgr[x].B = (bgr[x].B < 0 ) ? 0 : bgr[x].B;
                 bgr[x].G = (bgr[x].G < 0 ) ? 0 : bgr[x].G;
                 bgr[x].R = (bgr[x].R < 0 ) ? 0 : bgr[x].R;
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
        uint32 imsize = rgbIm.getRowColSize();
        uint32 ch  = rgbIm.getChannel();
        assert(ch==3);

        Mat A(rgbIm.getDatType(), row, col, 1);
        bgr_g<_T> *bgr = (bgr_g<_T> *) rgbIm.getDataPtr();
        _T *gray = A.getDataPtr<_T>();

        uint32 k;
        if( HowToGray==0){
            for( k=0 ; k < imsize; k++)
                gray[k] = (299*bgr[k].R + 587*bgr[k].G + 114*bgr[k].B)/1000;
        }else if( HowToGray==1){
            for( k=0 ; k < imsize; k++)
                gray[k] = (2126*bgr[k].R + 7152*bgr[k].G + 722*bgr[k].B)/10000;
        }else if( HowToGray==2){
            for( k=0 ; k < imsize; k++)
                gray[k] = (3333*bgr[k].R + 3334*bgr[k].G + 3333*bgr[k].B)/10000;
        }
        return A;
    }

    template <typename _T> inline Mat _conv_rgb2xyz(const Mat& rgbIm){
        /*
         * rgb -> XYZ
         * [0.4124564, 0.3575761, 0.1804375] [ r ]
         * [0.2126729, 0.7151522, 0.0721750]*[ g ]
         * [0.0193339, 0.1191920, 0.9503041] [ b ]
         */
        uint32 row    = rgbIm.getRow();
        uint32 col    = rgbIm.getCol();
        uint32 chsize = rgbIm.getChannel();
        if( chsize != 3) return Mat();

        Mat A(rgbIm.getDatType(), row, col, chsize);
        bgr_g<_T> *bgr = (bgr_g<_T> *)rgbIm.getDataPtr();
        xyz_g<_T> *xyz = (xyz_g<_T> *)A.getDataPtr();
        uint32 i;
        for( i=0; i < rgbIm.getRowColSize() ; ++i ){
            xyz[i].X = (412456*bgr[i].R + 357576*bgr[i].G + 180438*bgr[i].B)/1000000; // X
            xyz[i].Y = (212673*bgr[i].R + 715152*bgr[i].G +  72175*bgr[i].B)/1000000; // Y
            xyz[i].Z = ( 19334*bgr[i].R + 119192*bgr[i].G + 950304*bgr[i].B)/1000000; // Z
        }
        return A;
    }
    template <typename _T> inline Mat _conv_xyz2rgb(const Mat& xyzIm){
        /*
         * XYZ -> rgb
         * [ 3.2404542, -1.5371385, -0.4985314]  [ X ]
         * [-0.9692660,  1.8760108,  0.0415560] *[ Y ]
         * [ 0.0556434, -0.2040259,  1.0572252]  [ Z ]
         */
        uint32 row    = xyzIm.getRow();
        uint32 col    = xyzIm.getCol();
        uint32 chsize = xyzIm.getChannel();
        if( chsize != 3) return Mat();

        Mat A(xyzIm.getDatType(), row, col, chsize);
        xyz_g<_T> *xyz = (xyz_g<_T> *)xyzIm.getDataPtr();
        bgr_g<_T> *bgr = (bgr_g<_T> *)A.getDataPtr();
        uint32 i;
        for(i=0; i < xyzIm.getRowColSize(); ++i) {
            bgr[i].R = ( 3240454*xyz[i].X - 1537139*xyz[i].Y -  498531*xyz[i].Z)/1000000; // R
            bgr[i].G = ( -969266*xyz[i].X + 1876011*xyz[i].Y +   41556*xyz[i].Z)/1000000; // G
            bgr[i].B = (   55643*xyz[i].X -  204026*xyz[i].Y + 1057225*xyz[i].Z)/1000000; // B

            bgr[i].R = bgr[i].R < 0 ? 0 : bgr[i].R;
            bgr[i].G = bgr[i].G < 0 ? 0 : bgr[i].G;
            bgr[i].B = bgr[i].B < 0 ? 0 : bgr[i].B;
        }
        return A;
    }

    template <typename _T> inline Mat _conv_rgb2Yxy(const Mat& rgbIm){
        uint32 row    = rgbIm.getRow();
        uint32 col    = rgbIm.getCol();
        uint32 chsize = rgbIm.getChannel();
        if( chsize != 3) return Mat();

        Mat A(rgbIm.getDatType(), row, col, chsize);
        bgr_g<_T> *bgr = (bgr_g<_T> *)rgbIm.getDataPtr();
        Yxy_g<_T> *Yxy = (Yxy_g<_T> *)A.getDataPtr();

        _T X, Y, Z, W, x, y;
        uint32 i;
        for( i=0 ; i < rgbIm.getRowColSize(); ++i ){
            X = (412456*bgr[i].R + 357576*bgr[i].G + 180438*bgr[i].B)/1000000; // X
            Y = (212673*bgr[i].R + 715152*bgr[i].G +  72175*bgr[i].B)/1000000; // Y
            Z = ( 19334*bgr[i].R + 119192*bgr[i].G + 950304*bgr[i].B)/1000000; // Z
            W = X + Y + Z;
            if( W <= 0.0) {
                x = 0.0;
                y = 0.0;
            }else{
                x = X/W;
                y = Y/W;
            }
            Yxy[i].Y = Y;
            Yxy[i].x = x;
            Yxy[i].y = y;
        }
        return A;
    }

    template <typename _T> inline Mat _conv_Yxy2rgb(const Mat& YxyIm){
        uint32 row    = YxyIm.getRow();
        uint32 col    = YxyIm.getCol();
        uint32 chsize = YxyIm.getChannel();
        if( chsize != 3) return Mat();

        Mat A(YxyIm.getDatType(), row, col, chsize);
        Yxy_g<_T> *Yxy = (Yxy_g<_T> *)YxyIm.getDataPtr();
        bgr_g<_T> *bgr = (bgr_g<_T> *)A.getDataPtr();

        _T  X,Z,W;
        uint32 i;
        for(i=0 ; i < YxyIm.getRowColSize(); ++i){
            W = Yxy[i].Y/Yxy[i].y;
            X = Yxy[i].x * W;
            Z = W-Yxy[i].Y-X;
            bgr[i].R = ( 3240454*X - 1537139*Yxy[i].Y -  498531*Z)/1000000; // R
            bgr[i].G = ( -969266*X + 1876011*Yxy[i].Y +   41556*Z)/1000000; // G
            bgr[i].B = (   55643*X -  204026*Yxy[i].Y + 1057225*Z)/1000000; // B

            bgr[i].R = bgr[i].R < 0 ? 0 : bgr[i].R;
            bgr[i].G = bgr[i].G < 0 ? 0 : bgr[i].G;
            bgr[i].B = bgr[i].B < 0 ? 0 : bgr[i].B;
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

template <typename _T> inline Mat _bicubicIntp(const Mat& src, const int32 s){
    int32 mr = src.getRow();
    int32 mc = src.getCol();
    int32 mch = src.getChannel();
    DTYP  mdt = src.getDatType();

    int32 dr = mr*s;
    int32 dc = mc*s;
    int32 padw = 2;
    //Mat for boundary padding
    Mat A =	copy_padding(src,padw);
    int32 ar = A.getRow();
    int32 ac = A.getCol();
    int32 nc = ac*s;
    int32 nr = ar*s;
    Mat B = Mat::zeros(nr,nc,mch,mdt);

    // bicubic interpolation
    int32 y,x, oy, ox;
    _T *srcDat_p= A.getDataPtr<_T>();
    _T *desDat_p= B.getDataPtr<_T>();
    float xt,yt;
    float a_1[3],a0[3],a1[3],a2[3];
    float b_1[3],b0[3],b1[3],b2[3];
    int32 och_offset1 = ar*ac;
    int32 desCh_offset1 = nr*nc;
    int32 stp,k, orow;
    int32 yRowByteStep;
    int32 xOff = padw*s, yOff = padw*s;
    for(y=yOff, yRowByteStep=nc*yOff ; y < nr-yOff; ++y, yRowByteStep+=nc){
        oy   = y/s;
        orow = oy * ac;
        for(x=xOff ; x < nc-xOff; ++x){
            ox = x/s;
            xt = (float)(x-ox*s)/s;
            // f(x,-1)
            for( k=0, stp=orow-ac+ox; k < mch; ++k, stp+= och_offset1){
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
            for( k=0, stp=orow+ac+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b1[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }
            // f(x,2)
            for( k=0, stp=orow+ac*2+ox; k < mch; ++k, stp+= och_offset1){
                a_1[k]= srcDat_p[stp-1];
                a0[k] = srcDat_p[stp  ];
                a1[k] = srcDat_p[stp+1];
                a2[k] = srcDat_p[stp+2];
                b2[k] = cubic1d(a_1[k],a0[k],a1[k],a2[k],xt);
            }

            yt =(float)(y-oy*s)/s;
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

} // end of imgproc namespace
} // end of jmat namespace
#endif // JBIMGPROC_H

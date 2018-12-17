#ifndef JBIMGPROC_H
#define JBIMGPROC_H
#include "jbMat.h"

namespace jmat {
namespace imgproc {

    const double bt601_r2y[3][3]={{ 0.299    ,  0.587    ,  0.114 },
                                  {-0.168736 , -0.331264 ,  0.5   },
                                  { 0.5      , -0.418688 , -0.081312} };
    const double bt709_r2y[3][3]={{ 0.2126  ,  0.7152  , 0.0722  },
                                  {-0.11457 , -0.38543 , 0.5     },
                                  { 0.5     , -0.45415 ,-0.04585 }};

    const double bt601_y2r[3][3]={{ 1.0000 ,  -0.0000 ,   1.4020 },
                                  { 1.0000 ,  -0.3441 ,  -0.7141 },
                                  { 1.0000 ,   1.7720 ,   0.0000 }};

    const double bt709_y2r[3][3]={{ 1.0000 ,   0.0000 ,   1.5748 },
                                  { 1.0000 ,  -0.1873 ,  -0.4681 },
                                  { 1.0000 ,   1.8556 ,  -0.0000 }};

    Mat rgb2ycc(const Mat& rgbIm, const int32 sel_eq = 0);
    Mat ycc2rgb(const Mat& yccIm, const int32 sel_eq = 0);
    Mat rgb2gray(const Mat& rgbIm, const int32 HowToGray = 0);
    Mat histoPmf(const Mat& src,const int32 bins, const int32 step);
    Mat histoCmf(const Mat& src,const int32 bins, const int32 step);
    Mat clip_HistoPmf(const Mat& src, const int32 clipVal,const int32 bins, const int32 step);
    Mat clip_HistoCmf(const Mat& src, const int32 clipVal,const int32 bins, const int32 step);
    Mat clip_HistoEqual(const Mat& src, const Mat& histCmf, const int32 step);


    template <typename _T> Mat _rgb2ycc(const Mat& rgbIm, const int32 sel_eq );
    template <typename _T> Mat _ycc2rgb(const Mat& yccIm, const int32 sel_eq );
    template <typename _T> Mat _rgb2gray(const Mat& rgbIm, const int32 HowToGray);



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
        if( chsize != 3 ) {
            fprintf(stdout,"rgbIm is not three channel image\n");
            return Mat();
        }

        Mat  A(rgbIm.getDatType(), row, col, chsize);
        _T* srcDat_pt = rgbIm.getDataPtr<_T>();
        _T* tarDat_pt = A.getDataPtr<_T>();

        uint32 ch_offset1 = imsize;
        uint32 ch_offset2 = imsize << 1;

        uint32 x;
        if( sel_eq == 0){
            for(x= 0 ; x < imsize; x++ ){
                tarDat_pt[x           ] = bt601_r2y[0][0] * srcDat_pt[x  ] + bt601_r2y[0][1] * srcDat_pt[x+ch_offset1] + bt601_r2y[0][2] * srcDat_pt[x+ch_offset2];
                tarDat_pt[x+ch_offset1] = bt601_r2y[1][0] * srcDat_pt[x  ] + bt601_r2y[1][1] * srcDat_pt[x+ch_offset1] + bt601_r2y[1][2] * srcDat_pt[x+ch_offset2];
                tarDat_pt[x+ch_offset2] = bt601_r2y[2][0] * srcDat_pt[x  ] + bt601_r2y[2][1] * srcDat_pt[x+ch_offset1] + bt601_r2y[2][2] * srcDat_pt[x+ch_offset2];
            }
        }else if( sel_eq == 1){
            for(x= 0 ; x < imsize; x++ ){
                tarDat_pt[x           ] = bt709_r2y[0][0] * srcDat_pt[x  ] + bt709_r2y[0][1] * srcDat_pt[x+ch_offset1] + bt709_r2y[0][2] * srcDat_pt[x+ch_offset2];
                tarDat_pt[x+ch_offset1] = bt709_r2y[1][0] * srcDat_pt[x  ] + bt709_r2y[1][1] * srcDat_pt[x+ch_offset1] + bt709_r2y[1][2] * srcDat_pt[x+ch_offset2];
                tarDat_pt[x+ch_offset2] = bt709_r2y[2][0] * srcDat_pt[x  ] + bt709_r2y[2][1] * srcDat_pt[x+ch_offset1] + bt709_r2y[2][2] * srcDat_pt[x+ch_offset2];
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
         if( chsize != 3 ) {
             fprintf(stdout,"rgbIm is not three channel image\n");
             return Mat();
         }

         Mat A(yccIm.getDatType(), row, col, chsize);
         _T *srcDat_pt = yccIm.getDataPtr<_T>();
         _T *tarDat_pt = A.getDataPtr<_T>();

         uint32 ch_offset1 = imsize ;
         uint32 ch_offset2 = imsize << 1;
         uint32 x;
         _T tmp1, tmp2, tmp3;
         if( sel_eq == 0){
             for(x= 0 ; x < imsize; x++){
                 tmp1 = bt601_y2r[0][0] * srcDat_pt[x  ] + bt601_y2r[0][1] * srcDat_pt[x+ch_offset1] + bt601_y2r[0][2] * srcDat_pt[x+ch_offset2];
                 tmp2 = bt601_y2r[1][0] * srcDat_pt[x  ] + bt601_y2r[1][1] * srcDat_pt[x+ch_offset1] + bt601_y2r[1][2] * srcDat_pt[x+ch_offset2];
                 tmp3 = bt601_y2r[2][0] * srcDat_pt[x  ] + bt601_y2r[2][1] * srcDat_pt[x+ch_offset1] + bt601_y2r[2][2] * srcDat_pt[x+ch_offset2];

                 tarDat_pt[x           ] = (tmp1 < 0 ) ? 0 : tmp1;
                 tarDat_pt[x+ch_offset1] = (tmp2 < 0 ) ? 0 : tmp2;
                 tarDat_pt[x+ch_offset2] = (tmp3 < 0 ) ? 0 : tmp3;
             }
         }else if( sel_eq == 1){
             for(x= 0 ; x < imsize; x++){
                 tmp1 = bt709_y2r[0][0] * srcDat_pt[x  ] + bt709_y2r[0][1] * srcDat_pt[x+ch_offset1] + bt709_y2r[0][2] * srcDat_pt[x+ch_offset2];
                 tmp2 = bt709_y2r[1][0] * srcDat_pt[x  ] + bt709_y2r[1][1] * srcDat_pt[x+ch_offset1] + bt709_y2r[1][2] * srcDat_pt[x+ch_offset2];
                 tmp3 = bt709_y2r[2][0] * srcDat_pt[x  ] + bt709_y2r[2][1] * srcDat_pt[x+ch_offset1] + bt709_y2r[2][2] * srcDat_pt[x+ch_offset2];

                 tarDat_pt[x           ] = (tmp1 < 0 ) ? 0 : tmp1;
                 tarDat_pt[x+ch_offset1] = (tmp2 < 0 ) ? 0 : tmp2;
                 tarDat_pt[x+ch_offset2] = (tmp3 < 0 ) ? 0 : tmp3;

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
     *  Y = 0.333 * r + 0.334 * g + 0.333 * b
     */

        uint32 row = rgbIm.getRow();
        uint32 col = rgbIm.getCol();
        uint32 imsize = row * col;
        uint32 chsize = rgbIm.getChannel();
        if( chsize != 3 ) {
            fprintf(stdout,"rgbIm is not three channel image\n");
            return Mat();
        }

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
                tarDat_pt[k  ] =  0.333 * srcDat_pt[x  ] + 0.334 * srcDat_pt[x+ch_offset1] + 0.333 * srcDat_pt[x+ch_offset2];
        }

        return A;
    }

    template <typename _T> inline Mat _histoPmf(const Mat& src, const int32 bins , const int32 step){

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

    template <typename _T> inline Mat _histoCmf(const Mat& src, const int32 bins, const int32 step){

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

    template <typename _T> Mat _clip_HistoEqual(const Mat& src, const Mat& histCmf, const int32 step){

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
        for(int32 i=0; i < src.getLength(); i++){
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
} // end of imgproc namespace
} // end of jmat namespace
#endif // JBIMGPROC_H

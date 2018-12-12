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

    Mat rgb2ycc(const Mat& rgbIm, const int sel_eq = 0);
    Mat ycc2rgb(const Mat& yccIm, const int sel_eq = 0);
    Mat rgb2gray(const Mat& rgbIm, const int HowToGray = 0);
    Mat histoPmf(const Mat& src);
    Mat histoCmf(const Mat& src);
    Mat clip_HistoPmf(const Mat& src, const unsigned int clipVal);
    Mat clip_HistoCmf(const Mat& src, const unsigned int clipVal);
    Mat clip_HistoEqual(const Mat& src, const Mat& histCmf);


    template <typename _T> Mat _rgb2ycc(const Mat& rgbIm, const int sel_eq );
    template <typename _T> Mat _ycc2rgb(const Mat& yccIm, const int sel_eq );
    template <typename _T> Mat _rgb2gray(const Mat& rgbIm, const int HowToGray);



    template <typename _T> Mat _rgb2ycc(const Mat& rgbIm, const int sel_eq ){
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

        uint row = rgbIm.getRow();
        uint col = rgbIm.getCol();
        uint imsize = row * col;
        uint chsize = rgbIm.getChannel();
        if( chsize != 3 ) {
            fprintf(stdout,"rgbIm is not three channel image\n");
            return Mat();
        }

        Mat  A(rgbIm.getDatType(), row, col, chsize);
        _T* srcDat_pt = rgbIm.getDataPtr<_T>();
        _T* tarDat_pt = A.getDataPtr<_T>();

        uint ch_offset1 = imsize;
        uint ch_offset2 = imsize << 1;

        uint x;
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
    template <typename _T> Mat _ycc2rgb(const Mat& yccIm, const int sel_eq ){
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

         uint row = yccIm.getRow();
         uint col = yccIm.getCol();
         uint imsize = row * col;
         uint chsize = yccIm.getChannel();
         if( chsize != 3 ) {
             fprintf(stdout,"rgbIm is not three channel image\n");
             return Mat();
         }

         Mat A(yccIm.getDatType(), row, col, chsize);
         _T *srcDat_pt = yccIm.getDataPtr<_T>();
         _T *tarDat_pt = A.getDataPtr<_T>();

         uint ch_offset1 = imsize ;
         uint ch_offset2 = imsize << 1;
         uint x;
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


    template <typename _T> Mat _rgb2gray(const Mat& rgbIm, const int HowToGray){
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

        uint row = rgbIm.getRow();
        uint col = rgbIm.getCol();
        uint imsize = row * col;
        uint chsize = rgbIm.getChannel();
        if( chsize != 3 ) {
            fprintf(stdout,"rgbIm is not three channel image\n");
            return Mat();
        }

        Mat A(rgbIm.getDatType(), row, col, 1);
        _T *srcDat_pt = rgbIm.getDataPtr<_T>();
        _T *tarDat_pt = A.getDataPtr<_T>();

        uint ch_offset1 = imsize;
        uint ch_offset2 = imsize << 1;
        uint x,k;

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
} // end of imgproc namespace
} // end of jmat namespace
#endif // JBIMGPROC_H

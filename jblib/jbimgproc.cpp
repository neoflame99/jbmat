#include "jbimgproc.h"
#include "jbimgproc.h"


namespace jmat {

namespace imgproc {

Mat rgb2ycc(const Mat& rgbIm, const int32 sel_eq){
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

    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _rgb2ycc<double>(rgbIm, sel_eq);
    case DTYP::FLOAT  : return _rgb2ycc<float >(rgbIm, sel_eq);
    case DTYP::INT    : return _rgb2ycc<int32 >(rgbIm, sel_eq);
    case DTYP::UCHAR  : return _rgb2ycc<uchar >(rgbIm, sel_eq);
    default           : {
        fprintf(stderr, " Unsupported DTYP in rgb2ycc func.\n");
        return Mat();
        }
    }
}


Mat ycc2rgb(const Mat& rgbIm, const int32 sel_eq){
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

    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _ycc2rgb<double>(rgbIm, sel_eq);
    case DTYP::FLOAT  : return _ycc2rgb<float >(rgbIm, sel_eq);
    case DTYP::INT    : return _ycc2rgb<int32 >(rgbIm, sel_eq);
    case DTYP::UCHAR  : return _ycc2rgb<uchar >(rgbIm, sel_eq);
    default           : {
        fprintf(stderr, " Unsupported DTYP in ycc2rgb func.\n");
        return Mat();
        }
    }
}


Mat rgb2gray(const Mat& rgbIm, const int32 HowToGray){
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
    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _rgb2gray<double>(rgbIm, HowToGray);
    case DTYP::FLOAT  : return _rgb2gray<float >(rgbIm, HowToGray);
    case DTYP::INT    : return _rgb2gray<int32 >(rgbIm, HowToGray);
    case DTYP::UCHAR  : return _rgb2gray<uchar >(rgbIm, HowToGray);
    default           : {
        fprintf(stderr, " Unsupported DTYP in rgb2gray func.\n");
        return Mat();
        }
    }
}

Mat rgb2xyz(const Mat& rgbIm){
    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _conv_rgb2xyz<double>(rgbIm );
    case DTYP::FLOAT  : return _conv_rgb2xyz<float >(rgbIm );
    case DTYP::INT    : return _conv_rgb2xyz<int32 >(rgbIm );
    case DTYP::UCHAR  : return _conv_rgb2xyz<uchar >(rgbIm );
    default           : {
        fprintf(stderr, " Unsupported DTYP in rgb2xyz func.\n");
        return Mat();
        }
    }
}

Mat xyz2rgb(const Mat& xyzIm){
    switch(xyzIm.getDatType()){
    case DTYP::DOUBLE : return _conv_xyz2rgb<double>(xyzIm );
    case DTYP::FLOAT  : return _conv_xyz2rgb<float >(xyzIm );
    case DTYP::INT    : return _conv_xyz2rgb<int32 >(xyzIm );
    case DTYP::UCHAR  : return _conv_xyz2rgb<uchar >(xyzIm );
    default           : {
        fprintf(stderr, " Unsupported DTYP in xyz2rgb func.\n");
        return Mat();
        }
    }
}

Mat rgb2Yxy(const Mat& rgbIm){
    switch(rgbIm.getDatType()){
    case DTYP::DOUBLE : return _conv_rgb2Yxy<double>(rgbIm );
    case DTYP::FLOAT  : return _conv_rgb2Yxy<float >(rgbIm );
    case DTYP::INT    : return _conv_rgb2Yxy<int32 >(rgbIm );
    case DTYP::UCHAR  : return _conv_rgb2Yxy<uchar >(rgbIm );
    default           : {
        fprintf(stderr, " Unsupported DTYP in rgb2Yxy func.\n");
        return Mat();
        }
    }
}

Mat Yxy2rgb(const Mat& YxyIm){
    switch(YxyIm.getDatType()){
    case DTYP::DOUBLE : return _conv_Yxy2rgb<double>(YxyIm );
    case DTYP::FLOAT  : return _conv_Yxy2rgb<float >(YxyIm );
    case DTYP::INT    : return _conv_Yxy2rgb<int32 >(YxyIm );
    case DTYP::UCHAR  : return _conv_Yxy2rgb<uchar >(YxyIm );
    default           : {
        fprintf(stderr, " Unsupported DTYP in Yxy2rgb func.\n");
        return Mat();
        }
    }
}


Mat histoPmf(const Mat& src, const int32 bins, const double step){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _histoPmf<double>(src, bins, step);
    case DTYP::FLOAT  : return _histoPmf<float >(src, bins, step);
    case DTYP::INT    : return _histoPmf<int32 >(src, bins, step);
    case DTYP::UCHAR  : return _histoPmf<uchar >(src, bins, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in histoPmf func.\n");
        return Mat();
        }
    }
}

Mat histoCmf(const Mat& src, const int32 bins, const double step){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"histoCmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _histoCmf<double>(src, bins, step);
    case DTYP::FLOAT  : return _histoCmf<float >(src, bins, step);
    case DTYP::INT    : return _histoCmf<int32 >(src, bins, step);
    case DTYP::UCHAR  : return _histoCmf<uchar >(src, bins, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in histoCmf func.\n");
        return Mat();
        }
    }
}


Mat clip_HistoPmf(const Mat& src,const int32 clipVal,const int32 bins, const int32 step){
    uint32 ch  = src.getChannel();
    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoPmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"clip_histoPmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"clip_histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( step < 1) {
        fprintf(stderr,"clip_histoPmf : 'step' should be larger than or equal 1 \n");
        return Mat();
    }

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _clip_HistoPmf<double>(src, clipVal, bins, step);
    case DTYP::FLOAT  : return _clip_HistoPmf<float >(src, clipVal, bins, step);
    case DTYP::INT    : return _clip_HistoPmf<int32 >(src, clipVal, bins, step);
    case DTYP::UCHAR  : return _clip_HistoPmf<uchar >(src, clipVal, bins, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in clip_HistoPmf func.\n");
        return Mat();
        }
    }
}

Mat clip_HistoCmf(const Mat& src,const int32 clipVal,const int32 bins, const int32 step){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stderr,"clip_histoCmf : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stderr,"clip_histoCmf : src is not 1 channel matrix\n");
        return Mat();
    }else if( bins < 1) {
        fprintf(stderr,"clip_histoPmf : 'bins' should be larger than or equal 1 \n");
        return Mat();
    }else if( step < 1) {
        fprintf(stderr,"clip_histoPmf : 'step' should be larger than or equal 1 \n");
        return Mat();
    }

    switch(src.getDatType()){
    case DTYP::DOUBLE : return _clip_HistoCmf<double>(src, clipVal, bins, step);
    case DTYP::FLOAT  : return _clip_HistoCmf<float >(src, clipVal, bins, step);
    case DTYP::INT    : return _clip_HistoCmf<int32 >(src, clipVal, bins, step);
    case DTYP::UCHAR  : return _clip_HistoCmf<uchar >(src, clipVal, bins, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in clip_HistoCmf func.\n");
        return Mat();
        }
    }
}

Mat clip_HistoEqual(const Mat& src,const Mat& histCmf, const int32 step){
    uint32 ch  = src.getChannel();

    if( src.isEmpty() ){
        fprintf(stdout,"histoEqual : src argument is empty matrix\n");
        return Mat();
    }else if( ch != 1) {
        fprintf(stdout,"histoEqual : src is not 1 channel matrix\n");
        return Mat();
    }else if( histCmf.isEmpty()){
        fprintf(stdout,"histCmf is empty \n");
        return Mat();
    }
    switch(src.getDatType()){
    case DTYP::DOUBLE : return _clip_HistoEqual<double>(src, histCmf, step);
    case DTYP::FLOAT  : return _clip_HistoEqual<float >(src, histCmf, step);
    case DTYP::INT    : return _clip_HistoEqual<int32 >(src, histCmf, step);
    case DTYP::UCHAR  : return _clip_HistoEqual<uchar >(src, histCmf, step);
    default           : {
        fprintf(stderr, " Unsupported DTYP in clip_HistoEqual func.\n");
        return Mat();
        }
    }
}


Mat gaussMaskGen (const double sigma, const double factor, const uint32 ch ){
    // mask size sigma*factor*2+1
    uint32 hp = static_cast<uint32>( sigma*factor);
    uint32 sz = static_cast<uint32>( (hp<<1) +1 );

    Mat mask = Mat::zeros(sz, sz, ch, DTYP::DOUBLE);
    uint32 y, x, cc;
    int32 yp, xp;
    for( cc =0 ; cc < ch ; ++cc ){
        for (y=0, yp=-static_cast<int32>(hp) ; y < sz ; ++y, ++yp){
            for(x=0, xp=-static_cast<int32>(hp) ; x < sz; ++x, ++xp)
                mask.at<double>(y, x, cc) = exp(-((xp*xp + yp*yp)/(2*sigma*sigma)));
        }
    }
    Mat sum = mask.sum();
    for( cc =0 ; cc < ch ; ++cc ){
        for (y=0, yp=-static_cast<int32>(hp) ; y < sz ; ++y, ++yp){
            for(x=0, xp=-static_cast<int32>(hp) ; x < sz; ++x, ++xp)
                mask.at<double>(y, x, cc) /= sum.at<double>(cc);
        }
    }

    return mask;
}

Mat boxMaskGen( const uint32 sz, const uint32 ch){
    Mat mask = Mat::ones(sz,sz,ch,DTYP::DOUBLE);
    mask /= (sz*sz);
    return mask;
}

Mat nakaSigTonemap( Mat& src, Mat& localmean, const double globalmean, const double Imax){
    DTYP srcDtype  = src.getDatType();
    DTYP lmeanDtype = localmean.getDatType();

    if(srcDtype != lmeanDtype){
        fprintf(stderr,"the data types between src and localmean aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _nakaSigTm<double>( src, localmean, globalmean, Imax );
        break;
    case DTYP::FLOAT  :
        A   = _nakaSigTm<float >( src, localmean, globalmean, Imax);
        break;
    case DTYP::INT    :
        A   = _nakaSigTm<int32 >( src, localmean, globalmean, Imax);
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSigtonemap\n ");
    }

    return A;
}

Mat nakaSig3MeanTonemap( Mat& src, Mat& s_localmean, Mat& l_localmean, const double globalmean, const double Imax){
    DTYP srcDtype  = src.getDatType();
    DTYP smeanDtype = s_localmean.getDatType();
    DTYP lmeanDtype = l_localmean.getDatType();

    if(srcDtype != smeanDtype){
        fprintf(stderr,"the data types between src and s_localmean aren't the same!\n ");
        return Mat();
    }else if( srcDtype != lmeanDtype){
        fprintf(stderr,"the data types between src and l_localmean aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _nakaSig3MeanTm<double>( src, s_localmean, l_localmean, globalmean, Imax );
        break;
    case DTYP::FLOAT  :
        A   = _nakaSig3MeanTm<float >( src, s_localmean, l_localmean, globalmean, Imax);
        break;
    case DTYP::INT    :
        A   = _nakaSig3MeanTm<int32 >( src, s_localmean, l_localmean, globalmean, Imax);
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSig3MeanTonemap\n ");
    }

    return A;
}

Mat logRetinexTonemap( Mat& src, Mat& surround){
    DTYP srcDtype  = src.getDatType();
    DTYP surrDtype = surround.getDatType();

    if(srcDtype != surrDtype){
        fprintf(stderr,"the data types between src and surround aren't the same!\n ");
        return Mat();
    }

    Mat A;
    switch( srcDtype ){
    case DTYP::DOUBLE :
        A   = _logRetinexTm<double>( src, surround );
        break;
    case DTYP::FLOAT  :
        A   = _logRetinexTm<float >( src, surround );
        break;
    case DTYP::INT    :
        A   = _logRetinexTm<int32 >( src, surround );
        break;
    default:
        fprintf(stderr,"Unsupproted data type in nakaSigtonemap\n ");
    }

    return A;
}

Mat gamma( const Mat& src, const double gmvalue){
    DTYP srcDtype = src.getDatType();

    switch( srcDtype){
    case DTYP::DOUBLE : return _gamma<double>(src, gmvalue);
    case DTYP::FLOAT  : return _gamma<float >(src, gmvalue);
    case DTYP::INT    : return _gamma<int32 >(src, gmvalue);
    case DTYP::UCHAR  : return _gamma<uchar >(src, gmvalue);
    default: fprintf(stderr, "Unsuppretd data type in gamma func.\n"); return Mat();
    }
}

void fft(){
    uint32 k,dk, i;

    uint32 n= 3;
    uint32 n_minus_1 = n -1;
    uint32 t, t1, t2;

    uint32 mask;
    double *dat_r = new double[(1<<n)];
    double *dat_i = new double[(1<<n)];
    double dattmp;
    uint32 tNum = 1 << n;
    FILE *f1 = fopen("../a.txt","w");
    FILE *f2 = fopen("../b.txt","w");
    for(i=0; i < tNum; ++i){
        dat_r[i] = i;
        dat_i[i] = 0.;
    }

    fprintf(stdout,"Suffle\n");
    for( i=0 ; i < tNum ; i++){ // increaing index upto a half of 2^n
        t=0;
        // bit reverse of index
        for(k= 0; k <= n/2 ; k++){
            dk   = k << 1;    // k*2
            mask = 1 << k;
            t2   = i & mask;
            t1   = t2 << (n_minus_1-dk);
            t   |= t1;

            mask = 1 << (n_minus_1-k);
            t2   = i & mask;
            t1   = t2 >> (n_minus_1-dk);
            t   |= t1;
        }
        fprintf(f1,"%3d: %3d(%X)\n",i,t,t);
        if( t <= i) continue;
        dattmp   = dat_r[i];
        dat_r[i] = dat_r[t];
        dat_r[t] = dattmp;
    }
    for( i =0; i < tNum ; ++i){
        fprintf(f2,"%3d: %.4f\n",i,dat_r[i]);
    }

    delete [] dat_r;
    delete [] dat_i;
}

void fft_butterfly(double ra[], double ia[], uint32 step, uint32 len){
//    double f;
//    f = 2*3.141592;
//    ra[i]*cos(f/step)
}

} // end of imgproc namespace
} // end of jmat namespace

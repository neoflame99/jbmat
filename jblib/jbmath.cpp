#include "jbmath.h"
#include <stdio.h>
#include <float.h>
//#include <iostream>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

namespace jmat {


Mat mul(const Mat& mA,const Mat& mB){

    if(mA.getCol() != mB.getRow() ){
        fprintf(stderr," The number of columns of mA and the number of rows of mB are not match!\n");
        return Mat();
    }else if( mA.getChannel() != mB.getChannel()){
        fprintf(stderr," The numbers of channels of both mA and mB are not match!\n");
        return Mat();
    }

    uint32 aRow = mA.getRow();
    uint32 bCol = mB.getCol();
    uint32 ch   = mA.getChannel();

    Mat mO;
    DTYP maDtype = mA.getDatType();
    DTYP mbDtype = mB.getDatType();
    bool result_prod = false;

    if(maDtype == DTYP::CMPLX){
        switch(mbDtype){
        case DTYP::CMPLX :
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<cmplx ,cmplx ,cmplx>(mA, mB, mO); break;
        case DTYP::DOUBLE:
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<cmplx ,double,cmplx>(mA, mB, mO); break;
        case DTYP::FLOAT :
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<cmplx ,float ,cmplx>(mA, mB, mO); break;
        case DTYP::INT   :
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<cmplx ,int32 ,cmplx>(mA, mB, mO); break;
        case DTYP::UCHAR :
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<cmplx ,uchar ,cmplx>(mA, mB, mO); break;
        }
    }else if(maDtype == DTYP::DOUBLE){
        switch(mbDtype){
        case DTYP::CMPLX :
            mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<double, cmplx ,cmplx>(mA, mB, mO); break;
        case DTYP::DOUBLE:
            mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<double,double,double>(mA, mB, mO); break;
        case DTYP::FLOAT :
            mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<double,float ,double>(mA, mB, mO); break;
        case DTYP::INT   :
            mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<double,int32 ,double>(mA, mB, mO); break;
        case DTYP::UCHAR :
            mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<double,uchar ,double>(mA, mB, mO); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<float , cmplx ,cmplx>(mA, mB, mO); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<float ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _mul<float ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _mul<float ,int32   ,float >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _mul<float ,uchar ,float >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<int32 ,cmplx ,cmplx >(mA, mB, mO); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<int32 ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _mul<int32 ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _mul<int32 ,int32 ,int32 >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _mul<int32 ,uchar ,int32 >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, aRow, bCol, ch);
            result_prod = _mul<uchar ,cmplx ,cmplx >(mA, mB, mO); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _mul<uchar ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _mul<uchar ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _mul<uchar ,int32 ,int32 >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::UCHAR , aRow, bCol, ch);
            result_prod = _mul<uchar ,uchar ,uchar >(mA, mB, mO ); break;
        }
    }

    return (result_prod) ? mO : Mat();
}
Mat dot(const Mat& mA,const Mat& mB,uint32 dim){

    if(mA.getLength() != mB.getLength()){
        if(mA.getCol() != mB.getCol() )
            fprintf(stderr," The number of columns of mA and the number of rows of mB are not match!\n");
        else if( mA.getRow() != mB.getRow())
            fprintf(stderr," The numbers of rows of both mA and mB are not match!\n");
        else if( mA.getChannel() != mB.getChannel())
            fprintf(stderr," The numbers of channels of both mA and mB are not match!\n");

        return Mat();
    }


    uint32 aRow = mA.getRow();
    uint32 aCol = mA.getCol();
    uint32 bCol = mB.getCol();
    uint32 bRow = mB.getRow();
    uint32 ch   = mA.getChannel();

    uint32 oRow, oCol;
    if((aRow ==1 || aCol==1) && (bRow==1 || bCol==1)){
        /* MA and MB are row or column vectors, MO should be scalar or channel array*/
        oRow = 1;
        oCol = 1;
    }else if(dim==0){
        /* MA and MB are array, row-wise dot product such that MO is to be a column vector */
        oRow = aRow;
        oCol = 1;
    }else{
        /* MA and MB are array, column-wise dot product such that MO is to be a row vector */
        oRow = 1;
        oCol = aCol;
    }

    Mat mO;
    DTYP maDtype = mA.getDatType();
    DTYP mbDtype = mB.getDatType();
    bool result_prod = false;

    if(maDtype == DTYP::CMPLX){
        switch(mbDtype){
        case DTYP::CMPLX :
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<cmplx ,cmplx ,cmplx>(mA, mB, mO, dim); break;
        case DTYP::DOUBLE:
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<cmplx ,double,cmplx>(mA, mB, mO, dim); break;
        case DTYP::FLOAT :
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<cmplx ,float ,cmplx>(mA, mB, mO, dim); break;
        case DTYP::INT   :
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<cmplx ,int32 ,cmplx>(mA, mB, mO, dim); break;
        case DTYP::UCHAR :
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<cmplx ,uchar ,cmplx>(mA, mB, mO, dim); break;
        }
    }else if(maDtype == DTYP::DOUBLE){
        switch(mbDtype){
        case DTYP::CMPLX :
            mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<double, cmplx ,cmplx>(mA, mB, mO, dim); break;
        case DTYP::DOUBLE:
            mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<double,double,double>(mA, mB, mO, dim); break;
        case DTYP::FLOAT :
            mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<double,float ,double>(mA, mB, mO, dim); break;
        case DTYP::INT   :
            mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<double,int32 ,double>(mA, mB, mO, dim); break;
        case DTYP::UCHAR :
            mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<double,uchar ,double>(mA, mB, mO, dim); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<float , cmplx ,cmplx>(mA, mB, mO, dim); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<float ,double,double>(mA, mB, mO, dim); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , oRow, oCol, ch);
            result_prod = _dot_prod<float ,float ,float >(mA, mB, mO, dim); break;
        case DTYP::INT   : mO = Mat(DTYP::FLOAT , oRow, oCol, ch);
            result_prod = _dot_prod<float ,int32 ,float >(mA, mB, mO, dim); break;
        case DTYP::UCHAR : mO = Mat(DTYP::FLOAT , oRow, oCol, ch);
            result_prod = _dot_prod<float ,uchar ,float >(mA, mB, mO, dim); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<int32 ,cmplx ,cmplx >(mA, mB, mO, dim); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<int32 ,double,double>(mA, mB, mO, dim); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , oRow, oCol, ch);
            result_prod = _dot_prod<int32 ,float ,float >(mA, mB, mO, dim); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , oRow, oCol, ch);
            result_prod = _dot_prod<int32 ,int32 ,int32 >(mA, mB, mO, dim); break;
        case DTYP::UCHAR : mO = Mat(DTYP::INT   , oRow, oCol, ch);
            result_prod = _dot_prod<int32 ,uchar ,int32 >(mA, mB, mO, dim); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::CMPLX : mO = Mat(DTYP::CMPLX, oRow, oCol, ch);
            result_prod = _dot_prod<uchar ,cmplx ,cmplx >(mA, mB, mO, dim); break;
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, oRow, oCol, ch);
            result_prod = _dot_prod<uchar ,double,double>(mA, mB, mO, dim); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , oRow, oCol, ch);
            result_prod = _dot_prod<uchar ,float ,float >(mA, mB, mO, dim); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , oRow, oCol, ch);
            result_prod = _dot_prod<uchar ,int32 ,int32 >(mA, mB, mO, dim); break;
        case DTYP::UCHAR : mO = Mat(DTYP::UCHAR , oRow, oCol, ch);
            result_prod = _dot_prod<uchar ,uchar ,uchar >(mA, mB, mO, dim); break;
        }
    }

    return (result_prod) ? mO : Mat();
}
Mat triu(const Mat& mA)
{

    Mat utri = mA.copy();

    switch(utri.getDatType()){
    case DTYP::CMPLX  : _triu<cmplx >(utri); break;
    case DTYP::DOUBLE : _triu<double>(utri); break;
    case DTYP::FLOAT  : _triu<float >(utri); break;
    case DTYP::INT    : _triu<int32 >(utri); break;
    case DTYP::UCHAR  : _triu<uchar >(utri); break;
    }
    return utri;
}

Mat tril(const Mat& srcmat)
{

    Mat ltri = srcmat.copy();
    switch(ltri.getDatType()){
    case DTYP::CMPLX  : _tril<cmplx >(ltri); break;
    case DTYP::DOUBLE : _tril<double>(ltri); break;
    case DTYP::FLOAT  : _tril<float >(ltri); break;
    case DTYP::INT    : _tril<int32 >(ltri); break;
    case DTYP::UCHAR  : _tril<uchar >(ltri); break;
    }

    return ltri;
}

Mat augment(const Mat& srcmat)
{

    uint32 rows = srcmat.getRow();
    uint32 cols = srcmat.getCol();
    uint32 ch   = srcmat.getChannel();
    uint32 pivmax = (rows < cols)? rows : cols;
    uint32 augmentCols = cols + pivmax;

    DTYP datType = srcmat.getDatType();    
    Mat augm ;

    switch(datType){
    case DTYP::CMPLX  : augm = Mat(_augment<cmplx >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::DOUBLE : augm = Mat(_augment<double>(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::FLOAT  : augm = Mat(_augment<float >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::INT    : augm = Mat(_augment<int32 >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::UCHAR  : augm = Mat(_augment<uchar >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    }

    return augm;
}

Mat inverse(const Mat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    Mat invmat;
    switch(srcDtype){
    case DTYP::CMPLX  : invmat = _inverse<cmplx >(srcmat); break;
    case DTYP::DOUBLE : invmat = _inverse<double>(srcmat); break;
    case DTYP::FLOAT  : invmat = _inverse<float >(srcmat); break;
    default           : assert(false && "data type of srcmat into inverse is neither DOUBLE nor FLOAT");
                        fprintf(stderr,"data type of srcmat into inverse is neither DOUBLE nor FLOAT\n");
                        invmat = Mat();
    }

    return invmat;
}

Mat tranpose(const Mat& mA){

    if(mA.isEmpty()) return Mat();

    Mat t = mA.copy();
    t.transpose();
    return t;
}


Mat conv2d(const Mat& mA, const Mat& mB, std::string opt_conv, std::string opt_out ){

    if( mA.isEmpty()) {
        fprintf(stderr, "mA is empty\n ");
        return Mat();
    }else if( mB.isEmpty()){
        fprintf(stderr, "mB is empty\n ");
        return Mat();
    }else if( mA.getChannel() != mB.getChannel() ){
        fprintf(stderr, " channels of mA and mB are not the same! \n");
        return Mat();
    }

    uint32 ch    = mA.getChannel();
    bool fullout = false;    
    if( opt_out.compare("full")==0 )
        fullout = true;

    DTYP aDtype = mA.getDatType();
    DTYP bDtype = mB.getDatType();
    DTYP oDtype = (aDtype > bDtype) ? aDtype : bDtype;

    Mat mO;
    if(fullout){
        mO = Mat::zeros( mA.getRow() + mB.getRow() -1, mA.getCol() + mB.getCol() -1, ch, oDtype);
    }else{  // 'same'
        mO = Mat::zeros( mA.getRow(), mA.getCol(), mA.getChannel(), oDtype);
    }

    if(aDtype == DTYP::CMPLX ){
        switch(bDtype){
        case DTYP::CMPLX  : _conv2d<cmplx , cmplx , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::DOUBLE : _conv2d<double, double, cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<double, float , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<double, int32 , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<double, uchar , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::DOUBLE){
        switch(bDtype){
        case DTYP::CMPLX  : _conv2d<double, cmplx , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::DOUBLE : _conv2d<double, double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<double, float , double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<double, int32 , double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<double, uchar , double>(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::FLOAT){
        switch(bDtype){
        case DTYP::CMPLX  : _conv2d<float , cmplx , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::DOUBLE : _conv2d<float , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<float , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<float , int32 , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<float , uchar , float >(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::INT){
        switch(bDtype){
        case DTYP::CMPLX  : _conv2d<int32 , cmplx , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::DOUBLE : _conv2d<int32 , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<int32 , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<int32 , int32 , int32 >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<int32 , uchar , int32 >(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::UCHAR){
        switch(bDtype){
        case DTYP::CMPLX  : _conv2d<uchar , cmplx , cmplx >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::DOUBLE : _conv2d<uchar , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<uchar , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<uchar , int32 , int32 >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<uchar , uchar , uchar >(mA, mB, mO, fullout, opt_conv); break;
        }
    }

    return mO;
}

}

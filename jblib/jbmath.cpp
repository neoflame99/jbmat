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


Mat dot(const Mat& mA,const Mat& mB){

    if(mA.getCol() != mB.getRow() ){
        fprintf(stderr," The number of columns of mA and the number of rows of mB are not match!\n");
        return Mat();
    }else if( mA.getChannel() != mB.getChannel()){
        fprintf(stderr," The numbers of channels of both mA and mB are not match!\n");
        return Mat();
    }

    uint aRow = mA.getRow();
    uint bCol = mB.getCol();
    uint ch   = mA.getChannel();

    Mat mO;
    DTYP maDtype = mA.getDatType();
    DTYP mbDtype = mB.getDatType();
    bool result_prod = false;

    if(maDtype == DTYP::DOUBLE){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,double,double>(mA, mB, mO); break;
        case DTYP::FLOAT : mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,float ,double>(mA, mB, mO); break;
        case DTYP::INT   : mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,int   ,double>(mA, mB, mO); break;
        case DTYP::UCHAR : mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,uchar ,double>(mA, mB, mO); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<float ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,int   ,float >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,uchar ,float >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<int   ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,int   ,int   >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,uchar ,int   >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = Mat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = Mat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = Mat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,int   ,int   >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = Mat(DTYP::UCHAR , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,uchar ,uchar >(mA, mB, mO ); break;
        }
    }

    return (result_prod) ? mO : Mat();
}

Mat triu(const Mat& mA)
{

    Mat utri = mA.copy();

    switch(utri.getDatType()){
    case DTYP::DOUBLE : _triu<double>(utri); break;
    case DTYP::FLOAT  : _triu<float >(utri); break;
    case DTYP::INT    : _triu<int   >(utri); break;
    case DTYP::UCHAR  : _triu<uchar >(utri); break;
    }
    return utri;
}

Mat tril(const Mat& srcmat)
{

    Mat ltri = srcmat.copy();
    switch(ltri.getDatType()){
    case DTYP::DOUBLE : _tril<double>(ltri); break;
    case DTYP::FLOAT  : _tril<float >(ltri); break;
    case DTYP::INT    : _tril<int   >(ltri); break;
    case DTYP::UCHAR  : _tril<uchar >(ltri); break;
    }

    return ltri;
}

Mat augment(const Mat& srcmat)
{

    uint rows = srcmat.getRow();
    uint cols = srcmat.getCol();
    uint ch   = srcmat.getChannel();
    uint pivmax = (rows < cols)? rows : cols;
    uint augmentCols = cols + pivmax;

    DTYP datType = srcmat.getDatType();    
    Mat augm ;

    switch(datType){
    case DTYP::DOUBLE : augm = Mat(_augment<double>(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::FLOAT  : augm = Mat(_augment<float >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::INT    : augm = Mat(_augment<int   >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::UCHAR  : augm = Mat(_augment<uchar >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    }

    return augm;
}

Mat inverse(const Mat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    Mat invmat;
    switch(srcDtype){
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

    uint ch      = mA.getChannel();
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

    if(aDtype == DTYP::DOUBLE){
        switch(bDtype){
        case DTYP::DOUBLE : _conv2d<double, double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<double, float , double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<double, int   , double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<double, uchar , double>(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::FLOAT){
        switch(bDtype){
        case DTYP::DOUBLE : _conv2d<float , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<float , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<float , int   , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<float , uchar , float >(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::INT){
        switch(bDtype){
        case DTYP::DOUBLE : _conv2d<int   , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<int   , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<int   , int   , int   >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<int   , uchar , int   >(mA, mB, mO, fullout, opt_conv); break;
        }
    }else if(aDtype == DTYP::UCHAR){
        switch(bDtype){
        case DTYP::DOUBLE : _conv2d<uchar , double, double>(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::FLOAT  : _conv2d<uchar , float , float >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::INT    : _conv2d<uchar , int   , int   >(mA, mB, mO, fullout, opt_conv); break;
        case DTYP::UCHAR  : _conv2d<uchar , uchar , uchar >(mA, mB, mO, fullout, opt_conv); break;
        }
    }

    return mO;
}

}

#include "jbmath.h"
#include <stdio.h>
#include <float.h>
//#include <iostream>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

jbMat jbmath::dot(const jbMat& mA,const jbMat& mB){

    if(mA.getCol() != mB.getRow() ){
        fprintf(stderr," The number of columns of mA and the number of rows of mB are not match!\n");
        return jbMat();
    }else if( mA.getChannel() != mB.getChannel()){
        fprintf(stderr," The numbers of channels of both mA and mB are not match!\n");
        return jbMat();
    }

    uint aRow = mA.getRow();
    uint bCol = mB.getCol();
    uint ch   = mA.getChannel();

    jbMat mO;
    DTYP maDtype = mA.getDatType();
    DTYP mbDtype = mB.getDatType();
    bool result_prod = false;

    if(maDtype == DTYP::DOUBLE){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,double,double>(mA, mB, mO); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,float ,double>(mA, mB, mO); break;
        case DTYP::INT   : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,int   ,double>(mA, mB, mO); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,uchar ,double>(mA, mB, mO); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<float ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,int   ,float >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,uchar ,float >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<int   ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,int   ,int   >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,uchar ,int   >(mA, mB, mO ); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,double,double>(mA, mB, mO ); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,float ,float >(mA, mB, mO ); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,int   ,int   >(mA, mB, mO ); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::UCHAR , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,uchar ,uchar >(mA, mB, mO ); break;
        }
    }

    return (result_prod) ? mO : jbMat();
}

jbMat jbmath::triu(const jbMat& mA)
{

    jbMat utri = mA.copy();

    switch(utri.getDatType()){
    case DTYP::DOUBLE : _triu<double>(utri); break;
    case DTYP::FLOAT  : _triu<float >(utri); break;
    case DTYP::INT    : _triu<int   >(utri); break;
    case DTYP::UCHAR  : _triu<uchar >(utri); break;
    }
    return utri;
}

jbMat jbmath::tril(const jbMat& srcmat)
{

    jbMat ltri = srcmat.copy();
    switch(ltri.getDatType()){
    case DTYP::DOUBLE : _tril<double>(ltri); break;
    case DTYP::FLOAT  : _tril<float >(ltri); break;
    case DTYP::INT    : _tril<int   >(ltri); break;
    case DTYP::UCHAR  : _tril<uchar >(ltri); break;
    }

    return ltri;
}

jbMat jbmath::augment(const jbMat& srcmat)
{

    uint rows = srcmat.getRow();
    uint cols = srcmat.getCol();
    uint ch   = srcmat.getChannel();
    uint pivmax = (rows < cols)? rows : cols;
    uint augmentCols = cols + pivmax;

    DTYP datType = srcmat.getDatType();    
    jbMat augm ;

    switch(datType){
    case DTYP::DOUBLE : augm = jbMat(_augment<double>(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::FLOAT  : augm = jbMat(_augment<float >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::INT    : augm = jbMat(_augment<int   >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::UCHAR  : augm = jbMat(_augment<uchar >(srcmat, augmentCols), datType, rows, augmentCols, ch); break;
    }

    return augm;
}

jbMat jbmath::inverse(const jbMat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    jbMat invmat;
    switch(srcDtype){
    case DTYP::DOUBLE : invmat = _inverse<double>(srcmat); break;
    case DTYP::FLOAT  : invmat = _inverse<float >(srcmat); break;
    default           : assert(false && "data type of srcmat into inverse is neither DOUBLE nor FLOAT");
                        fprintf(stderr,"data type of srcmat into inverse is neither DOUBLE nor FLOAT\n");
                        invmat = jbMat();
    }

    return invmat;
}

jbMat jbmath::tranpose(const jbMat& mA){

    if(mA.isEmpty()) return jbMat();

    jbMat t = mA.copy();
    t.transpose();
    return t;
}


jbMat jbmath::conv2d(const jbMat& mA, const jbMat& mB, std::string opt_conv, std::string opt_out ){

    if( mA.isEmpty()) {
        fprintf(stderr, "mA is empty\n ");
        return jbMat();
    }else if( mB.isEmpty()){
        fprintf(stderr, "mB is empty\n ");
        return jbMat();
    }else if( mA.getChannel() != mB.getChannel() ){
        fprintf(stderr, " channels of mA and mB are not the same! \n");
        return jbMat();
    }

    uint ch      = mA.getChannel();
    bool fullout = false;    
    if( opt_out.compare("full")==0 )
        fullout = true;

    DTYP aDtype = mA.getDatType();
    DTYP bDtype = mB.getDatType();
    DTYP oDtype = (aDtype > bDtype) ? aDtype : bDtype;

    jbMat mO;
    if(fullout){
        mO = jbMat::zeros( mA.getRow() + mB.getRow() -1, mA.getCol() + mB.getCol() -1, ch, oDtype);
    }else{  // 'same'
        mO = jbMat::zeros( mA.getRow(), mA.getCol(), mA.getChannel(), oDtype);
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


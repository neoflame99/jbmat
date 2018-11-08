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
            result_prod = _dot_prod<double,double,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,float,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::INT   : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,int   ,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,uchar ,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<float ,double,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,float ,float >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::INT   : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,int   ,float >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,uchar ,float >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<int   ,double,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,float ,float >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,int   ,int   >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,uchar ,int   >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,double,double>(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,float ,float >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,int   ,int   >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::UCHAR , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,uchar ,uchar >(mA.getRawMat(), mB.getRawMat(), mO.getRawMat()); break;
        }
    }

    return (result_prod) ? mO : jbMat();
}

jbMat jbmath::triu(const jbMat& mA)
{

    jbMat utri = mA.copy();

    switch(utri.getDatType()){
    case DTYP::DOUBLE : _triu<double>(utri.getRawMat()); break;
    case DTYP::FLOAT  : _triu<float >(utri.getRawMat()); break;
    case DTYP::INT    : _triu<int   >(utri.getRawMat()); break;
    case DTYP::UCHAR  : _triu<uchar >(utri.getRawMat()); break;
    }
    return utri;
}

jbMat jbmath::tril(const jbMat& srcmat)
{

    jbMat ltri = srcmat.copy();
    switch(ltri.getDatType()){
    case DTYP::DOUBLE : _tril<double>(ltri.getRawMat()); break;
    case DTYP::FLOAT  : _tril<float >(ltri.getRawMat()); break;
    case DTYP::INT    : _tril<int   >(ltri.getRawMat()); break;
    case DTYP::UCHAR  : _tril<uchar >(ltri.getRawMat()); break;
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
    case DTYP::DOUBLE : augm = jbMat(_augment<double>(srcmat.getRawMat(), augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::FLOAT  : augm = jbMat(_augment<float >(srcmat.getRawMat(), augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::INT    : augm = jbMat(_augment<int   >(srcmat.getRawMat(), augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::UCHAR  : augm = jbMat(_augment<uchar >(srcmat.getRawMat(), augmentCols), datType, rows, augmentCols, ch); break;
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

jbMat jbmath::tranpose(const jbMat &mA){

    if(mA.isEmpty()) return jbMat();

    jbMat t = mA.copy();
    t.transpose();
    return t;
}

/*
jbMat jbMath::conv2d(const jbMat& mA, const jbMat& mB, std::string opt_conv, std::string opt_out ){

    if(mA.isEmpty()) {
        fprintf(stderr, "mA is empty\n ");
        return jbMat();
    }else if(mB.isEmpty()){
        fprintf(stderr, "mB is empty\n ");
        return jbMat();
    }else if( mA.getChannel() != mB.getChannel() ){
        fprintf(stderr, " channels of mA and mB are not the same! \n");
        return jbMat();
    }

    int tX , tY;
    int xdummy , ydummy;
    int mbHcol = mB.getCol()/2;
    int mbHrow = mB.getRow()/2;
    int ch = mA.getChannel();
    bool fullout = false;
    jbMat tmpO;

    if( opt_out.compare("full")==0 )
        fullout = true;

    if(fullout){
        xdummy = mB.getCol() -1;
        ydummy = mB.getRow() -1;
        tX = mA.getCol() + xdummy*2;
        tY = mA.getRow() + ydummy*2;

    }else{
        tX = mA.getCol() + mB.getCol() -1;
        tY = mA.getRow() + mB.getRow() -1;
    }

    jbMat tmpA(tY, tX, ch);
    if(fullout){
        tmpO = jbMat(mA.getRow() + mB.getRow() -1, mA.getCol() + mB.getCol() -1, ch);
    }else{  // 'same'
        tmpO = mA.copy();
    }

    if( fullout ){
        for(int ich=0; ich < ch; ich++){
            for(int y =0; y < tY; y++){
                for(int x=0; x < tX; x++){
                    if( (y >= xdummy && y < tY-ydummy) && ( x >= xdummy && x < tX-ydummy) )
                        tmpA(y,x,ich) = mA(y-ydummy,x-xdummy,ich);
                    else
                        tmpA(y,x,ich) = 0;
                }
            }
        }
    }else{
        for(int ich=0; ich < ch; ich++){
            for(int y =0; y < tY; y++){
                for(int x=0; x < tX; x++){
                    if(opt_conv.substr(0,4).compare("symm")==0){

                        if( y < mbHrow && x < mbHcol )         // left top corner
                            tmpA(y,x,ich) = mA(mbHrow-y-1, mbHcol-x-1, ich);
                        else if( y < mbHrow && x >= tX-mbHcol) // right top corner
                            tmpA(y,x,ich) = mA(mbHrow-y-1, mA.getCol()-x-mbHcol+tX-1, ich);
                        else if( y < mbHrow)                   // top
                            tmpA(y,x,ich) = mA(mbHrow-y-1, x-mbHcol, ich);
                        else if( y >= tY-mbHrow && x < mbHcol) // left bottom corner
                            tmpA(y,x,ich) = mA(mA.getRow()-y-mbHrow+tY-1, mbHcol-x-1,ich);
                        else if( y >= tY-mbHrow && x >= tX-mbHcol) // right bottom corner
                            tmpA(y,x,ich) = mA(mA.getRow()-y-mbHrow+tY-1, mA.getCol()-x-mbHcol+tX-1, ich);
                        else if( y >= tY-mbHrow)               // bottom
                            tmpA(y,x,ich) = mA(mA.getRow()-y-mbHrow+tY-1, x-mbHcol, ich);
                        else if( x < mbHcol )                  // left
                            tmpA(y,x,ich) = mA(y-mbHrow, mbHcol-x-1, ich);
                        else if( x >= tX - mbHcol)             // right
                            tmpA(y,x,ich) = mA(y-mbHrow, mA.getCol()-x-mbHcol+tX-1, ich);
                        else                                   // main
                            tmpA(y,x,ich) = mA(y-mbHrow,x-mbHcol,ich);
                    }else if(opt_conv.substr(0,4).compare("circ")==0){

                        if( y < mbHrow && x < mbHcol )         // left top corner
                            tmpA(y,x,ich) = mA(mA.getRow()-mbHrow+y, mA.getCol()-mbHcol+x, ich);
                        else if( y < mbHrow && x >= tX-mbHcol) // right top corner
                            tmpA(y,x,ich) = mA(mA.getRow()-mbHrow+y, mbHcol-tX+x, ich);
                        else if( y < mbHrow)                   // top
                            tmpA(y,x,ich) = mA(mA.getRow()-mbHrow+y, x-mbHcol, ich);
                        else if( y >= tY-mbHrow && x < mbHcol) // left bottom corner
                            tmpA(y,x,ich) = mA(mbHrow-tY+y, mA.getCol()-mbHcol+x,ich);
                        else if( y >= tY-mbHrow && x >= tX-mbHcol) // right bottom corner
                            tmpA(y,x,ich) = mA(mbHrow-tY+y, mbHcol-tX+x, ich);
                        else if( y >= tY-mbHrow)               // bottom
                            tmpA(y,x,ich) = mA(mbHrow-tY+y, x-mbHcol, ich);
                        else if( x < mbHcol )                  // left
                            tmpA(y,x,ich) = mA(y-mbHrow, mA.getCol()-mbHcol+x, ich);
                        else if( x >= tX - mbHcol)             // right
                            tmpA(y,x,ich) = mA(y-mbHrow, mbHcol-tX+x, ich);
                        else                                   // main
                            tmpA(y,x,ich) = mA(y-mbHrow,x-mbHcol,ich);
                    }else{

                        if( (y >= mbHrow && y < tY-mbHrow) && ( x >= mbHcol && x < tX-mbHcol) )
                            tmpA(y,x,ich) = mA(y-mbHrow,x-mbHcol,ich);
                        else
                            tmpA(y,x,ich) = 0;
                    }
                }
            }
        }
    }
   // tmpA.printMat();

    // convolution
    double sum;

    if(fullout){
        for(int ich=0; ich < ch ; ich++){
            for(int y=ydummy; y < tY; y++){
                for(int x=xdummy; x < tX; x++){
                    sum = 0.0;
                    for(int m= mB.getRow()-1; m >= 0; m--){
                        for(int n= mB.getCol()-1; n >= 0; n--){
                            sum += (tmpA(y-m,x-n,ich)* mB(m,n,ich));
                            //a = tmpA(y-mbHrow+m,x-mbHcol+n,ich);
                            //b = mB(m,n,ich);
                            //sum += (a+b);
                        }
                    }
                    tmpO(y-ydummy,x-xdummy,ich) = sum;
                }
            }
        }
    }else{
        for(int ich=0; ich < ch ; ich++){
            for(int y=mbHrow; y < tY-mbHrow; y++){
                for(int x=mbHcol; x < tX-mbHcol; x++){
                    sum = 0.0;
                    for(int m= mB.getRow()-1; m >= 0; m--){
                        for(int n= mB.getCol()-1; n >= 0; n--){
                            //sum += (tmpA(y-m,x-n,ich)* mB(-m+mbHrow,-n+mbHcol,ich));
                            sum += (tmpA(y+mbHrow-m,x+mbHcol-n,ich)* mB(m,n,ich));
                            //a = tmpA(y-mbHrow+m,x-mbHcol+n,ich);
                            //b = mB(m,n,ich);
                            //sum += (a+b);
                        }
                    }
                    tmpO(y-mbHrow,x-mbHcol,ich) = sum;
                }
            }
        }
    }

    return tmpO;
}
*/

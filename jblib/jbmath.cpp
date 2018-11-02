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
    uint aCol = mA.getCol();
    uint bCol = mB.getCol();
    uint ch   = mA.getChannel();

    jbMat mO;
    DTYP maDtype = mA.getDatType();
    DTYP mbDtype = mB.getDatType();
    bool result_prod = false;

    if(maDtype == DTYP::DOUBLE){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,double,double>(mA.getDataPtr<double>(), mB.getDataPtr<double>(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,float ,double>(mA.getDataPtr<double>(), mB.getDataPtr<float >(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::INT   : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,int   ,double>(mA.getDataPtr<double>(), mB.getDataPtr<int   >(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<double,uchar ,double>(mA.getDataPtr<double>(), mB.getDataPtr<uchar >(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        }
    }else if(maDtype == DTYP::FLOAT){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<float ,double,double>(mA.getDataPtr<float >(), mB.getDataPtr<double>(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,float ,float >(mA.getDataPtr<float >(), mB.getDataPtr<float >(), mO.getDataPtr<float >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::INT   : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,int   ,float >(mA.getDataPtr<float >(), mB.getDataPtr<int   >(), mO.getDataPtr<float >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<float ,uchar ,float >(mA.getDataPtr<float >(), mB.getDataPtr<uchar >(), mO.getDataPtr<float >(), aRow, aCol, ch, aCol, bCol, ch); break;
        }
    }else if(maDtype == DTYP::INT ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<int   ,double,double>(mA.getDataPtr<int   >(), mB.getDataPtr<double>(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,float ,float >(mA.getDataPtr<int   >(), mB.getDataPtr<float >(), mO.getDataPtr<float >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,int   ,int   >(mA.getDataPtr<int   >(), mB.getDataPtr<int   >(), mO.getDataPtr<int   >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<int   ,uchar ,int   >(mA.getDataPtr<int   >(), mB.getDataPtr<uchar >(), mO.getDataPtr<int   >(), aRow, aCol, ch, aCol, bCol, ch); break;
        }
    }else if(maDtype == DTYP::UCHAR ){
        switch(mbDtype){
        case DTYP::DOUBLE: mO = jbMat(DTYP::DOUBLE, aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,double,double>(mA.getDataPtr<uchar >(), mB.getDataPtr<double>(), mO.getDataPtr<double>(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::FLOAT : mO = jbMat(DTYP::FLOAT , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,float ,float >(mA.getDataPtr<uchar >(), mB.getDataPtr<float >(), mO.getDataPtr<float >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::INT   : mO = jbMat(DTYP::INT   , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,int   ,int   >(mA.getDataPtr<uchar >(), mB.getDataPtr<int   >(), mO.getDataPtr<int   >(), aRow, aCol, ch, aCol, bCol, ch); break;
        case DTYP::UCHAR : mO = jbMat(DTYP::UCHAR , aRow, bCol, ch);
            result_prod = _dot_prod<uchar ,uchar ,uchar >(mA.getDataPtr<uchar >(), mB.getDataPtr<uchar >(), mO.getDataPtr<uchar >(), aRow, aCol, ch, aCol, bCol, ch); break;
        }
    }

    return (result_prod) ? mO : jbMat();
}

jbMat jbmath::triu(const jbMat& mA)
{

    uint rows = mA.getRow();
    uint cols = mA.getCol();
    uint ch   = mA.getChannel();

    jbMat utri = mA.copy();

    switch(utri.getDatType()){
    case DTYP::DOUBLE : _triu(utri.getDataPtr<double>(), rows, cols, ch ); break;
    case DTYP::FLOAT  : _triu(utri.getDataPtr<float >(), rows, cols, ch ); break;
    case DTYP::INT    : _triu(utri.getDataPtr<int   >(), rows, cols, ch ); break;
    case DTYP::UCHAR  : _triu(utri.getDataPtr<uchar >(), rows, cols, ch ); break;
    }

    return utri;
}

jbMat jbmath::tril(const jbMat& srcmat)
{
    uint rows = srcmat.getRow();
    uint cols = srcmat.getCol();
    uint ch   = srcmat.getChannel();

    jbMat ltri = srcmat.copy();
    switch(ltri.getDatType()){
    case DTYP::DOUBLE : _tril(ltri.getDataPtr<double>(), rows, cols, ch ); break;
    case DTYP::FLOAT  : _tril(ltri.getDataPtr<float >(), rows, cols, ch ); break;
    case DTYP::INT    : _tril(ltri.getDataPtr<int   >(), rows, cols, ch ); break;
    case DTYP::UCHAR  : _tril(ltri.getDataPtr<uchar >(), rows, cols, ch ); break;
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
    //jbMat augm(datType, rows, augmentCols, ch);
    jbMat augm ;

    switch(datType){
    case DTYP::DOUBLE : augm = jbMat(_augment(srcmat.getDataPtr<double>(), rows, cols, ch, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::FLOAT  : augm = jbMat(_augment(srcmat.getDataPtr<float >(), rows, cols, ch, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::INT    : augm = jbMat(_augment(srcmat.getDataPtr<int   >(), rows, cols, ch, augmentCols), datType, rows, augmentCols, ch); break;
    case DTYP::UCHAR  : augm = jbMat(_augment(srcmat.getDataPtr<uchar >(), rows, cols, ch, augmentCols), datType, rows, augmentCols, ch); break;
    }

    return augm;
}
/*
jbMat jbMath::inverse(const jbMat& srcmat){
    int rows = srcmat.getRow();
    int cols = srcmat.getCol();
    if(rows != cols) {
        std::cout << "The inverse matrix cannot be computed because source matrix is not square!";
        return jbMat();
    }
    int pvmax= rows;

    jbMat mataug = augmentMatrix(srcmat);
    jbMat utri = triu(mataug);
    jbMat ltri = tril(utri);
    int ltri_col = ltri.getCol();
    int pv,j, pvr;

    double pivot;
    double* mat = ltri.getMat().get();
    for( pv=0 ; pv<pvmax ; pv++){
        pvr   = pv * ltri_col;
        pivot = mat[pvr+pv];
        if(pivot==0.0) continue;
        for(j=0 ; j < ltri_col ; j++)
            mat[pvr+j] /= pivot;
    }

    jbMat invmat(rows,cols);
    int cr;
    for(int i=0; i < rows ; i++){
        pvr = i * ltri_col + cols;
        cr  = i * cols;
        for( j=0 ; j < cols ; j++){
            invmat[cr+j] = ltri[pvr+j];
        }
    }

#ifdef _JB_DEBUG_
    std::cout << "Inverse Result: ";
    printMat(invmat);
#endif
    return invmat;
}

jbMat jbMath::tranpose(const jbMat &mA){

    if(mA.isEmpty()) return jbMat();

    jbMat t = mA;
    t.transpose();
#ifdef _JB_DEBUG_
    std::cout << "Transpose : ";
    printMat(t);
#endif
    return t;
}

void jbMath::printMat(const jbMat& Mat)
{
    const int bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    int i,j;
//    int rows = Mat.getRow();
    int cols = Mat.getCol();
    int ch   = Mat.getChannel();
    int len  = Mat.getLength();
    double *mat = Mat.getMat().get();
    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double = DBL_EPSILON ;

    double val;
    int k;
    for( k=0; k < ch; k++){
        fprintf(stdout,"channel: %d \n",k);
        for( i=0; i< len ; i+= cols*ch){
            snprintf(buf,bufsz,"[");
            for( j=0; j< cols*ch; j+=ch){
                val = mat[i+j+k];
                if( val >= neg_max_double && val <= pos_min_double)
                    val = 0.0;
                snprintf(tmp,bufsz," %.4f ",val);
                strncat(buf,tmp,bufsz);
            }

            strncat(buf,"]",1);
            fprintf(stdout,"%s\n",buf);
        }
    }
}


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

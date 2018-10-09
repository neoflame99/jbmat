#include "jbmath.h"
#include <stdio.h>
#include <float.h>
#include <iostream>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

jbMath::jbMath()
{
}

jbMath::~jbMath()
{
}

jbMat jbMath::mulMatrix(const jbMat& mA,const jbMat& mB){

    if(mA.getCol() != mB.getRow() ){
        fprintf(stderr," The number of columns of mA and the number of rows of mB are not match!\n");
        return jbMat();
    }else if( mA.getChannel() != mB.getChannel()){
        fprintf(stderr," The numbers of channels of both mA and mB are not match!\n");
        return jbMat();
    }

    int aRow = mA.getRow();
    int aCol = mA.getCol();
    int bCol = mB.getCol();
    int ch   = mA.getChannel();
    int i,j;
    double *matA;
    double *matB;
    double *matO;

    jbMat mO(aRow,bCol, ch );

    matA = mA.getMat().get();
    matB = mB.getMat().get();
    matO = mO.getMat().get();

    int k, lr, lc, rc, rra, rrb;
    double a,b,c;
    for(int m=0; m < ch; m++ ){
        for( i = 0 , lr =0, rra=0 ; i < aRow ; i++, lr +=bCol*ch, rra+=aCol*ch){
            for( j = 0, lc=0; j< bCol; j++, lc+= ch ){
                matO[lr+lc+m] = 0.0;
                c=0;
                for( k=0, rrb=0, rc =0; k < aCol; k++, rrb+=ch*bCol, rc += ch){
                    a = matA[rra+rc+m];
                    b = matB[rrb+lc+m];
                    c += a*b;
                    matO[lr+lc+m] += matA[rra+rc+m] * matB[rrb+lc+m];
                }
            }
        }
    }
#ifdef _JB_DEBUG_
    std::cout << " Matrix multimplication: ";
    printMat(mO);
#endif
    return mO;

}


jbMat jbMath::triu(const jbMat& mA)
{
    int pivtmax = (mA.getRow() < mA.getCol() ) ? mA.getRow() : mA.getCol();
    int rows = mA.getRow();
    int cols = mA.getCol();
    int ch   = mA.getChannel();

    jbMat utri = mA;
    double* utri_ma = utri.getMat().get();
    int   pv, i,j, cr, pvr;
    double fact;    

    for( int cc=0; cc < ch; cc++){
        for( pv=0; pv < pivtmax-1 ; pv++){
            //std::cout << "pv = " << pv << " ";
            pvr = pv* cols + cc;
            for(i = pv+1 ; i < rows ; i++){
                cr = i* cols + cc;
                fact = utri_ma[cr+pv] / utri_ma[pvr+pv];
                std::cout << "cr= " << cr <<", i = " << i << ", pv= " << pv << " pvr= " << pvr <<" fact = " << fact << " ";
                std::cout << "utri_ma[cr+pv] = " << utri_ma[cr+pv] <<" utri_ma[pvr+pv] = " << utri_ma[pvr+pv] << "\n";
                for(j=0 ; j<cols ; j++){
                    if(utri_ma[pvr+pv] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    utri_ma[cr+j] = (pv==j) ? 0 : utri_ma[cr+j] - utri_ma[pvr+j] * fact;
                }
            }
        }
    }

#ifdef _JB_DEBUG_
    std::cout << "Upper Triangular Result: ";
    printMat(utri);
#endif

    return utri;
}

jbMat jbMath::tril(const jbMat& srcmat)
{
    int rows = srcmat.getRow();
    int cols = srcmat.getCol();
    int ch   = srcmat.getChannel();
    int pivtmax = (rows < cols) ? rows : cols;

    jbMat ltri = srcmat;
    int   pv, i,j, cr, pvr;
    double fact;

    for(int cc=0; cc < ch; cc++){
        for(pv=pivtmax-1 ; pv>0 ; pv--){
            pvr = pv*cols*ch;
            for(i=pv-1 ; i >= 0; i--){
                cr = i* cols*ch;
                fact = ltri[cr+pv*ch] / ltri[pvr+pv*ch];
                //fprintf(stdout,"pv=%d, i=%d, fact=%f",pv,i,fact);
                for(j=0 ; j < cols; j++){
                    if(ltri[pvr+pv*ch] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    ltri[cr+j*ch] = (pv==j) ? 0 : ltri[cr+j*ch] - ltri[pvr+j*ch] * fact;
                }
            }
        }
    }

#ifdef _JB_DEBUG_
    std::cout << "Lower Triangular Result: ";
    printMat(ltri);
#endif

    return ltri;
}

jbMat jbMath::augmentMatrix(const jbMat& srcmat)
{
    int i,j;
    int rows = srcmat.getRow();
    int cols = srcmat.getCol();
    int ch   = srcmat.getChannel();
    int pivmax = (rows < cols)? rows : cols;
    int augmentCols = cols + pivmax;

    jbMat augm(rows,augmentCols,ch);
    int cr,scr;

    //--augmenting
    for( i=0 ; i<rows ; i++){
        scr = i*cols;
        cr  = i*augmentCols;
        for( j=0; j < augmentCols ; j++){
            if(j<cols)
                augm[cr+j] = srcmat[scr+j] ;
            else if(i==j-cols)
                augm[cr+j] =1.0;
            else
                augm[cr+j] =0.0;
        }
    }

#ifdef _JB_DEBUG_
    std::cout << "Augmenting Result: ";
    printMat(augm);
#endif

    return augm;
}

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

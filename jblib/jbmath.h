#ifndef JBMATH_H
#define JBMATH_H
#include <stddef.h>
#include "jbMat.h"
#include <string>
//#include <memory>
#include <iostream>
#include <assert.h>

namespace jmat {

    Mat mul(const Mat& mA,const Mat& mB);
    Mat triu(const Mat& mA);
    Mat tril(const Mat& mA);
    Mat augment(const Mat& mA);
    Mat inverse(const Mat& mA);
    Mat tranpose(const Mat& mA);

    Mat conv2d(const Mat& mA, const Mat& mB, std::string opt_conv="" , std::string opt_out="");

    template <typename _Ta, typename _Tb, typename _To > bool _mul(const Mat& rMA,const Mat& rMB, Mat& rMO);
    template <typename _Ta, typename _Tb, typename _To > bool _dot_prod(const Mat& rMA,const Mat& rMB, Mat& rMO,uint32 dim);
    template <typename _T > void _triu(Mat& utri_mat);
    template <typename _T > void _tril(Mat& ltri_mat);
    template <typename _T > std::shared_ptr<uchar> _augment(const Mat&, const uint32 augCols);
    template <typename _T> Mat _inverse(const Mat& srcmat);
    template <typename _Ta, typename _Tb, typename _To> void _conv2d(const Mat& mA, const Mat& mB, Mat& mO, const bool fullout,const std::string& opt_conv);


template <typename _Ta, typename _Tb, typename _To>
bool _mul(const Mat& rMA,const Mat& rMB, Mat& rMO){
    uint32 Ar  = rMA.getRow();
    uint32 Ac  = rMA.getCol();
    uint32 Ach = rMA.getChannel();
    uint32 Br  = rMB.getRow();
    uint32 Bc  = rMB.getCol();
    uint32 Bch = rMB.getChannel();
    uint32 Or  = rMO.getRow();
    uint32 Oc  = rMO.getCol();
    uint32 Och = rMO.getChannel();

    _Ta* MA  = rMA.getDataPtr<_Ta>();
    _Tb* MB  = rMB.getDataPtr<_Tb>();
    _To* MO  = rMO.getDataPtr<_To>();

    assert(Ac==Br);
    assert(Ach==Bch);
    assert(Or==Ar);
    assert(Oc==Bc);
    assert(Ach==Och);
    assert(MA!=nullptr);
    assert(MB!=nullptr);
    assert(MO!=nullptr);
/*
    if( Ac != Br || Ach != Bch ){
        fprintf(stderr, "sizes of ma and mb into _mul are not match!\n");
        return false;
    }else if( rMO.getRow() != Ar || rMO.getCol() != Bc || rMO.getChannel() != Ach){
        fprintf(stderr, "sizes of mo is not enough!\n");
        return false;
    }else if( MA == nullptr || MB == nullptr || MO == nullptr){
        fprintf(stderr, "one or more of ma, mb or mo into _mul are NULL!\n");
        return false;
    }
*/
    uint32 i,j;
    uint32 k, lr, lc, rc, rra, rrb;
    _Ta av;
    _Tb bv;
    _To cv;

    uint32 mo_chidx, ma_chidx, mb_chidx;
    uint32 Arc = Ar * Ac;
    uint32 Brc = Br * Bc;
    uint32 Orc = Ar * Bc;
    uint32 m;
    mo_chidx=0; ma_chidx=0; mb_chidx=0;
    for( m=0; m < Ach ; m++ ){
        for( i = 0 , lr =0, rra=0 ; i < Ar ; i++, lr += Bc, rra+= Ac){
            for( j = 0, lc=0; j< Bc ; j++, lc++ ){
                MO[lr+lc+mo_chidx] = 0;
                cv=0;
                for( k=0, rrb=0, rc =0; k < Ac; k++, rrb += Bc, rc++){
                    av = MA[rra+rc+ma_chidx];
                    bv = MB[rrb+lc+mb_chidx];
                    cv = av*bv;
                    MO[lr+lc+mo_chidx] += cv;
                    //MO[lr+lc+m] += MA[rra+rc+m] * MB[rrb+lc+m];
                }
            }
        }
        mo_chidx += Orc; ma_chidx += Arc; mb_chidx += Brc;
    }
    return true;
}
template <typename _Ta, typename _Tb, typename _To>
bool _dot_prod(const Mat& rMA,const Mat& rMB, Mat& rMO, uint32 dim){
    uint32 Ar  = rMA.getRow();
    uint32 Ac  = rMA.getCol();
    uint32 Ach = rMA.getChannel();
    uint32 Br  = rMB.getRow();
    uint32 Bc  = rMB.getCol();
    uint32 Bch = rMB.getChannel();
    uint32 Or  = rMO.getRow();
    uint32 Oc  = rMO.getCol();
    uint32 Och = rMO.getChannel();

    _Ta* MA  = rMA.getDataPtr<_Ta>();
    _Tb* MB  = rMB.getDataPtr<_Tb>();
    _To* MO  = rMO.getDataPtr<_To>();

    uint32 i,j;
    uint32 Arc = Ar * Ac;
    uint32 Brc = Br * Bc;
    uint32 m;
    _To sum;
    assert(Ach == Bch);
    assert(Ach == Och);
    if ( (Ar==1 || Ac==1) && (Br==1 || Bc==1)){
        /* MA and MB are row or column vectors, MO should be scalar or channel array*/
        assert(Arc == Brc);
        assert(Or  == 1  );
        assert(Oc  == 1  );
        for (m = 0; m < Ach; m++){
            sum = 0;
            for(i = 0; i < Arc ; i++)
                sum += MA[i] * MB[i];

            MO[m] = sum;
        }
    }else if(dim==0){
        /* MA and MB are array, row-wise dot product such that MO is to be a column vector */
        assert(Ar == Br);
        assert(Ac == Bc);
        assert(Oc == 1 );
        assert(Or == Ar);
        for(m = 0 ; m < Ach; m++){
            for(i=0; i < Ar; i++){
                sum = 0;
                for(j=0; j < Ac; j++)
                    sum += MA(i,j,m) * MB(i,j,m);
                MO(i,0,m) = sum;
            }
        }
    }else {
        /* MA and MB are array, column-wise dot product such that MO is to be a row vector */
        assert(Ar == Br);
        assert(Ac == Bc);
        assert(Oc == Ac);
        assert(Or == 1 );
        for(m = 0 ; m < Ach; m++){
            for(i=0; i < Ac; i++){
                sum = 0;
                for(j=0; j < Ar; j++){
                    sum += MA(j,i,m) * MB(j,i,m);
                }
                MO(0,i,m) = sum;
            }
        }
    }

    return true;
}
template <typename _T> void _triu( Mat& utri_mat){
    uint32 rows = utri_mat.getRow();
    uint32 cols = utri_mat.getCol();
    uint32 ch   = utri_mat.getChannel();

    uint32 pivtmax = (rows < cols ) ? rows : cols;
    uint32 pv, i,j, cr, pvr;
    uint32 rcstep = rows*cols;
    uint32 tlen   = rcstep*ch;
    double fact;
    _T* mat = utri_mat.getDataPtr<_T>();
    // do triu
    for( uint32 cc=0; cc < tlen ; cc+=rcstep){
        for( pv=0; pv < pivtmax-1 ; pv++){
            //std::cout << "pv = " << pv << " ";
            pvr = pv* cols + cc;
            for( i = pv+1 ; i < rows ; i++){
                cr = i* cols + cc;
                fact = mat[cr+pv] / mat[pvr+pv];
                //std::cout << "cr= " << cr <<", i = " << i << ", pv= " << pv << " pvr= " << pvr <<" fact = " << fact << " ";
                //std::cout << "utri_ma[cr+pv] = " << utri_ma[cr+pv] <<" utri_ma[pvr+pv] = " << utri_ma[pvr+pv] << "\n";
                for(j=0 ; j<cols ; j++){
                    if(mat[pvr+pv] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    mat[cr+j] = (pv==j) ? 0 : mat[cr+j] - mat[pvr+j] * fact;
                }
            }
        }
    }
}
template <> inline void _triu<cmplx>( Mat& utri_mat){
    uint32 rows = utri_mat.getRow();
    uint32 cols = utri_mat.getCol();
    uint32 ch   = utri_mat.getChannel();

    uint32 pivtmax = (rows < cols ) ? rows : cols;
    uint32 pv, i,j, cr, pvr;
    uint32 rcstep = rows*cols;
    uint32 tlen   = rcstep*ch;
    cmplx fact;
    cmplx zero(0,0);
    cmplx* mat = utri_mat.getDataPtr<cmplx>();
    // do triu
    for( uint32 cc=0; cc < tlen ; cc+=rcstep){
        for( pv=0; pv < pivtmax-1 ; pv++){
            //std::cout << "pv = " << pv << " ";
            pvr = pv* cols + cc;
            for( i = pv+1 ; i < rows ; i++){
                cr = i* cols + cc;
                fact = mat[cr+pv] / mat[pvr+pv];
                //std::cout << "cr= " << cr <<", i = " << i << ", pv= " << pv << " pvr= " << pvr <<" fact = " << fact << " ";
                //std::cout << "utri_ma[cr+pv] = " << utri_ma[cr+pv] <<" utri_ma[pvr+pv] = " << utri_ma[pvr+pv] << "\n";
                for(j=0 ; j<cols ; j++){
                    if(mat[pvr+pv] == zero){ std::cout << "Singular Matrix!"; break; }
                    mat[cr+j] = (pv==j) ? zero : mat[cr+j] - mat[pvr+j] * fact;
                }
            }
        }
    }
}

template <typename _T> void _tril( Mat& ltri_mat){
    uint32 rows = ltri_mat.getRow();
    uint32 cols = ltri_mat.getCol();
    uint32 ch   = ltri_mat.getChannel();

    uint32 pivtmax = (rows < cols) ? rows : cols;
    uint32 rcstep = rows*cols;
    uint32 tlen = rcstep*ch;
    uint32   pv, j, cr, pvr;
    int32 i;
    double fact;
    _T* mat = ltri_mat.getDataPtr<_T>();
    // do tril
    for(uint32 cc=0; cc < tlen ; cc+= rcstep){
        for(pv=pivtmax-1 ; pv>0 ; pv--){
            pvr = pv*cols + cc;
            for(i=pv-1 ; i >= 0; i--){
                cr = i* cols + cc;
                fact = mat[cr+pv] / mat[pvr+pv];
                //fprintf(stdout,"pv=%d, i=%d, fact=%f",pv,i,fact);
                for(j=0 ; j < cols; j++){
                    if(mat[pvr+pv] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    mat[cr+j] = (pv==j) ? 0 : mat[cr+j] - mat[pvr+j] * fact;
                }
            }
        }
    }
}
template <> inline void _tril<cmplx>( Mat& ltri_mat){
    uint32 rows = ltri_mat.getRow();
    uint32 cols = ltri_mat.getCol();
    uint32 ch   = ltri_mat.getChannel();

    uint32 pivtmax = (rows < cols) ? rows : cols;
    uint32 rcstep = rows*cols;
    uint32 tlen = rcstep*ch;
    uint32   pv, j, cr, pvr;
    int32 i;
    cmplx fact;
    cmplx zero(0,0);
    cmplx* mat = ltri_mat.getDataPtr<cmplx>();
    // do tril
    for(uint32 cc=0; cc < tlen ; cc+= rcstep){
        for(pv=pivtmax-1 ; pv>0 ; pv--){
            pvr = pv*cols + cc;
            for(i=pv-1 ; i >= 0; i--){
                cr = i* cols + cc;
                fact = mat[cr+pv] / mat[pvr+pv];
                //fprintf(stdout,"pv=%d, i=%d, fact=%f",pv,i,fact);
                for(j=0 ; j < cols; j++){
                    if(mat[pvr+pv] == zero){ std::cout << "Singular Matrix!"; break; }
                    mat[cr+j] = (pv==j) ? zero : mat[cr+j] - mat[pvr+j] * fact;
                }
            }
        }
    }
}
template <typename _T> std::shared_ptr<uchar> _augment(const Mat& srcmat, const uint32 augCols){
    uint32 rows = srcmat.getRow();
    uint32 cols = srcmat.getCol();
    uint32 ch   = srcmat.getChannel();

    uint32 pivmax = (rows < cols)? rows : cols;
    uint32 augmentCols = cols + pivmax;

    assert( augmentCols == augCols );

    uint32 len     = rows*augCols*ch;
    uint32 bytelen = len*sizeof(_T);

    std::shared_ptr<uchar> augm = std::shared_ptr<uchar>(new uchar[bytelen], std::default_delete<uchar[]>());
    _T* augm_ma = (_T*)augm.get();
    _T* mA      = srcmat.getDataPtr<_T>();

    uint32 i,j;
    uint32 cr,scr;
    uint32 cis, cio;
    uint32 rcstep_o  = rows*augCols;
    uint32 rcstep_s  = rows*cols;
    //--augmenting
    for( cis=0, cio=0; cio < len ; cis += rcstep_s, cio += rcstep_o){
        for( i=0 ; i < rows ; i++){
            scr = i*cols + cis;
            cr  = i*augmentCols + cio;
            for( j=0; j < augmentCols ; j++){
                if(j < cols)
                    augm_ma[cr+j] = mA[scr+j] ;
                else if( i == j-cols)
                    augm_ma[cr+j] =1.0;
                else
                    augm_ma[cr+j] =0.0;
            }
        }
    }

    return augm;
}

template <typename _T> Mat _inverse(const Mat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    if(!(srcDtype == DTYP::DOUBLE || srcDtype == DTYP::FLOAT)){
        assert(false && "data type of srcmat into inverse is neither DOUBLE nor FLOAT");
        fprintf(stderr,"data type of srcmat into inverse is neither DOUBLE nor FLOAT\n");
        return Mat();
    }
    uint32 rows = srcmat.getRow();
    uint32 cols = srcmat.getCol();
    uint32 chs  = srcmat.getChannel();

    if(rows != cols) {
        std::cout << "The inverse matrix cannot be computed because source matrix is not square!";
        return Mat();
    }

    Mat mataug = augment(srcmat);
    Mat utri   = triu(mataug);
    Mat ltri   = tril(utri);
    uint32 ltri_col = ltri.getCol();
    uint32 pv,j, pvr;

    double pivot;
    Mat invmat;
    uint32 rc = ltri.getRow() * ltri.getCol();
    uint32 rc2 = rows*cols;
    uint32 pvmax= rows;
    uint32 ci, cii;
    uint32 invchr_off , ltrichr_off;

    invmat = Mat(srcDtype, rows,cols,chs);
    _T* invmat_pt = invmat.getDataPtr<_T>();
    _T* ltri_pt   = ltri.getDataPtr<_T>();

    for( ci = 0; ci < ltri.getLength(); ci +=rc){
        for( pv=0 ; pv < pvmax ; pv++){
            pvr   = pv * ltri_col;
            ltrichr_off = ci + pvr;
            pivot = ltri_pt[ltrichr_off + pv];
            if(pivot==0.0) continue;
            for(j=0 ; j < ltri_col ; j++)
                ltri_pt[ltrichr_off + j] /= pivot;
        }
    }
    uint32 cr;
    for( ci = 0, cii=0 ; ci < ltri.getLength(); ci += rc, cii += rc2){
        for(uint32 i=0; i < rows ; i++){
            pvr = i * ltri_col + cols;
            cr  = i * cols;
            invchr_off = cii+cr;
            ltrichr_off = ci+pvr;
            for( j=0 ; j < cols ; j++){
                invmat_pt[invchr_off+j] = ltri_pt[ltrichr_off+j];
            }
        }
    }

    return invmat;
}
template <> inline Mat _inverse<cmplx>(const Mat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    if(!(srcDtype == DTYP::CMPLX || srcDtype == DTYP::CMPLX)){
        assert(false && "data type of srcmat into inverse is not CMPLX");
        fprintf(stderr,"data type of srcmat into inverse is not CMPLX\n");
        return Mat();
    }
    uint32 rows = srcmat.getRow();
    uint32 cols = srcmat.getCol();
    uint32 chs  = srcmat.getChannel();

    if(rows != cols) {
        std::cout << "The inverse matrix cannot be computed because source matrix is not square!";
        return Mat();
    }

    Mat mataug = augment(srcmat);
    Mat utri   = triu(mataug);
    Mat ltri   = tril(utri);
    uint32 ltri_col = ltri.getCol();
    uint32 pv,j, pvr;

    cmplx pivot;
    Mat invmat;
    uint32 rc = ltri.getRow() * ltri.getCol();
    uint32 rc2 = rows*cols;
    uint32 pvmax= rows;
    uint32 ci, cii;
    uint32 invchr_off , ltrichr_off;

    invmat = Mat(srcDtype, rows,cols,chs);
    cmplx* invmat_pt = invmat.getDataPtr<cmplx>();
    cmplx* ltri_pt   = ltri.getDataPtr<cmplx>();

    for( ci = 0; ci < ltri.getLength(); ci +=rc){
        for( pv=0 ; pv < pvmax ; pv++){
            pvr   = pv * ltri_col;
            ltrichr_off = ci + pvr;
            pivot = ltri_pt[ltrichr_off + pv];
            if(pivot.re==0.0 ) continue;
            for(j=0 ; j < ltri_col ; j++)
                ltri_pt[ltrichr_off + j] /= pivot;
        }
    }
    uint32 cr;
    for( ci = 0, cii=0 ; ci < ltri.getLength(); ci += rc, cii += rc2){
        for(uint32 i=0; i < rows ; i++){
            pvr = i * ltri_col + cols;
            cr  = i * cols;
            invchr_off = cii+cr;
            ltrichr_off = ci+pvr;
            for( j=0 ; j < cols ; j++){
                invmat_pt[invchr_off+j] = ltri_pt[ltrichr_off+j];
            }
        }
    }

    return invmat;
}
template <typename _Ta, typename _Tb, typename _To>
void _conv2d(const Mat& mA, const Mat& mB, Mat& mO, const bool fullout, const std::string& opt_conv ){

    uint32 tX , tY;
    uint32 xdummy , ydummy;
    uint32 mbHcol = mB.getCol()/2;
    uint32 mbHrow = mB.getRow()/2;
    uint32 ch     = mA.getChannel();

    if(fullout){
        xdummy = mB.getCol() -1;
        ydummy = mB.getRow() -1;
        tX = mA.getCol() + xdummy*2;
        tY = mA.getRow() + ydummy*2;
    }else{
        tX = mA.getCol() + mB.getCol() -1;
        tY = mA.getRow() + mB.getRow() -1;
    }

    //Mat tmpA(mA.getDatType(),tY, tX, ch);
    Mat tmpA = Mat::zeros(tY, tX, ch, mA.getDatType());
    uint32 ich, y, x;
    if( fullout ){
        for( ich=0; ich < ch; ich++){
            for( y =0; y < tY; y++){
                for( x=0; x < tX; x++){
                    if( (y >= xdummy && y < tY-ydummy) && ( x >= xdummy && x < tX-ydummy) )
                        tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-ydummy,x-xdummy,ich);
//                    else
//                        tmpA.at<_Ta>(y,x,ich) = 0;
                }
            }
        }
    }else{
        for( ich=0; ich < ch; ich++){
            for( y =0; y < tY; y++){
                for( x=0; x < tX; x++){
                    if(opt_conv.substr(0,4).compare("symm")==0){

                        if( y < mbHrow && x < mbHcol )         // left top corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-y-1, mbHcol-x-1, ich);
                        else if( y < mbHrow && x >= tX-mbHcol) // right top corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-y-1, mA.getCol()-x-mbHcol+tX-1, ich);
                        else if( y < mbHrow)                   // top
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-y-1, x-mbHcol, ich);
                        else if( y >= tY-mbHrow && x < mbHcol) // left bottom corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-y-mbHrow+tY-1, mbHcol-x-1,ich);
                        else if( y >= tY-mbHrow && x >= tX-mbHcol) // right bottom corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-y-mbHrow+tY-1, mA.getCol()-x-mbHcol+tX-1, ich);
                        else if( y >= tY-mbHrow)               // bottom
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-y-mbHrow+tY-1, x-mbHcol, ich);
                        else if( x < mbHcol )                  // left
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow, mbHcol-x-1, ich);
                        else if( x >= tX - mbHcol)             // right
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow, mA.getCol()-x-mbHcol+tX-1, ich);
                        else                                   // main
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow,x-mbHcol,ich);
                    }else if(opt_conv.substr(0,4).compare("circ")==0){

                        if( y < mbHrow && x < mbHcol )         // left top corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-mbHrow+y, mA.getCol()-mbHcol+x, ich);
                        else if( y < mbHrow && x >= tX-mbHcol) // right top corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-mbHrow+y, mbHcol-tX+x, ich);
                        else if( y < mbHrow)                   // top
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mA.getRow()-mbHrow+y, x-mbHcol, ich);
                        else if( y >= tY-mbHrow && x < mbHcol) // left bottom corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-tY+y, mA.getCol()-mbHcol+x,ich);
                        else if( y >= tY-mbHrow && x >= tX-mbHcol) // right bottom corner
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-tY+y, mbHcol-tX+x, ich);
                        else if( y >= tY-mbHrow)               // bottom
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(mbHrow-tY+y, x-mbHcol, ich);
                        else if( x < mbHcol )                  // left
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow, mA.getCol()-mbHcol+x, ich);
                        else if( x >= tX - mbHcol)             // right
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow, mbHcol-tX+x, ich);
                        else                                   // main
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow,x-mbHcol,ich);
                    }else{

                        if( (y >= mbHrow && y < tY-mbHrow) && ( x >= mbHcol && x < tX-mbHcol) )
                            tmpA.at<_Ta>(y,x,ich) = mA.at<_Ta>(y-mbHrow,x-mbHcol,ich);
//                        else
//                            tmpA.at<_Ta>(y,x,ich) = 0;
                    }
                }
            }
        }
    }
    //tmpA.printMat("tmpA");
    //-- convolution
    _To sum;
    if(fullout){
        for( ich=0; ich < ch ; ich++){
            for( y=ydummy; y < tY; y++){
                for( x=xdummy; x < tX; x++){
                    sum = 0;
                    for(int32 m= mB.getRow()-1; m >= 0; m--){
                        for(int32 n= mB.getCol()-1; n >= 0; n--){
                            sum += (tmpA.at<_Ta>(y-m,x-n,ich)* mB.at<_Tb>(m,n,ich));
                            //a = tmpA.at<_Ta>(y-mbHrow+m,x-mbHcol+n,ich);
                            //b = mB.at<_Tb>(m,n,ich);
                            //sum += (a+b);
                        }
                    }
                    mO.at<_To>(y-ydummy,x-xdummy,ich) = sum;
                }
            }
        }
    }else{
        for( ich=0; ich < ch ; ich++){
            for( y=mbHrow; y < tY-mbHrow; y++){
                for( x=mbHcol; x < tX-mbHcol; x++){
                    sum = 0.0;
                    for(int32 m= mB.getRow()-1; m >= 0; m--){
                        for(int32 n= mB.getCol()-1; n >= 0; n--){
                            //sum += (tmpA(y-m,x-n,ich)* mB(-m+mbHrow,-n+mbHcol,ich));
                            sum += (tmpA.at<_Ta>(y+mbHrow-m,x+mbHcol-n,ich)* mB.at<_Tb>(m,n,ich));
                            //a = tmpA.at<_Ta>(y-mbHrow+m,x-mbHcol+n,ich);
                            //b = mB.at<_Tb>(m,n,ich);
                            //sum += (a+b);
                        }
                    }
                    mO.at<_To>(y-mbHrow,x-mbHcol,ich) = sum;
                }
            }
        }
    }
}

}
#endif // JBMATH_H

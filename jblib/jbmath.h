#ifndef JBMATH_H
#define JBMATH_H
#include <stddef.h>
#include "jbMat.h"
#include <string>
//#include <memory>
#include <iostream>
#include <assert.h>

namespace  jbmath {

    jbMat dot(const jbMat& mA,const jbMat& mB);
    jbMat triu(const jbMat& mA);
    jbMat tril(const jbMat& mA);
    jbMat augment(const jbMat& mA);
    jbMat inverse(const jbMat& mA);
    jbMat tranpose(const jbMat& mA);

    jbMat conv2d(const jbMat& mA, const jbMat& mB, std::string opt_conv="" , std::string opt_out="");


    template <typename _Ta, typename _Tb, typename _To > bool _dot_prod(const rawMat rMA,const rawMat rMB, rawMat rMO);
    template <typename _T > void _triu( _T* utri_ma, const uint rows, const uint cols, const uint ch);
    template <typename _T > void _tril( _T* ltri_ma, const uint rows, const uint cols, const uint ch);
    template <typename _T > std::shared_ptr<uchar> _augment(const _T* mA, const uint rows, const uint cols, const uint ch, const uint augCols);
    template <typename _T> jbMat _inverse(const jbMat& srcmat);
}

template <typename _Ta, typename _Tb, typename _To>
bool jbmath::_dot_prod(const rawMat rMA,const rawMat rMB, rawMat rMO){
    uint Ar  = rMA.rows;
    uint Ac  = rMA.cols;
    uint Ach = rMA.chennels;
    uint Br  = rMB.rows;
    uint Bc  = rMB.cols;
    uint Bch = rMB.chennels;
    _Ta* MA  = reinterpret_cast<_Ta*>(rMA.dat_ptr); //(_Ta*) rMA.dat_ptr;
    _Tb* MB  = reinterpret_cast<_Tb*>(rMB.dat_ptr); //(_Tb*) rMB.dat_ptr;
    _To* MO  = reinterpret_cast<_To*>(rMO.dat_ptr); //(_To*) rMO.dat_ptr;

    if( Ac != Br || Ach != Bch ){
        fprintf(stderr, "sizes of ma and mb into _dot_prod_ are not match!\n");
        return false;
    }else if( rMO.rows != Ar || rMO.cols != Bc || rMO.chennels != Ach){
        fprintf(stderr, "sizes of mo is not enough!\n");
        return false;
    }else if( MA == nullptr || MB == nullptr || MO == nullptr){
        fprintf(stderr, "one or more of ma, mb or mo into _dot_prod_ are NULL!\n");
        return false;
    }

    uint i,j;
    uint k, lr, lc, rc, rra, rrb;
    _Ta av;
    _Tb bv;
    _To cv;

    uint mo_chidx, ma_chidx, mb_chidx;
    uint Arc = Ar * Ac;
    uint Brc = Br * Bc;
    uint Orc = Ar * Bc;
    uint m;
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
template <typename _T> void jbmath::_triu( _T* utri_ma, const uint rows, const uint cols, const uint ch){
    uint pivtmax = (rows < cols ) ? rows : cols;

    uint pv, i,j, cr, pvr;
    double fact;
    uint rcstep = rows*cols;
    uint tlen = rcstep*ch;
    // do triu
    for( uint cc=0; cc < tlen ; cc+=rcstep){
        for( pv=0; pv < pivtmax-1 ; pv++){
            //std::cout << "pv = " << pv << " ";
            pvr = pv* cols + cc;
            for( i = pv+1 ; i < rows ; i++){
                cr = i* cols + cc;
                fact = utri_ma[cr+pv] / utri_ma[pvr+pv];
                //std::cout << "cr= " << cr <<", i = " << i << ", pv= " << pv << " pvr= " << pvr <<" fact = " << fact << " ";
                //std::cout << "utri_ma[cr+pv] = " << utri_ma[cr+pv] <<" utri_ma[pvr+pv] = " << utri_ma[pvr+pv] << "\n";
                for(j=0 ; j<cols ; j++){
                    if(utri_ma[pvr+pv] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    utri_ma[cr+j] = (pv==j) ? 0 : utri_ma[cr+j] - utri_ma[pvr+j] * fact;
                }
            }
        }
    }
}

template <typename _T> void jbmath::_tril( _T* ltri_ma, const uint rows, const uint cols, const uint ch){

    uint pivtmax = (rows < cols) ? rows : cols;
    uint rcstep = rows*cols;
    uint tlen = rcstep*ch;
    uint   pv, j, cr, pvr;
    int i;
    double fact;

    // do tril
    for(uint cc=0; cc < tlen ; cc+= rcstep){
        for(pv=pivtmax-1 ; pv>0 ; pv--){
            pvr = pv*cols + cc;
            for(i=pv-1 ; i >= 0; i--){
                cr = i* cols + cc;
                fact = ltri_ma[cr+pv] / ltri_ma[pvr+pv];
                //fprintf(stdout,"pv=%d, i=%d, fact=%f",pv,i,fact);
                for(j=0 ; j < cols; j++){
                    if(ltri_ma[pvr+pv] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    ltri_ma[cr+j] = (pv==j) ? 0 : ltri_ma[cr+j] - ltri_ma[pvr+j] * fact;
                }
            }
        }
    }
}

template <typename _T> std::shared_ptr<uchar> jbmath::_augment(const _T* mA, const uint rows, const uint cols, const uint ch, const uint augCols){

    uint pivmax = (rows < cols)? rows : cols;
    uint augmentCols = cols + pivmax;

    assert( augmentCols == augCols );

    uint len     = rows*augCols*ch;
    uint bytelen = len*sizeof(_T);

    std::shared_ptr<uchar> augm = std::shared_ptr<uchar>(new uchar[bytelen], std::default_delete<uchar[]>());
    _T* augm_ma = (_T*)augm.get();


    uint i,j;
    uint cr,scr;
    uint cis, cio;
    uint rcstep_o  = rows*augCols;
    uint rcstep_s  = rows*cols;
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

template <typename _T>
jbMat jbmath::_inverse(const jbMat& srcmat){
    DTYP srcDtype = srcmat.getDatType();

    if(!(srcDtype == DTYP::DOUBLE || srcDtype == DTYP::FLOAT)){
        assert(false && "data type of srcmat into inverse is neither DOUBLE nor FLOAT");
        fprintf(stderr,"data type of srcmat into inverse is neither DOUBLE nor FLOAT\n");
        return jbMat();
    }
    uint rows = srcmat.getRow();
    uint cols = srcmat.getCol();
    uint chs  = srcmat.getChannel();

    if(rows != cols) {
        std::cout << "The inverse matrix cannot be computed because source matrix is not square!";
        return jbMat();
    }


    jbMat mataug = augment(srcmat);
    jbMat utri   = triu(mataug);
    jbMat ltri   = tril(utri);        
    uint ltri_col = ltri.getCol();
    uint pv,j, pvr;

    double pivot;
    jbMat invmat;
    uint rc = ltri.getRow() * ltri.getCol();
    uint rc2 = rows*cols;
    uint pvmax= rows;
    uint ci, cii;
    uint invchr_off , ltrichr_off;

    invmat = jbMat(srcDtype, rows,cols,chs);
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
    uint cr;
    for( ci = 0, cii=0 ; ci < ltri.getLength(); ci += rc, cii += rc2){
        for(uint i=0; i < rows ; i++){
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
#endif // JBMATH_H

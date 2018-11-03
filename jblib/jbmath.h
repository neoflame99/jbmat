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
    void printMat(const jbMat& Mat);
    jbMat conv2d(const jbMat& mA, const jbMat& mB, std::string opt_conv="" , std::string opt_out="");


    template <typename _Ta, typename _Tb, typename _To > bool _dot_prod(const _Ta* MA,const _Tb* MB, _To* MO,const uint Ar,const uint Ac,const uint Ach,const uint Br,const uint Bc,const uint Bch);
    template <typename _T > void _triu( _T* utri_ma, const uint rows, const uint cols, const uint ch);
    template <typename _T > void _tril( _T* ltri_ma, const uint rows, const uint cols, const uint ch);
    template <typename _T > std::shared_ptr<uchar> _augment(const _T* mA, const uint rows, const uint cols, const uint ch, const uint augCols);
    template <typename _T> jbMat _inverse(const jbMat& srcmat);
}

template <typename _Ta, typename _Tb, typename _To>
bool jbmath::_dot_prod(const _Ta* MA,const _Tb* MB, _To* MO, const uint Ar,const uint Ac,const uint Ach,const uint Br,const uint Bc,const uint Bch){
    if( Ac != Br || Ach != Bch ){
        fprintf(stderr, "sizes of ma and mb into _dot_prod_ are not match!\n");
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
    for(uint m=0; m < Ach; m++ ){
        for( i = 0 , lr =0, rra=0 ; i < Ar ; i++, lr += Bc*Ach, rra+= Ac*Ach){
            for( j = 0, lc=0; j< Bc ; j++, lc+= Ach ){
                MO[lr+lc+m] = 0;
                cv=0;
                for( k=0, rrb=0, rc =0; k < Ac; k++, rrb += Ach*Bc, rc += Ach){
                    av = MA[rra+rc+m];
                    bv = MB[rrb+lc+m];
                    cv = av*bv;
                    MO[lr+lc+m] += cv;
                    //MO[lr+lc+m] += MA[rra+rc+m] * MB[rrb+lc+m];
                }
            }
        }
    }
    return true;
}

template <typename _T> void jbmath::_triu( _T* utri_ma, const uint rows, const uint cols, const uint ch){
    uint pivtmax = (rows < cols ) ? rows : cols;

    uint   pv, i,j, cr, pvr;
    double fact;

    // do triu
    for( uint cc=0; cc < ch; cc++){
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

    uint   pv, i,j, cr, pvr;
    double fact;

    // do tril
    for(uint cc=0; cc < ch; cc++){
        for(pv=pivtmax-1 ; pv>0 ; pv--){
            pvr = pv*cols*ch;
            for(i=pv-1 ; i >= 0; i--){
                cr = i* cols*ch;
                fact = ltri_ma[cr+pv*ch] / ltri_ma[pvr+pv*ch];
                //fprintf(stdout,"pv=%d, i=%d, fact=%f",pv,i,fact);
                for(j=0 ; j < cols; j++){
                    if(ltri_ma[pvr+pv*ch] == 0.0){ std::cout << "Singular Matrix!"; break; }
                    ltri_ma[cr+j*ch] = (pv==j) ? 0 : ltri_ma[cr+j*ch] - ltri_ma[pvr+j*ch] * fact;
                }
            }
        }
    }
}

template <typename _T> std::shared_ptr<uchar> jbmath::_augment(const _T* mA, const uint rows, const uint cols, const uint ch, const uint augCols){

    uint pivmax = (rows < cols)? rows : cols;
    uint augmentCols = cols + pivmax;

    assert( augmentCols == augCols );

    uint len     = rows*cols*ch;
    uint bytelen = len*sizeof(_T);

    std::shared_ptr<uchar> augm = std::shared_ptr<uchar>(new uchar[bytelen], std::default_delete<uchar[]>());
    _T* augm_ma = (_T*)augm.get();

    uint i,j;
    uint cr,scr;
    // data copy
    for( i = 0 ; i < len; i++ )
        augm_ma[i] = mA[i];

    //--augmenting
    for( i=0 ; i < rows ; i++){
        scr = i*cols;
        cr  = i*augmentCols;
        for( j=0; j < augmentCols ; j++){
            if(j < cols)
                augm_ma[cr+j] = mA[scr+j] ;
            else if( i == j-cols)
                augm_ma[cr+j] =1.0;
            else
                augm_ma[cr+j] =0.0;
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

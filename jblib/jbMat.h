/*
 * Copyright (C) 2020. Jong B. Choi
 * License : MIT License
 * contact : neoflame99@naver.com
 *
 */

#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>

#include <initializer_list>
#include <vector>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <string>
#include <stdexcept>
#include "types.h"
#include <iostream>

namespace jmat {

typedef union _elemptr{
    uchar*  uch_ptr = nullptr;
    int*    int_ptr ;
    float*  f32_ptr ;
    double* f64_ptr ;
    cmplx*  cmx_ptr ;
}elemptr;

typedef struct _matRect{
    int32 sR, sC;
    int32 eR, eC;
    _matRect() = default;
    _matRect(int32 sr, int32 sc, int32 er, int32 ec ):sR(sr),sC(sc),eR(er),eC(ec){}
    //_matRect(std::initializer_list<int32> list){
    //    std::vector<int32> v(4);
    //    v.insert(v.begin(),list.begin(),list.end());
    //    sR = v.at(0);
    //    sC = v.at(1);
    //    eR = v.at(2);
    //    eC = v.at(3);
    // };
    inline void set(int32 sr=0, int32 sc=0, int32 er=0, int32 ec=0){
       sR= sr; sC=sc; eR= er; eC=ec;
    }
    bool inline operator== (const _matRect& other){
        return (this->sR == other.sR && this->sC == other.sC && this->eR == other.eR && this->eC == other.eC );
    }
    _matRect& operator<< (const int32 v){
        sR = sC;
        sC = eR;
        eR = eC;
        eC = v;
        return *this;
    }
}matRect;

class Mat;

class Mat{
private: // member fields
    shr_ptr mA;
    uchar   *dat_ptr;
    elemptr elptr;
    DTYP    datT;

    uint32 length;
    uint32 lenRowCol;
    uint32 row, col, Nch;
    uint32 byteStep;
    uint32 byteLen;
    uint32 stepCol;
    uint32 stepRow;

    std::string obj_name;

private: // initializing methods
    inline void sync_data_ptr(){ elptr.uch_ptr = dat_ptr = mA.get(); }
    inline void init(const uint32 r, const uint32 c, const uint32 ch,const DTYP dt,const bool do_alloc=true);
    inline void alloc(const uint32 len);

public: // constructors and destructor
    Mat();
    Mat(const DTYP dt, const uint32 r, const uint32 c, const uint32 ch );
    Mat(const DTYP dt, const uint32 r, const uint32 c, const uint32 ch, const std::string name);
    Mat(const DTYP dt, const uint32 rc);
    Mat(const shr_ptr ma, const DTYP dt, const uint32 r, const uint32 c, const uint32 ch );
    Mat(const Mat& mat); // copy constructor
    Mat(Mat&& mat);      // move constructor
    Mat(const std::initializer_list<double> list );
    Mat(const std::initializer_list<int32 > list );
    Mat(const std::initializer_list<float > list );
    Mat(const std::initializer_list<uchar > list );
    Mat(const std::initializer_list<cmplx > list );
    ~Mat();

public:
    //-- overloading operators
    Mat  operator-() const;  // unary minus

    Mat  operator+(const Mat& other ) const;
    Mat  operator-(const Mat& other ) const;
    Mat  operator*(const Mat& other ) const;
    Mat  operator/(const Mat& other ) const;

    friend Mat operator+(const Mat& A, const double scalar);
    friend Mat operator-(const Mat& A, const double scalar);
    friend Mat operator*(const Mat& A, const double scalar);
    friend Mat operator/(const Mat& A, const double scalar);

    friend Mat operator+(const double scalar, const Mat& A);
    friend Mat operator-(const double scalar, const Mat& A);
    friend Mat operator*(const double scalar, const Mat& A);
    friend Mat operator/(const double scalar, const Mat& A);

    Mat& operator= (const Mat& other ); // copy assignment
    Mat& operator= (Mat&& other );      // move assignment
    Mat& operator+=(const Mat& other );
    Mat& operator+=(const double scalar);
    Mat& operator-=(const Mat& other );
    Mat& operator-=(const double scalar);
    Mat& operator*=(const Mat& other );
    Mat& operator*=(const double scalar);
    Mat& operator/=(const Mat& other );
    Mat& operator/=(const double scalar);

    template <typename _T> inline _T& operator[](const uint32 i);
    template <typename _T> inline _T  operator[](const uint32 i) const;

    Mat& plusMat (const Mat& other);
    Mat& minusMat(const Mat& other);
    Mat& mulMat  (const Mat& other);
    Mat& divMat  (const Mat& other);

    Mat& plusScalar(const double scalar);
    Mat& plusScalar(const  float scalar);
    Mat& plusScalar(const  int32 scalar);
    Mat& plusScalar(const  uchar scalar);
    Mat& minusScalar(const double scalar);
    Mat& minusScalar(const  float scalar);
    Mat& minusScalar(const  int32 scalar);
    Mat& minusScalar(const  uchar scalar);
    Mat& mulScalar(const double scalar);
    Mat& mulScalar(const  float scalar);
    Mat& mulScalar(const  int32 scalar);
    Mat& mulScalar(const  uchar scalar);
    Mat& divByScalar(const double scalar);
    Mat& divByScalar(const  float scalar);
    Mat& divByScalar(const  int32 scalar);
    Mat& divByScalar(const  uchar scalar);
    Mat& divScalar(const double scalar);
    Mat& divScalar(const  float scalar);
    Mat& divScalar(const  int32 scalar);
    Mat& divScalar(const  uchar scalar);

    void setRowCol(const uint32 r,const uint32 c,const uint32 ch=1);
    void setChannel(const Mat& src, const uint32 srcCh=0, const uint32 tarCh=0,const uint32 Channels=1);
    void setName(const std::string name);

    Mat   copy() const;
    Mat   copySubMat(const uint32 startRow, const uint32 endRow, const uint32 startCol, const uint32 endCol) const;

    inline bool    isEmpty() const { return ((length <= 0) ? true : false); }
    inline shr_ptr getMat () const { return mA; }
    inline uint32  getLength() const{ return length; }
    inline uint32  getRow() const { return row; }
    inline uint32  getCol() const { return col; }
    inline uint32  getRowColSize() const {return lenRowCol; }
    inline uint32  getChannel() const { return Nch; }
    inline DTYP    getDatType() const { return datT; }
    inline uint32  getByteStep() const { return byteStep; }
    inline uint32  getByteLen() const { return byteLen; }

    int32 reshape(uint32 r, uint32 c, uint32 ch=1);
    void  transpose();
    void  changeDType(const DTYP dt);
    void  printMat() ;
    void  printMat(const std::string objname);
    Mat   max() const;
    Mat   min() const;
    Mat   sum() const;
    Mat   mean() const;
    Mat   std() const;
    Mat   var() const;
    Mat   sqrtm() const;

    inline elemptr getRowElptr(uint32 r=0) const;
    inline elemptr getElptr() const { return elptr; }
private: // other private methods

public : // static methods
    static Mat  ones (uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static Mat  zeros(uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static bool sliceCopyMat(const Mat& src, const matRect& srcSlice,const Mat& des, const matRect& desSlice );
    static Mat  repeat(const Mat& src, const uint32 rp_r, const uint32 rp_c, const uint32 rp_ch);
    template <typename _T> inline static Mat _repeat(const Mat& src, const uint32 r, const uint32 c, const uint32 ch);
    static Mat extractChannel(const Mat& src, const uint32 ch=0);

public : // public template methods
    template <typename _T> _T& at(uint32 i=0) const;
    template <typename _T> _T& at(uint32 r, uint32 c, uint32 nch=0) const;
    template <typename _T=uchar> _T* getDataPtr() const;
    template <typename _T=uchar> _T* getRowPtr(const uint32 r=0) const;
    template <typename _T> Mat _max() ;
    template <typename _T> Mat _min() ;
    template <typename _T> Mat _mean();
    template <typename _T> Mat _std() ;
    template <typename _T> Mat _sum() ;


private: // private template methods
    template <typename _T> void _print(_T* mdat);
    template <typename _Tslf, typename _Totr> void _plus_mat       (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _minus_mat      (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _multiply_mat   (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _divide_mat     (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _plus_scalar    (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _minus_scalar   (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _multiply_scalar(_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _dividing_scalar  (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _divided_by_scalar(_Tslf* self, _Totr scalar, uint32 len );
    template <typename _O, typename _T, typename _R> inline _O _dividing (_T a, _R b);

};

// initializing
inline void Mat::init(const uint32 r,const uint32 c,const uint32 ch,const DTYP dt,const bool do_alloc){
    row = r;
    col = c;
    Nch = ch;
    stepCol   = Nch;
    stepRow   = col*stepCol;
    length    = row*stepRow;
    lenRowCol = row*col;

    datT = dt;
    switch(datT){
    case DTYP::UCHAR  : byteStep = 1; break;
    case DTYP::INT    : byteStep = 4; break;
    case DTYP::FLOAT  : byteStep = 4; break;
    case DTYP::DOUBLE : byteStep = 8; break;
    case DTYP::CMPLX  : byteStep = sizeof(cmplx); break;
    }
    byteLen = length*byteStep;
    if(do_alloc)
        alloc(byteLen);
    sync_data_ptr();
}

//-- shallow copy version & using shared_ptr
inline void Mat::alloc(const uint32 len){
    mA = len==0 ? nullptr : shr_ptr (new uchar[len], std::default_delete<uchar[]>());
}

inline elemptr Mat::getRowElptr(uint32 r) const{
    assert( r < row);
    U64 offset = r*stepRow;
    elemptr ptrs;
    switch(datT){
    case DTYP::DOUBLE: ptrs.f64_ptr = elptr.f64_ptr+offset; break;
    case DTYP::FLOAT : ptrs.f32_ptr = elptr.f32_ptr+offset; break;
    case DTYP::INT   : ptrs.int_ptr = elptr.int_ptr+offset; break;
    case DTYP::UCHAR : ptrs.uch_ptr = elptr.uch_ptr+offset; break;
    case DTYP::CMPLX : ptrs.cmx_ptr = elptr.cmx_ptr+offset; break;
    default          : ptrs.uch_ptr = nullptr;
    }
    return ptrs;
}

template <typename _T> inline _T& Mat::at(uint32 i) const {
    assert(!isEmpty() && i < length);
    return ((_T *)dat_ptr)[i];
}
template <typename _T> inline _T& Mat::at(uint32 r, uint32 c, uint32 ch) const {
    uint32 i = r*stepRow + c*stepCol + ch;
    assert(!isEmpty() && i < length);

    return ((_T *)dat_ptr)[i];
}
template <typename _T> inline _T* Mat::getDataPtr() const {
    return (_T*)dat_ptr;
}
template <> inline cmplx* Mat::getDataPtr() const{
    return reinterpret_cast<cmplx*>(dat_ptr);
}
template <typename _T> inline _T* Mat::getRowPtr(const uint32 r) const{
    assert( r < row);
    if( r >= row ) throw std::out_of_range("r is out of range ");
    _T* ptr = ((_T*)dat_ptr)+(r*stepRow);
    return ptr;
}

template <typename _T> inline _T& Mat::operator[](const uint32 idx){
    assert( dat_ptr!=nullptr );
    assert( idx < byteLen );
    if( idx >= byteLen ) throw std::out_of_range("idx is out of range ");
    return ((_T*)dat_ptr)[idx];
}
template <typename _T> inline _T Mat::operator[](const uint32 idx) const {
    assert( dat_ptr!=nullptr );
    assert( idx < byteLen );
    if( idx >= byteLen ) throw std::out_of_range("idx is out of range ");
    return ((_T*)dat_ptr)[idx];
}

// array arithmetic methods
template <typename _Tslf, typename _Totr> void Mat::_plus_mat(_Tslf* self, _Totr* other, uint32 len){
    for(uint32 k=0; k < len; k++ )
        self[k] += other[k];
}
template <typename _Tslf, typename _Totr> void Mat::_minus_mat(_Tslf* self, _Totr* other, uint32 len){
    for(uint32 k=0; k < len; k++ )
        self[k] -= other[k];
}
template <typename _Tslf, typename _Totr> void Mat::_multiply_mat(_Tslf* self, _Totr* other, uint32 len){
    for(uint32 k=0; k < len; k++ )
        self[k] *= other[k];
}
template <typename _Tslf, typename _Totr> void Mat::_divide_mat(_Tslf* self, _Totr* other, uint32 len){
    for(uint32 k=0; k < len; k++ )
        self[k] /= other[k];
}
template <typename _Tslf, typename _Totr> void Mat::_plus_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] += scalar;
}
template <typename _Tslf, typename _Totr> void Mat::_minus_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] -= scalar;
}
template <typename _Tslf, typename _Totr> void Mat::_multiply_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] *= scalar;
}
template <typename _Tslf, typename _Totr> void Mat::_divided_by_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] /= scalar;
}
template <typename _Tslf, typename _Totr> void Mat::_dividing_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] = scalar / self[k];
}
/*
template <typename _T> inline void Mat::_print(_T* mdat){
    const int32 bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    uint32 i,j,k;

    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;

    if(datT==DTYP::DOUBLE || datT==DTYP::FLOAT){
        double val;
        for( i = 0; i < length; i += stepRow){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < stepRow; j+= stepCol ){ // columns
                snprintf(tmp, bufsz," (");
                strncat(buf, tmp, bufsz);
                for( k =0; k < Nch; ++k){
                    val = mdat[i+j+k];
                    if( val >= neg_max_double && val <= pos_min_double)
                        val = 0.0;

                    snprintf(tmp, bufsz,"% 6.3g",val);
                    strncat(buf, tmp, bufsz);
                    if( k < Nch-1)
                        strncat(buf, ",", 2);
                }
                snprintf(tmp, bufsz,") ");
                strncat(buf, tmp, bufsz);
            }
            strncat(buf,"]",2);
            fprintf(stdout,"%s\n",buf);
        }
    }else{
        int32 val;
        for( i = 0; i < length; i += stepRow){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < stepRow; j+= stepCol){ // columns
                snprintf(tmp,bufsz," (");
                strncat(buf,tmp,bufsz);
                for( k = 0 ; k < Nch; k++){
                    val = mdat[i+j+k];
                    if(val > 100000)
                        snprintf(tmp,bufsz,"% 5g",(double)val);
                    else
                        snprintf(tmp,bufsz,"% 5d",val);

                    strncat(buf,tmp,bufsz);
                    if( k < Nch-1)
                        strncat(buf, ",", 2);
                }
                snprintf(tmp,bufsz,") ");
                strncat(buf,tmp,bufsz);
            }
            strncat(buf,"]",2);
            fprintf(stdout,"%s\n",buf);
        }
    }
}
*/
template <typename _T> inline void Mat::_print(_T* mdat){
    uint32 i,j,k;
    std::stringstream ss;


    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;

    if(datT==DTYP::DOUBLE || datT==DTYP::FLOAT){
        ss.precision(3);
        ss.width(6);
        ss.setf(std::ios::scientific);
        double val;
        for( i = 0; i < length; i += stepRow){ // rows
            ss << "[";
            for( j=0; j < stepRow; j+= stepCol ){ // columns
                ss << " (";
                for( k =0; k < Nch; ++k){
                    val = mdat[i+j+k];
                    if( val >= neg_max_double && val <= pos_min_double)
                        val = 0.0;

                    ss << val ;
                    if( k < Nch-1) {
                        ss << ",";
                    }
                }
                ss << ") ";
            }
            ss << "]";
            std::cout << ss.str();
        }
    }else{        
        ss.width(5);
        int32 val;
        for( i = 0; i < length; i += stepRow){ // rows
            ss << "[";
            for( j=0; j < stepRow; j+= stepCol){ // columns
                ss << " (";
                for( k = 0 ; k < Nch; k++){
                    val = mdat[i+j+k];
                    if(val > 100000){
                        ss.setf(std::ios::scientific);
                        ss << (double)val;
                    }else{
                        ss.flags ( std::ios::right | std::ios::dec | std::ios::showbase );
                        ss << val;
                    }
                    if( k < Nch-1){
                        ss << ",";
                    }
                }
                ss << ") ";
            }
            ss << "] ";
            std::cout << ss.str();
        }
    }
}
template <> inline void Mat::_print<cmplx>(cmplx* mdat){
    const int32 bufsz = 32;
    char tmp[bufsz];
    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;
    uint32 i,j, k;
    std::stringstream ss;
    std::string sign;

    ss.precision(2);
    ss.width(5);
    ss.setf(std::ios::scientific);
    if(datT==DTYP::CMPLX){
        cmplx val;
        for( i = 0; i < lenRowCol; i += col){ // rows
            ss << "[";
            for( j=0; j < col; j++){          // columns
                ss << " (";
                for( k = 0 ; k < Nch; k++){
                    val = mdat[i+j+k];
                    if( val.re >= neg_max_double && val.re <= pos_min_double)
                        val.re = 0.0;
                    if( val.im >= neg_max_double && val.im <= pos_min_double)
                        val.im = 0.0;

                    snprintf(tmp,bufsz,"% 5.2g %+5.2gi",val.re, val.im);
                    ss << tmp;
                    if( k < Nch-1){
                        ss << ",";
                    }
                }
                ss << ") ";
            }
            ss << "]";
            std::cout << ss.str();
        }
    }
}

template <typename _T> Mat Mat::_max() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch = getChannel();
    uint32 k, m, n;
    Mat A(getDatType(),1,1,ch);

    for(k=0; k < ch; ++k)
        A.at<_T>(k) = datPtr[k];

    for(m = ch ; m < length; m+=ch ){
        for(k=0, n=m ; k < ch; ++k, ++n){
            if( A.at<_T>(k) < datPtr[n])
                A.at<_T>(k) = datPtr[n];
        }
    }
    return A;
}
template <> inline Mat Mat::_max<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch = getChannel();
    uint32 k, m, n;
    cmplx  tmp;
    double *large_mag_ch;
    double  large_mag, tmp_mag;

    Mat A(getDatType(),1,1,ch);

    large_mag_ch = new double[ch];
    for(k=0; k < ch; ++k){
        A.at<cmplx>(k) = datPtr[k];
        large_mag_ch[k] = datPtr[k].square();
    }

    for(m = ch; m < length; m+=ch){
        for(k=0, n=m; k < ch; ++k, ++n){
            tmp       = datPtr[n];
            tmp_mag   = tmp.square();
            large_mag = large_mag_ch[k];
            if(large_mag < tmp_mag){
                A.at<cmplx>(k) = tmp;
                large_mag_ch[k]= tmp_mag;
            }
        }
    }
    delete [] large_mag_ch;

    return A;
}
template <typename _T> Mat Mat::_min() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch = getChannel();
    uint32 k, m, n;
    Mat A(getDatType(),1,1,ch);

    for(k=0; k < ch; ++k)
        A.at<_T>(k) = datPtr[k];

    for(m = ch ; m < length; m+=ch ){
        for(k=0, n=m ; k < ch; ++k, ++n){
            if( A.at<_T>(k) > datPtr[n])
                A.at<_T>(k) = datPtr[n];
        }
    }
    return A;
}
template <> inline Mat Mat::_min<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch = getChannel();
    uint32 k, m, n;
    cmplx  tmp;
    double *large_mag_ch;
    double  large_mag, tmp_mag;

    Mat A(getDatType(),1,1,ch);

    large_mag_ch = new double[ch];
    for(k=0; k < ch; ++k){
        A.at<cmplx>(k) = datPtr[k];
        large_mag_ch[k] = datPtr[k].square();
    }

    for(m = ch; m < length; m+=ch){
        for(k=0, n=m; k < ch; ++k, ++n){
            tmp       = datPtr[n];
            tmp_mag   = tmp.square();
            large_mag = large_mag_ch[k];
            if(large_mag > tmp_mag){
                A.at<cmplx>(k) = tmp;
                large_mag_ch[k]= tmp_mag;
            }
        }
    }
    delete [] large_mag_ch;

    return A;
}
template <typename _T> Mat Mat::_mean() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch = getChannel();
    uint32 k, m, n;
    Mat A = Mat::zeros(1,1,ch,DTYP::DOUBLE);
    n = 0;
    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n)
            A.at<double>(k) += double(datPtr[n]);
    }
    A /= row*col;
    return A;
}
template <> inline Mat Mat::_mean<cmplx>() {
//reference https://en.wikipedia.org/wiki/Complex_random_variable#Expectation
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch    = getChannel();

    Mat A = Mat::zeros(1,1,ch,DTYP::CMPLX);
    uint32 k, m, n;
    n = 0;
    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n)
            A.at<cmplx>(k) += cmplx(datPtr[n]);
    }
    A /= getRowColSize();
    return A;
}
template <typename _T> Mat Mat::_sum() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch  = getChannel();
    Mat A = Mat::zeros(1,1,ch,DTYP::DOUBLE);
    uint32 k, m, n;
    n = 0;
    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n)
            A.at<double>(k) += double(datPtr[n]);
    }
    return A;
}
template <> inline Mat Mat::_sum<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch  = getChannel();
    Mat A = Mat::zeros(1,1,ch,DTYP::CMPLX);
    uint32 k, m, n;
    n = 0;
    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n)
            A.at<cmplx>(k) += cmplx(datPtr[n]);
    }
    return A;
}

template <typename _T> Mat Mat::_std() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch    = getChannel();

    Mat avg = _mean<_T>();
    Mat A   = Mat::zeros(1,1,ch,DTYP::DOUBLE);
    double diff;
    uint32 k, m, n;
    n = 0;
    uint32 Div = getRowColSize() -1;
    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n){
            diff = double(datPtr[n]) - avg.at<double>(k);
            A.at<double>(k) += diff*diff;
        }
    }

    for(k=0; k < ch; ++k){
        double v = A.at<double>(k)/Div;
        A.at<double>(k) = sqrt( v );
    }
    return A;
}
template <> inline Mat Mat::_std<cmplx>() {
//reference https://en.wikipedia.org/wiki/Complex_random_variable#Expectation

    cmplx* datPtr = this->getDataPtr<cmplx>();
    uint32 ch    = getChannel();

    uint32 k, m, n;
    Mat avg = _mean<cmplx>();
    Mat A   = Mat::zeros(1,1,ch,DTYP::DOUBLE);
    cmplx diff;
    n = 0;
    uint32 Div = getRowColSize() -1;
    // STD of complex numbers
    // E[ |Z - E[Z]|^2 ] = E[|Z|^2] - |E[Z]|^2 --> real valued Mat

    for(m = 0 ; m < length; m+=ch){
        for(k=0; k < ch; ++k, ++n){
            diff = datPtr[n] - avg.elptr.cmx_ptr[k];
            A.at<double>(k) += diff.square();
        }
    }
    A /= Div; // Varaince
    for(k=0; k < ch; ++k)  // standard deviation
        A.at<double>(k) = sqrt(A.at<double>(k));

    return A;
}

template <typename _T> inline Mat Mat::_repeat(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    uint32 sr = src.getRow();
    uint32 sc = src.getCol();
    uint32 sch= src.getChannel();
    uint32 nr = sr * rr;
    uint32 nc = sc * rc;
    uint32 nch= sch*rch;

    uint32 x, y, z, i, k, n;
    uint32 ch_xpnd_sz = nch*sc;
    uint32 cl_xpnd_sz = nch*nc;
    Mat des(src.getDatType(), nr, nc, nch);
    _T* sRow_ptr;
    _T* tRow_ptr;
    _T* ch_xpnd;

    for( y=0; y < sr ; ++y){
        sRow_ptr = src.getRowPtr<_T>(y);
        tRow_ptr = des.getRowPtr<_T>(y);
        ch_xpnd  = tRow_ptr;
        // make channel expanding array
        for(i=0; i < rch; ++i){
            for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                for(k=0 ; k < sch; ++k, ++n)
                    ch_xpnd[z+k] = sRow_ptr[n];
            }
        }
        // replicating the expanded channel data into new column of des Mat.
        for(i=1, tRow_ptr += ch_xpnd_sz; i < rc; ++i, tRow_ptr += ch_xpnd_sz){
            memcpy(tRow_ptr, ch_xpnd, src.byteStep*ch_xpnd_sz);
        }
    }
    // replicating the expanded column Mat into remaining rows of des Mat.
    for(y=0; y < sr ; ++y){
        sRow_ptr = des.getRowPtr<_T>(y);
        for( i=1, n=sr+y; i < rr; ++i, n+=sr){
            tRow_ptr = des.getRowPtr<_T>(n);
            memcpy(tRow_ptr, sRow_ptr, src.byteStep*cl_xpnd_sz);
        }
    }
    return des;
}
/*
template <> inline Mat Mat::_repeat<cmplx>(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    uint32 sr = src.getRow();
    uint32 sc = src.getCol();
    uint32 sch= src.getChannel();
    uint32 nr = sr * rr;
    uint32 nc = sc * rc;
    uint32 nch= sch*rch;

    uint32 x, y, z, i, k, n;
    uint32 ch_xpnd_sz = nch*sc;
    uint32 cl_xpnd_sz = nch*nc;
    Mat des(src.getDatType(), nr, nc, nch);
    cmplx* sRow_ptr;
    cmplx* tRow_ptr;
    cmplx* ch_xpnd;

    for( y=0; y < sr ; ++y){
        sRow_ptr = src.getRowPtr<cmplx>(y);
        tRow_ptr = des.getRowPtr<cmplx>(y);
        ch_xpnd  = tRow_ptr;
        // make channel expanding array
        for(i=0; i < rch; ++i){
            for( x=0, z=i*sch, n=0; x < sc ; x+=sch, z+=nch){
                for(k=0 ; k < sch; ++k, ++n)
                    ch_xpnd[z+k] = sRow_ptr[n];
            }
        }
        // replicating the expanded channel data into new column of des Mat.
        for(i=1, x = ch_xpnd_sz; i < rc; ++i){
            for(k=0 ; k < ch_xpnd_sz; ++k, ++x)
                tRow_ptr[x] = ch_xpnd[k];
        }
        //for(i=1, tRow_ptr += ch_xpnd_sz; i < rc; ++i, tRow_ptr += ch_xpnd_sz){
        //    memcpy(tRow_ptr, ch_xpnd, sizeof(cmplx)*ch_xpnd_sz);
        //    tRow_ptr += ch_xpnd_sz;   // move the pointer in a new column
        //}
    }
    // replicating the expanded column Mat into remaining rows of des Mat.
    for(y=0; y < sr ; ++y){
        sRow_ptr = des.getRowPtr<cmplx>(y);
        for( i=1, n=sr+y; i < rr; ++i, n+=sr){
            tRow_ptr = des.getRowPtr<cmplx>(n);
            for(k=0; k < cl_xpnd_sz; ++k)
                tRow_ptr[k] = sRow_ptr[k];
            //memcpy(tRow_ptr, sRow_ptr, sizeof(cmplx)*cl_xpnd_sz);
        }
    }
    return des;
}*/

} // namespace jmat

#endif // JBMAT_H

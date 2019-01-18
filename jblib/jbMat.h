#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>

#include <initializer_list>
#include <vector>
#include <float.h>
#include <assert.h>
#include "types.h"
#include <math.h>

#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

namespace jmat {

class Mat;

class Mat{
private: // member fields
    shr_ptr mA;
    uchar *dat_ptr;
    DTYP datT;

    uint32 length;
    uint32 lenRowCol;
    void alloc(uint32 len);
    uint32 row, col, Nch;
    uint32 byteStep;
    uint32 byteLen;

    std::string obj_name;

private:
    void sync_data_ptr(){
        if( mA.get() != dat_ptr)
            dat_ptr = mA.get();
    }
    void init(uint32 r, uint32 c, uint32 ch, DTYP dt);
    void initName();


public: // constructors and destructor

    Mat(DTYP dt = DTYP::DOUBLE);
    Mat(DTYP dt, uint32 r, uint32 c, uint32 ch );
    Mat(DTYP dt, uint32 r, uint32 c, uint32 ch, std::string name);
    Mat(DTYP dt, uint32 rc);
    Mat(shr_ptr ma, DTYP dt, uint32 r, uint32 c, uint32 ch );
    Mat(const Mat& mat);
    Mat( std::initializer_list<double> list );
    Mat( std::initializer_list<int32> list );
    Mat( std::initializer_list<float> list );
    ~Mat();

public:
    //-- overloading operators
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

    Mat& operator= (const Mat& other );
    Mat& operator+=(const Mat& other );
    Mat& operator+=(const double scalar);
    Mat& operator-=(const Mat& other );
    Mat& operator-=(const double scalar);
    Mat& operator*=(const Mat& other );
    Mat& operator*=(const double scalar);
    Mat& operator/=(const Mat& other );
    Mat& operator/=(const double scalar);

    Mat& plusMat    (const Mat& other);
    Mat& minusMat   (const Mat& other);
    Mat& multiplyMat(const Mat& other);
    Mat& divideMat  (const Mat& other);

    Mat& plusScalar(const double scalar);
    Mat& plusScalar(const  float scalar);
    Mat& plusScalar(const  int32 scalar);
    Mat& plusScalar(const  uchar scalar);
    Mat& minusScalar(const double scalar);
    Mat& minusScalar(const  float scalar);
    Mat& minusScalar(const  int32 scalar);
    Mat& minusScalar(const  uchar scalar);
    Mat& multiplyScalar(const double scalar);
    Mat& multiplyScalar(const  float scalar);
    Mat& multiplyScalar(const  int32 scalar);
    Mat& multiplyScalar(const  uchar scalar);
    Mat& dividedByScalar(const double scalar);
    Mat& dividedByScalar(const  float scalar);
    Mat& dividedByScalar(const  int32 scalar);
    Mat& dividedByScalar(const  uchar scalar);
    Mat& dividingScalar(const double scalar);
    Mat& dividingScalar(const  float scalar);
    Mat& dividingScalar(const  int32 scalar);
    Mat& dividingScalar(const  uchar scalar);

    void setRowCol(uint32 r, uint32 c, uint32 ch=1);
    void setChannelN(const Mat& src, const uint32 srcfromCh=0,const uint32 Channels=1, const uint32 tarToCh=0);
    void setName(std::string name);

    Mat   copy() const;
    Mat   copyChannelN(const uint32 NoCh=0) const;
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

    int32   reshape(uint32 r, uint32 c, uint32 ch=1);
    void    transpose();
    void    changeDType(const DTYP dt);
    void    printMat() ;
    void    printMat(const std::string objname);

public : // static methods
    static Mat ones (uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static Mat zeros(uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static int32 instant_count;

public : // public template methods
    template <typename _T> _T& at(uint32 i) const;
    template <typename _T> _T& at(uint32 r, uint32 c, uint32 nch=0) const;
    template <typename _T> _T* getDataPtr() const;
    template <typename _T> Mat max() ;
    template <typename _T> Mat min() ;
    template <typename _T> Mat mean();
    template <typename _T> Mat std() ;

private: // private template methods
    template <typename _T> void _print(_T* mdat);
    template <typename _Tsrc, typename _Ttar> void _type_change();
    template <typename _Tslf, typename _Totr> void _plus_mat       (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _minus_mat      (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _multiply_mat   (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _divide_mat     (_Tslf* self, _Totr* other, uint32 len );
    template <typename _Tslf, typename _Totr> void _plus_scalar    (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _minus_scalar   (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _multiply_scalar(_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _dividing_scalar  (_Tslf* self, _Totr scalar, uint32 len );
    template <typename _Tslf, typename _Totr> void _divided_by_scalar(_Tslf* self, _Totr scalar, uint32 len );
};


template <typename _T> inline _T& Mat::at(uint32 i) const {

    assert(!isEmpty() && i < length);

    return ((_T *)dat_ptr)[i];
}
template <typename _T> inline _T& Mat::at(uint32 r, uint32 c, uint32 ch) const {

    uint32 i = ch*lenRowCol + r*col + c;
    assert(!isEmpty() && i < length);

    return ((_T *)dat_ptr)[i];
}

template <typename _T> inline _T* Mat::getDataPtr() const {
    return (_T *)dat_ptr;
}

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
        self[k] = self[k] / scalar;
}
template <typename _Tslf, typename _Totr> void Mat::_dividing_scalar(_Tslf* self, _Totr scalar, uint32 len){
    for(uint32 k=0; k < len ; k++)
        self[k] = scalar / self[k];
}

template <typename _Tsrc, typename _Ttar> void Mat::_type_change(){
    _Tsrc* src_pt = (_Tsrc *)mA.get();
    _Ttar* tar_pt;
    try{
        tar_pt = new _Ttar[static_cast<unsigned long>(length)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Transpose Error: %s\n",ex.what());
        return;
    }

    uint32 k;
    for( k=0; k < length; k++){
        tar_pt[k] = (_Ttar) (src_pt[k]);
    }
    mA.reset((uchar *)tar_pt,std::default_delete<uchar[]>());
    sync_data_ptr();

    byteStep = sizeof(_Ttar);
    byteLen  = length * byteStep;
}

template <typename _T> void Mat::_print(_T* mdat){
    const int32 bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    uint32 i,j;

    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;

    uint32 k;
    uint32 ch_offset;

    if(datT==DTYP::DOUBLE || datT==DTYP::FLOAT){
        double val;
        for( k = 0 ; k < Nch; k++){
            fprintf(stdout,"channel: %d \n",k);
            ch_offset = k*lenRowCol;
            for( i = 0; i < lenRowCol; i += col){ // rows
                snprintf(buf,bufsz,"[");
                for( j=0; j < col; j++){          // columns
                    val = mdat[i+j+ch_offset];
                    if( val >= neg_max_double && val <= pos_min_double)
                        val = 0.0;
                    snprintf(tmp,bufsz," %10.4f ",val);
                    strncat(buf,tmp,bufsz);
                }
                strncat(buf,"]",1);
                fprintf(stdout,"%s\n",buf);
            }
        }
    }else{        
        int32 val;
        for( k = 0 ; k < Nch; k++){
            fprintf(stdout,"channel: %d \n",k);
            ch_offset = k*lenRowCol;
            for( i = 0; i < lenRowCol; i += col){ // rows
                snprintf(buf,bufsz,"[");
                for( j=0; j < col; j++){          // columns
                    val = mdat[i+j+ch_offset];
                    snprintf(tmp,bufsz," %7d ",val);
                    strncat(buf,tmp,bufsz);
                }
                strncat(buf,"]",1);
                fprintf(stdout,"%s\n",buf);
            }
        }
    }
}

template <typename _T> Mat Mat::max() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch  = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(getDatType(),1,ch,1);
    _T large;
    uint32 k, m, n;
    n=0;
    for(k=0 ; k < ch; k++){
        large = datPtr[n++];
        for(m = 1 ; m < rclen; m++, n++ ){
            if( large < datPtr[n])
                large = datPtr[n];
        }
        A.at<_T>(k) = large;
    }
    return A;
}
template <> inline Mat Mat::max<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(getDatType(),1,ch,1);
    cmplx large;
    cmplx tmp;
    double large_mag, tmp_mag;
    uint32 k, m, n;
    n=0;
    for(k=0 ; k < ch; k++){
        large = datPtr[n++];
        large_mag = large.re*large.re + large.im*large.im;
        for(m = 1 ; m < rclen; m++, n++ ){
            tmp = datPtr[n];
            tmp_mag = tmp.re*tmp.re + tmp.im*tmp.im;
            if( large_mag < tmp_mag){
                large = datPtr[n];
                large_mag = tmp_mag;
            }
        }
        A.at<cmplx>(k) = large;
    }
    return A;
}
template <typename _T> Mat Mat::min() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch  = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(getDatType(),1,ch,1);
    _T less;
    uint32 k, m, n;
    n=0;
    for(k=0 ; k < ch; k++){
        less = datPtr[n++];
        for(m = 1 ; m < rclen; m++, n++ ){
            if( less > datPtr[n])
                less = datPtr[n];
        }
        A.at<_T>(k) = less;
    }
    return A;
}
template <> inline Mat Mat::min<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(getDatType(),1,ch,1);
    cmplx less;
    cmplx tmp;
    double less_mag, tmp_mag;
    uint32 k, m, n;
    n=0;
    for(k=0 ; k < ch; k++){
        less = datPtr[n++];
        less_mag = less.re*less.re + less.im*less.im;
        for(m = 1 ; m < rclen; m++, n++ ){
            tmp = datPtr[n];
            tmp_mag = tmp.re*tmp.re + tmp.im*tmp.im;
            if( less_mag > tmp_mag){
                less     = datPtr[n];
                less_mag = tmp_mag;
            }
        }
        A.at<cmplx>(k) = less;
    }
    return A;
}
template <typename _T> Mat Mat::mean() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(DTYP::DOUBLE,1,ch,1);
    double sum;
    uint32 k, m, n;
    n = 0;
    for(k=0 ; k < ch; k++){
        sum = 0;
        for(m = 0 ; m < rclen; m++){
            sum += double(datPtr[n++]);
        }
        A.at<double>(k) = sum/rclen;
    }
    return A;
}
template <> inline Mat Mat::mean<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat A(DTYP::CMPLX,1,ch,1);
    cmplx sum;
    cmplx tmp;
    uint32 k, m, n;
    n = 0;
    for(k=0 ; k < ch; k++){
        sum.zero();
        for(m = 0 ; m < rclen; m++){
            sum += cmplx(datPtr[n++]);
        }
        A.at<cmplx>(k) = sum/rclen;
    }
    return A;
}

template <typename _T> Mat Mat::std() {
    _T* datPtr = this->getDataPtr<_T>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat avg = mean<_T>();
    Mat A(DTYP::DOUBLE,1,ch,1);
    double sqsum;
    double diff;
    double ch_avg;
    uint32 k, m, n;
    n = 0;
    uint32 Div = rclen -1;
    for(k=0 ; k < ch; k++){
        sqsum  = 0;
        ch_avg = avg.at<double>(k);
        for(m = 0 ; m < rclen; m++){
           diff = double(datPtr[n++]) - ch_avg;
           sqsum += diff*diff;
        }
        A.at<double>(k) = sqrt(sqsum / Div);
    }
    return A;
}
template <> inline Mat Mat::std<cmplx>() {
    cmplx* datPtr = this->getDataPtr<cmplx>();

    uint32 ch    = getChannel();
    uint32 rclen = getRowColSize();

    Mat avg = mean<cmplx>();
    Mat A(DTYP::CMPLX,1,ch,1);
    cmplx sqsum;
    cmplx diff;
    cmplx ch_avg;
    uint32 k, m, n;
    n = 0;
    uint32 Div = rclen -1;
    for(k=0 ; k < ch; k++){
        sqsum  = 0;
        ch_avg = avg.at<cmplx>(k);
        for(m = 0 ; m < rclen; m++){
           diff = cmplx(datPtr[n++]) - ch_avg;
           sqsum += diff*diff;
        }
        sqsum /= Div;
        A.at<cmplx>(k) = sqsum;
    }
    return A;
}

} // namespace jmat

#endif // JBMAT_H

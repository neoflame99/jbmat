#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>

#include <initializer_list>
#include <vector>
#include <float.h>
#include <assert.h>
#include "types.h"
#include <math.h>
#include <string.h>

#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

namespace jmat {

typedef union _elemptr{
    uchar*  uch_ptr;
    int*    int_ptr;
    float*  f32_ptr;
    double* f64_ptr;
    cmplx*  cmx_ptr;
}elemptr;

typedef struct _matRect{
    int32 sR, sC;
    int32 eR, eC;
    _matRect(int32 sr=0, int32 sc=0, int32 er=0, int32 ec=0 ):sR(sr),sC(sc),eR(er),eC(ec){};
    _matRect(std::initializer_list<int32> list){
        std::vector<int32> v;
        v.push_back(0);
        v.push_back(0);
        v.push_back(0);
        v.push_back(0);
        v.insert(v.begin(),list.begin(),list.end());
        sR = v.at(0);
        sC = v.at(1);
        eR = v.at(2);
        eC = v.at(3);
     };
    inline void set(int32 sr=0, int32 sc=0, int32 er=0, int32 ec=0){
       sR= sr; sC=sc; eR= er; eC=ec;
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

private:
    inline void sync_data_ptr(){ elptr.uch_ptr = dat_ptr = mA.get(); }
    void alloc(uint32 len);
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
    Mat( std::initializer_list<int32 > list );
    Mat( std::initializer_list<float > list );
    Mat( std::initializer_list<uchar > list );
    Mat( std::initializer_list<cmplx > list );
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

    Mat& operator= (Mat&& other );
    Mat& operator+=(const Mat& other );
    Mat& operator+=(const double scalar);
    Mat& operator-=(const Mat& other );
    Mat& operator-=(const double scalar);
    Mat& operator*=(const Mat& other );
    Mat& operator*=(const double scalar);
    Mat& operator/=(const Mat& other );
    Mat& operator/=(const double scalar);

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
    inline uint32  getByteLen() const { return byteLen; }

    int32   reshape(uint32 r, uint32 c, uint32 ch=1);
    void    transpose();
    void    changeDType(const DTYP dt);
    void    printMat() ;
    void    printMat(const std::string objname);

public : // static methods
    static Mat ones (uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static Mat zeros(uint32 r, uint32 c, uint32 ch= 1, DTYP dt = DTYP::DOUBLE);
    static int32 instant_count;
    static Mat repeat(const Mat& src, const uint32 rp_r, const uint32 rp_c, const uint32 rp_ch);
    static int32 sliceCopyMat(const Mat& src, const matRect& srcSlice,const Mat& des, const matRect& desSlice );

public : // public template methods
    template <typename _T> _T& at(uint32 i=0) const;
    template <typename _T> _T& at(uint32 r, uint32 c, uint32 nch=0) const;
    template <typename _T=uchar> _T* getDataPtr() const;
    template <typename _T> Mat _max() ;
    template <typename _T> Mat _min() ;
    template <typename _T> Mat _mean();
    template <typename _T> Mat _std() ;
    template <typename _T> Mat _sum() ;
    Mat max();
    Mat min() ;
    Mat mean();
    Mat std() ;
    Mat sum() ;

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
    template <typename _T> static Mat _repeat(const Mat& src, const uint32 r, const uint32 c, const uint32 ch);
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
    //return static_cast<_T*>(dat_ptr);
    return (_T*)dat_ptr;
}
template <> inline cmplx* Mat::getDataPtr() const{
    return reinterpret_cast<cmplx*>(dat_ptr);
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
template <> inline void Mat::_print<cmplx>(cmplx* mdat){
    const int32 bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;
    uint32 i,j, k;

    if(datT==DTYP::CMPLX){
        cmplx val;
        for( i = 0; i < lenRowCol; i += col){ // rows
            snprintf(buf,bufsz,"[");
            for( j=0; j < col; j++){          // columns
                snprintf(tmp,bufsz," (");
                strncat(buf,tmp,bufsz);
                for( k = 0 ; k < Nch; k++){
                    val = mdat[i+j+k];
                    if( val.re >= neg_max_double && val.re <= pos_min_double)
                        val.re = 0.0;
                    if( val.im >= neg_max_double && val.im <= pos_min_double)
                        val.im = 0.0;

                    snprintf(tmp,bufsz,"% 5.2g %+5.2gi",val.re, val.im);
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
    cmplx  large, tmp;
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
    cmplx  large, tmp;
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
    cmplx sum;
    cmplx tmp;
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
    for(k=0; k < ch; ++k)
        A.at<double>(k) = sqrt(A.at<double>(k) / Div);
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
        A.at<double>(k) = std::sqrt(A.at<double>(k));

    return A;
}


template <typename _T> Mat Mat::_repeat(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    // src must have only 1 channel
    uint32 sr = src.getRow();
    uint32 sc = src.getCol();
    uint32 nr = sr * rr;
    uint32 nc = sc * rc;
    uint32 nch= rch;
    DTYP srcDtype = src.getDatType();

    Mat A(srcDtype, nr, nc, nch);
    _T* aDat_pt   = A.getDataPtr<_T>();
    _T* srcDat_pt = src.getDataPtr<_T>();

    uint32 x,y,z,i, sx, sy, sy_sc, k;
    i= 0;
    for ( z =0 ; z < nch ; ++z){
        for( y=0, sy=0; y < nr ; ++y, sy = y % sr){
            sy_sc = sy*sc;
            for( x=0, sx=0; x < nc ; ++x, sx = x % sc){
                k = sy_sc + sx;
                aDat_pt[i++] = srcDat_pt[k];
            }
        }
    }
    return A;
}

template <> inline Mat Mat::_repeat<cmplx>(const Mat& src, const uint32 rr, const uint32 rc, const uint32 rch){
    // src must have only 1 channel
    uint32 sr = src.getRow();
    uint32 sc = src.getCol();
    uint32 nr = sr * rr;
    uint32 nc = sc * rc;
    uint32 nch= rch;
    DTYP srcDtype = src.getDatType();

    Mat A(srcDtype, nr, nc, nch);
    cmplx* aDat_pt   = A.getDataPtr<cmplx>();
    cmplx* srcDat_pt = src.getDataPtr<cmplx>();

    uint32 x,y,z,i, sx, sy, sy_sc, k;
    i= 0;

    for ( z =0 ; z < nch ; ++z){
        for( y=0, sy=0; y < nr ; ++y, sy = y % sr){
            sy_sc = sy*sc;
            for( x=0, sx=0; x < nc ; ++x, sx = x % sc){
                k = sy_sc + sx;
                aDat_pt[i++] = srcDat_pt[k];
            }
        }
    }
    return A;
}

} // namespace jmat

#endif // JBMAT_H

#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>
#include <memory>
#include <initializer_list>
#include <vector>
#include <float.h>

#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

//-- shallow copy version & using shared_ptr
typedef unsigned char uchar;
typedef unsigned int uint;

typedef std::shared_ptr<uchar> shr_ptr;
enum DTYP {UCHAR, INT, FLOAT, DOUBLE};

class jbMat{
private:    
    shr_ptr mA;
    uchar *dat_ptr;
    DTYP datT;

    uint length;
    uint lenRowCol;
    void alloc(uint len);
    uint row, col, Nch;
    uint byteStep;
    uint byteLen;

    std::string obj_name;

private:
    void sync_data_ptr(){
        if( mA.get() != dat_ptr)
            dat_ptr = mA.get();
    }
    void init(uint r, uint c, uint ch, DTYP dt);
    void initName();


public:

    jbMat(DTYP dt = DTYP::DOUBLE);
    jbMat(DTYP dt, uint r, uint c, uint ch );
    jbMat(DTYP dt, uint r, uint c, uint ch, std::string name);
    jbMat(DTYP dt, uint rc);
    jbMat(const jbMat& mat);
    jbMat( std::initializer_list<double> list );
    jbMat( std::initializer_list<int> list );
    jbMat( std::initializer_list<float> list );
    ~jbMat();

    shr_ptr getMat() const { return mA; }
    bool isEmpty() const { return ((length <= 0) ? true : false); }
    void setRowCol(uint r, uint c, uint ch=1);
    jbMat copy() const;
    jbMat copyChannelN(const uint NoCh=0) const;
    jbMat copySubMat(const uint startRow, const uint endRow, const uint startCol, const uint endCol) const;
    void setChannelN(const jbMat& src, const uint srcfromCh=0,const uint Channels=1, const uint tarToCh=0);
    void setName(std::string name);

    static jbMat ones (uint r, uint c, uint ch= 1, DTYP dt = DTYP::DOUBLE);
    static jbMat zeros(uint r, uint c, uint ch= 1, DTYP dt = DTYP::DOUBLE);
    static int instant_count;


    //-- overloading operators
    jbMat operator+(const jbMat& other);
    jbMat operator+(const double scalar);
    jbMat& operator+=(const jbMat& other);
    jbMat& operator+=(const double scalar);
    jbMat operator-(const jbMat& other);
    jbMat operator-(const double scalar);
    jbMat& operator-=(const jbMat& other);
    jbMat& operator-=(const double scalar);
    jbMat operator*(const jbMat& other);
    jbMat operator*(const double scalar);
    jbMat& operator*=(const jbMat& other);
    jbMat& operator*=(const double scalar);
    jbMat operator/(const jbMat& other);
    jbMat operator/(const double scalar);
    jbMat& operator/=(const jbMat& other);
    jbMat& operator/=(const double scalar);
    //jbMat& operator=(jbMat other);
    jbMat& operator=(const jbMat& other);


    template <typename _T> _T& at(uint i) const;
    template <typename _T> _T& at(uint r, uint c, uint nch=0) const;
    template <typename _T> _T* getDataPtr() const;

    //double& operator[](int i) const;
    //double& operator()(int i) const;
    //double& operator()(int r, int c) const;
    //double& operator()(int r, int c, int ch) const;

    jbMat& plusMat    (const jbMat& other);
    jbMat& minusMat   (const jbMat& other);
    jbMat& multiplyMat(const jbMat& other);
    jbMat& divideMat  (const jbMat& other);

    jbMat& plusScalar(const double scalar);
    jbMat& plusScalar(const  float scalar);
    jbMat& plusScalar(const    int scalar);
    jbMat& plusScalar(const  uchar scalar);
    jbMat& minusScalar(const double scalar);
    jbMat& minusScalar(const  float scalar);
    jbMat& minusScalar(const    int scalar);
    jbMat& minusScalar(const  uchar scalar);
    jbMat& multiplyScalar(const double scalar);
    jbMat& multiplyScalar(const  float scalar);
    jbMat& multiplyScalar(const    int scalar);
    jbMat& multiplyScalar(const  uchar scalar);
    jbMat& divideScalar(const double scalar);
    jbMat& divideScalar(const  float scalar);
    jbMat& divideScalar(const    int scalar);
    jbMat& divideScalar(const  uchar scalar);


    uint getLength() const{ return length; }
    uint getRow() const { return row; }
    uint getCol() const { return col; }
    uint getChannel() const { return Nch; }
    DTYP getDatType() const { return datT; }
    uint getByteStep() const { return byteStep; }

    uint reshape(uint r, uint c, uint ch=1);
    void transpose();
    void changeDType(const DTYP dt);
    void printMat() ;
    void printMat(const std::string objname);

private:
    template <typename _T> void _print(_T* mdat);
    template <typename _T> void _plus_mat    (_T* self, _T* other, uint len);
    template <typename _T> void _minus_mat   (_T* self, _T* other, uint len);
    template <typename _T> void _multiply_mat(_T* self, _T* other, uint len);
    template <typename _T> void _divide_mat  (_T* self, _T* other, uint len);
    template <typename _T1, typename _T2> void _plus_scalar(_T1* self, _T2 scalar, uint len );
    template <typename _T1, typename _T2> void _minus_scalar(_T1* self, _T2 scalar, uint len );
    template <typename _T1, typename _T2> void _multiply_scalar(_T1* self, _T2 scalar, uint len );
    template <typename _T1, typename _T2> void _divide_scalar(_T1* self, _T2 scalar, uint len );
    template <typename _Ts, typename _Tt> void _type_change();
};


template <typename _T> _T& jbMat::at(uint i) const {
    if(isEmpty()) return *((_T *)nullptr); //*((_T *)mA.get()); //*mA;

    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    return ((_T *)dat_ptr)[i]; //return ((_T *)mA.get())[i];
}
template <typename _T> _T& jbMat::at(uint r, uint c, uint ch) const {
    if(isEmpty()) return *((_T *)nullptr);
    uint i = ch*lenRowCol + r*col + c;
    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    return ((_T *)dat_ptr)[i];
}

template <typename _T> _T* jbMat::getDataPtr() const {
    return (_T *)dat_ptr;
}

template <typename _T> void jbMat::_plus_mat(_T* self, _T* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] += other[k];
}
template <typename _T> void jbMat::_minus_mat(_T* self, _T* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] -= other[k];
}
template <typename _T> void jbMat::_multiply_mat(_T* self, _T* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] *= other[k];
}
template <typename _T> void jbMat::_divide_mat(_T* self, _T* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] /= other[k];
}

template <typename _T1, typename _T2> void jbMat::_plus_scalar(_T1* self, _T2 scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] += scalar;
}
template <typename _T1, typename _T2> void jbMat::_minus_scalar(_T1* self, _T2 scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] -= scalar;
}
template <typename _T1, typename _T2> void jbMat::_multiply_scalar(_T1* self, _T2 scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] *= scalar;
}
template <typename _T1, typename _T2> void jbMat::_divide_scalar(_T1* self, _T2 scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] /= scalar;
}

template <typename _T> void jbMat::_print(_T* mdat){
    const int bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    int i,j;

    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;

    double val;
    uint k;
    uint ch_offset;

    if(datT==DTYP::DOUBLE || datT==DTYP::FLOAT){
        for( k = 0 ; k < Nch; k++){
            fprintf(stdout,"channel: %d \n",k);
            ch_offset = k*lenRowCol;
            for( i = 0; i < lenRowCol; i += col){ // rows
                snprintf(buf,bufsz,"[");
                for( j=0; j < col; j++){          // columns
                    val = mdat[i+j+ch_offset];
                    if( val >= neg_max_double && val <= pos_min_double)
                        val = 0.0;
                    snprintf(tmp,bufsz," %9.3f ",val);
                    strncat(buf,tmp,bufsz);
                }
                strncat(buf,"]",1);
                fprintf(stdout,"%s\n",buf);
            }
        }
    }else{
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

template <typename _Ts, typename _Tt> void jbMat::_type_change(){
    _Ts* src_pt = (_Ts *)mA.get();
    _Tt* tar_pt;
    try{
        tar_pt = new _Tt[static_cast<unsigned long>(length)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Transpose Error: %s\n",ex.what());
        return;
    }

    uint k;
    for(k=0; k < length; k++)
        tar_pt[k] = src_pt[k];

    mA.reset((uchar *)tar_pt,std::default_delete<uchar[]>());
    sync_data_ptr();

    byteStep = sizeof(_Tt);
    byteLen  = length * byteStep;
}
#endif // JBMAT_H

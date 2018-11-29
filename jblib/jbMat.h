#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>
#include <memory>
#include <initializer_list>
#include <vector>
#include <float.h>
#include <assert.h>

#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

//-- shallow copy version & using shared_ptr
typedef unsigned char uchar;
typedef unsigned int  uint;

typedef std::shared_ptr<uchar> shr_ptr;
enum    DTYP {UCHAR=0 , INT=1 , FLOAT=2, DOUBLE=3};
struct rawMat{
    uchar* dat_ptr;
    uint rows;
    uint cols;
    uint channels;
    DTYP dtype;
};
class jbMat{
private: // member fields
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


public: // constructors and destructor

    jbMat(DTYP dt = DTYP::DOUBLE);
    jbMat(DTYP dt, uint r, uint c, uint ch );
    jbMat(DTYP dt, uint r, uint c, uint ch, std::string name);
    jbMat(DTYP dt, uint rc);
    jbMat(shr_ptr ma, DTYP dt, uint r, uint c, uint ch );
    jbMat(const jbMat& mat);
    jbMat( std::initializer_list<double> list );
    jbMat( std::initializer_list<int> list );
    jbMat( std::initializer_list<float> list );
    ~jbMat();

public:
    //-- overloading operators
    jbMat  operator+(const jbMat& other ) const;
    jbMat  operator+(const double scalar) const;
    jbMat  operator-(const jbMat& other ) const;
    jbMat  operator-(const double scalar) const;
    jbMat  operator*(const jbMat& other ) const;
    jbMat  operator*(const double scalar) const;
    jbMat  operator/(const jbMat& other ) const;
    jbMat  operator/(const double scalar) const;

    jbMat& operator= (const jbMat& other );
    jbMat& operator+=(const jbMat& other );
    jbMat& operator+=(const double scalar);
    jbMat& operator-=(const jbMat& other );
    jbMat& operator-=(const double scalar);
    jbMat& operator*=(const jbMat& other );
    jbMat& operator*=(const double scalar);
    jbMat& operator/=(const jbMat& other );
    jbMat& operator/=(const double scalar);

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

    void    setRowCol(uint r, uint c, uint ch=1);
    void    setChannelN(const jbMat& src, const uint srcfromCh=0,const uint Channels=1, const uint tarToCh=0);
    void    setName(std::string name);

    jbMat   copy() const;
    jbMat   copyChannelN(const uint NoCh=0) const;
    jbMat   copySubMat(const uint startRow, const uint endRow, const uint startCol, const uint endCol) const;

    rawMat  getRawMat() const;
    inline bool isEmpty() const { return ((length <= 0) ? true : false); }
    inline shr_ptr getMat () const { return mA; }
    inline uint    getLength() const{ return length; }
    inline uint    getRow() const { return row; }
    inline uint    getCol() const { return col; }
    inline uint    getChannel() const { return Nch; }
    inline DTYP    getDatType() const { return datT; }
    inline uint    getByteStep() const { return byteStep; }
    uint    reshape(uint r, uint c, uint ch=1);
    void    transpose();
    void    changeDType(const DTYP dt);
    void    printMat() ;
    void    printMat(const std::string objname);

public : // static methods
    static jbMat ones (uint r, uint c, uint ch= 1, DTYP dt = DTYP::DOUBLE);
    static jbMat zeros(uint r, uint c, uint ch= 1, DTYP dt = DTYP::DOUBLE);
    static int instant_count;

public : // public template methods
    template <typename _T> inline _T& at(uint i) const;
    template <typename _T> inline _T& at(uint r, uint c, uint nch=0) const;
    template <typename _T> inline _T* getDataPtr() const;

private: // private template methods
    template <typename _T> void _print(_T* mdat);
    template <typename _Tsrc, typename _Ttar> void _type_change();
    template <typename _Tslf, typename _Totr> void _plus_mat       (_Tslf* self, _Totr* other, uint len );
    template <typename _Tslf, typename _Totr> void _minus_mat      (_Tslf* self, _Totr* other, uint len );
    template <typename _Tslf, typename _Totr> void _multiply_mat   (_Tslf* self, _Totr* other, uint len );
    template <typename _Tslf, typename _Totr> void _divide_mat     (_Tslf* self, _Totr* other, uint len );
    template <typename _Tslf, typename _Totr> void _plus_scalar    (_Tslf* self, _Totr scalar, uint len );
    template <typename _Tslf, typename _Totr> void _minus_scalar   (_Tslf* self, _Totr scalar, uint len );
    template <typename _Tslf, typename _Totr> void _multiply_scalar(_Tslf* self, _Totr scalar, uint len );
    template <typename _Tslf, typename _Totr> void _divide_scalar  (_Tslf* self, _Totr scalar, uint len );

};


template <typename _T> inline _T& jbMat::at(uint i) const {
    if(isEmpty()) return *((_T *)nullptr); //*((_T *)mA.get()); //*mA;

    assert(i < length);
    if(i >= length){
        //fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    return ((_T *)dat_ptr)[i]; //return ((_T *)mA.get())[i];
}
template <typename _T> inline _T& jbMat::at(uint r, uint c, uint ch) const {
    if(isEmpty()) return *((_T *)nullptr);
    uint i = ch*lenRowCol + r*col + c;
    assert(i < length);
    if(i >= length){
        fprintf(stderr,"The Index of jbMat is out of bound\n");
        i = length-1;
    }
    return ((_T *)dat_ptr)[i];
}

template <typename _T> inline _T* jbMat::getDataPtr() const {
    return (_T *)dat_ptr;
}

template <typename _Tslf, typename _Totr> void jbMat::_plus_mat(_Tslf* self, _Totr* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] += other[k];
}
template <typename _Tslf, typename _Totr> void jbMat::_minus_mat(_Tslf* self, _Totr* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] -= other[k];
}
template <typename _Tslf, typename _Totr> void jbMat::_multiply_mat(_Tslf* self, _Totr* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] *= other[k];
}
template <typename _Tslf, typename _Totr> void jbMat::_divide_mat(_Tslf* self, _Totr* other, uint len){
    for(uint k=0; k < len; k++ )
        self[k] /= other[k];
}

template <typename _Tslf, typename _Totr> void jbMat::_plus_scalar(_Tslf* self, _Totr scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] += scalar;
}
template <typename _Tslf, typename _Totr> void jbMat::_minus_scalar(_Tslf* self, _Totr scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] -= scalar;
}
template <typename _Tslf, typename _Totr> void jbMat::_multiply_scalar(_Tslf* self, _Totr scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] *= scalar;
}
template <typename _Tslf, typename _Totr> void jbMat::_divide_scalar(_Tslf* self, _Totr scalar, uint len){
    for(uint k=0; k < len ; k++)
        self[k] /= scalar;
}

template <typename _Tsrc, typename _Ttar> void jbMat::_type_change(){
    _Tsrc* src_pt = (_Tsrc *)mA.get();
    _Ttar* tar_pt;
    try{
        tar_pt = new _Ttar[static_cast<unsigned long>(length)];
    }catch(std::bad_alloc& ex){
        fprintf(stderr,"Transpose Error: %s\n",ex.what());
        return;
    }

    uint k;
    for( k=0; k < length; k++){
        tar_pt[k] = (_Ttar) (src_pt[k]);
    }
    mA.reset((uchar *)tar_pt,std::default_delete<uchar[]>());
    sync_data_ptr();

    byteStep = sizeof(_Ttar);
    byteLen  = length * byteStep;
}

template <typename _T> void jbMat::_print(_T* mdat){
    const int bufsz = 2049;
    char buf[bufsz]="\0";
    char tmp[bufsz];
    uint i,j;

    const double neg_max_double = -DBL_EPSILON ;
    const double pos_min_double =  DBL_EPSILON ;

    uint k;
    uint ch_offset;

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
        int val;
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


#endif // JBMAT_H

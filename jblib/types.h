#ifndef TYPES_H
#define TYPES_H
#include <memory>
#include <cstdint>
#include <float.h>

#define U08 uint8_t
#define I32 int32_t
#define U32 uint32_t
#define I64 int64_t
#define U64 uint64_t
#define F32 float
#define F64 double

namespace jmat {    
    typedef unsigned char  uchar;
    typedef uint8_t        uint8;
    typedef int8_t          int8;
    typedef uint32_t      uint32;
    typedef int32_t        int32;
    typedef uint64_t      uint64;
    typedef int64_t        int64;
    typedef std::shared_ptr<uchar> shr_ptr;
    enum class DTYP {UCHAR=0 , INT=1 , FLOAT=2, DOUBLE=3, CMPLX=4};
    struct _complex;
    typedef _complex  cmplx;

    struct _complex{
        double re;
        double im;
        // constructor
        _complex(double r):re(r),im(0){} // type conversion constructor
        _complex(double r, double i):re(r),im(i){}
        _complex():re(0),im(0){}
        _complex(const _complex& ot):re(ot.re),im(ot.im){} // copy constructor
        _complex(const _complex&& ot):re(ot.re),im(ot.im){} // move constructor

        // operator override
        friend inline _complex operator+(_complex& lhs, double rhs);
        friend inline _complex operator+(double lhs, _complex& rhs);
        //friend inline _complex operator+(const _complex& lhs, const _complex& rhs);
        inline _complex operator+(const _complex& rhs);

        friend inline _complex operator-(_complex& lhs, double rhs);
        friend inline _complex operator-(double lhs, _complex& rhs);
        inline _complex operator-(const _complex& rhs);

        friend inline _complex operator*(_complex& lhs, double rhs);
        friend inline _complex operator*(double lhs, _complex& rhs);
        friend inline _complex operator*(const _complex& lhs, const _complex& rhs);
        inline _complex operator*(const _complex& rhs);

        friend inline _complex operator/(_complex& lhs, double rhs);
        friend inline _complex operator/(double lhs, _complex& rhs);
        friend inline _complex operator/(const _complex& lhs, const _complex& rhs);

        friend inline double operator+=(const double lhs, const _complex& rhs);
        inline _complex& operator+=(const _complex& rhs);
        inline _complex& operator+=(const double rhs);

        friend inline double operator-=(const double lhs, const _complex& rhs);
        inline _complex& operator-=(const _complex& rhs);
        inline _complex& operator-=(const double rhs);

        friend inline double operator*=(const double lhs, const _complex& rhs);
        inline _complex& operator*=(const _complex& rhs);
        inline _complex& operator*=(const double rhs);

        friend inline double operator/=(const double lhs, const _complex& rhs);
        inline _complex& operator/=(const _complex& rhs);
        inline _complex& operator/=(const double rhs);


        inline _complex& operator=(const _complex& rhs) {  re = rhs.re; im = rhs.im; return *this; }
        inline _complex& operator=(const double rhs) { re=rhs; im=0; return *this;}
        inline _complex& operator=(const _complex&& rhs) { *this=rhs; return *this;} // move assignment

        friend bool operator==(const _complex& lhs,const double rhs);
        friend bool operator==(const double lhs, const _complex& rhs);
        inline bool operator==(const _complex& rhs) { return (re<=rhs.re+DBL_EPSILON && re >=rhs.re-DBL_EPSILON && im<=rhs.im+DBL_EPSILON && im >=rhs.im-DBL_EPSILON)? true : false; }

        inline void set_val(double r=0, double i=0) { re= r; im = i; }
        inline void zero() {set_val(0,0);}
        //inline void conj() {im=-im; }
        static inline _complex conj(const _complex& A) { return _complex(A.re, -A.im); }
        inline _complex conj() { return _complex(re, -im); }

    };

    inline _complex operator +(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re += rhs;
        return A;
    }
    inline _complex operator +(double lhs , _complex& rhs){
        _complex A ;
        A.re = lhs + rhs.re;
        return A;
    }
    inline _complex _complex::operator +(const _complex& rhs){
        _complex A ;
        A.re = re + rhs.re;
        A.im = im + rhs.im;
        return A;
    }
    inline _complex operator -(_complex& lhs, double rhs){
        _complex A;
        A.re = lhs.re - rhs;
        return A;
    }
    inline _complex operator -(double lhs , _complex& rhs){
        _complex A;
        A.re = lhs - rhs.re;
        A.im = -rhs.im;
        return A;
    }
    inline _complex _complex::operator -(const _complex& rhs){
        _complex A;
        A.re = re - rhs.re;
        A.im = im - rhs.im;
        return A;
    }

    inline _complex operator *(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re *= rhs;
        A.im *= rhs;
        return A;
    }
    inline _complex operator *(double lhs , _complex& rhs){
        _complex A = rhs;
        A.re *= lhs;
        A.im *= lhs;
        return A;
    }
    inline _complex operator *(const _complex& lhs ,const _complex& rhs){
        _complex A;
        A.re = lhs.re * rhs.re - lhs.im * rhs.im;
        A.im = lhs.re * rhs.im + lhs.im * rhs.re;
        return A;
    }
    inline _complex _complex::operator *(const _complex& rhs){
        _complex A = rhs;
        A.re = rhs.re * re - rhs.im * im;
        A.im = rhs.re * im + rhs.im * re;
        return A;
    }
    inline _complex operator /(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re /= rhs;
        A.im /= rhs;
        return A;
    }
    inline _complex operator /(double lhs , _complex& rhs){
        _complex A(rhs.re, -rhs.im);

        return ( lhs * A) /( rhs * A).re;
    }
    inline _complex operator /(const _complex& lhs, const _complex& rhs){
        _complex C(rhs.re, -rhs.im);
        _complex N = lhs * C;
        double   d = (rhs * C).re;
        return N / d;
    }

    inline double operator+=(const double lhs, const _complex& rhs){
        return lhs + rhs.re;
    }
    inline _complex& _complex::operator+=(const _complex& rhs){
        re += rhs.re;
        im += rhs.im;
        return *this;
    }
    inline _complex& _complex::operator+=(const double rhs){
        re += rhs;
        return *this;
    }
    inline double operator-=(const double lhs, const _complex& rhs){
        return lhs - rhs.re;
    }
    inline _complex& _complex::operator -=(const _complex& rhs){
        re -= rhs.re;
        im -= rhs.im;
        return (*this);
    }
    inline _complex& _complex::operator -=(const double rhs){
        re -= rhs;
        return (*this);
    }
    inline double operator*=(const double lhs, const _complex& rhs){
        return lhs * rhs.re;
    }
    inline _complex& _complex::operator *=(const _complex& rhs){
        *this = *this * rhs;

        return (*this);
    }
    inline _complex& _complex::operator *=(const double rhs){
        re *= rhs;
        im *= rhs;
        return (*this);
    }
    inline double operator/=(const double lhs, const _complex& rhs){
        return lhs / rhs.re;
    }
    inline _complex& _complex::operator /=(const _complex& rhs){
        *this = *this / rhs;

        return (*this);
    }
    inline _complex& _complex::operator /=(const double rhs){
        re /= rhs;
        im /= rhs;
        return (*this);
    }

    inline bool operator==(const _complex& lhs, const double rhs){
        return (lhs.re <= rhs+DBL_EPSILON && lhs.re >= rhs-DBL_EPSILON && lhs.im == 0.0) ? true : false;
    }
    inline bool operator==(const double lhs, const _complex& rhs){
        return (rhs.re <= lhs+DBL_EPSILON && rhs.re >= lhs-DBL_EPSILON && rhs.im == 0.0) ? true : false;
    }

}

namespace jmat {
static const double bt601_r2y[3][3]={{ 0.299    ,  0.587    ,  0.114 },
                                    {-0.168736 , -0.331264 ,  0.5   },
                                    { 0.5      , -0.418688 , -0.081312} };
static const double bt709_r2y[3][3]={{ 0.2126  ,  0.7152  , 0.0722  },
                                    {-0.11457 , -0.38543 , 0.5     },
                                    { 0.5     , -0.45415 ,-0.04585 }};

static const double bt601_y2r[3][3]={{ 1.0000 ,  -0.0000 ,   1.4020 },
                                    { 1.0000 ,  -0.3441 ,  -0.7141 },
                                    { 1.0000 ,   1.7720 ,   0.0000 }};

static const double bt709_y2r[3][3]={{ 1.0000 ,   0.0000 ,   1.5748 },
                                    { 1.0000 ,  -0.1873 ,  -0.4681 },
                                    { 1.0000 ,   1.8556 ,  -0.0000 }};
static const double rgb2xyz_bt709[3][3] = { {0.4124564, 0.3575761, 0.1804375},
                                           {0.2126729, 0.7151522, 0.0721750},
                                           {0.0193339, 0.1191920, 0.9503041} };
static const double xyz2rgb_bt709[3][3] = { { 3.2404542, -1.5371385, -0.4985314},
                                           {-0.9692660,  1.8760108,  0.0415560},
                                           { 0.0556434, -0.2040259,  1.0572252} };
}
#endif // TYPES_H

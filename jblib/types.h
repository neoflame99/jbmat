#ifndef TYPES_H
#define TYPES_H
#include <memory>
#include <cstdint>

namespace jmat {    
    typedef unsigned char  uchar;
    typedef uint8_t        uint8;
    typedef int8_t          int8;
    typedef uint32_t      uint32;
    typedef int32_t        int32;
    typedef uint64_t      uint64;
    typedef int64_t        int64;
    typedef std::shared_ptr<uchar> shr_ptr;
    enum    DTYP {UCHAR=0 , INT=1 , FLOAT=2, DOUBLE=3, CMPLX=4};    
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
        inline _complex operator+(const _complex& rhs);


        friend inline _complex operator-(_complex& lhs, double rhs);
        friend inline _complex operator-(double lhs, _complex& rhs);
        inline _complex operator-(const _complex& rhs);

        friend inline _complex operator*(_complex& lhs, double rhs);
        friend inline _complex operator*(double lhs, _complex& rhs);
        inline _complex operator*(const _complex& rhs);

        friend inline _complex operator/(_complex& lhs, double rhs);
        friend inline _complex operator/(double lhs, _complex& rhs);
        inline _complex operator/(const _complex& rhs);


        friend inline double operator+=(double lhs, _complex& rhs);
        inline _complex& operator+=(const _complex& rhs);
        inline _complex& operator+=(const double rhs);

        friend inline double operator-=(double lhs, _complex& rhs);
        inline _complex& operator-=(const _complex& rhs);
        inline _complex& operator-=(const double rhs);

        friend inline double operator*=(double lhs, _complex& rhs);
        inline _complex& operator*=(const _complex& rhs);
        inline _complex& operator*=(const double rhs);

        friend inline double operator/=(double lhs, _complex& rhs);
        inline _complex& operator/=(const _complex& rhs);
        inline _complex& operator/=(const double rhs);


        inline _complex& operator=(const _complex& rhs) {  re = rhs.re; im = rhs.im; return *this; }
        inline _complex& operator=(const double rhs) { re=rhs; im=0; return *this;}
        inline _complex& operator=(const _complex&& rhs) { *this=rhs; return *this;} // move assignment

        friend bool operator==(const _complex& lhs,const double rhs);
        friend bool operator==(const double lhs, const _complex& rhs);
        inline bool operator==(const _complex& rhs) { return (re==rhs.re && im==rhs.im)? true : false; }

        inline void set_val(double r=0, double i=0) { re= r; im = i; }
        inline void zero() {set_val(0,0);}
        inline void conj() {im=-im; }
    };

    inline _complex operator +(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re += rhs;
        return A;
    }
    inline _complex operator +(double lhs , _complex& rhs){
        _complex A = rhs;
        A.re += lhs;
        return A;
    }
    inline _complex _complex::operator +(const _complex& rhs){
        _complex A = rhs;
        A.re += re;
        A.im += im;
        return A;
    }


    inline _complex operator -(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re -= rhs;
        return A;
    }
    inline _complex operator -(double lhs , _complex& rhs){
        _complex A = rhs;
        A.re -= lhs;
        return A;
    }
    inline _complex _complex::operator -(const _complex& rhs){
        _complex A = rhs;
        A.re -= re;
        A.im -= im;
        return A;
    }

    inline _complex operator *(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re *= rhs;
        return A;
    }
    inline _complex operator *(double lhs , _complex& rhs){
        _complex A = rhs;
        A.re *= lhs;
        return A;
    }
    inline _complex _complex::operator *(const _complex& rhs){
        _complex A = rhs;
        A.re *= re;
        A.im *= im;
        return A;
    }

    inline _complex operator /(_complex& lhs, double rhs){
        _complex A = lhs;
        A.re /= rhs;
        return A;
    }
    inline _complex operator /(double lhs , _complex& rhs){
        _complex A = rhs;
        A.re /= lhs;
        return A;
    }
    inline _complex _complex::operator /(const _complex& rhs){
        _complex A = rhs;
        A.re /= re;
        A.im /= im;
        return A;
    }

    inline double operator+=(double lhs, _complex& rhs){
        lhs += rhs.re;
        return lhs;
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

    inline double operator-=(double lhs, _complex& rhs){
        lhs -= rhs.re;
        return lhs;
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

    inline double operator*=(double lhs, _complex& rhs){
        lhs *= rhs.re;
        return lhs;
    }
    inline _complex& _complex::operator *=(const _complex& rhs){
        re *= rhs.re;
        im *= rhs.im;
        return (*this);
    }
    inline _complex& _complex::operator *=(const double rhs){
        re *= rhs;
        return (*this);
    }

    inline double operator/=(double lhs, _complex& rhs){
        lhs /= rhs.re;
        return lhs;
    }
    inline _complex& _complex::operator /=(const _complex& rhs){
        re /= rhs.re;
        im /= rhs.im;
        return (*this);
    }
    inline _complex& _complex::operator /=(const double rhs){
        re /= rhs;
        return (*this);
    }

    inline bool operator==(const _complex& lhs, const double rhs){
        return (lhs.re == rhs && lhs.im == 0.0) ? true : false;
    }
    inline bool operator==(const double lhs, const _complex& rhs){
        return (rhs.re == lhs && rhs.im == 0.0) ? true : false;
    }

}
#endif // TYPES_H

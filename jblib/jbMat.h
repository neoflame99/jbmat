#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>
#include <memory>
#include <initializer_list>
#include <vector>

//-- shallow copy version & using shared_ptr
typedef std::shared_ptr<double> ptr_double;
class jbMat{
private:    
    ptr_double mA;

    int length;
    int lenRowCol;
    void alloc(int len);
    int row, col, Nch;

public:

    jbMat();
    jbMat(int r, int c, int ch=1);
    jbMat(int rc );
    jbMat(const jbMat& mat);
    jbMat( std::initializer_list<double> list );
    jbMat( std::initializer_list<int> list );
    jbMat( std::initializer_list<float> list );
    ~jbMat();
    //double *getMat() const { return mA; }
    ptr_double getMat() const { return mA; }
    bool isEmpty() const { return ((length <= 0) ? true : false); }
    void setRowCol(int r, int c, int ch=1);
    jbMat copy() const;
//    jbMat copyChannel(const int NoCh=-1);
//    jbMat copySubMat(const int r_start, const int r_end, const int c_start, const int c_end, const int ch_start, const int ch_end);
    static jbMat ones(int r, int c, int ch= 1);
    static jbMat zeros(int r, int c, int ch= 1);


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
    double& operator[](int i) const;
    double& operator()(int i) const;
    double& operator()(int r, int c) const;
    double& operator()(int r, int c, int ch) const;
    int getLength() const{ return length; }
    int getRow() const { return row; }
    int getCol() const { return col; }
    int getChannel() const { return Nch; }
    int reshape(int r, int c, int ch=1);
    void transpose();
    void printMat() const;
    void printMat(const std::string objname) const;
};

#endif // JBMAT_H

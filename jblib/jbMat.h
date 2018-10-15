#ifndef JBMAT_H
#define JBMAT_H

#include <stddef.h>
#include <memory>
#include <initializer_list>
#include <vector>
#ifdef _MACOS_
    #include <string>
#else
    #include <string.h>
#endif

//-- shallow copy version & using shared_ptr
typedef std::shared_ptr<double> ptr_double;
typedef unsigned int uint;

class jbMat{
private:    
    ptr_double mA;
    double *dat_ptr;

    int length;
    int lenRowCol;
    void alloc(int len);
    int row, col, Nch;

    std::string obj_name;

private:
    void sync_data_ptr(){
        if( mA.get() != dat_ptr)
            dat_ptr = mA.get();
    }
    void init(int r, int c, int ch);
    void initName();
public:


public:

    jbMat();
    jbMat(int r, int c, int ch=1);
    jbMat(int r, int c, int ch, std::string name);
    jbMat(int rc );
    jbMat(const jbMat& mat);
    jbMat( std::initializer_list<double> list );
    jbMat( std::initializer_list<int> list );
    jbMat( std::initializer_list<float> list );
    ~jbMat();

    ptr_double getMat() const { return mA; }
    bool isEmpty() const { return ((length <= 0) ? true : false); }
    void setRowCol(int r, int c, int ch=1);
    jbMat copy() const;
    jbMat copyChannelN(const uint NoCh=0) const;
    jbMat copySubMat(const uint startRow, const uint endRow, const uint startCol, const uint endCol) const;
    void setChannelN(const jbMat& src, const uint srcCh=0, const uint tarCh=0);
    void setChannelN(const jbMat& src, const uint srcfromCh=0,const uint Channels=1, const uint tarToCh=0);
    void setName(std::string name);

    static jbMat ones(int r, int c, int ch= 1);
    static jbMat zeros(int r, int c, int ch= 1);
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

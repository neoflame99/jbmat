#ifndef JBMATH_H
#define JBMATH_H
#include <stddef.h>
#include "jbMat.h"
#include <string>

class jbMath
{
public:
    jbMath();
    ~jbMath();
    static jbMat mulMatrix(const jbMat& mA,const jbMat& mB);
    static jbMat triu(const jbMat& mA);
    static jbMat tril(const jbMat& mA);
    static jbMat augmentMatrix(const jbMat& mA);
    static jbMat inverse(const jbMat& mA);
    static jbMat tranpose(const jbMat& mA);
    static void printMat(const jbMat& Mat);
    static jbMat conv2d(const jbMat& mA, const jbMat& mB, std::string opt_conv="" , std::string opt_out="");
};

#endif // JBMATH_H

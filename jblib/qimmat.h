#ifndef QIMMAT_H
#define QIMMAT_H

#include <QImage>
#include "jbMat.h"

class QimMat
{
public:
    QimMat();
    static jbMat qim2jbmat(const QImage& src);
    static QImage jbmat2qim(const jbMat& src);
};

#endif // QIMMAT_H

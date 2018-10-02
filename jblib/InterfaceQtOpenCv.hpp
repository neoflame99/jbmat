#ifndef INTERFACEQTOPENCV_HPP
#define INTERFACEQTOPENCV_HPP

//#define  _WIN32_

#include <QImage>
#include <QPixmap>
#include <QDebug>


#include <opencv2/opencv.hpp>



namespace IfQtOpencv{

    inline QImage cvMat2QImage(const cv::Mat &inMat, bool byCopy = true){
        switch ( inMat.type() ){
        case CV_8UC4:{ // 8-bit, 4 channel
            QImage image(inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB32 );
            if(byCopy){
                qDebug() << "cvmat2qimage_cpy: CV_8UC4";
                image = image.copy();
            }else{
                qDebug() << "cvmat2qimage_ref: CV_8UC4";
            }
            return image;
        }
        case CV_8UC3:{ // 8-bit, 3 channel
            QImage image( inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB888 );
            if(byCopy){
                qDebug() << "cvmat2qimage_cpy: CV_8UC3";
                image = image.copy();
            }else{
                qDebug() << "cvmat2qimage_ref: CV_8UC3";
            }
            return image.rgbSwapped();
        }
        case CV_8UC1:{ // 8-bit, 1 channel
            static QVector<QRgb>  sColorTable;
            // only create our color table once
            if ( sColorTable.isEmpty() ){
                for ( int i = 0; i < 256; ++i )
                    sColorTable.push_back( qRgb( i, i, i ) );
            }

            QImage image( inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_Indexed8 );
            if(byCopy){
                qDebug() << "cvmat2qimage_cpy: CV_8UC1";
                image = image.copy();
            }else{
                qDebug() << "cvmat2qimage_ref: CV_8UC1";
            }
            image.setColorTable( sColorTable );
            return image;
        }
        default:
            qWarning() << " cv::Mat image type not handled in switch:" << inMat.type();

            break;
        }

       return QImage();
    } // end of cvMat2QImage(...)

    inline QPixmap cvMat2QPixmap(const cv::Mat &MatIn, bool byCopy=true){
        return QPixmap::fromImage(cvMat2QImage(MatIn, byCopy));
    }


    inline cv::Mat QImage2CvMat( const QImage &inImage, bool inCloneImageData = true ){
       switch ( inImage.format() ){
          case QImage::Format_RGB32:
          case QImage::Format_ARGB32: { // 8-bit, 4 channel
             cv::Mat  mat( inImage.height(), inImage.width(), CV_8UC4, const_cast<uchar*>(inImage.bits()), inImage.bytesPerLine() );

             return (inCloneImageData ? mat.clone() : mat);
          }
          case QImage::Format_RGB888: { // 8-bit, 3 channel
             if ( !inCloneImageData )
                qWarning() << " Conversion requires cloning since we use a temporary QImage";

             QImage   swapped = inImage.rgbSwapped();

             return cv::Mat( swapped.height(), swapped.width(), CV_8UC3, const_cast<uchar*>(swapped.bits()), swapped.bytesPerLine() ).clone();
          }
          case QImage::Format_Indexed8:{ // 8-bit, 1 channel
             cv::Mat  mat( inImage.height(), inImage.width(), CV_8UC1, const_cast<uchar*>(inImage.bits()), inImage.bytesPerLine() );

             return (inCloneImageData ? mat.clone() : mat);
          }

          default:
             qWarning() << " QImage format not handled in switch:" << inImage.format();
             break;
       }

       return cv::Mat();
    }

    // If inPixmap exists for the lifetime of the resulting cv::Mat, pass false to inCloneImageData to share inPixmap's data
    // with the cv::Mat directly
    //    NOTE: Format_RGB888 is an exception since we need to use a local QImage and thus must clone the data regardless
    inline cv::Mat QPixmap2CvMat( const QPixmap &inPixmap, bool inCloneImageData = true )
    {
       return QImage2CvMat( inPixmap.toImage(), inCloneImageData );
    }


}



#endif // INTERFACEQTOPENCV_HPP

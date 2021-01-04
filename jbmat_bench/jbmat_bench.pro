QT += gui

CONFIG += c++17 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp \
    ../jblib/jbMat.cpp \
    ../jblib/jbmath.cpp \
    ../jblib/qimmat.cpp \
    ../jblib/jbimgproc.cpp \
    ../jblib/jbBmp.cpp

HEADERS += \
    ../jblib/types.h \
    ../jblib/satcast.h \
    ../jblib/jbMat.h \
    ../jblib/jbmath.h \
    ../jblib/qimmat.h \
    ../jblib/jbimgproc.h \
    ../jblib/jbBmp.h

win32{
    DEFINES += _WIN_
    message("Windows")
}unix:macx{
    DEFINES += _MACOS_
    message("Mac OS")
    DEFINES += _MACOS_
    INCLUDEPATH += "/usr/local/include/opencv4"
    CONFIG(debug, debug|release){
        LIBS += -L"/usr/local/lib" \
            -lopencv_core \
            -lopencv_highgui \
            -lopencv_imgproc \
            -lopencv_imgcodecs
    }else{
        LIBS += -L"/usr/local/lib" \
            -lopencv_core \
            -lopencv_highgui \
            -lopencv_imgproc \
            -lopencv_imgcodecs
    }
}unix:!macx{
    DEFINES += _LINUX_
    message("LINUX OS")
}

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += static
LIBS += -lz
#LIBS += -lm
#LIBS += -lcrypt

SOURCES += main.cpp \
    ../Depends/fasta_parser.cpp \
    clsalgorithm.cpp \
    clsscaffoldfiller.cpp \
    ../Depends/smith-waterman.cpp \
    ../Depends/local_alignment.cpp \
    ../Depends/KmerUtils.cpp \
    ../Depends/bam_parse.cpp \
    ../Depends/fai_parser.cpp \
    ../Depends/needleman_wunsch/alignment.c \
    ../Depends/needleman_wunsch/alignment_scoring.c \
    ../Depends/needleman_wunsch/needleman_wunsch.c \
    clscorestructure.cpp
INCLUDEPATH += ../Depends/
INCLUDEPATH += ../Depends/needleman_wunsch/
INCLUDEPATH += ../Bamtools/include/needleman_wunsch/
INCLUDEPATH += /usr/include/c++/4.8/bits/

HEADERS += \
    ../Depends/kseq.h \
    ../Depends/fasta_parser.h \
    clsalgorithm.h \
    clsscaffoldfiller.h \
    ../Depends/smith-waterman.h \
    ../Depends/stdaln.h \
    ../Depends/local_alignment.h \
    ../Depends/KmerUtils.h \
    ../Depends/bam_parse.h \
    ../Depends/fai_parser.h \
    ../Depends/bam_alignment_record.h \
    ../Depends/bam_header_record.h \
    ../Depends/needleman_wunsch/alignment.h \
    ../Depends/needleman_wunsch/alignment_scoring.h \
    ../Depends/needleman_wunsch/needleman_wunsch.h \
    ../Depends/needleman_wunsch/uthash.h \
    clscorestructure.h


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../Bamtools/lib/release/ -lbamtools
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../Bamtools/lib/debug/ -lbamtools
else:unix: LIBS += -L$$PWD/../Bamtools/lib/ -lbamtools

INCLUDEPATH += $$PWD/../Bamtools/include
DEPENDPATH += $$PWD/../Bamtools/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../Bamtools/lib/release/libbamtools.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../Bamtools/lib/debug/libbamtools.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../Bamtools/lib/release/bamtools.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../Bamtools/lib/debug/bamtools.lib
else:unix: PRE_TARGETDEPS += $$PWD/../Bamtools/lib/libbamtools.a


#OTHER_FILES += \
#    ../Bamtools/lib/libbamtools.so.2.3.0


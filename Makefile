#--------------------------------------------------------------------------
# следующие 3 строки лучше не трогать
ifeq (,$(WORKDIR))
WORKDIR=..
endif

# Всё что делается делается вслух :) - желательно оставить как есть
VERBOSEMAKE=""
# Если ключик -g при сборки не нужен - если уверен в себе
#NDEBUG=""

# Компилятор для *.с файлов
#CC_LOCAL := gcc
# Компилятор для *.cc и *.cpp файлов
#CXX_LOCAL := g++
# Компилятор для *.f, *.for и *.F файлов
F77_LOCAL := g77
# Линкер для этого пакета по умолчанию
#LD_LOCAL := g++
# Раскомментировать, если необходимо задавать линкер для каждого
# бинарника в отдельности. Обычно не требуется.
#LD_MULTI := ""

# Дополнительные опции для компилятора C/C++
COPTS = -I. `root-config --cflags` #-ansi -pedantic
NOOPT = ""
# Дополнительные опции для компилятора Fortran
#FOPTS = -I.
FOPTS  = -g -fvxt -Wall -fno-automatic -finit-local-zero \
-fno-second-underscore -ffixed-line-length-120 -Wno-globals \
-DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC -DCERNLIB_BLDLIB \
-DCOMPDATE="'$(COMPDATE)'" -I$(CERN)/pro/include -I$(CERN)/pro/include/geant321

# Дополнительные опции для линкера
LDOPTS = -Xlinker -rpath -Xlinker `root-config --libdir`

#Если определена переменная ONLYBINARY, то библиотека в пакете отсутствует
ONLYBINARY=""

# Если необходимо добавить CERNLIB, то лучше воспользоваться этим
# ключиком. Церновские библиотеки подшиваются в самом конце. По
# умалчанию подшиваются следущие библиотеки jetset74 mathlib graflib
# geant321 grafX11 packlib
CERNLIBRARY = ""

# Список библиотек, если вам не нравится стандартный набор, который
# расширяем по требованию. В любом случае необходимо раскомментировать
# CERNLIBRARY
#CERNLIBS =
CERNLIBS = jetset74 mathlib graflib geant lapack3 blas packlib

# где вываливать исполняемые файлы
BINDIR := ./

# Дополнительные либы (вставляются после либ)
LIB_LOCAL= `root-config --libs` -lMinuit -lpq -lcrypt -lbz2 -ldl -lg2c

# Определим, какие программы мы будем собирать
BINARIES = analysis_j-psi_pi_pi_pi0

# укажем, из каких модулей этого пакета они состоят
# (эти модули не будут включены в библиотеку)
# и какие библиотеки надо подключить при сборке
analysis_j-psi_pi_pi_pi0_MODULES := analysis_j-psi_pi_pi_pi0
analysis_j-psi_pi_pi_pi0_LIBS := KaFramework KrAtc KDisplay VDDCRec KrVDDCMu KrMu \
KrdEdxPId KrDCCalibdEdx DchdEdxDataRoot VDDCRec KrToF KsToF KEmcRec LKrTools \
VDDCRec KsTrg KdConvert KrObjects KdDCSim FitTools DchGeom ReadNat KDB AppFramework KrKRec KrDONLP2

# следующую строку лучше не трогать
include $(WORKDIR)/KcReleaseTools/rules.mk


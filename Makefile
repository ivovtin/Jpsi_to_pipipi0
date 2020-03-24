#--------------------------------------------------------------------------
# ��������� 3 ������ ����� �� �������
ifeq (,$(WORKDIR))
WORKDIR=..
endif

# �ӣ ��� �������� �������� ����� :) - ���������� �������� ��� ����
VERBOSEMAKE=""
# ���� ������ -g ��� ������ �� ����� - ���� ������ � ����
#NDEBUG=""

# ���������� ��� *.� ������
#CC_LOCAL := gcc
# ���������� ��� *.cc � *.cpp ������
#CXX_LOCAL := g++
# ���������� ��� *.f, *.for � *.F ������
F77_LOCAL := g77
# ������ ��� ����� ������ �� ���������
#LD_LOCAL := g++
# �����������������, ���� ���������� �������� ������ ��� �������
# ��������� � �����������. ������ �� ���������.
#LD_MULTI := ""

# �������������� ����� ��� ����������� C/C++
COPTS = -I. `root-config --cflags` #-ansi -pedantic
NOOPT = ""
# �������������� ����� ��� ����������� Fortran
#FOPTS = -I.
FOPTS  = -g -fvxt -Wall -fno-automatic -finit-local-zero \
-fno-second-underscore -ffixed-line-length-120 -Wno-globals \
-DCERNLIB_LINUX -DCERNLIB_UNIX -DCERNLIB_LNX -DCERNLIB_QMGLIBC -DCERNLIB_BLDLIB \
-DCOMPDATE="'$(COMPDATE)'" -I$(CERN)/pro/include -I$(CERN)/pro/include/geant321

# �������������� ����� ��� �������
LDOPTS = -Xlinker -rpath -Xlinker `root-config --libdir`

#���� ���������� ���������� ONLYBINARY, �� ���������� � ������ �����������
ONLYBINARY=""

# ���� ���������� �������� CERNLIB, �� ����� ��������������� ����
# ��������. ���������� ���������� ����������� � ����� �����. ��
# ��������� ����������� �������� ���������� jetset74 mathlib graflib
# geant321 grafX11 packlib
CERNLIBRARY = ""

# ������ ���������, ���� ��� �� �������� ����������� �����, �������
# ��������� �� ����������. � ����� ������ ���������� �����������������
# CERNLIBRARY
#CERNLIBS =
CERNLIBS = jetset74 mathlib graflib geant lapack3 blas packlib

# ��� ���������� ����������� �����
BINDIR := ./

# �������������� ���� (����������� ����� ���)
LIB_LOCAL= `root-config --libs` -lMinuit -lpq -lcrypt -lbz2 -ldl -lg2c

# ���������, ����� ��������� �� ����� ��������
BINARIES = analysis_j-psi_pi_pi_pi0

# ������, �� ����� ������� ����� ������ ��� �������
# (��� ������ �� ����� �������� � ����������)
# � ����� ���������� ���� ���������� ��� ������
analysis_j-psi_pi_pi_pi0_MODULES := analysis_j-psi_pi_pi_pi0
analysis_j-psi_pi_pi_pi0_LIBS := KaFramework KrAtc KDisplay VDDCRec KrVDDCMu KrMu \
KrdEdxPId KrDCCalibdEdx DchdEdxDataRoot VDDCRec KrToF KsToF KEmcRec LKrTools \
VDDCRec KsTrg KdConvert KrObjects KdDCSim FitTools DchGeom ReadNat KDB AppFramework KrKRec KrDONLP2

# ��������� ������ ����� �� �������
include $(WORKDIR)/KcReleaseTools/rules.mk


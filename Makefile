RM        = rm -f 
SHELL     =/bin/sh

OBJF    = qfs_func.o qfs_sigs.o inelas_xs.o F1F209_FAST.o sigrad.o polsig.o haprad_utils.o grsv2000.o DSSV-2008.o fDSS.o g2_f.o getg.o getgn.o g1g2types.o xs_pol.o g2_nlo_ww.o xyrad2d.o

cc_OBJS = g1g2_cf.o g1g2_nlo.o g1g2_nlo_cf.o g1g2solve_cf.o g1g2_dssv_ww.o g1g2_dssv_cf.o g1g2_grsv_ww.o g1g2_grsv_cf.o

my_objs = $(OBJF) $(cc_OBJS)

MYOS := $(subst -,,$(shell uname))

#CERNLIBS = -Wl,-static -lgeant321 -lpawlib-lesstif -lpawlib -lmathlib -lgraflib -lpacklib-lesstif -lgrafX11 -lpacklib -lkernlib -Wl,-dy -lXbae -lXm -lXaw -llapack -lm -lXt -lX11 -lnsl -lcrypt -ldl

CERNLIBS = -Wl,-static -lpdflib804 -lgeant321 -lpawlib-lesstif -lpawlib -lmathlib -lgraflib -lpacklib-lesstif -lgrafX11 -lpacklib -lkernlib -Wl,-dy -lXaw -llapack -lm -lXt -lX11 -lnsl -lcrypt -ldl

#For use with gfortran compiler
# -fno-automatic - all program storage treated as static
ifeq ($(MYOS),Linux)
  LIBROOT = CTP/O.Linux/Linux/lib
#  CERN_ROOT = /site/cernlib/x86_64_rhel6/2005
# JLab
  #CERN_ROOT = /apps/cernlib/i386_fc8/2005
# 32 bit, standard Fedora distributuion
#  CERN_ROOT = /usr/lib/cernlib/2006
# 64 bit, standard Fedora distributuion
#  CERN_ROOT =  /usr/lib64/cernlib/2006#
  #CERN_ROOT = /w/halla-1/g2p/software/cernlib/2006#
  
 # CERN_ROOT=/usr/lib64/cernlib/2006
#  FFLAGSA=-O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
  FFLAGSA=-O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fno-align-commons -lstdc++
#  INCLUDES=-I. -Ibigbite
  INCLUDES=-I.
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
#  OTHERLIBS =-L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib64
  OTHERLIBS =-L/usr/lib64
  FC  := gfortran
  FF := gfortran


ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
LHAPDFLIBS    = $(shell lhapdf-config --cflags --ldflags)
CXX           = g++
CXXFLAGS      =  -lgfortran -lHistPainter -Wall -frtti -fexceptions -fPIC \
                   -DLINUXVERS -lgsl -lgslcblas -lm -I$(ROOTSYS)/include -I/var/phy/project/mepg/tl190/LHAPDF/include -I/var/phy/project/mepg/tl190/LHAPDF/boost_1_59_0 -L/var/phy/project/mepg/tl190/LHAPDF/lib -lLHAPDF  -O

LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

ALL_LIBS =  $(GLIBS) $(LIBS) $(LHAPDFLIBS) $(OTHERLIBS)
endif

%.o: %.f
	$(FF) $(FFLAGS) -c $< -o $@


%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

none: test_polrad_radiate_xyrad2d

all: test_polrad_radiate_xyrad2d

test_polrad_radiate_xyrad2d: test_polrad_radiate_xyrad2d.o $(my_objs) Makefile
	$(FF) $(OSF_SHARED) -o $@ $(FFLAGS) $(my_objs) test_polrad_radiate_xyrad2d.o $(ALL_LIBS)

clean:
	$(RM) *.o  test_polrad_radiate_xyrad2d


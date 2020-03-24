ifeq (,$(CC_LOCAL))
CC           := gcc
else
CC           := $(CC_LOCAL)
endif
ifeq (,$(CXX_LOCAL))
CXX           := g++-4.4
else
CXX           := $(CXX_LOCAL)
endif

CPPFLAGS+=$(shell root-config --cflags) -I.
#LDFLAGS+=$(shell root-config --libs) -Wl,-rpath,$(shell root-config --libdir)
LDFLAGS+=$(shell root-config --libs --glibs) -lGui -lMinuit -lNew -Wl,-rpath -Wl,$(shell root-config --libdir)

.PHONY: all clean

all: analysis_jpsito3pi

%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $^

analysis_jpsito3pi: analysis_jpsito3pi.o
	$(CXX) -o $@ $(LDFLAGS) $^
clean:
	rm -f analysis_jpsito3pi.o analysis_jpsito3pi


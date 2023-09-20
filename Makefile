ROOTCFLAGS     = $(shell root-config --cflags)
ROOTLIBS       = $(shell root-config --libs)
ROOTGLIBS      = $(shell root-config --glibs) 

INCLUDES       = -I./include
CXX            = g++
CXXFLAGS       = -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  $(INCLUDES) 
CXXFLAGS      += $(ROOTCFLAGS)

LD             = g++
LDFLAGS        = 

SOFLAGS        = -O -shared  -fPIC 
LIBS           = $(ROOTLIBS) 

GLIBS         = $(ROOTGLIBS) -lMinuit -lMinuit2 -lTreePlayer -lGenVector

SRCS = src/lepton_candidate.cc src/photon_candidate.cc src/proton_candidate.cc src/analysis.C
OBJS =  $(patsubst %.C,%.o,$(SRCS:.cc=.o))

LIB=lib/main.so


.SUFFIXES: .cc,.C,.hh,.h

# Rules ====================================
all: $(LIB)  bin/RunAll

lib : $(LIB)
$(LIB): $(OBJS)
	@echo "Creating library $(LIB)"
	mkdir -p lib
	$(LD) $(LDFLAGS) $(GLIBS) $(SOFLAGS) $(OBJS) -o $(LIB)
	@echo "$(LIB) successfully compiled!"

bin/RunAll : src/main.C $(LIB)	
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -ldl $(LDFLAGS) -o $@ $^ $(GLIBS)

clean:
	$(RM) $(OBJS)	
	$(RM) $(LIB)
	$(RM) bin/RunAll

purge:
	$(RM) $(OBJS)

deps: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)

CPPFLAGS += -g

TARGET = CreateAcceptanceSMS
TARGET2 = DM_Signal_SYS_Final
TARGET3 = EvalSys
TARGET4 = PlotDMsignal

SRC = src/CreateAcceptanceSMS.cc src/hlt.cc src/HelperFunctions.cc
SRC2 = src/DM_Signal_SYS_Final.cc src/hlt.cc
SRC3 = src/EvaluateSystematics.cc src/DM_1DRatio.cc src/DM_2DRatio.cc src/DM_Base.cc
SRC4 = src/PLOT_DM_SIGNAL.cc hlt.cc src/DM_1DRatio.cc src/DM_Base.cc

OBJ = $(SRC:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)

all : $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET3) : $(OBJ3)
	$(LD) $(CPPFLAGS) -o $(TARGET3) $(OBJ3) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET4) : $(OBJ4)
	$(LD) $(CPPFLAGS) -o $(TARGET4) $(OBJ4) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) *~


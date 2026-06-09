CXX = g++

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXXFLAGS = -Wall -Wextra -std=c++17 -Iinclude $(ROOTCFLAGS)


PLOTUTILS_INC = /exp/minerva/app/users/mmehmood/MAT_AL9/opt/include
PLOTUTILS_LIB = /exp/minerva/app/users/mmehmood/MAT_AL9/opt/lib
CXXFLAGS += -I$(PLOTUTILS_INC)
LDFLAGS  += -L$(PLOTUTILS_LIB)
LDLIBS   += -lMAT -lMAT-MINERvA

SRC = \
	src/NuECrossSection.cxx \
        src/MatrixBuilder.cxx \
        src/ElectronEnergyDistribution.cxx \
	app/make_xsec_matrices.cxx

TARGET = make_xsec_matrices

all:
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS) $(LDLIBS) $(ROOTLIBS)


clean:
	rm -f $(TARGET)

.SUFFIXES: .cpp .o

CXX = icpc

LIBDIR = /opt/vfem
ARPACK_DIR = $(LIBDIR)/ARPACK
SUPERLU_DIR = $(LIBDIR)/SuperLU
ARPACKPP_DIR = $(LIBDIR)/arpack++

INCS = -I$(LIBDIR)/yaml-cpp/include -I$(ARPACKPP_DIR)/examples/areig -I$(ARPACKPP_DIR)/include
LIBS = -L$(LIBDIR)/yaml-cpp/build -lyaml-cpp -lm $(SUPERLU_DIR)/superlu_linux.a $(SUPERLU_DIR)/blas_linux.a $(ARPACK_DIR)/libarpack_linux.a -lifport -lifcore

#CXXFLAGS = -Wno-write-strings -O2 $(INCS)
#CXXFLAGS = -O2 $(INCS)
CXXFLAGS = -g $(INCS)

OBJS = vfem.o calculation.o postprocess.o matrix.o

TARGET = ../run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJS) $(TARGET)

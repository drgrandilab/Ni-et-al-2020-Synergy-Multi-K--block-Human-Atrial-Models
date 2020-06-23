# Makefile
CXX= icpc
CFLAGS =  -Ilib -O3 -qopenmp -std=c++11 #-lz
#CFLAGS = -Wall -Ilib -fopenmp
#CFLAGS = -g -Ilib -ipo -O3 -openmp
CXXFLAGS = -Ilib -O3 -qopenmp -std=c++11 #-lz
LDFLAGS=-lz
LDLIBS=-lz


ONE_D_NCZ_VEC_ERP:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
ONE_D_NCZ_VEC:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
ONE_D_NCZ_VEC_Ini_Cond:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))


NCZ_Model_for_IC:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
# ONE_D_Updated_CNZ_C_VEC_IKurB:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
# ONE_D_Updated_ORd_INaB:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))
clean:
	rm ONE_D_TNNP_C_VEC_INICOND  lib/*.o


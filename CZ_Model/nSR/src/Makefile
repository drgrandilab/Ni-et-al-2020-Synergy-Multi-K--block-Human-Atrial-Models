# Makefile

CXX=g++

CFLAGS =  -Ilib -O3 -std=c++11   -g
CXXFLAGS = -Ilib -O3 -std=c++11  -g

# LDFLAGS=  -lz 

NCZ_Model:$(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))

clean:
	rm lib/*.o main


# 
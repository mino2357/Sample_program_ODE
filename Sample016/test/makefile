CXX = g++
#CXX = /usr/local/Cellar/gcc/7.2.0_1/bin/g++-7
TARGET = main_sim
CXXFLAGS = -std=c++1z -Wall -Wextra -O2 -march=native -mtune=native 
LDLFLAGS = -lstdc++
SRCS = main_sim.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

run: all
	./$(TARGET)

clean:
	rm $(TARGET)

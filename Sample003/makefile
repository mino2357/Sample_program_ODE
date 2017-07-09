CXX = /usr/local/bin/clang++ 
TARGET = main
CXXFLAGS = -std=c++1z -Wall -Wextra -O2 #-I/usr/local/include
LDFLAGS = #-L/usr/local/lib
LDLIBS = #-lm -lpng -lSiv3D -lboost_filesystem -lboost_system -lfreetype -lturbojpeg -lz -lGLEW -lharfbuzz -lc++
SRCS = main.cpp 
OBJS := $(SRCS:.cpp=.o)

all:$(TARGET)

# $(TARGET): $(OBJS)
#		$(CXX) $(LDLIBS) -o $@ $(OBJS)

run: all
	./$(TARGET)

clean:
	rm $(TARGET)

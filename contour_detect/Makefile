CC = gcc
CXX = g++
CFLAGS+=	-Wall -g -fopenmp 
#CFLAGS+=	-DHAVE_OPENMP 
CFLAGS+=	`pkg-config --cflags opencv`
LDFLAGS+=	`pkg-config --libs opencv`
CXXFLAGS:=	$(CFLAGS) 
FILES:=	*.cpp

all: contour canny contour_cpp
contour:  contour.c
		$(CC) $(CXXFLAGS)  $? -o  $@ $(LDFLAGS)

contour_cpp:  contour.cpp
		$(CXX) $(CXXFLAGS)  $? -o  $@ $(LDFLAGS)

canny:  canny.cpp
		$(CXX) $(CXXFLAGS)  $? -o  $@ $(LDFLAGS)

clean:
		rm -rf contour canny contour_cpp

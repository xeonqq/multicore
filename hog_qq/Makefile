CXX = g++
CFLAGS+=	-Wall -g -fopenmp 
#CFLAGS+=	-DHAVE_OPENMP 
CFLAGS+=	`pkg-config --cflags opencv`
LDFLAGS+=	`pkg-config --libs opencv`
CXXFLAGS:=	$(CFLAGS) 
FILES:=	*.cpp

all: hog_opencv hog_naive
hog_opencv:  hog_opencv.cpp
		$(CXX) $(CXXFLAGS)  hog_opencv.cpp -o  hog_opencv $(LDFLAGS)

hog_naive:  hog_naive.cpp
		$(CXX) $(CXXFLAGS)  hog_naive.cpp -o  hog_naive $(LDFLAGS)


clean:
		rm -rf hog_opencv hog_naive

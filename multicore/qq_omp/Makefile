CFLAGS+=	-fopenmp

CXXFLAGS:=	$(CFLAGS)

FILES:=	src/*.c 

all: hello




hello: $(FILES)
		$(CC) $(CXXFLAGS) $(FILES) -o hello  


clean:
		rm -rf hello

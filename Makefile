CFLAGS=-g -c -O0 -Wall -Wextra -Wwrite-strings
prog3: prog3.o image.o
	g++ -g -oprog3 prog3.o image.o
prog3.o: prog3.cpp image.h
	g++ ${CFLAGS} -oprog3.o prog3.cpp
image.o: image.cpp image.h
	g++ ${CFLAGS} -oimage.o image.cpp
clean:
	rm *.o
	rm prog3
	rm *.bmp
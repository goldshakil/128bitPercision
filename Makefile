# Makefile for Long Double

TARGET=ldfp
SOURCE=ldfp.c main.c
CFLAGS=-O3 -Wall
STUDENTID=201812345

$(TARGET):$(SOURCE)
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -rf *.o $(TARGET) 


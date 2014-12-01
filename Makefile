CC=gcc
CFLAGS=-Wall -Werror -Wextra -ansi -pedantic -std=c99
LDFLAGS=
LIBS=-lm
SOURCES=mesh.c crs_matrix.c fem.c norm.c exercise5.c implement_me.c list.c
OBJECTS=$(SOURCES:.c=.o)
TARGET=exercise5

.PHONY=all debug clean

all: $(SOURCES) $(TARGET)

debug: CC += -DDEBUG -DPRINT_DEBUG -g
debug: $(SOURCES) $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

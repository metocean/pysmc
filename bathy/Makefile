CC = gcc
LFLAGS = -lnetcdf `pkg-config --cflags --libs glib-2.0 gmodule-2.0` -lm
CFLAGS = -g $(LFLAGS)
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES = etopo_interp_centred.o
HEADERS =  
SRCS_C = etopo_interp.c
SRCS = $(SRCS_C) $(HEADERS)

all: etopo_interp_centred.exe

etopo_interp_centred.exe: $(OBJFILES)
	$(CC) $(OBJFILES) $(LFLAGS) -o etopo_interp_centred.exe

etopo_interp_centred.o: 
	$(COMPILE) etopo_interp_centred.c


clean:
	rm -f *.o

TAGS: $(SRCS)
	etags $(SRCS)
tags: $(SRCS)
	etags $(SRCS)

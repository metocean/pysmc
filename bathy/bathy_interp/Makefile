CC = gcc
LFLAGS = -lnetcdf `pkg-config --cflags --libs glib-2.0 gmodule-2.0` -lm `pkg-config --cflags --libs python3`
CFLAGS = -g $(LFLAGS)
COMPILE = $(CC) $(CFLAGS) -c
OBJFILES = etopo_interp_centred.o
HEADERS =  
SRCS_C = etopo_interp.c
SRCS = $(SRCS_C) $(HEADERS)


# defaults
resolution=1
bathy_src=etopo1

all: etopo_interp_centred.exe

etopo_interp_centred.exe: $(OBJFILES)
	$(CC) $(OBJFILES) $(LFLAGS) -o etopo_interp_centred.exe

etopo_interp_centred.o: etopo_interp_centred.c
	$(COMPILE) etopo_interp_centred.c

clean:
	rm -f *.o

run: etopo_interp_centred.exe
	./etopo_interp_centred.exe 0 360 2 -70 70 2 in.nc out.nc

glob: etopo_interp_centred.exe
	./etopo_interp_centred.exe 0 360 $(resolution) -80 80 $(resolution) /source/gridgen/noaa/reference_data/$(bathy_src).nc glob-$(bathy_src)-$(resolution).nc  

nz: etopo_interp_centred.exe
	./etopo_interp_centred.exe 165 180 $(resolution) -50 -30 $(resolution) /source/gridgen/noaa/reference_data/$(bathy_src).nc nz-$(bathy_src)-$(resolution).nc  
	ncview nz-$(bathy_src)-$(resolution).nc


TAGS: $(SRCS)
	etags $(SRCS)
tags: $(SRCS)
	etags $(SRCS)

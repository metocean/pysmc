from bathy_interp import bathy_interp

filename_in = "/net/homes/home/tdurrant/Public/etopo1.nc"
filename_out = "/tmp/NZ.nc"

res = bathy_interp(filename_in,
	           filename_out,
	           155, 180, 1.0,
	           30, 70, 1.0,
                   -32767)

assert(res == filename_out)

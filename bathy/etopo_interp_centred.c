
/*Interpolate etopo bathymetry to a coarser grid. */
/*Uses the same mehod as gridgen matlab code*/

#include <glib.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>

#define max(a, b) (((a) > (b)) ? (a) : (b)) 
#define min(a, b) (((a) < (b)) ? (a) : (b)) 


int main(int argc, char *argv[])
{
  int i, j;
  int ncid;
  int status;
  double * lons_in, * lats_in;
  int lonid, latid, zid, ndim, natt, countid;
  size_t nlons, nlats;
  nc_type type;
  int dimids[NC_MAX_VAR_DIMS];
  double fill_value_in;
  double cut_off = 0.1;

  // Here output grid params
  /*double lon0 = 0, lon1 = 360, dlon = 0.0625, lon;*/
  /*double lat0 = -75, lat1 = 75, dlat = 0.0625, lat;*/

  if( argc < 9 ) {
     printf("Too few arguments supplied.\n", argc);
     printf("eight arguments expected.\n \
     lon0, lon1, dlon, lat0, lat1, dlat, input_netcdf, output_netcdf \n");
     return(0);
  }
  if( argc > 9 ) {
     printf("Too many arguments supplied.\n", argc);
     printf("eight arguments expected.\n \
     lon0, lon1, dlon, lat0, lat1, dlat, input_netcdf, output_netcdf \n");
     return(0);
  }

  double lon0 = atof(argv[1]), lon1 = atof(argv[2]), dlon = atof(argv[3]), lon;
  double lat0 = atof(argv[4]), lat1 = atof(argv[5]), dlat = atof(argv[6]), lat;
  char * fnmin = argv[7];
  char * fnmout = argv[8];
  double ddlon = dlon/2.;
  double ddlat = dlat/2.;
  short fill_value = -32767;
  printf(" -- Grid Setup --\n");
  printf("    lon0, lon1, dlon:  %4.2f, %4.2f, %4.2f \n", lon0, lon1, dlon);
  printf("    lat0, lat1, dlat:  %4.2f, %4.2f, %4.2f \n", lat0, lat1, dlat);
  printf("    input :  %s\n", fnmin);
  printf("    output:  %s\n", fnmout);
  printf(" ------------------\n");

  // Here input file
  printf("-- Opening parent bathy %s \n", fnmin);
  status = nc_open(fnmin,
           NC_NOWRITE, &ncid);
  g_assert(status == NC_NOERR);
  
  int ncout;
  status = nc_create(fnmout, NC_NETCDF4, &ncout);
  g_assert(status == NC_NOERR);
  
  status = nc_inq_varid (ncid, "lon", &lonid);
  g_assert(status == NC_NOERR);
  status = nc_inq_varid (ncid, "lat", &latid);
  g_assert(status == NC_NOERR);
  status = nc_inq_varid (ncid, "z", &zid);
  g_assert(status == NC_NOERR);

  int londimid, latdimid;
  status = nc_inq_dimid (ncid, "lon", &londimid);
  g_assert(status == NC_NOERR);
  status = nc_inq_dimid (ncid, "lat", &latdimid);
  g_assert(status == NC_NOERR);

  status = nc_inq_dimlen(ncid, londimid, &nlons);
  g_assert(status == NC_NOERR);
  status = nc_inq_dimlen(ncid, latdimid, &nlats);
  g_assert(status == NC_NOERR);

  lons_in = malloc(nlons*sizeof(double));
  status = nc_get_var_double(ncid, lonid, lons_in);
  g_assert(status == NC_NOERR);
  lats_in = malloc(nlats*sizeof(double));
  status = nc_get_var_double(ncid, latid, lats_in);
  g_assert(status == NC_NOERR);
  
  status = nc_get_att_double(ncid, zid, "_FillValue", &fill_value_in);
  g_assert(status == NC_NOERR);


  // Check max grid extent
  lon1 = min(360-dlon, lon1);

  int nlat_out=0, nlon_out=0;
  for ( lat = lat0; lat+dlat <= lat1+dlat; lat+=dlat ){
    nlat_out++;
       // fprintf(stderr, "%i %f\n", nlat_out, lat);<]
  }
  for ( lon = lon0; lon+dlon <= lon1+dlon; lon+=dlon )
    nlon_out++;
  double lats_out[nlat_out];
  double lons_out[nlon_out];
  i = 0;
  for ( lat = lat0; lat+dlat <= lat1+dlat; lat+=dlat ){
    lats_out[i++] = lat;
    //    fprintf(stderr, "%i %f\n", i, lat);
  }
  i = 0;
  for ( lon = lon0; lon+dlon <= lon1+dlon; lon+=dlon )
    lons_out[i++] = lon;
  
  //fprintf(stderr, "%i %i\n",nlon_out, nlat_out);
  int latdimid_out, londimid_out, latvarid_out, lonvarid_out, zvarid_out, countvarid_out;
  status = nc_def_dim (ncout, "lat", nlat_out, &latdimid_out);
  g_assert(status == NC_NOERR);
  status = nc_def_dim (ncout, "lon", nlon_out, &londimid_out);
  g_assert(status == NC_NOERR);
  dimids[0] = latdimid_out;
  status = nc_def_var (ncout, "lat", NC_FLOAT, 1, dimids, &latvarid_out);
  g_assert(status == NC_NOERR);
  dimids[0] = londimid_out;
  status = nc_def_var (ncout, "lon", NC_FLOAT, 1, dimids, &lonvarid_out);
  g_assert(status == NC_NOERR);
  dimids[0] = latdimid_out;
  dimids[1] = londimid_out;
  status = nc_def_var (ncout, "z", NC_SHORT, 2, dimids, &zvarid_out);
  g_assert(status == NC_NOERR);
  status = nc_def_var (ncout, "count", NC_INT, 2, dimids, &countvarid_out);
  g_assert(status == NC_NOERR);
  status = nc_put_att_short (ncout, zvarid_out, "_FillValue", NC_SHORT,
                 1, &fill_value);
  g_assert(status == NC_NOERR);
  status = nc_put_att_int (ncout, countvarid_out, "_FillValue", NC_INT,
                 1, &fill_value);
  g_assert(status == NC_NOERR);
  status = nc_enddef(ncout);
  g_assert(status == NC_NOERR);
  
  status = nc_put_var_double (ncout, latvarid_out, lats_out);
  g_assert(status == NC_NOERR);
  status = nc_put_var_double (ncout, lonvarid_out, lons_out);
  g_assert(status == NC_NOERR);
  
  //return 1;
  //dimids[0] = 
  //status = nc_def_var (*ncid, "lon", NC_FLOAT, 1, dimids, &var->id);
  
  size_t startlon[nlon_out], startlat[nlat_out];
  size_t countlon[nlon_out], countlat[nlat_out];
  int ilat, ilon;
  for ( i = 0; i < nlon_out; i++ ) {
    countlon[i] = 0;
    startlon[i] = -999;
  }
  for ( i = 0; i < nlat_out; i++ ) {
    countlat[i] = 0;
    startlat[i] = -999;
  }
  
  for ( ilon = 0; ilon < nlons; ilon++ ) {
    if (  *(lons_in+ilon) < lon0-ddlon || *(lons_in+ilon) > lon1+dlon+ddlon )
      continue;
    i = ((*(lons_in+ilon)-lon0)+ddlon)/dlon;
    if ( lon0-ddlon + i*dlon == *(lons_in+ilon) && ilon != 0 ) {
      if ( startlon[i-1] == -999 )
        startlon[i-1] = ilon;
        countlon[i-1]++;
      //fprintf(stderr, "%f %i ++\n",*(lons_in+ilon), i-1);
    }

    if ( i < nlon_out ) {
      if ( startlon[i] == -999 )
        startlon[i] = ilon;
        countlon[i]++;
      //fprintf(stdout, "%f %i\n",*(lons_in+ilon), i);
    }
  }


  for ( ilat = 0; ilat < nlats; ilat++ ) {
    if (  *(lats_in+ilat) < lat0-ddlat || *(lats_in+ilat) > lat1+dlat+ddlat )
      continue;
    i = (*(lats_in+ilat)-lat0+ddlat)/dlat;
    if ( lat0-ddlat + i*dlat == *(lats_in+ilat) && ilat != 0 ) {
       // fprintf(stderr, "%f %i ++\n",*(lats_in+ilat), i-1);
      if ( startlat[i-1] == -999 )
        startlat[i-1] = ilat;
        countlat[i-1]++;
        // fprintf(stderr, "%f %i ++\n",*(lats_in+ilat), i-1);<]
    }

    if ( i < nlat_out ) {
      if ( startlat[i] == -999 )
        startlat[i] = ilat;
        countlat[i]++;
        // fprintf(stderr, "%f %i\n",*(lats_in+ilat), i);<]
    }
  }


  //  for ( ilon = 0; ilon < nlon_out; ilon++ )
  //  fprintf(stderr, "%i\n", countlon[ilon]);
  //return 1;

  printf("-- looping through cells\n");

  
  size_t start[2];
  size_t count[2];
  double sum;
  int num, total;
  size_t index[2];
  for ( ilat = 0; ilat < nlat_out; ilat++ ) {
    start[0] = startlat[ilat];
    count[0] = countlat[ilat];
    
    for ( ilon = 0; ilon < nlon_out; ilon++ ) {
      start[1] = startlon[ilon];
      count[1] = countlon[ilon];

      index[0] = ilat;
      index[1] = ilon;

      if ( start[0] == -999 || start[1] == -999 || count[0] == 0 || count[1] == 0 ) {
    status = nc_put_var1_short(ncout, zvarid_out, index,
                    &fill_value);
    g_assert(status == NC_NOERR);
    continue;
      }
    
      
      double data[count[0]][count[1]];
      status = nc_get_vara_double(ncid, zid, start, count, *data);
      g_assert(status == NC_NOERR);
      
      sum = 0.;
      num = 0;
      total = 0;

      for ( i = 0; i < count[0]; i++ ) {
        for ( j = 0; j < count[1]; j++ ) {
          total++;
        // Here criteria to change //
          if ( data[i][j] < 0 && data[i][j] != fill_value_in ) {
          sum += data[i][j];
          num++;
      }
    }
      }
      if ( num == 0 )
    sum = fill_value;
      else if ( num/(float) total < (1-cut_off) )
    sum = fill_value;
      else
    sum /= num;

      if ( ilon == 0 && (((int)(100*ilat/nlat_out) % 10) == 0) ) {
            printf("---- Processed %i of %i rows... \n", ilat, nlat_out);
      }

      status = nc_put_var1_double(ncout, zvarid_out, index,
                    &sum);
      g_assert(status == NC_NOERR);
      status = nc_put_var1_int(ncout, countvarid_out, index,
                    &num);
      g_assert(status == NC_NOERR);
    }
  }
  
  nc_close(ncout);
  nc_close(ncid);

  printf("-- Closed output %s \n", fnmout);
  printf("-- Finished successfully", fnmout);

  free(lons_in);
  free(lats_in);
  return 1;
}

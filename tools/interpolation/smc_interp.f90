!----------------------------------------------------------------------
! Program intepolate SMC gridded model data (in netCDF format) to a new
! grid using a pre-calculated remap file generated using SCRIP.
!
! Chris Bunney
! 24-May-2017
!----------------------------------------------------------------------
      PROGRAM SMC_INTERP

      USE NETCDF
      USE remap

      IMPLICIT NONE

      INTEGER, PARAMETER        :: MAX_FLDS = 100

      ! Filenames:
      CHARACTER(LEN=256)        :: remap_filename, smc_filename, out_filename

      ! Unit numbers for remap, SMC file and output file/
      INTEGER                   :: ncsmc, ncout

      ! Storage for SMC field data:
      REAL, ALLOCATABLE         :: smcfld(:)
      DOUBLE PRECISION, ALLOCATABLE :: times(:)

      ! General
      INTEGER                   :: err, dimid, varid, smc_varid
      INTEGER                   :: outdims(3), dimvars(3)
      REAL                      :: lsm_thresh = 0.5, MDI, fac,          &
                                   fracsea, outmdi=-32767
      LOGICAL                   :: lsm_invariate = .TRUE.

      ! SMC input grid:
      INTEGER                   :: nsea, ntimes, itime, di, ix, iy,     &
                                   n, i, ifld, nflds

      CHARACTER(LEN=20)         :: flds(MAX_FLDS)

      ! Output grid def:
      REAL                      :: dst_x1, dst_y1, dst_dx, dst_dy,     &
                                   dst_plat, dst_plon
      REAL, ALLOCATABLE         :: dst_x(:), dst_y(:)

      NAMELIST /inputs/                                                &
            remap_filename, smc_filename, out_filename,                &
            flds, fracsea

      !-----------------
      ! TESTING! 
      !-----------------
      !remap_filename='/home/h06/frey/scriptest/rmp_gblsmc_to_ukrot_conserv.nc'
!      remap_filename='/home/h06/frey/scriptest/remap/rmp_gblsmc_to_eurortd_cnsrv.nc'
!      smc_filename = '/scratch/frwave/wave_rolling_archive/gbl/gbl_2017052500.nc1'
!      out_filename='out.nc'

      dst_nx=549
      dst_ny=492
      dst_x1=-11.650
      dst_y1=-19.940
      dst_dx=0.08
      dst_dy=0.08
      dst_plon=177.50
      dst_plat=37.55

      fracsea=0.5
      flds(:) = ''
      nflds = 0

      ! Read namelist
      OPEN(10, file='interp.nl')
      READ(10, nml=inputs)
      CLOSE(10)

      ! Count number of fields:
      DO i=1,MAX_FLDS
         IF(flds(i) .EQ. '') EXIT
         nflds = nflds + 1
      ENDDO
      WRITE(*,'("Read in ",I3," field names to process")') nflds

      !-----------------
      ! END TESTING! 
      !-----------------


      !--------------------------------------------------/-
      ! Read the remap file:    
      CALL read_remap_file(remap_filename)


      !--------------------------------------------------/-
      ! Open input netCDF SMC file
      err = NF90_OPEN(smc_filename, NF90_NOWRITE, ncsmc)
      IF(err /= NF90_NOERR) THEN
         PRINT*,'Error opening SMC file: '//smc_filename
      ENDIF
     
      ! Get dimension sizes:
      err = NF90_INQ_DIMID(ncsmc, 'seapoint', dimid)
      if(err /= NF90_NOERR) then
         PRINT*,'failed to get seapoint id'
         stop
      endif
      err = NF90_INQUIRE_DIMENSION(ncsmc, dimid, len=nsea)
      if(err /= NF90_NOERR) then
         PRINT*,'failed to get seapoint size'
         stop
      endif

      err = NF90_INQ_DIMID(ncsmc, 'time', dimid)
      err = NF90_INQUIRE_DIMENSION(ncsmc, dimid, len=ntimes)

      ! Allocate input field storage:
      ALLOCATE(smcfld(nsea))


      !--------------------------------------------------/-
      !! Create output nc
      err = NF90_CREATE(out_filename, IOR(NF90_CLOBBER, NF90_HDF5), ncout)
      IF(err /= NF90_NOERR) THEN
         PRINT*,'Error opening output file: '//smc_filename
      ENDIF

      err = NF90_DEF_DIM(ncout, 'time', ntimes, outdims(3))
      err = NF90_DEF_DIM(ncout, 'latitude', dst_ny, outdims(2))
      err = NF90_DEF_DIM(ncout, 'longitude', dst_nx, outdims(1))

      ! Dimension vars (copy attributes from SMC file):
      err = NF90_DEF_VAR(ncout, 'time', NF90_DOUBLE, outdims(3), dimvars(3))
      CALL copy_attrs(ncsmc, ncout, 'time')

      err = NF90_DEF_VAR(ncout, 'latitude', NF90_FLOAT, outdims(2), dimvars(2))
      CALL copy_attrs(ncsmc, ncout, 'latitude')

      err = NF90_DEF_VAR(ncout, 'longitude', NF90_FLOAT, outdims(1), dimvars(1))
      CALL copy_attrs(ncsmc, ncout, 'longitude')

      err = NF90_ENDDEF(ncout)


      ! Put dimension variable data:
      ALLOCATE(dst_x(dst_nx), dst_y(dst_ny))
      DO i=1,dst_nx
         dst_x(i) = dst_x1 + (i-1) * dst_dx
      ENDDO
      DO i=1,dst_ny
         dst_y(i) = dst_y1 + (i-1) * dst_dy
      ENDDO
      err = NF90_PUT_VAR(ncout, dimvars(1), dst_x)
      err = NF90_PUT_VAR(ncout, dimvars(2), dst_y)

      ! Copy times:
      ALLOCATE(times(ntimes))
      err = NF90_INQ_VARID(ncsmc, 'time', varid)
      err = NF90_GET_VAR(ncsmc, varid, times)
      err = NF90_PUT_VAR(ncout, dimvars(3), times)
      DEALLOCATE(times)


      ! -----------------------------------------------------------
      ! If coast mask is invariate, then calculate now:
      ALLOCATE(lsm(dst_nx,dst_ny))
      lsm(:,:) = .TRUE. ! false is sea, true is land
      IF(lsm_invariate) THEN
         DO n = 1,nlinks
            di = dst_addr(n)
            ix = MOD(di - 1,dst_nx) + 1
            iy = INT((di - 1)/ dst_nx) + 1
            IF(dst_frac(di) .GE. lsm_thresh) lsm(ix,iy) = .FALSE.
         ENDDO
      ENDIF



      ! ----------------------------------------------------------
      ! Loop over fields:
      DO ifld = 1,nflds

         ! ----------------------------------------------------------
         ! Interoplate!
         ! Get handle to input SMC variable:
         err = NF90_INQ_VARID(ncsmc, TRIM(flds(ifld)), smc_varid)
         err = NF90_GET_ATT(ncsmc, smc_varid, '_FillValue', MDI)
         err = NF90_GET_ATT(ncsmc, smc_varid, 'scale_factor', fac)

         ! Create new output variable
         err = NF90_REDEF(ncout)
         err = NF90_DEF_VAR(ncout, TRIM(flds(ifld)), NF90_FLOAT, outdims, varid)
         CALL copy_attrs(ncsmc, ncout, TRIM(flds(ifld)))
         err = NF90_PUT_ATT(ncout, varid, '_FillValue', outmdi)
         err = NF90_PUT_ATT(ncout, varid, 'scale_factor', 1.0) ! TODO: Scaling turned off!!
         err = NF90_ENDDEF(ncout)

         ! Loop over times:
         DO itime = 1,ntimes
   
            WRITE(*,'("Interpolating field: ",A10,"; timestep:",I3 )')  &
                    flds(ifld), itime
   
            ! Get input field for this timestep:
            err = NF90_GET_VAR(ncsmc, smc_varid, smcfld, start=(/1,itime/), count=(/nsea,1/))
            if(err /= NF90_NOERR) then
               PRINT*,'failed to load ',flds(ifld)
               stop
            endif
   
            ! Scale field (at non-mdi values) if scale_factor set:
            IF(fac .NE. 1.0) THEN
               WHERE(smcfld .NE. mdi) smcfld = smcfld * fac
            ENDIF
   
            ! Interpolate field:
            CALL interpolate(smcfld, mdi)
            
            ! Write field to netcdf output file:
            err = NF90_PUT_VAR(ncout, varid, outfld, start=(/1,1,itime/), count=(/dst_nx,dst_ny,1/))
   
         ENDDO ! Time loop
      ENDDO ! field loop

      err = NF90_CLOSE(ncout)
      err = NF90_CLOSE(ncsmc)


      CONTAINS 

      SUBROUTINE copy_attrs(nc1, nc2, varname)

      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: nc1, nc2
      CHARACTER(LEN=*), INTENT(IN)  :: varname

      ! Locals
      INTEGER                       :: err, n
      INTEGER                       :: varid1, varid2
      CHARACTER(LEN=128)            :: attname


      err = NF90_INQ_VARID(nc1, TRIM(varname), varid1)
      err = NF90_INQ_VARID(nc2, TRIM(varname), varid2)

      err = NF90_INQUIRE_VARIABLE(nc1, varid1, nAtts=n)
      print*,'Variable has ',n,'attrs'

      DO i=1,n
         err = NF90_INQ_ATTNAME(nc1, varid1, i, attname)
         PRINT*,"Copying attribute",attname
         err = NF90_COPY_ATT(nc1, varid1, TRIM(attname), nc2, varid2)
      ENDDO

      END SUBROUTINE copy_attrs


      END PROGRAM SMC_INTERP

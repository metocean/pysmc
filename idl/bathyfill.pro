;;  Fill isolated water points as land to avoid single sea point 
;;  in wave model grid.  Check all 8 points around a single sea
;;  point and fill it as land if more than mnpt points are land.
;;                               JGLi01Sep2014
; First created:    1 Sep 2014   Jian-Guo Li
; Last modified:    1 Sep 2014   Jian-Guo Li
;+
; name:     bathyfill  
;
; purpose:  Reads bathy file and fill single sea points.
;
; usage:    bathyfill, bathy, ncol, nrow, MinPnt=mnpt, LndSgn=lnsg, Global=glob
;
; input:    bathy--- 2-D bathymetry array.
;           ncol --- number of column in bathy.
;           nrow --- number of rows in bathy.  
;           mnpt --- minimum land points around a sea point to invoke filling.
;                    default 6 out of maximum 8.
;           lnsg --- sign of land point bathy vaelue (default positive).
;           glob --- 1 for global bathy to use periodical NCol loop. 
; output:   bathy--- same as input but filled single sea points as land.
;-
 
PRO  BATHYFILL, Bathy, NCol, NRow, MinPnt=mnpt, LndSgn=lnsg, Global=glob

;; Check Input NCol, NRow, should be greater or equal 3
     IF( (NCol LT 3) OR (NRow LT 3) ) THEN Begin
         Print, "NCol and NRow have to be greater than 3.", NCol, NRow
         Return
     ENDIF

;; Default settings if not provided by Keywords
  IF( NOT  PARAM_PRESENT(mnpt) ) THEN mnpt=6 
  IF( NOT  PARAM_PRESENT(lnsg) ) THEN lnsg=1 
  IF( NOT  PARAM_PRESENT(glob) ) THEN glob=0     ;; Assume regional

;; Check minimum land points around sea point, warning if less than 4
     IF( (mnpt LT 4) ) THEN Begin
         Print, "Minimum land points (out of 8) around a sea point too low ", mnpt
         Print, "Recommend minimum 4 points out of 8.  Default 6."
         Return
     ENDIF

;; Temporary array filin for input bathy data
;; Ensure land bathy is always positive.
      filin = bathy*lnsg

;; Regional model loop over inner cells only
      NR = NRow - 1 
      NC = NCol - 1 
;; Global bathy loop over last colum then first one
;; by wrapping over NCol.
      IF( glob GT 0 ) then NC = NCol + 1 

;;    Display input settings
      Help, Bathy, Filin
      Print, "NCol, NRow, NC, NR, mnpt, lnsg, glob"
      Print,  NCol, NRow, NC, NR, mnpt, lnsg, glob
 
;; Fill single grid width river or lake and straighten coastal line
   NRound=0
   NSS=10L
   WHILE (NSS GT 1L ) DO BEGIN
   NS=0L
     FOR ii=1, NC-1 DO BEGIN
         i   = (ii  ) mod NCol 
         ip1 = (ii+1) mod NCol
         im1 = (ii-1) mod NCol

     IF(ii GE NCol-1 AND NRound EQ 0) THEN Print, ii, i, ip1, im1

     FOR j=1, NR-1 DO BEGIN
        IF( filin(i,j) LT 0.0 ) THEN BEGIN
        NIJ=0
        SIJ=0.0
           IF( filin(im1,j) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(im1,j)
           ENDIF
           IF( filin(im1,j-1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(im1,j-1)
           ENDIF
           IF( filin(im1,j+1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(im1,j+1)
           ENDIF
           IF( filin(ip1,j) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(ip1,j)
           ENDIF
           IF( filin(ip1,j-1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(ip1,j-1)
           ENDIF
           IF( filin(ip1,j+1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(ip1,j+1)
           ENDIF
           IF( filin(i,j+1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(i,j+1)
           ENDIF
           IF( filin(i,j-1) GE 0.0 ) THEN BEGIN
               NIJ=NIJ+1
               SIJ=SIJ+filin(i,j-1)
           ENDIF

           IF(NIJ GE mnpt) THEN BEGIN
             newdepth = SIJ/FLOAT(NIJ)
             NS=NS+1L
             IF(NS MOD 500 EQ 0 ) THEN PRINT,   " i,j,N,Sea-Lnd =", i, j, NIJ, filin(i,j), newdepth
             filin(i,j) = newdepth
           ENDIF
        ENDIF
     ENDFOR
     ENDFOR
       NSS=NS
       NRound = NRound + 1
     PRINT,    " ********* Round", NRound, " filled sea points =", NSS
  ENDWHILE

;; Reset Bathy with filled data
   Bathy = filin*lnsg

   PRINT,    " All filling done and bathy updated." 
 
 RETURN

 END



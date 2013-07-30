;PROCEDURE MAGSYM6.PRO
;A PROCEDURE TO DISPLAY RESULTS FROM MAGNETOCONVECTION & DYNAMO CALCULATIONS
;USING THE CODE MAG BY G.A.GLATZMEIER WITH MODIFICATIONS BY U.CHRISTENSEN
;SOME LONGITUDINAL SYMMETRY MAY BE ASSUMED IN THE CALCULATION.
;ADAPTED FOR PARALLEL OUTPUT OF GRAPHIC DATA BY C.KUTZNER
;VERSION6 DISPLAYS INNER CORE DATA AS WELL.
 
       VERSION="Magsym 6  (September 4, 2000)"

;DETERMINE OPERATING SYSTEM
       CASE !VERSION.OS_FAMILY OF
          "MacOS" :  BEGIN
             PLATFORM="MAC"
             SCWINDOW=0.95
             END
          "unix" : BEGIN
             PLATFORM="X"
             SCWINDOW=1.5
             END
          "Windows" : BEGIN
             PLATFORM="WIN" 
             SCWINDOW=1.25
             END
        ENDCASE
     
;ACCESS RGB-VALUES OF COLOR TABLE
        COMMON COLORS,R_ORIG,G_ORIG,B_ORIG,R_CURR,G_CURR,B_CURR

;SIZE OF THE STABLE LAYER
        VCOND=0.8

;DEFINE THE VALUE OF PI:
	PI= 3.1415926

LABELREAD: CLOSE, 1  ;close old file if any
;READ IN DATA FILENAME.
	INFILE=''
	PRINT, FORMAT='("Enter data filename:")'
	READ, INFILE

;OPEN DATA FILE AND READ HEADER

        GRAPHOUTVERSION=1  ; default version
        NN_IC=0
        ;PRINT, FORMAT='$, "Format? Old=0, New=1")
        ;READ, GRAPHOUTVERSION
        ;
        ;---BINARY DATA FROM PARALLEL OUTPUT ROUTINE----------------------------
        ;

        IF ( STRMID(INFILE,0,1) EQ 'G' ) THEN BEGIN

          GRAPHOUTVERSION=3   ; GRAPHOUTVERSION >=3, this eventually gets
          ;                     overwritten if the version is even newer
          OPENU, 1, INFILE, /F77_UNFORMATTED

          ; determine file format (GRAPHOUTVERSION) and read run identifyer RUNID
          INFOSTRING='12345678901234567890'  ; 20 byte
          READU,1,INFOSTRING
	  PRINT,'Infostring is:',INFOSTRING
          IF STRPOS(INFOSTRING, 'Graphout_Version_5' ) NE -1 THEN $
             GRAPHOUTVERSION=5 $
	  ELSE IF STRPOS(INFOSTRING, 'Graphout_Version_7' ) NE -1 THEN $
	     GRAPHOUTVERSION=7

	  print, 'Graphout version is:',graphoutversion

          ; runid is a 64 byte string, this has to be exact for binary data.
          RUNID='1234567890123456789012345678901234567890123456789012345678901234'
          READU, 1, RUNID
	  PRINT, 'Runid is:',RUNID

; now read the header data:
	  IF ( GRAPHOUTVERSION EQ 7 ) THEN BEGIN

            READU, 1, TIME,NNf,NIf,NJf,NN_ICf, $
                      NSYMf,NFSIZEf,RA,EK,PR,PM, $
                      RADRATIO,SIGMA_RATIO
	     NGRAD=1
             NGCOLAT=1
             NGLON=1
	     NGRAD_IC=1
	     NN=FIX(NNf)
             NI=FIX(NIf)
             NJ=FIX(NJf)
             NN_IC=FIX(NN_ICf)
             NSYM=FIX(NSYMf)
             NFSIZE=FIX(NFSIZEf)

	  ENDIF ELSE BEGIN

	    PRINT,'Not supported any more!'
	    STOP

	  ENDELSE

        END ELSE BEGIN

	  PRINT,' ASCII NOT SUPPORTED'

	ENDELSE

        PRINT        
        PRINT, RUNID

 	PRINT, FORMAT='("    Ra =", E11.3, "  n_r     =    ", I3 )', RA , nn
        PRINT, FORMAT='("    Ek =", E11.3, "  n_theta =    ", I3 )', EK , nj
        PRINT, FORMAT='("    Pr =", E11.3, "  sym     =    ", I3 )', Pr , NSYM
        PRINT, FORMAT='("    Pm =", E11.3, " radratio = ", F8.6)', Pm, RADRATIO
 	PRINT, FORMAT='(" Sratio=", E11.3, "  n_r_ic  =    ", I3 )', sigma_ratio , nn_ic
        PRINT, FORMAT='("  TIME =", E11.3)', time

        IMAG=1
        IF PM LE 0.0 AND GRAPHOUTVERSION LT 4 THEN BEGIN
          PRINT, FORMAT='($, "magnetic field data in input file? yes=1 no=0")'
          READ,IMAG
        END
        IF ( Pm EQ 0 ) THEN IMAG=0

        PRINT, FORMAT='($, "color=1,  black & white=0, read another file=2, quit=-1 ?")'
        READ, LCOLOR
        IF LCOLOR EQ  2 THEN GOTO, LABELREAD ELSE $
        IF LCOLOR EQ -1 THEN GOTO, LABEL99
        PRINT, FORMAT='($, "draw zero contour? yes=1 no=0")'
        READ,IZERO
        ISINGLE=0
        IF IZERO LT 0 THEN ISINGLE=1
        IZERO=ABS(IZERO)

;DEFINES ARRAY SIZES
;	NR=(NN+NGRAD-1)/NGRAD                   ;# RADIAL GRID POINTS
	if ( graphoutversion le 5 ) then begin
           NR=CEIL(NN/NGRAD)                       ;# radial grid points in outer core
           NR_IC=CEIL(NN_IC/NGRAD_IC)              ;# radial grid points in inner core
           IF GRAPHOUTVERSION GE 2 THEN $          ;# THETA  GRID POINTS
              NT=2*CEIL(NI/(2.*ngcolat)) $
           ELSE $
 	      NT=(NI+NGCOLAT-1)/NGCOLAT
           NPFULL=CEIL(NJ/NGLON)
	end else begin
	   nr=nn
	   nr_ic=nn_ic
	   nt=ni
	   npfull=nj
	end 

        NP=NPFULL/NSYM                          ;# PHI-POINTS IN DATA FILE
	PI2NP = 2*PI/NPFULL
        PI4NP = 4*PI/NPFULL

;DECLARE COORDINATE ARRAYS
	THETA = FLTARR(NT)
        STHET = FLTARR(NT)
        CtHET = FLTARR(NT)
        STHET2= FLTARR(NPFULL,NT)
        CTHET2= FLTARR(NPFULL,NT)
	LAT   = FLTARR(NT)
        XIC   = FLTARR(NT)
        YIC   = FLTARR(NT)
        PHI   = FLTARR(NPFULL)
	LON   = FLTARR(NPFULL)
        LONWRAP  = FLTARR(NPFULL+1)
        THZON = FLTARR(NT,NR+NR_IC)
        THPLT = FLTARR(NT,NR+NR_IC)
        TH3D = FLTARR(NPFULL,NT,NR+NR_IC)
        PHIEQ = FLTARR(NPFULL,NR+NR_IC)
        PHI3D = FLTARR(NPFULL,NT,NR+NR_IC)

;read in array THETA of colatitudes of data points (radians):
        IF (GRAPHOUTVERSION MOD 2) EQ 1 THEN $
          READU,1,THETA $
        ELSE $
          READF,1,THETA
;create 2D array THZON and 3D array TH3D
        FOR IT=0,NT-1 DO BEGIN
            THZON(IT,*)=THETA(IT)+PI/2
            TH3D(*,IT,*)=THETA(IT)
        ENDFOR
;create 2D array THPLT for plotting
        THPLT=THZON & THPLT(0,*)=PI/2 & THPLT(NT-1,*)=3*PI/2

;create arrays of sin(theta) and cos(theta)
        STHET=SIN(THETA)
        CTHET=COS(THETA)

;converts THETA to an array LAT of latitudes 90 to -90 degs North to South:
	LAT=90-180*THETA/PI
;create array PHI of longitudes of data points (in radians)
	PHI=PI2NP*FINDGEN(NPFULL)
;create phi array for three-d and for equatorial contours
        FOR IP=0,NPFULL-1 DO BEGIN
            PHIEQ(IP,*)=PHI(IP)
            STHET2(IP,*)=STHET
            CTHET2(IP,*)=CTHET
            PHI3D(IP,*,*)=PHI(IP)
        ENDFOR
        PHIEQ((NPFULL-1),0)=PHIEQ(0,0) ;DEVICE TO CLOSE OUTER CIRCLE

;converts PHI to array LON from -180 to 180 degrees
	LON=180*PHI/PI - 180.0
        LONWRAP(0:NPFULL-1)=LON & LONWRAP(NPFULL)=180

;DIMENSION STATEMENTS FOR 2D AND 3D ARRAYS
	R=FLTARR(NR+NR_IC)          ;array of radii to spherical 2D surfaces
        REQ=FLTARR(NPFULL,NR+NR_IC) ;2D array of radii for equatorial contouring
        TEQ=FLTARR(NPFULL,NR)       ;2D array of temps for equatorial contouring
        WZE=FLTARR(NPFULL,NR)       ;2D array of vort for equatorial contouring
        BZE=FLTARR(NPFULL,NR+NR_IC) ;2D array of field for equatorial contouring
        RZON=FLTARR(NT,NR+NR_IC)    ;2D array of radii for meridional contouring
        TZON=FLTARR(NT,NR)          ;2D array of temps for meridional contouring
     dTdtZON=FLTARR(NT,NR)          ;2D array of dTdt for meridional contouring
     dTdtMER=FLTARR(NPFULL,NR)      ;2D array of dTdt for meridional contouring
        TMER=FLTARR(NPFULL,NR)      ;2D array of dTdt for meridional contouring
       TMER2=FLTARR(NPFULL,NR)      ;2D array of dTdt for meridional contouring
    dTdtMER2=FLTARR(NT,NR)      ;2D array of dTdt for meridional contouring
        BZON=FLTARR(NT,NR+NR_IC)    ;2D array of Bphi for meridional contouring
        VZON=FLTARR(NT,NR)          ;2D array of Vphi for meridional contouring
    dVPdzZON=FLTARR(NT,NR)          ;2D array of dVPdz for meridional contouring
    dVPdzMER=FLTARR(NPFULL,NR)      ;2D array of dVPdz for meridional contouring
   dVPdzMER2=FLTARR(NT,NR)      ;2D array of dVPdz for meridional contouring
        JZON=FLTARR(NT,NR+NR_IC)    ;2D array of Jphi for meridional contouring
       BRZON=FLTARR(NT,NR+NR_IC)    ;2D array of Br for meridional contouring
       BTZON=FLTARR(NT,NR+NR_IC)    ;2D array of Btheta for meridional contouring
       VRZON=FLTARR(NT,NR)          ;2D array of Vr for meridional contouring
       VTZON=FLTARR(NT,NR)          ;2D array of Vt for meridional contouring
       OMZON=FLTARR(NT,NR)          ;2D array of omega-effect; meridional contouring
       HBZON=FLTARR(NT,NR)          ;2D array of alpha-effect for Bpol; meridional
        HZON=FLTARR(NT,NR)          ;2D array of zonally averaged helicity
        WZON=FLTARR(NT,NR)          ;2D array of Z-Vorticity for meridional contour.
        XZON=FLTARR(NT,NR+NR_IC)    ;2D array of Xdist for meridional contouring
      ROTZON=FLTARR(NT,NR)          ;2D array of omega for meridional contouring
          B2=FLTARR(NT,NR)          ;2D array of zonal average of rms B^2
           A=FLTARR(NT,NR+NR_IC)    ;2D array of poloidal magnetic potential
         PSI=FLTARR(NT,NR)          ;2D array of meridional streamfunction
         VXM=FLTARR(NT,NR)          ;2D array of meridional velocity
         VZM=FLTARR(NT,NR)          ;2D array of meridional velocity
       VXEQ=FLTARR(NPFULL,NR)       ;2D array of vx in equatorial plane
       VYEQ=FLTARR(NPFULL,NR)       ;2D array of vy in equatorial plane
     WORK2D=FLTARR(NT,NR)           ;2d work array for zonal plots
; arrays for full 3D data
      WORK=FLTARR(NPFULL,NT,NR)       ;3d work array
      R3D =FLTARR(NPFULL,NT,NR+NR_IC) ;3d radii        in 3D shell
	T =FLTARR(NPFULL,NT,NR)       ;temperature     in 3D shell
	H =FLTARR(NPFULL,NT,NR)       ;helicity        in 3D shell
	VR=FLTARR(NPFULL,NT,NR)       ;radial velocity in 3D shell
	VT=FLTARR(NPFULL,NT,NR)       ;theta  velocity in 3D shell
	VP=FLTARR(NPFULL,NT,NR)       ;phi    velocity in 3D shell
     dVPdz=FLTARR(NPFULL,NT,NR)       ;partial derivative of VP in order to z
      dTdt=FLTARR(NPFULL,NT,NR)       ;partial derivative of T in order to theta
	WR=FLTARR(NPFULL,NT,NR)       ;radial vorticity in 3D shell
	WT=FLTARR(NPFULL,NT,NR)       ;theta  vorticity in 3D shell
	WP=FLTARR(NPFULL,NT,NR)       ;phi    vorticity in 3D shell
	B2=FLTARR(NPFULL,NT,NR+NR_IC) ;rms    B^2      in 3D shell
	BR=FLTARR(NPFULL,NT,NR+NR_IC) ;radial field    in 3D shell
	BT=FLTARR(NPFULL,NT,NR+NR_IC) ;theta  field    in 3D shell
	BP=FLTARR(NPFULL,NT,NR+NR_IC) ;phi    field    in 3D shell
	BX=FLTARR(NPFULL,NT,NR+NR_IC) ;radial field    in 3D shell
	BY=FLTARR(NPFULL,NT,NR+NR_IC) ;radial field    in 3D shell
	BZ=FLTARR(NPFULL,NT,NR+NR_IC) ;radial field    in 3D shell
	BS3D=FLTARR(NPFULL,NT,NR+NR_IC) ;radial field    in 3D shell
	BMagR=FLTARR(NR+NR_IC)        ;magnetic energy over radius (1D)
        OMEG=FLTARR(NPFULL,NT,NR)     ;omega field in 3D shell
        WRAP=FLTARR(NPFULL+1,NT)      ; augmented array on spherical surface

        TMEAN = FLTARR(NR)            ; mean temperature as fct of radius
        AERR  = FLTARR(NR+NR_IC)      ; integration error for field lines
        PERR  = FLTARR(NR)            ; integration error for stream lines

        NLV=17 + IZERO                ;contouring levels
        NLVV=17 + IZERO               ;contouring levels Vr
        NLVA=29                       ;meridional current contours
        NLVP=13                       ;meridional stream line contours

        TV=FLTARR(NLV)                ;temperature contouring levels
        VV=FLTARR(NLVV)               ;velocity contouring levels
        BV=FLTARR(NLV)                ;magn. field contouring levels
        HV=FLTARR(NLV)                ;helicity contouring levels

        NARROW=28                     ; # of points for drawing arrows
        ARRST=[2./NARROW,2./NARROW]
        RFAC=NARROW/(NARROW+2.0)
        NSTYLE=3                      ; Line style for BR
        NLS   =3                      ; Line style for values < 0

;*******************************************************************************
;LOOP READS IN DATA ARRAYS FROM FORTRAN DATAFILE INTO IDL ARRAYS
;DIFFERENT LOOPS FOR DIFFERENT GRAPHOUT VERSIONS
;*******************************************************************************

; #1: UNSORTED GRAPHIC DATA FROM PARALLEL CODE VERSION *************************
      IF GRAPHOUTVERSION GE 2 THEN BEGIN 

; define arrays that contain data from input file:
	T1 =FLTARR(NP,NT,NR)       ;temperature     in 1/minc part of 3D shell      
	VR1=FLTARR(NP,NT,NR)       ;radial velocity in 1/minc part of 3D shell
	VT1=FLTARR(NP,NT,NR)       ;theta  velocity in 1/minc part of 3D shell
	VP1=FLTARR(NP,NT,NR)       ;phi    velocity in 1/minc part of 3D shell
	BR1=FLTARR(NP,NT,NR+NR_IC) ;radial field    in 1/minc part of 3D shell
	BT1=FLTARR(NP,NT,NR+NR_IC) ;theta  field    in 1/minc part of 3D shell
	BP1=FLTARR(NP,NT,NR+NR_IC) ;phi    field    in 1/minc part of 3D shell
	BR2=FLTARR(NP,NT,NR+NR_IC) ;radial field    in 1/minc part of 3D shell
	BT2=FLTARR(NP,NT,NR+NR_IC) ;theta  field    in 1/minc part of 3D shell
	BP2=FLTARR(NP,NT,NR+NR_IC) ;phi    field    in 1/minc part of 3D shell

;       PRINT, FORMAT='($, "reading radial level (i1-i2) ... OC: ")'
        PRINT, FORMAT='($, "reading data OC: ")'
	DUMMY=FLTARR(NP)
        DUMMYI=FLTARR(4)

        FOR IR=1,NR*NFSIZE DO BEGIN ;*** LOOP OVER # OF DATA BLOCKS ************
            IF GRAPHOUTVERSION MOD 2 EQ 0 THEN $ ; version=2 or 4 or 6: ascii-data
                                                   ;         3 or 5 or 7: binary-data
                READF,1, KC,R1,i1,i2 $
                ELSE BEGIN
                    READU,1, DUMMYI
                    KC=DUMMYI(0)
                    R1=DUMMYI(1)
                    i1=DUMMYI(2)
                    i2=DUMMYI(3)
                ENDELSE

;	  PRINT,FORMAT='($, I2, " (", I3, "-", I3, "), ")', KC, i1, i2
                PRINT,FORMAT='($,".")'                  
                R(KC)=R1
                REQ(*,KC)=R1
                RZON(*,KC)=R1
                R3D(*,*,KC)=R1

     	    FOR i=i1,i2 DO BEGIN                                 ;READ TEMPERATURE
              READU,1,DUMMY
              T1(*,i-1,KC)=DUMMY
            ENDFOR
       	    FOR i=i1,i2 DO BEGIN                                  ;READ VELOCITIES
              READU,1,DUMMY
              VR1(*,i-1,KC)=DUMMY
            ENDFOR
       	    FOR i=i1,i2 DO BEGIN
               READU,1,DUMMY
               VT1(*,i-1,KC)=DUMMY
               ENDFOR
             FOR i=i1,i2 DO BEGIN
               READU,1,DUMMY
               VP1(*,i-1,KC)=DUMMY
             ENDFOR

            IF IMAG GT 0 THEN BEGIN              ;READ MAGNETIC FIELD (IF PRESENT)
      	      FOR i=i1,i2 DO BEGIN
                READU,1,DUMMY
                BR1(*,i-1,KC)=DUMMY
 	        ENDFOR
      	      FOR i=i1,i2 DO BEGIN
                READU,1,DUMMY
                BT1(*,i-1,KC)=DUMMY
              ENDFOR
     	      FOR i=i1,i2 DO BEGIN
                READU,1,DUMMY
                BP1(*,i-1,KC)=DUMMY
              ENDFOR
            ENDIF

        ENDFOR ;*** END OF LOOP OVER # OF DATA BLOCKS **************************


;THIS PART IS FOR THE INNER-CORE DATA (ONLY MAGNETIC FIELD) ********************
 
        IF GRAPHOUTVERSION GE 4 AND NN_IC GT 0 AND IMAG THEN BEGIN

          PRINT
          PRINT, FORMAT='($, "IC: ")'
	  FOR IR=1,NR_IC  DO BEGIN  ;*** LOOP OVER # OF DATA BLOCKS IC *********
            READU,1, KCf,R1,i1f,i2f
            KC=FIX(KCf)
            i1=FIX(i1f)
            i2=FIX(i2f)

 	    PRINT,FORMAT='($,".")'           
       
	    R(KC)=R1
	    REQ(*,KC)=R1
	    RZON(*,KC)=R1
	    R3D(*,*,KC)=R1

     	    FOR i=i1,i2 DO BEGIN
               READU,1,DUMMY
               BR1(*,i-1,KC)=DUMMY
 	    ENDFOR
      	    FOR i=i1,i2 DO BEGIN
               READU,1,DUMMY
               BT1(*,i-1,KC)=DUMMY
            ENDFOR
     	    FOR i=i1,i2 DO BEGIN
               READU,1,DUMMY
               BP1(*,i-1,KC)=DUMMY
            ENDFOR
          ENDFOR
        ENDIF
;*******************************************************************************


;THE PARALLEL VERSION OF MAG-CODE WRITES THETA GRID POINTS ALTERNATELY FOR
;NORTHERN AND SOUTHERN HEMISPHERE. NOW THE DATA HAS TO BE SORTED:

       PRINT
       PRINT, FORMAT='($, "sorting OC: ")'
       DUMMY2=FLTARR(NP,NT)

       FOR IR=0, NR-1 DO BEGIN ;*** LOOP OVER RADIAL LEVELS (OC) ***************
         PRINT,FORMAT='($,".")'                  

         DUMMY2=T1(*,*,IR)
         FOR i=0,NT/2-1 DO BEGIN                               ;SORT TEMPERATURE
           T1(*,i,IR)=DUMMY2(*,2*i)
           T1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
         ENDFOR 

         DUMMY2=VR1(*,*,IR)
         FOR i=0,NT/2-1 DO BEGIN                                  ;SORT VELOCITY  
           VR1(*,i,IR)=DUMMY2(*,2*i)
           VR1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
         ENDFOR         
         DUMMY2=VP1(*,*,IR)
         FOR i=0,NT/2-1 DO BEGIN
           VP1(*,i,IR)=DUMMY2(*,2*i)
           VP1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
         ENDFOR  
         DUMMY2=VT1(*,*,IR)
         FOR i=0,NT/2-1 DO BEGIN
           VT1(*,i,IR)=DUMMY2(*,2*i)
           VT1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
         ENDFOR 
 
 	 IF IMAG GT 0 THEN BEGIN               ;SORT MAGNETIC FIELD (IF PRESENT)
           DUMMY2=BR1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BR1(*,i,IR)=DUMMY2(*,2*i)
             BR1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
           ENDFOR         
           DUMMY2=BP1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BP1(*,i,IR)=DUMMY2(*,2*i)
             BP1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
           ENDFOR  
           DUMMY2=BT1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BT1(*,i,IR)=DUMMY2(*,2*i)
             BT1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
           ENDFOR  
	 ENDIF
       ENDFOR ;*** END OF LOOP OVER RADIAL LEVELS (OC) *************************

       IF GRAPHOUTVERSION GE 4 AND NN_IC GT 0 THEN BEGIN
         PRINT,FORMAT='($," IC: ")'
         FOR IR=NR, NR+NR_IC-1 DO BEGIN   ;*** LOOP OVER RADIAL LEVELS (IC) ****
           PRINT,FORMAT='($,".")'
                  
           DUMMY2=BR1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BR1(*,i,IR)=DUMMY2(*,2*i)
             BR1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
             ENDFOR         

           DUMMY2=BP1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BP1(*,i,IR)=DUMMY2(*,2*i)
             BP1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
             ENDFOR  

           DUMMY2=BT1(*,*,IR)
           FOR i=0,NT/2-1 DO BEGIN
             BT1(*,i,IR)=DUMMY2(*,2*i)
             BT1(*,i+NT/2,IR)=DUMMY2(*,NT-1-2*i)
             ENDFOR

           ENDFOR ;* END OF LOOP OVER RADIAL LEVELS (IC) ***********************
         ENDIF

;SET LOWEST T TO ZERO
         i=MIN(T1)
         T1=T1-i

;TAKE NSYM-FOLD SYMMETRY INTO ACCOUNT:
       IF NSYM EQ 1 THEN BEGIN
         T= T1
         VR=VR1
  	 VP=VP1
         VT=VT1
	 IF IMAG GT 0 THEN BEGIN
           BR=BR1
	   BP=BP1
           BT=BT1
         ENDIF
       END $
       ELSE BEGIN
         PRINT
         PRINT, FORMAT='($, "symmetry ")'
         FOR i=0,(NSYM-1)*NP,NP DO BEGIN
           PRINT, FORMAT='($, ".")'
           T (i:NP-1+i,*,*)=T1 (0:NP-1,*,*)
           VR(i:NP-1+i,*,*)=VR1(0:NP-1,*,*)
           VP(i:NP-1+i,*,*)=VP1(0:NP-1,*,*)
           VT(i:NP-1+i,*,*)=VT1(0:NP-1,*,*)
 	   IF IMAG GT 0 THEN BEGIN
             BR(i:NP-1+i,*,*)=BR1(0:NP-1,*,*)
             BP(i:NP-1+i,*,*)=BP1(0:NP-1,*,*)
             BT(i:NP-1+i,*,*)=BT1(0:NP-1,*,*)
           ENDIF
         ENDFOR
       ENDELSE

;MEAN TEMPERATURE (APPROXIMATE) AND OUTPUT OF RADII
        PRINT
        PRINT, "  level   radius    mean T"
        FOR IR=0, NR-1 DO BEGIN
          TM=0.0
          SM=0.0
          FOR IT=0,NT-1 DO BEGIN
            TM=TM+TOTAL(T1(*,IT,IR))*STHET(IT)
            SM=SM+STHET(IT)
            ENDFOR
          TMEAN(IR)=TM/(SM*NP)
          PRINT, FORMAT='(I5, F11.5, F11.5)', IR, R[IR], TMEAN(IR)
          ENDFOR 
        FOR IR=NR, NR+NR_IC-1 DO BEGIN
          IF IR EQ NR THEN PRINT, " --- inner core ---"
          PRINT, FORMAT='(I5, F11.5)', IR, R[IR]          
          END     
  
      END $

; #2: SORTED GRAPHIC DATA FROM SEQUENTIAL CODE VERSION *************************

      ELSE BEGIN

; define arrays that contain data from input file:
	T1 =FLTARR(NP,NT)        ;temperature     on spherical 2D surface 
	VR1=FLTARR(NP,NT)        ;radial velocity on spherical 2D surface 
	VT1=FLTARR(NP,NT)        ;theta  velocity on spherical 2D surface 
	VP1=FLTARR(NP,NT)        ;phi    velocity on spherical 2D surface 
	BR1=FLTARR(NP,NT)        ;radial field    on spherical 2D surface  
	BT1=FLTARR(NP,NT)        ;theta  field    on spherical 2D surface 
	BP1=FLTARR(NP,NT)        ;phi    field    on spherical 2D surface

	FOR IR=0,NR-1 DO BEGIN ;*** LOOP OVER RADIAL LEVELS ********************
	  READF,1, format='(i4,e13.5)', KC, R1 
          PRINT, IR,R1
          R(IR)=R1
          REQ(*,IR)=R1
          RZON(*,IR)=R1
          R3D(*,*,IR)=R1
	  READF,1,T1
	  READF,1,VR1
	  READF,1,VT1
	  READF,1,VP1
          READF,1,BR1
	  READF,1,BT1
	  READF,1,BP1

;         MEAN TEMPERATURE (APPROXIMATE)
          TM=0.0
          SM=0.0
          FOR IT=0,NT-1 DO BEGIN
            TM=TM+TOTAL(T1(*,IT))*STHET(IT)
            SM=SM+STHET(IT)
          ENDFOR
          TMEAN(IR)=TM/(SM*NP)

;         MAKES 3D SHELL ARRAYS FROM 
;         2D SPHERICAL SURFACE ARRAYS:
  	  FOR J=0,NP-1 DO BEGIN
            FOR NS=1,NSYM DO BEGIN
              JNS=(NS-1)*NP+J
       	      T(JNS,*,IR) = T1(J,*)
    	      VR(JNS,*,IR)= VR1(J,*)
    	      VT(JNS,*,IR)= VT1(J,*)
    	      VP(JNS,*,IR)= VP1(J,*)
    	      BR(JNS,*,IR)= BR1(J,*)
    	      BT(JNS,*,IR)= BT1(J,*)
    	      BP(JNS,*,IR)= BP1(J,*)   
            ENDFOR
          ENDFOR        

	ENDFOR ;*** END OF LOOP OVER RADIAL LEVELS ****************************

      ENDELSE
      CLOSE,1   ;CLOSE DATA FILE

;NOW WE GOT THE DATA **********************************************************

;NON-DIMENSIONALIZES RADII VALUES SUCH THAT R_CMB = 1.0
	RADTOP=R(0)
	R=R/RADTOP
        RSCALE=RADTOP-R(NR-1)
        REQ=REQ/RADTOP
        RZON=RZON/RADTOP
        R3D=R3D/RADTOP

;FORM X-DISTANCE ARRAY
	XZON=RZON*(ABS(COS(THZON)))

;FORM INNER CORE POLYGON ARRAYS
       YIC=(RZON(*,NR-1)-0.004)*SIN(THZON(*,NR-1))
       XIC=(RZON(*,NR-1)-0.004)*ABS(COS(THZON(*,NR-1)))
       XIC(0)= -0.004 & XIC(NT-1)=-0.004

;FORM AZIMUTHAL AVERAGES
	TZON=TOTAL(T,1)/NPFULL
	BZON=TOTAL(BP,1)/NPFULL
        BZON(*,0)=.0
        BZON(*,NR+NR_IC-1)=.0
	BRZON=TOTAL(BR,1)/NPFULL
	BTZON=TOTAL(BT,1)/NPFULL
	VZON=TOTAL(VP,1)/NPFULL
	VRZON=TOTAL(VR,1)/NPFULL
	VTZON=TOTAL(VT,1)/NPFULL
	ROTZON=VZON/XZON

;FORM PARTIAL DERIVATIVES FOR THE THERMAL WIND PLOT (T =FLTARR(NPFULL,NT,NR))

               dVPdz = COS(TH3D)*(SHIFT(VP,0,0,-1)-SHIFT(VP,0,0,1))/$
                       (SHIFT(R3D,0,0,-1)-SHIFT(R3D,0,0,1))$
                     - (SIN(TH3D)/R3D)*(SHIFT(VP,0,-1,0)-SHIFT(VP,0,1,0))/$
                       (SHIFT(TH3D,0,-1,0)-SHIFT(TH3D,0,1,0))
               dVPdzZON=TOTAL(dVPdz,1)/NPFULL
               dVPdzMER=TOTAL(dVPdz,2)/NT
               dVPdzMER2=(dVPdzMER(0:NPFULL/2-1,*)+dVPdzMER(NPFULL/2:NPFULL-1,*))/2

               dTdt = ((SHIFT(T,0,-1,0)-SHIFT(T,0,1,0))/$
	              (SHIFT(TH3D,0,-1,0)-SHIFT(TH3D,0,1,0)));*$
;                      ((RA*EK)/(2*Pr))
;               dTdt = dTdt/R3D
               dTdtZON=TOTAL(dTdt,1)/NPFULL
               dTdtMER=TOTAL(dTdt,2)/NT
               dTdtMER2=(dTdtMER(0:NPFULL/2-1,*)+dTdtMER(NPFULL/2:NPFULL-1,*))/2
               TMER=TOTAL(T,2)/NT
               TMER2=(TMER(0:NPFULL/2-1,*)+TMER(NPFULL/2:NPFULL-1,*))/2

;FORM VORTICITY VECTOR COMPONENTS

	       WP = (SHIFT(VT,0,0,-1)-SHIFT(VT,0,0,1))/$
	            (SHIFT(R3D,0,0,-1)-SHIFT(R3D,0,0,1))$
	          + VT/R3D $
	          - (SHIFT(VR,0,-1,0)-SHIFT(VR,0,1,0))/$
	            (R3D*(SHIFT(TH3D,0,-1,0)-SHIFT(TH3D,0,1,0)))

	       WT = (SHIFT(VR,-1,0,0)-SHIFT(VR,1,0,0))/$
	            (R3D*SIN(TH3D)*PI4NP)$
	          - (SHIFT(VP,0,0,-1)-SHIFT(VP,0,0,1))/$
	            (SHIFT(R3D,0,0,-1)-SHIFT(R3D,0,0,1))$
	          - VP/R3D
	       WR = (SHIFT(VP,0,-1,0)-SHIFT(VP,0,1,0))/$
	            (R3D*(SHIFT(TH3D,0,-1,0)-SHIFT(TH3D,0,1,0)))$
	          + COS(TH3D)*VP/(R3D*SIN(TH3D)) $
	          - (SHIFT(VT,-1,0,0)-SHIFT(VT,1,0,0))/$
	            (R3D*SIN(TH3D)*PI4NP)

	WR(*,0,*)=WR(*,1,*) & WR(*,NT-1,*)=WR(*,NT-2,*)
        WT(*,*,0)=WT(*,*,1) & WT(*,*,NR-1)=WT(*,*,NR-2)
        WP(*,*,0)=WP(*,*,1) & WP(*,*,NR-1)=WP(*,*,NR-2)

;real BR, BT, BP:

	BR2=BR*R3D*R3D
	BT2=BT*R3D*sin(TH3D)
	BP2=BP*R3D*sin(TH3D)

;FORM HELICITY
	H = WR*VR + WT*VT + WP*VP

;ZONAL AVERAGE
        HZON=TOTAL(H,1)/NPFULL
        HBZON=-TOTAL((WR*VR+WT*VT)*BP,1)/NPFULL

;FIELD LINES OF ZONALLY AVERAGED MAGNETIC AND VELOCITY FIELD
        A(0,*)=-0.5*XZON(0,*)*BRZON(0,*)*THETA(0)
        PSI(0,*)=-0.5*XZON(0,*)*VRZON(0,*)*THETA(0)

	FOR IT=1,NT-1 DO BEGIN
        	A(IT,*)=   A(IT-1,*)$
                 - (XZON(IT,*)*BRZON(IT,*) + XZON(IT-1,*)*BRZON(IT-1,*))$
        	   *(THZON(IT,*)-THZON(IT-1,*))

        	PSI(IT,*)=PSI(IT-1,*) - (XZON(IT,*)*VRZON(IT,*)$
        	                       + XZON(IT-1,*)*VRZON(IT-1,*))$
                	    *(THZON(IT,*)-THZON(IT-1,*))
	ENDFOR

; ERROR DETERMINATION
        AERR=A(NT-1,*)-0.5*(XZON(NT-1,*)*BRZON(NT-1,*))*(PI-THETA(NT-1))
        PERR=PSI(NT-1,*)-0.5*(XZON(NT-1,*)*VRZON(NT-1,*))*(PI-THETA(NT-1))
; ERROR DISTRIBUTED LINEARLY
        FOR IT=0,NT-1 DO BEGIN
           A(IT,*)=A(IT,*)-AERR(*)*THETA(IT)/PI
           PSI(IT,*)=PSI(IT,*)-PERR(*)*THETA(IT)/PI
        ENDFOR

        A=0.5*A/COS(THZON)
        PSI=0.5*PSI/COS(THZON)


;FIND EXTREME VALUES
        EPS=0.0000001
        TMAX=MAX(ABS(T))+EPS
        VMAX=MAX(ABS(VR))+EPS
        WMAX=MAX(ABS(WR))+EPS
        HMAX=MAX(ABS(H))+EPS
        BRMAX=MAX(ABS(BR))+EPS
        print,MAX(ABS(BP))
        print,EPS
        BPMAX=MAX(ABS(BP))+EPS
        HBMAX=MAX(ABS(HBZON))+EPS
        AMAX=MAX(ABS(A))+EPS
        PMAX=MAX(ABS(PSI))+EPS
;
;CALCULATE TEMPERATURE GRADIENTS AT IR=0 AND IR=NR-1 (HEATFLOW)
;
        TMAXOUT=MAX(T(*,*,1)-T(*,*,0))
        TMAXIN =MAX(T(*,*,NR-1)-T(*,*,NR-2))
        HEATFLOW_CMB=(T(*,*,1)-T(*,*,0))/TMAXOUT
;        openw, 44, 'test'
;	printf, 44, heatflow_cmb
        HEATFLOW_ICB=(T(*,*,NR-1)-T(*,*,NR-2))/TMAXIN

;
        VSTEP=2*VMAX/(NLVV-2)
        WSTEP=2*WMAX/(NLV-2)
        HSTEP=2*HMAX/(NLVV-2)
        BRSTEP=2*BRMAX/(NLV-2)
        BPSTEP=2*BPMAX/(NLV-2)
        HBSTEP=2*HBMAX/(NLV-2)
;TSTEP based on mean temperature at ICB, works better for flux
;boundary conditions [old: TSTEP=TMAX/(NLV-2)]
        TSTEP=(TMEAN(NR-1)+EPS)/(NLV-2)
        ASTEP=2.*AMAX/(NLVA-1)
        PSTEP=2.*PMAX/(NLVP-1)
        OMSTEP=0.0
        JSTEP=0.0
        TADD=0.0

;ASSIGN CONTOURING LEVELS
        TV=(FINDGEN(NLV)-1)*TSTEP - 0.5
        TV(0)=-100.
        TV(1)=-0.5+EPS                     ; make sure that contours are
        TV(NLV-1)=TMAX-0.5-EPS             ; drawn for T=0 and T=1
        VV=(FINDGEN(NLVV)-1)*VSTEP-VMAX
        VV(0)=-100.*VMAX
        BRV=(FINDGEN(NLV)-1)*BRSTEP-BRMAX
        BRV(0)=-100.*BRMAX
        BPV=(FINDGEN(NLV)-1)*BPSTEP-BPMAX
        BPV(0)=-100.*BPMAX
        HV=(FINDGEN(NLVV)-1)*HSTEP-HMAX
        HV(0)=-100.*HMAX
        HBV=(FINDGEN(NLV)-1)*HBSTEP-HBMAX
        HBV(0)=-100.*HBMAX
        WV=(FINDGEN(NLV)-1)*WSTEP-WMAX
        WV(0)=-100.*WMAX
        AV=(FINDGEN(NLVA)+0.0)*ASTEP-AMAX
        PV=(FINDGEN(NLVP)+0.5)*PSTEP-PMAX

;ASSIGN CONTOUR AMPLIFICATION FACTORS
        TFAC=1.
        BFAC=1.
        VFAC=1.
        HFAC=1.
        ARRV=1.  
	ARRB=1.

;HEADLINES
        ANNTEXT=''
        ANNTEXT1=STRTRIM(RUNID,2)
        ANNTEXT2=STRING(FORMAT='("Ra=",E8.2,", Ek=",E8.2,", Pr=",E8.2,", Pm=",E8.2)', RA,EK,PR,PM)
        ANNTEXT3=STRING("time=", STRTRIM(STRING(TIME),2))

;COLORS AND PLOTTING WINDOW INFORMATION
        IF LCOLOR GT 0 THEN CTABLE=39 ELSE CTABLE=0 ; 39: rainbow, 0: b/w
        DEVICE, DECOMPOSE=0
	LOADCT,CTABLE                   ;LOADS COLOR TABLE NUMBER

IF LCOLOR GT 0 THEN BEGIN
; SELF-DEFINED COLOR TABLE
; BLACK COLOR AT START OF SPECTRUM
        R_CURR(000)= 00 & G_CURR(000)= 00 & B_CURR(000)= 00
; WHITE COLOR AT END OF SPECTRUM
        R_CURR( 23)=255 & G_CURR( 23)=255 & B_CURR( 23)=255
        R_CURR( 24)=255 & G_CURR( 24)=255 & B_CURR( 24)=255
; DARK GREY
        R_CURR( 01)= 90 & G_CURR( 01)= 90 & B_CURR( 01)= 90
; LIGHT GREY
        R_CURR( 02)=190 & G_CURR( 02)=190 & B_CURR( 02)=190

        R_CURR( 5)= 00 & G_CURR( 5)=  0 & B_CURR( 5)=  0
        R_CURR( 6)= 55 & G_CURR( 6)=  0 & B_CURR( 6)=195
        R_CURR( 7)= 20 & G_CURR( 7)=  0 & B_CURR( 7)=245
        R_CURR( 8)=  0 & G_CURR( 8)= 60 & B_CURR( 8)=255
        R_CURR( 9)=  0 & G_CURR( 9)=110 & B_CURR( 9)=255
        R_CURR(10)=  0 & G_CURR(10)=150 & B_CURR(10)=255
        R_CURR(11)= 20 & G_CURR(11)=200 & B_CURR(11)=255
        R_CURR(12)= 70 & G_CURR(12)=240 & B_CURR(12)=255
        R_CURR(13)=150 & G_CURR(13)=255 & B_CURR(13)=240
        R_CURR(14)=225 & G_CURR(14)=255 & B_CURR(14)=225

        R_CURR(15)=255 & G_CURR(15)=255 & B_CURR(15)=185
        R_CURR(16)=255 & G_CURR(16)=255 & B_CURR(16)= 75
        R_CURR(17)=255 & G_CURR(17)=225 & B_CURR(17)=  0
        R_CURR(18)=255 & G_CURR(18)=175 & B_CURR(18)=  0
        R_CURR(19)=255 & G_CURR(19)=125 & B_CURR(19)=  0
        R_CURR(20)=255 & G_CURR(20)= 75 & B_CURR(20)=  0
        R_CURR(21)=255 & G_CURR(21)=  0 & B_CURR(21)=  0
        R_CURR(22)=230 & G_CURR(22)=  0 & B_CURR(22)= 35

        TVLCT,R_CURR,G_CURR,B_CURR

  VCOLORS=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,$
           22,7]
  COLINV=[22,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,22]
  AACOLOR=[6,6,7,7,8,8,9,9,10,10,11,11,12,12,13, $
           16,17,17,18,18,19,19,20,20,21,21,22,22,22]
  LIMIT=-1.E20
  BWCOL=0
  BWCOLI=0
  ICCOL=2
  LCINV=1
END $
;                grey-scale colors
ELSE BEGIN
  VCOLORS=[135,150,165,180,195,210,225,240,255,255,240,225,210,195,180,165,150,135]
  AACOLOR=[0]
  LIMIT=0.0
  BWCOL=0
  BWCOLI=255
  ICCOL=255
  LCINV=0
END

; Define window size and plotting frames
  XDIM=18.0
   IF ISINGLE GT 0 THEN XDIM=12.0
  YDIM=24.0
  IGIFPS=0

;** WINDOW SIZE IN PIXELS (ALSO SIZE OF GIF FILE)
  XWINDOW=405
   IF ISINGLE GT 0 THEN XWINDOW=300
  YWINDOW=540
  XWIND=XWINDOW*SCWINDOW
  YWIND=YWINDOW*SCWINDOW
  CSZ=1.00 & CSB=1.5 & CSS=0.720       ; CHARACTER SIZES
  SZF=1.50                              ; CHARACTER SIZE FACTOR


;*** Normalized coordinates
  XY1=[ 0.5/XDIM,12./YDIM, 8.5/XDIM,20./YDIM] & XT1=4.5/XDIM & YT1=20.2/YDIM
  XY2=[ 9.5/XDIM,12./YDIM,17.5/XDIM,20./YDIM] & XT2=13.5/XDIM & YT2=20.2/YDIM
  XY3=[ 0.5/XDIM, 2./YDIM, 8.5/XDIM,10./YDIM] & XT3=4.5/XDIM & YT3=10.2/YDIM
  XY4=[ 9.5/XDIM, 2./YDIM,17.5/XDIM,10./YDIM] & XT4=13.5/XDIM & YT4=10.2/YDIM
  XS0=[ 1./XDIM, 2./YDIM, 11.0/XDIM,22./YDIM]
  XS1=[ 3./XDIM,11./YDIM, 7.5/XDIM,20./YDIM] & XTS1=5.0/XDIM & YTS1=20.2/YDIM
  XS2=[11./XDIM,11./YDIM,15.5/XDIM,20./YDIM] & XTS2=13./XDIM & YTS2=YTS1
  XS3=[ 3./XDIM, 0.5/YDIM, 7.5/XDIM,9.5/YDIM] & XTS3=XTS1 & YTS3=9.7/YDIM
  XS4=[11./XDIM, 0.5/YDIM,15.5/XDIM,9.5/YDIM] & XTS4=XTS2 & YTS4=YTS3
  XC1=[1.5/XDIM,12./YDIM,17./XDIM,20./YDIM]
  XC2=[1.5/XDIM, 1.5/YDIM,17./XDIM,9.5/YDIM]
  XC0=[1.5/XDIM, 1.5/YDIM,17./XDIM,20./YDIM]
  XSS=[1.5/XDIM, 4.0/YDIM,16.5/XDIM,19./YDIM]

;*** SET DEVICE AND OPEN GRAPHIC WINDOW
SET_PLOT,PLATFORM
!P.FONT=0
;SCHRIFT_NORMAL='-adobe-helvetica-medium-r-normal--14-100-100-100-p-76-iso8859-1'
;SCHRIFT_GROSS='-adobe-helvetica-bold-r-normal--17-120-100-100-p-92-iso8859-1'
;DEVICE, SET_FONT = SCHRIFT_NORMAL
DEVICE,DECOMPOSE=0,RETAIN=2
!P.BACKGROUND=255 ;background color=white
WINDOWCOUNT=0


WINDOW,WINDOWCOUNT,xsize=XWIND,ysize=YWIND,title=VERSION
ERASE


;******************************************************************************
;**************  HERE THE PLOTTING MENU STARTS  **************************
;******************************************************************************

LABEL0: 
        !P.BACKGROUND=255 ;background color=white
        !P.THICK=1.75
        !P.CHARTHICK=1
        !P.NOERASE=1
        !P.MULTI=[0,0,2]
        !P.FONT=0

        print,"VV:"
        print,VV
        print,"TV:"
        print,TV
        print,"BPV:"
        print,BPV
        print,"WV:"
        print,WV

        IF IGIFPS EQ 1 THEN GOTO,LABELOUT
        PRINT, "map=1, closeup=2, equator=3, slice=4, lon. average=5/6, Bcmb=7, full map=8"
        PRINT, "change contour steps=10, window control=11, ", $
               "arrowlength=12, annotations=13"
        PRINT, FORMAT='($, "TW plot=15, show contour steps=14, print PS=21, print GIF=22, quit=-1")'
        READ, IOPTION

        CASE IOPTION OF
       -1: GOTO, LABEL99
        0: GOTO, LABEL0
        1: GOTO, LABEL1
        2: GOTO, LABEL2
        3: GOTO, LABEL3
        4: GOTO, LABEL4
        5: GOTO, LABEL5
        6: GOTO, LABEL6
        7: GOTO, LABEL7
        8: GOTO, LABEL8
        10: GOTO, LABEL10
        11: GOTO, LABEL11
        12: GOTO, LABEL12
        13: GOTO, LABEL13
        14: GOTO, LABEL14
        15: GOTO, LABEL15
        21: GOTO, LABELOUT
        22: GOTO, LABELOUT
        ELSE: GOTO, LABEL0
        ENDCASE


LABEL10:  PRINT, FORMAT='("current factors for T,v,B,H&W",4F7.2)',$
                TFAC,VFAC,BFAC,HFAC
          PRINT, FORMAT='($, "new factors?")'
          READ, TFAC,VFAC,BFAC,HFAC
          GOTO, LABEL0

LABEL11:  PRINT, FORMAT='($, "scale window=1, add window=2")'
          READ, IOPTION
          IF IOPTION EQ 1 THEN BEGIN
            PRINT, FORMAT='("window scale factor?  current=",F7.3)',SCWINDOW
            READ, SCWINDOW
            XWIND=XWINDOW*SCWINDOW & YWIND=YWINDOW*SCWINDOW
            CSZ=SCWINDOW*0.8 & CSB=1.12*SCWINDOW & CSS=0.60*SCWINDOW
            WINDOW,WINDOWCOUNT,xsize=XWIND,ysize=YWIND,title=VERSION
            ERASE
          ENDIF
          IF IOPTION EQ 2 THEN BEGIN
            WINDOWCOUNT=WINDOWCOUNT+1
            WINDOW,WINDOWCOUNT,xsize=XWIND,ysize=YWIND,title=VERSION
            ERASE
          ENDIF
         
          GOTO, LABEL0

LABEL12:  PRINT, FORMAT='("arrow length v / B ?  current=",2F7.3)',ARRV,ARRB
          READ, ARRV,ARRB
          GOTO, LABEL0

LABEL13:  PRINT, "new annotation text: enter line-# (1-3) [blank] text"
          READ, ITEXTLINE,ANNTEXT
          CASE ITEXTLINE OF
          1: ANNTEXT1=ANNTEXT
          2: ANNTEXT2=ANNTEXT
          3: ANNTEXT3=ANNTEXT
          ELSE:  GOTO,LABEL0
          ENDCASE
          GOTO, LABEL0

LABEL14:  PRINT," current contouring interval settings: "
          PRINT," radial velocity=   ",VSTEP/VFAC
          PRINT," azimuthal velocity=",VSTEP/HFAC
          PRINT," radial field=      ",BRSTEP/BFAC
          PRINT," azimuthal field=   ",BPSTEP/BFAC
          PRINT," vorticity=         ",WSTEP*RSCALE/HFAC
          PRINT," helicity=          ",HSTEP*RSCALE/HFAC
          PRINT," angular velocity=  ",3*VSTEP/VFAC
          PRINT," temperature (maps)=",TSTEP/TFAC,"  add= ",TADD
          PRINT," max temp. grad. =  ",TMAXOUT/(R(0)-R(1)),$
                                       TMAXIN/(R(NR-2)-R(NR-1))
          PRINT," toroidal current=  ",JSTEP*RSCALE/TFAC
          PRINT," omega=             ",OMSTEP/VFAC
          PRINT," helicity * B_phi=  ",HBSTEP*RSCALE/HFAC
          PRINT," psi=               ",PSTEP*ARRV

          PRINT," "
          PRINT," MAXIMA: TOTAL/ZONAL"
          PRINT,"    Radial velocity   =",VMAX," / ",MAX(ABS(VRZON))
          PRINT,"    Azimuthal velocity=",MAX(ABS(VP))," / ",MAX(ABS(VZON))
          PRINT,"    Radial    field   =",BRMAX," / ",MAX(ABS(BRZON))
          PRINT,"    Azimuthal field   =",BPMAX," / ",MAX(ABS(BZON))
          PRINT,"    Radial Vorticity  =",WMAX
          PRINT,"    Helicity          =",HMAX
          PRINT,"    Zonal alpha       =",HBMAX

          GOTO, LABEL0

; ********************************************************************************
; ***************  PLOTTING OF GLOBAL MAP PROJECTION  *********************
; ********************************************************************************

LABEL1: PRINT, FORMAT='($, "radial level, level for B&T, inclination, longitude?")'
        READ,NRAD,NRB,INCL,LONSHIFT
        IF (NRAD LT 0 OR NRAD GE NR) OR (NRB LT 0 OR NRB GE NR+NR_IC) THEN BEGIN
          PRINT, "'radial level' or 'level for B&T' out of range"
          GOTO, LABEL1
          ENDIF
        IF INCL GT 90.0 OR INCL LT -90 THEN BEGIN
          PRINT, "'inclination' out of range"
          GOTO, LABEL1
          ENDIF
        IHEATFLOW=0
        IF NRB EQ 0 OR NRB EQ (NR-1) THEN BEGIN
          PRINT, FORMAT='($, "temperature=0   heat flow=1 ?")'
          READ, IHEATFLOW
          ENDIF
        PRINT, FORMAT='($, "helicity=0  z-vorticity=1 ?")'
        READ, IVORT
;       PRINT, "RADIAL VALUES=0  Z-VALUES=1    ?"
;       READ,IRADZ
        IRADZ=0
        IPAGE=1

LABEL1A: ERASE

;DEFINES IDL MAPPING PROJECTION

;MAP TEMPERATURE (OR HEATFLOW)
!P.POSITION=XY1
        IF NRB GE NR THEN NRT=NRAD ELSE NRT=NRB  ;make shure that T is
	;                                         mapped in the outer core
        RADSTRING=STRING(FORMAT='("r=",F5.3)', R(NRT))
        IF IHEATFLOW GT 0 THEN BEGIN
          PTITLE='heat flow at ' + RADSTRING
          IF NRT EQ 0    THEN BEGIN
            WRAP(0:NPFULL-1,*)=HEATFLOW_CMB
            WRAP(NPFULL,*)=HEATFLOW_CMB(0,*)
            ENDIF
          IF NRT EQ NR-1 THEN BEGIN
            WRAP(0:NPFULL-1,*)=HEATFLOW_ICB
            WRAP(NPFULL,*)=HEATFLOW_ICB(0,*)
            ENDIF
          END $
        ELSE BEGIN
          PTITLE='temperature at ' + RADSTRING
          WRAP(0:NPFULL-1,*)=T(*,*,NRT) & WRAP(NPFULL,*)=T(0,*,NRT)
        ENDELSE

        MAP_SET,INCL,LONSHIFT,0,/ORTHOGRAPHIC,/NOBORDER,/NOERASE
        TADD=0.5/TFAC+TMEAN(NRT)*(TFAC-1.)/TFAC

      IF LCOLOR EQ 0 THEN $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=TV/TFAC+TADD,$
         COLOR=0,/OVERPLOT $
      ELSE $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=TV/TFAC+TADD,$
         /CELL_FILL,$
         C_COLORS=LCINV*VCOLORS,/OVERPLOT
 	MAP_GRID,COLOR=0,LATDEL=30,LONDEL=45,GLINETHICK=1

XYOUTS,XT1,YT1,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL

;MAP RADIAL VELOCITY
!P.POSITION=XY2
       RADSTRING=STRING(FORMAT='("r=",F5.3)', R(NRAD))
       IF IRADZ LT 1 THEN BEGIN
        WRAP(0:NPFULL-1,*)=VR(*,*,NRAD) & WRAP(NPFULL,*)=VR(0,*,NRAD)
        PTITLE='radial velocity at ' + RADSTRING
       END $
       ELSE BEGIN
        WRAP(0:NPFULL-1,*)=VR(*,*,NRAD)*CTHET2-VT(*,*,NRAD)*STHET2
        WRAP(NPFULL,*)=WRAP(0,*)
        PTITLE='z-velocity / ' + RADSTRING
       END
 	MAP_SET,INCL,LONSHIFT,0,/ORTHOGRAPHIC,/ADVANCE,/NOBORDER,$
                    /NOERASE
       IF LCOLOR EQ 0 THEN $
  	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=VV/VFAC,/OVERPLOT,$
         COLOR=0,C_LINESTYLE=NLS*(VV LT LIMIT) $
       ELSE $
  	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=VV/VFAC,/OVERPLOT,$
         /CELL_FILL,$
         C_COLORS=LCINV*VCOLORS,C_LINESTYLE=NLS*(VV LT LIMIT)
 	MAP_GRID,COLOR=0,LATDEL=30,LONDEL=45,GLINETHICK=1
XYOUTS,XT2,YT2,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL

;MAP RADIAL FIELD
!P.POSITION=XY3
       RADSTRING=STRING(FORMAT='("r=",F5.3)', R(NRB))
       IF IRADZ LT 1 THEN BEGIN
        WRAP(0:NPFULL-1,*)=BR(*,*,NRB) & WRAP(NPFULL,*)=BR(0,*,NRB)
        PTITLE='radial field at ' + RADSTRING
       END $
       ELSE BEGIN
        WRAP(0:NPFULL-1,*)=BR(*,*,NRB)*CTHET2-BT(*,*,NRB)*STHET2
        WRAP(NPFULL,*)=WRAP(0,*)
        PTITLE='z-field at ' + RADSTRING
       END
 	MAP_SET,INCL,LONSHIFT,0,/ORTHOGRAPHIC,/ADVANCE,/NOBORDER,/NOERASE
       IF LCOLOR EQ 0 THEN $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=BRV/BFAC,/OVERPLOT,$
         COLOR=0,C_LINESTYLE=NLS*(BRV LT LIMIT) $
       ELSE $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=BRV/BFAC,/OVERPLOT,$
         /CELL_FILL,$
         C_COLORS=LCINV*VCOLORS,C_LINESTYLE=NLS*(BRV LT LIMIT)
       MAP_GRID,COLOR=0,LATDEL=30,LONDEL=45,GLINETHICK=1
       XYOUTS,XT3,YT3,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL

;MAP HELICITY / Z-VORTICITY
!P.POSITION=XY4
      RADSTRING=STRING(FORMAT='("r=",F5.3)', R(NRAD))
      IF IVORT LT 1 THEN BEGIN
        PTITLE='helicity at '  + RADSTRING
        WRAP(0:NPFULL-1,*)=H(*,*,NRAD) & WRAP(NPFULL,*)=H(0,*,NRAD)
	MAP_SET,INCL,LONSHIFT,0,/ORTHOGRAPHIC,/ADVANCE,/NOBORDER,/NOERASE
       IF LCOLOR EQ 0 THEN $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=HV/HFAC,/OVERPLOT,$
          COLOR=0,C_LINESTYLE=NLS*(HV LT LIMIT) $
       ELSE $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=HV/HFAC,/OVERPLOT,$
          /CELL_FILL,$
          C_COLORS=LCINV*VCOLORS,C_LINESTYLE=NLS*(HV LT LIMIT)
      END $
      ELSE BEGIN
        PTITLE='z-vorticity at ' + RADSTRING
        WRAP(0:NPFULL-1,*)=WR(*,*,NRAD)*CTHET2(*,*)-WT(*,*,NRAD)*STHET2(*,*)
            WRAP(NPFULL,*)=WRAP(0,*)
	MAP_SET,INCL,LONSHIFT,0,/ORTHOGRAPHIC,/ADVANCE,/NOBORDER,/NOERASE
       IF LCOLOR EQ 0 THEN $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=WV/HFAC,/OVERPLOT,$
          COLOR=0,C_LINESTYLE=NLS*(WV LT LIMIT) $
       ELSE $
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=WV/HFAC,/OVERPLOT,$
          /CELL_FILL,$
          C_COLORS=LCINV*VCOLORS
      END

      XYOUTS,XT4,YT4,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL
      MAP_GRID,COLOR=0,LATDEL=30,LONDEL=45,GLINETHICK=1

      GOTO, LABELTXT


; *********************************************************************************
; *********** VELOCITIES AND MAGN. FIELD ON CLOSEUP MAP *******************
; *********************************************************************************

LABEL2: PRINT,"NRAD LONMIN/MAX LATMIN/MAX STEP IV(V=1 B=2 V&B=0)"
        READ,NRAD,LONMIN,LONMAX,LATMIN,LATMAX,IASTEP,IV
        IPAGE=2
;DETERMINE ARRAY BOUNDS

IPMIN=FIX(NPFULL*(180+LONMIN)/360)  & IPMN1=IPMIN-1
IPMAX=FIX(NPFULL*(180+LONMAX)/360)  & IPMX1=IPMAX+1
ITMIN=FIX(NT*(90-LATMAX)/180-0.5)   & ITMN1=ITMIN-1
ITMAX=FIX(NT*(90-LATMIN)/180-0.5)   & ITMX1=ITMAX+1
IF IPMN1 LT 0 OR IPMX1 GE NPFULL OR NRAD LT 0 THEN BEGIN
  PRINT,'bad range'
  GOTO, LABEL2
ENDIF

IPAMAX=(IPMAX-IPMIN+1)/IASTEP-1
ITAMAX=(ITMAX-ITMIN+1)/IASTEP-1
IAP=INDGEN(IPAMAX)*IASTEP+IPMIN
IAT=INDGEN(ITAMAX)*IASTEP+ITMIN
XTEXT=LON((IPMIN+IPMAX)/2)
YTEXT=LAT(ITMIN)+0.12*(LATMAX-LATMIN)

LABEL2A:  ERASE

IF IV EQ 2 THEN GOTO,LABEL2B
IF IV EQ 1 THEN XC=XC0 ELSE XC=XC1
!P.POSITION=XC
IF LCOLOR LT 1 THEN $
CONTOUR,VR(IPMN1:IPMX1,ITMN1:ITMX1,NRAD),LON(IPMN1:IPMX1),LAT(ITMN1:ITMX1),$
            LEVELS=VV/VFAC,COLOR=064,XSTYLE=5,YSTYLE=5,C_LINESTYLE=NLS*(VV LT 0) $
ELSE $
CONTOUR,VR(IPMN1:IPMX1,ITMN1:ITMX1,NRAD),LON(IPMN1:IPMX1),LAT(ITMN1:ITMX1),$
            LEVELS=VV/VFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=5,YSTYLE=5


VELOVECT,VP(IAP,IAT,NRAD),-VT(IAP,IAT,NRAD),$
         LON(IAP),LAT(IAT),$
         LENGTH=ARRV,COLOR=0
XYOUTS,9.25/XDIM,20.3/YDIM,'velocity',SIZE=CSZ*SZF,/NORMAL,$
         ALIGNMENT=0.5,COLOR=0
XYOUTS,9.25/XDIM,XC(1)-1.2/YDIM,'longitude',SIZE=CSS*SZF,/NORMAL,$
         ALIGNMENT=0.5,COLOR=0
XYOUTS,0.60/XDIM,(XC(1)+XC(3))/2.,'latitude',SIZE=CSS*SZF,/NORMAL,$
         ALIGNMENT=0.5,ORIENTATION=-90,COLOR=0
IF IV GT 0 THEN GOTO,LABEL2C

LABEL2B:  IF IV EQ 0 THEN XC=XC2 ELSE XC=XC0
!P.POSITION=XC

IF LCOLOR LT 1 THEN $
CONTOUR,BR(IPMN1:IPMX1,ITMN1:ITMX1,NRAD),LON(IPMN1:IPMX1),LAT(ITMN1:ITMX1),$
          LEVELS=BRV/BFAC,COLOR=064,XSTYLE=5,YSTYLE=5,C_LINESTYLE=NLS*(BRV LT 0) $
ELSE $
CONTOUR,BR(IPMN1:IPMX1,ITMN1:ITMX1,NRAD),LON(IPMN1:IPMX1),LAT(ITMN1:ITMX1),$
                    LEVELS=BRV/BFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=5,YSTYLE=5


VELOVECT,BP(IAP,IAT,NRAD),-BT(IAP,IAT,NRAD),$
         LON(IAP),LAT(IAT),$
         LENGTH=ARRB,COLOR=0
         XYOUTS,9.25/XDIM,XC(3)+0.3/YDIM,'magnetic field',SIZE=CSZ*SZF,$
         /NORMAL,ALIGNMENT=0.5,COLOR=0
XYOUTS,9.25/XDIM,XC(1)-1.2/YDIM,'longitude',SIZE=CSS*SZF,$
         /NORMAL,ALIGNMENT=0.5,COLOR=0
XYOUTS,0.60/XDIM,(XC(1)+XC(3))/2.,'latitude',SIZE=CSS*SZF,$
         /NORMAL,ALIGNMENT=0.5,ORIENTATION=-90,COLOR=0

LABEL2C:  
        DEVICE, SET_FONT=SCHRIFT_GROSS
        XYOUTS,0.5,22.7/YDIM,ANNTEXT1,/NORM,SIZE=CSB*SZF,$
          ALIGNMENT=0.5,COLOR=0

        DEVICE, SET_FONT=SCHRIFT_NORMAL
        XYOUTS,0.5,22.0/YDIM,ANNTEXT2,/NORM,SIZE=CSZ*SZF,$
          ALIGNMENT=0.5,COLOR=0
        XYOUTS,0.5,21.5/YDIM,ANNTEXT3,/NORM,SIZE=CSZ*SZF,$
          ALIGNMENT=0.5,COLOR=0

        IF IGIFPS LT 1 THEN BEGIN
          VHMAX=SQRT(MAX(VP(IPMIN:IPMAX,ITMIN:ITMAX,NRAD)^2 $
                        +VT(IPMIN:IPMAX,ITMIN:ITMAX,NRAD)^2))
          BHMAX=SQRT(MAX(BP(IPMIN:IPMAX,ITMIN:ITMAX,NRAD)^2 $
                        +BT(IPMIN:IPMAX,ITMIN:ITMAX,NRAD)^2))
          PRINT,FORMAT='("MAX VH=",F8.3,"  BH=",F8.3)',VHMAX,BHMAX
        ENDIF
        GOTO, LABEL0

; ********************************************************************************
; **************** PLOT EQUATORIAL SECTIONS *********************************
; ********************************************************************************

LABEL3: PRINT, FORMAT='($, "velocity arrows=0, field arrows=1, TW=2, B components=3, B^2=4, B^2 (log)=5")'
        READ, IOPT1
        IPAGE=3
        TEQ(*,*)=(T(*,NT/2,*)+T(*,(NT/2)-1,*))/2  ;FORM EQUATORIAL TEMPS

;        dVPdzMER(*,*)=(dVPdz(*,NT/2,*)+dVPdz(*,(NT/2)-1,*))/2  ;FORM EQUATORIAL TEMPS
;        dTdtMER(*,*)=(dTdt(*,NT/2,*)+dTdt(*,(NT/2)-1,*))/2  ;FORM EQUATORIAL TEMPS
        dTdtMER=TMER
;        dVPdzMER2=(dVPdzMER(0:NPFULL/2-1,*)+dVPdzMER(NPFULL/2:NPFULL-1,*))/2
;        dTdtMER2=(dTdtMER(0:NPFULL/2-1,*)+dTdtMER(NPFULL/2:NPFULL-1,*))/2
        TMER2=(TMER(0:NPFULL/2-1,*)+TMER(NPFULL/2:NPFULL-1,*))/2
        dTdtMER2=TMER2


        ERASE
        IF IOPT1 NE 2 THEN GOTO, LABEL3A

        !P.POSITION=XY1   &  IFR=10  & COVER_IC=1
        POLAR_CONTOUR,dVPdzMER,PHIEQ(*,0:NR-1),REQ(*,0:NR-1), $
          LEVELS=VV/VFAC,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,$
          XSTYLE=4,YSTYLE=4
        
        XYOUTS,XT1,YT1,'Meridional average of dVPdz',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
          /NORMAL,COLOR=BWCOL
        GOTO, LABELCOV

L15EQ:  !P.POSITION=XY2   &  IFR=11  & COVER_IC=1
        POLAR_CONTOUR,dTdtMER,PHIEQ(*,0:NR-1),REQ(*,0:NR-1), $
          LEVELS=TV/TFAC,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,$
          XSTYLE=4,YSTYLE=4
;         print,TV
;         print,(TV+0.5)/1000
    
        XYOUTS,XT2,YT2,'Meridional average of dTdt',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
          /NORMAL,COLOR=BWCOL
        GOTO, LABELCOV

L15EQ1: !P.POSITION=XS3    &  ISF=72
        IF LCOLOR LT 1 THEN POLAR_CONTOUR,dVPdzMER2,THPLT(*,0:NR-1),RZON(*,0:NR-1), COLOR=0,$
          LEVELS=VV/HFAC,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0) $
        ELSE $
          POLAR_CONTOUR,dVPdzMER2,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
          LEVELS=VV/VFAC,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
        COVER_IC=1
        GOTO, LABELSCOV
        
L15EQ2: !P.POSITION=XS4     &  ISF=44
        POLAR_CONTOUR,dTdtMER2(*,0:NR-1),THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=TV/TFAC,$
          C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
        
        COVER_IC=1
        GOTO, LABELSCOV


LABEL3A:  ERASE

!P.POSITION=XY1   &  IFR=1  & COVER_IC=1
        POLAR_CONTOUR,TEQ,PHIEQ(*,0:NR-1),REQ(*,0:NR-1), $
                      LEVELS=TV+0.5,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,$
                      XSTYLE=4,YSTYLE=4

XYOUTS,XT1,YT1,'temperature',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
         /NORMAL,COLOR=BWCOL
        GOTO, LABELCOV

L2:    
        IF IOPT1 EQ 3 THEN BEGIN
            !P.POSITION=XY2  &  IFR=2  & COVER_IC=1
            BZE(*,*)=-(BR2(*,NT/2,*)+BR2(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL B_R
            BRL=( (FINDGEN(NLV)-1)*(1*(MAX(BZE)-MIN(BZE))/(NLV-2)) ) - abs(MIN(BZE))
            BRL(0)=-100.*MAX(abs(BZE))
            print,'levels BR: '
            print,BRL
            IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BRL/BFAC,$
              COLOR=0, XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(WV LT 0) $
            ELSE $
              POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BRL/BFAC,C_COLORS=VCOLORS,$
              /FILL,XSTYLE=4,YSTYLE=4
            
            XYOUTS,XT2,YT2,'Br',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
              /NORMAL,COLOR=BWCOL
        ENDIF ELSE IF IOPT1 EQ 4 THEN BEGIN
            !P.POSITION=XY2 &  IFR=2
            IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
;            BX(*,*,*)=BR2*SIN(pi/2-TH3D)*COS(PHI3D)+BT2*COS(pi/2-TH3D)*COS(PHI3D)-BP2*SIN(PHI3D)
;            BY(*,*,*)=BR2*SIN(pi/2-TH3D)*SIN(PHI3D)+BT2*COS(pi/2-TH3D)*SIN(PHI3D)+BP2*COS(PHI3D)
;            BZ(*,*,*)=BR2*COS(pi/2-TH3D)-BP2*SIN(TH3D)
;;            BX(*,*,*)=BR*SIN(pi/2-TH3D)*COS(pi/2-PHI3D)+BT*COS(pi/2-TH3D)*COS(pi/2-PHI3D)-BP*SIN(pi/2-PHI3D)
;;            BY(*,*,*)=BR*SIN(pi/2-TH3D)*SIN(pi/2-PHI3D)+BT*COS(pi/2-TH3D)*SIN(pi/2-PHI3D)+BP*COS(pi/2-PHI3D)
;;            BZ(*,*,*)=BR*COS(pi/2-TH3D)-BP*SIN(pi/2-TH3D)
            BS3D(*,*,*)=SQRT(BR*BR*BR*BR+BT*BT*BT*BT+BP*BP*BP*BP);BX*BX+BY*BY+BZ*BZ;SQRT(BR*BR+BT*BT+BP*BP)
;            BZE(*,*)=(BS3D(*,NT/2,*)+BS3D(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL Babs
            BZE(*,*)=total(BS3D(*,*,*),2)/NT
;;            BZE(*,*)=alog(BZE(*,*))
;;            BZE(WHERE(BZE EQ '-Inf'))=0.0
            print,'min = ',min(BZE)
            print,'max = ',max(BZE)
            BS3=( (FINDGEN(NLV)-1)*(1*(MAX(BZE)-MIN(BZE))/(NLV-2)) ) - abs(MIN(BZE))
            BS3(0)=-100.*MAX(abs(BZE))
            print,'levels B^2: '
            print,BS3
        IF LCOLOR LT 1 THEN   POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BS3/BFAC,$
              COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BRV LT 0)  $
           ELSE $
              POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BS3/BFAC,C_COLORS=VCOLORS,/FILL,$
              XSTYLE=4,YSTYLE=4

            XYOUTS,XT2,YT2,'B^2',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
              /NORMAL,COLOR=BWCOL
        ENDIF ELSE IF IOPT1 EQ 5 THEN BEGIN
            !P.POSITION=XY2 &  IFR=2
            IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
            BX(*,*,*)=BR2*SIN(pi/2-TH3D)*COS(PHI3D)+BT2*COS(pi/2-TH3D)*COS(PHI3D)-BP2*SIN(PHI3D)
            BY(*,*,*)=BR2*SIN(pi/2-TH3D)*SIN(PHI3D)+BT2*COS(pi/2-TH3D)*SIN(PHI3D)+BP2*COS(PHI3D)
            BZ(*,*,*)=BR2*COS(pi/2-TH3D)-BP2*SIN(TH3D)
            BS3D(*,*,*)=BX*BX+BY*BY+BZ*BZ
            BZE(*,*)=(BS3D(*,NT/2,*)+BS3D(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL Babs
            BZE(*,*)=alog(BZE(*,*))
            BZE(WHERE(BZE EQ '-Inf'))=0.0
            print,'min = ',min(BZE)
            print,'max = ',max(BZE)
            BS3=( (FINDGEN(NLV)-1)*(1*(MAX(BZE)-MIN(BZE))/(NLV-2)) ) - abs(MIN(BZE))
            BS3(0)=-100.*MAX(abs(BZE))
            print,'levels B^2: '
            print,BS3
        IF LCOLOR LT 1 THEN   POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BS3/BFAC,$
              COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BRV LT 0)  $
           ELSE $
              POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BS3/BFAC,C_COLORS=VCOLORS,/FILL,$
              XSTYLE=4,YSTYLE=4

            XYOUTS,XT2,YT2,'B^2 (log)',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
              /NORMAL,COLOR=BWCOL
        ENDIF ELSE BEGIN
            !P.POSITION=XY2 &  IFR=2
            IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
            BZE(*,*)=-(BT(*,NT/2,*)+BT(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL B_Z
            IF LCOLOR LT 1 THEN   POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BRV/BFAC,$
              COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BRV LT 0)  $
            ELSE $
              POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BRV/BFAC,C_COLORS=VCOLORS,/FILL,$
              XSTYLE=4,YSTYLE=4

            XYOUTS,XT2,YT2,'z-field',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
              /NORMAL,COLOR=BWCOL
        ENDELSE
        GOTO, LABELCOV
        

L3:
        IF IOPT1 EQ 3 THEN BEGIN
            !P.POSITION=XY3  &  IFR=3  & COVER_IC=1
            BZE(*,*)=-(BT2(*,NT/2,*)+BT2(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL B_T
            BTL=( (FINDGEN(NLV)-1)*(1*(MAX(BZE)-MIN(BZE))/(NLV-2)) ) - abs(MIN(BZE))
            BTL(0)=-100.*MAX(abs(BZE))
            print,'levels BT: '
            print,BTL
           IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BTL/BFAC,$
              COLOR=0, XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(WV LT 0) $
            ELSE $
              POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BTL/BFAC,C_COLORS=VCOLORS,$
              /FILL,XSTYLE=4,YSTYLE=4
            
            XYOUTS,XT3,YT3,'Btheta',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
              /NORMAL,COLOR=BWCOL
            GOTO, LABELCOV
        ENDIF ELSE BEGIN
            !P.POSITION=XY3   &  IFR=3  & COVER_IC=0
            IF ((IOPT1 LT 1) OR (IOPT1 EQ 4) OR (IOPT1 EQ 5)) THEN BEGIN
;FORM EQUATORIAL VELOCITY PLOT
                VXEQ=0.5*(VR(*,NT/2-1,*)+VR(*,NT/2,*))*COS(PHIEQ(*,0:NR-1))$
                  -0.5*(VP(*,NT/2-1,*)+VP(*,NT/2,*))*SIN(PHIEQ(*,0:NR-1))
                VYEQ=0.5*(VR(*,NT/2-1,*)+VR(*,NT/2,*))*SIN(PHIEQ(*,0:NR-1))$
                  +0.5*(VP(*,NT/2-1,*)+VP(*,NT/2,*))*COS(PHIEQ(*,0:NR-1))
                PTITLE='equatorial velocity'
                ARRLEN=ARRV
                UX=POLAR_SURFACE(VXEQ,REQ(*,0:NR-1),PHIEQ(*,0:NR-1),SPACING=ARRST,BOUNDS=[-1,-1,1,1])
                UY=POLAR_SURFACE(VYEQ,REQ(*,0:NR-1),PHIEQ(*,0:NR-1),SPACING=ARRST,BOUNDS=[-1,-1,1,1])
                                ; Set values in inner core to zero
                RINCORE2=R(NR-1)*R(NR-1)
                FOR IX=0,28 DO BEGIN
                    X2=(IX-14)*(IX-14)/196.
                    FOR IY=0,32 DO BEGIN
                        Y2=(IY-14)*(IY-14)/196.
                        IF X2+Y2 LT RINCORE2 THEN BEGIN
                            UX(IX,IY)=0.0
                            UY(IX,IY)=0.0
                        END
                    END
                END  
            END $
            ELSE BEGIN
;FORM EQUATORIAL FIELD PLOT
                IF IOPT1 eq 1 THEN BEGIN
                VXEQ=0.5*(BR(*,NT/2-1,*)-BR(*,NT/2,*))*COS(PHIEQ)$
                  -0.5*(BP(*,NT/2-1,*)-BP(*,NT/2,*))*SIN(PHIEQ)/REQ
                VYEQ=0.5*(BR(*,NT/2-1,*)-BR(*,NT/2,*))*SIN(PHIEQ)$
                  +0.5*(BP(*,NT/2-1,*)-BP(*,NT/2,*))*COS(PHIEQ)/REQ
                PTITLE='dB.h/dz'
                ARRLEN=ARRB
                UX=POLAR_SURFACE(VXEQ,REQ,PHIEQ,SPACING=ARRST,BOUNDS=[-1,-1,1,1])
                UY=POLAR_SURFACE(VYEQ,REQ,PHIEQ,SPACING=ARRST,BOUNDS=[-1,-1,1,1])
            endif
        END
    ENDELSE
        
        IF IOPT1 ne 3 THEN BEGIN
            IF IGIFPS LT 1 THEN BEGIN
                VHMAX=SQRT(MAX(UX^2+UY^2))
                IF IOPT1 EQ 1 THEN VHMAX=VHMAX*RSCALE/(THETA(NT/2)-THETA(NT/2-1))
                PRINT,FORMAT='("max VH or BH/DZ= ",F8.3)',VHMAX
            ENDIF
            VELOVECT,UX,UY,COLOR=BWCOL,LENGTH=ARRLEN,XSTYLE=5,YSTYLE=5

            XYOUTS,XT3,YT3,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL
            PLOT,/POLAR,RFAC*REQ(*,0),PHIEQ(*,0),COLOR=BWCOL,/NOERASE,$
              XRANGE=[-1,1],YRANGE=[-1,1],XSTYLE=5,YSTYLE=5
            OPLOT,/POLAR,RFAC*REQ(*,NR-1),PHIEQ(*,NR-1),COLOR=BWCOL
        ENDIF
        

L4:
    IF IOPT1 EQ 3 THEN BEGIN
        !P.POSITION=XY4  &  IFR=4  & COVER_IC=1
        BZE(*,*)=-(BP2(*,NT/2,*)+BP2(*,(NT/2)-1,*))/2 ;FORM EQUATORIAL B_P
        BPL=( (FINDGEN(NLV)-1)*(1*(MAX(BZE)-MIN(BZE))/(NLV-2)) ) - abs(MIN(BZE))
        BPL(0)=-100.*MAX(abs(BZE))
        print,'levels BP: '
        print,BPL
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BPL/BFAC,$
          COLOR=0, XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(WV LT 0) $
        ELSE $
          POLAR_CONTOUR,BZE,PHIEQ,REQ,LEVELS=BPL/BFAC,C_COLORS=VCOLORS,$
          /FILL,XSTYLE=4,YSTYLE=4
        
        XYOUTS,XT4,YT4,'Bphi',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
          /NORMAL,COLOR=BWCOL
    ENDIF ELSE BEGIN
        !P.POSITION=XY4  &  IFR=4  & COVER_IC=1
        WZE(*,*)=-(WT(*,NT/2,*)+WT(*,(NT/2)-1,*))/2 ;EQUATORIAL W_Z
        WZE=SMOOTH(WZE,3)
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,WZE,PHIEQ(*,0:NR-1),REQ(*,0:NR-1),LEVELS=WV/HFAC,$
          COLOR=0, XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(WV LT 0) $
        ELSE $
          POLAR_CONTOUR,WZE,PHIEQ(*,0:NR-1),REQ(*,0:NR-1),LEVELS=WV/HFAC,C_COLORS=VCOLORS,$
          /FILL,XSTYLE=4,YSTYLE=4

        XYOUTS,XT4,YT4,'z-vorticity',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
          /NORMAL,COLOR=BWCOL
    ENDELSE
;;    BMagR=TOTAL(TOTAL(BS3D,1),2)/(NPFULL*NT)
;    BMagR=TOTAL(BZE(*,*),1)/NPFULL
;    window,1 & plot,rzon(0,*),BMagR, color = 45, background = 200


; PLOT INNER AND OUTER CIRCLE AND COVER INNER CORE

LABELCOV:    OPLOT,/POLAR,REQ(*,0),PHIEQ(*,0),COLOR=BWCOL ; DRAW OUTER EQUATOR
        IF COVER_IC GT 0 THEN BEGIN
          POLYFILL,XIC,YIC,COLOR=!P.BACKGROUND & POLYFILL,-XIC,YIC,COLOR=!P.BACKGROUND
          END
        OPLOT,/POLAR,REQ(*,NR-1),PHIEQ(*,NR-1),COLOR=BWCOL ;DRAW INNER EQUATOR
       OPLOT,/POLAR,VCOND*REQ(*,0),PHIEQ(*,0),COLOR=BWCOL,LINESTYLE=2 ;DRAW STABLE LAYER
        CASE IFR OF
        1: GOTO, L2
        2: GOTO, L3
        3: GOTO, L4
        10: GOTO, L15EQ
        11: GOTO, L15EQ1
        12: GOTO, L15EQ2
        4: GOTO, LABELTXT
    ENDCASE


;*****************************************************************************
;******************** PLOT THERMAL WIND COMPARISON ***************************
;*****************************************************************************

LABEL15:  IPAGE=15

LABEL15A: ERASE
          
          !P.POSITION=XS1     &  ISF=70
          POLAR_CONTOUR,dTdtZON(*,0:NR-1),THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=TV/TFAC,$
            C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
          
          XYOUTS,XTS1,YTS1,'dT/dtheta',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
          COVER_IC=1
          GOTO, LABELSCOV

L15:      !P.POSITION=XS2     &  ISF=44
          IF ISINGLE GT 0 THEN !P.POSITION=XS0
          IF LCOLOR LT 1 THEN POLAR_CONTOUR,dVPdzZON,THPLT(*,0:NR-1),RZON(*,0:NR-1), COLOR=0,$
            LEVELS=VV/HFAC,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0) $
          ELSE $
            POLAR_CONTOUR,dVPdzZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
            LEVELS=VV/VFAC,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
          IF ISINGLE LT 1 THEN $
            XYOUTS,XTS2,YTS2,'dVP/dz',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL
          COVER_IC=1
          GOTO, LABELSCOV


          GOTO, LABEL0


;*****************************************************************************
;******************** PLOT MERIDIONAL SLICES ******************************
;*****************************************************************************

LABEL4:	PRINT, FORMAT='($, "slice longitude in degrees")'
	READ,PHIMS
        PHIMS=PI*PHIMS/180
        IPMS=FIX(PHIMS/PI2NP)
        PRINT, FORMAT='("angle=", F10.5, "°")', IPMS*360./NPFULL
        IPMS=IPMS+NPFULL/2 & IF IPMS GE NPFULL THEN IPMS=IPMS-NPFULL
        PRINT, FORMAT='($, "merid. v=0 / B=1 / TW=2, z-vorticity=0 / v_phi=1")'
        READ, IOPT1,IOPT2
        IF IOPT1 EQ 2 THEN GOTO, J4
        PRINT, FORMAT='($, "T=0 / B_phi=1; helicity=0 / ang. velocity=1")'
        READ, IOPT3,IOPT4
J4:     IPAGE=4

LONMS=180*PHIMS/PI - 180

;PLOT MERIDIONAL SLICES

LABEL4A:   ERASE


L40:    IF IOPT1 NE 2 THEN GOTO, L41
          !P.POSITION=XS1     &  ISF=40
          POLAR_CONTOUR,dTdt(IPMS,*,0:NR-1),THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=TV/TFAC,$
            C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
          
          XYOUTS,XTS1,YTS1,'dT/dtheta',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
          COVER_IC=1
          GOTO, LABELSCOV

L415:     !P.POSITION=XS2     &  ISF=44
          IF ISINGLE GT 0 THEN !P.POSITION=XS0
          IF LCOLOR LT 1 THEN POLAR_CONTOUR,dVPdzZON,THPLT(*,0:NR-1),RZON(*,0:NR-1), COLOR=0,$
            LEVELS=VV/HFAC,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0) $
          ELSE $
            POLAR_CONTOUR,dVPdz(IPMS,*,*),THPLT(*,0:NR-1),RZON(*,0:NR-1),$
            LEVELS=VV/VFAC,C_COLORS=LCINV*VCOLORS,FILL=LCOLOR,XSTYLE=4,YSTYLE=4
          IF ISINGLE LT 1 THEN $
            XYOUTS,XTS2,YTS2,'dVP/dtheta',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL
          COVER_IC=1
          GOTO, LABELSCOV


L41:    IF IOPT1 LT 0 THEN GOTO, L42
      !P.POSITION=XS1    &  ISF=41
        IF IOPT1 EQ 0 THEN BEGIN
          VXM=VR(IPMS,*,*)*SIN(THZON(*,0:NR-1))+VT(IPMS,*,*)*COS(THZON(*,0:NR-1))
          VZM=VR(IPMS,*,*)*COS(THZON(*,0:NR-1))-VT(IPMS,*,*)*SIN(THZON(*,0:NR-1))
          PTITLE='meridional velocity'
          COVER_IC=1  ; set COVER_IC to 1 so that the inner core is covered in LABELSCOV 
          ARRLEN=ARRV
          UX=POLAR_SURFACE(VXM,RZON(*,0:NR-1),THPLT(*,0:NR-1),SPACING=ARRST)
          UZ=POLAR_SURFACE(VZM,RZON(*,0:NR-1),THPLT(*,0:NR-1),SPACING=ARRST)
          END $
        ELSE BEGIN
          VXM=BR(IPMS,*,*)*SIN(THZON)+BT(IPMS,*,*)*COS(THZON)
          VZM=BR(IPMS,*,*)*COS(THZON)-BT(IPMS,*,*)*SIN(THZON)
          PTITLE='meridional field'
          COVER_IC=0
          ARRLEN=ARRB
          UX=POLAR_SURFACE(VXM,RZON,THPLT,SPACING=ARRST)
          UZ=POLAR_SURFACE(VZM,RZON,THPLT,SPACING=ARRST)
        END

          VM2=VXM*VXM+VZM*VZM
          VMAX=SQRT(MAX(VM2))



        IF IGIFPS LT 1 THEN BEGIN
          VHMAX=SQRT(MAX(UX^2+UZ^2))
          PRINT,FORMAT='("MAX V OR B= ",F8.3)',VHMAX
        ENDIF

          VELOVECT,UZ,UX,COLOR=BWCOL,LENGTH=ARRLEN,XSTYLE=5,YSTYLE=5

          XYOUTS,XTS1,YTS1,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL

L42: IF IOPT2 LT 0 THEN GOTO, L43
     !P.POSITION=XS2    &  ISF=42

        IF IOPT2 GT 0 THEN BEGIN
          WZON=VP(IPMS,*,*)/XZON  & VLV=3*VV/VFAC & PTITLE='angular velocity'
        END ELSE BEGIN
           WZON=WR(IPMS,*,*)*SIN(THZON(*,0:NR-1))+WT(IPMS,*,*)*COS(THZON(*,0:NR-1))
           VLV=WV/HFAC    & PTITLE='z-vorticity'
        END
           IF LCOLOR LT 1 THEN  POLAR_CONTOUR,WZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
              LEVELS=VLV,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VLV LT 0) $
           ELSE $
           POLAR_CONTOUR,WZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                      LEVELS=VLV,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
          IF IOPT1 GE 0 THEN $
           XYOUTS,XTS2,YTS2,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        COVER_IC=1
        GOTO, LABELSCOV

L43:  IF IOPT3 LT 0 THEN GOTO, L44
      !P.POSITION=XS3      &  ISF=43
        IF IOPT3 LT 1 THEN BEGIN
         IF LCOLOR LT 1 THEN $
          POLAR_CONTOUR,T(IPMS,*,*),THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                      LEVELS=TV+0.5,COLOR=0,$
                      XSTYLE=4,YSTYLE=4  $
         ELSE $
          POLAR_CONTOUR,T(IPMS,*,*),THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                      LEVELS=TV+0.5,C_COLORS=LCINV*VCOLORS,$
                      /FILL, XSTYLE=4,YSTYLE=4
          XYOUTS,XTS3,YTS3,'temperature',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
          COVER_IC=1
        END $
        ELSE BEGIN
          IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BP(IPMS,*,*),THPLT,RZON,$
             LEVELS=BPV/BFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BPV LT 0)  $
          ELSE $
          POLAR_CONTOUR,BP(IPMS,*,*),THPLT,RZON,$
                      LEVELS=BPV/BFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
          XYOUTS,XTS3,YTS3,'azimuthal field',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
            IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
        END
        GOTO, LABELSCOV


L44: IF IOPT4 LT 0 THEN GOTO, LABELTXT
     !P.POSITION=XS4   &  ISF=44
      IF IOPT4 LT 1 THEN BEGIN
        WZON=H(IPMS,*,*) & VLV=HV/HFAC & PTITLE='helicity'
      END ELSE BEGIN
        WZON=VP(IPMS,*,*) & VLV=VV/VFAC & PTITLE='azimuthal velocity'
      END
      IF LCOLOR LT 1 THEN  POLAR_CONTOUR,WZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
         LEVELS=VLV,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VLV LT 0)  $
      ELSE $
      POLAR_CONTOUR,WZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                  LEVELS=VLV,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        XYOUTS,XTS4,YTS4,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
           /NORMAL,COLOR=BWCOL
        COVER_IC=1
        GOTO, LABELSCOV


;**********************************************************************
;************************ PLOT ZONAL AVERAGES *************************
;**********************************************************************

LABEL5: PRINT, FORMAT='($, "angular velocity=0, poloidal field=1, log rms B^2=2, rms B^2=3")'
        READ, IOPT4
        IPAGE=5
LABEL5A: ERASE

       IF ISINGLE GT 0 THEN GOTO, L52

L51:  !P.POSITION=XS1     &  ISF=51
        POLAR_CONTOUR,TZON(*,0:NR-1),THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=TV+0.5,C_COLORS=LCINV*VCOLORS,$
                      FILL=LCOLOR,XSTYLE=4,YSTYLE=4

        XYOUTS,XTS1,YTS1,'temperature',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        COVER_IC=1
        GOTO, LABELSCOV

L52:  !P.POSITION=XS2     &   ISF=52
      PRINT,"ARGH"
      PRINT,"ISINGLE",ISINGLE
      PRINT,"LCOLOR",LCOLOR
      PRINT,"ARRV",ARRV
      IF ISINGLE GT 0 THEN !P.POSITION=XS0
      IF LCOLOR LT 1 THEN POLAR_CONTOUR,VZON,THPLT(*,0:NR-1),RZON(*,0:NR-1), COLOR=0,$
        LEVELS=VV/VFAC,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0) $
      ELSE $
        POLAR_CONTOUR,VZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
        LEVELS=VV/VFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
      IF ISINGLE LT 1 THEN $
        XYOUTS,XTS2,YTS2,'azimuthal flow',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
        /NORMAL,COLOR=BWCOL
      COVER_IC=1
      GOTO, LABELSCOV

L53:   IF ISINGLE GT 0 THEN GOTO, LABELTXT
      !P.POSITION=XS3   &  ISF=53
        POLAR_CONTOUR,BZON,THPLT,RZON, $
              LEVELS=BPV/BFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BZON,THPLT,RZON, $
            LEVELS=BPV/BFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BPV LT 0)
        XYOUTS,XTS3,YTS3,'azimuthal field',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
        GOTO, LABELSCOV

L54:  !P.POSITION=XS4    &  ISF=54
      IF IOPT4 EQ 0 THEN BEGIN
        POLAR_CONTOUR,ROTZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                LEVELS=3*VV/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,ROTZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
           LEVELS=3*VV/HFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0)
        XYOUTS,XTS4,YTS4,'angular velocity',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        COVER_IC=1
      END $
      ELSE IF IOPT4 EQ 2 THEN BEGIN
;        BX(*,*,*)=BR2*SIN(pi/2-TH3D)*COS(PHI3D)+BT2*COS(pi/2-TH3D)*COS(PHI3D)-BP2*SIN(PHI3D)
;        BY(*,*,*)=BR2*SIN(pi/2-TH3D)*SIN(PHI3D)+BT2*COS(pi/2-TH3D)*SIN(PHI3D)+BP2*COS(PHI3D)
;        BZ(*,*,*)=BR2*COS(pi/2-TH3D)-BP2*SIN(TH3D)
;        B2=SQRT(BX*BX*BX*BX+BY*BY*BY*BY+BZ*BZ*BZ*BZ)
        B2=SQRT(BR*BR*BR*BR+BT*BT*BT*BT+BP*BP*BP*BP)
	B2ZON=TOTAL(B2,1)/NPFULL
        B2ZON(*,0)=.0
        B2ZON(*,NR+NR_IC-1)=.0
        BS3=( (FINDGEN(NLV)-1)*(1*(MAX(B2ZON)-MIN(B2ZON))/(NLV-2)) )-abs(MIN(B2ZON))
        BS3(0)=-100.*MAX(abs(B2ZON))

        POLAR_CONTOUR,alog(B2ZON),THPLT,RZON, $
              LEVELS=BS3/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,alog(B2ZON),THPLT,RZON, $
            LEVELS=BS3/HFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BPV LT 0)
        XYOUTS,XTS4,YTS4,'zonal average of B^2',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
        GOTO, LABELSCOV
      END $
      ELSE IF IOPT4 EQ 3 THEN BEGIN
        B2=SQRT(BR*BR*BR*BR+BT*BT*BT*BT+BP*BP*BP*BP)
	B2ZON=TOTAL(B2,1)/NPFULL
        B2ZON(*,0)=.0
        B2ZON(*,NR+NR_IC-1)=.0
        BS3=( (FINDGEN(NLV)-1)*(1*(MAX(B2ZON)-MIN(B2ZON))/(NLV-2)) )-abs(MIN(B2ZON))
        BS3(0)=-100.*MAX(abs(B2ZON))

        POLAR_CONTOUR,B2ZON,THPLT,RZON, $
              LEVELS=BS3/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,B2ZON,THPLT,RZON, $
            LEVELS=BS3/HFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BPV LT 0)
        XYOUTS,XTS4,YTS4,'zonal average of B^2',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
        GOTO, LABELSCOV
      END $
      ELSE BEGIN
        POLAR_CONTOUR,A,THPLT,RZON,$
                      LEVELS=AV,C_COLORS=AACOLOR,XSTYLE=4,YSTYLE=4
        XYOUTS,XTS4,YTS4,'poloidal field lines',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
      END

        GOTO, LABELSCOV

;**********************************************************************************
;***************** SECOND SET OF ZONAL AVERAGES ****************************
;**********************************************************************************

LABEL6: PRINT, FORMAT='($, "  velocities=0, helicity=1 / helicity*B_phi=2")'
        READ, IOPT1
        PRINT, FORMAT='($, " field lines=0 / omega=1 / stream lines=2")'
        READ, IOPT3
        IPAGE=6
LABEL6A: ERASE

;TAKE CURL OF B_ZON TO GET JZON
        JZON=(SHIFT(BTZON,0,-1)-SHIFT(BTZON,0,1))/$
    	     (SHIFT(RZON,0,-1)-SHIFT(RZON,0,1))$
		+ (BTZON/RZON) $
		-(SHIFT(BRZON,-1,0)-SHIFT(BRZON,1,0))/$
       		 (RZON*(SHIFT(THZON,-1,0)-SHIFT(THZON,1,0)))
	JZON(0,*)=0 & JZON(NT-1,*)=0    ;CORRECT ENDPOINTS
;;;	JZON(*,0)=JZON(*,1) & JZON(*,NR+NR_IC-1)=JZON(*,NR+NR_IC-2)
	JZON(*,NR+NR_IC-2)=JZON(*,NR+NR_IC-3) * ( 1 - $
                           (RZON(*,NR+NR_IC-3)-RZON(*,NR+NR_IC-2)) / $
                            RZON(*,NR+NR_IC-3)  )
	JZON(*,0)=.0 & JZON(*,NR+NR_IC-1)=.0

;	FOR IP=0,NT-1 DO BEGIN
;	FOR N=0,NR+NR_IC-1 DO BEGIN
;	    PRINT, N,IP,JZON(IP,N)
;       ENDFOR
;	ENDFOR
;	PRINT, BTZON(92,30),BTZON(92,31),BTZON(92,32),BTZON(92,33)
;	PRINT, RZON(92,30),RZON(92,31),RZON(92,32),RZON(92,33)
;	PRINT, BRZON(91,31),BRZON(93,31),BRZON(91,32),BRZON(93,32)
;	PRINT, THZON(91,31),THZON(93,31),THZON(91,32),THZON(93,32)
	

;CALCULATE OMZON, THE TOROIDAL FIELD OMEGA-EFFECT SOURCE TERM
        IF IOPT3 EQ 1 THEN BEGIN

;             FOR IP=0,NPFULL-1 DO  $
;               WORK(IP,*,*)=VP(IP,*,*)/XZON      ;Angular velocity

;               OMEG=BR*(SHIFT(WORK,0,0,-1)-SHIFT(WORK,0,0,1))/$
;                  (SHIFT(R3D, 0,0,-1)-SHIFT(R3D, 0,0,1))    $
;                +BT*(SHIFT(WORK,0,1,0) -SHIFT(WORK,0,-1,0))/$
;                (R3D*(SHIFT(TH3D,0,1,0) -SHIFT(TH3D,0,-1,0)))

;             FOR IP=0,NPFULL-1 DO  $
;               OMEG(IP,*,*)=OMEG(IP,*,*)*XZON

;             OMZON=TOTAL(OMEG,1)/NPFULL


            WORK2D=VZON/XZON    ;Angular velocity
            OMZON=BRZON*(SHIFT(WORK2D,0,-1)-SHIFT(WORK2D,0,1))/$
              (SHIFT(RZON,0,-1)-SHIFT(RZON,0,1))    $
              +BTZON*(SHIFT(WORK2D,1,0) -SHIFT(WORK2D,-1,0))/$
              (RZON*(SHIFT(THZON,1,0) -SHIFT(THZON,-1,0)))
            

        	OMZON(0,*)=0 & OMZON(NT-1,*)=0    ;CORRECT ENDPOINTS
;;;     	OMZON(*,0)=OMZON(*,1) & OMZON(*,NR-1)=OMZON(*,NR-2)
        	OMZON(*,0)=0. & OMZON(*,NR-1)=.0

        ENDIF
        IF IOPT3 EQ -1 THEN BEGIN
            OMZON=BRZON*(SHIFT(ROTZON,0,-1)-SHIFT(ROTZON,0,1))/$
                  (SHIFT(RZON,0,-1)-SHIFT(RZON,0,1))           $
                 +BTZON*(SHIFT(ROTZON,1,0)-SHIFT(ROTZON,-1,0))/ $
                  (RZON*(SHIFT(THZON,1,0)-SHIFT(THZON,-1,0)))

            OMZON=OMZON*XZON

        	OMZON(0,*)=0 & OMZON(NT-1,*)=0    ;CORRECT ENDPOINTS
        	OMZON(*,0)=OMZON(*,1) & OMZON(*,NR-1)=OMZON(*,NR-2)

        ENDIF

        JMAX=MAX(ABS(JZON))+EPS
        IF NR GE 9 THEN OMMAX=MAX(ABS(OMZON(*,4:NR-5)))+EPS ELSE OMMAX=MAX(ABS(OMZON))+EPS
        JSTEP=2*JMAX/(NLV-2)
        OMSTEP=2*OMMAX/(NLV-2)
;ASSIGN CONTOURING LEVELS
        JV=(FINDGEN(NLV)-1)*JSTEP-JMAX
        JV(0)=-100.*JMAX
        OMV=(FINDGEN(NLV)-1)*OMSTEP-OMMAX
        OMV(0)=-100.*OMMAX

L61:   !P.POSITION=XS1  &  ISF=61
        IF IOPT1 LE 0 THEN BEGIN
         IF LCOLOR LT 1 THEN BEGIN
            POLAR_CONTOUR,VZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
            LEVELS=VV/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
            POLAR_CONTOUR,VZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),COLOR=0,/OVERPLOT,$
            LEVELS=VV/HFAC,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(VV LT 0)
         END $
         ELSE $
           POLAR_CONTOUR,VZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
            LEVELS=VV/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
         IF ARRV GT 0 AND LCOLOR GT 0 THEN $
           POLAR_CONTOUR,PSI,THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=PV*ARRV, $
              COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NSTYLE*(PV LT 0)

          IF IOPT1 GE 0 THEN $
          XYOUTS,XTS1,YTS1,'velocity',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
          COVER_IC=1
        END

        IF IOPT1 EQ 1 THEN BEGIN
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,HZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
             LEVELS=HV/HFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(HV LT 0)  $
        ELSE $
          POLAR_CONTOUR,HZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
             LEVELS=HV/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        XYOUTS,XTS1,YTS1,'zonal helicity',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
           COVER_IC=1
        END

        IF IOPT1 EQ 2 THEN BEGIN
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,HBZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
              LEVELS=HBV/HFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(HBV LT 0)  $
        ELSE $
          POLAR_CONTOUR,HBZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
                LEVELS=HBV/HFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        IF IOPT1 GE 0 THEN $
        XYOUTS,XTS1,YTS1,'-helicity x Bphi',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
           COVER_IC=1
        END
        GOTO, LABELSCOV


L62:  IF IOPT1 LT 0 THEN GOTO, LABELTXT
      !P.POSITION=XS2   &  ISF=62
     IF LCOLOR LT 1 THEN  POLAR_CONTOUR,BZON,THPLT,RZON,$
       LEVELS=BPV/BFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(BPV LT 0)  $
     ELSE $
       POLAR_CONTOUR,BZON,THPLT,RZON,$
        LEVELS=BPV/BFAC,C_COLORS=VCOLORS,XSTYLE=4,YSTYLE=4,/FILL
        XYOUTS,XTS2,YTS2,'azimuthal field',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
         IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
     GOTO,LABELSCOV

L63: !P.POSITION=XS3  &  ISF=63
        IF LCOLOR LT 1 THEN POLAR_CONTOUR,JZON,THPLT,RZON,$
              LEVELS=JV/TFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(JV LT 0) $
        ELSE BEGIN
          POLAR_CONTOUR,JZON,THPLT,RZON,$
                 LEVELS=JV/TFAC,C_COLOR=VCOLORS,XSTYLE=4,YSTYLE=4,/FILL
          POLAR_CONTOUR,A,THPLT,RZON,LEVELS=AV,COLOR=0,XSTYLE=4,YSTYLE=4
        END
        XYOUTS,XTS3,YTS3,'toroidal current',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
           IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
         GOTO, LABELSCOV



L64:   !P.POSITION=XS4  &  ISF=64
        IF IOPT3 EQ 0 THEN BEGIN
        POLAR_CONTOUR,A,THPLT,RZON,LEVELS=AV,C_COLORS=AACOLOR,XSTYLE=4,YSTYLE=4
        PTITLE='poloidal fieldlines'
        IF NN_IC GT 0 THEN COVER_IC=0 ELSE COVER_IC=1
        END
        IF IOPT3 EQ 1 OR IOPT3 EQ -1 THEN BEGIN
        IF LCOLOR LT 1 THEN  POLAR_CONTOUR,OMZON,THPLT,RZON,$
         LEVELS=OMV/VFAC,COLOR=0,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(OMV LT 0) $
        ELSE $
         POLAR_CONTOUR,OMZON,THPLT(*,0:NR-1),RZON(*,0:NR-1),$
          LEVELS=OMV/VFAC,C_COLORS=VCOLORS,/FILL,XSTYLE=4,YSTYLE=4
        PTITLE='mean omega'
        COVER_IC=1
        END
        IF IOPT3 GT 1 THEN BEGIN
        POLAR_CONTOUR,PSI,THPLT(*,0:NR-1),RZON(*,0:NR-1),LEVELS=PV, $
           C_COLORS=LCINV*VCOLORS,XSTYLE=4,YSTYLE=4,C_LINESTYLE=NLS*(PV LT LIMIT)
         PTITLE='meridional streamlines'
        COVER_IC=1
        END

        XYOUTS,XTS4,YTS4,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,$
            /NORMAL,COLOR=BWCOL
        GOTO, LABELSCOV


;*******************************************************************************
;*******  COVER UP INNER CORE AND DRAW CIRCLES FOR SLICE PLOTS  ****************
;*******************************************************************************

LABELSCOV:  
        IF COVER_IC EQ 1  THEN POLYFILL,-XIC,YIC,COLOR=!P.BACKGROUND
        OPLOT,/POLAR,RZON(*,0),THPLT(*,0),COLOR=BWCOL ; DRAW OUTER EQUATOR
        OPLOT,/POLAR,VCOND*RZON(*,0),THPLT(*,0),COLOR=BWCOL,LINESTYLE=2 ;DRAW STABLE LAYER
        OPLOT,/POLAR,RZON(*,NR-1),THPLT(*,NR-1),COLOR=BWCOL ;DRAW INNER EQUATOR
        OPLOT,/POLAR,RZON(0,*),THPLT(0,*),COLOR=BWCOL
        OPLOT,/POLAR,RZON(NT-1,*),THPLT(NT-1,*),COLOR=BWCOL
        CASE ISF OF
        40: GOTO, L415
        41: GOTO, L42
        42: GOTO, L43
        43: GOTO, L44
        44: GOTO, LABELTXT
        51: GOTO, L52
        52: GOTO, L53
        53: GOTO, L54
        54: GOTO, LABELTXT
        61: GOTO, L62
        62: GOTO, L63
        63: GOTO, L64
        64: GOTO, LABELTXT
        70: GOTO, L15
        71: GOTO, L15EQ1
        72: GOTO, L15EQ2
        ENDCASE


;*******************************************************************************
;***************** PLOTTING OF GLOBAL MAP PROJECTION  **************************
;*******************************************************************************

LABEL7: PRINT, FORMAT='($, "B_r at CMB: inclination, longitude?")'
        READ,INCL,LONSHIFT
        IPAGE=7

LABEL7A:  ERASE

;MAP RADIAL FIELD
!P.POSITION=XSS
        WRAP(0:NPFULL-1,*)=BR(*,*,0) & WRAP(NPFULL,*)=BR(0,*,0)
        PTITLE='radial field'
	MAP_SET, INCL,LONSHIFT,0,/ORTHOGRAPHIC,/ADVANCE,/NOBORDER
	CONTOUR,WRAP,LONWRAP,LAT,LEVELS=BRV/BFAC,/OVERPLOT,$
         /CELL_FILL,$
         C_COLORS=LCINV*VCOLORS,C_LINESTYLE=NLS*(BRV LT LIMIT)
;XYOUTS,XT3,YT3,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL
	MAP_GRID,COLOR=BWCOL,GLINETHICK=1

        GOTO, LABELTXT


;*******************************************************************************
;*************** PLOT IN AITOFF PROJECTION *************************************
;*******************************************************************************

LABEL8: 
        IPAGE=8
        ERASE

LABEL81:
        PRINT, FORMAT='($, "(B=0  VR=1 VP=2 VT=3), radial level, longitude? [upper plot]")'
        READ, IOPTBV1, NRAD1, LONSHIFT1
        IF NRAD1 LT 0 OR NRAD1 GT NR+NR_IC-1 $
        OR (IOPTBV1 EQ 1 AND NRAD1 GT NR-1) THEN BEGIN
          PRINT, "radial level out of range"
          GOTO, LABEL81
        ENDIF

LABEL82:
        PRINT, FORMAT='($, "(B=0  VR=1 VP=2 VT=3), radial level, longitude? [lower plot]")'
        READ, IOPTBV2, NRAD2, LONSHIFT2
        IF NRAD2 LT 0 OR NRAD2 GT NR+NR_IC-1 $
        OR (IOPTBV2 EQ 1 AND NRAD2 GT NR-1) THEN BEGIN
          PRINT, "radial level out of range"
          GOTO, LABEL82
        ENDIF

LABEL8A:
        FOR I=1,2 DO BEGIN ;=========
          IF I EQ 1 THEN BEGIN
            IOPTBV=IOPTBV1
            NRAD=NRAD1
            LONSHIFT=LONSHIFT1
            !P.POSITION=[0.8/XDIM,10.5/YDIM,17.2/XDIM,18.5/YDIM]
            YTXT=19./YDIM
          END ELSE BEGIN
            IOPTBV=IOPTBV2
            NRAD=NRAD2
            LONSHIFT=LONSHIFT2
            !P.POSITION=[0.8/XDIM, 0.5/YDIM,17.2/XDIM, 8.5/YDIM]
            YTXT= 9./YDIM
          END

	  MAP_SET, 0, LONSHIFT, 0, /AITOFF, /NOBORDER, /NOERASE        

          IF IOPTBV EQ 0 THEN BEGIN 
          ;*** MAP RADIAL FIELD
            PTITLE=STRING(FORMAT='("radial field at r=", F5.3)', R(NRAD))
            WRAP(0:NPFULL-1,*)=BR(*,*,NRAD) & WRAP(NPFULL,*)=BR(0,*,NRAD)
            IF LCOLOR GT 0 THEN $
	      CONTOUR,WRAP,LONWRAP,LAT,LEVELS=BRV/BFAC,/OVERPLOT,$
              /CELL_FILL,/NOERASE,C_COLORS=LCINV*VCOLORS 
            IF LCOLOR EQ 0 THEN $
              CONTOUR,WRAP,LONWRAP,LAT,LEVELS=BRV/BFAC,/OVERPLOT,$
              /NOERASE,COLOR=0,C_LINESTYLE=NLS*(BRV LT LIMIT) 
          END ELSE BEGIN 
          ;*** MAP RADIAL VELOCITY
              CASE IOPTBV OF
                  1: BEGIN
                      PTITLE=STRING(FORMAT='("radial velocity at r=", F5.3)', R(NRAD))
                      WRAP(0:NPFULL-1,*)=VR(*,*,NRAD) & WRAP(NPFULL,*)=VR(0,*,NRAD)
                  END
                  2: BEGIN
                      PTITLE=STRING(FORMAT='("zonal velocity at r=", F5.3)', R(NRAD))
                      WRAP(0:NPFULL-1,*)=VP(*,*,NRAD) & WRAP(NPFULL,*)=VP(0,*,NRAD)
                  END
                  3: BEGIN
                      PTITLE=STRING(FORMAT='("meridional velocity at r=", F5.3)', R(NRAD))
                      WRAP(0:NPFULL-1,*)=VT(*,*,NRAD) & WRAP(NPFULL,*)=VT(0,*,NRAD)
                  END
                  ELSE: GOTO, LABEL81
              ENDCASE
;;;            PTITLE=STRING(FORMAT='("meridional velocity at r=", F5.3)', R(NRAD))
;;            PTITLE=STRING(FORMAT='("zonal velocity at r=", F5.3)', R(NRAD))
;            PTITLE=STRING(FORMAT='("radial velocity at r=", F5.3)', R(NRAD))
;            WRAP(0:NPFULL-1,*)=VR(*,*,NRAD) & WRAP(NPFULL,*)=VR(0,*,NRAD)
;;            WRAP(0:NPFULL-1,*)=VP(*,*,NRAD) & WRAP(NPFULL,*)=VP(0,*,NRAD)
;;;            WRAP(0:NPFULL-1,*)=VT(*,*,NRAD) & WRAP(NPFULL,*)=VT(0,*,NRAD)
            IF LCOLOR GT 0 THEN $
  	      CONTOUR,WRAP,LONWRAP,LAT,LEVELS=VV/VFAC,/OVERPLOT,$
              /CELL_FILL,/NOERASE,C_COLORS=LCINV*VCOLORS 
            IF LCOLOR EQ 0 THEN $
              CONTOUR,WRAP,LONWRAP,LAT,LEVELS=VV/VFAC,/OVERPLOT,$
              /NOERASE,COLOR=0,C_LINESTYLE=NLS*(VV LT LIMIT)
          END

 	  MAP_GRID,COLOR=1,LATDEL=30,LONDEL=45,GLINETHICK=1
          XYOUTS,9.0/XDIM,YTXT,PTITLE,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL,COLOR=BWCOL

        END ;=========

        GOTO, LABELTXT

;******************************************************************************
;********** ADD TITLES ON TOP OF PAGE *****************************************
;******************************************************************************

LABELTXT: $
;        DEVICE, SET_FONT=SCHRIFT_GROSS
        XYOUTS,0.5,22.7/YDIM,ANNTEXT1,/NORM,SIZE=CSB*SZF,ALIGNMENT=0.5,COLOR=BWCOL
;        DEVICE, SET_FONT=SCHRIFT_NORMAL
        XYOUTS,0.5,22.0/YDIM,ANNTEXT2,/NORM,SIZE=CSZ*SZF,ALIGNMENT=0.5,COLOR=BWCOL
        XYOUTS,0.5,21.5/YDIM,ANNTEXT3,/NORM,SIZE=CSZ*SZF,ALIGNMENT=0.5,COLOR=BWCOL
        GOTO,LABEL0

;**************************************************************************
;*********  OPTION FOR CREATING PS-FILE OR GIF_FILE *****************
;**************************************************************************
LABELOUT:   IF IGIFPS EQ 1 THEN BEGIN
             DEVICE,/CLOSE
             SET_PLOT, 'X'
             IGIFPS=0
             SZF=1.0
             GOTO, LABEL0
            ENDIF

            IGIFPS=IOPTION-20
            IF IGIFPS LT 1 OR IGIFPS GT 2 THEN GOTO,LABEL0
            OUTFILE=''
            PRINT, FORMAT='($, "output file name?")'
            READ, OUTFILE

            IF IGIFPS EQ 2 THEN  BEGIN
              WRITE_GIF,OUTFILE,TVRD()
              IGIFPS=0
              GOTO, LABEL0
            END  $
            ELSE BEGIN
             SET_PLOT,'PS'
             IF LCOLOR GT 0 THEN $
               DEVICE,FILENAME=OUTFILE,/COLOR $
             ELSE $
               DEVICE,FILENAME=OUTFILE
             DEVICE,XSIZE=XDIM,YSIZE=YDIM,XOFFSET=1.8,YOFFSET=2.0
             SZF=.8*SCWINDOW
             !P.FONT=1
             CASE IPAGE OF
               1: GOTO,LABEL1A
               2: GOTO,LABEL2A
               3: GOTO,LABEL3A
               4: GOTO,LABEL4A
               5: GOTO,LABEL5A
               6: GOTO,LABEL6A
               7: GOTO,LABEL7A
               8: GOTO,LABEL8A
               15: GOTO,LABEL15A
               ELSE: GOTO,LABEL0
            ENDCASE
          END
            GOTO, LABEL0

;***********************************************************************

LABEL99:        !P.MULTI=0    ;RETURN TO ONE PLOT PER PAGE
;ENDS PROCEDURE MAGSYM5.PRO
	END


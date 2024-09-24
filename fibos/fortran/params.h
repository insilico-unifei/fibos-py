c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *c
c* SIMS is written by Yury N Vorobjev, Computational Structural Biology Group,*c
c* Department of Biochemistry and Biophysics,                                 *c
c* University of North Carolina at Chapel Hill, Chapel Hill, NC 27599, USA    *c
c* e-mail: vorobjev@femto.med.unc.edu                                         *c
c* Permanent adress: Novosibirsk Institute of Bioorganic Chemistry,           *c
c* 8 Lavrentjeva Ave., Novosibirsk 630090, Russia                             *c
c* Copyright 1997. All rights reserved.                                       *c
c* SIMS method description: Biophysical J. 73:722-732, (1997)                 *c
c* SIMS: computation of a Smooth Invariant Molecular Surface.                 *c
c* Yury N Vorobjev and Jan Hermans                                            *c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *c
c f77

        implicit none

        integer MAXATM
c       PARAMETER (MAXATM=4000)  !Regular-maximal number of atoms in molecule
         PARAMETER (MAXATM=12500)  !Large-maximal number of atoms in molecule
 
c parameters set up permanently: 
	integer nres_max
	PARAMETER (nres_max = MAXATM/5)
        integer MAXNBR
        PARAMETER (MAXNBR=150)    
        integer MAXSPH
        PARAMETER (MAXSPH=1000)
        integer MAXCIR
        PARAMETER (MAXCIR=1000)
        integer maxarc
        PARAMETER (maxarc=100)
        integer MAXYON
        PARAMETER (MAXYON=MAXATM/3)
        integer MAXCUB
        PARAMETER (MAXCUB=40)
	integer MAXCHIL
	parameter(MAXCHIL = 1)
	integer MAXDOT
        PARAMETER (MAXDOT = 40*MAXATM)  
	integer maxprob
	parameter (maxprob= MAXATM)   !max number for probe positions
        integer ndotccmx            !max dots on one concave face
        parameter(ndotccmx=MAXSPH/2)
        integer narcmx
        parameter(narcmx=maxarc)
        integer nrotmx
        parameter (nrotmx = 50)
        integer ndotsdmx
        parameter(ndotsdmx=narcmx*nrotmx)
        integer ndottypemx
        parameter(ndottypemx=5)   
        integer ndot_smoothmx
        parameter(ndot_smoothmx = MAXDOT/20)
	integer ndotsmp2_max
	parameter (ndotsmp2_max=MAXDOT/200)
	integer nsmp2_max
	 parameter (nsmp2_max = ndotsmp2_max/2)

         integer kanalxyz
         parameter (kanalxyz=11)
         integer kanalin
         parameter (kanalin=12)
         integer kanalp
         parameter (kanalp=6)
         integer kanalrad
         parameter (kanalrad=13)
         integer kanalq
         parameter (kanalq=14)
         integer kanalpdbout
         parameter(kanalpdbout=15)

         integer kanalx,kanals 
         parameter(kanals=22)
         parameter(kanalx=23)  
         integer kanalz
         parameter(kanalz=24)
c------------------------------------------------------------------------------
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *c
c* SIMS is written by Yury N Vorobjev, Computational Structural Biology Group,*c
c* Department of Biochemistry and Biophysics,                                 *c
c* University of North Carolina at Chapel Hill, Chapel Hill, NC 27599, USA    *c
c* e-mail: vorobjev@femto.med.unc.edu                                         *c
c* Permanent adress: Novosibirsk Institute of Bioorganic Chemistry,           *c
c* 8 Lavrentjeva Ave., Novosibirsk 630090, Russia                             *c
c* Copyright 1997. All rights reserved.                                       *c
c* SIMS method description: Biophysical J. 73:722-732, (1997)                 *c
c* SIMS: computation of a Smooth Invariant Molecular Surface.                 *c
c* Yury N Vorobjev and Jan Hermans                                            *c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *c
c f77
c
c  input_sims.h  : input data, OPT-parameters, timing
c

        real*8 time_zero,time_total,time_inc,time_last
	common/timecon/time_zero,time_total,time_inc,time_last 
        integer outdetl
        common/outstl/outdetl

	real*8 dotden_l,dotden_h,rp_rpr,rad_sm,ster_ZOOM
	common/surfprm/dotden_l,dotden_h,rp_rpr,rad_sm,ster_ZOOM

        logical OPT_dot_surface,OPT_dot_file,OPT_pdb_surface
        logical OPT_sterdot_file,OPT_dot_midas,OPT_VMDpdb_surf
	logical OPT_date, OPT_refcut
	logical OPT_pdbdotkin, OPT_dotnrmkin,OPT_pdbrext
        common/OPT_sims/OPT_dot_surface,OPT_dot_file,OPT_pdb_surface,
     &             OPT_sterdot_file,OPT_dot_midas,OPT_VMDpdb_surf,
     &             OPT_date,OPT_refcut,OPT_pdbdotkin, OPT_dotnrmkin,
     &             OPT_pdbrext
c.
c...........................................................................

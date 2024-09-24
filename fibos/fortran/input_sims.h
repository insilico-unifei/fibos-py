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

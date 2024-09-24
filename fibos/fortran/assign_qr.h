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
c parameters for assign_qr  subroutine
c
        integer nclist,nrlist,nrmax,ncmax,natmx,ncrgmx

	parameter (nclist = 1000)
	parameter (nrlist = 1000)
	parameter (nrmax = 1000)
	parameter (ncmax = 1000)
	parameter (natmx = maxatm)
	parameter (ncrgmx = 10000)
      
      dimension atnam(nrmax)		
      dimension rnam(nrmax)		
      real*8 radt(nrmax)		
      integer irlink
      dimension irlink(nrlist)	
      integer irnumb
      dimension irnumb(nrlist)	
      dimension catnam(ncmax)		
      dimension cchn(ncmax)		
      dimension crnam(ncmax)		
      dimension crnum(ncmax)		
      real*8 chrgvt(ncmax)		
      integer iclink
      dimension iclink(nclist)	
      integer icnumb
      dimension icnumb(nclist)	
      real*8 atcrd(3,natmx),atrad(natmx),atcrg(natmx)	
      dimension atmcrg(4,ncrgmx)
      real*4 atmcrg

        character*1 chn,schn,cchn
        character*3 crnam,rnam,sres,res
        character*4 rnum,snum,crnum
        character*6 atnam,catnam
        character*6 atm     
c---------------------------------------------------
	common
     &	/linkk/  irlink,irnumb,iclink,icnumb,irtot,ictot
     &	/name/  atnam,rnam,catnam,cchn,crnam,crnum
     &	/value/ radt,chrgvt
c---------------------------------------------------

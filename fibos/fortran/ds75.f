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

c:sims-main
c
        subroutine runSIMS (method)

        include "surf-sims.h"
	include "input_sims.h"
	character*3 atnamel(maxatm),resnamel(maxatm)
	integer natoml,iatnumbl(maxatm),resnumbl(maxatm)
	real*8 atmxyz(3,maxatm),atmrad(maxatm),atmcha(maxatm)
        real*8 dotcrd(3,maxdot), dotnrm(3,maxdot),dotarea(maxdot)
        integer dot_IATNUM(maxdot),dot_ISHAPE(maxdot)
	integer smdotn_at(maxatm),smdotstartn_at(maxatm)
        integer  ndots
        integer method
        real*8 d_probe,dotden
	character*8 hour
	character*9 day
	integer iday,iyear
	character*3 cmonth
        character*80 structname

        integer i,j,k
	character*12 crd_input
	logical fext
	logical CONTROL
        integer natoml_pr,natoml_het
        integer atom_type(maxatm)
        integer timeindex,iiprint0,iiprint1
        logical atom,hetatom
        logical OPT_pdb
	character*6 head,head_l(maxatm)
        character*80 liner,line1,line2,line3
        character*2 opt_ass

c default
        iiprint0=1
        OPT_pdb = .FALSE.
        opt_ass= 'R'

c       write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
c       write(kanalp,*)'* * * SMOOTH INVARIANT MOLECULAR SURFACE * * *'
c       write(kanalp,*)'* * * Yury Vorobjev, 1997:               * * *'
c       write(kanalp,*)'* * * vorobjev@femto.med.unc.edu         * * *'
c       write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
c
c      call date(day)
c      call time(hour)
c      write(kanalp,'(a12,a9,a3,a8)')'SIMS: start:',day,' : ',hour

c input parameters
c      call input_sims

       OPT_refcut = .false.
       OPT_dot_surface =.false.
       OPT_dot_file=.false.
       OPT_dot_midas = .true.
       OPT_pdb_surface = .false.
       OPT_VMDpdb_surf = .false.
       OPT_date = .false.
       OPT_pdbdotkin = .false.
       OPT_dotnrmkin = .false.
       OPT_pdbrext = .false.
       OPT_sterdot_file = .false.


c      if(OPT_date)then
clastday='1-Jan-00'
c      if(day(3:3).ne.'-') day=' '//day

c      read(day,'(i2,1x,a3,1x,i2)')iday,cmonth,iyear
c       write(kanalp,'(i2,a1,a3,a1,i2)')iday,'-',cmonth,'-',iyear

c      if(iyear.eq.99.and.
c    &  (cmonth.eq.'Oct'.or.cmonth.eq.'Nov'.or.cmonth.eq.'Dec'))then
c      write(kanalp,'(a44)')'**WARNING
c      write(kanalp,'(a60)')
c    & 'SIMS: permission will expire soon, contact Yury Vorobjev: '
c      write(kanalp,'(a36)')
c    & '      vorobjev@femto.med.unc.edu or'
c      write(kanalp,'(a36)')
c    & '      vorobjev@milli.med.unc.edu or'
c      write(kanalp,'(a36)')
c    & '      vorobjev@niboch.nsc.ru       '
c      end if

c      if(abs(iyear-99).gt.5)then
c      write(kanalp,'(a52)')
c    & 'SIMS: permission is expired, contact Yury Vorobjev: '
c      write(kanalp,'(a36)')
c    & '      vorobjev@femto.med.unc.edu or'
c      write(kanalp,'(a36)')
c    & '      vorobjev@milli.med.unc.edu or'
c      write(kanalp,'(a36)')
c    & '      vorobjev@niboch.nsc.ru       '

c      stop
c      end if
c      end if

c      if(outdetl.ge.iiprint0)CONTROL=.true.
       CONTROL=.false.
c.......................................................................
c READ in xyz + rad + q data
       crd_input='part_i.pdb'
       inquire(file=crd_input, exist=fext)
       if(.not. fext) then
       write(kanalp,*) 'part_i.pdb file does not exist '
       stop
       endif
c
       open(unit=kanalxyz,file=crd_input,form='formatted',
     &     status='old' )


c READ PDB file and assign to (local) variables
c       rewind kanalxyz
c       read(kanalxyz,'(a80)')liner
c       if(liner(1:5).eq.'TITLE')then
c       structname = liner(7:80)
c       else
c       structname = 'UNTITLED MOLECULE'
c       if(CONTROL)then
c       write(kanalp,*)
c       write(kanalp,*)'WARNING
c       write(kanalp,*) '(make title line in pdb file: TITLE:MyTitle)'
c       write(kanalp,*)
c       end if
c       end if

c       if(CONTROL)then
c       write(kanalp,'(a6,a70)')'TITLE:',structname
c       end if

           natoml_pr=0
           natoml_het=0
           k=0
c read xyz
          rewind kanalxyz

           do i=1,maxatm
           read(kanalxyz,'(a80)',end=401)liner

           head=liner(1:6)
           if(head.eq.'ATOM  '.or.head.eq.'HETATM')then
           k=k+1
           if(head.eq.'ATOM  ')then
           atom=.true.
           hetatom=.false.
           natoml_pr=natoml_pr+1
           else
           atom=.false.
           end if

           if(head.eq.'HETATM')then
           hetatom=.true.
           atom=.false.
           natoml_het=natoml_het+1
           else
           hetatom=.false.
           end if

           if(atom) atom_type(k)=1
           if(hetatom)atom_type(k)=0

           if(k.gt.maxatm)then
           write(kanalp,*)'ERROR: Too much atoms in PDB file'
           write(kanalp,*)'       Allowed maxatm = ',maxatm
           write(kanalp,*)'       increase parameter (maxatm)...'
           stop
           end if

           OPT_pdbrext = .false.

	   if(.not.OPT_pdbrext)then
           read(liner,7071)
     &     head_l(k),
     &     iatnumbl(k),atnamel(k),resnamel(k),resnumbl(k),
     &     (atmxyz(j,k),j=1,3)
	   end if

	   if(OPT_pdbrext)then
           read(liner,7081)
     &     head_l(k),
     &     iatnumbl(k),atnamel(k),resnamel(k),resnumbl(k),
     &     (atmxyz(j,k),j=1,3),atmrad(k)
	   end if

            end if
	    end do

401         natoml=k

            if(natoml.ne.(natoml_pr+natoml_het))then
            write(kanalp,*)'ERROR: BAD atom COUNT'
            stop
            end if

            if(CONTROL)then
            write(kanalp,*)
            write(kanalp,*)'FINISH reading *.pdb '
            write(kanalp,*)'Total atoms=',natoml
            write(kanalp,*)
            end if

            close(kanalxyz)

c assign Rad
	    if(.not.OPT_pdbrext)then
            call assign_rq(opt_ass,
     &                  natoml,atnamel,resnamel,atmrad,atmcha)

           if(CONTROL)then
            write(kanalp,*)'Finish assignment of Atomic Radii...'
            write(kanalp,*)
           end if
	   else
	   write(kanalp,*)'Atomic Radii are taken from PDB ...'
	   end if

           if(outdetl.ge.iiprint0)then
           write(kanalp,*)
           write(kanalp,*)
     &     'REMARK : initial PDB file with assigned AtRadii...'
           do k=1,natoml
           write(kanalp,7081)
     &     head_l(k),
     &     iatnumbl(k),atnamel(k),resnamel(k),resnumbl(k),
     &     (atmxyz(j,k),j=1,3), atmrad(k)
           end do
           write(kanalp,*)
	   end if


c do SIMS calculation:
              IF(method .eq. 1)then
                dotden_h=5.0
              END IF
              IF(method .eq. 2)then
                dotden_h=5.0
              END IF
              dotden=dotden_h
              rp_rpr=1.4
              d_probe=rp_rpr*2.0d0

              call surf_sims(atmxyz,
     &                 atnamel,resnumbl,resnamel,
     &                 atmrad,
     &                 natoml,d_probe,dotden,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,ndots,
     &                 smdotn_at,smdotstartn_at, method)


        if(outdetl.ge.0)then
        timeindex=1
c       call timit(timeindex)
        end if

c.........................................................................
      close(kanalp)

7071    format(a6,1x,i4,2x,a3,1x,a3,2x,i4,4x,3d8.3)
7081    format(a6,1x,i4,2x,a3,1x,a3,2x,i4,4x,3f8.3,f6.3)
c        stop
	end
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
c.  Double precission input/output  arguments

       subroutine surf_sims(coords,
     &                 atnamel,resnumbl,resnamel,
     &                 atom_rad,
     &                 natoms,d_probe,dotden,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,ndots,
     &                 dot_num_atom,dot_startn_atom, method)

c INPUT:
c       coords(3,*) coord of atoms
c       atnamel(*)  atom name
c       resnumbl(*) residue number of atom(*)
c       resnamel(*) residue name of atom(*)
c       atom_rad(*) atomic radii
c       natoms      number of atoms
c       d_probe     probe diameter
c       dotden      dot density per A^2
C       outsurf     0/1 - flag to write dot-file
c
c RESULT:
c       dotcrd(3,*) dot coords
c       dotnrm(3,*) normal vectors
c       dotarea(*)  dot area
c       dot_iatnum(*) atom number to which dot(*) belongs too
c       dot_ishape(*) face type (1,2,3 - contact, saddle, reentrant concave)
c                     to which dot(*) belongs too
c       ndots         number of dots on sufrace
c
c
        include'surf-sims.h'
        include'input_sims.h'

        real*8 atom_rad(*),coords(3,*),dotarea(*),
     &       dotcrd(3,*),dotden,dotnrm(3,*),d_probe

        character*3 atnamel(*),resnamel(*)
        integer resnumbl(*)
        integer method
	integer natoms,ndots
        integer*4 dot_IATNUM(maxdot),dot_ISHAPE(maxdot)

        integer*4 dot_num_atom(maxatm),dot_startn_atom(maxatm)

c
c local VARIABLES:
c
        integer id,iatnn,nnnd
    	integer ndotrem,normbad
        integer iprint0,iprint1
        integer kanalpl
	character*4 satom
	character*18 atomref
	character*6 dotsch
	real*8 badsurf_mx,badsurf,surf,surf1,surf0
	real*8 areaRen
	real*8 athresh,athreshmx(ndottypemx)
	real*8 showrad,shdist

	integer ndotconc_ijk(maxnbr)
	integer maxnbr1

	integer nprobpos,nprobpos_t
	integer probmax1,maxdot1

	    integer start_dotn(maxprob),stop_dotn(maxprob)
	    integer ndotccsd, ndotccsd_t

        logical*1 WDYON(maxdot)

    	logical*1 WDYONRm(maxdot)

        integer WDYONSm(2,maxdot)

        integer ndot_cc, ndot_sd
        real*8 uvx(3)
        real*8 rpr,rpr2
	    real*8 dot_vij,dot_vik,dot_vjk
        real*8 center(3)
        real*8 vpi_m(3),vpj_m(3),vpk_m(3)
        real*8 curv_triang_area
        real*8 area_cc
        real*8 dens
        real*8 vdot_cc(3,ndotccmx),vndot_cc(3,ndotccmx)
        real*8 area_dot_cc(ndotccmx)

        logical dotcc01_run
        logical OPT_CCdot_ProbProj
	    logical OPT_CCdot_TriDiv
	    logical OPT_sccdens
     	real*8 OPT_sccdenK,densSCC
        real*8 OPT_triTetMx
        integer n_cc_tridiv,n_cc_PrbPrj

        real*8 CO(3,maxatm)
        real*8 rad(maxatm)
        real*8 atom_srf(maxatm)
        integer IAS(maxatm)
        integer*2 MOLNUM(maxatm)
        integer*2 IAT(maxatm)
        logical*1 SRS(maxatm)
        integer*2 ICO(3,maxatm),ICUPTR(maxatm)
        integer*2 icol
        real*8 COMIN(3)
        integer*2 ICUBE(MAXCUB,MAXCUB,MAXCUB)
        logical*1 SCUBE(MAXCUB,MAXCUB,MAXCUB)
        logical*1 SSCUBE(MAXCUB,MAXCUB,MAXCUB)
C

        real*8 ERNBR(maxnbr)
        real*8 DISNBR(maxnbr)
        real*8 CNBR(3,maxnbr),RNBR(maxnbr)
        logical*1 SNBR(maxnbr),MNBR(maxnbr)
        integer*2 MOLNBR(maxnbr)
        integer*2 INBR(maxnbr)
        integer*2 LKNBR(maxnbr)
        integer*2 ITNL(maxnbr)
        real*8 UP(3,maxsph),TETP(maxsph)
        real*8 UA(3,maxchil,maxsph)
	real*8 dotchild(3,maxchil)
	real*8 dotchildr(3,maxchil)
	logical dotchsf(maxchil)
	logical dotchrf(maxchil)
	real*8 ARUA(maxsph)
	real*8 ARUP(maxsph)
	integer NCHIUA(maxsph)
	integer NCHmx
	integer NUV

        real*8 CI(3), CJ(3), CK(3)
        logical SI,SJ,SK,SNS

	    real*8 Q1(3),UQ1(3),UT1(3),TVP(3)
        real*8 VIJ(3),UIJ(3),Q(3),T(3),CIJK(3),VIJK(3),UIJK(3)
        real*8 BIJ(3),AIJ(3),BIJK(3),AIJK(3,2),AIJP(3,2),A(3)
        real*8 P(3),PIJ(3),PIJP(3,2),PIJK(3,2),PIPT(3)
	    real*8 PrPIJK(3)
        real*8 VPI(3),VPJ(3),VPK(3)
        real*8 VPS0(3,maxarc),VBS0(3,maxarc,2),VBS(3),ARCA(maxarc)
        real*8 VECTOR(3),centerm(3)
	    real*8 RId,RJd,HIJd,DIJd
	    real*8 RI2
	    real*8 v0_arc(3,narcmx)
        real*8 dot_xyz(3,ndotsdmx)
        real*8 dot_vn(3,ndotsdmx),dot_area(ndotsdmx)
        integer dot_type(ndotsdmx),dot_atom(ndotsdmx)

	logical dot_cross(maxarc)
	real*8 H05(3,3),GH05GT(3,3)
	real*8 rotate_angl
	real*8 arch,archI,archJ
	real*8 storanal
	real*8 areaARC,areaARC_a
	real*8 areaSDnorm
	real*8 areaSADDLE_n,areaSADDLE_a
	real*8 areaCONCAVE,areaCONC_ij
	real*8 aa1,aa2,aa3,aa4,vbkat(3)
        real*8 volume
	real*8 saddle_angl_min
	real*8 cshift_k(3),ashift
	real*8 Ar_singularSD,Ar_singularCC
	real*8 ar_singCC_ij

	integer NARCT,narci,narcj
	integer ndotSADDLE
	integer ndotCONCAVE,ndotCONC_ij
	integer N_singularSD
	integer n_singulTor
	integer N_singularCC
	integer n_singCC_ij
	integer iconcf
	integer iavert,javert,kavert
	integer nijs_face
	integer nsf,k1,k2,k3,isf,it,ich
	integer nn1,nn2,nn3
	integer putsmall,putsmall_smp
	integer arc_cross
        integer sadlsmtypei,sadlsmtypej
	integer nsmp2tot
	integer ismp2
	integer nsmp2mx,ndotsmp2mx
	integer ndotsmp2
	integer dotsmp2_nsm(ndotsmp2_max)
	integer dotsmp2_dgl(ndotsmp2_max)
	real*8  smp2xyz(3,2)
	real*8  smp2txyz(3,nsmp2_max)

	logical free_tor
	logical singularSD
	logical singularCC
	logical delete_pp
	logical doshift_k

	integer ntor_i(maxatm)
	integer Jattor_i(maxnbr,maxatm)
	integer MAXntor_i
	integer frame_cx
	real*8  dot_aw(maxsph)
	real*8  dotCO(3,maxsph)
	real*8 sizeDot,sizeDot4
	real*8 sizeBE
	real*8 Dstor_i(maxnbr,maxatm)
	real*8 delta(maxnbr)
	real*8 range_cx
	real*8 sintet
	real*8 dotvUIJ
	integer OPT_sizeBE
	integer itormx
	logical lastcut


         integer nconc_glob_mx
         integer nconc_ij_mx
	     integer ndotCCmax
         parameter (nconc_glob_mx=maxprob)
         parameter (nconc_ij_mx=maxnbr)
	     parameter (ndotCCmax=256)
         integer  nconc_glob,nconc_ij
         integer  nvertex(maxatm)
         integer  atvert_conc(maxnbr,maxatm)

         integer  conc_gl_ijk(6,nconc_glob_mx)

         integer nconc_glb_nprpos(nconc_glob_mx)

         integer  conc_ij_k(2,nconc_ij_mx)

	 integer  concf_ijf(nconc_ij_mx)

         real*8   cprobe_ijk(3,nconc_glob_mx)
         logical*1 probe_WYONRm(nconc_glob_mx)

         real*8 cvertex(3,maxnbr,maxatm)

         real*8 cvert_probe(3,3,maxprob)
         integer prob_rotConn(3,maxprob)
         real*8 prob_rotAngl(3,maxprob)
         real*8 prob_rotAltt(3,maxprob)
         integer  katt,kpos,icc
         integer OPT_printcc,OPT_printsd
	 integer OPT_printcx
	 integer OPT_printrem
	 logical OPT_overlap_RE_remove
	 logical OPT_removeCC

       integer nbr_mx,atm_mx,nconc_mx
       integer icfIJ,kref
       real*8 AIJarcS(3,maxnbr)
       real*8 GarcS(3,3,maxnbr)
       real*8 RotAng(maxnbr)
       integer SignRotS(maxnbr)

	integer nprob_mx
	integer ndot_mx
	integer nconcij_mx

          integer ndot_smooth
          integer ndot_smmx
          integer dot_smooth_gl(ndot_smoothmx)
          integer prpr1_smooth(2,ndot_smoothmx)

        real*8 G(3,3),H(3,3),GHGT(3,3),POW(3,3)

        logical*1 YONdot,YONdotSm
        integer*2 IMOL
        integer*4 NP

        logical BOTH,PAIR(2),AYON(maxarc)
        logical FOUND
        logical*1 YONPRB
        logical*1 BURY


        integer  nPY_gl(MAXYON)


        integer*4 NIAS(3)
        integer*4 NSHAPE(ndottypemx),NLOST(ndottypemx)
    	real*8 ARLOST(ndottypemx)

        integer*2 IATNUM(3),ISHAPE
        integer*2 IB
        character*80 filnam
    	character*10 number
    	real*8 OUTCO(3),OUTVEC(3)

        logical COLLID
        logical BURIED

        real*8 D,RP,WIDTH,D2,AREAC,AREAR
        real*8 RI,SUNI,radMAX,FI,SUMI,HIJ
        real*8 VECT1,VECT2,VECT3,RJ
        real*8 DIJ,F1,F2,RK,DK,RIJK,DIJK
        real*8 HIJK,AREA,SIGN,DSI,DSJ,DSK
        real*8 RK2,RIJ,AVH,ANGLE,DUIJ,HT,DP
    	real*8 ANGLE05
        real*8 DP2,RP2,X,Y,Z,SUM,F
    	real*8 FIJ,FJI
    	real*8 RNROT
    	real*8 DPt,RP2t,DP2t

    	logical compare_ijk,compare_ij
        real*8 DET,DOT,DIST,DIST2,ANORM
        integer*4 IBURY,IABLS,ISPH,K,N,I
        integer*4 J,IATOM,IPTR,J1,I1,MAXNB,I2
        integer*4 NY,NUP,NCIRC,ICI,ICJ,ICK,NNBR
        integer*4 NYsym,NYnsym,NYRm
        integer*4 NIMOL,JCK,JCJ,JCI
        integer*4 JATOM,JMOLD,IUSE,JMIN,JNBR
        integer*4 JMINNBR,LFK,L,L1,L2,MUTUAL,KNBR
        integer*4 KATOM,IP,NROT,NARC,IROT,IDX,IY
        integer*4 NV,IFRLST,JP,IHASH,IPREV,IFREE
        integer*4 NEAT,NYEAT,IPT,kolor,JMINBR,LKF

        integer ndot
    	integer icyv,icyv2,newtypeyv
    	integer idot,idot_g
    	integer*4 ndot_type(ndottypemx),dndot_type(ndottypemx)
    	real*8 surf_type(ndottypemx),dsurf_type(ndottypemx)
    	real*8 atomradyv

        real*8 dst,hst,HHst,xstL,xstR,ystL,ystR
        real*8 dxstL,dxstR
        integer kanalz1,kanalz2,kanalz3,kanalz4,kanalz5

        integer iiprint0,iiprint1,iiprint2,iiprint1m
C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	 real*8 PI2,null,onehalf,one,two,three,four,big,small
	 common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small
         real*8 nulls

         real*8 toler_bur
	     real*8 ktoler_yon
         real*8 toler2_yon
         real*8 ktoler_cross,ktoler_hole

		 real*8 hmxsmth
         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
	     common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &	 toler_cx,toler_ovr,toler_cross

         logical CONTROL
         logical fext
	     character*16 filename


          CONTROL = .false.
          kanalpl = kanalp

         if(kanalpl.ne.6)then
         open(unit=kanalpl,file='SIMSres.out',form='formatted',
     &                     status='unknown')
          end if

	 null=0.0d0
         nulls = null
	 small = 1.0d-6
	 onehalf=0.5d0
	 one = 1.0d0
	 two = 2.0d0
	 three = 3.0d0
	 four = 4.0d0
	 big = 1.0d6
         PI2 = DACOS(-one)

        IBURY = 0
        IABLS = 0
        D=dotden
        RP=d_probe*0.5d0

	 OPT_CCdot_TriDiv = .false.
	 dotcc01_run=.false.
	 OPT_CCdot_ProbProj = .not.OPT_CCdot_TriDiv
	 OPT_triTetMx = 0.00d0
         OPT_triTetMx = dcos(PI2*OPT_triTetMx)*RP**2

	 OPT_sccdenK = 1.0d0
	 OPT_sccdens = .false.
         if(OPT_sccdens)then
	 densSCC = OPT_sccdenK*dotden
	 else
	 densSCC = dotden
	 end if

c	OPT_printcc=1
	OPT_printcc=0

c	OPT_printsd=1
	OPT_printsd=0

c	OPT_printcx=1
	OPT_printcx=0

c	OPT_printrem=1
	OPT_printrem=0

	OPT_overlap_RE_remove=.true.
c	OPT_overlap_RE_remove=.false.

c	OPT_removeCC=.false.
	OPT_removeCC=.true.

	OPT_sizeBE=0
        if(OPT_refcut) OPT_sizeBE=1

        iiprint1m=3
    	iiprint0=0
    	iiprint1=1
    	iiprint2=2

    	badsurf_mx=0.005d0

        n_cc_tridiv = 0
        n_cc_PrbPrj = 0
	nprobpos = 0
	nprobpos_t = 0
	ndotccsd = 0
	ndotccsd_t = 0
        ndot_smooth = 0

	do i=1,maxprob
	    start_dotn(i)=0
        stop_dotn(i) = 0
        probe_WYONRm(i)=.false.
        do k=1,3
         prob_rotConn(k,i)= 0
         prob_rotAngl(k,i)= 0.0d0
         prob_rotAltt(k,i) = 0.0d0
        end do
	end do

	do i=1,maxatm
	nvertex(i)=0
	end do

	do i=1,maxdot
	    WDYONRm(i)=.false.
        WDYON(i) = .false.
        WDYONSm(1,i) = 0
        WDYONSm(2,i) = 0
	end do

        do i=1,ndot_smoothmx
        dot_smooth_gl(i) = 0
        prpr1_smooth(1,i) = 0
        prpr1_smooth(2,i) = 0
        end do
		nsmp2tot=0
		ndotsmp2=0
		sadlsmtypei=4
		sadlsmtypej=5
		do i=1,nsmp2_max
		do k=1,3
		smp2txyz(k,i)=0.0d0
		end do
		end do

	putsmall=1

        putsmall_smp = 1

	areaCONCAVE=0.0d0
	areaSADDLE_a=0.0d0
	areaSADDLE_n=0.0d0
	MAXntor_i = 0
	frame_cx = 0
	ndotCONCAVE=0
	ndotSADDLE=0
	probmax1 = maxprob
	maxnbr1=maxnbr
	maxdot1=maxdot
	NCHmx = maxchil

  	toler_nb=1.0d-8

        toler_bur = 0.0d0

        toler_pr=1.0d-6

  	toler_d=1.0d-6

	toler_cx = 0.0d0

        ktoler_yon = 0.0d0
	toler_yon= 0.0d0
	toler2_yon=toler_yon**2/(RP*two)

        toler_ovr = toler_yon

	ktoler_hole = 2.0d0

          toler_cross=0.0d0
ctt
c        write(kanalpl,'(a23,f8.3)')'SIMS: SolvProbRad(A)  ',RP,
c    &  'SIMS: SmoothingProbRad',rad_sm, 'SIMS: DotDens  per A^2',dotden


	saddle_angl_min = 0.05d0
        ashift=20.0d0*toler_pr

	singularSD=.false.
	singularCC=.false.
	N_singularCC=0
	N_singularSD=0
	n_singulTor=0
	Ar_singularCC=0.0d0
	Ar_singularSD=0.0d0
        radMAX = 0.0d0

        if(natoms.gt.maxatm)then
        write(kanalp,*)'FATAL ERROR: natoms > maxatm !!!'
        stop
        end if

	do n=1,natoms
        rad(n)=atom_rad(n)
        if(radMAX.lt.rad(n))radMAX=rad(n)
        end do

        WIDTH = 2 * (radMAX + RP)
        do k = 1,3
           NIAS(k) = 0
           COMIN(K) = big
	   centerm(k)= 0.0d0
        end do

	do  N=1,natoms
	IAS(N)=2
	MOLNUM(N)=1
        atom_srf(N)=0.0d0
	do K=1,3
	CO(K,N)=coords(K,N)
	centerm(k)=centerm(k)+coords(K,N)
	enddo

        do K = 1,3
           IF (CO(K,N) .LT. COMIN(K)) COMIN(K) = CO(K,N)
           IF (IAS(N) .EQ. K-1) NIAS(K) = NIAS(K) + 1
        end do
        end do

        do k=1,3
	centerm(k)=centerm(k)/natoms
	end do

C     SET UP CUBE ARRAYS
C     FIRST THE integer COORDINATE ARRAYS
        DO I = 1,natoms

           DO K = 1,3
              icol = (CO(K,I)-COMIN(K))/WIDTH
              ICO(K,I) = icol + 1

        IF (ICO(K,I) .LT. 1)then
        write(*,*)'CUBE COORDINATE TOO SMALL'
        write(*,*)'Total Atoms=',natoms

	do J=1,natoms
        write(*,'(a4,i4,a8,3f15.10)')'at=',J,
     & 'CO(K,J)=',CO(1,J),CO(2,J),CO(3,J)
        end do
        write(*,'(a12,3f15.7,a7,f10.6)')'COMIN(K)=',COMIN,
     & 'WIDTH=',WIDTH
        STOP
        end if

              IF (ICO(K,I) .GT. MAXCUB)
     1         STOP 'CUBE COORDINATE TOO LARGE'

         end do
         end do

C INITIALIZE HEAD POINTER AND SRN=2 ARRAYS

        DO  K = 1,MAXCUB
           DO J = 1,MAXCUB
              DO I = 1,MAXCUB

                 ICUBE(I,J,K) = 0
                 SCUBE(I,J,K) = .FALSE.
                 SSCUBE(I,J,K) = .FALSE.
               end do
            end do
          end do

C INITIALIZE LINKED LIST POINTERS

        DO I = 1,natoms
           ICUPTR(I) = 0
        end do

C SET UP HEAD AND LATER POINTERS FOR EACH ATOM

        DO 1250 IATOM = 1,natoms

C SKIP ATOMS WITH SURFACE REQUEST NUMBERS OF ZERO

           IF (IAS(IATOM) .EQ. 0) GO TO 1250

           I = ICO(1,IATOM)
           J = ICO(2,IATOM)
           K = ICO(3,IATOM)

           IF (ICUBE(I,J,K) .LE. 0) THEN

C     FIRST ATOM IN THIS CUBE
              ICUBE(I,J,K) = IATOM
           ELSE

C     ADD TO END OF LINKED LIST
              IPTR = ICUBE(I,J,K)
1100          CONTINUE

C CHECK FOR DUPLICATE COORDINATES

              D2 = DIST2(CO(1,IATOM),CO(1,IPTR))

              IF (D2 .LE. null) THEN

                 IAS(IATOM) = 0

                 WRITE (kanalpl,1150) IATOM,IPTR

1150             FORMAT(1X,'ATOM',I5,' DROPPED (SAME CO AS ',I5,')')

                 GO TO 1250

              END IF

              IF (ICUPTR(IPTR) .LE. 0) GO TO 1200

C MOVE ON DOWN THE LIST

              IPTR = ICUPTR(IPTR)

              GO TO 1100

1200          CONTINUE

C STORE ATOM NUMBER

              ICUPTR(IPTR) = IATOM

           END IF

C CHECK FOR SURFACED ATOM

           IF (IAS(IATOM) .EQ. 2) SCUBE(I,J,K) = .TRUE.

1250    CONTINUE

C CHECK FOR 3 X 3 X 3 WITH SOME SRN = 2

        DO 1550 K = 1,MAXCUB
           DO 1500 J = 1,MAXCUB
              DO 1450 I = 1,MAXCUB

                 IF (ICUBE(I,J,K) .EQ. 0) GO TO 1450

C CHECK WHETHER THIS CUBE OR ANY ADJACENT CUBE HAS SRN = 2

                 DO 1400 K1 = K-1,K+1

                    IF (K1 .LT. 1 .OR. K1 .GT. MAXCUB) GO TO 1400

                    DO 1350 J1 = J-1,J+1

                       IF (J1 .LT. 1 .OR. J1 .GT. MAXCUB) GO TO 1350

                       DO 1300 I1 = I-1,I+1

                       IF (I1 .LT. 1 .OR. I1 .GT. MAXCUB) GO TO 1300

                          IF (SCUBE(I1,J1,K1)) SSCUBE(I,J,K) = .TRUE.

1300                   CONTINUE

1350                CONTINUE

1400             CONTINUE

1450          CONTINUE

1500       CONTINUE

1550    CONTINUE

C INITIALIZATION
C MAXIMUM NUMBER OF NEIGHBORS ANY ATOM HAS

        MAXNB = 0

C NUMBERS OF SURFACE POINTS
        DO  K = 1,ndottypemx
           NSHAPE(K) = 0
           NLOST(K) = 0
	   ARLOST(K)=null
         end do

        NY = 0
        NYsym = 0
        NYnsym = 0

C CONTACT AND REENTRANT AREAS

        AREAC = null
        AREAR = null

C STOP IF DENSITY IS NOT POSITIVE

        IF (D .LE. null) STOP 'NON-POSITIVE DENSITY'

C INITIALIZE SOME REENTRANT SURFACE TO FALSE FOR EACH ATOM

        do IATOM = 1,natoms
           SRS(IATOM) = .FALSE.
        end do

C ...................................................................


	NUP = (4.0d0*PI2*RP**2 )*densSCC

        IF (NUP .GT. maxsph) NUP = maxsph
        IF (NUP .LT. 12) NUP = 12
	call GENUN01(RP,UP,ARUP,TETP,NUP)

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        if(CONTROL)write(*,*)'MS: precalc finish'

    	idot=0
        idot_g=0
    	nconc_glob=0

	 do i=1,maxatm
	 ntor_i(i)=0
	 do k=1,maxnbr
	 atvert_conc(k,i)=0
	 end do
	 end do
c  * * * * * * * *  * *  * * *  *   * * * * * * * * *  * * * * * ** * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C BIG LOOP FOR EACH ATOM
        do IATOM = 1, natoms

        if(CONTROL)write(*,*)'MS: bigLoop iatom=',IATOM

C SKIP IGNORED ATOMS

           IF (IAS(IATOM) .EQ. 0) GO TO 5650


           ICI = ICO(1,IATOM)
           ICJ = ICO(2,IATOM)
           ICK = ICO(3,IATOM)

C SKIP IATOM IF ITS CUBE AND ADJOINING CUBES CONTAIN ONLY BLOCKERS

           IF (.NOT. SSCUBE(ICI,ICJ,ICK)) GO TO 5650

C TRANSFER VALUES FROM LARGE ARRAYS TO IATOM VARIABLES

           RI = rad(IATOM)

           SI = IAS(IATOM) .EQ. 2

           DO  K = 1,3
              CI(K) = CO(K,IATOM)
           end do

           IMOL = MOLNUM(IATOM)

C GATHER THE NEIGHBORING ATOMS OF IATOM

C INITIALIZE NUMBER OF NEIGHBORS, AND NUMBER OF NEIGHBORS IN THE

C SAME MOLECULE AS ATOM I

           NNBR = 0
           NIMOL = 0

C INITIALIZE SRN = 2 FOR SOME NEIGHBOR TO FALSE

           SNS = .FALSE.

C SAVE A LITTLE TIME FOR DISTANCE CHECK

           SUMI = two * RP + RI

C CHECK IATOM CUBE AND ADJACENT CUBES FOR NEIGHBORING ATOMS

           DO 2550 JCK = ICK-1,ICK+1

              IF (JCK .LT. 1 .OR. JCK .GT. MAXCUB) GO TO 2550

              DO 2500 JCJ = ICJ-1,ICJ+1

                 IF (JCJ .LT. 1 .OR. JCJ .GT. MAXCUB) GO TO 2500

                 DO 2450 JCI = ICI-1,ICI+1

                    IF (JCI .LT. 1 .OR. JCI .GT. MAXCUB) GO TO 2450

                    JATOM = ICUBE(JCI,JCJ,JCK)

2300                CONTINUE

C CHECK FOR END OF LINKED LIST FOR THIS CUBE

                    IF (JATOM .LE. 0) GO TO 2400

C DISTANCE CHECK

                    SUM = SUMI + rad(JATOM) + toler_nb
c

                    VECT1 = DABS(CO(1,JATOM) - CI(1))

                    IF (VECT1 .GE. SUM) GO TO 2350

                    VECT2 = DABS(CO(2,JATOM) - CI(2))

                    IF (VECT2 .GE. SUM) GO TO 2350

                    VECT3 = DABS(CO(3,JATOM) - CI(3))

                    IF (VECT3 .GE. SUM) GO TO 2350

                    D2 = VECT1 ** 2 + VECT2 ** 2 + VECT3 ** 2

                    IF (D2 .GE. SUM ** 2) GO TO 2350

C IATOM IS NOT ITS OWN NEIGHBOR

                    IF (IATOM .EQ. JATOM) GO TO 2350

	   aa1=dsqrt(D2)
	   if(aa1+RI.le.rad(JATOM)+toler_bur)then
	   if((OPT_printcx.ge.iiprint1).or.CONTROL)then
	   write(*,*)'Iatom ',IATOM, ' is inside of Jatom',JATOM
	   end if
	   goto 5650
	   end if

	   if(aa1+rad(JATOM).le.RI+toler_bur) then
	   if((OPT_printcx.ge.iiprint1).or.CONTROL)then
	   write(*,*)'Jatom ',JATOM, ' is inside of Iatom',IATOM
	   end if
	   go to 2350
	   end if

C WE HAVE A NEW NEIGHBOR
                    NNBR = NNBR + 1

C CHECK FOR NEIGHBOR OVERFLOW
           if (NNBR .GT. maxnbr)then
           write(kanalpl,*)'FATAL ERROR: maxnbr too small'
        	  stop
           end if

C SAVE ATOM NUMBER IN TEMPORARY ARRAY

                    ITNL(NNBR) = JATOM

C CHECK WHETHER SURFACED NEIGHBOR IN SAME MOLECULE

                IF (IAS(JATOM) .EQ. 2 .AND. MOLNUM(JATOM) .EQ. IMOL)

     1       SNS = .TRUE.

C COUNT THE NUMBER OF ATOMS IN THE SAME MOLECULE AS IATOM

                    IF (IMOL .EQ. MOLNUM(JATOM)) NIMOL = NIMOL + 1

2350                CONTINUE

C GET NUMBER OF NEXT ATOM IN CUBE

                    JATOM = ICUPTR(JATOM)

                    GO TO 2300

2400                CONTINUE

2450             CONTINUE

2500          CONTINUE

2550       CONTINUE

C KEEP TRACK OF MAXIMUM NUMBER OF NEIGHBORS
C FOR ARRAY-DIMENSIONING PURPOSES

           IF (NNBR .GT. MAXNB) MAXNB = NNBR

C NO SURFACE FOR ATOM I IF BURIED ONLY FLAG SET AND
C THERE ARE NO NEIGHBORS FROM THE OTHER MOLECULES

           IF (IBURY .EQ. 1 .AND. NIMOL .EQ. NNBR) GO TO 5650

C NO SURFACE IF IATOM AND ALL NEIGHBOR ATOMS
C IN THE SAME MOLECULE HAVE SURFACE REQUEST NUMBERS < 2

           IF (.NOT. SI .AND. .NOT. SNS) GO TO 5650

C SET UP NEIGHBORS ARRAYS WITH JATOM IN INCREASING ORDER
C INITIALIZE MINIMUM NEIGHBOR ATOM NUMBER

           JMOLD = 0
           if(NNBR.GE.1)then
           DO 2700 IUSE = 1,NNBR

              JMIN = natoms + 1

              DO 2600 JNBR = 1,NNBR

C DON'T USE ONES ALREADY SORTED

                 IF (ITNL(JNBR) .LE. JMOLD) GO TO 2600

                 IF (ITNL(JNBR) .LT. JMIN) THEN

                    JMIN = ITNL(JNBR)

                    JMINBR = JNBR

                 END IF

2600          CONTINUE

              JMOLD = JMIN

              JNBR = JMINBR

              JATOM = ITNL(JNBR)

C TRANSFER ATOM NUMBER, COORDINATES, radIUS, SURFACE REQUEST NUMBER,
C MOLECULE NUMBER, EXPANDED radIUS, DISTANCE FROM IATOM

              INBR(IUSE) = JATOM

              DO K = 1,3
                 CNBR(K,IUSE) = CO(K,JATOM)
              end do

              RNBR(IUSE) = rad(JATOM)

              SNBR(IUSE) = IAS(JATOM) .EQ. 2

              MOLNBR(IUSE) = MOLNUM(JATOM)

              ERNBR(IUSE) = RNBR(IUSE) + RP

              DISNBR(IUSE) = DIST2(CI,CNBR(1,IUSE))

C INITIALIZE LINK TO NEXT FARTHEST OUT NEIGHBOR

              LKNBR(IUSE) = 0

2700       CONTINUE
           end if

C SET UP A LINKED LIST OF NEIGHBORS IN ORDER OF
C INCREASING DISTANCE FROM IATOM
C INITIALIZE POINTER TO FIRST NEIGHBOR TO 0

           LKF = 0

C LOOK FOR NEIGHBOR IN SAME MOLECULE
C WE WANT ONLY ATOMS IN SAME MOLECULE FOR COLLISION CHECK

           DO 2750 L = 1,NNBR

              IF (IMOL .NE. MOLNBR(L)) GO TO 2750

              LKF = L

              GO TO 2800

2750       CONTINUE

           IF (LKF .EQ. 0) GO TO 3000

2800       CONTINUE

C PUT REMAINING NEIGHBORS IN LINKED LIST AT PROPER POSITION

           DO 2950 L = LKF+1,NNBR

              IF (IMOL .NE. MOLNBR(L)) GO TO 2950

              L1 = 0

              L2 = LKF

2850          CONTINUE

              IF (DISNBR(L) .LT. DISNBR(L2)) GO TO 2900

              L1 = L2

              L2 = LKNBR(L2)

              IF (L2 .NE. 0) GO TO 2850

2900          CONTINUE

C ADD TO LIST

              IF (L1 .EQ. 0) THEN

                 LKF = L

                 LKNBR(L) = L2

              ELSE

                 LKNBR(L1) = L

                 LKNBR(L) = L2

              END IF

2950       CONTINUE

3000       CONTINUE

c          if(CONTROL)write(*,*)'MS:Finish neighbList of IATOM'

           IF (RP .EQ. null) GO TO 5200

           IF (NIMOL .LE. 0) GO TO 5200


C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

           if(NNBR.GE.1)then
           DO 5150 JNBR = 1,NNBR

              JATOM = INBR(JNBR)

              if(CONTROL)
     &     write(*,*)'MS:MEDIUM LOOP JNBR,JATOM=',JNBR,JATOM


              IF (JATOM .LE. IATOM) GO TO 5150

              IF (IMOL .NE. MOLNBR(JNBR)) GO TO 5150

              RJ = RNBR(JNBR)
              SJ = SNBR(JNBR)

              DO  K = 1,3
                 CJ(K) = CNBR(K,JNBR)
              end do

              DO K = 1,3
                 VIJ(K) = CJ(K) - CI(K)
              end do

	      DIJ = ANORM(VIJ)

	       IF (DIJ .LE. null) THEN

                 WRITE (kanalpl,3150) IATOM,JATOM

3150             FORMAT(1X,'ATOMS',2I5,' HAVE THE SAME CENTER')

                 GO TO 5150

              END IF

ctt
ctt             write(kanalp,*)'aft3150: VIJ=',VIJ

              CALL VNORM(VIJ,UIJ)
ctt
ctt              write(kanalp,*)'aft3150: UIJ=',UIJ
ctt

              CALL VPERP(UIJ,Q)

              CALL CROSS(UIJ,Q,T)


	      aa1 = onehalf*((RI+RP)**2-(RJ+RP)**2)/DIJ**2

	      F = onehalf + aa1
	      FIJ = F*DIJ
	      FJI = (onehalf - aa1)*DIJ


              DO K = 1,3
              BIJ(K) = CI(K) + F * VIJ(K)
              end do

ctt
ctt        write(kanalp,'(a12,3f6.3,1x,4f6.3)')'aft3200:BIJ', BIJ,CI
ctt        write(kanalp,'(a14,3f6.3,1x,f6.3)')'aft3200:VIJ,F',VIJ,F

              F1 = (RI+RJ+2*RP)**2 - DIJ**2

              IF (F1 .LE. null ) GO TO 5150
              F2 = DIJ**2 - (RI-RJ)**2

              IF (F2 .LE. null) GO TO 5150


              HIJ = DSQRT(F1*F2) / (two*DIJ)

cttsep9815
c              write(kanalp,'(a20,2i3,1x,f8.4)')'IAt,Jat, HIJ =',
c     &     IATOM, JATOM, HIJ


          if(HIJ.lt.null)then
          write(kanalpl,*)
     & 'WARNING!  HEIGHT (radIUS OF SADDLE CIRCLE) ',
     & ' too SMALL=', HIJ
          end if

              DO K = 1,3
                 AIJ(K) = HIJ * Q(K)
              end do


              MUTUAL = 0
              DO 3300 KNBR = 1, NNBR
                 D2 = DIST2(CJ,CNBR(1,KNBR))
                 MNBR(KNBR) = D2 .LT. (two*RP+RJ+RNBR(KNBR))**2

     1    .AND. KNBR .NE. JNBR

                 IF (MNBR(KNBR)) MUTUAL = MUTUAL + 1

3300          CONTINUE

C . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

              ISHAPE = 3


	       doshift_k=.false.
	       KNBR=1
               KATOM = INBR(KNBR)

	       do k=1,3
	       cshift_k(k)=null
	       end do

4020	       nconc_ij=0
	       areaCONC_ij=null
	       ndotCONC_ij=0
	       n_singCC_ij = 0
	       ar_singCC_ij = 0
		if(doshift_k)then
		 aa1=IATOM**2+JATOM**2+KATOM**2
		 cshift_k(1)=ashift*IATOM**2/aa1
		 cshift_k(2)=ashift*JATOM**2/aa1
		 cshift_k(3)=ashift*KATOM**2/aa1

ctt		    if(OPT_printcc.ge.iiprint1)then
                write(kanalpl,*)'! ! ! ! ! ! ! ! ! ! ! ! ! !'
		write(kanalpl,*)'Break Sym: DO shift for KATOM=',KATOM
ctt		write(*,*)'cshift_k(k)=',cshift_k
ctt		    end if

		do k=1,3
		CK(k)=CK(k)+cshift_k(k)
		CNBR(k,KNBR)=CK(k)
		CO(k,KATOM) = CK(k)
		end do

		doshift_k=.false.
		end if

		  if(OPT_printcc.ge.iiprint1)then
           write(kanalpl,*)' * ! * ! * ! * ! * ! * ! * ! *'
	       write(kanalpl,*)
     &    ' START TRIPLETS IATOM< JATOM=', IATOM,JATOM
	       write(kanalpl,*)
		  end if

	       KNBR=0

4200           KNBR=KNBR+1

	       if(KNBR.GT.NNBR)goto 4202

                 IF (.NOT. MNBR(KNBR)) GO TO 4200

                 KATOM = INBR(KNBR)

                 IF (KATOM .LE. JATOM) GO TO 4200
                 SK = SNBR(KNBR)

                 IF (.NOT. (SI .OR. SJ .OR. SK)) GO TO 4200

                 IF (IMOL .NE. MOLNBR(KNBR)) GO TO 4200

                 RK = RNBR(KNBR)

                 DO K = 1,3
                    CK(K) = CNBR(K,KNBR)
                 end do

		    if(OPT_printcc.ge.iiprint2)then
		 write(kanalpl,*)'XYZ I=',CI
		 write(kanalpl,*)'XYZ J=',CJ
		 write(kanalpl,*)'XYZ K=',CK
                    end if


           DK = UIJ(1) * (BIJ(1)-CK(1)) + UIJ(2) * (BIJ(2)-CK(2)) +

     1UIJ(3) * (BIJ(3)-CK(3))

                 RIJK = (RK+RP) ** 2 - DK ** 2

                 IF (RIJK .LE. null) GO TO 4200
                 RIJK = DSQRT(RIJK)

                 DO K = 1,3
                    CIJK(K) = CK(K) + DK * UIJ(K)
                 end do

                 DO K = 1,3
                    VIJK(K) = CIJK(K) - BIJ(K)
                 end do

                 DIJK = ANORM(VIJK)

                 IF (DIJK .LE. null) THEN

                    WRITE (kanalpl,3500) IATOM,JATOM,KATOM

3500                FORMAT(1X,'ATOMS',3I5,' HAVE CONCENTRIC CIRCLES')

                    GO TO 4200

                 END IF

                 F = onehalf * (one+(HIJ**2-RIJK**2)/DIJK**2)

                 DO  K = 1,3
                    BIJK(K) = BIJ(K) + F * VIJK(K)
                  end do

                 F1 = (HIJ+RIJK)**2 - DIJK**2

                 IF (F1 .LE. null) GO TO 4200

                 F2 = DIJK**2 - (HIJ-RIJK)**2

                 IF (F2 .LE. null) GO TO 4200

                 HIJK = DSQRT(F1*F2) / (two*DIJK)
				 hmxsmth = HIJK*rad_sm/(RP+rad_sm)

                 CALL VNORM(VIJK,UIJK)

                 CALL CROSS(UIJ,UIJK,AIJK)


                 DO K = 1,3
                    AIJK(K,1) = HIJK * AIJK(K,1)
                    AIJK(K,2) = - AIJK(K,1)
                 end do

                 DO 3700 IP = 1,2

                    DO K = 1,3
                       PIJK(K,IP) = BIJK(K) + AIJK(K,IP)
                    end do

            PAIR(IP) = .NOT. COLLID(PIJK(1,IP),RP,CNBR,ERNBR,MNBR,

     1NNBR,maxnbr,ISHAPE,JNBR,KNBR,MOLNBR,IMOL,LKF,LKNBR)

	     if(OPT_printcc.ge.iiprint1)then
	    write(kanalpl,*)'start IP,PAIR(IP)=',IP,PAIR(IP)
             end if

3700             CONTINUE


                 IF (.NOT. PAIR(1) .AND. .NOT. PAIR(2)) GO TO 4200

                 BOTH = PAIR(1) .AND. PAIR(2)

                 SRS(IATOM) = .TRUE.
                 SRS(JATOM) = .TRUE.
                 SRS(KATOM) = .TRUE.

                 DO 4150 IP = 1,2

                    IF (.NOT. PAIR(IP)) GO TO 4150

                    BURY = .FALSE.

           IF (IBURY .GT. 0) BURY = BURIED(PIJK(1,IP),RP,CNBR,RNB

     1R,MNBR,NNBR,maxnbr,ISHAPE,JNBR,KNBR,MOLNBR,IMOL)

                 IF (IBURY .EQ. 1 .AND. .NOT. BURY) GO TO 4150

                 YONPRB = HIJK .LE. RP+toler_cross + toler_yon

                    DO K = 1,3

                       VPI(K) = (CI(K) - PIJK(K,IP)) * RP / (RI + RP)
                       VPJ(K) = (CJ(K) - PIJK(K,IP)) * RP / (RJ + RP)
                       VPK(K) = (CK(K) - PIJK(K,IP)) * RP / (RK + RP)

                    end do

	    nconc_glob=nconc_glob+1
	    nconc_ij=nconc_ij+1
  	    nprobpos=nprobpos+1

               if(nprobpos.ge.maxprob)then
		write(kanalpl,*)'SIMS::ERROR:maxprob too small',maxprob
		stop
		end if
	  if(nconc_glob.gt.nconc_glob_mx)then
          write(kanalpl,*)'ERROR:Parameter nconc_glob_mx too SMALL'
	  stop
	  end if

	   if(nconc_ij.gt.nconc_ij_mx)then
	   write(kanalpl,*)'ERROR:Parameter nconc_ij_mx too SMALL'
	   stop
	   end if

          conc_ij_k(1,nconc_ij)=KATOM
	      conc_ij_k(2,nconc_ij)=nconc_glob


	    nvertex(IATOM)=nvertex(IATOM)+1
	    nvertex(JATOM)=nvertex(JATOM)+1
	    nvertex(KATOM)=nvertex(KATOM)+1

	    if(OPT_printcc.ge.iiprint1)then
	    write(kanalpl,*)'CONCAVE face Nglob=',nconc_glob
	    write(kanalpl,*)'concave face ij_k =',nconc_ij
	    write(kanalpl,*)'IATOM,JATOM,KATOM=',IATOM,JATOM,KATOM
	    write(kanalpl,*)'nvertex(IATOM)=',nvertex(IATOM)
	    write(kanalpl,*)'nvertex(JATOM)=',nvertex(JATOM)
	    write(kanalpl,*)'nvertex(KATOM)=',nvertex(KATOM)
	    write(kanalpl,'(a3,i3,a15,3f8.4)')
     &		' IP=',IP,' PROBE xyz =', (PIJK(k,IP),k=1,3)
            end if

	    if(nvertex(IATOM).gt.maxnbr.or.
     &         nvertex(JATOM).gt.maxnbr.or.
     &         nvertex(KATOM).gt.maxnbr)then
	    write(kanalpl,*)
     & 'NUMBER OF VERTEX on ATOM is LARGE then maxnbr'
	    stop
	    end if

	  do k=1,3
	  cprobe_ijk(k,nconc_glob)=PIJK(k,IP)
          end do
ctt
         if(OPT_printcc.ge.iiprint1)then
          write(kanalpl,'(a24,3i5,/,3f18.13)')
     & 'IATOM,JATOM,KATOM,PIJK ',
     &   IATOM,JATOM,KATOM,(PIJK(k,IP),k=1,3)
         end if

	  delete_pp=.false.
	  if(nconc_ij.gt.1)then
          do icc=1,nconc_ij-1
	  aa1=big

	  if(KATOM.ne.conc_ij_k(1,icc))then
	  aa1=null

	  do k=1,3
	  aa1=aa1+dabs(PIJK(k,IP)-cprobe_ijk(k,conc_ij_k(2,icc)))
	  end do
	  end if
	  if(aa1.lt.toler_pr)delete_pp=.true.
          end do
          end if

	  if(nconc_glob.gt.1)then
	  if(.not.delete_pp)then

          do icc=1,nconc_glob-1
	  nn1=conc_gl_ijk(1,icc)
	  nn2=conc_gl_ijk(2,icc)
	  nn3=conc_gl_ijk(3,icc)

	  if(compare_ij(IATOM,JATOM,nn1,nn2,nn3))then
	  aa1=null
	  do k=1,3
	  aa1=aa1+dabs(PIJK(k,IP)-cprobe_ijk(k,icc))
	  end do
	  if(aa1.lt.toler_pr)delete_pp=.true.
	  end if

          end do

	  end if
	  end if


	  if(delete_pp)then

            nprobpos=nprobpos - nconc_ij
            nconc_glob=nconc_glob - nconc_ij
	    nvertex(IATOM)=nvertex(IATOM) - nconc_ij
	    nvertex(JATOM)=nvertex(JATOM) - nconc_ij
	    do icc=1,nconc_ij
	    k1=conc_ij_k(1,icc)
	    nvertex(k1)=nvertex(k1) - 1
	    end do

	  nconc_ij=0
	  doshift_k=.true.

	   if(OPT_printcc.ge.iiprint1)then
	  write(kanalpl,*)'DO shift:'
	  write(kanalpl,*)'IP,IATOM,JATOM,KATOM=',IP,IATOM,JATOM,KATOM
	  write(kanalpl,*)'UPDATED nvertex = ',(nvertex(k),k=1,natoms)
	  end if

	  goto 4020
	  end if
	  conc_gl_ijk(1,nconc_glob)=IATOM
	  conc_gl_ijk(2,nconc_glob)=JATOM
	  conc_gl_ijk(3,nconc_glob)=KATOM
	  conc_gl_ijk(4,nconc_glob)=nvertex(IATOM)
	  conc_gl_ijk(5,nconc_glob)=nvertex(JATOM)
	  conc_gl_ijk(6,nconc_glob)=nvertex(KATOM)

	  nconc_glb_nprpos(nconc_glob)=nprobpos

	  atvert_conc(nvertex(IATOM),IATOM)=nconc_glob
          atvert_conc(nvertex(JATOM),JATOM)=nconc_glob
	  atvert_conc(nvertex(KATOM),KATOM)=nconc_glob


c VERTEX coordinates
            do k=1,3

	    cvertex(k,nvertex(IATOM),IATOM)= VPI(k)
	    cvertex(k,nvertex(JATOM),JATOM)= VPJ(k)
	    cvertex(k,nvertex(KATOM),KATOM)= VPK(k)

            cvert_probe(k,1,nprobpos)=VPI(k)
            cvert_probe(k,2,nprobpos)=VPJ(k)
            cvert_probe(k,3,nprobpos)=VPK(k)
	    end do

	    do k=1,3
	    center(k)=null
	    vpi_m(k)=VPI(k)
	    vpj_m(k)=VPJ(k)
	    vpk_m(k)=VPK(k)
	    PrPIJK(k)=PIJK(k,IP)
	    end do

        rpr=RP
	    rpr2=rpr*rpr
	    dens = densSCC

	    SIGN = DET(VPI,VPJ,VPK)
	    dot_vij = DOT(VPI,VPJ)
	    dot_vik = DOT(VPI,VPK)
	    dot_vjk = DOT(VPJ,VPK)

	   if(OPT_CCdot_TriDiv)then
           if(dot_vij.lt.OPT_triTetMx
     &        .or.dot_vik.lt.OPT_triTetMx
     &         .or.dot_vjk.lt.OPT_triTetMx)then
           OPT_CCdot_TriDiv = .false.
	    dotcc01_run=.false.
           end if
           OPT_CCdot_ProbProj = .not.OPT_CCdot_TriDiv
	   end if


	    if(SIGN.eq.null)then
	    write(kanalpl,*)
     & 'Warn: Concave vect VPI,VPJ,VPK are in plane !!'
	    end if

	    if(dot_vij.eq.-rpr2)then
	    write(kanalpl,*)
     &  'WARN: Concave vect VPI,VPJ are in line !!'
	    end if
	    if(dot_vik.eq.-rpr2)then
	    write(kanalpl,*)
     &  'WARN: Concave vect VPI,VPK are in line !!'
	    end if
	    if(dot_vjk.eq.-rpr2)then
	    write(kanalpl,*)
     &  'WARN: Concave vect VPJ,VPK are in line !!'
	    end if

	    if(OPT_CCdot_triDiv)then
            n_cc_tridiv=n_cc_tridiv + 1

	    call gener_conc_dot01(vpi_m,vpj_m,vpk_m,rpr,dens,
     &      area_cc,vdot_cc,area_dot_cc,vndot_cc,ndot_cc,dotcc01_run)

	    if(ndot_cc.gt.ndotccmx)then
	    write(kanalpl,*)'ERROR: param ndotccmx is small',
     &	    ' increase to',ndot_cc
	    stop
	    end if

	    end if

	    if(OPT_CCdot_ProbProj.or.(.not.dotcc01_run))then
            n_cc_PrbPrj = n_cc_PrbPrj + 1

            call gener_conc_dot02(UP,ARUP,TETP,PrPIJK,
     &            vpi_m,vpj_m,vpk_m,rpr,dens,
     &            area_cc,vdot_cc,area_dot_cc,vndot_cc,NUP,ndot_cc)


	    end if

        if(ndot_cc.gt.ndotccmx)then
        write(*,*)'MAIN-sims ERROR: ndotccmx too small!!'
        stop
        end if


	    ndotCONC_ij=ndotCONC_ij + ndot_cc

             if(OPT_printcc.ge.iiprint1)then
	    write(kanalpl,*)'conc_dot generated ',ndot_cc
             end if

	    areaCONC_ij=areaCONC_ij + area_cc

	    if(OPT_printcc.ge.iiprint1)then
ctt           iatnn = 1
ctt           if(IATOM.eq.iatnn.or.JATOM.eq.iatnn.or.KATOM.eq.iatnn)then
	    write(kanalp,*)'CONCAVE face N=',nconc_glob
	    write(kanalp,'(a22,3i5)')'IATOM,JATOM,KATOM=',
     &      IATOM,JATOM,KATOM
            write(kanalp,'(a6,3f10.6)')'VPI=',VPI
            write(kanalp,'(a6,3f10.6)')'VPJ=',VPJ
            write(kanalp,'(a6,3f10.6)')'VPK=',VPK
	    write(kanalpl,*)' IP=',IP,' PAIR(iP)=',PAIR(IP)
	    write(kanalp,*)'BOTH = ', BOTH, ' HIJK:', HIJK,
     &   		' YONPRB:',YONPRB
	    write(kanalp,'(a12,3f10.6)')
     &      'PROBE xyz =', (PIJK(k,IP),k=1,3)
	    write(kanalp,*)'curv_triang_area=',area_cc
ctt            end if
            end if

             NP = 1
             start_dotn(nprobpos)=idot_g + 1

	       do 4000 ndot=1,ndot_cc

                   aa1=null
		   do k=1,3
		   aa1=aa1+vdot_cc(k,ndot)*AIJK(k,IP)
           uvx(k)=vdot_cc(k,ndot) + PIJK(k,IP)
           end do


           aa1 = - aa1/HIJK
           YONdot = aa1.GE.(HIJK - toler_cross)

cttcsep9818
ctt        write(kanalp,'(a26,3i3,1x,f8.4,1x,1L)')
ctt     &      		        'IAT,JAT,KAT,HIJK,YONdot:',
ctt     &                      IATOM,JATOM,KATOM,HIJK,YONdot
ctt        write(kanalp,'(a23,/,4D16.8,1x,L1)')'toler_y,aa1,aa2,aa3,YONdot',
ctt     &  toler_yon,aa1,aa2,aa3,YONdot
ctt

	       if(YONdot .AND. BOTH)then
	       singularCC=.true.
	       n_singCC_ij = n_singCC_ij + 1
	       ar_singCC_ij = ar_singCC_ij + area_dot_cc(ndot)

	       if(OPT_printcc.ge.iiprint1)then
	       write(kanalpl,*)'OVERLAPPING CONCAVE dot is REMOVED'
	       write(kanalpl,*)'number I=',ndot
	       if(OPT_printcc.ge.2)then
	       write(kanalpl,*)'coordinates dot'
           write(kanalpl,'(3f8.3)')(vdot_cc(k,ndot),k=1,3)
	       end if
	       end if
	       GO TO 4000
	       end if

            NP =NP + 1
            idot_g = idot_g +1
            if(idot_g.gt.maxdot)then
            write(kanalpl,*)'ERROR:MAXDOT is too SMALL:increase value'
            stop
            end if

            ndotccsd=ndotccsd + 1
            NSHAPE(3) = NSHAPE(3) + 1

            do k=1,3
            dotcrd(k,idot_g) = uvx(k)
            dotnrm(k,idot_g) = vndot_cc(k,ndot)
            end do

            dotarea(idot_g) = area_dot_cc(ndot)
            dot_ISHAPE(idot_g) = ISHAPE

            WDYON(idot_g) = YONdot


csep9825           YONdotSm =
csep9825     &     (aa1.GE.HIJK - toler_yon - toler_cross).AND.(.not.YONdot)

		   YONdotSm =
     &     (aa1.GE.HIJK - hmxsmth).AND.(.not.YONdot).and.
     &     (HIJK-toler_cross.LT.RP)

	        if(YONdotSm .AND. BOTH)then
            ndot_smooth = ndot_smooth + 1
            WDYONSm(1,idot_g) = WDYONSm(1,idot_g) +1

            if(ndot_smooth.gt.ndot_smoothmx)then
            write(kanalpl,*)'ERROR'
            stop
            end  if

            dot_smooth_gl(ndot_smooth) = idot_g

            prpr1_smooth(1,ndot_smooth) = nprobpos

            if(IP.eq.1)then
            prpr1_smooth(2,ndot_smooth) = nprobpos + 1
            else
            prpr1_smooth(2,ndot_smooth) = nprobpos - 1
            end if

            end if


                     DSI = 0.0d0
                     DSJ = 0.0d0
                     DSK = 0.0d0
               do k=1,3
               DSI = DSI + (uvx(k) - CI(k))**2
               DSJ = DSJ + (uvx(k) - CJ(k))**2
               DSK = DSK + (uvx(k) - CK(k))**2
               end do

               DSI = dsqrt(DSI) - RI
               DSJ = dsqrt(DSJ) - RJ
               DSK = dsqrt(DSK) - RK

                  IF (DSI .LE. DSJ .AND. DSI .LE. DSK) GO TO 3850
                  IF (DSJ .LE. DSI .AND. DSJ .LE. DSK) GO TO 3900

                       IF (.NOT. SK) GO TO 4000

                       dot_IATNUM(idot_g) = KATOM

                       GO TO 3950

3850                   CONTINUE

                       IF (.NOT. SI) GO TO 4000

                       dot_IATNUM(idot_g) = IATOM

                       GO TO 3950

3900                   CONTINUE

                       IF (.NOT. SJ) GO TO 4000

                       dot_IATNUM(idot_g) = JATOM

3950                   CONTINUE

4000                CONTINUE

          NP = NP - 1
	      ndotconc_ijk(nconc_ij) = NP
	  stop_dotn(nprobpos)=start_dotn(nprobpos) +
     &                        ndotconc_ijk(nconc_ij) - 1

	  if(OPT_printcc.ge.iiprint1)then
          write(kanalpl,*)
     &  ' IATOM,JATOM,KATOM=', IATOM,JATOM,KATOM
	      write(kanalpl,*)'This CC face nconc_ij',nconc_ij,
     &  ' ProbePos=',nprobpos, ' has Ndot=',
     &    ndot_cc

caug98
         write(kanalpl,'(a5,a25,a6,a6,a6,a5,a5)')
     & 'idot', ' xyz ', ' are ', 'WDYON',' WDYSm ','DTyp',
     & ' iat '
	 do k1=start_dotn(nprobpos),stop_dotn(nprobpos)
	    write(kanalpl,'(i5,2x,3f8.3,f6.3,2x,L1,3x,i5,i5,i5)')
     &  k1,(dotcrd(k2,k1),k2=1,3),dotarea(k1),WDYON(k1),
     &  WDYONSm(1,k1),dot_ISHAPE(k1),dot_IATNUM(k1)
        end do

        write(kanalpl,*)'********************************'
        end if

        if(OPT_printcc.ge.iiprint1)then
           write(kanalpl,'(a40,4i6,a18,i6)')
     &  'end of CC:IATOM,JATOM KATOM: IP:',
     &  IATOM,JATOM,KATOM,IP, ' Nsmothed dots:',ndot_smooth
        end if


                    IF (YONPRB) THEN

		  IF (NY .GE. MAXYON)then
           write(kanalpl,*)'FATAL ERROR MAXYON too small !!'
         	    stop
	          end if

                    if(.not.BOTH)then
                       NY = NY + 1
                       nPY_gl(NY) = nprobpos

                    probe_WYONRm(nconc_glob)=.true.

                     end if

                       if(BOTH)NYsym = NYsym + 1

                    END IF

4150             CONTINUE

	       if(KNBR.le.NNBR-1)goto 4200

4202           continue


	       if(nconc_ij.ge.1)then

	      areaCONCAVE=areaCONCAVE + areaCONC_ij
	      ndotCONCAVE=ndotCONCAVE + ndotCONC_ij
	      N_singularCC=N_singularCC + n_singCC_ij
	      Ar_singularCC=Ar_singularCC + ar_singCC_ij

	       end if

C**************************************************************************
C SADDLE-SHAPED REENTRANT

              ISHAPE = 2

C

C CHECK FOR NEITHER ATOM TO BE SURFACES

              IF (.NOT. (SI .OR. SJ)) GO TO 5150

C

C SPECIAL CHECK FOR BURIED TORI

C

C IF BOTH ATOMS ARE MARKED TO BE SURFACE,

C BUT NEITHER ATOM HAS ANY REENTRANT SURFACE SO FAR

C (AFTER TRIANGLES WITH ALL KATOMS HAVE BEEN CHECKED)

C AND IF THERE IS SOME MUTUAL NEIGHBOR IN THE SAME MOLECULE

C CLOSE ENOUGH SO THAT THE torus CANNOT BE FREE,

C THEN WE KNOW THAT THIS MUST BE A BURIED torus

	      free_tor=.true.
              iconcf=0

	      if(OPT_printsd.ge.iiprint1)then
	      write(kanalpl,*)'***********************************'
	      write(kanalpl,*)'SADDLE for IATOM,JATOM=',IATOM,JATOM
	      write(kanalpl,*)'nvertex(IATOM)=',nvertex(IATOM),
     &   ' nvertex(JATOM)=',nvertex(JATOM)
	      end if

	      if(nvertex(IATOM).ge.1.and.nvertex(JATOM).ge.1)then
	      do i=1,nvertex(IATOM)
	      iconcf=atvert_conc(i,IATOM)

	      do j=1,nvertex(JATOM)
	      if(iconcf.eq.atvert_conc(j,JATOM))then
	      free_tor=.FALSE.
              iavert=i
	      javert=j
	      go to 4240
	      end if
	      end do
	      end do
	      end if

4240          continue

	      if(OPT_printsd.ge.iiprint1)then
	      write(kanalpl,*)'SADDLE: Iatom, Jatom=',IATOM,JATOM
	      write(kanalpl,*)'This FREE torus status=',free_tor
	      if(.not.free_tor)then
	      write(kanalpl,*)'This SADDLE come from CONC face=',iconcf
	      end if
	      end if
	      ntor_i(IATOM)=ntor_i(IATOM)+1
	      ntor_i(JATOM)=ntor_i(JATOM)+1
	      if(ntor_i(IATOM).gt.maxnbr.or.
     &        ntor_i(JATOM).gt.maxnbr)then
	      write(kanalpl,*)'ERROR'
	      write(kanalpl,*)
     &    'NUMBER of torus for IATOM is LARGE maxnbr!!'
	      stop
	      end if

	      Jattor_i(ntor_i(IATOM),IATOM)=JATOM
	      Jattor_i(ntor_i(JATOM),JATOM)=IATOM

	      Dstor_i(ntor_i(IATOM),IATOM)=RI*FIJ/(RI+RP)
	      Dstor_i(ntor_i(JATOM),JATOM)=RJ*FJI/(RJ+RP)


	     if(free_tor.AND. MUTUAL .GT. 0) THEN

                 DO 4250 KNBR = 1,NNBR

                    IF (.NOT. MNBR(KNBR)) GO TO 4250

                    IF (IMOL .NE. MOLNBR(KNBR)) GO TO 4250

                    D2 = DIST2(BIJ,CNBR(1,KNBR))

		    aa2=null
		    do k=1,3
		    vbkat(k)=CNBR(k,KNBR)-BIJ(k)
		    aa2=aa2+vbkat(k)**2
		    end do
		    aa1= DOT(vbkat,UIJ)

                    if(aa2.lt.toler_bur)then
		    aa3 = null
		    else
		    aa3=one - aa1**2/aa2
		    if(aa3.lt.null)aa3=null
		    aa3=dsqrt(D2*aa3)
		    end if

		    RK2=ERNBR(KNBR) ** 2 - HIJ ** 2 + two*aa3*HIJ

		 if(D2 .LT. RK2+toler_bur)then

		 if(OPT_printsd.ge.iiprint1)then
		 write(kanalpl,*)'THIS torus is BURRIED'
		 end if
        	  GO TO 5150

		 end if

4250             CONTINUE

              END IF
              if(free_tor)then

	      nijs_face=1

              RotAng(1) = 0.0d0
              RotAng(2) = two*PI2


              DO K = 1,3

                 G(K,1) = UIJ(K)
                 G(K,2) = Q(K)
                 G(K,3) = T(K)
		 AIJ(K) = HIJ*Q(K)

              end do

ctt
ctt              write(*,'(a13,9f6.3)')'FreeTr : G=',G

	      end if

	      if(.not. free_tor)then
c
	atm_mx=maxatm
	nbr_mx=maxnbr
	nconc_mx=nconc_glob_mx

      	call torus_rot_angl02(atm_mx,nbr_mx,nconc_mx,
     &          IATOM,JATOM,UIJ,HIJ,nvertex,atvert_conc,cvertex,
     &          cprobe_ijk,conc_gl_ijk,
     &          icfIJ,concf_ijf,RotAng,SignRotS,AIJarcS,GarcS)

	 nijs_face=icfIJ/2

         do i=1,nijs_face
         i1=2*i-1
         i2=i1+1
         if(prob_rotConn(1,concf_ijf(i1)).eq.0)then
         prob_rotConn(1,concf_ijf(i1)) =  concf_ijf(i2)
         prob_rotAngl(1,concf_ijf(i1)) = RotAng(i2)-RotAng(i1)
         prob_rotAltt(1,concf_ijf(i1)) = HIJ

         else

         if(prob_rotConn(2,concf_ijf(i1)).eq.0)then
         prob_rotConn(2,concf_ijf(i1)) = concf_ijf(i2)
         prob_rotAngl(2,concf_ijf(i1)) = RotAng(i2)-RotAng(i1)
         prob_rotAltt(2,concf_ijf(i1)) = HIJ

         else

         if(prob_rotConn(3,concf_ijf(i1)).eq.0)then
         prob_rotConn(3,concf_ijf(i1)) = concf_ijf(i2)
         prob_rotAngl(3,concf_ijf(i1)) = RotAng(i2)-RotAng(i1)
         prob_rotAltt(3,concf_ijf(i1)) = HIJ

         else

         write(kanalpl,*)'WARNING:prob_rotConn: probN > 3'
         end if

         end if

         end if
         if(prob_rotConn(1,concf_ijf(i2)).eq.0)then
         prob_rotConn(1,concf_ijf(i2))= concf_ijf(i1)
         prob_rotAngl(1,concf_ijf(i2))=RotAng(i1)-RotAng(i2)
         prob_rotAltt(1,concf_ijf(i2)) = HIJ

         else

         if(prob_rotConn(2,concf_ijf(i2)).eq.0)then
         prob_rotConn(2,concf_ijf(i2))=concf_ijf(i1)
         prob_rotAngl(2,concf_ijf(i2)) = RotAng(i1)-RotAng(i2)
         prob_rotAltt(2,concf_ijf(i2)) = HIJ

         else

         if(prob_rotConn(3,concf_ijf(i2)).eq.0)then
         prob_rotConn(3,concf_ijf(i2))=concf_ijf(i1)
         prob_rotAngl(3,concf_ijf(i2))=RotAng(i1)-RotAng(i2)
         prob_rotAltt(3,concf_ijf(i2)) = HIJ

         else
         write(kanalpl,*)'WARNING:prob_rotConn: probN > 3'

         end if

         end if

         end if

         end do

         end if
	      rpr=RP
	      RId=RI
	      RJd=RJ
	      HIJd=HIJ
	      DIJd=DIJ
	      dens = densSCC

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'call ARC_points: generator'
	    end if

            arc_cross = 0

         call arc_points02c(IATOM,JATOM,
     & RId,RJd,DIJd,HIJd,rpr,dens,rad_sm,
     & putsmall,putsmall_smp,NARCT,narci,narcj,
     & v0_arc,storanal,arc_cross,free_tor,
     & nijs_face,RotAng,smp2xyz,
     & dot_xyz,dot_vn,dot_area,dot_type,dot_atom,ndot_sd)

	    if(arc_cross.gt.0)n_singulTor = n_singulTor + 1

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'ARC: arc_cross:',arc_cross, ' generated'
	    end if

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'ARC: generated'
	    write(kanalpl,*)'RI,RJ,HIJ,RP=',RI,RJ,HIJ,RP
	    write(kanalpl,*)'Analyt AREA for this FULLtorus',storanal
	    write(kanalpl,*)'ARC_CROSS IJ axis=',arc_cross

	    write(kanalpl,*)'NARCT,narci,narcj=',NARCT,narci,narcj
	    end if

            if(OPT_printsd.ge.iiprint2)then
	    write(kanalpl,*)'ARC points in LOC syst:'
	    do i=1,NARCT
	    write(*,*)'xyz=',(v0_arc(k,i),k=1,3)
	    end do
	    end if

            if(NARCT .GT. maxarc)then
	    write(*,*)'FATAL ERROR'
	    stop
	    end if

	    if(NARCT.EQ.0)then
	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'SADDLE ARC: no points '
	    end if
	    goto 5150
	    end if

	    if(.not.free_tor)then
            nsf = 1
            do k1=1,3
	    AIJ(k1)= AIJarcS(k1,nsf)
	    do k2=1,3
	    G(k2,k1)=GarcS(k2,k1,nsf)
	    end do
	    end do
	    end if

	    if(arc_cross.eq.2)then
	    do ismp2=1,2
        nsmp2tot = nsmp2tot + 1
		if(nsmp2tot.gt.nsmp2_max)then
	    write(*,*)'FATAL ERROR'
	    stop
	    end if

		call MULTV(smp2xyz(1,ismp2),G,smp2txyz(1,nsmp2tot))
		do k=1,3
		smp2txyz(k,nsmp2tot) = smp2txyz(k,nsmp2tot) + BIJ(k)
		end do

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'SD-arc_cross=2:nsmp2tot:',nsmp2tot
		write(kanalpl,'(a18,3f8.3)')
     &		'initial smp2zyz:',(smp2xyz(k,ismp2),k=1,3)
		write(kanalpl,'(a18,3f8.3)')
     & 		'rotsmp2txyz:',(smp2txyz(k,nsmp2tot),k=1,3)
        end if

		end do
		end if

	       areaARC=null

            do ndot=1,ndot_sd
            idot_g = idot_g +1
            if(idot_g.gt.maxdot)then
            write(kanalpl,*)'ERROR:MAXDOT is too SMALL:increase value'
            stop
            end if

            NSHAPE(2) = NSHAPE(2) + 1
            ndotccsd=ndotccsd + 1

            call MULTV(dot_xyz(1,ndot),G,dotcrd(1,idot_g))
            do k=1,3
            dotcrd(k,idot_g) = dotcrd(k,idot_g) + BIJ(k)
            end do

            call MULTV(dot_vn(1,ndot),G,dotnrm(1,idot_g))

            dotarea(idot_g) = dot_area(ndot)
            dot_IATNUM(idot_g) = dot_atom(ndot)
            dot_ISHAPE(idot_g) = dot_type(ndot)
			if(arc_cross.eq.2)then
			if(dot_type(ndot).eq.sadlsmtypei)then
			ndotsmp2 = ndotsmp2 +1
			if(ndotsmp2.ge.ndotsmp2_max)then
			write(*,*)'FATAL ERROR'
			end if
			dotsmp2_nsm(ndotsmp2)=nsmp2tot - 1
			dotsmp2_dgl(ndotsmp2)=idot_g

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'SD-arc_cross=2:ndotsmp2:',ndotsmp2
		write(kanalpl,*)'dotsmp2_nsm(ndotsmp2):',dotsmp2_nsm(ndotsmp2)
		write(kanalpl,*)'dotsmp2_dgl(ndotsmp2):',dotsmp2_dgl(ndotsmp2)
		write(kanalpl,*)'dot_ISHAPE(idot_g):',dot_ISHAPE(idot_g)
		write(kanalpl,'(a10,3f8.3)')'dotcrd:',(dotcrd(k,idot_g),k=1,3)
        end if
			end if

        	if(dot_type(ndot).eq.sadlsmtypej)then
			ndotsmp2 = ndotsmp2 +1
			dotsmp2_nsm(ndotsmp2)=nsmp2tot
			dot_ISHAPE(idot_g) = sadlsmtypei
			dotsmp2_dgl(ndotsmp2)=idot_g

	    if(OPT_printsd.ge.iiprint1)then
	    write(kanalpl,*)'SD-arc_cross=2:ndotsmp2:',ndotsmp2
		write(kanalpl,*)'dotsmp2_nsm(ndotsmp2):',dotsmp2_nsm(ndotsmp2)
		write(kanalpl,*)'dotsmp2_dgl(ndotsmp2):',dotsmp2_dgl(ndotsmp2)
		write(kanalpl,*)'dot_ISHAPE(idot_g):',dot_ISHAPE(idot_g)
		write(kanalpl,'(a10,3f8.3)')'dotcrd:',(dotcrd(k,idot_g),k=1,3)
        end if

            end if
			end if

            end do

ctt
	    if(OPT_printsd.ge.iiprint1)then
        write(kanalpl,*)'ARC points in LOC syst :'
        write(kanalpl,*)'xyz=',(v0_arc(k,NARCT),k=1,3)
        write(kanalpl,'(a10,9f6.3)')'Rot mat G:', G
	    write(kanalpl,*)'ARC points in ROTATED LOC syst:'
	    write(kanalpl,*)'xyz=',(VBS0(k,NARCT,1),k=1,3)
        end if


                 SRS(IATOM) = .TRUE.
                 SRS(JATOM) = .TRUE.

5100          CONTINUE

	      ndotSADDLE=ndotSADDLE+ndot_sd
	      areaSADDLE_n=areaSADDLE_n+storanal

	      if(OPT_printsd.ge.iiprint1)then
	      write(kanalpl,*)'Tot area SADDLE_n=',areaSADDLE_n
	      write(kanalpl,*)'Tot dots SADDLE=',ndotSADDLE
	      write(kanalpl,*)'**********************************'
	      end if


    	nprob_mx=maxprob
    	ndot_mx = maxdot
    	nconcij_mx=nconc_ij_mx
    	nconc_mx = nconc_glob_mx
        ndot_smmx = ndot_smoothmx

    	if(.not.free_tor.and.OPT_removeCC.and.
     &	                      (arc_cross.gt.0))then

	call remove_ccdot02d(nprob_mx,ndot_mx,nconcij_mx,nconc_mx,
     &                  icfIJ,arc_cross,
     &                  RP,BIJ,UIJ,HIJ,rad_sm,dens,CI,CJ,RI,RJ,
     &                  cprobe_ijk,nconc_glb_nprpos,probe_WYONRm,
     &                  cvert_probe,
     &                  concf_ijf,RotAng,SignRotS,
     &                  start_dotn,stop_dotn,dotcrd,WDYON,WDYONRm,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl)

        if(OPT_printsd.ge.iiprint1)then
           write(kanalpl,'(a42,2i6,a18,i6)')
     &  'end of remove_ccdot02:IATOM,JATOM:',
     &  IATOM,JATOM,' smothed dots:',ndot_smooth
        end if

	end if

5150       CONTINUE
           end if

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


5200       CONTINUE


           ISHAPE = 1


           IF (.NOT. SI) GO TO 5650

	   if(OPT_printcx.ge.iiprint1)then
          write(kanalpl,*)'* * * * * * * * * * * '
	  write(kanalpl,*)'CONVEX : ATOM IATOM ',IATOM
	   end if


           IF (.NOT. BURY .AND. RP .GT. null .AND.

     1    NIMOL .GT. 0 .AND. .NOT. SRS(IATOM)) THEN

	   if(OPT_printcx.ge.iiprint1)then
	  write(kanalpl,*)
     & 'ATOM IATOM ',IATOM , ' COMPLETELY INACCESSIBLE'
	   end if

	  GO TO 5650

	  end if

	NUV = (4.0d0*PI2*RI**2 )*dotden

        IF (NUV .GT. maxsph) NUV = maxsph
        IF (NUV .LT. 12) NUV = 12
        IF (method .EQ. 1)then
       	    CALL GENUN02(RI,UA,ARUA,NCHIUA,NUV,NCHmx)
        END IF
        IF (method .EQ. 2)then
            CALL GENUN03(RI,UA,ARUA,NUV)
        END IF

           AREA = (four*PI2 * RI ** 2) / dfloat(NUV)

	   sizeDot = dsqrt(AREA)
	   sizeDot4 = 0.75d0*sizeDot
	   range_cx = one
	   RI2 = RI**2

	   do I=1,NUV
	   dot_aw(I)=one
	   end do

           IB = 0


           JNBR = 0
           KNBR = 0

	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'***********************************'
	   write(kanalpl,*)'CONVEX face IATOM=',IATOM
           end if
	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'OPT_sizeBE = 0',OPT_sizeBE
	   end if
	   if(ntor_i(IATOM).le.0)then
	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'CONVEX LOOP: atom',IATOM,' HAS no torus'
	   end if

       open(unit = kanalz, file = 'dotCO.txt', status = 'unknown')

	   do I=1,NUV
	   dot_aw(I)=ARUA(I)
	   do k=1,3
	   dotCO(k,I)=UA(k,1,I)
       write(kanalz,'(i5,5f8.3)')i,dotCO(k,I)
	   end do
	   end do
	   end if
       close(kanalz)

           if(ntor_i(IATOM).ge.1)then

	   if(MAXntor_i.lt.ntor_i(IATOM))MAXntor_i=ntor_i(IATOM)

	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'CONVEX:CUT dot from IATOM=',IATOM
	   write(kanalpl,*)'totalNumb of torus(IATOM)=',ntor_i(IATOM)
	   end if

	   if(OPT_printcx.ge.iiprint1)then
	   if(IATOM.lt.10)then
	   write(kanalpl,*)'IATOM',IATOM
	   write(kanalpl,*)'make torus with JATOM:',
     &	   (Jattor_i(k1,IATOM),k1=1,ntor_i(IATOM))
	   do k1=1,ntor_i(IATOM)
	   write(kanalpl,*)'Ds, ',Dstor_i(k1,IATOM)
	   end do
	   end if
	   write(kanalpl,*)'MAXntor_i=',MAXntor_i
	   end if

	   frame_cx = 0
	   JATOM=Jattor_i(1,IATOM)

	   do k=1,3
	   VIJ(k)=CO(k,JATOM) - CI(k)
	   end do

	   call VNORM(VIJ,UIJ)

	   if(ntor_i(IATOM).eq.1)then

	   call VPERP(UIJ,Q)
	   frame_cx=1
	   end if

	   if(ntor_i(IATOM).ge.2)then

	   KATOM=Jattor_i(2,IATOM)
	   do k=1,3
	   VIJ(k)=CO(k,KATOM)-CI(k)
	   end do

           aa1=DOT(VIJ,UIJ)
	   do k=1,3
	   VIJ(k)=VIJ(k) - aa1*UIJ(k)
	   end do

c cheak linearity
	   aa2 = ANORM(VIJ)
	   if(aa2.le.small)then

	   if(ntor_i(IATOM).ge.3)then

	   KATOM=Jattor_i(3,IATOM)
	   do k=1,3
	   VIJ(k)=CO(k,KATOM)-CI(k)
	   end do

           aa1=DOT(VIJ,UIJ)
	   do k=1,3
	   VIJ(k)=VIJ(k) - aa1*UIJ(k)
	   end do
	   call VNORM(VIJ,Q)
	   frame_cx=3

	   else
	   call VPERP(UIJ,Q)
	   frame_cx=1
	   end if

	   else
	   call VNORM(VIJ,Q)
	   frame_cx=3
	   end if

	   end if


	   call CROSS(UIJ,Q,T)


	   do k=1,3
	   G(k,3)=T(k)
	   G(k,2)=Q(k)
	   G(k,1)=UIJ(k)
	   end do

	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'frameIATOM',IATOM,' frame_axis',frame_cx

	   if(ntor_i(IATOM).lt.2)then
	   write(kanalpl,*)'One axis frame IATOM',IATOM,' ntor_i ',
     &	   ntor_i(IATOM)
	   end if

	   write(kanalpl,'(a14,3f9.6)')'rmatr:G(k,1)=',(G(k,1),k=1,3)
	   write(kanalpl,'(a14,3f9.6)')'rmatr:G(k,2)=',(G(k,2),k=1,3)
	   write(kanalpl,'(a14,3f9.6)')'rmatr:G(k,3)=',(G(k,3),k=1,3)
	   end if
c
c
	   do I=1,NUV
	   dot_aw(I)=one

           do k=1,3
	   dotCO(k,I)=null
	   end do

	   do k2=1,3
	   do k1=1,3
	   dotCO(k1,I)=dotCO(k1,I) + G(k1,k2)*UA(k2,1,I)
	   end do
	   end do

	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,'(a22,i4,1x,3f8.4)')
     &	   'initial I, dotCO=',I,(UA(k,1,I),k=1,3)
	   write(kanalpl,'(a22,i4,1x,3f8.4)')
     &	   'Rotated I, dotCO=',I,(dotCO(k,I),k=1,3)
           end if

	   end do

	   I=0

5320       I=I+1
	   if(I.gt.NUV)goto 5325

           if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'IATOM ',IATOM,' I dot=',I,'area=',ARUA(I)
             end if

           sizeDot = dsqrt(ARUA(I))
           sizeDot4 = 0.50d0*sizeDot

	   do ich=1,NCHmx
	   dotchsf(ich)=.true.
	   dotchrf(ich)=.false.

	   do k=1,3
	   dotchild(k,ich)=UA(k,ich,I)
	   dotchildr(k,ich)=null
	   end do
	   end do

	   dotchrf(1)=.true.
	   do k=1,3
	   dotchildr(k,1)=dotCO(k,I)
	   end do
	   lastcut = .false.

	   do it=1,ntor_i(IATOM)

	   if(it.eq.ntor_i(IATOM))lastcut=.true.

	   JATOM=Jattor_i(it,IATOM)
	   do k=1,3
	   VIJ(k)=CO(k,JATOM) - CI(k)
	   end do

	   call VNORM(VIJ,UIJ)

	   aa1=null
	   do k=1,3
	   aa1=aa1+UIJ(k)*dotCO(k,I)
	   end do

       dotvUIJ = aa1


	   if(OPT_sizeBE.eq.0)then
	   if(dotvUIJ.gt.Dstor_i(it,IATOM) - toler_d)then
	   dot_aw(I)=null
           if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'DOT removed '
	   end if

	   goto 5320
	   end if
           end if

	   if(OPT_sizeBE.eq.1)then
	   aa2 = one-(Dstor_i(it,IATOM)/RI)**2
	   if(aa2.lt.null)aa2=null
	   sintet = dsqrt(aa2)
	   sizeBE=sizeDot4*sintet

	   if(sizeBE.gt.null)then
           delta(it)=(dotvUIJ - Dstor_i(it,IATOM))/sizeBE
	   else
	   delta(it)=null
	   end if

	   if(delta(it) .ge. range_cx)then
	   dot_aw(I)=null
           if(OPT_printcx.ge.iiprint1)then
           write(kanalpl,*)'CUT: it:',it,' dot ',I, ' is deleted**!!'
	   end if
	   goto 5320
	   end if

	   if(delta(it) .le. - range_cx)then
           if(OPT_printcx.ge.iiprint1)then
           write(kanalpl,*)'CUT it:',it,' dot ',I, ' is accepted*!!'
	   end if

	   if(.not.lastcut) goto 5322
	   end if

           if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'it, delta(it):',it,delta(it)
	   write(kanalpl,*)' DsDot=',aa1
	   write(kanalpl,*)' sizeDot4,sizeBE=', sizeDot4,sizeBE
           end if

	  call cxdot_refcut05(G,dotchild,dotchildr,dotchrf,dotchsf,
     &	  dotCO(1,I),UIJ,Dstor_i(it,IATOM),
     &               RI,dot_aw(I),NCHmx,lastcut)


	   end if

5322       continue
	   end do
ctt
	   if(OPT_printcx.ge.iiprint1)then
	   write(*,*)'IATOM=',IATOM,'ACCEP Ndot,dot_aw=',I,dot_aw(I)
           end if

	   goto 5320

5325       continue

	   do I=1,NUV
	   dot_aw(I)=dot_aw(I)*ARUA(I)
	   end do

	   if(OPT_printcx.ge.iiprint1)then
           do I=1,NUV
	   write(*,*)'final dot I area=',I,dot_aw(I)
           end do
	   end if

	   end if

           do I = 1,NUV

	      if(dot_aw(I).le.toler_cx)goto 5600


              NSHAPE(1) = NSHAPE(1) + 1
              idot_g = idot_g + 1
            if(idot_g.gt.maxdot)then
            write(kanalpl,*)'ERROR:MAXDOT is too SMALL:increase value'
            stop
            end if


              AREAC = AREAC + dot_aw(I)

              do K = 1,3
	          dotcrd(K,idot_g) = CI(K) + dotCO(K,I)
              dotnrm(K,idot_g) = dotCO(K,I)/RI
              end do

	      dot_IATNUM(idot_g)=IATOM
	      dot_ISHAPE(idot_g)=ISHAPE
	      dotarea(idot_g)=dot_aw(I)

5600         end do

	   if(OPT_printcx.ge.iiprint1)then
	   write(kanalpl,*)'END CONVEX Surf IATOM=', IATOM
	   end if

5650       end do


	   nprobpos_t = nprobpos
	   ndotccsd_t = ndotccsd

           NYRm = 0
           do i = 1,nprobpos_t
           if(probe_WYONRm(i))then
           NYRm = NYRm + 1

ctt
           if(OPT_printrem.ge.iiprint1)then
           write(kanalpl,*)
           write(kanalpl,'(a25,i6,1x,3f8.3,a20,1x,L1)')
     &     'DeepProb I, xyz =', i,(cprobe_ijk(k,i),k=1,3),
     &     ' probe_WYONRm(i):',probe_WYONRm(i)
           write(kanalpl,'(a40,i6,2x,3i6)')
     &     'DeepPprob_rotConn:nP,n1,n2,n3 ',
     &     i, (prob_rotConn(k,i),k=1,3)
           write(kanalpl,'(a40,i6,2x,3f8.4)')
     &     'DeepPprob_rotAngl:ip,a1,a2,a3 ',
     &     i, (prob_rotAngl(k,i),k=1,3)
           write(kanalpl,'(a40,i6,2x,3f8.4)')
     &     'DeepPprob_rotAltt:ip,h1,h2,h3 ',
     &     i, (prob_rotAltt(k,i),k=1,3)
           end if

           else

           if(OPT_printrem.ge.iiprint1)then
           write(kanalpl,*)
           write(kanalpl,'(a25,i6,1x,3f10.5)')
     &     'ToruProb I, xyz =', i,(cprobe_ijk(k,i),k=1,3)
           write(kanalpl,'(a40,i6,2x,3i6)')
     &    'TORPprob_rotConn: nP,n1,n2,n3 ',
     &     i, (prob_rotConn(k,i),k=1,3)
           write(kanalpl,'(a40,i6,2x,3f8.4)')
     &     'TORPPprob_rotAngl:ip,a1,a2,a3 ',
     &     i, (prob_rotAngl(k,i),k=1,3)
           write(kanalpl,'(a40,i6,2x,3f8.4)')
     &     'TORPprob_rotAltt:ip,h1,h2,h3 ',
     &     i, (prob_rotAltt(k,i),k=1,3)
           end if

           end if

           end do

	   if(OPT_printrem.ge.iiprint1)then
	   write(kanalpl,'(a22,i6,a36,i5)')
     &     'SIMS: Total nprobpos ',nprobpos,
     &     'Numb of IndefiniteDeepProbe ',NYRm
	   end if

	   if(OPT_printrem.gt.iiprint1)then

	   do k=1,nprobpos_t
           write(kanalpl,*)' probe-pos', k
           write(kanalpl,*)'start,stop_dotn = ',start_dotn(k),
     & stop_dotn(k)

	   write(kanalpl,*)'dots:'

	   do k1=start_dotn(k),stop_dotn(k)
	   write(kanalpl,'(i5,2x,3f14.10,f7.5,2x,L2)')
     &  k1,(dotcrd(k2,k1),k2=1,3),dotarea(k1),WDYONRm(k1)
	   end do
	   end do
	   end if

	    nprob_mx=maxprob
	    ndot_mx = maxdot
	    nconc_mx = nconc_glob_mx
        ndot_smmx = ndot_smoothmx

            call remove_ccdot03b(nprob_mx,ndot_mx,nconc_mx,nprobpos,
     &                  RP,rad_sm,dens,
     &                  conc_gl_ijk,nconc_glb_nprpos,
     &                  cprobe_ijk,cvert_probe,probe_WYONRm,
     &                  start_dotn,stop_dotn,dotcrd,WDYON,WDYONRm,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl)


           call smooth_ccdot02c(nprob_mx,ndot_mx,
     &                  RP,rad_sm,
     &                  cprobe_ijk,
     &                  dotcrd,dotnrm,dotarea,dot_ISHAPE,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl,WDYONRm)

		   if(nsmp2tot.ge.2)then
     		 nsmp2mx = nsmp2_max
	    	 ndotsmp2mx = ndotsmp2_max
			call remove_sdsmp2(ndot_mx,nsmp2mx,ndotsmp2mx,
     &                  nsmp2tot,ndotsmp2,rad_sm,
     &                  dotsmp2_dgl,dotsmp2_nsm,smp2txyz,
     &                  dotcrd,dotnrm,dotarea,dot_ISHAPE,
     &                  WDYONRm)

	    	 end if

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          ndots=idot_g
	      surf1 = null
	      surf0 = null
              surf  = null
       do k=1,3
       VECTOR(k)=0.0d0
       end do

       do id = 1,ndots
       if(.not.WDYONrm(id))then
       surf0=surf0 + dotarea(id)
       do k=1,3
        VECTOR(k)=VECTOR(k) + dotcrd(k,id)*dotarea(id)
       end do
       end if
       end do

        if(surf0.gt.0.0d0)then
        do k=1,3
        VECTOR(k)=VECTOR(k)/surf0
        end do
        end if

        volume=0.0d0
        do id = 1,ndots
       if(.not.WDYONrm(id))then
        aa1=0.0d0
        do k=1,3
         aa1=aa1 + (dotcrd(k,id)-VECTOR(k))*dotnrm(k,id)
        end do
        volume=volume + aa1*dotarea(id)
	end if
        end do
        volume=volume/3.0d0

       do k=1,ndottypemx
       surf_type(k)=null
       dsurf_type(k)=null
       ndot_type(k)=0
       dndot_type(k)=0
       end do

       athreshmx(1) = 0.0d0
       athreshmx(2) = 0.0d0
       athreshmx(3) = 0.0d0
       athreshmx(4) = 0.0d0
       athreshmx(5) = 0.0d0

	ndotrem=0
        nnnd=ndots
        ndots=0

	if(outdetl.gt.iiprint1)then
        write(kanalpl,*)
	write(kanalpl,*)'REMOVED small area DOTS'
	write(kanalpl,*)'       id  iat type  area               xyz',
     &  '                 nvect'
	end if

        do id=1,nnnd
        k=dot_ISHAPE(id)
       if(WDYONRm(id))then
       NLOST(k) = NLOST(k) + 1
       ARLOST(k) = ARLOST(k) + dotarea(id)
       end if

       if(.not.WDYONRm(id))then

       if(dotarea(id).le.athreshmx(k)/dotden)then
       dsurf_type(k)=dsurf_type(k) + dotarea(id)
       dndot_type(k)=dndot_type(k) + 1
       surf1=surf1+dotarea(id)
       ndotrem=ndotrem + 1
        if(outdetl.gt.iiprint1)then
	write(kanalpl,93)id,dot_IATNUM(id),dot_ISHAPE(id),dotarea(id),
     &  (dotcrd(k1,id),k1=1,3),(dotnrm(k1,id),k1=1,3)
        end if

        else
        ndots=ndots+1
        dot_ISHAPE(ndots)=dot_ISHAPE(id)
        dot_IATNUM(ndots)=dot_IATNUM(id)
	dotarea(ndots)=dotarea(id)
	surf = surf + dotarea(ndots)

	surf_type(dot_ISHAPE(ndots))=surf_type(dot_ISHAPE(ndots))+
     &  dotarea(ndots)
	ndot_type(dot_ISHAPE(ndots))=ndot_type(dot_ISHAPE(ndots))+1

        atom_srf(dot_IATNUM(ndots))=atom_srf(dot_IATNUM(ndots))
     &  + dotarea(ndots)

        do k1=1,3
	dotcrd(k1,ndots)=dotcrd(k1,id)
	dotnrm(k1,ndots)=dotnrm(k1,id)
        end do

	end if
        end if

	end do

        areaRen=surf0/surf
	do id=1,ndots
        dotarea(id)=dotarea(id)*areaRen
	end do


        call order_dots(natoms,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,
     &                 dot_num_atom,dot_startn_atom,ndots)


	call surf_kin(coords,
     &                 atnamel,resnumbl,resnamel,
     &                 natoms,atom_rad,dotden,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,
     &                 dot_num_atom,dot_startn_atom,ndots)




C * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * * * * * * *
	if(OPT_dot_file)then
	open(unit=kanalz,file='SIMS_dot.xyz',status='unknown')
	write(kanalz,'(a21)')'#SIMS: RESULT DOTS'
	write(kanalz,'(a54,a17)')
     & '#DOT    id   iat type  area    x       y       z      ',
     & 'nvx    nvy    nvz'

	do id = 1, ndots

	write(kanalz,'(a4,1x,i5,1x,i4,1x,i3,1x,f6.3,1x,3f8.3,1x,3f7.3)')
     &   'DOT:',id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  (dotcrd(k1,id),k1=1,3)
     & ,(dotnrm(k1,id),k1=1,3)
	end do

        close(kanalz)
        end if

	if(OPT_dot_midas)then
	fext=.false.
	filename='part_i.ms'
cinquire(file=filename, exist=fext)
cif(fext)then
cfilename='SIMS_midas_2.xyz'
cend if
	open(unit=kanalz,file=filename,status='unknown')
cwrite(kanalz,'(a34)')'#SIMS: RESULT DOTS in Midas format'

	normbad = 0
	aa2=0.0d0
        do i = 1,natoms
        write(kanalz,'(a3,i5,2x,a3,f8.3,1x,f8.3,1x,f8.3,1x,a1)')
     &	resnamel(i),resnumbl(i),atnamel(i),(coords(k,i),k=1,3),'A'

        do id=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1
        write(kanalz,'(a3,i5,2x,a3,f8.3,1x,f8.3,1x,f8.3,1x,
     &     a2,i1, 1x,f7.5,1x,f6.3,1x,f6.3,1x,f6.3)')
     &  resnamel(i),resnumbl(i),atnamel(i),(dotcrd(k,id),k=1,3),
     &  'ST',dot_ISHAPE(id), dotarea(id),(dotnrm(k,id),k=1,3)

	aa2=aa2+dotarea(id)
        aa1=0.0d0
        do k1=1,3
	aa1=aa1 + dotnrm(k1,id)*dotnrm(k1,id)
        end do
        if(abs(1.0d0-aa1).ge.0.0001)then
        normbad=normbad+1
        write(kanalz,'(a25,f8.5,a6,i6)')
     &   'WARNING: BAD NORM: norm:',aa1,' nbad:',normbad
	 end if

        end do
        end do
c       write(kanalz,'(a7, f12.4)')'#Area: ',aa2
        close(kanalz)
        end if

	if(OPT_sterdot_file)then
        kanalz1 = kanalz + 1
        kanalz2 = kanalz + 2
        kanalz3 = kanalz + 3
        kanalz4 = kanalz + 4
	kanalz5 = kanalz + 5

	open(unit=kanalz,file='SIMS_sterdot.xyz',status='unknown')
	open(unit=kanalz1,file='SIMS_sterdot1.xyz',status='unknown')
	open(unit=kanalz2,file='SIMS_sterdot2.xyz',status='unknown')
	open(unit=kanalz3,file='SIMS_sterdot3.xyz',status='unknown')
	open(unit=kanalz4,file='SIMS_sterdot4.xyz',status='unknown')
	open(unit=kanalz5,file='SIMS_steratm.xyz',status='unknown')
        HHst = 25.0/ster_ZOOM
        hst = 12.0
        dst = 3.4
        dxstL=0.0
        dxstR=0.0

        write(kanalz,*)
	write(kanalz,'(a21)')'#SIMS: RESULT DOTS STEREO-PAIR'
	write(kanalz,'(a61)')
     & '#DOT:    id   iat type  area    xL      yL      xR      yR'

	iatom=1
	showrad=6.0d0**2

        do i = 1,natoms
	shdist=0.0d0
        do k =1,3
        uvx(k) = coords(k,i) - VECTOR(k)
	shdist=shdist+(coords(k,i)-coords(k,iatom))**2
        end do
	if(shdist.lt.showrad)then
        if(dot_num_atom(i).gt.0)then
        dxstL = (uvx(1) + dst)*hst/(HHst+hst+uvx(3))
        xstL = -dst +dxstL
        ystL = uvx(2)*dsqrt(hst*hst + (dst+xstL)**2)/
     &         dsqrt((HHst+hst+uvx(3))**2+(uvx(1)+dst)**2)


        dxstR = (uvx(1) - dst)*hst/(HHst+hst+uvx(3))
        xstR = dst + dxstR
        ystR = uvx(2)*dsqrt(hst*hst + (xstR-dst)**2)/
     &         dsqrt((HHst+hst+uvx(3))**2+(uvx(1)-dst)**2)

	if(dxstR.gt.-dst.and.dxstL.lt.dst)then
	write(kanalz5,'(1x,i5,1x,i5,2x,f5.2,1x,i5,1x,2f8.2,1x,2f8.2)')
     &   i,iatom,dsqrt(shdist),
     &   dot_num_atom(i),
     &  xstL, ystL,xstR,ystR
	 end if

        do id=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1
	shdist=0.0d0
        do k =1,3
        uvx(k) = dotcrd(k,id) - VECTOR(k)
	shdist=shdist+(dotcrd(k,id)-coords(k,iatom))**2
        end do
        dxstL = (uvx(1) + dst)*hst/(HHst+hst+uvx(3))
        xstL = -dst +dxstL
        ystL = uvx(2)*dsqrt(hst*hst + (dst+xstL)**2)/
     &         dsqrt((HHst+hst+uvx(3))**2+(uvx(1)+dst)**2)


        dxstR = (uvx(1) - dst)*hst/(HHst+hst+uvx(3))
        xstR = dst + dxstR
        ystR = uvx(2)*dsqrt(hst*hst + (xstR-dst)**2)/
     &         dsqrt((HHst+hst+uvx(3))**2+(uvx(1)-dst)**2)

        if(dxstR.gt.-dst.and.dxstL.lt.dst)then
	write(kanalz,'(1x,i5,1x,i4,1x,i3,1x,f6.3,1x,2f8.3,1x,2f8.3)')
     &   id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  xstL, ystL,xstR,ystR

        if(dot_ISHAPE(id).eq.1)then
	write(kanalz1,'(1x,i5,1x,i4,1x,i3,1x,f6.3,1x,2f8.3,1x,2f8.3)')
     &  id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  xstL, ystL,xstR,ystR
        end if

        if(dot_ISHAPE(id).eq.2)then
	write(kanalz2,'(1x,i5,1x,i4,1x,i3,1x,f6.3,1x,2f8.3,1x,2f8.3)')
     &  id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  xstL, ystL,xstR,ystR
        end if

        if(dot_ISHAPE(id).eq.3)then
	write(kanalz3,'(1x,i5,1x,i4,1x,i3,1x,f6.3,1x,2f8.3,1x,2f8.3)')
     &  id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  xstL, ystL,xstR,ystR
        end if

        if(dot_ISHAPE(id).eq.4)then
	write(kanalz4,'(1x,i5,1x,i4,1x,i3,1x,f6.3,1x,2f8.3,1x,2f8.3)')
     &  id,dot_IATNUM(id),dot_ISHAPE(id),
     &   dotarea(id),
     &  xstL, ystL,xstR,ystR
        end if

        end if
	end do
	end if
	end if
	end do

        close(kanalz)
        close(kanalz1)
        close(kanalz2)
        close(kanalz3)
        close(kanalz4)
	    close(kanalz5)
        end if
C
C
CSurface file for insigthII
C
	if(OPT_dot_surface)then
	open(unit=kanals,file='SIMS_dot.surf',status='unknown')
	dotsch='DOTS  '
	kref=2
	atomref='ATOM_REF '//resnamel(kref)//':'//'2'//':'//'N'

c. xyz of a reference atom for outdot.surf file
        do k=1,natoms
        if(atnamel(k).eq.'N  '.and.resnumbl(k).eq.kref)then
        atomref='ATOM_REF '//resnamel(k)//':'//'2'//':'//'N  '
        write(kanals,'(a6)')dotsch
        write(kanals,'(a18,3f8.3)')atomref,(coords(j,k),j=1,3)
        end if
        end do

c DOT Surface file for insigthII

       do id=1,ndots
       if(dot_ISHAPE(id).eq.1)kolor=60
       if(dot_ISHAPE(id).eq.2)kolor=120
       if(dot_ISHAPE(id).eq.3)kolor=20
       if(dotarea(id).lt.0.0D0)kolor=360
       write(kanals,'(3f12.6,i5)') (dotcrd(k,id),
     & k=1,3),kolor
       end do
       close(kanals)
       end if

c MOLXYZ.PDB file consisten with SURFACE
	if(OPT_dot_surface.or.OPT_pdb_surface)then
	open(unit=kanalx,file='SIMS_mol.pdb',status='unknown')
        satom='ATOM'
	write(kanalx,'(a64)')
     &'REMARK: COORDs consistent with dot SURFACE, AtRad, AtomSurf, A^2'
        write(kanalx,'(a58,a11)')
     &'REMARK                            X       Y       Z    Rat',
     &'   SimsArea'
        do k=1,natoms
	write (kanalx,7071)satom,k,atnamel(k),resnamel(k),
     &      resnumbl(k),
     &      (coords(j,k),j=1,3), atom_rad(k),atom_srf(k)
         end do
       close(kanalx)
       end if

cSIMS02:SIMS_VMDdot.pdb file output
	if(OPT_VMDpdb_surf)then
	open(unit=kanalx,file='SIMS_VMDdot.pdb',status='unknown')
	satom='ATOM'
	atnamel(1)='XE'
	resnamel(1)='XE'
	write(kanalx,'(a42)')
     & 'REMARK: dot SURFACE in pdb format for VMD '
        write(kanalx,'(a59,a11)')
     & 'REMARK                            X       Y       Z    Adot'
        do id=1,ndots
        write (kanalx,7071)satom,id,atnamel(1),resnamel(1),
     &      id,
     &  (dotcrd(k1,id),k1=1,3),
     &  dotarea(id)
        end do
        close(kanalx)
        end if

	if(outdetl.ge.iiprint1m)then
        write(kanalpl,*)
	write(kanalpl,'(a13,i6)')'SIMS: natom:',natoms
	write(kanalpl,'(a22,3f10.6)')'SIMS:Molecule Center:',centerm
        write(kanalpl,'(a22,3f10.6)')'SIMS:Surface  Center:',VECTOR
        write(kanalpl,'(a31,i6)')'SIMS:Total number of Dots    :',ndots
	write(kanalpl,'(a31,f10.3)')
     &  'SIMS:Total Surface area, A^2 :',surf0
	write(kanalpl,'(a31,f12.3)')
     &  'SIMS:Solvet-Excluded Vol A^3 :',volume

        if(outdetl.ge.iiprint1)then
        write(kanalpl,*)
	write(kanalpl,'(a35,4(1x,i5,7x))')
     &  'SIMS:Deleted Ndot    types(1,2,3,4) ',
     &	dndot_type
	write(kanalpl,'(a35,4f12.5)')
     &	'SIMS:Deleted DotArea types(1,2,3,4) ', dsurf_type
	write(kanalpl,'(a35,2f12.5)')
     &	'SIMS:Deleted DotAreaTot  areaReNorm:', surf1,areaRen
	write(kanalpl,*)
	write(kanalpl,*)'SIMS: method of CONVEX dot cut:OPT_sizeBE',
     &  OPT_sizeBE
        end if

        write(kanalpl,*)
	write(kanalpl,'(a30,5f10.3)')
     & 'SIMS:DotArea types(1,2,3,4,5):',surf_type
        write(kanalpl,'(a30,5(1x,i5,4x))')
     & 'SIMS:  Ndot  types(1,2,3,4,5):',ndot_type
        write(kanalpl,'(a30,i5)')
     & 'SIMS:         N_singularTorus:',n_singulTor

       if(outdetl.ge.iiprint1)then
        write(kanalpl,*)
	write(kanalpl,'(a41,2f12.5)')
     & 'SIMS: Area of SADDLE faces (anlyt,numr) ',
     &  areaSADDLE_a,areaSADDLE_n
	write(kanalpl,'(a40,f12.5)')
     &  'SIMS: ANALYTICAL area of CONCAVE faces:',
     &  areaCONCAVE
        write(kanalpl,'(a19,i6,a13,i6)')
     &  'SIMS: N_CC_triDiv=', n_cc_tridiv,' N_CC_PrbPrj=',n_cc_PrbPrj
        end if

caug98
        write(kanalpl,*)
	write(kanalpl,'(a26,i5,a10,f9.4)')
     &  'SIMS:DotTypeRemoved: CX :', NLOST(1),' AreaCX: ',ARLOST(1)
        write(kanalpl,'(a26,i5,a10,f9.4)')
     &  'SIMS:DotTypeRemoved: SD :', NLOST(2),' AreaSD: ',ARLOST(2)
        write(kanalpl,'(a26,i5,a10,f9.4)')
     &  'SIMS:DotTypeRemoved: CC :', NLOST(3),' AreaCC: ',
     &  ARLOST(3)
        write(kanalpl,'(a26,i5,a10,f9.4)')
     &  'SIMS:DotTypeRemoved: ST4:', NLOST(4),' AreaST: ',ARLOST(4)
         write(kanalpl,'(a26,i5,a10,f9.4)')
     &  'SIMS:DotTypeRemoved: ST5:', NLOST(5),' AreaST: ',ARLOST(5)


	write(kanalpl,*)
	write(kanalpl,'(a27,i5,a21,f9.4)')
     &  'SIMS:SymSinglDotCCRemoved:',N_singularCC,
     &  'AreaSymSingul_CC:  ',AR_singularCC

caug98
        write(kanalpl,*)
        write(kanalpl,'(a32,2i5)')
     &  'SIMS:Deep probe: NYnsym, NYsym:',
     &  NY,NYsym

	badsurf=surf1/surf

	if(badsurf.gt.badsurf_mx)then
	write(kanalpl,*)
	write(kanalpl,*)
     & '!!! WARNINIG !!! too much surf is removed %reltval=',badsurf
        end if
        end if

c       write(kanalpl,*)
c       write(kanalpl,*)'SIMS: normalEND of SIMS ...'

93      format(i5,1x,i4,i2,2x,f8.5,2x,3f8.3,2x,3f7.3,1x,f6.4)
7071    format(a4,3x,i4,2x,a3,1x,a3,2x,i4,4x,3f8.3,f6.3,f8.3)


	RETURN
        END
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
        real*8 FUNCTION DIST(A,B)
        implicit none
        real*8 A(3)
        real*8 B(3)

        DIST = DSQRT((A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2)

        return
        END

        real*8 FUNCTION DIST2(A,B)

        implicit none
        real*8 A(3)
        real*8 B(3)

        DIST2 = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2

        return
        end

        real*8 FUNCTION ANORM(A)

        implicit none
        REAL*8 A(3)

        ANORM = DSQRT(A(1)**2 + A(2)**2 + A(3)**2)

        return
        end

         real*8 FUNCTION DOT(A,B)

        implicit none
        REAL*8 A(3)
        REAL*8 B(3)

        DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)

           RETURN
           END


        SUBROUTINE CROSS(A,B,C)

        implicit none
        REAL*8 A(3)
        REAL*8 B(3)
        REAL*8 C(3)

        C(1) = A(2) * B(3) - A(3) * B(2)
        C(2) = A(3) * B(1) - A(1) * B(3)
        C(3) = A(1) * B(2) - A(2) * B(1)

        RETURN
        END

        SUBROUTINE MULTV(V,A,W)

        implicit none
        REAL*8 A(3,3)
        REAL*8 V(3)
        REAL*8 W(3),U(3)
	integer I

        DO I = 1, 3
           W(I) = A(I,1)*V(1) + A(I,2)*V(2) + A(I,3)*V(3)
        end do

        RETURN
        END

        SUBROUTINE VNORM(A,B)

        implicit none
        REAL*8 A(3),B(3)
	real*8 V
	integer K

        V = DSQRT(A(1)**2 + A(2)**2 + A(3)**2)

        DO  K = 1,3
           B(K) = A(K) / V
        end do

        RETURN
        END

        SUBROUTINE VPERP(A,B)
	implicit none
        REAL*8 A(3)
        REAL*8 B(3)
        REAL*8 P(3)

	real*8 SMALL,DT
	integer M,K

        SMALL = 10000.0d0
        M = 0
        DO K = 1,3
           IF (DABS(A(K)) .GE. SMALL) GO TO 50
           SMALL = DABS(A(K))
           M = K
50       end do

        DO K = 1,3
           B(K) = 0.0d0
           IF (K .EQ. M) B(K) = 1.0d0
         end do

        DT = A(M) / (A(1)**2 + A(2)**2 + A(3)**2)

        DO K = 1, 3
           P(K) = DT * A(K)
           B(K) = B(K) - P(K)
        end do

        CALL VNORM(B,B)

        RETURN
        END

        SUBROUTINE CAT(A,B)

        implicit none
        REAL*8 A(3,3)
        REAL*8 B(3,3)
        REAL*8 TEMP(3,3)
	integer I,J

        DO I = 1,3

           DO J = 1,3

           TEMP(I,J) = A(I,1)*B(1,J) + A(I,2)*B(2,J) + A(I,3)*B(3,J)

           end do
        end do

        DO I = 1,3
           DO J = 1,3

              A(I,J) = TEMP(I,J)

           end do
        end do

        RETURN
        END

        SUBROUTINE CONJ(H,G,GHGT)

        implicit none
        REAL*8 G(3,3)
        REAL*8 H(3,3)
        REAL*8 GHGT(3,3)
        REAL*8 GT(3,3)

	integer K,L

        CALL IMATX(GHGT)
        CALL CAT(GHGT,G)
        CALL CAT(GHGT,H)

        DO K = 1,3
           DO L = 1,3
              GT(K,L) = G(L,K)
           end do
         end do

        CALL CAT(GHGT,GT)

        RETURN
        END

        SUBROUTINE IMATX(A)

        implicit none
        REAL*8 A(3,3)
	integer I,J

        DO I = 1,3

           DO J = 1,3

              A(I,J) = 0.0d0

           end do

           A(I,I) = 1.0d0

        end do

        RETURN
        END

        real*8 FUNCTION DET(A,B,C)

        implicit none
	real*8 DOT
        REAL*8 A(3)
        REAL*8 B(3)
        REAL*8 C(3)
        REAL*8 AB(3)

        CALL CROSS(A,B,AB)

        DET = DOT(AB,C)

        RETURN
        END

        LOGICAL FUNCTION COLLID(P,RP,CNBR,ERNBR,MNBR,NNBR,MAXNBR,ISHAPE,

     1   JNBR,KNBR,MOLNBR,IMOL,LKF,LKNBR)


	implicit none

	integer*4 NNBR,JNBR,KNBR,MAXNBR,LKF
	real*8 toler
        REAL*8 P(3),RP
        REAL*8 CNBR(3,MAXNBR)
        REAL*8 ERNBR(MAXNBR)
        LOGICAL*1 MNBR(MAXNBR)
        INTEGER*2 MOLNBR(MAXNBR)
        INTEGER*2 LKNBR(MAXNBR)
        INTEGER*2 IMOL,ISHAPE
	integer I
	real*8 VECT1,VECT2,VECT3,SR2,DD2

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross


	toler= toler_nb

        I = LKF

        GO TO 100

50      CONTINUE

        I = LKNBR(I)

100     CONTINUE

        IF (I .EQ. 0) GO TO 150

        VECT1 = DABS(P(1) - CNBR(1,I)) + toler

        IF (VECT1 .GE. ERNBR(I)) GO TO 50

        VECT2 = DABS(P(2) - CNBR(2,I)) + toler

        IF (VECT2 .GE. ERNBR(I)) GO TO 50

        VECT3 = DABS(P(3) - CNBR(3,I)) + toler

        IF (VECT3 .GE. ERNBR(I)) GO TO 50

        IF (I .EQ. JNBR .OR. I .EQ. KNBR) GO TO 50

        SR2 = ERNBR(I) ** 2

        DD2 = VECT1 ** 2 + VECT2 ** 2 + VECT3 ** 2

        IF (DD2 .GE. SR2) GO TO 50

        COLLID = .TRUE.

        RETURN

150     CONTINUE

        COLLID = .FALSE.

        RETURN

        END


        LOGICAL FUNCTION BURIED(P,RP,CNBR,RNBR,MNBR,NNBR,MAXNBR,ISHAPE,

     1   JNBR,KNBR,MOLNBR,IMOL)



	implicit none

	integer*4 NNBR,JNBR,KNBR,MAXNBR

        REAL*8 P(3),RP

        REAL*8 CNBR(3,MAXNBR)

        REAL*8 RNBR(MAXNBR)

        LOGICAL*1 MNBR(MAXNBR)

        INTEGER*2 MOLNBR(MAXNBR)

        INTEGER*2 IMOL,ISHAPE
	integer I
	real*8 VECT1,VECT2,VECT3,SR2,DD2,SUMRAD


        IF (NNBR .LE. 0) GO TO 100


        DO 50 I = 1, NNBR

           IF (IMOL .EQ. MOLNBR(I)) GO TO 50

           IF (ISHAPE .GT. 1 .AND. I .EQ. JNBR) GO TO 50

           IF (ISHAPE .EQ. 3 .AND. (I .EQ. KNBR .OR. .NOT. MNBR(I)))

     1   GO TO 50

           SUMRAD = RP + RNBR(I)

           VECT1 = DABS(P(1) - CNBR(1,I))

           IF (VECT1 .GE. SUMRAD) GO TO 50

           VECT2 = DABS(P(2) - CNBR(2,I))

           IF (VECT2 .GE. SUMRAD) GO TO 50

           VECT3 = DABS(P(3) - CNBR(3,I))

           IF (VECT3 .GE. SUMRAD) GO TO 50

           SR2 = SUMRAD ** 2

           DD2 = VECT1 ** 2 + VECT2 ** 2 + VECT3 ** 2

           IF (DD2 .LT. SR2) GO TO 150

50      CONTINUE

100     CONTINUE

        BURIED = .FALSE.

        GO TO 200

150     CONTINUE

        BURIED = .TRUE.

200     CONTINUE

        RETURN

        END


        SUBROUTINE GENUN(U,AR,N)

	implicit none
	integer N
        REAL*8 U(3,N)
	real*8 AR(N)

	real*8 PI,FI,FJ,Z,XY,X,Y
	real*8 twoPI
	real*8 dNHOR,dNVERT
	real*8 aat,aaf
	real*8 area
	real*8 dtet,dtet2,sdtet,cdtet
	integer NEQUAT,NVERT,I,NHOR,J,NU
       real*8 PI2,null,onehalf,one,two,three,four,big,small
       common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

	PI = PI2
	twoPI=two*PI

	area = null

        NEQUAT = DSQRT(dfloat(N ) * PI)

        NVERT = onehalf * NEQUAT
        IF (NVERT .LT. 2) NVERT = 2

        NU = 0

	dNVERT=dfloat(NVERT)

	dtet = PI/dNVERT
	dtet2= dtet*onehalf
	sdtet =two* dsin(dtet2)
	cdtet = dcos(dtet2)

        DO 100 I = 0,NVERT

	   FI = dtet*I

           Z = DCOS(FI)

           XY = DSIN(FI)

	   if(I.eq.0.or.I.eq.NVERT)then
	   aat=1.0d0 - cdtet
	   else
	   aat=XY*sdtet
	   end if

           NHOR = NEQUAT * XY

           IF (NHOR .LT. 1) NHOR = 1

	   dNHOR=dfloat(NHOR)

	   aaf=twoPI/dNHOR

           DO 50 J = 0,NHOR-1

              FJ = aaf*J

              X = DCOS(FJ) * XY

              Y = DSIN(FJ) * XY

              IF (NU .GE. N) GO TO 150

              NU = NU + 1

              U(1,NU) = X

              U(2,NU) = Y

              U(3,NU) = Z

	      AR(NU) = aaf*aat

	      area=area + AR(NU)

50         CONTINUE

100     CONTINUE

150     CONTINUE

        N = NU

        RETURN

        END
c
        SUBROUTINE GENUN01(RI,U,AR,TET,N)

	implicit none

	integer N, kanalz
	real*8 RI
        REAL*8 U(3,*)
	real*8 AR(*),TET(*)

	real*8 PI,FI,FJ,Z,XY,X,Y
	real*8 twoPI
	real*8 dNHOR,dNVERT
	real*8 aat,aaf,RI2
	real*8 area
	real*8 dtet,dtet2,sdtet,cdtet
	real*8 OPT_nhorkk

	integer NEQUAT,NVERT,I,NHOR,J,NU
	integer Nlarg
	integer NHORmin
	logical CONTROL

       real*8 PI2,null,onehalf,one,two,three,four,big,small
       common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

        data NHORmin/6/

	CONTROL = .false.
	OPT_nhorkk=1.3d0
	Nlarg=N*OPT_nhorkk
	OPT_nhorkk=0.5d0*OPT_nhorkk

	PI = PI2
	twoPI=two*PI
	RI2 = RI*RI

	area = null
        NEQUAT = DSQRT(dfloat(N ) * PI)
        NVERT = onehalf * NEQUAT
	dNVERT=dfloat(NVERT)

	I=nvert/2
	nvert = I + I
        IF (NVERT .LT. 2) NVERT = 2

	if(CONTROL)then
	write(*,*)'GENUN01: Nvert:',nvert
	endif

        NU = 0

	dtet = PI/dNVERT
	dtet2= dtet*onehalf
	sdtet =two* dsin(dtet2)
	cdtet = dcos(dtet2)

        DO 100 I = 0,NVERT

	   FI = dtet*I
           Z = RI*DCOS(FI)
           XY = DSIN(FI)

	   if(I.eq.0.or.I.eq.NVERT)then
	   aat=1.0d0 - cdtet
	   else
	   aat=XY*sdtet
	   end if

           NHOR = NEQUAT * XY
	   NHOR=int(OPT_nhorkk*NHOR + 0.5d0)*2
	   IF (NHOR .LT. NHORmin) NHOR = NHORmin
           if(I.eq.0.or.I.eq.NVERT) NHOR = 1

	if(CONTROL)then
	write(*,*)'GENUN01: iVert:',i,' Nhor:',NHOR
	endif

           dNHOR=dfloat(NHOR)
	   aaf=twoPI/dNHOR

	   XY = XY*RI

           DO 50 J = 0,NHOR-1

              FJ = aaf*J
              X = DCOS(FJ) * XY
              Y = DSIN(FJ) * XY

         IF (NU .GE. Nlarg)then
	 write(*,*)'GENUN01: probeSphe is not completed'
	 GO TO 150
	 end if

              NU = NU + 1
              U(1,NU) = X
              U(2,NU) = Y
              U(3,NU) = Z
	      AR(NU) = aaf*aat*RI2
	      TET(NU) = FI
	      area=area + AR(NU)

50      CONTINUE
100     CONTINUE
150     CONTINUE

        N = NU

        aaf = 4.0d0*PI*RI2/area
        open(unit = kanalz, file = 'teste.txt', status = 'unknown')
        do i=1,NU
C        AR(i) = aaf*AR(i)
C        TESTE = TESTE+AR(i)
        write(kanalz,'(i5,5f8.3)')i,U(1,i),U(2,i),U(3,i),AR(i),TET(i)
        end do
        close(kanalz)


	if(CONTROL)then
	write(*,'(a12,i5,a12,f16.14)')
     &	'GENUN01:N =',NU,' area norm=',aaf
	write(*,*)'dotN   x   y   z area TETA'
	do i=1,NU
c	write(*,'(i5,5f8.3)')i,U(1,i),U(2,i),U(3,i),AR(i),TET(i)!alterado
	enddo
	endif

	do i=1,NU
	AR(i) = aaf*AR(i)
	end do

        RETURN
        END
c
        SUBROUTINE GENUN02(RI,U,AR,NCHI,N,NC)

	implicit none
	integer N,NC, kanalz

	REAL*8 RI
        REAL*8 U(3,NC,N)
	real*8 AR(N)
	integer NCHI(N)

	real*8 PI,FI,FJ,Z,XY,X,Y
	real*8 twoPI
	real*8 dNHOR,dNVERT
	real*8 aat,aaf
	real*8 area,RI2
	real*8 dtet,dtet2,sdtet,cdtet

       real*8 dfich9(9),dtetch9(9)
       real*8 dfich09(9),dtetch09(9)
       real*8 dfich5(5),dtetch5(5)
       real*8 dfich05(5),dtetch05(5)
       real*8 dfich(9),dtetch(9)
       real*8 dfich0(9),dtetch0(9)
       real*8 chstep9,chstep5
       integer ncloc
       real*8 chstep,chsfi
       real*8 chstet,chstet1,chstet2
       real*8 aa,bb,tetch,fich
       real*8 aafc

       integer NEQUAT,NVERT,I,NHOR,J,NU
       integer NHORmin
       integer ich,ich1,k,k1,k2
       integer print,iprint0,iprint1

	logical CONTROL

       real*8 PI2,null,onehalf,one,two,three,four,big,small
       common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small
c
       data dfich9/0.0d0,1.0d0,0.0d0,-1.0d0,-1.0d0,-1.0d0,
     &       0.0d0,1.0d0,1.0d0/

       data dtetch9/0.0d0,-1.0d0,-1.0d0,-1.0d0,0.0d0,1.0d0,
     & 1.0d0,1.0d0,0.0d0/

       data dtetch09/0.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,
     & 1.0d0,1.0d0,1.0d0/

       data dfich09/0.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,
     & 6.0d0,7.0d0/

       data chstep9/0.33333d0/

c
       data dfich5/0.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0/
       data dtetch5/0.0d0,-1.0d0,-1.0d0,1.0d0,1.0d0/
       data dtetch05/0.0d0,1.0d0,1.0d0,1.0d0,1.0d0/
       data dfich05/0.0d0,1.0d0,3.0d0,5.0d0,7.0d0/
       data chstep5/0.250d0/

       data NHORmin/3/

       iprint0 = 0
       iprint1 = 1
       print  =  -1

       CONTROL = .false.

	if(print.ge.iprint0)then
	write(*,*)'GENUN02 start:'
	end if
c
c ASSIGN data
c
	ncloc = NC
	if(NC.eq.9)then
	do k=1,ncloc
        dfich(k)=dfich9(k)
	dtetch(k)=dtetch9(k)
	dfich0(k)=dfich09(k)
	dtetch0(k)=dtetch09(k)
	end do
	chstep = chstep9
	end if

	if(NC.eq.5)then
	do k=1,ncloc
        dfich(k)=dfich5(k)
	dtetch(k)=dtetch5(k)
	dfich0(k)=dfich05(k)
	dtetch0(k)=dtetch05(k)
	end do
	chstep = chstep5
	end if
c
	PI = PI2
        twoPI=two*PI
	RI2 = RI**2

	area = null

        NEQUAT = DSQRT(dfloat(N ) * PI)
        NVERT = onehalf * NEQUAT
        IF (NVERT .LT. 2) NVERT = 2

        NU = 0

	dNVERT=dfloat(NVERT)

	dtet = PI/dNVERT
	dtet2= dtet*onehalf
	chstet1 = chstep*dtet
	chstet2 = chstet1*two

	sdtet =two* dsin(dtet2)
	cdtet = dcos(dtet2)

        DO 100 I = 0,NVERT

	   FI = dtet*I

           Z = RI*DCOS(FI)

           XY = DSIN(FI)

	   if(I.eq.0.or.I.eq.NVERT)then
	   aat=1.0d0 - cdtet
	   else
	   aat=XY*sdtet
	   end if

           NHOR = NEQUAT * XY

           IF (NHOR .LT. NHORmin) NHOR = NHORmin
	   if(I.eq.0.or.I.eq.NVERT) NHOR = 1

	   dNHOR=dfloat(NHOR)

	   aaf=twoPI/dNHOR

	   if(NHOR.eq.1)then
	   chsfi = twoPI*0.125d0
	   chstet = dtet*0.66667d0
           else
	   chsfi = aaf*chstep
	   chstet = chstet1
	   end if
        open(unit = kanalz, file = 'genun2.txt', status = 'unknown')

	   XY = XY*RI

           DO 50 J = 0,NHOR-1

              FJ = aaf*J

              X = DCOS(FJ) * XY

              Y = DSIN(FJ) * XY

              IF (NU .GE. N)then
	      write(*,*)'GENUN02:ERROR: sphera is not completed'
	      GO TO 150
	      end if

              NU = NU + 1

	      NCHI(NU) = 1

              U(1,1,NU) = X

              U(2,1,NU) = Y

              U(3,1,NU) = Z
                        write(kanalz,'(i5,5f10.3)')i,U(1,1,NU),U(2,1,NU),U(3,1,NU)

	      AR(NU) = RI2*aaf*aat

	      area = area + AR(NU)

	      if(ncloc.ge.5)then
	      do ich = 2,ncloc
	      if(NHOR.eq.1)then
	      if(I.eq.0) tetch = FI + dtetch0(ich)*chstet
	      if(I.eq.NVERT) tetch = FI - dtetch0(ich)*chstet

	      fich  = dfich0(ich)*chsfi
	      aa = RI*dsin(tetch)
	      else
	      tetch = FI + dtetch(ich)*chstet
              aa = RI*dsin(tetch)

	      if(XY.gt.small)then
	      aafc=aa/XY
	      else
	      aafc=one
	      end if

	      fich  = FJ + dfich(ich)*chsfi*aafc
	      end if

	      U(1,ich,NU)=aa*dcos(fich)
	      U(2,ich,NU)=aa*dsin(fich)
	      U(3,ich,NU)=RI*dcos(tetch)
          write(kanalz,'(i5,5f10.3)')i,U(1,ich,NU),U(2,ich,NU)

	      NCHI(NU)=NCHI(NU) + 1
	      end do
	      end if


50      CONTINUE

100     CONTINUE

150     CONTINUE

        N = NU

        aaf = 4.0d0*PI*RI2/area
	if(CONTROL)then
	write(*,'(a12,i5,a12,f16.14)')
     &	'GENUN02:N =',NU,' area norm=',aaf
	end if

	do i=1,NU
	AR(i) = aaf*AR(i)
	end do


        close(kanalz)

	if(print.ge.iprint1)then
	write(*,*)'coordinates'
	do k=1,N
	write(*,'(a10,i5,3f10.6)')'parent ',k,(U(k1,1,k),k1=1,3)
	do k2=2,ncloc
	write(*,'(a10,i5,3f10.6)')'childr ',k2,(U(k1,k2,k),k1=1,3)
	end do
	end do
	end if
	if(print.ge.iprint0)then
	write(*,*)'GENUN02 finish: Nvect', N
	end if

        RETURN
        END

        SUBROUTINE GENUN03(RI,U,AR,N)
            implicit none
            integer N
            REAL*8 RI, valor, goldenratio, ang
            REAL*8 U(3,N)
            real*8 AR(N)
            real pi
            parameter (pi = 3.1415927)
            real i(N)
            real z(N)
            real theta(N)
            real phi(N)
            real out(3,N)
            logical out_xyz
            integer k
            real phi_1, phi_2
            real area_dot
            integer kanalz
            out_xyz = .true.
            area_dot = (4*pi*RI**2)/N

               if(N<1) then
                stop
               end if
               valor = 5
               goldenratio = (1+sqrt(valor))/2
               do 100 k = 1,N
                    i(k) = k - 0.5
100         continue
            do 101 k = 1,N
                z(k) = 1-(2*i(k))/N
101          continue
            do 102 k = 1,N
C            ang = MAX(-1.0, MIN(1.0,z(k)))
                if(z(k)<-1) then
                 ang = -1
                else if(z(k)>1) then
                 ang = 1
                else
                    ang = z(k)
                end if
                theta(k) = acos(ang)
102         continue
            do 103 k = 1,N
                phi_1 = (2*pi*i(k))/goldenratio
                phi_2 = 2*pi
                phi(k) = mod(phi_1, phi_2)
103         continue

            if (out_xyz) then
                open(unit=kanalz,file="genun03.ms",status='unknown')
                do 104 k = 1,N
                    U(1,k) = RI*sin(theta(k))*cos(phi(k))
                    U(2,k) = RI*sin(theta(k))*sin(phi(k))
                    U(3,k) = RI*z(k)
                    write(kanalz,'(i3, 3f8.3)')k,U(1,k),U(2,k),U(3,k)
104             continue
           end if
           continue

           do 105 k = 1,N
                AR(k) = area_dot
105         continue

             close(kanalz)


        RETURN
        END


c------------------------------------------------------------------
c          sims:  subroutines
c
c------------------------------------------------------------------
	logical function compare_ij(i,j,n,m,l)

	implicit none
	integer i,j,n,m,l

	compare_ij=.false.

	if(.not.(i.eq.n.or.i.eq.m.or.i.eq.l))then
	return
	end if
	if(.not.(j.eq.n.or.j.eq.m.or.j.eq.l))then
	return
	end if

	compare_ij=.true.
	return
	end

	logical function compare_ijk(i,j,k,n,m,l)

	implicit none
	integer i,j,k,n,m,l

	compare_ijk=.false.

	if(.not.(n.eq.i.or.n.eq.j.or.n.eq.k))then
	return
	end if
	if(.not.(m.eq.i.or.m.eq.j.or.m.eq.k))then
	return
	end if
	if(.not.(l.eq.i.or.l.eq.j.or.l.eq.k))then
	return
	end if

	compare_ijk=.true.
	return
	end

c
        real*8 function triang_area(a,b,c)
        implicit none
        real*8 a(3),b(3),c(3)
        real*8 aa,bb,ab,ak,bk,ar
        integer k
        aa=0.0D0
        bb=0.0D0
        ab=0.0D0
        do k=1,3
          ak=b(k)-a(k)
          bk=c(k)-a(k)
          aa=aa+ak*ak
          bb=bb+bk*bk
          ab=ab+ak*bk
        enddo

	triang_area=0.0d0

	ar=aa*bb-ab*ab
	if(ar.ge.0.0d0)triang_area=0.5D0*dsqrt(ar)
        return
        end
c------------------------------------------------------------------------
        real*8 function curv_triang_area(a,b,c,center,r)
c
c
        implicit none
        real*8 a(3),b(3),c(3)
	real*8 center(3),r
        real*8 tol /1.0D-4/, one /1.0D0/
        real*8 triang_area
        real*8 aa,bb,cc,ak,bk,cost12,cost13,cost23,r2,
     &   sint12,sint13,sint23,cosa,cosb,cosg,sina,sinb,sing,cossurf
        integer k

c       print *,'a:',a
c       print *,'b:',b
c       print *,'c:',c

        cost12=0.0D0
        cost13=0.0D0
        cost23=0.0D0
        r2=r*r
        do k=1,3
          aa=a(k)-center(k)
          bb=b(k)-center(k)
          cc=c(k)-center(k)
          cost12=cost12+aa*bb
          cost13=cost13+aa*cc
          cost23=cost23+bb*cc
        enddo
        cost12=cost12/r2
        cost13=cost13/r2
        cost23=cost23/r2
        sint12=1.0D0-cost12*cost12
        sint13=1.0D0-cost13*cost13
        sint23=1.0D0-cost23*cost23
        if (sint12.lt.tol .or. sint13.lt.tol .or. sint23.lt.tol) then
          curv_triang_area=triang_area(a,b,c)
          return
        endif
        sint12=dsqrt(sint12)
        sint13=dsqrt(sint13)
        sint23=dsqrt(sint23)
c        print *,'cost12=',cost12,' cost13=',cost13,' cost23=',cost23
        cosa=(cost23-cost12*cost13)/(sint12*sint13)
        cosb=(cost13-cost12*cost23)/(sint12*sint23)
        cosg=(cost12-cost13*cost23)/(sint13*sint23)
        sina=1.0D0-cosa*cosa
        sinb=1.0D0-cosb*cosb
        sing=1.0D0-cosg*cosg
        if (sina.lt.tol .or. sinb.lt.tol .or. sing.lt.tol) then
          curv_triang_area=triang_area(a,b,c)
          return
        endif
        sina=dsqrt(sina)
        sinb=dsqrt(sinb)
        sing=dsqrt(sing)

c        print *,'cosa=',cosa,' cosb=',cosb,' cosg=',cosg
c        print *,'cos(area/R**2)=',-cosa*cosb*cosg+sina*sinb*cosg
c     &                          +sina*cosb*sing+cosa*sinb*sing

        cossurf=-cosa*cosb*cosg+sina*sinb*cosg
     &                          +sina*cosb*sing+cosa*sinb*sing
        if (dabs(cossurf).gt.1.0D0) then
          curv_triang_area=0.0D0
        else
          curv_triang_area=r2*dacos(cossurf)
        endif
        return
        end

c

	 subroutine gener_conc_dot01(v1,v2,v3,rpr,dens,
     &	 area_cc,vdot,area_dot,vndot,ndot,dot01_run)

         implicit none
         integer ndotccmx
         parameter(ndotccmx=128)
         integer ndot
         real*8 v1(3),v2(3),v3(3),rpr,dens
         real*8 vdot(3,ndotccmx),area_dot(ndotccmx),vndot(3,ndotccmx)
         logical dot01_run

         real*8 area_cc,center(3)
	 real*8 aa,bb,cc,r2,empk,aa1
         real*8 cost(3),tedgemx,costm,dens_l
         real*8 curv_triang_area

	 logical control

         integer tmax,i,k,j,klev,nd,klevmax,nd_l
         parameter(tmax=4)
         integer ntriangl(tmax)
         data ntriangl/1,4,16,64/
	 data center/0.0d0,0.0d0,0.0d0/

	 control = .false.

         dot01_run=.true.
         klevmax=4
	 area_cc=curv_triang_area(v1,v2,v3,center,rpr)

	nd= dint(area_cc*dens + 0.5d0)

        empk = 2.00d0
        r2=rpr*rpr
        do k=1,3
           cost(k)=0.0d0
        end do

        do k=1,3
          aa=v1(k)-center(k)
          bb=v2(k)-center(k)
          cc=v3(k)-center(k)
          cost(1)=cost(1)+aa*bb
          cost(2)=cost(2)+aa*cc
          cost(3)=cost(3)+bb*cc
        enddo

        costm=1.0d6
        do k=1,3
        if(costm.gt.cost(k))costm=cost(k)
        end do

        tedgemx=rpr*dacos(costm/r2)
        dens_l = empk*dsqrt(1.25d0/dens)
        nd_l = dint(tedgemx/dens_l)
        if(nd_l.lt.1)nd_l=1
        nd_l = 4**(nd_l - 1)

        if(nd_l.gt.nd)nd=nd_l


	klev = 1
	if(nd.le.ntriangl(2)/2)klev=1
	if(nd.gt.ntriangl(2)/2.and.nd.le.ntriangl(3)/2)klev=2
	if(nd.gt.ntriangl(3)/2.and.nd.le.ntriangl(4)/2)klev=3
	if(nd.gt.ntriangl(4)/2)klev=4
        if(klev.gt.klevmax)klev=klevmax

cmar19        if(nd.gt.1.25*ntriangl(4))then
        if(nd.gt.ntriangl(4))then
	ndot=0
        dot01_run = .false.
        goto 2000
        end if


        if(klev.eq.1)then
	call gener_cc1(v1,v2,v3,rpr,area_cc,ndot,vdot,area_dot,vndot)
        goto 1010
	end if

	if(klev.eq.2)then
	call gener_cc2(v1,v2,v3,rpr,area_cc,ndot,vdot,area_dot,vndot)
	goto 1010
	end if

	if(klev.eq.3)then
       call gener_cc3(v1,v2,v3,rpr,area_cc,ndot,vdot,area_dot,vndot)
       goto 1010
       end if

       if(klev.eq.4)then
       call gener_cc4(v1,v2,v3,rpr,area_cc,ndot,vdot,area_dot,vndot)
       goto 1010
       end if

1010   continue

       if(control)then
       aa=0.0d0
       bb=0.0d0
       aa1=area_cc/dfloat(ndot)
       do i=1,ndot
       aa=aa + area_dot(i)
       bb=bb + (area_dot(i)-aa1)**2
       end do
       bb=dsqrt(bb/dfloat(ndot))
       bb=bb/aa1
       write(*,*)'gener_cc_dot:ndot=',ndot
       write(*,*)'Area RMS for Dots bb(relative)=',bb
       write(*,*)'AREA_CC control: area_cc, sum=',area_cc,aa
       end if

2000       return
       end

c
       subroutine gener_cc1(v1,v2,v3,rpr,area,ndot,vdot,area_dot,vndot)

       implicit none
         integer ndotccmx
         parameter(ndotccmx=128)
         integer ndot
         real*8 v1(3),v2(3),v3(3),rpr,area
         real*8 vdot(3,ndotccmx),area_dot(ndotccmx),vndot(3,ndotccmx)

       real*8 s,rs
       integer i,k

       ndot=1
       s=0.0d0
       do k=1,3
       vdot(k,ndot)=v1(k)+v2(k)+v3(k)
       s=s + vdot(k,ndot)**2
       end do

       s=dsqrt(s)
       do k=1,3
       vndot(k,ndot)=-vdot(k,ndot)/s
       vdot(k,ndot)=-vndot(k,ndot)*rpr
       end do

       area_dot(ndot)=area

       return
       end
c
       subroutine gener_cc2(v1,v2,v3,rpr,area,ndot,vdot,area_dot,vndot)


       implicit none
         integer ndotccmx
         parameter(ndotccmx=128)
         integer ndot
         real*8 v1(3),v2(3),v3(3),rpr,area
         real*8 vdot(3,ndotccmx),area_dot(ndotccmx),vndot(3,ndotccmx)

c
       integer ntrimx
       parameter (ntrimx=4)
       integer nvertmx
       parameter (nvertmx=3+3)

c real function:
       real*8 curv_triang_area

       real*8 one3,s,rs,ar,vertex(3,nvertmx)
       real*8 center(3)
       integer vertex_parent(2,nvertmx)
       integer tri_parent(3,ntrimx)
       integer i,j,k,l
       logical OPT_AreaEachDot,OPT_norm
       data vertex_parent/1,1, 2,2, 3,3, 1,2, 2,3, 1,3/
       data tri_parent/4,5,6, 1,4,6, 3,5,6, 2,4,5/
       data center/0.0d0,0.0d0,0.0d0/

       OPT_AreaEachDot=.true.
       OPT_norm = .true.

       ndot=ntrimx
       ar=area/dfloat(ndot)
       one3=1.0d0/3.0d0

       do l=1,3
       vertex(l,1)=v1(l)
       vertex(l,2)=v2(l)
       vertex(l,3)=v3(l)
       end do

       do i=4,nvertmx
       s=0.0d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       if(OPT_norm)then
       s=s+vertex(l,i)**2
       end if
       end do
       if(OPT_norm)then
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end if
       end do

       do i=1,ntrimx
       s=0.0d0
       do k=1,3
       vdot(k,i)=vertex(k,tri_parent(1,i))+
     &           vertex(k,tri_parent(2,i))+
     &           vertex(k,tri_parent(3,i))
       s=s+vdot(k,i)**2
       end do

       s=dsqrt(s)
       do k=1,3
       vndot(k,ndot)= -vdot(k,ndot)/s
       vdot(k,ndot)=-vndot(k,ndot)*rpr
       end do



       if(.not.OPT_AreaEachDot)then
       area_dot(i)=ar
       else

      area_dot(i)=curv_triang_area(vertex(1,tri_parent(1,i)),
     &	                         vertex(1,tri_parent(2,i)),
     &                           vertex(1,tri_parent(3,i)),
     &                           center,rpr)

       end if

       end do

c
       return
       end
c

       subroutine gener_cc3(v1,v2,v3,rpr,area,ndot,vdot,area_dot,vndot)


       implicit none
         integer ndotccmx
         parameter(ndotccmx=128)
         integer ndot
         real*8 v1(3),v2(3),v3(3),rpr,area
         real*8 vdot(3,ndotccmx),area_dot(ndotccmx),vndot(3,ndotccmx)
       real*8 curv_triang_area

       integer  ntrimx
       parameter (ntrimx=16)
       integer nvertmx
       parameter (nvertmx=3+12)

       real*8 one3,s,rs,ar,vertex(3,nvertmx)
       real*8 center(3)
       integer vertex_parent(2,nvertmx)
       integer tri_parent(3,ntrimx)
       integer nvert1,nvert2,nvert3
       integer i,j,k,l
       logical OPT_AreaEachDot,OPT_norm

       data center/0.0d0,0.0d0,0.0d0/
       data vertex_parent/1,1, 2,2, 3,3,
     & 1,2, 2,3, 1,3,
     & 1,4, 1,6, 4,6, 4,5, 5,6, 3,5, 3,6, 2,4, 2,5/

       data tri_parent/9,10,11, 6,9,11, 4,9,10, 5,10,11,
     & 1,7,8, 7,8,9, 6,8,9, 4,7,9,
     & 3,12,13, 11,12,13, 5,11,12, 6,11,13,
     & 10,14,15, 4,10,14, 5,10,15, 2,14,15 /

       OPT_AreaEachDot=.true.
       OPT_norm = .true.

       ndot=ntrimx
       ar=area/dfloat(ndot)
       one3=1.0d0/3.0d0

       nvert1=3
       do l=1,3
       vertex(l,1)=v1(l)
       vertex(l,2)=v2(l)
       vertex(l,3)=v3(l)
       end do

       nvert2=6
       do i=nvert1+1,nvert2
       s=0.0d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       if(OPT_norm)then
       s=s+vertex(l,i)**2
       end if
       end do

       if(OPT_norm)then
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end if
       end do

       nvert3=nvertmx
       do i=nvert2+1,nvert3
       s=0.0d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       if(OPT_norm)then
       s=s+vertex(l,i)**2
       end if
       end do
        if(OPT_norm)then
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end if
       end do

       do i=1,ntrimx
       s=0.0d0
       do k=1,3
       vdot(k,i)=vertex(k,tri_parent(1,i))+
     &           vertex(k,tri_parent(2,i))+
     &           vertex(k,tri_parent(3,i))
       s=s+vdot(k,i)**2
       end do

       s=dsqrt(s)
       do k=1,3
       vndot(k,ndot)= -vdot(k,ndot)/s
       vdot(k,ndot)=-vndot(k,ndot)*rpr
       end do

       if(.not.OPT_AreaEachDot)then
       area_dot(i)=ar
       else

      area_dot(i)=curv_triang_area(vertex(1,tri_parent(1,i)),
     &	                           vertex(1,tri_parent(2,i)),
     &                             vertex(1,tri_parent(3,i)),
     &                             center,rpr)

       end if
       end do

c
       return
       end
c
       subroutine gener_cc4(v1,v2,v3,rpr,area,ndot,vdot,area_dot,vndot)

       implicit none
       integer ndotccmx
       parameter(ndotccmx=128)
       integer ndot
       real*8 v1(3),v2(3),v3(3),rpr,area
       real*8 vdot(3,ndotccmx),area_dot(ndotccmx),vndot(3,ndotccmx)
       real*8 curv_triang_area

       integer ntrimx
       parameter (ntrimx=64)
       integer nvertmx
       parameter (nvertmx=3+42)

       real*8 one3,s,rs,ar,vertex(3,nvertmx)
       real*8 center(3)
       integer vertex_parent(2,nvertmx)
       integer nvert1,nvert2,nvert3,nvert4
       integer tri_parent(3,ntrimx)
       integer i,j,k,l
       logical OPT_AreaEachDot

       data center/0.0d0,0.0d0,0.0d0/
       data vertex_parent/1,1,2,2,3,3,
     & 1,2,2,3,1,3,
     & 1,4,1,6,4,6,4,5,5,6,3,5,3,6,2,4,2,5,
     & 9,10, 10,11, 9,11, 6,9, 6,11,
     & 4,9, 4,10, 5,11, 5,10, 6,8,
     & 8,9, 7,9, 4,7, 7,8, 1,8,
     & 1,7, 6,13, 11,13, 12,13, 11,12,
     & 5,12, 3,12, 3,13, 4,14, 10,14,
     & 10,15, 14,15, 5,15,2,14, 2,15 /

       data tri_parent/16,17,18, 9,16,18, 11,17,18, 10,16,17,
     &                 18,19,20, 9,18,19, 11,18,20, 6,19,20,
     &                 21,22,16, 9,16,21, 10,16,22, 4,21,22,
     &                 17,23,24, 10,17,24, 11,17,23, 5,23,24,

     &                 19,25,26, 6,19,25, 9,19,26, 8,25,26,
     &                 9,21,27, 21,27,28, 4,21,28, 7,27,28,
     &                 9,26,27, 8,26,29, 26,27,29, 7,27,29,
     &                 8,29,30, 7,29,31, 29,30,31, 1,30,31,

     &                 6,20,32, 20,32,33, 13,32,33, 11,20,33,
     &                 13,33,34, 11,33,35, 33,34,35, 12,34,35,
     &                 11,23,35, 23,35,36, 12,35,36, 5,23,36,
     &                 13,34,38, 34,37,38, 12,34,37, 3,37,38,

     &                 22,39,40, 10,22,40, 4,22,39, 14,39,40,
     &                 10,40,41, 14,40,42, 15,41,42, 40,41,42,
     &                 10,24,41, 5,24,43,  24,41,43, 15,41,43,
     &                 42,44,45, 14,42,44, 15,42,45, 2,44,45/



       OPT_AreaEachDot=.true.

       nvert1=3
       nvert2=6
       nvert3=15
       nvert4=nvertmx

       ndot=ntrimx
       ar=area/dfloat(ndot)
       one3=1.0d0/3.0d0

       do l=1,nvert1
       vertex(l,1)=v1(l)
       vertex(l,2)=v2(l)
       vertex(l,3)=v3(l)
       end do

       do i=nvert1+1,nvert2
       s=0.d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       s=s+vertex(l,i)**2
       end do
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end do


       do i=nvert2+1,nvert3
       s=0.d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       s=s+vertex(l,i)**2
       end do
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end do



       do i=nvert3+1,nvert4
       s=0.d0
       do l=1,3
       vertex(l,i)=(vertex(l,vertex_parent(1,i))+
     &              vertex(l,vertex_parent(2,i)))*0.5d0
       s=s+vertex(l,i)**2
       end do
       rs=rpr/dsqrt(s)
       do l=1,3
       vertex(l,i)=vertex(l,i)*rs
       end do
       end do

       do i=1,ntrimx
       s=0.0d0
       do k=1,3
       vdot(k,i)=vertex(k,tri_parent(1,i))+
     &           vertex(k,tri_parent(2,i))+
     &           vertex(k,tri_parent(3,i))
       s=s+vdot(k,i)**2
       end do

       s=dsqrt(s)
       do k=1,3
       vndot(k,ndot)= -vdot(k,ndot)/s
       vdot(k,ndot)=-vndot(k,ndot)*rpr
       end do

       if(.not.OPT_AreaEachDot)then
       area_dot(i)=ar
       else

      area_dot(i)=curv_triang_area(vertex(1,tri_parent(1,i)),
     &	                           vertex(1,tri_parent(2,i)),
     &                             vertex(1,tri_parent(3,i)),
     &                             center,rpr)

       end if

       end do

c
       return
       end
c

	 subroutine gener_conc_dot02(UP,ARUP,TETP,PrPIJK,
     &	 v1,v2,v3,rpr,dens,
     &	 area_cc,vdot,area_dot,vndot,NUP,ndot)

         implicit none
	 real*8 UP(3,*),ARUP(*),TETP(*)
	 real*8 PrPIJK(3)
	 integer NUP
         integer ndot
         real*8 v1(3),v2(3),v3(3),rpr,dens
         real*8 vdot(3,*),area_dot(*),vndot(3,*)
	 real*8 area_cc

c local variables
	 real*8 vv(3,4),sv(3),vt(3,4)
	 real*8 G(3,3),x(3),y(3),z(3),u(3)
         real*8 v11(3),v21(3),v31(3),z1(3)
	 real*8 w12(3),w13(3),w23(3)
         real*8 center(3)
	 real*8 aa,bb,cc
	 real*8 s,s1,s2,s3,tetmx
	 real*8 s11,s21,s31
	 real*8 null,sign,big
	 real*8 area
         real*8 toler_a,toler_b,toler_c(3)
c real function:
         real*8 curv_triang_area
	 real*8 DOT

c
	 logical OPT_equator, OPT_polar
	 logical CONTROL
         integer i,k,j,l,lmin
	 data center/0.0d0,0.0d0,0.0d0/

	 OPT_polar=.true.
	 OPT_equator=.not.OPT_polar
	 CONTROL = .false.

	 null = 0.0d0
	 big = 1.0d10
         toler_a = 0.125d0/dsqrt(dens)
	 toler_b = 0.995d0
	 do k=1,3
	 toler_c(k)=0.0d0
         end do
	 if(CONTROL)then
	 write(*,*)'gener_CC02'
	 write(*,'(a6,3f10.6)')'v1=',v1
	 write(*,'(a6,3f10.6)')'v2=',v2
	 write(*,'(a6,3f10.6)')'v3=',v3
	 end if

	  do k=1,3
	  vv(k,1)=v1(k)
	  vv(k,2)=v2(k)
	  vv(k,3)=v3(k)
	  sv(k)=0.0d0
	  end do

        s=0.0d0
	s1=0.0d0
	s2=0.0d0
	s3=0.0d0
        do k=1,3
	  z(k)=v1(k)+v2(k)+v3(k)
	  s=s + z(k)**2
	  w12(k)=v2(k) - v1(k)
	  w13(k)=v3(k) - v1(k)
	  w23(k)=v3(k) - v2(k)
	  s1=s1+w12(k)**2
	  s2=s2+w13(k)**2
	  s3=s3+w23(k)**2
        end do

	s = dsqrt(s)
	s1= dsqrt(s1)
	s2= dsqrt(s2)
	s3= dsqrt(s3)
	do k=1,3
	z(k) = z(k)/s
        z1(k)=z(k)*rpr
	vv(k,4)=z1(k)
	w12(k)=w12(k)/s1
	w13(k)=w13(k)/s2
	w23(k)=w23(k)/s3
	end do

	tetmx=big
	lmin=1
	do l=1,3
	do k=1,3
        sv(l) = sv(l) + z(k)*vv(k,l)
	end do
	if(tetmx.gt.sv(l))then
	tetmx=sv(l)
	lmin=l
	end if
	end do
        tetmx = tetmx/rpr
	tetmx = dacos(tetmx)

        if(CONTROL)then
	write(*,*)'GENUN01:mostRemoteVert:',lmin,' tetmx = ', tetmx
	end if

	s=null
	do k=1,3
	x(k) = vv(k,lmin) - z(k)*sv(lmin)
	s = s + x(k)**2
	end do
	s = dsqrt(s)
	do k=1,3
	x(k) = x(k)/s
	end do

	call CROSS(z,x,y)

	if(OPT_polar)then
	do k=1,3
	G(k,1) = x(k)
	G(k,2) = y(k)
	G(k,3) = z(k)
	end do
	end if
	if(OPT_equator)then
	do k=1,3
	G(k,3) = y(k)
	G(k,2) = x(k)
	G(k,1) = z(k)
	end do
	end if

	if(CONTROL)then
	write(*,*)'Rotatiom Matr G:'
	 write(*,'(a6,3f10.6)')'G(k,1)=',(G(k,1), k=1,3)
	 write(*,'(a6,3f10.6)')'G(k,2)=',(G(k,2), k=1,3)
	 write(*,'(a6,3f10.6)')'G(k,3)=',(G(k,3), k=1,3)
        do l=1,4
        do j=1,3
        u(j)=0.0
        do k=1,3
        u(j)=u(j)+G(k,j)*vv(k,l)
        end do
        end do
        do k=1,3
        vt(k,l)=u(k)
        end do
        end do
	write(*,*)'Rotated vertexes: toLocal x,y,z'
	 write(*,'(a9,3f10.6)')'v1(k)=',(vt(k,1), k=1,3)
	 write(*,'(a9,3f10.6)')'v2(k)=',(vt(k,2), k=1,3)
	 write(*,'(a9,3f10.6)')'v3(k)=',(vt(k,3), k=1,3)
	 write(*,'(a9,3f10.6)')'z1(k)=',(vt(k,4), k=1,3)
	end if

	s11=0.0d0
	s21=0.0d0
	s31=0.0d0
        do k=1,3
        v11(k) = z1(k)-v1(k)
        v21(k) = z1(k)-v2(k)
        v31(k) = z1(k)-v3(k)
	s11=s11 + v11(k)**2
        s21=s21 + v21(k)**2
	s31=s31 + v31(k)**2
        end do

	s11=dsqrt(s11)
	s21=dsqrt(s21)
	s31=dsqrt(s31)

	do k=1,3
	v11(k)=v11(k)/s11
	v21(k)=v21(k)/s21
	v31(k)=v31(k)/s31
	end do

	s1=0.0d0
	s2=0.0d0
	s3=0.0d0
        do k=1,3
	s1=s1 + w12(k)*v11(k)
	s2=s2 + w12(k)*v21(k)
	s3=s3 + w13(k)*v31(k)
	end do
        s1 = dsqrt(1.0d0 - s1**2)
	    s2 = dsqrt(1.0d0 - s2**2)
	    s3 = dsqrt(1.0d0 - s3**2)

	if(s1.gt.0.0d0)then
	toler_c(1)=toler_a/s1
	if(s11.lt.toler_c(1))then
	toler_c(1)=toler_b*s11
	end if
	else
	toler_c(1)=toler_b*s11
        end if

	if(s2.gt.0.0d0)then
	toler_c(2)=toler_a/s2
	if(s21.lt.toler_c(2))then
	toler_c(2)=toler_b*s21
	end if
	else
	toler_c(2)=toler_b*s21
	end if

	if(s3.gt.0.0d0)then
	toler_c(3)=toler_a/s3
	if(s31.lt.toler_c(3))then
	toler_c(3)=toler_b*s31
	end if
	else
	toler_c(3)=toler_b*s31
	end if

        do k=1,3
        v11(k) = v1(k) + v11(k)*toler_c(1)
        v21(k) = v2(k) + v21(k)*toler_c(2)
        v31(k) = v3(k) + v31(k)*toler_c(3)
        end do
	if(CONTROL)then
	do k=1,3
	vv(k,1)=v11(k)
	vv(k,2)=v21(k)
	vv(k,3)=v31(k)
	end do

        do l=1,4
        do j=1,3
        u(j)=0.0
        do k=1,3
        u(j)=u(j)+G(k,j)*vv(k,l)
        end do
        end do
        do k=1,3
        vt(k,l)=u(k)
        end do
        end do
	write(*,*)'Rotated shiftedVertexes: toLocal x,y,z'
	 write(*,'(a9,3f10.6)')'v1(k)=',(vt(k,1), k=1,3)
	 write(*,'(a9,3f10.6)')'v2(k)=',(vt(k,2), k=1,3)
	 write(*,'(a9,3f10.6)')'v3(k)=',(vt(k,3), k=1,3)
	 write(*,'(a9,3f10.6)')'z1(k)=',(vt(k,4), k=1,3)
	end if

	call CROSS(v21,v31,x)
	call CROSS(v31,v11,y)
	call CROSS(v11,v21,z)

        sign = DOT(v11,x)

	ndot = 0
	area = null

	do i=1,NUP

	do k=1,3
	u(k)=null
	do l=1,3
	u(k)=u(k) + G(k,l)*UP(l,i)
	end do
	end do
        s1=DOT(u,x)*sign
        if(s1.gt.null)then
		 s2=DOT(u,y)*sign
		 if(s2.gt.null)then
		  s3=DOT(u,z)*sign
		   if(s3.gt.null)then
c
	 ndot = ndot +1
	 do k=1,3
	 vndot(k,ndot) = -u(k)/rpr
	 vdot(k,ndot) = u(k)
	 end do

	 area_dot(ndot)= ARUP(i)
	 area=area + area_dot(ndot)
	 if(CONTROL)then
	 write(*,*)'dot accepted i=',i,'tet(i)=',TETP(i)
	 write(*,*)'UPprsph, xyz                    area'
	write(*,'(i5,1x,3f10.6,2x,30x,2x,f10.7)')i,
     &	(UP(k,i),k=1,3),ARUP(i)
	write(*,*)'dotRot:ndot, xyz            vnorm   area'
	write(*,'(i5,1x,3f10.6,2x,3f10.6,2x,f10.7)')ndot,
     &	(vdot(k,ndot),k=1,3),(vndot(k,ndot),k=1,3),area_dot(ndot)
	 end if

	 end if
	 end if
	 end if


         end do


1001     continue

	 if(ndot.eq.0)then
	 write(*,*)'gener_conc02:ERROR: NO dots on CC face'
	 stop
	 end if

	 area_cc=curv_triang_area(v1,v2,v3,center,rpr)
	 if(CONTROL)then
	 write(*,*)'Analcurv_triang_area=',area_cc
	 write(*,*)'Num      triang_area=',area
	 write(*,*)'ndot =',ndot
	 end if

	 if(area.ne.null)then
         aa = area_cc/area
	 do i=1,ndot
         area_dot(i) = area_dot(i)*aa
	 end do
	 end if

	if(CONTROL)then
	write(*,*)'CCndot:',ndot,' ReNormAreaCoeff: aa=',aa
	write(*,*)'dot i, xyz            vnorm   area'
	do i=1,ndot
	write(*,'(i5,1x,3f10.6,2x,3f10.6,2x,f10.7)')i,
     &	(vdot(k,i),k=1,3),(vndot(k,i),k=1,3),area_dot(i)
	end do
	write(*,*)'end dotcc02 ***********'
	end if

	 return
	 end
       subroutine arc_points02(IATOM,JATOM,
     & ri,rj,dij,hij,rpr,dens,rad_sms,
     & putsm,putsm_smp,narc,ni,nj,
     & v0arc,stor,arc_cross,free_tor,
     & nijs_face,RotAng,
     & dot_xyz,dot_vn,dot_area,dot_type,dot_atom,ndot)


       implicit none
       integer narcmx,nrotmx,nbr_mx,ndotmx

       parameter(narcmx=100)
       parameter(nrotmx=50 )
       parameter(ndotmx=narcmx*nrotmx)
       parameter(nbr_mx = 180)

       integer putsm

       integer putsm_smp

       integer IATOM,JATOM,nijs_face
       integer narc,ni,nj
       real*8 RotAng(nbr_mx)
       real*8 v0arc(3,narcmx)
       real*8 ri,rj,dij,hij,rpr,dens
       real*8 rad_sms,rad_sm
       real*8 stor
       logical arc_cross,free_tor

       integer ndot
       real*8 dot_xyz(3,ndotmx),dot_vn(3,ndotmx),dot_area(ndotmx)
       integer dot_type(ndotmx),dot_atom(ndotmx)

       real*8 PI2,null,onehalf,one,two,three,four,big,small
       common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

       real*8 teti,tetj,dtet2
       real*8 dens2,dens2_p,dens2_r,dens2_sm
       real*8 di,dj,di2,dj2,tet,tetmx,dtet
       real*8 pi,cost,cosdt,rip,rjp,aa
       real*8 step,bb,rpr_s
       real*8 dfi,fi_rot,start_ang
       real*8 ditet_i,ditet_j
       real*8 sinditet_i,sinditet_j
       real*8 rhditet_i,rhditet_j
       real*8 sinditet_s,ditet_s
       real*8 d_arpi,d_arpj,d_arps
       real*8 costi,costj,tetai,tetaj,rpsm,costsm
       real*8 tetasm,sintsm,pi05
       real*8 psi_xyz(3),psj_xyz(3)
       real*8 stor_i,stor_j,smsph_i,stor_ijs
       real*8 area_i,area_j, area_si, area_sj
       real*8 arc_tet(narcmx)
       real*8 rot_ang,dsint,dcost,dsint2
       real*8 drot_max
       real*8 rpr2,rpr22,rad_sm2,rprh
       real*8 toler_crossloc

       integer arc_type(narcmx)
       integer sadl_ptype, sadl_smtype

       logical changeSg
       logical testpr
       logical one_side
       integer ni_s,nj_s,ndotl
       integer i,j,k,n
       integer n1,n2,n3
       integer isf,nsf,nr

c        testpr=.true.
       testpr=.false.
       one_side=.false.

        pi=PI2

        sadl_ptype = 2
        sadl_smtype = 4

	toler_crossloc = 0.0d0
	if(free_tor)toler_crossloc = toler_cross

        drot_max =  pi/1.5d0

        dens2=dsqrt(dens)
        dens2_p=rpr*dens2
        rpr2 = rpr*rpr
        rpr22=rpr2*2.0d0
        rprh = rpr*hij
cct
        rad_sm = rad_sms
        if(rad_sm.eq.0.0d0)rad_sm = 0.1d-9
        rad_sm2=rad_sm*rad_sm

	rip=ri+rpr
	rjp=rj+rpr
	stor=null


	if(ri.gt.rj)then
	cost=hij/rip
	changeSg=.false.
	else
	cost=hij/rjp
	changeSg=.true.
	end if

	tetmx=dacos(cost)

	cosdt=(rip**2 + rjp**2 - dij**2)/(2.0d0*rip*rjp)
	dtet=dacos(cosdt)

	if(testpr)then
	write(*,*)'arc_points:ATOMIJ:',IATOM,JATOM
	write(*,*)'arc_points: tetmx, dtet=',tetmx,dtet
        end if
         arc_cross=.false.
         if(hij.le.rpr + toler_crossloc)arc_cross=.true.

	 if(tetmx.ge.dtet)one_side=.true.
	 if(one_side)arc_cross=.false.
ctt
        if(testpr)then
        write(*,*)'arc_points: hij,toler_cross:',
     &  hij,toler_cross
	 write(*,*)'one_side:',one_side,' arc_cross:',arc_cross
        end if

	dtet2=dtet*0.5d0
	ni=0
	nj=0
	if(putsm.eq.2)then
	ni=DINT(dtet2*dens2_p+0.5d0)
	if(ni.eq.0)ni=1
	narc=2*ni
	nj=ni

	else

	narc=DINT(dtet*dens2_p+0.5d0)
        if(narc.eq.0.and.putsm.eq.1)narc=1
	nj=narc/2
	ni=narc-nj
	if(changeSg)then
	ni=narc/2
	nj=narc-ni
	end if

	end if

	if(narc.lt.1)goto 1001

         if(narc.gt.narcmx)then
         write(*,*)'ERROR : sub: nrc_points , narcmx is small=',narcmx
	 stop
         end if

	di=dtet/dfloat(narc)
	di2=di*0.5d0

        ditet_i=di2
        ditet_j=di2
        rhditet_i=rpr*hij*2.0d0*ditet_i
        rhditet_j=rhditet_i
        sinditet_i=dsin(ditet_i)
        sinditet_j=sinditet_i
        ditet_s=0.0d0

	tet=tetmx+di2

	if(changeSg)then
	tet=tetmx - dtet
	tet=tet-di2
	di=-di
	end if

        if(.not.arc_cross)then
        rpr_s = rpr
        if(changeSg)rpr_s = -rpr
	do i=1,narc
	tet=tet - di
	v0arc(1,i)= -rpr_s*dsin(tet)
	v0arc(2,i)= hij - rpr*dcos(tet)
        v0arc(3,i)= 0.0d0
        if(i.le.ni)then
        arc_type(i)=1
        else
        arc_type(i)=2
        end if
        arc_tet(i) = tet
	end do

        dsint2=dsin(tetmx-di*ni)
        stor_i=rhditet_i*ni-rpr2*(dsin(tetmx)-dsint2)
        stor_j=rhditet_j*nj-rpr2*(dsint2 -dsin(tetmx-dtet))

        if(changeSg)then
        di = -di
        dsint2=dsin(tetmx-narc*di+ni*di)
        stor_i=rhditet_i*ni-rpr2*(dsint2 - dsin(tetmx-narc*di))
        stor_j=rhditet_j*nj-rpr2*(dsin(tetmx) - dsint2)
        end if

        stor_ijs=stor_i + stor_j
        stor=stor_ijs
ctt
          if(testpr)then
          write(*,*)'torus: arc_cross=',arc_cross
          write(*,*)'putsmall = ',putsm
          write(*,*)'dot_density=',dens
          write(*,'(a50,3i5)')
     &    'arc_points: on v0arc: narc,ni,nj=',
     &    narc,ni,nj
          write(*,*)'nn    x       y       z'
          do i=1,narc
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          write(*,*)'analitical area(1radrot):stor-i,j=',stor_i,stor_j
          end if

        goto 1001
	end if


          if(arc_cross)then
          costi = hij/rip
          costj = hij/rjp
          tetai = dacos(costi)
          tetaj = dacos(costj)

          rpsm = rpr + rad_sm
          costsm = hij/rpsm
          tetasm = dacos(costsm)
          sintsm = dsin(tetasm)

          do k=1,3
          psi_xyz(k) = 0.0d0
          psj_xyz(k) = 0.0d0
          end do

          psj_xyz(1) = rpsm*sintsm
          psi_xyz(1) = -psj_xyz(1)
ctt
          if(testpr)then
          write(*,*)'hij,rad_smooth probe=',hij,rad_sm
          write(*,'(a30,3f10.6)')
     &    'costi,costj,costsm=',costi,costj,costsm
          write(*,'(a45,2f8.4)')
     &   'arc_points:singular arc: psi,psj_xyz',psi_xyz(1),psj_xyz(1)
          end if

          pi05=pi*0.5d0
          dens2_sm = dens2*rad_sm

          ni=DINT((tetai-tetasm)*dens2_p+0.5d0)
          if(ni.eq.0.and.putsm.eq.1)ni=1

          ni_s=DINT((pi05-tetasm)*dens2_sm  + 0.5d0)
          if(ni_s.eq.0.and.putsm_smp.eq.1)ni_s=1

          nj=DINT((tetaj-tetasm)*dens2_p+0.5d0)
          if(nj.eq.0.and.putsm.eq.1)nj=1
          nj_s=ni_s

          di = (tetai-tetasm)/dfloat(ni)
          di2=0.5*di
          ditet_i=di2
          rhditet_i=rpr*hij*2.0d0*ditet_i
          sinditet_i=dsin(ditet_i)

          tet=tetai+di2
          if(ni.ge.1)then
          do i = 1, ni
          tet = tet - di
          v0arc(1,i)= -rpr*dsin(tet)
          v0arc(2,i)= hij - rpr*dcos(tet)
          v0arc(3,i)= 0.0d0
          arc_tet(i) = tet
          arc_type(i)=1
          end do
          end if

          di = (pi05-tetasm)/dfloat(ni_s)
          di2=0.5*di
          ditet_s=di2
          sinditet_s=2.0d0*dsin(ditet_s)*rad_sm*rad_sm
          tet=pi05-tetasm + di2
          if(ni_s.ge.1)then
          do i = 1+ni, ni+ni_s
          tet = tet - di
          v0arc(1,i) = rad_sm*dcos(tet) + psi_xyz(1)
          v0arc(2,i)= rad_sm*dsin(tet)
          v0arc(3,i)= 0.0d0
          arc_type(i)=3
          arc_tet(i) = tet
          end do
          end if

          ni = ni + ni_s

          if(testpr)then
          write(*,*)'arc_points:singular arc:v0arc-iat: nisp,ni_sm=',
     &    (ni-ni_s),ni_s
          write(*,*)'nn    x       y       z'
          do i=1,ni
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          end if

          di = (tetaj-tetasm)/dfloat(nj)
          di2=0.5*di
          ditet_j=di2
          rhditet_j=rpr*hij*2.0d0*ditet_j
          sinditet_j=dsin(ditet_j)
          tet=tetaj+di2

          if(nj.ge.1)then
          do i = ni+1, nj+ni
          tet = tet - di
          v0arc(1,i)= rpr*dsin(tet)
          v0arc(2,i)= hij - rpr*dcos(tet)
          v0arc(3,i)= 0.0d0
          arc_type(i)=2
          arc_tet(i) = tet
          end do
          end if
          di = (pi05-tetasm)/dfloat(ni_s)
          di2=0.5*di
          tet=pi05-tetasm + di2

          if(nj_s.ge.1)then
          do i = 1+ni+nj, ni+nj+nj_s
          tet = tet - di
          v0arc(1,i) = psj_xyz(1) - rad_sm*dcos(tet)
          v0arc(2,i) = rad_sm*dsin(tet)
          v0arc(3,i) = 0.0d0
          arc_type(i) = 4
          arc_tet(i) = tet
          end do
          end if

          nj = nj + ni_s

          if(testpr)then
          write(*,*)'arc_points:singular arc:v0arc-iat: nisp,ni_sm=',
     &    (ni-ni_s),ni_s
          write(*,*)'nn    x       y       z'
          do i=1+ni,nj+ni
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          end if

          narc = ni + nj

         stor_i = rprh*(tetai-tetasm)-rpr2*(dsin(tetai)-sintsm)

         smsph_i=rad_sm2*(1.0d0 - sintsm)


         stor_j= rprh*(tetaj-tetasm)-rpr2*(dsin(tetaj)-sintsm)

         stor_ijs = stor_i + 2.0d0*smsph_i + stor_j

ctt
        if(testpr)then
        write(*,'(a60,3f8.5)')
     &  'singul arc:analyt area(per 1rad):stor_i,smsph_i,stor_j',
     &  stor_i,smsph_i,stor_j
        end if

	end if

1001    continue

        ndotl = 0
        stor = 0.0d0


        do isf = 1,nijs_face
        nsf=2*isf - 1

        rot_ang = RotAng(nsf+1)-RotAng(nsf)

ctt
        if(testpr)then
        write(*,'(a45,2i5,2x,f8.6)')'arc_point:nijs_face,isf,rot_ang=',
     &  nijs_face,isf,rot_ang
        end if

        stor = stor + stor_ijs*rot_ang
cttw
        start_ang = RotAng(nsf) - RotAng(1)

        if(testpr)then
        write(*,'(a30,f8.4)')'arc_point: start_ang =',start_ang
        end if

        area_i=0.0d0
        area_j=0.0d0
        area_si=0.0d0
        area_sj=0.0d0


        do i = 1,narc

        dens2_r = v0arc(2,i)*dens2
        nr = DINT(rot_ang*dens2_r + 0.5d0)
        if(nr.eq.0.and.putsm.eq.1)nr=1

        dfi = rot_ang/dfloat(nr)

        if(dfi.ge.drot_max)then
        dfi=drot_max
        nr = DINT(rot_ang/dfi + 0.5d0)
        end if
ctt
        if(testpr)then
        write(*,*)'arc_point:i=',i,' numb_rot:nr=',nr
        end if
ctt

        if(nr.ge.1)then
         dfi = rot_ang/dfloat(nr)
         fi_rot = start_ang - dfi*0.5d0

         if(arc_type(i).eq.1)then

          d_arpi=dfi*(rhditet_i - rpr22*dcos(arc_tet(i))*sinditet_i)
         end if

         if(arc_type(i).eq.2)then

          d_arpj=dfi*(rhditet_j - rpr22*dcos(arc_tet(i))*sinditet_j)
         end if
         if(arc_type(i).eq.3.or.arc_type(i).eq.4)then
         d_arps=dfi*dsin(arc_tet(i))*sinditet_s
         end if

         do j= 1,nr
         fi_rot = fi_rot + dfi
         ndotl=ndotl + 1

         if(ndotl.gt.ndotmx)then
         write(*,*)'ERROR'
         stop
         end if

         dcost=dcos(fi_rot)
         dsint=dsin(fi_rot)

         dot_xyz(1,ndotl)=v0arc(1,i)
         dot_xyz(2,ndotl)=v0arc(2,i)*dcost
         dot_xyz(3,ndotl)=v0arc(2,i)*dsint

         if(arc_type(i).eq.1)then
         dot_vn(1,ndotl) = -v0arc(1,i)/rpr
         dot_vn(2,ndotl) = (hij - v0arc(2,i))*dcost/rpr
         dot_vn(3,ndotl) = (hij - v0arc(2,i))*dsint/rpr
         dot_area(ndotl) = d_arpi
         dot_type(ndotl) = sadl_ptype
         dot_atom(ndotl) = IATOM
         area_i=area_i + d_arpi
         end if

         if(arc_type(i).eq.2)then
         dot_vn(1,ndotl) = -v0arc(1,i)/rpr
         dot_vn(2,ndotl) = (hij - v0arc(2,i))*dcost/rpr
         dot_vn(3,ndotl) = (hij - v0arc(2,i))*dsint/rpr
         dot_area(ndotl) = d_arpj
         dot_type(ndotl) = sadl_ptype
         dot_atom(ndotl) = JATOM
         area_j=area_j + d_arpj
         end if
         if(arc_type(i).eq.3)then
         area_si=area_si+ d_arps
         do k = 1,3
         dot_vn(k,ndotl) = (dot_xyz(k,ndotl)-psi_xyz(k))/rad_sm
         end do
         dot_area(ndotl) =  d_arps
         dot_type(ndotl) = sadl_smtype
         dot_atom(ndotl) = IATOM
         end if

         if(arc_type(i).eq.4)then
         area_sj=area_sj+ d_arps
         do k = 1,3
         dot_vn(k,ndotl) = (dot_xyz(k,ndotl)-psj_xyz(k))/rad_sm
         end do

         dot_area(ndotl) =  d_arps
         dot_type(ndotl) = sadl_smtype
         dot_atom(ndotl) = JATOM
         end if

         if(testpr)then
	 write(*,'(a5,i3,a3,i5,a5,i3,a5,f6.3,1x,a4,3f6.3,a5,3f6.3)')
     &   'rotj:',j,'id=',ndotl,'type:',dot_type(ndotl),
     &   'area:',dot_area(ndotl),
     &   ' xyz',(dot_xyz(k,ndotl),k=1,3),
     &   'nvec=', (dot_vn(k,ndotl),k=1,3)
         end if

         end do
         end if
         end do

ctt
         if(testpr)then
         write(*,*)'arc_point: dot in loc xyz: saddle isf=',isf
         write(*,*)'compare areas(anal vs num) arci,arcj,arcsm'
         write(*,'(a18,2f12.8)')'stor_i,area_i',stor_i*rot_ang,area_i
         write(*,'(a18,2f12.8)')'stor_j,area_j',stor_j*rot_ang,area_j
         write(*,'(a35,3f12.8)')'smsphi*rot_ang,area_si,srea_sj',
     &   smsph_i*rot_ang,area_si ,area_sj
         end if

         end do

         ndot = ndotl

	return
	end

c
      	subroutine torus_rot_angl02(atm_mx,nbr_mx,nconc_mx,
     &                  iat,jat,UIJ,HIJ,nvertex,atvert_conc,cvertex,
     &                  cprobe_ijk,conc_gl_ijk,
     &                  icfIJ,concfS_ij,RotAng,SignRotS,AIJarcS,GarcS)

       implicit none

       integer atm_mx,nbr_mx,nconc_mx
       integer iat,jat
       integer nvertex(atm_mx)
       integer atvert_conc(nbr_mx,atm_mx)
       real*8 UIJ(3),HIJ
       real*8 cvertex(3,nbr_mx,atm_mx)
       real*8 cprobe_ijk(3,nconc_mx)
       integer conc_gl_ijk(6,nconc_mx)

       integer nbr_mxx
       parameter(nbr_mxx=180)
       integer iconcf,icf
       integer concf_ij(nbr_mxx)
       real*8 AIJarc(3,nbr_mxx)
       real*8  Garc(3,3,nbr_mxx)
       integer SignRot(nbr_mxx)
       integer order(nbr_mxx)

       integer i,j,iperm
       integer k,katt,kpos,ki,kj
       integer iavert,javert,kavert
       integer IATOM,JATOM,KATOM
       integer OPT_printsd,iiprint0,iiprint1,iiprint2

       real*8 VPI(3),VPJ(3),VPK(3),TVP(3),PIJK(3)
       real*8 UQ1(3),UT1(3),vvv(3)

       real*8 SIGN,aa1,aa2
       real*8 PI
       real*8 toler

       logical exist_i,exist_j,exist_k

c
	 real*8 PI2,null,onehalf,one,two,three,four,big,small
	 common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

       integer icfIJ
       real*8 AIJarcS(3,nbr_mx)
       real*8 GarcS(3,3,nbr_mx)
       real*8 RotAng(nbr_mx)
       integer SignRotS(nbr_mx)
       integer concfS_ij(nbr_mx)

       logical OPT_printwarn
c
       real*8 ANORM
       real*8 DOT

	PI=PI2
	toler=1.0d-6
	iiprint0=0
	iiprint1=1
	iiprint2=2

c	OPT_printsd=1
        OPT_printsd=0

	OPT_printwarn = .true.


c initialize
       IATOM=iat
       JATOM=jat
       icf=0
       do i=1,nbr_mxx
       concf_ij(i)=0
       order(i)=0
       end do


       if(OPT_printsd.ge.iiprint1)then
       write(*,*)'*****************************************'
       write(*,*)'START SUB torus_angl_rot'
       write(*,*)'nvertex(IATOM)=',nvertex(IATOM),
     &    ' nvertex(JATOM)=',nvertex(JATOM)
       end if
c
       if(nvertex(IATOM).ge.1.and.nvertex(JATOM).ge.1)then
       do i=1,nvertex(IATOM)

       iconcf=atvert_conc(i,IATOM)

       do j=1,nvertex(JATOM)
       if(iconcf.eq.atvert_conc(j,JATOM))then

       icf=icf+1

c control
       if(icf.gt.nbr_mx.or.icf.gt.nbr_mxx)then
       write(*,*)'Sub torus_rot_angl:ERROR: nbr_mx,nbr_mxx SMALL'
       stop
       end if

       concf_ij(icf)=iconcf
       iavert=i
       javert=j


	       if(OPT_printsd.ge.iiprint1)then
	       write(*,*)
	       write(*,*)'   icf =  ', icf
	       write(*,*)'Vertex Iat,Jat=',iavert,javert
               end if

        exist_i=.false.
        exist_j=.false.
        exist_k=.false.

	do k=1,3
	if(conc_gl_ijk(k,iconcf).eq.IATOM)then
	ki=k
	exist_i=.true.
	end if
	end do

	do k=1,3
	if(conc_gl_ijk(k,iconcf).eq.JATOM.and.k.ne.ki)then
	kj=k
	exist_j=.true.
	end if
	end do

	do k=1,3
	if(k.ne.ki.and.k.ne.kj.and.exist_i.and.exist_j)then
	exist_k=.true.
	kpos=k
	KATOM=conc_gl_ijk(k,iconcf)
	end if
	end do

c check
	if(.not.exist_i.and.exist_j.and.exist_k)then
	write(*,*)'atom K is not found'
	stop
	end if


	kavert=conc_gl_ijk(3+kpos,iconcf)

        do k=1,3
	VPI(k)=cvertex(k,iavert,IATOM)
	VPJ(k)=cvertex(k,javert,JATOM)
	VPK(k)=cvertex(k,kavert,KATOM)
	PIJK(k)=cprobe_ijk(k,iconcf)
        UQ1(k)= -VPI(k)-VPJ(k)
	TVP(k)= UQ1(k)- VPK(k)
        end do

c
	      if(OPT_printsd.ge.iiprint1)then
	write(*,*)'ARC  from iconcf=', iconcf
	write(*,*)'IATOM,JATOM,KATOM=',IATOM,JATOM,KATOM
	write(*,*)'iavert,javert,kavert=',iavert,javert,kavert
	write(*,*)'VERTEX VPI=',VPI
	write(*,*)'VERTEX VPJ=',VPJ
	write(*,*)'VERTEX VPK=',VPK
	write(*,*)'PROBE PIJK=',PIJK
              end if
c
	       aa1=DOT(UIJ,UQ1)
	       do k=1,3
	       UQ1(k)=UQ1(k)-aa1*UIJ(k)
	       end do

	       aa1= ANORM(UQ1)
	       do k=1,3
	       UQ1(k)=UQ1(k)/aa1
	       AIJarc(k,icf)=HIJ*UQ1(k)
	       end do

	       call CROSS(UIJ,UQ1,UT1)

              do K = 1,3
                 Garc(K,1,icf) = UIJ(K)
                 Garc(K,2,icf) = UQ1(K)
                 Garc(K,3,icf) = UT1(K)
              end do
c
	      if(OPT_printsd.ge.iiprint1)then
	write(*,*)'FRAME : '
	write(*,*)'UIJ=',UIJ
	write(*,*)'UQ1=',VPJ
	write(*,*)'UT1=',UT1
	      end if

	      SIGN = DOT(TVP,UT1)

	      if(SIGN.lt.null)then
	      SignRot(icf)=-1
	      else
	      SignRot(icf)=1
	      end if

	      if(OPT_printsd.ge.iiprint1)then
	      write(*,*)'SIGN of ARC ROTATION=',SignRot(icf)
	      write(*,*)'TVP=',TVP
	      write(*,*)'UT1_icf=',UT1
	      write(*,*)'UT1_1  =',(Garc(K,3,1),K=1,3)
	      end if


	      RotAng(1)=null
              if(SignRot(1).eq.-1)RotAng(1) = two*PI

	      if(icf.gt.1)then

	      aa1=DOT(Garc(1,3,1),UT1)

	      if(aa1.gt.one)then
	      if(OPT_printwarn)then
	      write(*,*)'WARNING in torus_ang_rot:aa1>one',aa1
	      end if
	      aa1 = one
	      end if
	      if(aa1.lt.-one)then
	      if(OPT_printwarn)then
	      write(*,*)'WARNING in torus_ang_rot:aa1<-one',aa1
	      end if
	      aa1 = -one
	      end if

	      RotAng(icf)=dacos(aa1)

              call CROSS(Garc(1,2,1),UQ1,vvv)

	      aa2=DOT(UIJ,vvv)


ctt      aa2=aa2*SignRot(1)

	      if(aa2.lt. 0.0d0)
     &	      RotAng(icf)=two*PI-RotAng(icf)

               if(OPT_printsd.ge.iiprint1)then
	       write(*,*)'icf=', icf
	       write(*,*)'COS(rotangl)=',aa1,' RotAng=',RotAng(icf)
	       end if

	      end if

	      end if
	      end do
	      end do
	      end if

	      icfIJ=icf

	      if(icfIJ.ne.2*(icfIJ/2))then
	      write(*,*)'SUbr torus_rot_ang: ERROR BAD icfIJ=',icfIJ
	      stop
	      end if

	     do i=1,icfIJ
	     order(i)=i
	     end do

            if(OPT_printsd.ge.iiprint1)then
            write(*,*)'Initial RotAng '
	    do i=1,icfIJ
	    write(*,*)i,RotAng(i),order(i)
	    end do
            end if

5000        iperm=0

	     do i=2,icfIJ
	     if(RotAng(i).lt.RotAng(i-1))then
	     iperm=iperm+1
	     aa1=RotAng(i-1)
	     RotAng(i-1)=RotAng(i)
	     RotAng(i) = aa1
	     j=order(i-1)
	     order(i-1)=order(i)
	     order(i) = j
	     end if
	     end do

	     if(iperm.ge.1)goto 5000

	    do i=1,icfIJ
	    j=order(i)
            SignRotS(i)=SignRot(j)
              do K = 1,3
	         AIJarcS(K,i) = AIJarc(K,j)
                 GarcS(K,1,i) = Garc(K,1,j)
                 GarcS(K,2,i) = Garc(K,2,j)
                 GarcS(K,3,i) = Garc(K,3,j)
              end do
	      concfS_ij(i)=concf_ij(j)
	      end do
ctt
           if(OPT_printsd.ge.iiprint1)then
            write(*,*)'SORTED:i, RotAng(i),SignRot(i),order(i) '
	    do i=1,icfIJ
	    write(*,*)i,RotAng(i),SignRotS(i),order(i)
	    end do
	    end if

        if(OPT_printsd.ge.iiprint1)then
       write(*,*)'END SUB torus_angl_rot'
       write(*,*)'*****************************************'
       end if

            return
	    end

c
	  subroutine cxdot_refcut09(GG,dotch,dotchr,dotrot,dotstat,
     &	  dotCO,UIJ,Ds,RI,dotaw,nchmx)


           implicit none

	   integer nchmx
	   real*8 GG(3,3)
	   real*8 dotch(3,nchmx)
	   real*8 dotchr(3,nchmx)
	   logical dotstat(nchmx)
	   logical dotrot(nchmx)
	   logical lastcut
	   real*8 dotCO(3),UIJ(3)
	   real*8 X(3),Y(3)
	   real*8 dotchg(3)
	   real*8 Ds,dotaw,RI

	 real*8 PI2,null,onehalf,one,two,three,four,big,small
	 common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small
         real*8 nulls
         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

	   integer i,j,k,k1,k2,nch
	   integer nnch
	   real*8 dnch,aa1,aa2
	   real*8 shif1,shif2,shif3
           integer print, iprint0,iprint1,iprintm
	   logical filter1,special(5)
	   integer dofilter
	   logical filter2,filter3,filter4

	   print =  -1
	   iprint0 = 0
	   iprint1 = 1
	   iprintm = -1
	   nnch =nchmx
	   filter1=.true.
	   filter2 = .false.
	   filter3 = .false.
	   filter4 = .false.

	   dofilter = -1
	   dofilter = 0
	   lastcut = .true.

	   do k=1,5
	   special(k)=.false.
	   end do
c
	  if(print.ge.iprint0)then
	  write(*,*)'********************************************'
	  write(*,*)'in cxdot_refcut:'
	  end if
	  if(print.ge.iprint0)then
	  write(*,*)'dotCO=',dotCO
	  write(*,*)'UIJ =',UIJ
	  write(*,*)'Ds, RI ',Ds, RI, dotaw
          end if

c check his position

	   dotaw = null
	   do k=1,3
	   dotCO(k)= dotchr(k,1)
	   end do
	   do k=1,3
	   X(k)=null
	   end do

	   nch = 0

	   do i=1,nnch

	   if(.not.dotrot(i))then
	   do k=1,3
	   dotchr(k,i)=null
	   end do
           do k2=1,3
	   do k1=1,3
	   dotchr(k1,i)=dotchr(k1,i) + GG(k1,k2)*dotch(k2,i)
	   end do
	   end do

	   dotrot(i)=.true.
	   end if

	   if(dotstat(i))then
	   aa1 = null
	   do k = 1,3
	   aa1=aa1 + dotchr(k,i)*UIJ(k)
	   end do

	  if(print.ge.iprint1)then
          write(*,*)'INITIAL LocXYZ child i =',i
	  write(*,*)'dotchg=',(dotch(k,i),k=1,3)
	  write(*,*)'ROTATED child i =',i
	  write(*,*)'dotchg=',(dotchr(k,i),k=1,3)
	  end if

	   if(aa1.gt.Ds )then
	  if(print.ge.iprint1)then
	  write(*,*)'delete child i=',i, ' dotDs=',aa1
          end if

	   dotstat(i) = .false.

	   end if
	   end if

	   end do

c
c
	  if(lastcut)then


1005      continue

	  if(filter1)then
	  if(.not.dotstat(1) .and.(dotstat(2).and.dotstat(3)
     &       .and.dotstat(4).and.dotstat(5).and.dotstat(6)
     &       .and.dotstat(7).and.dotstat(8).and.dotstat(9)))then
	    do i = 2,9
	    dotstat(i)=.false.
	    end do
	    special(1)=.true.
	    end if

          if((.not.dotstat(1).and. .not.dotstat(9)).and.
     &       (dotstat(2).and.dotstat(3).and.dotstat(4)
     &       .and.dotstat(5).and.dotstat(6)
     &       .and.dotstat(7).and.dotstat(8)))then
	  dotstat(2)=.false.
	  dotstat(8)=.false.
	  special(2) = .true.
	  end if

	  if((.not.dotstat(1).and. .not.dotstat(3)).and.
     &       (dotstat(9).and.dotstat(2).and.dotstat(4)
     &       .and.dotstat(5).and.dotstat(6)
     &       .and.dotstat(7).and.dotstat(8)))then
	  dotstat(2)=.false.
	  dotstat(4)=.false.
	  special(3)=.true.
	  end if

	  if((.not.dotstat(1).and. .not.dotstat(7)).and.
     &       (dotstat(9).and.dotstat(2).and.dotstat(3)
     &       .and.dotstat(5).and.dotstat(6)
     &       .and.dotstat(4).and.dotstat(8)))then
	  dotstat(6)=.false.
	  dotstat(8)=.false.
	  special(4)=.true.
	  end if

	  if((.not.dotstat(1).and. .not.dotstat(5)).and.
     &       (dotstat(9).and.dotstat(2).and.dotstat(4)
     &       .and.dotstat(3).and.dotstat(6)
     &       .and.dotstat(7).and.dotstat(8)))then
	  dotstat(6)=.false.
	  dotstat(4)=.false.
	  special(5)=.true.
	  end if
	  end if

	  if(filter2)then
c b1:
	  if((.not.dotstat(1) .and. .not.dotstat(2) .and.
     &        .not.dotstat(6)) .and.
     &       (dotstat(3).or.dotstat(4).or.dotstat(5)).and.
     &       (dotstat(7).or.dotstat(8).or.dotstat(9)))then
	  write(*,*)' b1: 1 2 6 - F'
	  do k=1,9
	  dotstat(k)=.false.
	  end do
	  end if
c b2:
	  if((.not.dotstat(1) .and. .not.dotstat(4) .and.
     &        .not.dotstat(8)) .and.
     &       (dotstat(2).or.dotstat(3).or.dotstat(9)).and.
     &       (dotstat(5).or.dotstat(6).or.dotstat(7)))then
	  write(*,*)' b2: 1 4 8 - F'
	  do k=1,9
	  dotstat(k)=.false.
	  end do
	  end if
c b3:
	  if((.not.dotstat(1) .and. .not.dotstat(3) .and.
     &        .not.dotstat(7)) .and.
     &       (dotstat(4).or.dotstat(5).or.dotstat(6)).and.
     &       (dotstat(2).or.dotstat(8).or.dotstat(9)))then
	  write(*,*)' b1: 1 3 7 - F'
	  do k=1,9
	  dotstat(k)=.false.
	  end do
	  end if
c b4:
	  if((.not.dotstat(1) .and. .not.dotstat(5) .and.
     &        .not.dotstat(9)) .and.
     &       (dotstat(3).or.dotstat(4).or.dotstat(2)).and.
     &       (dotstat(7).or.dotstat(8).or.dotstat(6)))then
	  write(*,*)' b1: 1 5 9 - F'
	  do k=1,9
	  dotstat(k)=.false.
	  end do
	  end if

	  end if

	  if(filter3)then
c b5:
	  if((.not.dotstat(3) .and. .not.dotstat(5)).and.
     &       dotstat(4).and.
     &       (dotstat(1).or.dotstat(2).or.dotstat(6).or.
     &        dotstat(7).or.dotstat(8).or.dotstat(9)))then
	     dotstat(4)=.false.
	  end if
c b6:
	  if((.not.dotstat(7) .and. .not.dotstat(5)).and.
     &       dotstat(6).and.
     &       (dotstat(1).or.dotstat(4).or.dotstat(8).or.
     &        dotstat(2).or.dotstat(3).or.dotstat(9)))then
	     dotstat(6)=.false.
	  end if
c b7:
	  if((.not.dotstat(7) .and. .not.dotstat(9)).and.
     &       dotstat(8).and.
     &       (dotstat(1).or.dotstat(2).or.dotstat(6).or.
     &        dotstat(3).or.dotstat(4).or.dotstat(5)))then
	     dotstat(8)=.false.
	  end if
c b8:
	  if((.not.dotstat(3) .and. .not.dotstat(9)).and.
     &       dotstat(2).and.
     &       (dotstat(1).or.dotstat(4).or.dotstat(8).or.
     &        dotstat(7).or.dotstat(5).or.dotstat(6)))then
	     dotstat(2)=.false.
	  end if

	  end if

	  if(filter4)then
c b9:
	  if((.not.dotstat(1) .and. .not.dotstat(2)))then
	  dotstat(3)=.false.
	  dotstat(9)=.false.
	  end if
c b10:
	  if((.not.dotstat(1) .and. .not.dotstat(4)))then
	  dotstat(3)=.false.
	  dotstat(5)=.false.
	  end if
c b11:
	  if((.not.dotstat(1) .and. .not.dotstat(6)))then
	  dotstat(7)=.false.
	  dotstat(5)=.false.
	  end if
c b12:
	  if((.not.dotstat(1) .and. .not.dotstat(8)))then
	  dotstat(7)=.false.
	  dotstat(9)=.false.
	  end if

	  end if

	  dofilter = dofilter + 1
c	  if(dofilter.le.0)goto 1005


             nch = 0
	   do k=1,3
	   X(k)=null
	   end do

	   do i = 1,nnch
	   if(dotstat(i))then
	   nch=nch + 1
           do k=1,3
	   X(k)= X(k) + dotchr(k,i)
	   end do
	   end if
	   end do

c result
	   if(nch.gt.0)then
	   dotaw = dfloat(nch)/dfloat(nnch)

	   call VNORM(X,Y)

	   do k=1,3
	   dotCO(k)=RI*Y(k)
	   end do
	   end if

	   if(nch.eq.0)then
	   dotaw = null
	   do k=1,3
	   dotCO(k)= dotchr(k,1)
	   end do

	   end if
c
	  if(print.ge.iprintm)then
	  do k=1,5
	  if(special(k))then
	  write(*,*)'res cxdot_refcut:special',special,'  nch ',nch
	  end if
	  end do
	  end if

	  end if

	  if(print.ge.iprint0)then
	  write(*,*)'nch =',nch
	  write(*,*)'dotCO=',dotCO
	  write(*,*)'dotstat=',(dotstat(k),k=1,nnch),'  dotaw =',dotaw
          end if

          return
	  end
c
	  subroutine cxdot_refcut05(GG,dotch,dotchr,dotrot,dotstat,
     &	  dotCO,UIJ,Ds,RI,dotaw,nchmx,lastcut)

           implicit none

	   integer nchmx
	   real*8 GG(3,3)
	   real*8 dotch(3,nchmx)
	   real*8 dotchr(3,nchmx)
	   logical dotstat(nchmx)
	   logical dotrot(nchmx)
	   logical lastcut
	   real*8 dotCO(3),UIJ(3)
	   real*8 X(3),Y(3)
	   real*8 dotchg(3)
	   real*8 Ds,dotaw,RI

	 real*8 PI2,null,onehalf,one,two,three,four,big,small
	 common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small
         real*8 nulls

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

	   integer i,j,k,k1,k2,nch
	   integer nnch
	   real*8 dnch,aa1,aa2
	   real*8 shif1,shif2,shif3
           integer print, iprint0,iprint1,iprintm
	   logical filter1

	   print = -1
	   iprint0 = 0
	   iprint1 = 1
	   iprintm = -1
	   nnch =nchmx

	   filter1=.true.

c	   lastcut = .true.

	  if(print.ge.iprint0)then
	  write(*,*)'********************************************'
	  write(*,*)'in cxdot_refcut:'
	  end if
	  if(print.ge.iprint0)then
	  write(*,'(a10,3f10.6)')'dotCO=',dotCO
	  write(*,'(a10,3f10.6)')'UIJ =',UIJ
	  write(*,'(a12,3f10.6)')'Ds, RI, dotaw: ',Ds, RI, dotaw
          end if


	   do i=1,nnch

	   if(.not.dotrot(i))then
	   do k=1,3
	   dotchr(k,i)=null
	   end do
           do k2=1,3
	   do k1=1,3
	   dotchr(k1,i)=dotchr(k1,i) + GG(k1,k2)*dotch(k2,i)
	   end do
	   end do

	   dotrot(i)=.true.
	   end if

	   if(dotstat(i))then
	   aa1 = null
	   do k = 1,3
	   aa1=aa1 + dotchr(k,i)*UIJ(k)
	   end do

	  if(print.ge.iprint1)then
          write(*,'(a20,i5)')'INITIAL LocXYZ child i =',i
	  write(*,'(a12,3f10.6)')'dotch-loc:',(dotch(k,i),k=1,3)
	  write(*,'(a12,i5)')'ROTATED child i:',i
	  write(*,'(a12,3f10.6)')'dotch-rot:',(dotchr(k,i),k=1,3)
	  end if

	   if(aa1.gt.Ds)then
	  if(print.ge.iprint1)then
	  write(*,'(a12,i5,a14,2f10.6)')
     &	      'delete child i=',i,' dotDs, aa1:',Ds,aa1
          end if

	   dotstat(i) = .false.

	   end if
	   end if

	   end do

c1005      continue

	  if(lastcut)then

	  if(filter1)then

	  if(.not.dotstat(1).and.(dotstat(2).and.dotstat(3)
     &       .and.dotstat(4).and.dotstat(5)))then
	  do k = 2,5
	  dotstat(k)=.false.
	  end do
	  end if

	  if((.not.dotstat(2) .and. .not.dotstat(4)).and.
     &    (dotstat(3).and.dotstat(5)))then
	    dotstat(3)=.false.
	    dotstat(5)=.false.
	    end if

	  if((.not.dotstat(3) .and. .not.dotstat(5)).and.
     &    (dotstat(2).and.dotstat(4)))then
	    dotstat(2)=.false.
	    dotstat(4)=.false.
	    end if
	  end if

           nch = 0
	   do k=1,3
	   X(k)=null
	   end do

	   do i = 2,nnch

	   if(dotstat(i))then
	   nch=nch + 1
           do k=1,3
	   X(k)= X(k) + dotchr(k,i)
	   end do
	   end if
	   end do

	   if(nch.gt.0)then
	   dotaw = dfloat(nch)/dfloat(nnch - 1)
	   call VNORM(X,Y)
	   do k=1,3
	   dotCO(k)=RI*Y(k)
	   end do

	   else
	   dotaw = null
	   end if
c
	  if(print.ge.iprint0)then
	  write(*,*)'lastcut:Totalidotch: nch =',nch
	  write(*,'(a10,3f10.6)')'dotCO=',dotCO
	  write(*,'(a10,5L1,a10,f10.6)')
     &	      'dotstat=',(dotstat(k),k=1,nnch),'  dotaw =',dotaw
          end if
	  end if

          return
	  end
c
       subroutine inprizm(v1,v2,v3,u,inp)

       implicit none
       real*8 v1(3),v2(3),v3(3),u(3)
       logical*1 inp
       real*8 DOT
       real*8 s,v23(3),v31(3),v12(3)
	   real*8 toler
       logical*1 testpr

	   testpr=.false.

       if(testpr)then
	   write(*,'(a13,3(a4,3f7.3))')
     &	   'inprizM:',' v1:',v1,' v2:',v2,' v3:',v3
	   write(*,'(a6,3f7.3)')'u: ',u
	   end if

       inp = .false.

	   toler = -0.1d-06

       call CROSS(v2,v3,v23)
       s = DOT(v1,v23)

	   if(testpr)then
	   write(*,'(a6,3f7.3)')'v23: ',v23
	   write(*,*)'inprizM:s,s*DOT(u,v23):',s,s*DOT(u,v23)
	   end if

       if(s*DOT(u,v23).LE.toler)goto 100

       call CROSS(v3,v1,v31)
	   if(testpr)then
	   write(*,'(a6,3f7.3)')'v31: ',v31
	   write(*,*)'inprizM:s*DOT(u,v31):',s*DOT(u,v31)
	   end if
       if(s*DOT(u,v31).LE.toler)goto 100

       call CROSS(v1,v2,v12)
	   if(testpr)then
	   write(*,'(a6,3f7.3)')'v12: ',v12
	   write(*,*)'inprizM:s*DOT(u,v12):',s*DOT(u,v12)
	   end if
       if(s*DOT(u,v12).LE.toler)goto 100

       inp = .true.

	   if(testpr)then
	   write(*,*)'inprzm:',inp
	   end if

100    return
       end

         subroutine order_dots(natoms,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,
     &                 dot_num_atom,dot_startn_atom,ndots)

        include'surf-sims.h'
        include'input_sims.h'

        real*8 dotarea(maxdot),dotcrd(3,maxdot),dotnrm(3,maxdot)
        integer natoms,ndots
        integer*4 dot_IATNUM(maxdot),dot_ISHAPE(maxdot)
        integer*4 dot_num_atom(maxatm),dot_startn_atom(maxatm)

        integer*4 locint(maxdot),dotpos(maxatm)
        real*8 loc1real(maxdot),loc3real(3,maxdot)
        real atomarea(maxatm),atomarea2(maxatm),area_tot
        integer*4 i,j,k,ndin
        integer iprev,nsurfatom
        integer kanalpl
        logical start(maxatm)
        logical CONTROL

        kanalpl=kanalp
        CONTROL=.false.

        do i=1,maxatm
        dot_num_atom(i)=0
        dot_startn_atom(i)=0
        dotpos(i)=0
        atomarea(i)=0.0
        atomarea2(i)=0.0
        start(i)=.false.
        end do

        do j=1,ndots
        dot_num_atom(dot_IATNUM(j))=dot_num_atom(dot_IATNUM(j))+1
        atomarea(dot_IATNUM(j))=atomarea(dot_IATNUM(j))+dotarea(j)
        end do
        ndin=0
        do i=1,natoms
        if(dot_num_atom(i).ne.0)then
        dot_startn_atom(i)=ndin + 1
        ndin = ndin + dot_num_atom(i)
        end if
        end do

        do i=1,natoms
        dotpos(i)=dot_startn_atom(i)
        start(i)=.false.
        end do

        do i=1,ndots
        if(start(dot_IATNUM(i)))then
        dotpos(dot_IATNUM(i))=dotpos(dot_IATNUM(i))+1
        else
        start(dot_IATNUM(i))=.true.
        end if
        locint(dotpos(dot_IATNUM(i)))=dot_ISHAPE(i)
        loc1real(dotpos(dot_IATNUM(i)))=dotarea(i)
        do k=1,3
        loc3real(k,dotpos(dot_IATNUM(i)))=dotcrd(k,i)
        end do
        end do

        do i=1,ndots
        dot_ISHAPE(i)=locint(i)
        dotarea(i)=loc1real(i)
        do k=1,3
        dotcrd(k,i)=loc3real(k,i)
        end do
        end do

        do i=1,natoms
        dotpos(i)=dot_startn_atom(i)
        start(i)=.false.
        end do

        do i=1,ndots
        if(start(dot_IATNUM(i)))then
        dotpos(dot_IATNUM(i))=dotpos(dot_IATNUM(i))+1
        else
        start(dot_IATNUM(i))=.true.
        end if
        locint(dotpos(dot_IATNUM(i)))=dot_IATNUM(i)
        do k=1,3
        loc3real(k,dotpos(dot_IATNUM(i)))=dotnrm(k,i)
        end do
        end do
        do i=1,ndots
        dot_IATNUM(i)=locint(i)
        do k=1,3
        dotnrm(k,i)=loc3real(k,i)
        end do
        end do

        if(CONTROL)then
        write(kanalpl,*)'CONTROLtest: ndot_total =' , ndin
        write(kanalpl,*)'test ia  area1    area2   ndot   nstart'
        area_tot=0.0
        do i=1,natoms
        if(dot_startn_atom(i).ne.0)then
        do j=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1
        atomarea2(i)=atomarea2(i)+dotarea(j)
        end do
        area_tot=area_tot+atomarea2(i)
        end if
        write(kanalpl,'(i5,2x,f8.3,f8.3,i5,2x,i5)')
     &  i,atomarea(i),atomarea2(i),dot_num_atom(i),dot_startn_atom(i)
        end do
        write(kanalpl,'(a30,f10.4)')'CONTROLtest: area_tot=',area_tot
        end if
c       write(kanalpl,'(a31)')'SIMS: dot ordering is complete'
        return
        end
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
c.  Double precission input/output  arguments

c.creates *.kin file for *.pdb and surface dots with normal vectors
       subroutine surf_kin(coords,
     &                 atnamel,resnumbl,resnamel,
     &                 natoms,atom_rad,dotden,
     &                 dotcrd,dotnrm,dotarea,
     &                 dot_IATNUM,dot_ISHAPE,
     &         dot_num_atom,dot_startn_atom,ndots)

c INPUT:
c       coords(3,*) coord of atoms
c       atnamel(*)  atom name
c       resnumbl(*) residue number of atom(*)
c       resnamel(*) residue name of atom(*)
c
c RESULT:
c       SIMS-pdb-dot.kin  file
c
        include'surf-sims.h'
        include'input_sims.h'

        real*8 atom_rad(*),coords(3,*),dotarea(*),
     &       dotcrd(3,*),dotnrm(3,*),dotden

        character*3 atnamel(*),resnamel(*)
	integer resnumbl(*)
	integer natoms,ndots
	integer*4 dot_IATNUM(maxdot),dot_ISHAPE(maxdot)
        integer*4 dot_num_atom(maxatm),dot_startn_atom(maxatm)
c
c local VARIABLES:
c
	integer dtypemx,nvertmx
	parameter (dtypemx=5)
	parameter (nvertmx=7)
	integer ndotimage(dtypemx)
	real*8 pdotimage(2,nvertmx,dtypemx)
        real*8 vert(3,nvertmx)
        integer atom_connect(4,maxatm)
        integer id,iatnn,nnnd,normbad
        integer iprint0,iprint1
        integer kanalpl
	integer i,j,k,l,nnn,dtype
	integer kanalzl
	real*8 normfff,dotxyz(3)
	real*8 x(3),y(3)
	real*8 ddotsh,a
	character*4 satom
	character*1 linksimb,leftb,rightb,point,vect,comma
	character*1 simbol,space
	character*8 day
	character*9 hour
c       logical OPT_pdbdotkin
c	logical OPT_dotnrmkin
	logical CONTROL

	data ndotimage /5,5,4,6,7/

caug9826: 1=CX(square); 2=SD(square); 3=CC(triangle);
caug9826: 4=smoothed(SDcone)(pentagon);
caug9826: 5=smoothed(CCedge)(sixagone)
	data pdotimage
     &	/1.0,0.0, 0.0,-1.0, -1.0,0.0, 0.0,1.0, 1.0,0.0,
     &   0.0,0.0, 0.0,0.0,
     &  1.0,0.0, 0.0,-1.0, -1.0,0.0, 0.0,1.0, 1.0,0.0, 0.0,0.0,
     &  0.0,0.0,
     &  1.0,0.0, -0.45,-0.81, -0.45,0.81, 1.0,0.0, 0.0,0.0,
     &  0.0,0.0, 0.0,0.0,
     &  0.0,1.0, 0.88,0.47, 0.59,-0.81, -0.59,-0.81, -0.88,0.47,
     &  0.0,1.0, 0.0,0.0,
     &  0.0,1.0, 0.87,0.5, 0.87,-0.5, 0.0,-1.0, -0.87,-0.5,
     &  -0.87,0.5, 0.0,1.0/

        CONTROL = .false.
c        OPT_pdbdotkin=.true.
c	OPT_dotnrmkin=.true.
c	OPT_dotnrmkin=.false.

	normfff=1.0d0
	kanalzl=kanalz
	ddotsh = 0.25d0
	if(CONTROL)then
	write(*,*)'In surf_kin:1:'
         end if
	 call connect_aamg(coords,
     &                 atnamel,resnumbl,resnamel,
     &                 natoms,
     &                 atom_connect)

	if(OPT_pdbdotkin)then
	open(unit=kanalzl,file='SIMS_pdbdot.kin',status='unknown')
        write(kanalzl,'(a5)')'@text'
	write(kanalzl,*)'SMOOTH INVARIANT MOLECULAR SURFACE '
	write(kanalzl,*)'AUTHOR       Yury Vorobjev, 1998   '
c
c      call date(day)
c      call time(hour)
c      write(kanalzl,'(a12,a9,a3,a8)')'SIMS: start:',day,' : ',hour

       write(kanalzl,'(a11)')'@kinemage 1'
       write(kanalzl,'(a8)')'@caption'
       write(kanalzl,*)'coordinates from file: molec.cor'
       write(kanalzl,'(a9)')'@onewidth'
       write(kanalzl,'(a17)')'@group {molecule}'
       write(kanalzl,'(a21)')'@subgroup {mainchain}'
       write(kanalzl,'(a29)')'@vectorlist {mc} color= white'
	leftb='{'
	rightb='}'
        point='P'
	vect ='L'
	comma=','
	space=' '

	do i=1,natoms

	do l=1,4

	j=atom_connect(l,i)
	if(j.gt.i)then
c Point for atom i:
	write(kanalzl,'(a1,1x,a3,1x,a3,1x,i4,1x,a1,
     &                             1x,a1,1x,3(f8.3,a1,1x),a1)')
     &        leftb,atnamel(i),resnamel(i),resnumbl(i),rightb,
     &        point,coords(1,i),comma,
     &              coords(2,i),comma,
     &              coords(3,i),space,space

c Vector i-->l:
	write(kanalzl,'(a1,1x,a3,1x,a3,1x,i4,1x,a1,
     &                             1x,a1,1x,3(f8.3,a1,1x),a1)')
     &        leftb,atnamel(j),resnamel(j),resnumbl(j),rightb,
     &        vect,coords(1,j),comma,
     &             coords(2,j),comma,
     &             coords(3,j),space,space
	end if

	end do
	end do
cend  all atom

c add dotnrm
	if(OPT_dotnrmkin)then
	 write(kanalzl,'(a33)')'@vectorlist {surfn} color= yellow'

	do i=1,natoms

	do id=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1

c dots for atom i:
	write(kanalzl,'(a1,1x,a3,1x,a3,1x,i4,1x,i5,1x,a1,
     &                             1x,a1,1x,3(f8.3,a1,1x),a1)')
     &        leftb,atnamel(i),resnamel(i),resnumbl(i),id,rightb,
     &        point,dotcrd(1,id),comma,
     &              dotcrd(2,id),comma,
     &              dotcrd(3,id),space,space
c Vector i-->l:

	do k=1,3
	dotxyz(k)=dotcrd(k,id) + normfff*dotarea(id)*dotnrm(k,id)
	end do

	write(kanalzl,'(a1,1x,a3,1x,a3,1x,i4,1x,i5,1x,a1,
     &                             1x,a1,1x,3(f8.3,a1,1x),a1)')
     &        leftb,atnamel(i),resnamel(i),resnumbl(i),id,rightb,
     &        vect,dotxyz(1),comma,
     &             dotxyz(2),comma,
     &             dotxyz(3),space,space

	end do
	end do
	end if

c add dotsquare
	if(OPT_dotnrmkin)then
	 write(kanalzl,'(a33)')'@vectorlist {surfsq} color= green'

	do i=1,natoms

	do id=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1

	call VPERP(dotnrm(1,id),x)
        call CROSS(dotnrm(1,id),x,y)
	a = ddotsh*dsqrt(dotarea(id))

c dots for atom i:
	dtype = dot_ISHAPE(id)
	do l=1,ndotimage(dtype)
	do k=1,3
	vert(k,l)=dotcrd(k,id) +
     &     	a*pdotimage(1,l,dtype)*x(k) +
     &          a*pdotimage(2,l,dtype)*y(k)
	end do

	simbol=vect
	if(l.eq.1)simbol=point

	write(kanalzl,'(a1,1x,a3,1x,a3,1x,i4,1x,i5,1x,a1,
     &                             1x,a1,1x,3(f8.3,a1,1x),a1)')
     &        leftb,atnamel(i),resnamel(i),resnumbl(i),id,rightb,
     &        simbol,vert(1,l),comma,
     &               vert(2,l),comma,
     &               vert(3,l),space,space

        end do
c
	end do
	end do
	end if

c dots
c	write(kanalzl,'(a29)')'@dotlist {surf} color= yellow'
	write(kanalzl,'(a28)')'@dotlist {surf} color= green'
	do i=1,natoms
	do id=dot_startn_atom(i),dot_startn_atom(i)+dot_num_atom(i)-1

	do k=1,3
	dotxyz(k)=dotcrd(k,id)
	end do

	write(kanalzl,'(3(f8.3,a1,1x))')
     &        dotxyz(1),comma,dotxyz(2),comma,dotxyz(3),space

	end do
	end do

	close(kanalzl)

	end if

	RETURN
        END

         subroutine connect_aamg(atmxyz,
     &                 atnamel,resnumbl,resnamel,
     &                 natoms,
     &                 atom_connect)

        include'surf-sims.h'
	include'input_sims.h'
        real*8 atmxyz(3,*)
	character*3 atnamel(*),resnamel(*)
	integer resnumbl(*)
	integer natoms
	integer atom_connect(4,*)


c local variables
	  integer atom_neighnumb(maxatm)
	  integer nat_inres(nres_max)
	  integer stop_res(nres_max),start_res(nres_max)
	  integer ia,k,ja,ir,irn,jr,iatom,jatom,natjr
	  integer i,nneigh, iiprint0,iiprint1
          real*4 radvi,radvj,radvC,radvH,d2
          character*1 atnamei,atnamej
          integer nvalneigh_max
	  integer nres,strtres
	  integer nnatmx
          logical CONTROL
	  integer outdetll
	  integer kanalpl

	kanalpl = kanalp
	 outdetll = 0
c        CONTROL = .true.
        CONTROL = .false.
        iiprint0=0
        iiprint1=1

	 nvalneigh_max=4
	 nnatmx=natoms
c         nres=0
         strtres=resnumbl(1)
         nres=strtres-1
	 do i = 1,natoms
	  if(resnumbl(i).gt.nres)then
	   nres=resnumbl(i)
	   start_res(nres)=i
           stop_res(nres-1)=i-1
          endif
	 end do
         stop_res(nres)=natoms
	 stop_res(nres)=natoms
	if(nres.gt.nres_max)then
	write(kanalpl,*)'PARAMETR:ERROR:connect_aamg:
     &  nres_max TOO small',nres
        stop
	end if
	 if(CONTROL)then
	 write(*,*)'In connect_aamg: nres:', nres
	 end if

	radvC=0.95
	radvH=0.60
c
	do ia=1,nnatmx
	atom_neighnumb(ia)=0
	do k=1,nvalneigh_max
	atom_connect(k,ia)=0
	end do
	end do

c
	 do ia=1,nres_max
	 nat_inres(ia)=0
	 end do
         do ir=strtres,nres
         nat_inres(ir)=stop_res(ir)-start_res(ir)+1
         end do
        do ir=strtres,nres
	if(nat_inres(ir).ge.1)then
	do iatom=start_res(ir),stop_res(ir)
        if(atom_neighnumb(iatom).eq.nvalneigh_max)then
        goto 101
        end if

	nneigh=0
	radvi=radvC
	atnamei=atnamel(iatom)(1:1)
        if(atnamei.eq.'H')radvi=radvH
c
        irn=ir+1
	if(irn.gt.nres)irn=nres

	do jr=ir,irn
	natjr=nat_inres(jr)
        if(natjr.ge.1)then
        do jatom=start_res(jr),stop_res(jr)
	radvj=radvC
	atnamej=atnamel(jatom)(1:1)
        if(atnamej.eq.'H')radvj=radvH

	d2=0.0
	do k=1,3
	d2=d2+(atmxyz(k,iatom)-atmxyz(k,jatom))**2
	end do

	if(iatom.lt.jatom.and.d2.lt.(radvi+radvj)**2)then

	if(outdetll.ge.iiprint1)then
	write(kanalpl,*)'1out: ir, jr:',ir,jr
        write(kanalpl,*)'iatom=',iatom,'jatom=',jatom,'d2=',d2
c        write(kanalpl,*)'radvi,radvj=',radvi,radvj
        end if
        atom_neighnumb(iatom)=atom_neighnumb(iatom)+1
        atom_neighnumb(jatom)=atom_neighnumb(jatom)+1
        if(atom_neighnumb(iatom).gt.nvalneigh_max.or.
     &  atom_neighnumb(jatom).gt.nvalneigh_max)then
         write(kanalpl,*)'PARAM nvalneigh_max is too low'
        end if
	atom_connect(atom_neighnumb(iatom),iatom)=jatom
	atom_connect(atom_neighnumb(jatom),jatom)=iatom
	end if

	end do
        end if
	end do

101     continue
	end do
        end if
	end do
c
	if(outdetll.ge.iiprint1)then
	write(kanalpl,*)'connect_aamg:atom-atom bond matrix'
        write(kanalpl,*)'nres=', nres
	do ia=1,nnatmx
	write(kanalpl,'(a4,i4,2x,a5,4i4)')'ia=',ia,
     &	'list=',(atom_connect(k,ia),k=1,4)
	end do
	end if
	if(CONTROL)then
	write(*,*)'connect_aamg:atom-atom bond matrix'
        write(*,*)'nres=', nres
	do ia=1,nnatmx
	write(*,'(a4,i4,2x,a5,4i5)')'ia=',ia,
     &	'list=',(atom_connect(k,ia),k=1,4)
	end do
	write(*,*)'END connect_aamg:'
	end if

        return
	end
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
       subroutine smooth_ccdot02c(nprob_mx,ndot_mx,
     &                  RP,rad_sm,
     &                  cprobe_ijk,
     &                  dotcrd,dotnvec,dotarea,dot_type,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl,WDYONRm)

        include 'surf-sims.h'

       integer nprob_mx,ndot_mx
       real*8 RP,rad_sm

       real*8 cprobe_ijk(3,nprob_mx)

       real*8 dotcrd(3,ndot_mx)
       real*8 dotnvec(3,ndot_mx)
       real*8 dotarea(ndot_mx)
       integer dot_type(ndot_mx)

        logical*1 WDYONRm(ndot_mx)
        integer WDYONSm(2,ndot_mx)

        integer ndot_smmx
        integer prpr1_smooth(2,ndot_smmx)
        integer dot_smooth_gl(ndot_smmx)
        integer ndot_smooth

        real*8 PI2,null,onehalf,one,two,three,four,big,small
        common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

        integer wdyonsm_n(maxdot)
    	integer indxdclose(ndot_smoothmx)
    	logical*1 dotCdist(ndot_smoothmx)
        real*8 dot_sm_loc(3,ndot_smoothmx)
        real*8 dot_nv_loc(3,ndot_smoothmx)
    	real*8 dot_are_loc(ndot_smoothmx)

    	integer id,k,idgl
    	integer npi,npj
    	integer jd,jdgl,cd,cdgl
        integer smooth_type
    	integer ndclose
        real*8 PI22,prs,rpdrs,rsdrp
        real*8 pxyz(3),p1xyz(3),p1mp(3)
        real*8 dot(3),dotp(3),p1pp(3),dotpsm(3)
        real*8 vvp(3),dot1(3),psm(3),dnvec(3)
        real*8 cost,sint,dv,ddv,ddd
        real*8 Hpp,hsm,kar,aret,are1t,are1,are
        real*8 dist,distav,dd
		real*8 costw,costw1,sinrs
    	real*8 kdotclose
        logical*1 testpr
        integer nsmoothl
    	logical OPT_clastering

caug98
        testpr = .false.
    	OPT_clastering = .true.

        smooth_type = 5
        nsmoothl = 0
        aret=0.0d0
        are1t=0.0d0

        if(rad_sm.le.0.0d0.or.ndot_smooth.le.0)goto 1000

        PI22 = PI2*0.5d0
        prs = RP + rad_sm
        rpdrs = RP/rad_sm
        rsdrp = rad_sm/RP

        do id = 1, ndot_smooth
        idgl=dot_smooth_gl(id)
        wdyonsm_n(idgl) = 0
        dot_are_loc(id) = 0.0d0
        do k =1,3
        dot_sm_loc(k,id) = 0.0d0
        dot_nv_loc(k,id) = 0.0d0
        end do
        end do

        if(testpr)then
		write(*,*)
		write(*,*)'in smooth_ccdot02c: start*****************'
        write(*,'(a58)')
     & 'dotListToBeSmoothed: id,idg, dot_type,WDYONSm,npi,npj,WDYONRm'
        if(ndot_smooth.le.0) goto 1000
        do id = 1, ndot_smooth
        idgl=dot_smooth_gl(id)
        npi=prpr1_smooth(1,id)
        npj=prpr1_smooth(2,id)
        write(*,'(20x,6i6,3x,L1)')
     &  id,idgl,dot_type(idgl),WDYONSm(1,idgl),
     &  npi,npj,WDYONRm(idgl)
        end do
        end if

        if(ndot_smooth.le.0) goto 1000
        do id = 1, ndot_smooth

        idgl=dot_smooth_gl(id)

        if(WDYONSm(1,idgl).gt.0)then
        nsmoothl = nsmoothl + 1
        npi=prpr1_smooth(1,id)
        npj=prpr1_smooth(2,id)

        do k=1,3
        pxyz(k)=cprobe_ijk(k,npi)
        p1xyz(k)=cprobe_ijk(k,npj)
        p1mp(k)= p1xyz(k)-pxyz(k)
        dot(k) = dotcrd(k,idgl)
        dotp(k)= dot(k) - pxyz(k)
        p1pp(k)= 0.5d0*(p1xyz(k)+pxyz(k))
        end do

        if(testpr)then
	    write(*,'(a16,3f8.3,a4,3f8.3)')
     &  'smooth_cc: p:',pxyz, ' p1:',p1xyz
        end if

        dist=dsqrt(p1mp(1)**2+p1mp(2)**2+p1mp(3)**2)
		if(dist.le.0.0d0)then
        write(*,'(a40,3i6,f8.4)')'smooth:id, npi,npj,dist=zero',
     &  id,npi,npj,dist
        goto 2000
        end if

        p1mp(1)=p1mp(1)/dist
        p1mp(2)=p1mp(2)/dist
        p1mp(3)=p1mp(3)/dist

        dd=dotp(1)*p1mp(1)+dotp(2)*p1mp(2)+dotp(3)*p1mp(3)
        vvp(1)=dotp(1)-p1mp(1)*dd
        vvp(2)=dotp(2)-p1mp(2)*dd
        vvp(3)=dotp(3)-p1mp(3)*dd

        dv=dsqrt(vvp(1)**2+vvp(2)**2+vvp(3)**2)

        if(dv.le.0.0d0)then
       write(*,'(a40,3i6,f8.4)')'smooth:id, npi,npj,dv=zero',
     &  id,npi,npj,dv
        goto 2000
        end if

        dist=0.5d0*dist
        cost = dist/prs
        sint = dsqrt(1.0d0 - cost**2)
        Hpp  = RP*sint
        hsm = rad_sm*sint
        ddv=prs*sint/dv
    	costw = dist/RP

		kar=0.0d0
		if(dist.le.RP+toler_cross)then
		if(dist.lt.RP)then
		costw = dacos(costw)
		else
		costw = 0.0d0
		end if

		costw1=dacos(cost)-costw
		sinrs = dasin(cost)
		kar = rad_sm*sinrs/(RP*costw1)
		else
		kar=0.0d0
		end if

        psm(1)=p1pp(1)+vvp(1)*ddv
        psm(2)=p1pp(2)+vvp(2)*ddv
        psm(3)=p1pp(3)+vvp(3)*ddv

        dotpsm(1)=dot(1) - psm(1)
        dotpsm(2)=dot(2) - psm(2)
        dotpsm(3)=dot(3) - psm(3)

        ddd=dsqrt(dotpsm(1)**2+dotpsm(2)**2+dotpsm(3)**2)
        if(ddd.le.0.0d0)then
       write(*,'(a40,3i6,f8.4)')'smooth:id, npi,npj,ddd=zero',
     &  id,npi,npj,dv
        goto 2000
        end if

        ddd = rad_sm/ddd
        dot1(1)=psm(1)+ddd*dotpsm(1)
        dot1(2)=psm(2)+ddd*dotpsm(2)
        dot1(3)=psm(3)+ddd*dotpsm(3)

        dnvec(1)=(dot1(1) - psm(1))/rad_sm
        dnvec(2)=(dot1(2) - psm(2))/rad_sm
        dnvec(3)=(dot1(3) - psm(3))/rad_sm

        are=dotarea(idgl)
        if(WDYONSm(1,idgl).gt.1)are=dotarea(idgl)/WDYONSm(1,idgl)
        aret = aret + are
        are1=are*kar
        are1t =  are1t + are1

        dot_type(idgl) = smooth_type
        dot_are_loc(id)=are1

        do k=1,3
        dot_sm_loc(k,id)=dot1(k)
        dot_nv_loc(k,id)=dnvec(k)
        end do

        if(testpr)then
         write(*,'(a30,4i6)')
     &  'Res:smooth_dot:id,idgl,npi,npj=',id,idgl,npi,npj
         write(*,'(a5,2f8.3,a6,3f8.3)')' d,H:',dist,Hpp,' psm:',psm
         write(*,'(a8,2f6.3,a16,3f7.3,1x,3f7.3)')
     &  'dv,ddv:',dv,ddv,' p1pp(3),vvp(3)',p1pp,vvp
         write(*,'(a70)')
     &  'Res: dotxyz      dotsmxyz      dotnrm     knormarea,area'
         write(*,'(3f7.3,1x,3f7.3,1x,3f6.3,1x,2f6.3)')
     &  dot,dot1,dnvec,kar,are1
        write(*,*)
        end if


            end if
2000        end do

        do id = 1, ndot_smooth
	    idgl=dot_smooth_gl(id)
        if(WDYONSm(1,idgl).gt.0)then
        wdyonsm_n(idgl) = 0
        dotarea(idgl) = 0.0d0
        do k =1,3
        dotcrd(k,idgl) = 0.0d0
        dotnvec(k,idgl) = 0.0d0
        end do
	    end if
        end do

	   do id = 1, ndot_smooth
	   idgl=dot_smooth_gl(id)
        if(WDYONSm(1,idgl).gt.0)then
        wdyonsm_n(idgl)=wdyonsm_n(idgl)+1

        dotarea(idgl)=dotarea(idgl) + dot_are_loc(id)

        do k=1,3
        dotcrd(k,idgl)=dotcrd(k,idgl)
     &      + (dot_sm_loc(k,id) - dotcrd(k,idgl))/wdyonsm_n(idgl)
        dotnvec(k,idgl)=dotnvec(k,idgl)
     &      + (dot_nv_loc(k,id) - dotnvec(k,idgl))/wdyonsm_n(idgl)
        end do

        if(WDYONSm(1,idgl).gt.1)then
        dd=dsqrt(dotnvec(1,idgl)**2+dotnvec(2,idgl)**2+
     &  dotnvec(3,idgl)**2)
        dotnvec(3,idgl)=dotnvec(3,idgl)/dd
        dotnvec(2,idgl)=dotnvec(2,idgl)/dd
        dotnvec(1,idgl)=dotnvec(1,idgl)/dd
        end if

        if(testpr)then
        write(*,'(a24,2i6)')'smooth_dotFiN:id,idgl:',id,idgl
        write(*,'(a12,3f8.3,1x,3f8.3,1x,f6.3)')'dx,nv,are: ',
     &  (dotcrd(k,idgl),k=1,3),
     &  (dotnvec(k,idgl),k=1,3),dotarea(idgl)
        end if
	end if
        end do

	   call clust_smooth_dot(ndot_mx,ndot_smooth,
     &                  dotcrd,dotnvec,dotarea,dot_type,
     &                  dot_smooth_gl,WDYONRm)

        if(testpr)then
        write(*,'(a58)')
     & 'dotListBeSmoothed: id, idg, dot_type,WDYONSm,npi,npj'
        if(ndot_smooth.le.0) goto 1000
        do id = 1, ndot_smooth
        idgl=dot_smooth_gl(id)
        npi=prpr1_smooth(1,id)
        npj=prpr1_smooth(2,id)
        write(*,'(20x,6i6,1x,L1)')
     &  id,idgl,dot_type(idgl),WDYONSm(1,idgl),
     &  npi,npj,WDYONRm(idgl)
        end do
		write(*,*)'in smooth_ccdot02c: FINISH'
        end if

1000     continue

	 if(testpr)then
	 write(*,*)
         write(*,'(a38,i6)')
     &  'SIMS:smooth_dotCC: Ndots are smoothed',nsmoothl
        write(*,'(a30,2f12.5)')'SIMS:Singular->Smoothed area:',
     &  aret,are1t
	 end if

         return
         end
       subroutine remove_ccdot02d(nprob_mx,ndot_mx,nconcij_mx,nconc_mx,
     &                  icfIJ,arc_cross,
     &                  RP,BIJ,UIJ,HIJ,rad_sm,dens,Aic,Ajc,radi,radj,
     &                  cprobe_ijk,nconc_glb_nprpos,probe_WYONRm,
     &                  cvert_probe,
     &                  concf_ijf,RotAng,SignRot,
     &                  start_dotn,stop_dotn,dotcrdl,WDYON,WDYONRm,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl)


       implicit none
       integer nprob_mx,ndot_mx,nconc_mx,nconcij_mx
       real*8 BIJ(3)
       real*8 UIJ(3)
       real*8 Aic(3),Ajc(3),radi,radj
       real*8 RP,HIJ,rad_sm,dens
       real*8 cprobe_ijk(3,nconc_mx)
       integer nconc_glb_nprpos(nconc_mx)
       logical*1 probe_WYONRm(nconc_mx)
       logical*1 WDYON(ndot_mx)
	   real*8 cvert_probe(3,3,nprob_mx)

       integer concf_ijf(nconcij_mx)
       integer SignRot(nconcij_mx)
       real*8 RotAng(nconcij_mx)

       integer icfIJ, arc_cross
       integer start_dotn(nprob_mx),stop_dotn(nprob_mx)
       real*8 dotcrdl(3,ndot_mx)

        logical*1 WDYONRm(ndot_mx)
        integer  WDYONSm(2,ndot_mx)
        integer ndot_smooth
        integer ndot_smmx
        integer prpr1_smooth(2,ndot_smmx)
        integer dot_smooth_gl(ndot_smmx)
    	integer nedgCppglb,nedgCpploc

        real*8 PI2,null,onehalf,one,two,three,four,big,small
        common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

        real*8 toler_in
		real*8 toler_cos
        real*8 rpr2,rsm2
        real*8 angsg
    	integer ic,jc,id,k,k1
    	integer icgl,jcgl
    	integer nst(2),nfin(2),nprb(2)
    	integer ipp,ip2,ip
    	real*8 aa,dist,dpp1,dpp12
    	real*8 sizecorr
    	real*8 r110,r11rsm,sint0,sintr,r110s
    	real*8 Hsmin,Hrem
    	real*8 toler_crossloc
    	real*8 xyzd(3),pxyz(3),p1xyz(3),p1pvec(3),p2xyz(3)
    	real*8 diu,dju,aidi,ajdj,aibij,ajbij,dist2
    	real*8 sibij,aisij,ajsij,HIJ2,rpps,rpps2,bijo
		real*8 vpi(3),vpj(3),vpk(3),vpijk(3)
		real*8 dotpxyz(3),viti(3),vjtj(3)
		real*8 sic(3),sjc(3)
		real*8 psiv(3),psjv(3),p1siv(3),p1sjv(3)
		real*8 npsi(3),npsj(3)
		real*8 sij2c(3),nppss(3),dotp2xyz(3)
		real*8 sij21c(3),dots2xyz(3)
		real*8 sim,sjm,dsim,dsjm,sijm,dppss
		real*8 plshift
		real*8 sicos,sjcos,sicosm,sjcosm
		real*8 vitim,vjtjm,ai,aj,bi
		real*8 DOT
        logical*1 testpr
        logical*1 inprzm


        testpr = .false.

        if(testpr)then
		write(*,*)
        write(*,*)'removeCCdot02d:start*******************'
	    end if

       toler_in = 0.1d-6
       toler_crossloc = 0.1d-6
	   plshift = 1.0d0
	   plshift = plshift*rad_sm
       toler_cos=0.15d0

       rpr2 = RP*RP
       rsm2=rad_sm**2
       sizecorr=0.50d0
       toler_in = - 0.5d0*sizecorr/dsqrt(dens)
       HIJ2=HIJ*HIJ
       rpps=RP + rad_sm
       rpps2=rpps**2

       aibij = dsqrt((RP+radi)**2 - HIJ2)
       ajbij = dsqrt((RP+radj)**2 - HIJ2)
       sibij = dsqrt(rpps2 - HIJ2)
       aisij = aibij - sibij
       ajsij = ajbij - sibij
       aa =  sibij*rad_sm/(rad_sm+RP)
       aidi = aisij + aa
       ajdj = ajsij + aa
	   vitim = RP*aisij/(RP+radi)
	   vjtjm = RP*ajsij/(RP+radj)
	   if(testpr)then
	   write(*,'(a20)')
     & 'in remove_ccdot02:'
	   write(*,'(a30,5f8.3)')'RP,HIJ,rad_sm,radi,radj:',
     & RP,HIJ,rad_sm,radi,radj
	   write(*,'(a45,4f7.3)')
     &	 'arccross=2:geom:aibij,ajbij,aisij,ajsij:',
     &	 aibij,ajbij,aisij,ajsij
	   end if



       ic=icfIJ/2
       if(ic*2 .ne. icfIJ)then
       write(*,*)'remove_ccdot02b: ERROR: icfIJ# 2*k (not even)'
       stop
       end if
       nedgCpploc = ic
       do ipp=1,nedgCpploc
	    ip2=2*ipp
        if(ipp.eq.1)then
    	ic=1
    	jc=icfIJ
    	else
        ic=ip2-2
    	jc=ic+1
        end if

    	nedgCppglb = nedgCppglb + 1
        icgl=concf_ijf(ic)
    	jcgl=concf_ijf(jc)

        if(testpr)then
		write(*,*)'* * * * Edge loop * * * * * * '
        write(*,*)'remove_cc02b:icfIJ: ic,jc = ',icfIJ,ic,jc
        write(*,*)'remove_cc02b: icgl,jcgl = ',icgl,jcgl
        end if

       nprb(1)=nconc_glb_nprpos(icgl)
       nst(1) = start_dotn(nprb(1))
       nfin(1)= stop_dotn(nprb(1))

       nprb(2)=nconc_glb_nprpos(jcgl)
       nst(2) = start_dotn(nprb(2))
       nfin(2)= stop_dotn(nprb(2))


       dpp12=0.0d0
       do k=1,3
       pxyz(k)=cprobe_ijk(k,icgl)
       p1xyz(k) = cprobe_ijk(k,jcgl)
       p1pvec(k) = p1xyz(k) - pxyz(k)
       p2xyz(k) = (p1xyz(k) + pxyz(k))*0.5d0
       dpp12=dpp12 + p1pvec(k)**2
       end do
       bijo = rsm2
       dpp1 = dsqrt(dpp12)
       dpp12 = dpp1*0.5d0

    	do k=1,3
    	p1pvec(k) = p1pvec(k)/dpp1
    	end do
    	Hsmin = dpp12*RP/(RP + rad_sm)

       	Hrem = dpp12 - toler_crossloc

       if(testpr)then
       write(*,'(a39,3f9.3)')'remove_ccdot02c:dpp12,Hsmin,Hrem:',
     & dpp12,Hsmin,Hrem
       end if

       do ip=1,2

       if(ip.eq.2)then
       do k=1,3
       pxyz(k)=cprobe_ijk(k,jcgl)
       p1xyz(k) = cprobe_ijk(k,icgl)
       p1pvec(k) = -p1pvec(k)
       end do
       end if

        if(testpr)then
		write(*,*)'* * * PROBE loop* * * * * * * * * *'
		write(*,'(a22,2i5)')'removeCCdot02d:ipp,ip;',ipp,ip
		write(*,'(a5,3f7.3,a5,3f7.3,a5,3f7.3)')
     &  'Aic:',Aic,' Ajc:',Ajc, ' UIJ:',UIJ
        write(*,'(a10,3f8.3,1x,2(a5,3f8.3))')
     &	'probexyz:p',pxyz, ' p1:',p1xyz,' pp2:',p2xyz
        end if
	   if(arc_cross.eq.2)then

	   if(testpr)then
		ai=RP/(RP+radi)
		aj=RP/(RP+radj)
		bi=RP/(RP+rad_sm)
       do k=1,3
	   vpi(k)=cvert_probe(k,1,nprb(ip))
	   vpj(k)=cvert_probe(k,2,nprb(ip))
	   vpk(k)=cvert_probe(k,3,nprb(ip))
	   vpijk(k)=vpi(k)+vpj(k)+vpk(k)
		vpi(k)= (Aic(k)-pxyz(k))*ai
		viti(k)=(Aic(k)+aisij*UIJ(k)-pxyz(k))*bi
		vjtj(k)=(Ajc(k)-ajsij*UIJ(k)-pxyz(k))*bi
		vpj(k)= (Ajc(k)-pxyz(k))*aj
		vpk(k)=vpijk(k)-vpi(k)-vpj(k)
		vpi(k)=vpi(k)+pxyz(k)
		vpj(k)=vpj(k)+pxyz(k)
		vpk(k)=vpk(k)+pxyz(k)
		viti(k) = viti(k) + pxyz(k)
		vjtj(k) = vjtj(k) + pxyz(k)
	   end do

		 write(*,'(a35,3f8.3)')':arc_cross2:'
		 write(*,'(3(a8,3f8.3))')'vpi:',vpi,' vpj:',
     &	 vpj, ' vpk:',vpk
		 write(*,'(2(a8,3f8.3))')'viti:',viti,' vjtj:',vjtj
		 end if

		dppss=0.0d0
        sicos=0.0d0
		sim=0.0d0
		do k=1,3
		sic(k)= Aic(k)+aisij*UIJ(k)
		sjc(k)= Ajc(k)-ajsij*UIJ(k)
		sij2c(k)=0.5d0*(sic(k)+sjc(k))
		psiv(k)=pxyz(k)-sic(k)
		psjv(k)=pxyz(k)-sjc(k)
		p1siv(k)=p1xyz(k)-sic(k)
		p1sjv(k)=p1xyz(k)-sjc(k)
		nppss(k)=p2xyz(k) - sij2c(k)
        dppss=dppss + nppss(k)**2
        sicos = sicos + (pxyz(k)-sic(k))*UIJ(k)
		sim = sim + (pxyz(k)-sic(k))**2
		end do
		sicosm = sicos/dsqrt(sim)
		sicosm = sicosm + (1.0d0-sicosm)*toler_cos
		sjcosm = sicosm
		dppss = dsqrt(dppss)

		call CROSS(p1siv,psiv,npsi)
		call CROSS(psjv,p1sjv,npsj)

		sim = DOT(npsi,UIJ)
		sjm = DOT(npsj,UIJ)

         if(testpr)then
		 write(*,'(a20,3f8.3)')'***smooth probes:'
		 write(*,'(2(a8,3f8.3))')
     &   'sic:',sic,' sjc:',sjc,' sij2:',sij2c
		 write(*,'(2(a8,3f8.3),a18,2f6.2)')
     &	 ' npsi0:',npsi, ' npsj0:',npsj,
     &   ' signi, signj:', sim,sjm
		 write(*,'(a16,2f8.3)')'sicosm, sjcosm:',sicosm,sjcosm
		 end if

		if(sim.lt.0.0d0)then
		do k=1,3
		npsi(k)=-npsi(k)
		end do
		sim = -sim
		end if
        if(sjm.gt.0.0d0)then
		do k=1,3
        npsj(k)=-npsj(k)
		end do
		sjm = - sjm
		end if

		ai=0.0d0
		aj=0.0d0
		do k = 1,3
		ai=ai + npsi(k)**2
		aj=aj + npsj(k)**2
		end do
		ai = dsqrt(ai)
		aj = dsqrt(aj)

        do k=1,3
		npsi(k)=npsi(k)/ai
		npsj(k)=npsj(k)/aj
		nppss(k)=nppss(k)/dppss
		sij21c(k)=sij2c(k)-nppss(k)*plshift
		end do

         if(testpr)then
		 write(*,'(2(a8,3f8.3),/,a18,2f6.2)')
     &		 ' npsi:',npsi, ' npsj:',npsj,
     &   ' signi, signj:', sim,sjm
		 write(*,'(a10,f8.3,a8,3f8.3)')
     &	 'dppss:', dppss, ' nppss:', nppss
		 write(*,'(a10,3f8.3)')'sij21c:',sij21c
		 end if

	     end if

       if((nst(ip).ge.1).and.(nfin(ip).ge.1))then
       do id=nst(ip),nfin(ip)
        dist = 0.0d0
        do k= 1,3
		dotpxyz(k)=dotcrdl(k,id)-pxyz(k)
        dist = dist + p1pvec(k)*dotpxyz(k)
        end do

        if(testpr)then
        write(*,'(a20,i5,1x,a6,3f8.4,1x,a5,L1,a5,2i2)')
     &	'removeCCd02d(IN):id:',id,
     &  'dotcrd:',(dotcrdl(k,id),k=1,3),
     &  ' Wrm:',WDYONRm(id),' Wsm:',WDYONSm(1,id),WDYONSm(2,id)
		 write(*,'(a12,3f8.3,a8,f8.3)')
     &  'dotpxyz:',dotpxyz, ' dist:',dist
		 write(*,'(2(a12,3f8.3))')
     &  'dotp2xyz:',dotp2xyz, ' dots2xyz:', dots2xyz

         end if

        WDYONSm(2,id) = arc_cross

        if(dist.GE.Hrem)then
        if((WDYON(id).and.(.not.WDYONRm(id))))then

        WDYONRm(id)= .true.

        if(testpr)then
        write(*,'(a28,i5,a5,L1,a5,2i2)')
     &	'removeCCd02d:Remove0:id:',id,
     &  ' Wrm:',WDYONRm(id),' Wsm:',WDYONSm(1,id),WDYONSm(2,id)
        end if

	    end if

	    else
	    if((dist .GT. Hsmin).and.(rad_sm.GT.0.0d0))then

        WDYONSm(1,id) = WDYONSm(1,id) + 1
		ndot_smooth = ndot_smooth +1

        if(ndot_smooth.gt.ndot_smmx)then
        write(*,*)'ERROR:ndot_smoothmx is SMALL'
        stop
        end if

        dot_smooth_gl(ndot_smooth) = id
        prpr1_smooth(1,ndot_smooth) =  icgl
        prpr1_smooth(2,ndot_smooth) =  jcgl

        if(testpr)then
        write(*,'(a25,i5,1x,a5,L1,a5,2i2)')
     &	'removeCCd02d:Smooth-dot:id:',id,
     &  ' Wrm:',WDYONRm(id),' Wsm:',WDYONSm(1,id),WDYONSm(2,id)
        end if

		end if
        end if

        if(arc_cross.eq.2)then
		sim=0.0d0
		sjm=0.0d0
		sijm=0.0d0
		dsim=0.0d0
		dsjm=0.0d0
		sicos=0.0d0
		sjcos=0.0d0
    	do k=1,3
		dotp2xyz(k)=dotcrdl(k,id)-p2xyz(k)
		dots2xyz(k)=dotcrdl(k,id)-sij21c(k)
        sim=sim + dotp2xyz(k)*npsi(k)
		sjm=sjm + dotp2xyz(k)*npsj(k)
		sijm=sijm+dots2xyz(k)*nppss(k)
		dsim = dsim +(dotcrdl(k,id)-sic(k))**2
		sicos=sicos + (dotcrdl(k,id)-sic(k))*UIJ(k)
		dsjm = dsjm + (dotcrdl(k,id)-sjc(k))**2
		sjcos=sjcos - (dotcrdl(k,id)-sjc(k))*UIJ(k)
		end do
		sicos = sicos/dsqrt(sim)
		sjcos = sjcos/dsqrt(sjm)

		 if(testpr)then
		 write(*,'(a17,3f8.3)')' sim, sjm, sijm:', sim,sjm,sijm
		 write(*,'(a17,2f8.3)')'sicos, sjcos:',sicos, sjcos
		 end if


		if(sim.ge.toler_in .and. sjm.ge.toler_in
     &	     .and. sijm .ge. toler_in)then
        if(sjcos.ge.sjcosm .and. sicos.ge.sicosm)then

		WDYONRm(id)=.true.

		if(testpr)then
        write(*,'(a25,i5,1x,a5,L1,a5,2i2)')
     &	'removeCCd02d:REmove2:id:',id,
     &  ' Wrm:',WDYONRm(id),' Wsm:',WDYONSm(1,id),WDYONSm(2,id)
		 write(*,*)
        end if

		end if
	    end if
		end if

       end do
       end if
       end do
       end do

        if(testpr)then
        write(*,'(a40,i5)')
     &	'removeCCdot02d:END:ndot_smooth: ', ndot_smooth
	    write(*,'(a30)')'* * * * * * * * *  * * * * * * * * '
        end if

       return
	   end
       subroutine arc_points02c(IATOM,JATOM,
     & ri,rj,dij,hij,rpr,dens,rad_sms,
     & putsm,putsm_smp,narc,ni,nj,
     & v0arc,stor,arc_cross,free_tor,
     & nijs_face,RotAng,smp2xyz,
     & dot_xyz,dot_vn,dot_area,dot_type,dot_atom,ndot)


       implicit none
       integer narcmx,nrotmx,nbr_mx,ndotmx

       parameter(narcmx=100)
       parameter(nrotmx=50 )
       parameter(ndotmx=narcmx*nrotmx)
       parameter(nbr_mx = 180)

       integer putsm

       integer putsm_smp

       integer IATOM,JATOM,nijs_face
       integer narc,ni,nj
       real*8 RotAng(nbr_mx)
       real*8 v0arc(3,narcmx)
       real*8 ri,rj,dij,hij,rpr,dens
       real*8 rad_sms,rad_sm
       real*8 stor
       logical free_tor
       integer arc_cross

       integer ndot
	   real*8 smp2xyz(3,2)
       real*8 dot_xyz(3,ndotmx),dot_vn(3,ndotmx),dot_area(ndotmx)
       integer dot_type(ndotmx),dot_atom(ndotmx)

       real*8 PI2,null,onehalf,one,two,three,four,big,small
       common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross

       real*8 teti,tetj,dtet2
       real*8 dens2,dens2_p,dens2_r,dens2_sm
       real*8 di,dj,di2,dj2,tet,tetmx,dtet
       real*8 pi,cost,cosdt,rip,rjp,aa
       real*8 step,bb,rpr_s
       real*8 dfi,fi_rot,start_ang
       real*8 ditet_i,ditet_j
       real*8 sinditet_i,sinditet_j
       real*8 rhditet_i,rhditet_j
       real*8 sinditet_s,ditet_s
       real*8 d_arpi,d_arpj,d_arps
       real*8 costi,costj,tetai,tetaj,rpsm,costsm
       real*8 tetasm,sintsm,pi05
       real*8 psi_xyz(3),psj_xyz(3)
       real*8 stor_i,stor_j,smsph_i,stor_ijs
       real*8 area_i,area_j, area_si, area_sj
       real*8 arc_tet(narcmx)
       real*8 rot_ang,dsint,dcost,dsint2
       real*8 drot_max
       real*8 rpr2,rpr22,rad_sm2,rprh
       real*8 toler_crossloc
       real*8 start_ang0,rot_ang0

       integer arc_type(narcmx)
       integer sadl_ptype,sadl_smtypei,sadl_smtypej

       logical changeSg
       logical testpr
       logical one_side
       integer ni_s,nj_s,ndotl
       integer i,j,k,n
       integer n1,n2,n3
       integer isf,nsf,nr
       integer done34

       testpr = .false.
       one_side=.false.

        pi=PI2

        sadl_ptype = 2
        sadl_smtypei = 4
		sadl_smtypej = 4
		if(toler_cross.gt.0.0d0)sadl_smtypej = 5
    	done34 = 0

	   toler_crossloc = toler_cross

        drot_max =  pi/1.5d0

        dens2=dsqrt(dens)
        dens2_p=rpr*dens2
        rpr2 = rpr*rpr
        rpr22=rpr2*2.0d0
        rprh = rpr*hij
cct
        rad_sm = rad_sms
        if(rad_sm.eq.0.0d0)rad_sm = 0.1d-9
        rad_sm2=rad_sm*rad_sm

	rip=ri+rpr
	rjp=rj+rpr
	stor=null


	if(ri.gt.rj)then
	cost=hij/rip
	changeSg=.false.
	else
	cost=hij/rjp
	changeSg=.true.
	end if

	tetmx=dacos(cost)

	cosdt=(rip**2 + rjp**2 - dij**2)/(2.0d0*rip*rjp)
	dtet=dacos(cosdt)

	if(testpr)then
	write(*,*)'arc_points:ATOMIJ:',IATOM,JATOM
	write(*,*)'arc_points: tetmx, dtet=',tetmx,dtet
        end if

        arc_cross=0
	if(hij.lt.rpr + toler_crossloc)arc_cross=2
	if(hij.le.rpr)arc_cross=1
	if(hij.le.rpr + toler_crossloc)arc_cross=1

	 if(tetmx.ge.dtet)one_side=.true.
	 if(one_side)arc_cross=0
ctt
        if(testpr)then
        write(*,*)'arc_points: hij,toler_cross:',
     &  hij,toler_crossloc
	    write(*,*)'one_side:',one_side,' arc_cross:',arc_cross
        end if

	dtet2=dtet*0.5d0
	ni=0
	nj=0
	if(putsm.eq.2)then
	ni=DINT(dtet2*dens2_p+0.5d0)
	if(ni.eq.0)ni=1
	narc=2*ni
	nj=ni

	else

	narc=DINT(dtet*dens2_p+0.5d0)
        if(narc.eq.0.and.putsm.eq.1)narc=1
	nj=narc/2
	ni=narc-nj
	if(changeSg)then
	ni=narc/2
	nj=narc-ni
	end if

	end if

	if(narc.lt.1)goto 1001

         if(narc.gt.narcmx)then
         write(*,*)'ERROR : sub: nrc_points , narcmx is small=',narcmx
	 stop
         end if

	di=dtet/dfloat(narc)
	di2=di*0.5d0

        ditet_i=di2
        ditet_j=di2
        rhditet_i=rpr*hij*2.0d0*ditet_i
        rhditet_j=rhditet_i
        sinditet_i=dsin(ditet_i)
        sinditet_j=sinditet_i
        ditet_s=0.0d0

	tet=tetmx+di2

	if(changeSg)then
	tet=tetmx - dtet
	tet=tet-di2
	di=-di
	end if

	if(arc_cross.eq.0)then
        rpr_s = rpr
        if(changeSg)rpr_s = -rpr
	    do i=1,narc
	    tet=tet - di
	    v0arc(1,i)= -rpr_s*dsin(tet)
	    v0arc(2,i)= hij - rpr*dcos(tet)
        v0arc(3,i)= 0.0d0
        if(i.le.ni)then
        arc_type(i)=1
        else
        arc_type(i)=2
        end if
        arc_tet(i) = tet
	    end do

        dsint2=dsin(tetmx-di*ni)
        stor_i=rhditet_i*ni-rpr2*(dsin(tetmx)-dsint2)
        stor_j=rhditet_j*nj-rpr2*(dsint2 -dsin(tetmx-dtet))

        if(changeSg)then
        di = -di
        dsint2=dsin(tetmx-narc*di+ni*di)
        stor_i=rhditet_i*ni-rpr2*(dsint2 - dsin(tetmx-narc*di))
        stor_j=rhditet_j*nj-rpr2*(dsin(tetmx) - dsint2)
        end if

        stor_ijs=stor_i + stor_j
        stor=stor_ijs
ctt
          if(testpr)then
          write(*,*)'torus: arc_cross=',arc_cross
          write(*,*)'putsmall = ',putsm
          write(*,*)'dot_density=',dens
          write(*,'(a50,3i5)')
     &    'arc_points: on v0arc: narc,ni,nj=',
     &    narc,ni,nj
          write(*,*)'nn    x       y       z'
          do i=1,narc
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          write(*,*)'analitical area(1radrot):stor-i,j=',stor_i,stor_j
          end if

        goto 1001
	end if


caug9831          if(arc_cross)then
	      if(arc_cross.gt.0)then
          costi = hij/rip
          costj = hij/rjp
          tetai = dacos(costi)
          tetaj = dacos(costj)

          rpsm = rpr + rad_sm
          costsm = hij/rpsm
          tetasm = dacos(costsm)
          sintsm = dsin(tetasm)

          do k=1,3
          psi_xyz(k) = 0.0d0
          psj_xyz(k) = 0.0d0
          end do

          psj_xyz(1) = rpsm*sintsm
          psi_xyz(1) = -psj_xyz(1)

		  do k=1,3
          smp2xyz(k,1)=psi_xyz(k)
		  smp2xyz(k,2)=psj_xyz(k)
		  end do

ctt
          if(testpr)then
          write(*,*)'hij,rad_smooth probe=',hij,rad_sm
          write(*,'(a30,3f10.6)')
     &    'costi,costj,costsm=',costi,costj,costsm
          write(*,'(a45,2f8.4)')
     &   'arc_points:singular arc: psi,psj_xyz',psi_xyz(1),psj_xyz(1)
          end if

          pi05=pi*0.5d0
          dens2_sm = dens2*rad_sm

          ni=DINT((tetai-tetasm)*dens2_p+0.5d0)
          if(ni.eq.0.and.putsm.eq.1)ni=1

          ni_s=DINT((pi05-tetasm)*dens2_sm  + 0.5d0)
          if(ni_s.eq.0.and.putsm_smp.eq.1)ni_s=1

          nj=DINT((tetaj-tetasm)*dens2_p+0.5d0)
          if(nj.eq.0.and.putsm.eq.1)nj=1
          nj_s=ni_s

          di = (tetai-tetasm)/dfloat(ni)
          di2=0.5*di
          ditet_i=di2
          rhditet_i=rpr*hij*2.0d0*ditet_i
          sinditet_i=dsin(ditet_i)

          tet=tetai+di2
          if(ni.ge.1)then
          do i = 1, ni
          tet = tet - di
          v0arc(1,i)= -rpr*dsin(tet)
          v0arc(2,i)= hij - rpr*dcos(tet)
          v0arc(3,i)= 0.0d0
          arc_tet(i) = tet
          arc_type(i)=1
          end do
          end if

          di = (pi05-tetasm)/dfloat(ni_s)
          di2=0.5*di
          ditet_s=di2
          sinditet_s=2.0d0*dsin(ditet_s)*rad_sm*rad_sm
          tet=pi05-tetasm + di2
          if(ni_s.ge.1)then
          do i = 1+ni, ni+ni_s
          tet = tet - di
          v0arc(1,i) = rad_sm*dcos(tet) + psi_xyz(1)
          v0arc(2,i)= rad_sm*dsin(tet)
          v0arc(3,i)= 0.0d0
          arc_type(i)=3
          arc_tet(i) = tet
          end do
          end if

          ni = ni + ni_s

          if(testpr)then
          write(*,*)'arc_points:singular arc:v0arc-iat: nisp,ni_sm=',
     &    (ni-ni_s),ni_s
          write(*,*)'nn    x       y       z'
          do i=1,ni
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          end if

          di = (tetaj-tetasm)/dfloat(nj)
          di2=0.5*di
          ditet_j=di2
          rhditet_j=rpr*hij*2.0d0*ditet_j
          sinditet_j=dsin(ditet_j)
          tet=tetaj+di2

          if(nj.ge.1)then
          do i = ni+1, nj+ni
          tet = tet - di
          v0arc(1,i)= rpr*dsin(tet)
          v0arc(2,i)= hij - rpr*dcos(tet)
          v0arc(3,i)= 0.0d0
          arc_type(i)=2
          arc_tet(i) = tet
          end do
          end if
          di = (pi05-tetasm)/dfloat(ni_s)
          di2=0.5*di
          tet=pi05-tetasm + di2

          if(nj_s.ge.1)then
          do i = 1+ni+nj, ni+nj+nj_s
          tet = tet - di
          v0arc(1,i) = psj_xyz(1) - rad_sm*dcos(tet)
          v0arc(2,i) = rad_sm*dsin(tet)
          v0arc(3,i) = 0.0d0
          arc_type(i) = 4
          arc_tet(i) = tet
          end do
          end if

          nj = nj + ni_s

          if(testpr)then
          write(*,*)'arc_points:singular arc:v0arc-iat: nisp,ni_sm=',
     &    (ni-ni_s),ni_s
          write(*,*)'nn    x       y       z'
          do i=1+ni,nj+ni
          write(*,'(i5,3f8.3)')i,(v0arc(k,i),k=1,3)
          end do
          end if

          narc = ni + nj

         stor_i = rprh*(tetai-tetasm)-rpr2*(dsin(tetai)-sintsm)

         smsph_i=rad_sm2*(1.0d0 - sintsm)


         stor_j= rprh*(tetaj-tetasm)-rpr2*(dsin(tetaj)-sintsm)

         stor_ijs = stor_i + 2.0d0*smsph_i + stor_j

ctt
        if(testpr)then
        write(*,'(a60,3f8.5)')
     &  'singul arc:analyt area(per 1rad):stor_i,smsph_i,stor_j',
     &  stor_i,smsph_i,stor_j
        end if

	end if

1001    continue

        ndotl = 0
        stor = 0.0d0
        done34=0

        do isf = 1,nijs_face
        nsf=2*isf - 1
	    if(arc_cross.eq.2)done34=done34+1

        rot_ang0 = RotAng(nsf+1)-RotAng(nsf)
caug9831        rot_ang = RotAng(nsf+1)-RotAng(nsf)

ctt
        if(testpr)then
        write(*,'(a45,2i5,2x,f8.6)')'arc_point:nijs_face,isf,rot_ang0=',
     &  nijs_face,isf,rot_ang0
        end if

        stor = stor + stor_ijs*rot_ang0
cttw
        start_ang0 = RotAng(nsf) - RotAng(1)
caug9831          start_ang = RotAng(nsf) - RotAng(1)

        if(testpr)then
        write(*,'(a30,f8.4)')'arc_point: start_ang =',start_ang
        end if

        area_i=0.0d0
        area_j=0.0d0
        area_si=0.0d0
        area_sj=0.0d0


        do i = 1,narc

	 if(arc_cross.eq.2.and.
     &	 (arc_type(i).eq.3.or.arc_type(i).eq.4))then
	 if(done34.eq.1)then
	 rot_ang = two*pi
	 start_ang = 0.0d0
	 else
	 goto 2002
	 end if
	 else
	 rot_ang = rot_ang0
	 start_ang = start_ang0
	 end if

        dens2_r = v0arc(2,i)*dens2
        nr = DINT(rot_ang*dens2_r + 0.5d0)
        if(nr.eq.0.and.putsm.eq.1)nr=1

        dfi = rot_ang/dfloat(nr)

        if(dfi.ge.drot_max)then
        dfi=drot_max
        nr = DINT(rot_ang/dfi + 0.5d0)
        end if
ctt
        if(testpr)then
        write(*,*)'arc_point:i=',i,' numb_rot:nr=',nr
        end if
ctt

        if(nr.ge.1)then
         dfi = rot_ang/dfloat(nr)
         fi_rot = start_ang - dfi*0.5d0

         if(arc_type(i).eq.1)then

          d_arpi=dfi*(rhditet_i - rpr22*dcos(arc_tet(i))*sinditet_i)
         end if

         if(arc_type(i).eq.2)then

          d_arpj=dfi*(rhditet_j - rpr22*dcos(arc_tet(i))*sinditet_j)
         end if
         if(arc_type(i).eq.3.or.arc_type(i).eq.4)then
         d_arps=dfi*dsin(arc_tet(i))*sinditet_s
         end if

         do j= 1,nr
         fi_rot = fi_rot + dfi
         ndotl=ndotl + 1

         if(ndotl.gt.ndotmx)then
         write(*,*)'ERROR'
         stop
         end if

         dcost=dcos(fi_rot)
         dsint=dsin(fi_rot)

         dot_xyz(1,ndotl)=v0arc(1,i)
         dot_xyz(2,ndotl)=v0arc(2,i)*dcost
         dot_xyz(3,ndotl)=v0arc(2,i)*dsint

         if(arc_type(i).eq.1)then
         dot_vn(1,ndotl) = -v0arc(1,i)/rpr
         dot_vn(2,ndotl) = (hij - v0arc(2,i))*dcost/rpr
         dot_vn(3,ndotl) = (hij - v0arc(2,i))*dsint/rpr
         dot_area(ndotl) = d_arpi
         dot_type(ndotl) = sadl_ptype
         dot_atom(ndotl) = IATOM
         area_i=area_i + d_arpi
         end if

         if(arc_type(i).eq.2)then
         dot_vn(1,ndotl) = -v0arc(1,i)/rpr
         dot_vn(2,ndotl) = (hij - v0arc(2,i))*dcost/rpr
         dot_vn(3,ndotl) = (hij - v0arc(2,i))*dsint/rpr
         dot_area(ndotl) = d_arpj
         dot_type(ndotl) = sadl_ptype
         dot_atom(ndotl) = JATOM
         area_j=area_j + d_arpj
         end if
         if(arc_type(i).eq.3)then
         area_si=area_si+ d_arps
         do k = 1,3
         dot_vn(k,ndotl) = (dot_xyz(k,ndotl)-psi_xyz(k))/rad_sm
         end do
         dot_area(ndotl) =  d_arps
         dot_type(ndotl) = sadl_smtypei
         dot_atom(ndotl) = IATOM
         end if

         if(arc_type(i).eq.4)then
         area_sj=area_sj+ d_arps
         do k = 1,3
         dot_vn(k,ndotl) = (dot_xyz(k,ndotl)-psj_xyz(k))/rad_sm
         end do

         dot_area(ndotl) =  d_arps
         dot_type(ndotl) = sadl_smtypej
         dot_atom(ndotl) = JATOM
         end if

         if(testpr)then
	 write(*,'(a5,i3,a3,i5,a5,i3,a5,f6.3,1x,a4,3f6.3,a5,3f6.3)')
     &   'rotj:',j,'id=',ndotl,'type:',dot_type(ndotl),
     &   'area:',dot_area(ndotl),
     &   ' xyz',(dot_xyz(k,ndotl),k=1,3),
     &   'nvec=', (dot_vn(k,ndotl),k=1,3)
         end if

         end do
         end if
2002     continue
         end do

ctt
         if(testpr)then
         write(*,*)'arc_point: dot in loc xyz: saddle isf=',isf
         write(*,*)'compare areas(anal vs num) arci,arcj,arcsm'
         write(*,'(a18,2f12.8)')'stor_i,area_i',stor_i*rot_ang,area_i
         write(*,'(a18,2f12.8)')'stor_j,area_j',stor_j*rot_ang,area_j
         write(*,'(a35,3f12.8)')'smsphi*rot_ang,area_si,srea_sj',
     &   smsph_i*rot_ang,area_si ,area_sj
         end if

         end do

         ndot = ndotl

	return
	end
       subroutine remove_ccdot03b(nprob_mx,ndot_mx,nconc_mx,nprobpost,
     &                  RP,rad_sm,dens,
     &                  conc_gl_ijk,nconc_glb_nprpos,
     &                  cprobe_ijk,cvert_probe,probe_WYONRm,
     &                  start_dotn,stop_dotn,dotcrdl,WDYON,WDYONRm,
     &                  ndot_smmx,WDYONSm,ndot_smooth,
     &                  prpr1_smooth,dot_smooth_gl)

       implicit none
       integer nprob_mx,ndot_mx,nconc_mx,nprobpost
       real*8 RP,rad_sm,dens

       real*8 cprobe_ijk(3,nconc_mx)
       real*8 cvert_probe(3,3,nprob_mx)
	   integer conc_gl_ijk(6,nprob_mx)
       integer nconc_glb_nprpos(nprob_mx)

       logical*1 WDYON(ndot_mx)
       logical*1 probe_WYONRm(nconc_mx)

       integer start_dotn(nprob_mx),stop_dotn(nprob_mx)

       real*8 dotcrdl(3,ndot_mx)

        logical*1 WDYONRm(ndot_mx)
        integer WDYONSm(2,ndot_mx)
        integer ndot_smooth
        integer ndot_smmx
        integer prpr1_smooth(2,ndot_smmx)
        integer dot_smooth_gl(ndot_smmx)

        real*8 PI2,null,onehalf,one,two,three,four,big,small
        common/CONSTANTI/PI2,null,onehalf,one,two,three,four,big,small

         real*8  toler_nb,toler_pr,toler_yon,toler_d,toler_cx,toler_ovr
         real*8  toler_cross
         common/TOLERANCE/toler_nb,toler_pr,toler_yon,toler_d,
     &   toler_cx,toler_ovr,toler_cross
         real*8 DOT

        integer nprobmx_loc
        parameter(nprobmx_loc=1000)
        real*8 toler_in
	    integer i,j,id,k,ip
	    integer nst,nfin,npi,npj
		integer npip,npjp,nstor
        integer NY
        integer prob_gl(nprobmx_loc)
        integer nrem,nsmoo
        logical*1 YON(nprobmx_loc)
        logical*1 inp1,inp2,inp3
        real*8 PY(3,nprobmx_loc)

        real*8 vpi1(3),vpi2(3),vpi3(3)
        real*8 vpj1(3),vpj2(3),vpj3(3)
        real*8 upj1(3),upj2(3),upj3(3)
        real*8 upi1(3),upi2(3),upi3(3)
        real*8 dotip1(3)
	    real*8 dist,dist2,distij
        real*8 distm,distm2,distum,distum2
        real*8 distu1,distu2,distu3
        real*8 rpr2,dpp1,Hsmmin,Hsm,dpp12
	    real*8 xyzd(3),pxyz(3),p1xyz(3),p1pvec(3)
        logical*1 testpr,OPT_checktor

        testpr = .false.

        OPT_checktor = .false.
        nrem = 0
        nsmoo = 0
		nstor=0

       toler_in = 0.1d-6
       rpr2=RP**2
       distm = 2.0d0*(RP+toler_cross)
       distm2 = distm**2
       distum = RP + 2.0d0*toler_cross
       distum2=distum**2

       if(testpr)then
	   write(*,*)
	   write(*,*)'in remove_ccdot03b: start ******************'
       write(*,'(a48,6f8.3)')
     & 'remove_ccdot03b:RP,rad_sm,distm,distm2,distum,distum2=',
     & RP,rad_sm,distm,distm2,distum,distum2
       end if

       NY = 0
       do i = 1, nprobpost
       if(probe_WYONRm(i))then
       NY = NY + 1
       if(NY.gt.nprobmx_loc)then
       write(*,*)'remove_ccdot03:ERROR: nprobmx_loc too small'
       stop
       end if

       prob_gl(NY)=i
       YON(NY)=.true.
       do k=1,3
       PY(k,NY) = cprobe_ijk(k,i)
       end do

       if(testpr)then
           write(*,'(a25,i6,1x,3f10.5)')
     &     'DeepProb I, xyz =', i,(cprobe_ijk(k,i),k=1,3)
       end if

       end if
       end do

       if(testpr)then
       write(*,'(a40,i5,2x,a20,i5)')
     & 'remove_cc03b:totalNprobePos=',nprobpost,
     & ' found NYprobe=', NY
       end if

       if(NY.le.1)goto 1000


        do i = 1, NY-1
        npi=prob_gl(i)
        npip=nconc_glb_nprpos(npi)

        if(testpr)then
        write(*,*)'remove_cc03b: i,YON(i) = ', i,YON(i)
        end if

        if(.not.YON(i))goto 100
C:SIMS02:

        do k=1,3
         vpi1(k)=cvert_probe(k,1,npi)
         vpi2(k)=cvert_probe(k,2,npi)
         vpi3(k)=cvert_probe(k,3,npi)
        end do

csep9802:
	 if(testpr)then
         write(*,'(a30,3f8.3)')'remove_ccdot03b:vpi1:',vpi1
         write(*,'(a30,3f8.3)')'remove_ccdot03b:vpi2:',vpi2
         write(*,'(a30,3f8.3)')'remove_ccdot03b:vpi3:',vpi3
	 end if

        j=i
200     j=j+1

        if(j.GT.NY)goto 100

         npj=prob_gl(j)
         npjp=nconc_glb_nprpos(npj)
         nstor=0
		 if(OPT_checktor)then
		 do k=1,3
		 if(conc_gl_ijk(1,npip).eq.conc_gl_ijk(k,npjp).or.
     &      conc_gl_ijk(2,npip).eq.conc_gl_ijk(k,npjp).or.
     &      conc_gl_ijk(3,npip).eq.conc_gl_ijk(k,npjp))nstor=nstor +1
		 end do

		 if(testpr)then
         write(*,*)'remove_cc03b: npi,npj=',npi,npj
         write(*,'(a10,3f10.4)')'Pi xyz=', PY(1,i),PY(2,i),PY(3,i)
         write(*,'(a10,3f10.4)')'Pj xyz=', PY(1,j),PY(2,j),PY(3,j)
     	 write(*,*)'Probes(gl) has eqvIJK:',nstor
         end if

		 if(nstor.ge.2) goto 200
		 end if

         p1pvec(1) = PY(1,j) - PY(1,i)
         if(dabs(p1pvec(1)).ge.distm)goto 200

         p1pvec(2) = PY(2,j) - PY(2,i)
         if(dabs(p1pvec(2)).ge.distm)goto 200

         p1pvec(3) = PY(3,j) - PY(3,i)
         if(dabs(p1pvec(3)).ge.distm)goto 200

         distij = p1pvec(1)**2 + p1pvec(2)**2 + p1pvec(3)**2
         if(distij.ge.distm2)goto 200

         if(testpr)then
         distij = dsqrt(distij)
         write(*,*)'remove_cc03b: probes can have intersection:'
		 write(*,*)'pri,prj,distij,distmx:', i,j,distij,distm
	     write(*,'(a10,3f8.3)')'p1pvec():',p1pvec
         end if

         do k=1,3
         vpj1(k)=cvert_probe(k,1,npj)
         vpj2(k)=cvert_probe(k,2,npj)
         vpj3(k)=cvert_probe(k,3,npj)

         upj1(k)=cvert_probe(k,1,npj)+p1pvec(k)
         upj2(k)=cvert_probe(k,2,npj)+p1pvec(k)
         upj3(k)=cvert_probe(k,3,npj)+p1pvec(k)
         end do

	     if(testpr)then
         write(*,'(a30,3f8.3)')'remove_ccdot03b:upj1:',upj1
         write(*,'(a30,3f8.3)')'remove_ccdot03b:upj2:',upj2
         write(*,'(a30,3f8.3)')'remove_ccdot03b:upj3:',upj3
	     end if

         call inprizm(vpi1,vpi2,vpi3,upj1,inp1)

         if(.not.inp1)then

	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vj1 is OUT:'
	  write(*,'(a20,f8.4,2x,L1)')'distu1, inp1 = ',distu1, inp1
	  end if

           call inprizm(vpi1,vpi2,vpi3,upj2,inp2)

            if(.not.inp2)then

	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vj2 is OUT:'
	  write(*,'(a20,f8.4,2x,L1)')'distu2, inp2 = ',distu2, inp2
	  end if

             call inprizm(vpi1,vpi2,vpi3,upj3,inp3)

             if(.not.inp3)then

	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vj3 is OUT:'
	  write(*,'(a20,f8.4,2x,L1)')'distu3, inp3 = ',distu3, inp3
	  end if

         call inprizm(upj1,upj2,upj3,vpi1,inp1)
         if(.not.inp1)then
	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vi1 is OUT:'
	  end if

           call inprizm(upj1,upj2,upj3,vpi2,inp2)
           if(.not.inp2)then
	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vi2 is OUT:'
	  end if

            call inprizm(upj1,upj2,upj3,vpi3,inp3)
             if(.not.inp3)then
	  if(testpr)then
	  write(*,'(a30)')'remove_ccdot03b: vi3 is OUT:'
	  end if
	     goto 200

	   end if
         end if
          end if
           end if
            end if
             end if

		   if(testpr)then
           write(*,*)'remove_cc03b: probe i,j have intersected tringles'
           end if

		do ip=1,2

		if(ip.eq.1)then
        nst = start_dotn(npi)
        nfin= stop_dotn(npi)
		do k=1,3
        pxyz(k) = PY(k,i)
        p1xyz(k) = PY(k,j)
		end do
		end if

		if(ip.eq.2)then
        nst = start_dotn(npj)
        nfin= stop_dotn(npj)
		do k=1,3
        pxyz(k) = PY(k,j)
        p1xyz(k) = PY(k,i)
		end do
		end if

		do k=1,3
        p1pvec(k) = p1xyz(k) - pxyz(k)
        end do

        dpp1 = p1pvec(1)**2 + p1pvec(2)**2 + p1pvec(3)**2
        dpp1 = dsqrt(dpp1)
        dpp12 = dpp1*0.5d0
        Hsmmin = dpp12*RP/(RP + rad_sm)

        if((nst.ge.1).and.(nfin.ge.1))then
        do id=nst,nfin
        if(.not.WDYONRm(id))then
        do k=1,3
        xyzd(k)=dotcrdl(k,id)
        end do

        Hsm = 0.0d0
        do k= 1,3
        Hsm = Hsm +(xyzd(k)-pxyz(k))*p1pvec(k)
        end do
        Hsm = Hsm/dpp1

        if(rad_sm.GT.0.0d0)then
        if((Hsm .GT. Hsmmin ).and.(Hsm.LT.dpp12))then
        WDYONSm(1,id) =  WDYONSm(1,id) + 1
        ndot_smooth = ndot_smooth +1
        nsmoo = nsmoo + 1

        if(testpr)then
        write(*,'(a30,3f10.5)')
     & 'Smoothing if:Hsmmin<Hsm<dpp12 =', Hsmmin,Hsm,dpp12
             write(*,'(a20,i3,i5,2x,3f8.4,1x,a16,2i3,1x,L1)')
     &  'pi-dot:ip,id,xyz:',
     &  ip,id,(dotcrdl(k,id),k=1,3),
     &  ' WSm(2),WRm:',WDYONSm(1,id),WDYONSm(2,id),WDYONRm(id)
        end if

        if(ndot_smooth.gt.ndot_smmx)then
        write(*,*)'ERROR:ndot_smoothmx is SMALL'
        stop
        end if

        dot_smooth_gl(ndot_smooth) = id

		if(ip.eq.1)then
        prpr1_smooth(1,ndot_smooth) =  npi
        prpr1_smooth(2,ndot_smooth) =  npj
		end if
		if(ip.eq.2)then
        prpr1_smooth(1,ndot_smooth) =  npj
        prpr1_smooth(2,ndot_smooth) =  npi
        end if

        end if
        end if

csep9814        if(WDYON(id).and.(.not.WDYONRm(id)))then
        if(WDYON(id) .and. .not.WDYONRm(id))then

		if(Hsm.ge.dpp12)then
        WDYONRm(id)= .true.
        nrem = nrem + 1

        if(testpr)then
        write(*,'(a30,2f10.5)')
     & 'REmovimg if:Hsm>dpp12:', Hsm,dpp12
        write(*,'(a20,i3,i5,2x,3f8.4,1x,a16,2i3,1x,L1)')
     &	'pi-dot:ip,id,xyz:',
     &  ip,id,(dotcrdl(k,id),k=1,3),
     &  ' WSm(2),WRm:',WDYONSm(1,id),WDYONSm(2,id),WDYONRm(id)
        end if

       end if
       end if
	   end if
       end do
       end if
	   end do


       YON(j) = .false.
       if(j.lt.NY)goto 200

       YON(i) = .false.
100    end do

       if(testpr)then
       write(*,'(a48,2i6)')
     & 'SIMS:remove_cc03b:N dots removed;to be smoothed ',
     &  nrem,nsmoo
	   write(*,*)'in remove_ccdot03b: finish'
       end if

1000   return
       end
       subroutine remove_sdsmp2(ndot_mx,nsmp2mx,ndotsmp2mx,
     &                  nsmp2tot,ndotsmp2,rad_sm,
     &                  dotsmp2_dgl,dotsmp2_nsm,smp2txyz,
     &                  dotcrd,dotnvec,dotarea,dot_type,
     &                  WDYONRm)


       integer ndot_mx,nsmp2mx,ndotsmp2mx
	   integer nsmp2tot,ndotsmp2
	   integer dotsmp2_dgl(ndotsmp2mx)
	   integer dotsmp2_nsm(ndotsmp2mx)
	   real*8 rad_sm
       real*8 smp2txyz(3,nsmp2mx)
	   real*8 dotcrd(3,ndot_mx)
	   real*8 dotnvec(3,ndot_mx)
	   real*8 dotarea(ndot_mx)
	   integer dot_type(ndot_mx)
	   logical*1 WDYONRm(ndot_mx)

		integer is,js,i,j,k
		integer id,idgl,irem
		real*8 rsm2,drm,drm0,drm1,kdist
		real*8 OPT_dhole
		real*8 dist,d
		logical testpr

		testpr = .false.
		irem=0
        OPT_dhole = 0.5d0
        rsm2=((2.0d0 +OPT_dhole)*rad_sm)**2
		drm0=rad_sm**2
		drm1=((1.0d0+OPT_dhole)*rad_sm)**2
		kdist = (drm1-drm0)/rsm2
        if(testpr)then
		write(*,*)'remove_sdsmp2: nsmp2tot:',nsmp2tot
		write(*,*)'remove_sdsmp2: ndotsmp2:',ndotsmp2
		end if
		do id = 1,ndotsmp2
		idgl = dotsmp2_dgl(id)
		if(.not.WDYONRm(idgl))then
		is = dotsmp2_nsm(id)

		if(testpr)then
		write(*,'(a34,3i5,1x,a8,3f7.3)')
     &	'id,dotsmp2_nsm(id),dotsmp2_dgl(id):',
     &  id,dotsmp2_nsm(id),idgl,
     &  ' dotcrd:',(dotcrd(k,idgl),k=1,3)
		write(*,'(a15,3f8.3)')
     &	'smp2txyz(is):',(smp2txyz(k,is),k=1,3)
		end if

    	do js = 1,nsmp2tot
		if(is.ne.js)then

		if(testpr)then
		write(*,'(a15,i5,1x,3f8.3)')
     &	'js, smp2txyz(js):',js,(smp2txyz(k,js),k=1,3)
		end if

         d = 0.0d0
		 dist = 0.0d0
		 do k =1,3
		 d = d + (dotcrd(k,idgl)-smp2txyz(k,js))**2
		 dist =dist + (smp2txyz(k,is)-smp2txyz(k,js))**2
		 end do
		 if(dist.lt.rsm2)then
		 drm = drm0 + kdist*dist
		 else
		 drm = drm1
		 end if
		 if(d.lt.drm)then
		 WDYONRm(idgl)=.true.
		 irem = irem +1
		 if(testpr)write(*,*)'dot id,idgl removed:',id,idgl
		 goto 1002
		 end if

		 end if
		 end do
		 end if
1002     continue
		 end do

		 if(testpr)then
		 write(*,*)'remove_smp2: dotremooved:', irem
		 end if

	   call clust_smooth_dot(ndot_mx,ndotsmp2,
     &                  dotcrd,dotnvec,dotarea,dot_type,
     &                  dotsmp2_dgl,WDYONRm)

		 return
		 end
	   subroutine clust_smooth_dot(ndot_mx,ndot_smooth,
     &                  dotcrd,dotnvec,dotarea,dot_type,
     &                  dot_smooth_gl,WDYONRm)


	   include 'surf-sims.h'
       integer ndot_mx

       real*8 dotcrd(3,ndot_mx)
       real*8 dotnvec(3,ndot_mx)
       real*8 dotarea(ndot_mx)
       integer dot_type(ndot_mx)

        logical*1 WDYONRm(ndot_mx)

		integer ndot_smooth
        integer dot_smooth_gl(ndot_smooth)

        integer wdyonsm_n(maxdot)
    	integer indxdclose(ndot_smoothmx)
    	logical*1 dotCdist(ndot_smoothmx)
        real*8 dot_sm_loc(3,ndot_smoothmx)
        real*8 dot_nv_loc(3,ndot_smoothmx)
    	real*8 dot_are_loc(ndot_smoothmx)

    	integer id,k,idgl
    	integer jd,jdgl,cd,cdgl
    	integer ndclose

		real*8 aret,are1t,are1,are
        real*8 dist,distav,dd
    	real*8 kdotclose
        logical*1 testpr
        integer nsmoothl
    	logical OPT_clastering

csep98
        testpr = .false.
    	OPT_clastering = .true.
        nsmoothl = 0

		if(testpr)then
		write(*,*)'n sub:clast_smooth_dot: start:'
        end if
	if(OPT_clastering)then
	if(ndot_smooth.ge.2)then
	kdotclose=0.50d0
	kdotclose = 0.5d0*kdotclose
	do id = 1, ndot_smooth
	dotCdist(id)=.false.
	end do

	do id = 1, ndot_smooth-1
	idgl=dot_smooth_gl(id)
	if(.not.dotCdist(id).and. (.not.WDYONRm(idgl)))then

        if(testpr)then
    	write(*,*)
        write(*,'(a24,2i6)')'Close:idgl:',idgl
        write(*,'(a12,3f8.3,1x,3f8.3,1x,f6.3)')'dx,nv,are: ',
     &  (dotcrd(k,idgl),k=1,3),
     &  (dotnvec(k,idgl),k=1,3),dotarea(idgl)
        end if

	ndclose=0
	do k=1,ndot_smooth
        indxdclose(k)=0
	end do

	distav =0.0d0
	are1=dotarea(idgl)
	dotCdist(id)=.true.

    	do jd = id+1,ndot_smooth
    	jdgl = dot_smooth_gl(jd)
    	if(.not.dotCdist(jd).and.(idgl.ne.jdgl)
     &     .and. (.not.WDYONRm(jdgl)))then

    	dist=0.0d0
	    do k=1,3
    	dist = dist + (dotcrd(k,idgl)-dotcrd(k,jdgl))**2
    	end do
    	are=dotarea(idgl)+dotarea(jdgl)

	    if(dist.lt.kdotclose*are)then
        dotCdist(jd)=.true.
	    ndclose = ndclose+1
	    distav = distav + dist
    	are1 = are1 + dotarea(jdgl)
        indxdclose(ndclose)=jdgl

        if(testpr)then
        write(*,'(a24,2i6)')'CdotLst:jdgl,ndclose:',jdgl,ndclose
    	write(*,'(a24,2f8.3)')'k*are, dist:',kdotclose*are,dist
        write(*,'(a12,3f8.3,1x,3f8.3,1x,f6.3)')'dx,nv,are: ',
     &  (dotcrd(k,jdgl),k=1,3),
     &  (dotnvec(k,jdgl),k=1,3),dotarea(jdgl)
        end if

    	end if

    	end if
    	end do

    	if(ndclose.ge.1)then
    	are1 = are1/(ndclose+1)
        distav=dsqrt(distav/ndclose/are1)
    	dotarea(idgl) = are1*(1.0d0 + distav)

        do cd=1,ndclose
    	dd = cd +1.0d0
        cdgl = indxdclose(cd)
    	WDYONRm(cdgl)=.true.
        do k=1,3
        dotcrd(k,idgl)=dotcrd(k,idgl)
     &      + (dotcrd(k,cdgl) - dotcrd(k,idgl))/dd
        dotnvec(k,idgl)=dotnvec(k,idgl)
     &      + (dotnvec(k,cdgl) - dotnvec(k,idgl))/dd
        end do
    	end do

    	dd=dsqrt(dotnvec(1,idgl)**2+dotnvec(2,idgl)**2+
     &  dotnvec(3,idgl)**2)
        dotnvec(3,idgl)=dotnvec(3,idgl)/dd
        dotnvec(2,idgl)=dotnvec(2,idgl)/dd
        dotnvec(1,idgl)=dotnvec(1,idgl)/dd

    	end if

        if(testpr)then
        write(*,'(a24,2i6)')'FinAvClose:idgl:',idgl
        write(*,'(a12,3f8.3,1x,3f8.3,1x,f6.3)')'dx,nv,are: ',
     &  (dotcrd(k,idgl),k=1,3),
     &  (dotnvec(k,idgl),k=1,3),dotarea(idgl)
        end if

    	end if
        end do
    	end if
    	end if

		if(testpr)then
		write(*,*)'n sub:clast_smooth_dot: finish'
        end if

		 return
		 end
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
c assign atomic radii from dictionary file:radii
c
        subroutine assign_rq(opt_ass,
     &                       natt_s,atnam_s,rnam_s,atrad_s,atcrg_s)
c
        include 'surf-sims.h'
	include 'assign_qr.h'
c..................................................................
       integer iopt,irefq
       integer natt_s,rnum_s(maxatm)
       character*2 opt_ass
       character*3 atnam_s(maxatm)
       character*3 rnam_s(maxatm)
       real*8 atrad_s(maxatm),atcrg_s(maxatm)
       real*8 atxyz_s(3,maxatm)
       real*8 aaqqd,chrgv,rad
       real*4 aaqq
c...................................................................
       character*60 comment
       character*10 nxtlbl
	character*9 bclab(5)
	character*80 line,filnam
	character*6 head
	character*3 tres
	character*6 tatm
	character*24 crdstr
	character*13 radstr
	logical iper,iconc,isite,iatout,ibios,isph
	logical iauto,donon,dmmy2

        character*9 radfil,crgfil
        character*10 pdbfili,pdbfilo
        character*4 num_to_char
        logical CONTROL
	character*7 char7
        integer i,j,k,iun11,iun12,iun19
        integer nchrec
        integer irtot,ictot,nent
        integer nrdrec,nqass
	integer iiprint0,iiprint1
        real qnet,qmin,qplus
        integer n,ifind,natom
        character*80 liner,line1,line2,line3

        CONTROL = .false.

c ASSIGN q and atRad
c initialization
        if(opt_ass.eq.'RQ')then
	do i=1,natt_s
	atrad_s(i)=0.0d0
	atcrg_s(i)=0.0d0
	end do
        end if

        if(opt_ass.eq.'R ')then
        do i=1,natt_s
        atrad_s(i)=0.0d0
        end do
        end if

        if(opt_ass.eq.'Q ')then
	do i=1,natt_s
	atcrg_s(i)=0.0d0
	end do
        end if

c.........................................................................
      radfil= 'radii'
      crgfil= 'atcrg.crg'
c
c initialize rad, charge linklists
c
	do 9000 i = 1,nrlist
	  irnumb(i) = 0
	  irlink(i) = 0
9000	continue
	do 9001 i = 1,nclist
	  icnumb(i) = 0
	  iclink(i) = 0
9001	continue

        iatout=.false.
        if(iatout) then
        iun19=kanalpdbout
	open(unit=iun19,file='molec_rq.pdb',status='unknown')
        end if

c
c read radius file
c
        if(opt_ass.eq.'R '.or.opt_ass.eq.'RQ')then
        iun11=kanalrad
	open(unit=iun11,file=radfil,status='old',err=901)
cwrite(kanalp,*)'reading radii from file'
cwrite(kanalp,*)radfil
        irtot=0
c
105	read(iun11,201,end=901)comment
	if(comment(1:1).eq.'!') then
	  write(kanalp,*)comment
	  goto 105
	end if

c......read file radii.........................................

        nrdrec = 0
100   continue
	  nrdrec = nrdrec + 1
	  if(nrdrec.gt.nrmax) then
	    write(kanalp,*)'assingRQ:maximum # of radius records exceeded'
	    write(kanalp,*)' increase nrmax'
	    stop
	  end if
	  read(iun11,200,err=904,end=300)atm,res,rad
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  atnam(nrdrec) = atm
	  rnam(nrdrec)  = res
	  radt(nrdrec)  = rad
	  call rent(atm,res,nrdrec)
	goto 100
c..................................................................
300   continue
	close(iun11)
	nrdrec = nrdrec - 1
c 	write(kanalp,*)'# of radius parameter records:',nrdrec
        end if

c
c read charge parameter file
c
        if(opt_ass.eq.'Q '.or.opt_ass.eq.'RQ')then
        iun12=kanalq
	open(unit=iun12,file=crgfil,status='old',err=902)
	write(kanalp,*)'reading charges from file'
	write(kanalp,*)crgfil
        ictot=0
c
106	read(iun12,201,end=901)comment
	if(comment(1:1).eq.'!') then
	  write(kanalp,*)comment
	  goto 106
	end if
	nchrec = 0
c................read file atcrg.crg.....................
101   continue
	  nchrec = nchrec + 1
	  if(nchrec.gt.ncmax) then
	    write(kanalp,*)' maximum # of charge records exceeded'
	    write(kanalp,*)' - increase ncmax'
	    stop
	  end if
	  read(iun12,202,err=905,end=301)atm,res,rnum,chn,chrgv
	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
	  catnam(nchrec) = atm
	  crnam(nchrec)  = res
	  crnum(nchrec)  = rnum
	  cchn(nchrec)   = chn
	  chrgvt(nchrec) = chrgv
	  call cent(atm,res,rnum,chn,nchrec)
	goto 101
301   continue
	close(iun12)
c..................................................................
	nchrec = nchrec - 1
	write(kanalp,*)'# of charge parameter records:',nchrec

        end if
c
c atom coordinate file assign radii and charges
c
        natom = 0
	qnet = 0.
	nqass = 0
	qplus = 0.
	qmin = 0.
103   continue

c................. assign rad  and q to atoms  .................

         head = 'ATOM'
         call up(head,6)
c
	  natom = natom + 1
          if(natom.gt.natt_s) goto 303

	  if(natom.gt.natmx)then
	    print *,'Max # of atoms exceeded: increase natmx: ',natmx
	    stop
	  end if

          atm = atnam_s(natom)
          res = rnam_s(natom)

          chn = ' '
          line(1:40)= '                                         '
          line(41:80)='                                         '
          line(1:6)=head
          line(12:16)=atm
          line(18:20)=res
	  line(23:26)='    '
          line(22:22)=chn

	  call up(atm,6)
	  call elb(atm,6)
	  call up(res,3)
	  call elb(res,3)
	  call up(rnum,4)
	  call elb(rnum,4)
	  call up(chn,1)
	  call elb(chn,1)
c
        if(opt_ass.eq.'R '.or.opt_ass.eq.'RQ')then

	  call rfind(atm,res,ifind,n)
	  if(ifind.eq.0) then
          tres = '   '
	    call rfind(atm,tres,ifind,n)
	    if(ifind.eq.0) then
		tatm = atm(1:1)//'     '
	      call rfind(tatm,tres,ifind,n)
	      if(ifind.eq.0) then
	        write(kanalp,*)'no radius record found for'
	        write(kanalp,*)line
                write(kanalp,*)natom,atnam_s(natom),rnam_s(natom)
	        stop
		end if
	    end if
	  end if
	  atrad(natom) = radt(n)
	  atrad_s(natom) = dble(radt(n))

          end if
c
        if(opt_ass.eq.'Q '.or.opt_ass.eq.'RQ')then

	  call cfind(atm,res,rnum,chn,ifind,n)
	  if(ifind.eq.0) then
	    schn = chn
	    chn = ' '
	    call cfind(atm,res,rnum,chn,ifind,n)
	    if(ifind.eq.0) then
		chn = schn
		snum = rnum
		rnum = '    '
	      call cfind(atm,res,rnum,chn,ifind,n)
	      if(ifind.eq.0) then
	        schn = chn
	        chn = ' '
	        call cfind(atm,res,rnum,chn,ifind,n)
		  if(ifind.eq.0) then
		    chn = schn
		    rnum = snum
		    sres = res
		    res = '   '
	          call cfind(atm,res,rnum,chn,ifind,n)
	          if(ifind.eq.0) then
	            schn = chn
	            chn = ' '
	            call cfind(atm,res,rnum,chn,ifind,n)
	            if(ifind.eq.0) then
		        chn = schn
		        snum = rnum
		        rnum = '    '
	              call cfind(atm,res,rnum,chn,ifind,n)
	              if(ifind.eq.0) then
	                schn = chn
	                chn = ' '
	                call cfind(atm,res,rnum,chn,ifind,n)
			  endif
			endif
		    endif
		  end if
		end if
	    end if
	  end if
c
c net charge etc
c
	  if(ifind.ne.0) then
	    chrgv = dble(chrgvt(n))
	    nqass = nqass + 1
	  else
	    chrgv = 0.0d0
	  end if

	  qnet = qnet + chrgv

	  atcrg(natom) = chrgv
	  atcrg_s(natom) = chrgv

	  if(chrgv.gt.0.0d0) then
	    qplus = qplus + chrgv
	  end if
	  if(chrgv.lt.0.0d0) then
	    qmin = qmin + chrgv

	  end if

          end if
c
	  if(iatout) then
          write(radstr,206)atrad(natom),atcrg(natom)
          line(55:67) = radstr
          write(line(31:54),207)(atxyz_s(j,natom),j=1,3)
          write(iun19,204)line
	  end if

          if(natom.eq.natt_s)goto 303

      go to 103
303   continue
	if(iatout) then
	  close(iun19)
	end if
c
	goto 999
901	write(kanalp,*) 'unexpected end or non-existence of radius file'
	write(kanalp,*)radfil
   	stop
902	write(kanalp,*) 'unexpected end or non-existence of charge file'
   	stop
904	write(kanalp,*) 'error in reading radius file'
   	stop
905	write(kanalp,*) 'error in reading charge file'
   	stop
999	continue

1099    continue

c
200   format(A6,A3,D8.3)
201   format(a)
202   format(A6,A3,A4,A1,D8.3)
204   format(a80)
205   format(3f8.3)
206   format(F6.2,F7.3)
207   format(3f8.3)

1009        return

	end
c------------------------------------------------------------------
      subroutine cent(atm,res,rnum,chn,nent)
c
        include 'surf-sims.h'
	include "assign_qr.h"

        integer irtot,ictot,nent,n,new
        integer ichash
c
	if(ictot.eq.nclist) then
	  write(6,*)'charge list full- increase nclist'
	  stop
	end if

	n = ichash(atm,res,rnum,chn)
	if(icnumb(n).ne.0) then
c
9000	   continue
	   if(iclink(n).eq.0)goto 9001
            n = iclink(n)
	   goto 9000
9001	   continue
c
         new = 1
9002	   continue
	   if(icnumb(new).eq.0)goto 9003
            new = new + 1
    	   goto 9002
9003     continue

c
c found one- addlink
c
         iclink(n) = new
         n = new
      end if

c
      icnumb(n) = nent
      iclink(n) = 0
	ictot = ictot + 1
      return
      end

c--------------------------------------------------------
      subroutine cfind(atm,res,rnum,chn,ifind,n)
c
c	find entry nres in hash table and check match with res
c
        include 'surf-sims.h'
	include "assign_qr.h"
        integer irtot,ictot,nent,n,new,ifind
        integer ichash

      n = ichash(atm,res,rnum,chn)
c	check for match
      ifind = 0
100   continue
c	while no match and not at end of linklist
      if(icnumb(n).eq.0) then
        ifind = 0
        return
      end if

      if((res.eq.crnam(icnumb(n))).and.(atm.eq.catnam(icnumb(n)))
     &  .and.(rnum.eq.crnum(icnumb(n))).and.(chn.eq.cchn(icnumb(n))))
     &  then
        n = icnumb(n)
        ifind = 1
        return
      else
        if(iclink(n).ne.0) then
          n = iclink(n)
        else
          ifind = 0
          return
        end if
      end if
      go to 100
      end

c----------------------------------------------------------
      subroutine elb(txt,len)
c
c eliminate leading blanks from a character string
c
      character*(*) txt
      character*80 save
      integer i,len

      do 9000 i=1,len
        if(txt(i:i).ne.' ') then
          nonb = i
          go to 100
        end if
9000  continue
      return
100   continue
      save = txt(nonb:len)
      txt = save
      return
      end
c------------------------------------------------------------------
      subroutine rent(atm,res,nent)
c
c enter character strings res,atm into hash table for radii
c by assigning it a number nent
c
        include 'surf-sims.h'
	include 'assign_qr.h'
        integer irtot,ictot,nent,n,new
        integer irhash
c
c check to see if there is room
c
	if(irtot.eq.nrlist) then
	  write(6,*)' radii list full- increase nrlist'
	  stop
	end if

	n = irhash(atm,res)
	if(irnumb(n).ne.0) then
c
9000	continue
	if(irlink(n).eq.0)goto 9001
            n = irlink(n)
	goto 9000
9001	continue
c
         new = 1
9002	continue
	if(irnumb(new).eq.0)goto 9003
            new = new + 1
	goto 9002
9003	continue
c
         irlink(n) = new
         n = new
      end if

c
      irnumb(n) = nent
      irlink(n) = 0
	irtot = irtot + 1
      return
      end
c--------------------------------------------------------
      subroutine rfind(atm,res,ifind,n)
c
c	find entry res in radius hash table and check match
c
        include 'surf-sims.h'
	include 'assign_qr.h'
        integer irtot,ictot,nent,n,new,ifind
        integer irhash

      n = irhash(atm,res)
c	check for match
      ifind = 0
100   continue
c	while no match and not at end of linklist
      if(irnumb(n).eq.0) then
        ifind = 0
        return
      end if

      if((res.eq.rnam(irnumb(n))).and.(atm.eq.atnam(irnumb(n)))) then
        n = irnumb(n)
        ifind = 1
        return
      else
        if(irlink(n).ne.0) then
          n = irlink(n)
        else
          ifind = 0
          return
        end if
      end if
      go to 100
      end

c------------------------------------------------------------
      subroutine up(txt,len)
c
c convert character string to upper case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do 9000 i=1,len
        if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
9000	continue

      txt = save
      return
      end

c----------------------------------------------------------------
      function ichash(atxt,rtxt,ntxt,ctxt)
c
c produce hash number from atom and residue name,rsidue number and chain
c name for charge assignment
c
       include 'surf-sims.h'
       include 'assign_qr.h'

       character*6 atxt
       character*3 rtxt
       character*4 ntxt
       character*1 ctxt
       character*38 string
       integer irhash,ichash
       integer n,i,j
       integer irtot,ictot

       data string /'* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

       n = 1
       do 9000 i = 1,3
        j = index(string,rtxt(i:i))
        n = 5*n + j
9000	continue
      do 9001 i = 1,6
        j = index(string,atxt(i:i))
        n = 5*n + j
9001	continue
      do 9002 i = 1,4
        j = index(string,ntxt(i:i))
        n = 5*n + j
9002	continue
      do 9003 i = 1,1
        j = index(string,ctxt(i:i))
        n = 5*n + j
9003	continue
	n = abs(n)
      ichash = mod(n,nclist) + 1
      return
      end

c----------------------------------------------------------------
      function irhash(atxt,rtxt)
c
c produce hash number from atom and residue name
c for radius assignment
c
        include 'surf-sims.h'
	include 'assign_qr.h'

      character*6 atxt
      character*3 rtxt
      character*38 string
      integer n,i,j,ictot,irtot
      integer ichash,irhash
      data string /'* 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      n = 1
      do 9000 i = 1,3
        j = index(string,rtxt(i:i))
        n = 5*n + j
9000	continue
      do 9001 i = 1,6
        j = index(string,atxt(i:i))
        n = 5*n + j
9001	continue
	n = abs(n)
      irhash = mod(n,nrlist) + 1
      return
      end


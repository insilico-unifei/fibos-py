c**********************************************************************
c       program to calculate the exposed, buried and pocket surfaces
c       for each residue in a protein
c
c       author  : n.pattabiraman.
c
c       version : 1.0 dated :  apr-20-1992
c       version : 1.1 dated :  nov-30-1993
c
c      version 2: p. fleming, 1994
c      version 3: p. fleming, 1995
c      rewritten as version 4: p. fleming, 1997
c
C* Copyright (c) Laboratory for the Structure of Matter,
C* Naval Research Laboratory, Washington DC 20375 and Geo-Centers Inc
C* Fort Washington, MD 20744.
C* and
C* Yale University, New Haven, CT 06520
C* All rights reserved.
C
C  The program either in full or in part should not be distributed
C  in any form or by any means without the prior written
C  permission of the author:
C
C* N. Pattabirman
C* Bldg 430/ Room 206
C* PRI-NCI/FCRDC
C* P. O. Box B
C* Frederick, MD 21702-1201
C* E-mail : pattabir@fcrfv1.ncifcrf.gov
C
C* Patrick Fleming
C* Yale University
C* Dept. of Molecular Biophysics and Biochemistry
C* P. O. Box 208114
C* 260 Whitney Ave.
C* New Haven, CT 06520-8114
C* fleming@csb.yale.edu
C
c**********************************************************************
       subroutine surfcal
       parameter ( maxat=20000,maxseg=2000)

       common/resids/ nseg,secnm(maxseg),secsq(maxseg),atnm(maxat),
     1               iats(maxseg),iate(maxseg),seqcd(maxseg),natom

       common/coord/cr(maxat,3),secsqe(maxseg),rescid(maxseg)
       common/anput/resinf(5),residen
       common/vdw/vdwrad(maxat),inter,icter

       dimension attype(50),vdwr(50),restyp(50)

       character*3 attype,restyp
       character *1 resext,secsqe,rescid
       character*1 chainid
       character*13 atiden
       character *4 atnm,atname
       character*60 infile
       character*80 line
       character*13 resinf,residen
       character*9 dumyc1,predumy
       integer resno,secsq,numat,seqcd,first,rayflag
       integer dummy
       integer kanal
       character *3 resn,secnm


 10    format(6x,5i5)
770    format(' type in the name of the coordinates file '/)
880    format(' type in the name of the output file '/)
881    format(' type in the residue range '/)
110    format(a)
115    format(i1)
1989   format(' the number of atoms is larger than ',i5/)
2989   format(' the number of residue is larger than ',i5/)

c Now read "part.inp" for info on where this residue is
c in the sequence. This is necessary for handling N and
c C termini.
c The format of the "part.inp" file for first, middle
c and last residues is as follows for a 6 residue peptide:

c first residue
c 1
c
c ALA    1
c VAL    2    N
c part_v.pdb
c part_i.ms

c middle residue
c 1
c ALA    4    C
c ALA    4    O
c VAL    5
c ALA    6    N
c part_v.pdb
c part_i.ms

c last residue
c 1
c VAL    5    C
c VAL    5    O
c ALA    6
c        7    N
c part_v.pdb
c part_i.ms

c read the first line of "part.inp"
c if 0, no ray file, if 1, write ray display file
        print *,path
        open(unit=kanal,file='part.inp',status='old')

        read(kanal,"(a)")rayflag
        rayflag = 1

	i = 1		!counter for reading part.inp
c now read the second line. this will be
c   the res and atoms n-terminal to the residue
c   of concern

987 	continue	!return here every line
	read(kanal,110) resinf(i) 		!i.e., THR   22    C

	if(resinf(i)(1:5).eq.'     ') then
	   if(i.eq.1)then	!this is res 1 and no i-1 exists
	      i = i + 1
	      resinf(i) = resinf(i-1) !to account for C and 0
	   endif
	endif
	i = i + 1
	if(i.lt.5) go to 987	!go back up and read next line
                                !until reach line 5, which is last
                                !line of residue information.

        residen = resinf(3) 	!this is current residue

c check if residue i is n-ter or c-ter

	inter = 0
 	icter = 0
	if(resinf(1)(1:5).eq. '     ') inter = 1 !n-terminal
	if(resinf(4)(1:5).eq. '     ') icter = 1 !c-terminal

c display on screen residue being calculated
        write(6,57)resinf(3)
57      format(5x,a13)

c readin the vdw radii
	call assvdw(attype,vdwr,nvdwt,restyp)

c        do 1234 dummy=1,nvdwt
c        print *, attype(dummy),restyp(dummy),vdwr(dummy)
c1234    continue

1000   continue

c read the name of the coords file containing all the
c  rest of the protein besides the residue of concern
       read(kanal,110) infile	!should be part_v.pdb
				!the rest of the pdb file
        close(kanal)

c open that file as 1
       open(unit=1,file=infile, status='old')
       rewind 1

c open file.srf to write out the *.srf info for this res
       open(unit=2,file='file.srf',status='unknown', position='append')

c initialize
       nseg=0
       natom=0
       lresno=-9999
       first = 0

100    read(1,110) line 	!read part_v.pdb
c      if(line(1:4).eq.'TER ') then
c        first = 0
c      endif
       if(line(1:4).eq.'ATOM') then

c to check the format for atom names in pdb format
c the atom name starts at 13 or 14 column
c this moves atom name to 13:15

         if(line(13:13).eq.' ') then
           line(13:13) = line(14:14)
           line(14:14) = line(15:15)
           line(15:15) = line(16:16)
           line(16:16) = ' '
         endif
c
c      to check for chain id
c
c        if(first .eq. 0) then
c          first = 1
c          chainid = ' '
c          if(line(22:22) .ne. ' ') then
c               chainid = line(22:22)
c          endif
c        endif
c        if(first .eq. 1) then
c          if(line(22:22) .ne. chainid) then
c              print *,'error in chain id'
c              print *, line(1:22)
c              stop
c          endif
c        endif
c
c      to check for alternate conformation
c
         if(line(17:17) .ne. ' ')then
           print *, ' alternate conformation is present.  check '
           print *,line(1:22)
	   stop
         endif

         natom=natom+1
         if(natom.gt.maxat) then
c          write(6,1989) maxat
           stop
         endif
c
c      assigning vdw radius
c
         vdwrad(natom) = 0.0
         do 1256 ira = 1, nvdwt
c first the 3 letter atoms
           if(attype(ira)(3:3).ne.' ') then
             if(line(13:15).eq.attype(ira)(1:3)
     &.and.(restyp(ira).eq.'   ')) then
                vdwrad(natom) = vdwr(ira)
             elseif((line(13:15).eq.attype(ira)(1:3))
     &.and.(line(18:20).eq.restyp(ira)))then
                vdwrad(natom) = vdwr(ira)
                go to 1257
             endif

c second the 2 letter atoms
           elseif(attype(ira)(2:2).ne.' ') then
             if(line(13:14).eq.attype(ira)(1:2)
     &.and.(restyp(ira).eq.'   ')) then
                vdwrad(natom) = vdwr(ira)
             elseif((line(13:15).eq.attype(ira)(1:3))
     &.and.(line(18:20).eq.restyp(ira)))then
                vdwrad(natom) = vdwr(ira)
                go to 1257
             endif

c last the 1 letter atoms
           else if(attype(ira)(2:2).eq.' ')then
             if(line(13:13).eq.attype(ira)(1:1)) then
               vdwrad(natom) = vdwr(ira)
             endif
           endif
1256     continue
1257     continue
c        print *, line(13:15), vdwrad(natom)
         if(vdwrad(natom).le.0.001)then
           print *, ' Atom not found in radii file  ',natom
           stop
         endif
c
c  assign coords to array

       read(line,120) atnm(natom),resn,resno,resext,(cr(natom,l),l=1,3)
c      write(6,120) atnm(natom),resn,resno,resext,(cr(natom,l),l=1,3)
120    format(12x,a4,1x,a3,2x,i4,a1,3x,3f8.3)

         if(resno.gt.99) then
           if(resext.ne.' '.and.chainid.ne.' ')then
             print *, ' you may have to reorder the residus numbers'
             stop
           endif
         endif

         if(resno.eq.lresno) go to 876
         nseg=nseg+1
         if(nseg.gt.maxseg) then
c          write(6,2989) maxseg
           stop
         endif
         secsq(nseg)=resno	!residue number
         secnm(nseg)=resn	!residue name, "ALA"
         secsqe(nseg) = resext
         iats(nseg)=natom
         rescid(nseg) = chainid
         lresno=resno
 876     continue
         iate(nseg)=natom
       else if (line(1:3).eq.'END') then
         go to 101
       end if
       go to 100
 101   continue

       close(1)	!close part_v.pdb

c Call the main subroutine
       call bsurf(rayflag,attype,vdwr,nvdwt,restyp)
        close(2)
c       stop
       end
c
c****************************************************************

       subroutine asorder(nuniq,ncount,parea,dtmin,iuni)

c      this subroutine reverses the order of array

       dimension nuniq(1200),ncount(1200),parea(1200),dtmin(1200)

	do 100 i = 1,iuni-1
	  do 101 j = i+1,iuni
	    if(ncount(i).le.ncount(j))then
              ncl  = ncount(i)
	      ncount(i) = ncount(j)
	      ncount(j) = ncl
	      nul = nuniq(i)
	      nuniq(i) = nuniq(j)
	      nuniq(j) = nul
	      area = parea(i)
	      parea(i) = parea(j)
	      parea(j) = area
	      dtm = dtmin(i)
	      dtmin(i) = dtmin(j)
	      dtmin(j) = dtm
	    endif
101	  continue
100	continue
c
	return
	end
c
c*****************************************************************

       subroutine bsurf(rayflag,attype,vdwr,nvdwt,restyp)

       parameter ( maxat=20000,maxseg=2000)
       common/resids/ nseg,secnm(maxseg),secsq(maxseg),atnm(maxat),
     1               iats(maxseg),iate(maxseg),seqcd(maxseg),natom

       common/coord/cr(maxat,3),secsqe(maxseg),rescid(maxseg)
       common/anput/resinf(5),residen
       common/vdw/vdwrad(maxat),inter,icter

       dimension attype(50),vdwr(50),restyp(50)
       character*3 attype,restyp

       dimension ipick(1200),retropt(1200)
       dimension dminpt(1200),iadmin(1200)
       dimension nuniq(1200),ncount(1200)
       dimension ptarea(1200),parea(1200),dtmin(1200)
       character*80 line,sline
       character *1 resext,secsqe,rescid

       character*13 atiden
       character *4 atnm,atname,dumyat,fratom,toatom
       character*60 infile
       character*13 resinf,residen
       integer secsq,seqcd,rayflag
       character *3 secnm

       dimension spc(3),ptc(3)

c open the residue surface file to read

c       read(5,110) infile 	!part_i.ms
        infile = "part_i.ms"
110    format(a)
       open(unit=3,file=infile, status='old')
       rewind 3		!part_i.ms

       open(7,file='buried.ms',status='unknown')
       open(8,file='exposed.ms',status='unknown')
       open(10,file='scratch.ms',status='unknown')
       if (rayflag .eq. 1) then
         open(11,file='raydist.lst',status='unknown')
c find end of output file
          do while (.true.)
             read (11,'(a)',end=775) line
          end do
775       continue
c reposition before EOF marker
          backspace(11)
       endif

       nofexp = 0
       nofbur = 0
       ntotpt = 0
       natmpt = 0	!number atom surface points
       ifirst = 0
       ipoint = 0
       ifounds = 0
       iend = 0 	!flag for end of residue
       count = 0
       iretro = 0
       bs_tot = 0.0
       bs_side = 0.0
       bs_back = 0.0
       davg_side = 0.0
       davg_back = 0.0
       nofbur_side = 0
       nofbur_back = 0
       ep_tot = 0.0
       ep_side = 0.0
       ep_back = 0.0
       dcutsq = 6.6 ** 2
       numray = 0	!number of rays each occlusion
       numray_side = 0	!number of rays each occlusion
       numray_back = 0	!number of rays each occlusion
       dstavg = 0.0
       dstavg_side = 0.0
       dstavg_back = 0.0
       dmin = 0.0

c ******************************************************************
c loop to read the surface file
c ******************************************************************

1001   continue		!main return point
       if(iend.eq.1) go to 1002     !finished-write residue summary
				    !this gets set just below 1005

       read(3,110,end=1005) line 	!part_i.ms
       go to 1006	!skip next command until end of part_i.ms

1005   continue		!here at end of part_i.ms
       iend = 1
       go to 1007

1006   continue
       if(line(41:41).eq.'A') then  !atom record in ms file
				    !else go to "surface pointsxxx"
         if(ifirst.eq.1) then	!have previous atom's dots in scratch.ms
1007       continue	!all this atoms surf records in scracth.ms
           rewind (10)		!scratch.ms
c ***************************************************************
c loop over all the atoms
c****************************************************************
c          if(natmpt.eq.0) go to 1903
	   iflati = 0
	   do 9008 iw = 1,natmpt
	     ipick(iw) = 0
9008	   continue
           do 1500 ia = 1, natom
             xdiff =  cr(ia,1) - xatc
             ydiff =  cr(ia,2) - yatc  !atc = curr atom coords
             zdiff =  cr(ia,3) - zatc
             distsq = xdiff ** 2 + ydiff ** 2 + zdiff ** 2

c if this atom greater than 6.5(?) A^2 distant exit
c            print *, distsq, dcutsq
             if(distsq.gt.dcutsq) go to 1500
             ipoint = 0

c ***************************************************************
c      loop to check the intersection of normals with atoms
c ***************************************************************

2000         continue
             ipoint = ipoint + 1	!point counter
             if(ipoint.gt.natmpt) go to 888
             read(10,110)sline
             read(sline,554) atiden,xpc,ypc,zpc,sarea,xnor,ynor,znor
             spc(1)=cr(ia,1)
             spc(2)=cr(ia,2)	!coords of possible occluding atoms
             spc(3)=cr(ia,3)
             ptc(1)=xpc
             ptc(2)=ypc		!coords of point
             ptc(3)=zpc
             sphrad = vdwrad(ia) !radius of occluding atom
             dcl = xnor
             dcm = ynor		!cos angles of normal
             dcn = znor
c check with neighboring atoms for intersection
             call intsect(spc,dcl,dcm,dcn,ptc,sphrad,
     1             iretro,iflag,dmin,dist6)
c dist6 passed to check on things

             if(iflag.eq.-1) go to 2500   !no intersection

c	scaling the dmin by the sum of van der waals radii

             dmin = dmin/2.8

c set flag for vdw overlap (1 equals overlap)
c            retropt(ipoint) = iretro

c find minimum ray length for each point (in case goes through)
	     if(ipick(ipoint) .eq. 0)then  !dot not included yet
	       dminpt(ipoint) = dmin
	       iadmin(ipoint) = ia 	!atm number occluding atom
	       ptarea(ipoint) = sarea
               retropt(ipoint) = iretro
	     endif
	     if(ipick(ipoint) .eq. 1) then  !dot included already
	       if(dmin.lt.dminpt(ipoint)) then !update closest contact
		 dminpt(ipoint) = dmin
		 iadmin(ipoint) = ia
		 ptarea(ipoint) = sarea
                 retropt(ipoint) = iretro
	       endif
	     endif

	     iflati = 1
             ipick(ipoint) = 1      !tag dot as included

2500         continue
             go to 2000             !check another normal
888          continue	!here if(ipoint.gt.natmpt)
                        !at end of scratch.ms
             rewind (10) !scratch.ms

1500       continue	!here if(distsq.gt.dcutsq)
           rewind(10)

c if no intersections, got here with iflati = 0
c add 1 to be consistent with exit from above
           if(iflati .eq.0) ipoint = natmpt + 1

c now sum areas for each patch to an occluded atom

c initialize
           davg = 0.0
           davg_back = 0.0
           davg_side = 0.0
           aep_tot = 0.0
           abs_tot = 0.0

           do 1901 ip = 1,ipoint-1
             read(10,110)sline

c read atom id area associated with point from ms file
c and add it to exposed surf area bin

             read(sline,1109) atname,areapt
1109         format(9x,a4,31x,f6.3)

c            write(6,1110) atname,areapt
1110         format('xx',a4,'xx',f6.3)
c but bail out if buried
             if(ipick(ip).eq.1) go to 1902 !if buried

             nofexp = nofexp + 1
             ep_tot = ep_tot + areapt	!residue exposed area
             aep_tot = aep_tot + areapt	!atom exposed area
c
c backbone atoms n,c,o,ca
c
             if(atname.eq.' N  ' .or. atname .eq.' C  ' .or.
     1       atname .eq. ' O ' .or. atname .eq. ' CA ' ) then
               ep_back = ep_back + areapt
             else
                 ep_side = ep_side + areapt
             endif
             ifounds = 1
             write(8,110)sline 	!exposed.ms
c this will be multiplied by 2.8 for display
             distvdw = 1.0

c if point is surface point write it to raydist.lst
c  this commented out to allow only drawing of occluded
c  below.
c            if (rayflag .eq. 1) then
c              if (sline(41:41) .eq. 'S')then
c       	 write(11,1123)sline(1:39),distvdw,sline(44:72)
1123               format(a39,f5.2,a29)
c              endif
c            endif
             go to 1901

1902         continue	!here if point is buried

c zero rays if vdw overlap
             if (retropt(ip) .eq. 1)then
                dminpt(ip) = 0.0
             endif

c If occluding atom number is same as this atom - stop
	     ipat = iadmin(ip)
	     do 878 ir = 1, nseg
	       if(ipat .ge. iats(ir).and.ipat.le.iate(ir)) go to 898
878	     continue
	     stop
898	     continue

c accumulate raylength
	     davg = davg + dminpt(ip)

c add surface area to buried bin

             bs_tot = bs_tot + areapt
             abs_tot  = abs_tot + areapt

c backbone atoms n,c,o,ca

             if(atname.eq.' N  ' .or. atname .eq.' C  ' .or.
     1       atname .eq. ' O  ' .or. atname .eq. ' CA ' ) then
	         bs_back = bs_back + areapt
                 davg_back = davg_back + dminpt(ip)
                 nofbur_back = nofbur_back + 1
             else
	         bs_side = bs_side + areapt
                 davg_side = davg_side + dminpt(ip)
                 nofbur_side = nofbur_side + 1
             endif
             write(7,110)sline 	!buried.ms
             distvdw = dminpt(ip)

c write ray to raydist.lst for display
c   this will be only occluded rays

             if (rayflag .eq. 1) then
               if (sline(41:41) .eq. 'S')then
  	           write(11,1122)sline(1:39),distvdw,sline(44:72),
     1             secnm(ir)(1:3),secsq(ir)
1122               format(a39,f5.2,a29,1x,a3,i4)
               endif
             endif
             nofbur = nofbur + 1

c return and read another sline until end reach ipoint value
1901       continue
c
c calc average of raylength for this occluded atom patch
           if(nofbur .eq.0)then
            davg=davg
           else
              davg = davg/nofbur
           endif

           if(nofbur_back .eq.0)then
            davg_back=davg
           else
              davg_back = davg_back/nofbur_back
           endif

           if(nofbur_side .eq.0)then
            davg_side=davg
           else
              davg_side = davg_side/nofbur_side
           endif

	   iuni = 0
c
           do 543 i = 1,ipoint-1	!do for each point
	     if(ipick(i).eq.0) go to 543  !no intersec this point
	     if(iuni.eq.0) then	!first time for this occluding atom
	       iuni = iuni + 1
	       nuniq(iuni) = iadmin(i)	!atom number occluding atom
	       ncount(iuni) = 1		!first occluded ray
	       parea(iuni) = ptarea(i)	!point associated area

	       dtmin(iuni) = dminpt(i)
	       go to 543
	     endif
	     do 443 j = 1,iuni	!subsequent times this occlud atm
c if this is still same occluding atom
	       if(iadmin(i).eq.nuniq(j))then
	         ncount(j) = ncount(j) + 1	!increment rays
	         parea(j) = parea(j) + ptarea(i)

c summing the dmin for all points hitting the same atom

	         dtmin(j) = dtmin(j) + dminpt(i)

	         go to 543
	       endif
443	     continue
	     iuni = iuni + 1
	     nuniq(iuni) = iadmin(i)
	     ncount(iuni) = 1
	     parea(iuni) = ptarea(i)
	     dtmin(iuni) = dminpt(i)
543	   continue

           if(abs_tot.gt.0.0) then
	     dstavg = dstavg + davg
	     numray = numray + 1

             if(atname.eq.' N  ' .or. atname .eq.' C  ' .or.
     1       atname .eq. ' O  ' .or. atname .eq. ' CA ' ) then
	        dstavg_back = dstavg_back + davg_back
                numray_back = numray_back + 1
             else
	        dstavg_side = dstavg_side + davg_side
                numray_side = numray_side + 1
             endif
             call asorder(nuniq,ncount,parea,dtmin,iuni)
c
	     satot = aep_tot + abs_tot
       	     write(2,3344) sline(10:13),aep_tot,abs_tot,satot,davg,
     1                  resinf(3)
3344         format('AVG for ATOM:',a4,f9.3,' es',f9.3,' os',
     1             f9.3,' ts',f9.3,' Rln',1x,a8)

             do 174 ip = 1, iuni

c if occluding atom number is same at this atom - stop
               ipat = nuniq(ip) !atom number occluding atom

	       do 276 ir = 1, nseg
		 if(ipat .ge. iats(ir).and.ipat.le.iate(ir)) go to 298
276	       continue
	       stop
298	       continue


	       fratom(1:4) = sline(10:13)
	       toatom(1:4) = atnm(ipat)
c        print*, xatc,yatc,zatc,sline(10:13)
c              print*, (cr(ipat,l),l=1,3),atnm(ipat)

             xdiff =  cr(ipat,1) - xatc
             ydiff =  cr(ipat,2) - yatc  !atc = curr atom coords
             zdiff =  cr(ipat,3) - zatc
             dxx = sqrt(xdiff ** 2 + ydiff ** 2 + zdiff ** 2)

	       icf = 0
	       do 143 ih = 1,4
	         if(fratom(ih:ih) .eq. ' ')go to 143
	         do 145 ic = ih,	4
	           icf = icf + 1
	           fratom(icf:icf) = fratom(ic:ic)
	           fratom(ic:ic) = ' '
145	         continue
	         go to 149
143	       continue
149 	       continue
	       do 146 ih = 1,4
	         if(fratom(ih:ih).eq.' ') fratom(ih:ih) = '_'
	         if(toatom(ih:ih).eq.' ') toatom(ih:ih) = '_'
146	       continue
	       write(2,7734) resinf(3),fratom,secnm(ir)(1:3)
     1             ,secsq(ir),toatom,ncount(ip),parea(ip)
     2             ,dtmin(ip)/ncount(ip),dxx
7734	       format('INF',1x,a8,'@',a4,'>',a3,i4,'@',a4
     1           ,i6,' pts',f8.3,' A2 ',f8.3,' Rlen',f6.2,' Dxx')
174	     continue
           endif	!.bs_tot.gt.0.0
1903       continue	!here if(natmpt.eq.0)

c here if ifirst = 0 (just read atom records)
         endif

         natmpt = 0
         ifirst = 1
         if(iend.eq.1) go to 1002 !when eof part_i.ms, write atm summary
         write(7,110)line	!buried.ms
         write(8,110)line	!exposed.ms

c reading the x,y and z coordinates from the 'A' record

         read(line,554) atiden,xatc,yatc,zatc

         ip = 0

c check  which column first letter of atom id is in

         do 231 id = 10,13
	   if(line(id:id) .eq. ' ') go to 231
	   ip = ip + 1
	   dumyat(ip:ip) = line(id:id)
c dumyat now is first letter of atom id
231    	 continue

c assigning vdw radius

         sprvdw = 0.0
         do 1258 ira = 1, nvdwt
c first the 3 letter atoms
           if(attype(ira)(3:3).ne.' ') then
             if(dumyat(1:3).eq.attype(ira)(1:3)
     &.and.(restyp(ira).eq.'   ')) then
                sprvdw = vdwr(ira)
             elseif((dumyat(1:3).eq.attype(ira)(1:3))
     &.and.(atiden(1:3).eq.restyp(ira)))then
                sprvdw = vdwr(ira)
                go to 1259
             endif

c second the 2 letter atoms
           elseif(attype(ira)(2:2).ne.' ') then
             if(dumyat(1:2).eq.attype(ira)(1:2)
     &.and.(restyp(ira).eq.'   ')) then
                sprvdw = vdwr(ira)
             elseif((dumyat(1:2).eq.attype(ira)(1:3))
     &.and.(atiden(1:3).eq.restyp(ira)))then
                sprvdw = vdwr(ira)
                go to 1259
             endif

c last the 1 letter atoms
           else if(attype(ira)(2:2).eq.' ')then
             if(dumyat(1:1).eq.attype(ira)(1:1)) then
               sprvdw = vdwr(ira)
             endif
           endif

1258     continue
1259     continue
         if(sprvdw.le.0.001)then
           print *, ' sphere vdwrad is zero '
           print *, line
           stop
         endif

         nofbur = 0
         nofbur_back = 0
         nofbur_side = 0
         rewind(10)	!scratch.ms
         go to 1001	!go back up (out of loop) and read next line
                        !of part_i.ms

       endif 		!from line(41:41) .eq. 'A'

c *******************************************************************
c surface pointsxxx
c
c format(a13,3(1x,f8.3),4x,4f7.3)
cALA    6  N     4.210    7.064   -2.933 SC0  0.219 -0.155  0.268  0.951
c  atiden         xpc      ypc      zpc       sarea   xnor   ynor   znor
c *******************************************************************
       ntotpt = ntotpt + 1
       read(line,554) atiden,xpc,ypc,zpc,sarea,xnor,ynor,znor

c      print *, atiden
c write surface records for this atom
         if(atiden(1:8).eq.resinf(3)(1:8))then
           natmpt = natmpt + 1
c write this line of part_i(v).ms to scratch.ms

           write(10,111) line	!scratch.ms
111        format(a80)
         endif
         go to 1001	!go back to main loop beginning

c *******************************************************************
c      end of loop to read the surface file
c *******************************************************************

1002   continue
       write(2,1213) ep_tot,resinf(3)
       write(2,1214) ep_back,resinf(3)
       write(2,1215) ep_side,resinf(3)
       write(2,4213) bs_tot,resinf(3)
       write(2,4214) bs_back,resinf(3)
       write(2,4215) bs_side,resinf(3)

c so that we won't divide by zero, but we want zero in output
       if(numray.eq.0)then
         dstavg = 0.0001
         numray = 1
       endif
       if(numray_back.eq.0)then
         dstavg_back = 0.0001
         numray_back = 1
       endif
       if(numray_side.eq.0)then
         dstavg_side = 0.0001
         numray_side = 1
       endif

       write(2,5215) dstavg/numray,resinf(3)
       write(2,5217) dstavg_back/numray_back,resinf(3)
       write(2,5219) dstavg_side/numray_side,resinf(3)

554    format(a13,3(1x,f8.3),4x,4f7.3)
1213   format('    ES_Total      :',f9.3,39x,a8)
1214   format('    ES_Backbone   :',f9.3,39x,a8)
1215   format('    ES_Side Chain :',f9.3,39x,a8)
4213   format('    OS_Total      :',f9.3,39x,a8)
4214   format('    OS_Backbone   :',f9.3,39x,a8)
4215   format('    OS_Side Chain :',f9.3,39x,a8)
5215   format('    Ave_Raylength :',f9.3,39x,a8)
5217   format('    MC_RAYLENGTH  :',f9.3,39x,a8)
5219   format('    SC_RAYLENGTH  :',f9.3,39x,a8)

       close(7)
       close(8)
       close(10)
       if (rayflag .eq. 1) then
         close(11)
       endif
       close(3)
       return
       end

c****************************************************************

        subroutine intsect(spc,dcl,dcm,dcn,ptc,sphrad,
     1          iretro,iflag,dmin,dist6)

	dimension spc(3),ptc(3),ptofin(2,3)
	xsp = spc(1)
        ysp = spc(2)	!coords of possible occluding atom
	zsp = spc(3)
	xpt = ptc(1)
	ypt = ptc(2) 	!coords of point
	zpt = ptc(3)

	a = 1.0
	b = dcl * (xpt - xsp)  !(cos_angle_x * (point->atom xdist))
        b = b + dcm * (ypt - ysp)
        b = b + dcn * (zpt - zsp)
	b = 2.0 * b
	c = (xpt**2 + xsp**2 + ypt**2 + ysp**2 + zpt**2 + zsp ** 2)
	c = c - 2 * ( xpt * xsp + ypt * ysp + zpt * zsp)
	c = c - sphrad**2
	bs4ac = b**2 - 4*a*c
	if(bs4ac.lt.0.0) then
	  iflag = -1	!no intersection
	  go to 999
	endif
	anuera = -b + sqrt(bs4ac)
	denom  = 2 *a
	root1 = anuera/denom
	anuera = -b - sqrt(bs4ac)
	root2 = anuera/denom
	ptofin(1,1) = xpt + root1*dcl
	ptofin(1,2) = ypt + root1*dcm
	ptofin(1,3) = zpt + root1*dcn
	ptofin(2,1) = xpt + root2*dcl
	ptofin(2,2) = ypt + root2*dcm
	ptofin(2,3) = zpt + root2*dcn
c
c	to check whether the vdw sphere is along the forward
c	or backward direction of the normal.
c	this was done by calculating the distances between
c	the surface point and the intersecting points (d1,d3) and between
c	the point at a unit distance from the surface point along
c	the normal and the intersecting points (d2, d4).
c	if d1 > d2 or d3 > d4 then the sphere is along the direction
c	of the normal
c
c
	xpton = xpt + dcl
	ypton = ypt + dcm
	zpton = zpt + dcn
c
c	intersecting point 1
c
	dist1 = (xpt - ptofin(1,1)) ** 2
	dist1 = dist1 + (ypt - ptofin(1,2)) ** 2
	dist1 = dist1 + (zpt - ptofin(1,3)) ** 2

	dist2 = (xpton - ptofin(1,1)) ** 2
	dist2 = dist2 + (ypton - ptofin(1,2)) ** 2
	dist2 = dist2 + (zpton - ptofin(1,3)) ** 2
	jflag = -1
	if(dist1.gt.dist2) then
	  jflag = 1
	endif
c
c	intersecting point 2
c
	dist3 = (xpt - ptofin(2,1)) ** 2
	dist3 = dist3 + (ypt - ptofin(2,2)) ** 2
	dist3 = dist3 + (zpt - ptofin(2,3)) ** 2

	dist4 = (xpton - ptofin(2,1)) ** 2
	dist4 = dist4 + (ypton - ptofin(2,2)) ** 2
	dist4 = dist4 + (zpton - ptofin(2,3)) ** 2
	kflag = -1
	if(dist3.gt.dist4) then
	  kflag = 1
	endif

	iflag = -1
	if(jflag .eq. 1 .or. kflag .eq. 1) then
	  iflag = 1
	endif

	dmin = dist1
	if(dmin.gt.dist3) dmin = dist3
	dmin = sqrt(dmin)

c now check for vdw overlap
        iretro = 0
        dist5 = (xpt - xsp) ** 2
        dist5 = dist5 + ((ypt - ysp) ** 2)
        dist5 = dist5 + ((zpt - zsp) ** 2)
        dist6 = sqrt(dist5)
c this is distance from point to occluding atom center
        if (dist6 .lt. sphrad) then
           iretro = 1
        endif
c now check for distance greater than water molecule

        if (dmin .gt. 2.8) then
            iflag = -1
        endif

999	continue

	return
	end

c********************************************************

	subroutine assvdw(attype,vdwr,nvdwt,restyp)

        integer i,ntype
        real vdwr
	character*3 attype,restyp
	dimension attype(50),vdwr(50),restyp(50)

	open(unit=12,file='radii',status='old')

	i = 0
	ntype = 50
100	continue

c if it wasn't found in list say so
        i = i + 1
cif(i.gt.ntype) then
c        print *, ' attype for vdw radii is greater than ',ntype
c        stop
c       endif

c read radii from list

	read(12,15,end=101) attype(i),restyp(i),vdwr(i)
15	format(a3,3x,a3,2x,f4.2)
	go to 100	!go back to next read

101	continue
	nvdwt = i-1	!number of different atom radii read
	close(12)
	return
	end

c ******************************************************************

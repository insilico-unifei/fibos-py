c*************************************************************
c	Program RESPAK
c
c       Version 75
c
c	Calculates a packing parameter for each residue
c	in a protein based on the packing value equal to
c	(sum[occ surf])*(1-<raylength>)	
c	This is normalized by dividing by the following
c	(totsurf(restype))
c	where totsurf(restype) is the surface area for the
c	residue, i, and is sum of ES_Total and OS_Total 
c	Requires a "prot.srf" file from the OS package.
c       This version uses the distributed OS package.
c
c       Copyright (c) Yale University, New Haven, CT 06520
c       All rights reserved.
c       The program either in full or in part should not be distributed
c       in any form or by any means without the prior written 
c       permission of the author:
c
c       Patrick Fleming
c       Yale University
c       Dept. of Molecular Biophysics and Biochemistry
c       P. O. Box 208114
c       260 Whitney Ave.
c       New Haven, CT 06520-8114
c       fleming@csb.yale.edu

c************************************************************

	subroutine respak

	integer i,j,k,lunit
	integer maxres

	parameter (maxres=10000)

	character fname*40,prompt(2)*22,line*80,resnam(maxres)*3
        character aa(24)*3

	character resnum(maxres)*4
	integer cntcres(maxres)

 	integer current


c The "packing value" = sum[occ surf*(1-raylength)]

        real respv(maxres)
        real pupv(maxres)
        real resos(maxres)
        real totos(maxres)
	real os,raylen
        real ES_Total(maxres)
        real OS_Total(maxres)
	real srfaa(maxres)
        real norpak(maxres)

	integer dummy
	integer frstres

c These used to identify restype for assigning total surface area
        data aa/'ALA','ARG','ASN','ASP','CYS','CYT',
     1  'GLN','GLU','GLY','HIS','ILE',
     2  'LEU','LYS','MET','PHE','PRO',
     3  'SER','THR','TRP','TYR','VAL',
     4  'UNK','HEX','PRG'/

c Get filenames and open files

c       data prompt /
c    &  'Name of the .srf file?','Name for output file? ' /

        lunit=1
c       call askfil(lunit,fname,'old',prompt(lunit))
c For the OS package hardwire the file names
        open(unit=lunit,file='prot.srf',status='unknown')

        lunit=2
c       call askfil(lunit,fname,'unk',prompt(lunit))
c For the OS package hardwire the file names
        open(unit=lunit,file='prot.pak',status='unknown')

c Initialize arrays

        do i=1,maxres
          respv(i)=0.0
        end do
        do i=1,maxres
          pupv(i)=0.0
        end do
        do i=1,maxres
          resos(i)=0.0
        end do
        do i=1,maxres
          totos(i)=0.0
        end do
	do i=1,maxres
	  cntcres(i)=0
	end do

c Formats
100	  format(a80)
102	format(71x,i4)
104     format(20x,f9.3)
108     format(4x,'Resnum',2x,'Resname',7x,'OS',5x,'os*[1-raylen]'
     1   ,3x,'OSP')
110     format(5x,i4,5x,a3,5x,f7.2,5x,f7.2,5x,f7.3)
112     format(5x,i4,'      ?         0.0         0.0        0.0')

c Start main program

	i=0
	j=1
	current = 1

c Get number of first residue for later register
	read(1,102)frstres

	i=frstres

	do while (.true.)
101	  continue
	  read(1,100,end=909)line
	  if(line(1:3) .eq. 'INF')then
	    resnam(i) = line(5:7)
	    resnum(i) = line(9:12)
	    backspace (unit = 1)
 	    call rescalc(cntcres,i,j,respv,current,resos)
            totos(i)=totos(i)+resos(i)
 	    pupv(i)=pupv(i)+respv(i)
 	  else if(line(1:3) .eq. 'AVG') then
 	    goto 101

          else if(line(5:12) .eq. 'ES_Total') then
            backspace (unit = 1)
            read(1,104)ES_Total(i)

          else if(line(5:12) .eq. 'OS_Total') then
            backspace (unit = 1)
            read(1,104)OS_Total(i)
            srfaa(i) = ES_Total(i) + OS_Total(i)

 	  else	if(line(5:10) .eq. 'SC_RAY') then
c
	    i=i+1

 	  end if

 	end do

909	continue

c Now calculate normalized packing
c norpak=[pupv/(total surf(i)]

        do k=frstres,(i-1)
          norpak(k)=pupv(k)/(srfaa(k)) 
        end do

c Write out packing to file, "prot.pak"
        write(2,108)
 	do k=frstres,(i-1)
            dummy = 0
          do j=1,24
            if(resnam(k).eq.aa(j))then
              write(2,110)k,resnam(k),totos(k),pupv(k),norpak(k)
              dummy = 1
            end if
          end do
            if (dummy .eq. 0) then
c
c
              write(2,112)k
            end if
 	end do
 	CLOSE(1)
 	CLOSE(2)
	end subroutine respak
c	end subroutine respak

c*************************************************************
c
        subroutine askfil (lunit, fname, age, prompt)
c
c*************************************************************

c Argument declarations
	integer lunit
	character*3 age
	character*40 fname
	character*(*) prompt

	write(6,100) prompt
	read(5,200) fname

	if (age .eq. 'old') then
	  open(unit=lunit,file=fname,status='old')
	else if (age .eq. 'new') then			
	  open(unit=lunit,file=fname,status='new')
	else if (age .eq. 'unk') then
	  open(unit=lunit,file=fname,status='unknown')
	end if

100	format(5x,a,'---->',$)
200	format(a40)

	return
	end

c**************************************************************
c
	subroutine search(array, target, found)
c
c Searches contact residue array to see if current occluding
c residue (number) has already been tagged.
c Returns logical true or false
c**************************************************************

c Argument Declarations
	integer array(*)
	integer target
	logical found

c Local Declarations
	integer i
	
c Compare each element to target
	i=1
	do while(.not. found)
	  if (array(i) .eq. target) then
	    found = .true.
	  else
	    i=i+1
	  end if
	end do

	return
	end

c************************************************************
c
 	subroutine rescalc(cntcres,i,j,respv,current,resos)
c
c Calculates residue packing value (respv) for a residue
c************************************************************

        parameter (maxres=10000)
c Argument Declarations
        integer cntcres(1000)
        real respv(maxres)
        real resos(maxres)
 	integer current,i

c Local Declarations
        real os,raylen
 	logical found

        read(1,110)cntcres(i),os,raylen
110     format(21x,i4,18x,f5.3,7x,f5.3)
c       print*, os, raylen

c If raylen is greater than 1.0 make it 1.0
	if (raylen .gt. 1.0) then
	  raylen = 1.0
	end if

c Find out if contact residue has been tagged
c (Done only if working on current residue, not occluding res

c	if (current .eq. 1) then
c         found= .false.  
c         call search(cntcres,cntcres(i),found)
c         if(found) then
c           j=j+1
c         end if
c	end if

c Calculate residue packing value
        resos(i)=os
        respv(i)=(os*(1-raylen))
c       print*, respv(i)
c       print*, ' '


	return
	end

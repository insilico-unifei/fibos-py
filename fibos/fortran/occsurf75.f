 
        subroutine occsurf

c Calculates the weighted surface area, Sw(Qsr),
c  the weighted protein occlude surfac, Pw,
c  the ideal weighted protein occluded surface, Pi,
c  and the normalized protein surface ratio, Psr.
c
c The input file is called "prot.srf" and is the output of 
c "surfcal.f" 
c
c Needs "tot_aa_hsg.lis" in the CWD.
 
	common/aaa/ sfeaa(20,10),aaname(20),a1c(20),sfe(20)
	character*1 a1c,onelc
	character*3 aaname,strng
	character*9 qdate
	character*8 qtime
	character*40 filenam,plotfl,outfl
	character*80 line
	character*13 residen


10	format(a)
20	format(20x,f8.3)
40     format(4X,'Total_Res:    ',I6)
42     format(4X,'Pw:        ',F9.2)
44     format(4X,'Pi:        ',F9.2)
46     format(4X,'Psr:       ',F9.2)
 
c Get the distributions of aa occluded surface ratio
c  and put them in array sfeaa
	call aadist
 	
100	continue
	pren_tot = 0.0
	pridl_tot = 0.0
 
c Read the prot.srf file to get the first residue name
c  and number

 	open(unit=1,file='prot.srf',status='old')
	read(1,10) line
	residen = line(68:75)
c This would be something like: "ARG    1"

c Initialize
	psact = 0.0
	pside = 0.0
	ptact = 0.0
	ptide = 0.0
	itotres = 0
	itotgly = 0
300	continue

c Read each line looking for the residue summary
	read(1,10,end=888)line

c Update when new residue is seen
	if(line(68:75).ne.residen)then
	   residen = line(68:75)
	endif	

c Here is residue summary section
	if(line(5:11).eq.'ES_Side') then
	   read(line,20)esside
	else if(line(5:11).eq.'OS_Side') then
	   read(line,20)bsside
	else if(line(5:12).eq.'OS_Total') then
	   read(line,20)bstotal
	else if(line(5:12).eq.'ES_Total') then
	   read(line,20)estotal
	else if(line(5:10).eq.'SC_RAY') then
 
c Increment the total # residues counter
	   itotres = itotres + 1
 
c Translate aa name to one letter code
           strng = residen(1:3)
           call aa1to3(strng,onelc)

	   ts = bstotal + estotal	!Calc total surface area
	   if(ts.le.0.0) then
	     print *, ' Error: residue has zero amino acid surface ', 
     1       residen
	     go to 300
	   endif
c Now get the weighted surface area for the resdue, i and the
c  and the ideal weighted surface area for that residue type
  
           call getsfe(onelc,sfdeaa,sfdeaai,bstotal,estotal)

c Sum for calculation of protein as a whole
	   ptact = ptact + sfdeaa	!Actual value for i
	   ptide = ptide + sfdeaai 	!Ideal value for type

c Write the prot.eval file
           write(6,'(a,5x,F9.2,5x,a)')'Sw:  ',sfdeaa,residen(1:8)
           write(6,'(a,5x,F9.2,5x,a)')'St:  ',sfdeaai,residen(1:8)
           write(6,'(a,7x,F9.2,5x,a)')'Ri:',
     +               sfdeaa/sfdeaai,residen(1:8)
	endif

c Go back and read another residue
	go to 300

c Arrive here when eof reached on prot.srf file
888	continue
	
	CLOSE(1)
        inquire(1)

c Write the protein summary in prot.eval
        write(6,40) itotres
        write(6,42) ptact
        write(6,44) ptide
        write(6,46) ptact/ptide
	end
c****************************************************************
	subroutine aa1to3(strng,onelc)

c Translates 3 letter code to 1 letter code

	character*3 strng
	character*3 aa3c(21)
	character*1 aa1c(21)
C
	character*1 onelc
	natype = 21
	data aa3c/'ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS',
     1            'ILE','LEU','LYS','MET','PRO','SER','THR','VAL','TYR',
     2            'PHE','TRP','CYT'/
	data aa1c/'A','R','N','D','C','E','Q','G','H','I','L','K','M','P',
     1            'S','T','V','Y','F','W','C'/
	do 10 i = 1,natype
	if(strng(1:3).eq.aa3c(i)(1:3))then
	   ival = i
	   go to 11
	endif
10	continue
	    print *, ' Error in Amino acid ident ', strng(1:3)
	    stop
11	continue
	onelc = aa1c(ival)
	return
	end
c**************************************************************
	subroutine aadist

c Reads the data file containing the distributions of occluded
c  surface ratios for the 20 amino acids.
c Calculates 
	common/aaa/ sfeaa(20,10),aaname(20),a1c(20),sfe(20)
	character *3 aaname	
	character *1 onelc,a1c
	dimension nadist(20,10)

10	format(a3,f8.3)
30	format(I4)

c Open the distribution file in the CWD
 	open(unit=8,file='tot_aa_hsg.lis',status='old')
	
	i = 0
100	continue
	i = i + 1
	read(8,10,end=999) aaname(i),sfe(i)
c This is the first line of each aa entry 
c  It should be something like: "ALA  68.75"
c  sfe is the surface area of the residue

c Translate the 3 letter code to 1 letter code
	call aa1to3(aaname(i),onelc)
	a1c(i) = onelc
	max = -9999

c Now read the distribution of surface ratios.
c  These are in the form of a histogram with
c  ten bars, each bar give the incidence of that
c  occluded surface ratio, i.e, 0.0-0.09,0.1-0.19,
c  etc. 
	do 20 k = 1, 10
	  read(8,30) nadist(i,k)
	  if(nadist(i,k).gt.max)then
	    max = nadist(i,k)
	  endif
20	continue
 
c Now normalize this distribution to the total 
c  surface area of the residue and put it in 
c  the array sfeaa, called the weighted surf area
	do 40 kk = 1, 10
	   anum  = nadist(i,kk)/(1.0*max)
	   sfeaa(i,kk) = anum * sfe(i)
40	continue

c Go back and read the distribution for the next aa
	go to 100

999	continue	!Here after reading distributions
	return
	end
c**************************************************************
	subroutine getsfe(onelc,sfdeaa,sfdeaai,bstotal,estotal)

c This 
	common/aaa/ sfeaa(20,10),aaname(20),a1c(20),sfe(20)

	character*1 onelc
	character*1 a1c
	character*3 aaname
	
	do 10 i = 1, 20
	   if(onelc.eq.a1c(i)) then
	      iaa = i
	      go to 20
	   endif
10	continue
	print *, ' Error Does not match in getsfe routine '
	stop
20	continue
	sfdeaai = sfe(iaa)	!Total surface area of res
	sr = bstotal/(bstotal+estotal)
	iptr = 0
	if(sr .ge. 0.0 .and. sr .le. 0.1) iptr = 1
	if(sr .gt. 0.1 .and. sr .le. 0.2) iptr = 2
	if(sr .gt. 0.2 .and. sr .le. 0.3) iptr = 3
	if(sr .gt. 0.3 .and. sr .le. 0.4) iptr = 4
	if(sr .gt. 0.4 .and. sr .le. 0.5) iptr = 5
	if(sr .gt. 0.5 .and. sr .le. 0.6) iptr = 6
	if(sr .gt. 0.6 .and. sr .le. 0.7) iptr = 7
	if(sr .gt. 0.7 .and. sr .le. 0.8) iptr = 8
	if(sr .gt. 0.8 .and. sr .le. 0.9) iptr = 9
	if(sr .gt. 0.9 .and. sr .le. 1.0) iptr = 10
	if(iptr.eq.0)then
	   print *,' Error in iptr = ',iptr,sr,bstotal,estotal,iaa,
     &           a1c(iaa)
	   stop
	endif

c Get the weighted surface area for that relative occluded 
c  surface ratio
	sfdeaa = sfeaa(iaa,iptr)

	return
	end

C2345678 1 2345678 2 2345678 3 2345678 4 2345678 5 2345678 6 2345678 7 2
        subroutine renum

c Written by P.J. Fleming, March, 1994

c This progran renumbers the residues in a clean PDB file.
c It uses the existing number of the first residue as
c the first number and increments from there.

c Declarations
       character line*80,kywrd*4
       character pada*13,atm*4,res1*3,res2*3,padb*54
       character fname*30, prompt*30
       integer num,onenum

       fname = "temp.cln"

c       prompt = "temp1.cln"

c       data prompt / 'PDB filename to be renumbered?' /

c Open the current input and output files.
       call ask(1,fname,'old',prompt)
       open (unit=2,file='new.pdb',status='unknown')       
c       write(6,599)
c599    format('     Renumbered file is called new.pdb')

c Read the current first residue number
       num=0
5      read(1,'(a80)',END=90)line 
       if (line(1:4) .eq. 'ATOM') then
          backspace (unit=1)
          read(1,'(a12,1x,a4,a3,2x,i4,a54)')pada,atm,res1,onenum,padb
c         print*,res1
          num=onenum     !set num = to actual first res number
c Write first record to output
          write(2,'(a12,1x,a4,a3,2x,i4,a54)')pada,atm,res1,num,padb
       else    
          go to 5
       end if
       do while (.true.)
          read(1,'(a12,1x,a4,a3,2x,i4,a54)',end=95)pada,atm,res2,
     +          onenum,padb
c         print*,res1,res2
          if ((res1 .eq. res2) .and. (atm .ne. 'N   ') .and.
     1       (atm .ne. 'P   ') .and.
     2       (res2 .ne. 'WAT') .and. (res2 .ne. 'HOH') .and.
     3       (res2 .ne. 'ZN ') .and. (res2 .ne. 'PO4') .and.
     4       (res2 .ne. 'SO4') .and. (res2 .ne. 'HEM')) then 
             write(2,'(a12,1x,a4,a3,2x,i4,a54)')pada,atm,res2,num,padb
          else   !Now on to next residue
             if (pada(1:3) .eq. 'END') go to 95
	     if (pada(1:3) .eq. 'TER') go to 85
             num=num+1
             write(2,'(a12,1x,a4,a3,2x,i4,a54)')pada,atm,res2,num,padb
             res1=res2
85	  continue
          end if
       end do
90     print*,'ERROR: No ATOM record'
95     continue
       write(2,'(a3)')'END'
       close (unit=1)
       close (unit=2)
       end
*****************************************************************

        SUBROUTINE ask ( lunit, fname, age, prompt )

****************************************************************

        INTEGER lunit		!Unit number
        CHARACTER*3    age              ! either "OLD" or "NEW"
        CHARACTER*(*)  prompt, fname
C
C               THE SIZE OF THE MESSAGE BUFFER IS DEFINED AT RUN TIME
C               SINCE IT IS USED AS A PASSED LENGTH STRING. FNAME IS 
C               A LOCALLY DEFINED CONSTANT CURRENTLY SET AT 30.
C                                                             
C                                                             
c        WRITE(6,600) prompt
c        READ(5,500) FNAME
        FNAME = "temp.cln"
c        write(6,610)

        IF ( ( age .EQ. 'old' ) .OR. ( age .EQ. 'OLD' ) ) THEN
            OPEN ( UNIT = lunit, FILE = fname, STATUS = 'OLD')
        ELSE IF ( ( age .EQ. 'new' ) .OR. ( age .EQ. 'NEW' ) ) THEN
            OPEN ( UNIT = lunit, FILE = fname, STATUS = 'NEW')
        END IF
c        WRITE(6,630) fname, age

500     FORMAT ( A30 )
600     FORMAT ( 2x,A, '---->',$ )
610     FORMAT ( 1X )
630     FORMAT(5X,A,'Status=', A3)

        RETURN
        END


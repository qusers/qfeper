!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        By Masoud Kazemi 2013/08/09         !
!           revised 2014/06/02               !
!  top_parser module for writting fep files  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module top_parser

implicit none


!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 type decleration                                            !
!-----------------------------------------------------------------------------------------------------------------------------!

type  ::  atomtype_type
        integer(4)                :: code=0,atmno=0                                     !atom codes (a number that corresponds to vdwpar array index for mass, aii,bii, type name)
        character(len=10)         :: nam="0",indic="!"                                  !type name
        real(4)                   :: charge=0,mass=0,aii=0,bii=0,ci=0,ai=0              !type for different parameters of atoms
end type atomtype_type

type  ::  bondtype_type
        integer(4)                :: code=0,a1=0,a2=0,indic=0,line_counter=0            !bond codes,atom 1 and 2 of bond,indicator (used for mapping between different states)
        real(4)                   :: force=0,dis=0                                      !bond force and distance
end type bondtype_type

type  :: angltype_type
        integer(4)                :: code=0,a1=0,a2=0,a3=0,indic=0,line_counter=0       !angle codes,atom 1, 2 and 3 of angle,indicator (used for mapping between different states)
        real(4)                   :: force=0,dis=0                                      !angle force and angle (its called dis one can think as angular distance)
end type angltype_type

type  :: torstype_type
        integer(4)                :: nterm=0,line_counter(3)=0                          !for grouping tosions(used in q_lib)
        integer(4)                :: a1=0,a2=0,a3=0,a4=0                                !torsion codes,atom 1, 2,3 and 4 of torsion,indicator (used for mapping between different states)
        integer(4)                :: code(3)=0,indic(3)=0
        real(4)                   :: force(3)=0,dis(3)=0,multip(3)=0                    !type for different parameters of atoms
end type torstype_type

type  :: imprtype_type
        integer(4)                :: code=0,a1=0,a2=0,a3=0,a4=0,indic=0,line_counter=0 !bond codes,atom 1 and 2 of bond
        real(4)                   :: force=0,dis=0                                      !type for different parameters of atoms
end type imprtype_type

type ::  bondpar_type
        integer(4)                :: code
        real(4)                   :: force,dis
end type bondpar_type

type  ::  anglpar_type
        integer(4)                :: code
        real(4)                   :: force,dis
end type anglpar_type

type  ::  torspar_type
        integer(4)                :: code
        real(4)                   :: force,multip,dis
end type torspar_type

type  ::  imprpar_type
        integer(4)                :: code
        real(4)                   :: force,dis
end type imprpar_type

type  ::  vdwpar_type
        real(4)                   :: mass,aii,bii
        character(len=10)         :: nam
end type vdwpar_type

type  ::  top_type
        type(atomtype_type),allocatable,dimension(:)        :: atomtype
        type(bondtype_type),allocatable,dimension(:)        :: bondtype
        type(angltype_type),allocatable,dimension(:)        :: angltype
        type(torstype_type),allocatable,dimension(:)        :: torstype
        type(imprtype_type),allocatable,dimension(:)        :: imprtype
        type(bondpar_type) ,allocatable,dimension(:)        :: bondpar
        type(anglpar_type) ,allocatable,dimension(:)        :: anglpar
        type(torspar_type) ,allocatable,dimension(:)        :: torspar
        type(imprpar_type) ,allocatable,dimension(:)        :: imprpar
        type(vdwpar_type)  ,allocatable,dimension(:)        :: vdwpar
        integer(4)                                          :: atomtotal,bondtotal,angltotal,torstotal,imprtotal,qtomtotal
        integer(4)                                          :: bondpartotal,anglpartotal,torparmtotal,imprpartotal,vdwtotal
        logical                                             :: atomtype_flag,bondtype_flag,angltype_flag
        logical                                             :: torstype_flag,imprtype_flag
end type top_type


 logical                                                    :: top_flag
 character(len=1)                                           :: debug,split
 character(len=100)                                         :: qformat
 character(len=5)                                           :: nstate_char,option


public                ::        read_line

interface read_line
        module procedure read_line_int,read_line_real,read_line_char
end interface read_line


contains

!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 function top_reader                                         !
!-----------------------------------------------------------------------------------------------------------------------------!
!This subroutine read all needed parameter from topology file (qatom and nonqatoms)

subroutine top_reader(sfilename,top)
        character(len=*),intent(in)                 ::        sfilename
        type(top_type),intent(inout)                ::        top
!local
        character(len=132)                          ::        line
        integer(4)                                  ::        i,j,k,u,h,msg=0
        character(100)                              ::        header,formatt
!-----------------------------------------------------------------------------------------------------------------------------!
!initial control flags
        top%atomtype_flag = .false.
        top%bondtype_flag = .false.
        top%angltype_flag = .false.
        top%torstype_flag = .false.
        top%imprtype_flag = .false.

!open topology file
        u=freefile()
        open (file=sfilename, unit=u, form="formatted", access="sequential",status="old", iostat=msg, action="read")
        call errormsg(msg,"file could not be oppend",sfilename)
if(debug .ne. "p") print*,"open topology file ",sfilename

! find atoms code----------------------------------------------------------------------------------------------------
        formatt = '(a11,a25)'
        header = "No. of integer atom codes"
        call section(header,formatt,u,line,top%atomtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%atomtotal
        if ( .not. top%atomtype_flag ) call errormsg(msg," have no atom code. Nothing is there to do",sfilename)

if ( top%atomtype_flag ) then
!count how many line of atom code is there
if(debug .ne. "p") print*,"reading vdw codes"
        k=16*(CEILING(real(top%atomtotal)/16)) + 16

!read atom codes
        allocate (top%atomtype(k),stat=msg)
        call errormsg(msg,"atom type array could not be allocated")
        j=0
        do i=1, CEILING(real(top%atomtotal)/16)
                read (unit=u, fmt='(a)',iostat=msg)  line 

                if (msg < 0) call errormsg(1," seems corupted, In atom code section.",sfilename)
                read (unit=line,fmt=*,iostat=msg), (top%atomtype(j+h)%code,h=1,16)
                if (msg < 0) call read_line(top%atomtype(j+1:j+16)%code,line)
                top%atomtype(j+1:j+16)%atmno=[ j+1,j+2,j+3,j+4,j+5,j+6,j+7,j+8,j+9,j+10,j+11,j+12,j+13,j+14,j+15,j+16 ]
                j=j+16
        end do

! find atom charge
if(debug .ne. "p") print*,"reading charges"
        formatt = '(a11,a21)'
        header = "No. of atomic charges"
        call section(header,formatt,u,line,top%atomtype_flag)
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out charge section",sfilename)

!read atom charges
        j=0
        do i=1, CEILING(real(top%atomtotal)/10)
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In atom charge section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) (top%atomtype(j+h)%charge,h=1,10)
                if (msg < 0) call read_line(top%atomtype(j+1:j+10)%charge,line)
                j=j+10
        end do

! find how many vdw parameter
if(debug .ne. "p") print*,"reading total parameter"
        formatt = '(a11,a17)'
        header = "No. of atom types"
        call section(header,formatt,u,line,top%atomtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%vdwtotal
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out # of parameter section",sfilename)

        k=10*(CEILING(real(top%vdwtotal)/10)) + 10
        allocate (top%vdwpar(k),stat=msg)
        call errormsg(msg,"VdW parameter array could not be allocated")
        top%vdwpar(:)%mass =  0
        top%vdwpar(:)%aii  =  0
        top%vdwpar(:)%bii  =  0
        top%vdwpar(:)%nam  = "0"

!find atom masses
if(debug .ne. "p") print*,"reading masses"
        formatt = '(a1,a6)'
        header = "asses:"
        call section(header,formatt,u,line,top%atomtype_flag)
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out mass section",sfilename)

!read atom masses

        j=0
        do i=1, CEILING(real(top%vdwtotal)/10)
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In atom mass section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) (top%vdwpar(j+h)%mass,h=1,10)
                if (msg < 0) call read_line(top%vdwpar(j+1:j+10)%mass,line)
                j=j+10
        end do

!assign masses
if(debug .ne. "p") print*,"assigning masses"
        do i=1,top%atomtotal
                do j=1,top%vdwtotal
                        if ( abs(j - top%atomtype(i)%code) .lt. 0.01 ) then
                                top%atomtype(i)%mass = top%vdwpar(j)%mass
                        end if
                end do
        end do

!find atom vdw aii
if(debug .ne. "p") print*,"reading vdw aii"
        formatt = '(a1,a17)'
        header = "qrt (Aii) normal:"
        call section(header,formatt,u,line,top%atomtype_flag)

! search for charmm type
        if ( .not. top%atomtype_flag ) then 
                formatt = '(a1,a14)'
                header = "psilon normal:"
                call section(header,formatt,u,line,top%atomtype_flag)
                print*,top%atomtype_flag
        end if
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out aii section",sfilename)

!read atom vdw aii
        j=0
        do i=1, CEILING(real(top%vdwtotal)/8)
                read (unit=u, fmt='(a)',iostat=msg)  line
                if (msg < 0) call errormsg(1," seems corupted, In atom aii section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) (top%vdwpar(j+h)%aii,h=1,8)
                if (msg < 0) call read_line(top%vdwpar(j+1:j+8)%aii,line)
                j=j+8
        end do

!assign aii
if(debug .ne. "p") print*,"assigning vdw aii"
        do i=1,top%atomtotal
                do j=1,top%vdwtotal
                        if ( abs(j - top%atomtype(i)%code) .lt. 0.01 ) then
                        top%atomtype(i)%aii=top%vdwpar(j)%aii
                        end if
                end do
        end do

!find atom vdw bii
if(debug .ne. "p") print*,"reading vdw bii"
        formatt = '(a1,a17)'
        header = "qrt (Bii) normal:"
        call section(header,formatt,u,line,top%atomtype_flag)
! search for charmm type
        if ( .not. top%atomtype_flag ) then 
                formatt = '(a1,a9)'
                header = "* normal:"
                call section(header,formatt,u,line,top%atomtype_flag)
        end if
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out bii section",sfilename)

!read atom vdw bii
        j=0
        do i=1, CEILING(real(top%vdwtotal)/8)
                read (unit=u, fmt='(a)',iostat=msg)  line
                if (msg < 0) call errormsg(1," seems corupted, In atom bii section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) (top%vdwpar(j+h)%bii,h=1,8)
                if (msg < 0) call read_line(top%vdwpar(j+1:j+8)%bii,line)
                j=j+8
        end do

!assign bii
if(debug .ne. "p") print*,"assigning vdw bii"
        do i=1,top%atomtotal
                do j=1,top%vdwtotal
                        if ( abs(j - top%atomtype(i)%code) .lt. 0.01 ) then
                                top%atomtype(i)%bii = top%vdwpar(j)%bii
                        end if
                end do
        end do

!find atom vdw name
if(debug .ne. "p") print*,"reading vdw names"
        formatt = '(a11,a18)'
        header = "No. of atom types:"
        call section(header,formatt,u,line,top%atomtype_flag)
        if ( .not. top%atomtype_flag ) call errormsg(1," file have vdw with out name section",sfilename)

!read atom vdw name
        j=0
        do i=1, CEILING(real(top%vdwtotal)/8)
                read (unit=u, fmt='(a)',iostat=msg)  line
                if (msg < 0) call errormsg(1," seems corupted, In atom name section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) (top%vdwpar(j+h)%nam,h=1,8)
                if (msg < 0) call read_line(top%vdwpar(j+1:j+8)%nam,line)
                j=j+8
        end do

!assign name
if(debug .ne. "p") print*,"assigning vdw name"
        do i=1,top%atomtotal
                do j=1,top%vdwtotal
                        if ( abs(j - top%atomtype(i)%code) .lt. 0.01 ) then
                                top%atomtype(i)%nam   =  top%vdwpar(j)%nam
                                top%atomtype(i)%ci    =  0
                                top%atomtype(i)%ai    =  0
                                top%atomtype(i)%indic = "!"
                        end if
                end do
        end do
if(debug .ne. "p") print*,"deallocating parameter array"
deallocate(top%vdwpar)
end if!end atom type


!find bond type--------------------------------------------------------------------------------------------------------
        formatt = '(a19,a48)'
        header = "No. of bonds, no. of solute bonds. i - j - icode"
        call section(header,formatt,u,line,top%bondtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%bondtotal

if ( top%bondtype_flag ) then
!count how many line of bond types is there
if(debug .ne. "p") print*,"reading bonds code"
        k=5*(CEILING(real(top%bondtotal)/5))
        allocate (top%bondtype(k),stat=msg)
        call errormsg(msg,"bonds codes array could not be allocated")

!read the bond code 
        j=0
        do i=1, CEILING(real(top%bondtotal)/5)
                read (unit=u, fmt='(a)',iostat=msg)  line
                if (msg < 0) call errormsg(1," seems corupted, In bond type section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) top%bondtype(j+1)%a1, top%bondtype(j+1)%a2, top%bondtype(j+1)%code, &
                                                  top%bondtype(j+2)%a1, top%bondtype(j+2)%a2, top%bondtype(j+2)%code, &
                                                  top%bondtype(j+3)%a1, top%bondtype(j+3)%a2, top%bondtype(j+3)%code, &
                                                  top%bondtype(j+4)%a1, top%bondtype(j+4)%a2, top%bondtype(j+4)%code, &
                                                  top%bondtype(j+5)%a1, top%bondtype(j+5)%a2, top%bondtype(j+5)%code 
                j=j+5
        end do


!find bond parameter library
if(debug .ne. "p") print*,"reading bonds parameter"
        formatt = '(a11,a29)'
        header = "No. of bond codes. Parameters"
        call section(header,formatt,u,line,top%bondtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%bondpartotal
        if ( .not. top%bondtype_flag ) call errormsg(1," file have bond with no bond parameter section",sfilename)

        allocate (top%bondpar(top%bondpartotal),stat=msg)
        call errormsg(msg,"bonds parameter array could not be allocated")

!read the bond parameter
        do i=1, top%bondpartotal
                read (unit=u, fmt='(a)',iostat=msg)  line 
                read (unit=line,fmt=*,iostat=msg) top%bondpar(i)%code,top%bondpar(i)%force,top%bondpar(i)%dis
                if (msg < 0) call errormsg(1," seems corupted, In bond parameter section.",sfilename)
        end do

!assign parameter
if(debug .ne. "p") print*,"assigning bonds parameter"
        do i=1,top%bondtotal
                do j=1,top%bondpartotal

                        if ( abs(j - top%bondtype(i)%code) .lt. 0.01 ) then
                                        top%bondtype(i)%force = top%bondpar(j)%force
                                        top%bondtype(i)%dis   = top%bondpar(j)%dis
                        end if
                end do
        end do
if(debug .ne. "p") print*,"deallocating bonds parameter array"
deallocate(top%bondpar)
end if!end bond

!find angle type--------------------------------------------------------------------------------------------------------
        formatt = '(a19,a54)'
        header = "No. of angles, no. of solute angles. i - j - k - icode"
        call section(header,formatt,u,line,top%angltype_flag)
        read (unit=line,fmt=*,iostat=msg) top%angltotal

if ( top%angltype_flag ) then
!count how many line of angle types is there
if(debug .ne. "p") print*,"reading angles codes"
        k=3*(CEILING(real(top%angltotal)/3))
        allocate (top%angltype(k),stat=msg)
        call errormsg(msg,"angle codes array could not be allocated")

!read the angle type
        j=0
        do i=1, CEILING(real(top%angltotal)/3)
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In angle type section.",sfilename)
                read (unit=line,fmt=*,iostat=msg)                                                              &
                        top%angltype(j+1)%a1,top%angltype(j+1)%a2,top%angltype(j+1)%a3,top%angltype(j+1)%code, &
                        top%angltype(j+2)%a1,top%angltype(j+2)%a2,top%angltype(j+2)%a3,top%angltype(j+2)%code, &
                        top%angltype(j+3)%a1,top%angltype(j+3)%a2,top%angltype(j+3)%a3,top%angltype(j+3)%code
        j=j+3
end do

!find angle parameter library
if(debug .ne. "p") print*,"reading angles parameters"
        formatt = '(a11,a30)'
        header = "No. of angle codes. Parameters"
        call section(header,formatt,u,line,top%angltype_flag)
        read (unit=line,fmt=*,iostat=msg) top%anglpartotal
        if ( .not. top%angltype_flag ) call errormsg(1," file have angle with no angle parameter section",sfilename)

        allocate (top%anglpar(top%anglpartotal),stat=msg)
        call errormsg(msg,"angle parameter array could not be allocated")

!read angle parameter
        do i=1, top%anglpartotal
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In angle parameter section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) top%anglpar(i)%code,top%anglpar(i)%force,top%anglpar(i)%dis
        end do

!assign angle parameter
if(debug .ne. "p") print*,"assigning angles parameters"
        do i=1,top%angltotal
                do j=1,top%anglpartotal
                        if ( abs(j - top%angltype(i)%code) .lt. 0.01 ) then
                                        top%angltype(i)%force = top%anglpar(j)%force
                                        top%angltype(i)%dis   = top%anglpar(j)%dis
                        end if
                end do
        end do
if(debug .ne. "p") print*,"deallocate angles parameter array"
deallocate(top%anglpar)
end if !end angle

!find torsion type--------------------------------------------------------------------------------------------------------
        formatt = '(a19,a55)'
        header = "No. of torsions, solute torsions. i - j - k - l - icode"
        call section(header,formatt,u,line,top%torstype_flag)
        read (unit=line,fmt=*,iostat=msg) top%torstotal

if ( top%torstype_flag ) then
!count how many line of torsion types is there
if(debug .ne. "p") print*,"reading torsion codes"
        k=2*(CEILING(real(top%torstotal)/2))
        allocate (top%torstype(k),stat=msg)
        call errormsg(msg,"torsion codes array could not be allocated")

!read the torsion codes
        j=0
        do i=1, CEILING(real(top%torstotal)/2)
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In torsion type section.",sfilename)
                read (unit=line,fmt=*,iostat=msg)                                                                              &
                top%torstype(j+1)%a1,top%torstype(j+1)%a2,top%torstype(j+1)%a3,top%torstype(j+1)%a4,top%torstype(j+1)%code(1), &
                top%torstype(j+2)%a1,top%torstype(j+2)%a2,top%torstype(j+2)%a3,top%torstype(j+2)%a4,top%torstype(j+2)%code(1)
                j=j+2
        end do

!find torsion parameter library
if(debug .ne. "p") print*,"reading torsion parameters"
        formatt = '(a11,a32)'
        header = "No. of torsion codes. Parameters"
        call section(header,formatt,u,line,top%torstype_flag)
        read (unit=line,fmt=*,iostat=msg) top%torparmtotal
        if ( .not. top%torstype_flag ) call errormsg(1," file have torsion with no torsion parameter section",sfilename)

        allocate (top%torspar(top%torparmtotal),stat=msg)
        call errormsg(msg,"torsion parameter array could not be allocated")

!read torsion parameter
        do i=1, top%torparmtotal
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In torsion parameter section.",sfilename)
                read (unit=line,fmt=*,iostat=msg)                                                         &
                        top%torspar(i)%code,top%torspar(i)%force,top%torspar(i)%multip,top%torspar(i)%dis
        end do

!assign torsion parameter
if(debug .ne. "p") print*,"assigning torsion parameters"
        do i=1,top%torstotal
                do j=1,top%torparmtotal
                        if ( abs(j - top%torstype(i)%code(1) ) .lt. 0.01 ) then
                                        top%torstype(i)%force(1)  = top%torspar(j)%force
                                        top%torstype(i)%dis(1)    = top%torspar(j)%dis
                                        top%torstype(i)%multip(1) = top%torspar(j)%multip
                        end if
                end do
        end do
if(debug .ne. "p") print*,"deallocating torsion parameter array"
deallocate(top%torspar)
end if !end torsion

!find improper type------------------------------------------------------------------------------------------------------
        formatt = '(a19,a52)'
        header = "No. of impropers, solute impr. i - j - k - l - icode"
        call section(header,formatt,u,line,top%imprtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%imprtotal

if (top%imprtype_flag) then
!count how many line of improper types is there
if(debug .ne. "p") print*,"reading improper codes"
        k=2*(CEILING(real(top%imprtotal)/2))
        allocate (top%imprtype(k),stat=msg)
        call errormsg(msg,"array could not be allocated. too many torsions")

!read the improper codes
        j=0
        do i=1, CEILING(real(top%imprtotal)/2)
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In improper type section.",sfilename)
                read (unit=line,fmt=*,iostat=msg)                                                                           &
                top%imprtype(j+1)%a1,top%imprtype(j+1)%a2,top%imprtype(j+1)%a3,top%imprtype(j+1)%a4,top%imprtype(j+1)%code, &
                top%imprtype(j+2)%a1,top%imprtype(j+2)%a2,top%imprtype(j+2)%a3,top%imprtype(j+2)%a4,top%imprtype(j+2)%code
                j=j+2
        end do

!find improper parameter library
if(debug .ne. "p") print*,"reading improper parameters"
        formatt = '(a19,a63)'
        header = "No. of improper codes, type (1=harmonic,2=periodic). Parameters"
        call section(header,formatt,u,line,top%imprtype_flag)
        read (unit=line,fmt=*,iostat=msg) top%imprpartotal
        if ( .not. top%imprtype_flag ) call errormsg(1," file have improper with no improper parameter section",sfilename)

        allocate (top%imprpar(top%imprpartotal),stat=msg)
        call errormsg(msg,"improper parameter array could not be allocated")

!read improper parameter
        do i=1, top%imprpartotal
                read (unit=u, fmt='(a)',iostat=msg) line
                if (msg < 0) call errormsg(1," seems corupted, In improper parameter section.",sfilename)
                read (unit=line,fmt=*,iostat=msg) top%imprpar(i)%code,top%imprpar(i)%force,top%imprpar(i)%dis
        end do

!assign improper parameter
if(debug .ne. "p") print*,"assigning improper parameters"
        do i=1,top%imprtotal
                do j=1,top%imprpartotal
                        if ( abs(j - top%imprtype(i)%code) .lt. 0.01 ) then
                                        top%imprtype(i)%force = top%imprpar(j)%force
                                        top%imprtype(i)%dis   = top%imprpar(j)%dis
                        end if
                end do
        end do
if(debug .ne. "p") print*,"deallocating improper parameter array"
deallocate(top%imprpar)
end if !end improper
 close(u)

end subroutine top_reader
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             subroutine read_line_int                                        !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine read_line_int(var,line)
 character(*),intent(in)                                 ::        line
 integer(4),intent(inout)                                ::        var(:)
 integer(4)                                              ::        i,msg,h

        do i=1,size(var(:)) 
                read(line,*,iostat=msg),(var(h),h=1,i)
                if (msg < 0 ) return
        end do
end subroutine read_line_int
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             subroutine read_line_real                                       !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine read_line_real(var,line)
 character(*),intent(in)                          ::        line
 real(4),intent(inout)                            ::        var(:)
 integer(4)                                       ::        i,msg,h

        do i=1,size(var(:)) 
                read(line,*,iostat=msg),(var(h),h=1,i)
                if (msg < 0 ) return
        end do
end subroutine read_line_real
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             subroutine read_line_char                                       !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine read_line_char(var,line)
 character(*),intent(in)                           ::        line
 character(*),intent(inout)                        ::        var(:)
 integer(4)                                        ::        i,msg,h

        do i=1,size(var(:)) 
                read(line,*,iostat=msg),(var(h),h=1,i)
                if (msg < 0 ) return
        end do
end subroutine read_line_char
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 subroutine section                                          !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine section(header,formatt,u,line,flag)
        integer(4),intent(in)                          ::        u
        character(*),intent(in)                        ::        header,formatt
        character(*),intent(inout)                     ::        line
        logical,intent(inout)                          ::        flag
        integer(4)                                     ::        msg=0
        character(100)                                 ::        hit1,hit2
        flag = .false.
        msg = 0
        rewind u
        do while ( msg .ge. 0 )
                read (unit=u, fmt='(a)',iostat=msg)  line
                read ( line ,formatt), hit1,hit2
                if (trim(adjustl(hit2)) .eq. trim(adjustl(header))) then
                        flag = .true. 
                        print*,"Section [",trim(adjustl(header)),"] was found..."
                        return
                elseif (msg .lt. 0) then
                        flag = .false.
                        line = "0"
                        print*, "End of file reached no [",trim(adjustl(header)),"] was found..."
                end if
        end do

end subroutine section
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 subroutine errormsg                                         !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine errormsg(msg,emsg,filenam)
        integer(4)                ,intent(in)                   ::        msg
        character(len=*),intent(in)                             ::        emsg
        character(len=*),intent(in),optional                    ::        filenam
        if (present(filenam)) then
                if (msg > 0) then
                        print*, filenam, emsg,".   stoping program..."
                        stop
                end if
        else
                if (msg > 0) then
                        print*, emsg,".   stoping program..."
                        stop
                end if
        end if
end subroutine errormsg
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             subroutine option_read                                          !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine option_read
 integer                                                   ::      i
        option = trim(adjustl(option))
        do i=1,len_trim(adjustl(option))
                select case (option(i:i))
                        case ("p")
                                debug = "p"
                                print*, "extra print flag: ON"
                        case ("s")
                                split = "s"
                                print*, "Split flag:       ON"
                        case default
                                print*, "extra print flag: OFF"
                                print*, "Split flag:       OFF"
                end select
        end do

end subroutine option_read
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             subroutine instruction                                          !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine instruction
        print*,"Input file name not found"
        print*,"qfeper input ps"
        print*,"p and s are optional for printing more and spliting the fep states"
        print*,"qfeper h"
        print*,"Will print out the input format"
        stop
end subroutine instruction
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 subroutine help                                             !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine help
        print*,"Order of line is important-no empty line"
        print*,"line number     function"
        print*,"----------------------------------------"
        print*,"1               number of states"
        print*,"2               first topology file"
        print*,".               "
        print*,".               "
        print*,"n               last topology file"
        print*,"n+1   n+1   ... number of qatoms"
        print*,"n+2   n+2   ... topology number of qatom 1 in different topology files"
        print*,".      .        "
        print*,".      .        "
        print*,".      .        "
        print*,"n+j+1 n+j+1 ... topology number of qatom j in different topology files"
        print*,"----------------------------------------"
        print*,"first topology file is the reference state"
        stop
end subroutine help
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                      freefile : get the free unit number                                    !
!-----------------------------------------------------------------------------------------------------------------------------!
integer function freefile()

        integer(4)                                        ::        u
        logical                                           ::        used

        do u = 20, 999                
                inquire(unit=u, opened=used)
                if(.not. used)  then
                        freefile = u
                        return
                end if
        end do

        !if we get here then we're out of unit numbers
        write(*,20)
        stop 

20        format('ERROR: Failed to find an unused unit number')

end function freefile
!=============================================================================================================================!
end module top_parser


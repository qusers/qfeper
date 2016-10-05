!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        By Masoud Kazemi 2013/08/09         !
!           revised 2014/06/02               !
!   qfeper program for writting fep files    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!                    TODO                    
!the parm type is redundant. not urgent

!                    DONE
!the problem of torsions repetition solved.
!split fep file implemented.

program qfeper

use top_parser
use array_analysis

implicit none

!variable decleration-----------------------------------------------------------------------------------------------
 type(top_type),allocatable                             ::        top(:),qtop(:)
 character(len=50)                                      ::        ifilename,output
 character(len=50),allocatable                          ::        sfilename(:)
 integer(4)                                             ::        msg,nstate,nqatom,istate
 integer(4)                                             ::        u,i,j,k,h,tors_count,nsplit,bond_count,angl_count,impr_count
 integer(4),allocatable                                 ::        qstatemap(:,:),a_kapl(:,:),t_kapl(:,:),i_kapl(:,:)
 type(q_lib_type)                                       ::        q_lib
 type(q_lib_type)                                       ::        p_lib
 type(bondtype_type),allocatable                        ::        soft_pairs(:,:)

!Reading instruction file ------------------------------------------------------------------------------------------
        debug = "n"
        split = "n"
        call getarg(1,ifilename)
        call getarg(2,option)
        call option_read
        if(len_trim(ifilename) .eq. 0 ) call instruction
        if(trim(adjustl(ifilename)) .eq. "h" ) call help

        u = freefile()
        open (file=ifilename, unit=u, form="formatted", access="sequential", status="old", iostat=msg, action="read")
        call errormsg(msg," file cant be opend.",ifilename)

!read first record
        read (unit=u, fmt=*,iostat=msg)  nstate
        if ( nstate .eq. 0 .or. nstate .eq. 1 ) call errormsg(msg,"At least two states is needed.")
        allocate (sfilename(nstate),stat=msg) 
        call errormsg(msg,"State file name could not be allocated.")
        if (debug .eq. "p") write(*,*),"--------------------------------------------------------------------------"
!read state files
        do i=1,nstate
                read(unit=u, fmt=*,iostat=msg) sfilename(i)
                if(len_trim(sfilename(i)) .eq. 0 ) call errormsg(1,"No State file was specifies.")
        if (debug .eq. "p") write(*,'(a11,i2,a3,a50)'),"State file ",i,":  ",sfilename(i)
        end do
        if (debug .eq. "p") write(*,*),"--------------------------------------------------------------------------"
!read qatom equivalency
        read (unit=u, fmt=*,iostat=msg) nqatom
        if(nqatom .eq. 0 ) call errormsg(1,"No q atom was found.")
        allocate (qstatemap(nqatom,nstate),stat=msg)
        call errormsg(msg,"qatoms could not be allocated")
        if (debug .eq. "p") write(*,1) ('state',h, h=1,nstate)
        do i=1, nqatom
                read (unit=u,fmt=*,iostat=msg) (qstatemap(i,h),h=1,nstate)
                do j=1,nstate
                        if(qstatemap(i,j) .eq. 0 ) call errormsg(1,"qatom number cant be 0.")
                end do
        if (debug .eq. "p") write(*,2),(qstatemap(i,h), h=1,nstate)
        end do
1        format('qatom number in ' , 7(1x,a5,i2))
2        format( 19X,7(i4,3X))
        close (u)
        if (debug .eq. "p") write(*,*),"--------------------------------------------------------------------------"
!allocate and read the top ---------------------------------------------------------------------------------------
        allocate (top(nstate),stat=msg)
        call errormsg(msg,"topology variable could not be allocated.")
        do i=1,nstate
                call top_reader(sfilename(i),top(i)) 
        end do

! make q atom structure------------------------------------------------------------------------------------------
        call q_construct(top,qstatemap,qtop)
        deallocate(top)
! make q library-------------------------------------------------------------------------------------------------
        call q_library(qtop,q_lib)
        deallocate(qtop)
!make parameter libray-------------------------------------------------------------------------------------------
        call p_library(q_lib,p_lib)

!find bond breaking and take care of soft core and coupling------------------------------------------------------
        call c_library(q_lib,p_lib,soft_pairs,a_kapl,t_kapl,i_kapl)

!set the soft_pairs matrix to q number---------------------------------------------------------------------------
                do i=1,size(soft_pairs(:,1))
                        do j=1,size(qstatemap(:,1))
                                if (soft_pairs(i,1)%a1 .eq. qstatemap(j,1)) then
                                        soft_pairs(i,1)%a1 = j
                                end if
                                if (soft_pairs(i,1)%a2 .eq. qstatemap(j,1)) then
                                        soft_pairs(i,1)%a2 = j
                                end if
                        end do
                end do
        close(u)
!print extra on demand-------------------------------------------------------------------------------------------
if(debug .eq. "p") then
        write(nstate_char,'(i4)') nstate
        if(q_lib%atm_flag) then
        print*,"Q-atom-type"
        write(qformat,*)       "(i4,X,",trim(adjustl(nstate_char)),"(i4,X),X,",trim(adjustl(nstate_char)),&
                            "(a7,X),",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),&
                        "(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"

        do i=1,size(q_lib%atm(:,1))
        write(*,qformat),i,(q_lib%atm(i,h)%atmno,h=1,nstate),(q_lib%atm(i,h)%nam,h=1,nstate),(q_lib%atm(i,h)%charge,h=1,nstate) &
                ,(q_lib%atm(i,h)%aii,h=1,nstate),(q_lib%atm(i,h)%bii,h=1,nstate),(q_lib%atm(i,h)%mass,h=1,nstate)
        end do
        write(*,*),"--------------------------------------------------------------------------"
        end if

        if(q_lib%bond_flag) then
        print*,"Q-bond"
        write(qformat,*) "(i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
        do i=1,size(q_lib%bond(:,1))
        write(*,qformat),i,q_lib%bond(i,1)%a1,q_lib%bond(i,1)%a2,(q_lib%bond(i,h)%force,h=1,nstate), &
                                                                                        (q_lib%bond(i,h)%dis,h=1,nstate)
        end do
        write(*,*),"--------------------------------------------------------------------------"
        end if

        if(q_lib%angl_flag) then
        print*,"Q-angle"
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
        do i=1,size(q_lib%angl(:,1))
        write(*,qformat),i,q_lib%angl(i,1)%a1,q_lib%angl(i,1)%a2,q_lib%angl(i,1)%a3, &
                                                        (q_lib%angl(i,h)%force,h=1,nstate), (q_lib%angl(i,h)%dis,h=1,nstate)
        end do
        write(*,*),"--------------------------------------------------------------------------"
        end if

        if(q_lib%tors_flag) then
        print*,"Q-torsion"
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)), &
                                                                        "(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
        do i=1,size(q_lib%tors(:,1))
        write(*,qformat),i,q_lib%tors(i,1)%a1,q_lib%tors(i,1)%a2,q_lib%tors(i,1)%a3,q_lib%tors(i,1)%a4, &
         (q_lib%tors(i,h)%multip(1),h=1,nstate),(q_lib%tors(i,h)%force(1),h=1,nstate),(q_lib%tors(i,h)%dis(1),h=1,nstate)
        write(*,qformat),i,q_lib%tors(i,1)%a1,q_lib%tors(i,1)%a2,q_lib%tors(i,1)%a3,q_lib%tors(i,1)%a4, &
         (q_lib%tors(i,h)%multip(2),h=1,nstate),(q_lib%tors(i,h)%force(2),h=1,nstate),(q_lib%tors(i,h)%dis(2),h=1,nstate)
        write(*,qformat),i,q_lib%tors(i,1)%a1,q_lib%tors(i,1)%a2,q_lib%tors(i,1)%a3,q_lib%tors(i,1)%a4, &
         (q_lib%tors(i,h)%multip(3),h=1,nstate),(q_lib%tors(i,h)%force(3),h=1,nstate),(q_lib%tors(i,h)%dis(3),h=1,nstate)
        write(*,*),"#"
        end do
        write(*,*),"--------------------------------------------------------------------------"
        end if

        if(q_lib%impr_flag) then
        print*,"Q-improper"
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)), &
                                                                                                                "(f8.3,X))"
        do i=1,size(q_lib%impr(:,1))
        write(*,qformat),i,q_lib%impr(i,1)%a1,q_lib%impr(i,1)%a2,q_lib%impr(i,1)%a3,q_lib%impr(i,1)%a4, &
                                                        (q_lib%impr(i,h)%force,h=1,nstate),(q_lib%impr(i,h)%dis,h=1,nstate)
        end do
        write(*,*),"--------------------------------------------------------------------------"
        end if
end if

!writing fep file!------------------------------------------------------------------------------------------------------------
        u=freefile()
        istate = 0
        write(nstate_char,'(i4)') nstate
        write(qformat,*) "(",trim(adjustl(nstate_char)),"(i1),a4)"
        write(output,qformat) (istate + h,h=1,nstate),".fep"
        open (file=trim(adjustl(output)), unit=u, form="formatted", access="sequential",status="replace", &
                                                                                                iostat=msg, action="write")
        call errormsg(msg,"output could not be oppend")

!write the atom section(if here the atomtype_flag must be true)
        write (u,'(a5)'),"[FEP]"
        write (u,'(a7,i2)'),"states",nstate
        write (u,'(a9)'),"offset    "

if(q_lib%atm_flag) then
        write(qformat,*) "(3X,i4,3X,i4,a3,",trim(adjustl(nstate_char)),"(a7,1X))"
                write (u,'(/,a7)'),"[atoms]"
                do i=2,size(q_lib%atm(:,1))
                        write(u,qformat),i-1 ,q_lib%atm(i,1)%atmno,"  !",(q_lib%atm(i,h)%nam,h=1,nstate)
                end do

     write(qformat,*) "(3X,i4,3X,",trim(adjustl(nstate_char)),"(f9.4,2X)",",a3,i4,2X,",trim(adjustl(nstate_char)),"(a7,1X))"
                write (u,'(/,a16)'),"[change_charges]"
                do i=2,size(q_lib%atm(:,1))
        write(u,qformat),i-1 ,(q_lib%atm(i,h)%charge,h=1,nstate),"  !",q_lib%atm(i,1)%atmno,(q_lib%atm(i,h)%nam,h=1,nstate)
                end do

        write(qformat,*) "(1X,a7,1X,7(f8.3,2X),a3,i4)"
                write (u,'(/,a12)'),"[atom_types]"
                do i=2,size(p_lib%atm(:,1))
                        write(u,qformat),p_lib%atm(i,1)%nam,p_lib%atm(i,1)%aii,p_lib%atm(i,1)%bii,p_lib%atm(i,1)%ci &
                        ,p_lib%atm(i,1)%ai,p_lib%atm(i,1)%aii,p_lib%atm(i,1)%bii,p_lib%atm(i,1)%mass ,"  !",i-1 
                end do

        write(qformat,*) "(1X,i4,3X,",trim(adjustl(nstate_char)),"(a7,1X),a3,X,i4)"
                write (u,'(/,a14)'),"[change_atoms]"
                do i=2,size(q_lib%atm(:,1))
              write(u,qformat),i-1,(p_lib%atm(q_lib%atm(i,h)%code,1)%nam,h=1,nstate),"  !",q_lib%atm(i,1)%atmno
                end do
end if

if(q_lib%bond_flag) then
        write(qformat,*) "(1X,i4,3X,i4)"
                write (u,'(/,a12)'),"[soft_pairs]"
                do i=1,size(soft_pairs(:,1))
                        write(u,qformat),soft_pairs(i,1)%a1,soft_pairs(i,1)%a2
                end do
end if

        write (u,'(/,a16)'),"[excluded_pairs]"
        write (u,'(a20)'),"! i   j   s1  s2    "

if(q_lib%bond_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f7.2,3X,f8.3)"
        write (u,'(/,a12)'),"[bond_types]"
        do i=2,size(p_lib%bond(:,1))
        write(u,qformat),i-1 ,p_lib%bond(i,1)%force/8, 2.00 ,p_lib%bond(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a14)'),"[change_bonds]"
        do i=1,size(q_lib%bond(:,1))
        write(u,qformat),q_lib%bond(i,1)%a1,q_lib%bond(i,1)%a2,q_lib%bond(i,1:nstate)%code-1,"  !",i
        end do
end if

if(q_lib%angl_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3)"
        write (u,'(/,a13)'),"[angle_types]"
        do i=2,size(p_lib%angl(:,1))
        write(u,qformat),i-1 ,p_lib%angl(i,1)%force,p_lib%angl(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a15)'),"[change_angles]"
        do i=1,size(q_lib%angl(:,1))
        write(u,qformat),q_lib%angl(i,1)%a1,q_lib%angl(i,1)%a2,q_lib%angl(i,1)%a3,q_lib%angl(i,1:nstate)%code-1,"  !",i
        end do
end if

if(q_lib%tors_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3,3X,f8.3)"
        write (u,'(/,a15)'),"[torsion_types]"
        do i=2,size(p_lib%tors(:,1))
        write(u,qformat),i-1 ,p_lib%tors(i,1)%force(1),p_lib%tors(i,1)%multip(1),p_lib%tors(i,1)%dis(1)
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4,X,i4)"
        write (u,'(/,a17)'),"[change_torsions]"
        tors_count = 0
        do i=1,size(q_lib%tors(:,1))
        !find the MAX nterms
        do istate=1,nstate
                if (q_lib%tors(i,istate)%nterm .gt. q_lib%tors(i,1)%nterm) then
                        q_lib%tors(i,1)%nterm = q_lib%tors(i,istate)%nterm
                end if
        end do
        do j=1,q_lib%tors(i,1)%nterm
        tors_count = tors_count + 1 
        write(u,qformat),q_lib%tors(i,1)%a1,q_lib%tors(i,1)%a2,q_lib%tors(i,1)%a3,q_lib%tors(i,1)%a4 &
                                                        ,q_lib%tors(i,1:nstate)%code(j)-1,"  !",i,tors_count
        !use %indice array to index the torsion number in output for coupling 
        q_lib%tors(i,1:nstate)%line_counter(j) = tors_count
        end do
        end do
end if

if(q_lib%impr_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3)"
        write (u,'(/,a16)'),"[improper_types]"
        do i=2,size(p_lib%impr(:,1))
        write(u,qformat),i-1 ,p_lib%impr(i,1)%force,p_lib%impr(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a18)'),"[change_impropers]"
        do i=1,size(q_lib%impr(:,1))
        write(u,qformat),q_lib%impr(i,1)%a1,q_lib%impr(i,1)%a2,q_lib%impr(i,1)%a3,q_lib%impr(i,1)%a4 &
                                                                                ,q_lib%impr(i,1:nstate)%code-1,"  !",i
        end do
end if

if(q_lib%bond_flag) then
if(q_lib%angl_flag) then
        if (size(a_kapl(:,1)) .gt. 0) then
        write (u,'(/,a17)'),"[angle_couplings]"
        do i=1,size(a_kapl(:,1))
                        write(u,'(i3,3X)',advance="no"),a_kapl(i,1)
                        do k=2,size(a_kapl(i,:))
                                if (a_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),a_kapl(i,k)
                        end do
                        write(u,'(3X)')
        end do
        end if
end if

if(q_lib%tors_flag) then
        if (size(t_kapl(:,1)) .gt. 0) then
        !get the number of guessed couplings
        write (u,'(/,a20)'),"[torsion_couplings]"
        do i=1,size(t_kapl(:,1))
                !t_kapl 's first element is the torsion index in q_lib 
                do j=1,q_lib%tors(t_kapl(i,1),1)%nterm
                     write(u,'(i3,3X)',advance="no"),q_lib%tors(t_kapl(i,1),1)%line_counter(j)
                     do k=2,size(a_kapl(i,:))
                              if (t_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),t_kapl(i,k)
                     end do
                     write(u,'(3X)')
                end do
        end do
        end if
end if

if(q_lib%impr_flag) then
        if (size(i_kapl(:,1)) .gt. 0) then
        write (u,'(/,a22)'),"[improper_couplings]"
        do i=1,size(i_kapl(:,1))
             write(u,'(i3,3X)',advance="no"),i_kapl(i,1)
             do k=2,size(a_kapl(i,:))
                     if (i_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),i_kapl(i,k)
             end do
             write(u,'(3X)')
        end do
        end if

end if
end if
        write (u,'(/,a16)'),"[off_diagonals]"
        write (u,*),"!  statei   statej    i    j     Aij    Mij"

        write (u,'(/)')
        close(u)

!writing split fep files!------------------------------------------------------------------------------------------------------------
!in split mode the order of print might change unpreditebly. this cause problem in coupling matrix. adding a line_counter to varible 
!defenition to keep track of print order. %indic will be the p_library order, %line_counter will be the print order (coupling map is 
!based on q_lib order (indices)). the 0-0 terms should be there since they are defined in 1th .top, trmoving conditions

if (split .eq. "s") then
do nsplit = 1,nstate-1 !n state need n-1 split fep
        u=freefile()
        write(nstate_char,'(i4)') 2
        write(qformat,*) "(",trim(adjustl(nstate_char)),"(i1),a4)"
        write(output,qformat) (nsplit + h ,h=0,1),".fep"
        open (file=trim(adjustl(output)), unit=u, form="formatted", access="sequential",status="replace", &
                                                                                                iostat=msg, action="write")
        call errormsg(msg,"output could not be oppend")

!write the atom section(if here the atomtype_flag must be true)
        write (u,'(a5)'),"[FEP]"
        write (u,'(a7,i2)'),"states",2
        write (u,'(a9)'),"offset    "

if(q_lib%atm_flag) then
        write(qformat,*) "(3X,i4,3X,i4,a3,",trim(adjustl(nstate_char)),"(a7,1X))"
                write (u,'(/,a7)'),"[atoms]"
                do i=2,size(q_lib%atm(:,1))
                        write(u,qformat),i-1 ,q_lib%atm(i,1)%atmno,"  !",(q_lib%atm(i,nsplit + h)%nam,h=0,1)
                end do

     write(qformat,*) "(3X,i4,3X,",trim(adjustl(nstate_char)),"(f9.4,2X)",",a3,i4,2X,",trim(adjustl(nstate_char)),"(a7,1X))"
                write (u,'(/,a16)'),"[change_charges]"
                do i=2,size(q_lib%atm(:,1))
                write(u,qformat),i-1 ,(q_lib%atm(i,nsplit + h)%charge,h=0,1),"  !",q_lib%atm(i,1)%atmno,& 
                                                                                        (q_lib%atm(i,nsplit + h)%nam,h=0,1)
                end do

        write(qformat,*) "(1X,a7,1X,7(f8.3,2X),a3,i4)"
                write (u,'(/,a12)'),"[atom_types]"
                do i=2,size(p_lib%atm(:,1))
                        write(u,qformat),p_lib%atm(i,1)%nam,p_lib%atm(i,1)%aii,p_lib%atm(i,1)%bii,p_lib%atm(i,1)%ci &
                        ,p_lib%atm(i,1)%ai,p_lib%atm(i,1)%aii,p_lib%atm(i,1)%bii,p_lib%atm(i,1)%mass ,"  !",i-1 
                end do

        write(qformat,*) "(1X,i4,3X,",trim(adjustl(nstate_char)),"(a7,1X),a3,X,i4)"
                write (u,'(/,a14)'),"[change_atoms]"
                do i=2,size(q_lib%atm(:,1))
                write(u,qformat),i-1,(p_lib%atm(q_lib%atm(i,nsplit + h)%code,1)%nam,h=0,1),"  !",q_lib%atm(i,1)%atmno
                end do
end if

if(q_lib%bond_flag) then
        write(qformat,*) "(1X,i4,3X,i4)"
                write (u,'(/,a12)'),"[soft_pairs]"
                do i=1,size(soft_pairs(:,1))
                        !only declare soft_pairs if its 1-0 or 0-1 (breaking/forming in 2 consecetive states)
                        if ( abs(sum(soft_pairs(i,nsplit:nsplit+1)%code) - 1) .lt. 0.001 ) then 
                        write(u,qformat),soft_pairs(i,1)%a1,soft_pairs(i,1)%a2
                        end if
                end do
end if

        write (u,'(/,a16)'),"[excluded_pairs]"
        write (u,'(a20)'),"! i   j   s1  s2    "

if(q_lib%bond_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f7.2,3X,f8.3)"
        write (u,'(/,a12)'),"[bond_types]"
        do i=2,size(p_lib%bond(:,1))
        write(u,qformat),i-1 ,p_lib%bond(i,1)%force/8, 2.00 ,p_lib%bond(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a14)'),"[change_bonds]"
        bond_count = 0
        do i=1,size(q_lib%bond(:,1))
        q_lib%bond(i,1)%line_counter = 0                                        !reset for new fep file
!                if ( sum(q_lib%bond(i,nsplit:nsplit+1)%code-1) .gt. 0 ) then    !if the codes are not 0-0"the damn thing must be there :)"
                        bond_count = bond_count + 1
                        q_lib%bond(i,1)%line_counter = bond_count               !keep track which line its printed
                        q_lib%bond(i,1)%indic = i                               !keep track of q_lib order
        write(u,qformat),q_lib%bond(i,1)%a1,q_lib%bond(i,1)%a2,(q_lib%bond(i,nsplit + h)%code-1,h=0,1),"  !",bond_count
!                end if
        end do
end if

if(q_lib%angl_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3)"
        write (u,'(/,a13)'),"[angle_types]"
        do i=2,size(p_lib%angl(:,1))
        write(u,qformat),i-1 ,p_lib%angl(i,1)%force,p_lib%angl(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a15)'),"[change_angles]"
        angl_count = 0
        do i=1,size(q_lib%angl(:,1))
        q_lib%angl(i,1)%line_counter = 0                                        !reset for new fep file
!                if ( sum(q_lib%angl(i,nsplit:nsplit+1)%code-1) .gt. 0 ) then    !if the codes are not 0-0
                        angl_count = angl_count + 1
                        q_lib%angl(i,1)%line_counter = angl_count               !keep track which line its printed
        write(u,qformat),q_lib%angl(i,1)%a1,q_lib%angl(i,1)%a2,q_lib%angl(i,1)%a3,(q_lib%angl(i,nsplit + h)%code-1,h=0,1) &
                                                                                                            ,"  !",angl_count
!                end if
        end do
end if

if(q_lib%tors_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3,3X,f8.3)"
        write (u,'(/,a15)'),"[torsion_types]"
        do i=2,size(p_lib%tors(:,1))
        write(u,qformat),i-1 ,p_lib%tors(i,1)%force(1),p_lib%tors(i,1)%multip(1),p_lib%tors(i,1)%dis(1)
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4,X,i4)"
        write (u,'(/,a17)'),"[change_torsions]"
        tors_count = 0
        do i=1,size(q_lib%tors(:,1))
        !find the MAX nterms
        do istate=1,nstate
                if (q_lib%tors(i,istate)%nterm .gt. q_lib%tors(i,1)%nterm) then
                        q_lib%tors(i,1)%nterm = q_lib%tors(i,istate)%nterm
                end if
        end do

        !reset for new fep file for each tors triplet
        q_lib%tors(i,1)%line_counter(:) = 0
        !if not all codes of a triplet = 0 in 2 successive states
!        if ( (sum(q_lib%tors(i,nsplit)%code(:)-1) + sum(q_lib%tors(i,nsplit+1)%code(:)-1) ) .ne. 0 ) then 
                do j=1,q_lib%tors(i,1)%nterm
                        tors_count = tors_count + 1 
                        write(u,qformat),q_lib%tors(i,1)%a1,q_lib%tors(i,1)%a2,q_lib%tors(i,1)%a3,q_lib%tors(i,1)%a4 &
                                                        ,q_lib%tors(i,nsplit:nsplit+1)%code(j)-1,"  !",i,tors_count
                        !use %line_counter to index the torsion number in output for coupling 
                        q_lib%tors(i,1:nstate)%line_counter(j) = tors_count
                end do
!        end if
        end do
end if

if(q_lib%impr_flag) then
        write(qformat,*) "(i4,3X,f8.3,3X,f8.3)"
        write (u,'(/,a16)'),"[improper_types]"
        do i=2,size(p_lib%impr(:,1))
        write(u,qformat),i-1 ,p_lib%impr(i,1)%force,p_lib%impr(i,1)%dis
        end do

        write(qformat,*) "(i4,3X,i4,3X,i4,3X,i4,3X,",trim(adjustl(nstate_char)),"(i4,3X),a3,X,i4)"
        write (u,'(/,a18)'),"[change_impropers]"
        impr_count = 0
        do i=1,size(q_lib%impr(:,1))
        q_lib%impr(i,1)%line_counter = 0                                        !reset for new fep file
!                if ( sum(q_lib%impr(i,nsplit:nsplit+1)%code-1) .gt. 0 ) then    !if the codes are not 0-0
                        impr_count = impr_count + 1
                        q_lib%impr(i,1)%line_counter = impr_count               !keep track which line its printed
                write(u,qformat),q_lib%impr(i,1)%a1,q_lib%impr(i,1)%a2,q_lib%impr(i,1)%a3,q_lib%impr(i,1)%a4 &
                                                        ,(q_lib%impr(i,nsplit + h)%code-1,h=0,1),"  !",impr_count
!                end if
        end do
end if

if(q_lib%bond_flag) then
if(q_lib%angl_flag) then
        if (size(a_kapl(:,1)) .gt. 0) then
        write(qformat,*) "(",size(a_kapl(1,:)),"(i4,3X))"
        write (u,'(/,a17)'),"[angle_couplings]"
        do i=1,size(a_kapl(:,1))
                !if an angle is has %line_counter = 0 it means it not exist in this fep
                if ( q_lib%angl(a_kapl(i,1),1)%line_counter .ne. 0 ) then
                        write(u,'(i3,3X)',advance="no"),q_lib%angl(a_kapl(i,1),1)%line_counter
                        do k=2,size(a_kapl(i,:))
                         if (a_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),q_lib%bond(a_kapl(i,k),1)%line_counter
                        end do
                write(u,'(3X)')
                end if
        end do
        end if
end if

if(q_lib%tors_flag) then
        if (size(t_kapl(:,1)) .gt. 0) then
        !get the number of guessed couplings

        write (u,'(/,a20)'),"[torsion_couplings]"
        do i=1,size(t_kapl(1,:)) !over each tors
              if ( sum(q_lib%tors(t_kapl(i,1),1)%line_counter(:)) .ne. 0 ) then
                    do j=1,q_lib%tors(t_kapl(i,1),1)%nterm !over elements of each tors
                           write(u,'(i3,3X)',advance="no"),q_lib%tors(t_kapl(i,1),1)%line_counter(j)
                           do k=2,size(a_kapl(i,:))
                              if (t_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),q_lib%bond(t_kapl(i,k),1)%line_counter
                           end do
                           write(u,'(3X)')
                    end do
             end if
        end do
        end if
end if

if(q_lib%impr_flag) then
        if (size(i_kapl(:,1)) .gt. 0) then
        write (u,'(/,a22)'),"[improper_couplings]"
        do i=1,size(i_kapl(:,1))
                !if an impr is has %line_counter = 0 it means it not exist in this fep
                if ( q_lib%impr(i_kapl(i,1),1)%line_counter .ne. 0 ) then
                        write(u,'(i3,3X)',advance="no"),q_lib%impr(i_kapl(i,1),1)%line_counter
                        do k=2,size(a_kapl(i,:))
                         if (i_kapl(i,k) .ne. 0) write(u,'(i3,3X)',advance="no"),q_lib%bond(i_kapl(i,k),1)%line_counter
                        end do
                write(u,'(3X)')
                end if
        end do
        end if

end if
end if
        write (u,'(/,a16)'),"[off_diagonals]"
        write (u,*),"!  statei   statej    i    j     Aij    Mij"

        write (u,'(/)')
        close(u)
end do
end if

        deallocate(sfilename)

end program qfeper

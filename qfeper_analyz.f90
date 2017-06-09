!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        By Masoud Kazemi 2013/08/09         !
!           revised 2014/06/02               !
!Array_analysis module for writting fep files!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module array_analysis

use top_parser
implicit none

!-----------------------------------------------------------------------------------------------------------------------------!
!                                                                 type decleration                                            !
!-----------------------------------------------------------------------------------------------------------------------------!
!qtoms library
type q_lib_type
 type(atomtype_type),allocatable                        ::        atm(:,:)
 type(bondtype_type),allocatable                        ::        bond(:,:)
 type(angltype_type),allocatable                        ::        angl(:,:)
 type(torstype_type),allocatable                        ::        tors(:,:)
 type(imprtype_type),allocatable                        ::        impr(:,:)
 logical                                                ::        atm_flag
 logical                                                ::        bond_flag
 logical                                                ::        angl_flag
 logical                                                ::        tors_flag
 logical                                                ::        impr_flag
end type q_lib_type

public                 ::        operator(.same.),operator(.eqt.),q_construct
private                ::        bond_same_bond,angl_same_angl,tors_same_tors,atom_eqt_atom


interface operator(.same.)
        module procedure bond_same_bond,angl_same_angl,tors_same_tors,impr_same_impr
end interface operator(.same.)

interface operator(.eqt.)
        module procedure bond_eqt_bond,angl_eqt_angl,tors_eqt_tors,impr_eqt_impr,atom_eqt_atom
end interface operator(.eqt.)

interface operator(.incl.)
        module procedure angl_incl_bond,tors_incl_bond,impr_incl_bond
end interface operator(.incl.)


contains
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             q_construct                                                     !
!-----------------------------------------------------------------------------------------------------------------------------!
!remove non qtom and rename everything according to first state
subroutine q_construct(top,qstatemap,qtop)
        type(top_type),intent(inout)                                   :: top(:)
        integer(4),intent(in)                                          :: qstatemap(:,:)
        type(top_type),allocatable,intent(inout)                       :: qtop(:)
!local
        type(atomtype_type),allocatable                                :: qatmtype_tm(:)
        type(bondtype_type),allocatable                                :: bondtype_tm(:)
        type(angltype_type),allocatable                                :: angltype_tm(:) 
        type(torstype_type),allocatable                                :: torstype_tm(:)
        type(imprtype_type),allocatable                                :: imprtype_tm(:)
        integer(4)                                                     :: istate,i,j,k,nq,nstate
!----------------------------------------------------------------------------------------------------------------------------

        nstate = size(qstatemap(1,:))
        allocate(qtop(nstate))
!qatom types-----------------------------------------------------------------------------------------------------------------
if(debug .ne. "p") print*,"finding q atoms"
        do istate=1,nstate !over states
        qtop(istate)%atomtype_flag = top(istate)%atomtype_flag
        allocate(qatmtype_tm(top(istate)%atomtotal))

        nq = 0
        do i=1,size(qstatemap(:,istate)) !over qatoms
                do j=1,top(istate)%atomtotal!over top atoms
                        if ( abs(top(istate)%atomtype(j)%atmno - qstatemap(i,istate)) .lt. 0.01 ) then
                                nq = nq+1
                                qatmtype_tm(nq)%atmno  = qstatemap(i,1)
                                qatmtype_tm(nq)%code   = top(istate)%atomtype(j)%code
                                qatmtype_tm(nq)%nam    = top(istate)%atomtype(j)%nam
                                qatmtype_tm(nq)%indic  = top(istate)%atomtype(j)%indic
                                qatmtype_tm(nq)%charge = top(istate)%atomtype(j)%charge
                                qatmtype_tm(nq)%mass   = top(istate)%atomtype(j)%mass
                                qatmtype_tm(nq)%aii    = top(istate)%atomtype(j)%aii
                                qatmtype_tm(nq)%bii    = top(istate)%atomtype(j)%bii
                                qatmtype_tm(nq)%ci     = top(istate)%atomtype(j)%ci
                                qatmtype_tm(nq)%ai     = top(istate)%atomtype(j)%ai
                        end if
                end do
        end do
        allocate(qtop(istate)%atomtype(nq))
        qtop(istate)%atomtype(:) = qatmtype_tm(1:nq)
        qtop(istate)%atomtotal   = nq
        deallocate(qatmtype_tm,top(istate)%atomtype)

        end do !end state

!qbond types-----------------------------------------------------------------------------------------------------------------

        do istate=1,nstate !over states
        qtop(istate)%bondtype_flag = top(istate)%bondtype_flag
        if (top(istate)%bondtype_flag) then
        if(debug .ne. "p") print*,"finding q bonds"
        allocate(bondtype_tm(top(istate)%bondtotal))

        nq = 0
        do i=1,top(istate)%bondtotal !over bonds
                k  = 0 !for each bond
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%bondtype(i)%a1 - qstatemap(j,istate)) .lt. 0.01 ) then
                                bondtype_tm(nq+1)%a1 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%bondtype(i)%a2 - qstatemap(j,istate)) .lt. 0.01 ) then
                                bondtype_tm(nq+1)%a2 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                if ( abs(k - 2) .lt. 0.01) then
                        nq = nq+1
                        bondtype_tm(nq)%code  = top(istate)%bondtype(i)%code
                        bondtype_tm(nq)%indic = top(istate)%bondtype(i)%indic
                        bondtype_tm(nq)%force = top(istate)%bondtype(i)%force
                        bondtype_tm(nq)%dis   = top(istate)%bondtype(i)%dis
                end if
        end do
        allocate(qtop(istate)%bondtype(nq))
        qtop(istate)%bondtype(:) = bondtype_tm(1:nq)
        qtop(istate)%bondtotal   = nq
        deallocate(bondtype_tm,top(istate)%bondtype)

        end if
        end do !end states

!q angle types-----------------------------------------------------------------------------------------------------------------

        do istate=1,nstate !over states
        qtop(istate)%angltype_flag = top(istate)%angltype_flag
        if (top(istate)%angltype_flag) then
        if(debug .ne. "p") print*,"finding q angles"
        allocate(angltype_tm(top(istate)%angltotal))
        nq = 0
        do i=1,top(istate)%angltotal !over angle
                k  = 0 !for each angle
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%angltype(i)%a1 - qstatemap(j,istate)) .lt. 0.01 ) then
                                angltype_tm(nq+1)%a1 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%angltype(i)%a2 - qstatemap(j,istate)) .lt. 0.01 ) then
                                angltype_tm(nq+1)%a2 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%angltype(i)%a3 - qstatemap(j,istate)) .lt. 0.01 ) then
                                angltype_tm(nq+1)%a3 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                if ( abs(k - 3) .lt. 0.01) then
                        nq = nq+1
                        angltype_tm(nq)%code  = top(istate)%angltype(i)%code
                        angltype_tm(nq)%indic = top(istate)%angltype(i)%indic
                        angltype_tm(nq)%force = top(istate)%angltype(i)%force
                        angltype_tm(nq)%dis   = top(istate)%angltype(i)%dis
                end if
        end do
        allocate(qtop(istate)%angltype(nq))
        qtop(istate)%angltype(:) = angltype_tm(1:nq)
        qtop(istate)%angltotal   = nq
        deallocate(angltype_tm,top(istate)%angltype)

        end if
        end do !end states

!q torsion types-----------------------------------------------------------------------------------------------------------------

        do istate=1,nstate !over states
        qtop(istate)%torstype_flag = top(istate)%torstype_flag
        if (top(istate)%torstype_flag) then
        if(debug .ne. "p") print*,"finding q torsions"
        allocate(torstype_tm(top(istate)%torstotal))
        nq = 0
        do i=1,top(istate)%torstotal !over torsions
                k  = 0 !for each torsion
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%torstype(i)%a1 - qstatemap(j,istate)) .lt. 0.01 ) then
                                torstype_tm(nq+1)%a1 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%torstype(i)%a2 - qstatemap(j,istate)) .lt. 0.01 ) then
                                torstype_tm(nq+1)%a2 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%torstype(i)%a3 - qstatemap(j,istate)) .lt. 0.01 ) then
                                torstype_tm(nq+1)%a3 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%torstype(i)%a4 - qstatemap(j,istate)) .lt. 0.01 ) then
                                torstype_tm(nq+1)%a4 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                if ( abs(k - 4) .lt. 0.01) then
                        nq = nq+1
                        torstype_tm(nq)%code(1)   = top(istate)%torstype(i)%code(1)
                        torstype_tm(nq)%indic     = top(istate)%torstype(i)%indic
                        torstype_tm(nq)%force(1)  = top(istate)%torstype(i)%force(1)
                        torstype_tm(nq)%dis(1)    = top(istate)%torstype(i)%dis(1)
                        torstype_tm(nq)%multip(1) = top(istate)%torstype(i)%multip(1)
                end if 
        end do
        allocate(qtop(istate)%torstype(nq))
        qtop(istate)%torstype(:) = torstype_tm(1:nq)
        qtop(istate)%torstotal   = nq
        deallocate(torstype_tm,top(istate)%torstype)

        end if
        end do !end states

!q improper types-----------------------------------------------------------------------------------------------------------------

        do istate=1,nstate !over states
        qtop(istate)%imprtype_flag = top(istate)%imprtype_flag
        if (top(istate)%imprtype_flag) then
        if(debug .ne. "p") print*,"finding q impropers"
        allocate(imprtype_tm(top(istate)%imprtotal))
        nq = 0
        do i=1,top(istate)%imprtotal !over impropers
                k  = 0        !for each impropers
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%imprtype(i)%a1 - qstatemap(j,istate)) .lt. 0.01 ) then
                                imprtype_tm(nq+1)%a1 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%imprtype(i)%a2 - qstatemap(j,istate)) .lt. 0.01 ) then
                                imprtype_tm(nq+1)%a2 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%imprtype(i)%a3 - qstatemap(j,istate)) .lt. 0.01 ) then
                                imprtype_tm(nq+1)%a3 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                do j=1,size(qstatemap(:,istate)) !over qatoms
                        if ( abs(top(istate)%imprtype(i)%a4 - qstatemap(j,istate)) .lt. 0.01 ) then
                                imprtype_tm(nq+1)%a4 = qstatemap(j,1)
                                k=k+1
                        end if
                end do
                if ( abs(k - 4) .lt. 0.01) then
                        nq = nq+1
                        imprtype_tm(nq)%code   = top(istate)%imprtype(i)%code
                        imprtype_tm(nq)%indic  = top(istate)%imprtype(i)%indic
                        imprtype_tm(nq)%force  = top(istate)%imprtype(i)%force
                        imprtype_tm(nq)%dis    = top(istate)%imprtype(i)%dis
                end if 
        end do
        allocate(qtop(istate)%imprtype(nq))
        qtop(istate)%imprtype(:) = imprtype_tm(1:nq)
        qtop(istate)%imprtotal   = nq
        deallocate(imprtype_tm,top(istate)%imprtype)
        end if
        end do !end states

end subroutine q_construct
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             q_library                                                       !
!-----------------------------------------------------------------------------------------------------------------------------!
!find uniqe elements that are changing,get first element and itirate additivly to build the lib
subroutine q_library(qtop,q_lib)
 type(top_type),intent(inout)                                ::        qtop(:)
 type(q_lib_type),intent(inout)                              ::        q_lib
!local
 integer(4)                                                  ::        nstate,istate,i,j,h,length,stat,el_max
 type(q_lib_type)                                            ::        q_lib_tm,q_lib_tm2

!-----------------------------------------------------------------------------------------------------------------------------!

        nstate = size(qtop(:))
        q_lib%atm_flag  = .false.
        q_lib%bond_flag = .false.
        q_lib%angl_flag = .false.
        q_lib%tors_flag = .false.
        q_lib%impr_flag = .false.

!atom part of lib------------------------------------------------------------------------------------------------------------
!get max lentght of temporary matrix and update status flag 
        el_max = 1
        q_lib_tm%atm_flag = .false.
        do istate=1,nstate 
                el_max = (el_max + qtop(istate)%atomtotal) 
                if (qtop(istate)%atomtype_flag) q_lib_tm%atm_flag = .true.
        end do

!construct atoms library itiratively 
if (q_lib_tm%atm_flag) then
        if(debug .ne. "p") print*,"constructing  q atom library"
        allocate(q_lib_tm%atm(el_max,nstate))

!initialize the first element(just names)

        length = 1

!itirate to get uniqe names
        do istate=1,nstate !over state

        if (qtop(istate)%atomtype_flag) then !block undefined access
        do i=1,qtop(istate)%atomtotal !over qatoms
                stat = 0
                do j=1,length
                        if ( qtop(istate)%atomtype(i)%atmno .eq. q_lib_tm%atm(j,1)%atmno ) then
                        stat = stat + 1
                        end if
                end do
                if ( abs(stat - 0) .lt. 0.01 ) then
                        length = length + 1
                        q_lib_tm%atm(length,:)%atmno = qtop(istate)%atomtype(i)%atmno
                end if
        end do
        end if !block undefined access
        end do
        allocate(q_lib%atm(length,nstate))
        q_lib%atm(:,:) = q_lib_tm%atm(1:length,1:nstate)
        q_lib%atm_flag = q_lib_tm%atm_flag
!get parameter for lib at different states
!no further comparision for atoms even thogh they might not change good to have them all 

        do i=1,length !over lib
        do istate=1,nstate !over state
        if (qtop(istate)%atomtype_flag) then !block undefined access
                do j=1,qtop(istate)%atomtotal
                        if( q_lib%atm(i,istate)%atmno .eq. qtop(istate)%atomtype(j)%atmno) then
                                q_lib%atm(i,istate) = qtop(istate)%atomtype(j)
                        end if
                end do
        end if !block undefined access
        end do
        end do

deallocate(q_lib_tm%atm)
end if

!bond part of lib------------------------------------------------------------------------------------------------------------
!get max lentght of temporary matrix and update status flag 

        el_max = 1
        q_lib_tm%bond_flag = .false.
        do istate=1,nstate 
                el_max = (el_max + qtop(istate)%bondtotal)
                if (qtop(istate)%bondtype_flag) q_lib_tm%bond_flag = .true.
        end do

!construct atoms library itiratively 
if (q_lib_tm%bond_flag) then
        if(debug .ne. "p") print*,"constructing  q bond library"
        allocate(q_lib_tm%bond(el_max,nstate))

!initialize the first element(just names)

        length = 1

!itirate to get uniqe names
        do istate=1,nstate !over state

        if (qtop(istate)%bondtype_flag) then !block undefined access
        do i=1,qtop(istate)%bondtotal !over qatoms
                stat = 0
                do j=1,length
                        if ( qtop(istate)%bondtype(i) .same. q_lib_tm%bond(j,1) ) then
                        stat = stat + 1
                        end if
                end do
                if ( abs(stat - 0) .lt. 0.01 ) then
                        length = length + 1
                        q_lib_tm%bond(length,:)%a1 = qtop(istate)%bondtype(i)%a1
                        q_lib_tm%bond(length,:)%a2 = qtop(istate)%bondtype(i)%a2
                end if
        end do
        end if !block undefined access
        end do
        allocate(q_lib_tm2%bond(length,nstate))
        q_lib_tm2%bond(:,:) = q_lib_tm%bond(1:length,1:nstate)
        q_lib%bond_flag = q_lib_tm%bond_flag

!get parameter for lib at different states
        do i=1,length !over lib
        do istate=1,nstate !over state
        if (qtop(istate)%bondtype_flag) then !block undefined access
                do j=1,qtop(istate)%bondtotal
                        if( q_lib_tm2%bond(i,istate) .same. qtop(istate)%bondtype(j)) then
                                q_lib_tm2%bond(i,istate) = qtop(istate)%bondtype(j)
                        end if
                end do
        end if !block undefined access
        end do
        end do

if(debug .eq. "p") then
print*,"all-Q-bond"
write(nstate_char,'(i4)') nstate
write(qformat,*) "(i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
do i=1,length !over lib
write(*,qformat),i,q_lib_tm2%bond(i,1)%a1,q_lib_tm2%bond(i,1)%a2,(q_lib_tm2%bond(i,h)%force,h=1,nstate) &
                                                                                ,(q_lib_tm2%bond(i,h)%dis,h=1,nstate)
end do
write(*,*),"--------------------------------------------------------------------------"
end if

!remove the terms that dont change reletive to state one , first count then remove
        length = 0
        do i=1,size(q_lib_tm2%bond(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%bond(i,1) .eqt. q_lib_tm2%bond(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                end if
        end do

        allocate(q_lib%bond(length,nstate))


        length = 0
        do i=1,size(q_lib_tm2%bond(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%bond(i,1) .eqt. q_lib_tm2%bond(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                        q_lib%bond(length,:) = q_lib_tm2%bond(i,:)
                end if
        end do

deallocate(q_lib_tm2%bond,q_lib_tm%bond)
end if

!angle part of lib------------------------------------------------------------------------------------------------------------
!get max lentght of temporary matrix and update status flag 

        el_max = 1
        q_lib_tm%angl_flag = .false.
        do istate=1,nstate 
                el_max = (el_max + qtop(istate)%angltotal)
                if (qtop(istate)%angltype_flag) q_lib_tm%angl_flag = .true.
        end do

!construct atoms library itiratively 
if (q_lib_tm%angl_flag) then
        if(debug .ne. "p") print*,"constructing  q angle library"
        allocate(q_lib_tm%angl(el_max,nstate))

!initialize the first element(just names)
        length = 1

!itirate to get uniqe names
        do istate=1,nstate !over state

        if (qtop(istate)%angltype_flag) then !block undefined access
        do i=1,qtop(istate)%angltotal !over qatoms
                stat = 0
                do j=1,length
                        if ( qtop(istate)%angltype(i) .same. q_lib_tm%angl(j,1) ) then
                        stat = stat + 1
                        end if
                end do
                if ( abs(stat - 0) .lt. 0.01 ) then
                        length = length + 1
                        q_lib_tm%angl(length,:)%a1 = qtop(istate)%angltype(i)%a1
                        q_lib_tm%angl(length,:)%a2 = qtop(istate)%angltype(i)%a2
                        q_lib_tm%angl(length,:)%a3 = qtop(istate)%angltype(i)%a3
                end if
        end do
        end if !block undefined access
        end do
        allocate(q_lib_tm2%angl(length,nstate))
        q_lib_tm2%angl(:,:) = q_lib_tm%angl(1:length,1:nstate)
        q_lib%angl_flag = q_lib_tm%angl_flag

!get parameter for lib at different states
        do i=1,length !over lib
        do istate=1,nstate !over state
        if (qtop(istate)%angltype_flag) then !block undefined access
                do j=1,qtop(istate)%angltotal
                        if( q_lib_tm2%angl(i,istate) .same. qtop(istate)%angltype(j)) then
                                q_lib_tm2%angl(i,istate) = qtop(istate)%angltype(j)
                        end if
                end do
        end if !block undefined access
        end do
        end do

if(debug .eq. "p") then
print*,"all-Q-angles"
write(nstate_char,'(i4)') nstate
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
do i=1,length !over lib
write(*,qformat),i,q_lib_tm2%angl(i,1)%a1,q_lib_tm2%angl(i,1)%a2,q_lib_tm2%angl(i,1)%a3,    &
                                        (q_lib_tm2%angl(i,h)%force,h=1,nstate),(q_lib_tm2%angl(i,h)%dis,h=1,nstate)
end do
write(*,*),"--------------------------------------------------------------------------"
end if

!remove the terms that dont change reletive to state one , first count then remove
        length = 0
        do i=1,size(q_lib_tm2%angl(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%angl(i,1) .eqt. q_lib_tm2%angl(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                end if
        end do

        allocate(q_lib%angl(length,nstate))


        length = 0
        do i=1,size(q_lib_tm2%angl(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%angl(i,1) .eqt. q_lib_tm2%angl(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                        q_lib%angl(length,:) = q_lib_tm2%angl(i,:)
                end if
        end do

deallocate(q_lib_tm2%angl,q_lib_tm%angl)
end if

!torsion part of lib------------------------------------------------------------------------------------------------------------
!get max lentght of temporary matrix and update status flag 

        el_max = 1
        q_lib_tm%tors_flag = .false.
        do istate=1,nstate 
                el_max = (el_max + qtop(istate)%torstotal)
                if (qtop(istate)%torstype_flag) q_lib_tm%tors_flag = .true.
        end do

!construct atoms library itiratively 
if (q_lib_tm%tors_flag) then
        if(debug .ne. "p") print*,"constructing  q torsion library"
        allocate(q_lib_tm%tors(el_max,nstate))

!initialize the first element(just names)
        length = 1

!itirate to get uniqe names
        do istate=1,nstate !over state

        if (qtop(istate)%torstype_flag) then !block undefined access
        do i=1,qtop(istate)%torstotal !over qatoms
                stat = 0
                do j=1,length
                        if ( qtop(istate)%torstype(i) .same. q_lib_tm%tors(j,1) ) then
                        stat = stat + 1
                        end if
                end do
                if ( abs(stat - 0) .lt. 0.01 ) then
                        length = length + 1
                        q_lib_tm%tors(length,:)%a1 = qtop(istate)%torstype(i)%a1
                        q_lib_tm%tors(length,:)%a2 = qtop(istate)%torstype(i)%a2
                        q_lib_tm%tors(length,:)%a3 = qtop(istate)%torstype(i)%a3
                        q_lib_tm%tors(length,:)%a4 = qtop(istate)%torstype(i)%a4
                end if
        end do
        end if !block undefined access
        end do

        allocate(q_lib_tm2%tors(length,nstate))
        q_lib_tm2%tors(:,:) = q_lib_tm%tors(1:length,1:nstate)
        q_lib%tors_flag = q_lib_tm%tors_flag

!get parameter for lib at different states
!some extra work for tors.
! count how many times each uniqe name is printed in each state = nterms
! save tors parm according to their nterm
        do i=1,length !over lib
        do istate=1,nstate !over state
        if (qtop(istate)%torstype_flag) then !block undefined access
                stat = 0
                do j=1,qtop(istate)%torstotal
                        if( q_lib_tm2%tors(i,istate) .same. qtop(istate)%torstype(j)) then
                                stat = stat + 1
                                q_lib_tm2%tors(i,istate)%nterm = stat
                                if ( abs(stat - 1) .lt. 0.01) then
                                q_lib_tm2%tors(i,istate)%multip(1) = qtop(istate)%torstype(j)%multip(1)
                                q_lib_tm2%tors(i,istate)%force(1)  = qtop(istate)%torstype(j)%force(1)
                                q_lib_tm2%tors(i,istate)%dis(1)    = qtop(istate)%torstype(j)%dis(1)
                                elseif ( abs(stat - 2) .lt. 0.01) then
                                q_lib_tm2%tors(i,istate)%multip(2) = qtop(istate)%torstype(j)%multip(1)
                                q_lib_tm2%tors(i,istate)%force(2)  = qtop(istate)%torstype(j)%force(1)
                                q_lib_tm2%tors(i,istate)%dis(2)    = qtop(istate)%torstype(j)%dis(1)
                                elseif ( abs(stat - 3) .lt. 0.01) then
                                q_lib_tm2%tors(i,istate)%multip(3) = qtop(istate)%torstype(j)%multip(1)
                                q_lib_tm2%tors(i,istate)%force(3)  = qtop(istate)%torstype(j)%force(1)
                                q_lib_tm2%tors(i,istate)%dis(3)    = qtop(istate)%torstype(j)%dis(1)
                                end if
                        end if
                end do
        end if !block undefined access
        end do
        end do

if(debug .eq. "p") then
print*,"all-Q-torsion"
write(nstate_char,'(i4)') nstate
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)), &
                                                                        "(f8.3,X),",trim(adjustl(nstate_char)),"(f8.3,X))"
do i=1,length !over lib
write(*,qformat),i,q_lib_tm2%tors(i,1)%a1,q_lib_tm2%tors(i,1)%a2,q_lib_tm2%tors(i,1)%a3,q_lib_tm2%tors(i,1)%a4, &
(q_lib_tm2%tors(i,h)%multip(1),h=1,nstate),(q_lib_tm2%tors(i,h)%force(1),h=1,nstate),(q_lib_tm2%tors(i,h)%dis(1),h=1,nstate)
write(*,qformat),i,q_lib_tm2%tors(i,1)%a1,q_lib_tm2%tors(i,1)%a2,q_lib_tm2%tors(i,1)%a3,q_lib_tm2%tors(i,1)%a4, &
(q_lib_tm2%tors(i,h)%multip(2),h=1,nstate),(q_lib_tm2%tors(i,h)%force(2),h=1,nstate),(q_lib_tm2%tors(i,h)%dis(2),h=1,nstate)
write(*,qformat),i,q_lib_tm2%tors(i,1)%a1,q_lib_tm2%tors(i,1)%a2,q_lib_tm2%tors(i,1)%a3,q_lib_tm2%tors(i,1)%a4, &
(q_lib_tm2%tors(i,h)%multip(3),h=1,nstate),(q_lib_tm2%tors(i,h)%force(3),h=1,nstate),(q_lib_tm2%tors(i,h)%dis(3),h=1,nstate)
write(*,*),"#"
end do
write(*,*),"--------------------------------------------------------------------------"
end if

!remove the terms that dont change reletive to state one , first count then remove
!the equivalency of two torsions is a function of all 9 variables
        length = 0
        do i=1,size(q_lib_tm2%tors(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%tors(i,1) .eqt. q_lib_tm2%tors(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                end if
        end do

        allocate(q_lib%tors(length,nstate))

        length = 0
        do i=1,size(q_lib_tm2%tors(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%tors(i,1) .eqt. q_lib_tm2%tors(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                        q_lib%tors(length,:) = q_lib_tm2%tors(i,:)
                end if
        end do

deallocate(q_lib_tm2%tors,q_lib_tm%tors)
end if

!improper part of lib------------------------------------------------------------------------------------------------------------
!get max lentght of temporary matrix and update status flag 

        el_max = 1
        q_lib_tm%impr_flag = .false.
        do istate=1,nstate 
                el_max = (el_max + qtop(istate)%imprtotal)
                if (qtop(istate)%imprtype_flag) q_lib_tm%impr_flag = .true.
        end do

!construct atoms library itiratively 
if (q_lib_tm%impr_flag) then
        if(debug .ne. "p") print*,"constructing  q improper library"
        allocate(q_lib_tm%impr(el_max,nstate))

!initialize the first element(just names)
        length = 1

!itirate to get uniqe names
        do istate=1,nstate !over state

        if (qtop(istate)%imprtype_flag) then !block undefined access
        do i=1,qtop(istate)%imprtotal !over qatoms
                stat = 0
                do j=1,length
                        if ( qtop(istate)%imprtype(i) .same. q_lib_tm%impr(j,1) ) then
                        stat = stat + 1
                        end if
                end do
                if ( abs(stat - 0) .lt. 0.01 ) then
                        length = length + 1
                        q_lib_tm%impr(length,:)%a1 = qtop(istate)%imprtype(i)%a1
                        q_lib_tm%impr(length,:)%a2 = qtop(istate)%imprtype(i)%a2
                        q_lib_tm%impr(length,:)%a3 = qtop(istate)%imprtype(i)%a3
                        q_lib_tm%impr(length,:)%a4 = qtop(istate)%imprtype(i)%a4
                end if
        end do
        end if !block undefined access
        end do
        allocate(q_lib_tm2%impr(length,nstate))

        q_lib_tm2%impr(:,:) = q_lib_tm%impr(1:length,1:nstate)
        q_lib%impr_flag = q_lib_tm%impr_flag

!get parameter for lib at different states
        do i=1,length !over lib
        do istate=1,nstate !over state
        if (qtop(istate)%imprtype_flag) then !block undefined access
                do j=1,qtop(istate)%imprtotal
                        if( q_lib_tm2%impr(i,istate) .same. qtop(istate)%imprtype(j)) then
                                q_lib_tm2%impr(i,istate) = qtop(istate)%imprtype(j)
                        end if
                end do
        end if !block undefined access
        end do
        end do


if(debug .eq. "p") then
print*,"all-Q-improper"
write(nstate_char,'(i4)') nstate
        write(qformat,*) "(i4,X,i4,X,i4,X,i4,X,i4,X,",trim(adjustl(nstate_char)),"(f8.3,X),",trim(adjustl(nstate_char)), &
                                                                                                                "(f8.3,X))"

do i=1,length !over lib
write(*,qformat),i,q_lib_tm2%impr(i,1)%a1,q_lib_tm2%impr(i,1)%a2,q_lib_tm2%impr(i,1)%a3,q_lib_tm2%impr(i,1)%a4 &
                                            ,(q_lib_tm2%impr(i,h)%force,h=1,nstate),(q_lib_tm2%impr(i,h)%dis,h=1,nstate)
end do
write(*,*),"--------------------------------------------------------------------------"
end if

!remove the terms that dont change reletive to state one , first count then remove
        length = 0
        do i=1,size(q_lib_tm2%impr(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%impr(i,1) .eqt. q_lib_tm2%impr(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                end if
        end do

        allocate(q_lib%impr(length,nstate))

        length = 0
        do i=1,size(q_lib_tm2%impr(:,1)) !over qatoms
                stat = 0
                do istate=2,nstate
                        if (q_lib_tm2%impr(i,1) .eqt. q_lib_tm2%impr(i,istate) ) then
                        stat = stat + 1
                        end if
                end do
                if ( stat .lt. (nstate-1) ) then !change in at least one state
                        length = length+1
                        q_lib%impr(length,:) = q_lib_tm2%impr(i,:)
                end if
        end do

deallocate(q_lib_tm2%impr,q_lib_tm%impr)
end if


end subroutine q_library
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             p_library                                                       !
!-----------------------------------------------------------------------------------------------------------------------------!
!find uniqe elements q lib and make a parameter lib same strategies
subroutine p_library(q_lib,p_lib)
 type(q_lib_type),intent(in out)                            ::        q_lib
 type(q_lib_type),intent(out)                               ::        p_lib
!local
 integer(4)                                                 ::        nstate,istate,i,j,k,length,stat,el_max 
 type(q_lib_type)                                           ::        p_lib_tm 
!-----------------------------------------------------------------------------------------------------------------------------!

!make atom parmameter library
        if(q_lib%atm_flag) then
        if(debug .ne. "p") print*,"constructing  p atom library"
        nstate = size(q_lib%atm(1,:))
        el_max = (nstate*size(q_lib%atm(:,1))) + 1
        allocate(p_lib_tm%atm(el_max,1))

!find all uniqe parameters
        length = 1

        do i=1,size(q_lib%atm(:,1))        !over q_lib
                do istate=1,nstate        !over state
                        stat = 0
                        do j=1,length        !over p_lib
                                if ( q_lib%atm(i,istate) .eqt. p_lib_tm%atm(j,1) ) then
                                        stat = stat + 1
                                end if
                        end do
                        if ( abs(stat - 0 ) .lt. 0.01 ) then
                                length = length + 1 
                                p_lib_tm%atm(length,1) = q_lib%atm(i,istate)
                        end if
                end do
        end do

        allocate(p_lib%atm(length,1))
        p_lib%atm(:,1) = p_lib_tm%atm(1:length,1)

!assign the codes q_lib by p_lib
        do i=1,size(q_lib%atm(:,1)) ! over q_lib
                do istate=1,nstate
                        do j=1,size(p_lib%atm(:,1))
                                if ( q_lib%atm(i,istate) .eqt. p_lib%atm(j,1) ) then
                                        q_lib%atm(i,istate)%code = j
                                end if
                        end do
                end do
        end do
        deallocate(p_lib_tm%atm)
        end if

!make bond parmameter library
        if(q_lib%bond_flag) then
        if(debug .ne. "p") print*,"constructing  p bond library"
        nstate = size(q_lib%bond(1,:))
        el_max = (nstate*size(q_lib%bond(:,1))) + 1
        allocate(p_lib_tm%bond(el_max,1))

!find all uniqe parameters
        length = 1

        do i=1,size(q_lib%bond(:,1))        !over q_lib
                do istate=1,nstate        !over state
                        stat = 0
                        do j=1,length        !over p_lib
                                if ( q_lib%bond(i,istate) .eqt. p_lib_tm%bond(j,1) ) then
                                        stat = stat + 1
                                end if
                        end do
                        if ( abs(stat - 0 ) .lt. 0.01 ) then
                                length = length + 1 
                                p_lib_tm%bond(length,1) = q_lib%bond(i,istate)
                        end if
                end do
        end do

        allocate(p_lib%bond(length,1))
        p_lib%bond(:,1) = p_lib_tm%bond(1:length,1)

!assign the codes q_lib by p_lib
        do i=1,size(q_lib%bond(:,1)) ! over q_lib
                do istate=1,nstate
                        do j=1,size(p_lib%bond(:,1))
                                if ( q_lib%bond(i,istate) .eqt. p_lib%bond(j,1) ) then
                                        q_lib%bond(i,istate)%code = j
                                end if
                        end do
                end do
        end do
        deallocate(p_lib_tm%bond)
        end if

!make angle parmameter library
        if(q_lib%angl_flag) then
        if(debug .ne. "p") print*,"constructing  p angle library"
        nstate = size(q_lib%angl(1,:))
        el_max = (nstate*size(q_lib%angl(:,1))) + 1
        allocate(p_lib_tm%angl(el_max,1))

!find all uniqe parameters
        length = 1

        do i=1,size(q_lib%angl(:,1))        !over q_lib
                do istate=1,nstate        !over state
                        stat = 0
                        do j=1,length        !over p_lib
                                if ( q_lib%angl(i,istate) .eqt. p_lib_tm%angl(j,1) ) then
                                        stat = stat + 1
                                end if
                        end do
                        if ( abs(stat - 0 ) .lt. 0.01 ) then
                                length = length + 1 
                                p_lib_tm%angl(length,1) = q_lib%angl(i,istate)
                        end if
                end do
        end do

        allocate(p_lib%angl(length,1))
        p_lib%angl(:,1) = p_lib_tm%angl(1:length,1)

!assign the codes q_lib by p_lib
        do i=1,size(q_lib%angl(:,1)) ! over q_lib
                do istate=1,nstate
                        do j=1,size(p_lib%angl(:,1))
                                if ( q_lib%angl(i,istate) .eqt. p_lib%angl(j,1) ) then
                                        q_lib%angl(i,istate)%code = j
                                end if
                        end do
                end do
        end do
        deallocate(p_lib_tm%angl)
        end if

!****************************************************
!make torsion parmameter library-----
        if(q_lib%tors_flag) then
        if(debug .ne. "p") print*,"constructing  p torsion library"
        nstate = size(q_lib%tors(1,:))
        el_max = (3*nstate*size(q_lib%tors(:,1)))+1
        allocate(p_lib_tm%tors(el_max,1))
        if(debug .ne. "p") print*,"        finding uniqe torsions"
!find all uniqe parameters
        length = 1
        do i=1,size(q_lib%tors(:,1))        !over q_lib
                do istate=1,nstate        !over state
                do k=1,3                !over 3 components
                        stat = 0
                        do j=1,length        !over p_lib
                        if ( (abs(q_lib%tors(i,istate)%multip(k) - p_lib_tm%tors(j,1)%multip(1)) .lt. 0.01)  .and. &
                             (abs(q_lib%tors(i,istate)%force(k)  - p_lib_tm%tors(j,1)%force(1) ) .lt. 0.01)  .and. &
                             (abs(q_lib%tors(i,istate)%dis(k)    - p_lib_tm%tors(j,1)%dis(1)   ) .lt. 0.01) ) then
                                stat = stat + 1
                        end if
                        end do
                        if ( abs(stat - 0 ) .lt. 0.01 ) then
                                length = length + 1 
                                p_lib_tm%tors(length,1)%multip(1) = q_lib%tors(i,istate)%multip(k)
                                p_lib_tm%tors(length,1)%force(1)  = q_lib%tors(i,istate)%force(k)
                                p_lib_tm%tors(length,1)%dis(1)    = q_lib%tors(i,istate)%dis(k)
                        end if
                end do
                end do
        end do
        allocate( p_lib%tors(length,1) )
        p_lib%tors(:,1) = p_lib_tm%tors(1:length,1)
        if(debug .ne. "p") print*,"        assigning parameter to q torsions library "
!assign the codes q_lib by p_lib
        do i=1,size(q_lib%tors(:,1))        ! over q_lib
                do istate=1,nstate        !over state
                do k=1,3                !over 3 components
                        do j=1,size(p_lib%tors(:,1))
                        if ( (abs(q_lib%tors(i,istate)%multip(k) - p_lib_tm%tors(j,1)%multip(1)) .lt. 0.01)  .and. &
                             (abs(q_lib%tors(i,istate)%force(k)  - p_lib_tm%tors(j,1)%force(1))  .lt. 0.01)  .and. &
                             (abs(q_lib%tors(i,istate)%dis(k)    - p_lib_tm%tors(j,1)%dis(1))    .lt. 0.01) ) then

                                        q_lib%tors(i,istate)%code(k) = j
                                end if
                        end do
                end do
                end do
        end do
        deallocate(p_lib_tm%tors)
        end if
!*********************************************

!make improper parmameter library
        if(q_lib%impr_flag) then
        if(debug .ne. "p") print*,"constructing  p improper library"
        nstate = size(q_lib%impr(1,:))
        el_max = (nstate*size(q_lib%impr(:,1))) + 1
        allocate(p_lib_tm%impr(el_max,1))

!find all uniqe parameters
        length = 1
        do i=1,size(q_lib%impr(:,1))        !over q_lib
                do istate=1,nstate        !over state
                        stat = 0
                        do j=1,length        !over p_lib
                        if ( q_lib%impr(i,istate) .eqt. p_lib_tm%impr(j,1) ) then
                                        stat = stat + 1
                        end if
                        end do
                        if ( abs(stat - 0 ) .lt. 0.01 ) then
                        length = length + 1 
                        p_lib_tm%impr(length,1) = q_lib%impr(i,istate)
                        end if
                end do
        end do

        allocate(p_lib%impr(length,1))
        p_lib%impr(:,1) = p_lib_tm%impr(1:length,1)

!assign the codes q_lib by p_lib
        do i=1,size(q_lib%impr(:,1)) ! over q_lib
                do istate=1,nstate
                        do j=1,size(p_lib%impr(:,1))
                        if ( q_lib%impr(i,istate) .eqt. p_lib%impr(j,1) ) then
                                q_lib%impr(i,istate)%code = j
                        end if
                        end do
                end do
        end do
        deallocate(p_lib_tm%impr)
        end if

end subroutine p_library
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                               c_library                                                     !
!-----------------------------------------------------------------------------------------------------------------------------!
subroutine c_library(q_lib,p_lib,cbond,a_kapl,t_kapl,i_kapl)
!find uniqe elements q lib and make a parameter lib same strategies
 type(bondtype_type),allocatable,intent(inout)                 ::        cbond(:,:)
 type(q_lib_type),intent(inout)                                ::        q_lib,p_lib
 integer(4),allocatable,intent(inout)                          ::        a_kapl(:,:),t_kapl(:,:),i_kapl(:,:)
!local
 integer(4)                                                    ::        nstate,istate,i,j,length,stat
 integer(4),allocatable                                        ::        a_kapl_tm(:,:),t_kapl_tm(:,:),i_kapl_tm(:,:)
 type(bondtype_type),allocatable                               ::        cbond_tm(:,:)
 type(bondtype_type)                                           ::        dummy_bond
!-----------------------------------------------------------------------------------------------------------------------------!
!get the list of breaking or forming bonds
if(q_lib%bond_flag) then

        nstate = size(q_lib%bond(1,:))
        allocate(cbond_tm(size(q_lib%bond(:,1)),nstate))
        cbond_tm(:,:)%code = 1  !to keep track of softcores statewise(1 = yes, 0 = no)
        length = 0
        do i=1,size(q_lib%bond(:,1))
                stat = 0
                do istate=1,nstate
                        if (q_lib%bond(i,istate) .eqt. dummy_bond) then
                                stat = stat + 1
                                cbond_tm(length+1,istate)%code = 0
                        end if
                end do
                if ( stat .gt. 0) then
                        length = length + 1 
                        cbond_tm(length,1:nstate)%a1 = q_lib%bond(i,1:nstate)%a1
                        cbond_tm(length,1:nstate)%a2 = q_lib%bond(i,1:nstate)%a2
                        cbond_tm(length,1:nstate)%indic = i !keep track of array index for couplings
                end if
        end do
        allocate(cbond(length,nstate))
        cbond(:,:) = cbond_tm(1:length,1:nstate)
        deallocate (cbond_tm)

if (size(cbond) .gt. 0) then

!update the p_lib for soft_pairs variables
        do i=1,size(cbond(:,1))
                do istate=1,nstate
                        do j=1,size(q_lib%atm(:,1))
                                if ( abs(cbond(i,1)%a1 - q_lib%atm(j,1)%atmno) .lt. 0.01) then
                                        p_lib%atm((q_lib%atm(j,istate)%code),1)%ai = 1.58
                                end if
                                if ( abs(cbond(i,1)%a2 - q_lib%atm(j,1)%atmno) .lt. 0.01) then
                                        p_lib%atm((q_lib%atm(j,istate)%code),1)%ai = 1.58
                                end if
                        end do
                end do
        end do


!if an angle, torsion or improper contains any successive combination of bond atoms it would be suggested for couppling
!the implementation is done in this way to be in same sprit of common practice. in reallity its not needed in the first 
!place since all terms that include the breaking bond would be eliminated in coresbonding state according to parameter 
!file instruction

        if(q_lib%angl_flag) then
!find angles coupling.
        allocate(a_kapl_tm( size(q_lib%angl(:,1))*nstate , size(cbond(:,1)) + 1 ))
        a_kapl_tm(:,:) = 0
        length = 0                 !row count
        do i=1,size(q_lib%angl(:,1))
                stat = 1         !colum count
                do j=1,size(cbond(:,1))
                        if (q_lib%angl(i,1) .incl. cbond(j,1) ) then
                                stat   = stat + 1
                                a_kapl_tm(length+1,1) = i
                                a_kapl_tm(length+1,stat) = cbond(j,1)%indic
                        end if
                end do
                if ( stat .gt. 1) then
                                length = length + 1
                end if
        end do

        allocate(a_kapl( length , size(a_kapl_tm(1,:))))
        a_kapl(:,:) = a_kapl_tm(1:length,:)
        end if


        if(q_lib%tors_flag) then
!find torsion coupling.
        allocate(t_kapl_tm( size(q_lib%tors(:,1))*nstate , size(cbond(:,1)) + 1 ))
        t_kapl_tm(:,:) = 0
        length = 0                 !row count
        do i=1,size(q_lib%tors(:,1))
                stat = 1         !colum count
                do j=1,size(cbond(:,1))
                        if (q_lib%tors(i,1) .incl. cbond(j,1) ) then
                                stat   = stat + 1
                                t_kapl_tm(length+1,1) = i
                                t_kapl_tm(length+1,stat) = cbond(j,1)%indic
                        end if
                end do
                if ( stat .gt. 1) then
                                length = length + 1
                end if
        end do

        allocate(t_kapl( length , size(t_kapl_tm(1,:))))
        t_kapl(:,:) = t_kapl_tm(1:length,:)
        end if

        if(q_lib%impr_flag) then
!find the 0 1 improper
        allocate(i_kapl_tm( size(q_lib%impr(:,1))*nstate , size(cbond(:,1)) + 1 ))
        i_kapl_tm(:,:) = 0
        length = 0                 !row count
        do i=1,size(q_lib%impr(:,1))
                stat = 1           !colum count
                do j=1,size(cbond(:,1))
                        if (q_lib%impr(i,1) .incl. cbond(j,1) ) then
                                stat   = stat + 1
                                i_kapl_tm(length+1,1) = i
                                i_kapl_tm(length+1,stat) = cbond(j,1)%indic
                        end if
                end do
                if ( stat .gt. 1) then
                                length = length + 1
                end if
        end do

        allocate(i_kapl( length , size(i_kapl_tm(1,:))))
        i_kapl(:,:) = i_kapl_tm(1:length,:)
        end if



end if
end if
end subroutine c_library
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             angl_incl_bond                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function angl_incl_bond(a,b)

        type(angltype_type),intent(in)                :: a
        type(bondtype_type),intent(in)                :: b

        if ( ( ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a1 - b%a2) .lt. 0.01) .and. ( abs(a%a2 - b%a1) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a1) .lt. 0.01) .and. ( abs(a%a3 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a2) .lt. 0.01) .and. ( abs(a%a3 - b%a1) .lt. 0.01) ) ) then
                 angl_incl_bond = .true.
        else
                 angl_incl_bond = .false.
        end if
end function angl_incl_bond
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             tors_incl_bond                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function tors_incl_bond(a,b)

        type(torstype_type),intent(in)                :: a
        type(bondtype_type),intent(in)                :: b

        if ( ( ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a1 - b%a2) .lt. 0.01) .and. ( abs(a%a2 - b%a1) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a1) .lt. 0.01) .and. ( abs(a%a3 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a2) .lt. 0.01) .and. ( abs(a%a3 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a3 - b%a1) .lt. 0.01) .and. ( abs(a%a4 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a3 - b%a2) .lt. 0.01) .and. ( abs(a%a4 - b%a1) .lt. 0.01) ) ) then
                 tors_incl_bond = .true.
        else
                 tors_incl_bond = .false.
        end if
end function tors_incl_bond
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             impr_incl_bond                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function impr_incl_bond(a,b)

        type(imprtype_type),intent(in)                :: a
        type(bondtype_type),intent(in)                :: b

        if ( ( ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a1 - b%a2) .lt. 0.01) .and. ( abs(a%a2 - b%a1) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a1) .lt. 0.01) .and. ( abs(a%a3 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a2 - b%a2) .lt. 0.01) .and. ( abs(a%a3 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a3 - b%a1) .lt. 0.01) .and. ( abs(a%a4 - b%a2) .lt. 0.01) )   .or. &
             ( ( abs(a%a3 - b%a2) .lt. 0.01) .and. ( abs(a%a4 - b%a1) .lt. 0.01) ) ) then
                 impr_incl_bond = .true.
        else
                 impr_incl_bond = .false.
        end if
end function impr_incl_bond
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             atom_eqt_atom                                                   !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function atom_eqt_atom(a,b)

        type(atomtype_type),intent(in)                :: a,b

        if (( a%nam .eq. b%nam) .and. ( abs(a%mass - b%mass) .lt. 0.001 ) .and. ( abs(a%aii - b%aii) .lt. 0.001 ) .and. & 
        ( abs(a%bii - b%bii) .lt. 0.001 ) .and. ( abs(a%ci - b%ci) .lt. 0.001 ) .and. ( abs(a%ai - b%ai) .lt. 0.001 ) ) then

                 atom_eqt_atom = .true.
        else
                 atom_eqt_atom = .false.
        end if
end function atom_eqt_atom


!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             bond_same_bond                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function bond_same_bond(a,b)

        type(bondtype_type),intent(in)                :: a,b

        if ((  ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01)  ) .or. &
            (  ( abs(a%a1 - b%a2) .lt. 0.01) .and. ( abs(a%a2 - b%a1) .lt. 0.01)  )) then
                 bond_same_bond = .true.
        else
                 bond_same_bond = .false.
        end if
end function bond_same_bond
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             bond_eqt_bond                                                   !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function bond_eqt_bond(a,b)

        type(bondtype_type),intent(in)                :: a,b

        if ( ( abs(a%force - b%force) .lt. 0.001 ) .and. ( abs(a%dis - b%dis) .lt. 0.001 ) ) then
                 bond_eqt_bond = .true.
        else
                 bond_eqt_bond = .false.
        end if
end function bond_eqt_bond

!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             angl_same_angl                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function angl_same_angl(a,b)

        type(angltype_type),intent(in)                :: a,b

        if((  ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01) .and. ( abs(a%a3 - b%a3) .lt. 0.01)  ) .or. &
           (  ( abs(a%a1 - b%a3) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01) .and. ( abs(a%a3 - b%a1) .lt. 0.01)   ) ) then
                 angl_same_angl = .true.
        else
                 angl_same_angl = .false.
        end if
end function angl_same_angl
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             angl_eqt_angl                                                   !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function angl_eqt_angl(a,b)

        type(angltype_type),intent(in)                :: a,b

        if ( ( abs(a%force - b%force) .lt. 0.001 ) .and. ( abs(a%dis - b%dis) .lt. 0.001 ) ) then
                 angl_eqt_angl = .true.
        else
                 angl_eqt_angl = .false.
        end if
end function angl_eqt_angl
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             tors_same_tors                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function tors_same_tors(a,b)

        type(torstype_type),intent(in)                :: a,b

        if  ((( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01)     .and. &
              ( abs(a%a3 - b%a3) .lt. 0.01) .and. ( abs(a%a4 - b%a4) .lt. 0.01) )   .or.  &
             (( abs(a%a1 - b%a4) .lt. 0.01) .and. ( abs(a%a2 - b%a3) .lt. 0.01)    .and. &
              ( abs(a%a3 - b%a2) .lt. 0.01) .and. ( abs(a%a4 - b%a1) .lt. 0.01) ) )  then
                 tors_same_tors = .true.
        else
                 tors_same_tors = .false.
        end if
end function tors_same_tors
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             tors_eqt_tors                                                   !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function tors_eqt_tors(a,b)

        type(torstype_type),intent(in)                :: a,b

        if ( ( abs(a%force(1)  - b%force(1) ) .lt. 0.001 ) .and. ( abs(a%dis(1)    - b%dis(1)   ) .lt. 0.001 )   .and. &
             ( abs(a%multip(1) - b%multip(1)) .lt. 0.001 ) .and. ( abs(a%force(2)  - b%force(2) ) .lt. 0.001 )   .and. &
             ( abs(a%dis(2)    - b%dis(2)   ) .lt. 0.001 ) .and. ( abs(a%multip(2) - b%multip(2)) .lt. 0.001 )   .and. &
             ( abs(a%force(3)  - b%force(3) ) .lt. 0.001 ) .and. ( abs(a%dis(3)    - b%dis(3)   ) .lt. 0.001 )   .and. &
                                                                 ( abs(a%multip(3) - b%multip(3)) .lt. 0.001 )  ) then

                 tors_eqt_tors = .true.
        else
                 tors_eqt_tors = .false.
        end if
end function tors_eqt_tors
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             impr_same_impr                                                  !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function impr_same_impr(a,b)

        type(imprtype_type),intent(in)                :: a,b

        if((  ( abs(a%a1 - b%a1) .lt. 0.01) .and. ( abs(a%a2 - b%a2) .lt. 0.01)     .and.  &
              ( abs(a%a3 - b%a3) .lt. 0.01) .and. ( abs(a%a4 - b%a4) .lt. 0.01)  )   .or.  &
          (   ( abs(a%a1 - b%a4) .lt. 0.01) .and. ( abs(a%a2 - b%a3) .lt. 0.01)     .and.  &
              ( abs(a%a3 - b%a2) .lt. 0.01) .and. ( abs(a%a4 - b%a1) .lt. 0.01)  ) ) then
                 impr_same_impr = .true.
        else
                 impr_same_impr = .false.
        end if
end function impr_same_impr
!-----------------------------------------------------------------------------------------------------------------------------!
!                                                             impr_eqt_impr                                                   !
!-----------------------------------------------------------------------------------------------------------------------------!
logical function impr_eqt_impr(a,b)

        type(imprtype_type),intent(in)                :: a,b

        if ( ( abs(a%force - b%force) .lt. 0.001 ) .and. ( abs(a%dis - b%dis) .lt. 0.001 ) ) then
                 impr_eqt_impr = .true.
        else
                 impr_eqt_impr = .false.
        end if
end function impr_eqt_impr

end module array_analysis

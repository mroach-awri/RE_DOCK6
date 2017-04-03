!-------------------------------------------------------------------------------
!
! MIT License
! 
! Copyright (c) 2017 Michael Roach
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!-------------------------------------------------------------------------------
!
! Description: Module to read .mol2 files for conformation ensemble generator.
!
!-------------------------------------------------------------------------------

module mol2io_min
    
    use iso_fortran_env, only : ERROR_UNIT
    
    !---
    
    implicit none
    
    !---DATA TYPES---
    
    type :: mol2atom
        real :: xyz(3) ! atom coords
    end type mol2atom
    
    type :: mol2mol 
        character(len=30) :: headr(4) ! important stuff under @<TRIPOS>MOLECULE
        integer :: n_atoms ! needed for allocating and reading in the mol
        integer :: n_h_atoms = 0 ! number of non-hydrogen atoms
        integer :: n_bonds
        real :: gridScore
        character(len=110), allocatable :: atomline(:)
        character(len=50), allocatable :: bondline(:)
        type(mol2atom), allocatable :: atm(:)
    end type mol2mol
    
    contains
    
        subroutine readmol2(FH, m2m, io)
            
            !------
            
            character(len=150) :: strng
            character(len=50) :: arg, head
            integer :: i
            
            ! when reading in pre-mol comment lines (i.e. dock scores)
            character(len=20) :: cmmnt1, cmmnt2, cmmnt3
            
            ! when reading in the atom block
            integer :: jnk1, jnk4
            character(len=5) :: jnk2, jnk3
            character(len=10) :: jnk5
            real :: jnk6
            
            
            integer, intent(in) :: FH
            type(mol2mol), intent(out) :: m2m
            integer, intent(out) :: io
            
            !------
            
            ! iterate until start of molecule...
            do
                read(FH, '(A)', iostat=io) strng
                if (io /= 0) then
                    return
                end if
                if (trim(strng) == "@<TRIPOS>MOLECULE") then
                    exit
                end if
            end do
            
            ! line2, mol name
            read(FH, '(A)') strng
            read(strng, *) head
            m2m%headr(1) = trim(head)
            
            ! line3, the useful stuff
            read(FH, '(A)') strng
            read(strng, *) m2m%n_atoms, m2m%n_bonds
            m2m%headr(2) = trim(strng)
            
            ! we don't worry about anything else in the mol2 header area
            m2m%headr(3) = ""
            m2m%headr(4) = ""
            
            allocate(m2m%atomline(m2m%n_atoms))
            allocate(m2m%bondline(m2m%n_bonds))
            allocate(m2m%atm(m2m%n_atoms))
            
            ! iterate until we hit the atom block
            do
                read(FH, '(A)') strng
                if (trim(strng) == "@<TRIPOS>ATOM") then
                    exit
                end if
                if (trim(strng) /= "") then
                    read(strng,*) cmmnt1
                    if (trim(cmmnt1) == "#probe_cutoff:") then
                        m2m%headr(3) = trim(strng)
                    end if
                end if
            end do
            
            ! we read in all the values individually to standardize spacing
            do i=1, m2m%n_atoms
                read(FH, '(A)') strng
                
                read(strng,*) jnk1, jnk2, m2m%atm(i)%xyz(1), m2m%atm(i)%xyz(2), &
                    m2m%atm(i)%xyz(3), jnk3, jnk4, jnk5, jnk6
                write(m2m%atomline(i), FMT=10) jnk1, jnk2, m2m%atm(i)%xyz(1), m2m%atm(i)%xyz(2), &
                    m2m%atm(i)%xyz(3), jnk3, jnk4, jnk5, jnk6
10 FORMAT (I6, 1X, A7, 1X, 3F10.4, 1X, A3, 1X, I5, 1X, A6, 1X, F10.4)
                
                if (trim(jnk3) /= "H") then
                    m2m%n_h_atoms = m2m%n_h_atoms + 1
                end if
            end do
            
            ! there shouldn't be any blank lines between the atom block and the bonds block
            ! but just in case...
            do
                read(FH, '(A)') strng
                if (trim(strng) == "@<TRIPOS>BOND") then
                    exit
                end if
            end do
            
            !bonds block
            do i=1, m2m%n_bonds
                read(FH, '(A)') strng
                m2m%bondline(i) = trim(strng)
            end do
            
        end subroutine readmol2
        
        !---
        
        subroutine readmol2boxatoms(FH, xr, yr, zr, m2m, io)
            
            !------
            
            character(len=150) :: strng
            character(len=50) :: head
            integer :: i
            
            ! when reading in the atom block
            integer :: jnk1, allatoms
            character(len=5) :: jnk2, a_type
            
            ! hopefully the receptor shouldn't have more than 10000 atoms, we can always increase if needed
            type(mol2atom) :: bufferXYZ(10000)
            real :: xbuff, ybuff, zbuff
            
            integer, intent(in) :: FH
            real, intent(in) :: xr(2), yr(2), zr(2)
            type(mol2mol), intent(out) :: m2m
            integer, intent(out) :: io
            
            !------
            
            ! iterate until start of molecule
            do
                read(FH, '(A)', iostat=io) strng
                if (io /= 0) then
                    return
                end if
                if (trim(strng) == "@<TRIPOS>MOLECULE") then
                    exit
                end if
            end do
            
            ! line2, mol name, we don't use this for anything
            read(FH, '(A)') strng
            
            ! line3, the useful stuff
            read(FH, '(A)') strng
            read(strng, *) allatoms
            
            ! iterate until the start of the atom block
            do
                read(FH, '(A)') strng
                if (trim(strng) == "@<TRIPOS>ATOM") then
                    exit
                end if
            end do
            
            m2m%n_atoms = 0

            do i=1, allatoms
                read(FH, '(A)') strng
                read(strng,*) jnk1, jnk2, xbuff, ybuff, zbuff, a_type
                if (trim(a_type) /= "H") then
                    if ((xbuff > xr(1)) .AND. (xbuff < xr(2))) then
                        if ((ybuff > yr(1)) .AND. (ybuff < yr(2))) then
                            if ((zbuff > zr(1)) .AND. (zbuff < zr(2))) then
                                bufferXYZ((m2m%n_atoms + 1))%xyz(1) = xbuff
                                bufferXYZ((m2m%n_atoms + 1))%xyz(2) = ybuff
                                bufferXYZ((m2m%n_atoms + 1))%xyz(3) = zbuff
                                m2m%n_atoms = m2m%n_atoms + 1
                            end if
                        end if
                    end if
                end if
            end do
            
            allocate(m2m%atm(m2m%n_atoms))
            
            do i=1, m2m%n_atoms
                m2m%atm(i) = bufferXYZ(i)
            end do
            
        end subroutine readmol2boxatoms
        
        !---
        
        subroutine writeProbeRotation(FH, m2m, m2m_rot)
            
            ! script to write the probe rotation using the rotated probe's
            ! xyz coords but the master probe's atomblock lines
            integer, intent(in) :: FH
            type(mol2mol), intent(in) :: m2m
            type(mol2mol), intent(in) :: m2m_rot
            integer :: i
            
            
            write(FH, '(A17)') "@<TRIPOS>MOLECULE"
            
            do i=1, 4
                write(FH, '(A)') trim(m2m%headr(i))
            end do
            
            write(FH, '(A13)') "@<TRIPOS>ATOM"
            
            ! It might be more robust to re-read in all the individual values but the spacing is standardized
            do i=1, m2m%n_atoms
                write(FH, '(A16, 3F10.4, A54)') m2m%atomline(i)(1:15), m2m_rot%atm(i)%xyz(1), &
                    m2m_rot%atm(i)%xyz(2),m2m_rot%atm(i)%xyz(3), m2m%atomline(i)(45:103)
            end do
            
            write(FH, '(A13)') "@<TRIPOS>BOND"
            
            do i=1, m2m%n_bonds
                write(FH, '(A)') trim(m2m%bondline(i))
            end do
            
            write(FH,*) " "
            
        end subroutine writeProbeRotation
    
end module mol2io_min

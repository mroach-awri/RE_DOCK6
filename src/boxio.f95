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
! Description: Module to read UCSF Dock box .pdb format file for conformation 
!              ensemble generator.
!
!-------------------------------------------------------------------------------

module boxio
    
    implicit none
    
    type :: boxrange
        real :: xrange(2)
        real :: yrange(2)
        real :: zrange(2)
    end type boxrange
    
    contains 
        
        subroutine readbox(FH, boxpdb, io)
            
            implicit none
            
            integer, intent(IN) :: FH
            integer, intent(OUT) :: io
            type(boxrange), intent(OUT) :: boxpdb
            
            character(len=200) :: line
            character(len=10) :: rmk
            integer :: atomcount
            
            integer :: jnk1, jnk4
            character(len=4) :: jnk2, jnk3
            real :: x,y,z
            
            !---
            
            boxpdb%xrange(1:2) = (/ 9999.0, -9999.0 /)
            boxpdb%yrange(1:2) = (/ 9999.0, -9999.0 /)
            boxpdb%zrange(1:2) = (/ 9999.0, -9999.0 /)
            
            atomcount = 0
            
            do
                read(FH, '(A)') line
                if (io /= 0) then
                    return
                end if
                
                if (trim(line) /= "") then
                    read(line,*) rmk
                    if (trim(rmk) .eq. "ATOM") then
                        atomcount = atomcount + 1
                        
                        read(line,*) rmk, jnk1, jnk2, jnk3, jnk4, x, y, z
                        
                        if (x < boxpdb%xrange(1)) then
                            boxpdb%xrange(1) = x
                        end if
                        if (x > boxpdb%xrange(2)) then
                            boxpdb%xrange(2) = x
                        end if
                        if (y < boxpdb%yrange(1)) then
                            boxpdb%yrange(1) = y
                        end if
                        if (y > boxpdb%yrange(2)) then
                            boxpdb%yrange(2) = y
                        end if
                        if (z < boxpdb%zrange(1)) then
                            boxpdb%zrange(1) = z
                        end if
                        if (z > boxpdb%zrange(2)) then
                            boxpdb%zrange(2) = z
                        end if
                        
                    end if
                end if
                
                if (atomcount == 8) then
                    return
                end if
            end do
            
            
        end subroutine readbox
        
end module boxio

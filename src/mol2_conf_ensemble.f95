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
! Description: Stand-alone program, intended for use with RE_DOCK6. Will take a
!              receptor mol2 file, box.pdb file and a ligand mol2 file, and will
!              generate a conformation ensemble of each ligand to fit the box
!              and receptor structure.
!
!-------------------------------------------------------------------------------

program mol2_conformation_ensemble
    
    use iso_fortran_env, only : ERROR_UNIT
    use mol2io_min
    use boxio
    
    implicit none
    
    type(mol2mol) :: recMol2
    type(mol2mol) :: probeMol2
    type(mol2mol), allocatable :: probe_rots(:)
    
    type(boxrange) :: box
    
    ! OPTIONALLY PARSABLE, everything is still positional arguments though because lazy
    integer :: nsteps = 3
    real :: stepsize = 1.5
    real :: tether_dist = 5.0
    
    ! could possibly make these parsable in the future
    real :: atom_clash_dist = 1.5
    real :: cent_clash_dist = 2.0

    integer :: xsteps, ysteps, zsteps
    
    integer :: FHMOD = 100
    integer :: FHBOX = 200
    integer :: FHPRB = 300
    integer :: FHOUT = 400

    integer :: io
    integer :: outconfs = 0
    integer :: nconfs = 0
    integer :: i, ii, iii, iv, v, vi, rot
    real :: dtheta, theta, pi ! rotation stuff
    real, allocatable :: xyz1(:,:), xyz2(:,:), xyz3(:,:)
    real :: xyz_c(3), box_cent(3)
    real :: xdif, ydif, zdif, dist
    real, allocatable :: costh(:), sinth(:)
    logical :: clash, clash2, tether
    
    
    character(len=100) :: arg
    character(len=10) :: head
    
    call get_command_argument(1, arg)
    open(FHMOD, file=trim(arg), action="read", status="old", iostat=io)
    if(io/= 0) then
        write(ERROR_UNIT,*) "File: '", trim(arg), "' could not be opened for reading"
        write(ERROR_UNIT,*) "Usage: conf_ensmbl  receptor.mol2  box.pdb  &
            ligand(s).mol2  out_ligand(s).mol2  [n_rot]  [d_step] [teth_dist]"
        call exit(1)
    end if
    
    call get_command_argument(2, arg)
    open(FHBOX, file=trim(arg), action="read", status="old", iostat=io)
    if (io /= 0) then
        write(ERROR_UNIT,*) "File: '", trim(arg), "' could not be opened for reading"
        write(ERROR_UNIT,*) "Usage: conf_ensmbl  receptor.mol2  box.pdb  &
            ligand(s).mol2  out_ligand(s).mol2  [n_rot]  [d_step] [teth_dist]"
        call exit(1)
    end if
    
    call get_command_argument(3, arg)
    open(FHPRB, file=trim(arg), action="read", status="old", iostat=io)
    if (io /= 0) then
        write(ERROR_UNIT,*) "File: '", trim(arg), "' could not be opened for reading"
        write(ERROR_UNIT,*) "Usage: conf_ensmbl  receptor.mol2  box.pdb  &
            ligand(s).mol2  out_ligand(s).mol2  [n_rot]  [d_step] [teth_dist]"
        call exit(1)
    end if
    
    call get_command_argument(4, arg)
    open(FHOUT, file=trim(arg), action="write", iostat=io)
    if (io /= 0) then
        write(ERROR_UNIT,*) "File: '", trim(arg), "' could not be opened for writing"
        write(ERROR_UNIT,*) "Usage: conf_ensmbl  receptor.mol2  box.pdb  &
            ligand(s).mol2  out_ligand(s).mol2  [n_rot]  [d_step] [teth_dist]"
        call exit(1)
    end if
    
    call get_command_argument(5, arg)
    if (trim(arg) /= "") then
        read(arg,*) nsteps
    end if
    
    call get_command_argument(6, arg)
    if (trim(arg) /= "") then
        read(arg,*) stepsize
    end if
    
    call get_command_argument(7,arg)
    if(trim(arg) /= "") then
        read(arg,*) tether_dist
    end if
    
    allocate(costh(nsteps))
    allocate(sinth(nsteps))
    
    !---READ IN THE TARGET ATOM MODEL---    
    call readbox(FHBOX, box, io)
    if (io /= 0) then
        write(ERROR_UNIT,*) "failed to get coords from .pdb box file"
        write(ERROR_UNIT,*) "Usage: probe_ensmbl  receptor.mol2  in_box.pdb  &
            probe.mol2  out_probe_orients.mol2"
    end if
    
    call readmol2boxatoms(FHMOD, box%xrange, box%yrange, box%zrange, recMol2, io)
    if (io /= 0) then
        write(ERROR_UNIT,*) "Error reading receptor mol2 file"
        call exit(1)
    end if
    
    
!---LOOP TO ITERATE THROUGH THE PROBE LIBRARY, WILL STOP WHEN EOF REACHED---
    write(ERROR_UNIT,*) "Generating ligand conformation ensemble..."
    do
        !---READ IN THE PROBE---
        call readmol2(FHPRB, probeMol2, io)
        if (io /= 0) then
            exit
        end if
        
        nconfs = nconfs + 1
        
        !---GENERATE THE ROTATIONS---
        pi = 4. * atan(1.)
        theta = 2. * pi / nsteps
        
        do i=1,nsteps
            dtheta = theta * i
            costh(i) = cos(dtheta)
            sinth(i) = sin(dtheta)
        end do
        
        ! allocate variables
        allocate(xyz1(probeMol2%n_atoms,3))
        allocate(xyz2(probeMol2%n_atoms,3))
        allocate(xyz3(probeMol2%n_atoms,3))
        allocate(probe_rots(nsteps*nsteps*nsteps))
        
        do i=1, nsteps*nsteps*nsteps
            allocate(probe_rots(i)%atm(probeMol2%n_atoms))
        end do
        
        ! rotation step, nested loops are actually a bit faster than matrix multiplication
        rot = 1
        ! x-rotation
        do ii=1,nsteps
            xyz2(:,1) = probeMol2%atm(:)%xyz(1)
            xyz2(:,2) = costh(ii)*probeMol2%atm(:)%xyz(2) + sinth(ii)*probeMol2%atm(:)%xyz(3)
            xyz2(:,3) = -probeMol2%atm(:)%xyz(2)*sinth(ii) + probeMol2%atm(:)%xyz(3)*costh(ii)
            ! y-rotation
            do iii=1,nsteps
                xyz3(:,1) = xyz2(:,1)*costh(iii) + xyz2(:,3)*sinth(iii)
                xyz3(:,2) = xyz2(:,2)
                xyz3(:,3) = -xyz2(:,1)*sinth(iii) + xyz2(:,3)*costh(iii)
                ! z-rotation
                do iv=1,nsteps
                    probe_rots(rot)%atm(:)%xyz(1) = xyz3(:,1)*costh(iv) - xyz3(:,2)*sinth(iv)
                    probe_rots(rot)%atm(:)%xyz(2) = xyz3(:,1)*sinth(iv) + xyz3(:,2)*costh(iv)
                    probe_rots(rot)%atm(:)%xyz(3) = xyz3(:,3)

                    rot = rot + 1
                end do
            end do
        end do
        
        
        !---ITERATE ACCROSS THE GRID AND OUTPUT THE ROTATED PROBES---
            ! xyz_c = centroid for a grid iteration step
            
        xsteps = int((box%xrange(2) - box%xrange(1)) / stepsize)
        ysteps = int((box%yrange(2) - box%yrange(1)) / stepsize)
        zsteps = int((box%zrange(2) - box%zrange(1)) / stepsize)
        
        box_cent(1) = (box%xrange(2) + box%xrange(1)) / 2
        box_cent(2) = (box%yrange(2) + box%yrange(1)) / 2
        box_cent(3) = (box%zrange(2) + box%zrange(1)) / 2
        
        do i=1, xsteps
            xyz_c(1) = (i * stepsize) + box%xrange(1)
            do ii=1, ysteps
                xyz_c(2) = (ii * stepsize) + box%yrange(1)
                do iii=1, zsteps
                    xyz_c(3) = (iii * stepsize) + box%zrange(1) 
                    
                    tether = .FALSE.
                    clash = .FALSE.
                    
                    ! check if centroid is within tether dist of box center
                    call DistDiff(box_cent(1), box_cent(2), box_cent(3), xyz_c(1), xyz_c(2), xyz_c(3), dist)
                    if (dist > tether_dist) then
                        tether = .TRUE.
                    end if
                    
                    ! Check if the centroid clashes with any target atoms
                    if (tether .eqv. .FALSE.) then
                        do iv=1, recMol2%n_atoms
                            call DistDiff(xyz_c(1), xyz_c(2), xyz_c(3), recMol2%atm(iv)%xyz(1), &
                                    recMol2%atm(iv)%xyz(2), recMol2%atm(iv)%xyz(3), dist)
                            if (dist < cent_clash_dist) then 
                                clash = .TRUE.
                                exit
                            end if
                        end do
                    else if (tether .eqv. .TRUE.) then
                        clash = .TRUE.
                    end if
                    
                    if (clash .eqv. .FALSE.) then
                        do iv=1, (nsteps*nsteps*nsteps)
                            ! get centroid for orient
                            xdif = ( sum(probe_rots(iv)%atm(:)%xyz(1)) / probeMol2%n_atoms)
                            ydif = ( sum(probe_rots(iv)%atm(:)%xyz(2)) / probeMol2%n_atoms)
                            zdif = ( sum(probe_rots(iv)%atm(:)%xyz(3)) / probeMol2%n_atoms)
                            
                            ! get the offset between orient centroid and grid iteration index
                            xdif = (xdif - xyz_c(1))
                            ydif = (ydif - xyz_c(2))
                            zdif = (zdif - xyz_c(3))
                            
                            ! zero the orient to the grid iteration index and then
                            ! check for clashing atoms with receptor
                            probe_rots(iv)%atm(:)%xyz(1) = probe_rots(iv)%atm(:)%xyz(1) - xdif
                            probe_rots(iv)%atm(:)%xyz(2) = probe_rots(iv)%atm(:)%xyz(2) - ydif
                            probe_rots(iv)%atm(:)%xyz(3) = probe_rots(iv)%atm(:)%xyz(3) - zdif

                            clash2 = .FALSE.
                            
                            ! check for atoms in box and not clashing
                            do v=1, probeMol2%n_atoms
                                ! make sure atoms are within the box
                                if ( (probe_rots(iv)%atm(v)%xyz(1) < box%xrange(1)) .OR. &
                                            (probe_rots(iv)%atm(v)%xyz(1) > box%xrange(2)) ) then
                                    clash2 = .TRUE.
                                    exit
                                else if ( (probe_rots(iv)%atm(v)%xyz(2) < box%yrange(1)) .OR. &
                                            (probe_rots(iv)%atm(v)%xyz(2) > box%yrange(2)) ) then
                                    clash2 = .TRUE.
                                    exit
                                else if ( (probe_rots(iv)%atm(v)%xyz(3) < box%zrange(1)) .OR. &
                                            (probe_rots(iv)%atm(v)%xyz(3) > box%zrange(2)) ) then
                                    clash2 = .TRUE.
                                    exit
                                end if
                                
                                ! now check if they clash with the receptor
                                if (clash2 .EQV. .FALSE.) then
                                    do vi=1, recMol2%n_atoms
                                        call DistDiff(probe_rots(iv)%atm(v)%xyz(1), probe_rots(iv)%atm(v)%xyz(2), &
                                                    probe_rots(iv)%atm(v)%xyz(3), recMol2%atm(vi)%xyz(1), &
                                                    recMol2%atm(vi)%xyz(2), recMol2%atm(vi)%xyz(3), dist)
                                        if (dist < atom_clash_dist) then 
                                            clash2 = .TRUE.
                                            exit
                                        end if
                                    end do
                                end if
                            end do
                            
                            if (clash2 .eqv. .FALSE.)then
                                ! write to mol2 file
                                call writeProbeRotation(FHOUT, probeMol2, probe_rots(iv))
                                outconfs = outconfs + 1
                            end if
                        end do
                    end if
                end do
            end do
        end do
        
        ! deallocate variables
        deallocate(xyz1)
        deallocate(xyz2)
        deallocate(xyz3)
        deallocate(probe_rots)
    
    ! end loop for iterating over probe mol2 file
    end do
    write(ERROR_UNIT,*) "Generated ", outconfs, " total ligand orientations over", nconfs, "ligand conformers"
    
end program mol2_conformation_ensemble

!---SUBROUTINES---

subroutine DistDiff(x1,y1,z1,x2,y2,z2,dist)
  
  implicit none

  real :: x1, y1, z1, x2, y2, z2, dist

  dist = ((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)) + ((z2-z1)*(z2-z1))
  dist = SQRT(dist)

end subroutine DistDiff

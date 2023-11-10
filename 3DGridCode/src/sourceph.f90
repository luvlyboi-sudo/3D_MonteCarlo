module sourceph_mod
!! Module contains the routines to inialise a neutron, i.e different light sources.
    implicit none

    contains
        subroutine isotropic_point_src(packet, grid)
        !! set intial neutron position at (0.0, 0.0, 0.0) and sample neutron direction in an isotropic manner.

            use constants,    only : TWOPI, wp
            use gridset_mod,  only :cart_grid
            use neutron_class, only : neutron
            use random_mod,   only : ran2

            !> neutron object
            type(neutron),    intent(out) :: packet
            !> grid object
            type(cart_grid), intent(in)  :: grid

            !set packet position
            packet%pos%z = 0.0_wp
            packet%pos%x = 0.0_wp
            packet%pos%y = 0.0_wp

            ! set packet cosines
            packet%phi  = ran2()*twoPI
            packet%cosp = cos(packet%phi)
            packet%sinp = sin(packet%phi)
            packet%cost = 2._wp*ran2()-1._wp
            packet%sint = sqrt(1._wp - packet%cost**2)

            ! set direction vector
            packet%dir%x = packet%sint * packet%cosp  
            packet%dir%y = packet%sint * packet%sinp
            packet%dir%z = packet%cost

            ! set packet voxel
            packet%xcell=int(grid%nxg*(packet%pos%x+grid%dim%x)/(2._wp*grid%dim%x))+1
            packet%ycell=int(grid%nyg*(packet%pos%y+grid%dim%y)/(2._wp*grid%dim%y))+1
            packet%zcell=int(grid%nzg*(packet%pos%z+grid%dim%z)/(2._wp*grid%dim%z))+1

            packet%tflag = .false.

        end subroutine isotropic_point_src

        subroutine uniform_sphere_src(packet, grid, radius)
            use constants,    only : TWOPI, wp
            use gridset_mod,  only : cart_grid
            use neutron_class, only : neutron
            use random_mod,   only : ran2
        
            type(neutron),    intent(out) :: packet
            type(cart_grid), intent(in)  :: grid
            real(kind=wp),   intent(in)  :: radius
        
            real(kind=wp) :: x, y, z, r2
    
            do
                x = (ran2() * 2._wp - 1._wp) * radius
                y = (ran2() * 2._wp - 1._wp) * radius
                z = (ran2() * 2._wp - 1._wp) * radius
                r2 = x**2 + y**2 + z**2
                if (r2 <= radius**2) exit
            end do
        
            ! Set packet position
            packet%pos%x = x
            packet%pos%y = y
            packet%pos%z = z
        
            ! set packet cosines
            packet%phi  = ran2()*twoPI
            packet%cosp = cos(packet%phi)
            packet%sinp = sin(packet%phi)
            packet%cost = 2._wp*ran2()-1._wp
            packet%sint = sqrt(1._wp - packet%cost**2)

            ! set direction vector
            packet%dir%x = packet%sint * packet%cosp  
            packet%dir%y = packet%sint * packet%sinp
            packet%dir%z = packet%cost

            ! set packet voxel
            packet%xcell=int(grid%nxg*(packet%pos%x+grid%dim%x)/(2._wp*grid%dim%x))+1
            packet%ycell=int(grid%nyg*(packet%pos%y+grid%dim%y)/(2._wp*grid%dim%y))+1
            packet%zcell=int(grid%nzg*(packet%pos%z+grid%dim%z)/(2._wp*grid%dim%z))+1

            packet%tflag = .false.
        end subroutine uniform_sphere_src

        subroutine init_neutron_hemispheres(packet, grid, radius, separation, hemisphere)
            use constants,    only : TWOPI, wp, PI
            use gridset_mod,  only : cart_grid
            use neutron_class, only : neutron
            use random_mod,   only : ran2
            
            type(neutron),    intent(out) :: packet
            type(cart_grid), intent(in)  :: grid
            real(kind=wp),   intent(in)  :: radius, separation
            integer,         intent(in)  :: hemisphere  ! 1 for top, -1 for bottom
            
            real(kind=wp) :: theta, phi, sin_theta, z_offset, mag
            
            phi = ran2() * TWOPI
            theta = acos(1._wp - 2._wp * ran2())
            sin_theta = sin(theta)
            
            ! Determine the z_offset
            z_offset = hemisphere * separation / 2.0_wp
            
            packet%pos%x = radius * sin_theta * cos(phi)
            packet%pos%y = radius * sin_theta * sin(phi)
            packet%pos%z = radius * cos(theta) + z_offset
            
            ! Make the direction vector point towards the center of the sphere
            packet%dir%x = -packet%pos%x
            packet%dir%y = -packet%pos%y
            packet%dir%z = -packet%pos%z * hemisphere
            
            ! Normalize the direction vector
            mag = sqrt(packet%dir%x**2 + packet%dir%y**2 + packet%dir%z**2)
            packet%dir%x = packet%dir%x / mag
            packet%dir%y = packet%dir%y / mag
            packet%dir%z = packet%dir%z / mag
            
            packet%xcell = int(grid%nxg * (packet%pos%x + grid%dim%x) / (2._wp * grid%dim%x)) + 1
            packet%ycell = int(grid%nyg * (packet%pos%y + grid%dim%y) / (2._wp * grid%dim%y)) + 1
            packet%zcell = int(grid%nzg * (packet%pos%z + grid%dim%z) / (2._wp * grid%dim%z)) + 1
            
            packet%tflag = .false.
            
        end subroutine init_neutron_hemispheres
        
    
        
        
end module sourceph_mod

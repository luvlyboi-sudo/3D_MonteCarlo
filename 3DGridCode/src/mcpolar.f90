program mcpolar 

    !imports
    use constants,                only : resdir, wp
    use gridset_mod,              only : gridset, cart_grid
    use inttau2,                  only : tauint1
    use optical_properties_class, only : optical_properties, init_neutron_properties 
    use neutron_class,            only : neutron
    use random_mod,               only : ran2, init_seed
    use sourceph_mod,             only : isotropic_point_src, uniform_sphere_src,init_neutron_hemispheres
    use utils,                    only : set_directories, str
    use writer_mod,               only : writer

    implicit none

    !> variable that holds all information about the neutron to be simulated
    type(neutron) :: packet
    !> variable that holds the 3D grid information
    type(cart_grid)  :: grid
    !> optical properties variable
    type(optical_properties) :: opt_prop
    !> number of photons/neutrons to run in the simulation
    integer :: nphotons
    !> probabilites for absorb, scattering, fission
    real(kind=wp) :: cumulative_scatter, cumulative_absorb
    !> Storing the generated random number
    real(kind=wp) :: rand_val
    !> user defined seed
    integer :: seed
    !> temp variables related to I/O from param file
    integer :: nxg, nyg, nzg
    !> loop variable
    integer :: j
    !> file handle
    integer :: u
    !> temp variables related to I/O from param file
    real(kind=wp) :: xmax, ymax, zmax, mus, mua, muf, hgg


    type(neutron), allocatable :: neutron_bank(:)  !> The bank of neutrons to simulate 
    type(neutron), allocatable :: local_neutron_bank(:)  !> The bank of neutrons to simulate for the next gen
    !type(real(kind=wp)), allocatable :: keff(:) !> Array to store the critical coefficients of each generation  
    integer :: current_gen_num !> What generation is it 
    integer, parameter :: max_bank_size = 10000 !> Maximum number of neutrons in the bank before we terminate the program
    integer :: num_fission !> Number of neutrons generated in current generation
    integer :: num_fission_last !> Number of neutrons generated in the last generation 

    ! ########### [code to set up optical properties, grid, and scattering] ##############

    !set directory paths
    call set_directories()
    
    !set random seed
    seed = 42
    call init_seed(seed)
    
    !**** Read in parameters from the file input.params
    open(newunit=u,file=trim(resdir)//'input.params',status='old')
    read(u,*) nphotons ! This is really the number of neutrons, too lazy to change it everywhere 
    read(u,*) xmax
    read(u,*) ymax
    read(u,*) zmax
    read(u,*) nxg
    read(u,*) nyg
    read(u,*) nzg
    read(u,*) mus
    read(u,*) mua
    read(u,*) muf
    read(u,*) hgg
    close(u)
    
    call init_neutron_properties(mus, mua, muf, hgg, opt_prop) ! Set optical properties for fission
    call gridset(grid, opt_prop, nxg, nyg, nzg, xmax, ymax, zmax) ! Set up grid

    ! ################## Intialising Variables ##############
    current_gen_num = 1
    num_fission = 0 ! Update this later
    num_fission_last = nphotons ! For the first generation
    allocate(neutron_bank(max_bank_size)) ! Start by allocating nphotons spaces in the neutron bank, need to reallocate later
    allocate(local_neutron_bank(max_bank_size + 1000)) ! Ensuring we don't get a memory error so local bank must be greater than neutron_bank
    !allocate(keff(10)) 

    do j = 1, nphotons/2 !ensure nphotons is even :)

        !call isotropic_point_src(neutron_bank(j),grid) ! Populate the bank for first generation
        !call uniform_sphere_src(packet, grid, 1._wp)
        call init_neutron_hemispheres(neutron_bank(j),grid,10._wp,0._wp,1)
        print*, neutron_bank(j)%pos
        
    end do

    do j = 1, nphotons/2 !ensure nphotons is even :)

        call init_neutron_hemispheres(neutron_bank(j),grid,10._wp,0._wp,-1)
        print*, neutron_bank(j)%pos
        
    end do

    do while (.true.)

        ! ################## The Simulation #######################

        do j = 1, size(neutron_bank)
            
            packet = neutron_bank(j)
            
            call tauint1(packet, grid) ! Simulate the movement of neutron
            
            do while (.not. packet%tflag) ! While packet is alive
                
                print*, packet%pos
                ! Random number for determining if fission, scattering or absorption 
                rand_val = ran2()

                ! Calculate cumulative probabilities
                cumulative_scatter = opt_prop%mus / opt_prop%kappa
                cumulative_absorb = cumulative_scatter + opt_prop%mua / opt_prop%kappa

                if (rand_val < cumulative_scatter) then
                    ! Neutron is scattered
                    !print*,'Scattering'
                    call packet%scatter(opt_prop) 

                elseif (rand_val >= cumulative_scatter .and. rand_val < cumulative_absorb) then 
                    ! Neutron is absorbed 
                    !print*,'Absorbed'
                    exit

                else
                    !print*,'Fission'
                    ! Else the neutrons causes fission 
                    call packet%fission(opt_prop,local_neutron_bank,num_fission)
                    exit
                    
                end if

            end do

        end do
        
        print *, real(num_fission)/real(num_fission_last)

        ! Check if the number of fission neutrons exceeds the maximum bank size
        if (num_fission >= max_bank_size) then
            print *, "Stopping condition met: Number of fission neutrons exceeds maximum bank size."
            exit  ! Exit the simulation loop
        endif

        if (num_fission > 0) then
            ! Reallocate neutron_bank to the size of the new fission neutrons
            neutron_bank(1:num_fission) = local_neutron_bank(1:num_fission)
    
        else
            ! No fission occurred, end the simulation. This stopping condition is num_fission == 0
            print *, "Stopping condition met: Number of fission neutrons is 0"
            exit
        endif


        !Update for next generation
        !keff(current_gen_num) = real(num_fission)/real(num_fission_last)
        num_fission_last = num_fission
        num_fission = 0
        current_gen_num = current_gen_num + 1 ! Update the generation number

        if (current_gen_num > 200) then 
            exit
        endif

    end do 

end program mcpolar
    

    

    


    


program kep
    implicit none
    ! Precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: int64 = selected_int_kind(18)
    
    ! Physical parameters
    integer, parameter :: L = 200
    integer, parameter :: N = 100
    real(dp), parameter :: T_left = 4.0_dp
    real(dp), parameter :: T_right = 1.0_dp
    real(dp), parameter :: beta = 0.0_dp ! Energy exponent f(e) = e^beta
    
    ! !!! FIX: Ya no usamos rho_left/right para inyeccion, N es constante.
    ! Estas variables se pueden ignorar para el "Case 1" (thermal walls only)
    real(dp), parameter :: rho_left = 0.3_dp
    real(dp), parameter :: rho_right = 0.8_dp
    
    ! Simulation parameters
    integer(int64), parameter :: n_steps = 400000000000_int64
    integer(int64), parameter :: relaxation_steps = 9000000000_int64
    integer(int64), parameter :: measurement_steps = n_steps - relaxation_steps
    integer(int64), parameter :: sample_interval = 100000_int64
    integer(int64), parameter :: file_output_interval = 1000000000_int64
    
    ! Lattice arrays
    integer :: occupation(0:L+1)
    real(dp) :: energy(0:L+1)
    real(dp) :: temperature(0:L+1) ! Perfil teorico para inicializacion
    
    ! Averaged profiles (First Moment <x>)
    real(dp) :: avg_density(L)
    real(dp) :: avg_energy_density(L)
    real(dp) :: avg_temperature(L)
    
    ! !!! FIX: Arrays for Second Moments <x^2> for Error Calculation
    real(dp) :: avg_density_sq(L)
    real(dp) :: avg_energy_density_sq(L)
    
    ! Arrays for final error output
    real(dp) :: error_density(L)
    real(dp) :: error_energy(L)
    
    ! KMC variables
    real(dp) :: transition_rates(0:L)
    real(dp) :: total_rate
    real(dp) :: current_time
    
    ! Statistics
    integer(int64) :: sample_count ! Changed to int64 to be safe
    real(dp) :: total_energy
    integer :: n_particles
    
    ! Loop variables
    integer :: i, j, bond, iaux, iaux2
    integer(int64) :: step
    real(dp) :: dt, random_val, pair_energy, energy_share, local_temp
    
    ! Initialize random seed
    call initialize_random_seed()
    
    ! Initialize lattice
    call initialize_lattice()
    
    ! Initialize averaged profiles
    avg_density = 0.0_dp
    avg_energy_density = 0.0_dp
    avg_temperature = 0.0_dp
    
    ! !!! FIX: Initialize squared arrays
    avg_density_sq = 0.0_dp
    avg_energy_density_sq = 0.0_dp
    
    sample_count = 0
    
    write(*,*) '========================================'
    write(*,*) '            KEP SIMULATION (FIXED)      '
    write(*,*) '========================================'
    write(*,*) 'Lattice size L = ', L
    write(*,*) 'Number of particles N = ', N
    write(*,*) 'Left temperature T_L = ', T_left
    write(*,*) 'Right temperature T_R = ', T_right
    write(*,*) 'Total steps = ', n_steps
    write(*,*) 'Relaxation steps = ', relaxation_steps
    write(*,*) 'Measurement steps = ', measurement_steps
    write(*,*) '========================================'
    
    ! Run simulation
    do step = 1, n_steps
        
        ! Perform KMC step
        call perform_kmc_step()
        
        ! Measurements
        if (step > relaxation_steps) then
            if (mod(step - relaxation_steps, sample_interval) == 0) then
                call take_measurement()
            end if
        end if
        
        ! Write intermediate files
        if (step > relaxation_steps) then
            if (mod(step - relaxation_steps, file_output_interval) == 0) then
                call write_intermediate_results()
            end if
        end if
        
        ! Progress reporting
        if (mod(step, n_steps/10_int64) == 0) then
            if (step <= relaxation_steps) then
                write(*,*) 'Relaxation progress: ', int(100.0_dp * step / relaxation_steps), '%'
            else
                write(*,*) 'Measurement progress: ', int(100.0_dp * (step - relaxation_steps) / measurement_steps), '%'
            end if
        end if
    end do
    
    ! Calculate final averages and errors
    call calculate_final_stats()
    
    ! Print results
    call print_results()
    
    write(*,*) '========================================'
    write(*,*) 'SIMULATION COMPLETED SUCCESSFULLY!'
    write(*,*) '========================================'
    
contains

    subroutine initialize_random_seed()
        implicit none
        integer :: i, n, clock
        integer, allocatable :: seed_array(:)
        call random_seed(size=n)
        allocate(seed_array(n))
        call system_clock(count=clock)
        seed_array = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put=seed_array)
        deallocate(seed_array)
    end subroutine initialize_random_seed
    
    subroutine initialize_lattice()
        implicit none
        integer :: i, j, temp_occ
        real(dp) :: local_temp, temp_energy
        
        ! Initialize arrays
        occupation = 0
        energy = 0.0_dp
        temperature = 0.0_dp
        
        ! Set up theoretical linear temperature profile for initialization
        do i = 0, L+1
            local_temp = T_left + (T_right - T_left) * real(i, dp) / real(L+1, dp)
            temperature(i) = local_temp
        end do
        
        ! Place N particles
        do i = 1, N
            occupation(i) = 1
            energy(i) = -temperature(i) * log(max(get_random(), 1.0e-10_dp))
        end do

        ! Fisher-Yates shuffle (internal sites only)
        do i = L, 2, -1
            j = int(get_random() * i) + 1
            temp_occ = occupation(i)
            occupation(i) = occupation(j)
            occupation(j) = temp_occ
            temp_energy = energy(i)
            energy(i) = energy(j)
            energy(j) = temp_energy
        end do
        
        ! Initialize boundary conditions (Walls)
        call update_boundaries()
        
        ! Observables
        n_particles = sum(occupation(1:L))
        total_energy = sum(energy(1:L))
        current_time = 0.0_dp
        
        write(*,*) 'Lattice initialized. Particles: ', n_particles
    end subroutine initialize_lattice
    
    subroutine perform_kmc_step()
        implicit none
        
        ! 1. Update boundaries energies BEFORE calculating rates
        !    Esto asegura que la pared siempre tenga una energia fresca "termalizada"
        call update_boundaries()
        
        ! 2. Calculate transition rates
        call calculate_transition_rates()
        
        ! Safety check
        if (total_rate <= 1.0e-10_dp) return
        
        ! 3. Time step
        dt = -log(max(get_random(), 1.0e-10_dp)) / total_rate
        current_time = current_time + dt
        
        ! 4. Select bond
        random_val = get_random() * total_rate
        bond = 0
        do i = 0, L
            if (transition_rates(i) > random_val) then
                bond = i
                exit
            end if
        end do
        
        i = bond
        j = bond + 1
        
        ! 5. Perform transition
        iaux = occupation(i) + occupation(j)
        
        if (iaux == 2) then
            ! Collision (Energy exchange)
            ! Esto ocurre en el bulk (1+1) o en la pared (Wall+1)
            pair_energy = energy(i) + energy(j)
            energy(i) = get_random() * pair_energy
            energy(j) = pair_energy - energy(i)
            
        else if (iaux == 1) then
            ! Jump (Particle-Hole)
            ! !!! FIX: Esto solo puede ocurrir en el BULK (1..L).
            ! Las tasas de transicion para la pared saltando al vacio son 0 (ver calculate_transition_rates)
            ! por lo tanto aqui no hace falta un if extra, pero logicamente solo entrara aqui si no es una pared intentando saltar.
            pair_energy = energy(i) + energy(j)
            iaux2 = occupation(i)
            occupation(i) = occupation(j)
            occupation(j) = iaux2
            energy(i) = occupation(i) * pair_energy
            energy(j) = occupation(j) * pair_energy
        end if
        
        ! Update observables
        n_particles = sum(occupation(1:L))
        total_energy = sum(energy(1:L))
    end subroutine perform_kmc_step
    
    subroutine calculate_transition_rates()
        implicit none
        
        ! !!! FIX CRUCIAL:
        ! Bond 0 (Pared Izq <-> Sitio 1):
        ! Si sitio 1 esta vacio, tasa = 0 (Pared no puede saltar dentro).
        ! Si sitio 1 ocupado, tasa = f(E_wall + E_1) (Colision permitida).
        if (occupation(1) == 0) then
            transition_rates(0) = 0.0_dp
        else
            transition_rates(0) = (energy(0) + energy(1))**beta
        end if
        
        total_rate = transition_rates(0)
        
        ! Bulk bonds (1 to L-1) - Standard KEP logic
        do i = 1, L-1
            transition_rates(i) = total_rate + (energy(i) + energy(i+1))**beta
            total_rate = transition_rates(i)
        end do

        ! Bond L (Sitio L <-> Pared Der):
        ! Si sitio L esta vacio, tasa = 0 (Pared der no puede saltar dentro).
        if (occupation(L) == 0) then
            transition_rates(L) = total_rate ! Sumamos 0
        else
            transition_rates(L) = total_rate + (energy(L) + energy(L+1))**beta
        end if
        
        total_rate = transition_rates(L)
    end subroutine calculate_transition_rates
    
    subroutine update_boundaries()
        implicit none
        
        ! !!! FIX: Thermal Walls Implementation (Case 1)
        ! 1. Los sitios de la pared SIEMPRE estan ocupados (actuan como particulas fijas)
        occupation(0) = 1
        occupation(L+1) = 1
        
        ! 2. Su energia se resetea aleatoriamente siguiendo distribucion de Boltzmann (Exponencial)
        ! E ~ -T * ln(rand)
        energy(0) = -T_left * log(max(get_random(), 1.0e-10_dp))
        energy(L+1) = -T_right * log(max(get_random(), 1.0e-10_dp))
        
    end subroutine update_boundaries
    
    subroutine take_measurement()
        implicit none
        
        sample_count = sample_count + 1
        
        do i = 1, L
            ! Moment 1: <x>
            avg_density(i) = avg_density(i) + real(occupation(i), dp)
            avg_energy_density(i) = avg_energy_density(i) + energy(i)
            
            ! !!! FIX: Moment 2: <x^2> for Error calculation
            avg_density_sq(i) = avg_density_sq(i) + (real(occupation(i), dp)**2)
            avg_energy_density_sq(i) = avg_energy_density_sq(i) + (energy(i)**2)
        end do
    end subroutine take_measurement
    
    subroutine calculate_final_stats()
        implicit none
        real(dp) :: mean_rho, mean_e, var_rho, var_e, n_samp_real
        
        if (sample_count > 0) then
            n_samp_real = real(sample_count, dp)
            
            do i = 1, L
                ! 1. Calculate Means <x>
                mean_rho = avg_density(i) / n_samp_real
                mean_e = avg_energy_density(i) / n_samp_real
                
                ! 2. Calculate Standard Error = sqrt(Variance / N_samples)
                ! Variance = <x^2> - <x>^2
                
                ! Density Error
                var_rho = (avg_density_sq(i) / n_samp_real) - (mean_rho**2)
                if (var_rho < 0.0_dp) var_rho = 0.0_dp ! Numerical safety
                error_density(i) = sqrt(var_rho / n_samp_real)
                
                ! Energy Error
                var_e = (avg_energy_density_sq(i) / n_samp_real) - (mean_e**2)
                if (var_e < 0.0_dp) var_e = 0.0_dp
                error_energy(i) = sqrt(var_e / n_samp_real)
                
                ! Store final averages in the arrays for printing
                avg_density(i) = mean_rho
                avg_energy_density(i) = mean_e
                
                ! Calculate Temperature: e / rho
                if (mean_rho > 1.0e-10_dp) then
                    avg_temperature(i) = mean_e / mean_rho
                else
                    avg_temperature(i) = 0.0_dp
                end if
            end do
        end if
    end subroutine calculate_final_stats

    subroutine write_intermediate_results()
        implicit none
        integer :: i
        real(dp) :: mean_rho, mean_e, var_rho, var_e, n_samp_real
        real(dp) :: temp_density(L), temp_energy_density(L), temp_temperature(L)
        real(dp) :: temp_error_density(L), temp_error_energy(L)
        
        if (sample_count > 0) then
            n_samp_real = real(sample_count, dp)
            
            ! Calculate temporal averages and errors
            do i = 1, L
                ! 1. Calculate Means <x>
                mean_rho = avg_density(i) / n_samp_real
                mean_e = avg_energy_density(i) / n_samp_real
                
                ! 2. Calculate Standard Error = sqrt(Variance / N_samples)
                ! Variance = <x^2> - <x>^2
                
                ! Density Error
                var_rho = (avg_density_sq(i) / n_samp_real) - (mean_rho**2)
                if (var_rho < 0.0_dp) var_rho = 0.0_dp ! Numerical safety
                temp_error_density(i) = sqrt(var_rho / n_samp_real)
                
                ! Energy Error
                var_e = (avg_energy_density_sq(i) / n_samp_real) - (mean_e**2)
                if (var_e < 0.0_dp) var_e = 0.0_dp
                temp_error_energy(i) = sqrt(var_e / n_samp_real)
                
                ! Store temporal averages
                temp_density(i) = mean_rho
                temp_energy_density(i) = mean_e
                
                ! Calculate Temperature: e / rho
                if (mean_rho > 1.0e-10_dp) then
                    temp_temperature(i) = mean_e / mean_rho
                else
                    temp_temperature(i) = 0.0_dp
                end if
            end do
            
            ! Write to same files as print_results() (overwrite)
            open(unit=10, file='density_profile.dat', status='replace')
            write(10,*) '# Pos, Avg_Dens, Err_Dens'
            do i = 1, L
                write(10,*) i, temp_density(i), temp_error_density(i)
            end do
            close(10)
            
            open(unit=11, file='energy_density_profile.dat', status='replace')
            write(11,*) '# Pos, Avg_Energy, Err_Energy'
            do i = 1, L
                write(11,*) i, temp_energy_density(i), temp_error_energy(i)
            end do
            close(11)
            
            open(unit=12, file='temperature_profile.dat', status='replace')
            write(12,*) '# Pos, Temperature'
            do i = 1, L
                write(12,*) i, temp_temperature(i)
            end do
            close(12)
            
            write(*,*) 'Intermediate results written at step ', step, ' (samples: ', sample_count, ')'
        end if
    end subroutine write_intermediate_results
    
    subroutine print_results()
        implicit none
        integer :: i
        
        write(*,*) ''
        write(*,*) '========================================'
        write(*,*) 'FINAL RESULTS (With Errors)'
        write(*,*) '========================================'
        write(*,*) 'Samples: ', sample_count
        
        ! Files with Errors columns
        open(unit=10, file='density_profile.dat', status='replace')
        write(10,*) '# Pos, Avg_Dens, Err_Dens'
        do i = 1, L
            write(10,*) i, avg_density(i), error_density(i)
        end do
        close(10)
        
        open(unit=11, file='energy_density_profile.dat', status='replace')
        write(11,*) '# Pos, Avg_Energy, Err_Energy'
        do i = 1, L
            write(11,*) i, avg_energy_density(i), error_energy(i)
        end do
        close(11)
        
        open(unit=12, file='temperature_profile.dat', status='replace')
        write(12,*) '# Pos, Temperature' 
        ! Nota: El error de temperatura es mas complejo (propagacion de errores de cociente)
        ! Por ahora imprimimos el valor medio.
        do i = 1, L
            write(12,*) i, avg_temperature(i)
        end do
        close(12)
        
        write(*,*) 'Files written: density_profile.dat, energy_density_profile.dat, temperature_profile.dat'
    end subroutine print_results
    
    function get_random() result(r)
        implicit none
        real(dp) :: r
        call random_number(r)
    end function get_random

end program kep

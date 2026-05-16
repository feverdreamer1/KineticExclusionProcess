!===============================================================================
! SSEP_SIM_CORRECTED.F90
! Validacion SSEP Correcta (Usando Ghost Sites)
! - Beta = 0.0 (Sin dependencia de energia)
! - Rho_left = 0.8, Rho_right = 0.2 (Debe salir una linea recta)
!===============================================================================

program ssep_sim_corrected
    implicit none
    ! Precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: int64 = selected_int_kind(18)
    
    ! Physical parameters - SSEP LIMIT
    integer, parameter :: L = 200
    real(dp), parameter :: T_left = 5.0_dp
    real(dp), parameter :: T_right = 1.0_dp
    real(dp), parameter :: beta = 0.0_dp  ! CRUCIAL: beta=0 para SSEP
    
    ! Reservorios (Gradiente fuerte)
    real(dp), parameter :: rho_left = 0.8_dp
    real(dp), parameter :: rho_right = 0.2_dp
    
    ! Simulation parameters
    integer(int64), parameter :: n_steps = 400000000000_int64
    integer(int64), parameter :: relaxation_steps = 10000000000_int64
    integer(int64), parameter :: measurement_steps = n_steps - relaxation_steps
    integer(int64), parameter :: sample_interval = 100000_int64
    integer(int64), parameter :: file_output_interval = 1000000000_int64
    
    ! Lattice arrays
    integer :: occupation(0:L+1)
    real(dp) :: energy(0:L+1)
    
    ! Promedios y Errores
    real(dp) :: avg_density(L), avg_density_sq(L), error_density(L)
    
    ! KMC variables
    real(dp) :: transition_rates(0:L) ! Rates de los enlaces (Bond 0 es 0<->1)
    real(dp) :: total_rate
    real(dp) :: current_time
    
    ! Statistics
    integer(int64) :: sample_count
    
    ! Loop variables
    integer :: i, j, bond, iaux, iaux2
    integer(int64) :: step
    real(dp) :: dt, random_val, pair_energy, local_temp
    
    ! Inicializacion
    call initialize_random_seed()
    call initialize_lattice()
    
    ! Reset arrays
    avg_density = 0.0_dp
    avg_density_sq = 0.0_dp
    sample_count = 0
    
    write(*,*) '========================================'
    write(*,*) '       SSEP SIMULATION (CORRECTED)      '
    write(*,*) '========================================'
    write(*,*) 'L = ', L
    write(*,*) 'Rho_Left = ', rho_left, ' Rho_Right = ', rho_right
    write(*,*) 'Beta = ', beta
    write(*,*) '========================================'
    
    ! Run simulation
    do step = 1, n_steps
        
        call perform_kmc_step()
        
        ! Medidas
        if (step > relaxation_steps) then
            if (mod(step - relaxation_steps, sample_interval) == 0) then
                call take_measurement()
            end if
        end if
        
        ! Output intermedio
        if (step > relaxation_steps) then
            if (mod(step - relaxation_steps, file_output_interval) == 0) then
               call print_results() ! Escribimos sobre la marcha
               write(*,*) 'Saving at step:', step
            end if
        end if
        
        ! Progress
        if (mod(step, n_steps/10_int64) == 0) then
            write(*,*) 'Progress: ', int(100.0_dp * step / n_steps), '%'
        end if
    end do
    
    call calculate_final_stats()
    call print_results()
    
    write(*,*) 'DONE.'
    
contains

    subroutine initialize_random_seed()
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
        real(dp) :: p_init
        
        occupation = 0
        energy = 0.0_dp
        
        ! Rampa lineal inicial para acelerar convergencia
        do i = 1, L
            p_init = rho_left + (rho_right - rho_left) * real(i-1,dp)/real(L-1,dp)
            if (get_random() < p_init) then
                occupation(i) = 1
                ! En SSEP la energia no importa para el movimiento, pero la inicializamos
                energy(i) = -T_left * log(max(get_random(), 1.0e-10_dp)) 
            end if
        end do
        
        call update_boundaries()
    end subroutine initialize_lattice
    
    subroutine perform_kmc_step()
        implicit none
        
        ! 1. Actualizar SITIOS FANTASMA (0 y L+1)
        ! Esto simula el reservorio infinito.
        call update_boundaries()
        
        ! 2. Calcular rates de salto
        call calculate_transition_rates()
        
        if (total_rate <= 1.0e-10_dp) return
        
        ! 3. Tiempo
        dt = -log(max(get_random(), 1.0e-10_dp)) / total_rate
        current_time = current_time + dt
        
        ! 4. Elegir enlace
        random_val = get_random() * total_rate
        bond = 0
        do i = 0, L
            if (transition_rates(i) > random_val) then
                bond = i
                exit
            end if
        end do
        
        i = bond      ! Sitio Izq del enlace
        j = bond + 1  ! Sitio Der del enlace
        
        ! 5. Ejecutar movimiento
        iaux = occupation(i) + occupation(j)
        
        if (iaux == 2) then
            ! Colision (Solo intercambian energia, no se mueven)
            ! En SSEP puro beta=0 esto no afecta a la densidad, pero ocurre.
            pair_energy = energy(i) + energy(j)
            energy(i) = get_random() * pair_energy
            energy(j) = pair_energy - energy(i)
            
        else if (iaux == 1) then
            ! Salto (Particle-Hole)
            ! Esto incluye saltos desde/hacia los reservorios (0->1 o L->L+1)
            pair_energy = energy(i) + energy(j)
            iaux2 = occupation(i)
            occupation(i) = occupation(j)
            occupation(j) = iaux2
            energy(i) = occupation(i) * pair_energy
            energy(j) = occupation(j) * pair_energy
        end if
        
        ! Nota: No hace falta una rutina especial de inyeccion.
        ! Si bond=0 y occupation(0)=1 y occupation(1)=0, la particula salta de 0 a 1.
        ! Eso ES una inyeccion.
        
    end subroutine perform_kmc_step
    
    subroutine calculate_transition_rates()
        implicit none
        
        ! En SSEP (beta=0), rate = (Ei+Ej)^0 = 1.0
        ! Pero SOLO si al menos uno de los sitios tiene particula (si no, no hay evento).
        ! Sin embargo, KMP/KEP define rate = (Ei+Ej)^beta. Si ambos vacios, E=0, Rate=0.
        
        ! Bond 0 (Reservorio Izq <-> Sitio 1)
        transition_rates(0) = (energy(0) + energy(1))**beta
        total_rate = transition_rates(0)
        
        ! Bulk (1 <-> L)
        do i = 1, L-1
            transition_rates(i) = total_rate + (energy(i) + energy(i+1))**beta
            total_rate = transition_rates(i)
        end do

        ! Bond L (Sitio L <-> Reservorio Der)
        transition_rates(L) = total_rate + (energy(L) + energy(L+1))**beta
        total_rate = transition_rates(L)
        
        ! Nota importante:
        ! Aunque beta=0 haga que (...)**0 sea 1,
        ! si occupation(i)+occupation(j) == 0 (ambos vacios), energy es 0.
        ! 0^0 suele dar 1 en fortran, PERO logicamente no puede haber evento si no hay particulas.
        ! El truco es que si ambos son 0, no pasa nada en perform_kmc_step (iaux=0 no hace nada).
        ! Asi que podemos dejar el rate como 1, gastaremos CPU en elegir un evento nulo,
        ! o podemos optimizar poniendo rate=0 si ambos vacios.
        ! Para ser fieles al KEP general, lo dejamos con la formula de energia.
        ! Si energy=0 y beta=0, fortran suele dar 1. Si esto da problemas, añade un check.
        
    end subroutine calculate_transition_rates
    
    subroutine update_boundaries()
        implicit none
        
        ! Reservorio Izquierdo (Sitio Fantasma 0)
        if (get_random() < rho_left) then
            occupation(0) = 1
            energy(0) = -T_left * log(max(get_random(), 1.0e-10_dp))
        else
            occupation(0) = 0
            energy(0) = 0.0_dp
        end if
        
        ! Reservorio Derecho (Sitio Fantasma L+1)
        if (get_random() < rho_right) then
            occupation(L+1) = 1
            energy(L+1) = -T_right * log(max(get_random(), 1.0e-10_dp))
        else
            occupation(L+1) = 0
            energy(L+1) = 0.0_dp
        end if
        
    end subroutine update_boundaries
    
    subroutine take_measurement()
        implicit none
        real(dp) :: dens
        
        sample_count = sample_count + 1
        
        do i = 1, L
            dens = real(occupation(i), dp)
            avg_density(i) = avg_density(i) + dens
            avg_density_sq(i) = avg_density_sq(i) + (dens**2)
        end do
    end subroutine take_measurement
    
    subroutine calculate_final_stats()
        implicit none
        real(dp) :: mean_rho, var_rho, ns
        
        if (sample_count > 0) then
            ns = real(sample_count, dp)
            do i = 1, L
                mean_rho = avg_density(i) / ns
                var_rho = (avg_density_sq(i) / ns) - (mean_rho**2)
                if (var_rho < 0.0_dp) var_rho = 0.0_dp
                
                avg_density(i) = mean_rho
                error_density(i) = sqrt(var_rho / ns)
            end do
        end if
    end subroutine calculate_final_stats
    
    subroutine print_results()
        implicit none
        integer :: i
        open(unit=10, file='density_profile.dat', status='replace')
        write(10,*) '# Pos, Density, Error'
        do i = 1, L
            write(10,*) i, avg_density(i), error_density(i)
        end do
        close(10)
    end subroutine print_results
    
    function get_random() result(r)
        implicit none
        real(dp) :: r
        call random_number(r)
    end function get_random

end program ssep_sim_corrected
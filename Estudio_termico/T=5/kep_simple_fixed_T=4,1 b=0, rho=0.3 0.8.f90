program kep
    implicit none
    ! Precision
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: int64 = selected_int_kind(18)
    
    !---------------------------------------------------------------------------
    ! 1. PARAMETROS FISICOS (Ajustar segun necesidad)
    !---------------------------------------------------------------------------
    integer, parameter :: L = 200            ! Tamaño red
    
    ! Temperaturas
    real(dp), parameter :: T_left = 5.0_dp
    real(dp), parameter :: T_right = 1.0_dp
    
    ! Densidades de los reservorios
    real(dp), parameter :: rho_left = 0.5_dp
    real(dp), parameter :: rho_right = 0.5_dp
    
    ! Parametro Beta (0.0 = SSEP, 1.0 = KEP estandar)
    real(dp), parameter :: beta = 0.5_dp
    
    !---------------------------------------------------------------------------
    ! 2. PARAMETROS DE SIMULACION
    !---------------------------------------------------------------------------
    integer(int64), parameter :: n_steps = 400000000000_int64
    integer(int64), parameter :: relaxation_steps = 10000000000_int64
    
    ! Frecuencia de toma de muestras (para estadistica)
    integer(int64), parameter :: sample_interval = 100000_int64
    
    ! Frecuencia de escritura en disco (Output intermedio)
    ! Se escribiran los ficheros .dat cada vez que pasen estos pasos
    integer(int64), parameter :: file_output_interval = 1000000000_int64
    
    !---------------------------------------------------------------------------
    ! VARIABLES
    !---------------------------------------------------------------------------
    ! Arrays del sistema
    integer :: occupation(0:L+1)
    real(dp) :: energy(0:L+1)
    real(dp) :: temperature(0:L+1)
    
    ! Acumuladores para Promedios <x>
    real(dp) :: avg_density(L)
    real(dp) :: avg_energy_density(L)
    real(dp) :: avg_temperature(L)
    
    ! Acumuladores para Cuadrados <x^2> (Para el error)
    real(dp) :: avg_density_sq(L)
    real(dp) :: avg_energy_density_sq(L)
    
    ! Variables KMC
    real(dp) :: transition_rates(0:L)
    real(dp) :: total_rate
    real(dp) :: current_time
    
    ! Estadisticas
    integer(int64) :: sample_count
    
    ! Auxiliares
    integer :: i, j, bond, iaux, iaux2
    integer(int64) :: step
    real(dp) :: dt, random_val, pair_energy
    
    !---------------------------------------------------------------------------
    ! PROGRAMA PRINCIPAL
    !---------------------------------------------------------------------------
    call initialize_random_seed()
    call initialize_lattice()
    
    ! Limpiar acumuladores
    avg_density = 0.0_dp
    avg_energy_density = 0.0_dp
    avg_temperature = 0.0_dp
    avg_density_sq = 0.0_dp
    avg_energy_density_sq = 0.0_dp
    
    sample_count = 0
    
    write(*,*) '========================================'
    write(*,*) ' KEP SIMULATION (FULL VERSION)          '
    write(*,*) '========================================'
    write(*,*) 'L=', L, ' Beta=', beta
    write(*,*) 'Rho_L=', rho_left, ' Rho_R=', rho_right
    write(*,*) 'T_L=', T_left, ' T_R=', T_right
    write(*,*) 'Steps=', n_steps
    write(*,*) 'Output Interval=', file_output_interval
    write(*,*) '========================================'
    
    ! Bucle Principal
    do step = 1, n_steps
        
        call perform_kmc_step()
        
        if (step > relaxation_steps) then
            
            ! 1. Tomar medidas estadisticas
            if (mod(step - relaxation_steps, sample_interval) == 0) then
                call take_measurement()
            end if
            
            ! 2. Escribir resultados intermedios en disco
            if (mod(step - relaxation_steps, file_output_interval) == 0) then
                call write_intermediate_results()
                write(*,*) 'Auto-save at step: ', step
            end if
            
        end if
        
        ! Progreso en pantalla
        if (mod(step, n_steps/10_int64) == 0) then
            if (step <= relaxation_steps) then
                 write(*,*) 'Relaxing... ', int(100.0_dp * step / relaxation_steps), '%'
            else
                 write(*,*) 'Measuring... ', int(100.0_dp * (step-relaxation_steps) / (n_steps-relaxation_steps)), '%'
            end if
        end if
    end do
    
    ! Escritura final
    call write_intermediate_results()
    
    write(*,*) '========================================'
    write(*,*) 'SIMULATION FINISHED SUCCESSFULLY'
    write(*,*) '========================================'
    
contains

    !---------------------------------------------------------------------------
    ! SUBRUTINAS
    !---------------------------------------------------------------------------

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
        real(dp) :: p_target, t_target
        
        occupation = 0
        energy = 0.0_dp
        temperature = 0.0_dp
        
        ! Inicializacion suave (rampa lineal)
        do i = 1, L
            p_target = rho_left + (rho_right - rho_left) * real(i-1,dp)/real(L-1,dp)
            t_target = T_left + (T_right - T_left) * real(i-1,dp)/real(L-1,dp)
            
            if (get_random() < p_target) then
                occupation(i) = 1
                energy(i) = -t_target * log(max(get_random(), 1.0e-10_dp))
            end if
        end do

        call update_boundaries()
    end subroutine initialize_lattice
    
    subroutine perform_kmc_step()
        implicit none
        
        ! 1. Actualizar bordes (Reservorios Estocasticos)
        call update_boundaries()
        
        ! 2. Calcular rates
        ! Enlace 0 (Pared L <-> 1)
        transition_rates(0) = (energy(0) + energy(1))**beta
        total_rate = transition_rates(0)
        
        ! Bulk
        do i = 1, L-1
            transition_rates(i) = total_rate + (energy(i) + energy(i+1))**beta
            total_rate = transition_rates(i)
        end do

        ! Enlace L (L <-> Pared R)
        transition_rates(L) = total_rate + (energy(L) + energy(L+1))**beta
        total_rate = transition_rates(L)
        
        if (total_rate <= 1.0e-12_dp) return
        
        ! 3. Tiempo
        dt = -log(max(get_random(), 1.0e-12_dp)) / total_rate
        current_time = current_time + dt
        
        ! 4. Seleccionar evento
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
        
        ! 5. Ejecutar
        iaux = occupation(i) + occupation(j)
        
        if (iaux == 2) then
            ! Colision
            pair_energy = energy(i) + energy(j)
            energy(i) = get_random() * pair_energy
            energy(j) = pair_energy - energy(i)
        else if (iaux == 1) then
            ! Salto
            pair_energy = energy(i) + energy(j)
            iaux2 = occupation(i)
            occupation(i) = occupation(j)
            occupation(j) = iaux2
            energy(i) = occupation(i) * pair_energy
            energy(j) = occupation(j) * pair_energy
        end if
        
    end subroutine perform_kmc_step
    
    subroutine update_boundaries()
        implicit none
        ! Logica de Reservorios
        ! Pared Izquierda
        if (get_random() < rho_left) then
            occupation(0) = 1
            energy(0) = -T_left * log(max(get_random(), 1.0e-10_dp))
        else
            occupation(0) = 0
            energy(0) = 0.0_dp
        end if
        
        ! Pared Derecha
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
        real(dp) :: dens, ener
        
        sample_count = sample_count + 1
        
        do i = 1, L
            dens = real(occupation(i), dp)
            ener = energy(i)
            
            ! Acumular suma <x>
            avg_density(i) = avg_density(i) + dens
            avg_energy_density(i) = avg_energy_density(i) + ener
            
            ! Acumular suma cuadrados <x^2>
            avg_density_sq(i) = avg_density_sq(i) + (dens**2)
            avg_energy_density_sq(i) = avg_energy_density_sq(i) + (ener**2)
        end do
    end subroutine take_measurement

    !---------------------------------------------------------------------------
    ! ESCRITURA DE RESULTADOS (Intermedio y Final)
    !---------------------------------------------------------------------------
    subroutine write_intermediate_results()
        implicit none
        integer :: k
        real(dp) :: mean_rho, mean_e, var_rho, var_e, err_rho, err_e, temp_val
        real(dp) :: n_samp_real
        
        if (sample_count == 0) return
        
        n_samp_real = real(sample_count, dp)
        
        ! Abrimos ficheros (modo 'replace' sobrescribe el anterior)
        
        ! 1. Densidad
        open(unit=10, file='density_profile.dat', status='replace')
        write(10,*) '# Pos, Density, Err_Density'
        
        ! 2. Energia
        open(unit=11, file='energy_density_profile.dat', status='replace')
        write(11,*) '# Pos, Energy, Err_Energy'
        
        ! 3. Temperatura
        open(unit=12, file='temperature_profile.dat', status='replace')
        write(12,*) '# Pos, Temperature'
        
        do k = 1, L
            ! Calculamos medias y errores TEMPORALES (sin modificar acumuladores)
            mean_rho = avg_density(k) / n_samp_real
            mean_e   = avg_energy_density(k) / n_samp_real
            
            ! Varianza y Error Std
            var_rho = (avg_density_sq(k) / n_samp_real) - (mean_rho**2)
            var_e   = (avg_energy_density_sq(k) / n_samp_real) - (mean_e**2)
            
            if (var_rho < 0.0_dp) var_rho = 0.0_dp
            if (var_e < 0.0_dp) var_e = 0.0_dp
            
            err_rho = sqrt(var_rho / n_samp_real)
            err_e   = sqrt(var_e / n_samp_real)
            
            ! Temperatura local (e / rho)
            if (mean_rho > 1.0e-10_dp) then
                temp_val = mean_e / mean_rho
            else
                temp_val = 0.0_dp
            end if
            
            ! Escribir lineas
            write(10,*) k, mean_rho, err_rho
            write(11,*) k, mean_e, err_e
            write(12,*) k, temp_val
        end do
        
        close(10)
        close(11)
        close(12)
        
    end subroutine write_intermediate_results
    
    function get_random() result(r)
        implicit none
        real(dp) :: r
        call random_number(r)
    end function get_random

end program kep

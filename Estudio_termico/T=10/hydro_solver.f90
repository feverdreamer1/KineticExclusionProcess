program hydro_solver
  implicit none
  ! Solución del sistema estacionario usando FD + Newton (Jacobiano numérico)
  ! Unknowns: rho(1..N-1), T(1..N-1), c1, c2  --> total 2*(N-1)+2 variables

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer :: N, i, iter, maxiter, neq
  real(dp) :: beta, dx, tol, fd_eps, damping
  real(dp), allocatable :: x(:), rho(:), T(:)
  real(dp), allocatable :: U(:), F(:), dU(:,:)
  real(dp) :: rhoL, rhoR, TL, TR
  real(dp) :: normF
  logical :: converged

  ! Linear solver workspace
  real(dp), allocatable :: A(:,:), bvec(:)
  integer :: info

  ! -------------------------
  ! --- PARAMETROS (ajusta)
  ! -------------------------
  N = 400                ! numero de subintervalos (N+1 nodos); mayor N -> mayor resolucion
  beta = 0.5_dp     ! parametro beta del modelo
  rhoL = 0.5_dp
  rhoR = 0.5_dp
  TL = 10.0_dp
  TR = 1.0_dp

  maxiter = 60
  tol = 1.0e-9_dp
  fd_eps = 1.0e-8_dp
  damping = 0.5_dp   ! Damping factor for Newton iterations

  ! -------------------------
  dx = 1.0_dp / real(N,dp)
  allocate(x(0:N))
  do i = 0, N
     x(i) = real(i,dp) * dx
  end do

  ! allocate unknowns: interior nodes only
  neq = 2*(N-1) + 2    ! rho_1..rho_{N-1}, T_1..T_{N-1}, c1, c2
  allocate(U(neq), F(neq), dU(neq,neq), A(neq,neq), bvec(neq))

  ! Initial guess: linear interpolation for rho and T; c1,c2 initial 0
  do i = 1, N-1
     U(i) = rhoL + (rhoR - rhoL) * real(i,dp)/real(N,dp)        ! rho_i
     U(N-1 + i) = TL + (TR - TL) * real(i,dp)/real(N,dp)       ! T_i
  end do
  U(neq-1) = 0.0_dp    ! c1
  U(neq)   = 0.0_dp    ! c2

  converged = .false.

  do iter = 1, maxiter
     call residual(U, F, N, dx, beta, rhoL, rhoR, TL, TR)
     normF = vec_norm(F)
     write(*,'(A,I3,A,ES12.4)') 'Iter = ', iter, '  ||F|| = ', normF
     if (normF < tol) then
        converged = .true.
        exit
     end if

     ! Build Jacobian by forward finite differences: dU(:,j) columns
     call jacobian_fd(U, dU, N, dx, beta, rhoL, rhoR, TL, TR, fd_eps)

     ! Solve linear system dU * delta = -F  --> we reuse A and bvec
     A = dU
     bvec = -F

     call gauss_solve(A, bvec, neq, info)
     if (info /= 0) then
        write(*,*) 'Warning: gauss_solve returned info=', info
        write(*,*) 'Matrix appears singular at iteration ', iter
        ! Try with smaller damping
        damping = damping * 0.5_dp
        if (damping < 1.0e-6_dp) then
           write(*,*) 'Damping too small, stopping iteration'
           exit
        end if
        ! Apply update with reduced damping anyway
        U = U + damping * bvec
     else
        ! update unknowns with damping
        U = U + damping * bvec
     end if

  end do

  if (.not. converged) then
     write(*,*) 'No convergió en ', maxiter, ' iteraciones. ||F||=', normF
     write(*,*) 'Guardando solución no convergida de todos modos...'
  else
     write(*,*) 'Convergido en ', iter, ' iteraciones. ||F||=', normF
  end if

  ! Reconstruir perfiles completos
  allocate(rho(0:N), T(0:N))
  rho(0) = rhoL
  rho(N) = rhoR
  T(0) = TL
  T(N) = TR
  do i = 1, N-1
     rho(i) = U(i)
     T(i)   = U(N-1 + i)
  end do

  ! imprimir resultados en formato CSV
  open(unit=10, file='solution.dat', status='replace')
  write(10,'(A)') '# x, rho, T'
  write(10,'(A)') '# CSV format for plotting'
  write(10,'(A)') '# x,rho,T'
  do i = 0, N
     write(10,'(E16.8,A,E16.8,A,E16.8)') x(i), ',', rho(i), ',', T(i)
  end do
  close(10)
  write(*,*) 'Perfiles guardados en solution.dat (formato CSV)'
  
  ! Also create a .csv file for convenience
  open(unit=11, file='solution.csv', status='replace')
  write(11,'(A)') 'x,rho,T'
  do i = 0, N
     write(11,'(E16.8,A,E16.8,A,E16.8)') x(i), ',', rho(i), ',', T(i)
  end do
  close(11)
  write(*,*) 'Perfiles guardados tambien en solution.csv'
  write(*,*) 'Corrientes (c1,c2) = ', U(neq-1), U(neq)

contains

  subroutine residual(U, F, N, dx, beta, rhoL, rhoR, TL, TR)
    implicit none
    integer, intent(in) :: N
    real(dp), intent(in) :: dx, beta, rhoL, rhoR, TL, TR
    real(dp), intent(in) :: U(:)
    real(dp), intent(out) :: F(:)
    integer :: i, neq, idx_rho, idx_T
    real(dp) :: c1, c2
    real(dp) :: rho_im1, rho_i, rho_ip1, T_im1, T_i, T_ip1
    real(dp) :: D11, D12, D21, D22, denom
    real(dp) :: d_rho, d_T

    neq = size(U)
    ! indexing: rho_i (i=1..N-1) -> U(1..N-1)
    !           T_i   (i=1..N-1) -> U(N .. 2*(N-1))
    !           c1 -> U(neq-1), c2 -> U(neq)
    c1 = U(neq-1)
    c2 = U(neq)

    ! For interior nodes i=1..N-1: two equations per node
    do i = 1, N-1
       idx_rho = i
       idx_T   = (N-1) + i

       rho_i = U(idx_rho)
       T_i   = U(idx_T)

       if (i == 1) then
          rho_im1 = rhoL
          T_im1   = TL
       else
          rho_im1 = U(idx_rho - 1)
          T_im1   = U(idx_T - 1)
       end if
       if (i == N-1) then
          rho_ip1 = rhoR
          T_ip1   = TR
       else
          rho_ip1 = U(idx_rho + 1)
          T_ip1   = U(idx_T + 1)
       end if

       d_rho = (rho_ip1 - rho_im1) / (2.0_dp * dx)
       d_T   = (T_ip1   - T_im1)   / (2.0_dp * dx)

       call Dcoef(rho_i, T_i, beta, D11, D12, D21, D22)
       ! Equations: D11*d_rho + D12*d_T - c1 = 0
       F(2*(i-1)+1) = D11 * d_rho + D12 * d_T - c1
       !           D21*d_rho + D22*d_T - c2 = 0
       F(2*(i-1)+2) = D21 * d_rho + D22 * d_T - c2
    end do

    ! Two extra equations: enforce boundary consistency for c1,c2 by integrating definition
    ! Alternative approach: enforce global integral constraint -> here we add two equations that
    ! ensure that the discrete version at the left node reproduces c1,c2 computed with left derivatives:
    ! Use left-sided derivative at x=0 to define c1,c2 consistency:
    ! compute d_rho0 = (U(1)-rhoL)/(dx) approx forward; d_T0 = (U(N)-TL)/dx
    ! Note: U(1) is rho_1 (first interior point), U(N) is T_1 (first interior temperature)
    d_rho = (U(1) - rhoL) / dx
    d_T   = (U(N) - TL)   / dx
    call Dcoef(rhoL, TL, beta, D11, D12, D21, D22)
    F(neq-1) = D11 * d_rho + D12 * d_T - c1
    F(neq)   = D21 * d_rho + D22 * d_T - c2

  end subroutine residual

  subroutine Dcoef(rho, T, beta, D11, D12, D21, D22)
    implicit none
    real(dp), intent(in) :: rho, T, beta
    real(dp), intent(out) :: D11, D12, D21, D22
    real(dp) :: g1, g2

    ! Using expressions from el TFM:
    ! D11 = Gamma(beta+1) * T^beta
    ! D12 = beta * Gamma(beta+1) * rho*(1-rho) * T^(beta-1)
    ! D21 = Gamma(beta+2) * T^(beta+1)
    ! D22 = Gamma(beta+2) * T^beta * rho * ( (1+beta) + rho*(beta^2 - 7*beta -6)/12 )
    ! We use intrinsic gamma function (if your compiler lacks it, replace accordingly)
    g1 = gamma(beta + 1.0_dp)
    g2 = gamma(beta + 2.0_dp)

    D11 = g1 * T**beta
    if (T <= 0.0_dp) then
       D12 = 0.0_dp
       D21 = 0.0_dp
       D22 = 0.0_dp
    else
       D12 = beta * g1 * rho*(1.0_dp - rho) * T**(beta - 1.0_dp)
       D21 = g2 * T**(beta + 1.0_dp)
       D22 = g2 * T**beta * rho * ( (1.0_dp + beta) + rho * ((beta**2 - 7.0_dp*beta - 6.0_dp)/12.0_dp) )
    end if
  end subroutine Dcoef

  subroutine jacobian_fd(U, J, N, dx, beta, rhoL, rhoR, TL, TR, eps)
    implicit none
    real(dp), intent(in) :: U(:)
    real(dp), intent(out) :: J(:,:)
    integer, intent(in) :: N
    real(dp), intent(in) :: dx, beta, rhoL, rhoR, TL, TR, eps
    integer :: neq, jcol, i
    real(dp), allocatable :: F0(:), F1(:), Up(:)
    neq = size(U)
    allocate(F0(neq), F1(neq), Up(neq))
    call residual(U, F0, N, dx, beta, rhoL, rhoR, TL, TR)
    J = 0.0_dp
    do jcol = 1, neq
       Up = U
       Up(jcol) = U(jcol) + eps
       call residual(Up, F1, N, dx, beta, rhoL, rhoR, TL, TR)
       J(:,jcol) = (F1 - F0) / eps
    end do
    deallocate(F0, F1, Up)
  end subroutine jacobian_fd

  function vec_norm(v) result(nrm)
    implicit none
    real(dp), intent(in) :: v(:)
    real(dp) :: nrm
    integer :: i
    nrm = 0.0_dp
    do i = 1, size(v)
       nrm = nrm + v(i)**2
    end do
    nrm = sqrt(nrm)
  end function vec_norm

  subroutine gauss_solve(A, b, n, info)
    ! Simple Gaussian elimination with partial pivoting solving A*x=b
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    real(dp), intent(inout) :: b(n)
    integer, intent(out) :: info
    integer :: i,j,k,imax
    real(dp) :: maxa, tmp, factor, piv

    info = 0
    do k = 1, n-1
       ! pivot
       imax = k
       maxa = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) > maxa) then
             maxa = abs(A(i,k))
             imax = i
          end if
       end do
       if (maxa < 1.0e-16_dp) then
          info = k
          return
       end if
      if (imax /= k) then
         ! swap rows k and imax
         do j = 1, n
            tmp = A(k,j)
            A(k,j) = A(imax,j)
            A(imax,j) = tmp
         end do
         tmp = b(k); b(k) = b(imax); b(imax) = tmp
      end if
       piv = A(k,k)
       do i = k+1, n
          factor = A(i,k)/piv
          A(i,k) = 0.0_dp
          do j = k+1, n
             A(i,j) = A(i,j) - factor * A(k,j)
          end do
          b(i) = b(i) - factor * b(k)
       end do
    end do
    if (abs(A(n,n)) < 1.0e-16_dp) then
       info = n
       return
    end if
    ! Back substitution
    do i = n, 1, -1
       tmp = b(i)
       if (i < n) then
          do j = i+1, n
             tmp = tmp - A(i,j) * b(j)
          end do
       end if
       b(i) = tmp / A(i,i)
    end do
  end subroutine gauss_solve

end program hydro_solver
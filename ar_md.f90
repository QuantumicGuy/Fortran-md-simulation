program ar_md
    implicit none
    integer, parameter :: n = 1000
    integer, parameter :: steps = 100
    real*8, parameter :: dt = 0.001
    real*8, parameter :: box_size = 10.0
    real*8, parameter :: sigma = 3.405d0
    real*8, parameter :: epsilon = 0.2381d0
    real*8 :: r(3,n), v(3,n), f(3,n)
    integer :: i, j, k, step

    call initialize_positions(r, box_size, n)
    call initialize_velocities(v, n)

    do step = 1, steps
        call calculate_forces(r, f, n, box_size, sigma, epsilon)
        call update_positions_velocities(r, v, f, dt, n, box_size)
        if (mod(step, 10) == 0) then
            call print_energy(r, v, n, box_size, sigma, epsilon)
        end if
    end do

contains

    subroutine initialize_positions(r, box_size, n)
        real*8, intent(out) :: r(3,n)
        real*8, intent(in) :: box_size
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            r(1,i) = box_size * rand()
            r(2,i) = box_size * rand()
            r(3,i) = box_size * rand()
        end do
    end subroutine initialize_positions

    subroutine initialize_velocities(v, n)
        real*8, intent(out) :: v(3,n)
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            v(1,i) = 0.1 * (2.0 * rand() - 1.0)
            v(2,i) = 0.1 * (2.0 * rand() - 1.0)
            v(3,i) = 0.1 * (2.0 * rand() - 1.0)
        end do
    end subroutine initialize_velocities

    subroutine calculate_forces(r, f, n, box_size, sigma, epsilon)
        real*8, intent(in) :: r(3,n)
        real*8, intent(out) :: f(3,n)
        real*8, intent(in) :: box_size, sigma, epsilon
        real*8 :: dr(3), r2, r6, r12, force_factor
        integer :: i, j
        f = 0.0d0
        do i = 1, n - 1
            do j = i + 1, n
                dr = r(:,i) - r(:,j)
                dr = dr - box_size * anint(dr / box_size)
                r2 = dot_product(dr, dr)
                r6 = r2 * r2 * r2
                r12 = r6 * r6
                force_factor = 48.0d0 * epsilon * (r12 - 0.5d0 * r6) / (r2 * r6)
                f(:,i) = f(:,i) + force_factor * dr
                f(:,j) = f(:,j) - force_factor * dr
            end do
        end do
    end subroutine calculate_forces

    subroutine update_positions_velocities(r, v, f, dt, n, box_size)
        real*8, intent(inout) :: r(3,n), v(3,n)
        real*8, intent(in) :: f(3,n), dt, box_size
        integer, intent(in) :: n
        integer :: i
        do i = 1, n
            v(:,i) = v(:,i) + 0.5d0 * dt * f(:,i)
            r(:,i) = r(:,i) + dt * v(:,i)
            r(:,i) = r(:,i) - box_size * anint(r(:,i) / box_size)
            v(:,i) = v(:,i) + 0.5d0 * dt * f(:,i)
        end do
    end subroutine update_positions_velocities

    subroutine print_energy(r, v, n, box_size, sigma, epsilon)
        real*8, intent(in) :: r(3,n), v(3,n)
        real*8, intent(in) :: box_size, sigma, epsilon
        real*8 :: kinetic_energy, potential_energy, total_energy
        real*8 :: dr(3), r2, r6, r12
        integer :: i, j
        kinetic_energy = 0.0d0
        do i = 1, n
            kinetic_energy = kinetic_energy + 0.5d0 * dot_product(v(:,i), v(:,i))
        end do
        potential_energy = 0.0d0
        do i = 1, n - 1
            do j = i + 1, n
                dr = r(:,i) - r(:,j)
                dr = dr - box_size * anint(dr / box_size)
                r2 = dot_product(dr, dr)
                r6 = r2 * r2 * r2
                r12 = r6 * r6
                potential_energy = potential_energy + 4.0d0 * epsilon * (r12 - r6)
            end do
        end do
        total_energy = kinetic_energy + potential_energy
        print *, 'Step: ', step, ' Total Energy: ', total_energy
    end subroutine print_energy

end program ar_md
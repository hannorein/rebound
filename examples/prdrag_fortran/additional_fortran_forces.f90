module rebound
    use iso_c_binding
    real(c_double), bind(C, name="G") :: G
    integer(c_int), bind(C, name="N") :: N
    integer(c_int), bind(C, name="N_active") :: N_active
    type, bind(c) :: particle
        real(c_double) :: x, y, z, vx, vy, vz, ax, ay, az, m
    end type particle
end module rebound


subroutine tidal_forces(nbod, particles) bind(c, name='tidal_forces')
    use iso_c_binding
    use rebound
    implicit none
    integer(c_int) :: nbod
    type (particle), intent(inout) :: particles(nbod)

    !integer i, dragcoefficient
    !dragcoefficient = -1

    !print *, "This is in Fortran routine..."
    !print *, "x = ", nbod, particles(1)%x, particles(2)%x

    !do i = 2, nbod
        !particles(i)%ax = particles(i)%ax + dragcoefficient*particles(i)%vx
        !particles(i)%ay = particles(i)%ay + dragcoefficient*particles(i)%vy
        !particles(i)%az = particles(i)%az + dragcoefficient*particles(i)%vz
    !end do
end subroutine tidal_forces



subroutine force_radiation(nbod, particles) bind(c, name='force_radiation')
    use iso_c_binding
    use rebound
    implicit none
    integer(c_int) :: nbod
    type (particle), intent(inout) :: particles(nbod)


    real(c_double) :: prx, pry, prz, pr, prvx, prvy, prvz, c, rdot, F_r
    real(c_double) :: betaparticles

    integer i
    betaparticles = 0.01  ! beta parameter
    betaparticles = particles(1)%x

    do i = 3, N
        prx  = particles(i)%x - particles(1)%x
        pry  = particles(i)%y - particles(1)%y
        prz  = particles(i)%z - particles(1)%z
        pr   = sqrt(prx*prx + pry*pry + prz*prz)     ! distance relative to particles(1)
        prvx = particles(i)%vx - particles(1)%vx
        prvy = particles(i)%vx - particles(1)%vy
        prvz = particles(i)%vx - particles(1)%vz
        c      = 1.006491504759635e+04  ! speed of light in unit of G=1, M_sun=1, 1year=1
        rdot   = (prvx*prx + prvy*pry + prvz*prz)/pr  ! radial velocity relative to star
        F_r    = betaparticles*G*particles(1)%m/(pr*pr)
        F_r = 0

        ! Equation (5) of Burns, Lamy, Soter (1979)
        particles(i)%ax = particles(i)%ax + F_r*((1.-rdot/c)*prx/pr - prvx/c)
        particles(i)%ay = particles(i)%ay + F_r*((1.-rdot/c)*pry/pr - prvy/c)
        particles(i)%az = particles(i)%az + F_r*((1.-rdot/c)*prz/pr - prvz/c)
    end do
end subroutine force_radiation


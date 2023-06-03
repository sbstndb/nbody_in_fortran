module mod_compute
!---------------------------
        use iso_fortran_env, only : REAL64, REAL32
        use iso_fortran_env, only : INT32, INT64
        implicit none

        !! declare kinds
        !integer, parameter :: real_k = REAL64
        integer, parameter :: real_k = REAL32
        !integer, parameter :: int_k = INT64
        integer, parameter :: int_k = INT32

!----------------------------

        type particle_s
                real(real_k) :: x , y , z
                real(real_k) :: vx, vy, vz
        end type particle_s

contains
!-----------------------------

subroutine init(array , len)
        implicit none
        type(particle_s), allocatable, dimension(:), intent(inout) :: array
        integer(int_k), intent(in) :: len
                
        integer(int_k) :: i, j

        do i = 1 , len
                array(i)%x = 0.0
                array(i)%y = 0.0
                array(i)%z = 0.0

                array(i)%vx = 0.0
                array(i)%vy = 0.0
                array(i)%vz = 0.0
        end do
end subroutine init

subroutine move(part, n)
        implicit none
        type(particle_s), allocatable, dimension(:), intent(inout) :: part
        integer(int_k) :: n

        integer(int_k) :: i, j, k
        real(real_k), parameter :: softening = 1e-20
        real(real_k), parameter :: dt = 0.0001
        real(real_k) :: fx_t, fy_t, fz_t
        real(real_k) :: dx_t, dy_t, dz_t
        real(real_k) :: d_2, d_3_over_2

        do i = 1 , n
                fx_t = 0.0
                fy_t = 0.0
                fz_t = 0.0
                do j = 1 , n
                        dx_t = part(j)%x - part(i)%x
                        dy_t = part(j)%y - part(i)%y
                        dz_t = part(j)%z - part(i)%z
                        d_2 = dx_t*dx_t + dy_t*dy_t + dz_t*dz_t + softening
                        d_3_over_2 = d_2**(3.0/2.0)

                        fx_t = fx_t + dx_t/d_3_over_2
                        fy_t = fy_t + dy_t/d_3_over_2
                        fz_t = fz_t + dz_t/d_3_over_2
                end do
                part(i)%vx = part(i)%vx + dt*fx_t
                part(i)%vy = part(i)%vy + dt*fy_t
                part(i)%vz = part(i)%vz + dt*fz_t
        end do

        do i = 1 , n
                part(i)%x = part(i)%x + dt*part(i)%vx
                part(i)%y = part(i)%y + dt*part(i)%vy
                part(i)%z = part(i)%z + dt*part(i)%vz
        end do

end subroutine move

subroutine loop(part, n, step)
        implicit none
        type(particle_s), allocatable, dimension(:), intent(inout) :: part
        integer(int_k), intent(in) :: n , step       

        integer(int_k) :: i 

        do i = 1 , step
                print*, " -- iteration numéro " , i
                call move(part, n)
        end do
end subroutine loop        

subroutine show_p1(part)
        implicit none
        type(particle_s), allocatable, dimension(:), intent(inout) :: part

        print*, "La position de la particule 1 est  : " 
        print* , "   x : " , part(1)%x
        print* , "   y : " , part(1)%y
        print* , "   z : " , part(1)%z
end subroutine show_p1

subroutine compute()
        implicit none
        integer(int_k) :: n, step
        type(particle_s), allocatable, dimension(:) :: part

        n = 1000
        step = 100

        allocate(part(n))

        call init(part, n)
        call loop(part, n, step)
        call show_p1(part)

        deallocate(part)
end subroutine compute

!-----------------------------------
end module mod_compute

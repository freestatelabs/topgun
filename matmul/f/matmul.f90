program matmul

    integer :: i, m, n, p, nit, clock_start, clock_stop, clock_rate
    real :: e_time 
    real, allocatable :: A(:), B(:), C(:)

    nit = 1000
    m = 100 
    n = 100 
    p = 100 
    e_time = 0.0

    allocate(A(m*n))
    allocate(B(n*p))
    allocate(C(m*p))

    call random_seed() 
    call random_number(A) 
    call random_number(B)

    call system_clock(count_rate=clock_rate) !Find the time rate

    do i = 1,nit 
        C = 0.0 
        call system_clock(count=clock_start)     !Start Timer
        call mul(C, A, B, m, n, p)
        call system_clock(count=clock_stop)      ! Stop Timer

        e_time = e_time + real(clock_stop-clock_start)/real(clock_rate)
    end do

    e_time = e_time / Nit

    print *, "Time per iteration: ", 1e6*e_time, " us"

    deallocate(A)
    deallocate(B)
    deallocate(C)

end program matmul


subroutine mul(C, A, B, m, n, p) 
    integer, intent(in) :: m, n, p
    real, intent(inout) :: C(m*p) 
    real, intent(in) :: A(m*n), B(n*p)
    integer :: i, j, k

    do i = 1,m 
        do j = 1,p 
            do k = 1,n 
                C((j-1)*p+i) = C((j-1)*p+i) + A((k-1)*n + i) * B((j-1)*p+k)
            end do
        end do
    end do
end subroutine
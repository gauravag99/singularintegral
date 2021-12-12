program singintquad
    implicit none
    integer, parameter :: dp=8
    real(dp) :: lowlim, uplim
    integer :: steps
    complex(dp) :: integral

    lowlim = 0.0_dp
    uplim = 1.0_dp
    steps=200  !keep steps so that interval width ~ 1.0e-2 to 1.0e-3


    CALL integ13(nume,deno,lowlim,uplim,steps,integral)

    write(*,'(A,2ES10.2e2,"i")') "The integral value is", integral

    contains

    subroutine integ13(num,den,q,p,N,integral)
        !!!--------------------------
        !This subroutine performs integration by equidistant intervals for 'nearly' singular functions using quadratic approximations.
        ! func -> function to integrate, q-> lower bound, p-> upper bound, N -> number of intervals (including q and p)
        ! OUTPUT => integral
        ! Use num, den to declare numerator, denominator respectively of the integral A(y)/B(y)
        !!!---------------------------
        real(dp), intent(in) :: q,p
        integer, intent(in) :: N
        real(dp) :: h,x
        COMPLEX(dp) :: num,den, tp(3), int, a,b,c,d,e,f
        COMPLEX(dp), intent(out) :: integral
        integer :: i,k,t
        h=(p-q)/N              ! h is the interval width.
        x=(q+h/1000)           ! x tracks the interval. h/1000 extra factor is added to skew the program a bit from the sing.
        integral=0.
        k =0;t=0                ! dummy indices to keep track of which approximation gets triggered
        write(*,'(A,F7.3, " to", F7.3)') "Integral limits", q,p
        write(*,'(A,ES8.2e2,/,A)') 'Interval width = ', h, repeat("-",10)
        
        ! Main LOOP below
        do i=0,N-1
                !Variables named same as the writeup
                int=0.
                a =(2.) *(den(x+h)-2.*den(x+h/2) +den(x))
                b = -den(x+h) + 4* den(x+h/2) -3*den(x)
                c = den(x)
                d =(2.) *(num(x+h)-2.*num(x+h/2) + num(x))
                e = -num(x+h) - 3.*num(x) + 4.* num(x+h/2)
                f = num(x)

                if (abs(a) .le. 1.0e-6) then
                    ! if the quad coeff. is too low then switch to linear
                    b = den(x+h)-den(x)
                    c = den(x)
                    if (abs(b) .le. 1.0e-3) then
                    	! if linear approx also requires log expansion, upto 3rd order.
                        int=h/(2*c) * (2*d/3 +e + 2*f-(2*b*e+b*f)/(3*c)+2*b**2 *f/(3*c**2))
                        ! print*, 'linearapprox', int, x
                        k= k+1
                    else
                        ! int = h/b *( d/2 +e+ d*c/b + (d*c**2 / b**2 -e*c/b + f) * log(1+b/c) ) ! old one has some typo
                        int = h/(2*b**3) * (b*(-2*c*d+b*(d+2*e)) +2*( c**2 *d - b*c*e+b**2 *f)*log((1+b/c)) )
                        t=t+1
                        ! print*, 'linear', int, x
                    endif
                else
                    tp(1)= ( -b + sqrt(b**2 - 4*a*c) )/(2*a) !x1
                    tp(2)= ( -b - sqrt(b**2 - 4*a*c) )/(2*a)  !x2
                    if ((real(tp(1)) < q .and. real(tp(1)) > p) .or. (real(tp(2)) < q .and. real(tp(2)) > p)) then
                        print*, "Quadratic approx not possible decrease interval to force linear"
                        EXIT
                    endif
                    ! int = h/a * ( d + log(1-1/tp(2)) *(d*(-b)/a +e+(d*tp(1)**2+e*tp(1)+f)/(tp(2)-tp(1))) &
                    !     - log(1-1/tp(1)) * (d*tp(1)**2+e*tp(1)+f)/(tp(2)-tp(1))  )
                    ! print*, 'quad', int, x
                    tp(1) = sqrt(4*a*c - b**2)
                    int = h*((4*a**2*atan((b+2*a)/tp(1))-4*a**2*atan(b/tp(1)))*f+(-b*tp(1)* &
                    log(abs(c+b+a))+b*tp(1)*log(abs(c))+(2*b**2-4*a*c)*atan((b+2*a)/tp(1))+(4*a*c-2*b**2)*&
                    atan(b/tp(1))+2*a*tp(1))*d+e*a*tp(1)*log(abs(c+b+a))-e*a*tp(1)*&
                    log(abs(c))-2*e*a*b*atan((b+2*a)/tp(1))+2*e*a*b*atan(b/tp(1)))/(2*a**2*tp(1))
                endif

                integral = int+integral
                x=x+h
        end do
        write(*,'(A,I5,"/",I5)') "# log expansions in linear:", k,N
        write(*,'(A,I5,"/",I5)') "# total linear approximations:", k+t,N
        write(*,'(A,I5,"/",I5,/,A)') "# total quadratic approximations:", N-k-t,N, repeat("-",10)

    end subroutine

function nume(x)
    ! this is A(y) in integral(A(y)/B(y))
    implicit none
    COMPLEX(dp) :: nume,y
    real(dp) :: x
    y= complex(x,0.)
    ! nume= complex(abs(x),0.)
    ! nume = y-sqrt(y**2-1)
    !nume = sqrt(1-y**2)
    nume=1

    return
end function nume

function deno(x)
    ! this is B(y) in integral(A(y)/B(y))
    implicit none
    COMPLEX(dp) :: deno,y
    real(dp) :: x
    y = complex(x,0)
    ! deno = complex(0.9-x,1.0e-10)
    ! deno = sqrt(y**2 -1) + complex(0.,1.0e-10)
    !deno = 1-y + complex(0.,1.0e-10)
    deno = y**2

end function deno


end program singintquad

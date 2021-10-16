program xxblack_scholes
! 10/16/2021 04:13 PM calculates option price many times for timing purposes
! 10/16/2021 04:08 PM branched from xblack_scholes.f90 to xxblack_scholes.f90
! 10/16/2021 04:07 PM driver for call_price
use kind_mod         , only: dp
use black_scholes_mod, only: call_price
implicit none
integer, parameter :: niter = 10**8
integer            :: iter
real(kind=dp)      :: c,t1,t2
call cpu_time(t1)
! Example 15.6 p360 of Options, Futures, and other Derivatives (2015), 9th edition,
! by John C. Hull
do iter=1,niter
   c = call_price(s=42.0_dp,k=40.0_dp,r=0.1_dp,t=0.5_dp,vol=0.2_dp) ! Hull gets 4.76
end do
print*,"c =",c
call cpu_time(t2)
print*,"time elapsed = ",t2-t1
end program xxblack_scholes

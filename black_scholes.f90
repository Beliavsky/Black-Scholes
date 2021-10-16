module black_scholes_mod
use kind_mod, only: dp
implicit none
private
public :: call_price,alnorm,cumnorm
contains
pure elemental function call_price(s,k,r,t,vol) result(price)
! Black-Scholes price of a European call option
real(kind=dp), intent(in) :: s     ! stock price
real(kind=dp), intent(in) :: k     ! strike price
real(kind=dp), intent(in) :: r     ! annual interest rate -- 0.02 means 2%
real(kind=dp), intent(in) :: t     ! time to expiration in years
real(kind=dp), intent(in) :: vol   ! annualized volatility -- 0.30 means 30%
real(kind=dp)             :: price ! call price
real(kind=dp)             :: d1,d2,vol_sqrt_t
vol_sqrt_t = vol*sqrt(t)
d1 = (log(s/k) + (r + 0.5_dp*vol**2)*t)/vol_sqrt_t
d2 = d1 - vol_sqrt_t
price = s*cumnorm(d1) - k*exp(-r*t)*cumnorm(d2) 
end function call_price
!
!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 2001
pure elemental function alnorm(x, upper) result(fn_val)
real(dp), intent(in)  ::  x
logical , intent(in)  ::  upper
real(dp)              ::  fn_val
!  local variables
real(dp), parameter   ::  zero=0.0_dp, one=1.0_dp, half=0.5_dp, con=1.28_dp
real(dp)              ::  z, y
logical               ::  up
!  machine dependent constants
real(dp), parameter  ::  ltone = 7.0_dp, utzero = 18.66_dp
real(dp), parameter  ::  p = 0.398942280444_dp, q = 0.39990348504_dp,   &
                         r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
                         a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
                         b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
                         c1 = -3.8052e-8_dp, c2 = 3.98064794e-4_dp,     &
                         c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
                         c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
                         d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
                         d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
                         d5 = 30.789933034_dp
up = upper
z = x
if (z < zero) then
   up = .not. up
   z = -z
end if
if (z <= ltone  .or.  (up  .and.  z <= utzero)) then
   y = half*z*z
   if (z > con) then
      fn_val = r*exp( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
   else
      fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
   end if
else
   fn_val = zero
end if
if (.not. up) fn_val = one - fn_val
end function alnorm
!
pure elemental function cumnorm(x) result(fn_val)
! adapted from alnorm with upper = .false.
real(dp), intent(in) ::  x
real(dp)             ::  fn_val
!  local variables
real(dp), parameter  ::  zero=0.0_dp, one=1.0_dp, half=0.5_dp, con=1.28_dp
real(dp)             ::  z, y
logical              ::  up
!  machine dependent constants
real(dp), parameter  ::  ltone = 7.0_dp, utzero = 18.66_dp
real(dp), parameter  ::  p = 0.398942280444_dp, q = 0.39990348504_dp,   &
                         r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
                         a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
                         b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
                         c1 = -3.8052e-8_dp, c2 = 3.98064794e-4_dp,     &
                         c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
                         c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
                         d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
                         d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
                         d5 = 30.789933034_dp
z = x
if (z < zero) then
   up = .true.
   z = -z
else
   up = .false.
end if
if (z <= ltone  .or.  (up  .and.  z <= utzero)) then
   y = half*z*z
   if (z > con) then
      fn_val = r*exp( -y )/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
   else
      fn_val = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
   end if
else
   fn_val = zero
end if
if (.not. up) fn_val = one - fn_val
end function cumnorm
!
end module black_scholes_mod

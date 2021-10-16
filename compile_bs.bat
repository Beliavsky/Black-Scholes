@echo off
:: compile and run with gfortran and Intel Fortran
setlocal
if exist a.exe del a.exe
gfortran -O3 -std=f2018 kind.f90 black_scholes.f90 xxblack_scholes.f90
if exist a.exe a.exe
if exist kind.exe del kind.exe
ifort -O3 -stand:f18 kind.f90 black_scholes.f90 xxblack_scholes.f90
if exist kind.exe kind.exe

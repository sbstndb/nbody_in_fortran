ifort -c module.f90
ifort main.f90 -g -O module.o -o main

# polar_lid_driven_cavity
This repository contain the codes (written in C++) to compute the fluid flow inside a polar lid driven cavity (Ref. Fuchs and Tillmark, 1985, IJNMF). Finite Volume Method is used to discretize the governing equations.

# Reference
Fuchs, L. and Tillmark, N., 1985. Numerical and experimental study of driven flow in a polar cavity. International journal for numerical methods in fluids, 5(4), pp.311-329.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ driver.cpp -o output
```
To run this code

```bash
./output
```

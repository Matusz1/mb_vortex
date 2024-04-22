#ifndef _FOCK_INTERNAL_H_
#define _FOCK_INTERNAL_H_

#include <complex.h>
#include <stdio.h>
#include "fock.h"

typedef unsigned int uint;

double __bessel_energy(uint state_nr);
double complex __bessel_value(uint state_nr, double r, double phi);
double complex __bessel_value_xy(uint state_nr, double x, double y);
//double complex __bessel_value_xy_dx(uint state_nr, double x, double y);
//double complex __bessel_value_xy_dy(uint state_nr, double x, double y);
Vec2z __fock_bessel_grad(uint state_nr, double x, double y);
void __fock_bessel_fprint(FILE* f, uint state_nr);
double __bessel_matrix_element_delta_potential(const FockState* left, const FockState* right);

double __planewave_energy(uint state_nr);
double complex __planewave_value(uint state_nr, double x, double y);
//double complex __planewave_value_dx(uint state_nr, double x, double y);
//double complex __planewave_value_dy(uint state_nr, double x, double y);
Vec2z __fock_planewave_grad(uint state_nr, double x, double y);
void __fock_planewave_fprint(FILE* f, uint state_nr);
double __fock_planewave_matrix_element_delta_potential(const FockState* left, const FockState* right);

#endif /* _FOCK_INTERNAL_H_ */

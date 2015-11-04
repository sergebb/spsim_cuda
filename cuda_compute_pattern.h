/* 
 * File:   cuda_compute_pattern.h
 * Author: serge
 *
 * Created on October 28, 2014, 12:58 PM
 */

#ifndef CUDA_COMPUTE_PATTERN_H
#define	CUDA_COMPUTE_PATTERN_H

#include <stdio.h>

#include "config.h"

#include "molecule.h"
#ifdef	__cplusplus
extern "C" {
#endif

float *cuda_compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size,float atomsf[][9]);
float *cuda_calculate_pixel_solid_angle( int det_nx, int det_ny, float det_width, float det_height, float det_pixel_width, float det_pixel_height, float det_distance);

#ifdef	__cplusplus
}
#endif
#endif	/* CUDA_COMPUTE_PATTERN_H */


/*
 * Copyright (c) 2006 Filipe Maia
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


#ifndef _DIFFRACTION_H_
#define _DIFFRACTION_H_ 1

#include "config.h"

/*
   We're always gonna use 3D diffraction patterns

   the diffraction pattern can be sampled on a regualar grid
   or unstructured grid. In the first case d_s is used to calculate
   the scattering vector coordinates, and the grid order or row major
   (C style).

   In the second case HKL list is used to get the scattering vectors.
   HKL list should always be sorted in row major style (0,0,0; 0,0,1; 0,1,0; 0,1,1...)
*/

typedef struct{
  /* scattering intensity */
  float * I;
  /* scattering pattern */
  double * F_R;
  double * F_I;
  /* maximum scattering vector on each dimension */
  float max_s[3];
  /* scattering vector grid spacing on each dimension */
  float d_s[3];
  /* scattering vector list for unstructured grid */
  float * HKL_list;
  int HKL_list_size;
}Diffraction_Pattern;



/* Calculates the HKL indexes for all the pixels of a given detector */
float * get_HKL_list_for_detector(Options * opts, Experiment * exp, int * HKL_list_size);

/* Computes the diffraction pattern of a given structure on a given set of HKL points */
Diffraction_Pattern * compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size);

void calculate_pixel_solid_angle(Detector * det);

void calculate_thomson_correction(Detector * det);

void calculate_photons_per_pixel(Diffraction_Pattern * pattern, Detector * det, Experiment * experiment);

Diffraction_Pattern * load_pattern_from_file(Detector * det,char * filename,
					     float * HKL_list, int HKL_list_size);
#endif

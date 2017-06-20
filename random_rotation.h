/*
 * File:   random_rotation.h
 * Author: Sergey Bobkov
 *
 * Created on October 27, 2014, 2:16 PM
 */

#ifndef RANDOM_ROTATION_H
#define	RANDOM_ROTATION_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "config.h"
#include "molecule.h"

void init_random_orientation();
Molecule *rotate_molecule(Molecule *mol, double alpha, double beta, double gamma );
Molecule *uniform_rotate_molecule(Molecule *mol, int iteration, int total);
Molecule *uniform_rotate_molecule_2_angles(Molecule *mol, int iteration, int total);
Molecule *random_rotate_molecule(Molecule *mol);
Molecule *random_rotate_molecule_2_angles(Molecule *mol);
int recomended_total_number(int total_in_first_angle);

Molecule *add_molecule_with_shift(Molecule *orig_mol, Molecule *added_mol, double shift_x,double shift_y,double shift_z );

#ifdef	__cplusplus
}
#endif

#endif	/* RANDOM_ROTATION_H */


/*
 * File:   random_rotation.h
 * Author: Sergey Bobkov
 *
 * Created on April 8, 2015, 2:16 PM
 */

#ifndef RANDOM_POSITION_H
#define	RANDOM_POSITION_H
#include "config.h"
#include "molecule.h"
#include <time.h>

    double get_random_number(double standard_deviation);

#ifdef __cplusplus
extern "C" {
#endif
    void get_random_shift_molecule_in_px(float *shift_x,float *shift_y, Options *opts);
#ifdef __cplusplus
}
#endif

#endif	/* RANDOM_POSITION_H */

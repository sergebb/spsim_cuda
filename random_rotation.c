/*
 * Copyright (c) 2014 Bobkov Sergey
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

#include "random_rotation.h"

/**
 * R = ABC
 *
 *     |  cos G sin G   0 |
 * C = | -sin G cos G   0 |
 *     |    0     0     1 |
 *
 *     |  1     0     0   |
 * B = |  0   cos B sin B |
 *     |  0  -sin B cos B |
 *
 *     |  cos A sin A   0 |
 * A = | -sin A cos A   0 |
 *     |    0     0     1 |
 *
 *     |  cos A cos G - sin G sin A cos B      cos A sin G + sin A cos B cos G     sin A sin B |
 * R = | -sin A cos G - cos A cos B sin G     -sin A sin G + cos A cos B cos G     cos A sin B |
 *     |            sin G sin B                         -sin B cos G                cos B      |
 *
 * @param mol Original molecule to rotate
 * @param alpha first rotation around Z
 * @param beta Second rotation aroung Y
 * @param gamma Last rotation around Z
 * @return
 */
Molecule *rotate_molecule(Molecule *mol, double alpha, double beta, double gamma ){
    double sinA,sinB,sinG,cosA,cosB,cosG;
    double R[3][3];
    int i,d;

    Molecule *res = malloc(sizeof(Molecule));

    sinA = sin(alpha);
    cosA = cos(alpha);

    sinB = sin(beta);
    cosB = cos(beta);

    sinG = sin(gamma);
    cosG = cos(gamma);

    R[0][0] = cosA*cosG - sinG*sinA*cosB;
    R[0][1] = cosA*sinG + sinA*cosB*cosG;
    R[0][2] = sinA*sinB;

    R[1][0] = -1*sinA*cosG - cosA*cosB*sinG;
    R[1][1] = -1*sinA*sinG + cosA*cosB*cosG;
    R[1][2] = cosA*sinB;

    R[2][0] = sinG*sinB;
    R[2][1] = -1*sinB*cosG;
    R[2][2] = cosB;

    res->atomic_number = NULL;
    res->pos = NULL;
    res->natoms = mol->natoms;
    res->atomic_number = malloc(sizeof (int)*res->natoms);
    res->pos = malloc(sizeof (float)*res->natoms * 3);
    for (d = 0; d < 3; d++) {
            res->box[d] = mol->box[d];
    }
    for (i = 0; i < res->natoms; i++) {
        res->atomic_number[i] = mol->atomic_number[i];
        for (d = 0; d < 3; d++) {
            res->pos[3 * i + d] = R[d][0]*mol->pos[3 * i] +  R[d][1]*mol->pos[3 * i + 1] +  R[d][2]*mol->pos[3 * i + 2];
        }
    }

    return res;
}

Molecule *random_rotate_molecule(Molecule *mol){
    double alpha, beta, gamma;

    alpha = ((double)rand()/(double)RAND_MAX)*(2*M_PI) - M_PI; // (-Pi,Pi)
    beta  = ((double)rand()/(double)RAND_MAX)*(M_PI) - M_PI/2; // (-Pi/2,Pi/2)
    gamma = ((double)rand()/(double)RAND_MAX)*(2*M_PI) - M_PI; // (-Pi,Pi)

    //printf("%f\t%f\t%f\n",alpha,beta,gamma);

    return rotate_molecule( mol, alpha, beta, gamma );

}

void init_random_orientation(){
    srand(time(NULL));
    return;
}

Molecule *uniform_rotate_molecule(Molecule *mol, int iteration, int total){
    double alpha, beta, gamma;
    int N_alpha, N_beta, N_gamma;
    int i_alpha, i_beta, i_gamma;

    //N*(N/2)*N = Total
    //N = cubic_root(2*Total)

    N_alpha = (int) cbrt((double)2*total);
    N_beta = total/(N_alpha*N_alpha);
    N_gamma = (total)/(N_alpha*N_beta);

    i_alpha = iteration % N_alpha;
    i_beta = (iteration % (N_alpha*N_beta))/N_alpha;
    i_gamma = iteration/(N_alpha*N_beta);

//    printf("%d,%d,%d - (%d:%d N={%d,%d,%d})\n",i_alpha,i_beta,i_gamma,iteration,total,N_alpha,N_beta,N_gamma);

    alpha = (i_alpha/(double)N_alpha - 0.5)*M_PI*2; // (-Pi,Pi)
    beta  = (i_beta/(double)N_beta - 0.5)*M_PI; // (-Pi/2,Pi/2)
    gamma = (i_gamma/(double)N_gamma - 0.5)*M_PI*2; // (-Pi,Pi)

//    printf("%5.2f,%5.2f,%5.2f\n", alpha,beta,gamma);

    return rotate_molecule( mol, alpha, beta, gamma );
}

Molecule *uniform_rotate_molecule_2_angles(Molecule *mol, int iteration, int total){
    double alpha, beta, gamma;
    int N_alpha, N_beta, N_gamma;
    int i_beta, i_gamma;

    //N*(N/2)*N = Total
    //N = cubic_root(2*Total)

    N_alpha = (int) sqrt((double)2*total);
    N_beta = total/(N_alpha);
    N_gamma = (total)/(N_beta);

    i_beta =  iteration / N_gamma;
    i_gamma = iteration%(N_gamma);

//    printf("%d,%d,%d - (%d:%d N={%d,%d,%d})\n",i_alpha,i_beta,i_gamma,iteration,total,N_alpha,N_beta,N_gamma);

    alpha = 0; // Fixed
    beta  = (i_beta/(double)N_beta - 0.5)*M_PI; // (-Pi/2,Pi/2)
    gamma = (i_gamma/(double)N_gamma - 0.5)*M_PI*2; // (-Pi,Pi)

//    printf("%5.2f,%5.2f,%5.2f\n", alpha,beta,gamma);

    return rotate_molecule( mol, alpha, beta, gamma );
}

int recomended_total_number(int total_in_first_angle){
    int N = total_in_first_angle;
    return N*N*(N/2);
}




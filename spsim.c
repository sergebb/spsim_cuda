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

#include <stdlib.h>
#include "config.h"
#include "diffraction.h"
#include "molecule.h"
#include "io.h"
#include "mpi.h"
#include "noise.h"
#include "amplification.h"
#include "random_rotation.h"
#include "random_position.h"

int main(int argc, char ** argv) {
    Options * opts = set_defaults();
    Diffraction_Pattern * pattern;
    /* 20x100x100 nm box */
    float box[3] = {20e-9, 100e-9, 100e-9};
    float * HKL_list;
    int HKL_list_size = 0;
    Molecule * mol = NULL;
    Molecule * original_mol = NULL;
#ifdef MPI
    MPI_Init(&argc, &argv);
#endif
    read_options_file("spsim.conf", opts);
    write_options_file("spsim.confout", opts);
    if (opts->input_type == CHEM_FORMULA) {
        original_mol = get_Molecule_from_formula(opts->chem_formula, box);
    } else if (opts->input_type == PDB) {
        original_mol = get_Molecule_from_pdb(opts->pdb_filename);
    } else {
        fprintf(stderr, "input_type not specified!\n");
        exit(0);
    }
    char file_name[50];
    int i;
    if(!opts->random_position){
        HKL_list = get_HKL_list_for_detector(opts, opts->experiment, &HKL_list_size);
    }
    if(opts->random_orientation){
        init_random_orientation();
    }
    if (opts->gpu_compute) {
        gpu_calculate_pixel_solid_angle(opts->detector);
    } else {
        calculate_pixel_solid_angle(opts->detector);
    }
    if (opts->intermediate_files) {
        write_2D_array_to_vtk(opts->detector->solid_angle, opts->detector->nx, opts->detector->ny, "solid_angle.vtk");
    }

    for (i = 1; i <= opts->number_of_experiments; i++) {
        if (opts->random_orientation) {
            if (mol != NULL) free(mol);
            mol = random_rotate_molecule(original_mol);
        }else if(opts->uniform_rotation){
            if (mol != NULL) free(mol);
            mol = uniform_rotate_molecule_2_angles(original_mol,(i-1),opts->number_of_experiments);
        }else{
            mol = original_mol;
        }
        if(opts->random_position){
            HKL_list = get_HKL_list_for_detector(opts, opts->experiment, &HKL_list_size);
        }

        if (opts->sf_filename[0]) {
            pattern = load_pattern_from_file(opts->detector, opts->sf_filename, HKL_list, HKL_list_size);
        } else {
            if (opts->gpu_compute) {
                pattern = gpu_compute_pattern_on_list(mol, HKL_list, HKL_list_size);
            } else {
                pattern = compute_pattern_on_list(mol, HKL_list, HKL_list_size);
            }
        }
#ifdef MPI
        if (!is_mpi_master()) {
            return 0;
        }
#endif

        calculate_thomson_correction(opts->detector);
        calculate_photons_per_pixel(pattern, opts->detector, opts->experiment);
        generate_poisson_noise(opts->detector);
        calculate_electrons_per_pixel(opts->detector, opts->experiment);
//
        if (opts->intermediate_files) {
            sprintf(file_name, "%04d_scattering_factor.vtk", i);
            write_2D_array_to_vtk(pattern->I, opts->detector->nx, opts->detector->ny, file_name);
            sprintf(file_name, "%04d_thomson_correction.vtk", i);
            write_2D_array_to_vtk(opts->detector->thomson_correction, opts->detector->nx, opts->detector->ny, file_name);
            sprintf(file_name, "%04d_pattern.vtk", i);
            write_2D_array_to_vtk(opts->detector->photons_per_pixel, opts->detector->nx, opts->detector->ny, file_name);
            sprintf(file_name, "%04d_photon_count.vtk", i);
            write_2D_array_to_vtk(opts->detector->photon_count, opts->detector->nx, opts->detector->ny, file_name);
            sprintf(file_name, "%04d_electrons_per_pixel.vtk", i);
            write_2D_array_to_vtk(opts->detector->electrons_per_pixel, opts->detector->nx, opts->detector->ny, file_name);
        }

        calculate_real_detector_output(opts->detector, opts->experiment, opts->beamstop_radius);
        if(opts->binary_output){
            sprintf(file_name, "%04d_real_output.bin", i);
            write_2D_array_to_bin(opts->detector->real_output, opts->detector->nx / opts->detector->binning,
                    opts->detector->ny / opts->detector->binning, file_name);
        }else{
            sprintf(file_name, "%04d_real_output.vtk", i);
            write_2D_array_to_vtk(opts->detector->real_output, opts->detector->nx / opts->detector->binning,
                    opts->detector->ny / opts->detector->binning, file_name);
        }

        free(pattern->HKL_list);
        free(pattern->I);
        free(pattern);
    }
    free(HKL_list);
    free(opts->detector->real_output);
    free(opts->detector->thomson_correction);
    free(opts->detector->photons_per_pixel);
    free(opts->detector->photon_count);
    free(opts->detector->electrons_per_pixel);
    free(opts->detector->solid_angle);
    if(opts->random_orientation)free(mol);
    free(original_mol);
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

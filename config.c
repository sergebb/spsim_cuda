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


#include <string.h>
#include <libconfig.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include "config.h"
#include "mpi.h"


Options * set_defaults(){
  Options * opt = calloc(1,sizeof(Options));
  opt->chem_formula = calloc(1,sizeof(Chem_Formula));
  opt->experiment = calloc(1,sizeof(Experiment));
  opt->detector = calloc(1,sizeof(Detector));
  opt->gpu_compute = 0; //false
  opt->random_orientation = 0; //false
  opt->uniform_rotation = 0; //false
  opt->number_of_experiments = 1;
  opt->beamstop_radius = 0;
  opt->intermediate_files = 1;//true
  opt->binary_output = 0;//false
  opt->random_position = 0;
  opt->position_std = 1.0;
  opt->particle_num_min = opt->particle_num_max = 1;
  opt->particle_dist_min = opt->particle_dist_max = 0.0;
  return opt;
}

void read_options_file(char * filename, Options * res){
  config_t config;
  const char * tmp;
  config_setting_t * root;
  config_setting_t * s;
  config_setting_t * chem_formula;
  int i;
  char buffer[1024];

  config_init(&config);
  root = config_root_setting(&config);
  if(!config_read_file(&config,filename)){
    fprintf(stderr,"Error parsing %s on line %d:\n,%s\n",
	    filename,config_error_line(&config),config_error_text(&config));
    exit(1);
  }

  if(config_lookup(&config,"number_of_dimensions")){
    res->n_dims = config_lookup_int(&config,"number_of_dimensions");
  }
  if((tmp = config_lookup_string(&config,"input_type"))){
    if(strcmp(tmp,"chemical_formula") == 0){
      res->input_type = CHEM_FORMULA;
      res->chem_formula->n_elements = 0;
      chem_formula = config_setting_get_member(root, "chemical_formula");
      for(i = 0;i<100;i++){
	sprintf(buffer,"chemical_formula/atom%d",i);
	if(config_lookup_int(&config,buffer)){
	  res->chem_formula->n_elements++;
	  res->chem_formula->atomic_number = realloc(res->chem_formula->atomic_number,sizeof(int)*res->chem_formula->n_elements);
	  res->chem_formula->atomic_number[res->chem_formula->n_elements-1] = i;
	  res->chem_formula->quantity = realloc(res->chem_formula->quantity,sizeof(int)*res->chem_formula->n_elements);
	  res->chem_formula->quantity[res->chem_formula->n_elements-1] = config_lookup_int(&config,buffer);
	}
      }
    }else if(strcmp(tmp,"pdb") == 0){
      res->input_type = PDB;
      if(config_lookup(&config,"pdb_filename")){
	strcpy(res->pdb_filename,config_lookup_string(&config,"pdb_filename"));
      }else{
	fprintf(stderr,"Error: You need to specify \"pdb_filename\" when you use input_type = \"pdb\"\n");
	exit(0);
      }
    }
  }
  if(config_lookup(&config,"detector_distance")){
    res->detector->distance = config_lookup_float(&config,"detector_distance");
  }
  if(config_lookup(&config,"detector_width")){
    res->detector->width = config_lookup_float(&config,"detector_width");
  }
  if(config_lookup(&config,"detector_height")){
    res->detector->height = config_lookup_float(&config,"detector_height");
  }
  if(config_lookup(&config,"detector_pixel_width")){
    res->detector->pixel_width = config_lookup_float(&config,"detector_pixel_width");
  }
  if(config_lookup(&config,"detector_pixel_height")){
    res->detector->pixel_height = config_lookup_float(&config,"detector_pixel_height");
  }
  if(config_lookup(&config,"detector_quantum_efficiency")){
    res->detector->quantum_efficiency = config_lookup_float(&config,"detector_quantum_efficiency");
  }
  if(config_lookup(&config,"detector_electron_hole_production_energy")){
    res->detector->electron_hole_production_energy = config_lookup_float(&config,"detector_electron_hole_production_energy");
  }
  if(config_lookup(&config,"detector_readout_noise")){
    res->detector->readout_noise = config_lookup_float(&config,"detector_readout_noise");
  }
  if(config_lookup(&config,"detector_dark_current")){
    res->detector->dark_current = config_lookup_float(&config,"detector_dark_current");
  }
  if(config_lookup(&config,"detector_linear_full_well")){
    res->detector->linear_full_well = config_lookup_float(&config,"detector_linear_full_well");
  }
  if(config_lookup(&config,"detector_binning")){
    res->detector->binning = config_lookup_int(&config,"detector_binning");
  }
  if(config_lookup(&config,"detector_maximum_value")){
    res->detector->maximum_value = config_lookup_float(&config,"detector_maximum_value");
  }

  //Cuda support
  if(config_lookup(&config,"compute_using_gpu")){
    res->gpu_compute = config_lookup_bool(&config,"compute_using_gpu");
  }
  //Multiple experiments
  if(config_lookup(&config,"number_of_experiments")){
    res->number_of_experiments = config_lookup_int(&config,"number_of_experiments");
  }
  //Beamstop radius
  if(config_lookup(&config,"beamstop_radius")){
    res->beamstop_radius = config_lookup_int(&config,"beamstop_radius");
  }
  //Random orientation
  if(config_lookup(&config,"random_orientation")){
    res->random_orientation = config_lookup_bool(&config,"random_orientation");
  }
  //uniform rotation
  if(config_lookup(&config,"uniform_rotation")){
    res->uniform_rotation = config_lookup_bool(&config,"uniform_rotation");
  }
  //output as binary file
  if(config_lookup(&config,"binary_output")){
    res->binary_output = config_lookup_bool(&config,"binary_output");
  }
  //Print intermediate files
  if(config_lookup(&config,"intermediate_files")){
    res->intermediate_files = config_lookup_bool(&config,"intermediate_files");
  }
  //Random particle position
  if(config_lookup(&config,"random_position")){
    res->random_position = config_lookup_bool(&config,"random_position");
  }
  //Standard deviation of particle position (in pixels)
  if(config_lookup(&config,"position_std")){
    res->position_std = config_lookup_float(&config,"position_std");
  }

  if(config_lookup(&config,"experiment_wavelength")){
    res->experiment->wavelength = config_lookup_float(&config,"experiment_wavelength");
  }
  if(config_lookup(&config,"experiment_exposure_time")){
    res->experiment->exposure_time = config_lookup_float(&config,"experiment_exposure_time");
  }
  if(config_lookup(&config,"experiment_beam_intensity")){
    res->experiment->beam_intensity = config_lookup_float(&config,"experiment_beam_intensity");
  }
  if(config_lookup_string(&config,"precalculated_sf")){
    strcpy(res->sf_filename,config_lookup_string(&config,"precalculated_sf"));
  }
  if(config_lookup(&config,"particle_number")){
    s = config_lookup(&config,"particle_number");
    res->particle_num_min = config_setting_get_int_elem (s, 0);
    res->particle_num_max = config_setting_get_int_elem (s, 1);
    if(res->particle_num_max == 0){
      res->particle_num_max = res->particle_num_min;
    }
  }
  if(config_lookup(&config,"particle_distance")){
    s = config_lookup(&config,"particle_distance");
	res->particle_dist_min = config_setting_get_float_elem (s, 0);
	res->particle_dist_max = config_setting_get_float_elem (s, 1);
	if(res->particle_dist_max == 0){
	  res->particle_dist_max = res->particle_dist_min;
	}
  }
  res->detector->nx = rint(res->detector->width/res->detector->pixel_width);
  res->detector->ny = rint(res->detector->height/res->detector->pixel_height);
}

void write_options_file(char * filename, Options * res){
  config_t config;
  config_setting_t * root;
  config_setting_t * s;
  config_setting_t * e;
  config_setting_t * chem_formula;
  int i;
  char buffer[1024];
#ifdef MPI
  if(!is_mpi_master()){
    return;
  }
#endif
  config_init(&config);
  root = config_root_setting(&config);
  s = config_setting_add(root,"number_of_dimensions",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->n_dims);
  if(res->input_type == CHEM_FORMULA){
s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"chemical_formula");
    chem_formula = config_setting_add(root,"chemical_formula",CONFIG_TYPE_GROUP);
    for(i = 0;i<res->chem_formula->n_elements;i++){
      sprintf(buffer,"atom%d",res->chem_formula->atomic_number[i]);
      s = config_setting_add(chem_formula,buffer,CONFIG_TYPE_INT);
      config_setting_set_int(s,res->chem_formula->quantity[i]);
    }
  }else if(res->input_type == PDB){
    s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"pdb");
    s = config_setting_add(root,"pdb_filename",CONFIG_TYPE_STRING);
    config_setting_set_string(s,res->pdb_filename);
  }else if(!res->input_type){
    s = config_setting_add(root,"input_type",CONFIG_TYPE_STRING);
    config_setting_set_string(s,"none");
  }
  s = config_setting_add(root,"detector_distance",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->distance);
  s = config_setting_add(root,"detector_width",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->width);
  s = config_setting_add(root,"detector_height",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->height);
  s = config_setting_add(root,"detector_pixel_width",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->pixel_width);
  s = config_setting_add(root,"detector_pixel_height",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->pixel_height);
  s = config_setting_add(root,"detector_quantum_efficiency",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->quantum_efficiency);
  s = config_setting_add(root,"detector_electron_hole_production_energy",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->electron_hole_production_energy);
  s = config_setting_add(root,"detector_readout_noise",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->readout_noise);
  s = config_setting_add(root,"detector_dark_current",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->dark_current);
  s = config_setting_add(root,"detector_linear_full_well",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->linear_full_well);
  s = config_setting_add(root,"detector_binning",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->detector->binning);
  s = config_setting_add(root,"detector_maximum_value",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->detector->maximum_value);

  s = config_setting_add(root,"compute_using_gpu",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->gpu_compute);
  s = config_setting_add(root,"number_of_experiments",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->number_of_experiments);
  s = config_setting_add(root,"beamstop_radius",CONFIG_TYPE_INT);
  config_setting_set_int(s,res->beamstop_radius);
  s = config_setting_add(root,"random_orientation",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->random_orientation);
  s = config_setting_add(root,"uniform_rotation",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->uniform_rotation);
  s = config_setting_add(root,"intermediate_files",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->intermediate_files);
  s = config_setting_add(root,"random_position",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->random_position);
  s = config_setting_add(root,"position_std",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->position_std);
  s = config_setting_add(root,"binary_output",CONFIG_TYPE_BOOL);
  config_setting_set_bool(s,res->binary_output);

  s = config_setting_add(root,"experiment_wavelength",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->wavelength);
  s = config_setting_add(root,"experiment_exposure_time",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->exposure_time);
  s = config_setting_add(root,"experiment_beam_intensity",CONFIG_TYPE_FLOAT);
  config_setting_set_float(s,res->experiment->beam_intensity);

  s = config_setting_add(root,"particle_number",CONFIG_TYPE_ARRAY);
  e = config_setting_add(s, NULL, CONFIG_TYPE_INT);
  config_setting_set_int(e, res->particle_num_min);
  e = config_setting_add(s, NULL, CONFIG_TYPE_INT);
  config_setting_set_int(e, res->particle_num_max);

  s = config_setting_add(root,"particle_distance",CONFIG_TYPE_ARRAY);
  e = config_setting_add(s, NULL, CONFIG_TYPE_FLOAT);
  config_setting_set_int(e, res->particle_dist_min);
  e = config_setting_add(s, NULL, CONFIG_TYPE_FLOAT);
  config_setting_set_int(e, res->particle_dist_max);

  if(res->sf_filename[0]){
    s = config_setting_add(root,"precalculated_sf",CONFIG_TYPE_STRING);
    config_setting_set_string(s,res->sf_filename);
  }
  config_write_file(&config,filename);
}

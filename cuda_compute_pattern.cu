/*
 * File:   cuda_compute_pattern.cu
 * Author: Sergey Bobkov
 *
 * Created on October 28, 2014, 12:35 PM
 */

#include <cuda_runtime_api.h>
#include "cuda_compute_pattern.h"
#define BLOCK_DIM 256

#define CUDA_CHECK_ERROR(err)           \
	if (err != cudaSuccess) {          \
	printf("Cuda error: %s\n", cudaGetErrorString(err));    \
	printf("Error in file: %s, line: %i\n", __FILE__, __LINE__);  \
	}


texture<float, 2, cudaReadModeElementType> tex_atoms;

texture<int, cudaTextureType1D, cudaReadModeElementType> tex_atomic_number;

__device__ float scatt_factor(float d,int Z){
  int i;
  /* the 0.25 is there because the formula is (|H|/2)^2 */
  float res = tex2D(tex_atoms,8,Z);
  for(i = 0;i<4;i++){
    res+= tex2D(tex_atoms,i,Z)*exp((-tex2D(tex_atoms,i+4,Z))*d*d*0.25*1e-20);
  }
  return res;
}

__global__ void cuda_kernel_compute_pattern( int HKL_list_size, int natoms, float *dev_HKL_list, float * dev_pos, float *dev_i, int *dev_is_element_in_molecule) {
    int tNum = blockIdx.x * blockDim.x + threadIdx.x;
    double scattering_factor;
    double scattering_vector_length;
    double scattering_factor_cache[ELEMENTS];
    double dev_F_R,dev_F_I;

    if (tNum < HKL_list_size) {
        dev_F_R = 0;
        dev_F_I = 0;
        float* pixel_coord = dev_HKL_list+3*tNum;
        scattering_vector_length = sqrt(pixel_coord[0]*pixel_coord[0] + pixel_coord[1]*pixel_coord[1] + pixel_coord[2]*pixel_coord[2]);
//        scattering_vector_length = sqrt(dev_HLK_list[tNum]*dev_HLK_list[tNum] + dev_HLK_list[tNum + HKL_list_size]*dev_HLK_list[tNum + HKL_list_size] + dev_HLK_list[tNum + 2*HKL_list_size]*dev_HLK_list[tNum + 2*HKL_list_size]);
        for (int j = 0; j < ELEMENTS; j++) {
           if (dev_is_element_in_molecule[j]) {
                scattering_factor_cache[j] = scatt_factor(scattering_vector_length, j);
            }
        }
        for (int j = 0; j < natoms; j++) {
            scattering_factor = scattering_factor_cache[ tex1Dfetch(tex_atomic_number,j)];
            double sin_phase, cos_phase;
            sincos( 2 * M_PI * (pixel_coord[0] * dev_pos[3*j] + pixel_coord[1]* dev_pos[3*j + 1] + pixel_coord[2] * dev_pos[3*j + 2]), &sin_phase, &cos_phase );
            dev_F_R += scattering_factor * cos_phase;
            dev_F_I += scattering_factor * sin_phase;
        }
        dev_i[tNum] = (float)(dev_F_R * dev_F_R + dev_F_I * dev_F_I);
    }
}

//__global__ void cuda_kernel_mark_molecule( int natoms, int * dev_atomic_number, int *dev_is_element_in_molecule) {
//    int threadNum = blockIdx.x * blockDim.x + threadIdx.x;
//
//    if (threadNum < ELEMENTS) {
//        dev_is_element_in_molecule[threadNum] = 0;
//    }
//    __syncthreads();
//    if (threadNum < natoms) {
//        dev_is_element_in_molecule[dev_atomic_number[threadNum]] = 1;
//    }
//}

__global__ void cuda_kernel_init( int HKL_list_size, int natoms, int * dev_atomic_number, int *dev_is_element_in_molecule, double *dev_f, float *dev_i) {
    int threadNum = blockIdx.x * blockDim.x + threadIdx.x;

    if (threadNum < ELEMENTS) {
        dev_is_element_in_molecule[threadNum] = 0;
    }
    __syncthreads();

    if (threadNum < natoms ) {
        dev_is_element_in_molecule[dev_atomic_number[threadNum]] = 1;
    }
}

__global__ void cuda_kernel_solid_angle( float *det_solid_angle, int det_nx, int det_ny, float det_width,
        float det_height, float det_pixel_width, float det_pixel_height, float det_distance){

    int tNum = blockIdx.x * blockDim.x + threadIdx.x;
    double px,py;
    double corners[4][2];
    double corner_distance[4];
    double projected_pixel_sides[4];
    double r;
    int x = tNum / det_ny;
    int y = tNum % det_ny;
    int i;
    if (tNum < det_nx * det_ny) {
        px = ((x - (det_nx - 1.0) / 2.0) / det_nx) * det_width / 2;
        py = (((det_ny - 1.0) / 2.0 - y) / det_ny) * det_height / 2;
        r = sqrt(det_distance * det_distance + px * px + py * py);
        /* top left */
        corners[0][0] = px - det_pixel_width / 2;
        corners[0][1] = py + det_pixel_height / 2;
        corner_distance[0] = sqrt(det_distance * det_distance + corners[0][0] * corners[0][0] + corners[0][1] * corners[0][1]);
        /* top right */
        corners[1][0] = px + det_pixel_width / 2;
        corners[1][1] = py + det_pixel_height / 2;
        corner_distance[1] = sqrt(det_distance * det_distance + corners[1][0] * corners[1][0] + corners[1][1] * corners[1][1]);
        /* bottom right */
        corners[2][0] = px + det_pixel_width / 2;
        corners[2][1] = py - det_pixel_height / 2;
        corner_distance[2] = sqrt(det_distance * det_distance + corners[2][0] * corners[2][0] + corners[2][1] * corners[2][1]);
        /* bottom left */
        corners[3][0] = px - det_pixel_width / 2;
        corners[3][1] = py - det_pixel_height / 2;
        corner_distance[3] = sqrt(det_distance * det_distance + corners[3][0] * corners[3][0] + corners[3][1] * corners[3][1]);
        /* project on plane*/
        for (i = 0; i < 4; i++) {
            corners[i][0] *= r / corner_distance[i];
            corners[i][1] *= r / corner_distance[i];
        }
        /* top */
        projected_pixel_sides[0] = sqrt((corners[0][0] - corners[1][0])*(corners[0][0] - corners[1][0])+(corners[0][1] - corners[1][1])*(corners[0][1] - corners[1][1]));
        /* left */
        projected_pixel_sides[1] = sqrt((corners[0][0] - corners[3][0])*(corners[0][0] - corners[3][0])+(corners[0][1] - corners[3][1])*(corners[0][1] - corners[3][1]));
        det_solid_angle[tNum] = projected_pixel_sides[0] * projected_pixel_sides[1] / (r * r);

    }
}

/*
 *
 */
float *cuda_compute_pattern_on_list(Molecule * mol, float * HKL_list, int HKL_list_size,float atomsf[][9]){
    float *dev_HLK_list;
    int * dev_atomic_number;
    float * dev_pos;
    float *res_intencity;
//    double *res_f;
    float *dev_i;
    double *dev_f;
    cudaArray *dev_atomsf;

    cudaEvent_t start, stop;
    float time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);


//    float *tmp_HKL_list = new float[3*HKL_list_size];
//    for(int i=0;i<HKL_list_size;i++){
//        tmp_HKL_list[i] = HKL_list[3*i];
//        tmp_HKL_list[i+HKL_list_size] = HKL_list[3*i+1];
//        tmp_HKL_list[i+2*HKL_list_size] = HKL_list[3*i+2];
//    }
//    CUDA_CHECK_ERROR(cudaMemcpy(dev_HLK_list, tmp_HKL_list, 3*sizeof (float)*HKL_list_size, cudaMemcpyHostToDevice));
//    delete[] tmp_HKL_list;

//    float *tmp_pos = new float[3*mol->natoms];
//    for(int i=0;i<mol->natoms;i++){
//        tmp_pos[i] = mol->pos[3*i];
//        tmp_pos[i+mol->natoms] = mol->pos[3*i+1];
//        tmp_pos[i+2*mol->natoms] = mol->pos[3*i+2];
//    }
//    CUDA_CHECK_ERROR(cudaMemcpy(dev_pos, tmp_pos, 3*sizeof (float) * mol->natoms, cudaMemcpyHostToDevice));
//    delete[] tmp_pos;

    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_HLK_list), 3*sizeof (float) * HKL_list_size));
    CUDA_CHECK_ERROR(cudaMemcpy(dev_HLK_list, HKL_list, 3*sizeof (float)*HKL_list_size, cudaMemcpyHostToDevice));

    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_atomic_number), sizeof (int) * mol->natoms));
    CUDA_CHECK_ERROR(cudaMemcpy(dev_atomic_number, mol->atomic_number, sizeof (int) * mol->natoms, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaBindTexture((size_t)0,tex_atomic_number, dev_atomic_number, sizeof (int) * mol->natoms));

    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_pos), 3*sizeof (float) * mol->natoms));
    CUDA_CHECK_ERROR(cudaMemcpy(dev_pos, mol->pos, 3*sizeof (float) * mol->natoms, cudaMemcpyHostToDevice));

    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_f), 2*sizeof(double)*HKL_list_size));
    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_i), sizeof(float)*HKL_list_size));
//    CUDA_CHECK_ERROR(cudaMallocHost((void**) &(res_intencity), sizeof(float)*HKL_list_size));
//    CUDA_CHECK_ERROR(cudaMallocHost((void**) &(res_f), 2*sizeof(double)*HKL_list_size));
    res_intencity = (float*)malloc(sizeof(float)*HKL_list_size);

    int *dev_is_element_in_molecule;
    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_is_element_in_molecule), sizeof(int)*ELEMENTS));

    //Create Channel Descriptor. float is just for example. Change it to required data type.
    cudaChannelFormatDesc channel = cudaCreateChannelDesc<float>();

    //Allocate Memory
    CUDA_CHECK_ERROR(cudaMallocArray(&dev_atomsf, &channel, 9, ELEMENTS, cudaArrayDefault));
    CUDA_CHECK_ERROR(cudaMemcpyToArray(dev_atomsf, 0, 0, atomsf, sizeof(float)*9*ELEMENTS, cudaMemcpyHostToDevice));
    CUDA_CHECK_ERROR(cudaBindTextureToArray(tex_atoms, dev_atomsf, channel));

    dim3 gridSize = dim3(HKL_list_size % BLOCK_DIM == 0 ? HKL_list_size / BLOCK_DIM : HKL_list_size / BLOCK_DIM + 1, 1, 1);
    dim3 gridSizeAtoms = dim3(ELEMENTS % BLOCK_DIM == 0 ? ELEMENTS / BLOCK_DIM : ELEMENTS / BLOCK_DIM + 1, 1, 1);
    dim3 blockSize = dim3(BLOCK_DIM, 1, 1);

    cuda_kernel_init<<< gridSize, blockSize >>>( HKL_list_size, mol->natoms, dev_atomic_number, dev_is_element_in_molecule, dev_f, dev_i);
    cuda_kernel_compute_pattern<<< gridSize, blockSize >>>(HKL_list_size, mol->natoms, dev_HLK_list, dev_pos, dev_i, dev_is_element_in_molecule);

    CUDA_CHECK_ERROR(cudaDeviceSynchronize());

    CUDA_CHECK_ERROR(cudaMemcpy(res_intencity, dev_i, sizeof(float)*HKL_list_size, cudaMemcpyDeviceToHost));
//    CUDA_CHECK_ERROR(cudaMemcpy(res_f, dev_f, 2*sizeof(double)*HKL_list_size, cudaMemcpyDeviceToHost));

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
//    printf ("Time for the gpu compute: %f ms\n", time);

    CUDA_CHECK_ERROR(cudaFreeArray(dev_atomsf));

    CUDA_CHECK_ERROR(cudaFree(dev_HLK_list));
    CUDA_CHECK_ERROR(cudaFree(dev_atomic_number));
    CUDA_CHECK_ERROR(cudaFree(dev_pos));
    CUDA_CHECK_ERROR(cudaFree(dev_f));
    CUDA_CHECK_ERROR(cudaFree(dev_i));
    CUDA_CHECK_ERROR(cudaFree(dev_is_element_in_molecule));


    return res_intencity;

}


float *cuda_calculate_pixel_solid_angle( int det_nx, int det_ny, float det_width,
        float det_height, float det_pixel_width, float det_pixel_height, float det_distance){

    float *result;
    float *dev_solid_angle;

    CUDA_CHECK_ERROR(cudaMalloc((void**) &(dev_solid_angle), sizeof(float)*det_nx*det_ny));
    CUDA_CHECK_ERROR(cudaMallocHost((void**) &(result), sizeof(float)*det_nx*det_ny));
    result = (float*)malloc(sizeof(float)*det_nx*det_ny);

    dim3 gridSize = dim3(det_nx*det_ny % BLOCK_DIM == 0 ? det_nx*det_ny / BLOCK_DIM : det_nx*det_ny / BLOCK_DIM + 1, 1, 1);
    dim3 blockSize = dim3(BLOCK_DIM, 1, 1);

    cuda_kernel_solid_angle<<< gridSize, blockSize >>>(dev_solid_angle, det_nx, det_ny, det_width,
        det_height, det_pixel_width, det_pixel_height, det_distance);

    CUDA_CHECK_ERROR(cudaDeviceSynchronize());
    CUDA_CHECK_ERROR(cudaMemcpy(result, dev_solid_angle, sizeof(float)*det_nx*det_ny, cudaMemcpyDeviceToHost));

    CUDA_CHECK_ERROR(cudaFree(dev_solid_angle));

    return result;


}

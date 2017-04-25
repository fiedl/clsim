// These tests try to run the intersection code on the GPU.
//
// The example code is taken from:
// https://www.fixstars.com/en/opencl/book/OpenCLProgrammingBook/first-opencl-program/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define MEM_SIZE (10)
#define MAX_SOURCE_SIZE (0x100000)

int main()
{
  cl_device_id device_id = NULL;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;
  cl_mem memobj = NULL;
  cl_program program = NULL;
  cl_kernel kernel = NULL;
  cl_platform_id platform_id = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;

  float *results;
  results = (float *)malloc(MEM_SIZE * sizeof(float));

  FILE *fp;
  char fileName[] = "./intersection_opencl_tests_kernel.cl";
  char *source_str;
  size_t source_size;

  /* Load the source code containing the kernel*/
  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  /* Get Platform and Device Info */
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);

  /* Create OpenCL context */
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

  /* Create Command Queue */
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);

  /* Create Memory Buffer */
  memobj = clCreateBuffer(context, CL_MEM_READ_WRITE,MEM_SIZE * sizeof(float), NULL, &ret);

  /* Create Kernel Program from the source */
  program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
  (const size_t *)&source_size, &ret);

  /* Build Kernel Program */
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);

  // http://stackoverflow.com/a/9467325/2066546
  if (ret == CL_BUILD_PROGRAM_FAILURE) {
    // Determine the size of the log
    size_t log_size;
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    // Allocate memory for the log
    char *log = (char *) malloc(log_size);

    // Get the log
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

    // Print the log
    printf("%s\n", log);
  }

  /* Create OpenCL Kernel */
  kernel = clCreateKernel(program, "test_intersection_with_gpu", &ret);

  /* Set OpenCL Kernel Parameters */
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);

  /* Execute OpenCL Kernel */
  ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);

  /* Copy results from the memory buffer */
  ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
  MEM_SIZE * sizeof(float),results, 0, NULL, NULL);

  /* Display Result */
  printf("number_of_intersections: %e\n", results[0]);
  printf("intersection_x1: %e\n", results[1]);
  printf("intersection_y1: %e\n", results[2]);
  printf("intersection_x2: %e\n", results[3]);
  printf("intersection_y2: %e\n", results[4]);
  printf("intersection_s1: %e\n", results[5]);
  printf("intersection_s2: %e\n", results[6]);

  /* Finalization */
  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(memobj);
  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  free(source_str);

  return 0;
}

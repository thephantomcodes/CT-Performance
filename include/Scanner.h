#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

enum ProjectionDirection
{
  Forward,
  Backward
};

struct Scanner
{   
    double scanning_radius;
    double detector_length;
    int num_pixels;
    int num_views;
    int num_detectors;
    double phantom_radius;
    double field_of_view;
};

__global__ void project(Scanner scanner, double *img, double *sinogram, int view_begin, int view_end, ProjectionDirection projectionDirection);
// __global__ void rampFilter(int N, double *in);
__device__ void projectPoint(double src[2], double pt[2], double x);
__device__ void projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2]);
__device__ void intervalsIntersect(double interval1[2], double interval2[2]);
__device__ void rotatePoint(double point[2], double theta);

#endif  // PROJ_PARAMS_H

#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

enum class ProjectionDirection
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

void project(Scanner scanner, double *img, double *sinogram, int view_begin, int view_end, ProjectionDirection projectionDirection);
void rampFilter(int N, double *in);
double projectPoint(double src[2], double pt[2], double x);
void projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2]);
bool intervalsIntersect(double interval1[2], double interval2[2]);
void rotatePoint(double point[2], double theta);

#endif  // PROJ_PARAMS_H

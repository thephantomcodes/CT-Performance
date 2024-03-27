#include "Scanner.h"
#include <math.h>
#include <fftw3.h>

//__global__
double projectPoint(double src[2], double pt[2], double x)
{
  if(src[0] == pt[0]) return pt[0];
  
  double m = (src[1] - pt[1]) / (src[0] - pt[0]);
  double b = pt[1] - m*pt[0];
  return m*x + b;
}

//__global__
void projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2])
{
  interval[0] = projectPoint(src, pt1, x);
  interval[1] = projectPoint(src, pt2, x);
  if(interval[0] > interval[1])
    std::swap(interval[0], interval[1]);
}

//__global__
bool intervalsIntersect(double interval1[2], double interval2[2])
{
  return interval1[0] < interval2[1] && interval2[0] < interval1[1];
}

//__global__
void rotatePoint(double point[2], double theta)
{
  theta = 3.14159*theta/180.0;
  double point_copy[2] = {point[0], point[1]};
  double vec[2] = {cos(theta) , -sin(theta)};
  point[0] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
  vec[0] = sin(theta);
  vec[1] = cos(theta);
  point[1] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
}

//__global__
void project(Scanner scanner, double *img, double *sinogram, int view_begin, int view_end, ProjectionDirection projectionDirection)
{
  double det_len = scanner.detector_length/scanner.num_detectors;
  double det_begin = 0.5*scanner.detector_length;
  double col_begin = scanner.phantom_radius;
  double px_width = 2.0*scanner.phantom_radius/scanner.num_pixels;
  double rotation_delta = scanner.field_of_view/scanner.num_views;
  // double *col_sums = (double *)malloc(scanner.num_pixels*scanner.num_pixels * sizeof(double));
  // double *row_sums = (double *)malloc(scanner.num_detectors*scanner.num_views * sizeof(double));

  for(int v=view_begin; v<view_end; v++)
  {
    double theta = std::fmod(v*rotation_delta, 360.0);
    bool swap_indices = (theta > 45.0 && theta <= 135.0) || (theta > 225.0 && theta <= 315.0);
    theta -= ((double)swap_indices)*90.0;

    double src[] = {scanner.scanning_radius, 0.0};
    rotatePoint(src, theta);
    
    for(int d=0; d<scanner.num_detectors; d++)
    {
      double det1[] = {-scanner.scanning_radius, det_begin - d*det_len};
      double det2[] = {-scanner.scanning_radius, det_begin - (d+1)*det_len};
      rotatePoint(det1, theta);
      rotatePoint(det2, theta);
      
      double det_proj_interval[2];
      projectInterval(src, det1, det2, 0, det_proj_interval);
      
      for(int c=0; c<scanner.num_pixels; c++)
      {
        double px_x = (c+0.5)*px_width - col_begin;
        double cos_correction1 = src[0]/std::sqrt((det_proj_interval[0] - src[1])*(det_proj_interval[0] - src[1]) + src[0]*src[0]);
        double cos_correction2 = src[0]/std::sqrt((det_proj_interval[1] - src[1])*(det_proj_interval[1] - src[1]) + src[0]*src[0]);
        double cos_correction = std::abs(0.5*(cos_correction1 + cos_correction2));
        
        double pxb1 = (projectPoint(src, det1, px_x));
        double pxb2 = (projectPoint(src, det2, px_x));
        if(pxb1 > pxb2) std::swap(pxb1, pxb2);
        
        int px_bound1 = std::max((int)std::floor((pxb1 + col_begin)/px_width), 0);
        int px_bound2 = std::min((int)std::ceil((pxb2 + col_begin)/px_width), scanner.num_pixels);
        
        for(int r=px_bound1; r<px_bound2; r++)
        {
          double px1[] = {px_x, (r)*px_width - col_begin};
          double px2[] = {px_x, (r+1)*px_width - col_begin};
          double px_proj_interval[2];
          projectInterval(src, px1, px2, 0, px_proj_interval);
          
          if(intervalsIntersect(px_proj_interval, det_proj_interval))
          {
            double det_width = det_proj_interval[1] - det_proj_interval[0];
            double det_px_overlap = std::min(det_proj_interval[1], px_proj_interval[1]) - std::max(det_proj_interval[0], px_proj_interval[0]); 
            double weight = det_px_overlap / (det_width * cos_correction);
            int sinogram_index = (v+1)*scanner.num_detectors - d-1;
            int img_index = (r+1)*scanner.num_pixels - c-1;
            if(swap_indices)
              img_index = (c)*scanner.num_pixels + r;
            if(projectionDirection == ProjectionDirection::Forward)
              (sinogram)[sinogram_index] += weight * (img)[img_index];
            else
              (img)[img_index] += weight * (sinogram)[sinogram_index];
            
            #ifdef GEN_SART_WEIGHTS
            if(projectionDirection == ProjectionDirection::Forward)
            {
              row_sums[sinogram_index] += weight;
              col_sums[img_index] += weight;
            }
            #endif
          }
        }
      }
    }
  }
}

void rampFilter(int N, double *in)
{
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));
  fftw_plan fwd, inv;
  
  for(int i=0; i<N; i++)
  {
    fwd = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(fwd);
    for(int j=0; j<(N/2 + 1); j++)
    {
      double gain = (double)(j+1)/(N/2 + 1);
      out[j][0] *= gain;
      out[j][1] *= gain;
    }
    inv = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    fftw_execute(inv);
    fftw_destroy_plan(fwd);
    fftw_destroy_plan(inv);
    in += N;
  }
  
  fftw_free(out);
}

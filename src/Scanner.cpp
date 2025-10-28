#include "Scanner.h"
#include <iostream>
#include <cmath>
#include <cstring>
#include <fftw3.h>
#include <stdint.h>
#include <string>

#ifdef INSTR_RDTSC
#include <x86intrin.h>
#endif

#define PI (3.1415926535)
#define TO_RADIANS(X) (X*(PI/180.0))

namespace CT
{
	Scanner::Scanner(double scanning_radius_, double detector_length_, int num_pixels_, int num_views_, int num_detectors_, double phantom_radius_, double field_of_view_, double phase_)
    : scanning_radius(scanning_radius_)
    , detector_length(detector_length_)
    , num_pixels(num_pixels_)
    , num_views(num_views_)
    , num_detectors(num_detectors_)
    , phantom_radius(phantom_radius_)
    , field_of_view(field_of_view_)
    , phase(phase_)
  {
    det_len = detector_length/num_detectors;
    det_begin = 0.5*detector_length;
    col_begin = phantom_radius;
    px_width = 2.0*phantom_radius/num_pixels;
    rotation_delta = field_of_view/num_views;
    col_sums.resize(num_pixels*num_pixels);
    row_sums.resize(num_detectors*num_views);
  }
  
  void Scanner::PrintProjectionParameters()
  {
    std::cout << "scanning_radius: " << scanning_radius << '\n'
      << "detector_length: " << detector_length << '\n'
      << "num_pixels: " << num_pixels << "x" << num_pixels << '\n'
      << "num_views: " << num_views << '\n'
      << "num_detectors: " << num_detectors << '\n'
      << "phantom_radius: " << phantom_radius << '\n'
      << "field_of_view: " << field_of_view << std::endl;
  }

  double Scanner::projectPoint(double src[2], double pt[2], double x)
  {
    if(src[0] == pt[0]) return pt[0];

    #ifdef INSTR_RDTSC
    uint64_t st = __rdtsc();
    #endif

    double m = (src[1] - pt[1]) / (src[0] - pt[0]);
    double b = pt[1] - m*pt[0];

    #ifdef INSTR_RDTSC
    tick_counter[ProjectPoint] += __rdtsc()-st;
    call_counter[ProjectPoint]++;
    #endif

    return m*x + b;
  }

  void Scanner::projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2])
  {
    interval[0] = projectPoint(src, pt1, x);
    interval[1] = projectPoint(src, pt2, x);
    if(interval[0] > interval[1])
      std::swap(interval[0], interval[1]);
  }

  bool Scanner::intervalsIntersect(double interval1[2], double interval2[2])
  {
  #ifdef INSTR_RDTSC
    uint64_t st = __rdtsc();
  #endif
    bool b = interval1[0] < interval2[1] && interval2[0] < interval1[1];

  #ifdef INSTR_RDTSC
    tick_counter[IntervalsIntersect] += __rdtsc()-st;
    call_counter[IntervalsIntersect]++;
  #endif
    return b;
  }

  void Scanner::rotatePoint(double point[2], double theta)
  {
  #ifdef INSTR_RDTSC
    uint64_t st = __rdtsc();
  #endif

    theta = TO_RADIANS(theta);
    double point_copy[2] = {point[0], point[1]};
    double vec[2] = {cos(theta) , -sin(theta)};
    point[0] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
    vec[0] = sin(theta);
    vec[1] = cos(theta);
    point[1] = point_copy[0]*vec[0] + point_copy[1]*vec[1];

  #ifdef INSTR_RDTSC
    tick_counter[RotatePoint] += __rdtsc()-st;
    call_counter[RotatePoint]++;
  #endif
  }

  void Scanner::rotatePoint(double point[2], double cos_theta, double sin_theta)
  {
  #ifdef INSTR_RDTSC
    uint64_t st = __rdtsc();
  #endif

    // theta = TO_RADIANS(theta);
    double point_copy[2] = {point[0], point[1]};
    double vec[2] = {cos_theta , -sin_theta};
    point[0] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
    vec[0] = sin_theta;
    vec[1] = cos_theta;
    point[1] = point_copy[0]*vec[0] + point_copy[1]*vec[1];

  #ifdef INSTR_RDTSC
    tick_counter[RotatePoint] += __rdtsc()-st;
    call_counter[RotatePoint]++;
  #endif
  }

  void Scanner::project(std::vector<double> *img, std::vector<double> *sinogram, int view_begin, int view_end, ProjectionDirection projectionDirection)
  {
  #ifdef INSTR_RDTSC
    uint64_t st;
    uint64_t start = __rdtsc();
  #endif

    //for each view in thread
    for(int v=view_begin; v<view_end; v++)
    {
      //calculate theta and determine if rotation and index swap are needed
      double theta = std::fmod(v*rotation_delta + phase, 360.0);
      bool swap_indices = (theta > 45.0 && theta <= 135.0) || (theta > 225.0 && theta <= 315.0);
      theta -= ((double)swap_indices)*90.0;

      //rotate src location by theta
      double src[] = {scanning_radius, 0.0};
      // rotatePoint(src, theta);
      theta = TO_RADIANS(theta);
      double cos_theta = cos(theta);
      double sin_theta = sin(theta);
      rotatePoint(src, cos_theta, sin_theta);
      
      //for each detector
      for(int d=0; d<num_detectors; d++)
      {
        //calculate detector boundaries and rotate by theta
        double det1[] = {-scanning_radius, det_begin - d*det_len};
        double det2[] = {-scanning_radius, det_begin - (d+1)*det_len};
        // rotatePoint(det1, theta);
        // rotatePoint(det2, theta);
        rotatePoint(det1, cos_theta, sin_theta);
        rotatePoint(det2, cos_theta, sin_theta);
        
        //project detector boundaries to y-axis
        //det_proj_interval contains y-intercepts of each point
        double det_proj_interval[2];
        projectInterval(src, det1, det2, 0, det_proj_interval);

        #ifdef INSTR_RDTSC
          st = __rdtsc();
        #endif

        double cos_correction1 = src[0]/std::sqrt((det_proj_interval[0] - src[1])*(det_proj_interval[0] - src[1]) + src[0]*src[0]);
        double cos_correction2 = src[0]/std::sqrt((det_proj_interval[1] - src[1])*(det_proj_interval[1] - src[1]) + src[0]*src[0]);
        double cos_correction = std::abs(0.5*(cos_correction1 + cos_correction2));
          
        #ifdef INSTR_RDTSC
          tick_counter[CosCorrect] += __rdtsc()-st;
          call_counter[CosCorrect]++;
        #endif
        
        //for each column
        for(int c=0; c<num_pixels; c++)
        {
          double px_x = (c+0.5)*px_width - col_begin;
          double pxb1 = (projectPoint(src, det1, px_x));
          double pxb2 = (projectPoint(src, det2, px_x));
          if(pxb1 > pxb2) std::swap(pxb1, pxb2);
          
          int px_bound1 = std::max((int)std::floor((pxb1 + col_begin)/px_width), 0);
          int px_bound2 = std::min((int)std::ceil((pxb2 + col_begin)/px_width), num_pixels);
          
          for(int r=px_bound1; r<px_bound2; r++)
          {
            double px1[] = {px_x, (r)*px_width - col_begin};
            double px2[] = {px_x, (r+1)*px_width - col_begin};
            double px_proj_interval[2];
            projectInterval(src, px1, px2, 0, px_proj_interval);
            
            if(intervalsIntersect(px_proj_interval, det_proj_interval))
            {
              #ifdef INSTR_RDTSC
                st = __rdtsc();
              #endif

              double det_width = det_proj_interval[1] - det_proj_interval[0];
              double det_px_overlap = std::min(det_proj_interval[1], px_proj_interval[1]) - std::max(det_proj_interval[0], px_proj_interval[0]); 
              double weight = det_px_overlap / (det_width * cos_correction);
              int sinogram_index = (v+1)*num_detectors - d-1;
              int img_index = (r+1)*num_pixels - c-1;
              if(swap_indices)
                img_index = (c)*num_pixels + r;
              if(projectionDirection == ProjectionDirection::Forward)
                (*sinogram)[sinogram_index] += weight * (*img)[img_index];
              else
                (*img)[img_index] += weight * (*sinogram)[sinogram_index];

              #ifdef INSTR_RDTSC
                tick_counter[StoreProjValue] += __rdtsc()-st;
                call_counter[StoreProjValue]++;
              #endif

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

    #ifdef INSTR_RDTSC
      tick_counter[Total] += __rdtsc()-start;
      call_counter[Total]++;
      PrintTickCounts();
    #endif
  }

#ifdef INSTR_RDTSC
  void Scanner::PrintTickCounts()
  {
    std::ios_base::fmtflags f( std::cout.flags() );
    uint64_t sum = 0;
    for(uint8_t i=ProjectPoint; i<NumOps; i++)
    {
      sum += tick_counter[i];
      uint64_t ave_ticks = call_counter[i] ? tick_counter[i]/call_counter[i] : 0;
      std::cout << std::setprecision(1) << std::scientific
        << SubOpStrings[i] 
        << "\t" << (double)tick_counter[i] 
        << "\t" << (double)call_counter[i]
        << "\t" << (double)ave_ticks 
        << "\t" << (double)tick_counter[i]/tick_counter[Total] << "\n";
    }
    std::cout << "Sum " << (double)sum << " " << (double)tick_counter[Total]/sum << std::endl;
    std::cout.flags( f );
  }
#endif

  // Direct specification of ramp filter in frequency domain
  void Scanner::rampFilter(double *in)
  {
    fftw_complex *out;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (num_detectors/2 + 1));
    fftw_plan fwd, inv;
    
    for(int i=0; i<num_views; i++)
    {
      fwd = fftw_plan_dft_r2c_1d(num_detectors, in, out, FFTW_ESTIMATE);
      fftw_execute(fwd);
      for(int j=0; j<(num_detectors/2 + 1); j++)
      {
        double gain = (double)(j+1)/(num_detectors/2 + 1);
        out[j][0] *= gain;
        out[j][1] *= gain;
      }
      inv = fftw_plan_dft_c2r_1d(num_detectors, out, in, FFTW_ESTIMATE);
      fftw_execute(inv);
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      in += num_detectors;
    }
    
    fftw_free(out);
  }

  // Shepp-Logan specification of sinc-windowed ramp filter in time domain
  void Scanner::rampFilterSL(double *in)
  {
    fftw_complex *out, *ramp;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (num_detectors/2 + 1));
    ramp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (num_detectors/2 + 1));
    fftw_plan fwd, inv, ramp_plan;

    // Ramp filter in time domain
    double *filter;
    filter = (double *) malloc(sizeof(double) * num_detectors);
    const double pi_squared = PI * PI;
    double t_squared = (double)num_detectors*num_detectors;
    for(int i=0; i<num_detectors; i++)
    {
      filter[i] = -2 / (pi_squared * t_squared * (4*i*i - 1));
    }

    // Transform ramp to freq domain
    ramp_plan = fftw_plan_dft_r2c_1d(num_detectors, filter, ramp, FFTW_ESTIMATE);
    fftw_execute(ramp_plan);
    fftw_destroy_plan(ramp_plan);

    // Perform filtering
    for(int i=0; i<num_views; i++)
    {
      fwd = fftw_plan_dft_r2c_1d(num_detectors, in, out, FFTW_ESTIMATE);
      fftw_execute(fwd);
      for(int j=0; j<(num_detectors/2 + 1); j++)
      {
        out[j][0] = out[j][0]*ramp[j][0] - out[j][1]*ramp[j][1];
        out[j][1] = out[j][0]*ramp[j][1] + out[j][1]*ramp[j][0];
      }
      inv = fftw_plan_dft_c2r_1d(num_detectors, out, in, FFTW_ESTIMATE);
      fftw_execute(inv);
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      in += num_detectors;
    }

    fftw_free(out);
    fftw_free(ramp);
  }

  // Perform ramp filtering by factoring filter into differentiator and Hilbert operator
  void Scanner::rampFilterHilbert(double *in)
  {
    // Determine signal lengths
    int n = num_detectors;
    // Padded signal length for Tricomi inversion of Hilbert
    int m = 3*n-2;
    // Real FFT uses symmetry to use half the data for output
    int p = m/2 + 1;

    fftw_complex *out, *kernel_out;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * p);
    kernel_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * p);
    fftw_plan fwd, inv, kernel_plan;

    // Define Hilbert kernel
    double *kernel = (double *)calloc(m, sizeof(double));
    for(int i = 1; i<n; i+=2){
      kernel[i] = 2/(i* PI);
    }
    for(int i = m-n-1; i<m; i+=2){
      kernel[i] = 2/((i-m) * PI);
    }

    // Transform kernel to freq domain
    kernel_plan = fftw_plan_dft_r2c_1d(m, kernel, kernel_out, FFTW_ESTIMATE);
    fftw_execute(kernel_plan);

    // Create space for padded input signal, calloc fills with zeros
    double *in_padded = (double *)calloc(m, sizeof(double));

    // Ramp filter(in) = hilbert of d(in)/dn
    for(int i=0; i<num_views; i++)
    {
      // Compute first difference of input signal (one row of sinogram)
      for(int j=1; j<n; j++) 
        in_padded[j] = in[j+1] - in[j];

      // Transform input to freq domain
      fwd = fftw_plan_dft_r2c_1d(m, in_padded, out, FFTW_ESTIMATE);
      fftw_execute(fwd);

      // Perform Hilbert transform (convolution) in freq domain
      for(int j=0; j<p; j++)
      {
        // (a+bi) * (c + di) = (ac - bd) + (ad + cb)i
        out[j][0] = out[j][0]*kernel_out[j][0] - out[j][1]*kernel_out[j][1];
        out[j][1] = out[j][0]*kernel_out[j][1] + out[j][1]*kernel_out[j][0];
      }

      // Invert previous transform to time domain
      inv = fftw_plan_dft_c2r_1d(num_detectors, out, in_padded, FFTW_ESTIMATE);
      fftw_execute(inv);

      // Clean up memory and advance ptr to next sinogram row.
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      std::memcpy(in, in_padded, n);
      for(int i=0; i<n; i++)
        in[i] *= 1.0/(2.0 * PI);
      std::memset(in_padded, 0x00, m*sizeof(double));
      in += num_detectors;
    }

    free(kernel);
    free(in_padded);
    fftw_destroy_plan(kernel_plan);
    fftw_free(out);
    fftw_free(kernel_out);
  }
}
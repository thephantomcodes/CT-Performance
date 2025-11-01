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
    : m_scanning_radius(scanning_radius_)
    , m_detector_length(detector_length_)
    , m_num_pixels(num_pixels_)
    , m_num_views(num_views_)
    , m_num_detectors(num_detectors_)
    , m_phantom_radius(phantom_radius_)
    , m_field_of_view(field_of_view_)
    , m_phase(phase_)
  {
    m_det_len = m_detector_length/m_num_detectors;
    m_det_begin = 0.5*m_detector_length;
    m_col_begin = m_phantom_radius;
    m_px_width = 2.0*m_phantom_radius/m_num_pixels;
    m_rotation_delta = m_field_of_view/m_num_views;
    col_sums.resize(m_num_pixels*m_num_pixels);
    row_sums.resize(m_num_detectors*m_num_views);
  }
  
  void Scanner::PrintProjectionParameters()
  {
    std::cout << "m_scanning_radius: " << m_scanning_radius << '\n'
      << "m_detector_length: " << m_detector_length << '\n'
      << "m_num_pixels: " << m_num_pixels << "x" << m_num_pixels << '\n'
      << "m_num_views: " << m_num_views << '\n'
      << "m_num_detectors: " << m_num_detectors << '\n'
      << "m_phantom_radius: " << m_phantom_radius << '\n'
      << "m_field_of_view: " << m_field_of_view << std::endl;
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

  // void Scanner::rotatePoint(double point[2], double theta)
  // {
  // #ifdef INSTR_RDTSC
  //   uint64_t st = __rdtsc();
  // #endif

  //   theta = TO_RADIANS(theta);
  //   double point_copy[2] = {point[0], point[1]};
  //   double vec[2] = {cos(theta) , -sin(theta)};
  //   point[0] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
  //   vec[0] = sin(theta);
  //   vec[1] = cos(theta);
  //   point[1] = point_copy[0]*vec[0] + point_copy[1]*vec[1];

  // #ifdef INSTR_RDTSC
  //   tick_counter[RotatePoint] += __rdtsc()-st;
  //   call_counter[RotatePoint]++;
  // #endif
  // }

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
      //calculate v_theta and determine if rotation and index swap are needed
      double v_theta = std::fmod(v*m_rotation_delta + m_phase, 360.0);
      bool v_swap_indices = (v_theta > 45.0 && v_theta <= 135.0) || (v_theta > 225.0 && v_theta <= 315.0);
      v_theta -= ((double)v_swap_indices)*90.0;

      //rotate v_src location by v_theta
      double v_src[] = {m_scanning_radius, 0.0};
      // rotatePoint(v_src, v_theta);
      v_theta = TO_RADIANS(v_theta);
      double v_cos_theta = cos(v_theta);
      double v_sin_theta = sin(v_theta);
      rotatePoint(v_src, v_cos_theta, v_sin_theta);
      
      //for each detector
      for(int d=0; d<m_num_detectors; d++)
      {
        //calculate detector boundaries and rotate by v_theta
        double d_det1[] = {-m_scanning_radius, m_det_begin - d*m_det_len};
        double d_det2[] = {-m_scanning_radius, m_det_begin - (d+1)*m_det_len};
        // rotatePoint(d_det1, v_theta);
        // rotatePoint(d_det2, v_theta);
        rotatePoint(d_det1, v_cos_theta, v_sin_theta);
        rotatePoint(d_det2, v_cos_theta, v_sin_theta);
        
        //project detector boundaries to y-axis
        //d_det_proj_interval contains y-intercepts of each point
        double d_det_proj_interval[2];
        projectInterval(v_src, d_det1, d_det2, 0, d_det_proj_interval);

        #ifdef INSTR_RDTSC
          st = __rdtsc();
        #endif

        double d_cos_correction1 = v_src[0]/std::sqrt((d_det_proj_interval[0] - v_src[1])*(d_det_proj_interval[0] - v_src[1]) + v_src[0]*v_src[0]);
        double d_cos_correction2 = v_src[0]/std::sqrt((d_det_proj_interval[1] - v_src[1])*(d_det_proj_interval[1] - v_src[1]) + v_src[0]*v_src[0]);
        double d_cos_correction = std::abs(0.5*(d_cos_correction1 + d_cos_correction2));
          
        #ifdef INSTR_RDTSC
          tick_counter[CosCorrect] += __rdtsc()-st;
          call_counter[CosCorrect]++;
        #endif
        
        //for each column
        for(int c=0; c<m_num_pixels; c++)
        {
          double c_px_x = (c+0.5)*m_px_width - m_col_begin;
          double c_px_bounds[2] = {0};
          projectInterval(v_src, d_det1, d_det2, c_px_x, c_px_bounds);
          
          int c_r_idx_0 = std::max((int)std::floor((c_px_bounds[0] + m_col_begin)/m_px_width), 0);
          int c_r_idx_1 = std::min((int)std::ceil((c_px_bounds[1] + m_col_begin)/m_px_width), m_num_pixels);
          
          for(int r=c_r_idx_0; r<c_r_idx_1; r++)
          {
            double r_px1[] = {c_px_x, (r)*m_px_width - m_col_begin};
            double r_px2[] = {c_px_x, (r+1)*m_px_width - m_col_begin};
            double r_px_proj_interval[2];
            projectInterval(v_src, r_px1, r_px2, 0, r_px_proj_interval);
            
            // if(intervalsIntersect(r_px_proj_interval, d_det_proj_interval))
            //{
            
            #ifdef INSTR_RDTSC
              st = __rdtsc();
            #endif

            double r_det_width = d_det_proj_interval[1] - d_det_proj_interval[0];
            double r_det_px_overlap = std::min(d_det_proj_interval[1], r_px_proj_interval[1]) - std::max(d_det_proj_interval[0], r_px_proj_interval[0]); 
            double r_weight = r_det_px_overlap / (r_det_width * d_cos_correction);
            int r_sinogram_index = (v+1)*m_num_detectors - d-1;
            int r_img_index = (r+1)*m_num_pixels - c-1;
            if(v_swap_indices)
              r_img_index = (c)*m_num_pixels + r;
            if(projectionDirection == ProjectionDirection::Forward)
              (*sinogram)[r_sinogram_index] += r_weight * (*img)[r_img_index];
            else
              (*img)[r_img_index] += r_weight * (*sinogram)[r_sinogram_index];

            #ifdef INSTR_RDTSC
              tick_counter[StoreProjValue] += __rdtsc()-st;
              call_counter[StoreProjValue]++;
            #endif

            #ifdef GEN_SART_WEIGHTS
              if(projectionDirection == ProjectionDirection::Forward)
              {
                row_sums[r_sinogram_index] += r_weight;
                col_sums[r_img_index] += r_weight;
              }
            #endif

            //}
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
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (m_num_detectors/2 + 1));
    fftw_plan fwd, inv;
    
    for(int i=0; i<m_num_views; i++)
    {
      fwd = fftw_plan_dft_r2c_1d(m_num_detectors, in, out, FFTW_ESTIMATE);
      fftw_execute(fwd);
      for(int j=0; j<(m_num_detectors/2 + 1); j++)
      {
        double gain = (double)(j+1)/(m_num_detectors/2 + 1);
        out[j][0] *= gain;
        out[j][1] *= gain;
      }
      inv = fftw_plan_dft_c2r_1d(m_num_detectors, out, in, FFTW_ESTIMATE);
      fftw_execute(inv);
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      in += m_num_detectors;
    }
    
    fftw_free(out);
  }

  // Shepp-Logan specification of sinc-windowed ramp filter in time domain
  void Scanner::rampFilterSL(double *in)
  {
    fftw_complex *out, *ramp;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (m_num_detectors/2 + 1));
    ramp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (m_num_detectors/2 + 1));
    fftw_plan fwd, inv, ramp_plan;

    // Ramp filter in time domain
    double *filter;
    filter = (double *) malloc(sizeof(double) * m_num_detectors);
    const double pi_squared = PI * PI;
    double t_squared = (double)m_num_detectors*m_num_detectors;
    for(int i=0; i<m_num_detectors; i++)
    {
      filter[i] = -2 / (pi_squared * t_squared * (4*i*i - 1));
    }

    // Transform ramp to freq domain
    ramp_plan = fftw_plan_dft_r2c_1d(m_num_detectors, filter, ramp, FFTW_ESTIMATE);
    fftw_execute(ramp_plan);
    fftw_destroy_plan(ramp_plan);

    // Perform filtering
    for(int i=0; i<m_num_views; i++)
    {
      fwd = fftw_plan_dft_r2c_1d(m_num_detectors, in, out, FFTW_ESTIMATE);
      fftw_execute(fwd);
      for(int j=0; j<(m_num_detectors/2 + 1); j++)
      {
        out[j][0] = out[j][0]*ramp[j][0] - out[j][1]*ramp[j][1];
        out[j][1] = out[j][0]*ramp[j][1] + out[j][1]*ramp[j][0];
      }
      inv = fftw_plan_dft_c2r_1d(m_num_detectors, out, in, FFTW_ESTIMATE);
      fftw_execute(inv);
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      in += m_num_detectors;
    }

    fftw_free(out);
    fftw_free(ramp);
  }

  // Perform ramp filtering by factoring filter into differentiator and Hilbert operator
  void Scanner::rampFilterHilbert(double *in)
  {
    // Determine signal lengths
    int n = m_num_detectors;
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
    for(int i=0; i<m_num_views; i++)
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
      inv = fftw_plan_dft_c2r_1d(m_num_detectors, out, in_padded, FFTW_ESTIMATE);
      fftw_execute(inv);

      // Clean up memory and advance ptr to next sinogram row.
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      std::memcpy(in, in_padded, n);
      for(int i=0; i<n; i++)
        in[i] *= 1.0/(2.0 * PI);
      std::memset(in_padded, 0x00, m*sizeof(double));
      in += m_num_detectors;
    }

    free(kernel);
    free(in_padded);
    fftw_destroy_plan(kernel_plan);
    fftw_free(out);
    fftw_free(kernel_out);
  }
}
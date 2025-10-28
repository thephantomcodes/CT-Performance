#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>
#include <thread>
#include <stdint.h>
#include <string>

#undef GEN_SART_WEIGHTS

namespace CT
{
  enum class ProjectionDirection
  {
    Forward,
    Backward
  };

#ifdef INSTR_RDTSC
  enum SubOps
  {
    ProjectPoint = 0,
    ProjectInterval,
    IntervalsIntersect,
    RotatePoint,
    CosCorrect,
    StoreProjValue,
    Total,

    NumOps,
  };
#endif

  class Scanner
  {
    public:
			Scanner(double scanning_radius_
        , double detector_length_
        , int num_pixels_
        , int num_views_
        , int num_detectors_
        , double phantom_radius_
        , double field_of_view_
        , double phase_);
      
      const double scanning_radius;
      const double detector_length;
      const int num_pixels;
      const int num_views;
      const int num_detectors;
      const double phantom_radius;
      const double field_of_view;
      const double phase;
      
      double det_len;
      double det_begin;
      double col_begin;
      double px_width;
      double rotation_delta;
      std::vector<double> row_sums;
      std::vector<double> col_sums;

      void project(std::vector<double> *img, std::vector<double> *sinogram, int view_begin, int view_end, ProjectionDirection projectionDirection);
      void rampFilter(double *in);
      void rampFilterSL(double *in);
      void rampFilterHilbert(double *in);
      void PrintProjectionParameters();
    #ifdef INSTR_RDTSC
      void PrintTickCounts();
    #endif

    private:
      double projectPoint(double src[2], double pt[2], double x);
      void projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2]);
      bool intervalsIntersect(double interval1[2], double interval2[2]);
      void rotatePoint(double point[2], double theta);
      void rotatePoint(double point[2], double cos_theta, double sin_theta);

    #ifdef INSTR_RDTSC
      uint64_t tick_counter[NumOps] = {0};
      uint64_t call_counter[NumOps] = {0};
      std::string SubOpStrings[NumOps] = {
        "ProjectPoint",
        "ProjectInterval",
        "IntrvlIntersect",
        "RotatePoint",
        "CosCorrect",
        "StoreProjValue",
        "Total"
      };
    #endif
  };
}

#endif  // PROJ_PARAMS_H

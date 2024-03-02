#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

namespace Projection
{
  class ProjectionParameters
  {
    public:
			ProjectionParameters(double scanning_radius_, double detector_length_, int num_pixels_, int num_views_, int num_detectors_, double phantom_radius_, double field_of_view_, double phase_);
      
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
  };
  
  void PrintProjectionParameters(const ProjectionParameters &params);
}

#endif  // PROJ_PARAMS_H

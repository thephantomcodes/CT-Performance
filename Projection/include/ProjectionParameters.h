#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

namespace Projection
{
  class ProjectionParameters
  {
    public:
			ProjectionParameters(double scanning_radius_, double detector_length_, int num_pixels_, int num_views_, int num_detectors_, double phantom_radius_, double field_of_view_);
      
      const double scanning_radius;
      const double detector_length;
      const int num_pixels;
      const int num_views;
      const int num_detectors;
      const double phantom_radius;
      const double field_of_view;
  };
  
  void PrintProjectionParameters(const ProjectionParameters &params);
}

#endif  // PROJ_PARAMS_H

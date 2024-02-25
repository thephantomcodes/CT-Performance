#include "ProjectionParameters.h"
#include <iostream>
#include <cmath>

namespace Projection
{
	ProjectionParameters::ProjectionParameters(double scanning_radius_, double detector_length_, int num_pixels_, int num_views_, int num_detectors_, double phantom_radius_, double field_of_view_)
    : scanning_radius(scanning_radius_)
    , detector_length(detector_length_)
    , num_pixels(num_pixels_)
    , num_views(num_views_)
    , num_detectors(num_detectors_)
    , phantom_radius(phantom_radius_)
    , field_of_view(field_of_view_)
  {
  }
  
  void PrintProjectionParameters(const ProjectionParameters &params)
  {
    std::cout << "scanning_radius: " << params.scanning_radius << '\n'
      << "detector_length: " << params.detector_length << '\n'
      << "num_pixels: " << params.num_pixels << "x" << params.num_pixels << '\n'
      << "num_views: " << params.num_views << '\n'
      << "num_detectors: " << params.num_detectors << '\n'
      << "phantom_radius: " << params.phantom_radius << '\n'
      << "field_of_view: " << params.field_of_view << std::endl;
  }
}

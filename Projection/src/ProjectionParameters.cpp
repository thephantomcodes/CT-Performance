#include "ProjectionParameters.h"
#include <iostream>
#include <cmath>

namespace Projection
{
	const double radian_conversion = std::acos(-1)/180.0;

	PointValue::PointValue(double x_, double y_, double value_): x(x_), y(y_), value(value_) 
	{
	}	
	
	void RotatePointValue(PointValue *pointValue, double theta)
	{
		theta *= radian_conversion;
		double x = pointValue->x;
		double y = pointValue->y;
		pointValue->x = x*std::cos(theta) - y*std::sin(theta);
		pointValue->y = x*std::sin(theta) + y*std::cos(theta);
	}
	
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
  
  double Project_x(PointValue s, PointValue d, double x, double px_width)
  {
  	double m = (s.y - d.y)/(s.x - d.x);
  	return (m*(x - s.x) + s.y); //px_width;
  }
  
  double Project_y(PointValue s, PointValue d, double y, double px_width)
  {
  	double m_inv = (s.x - d.x)/(s.y - d.y);
  	return (m_inv*(y - s.y) + s.x); //px_width;
  }
}

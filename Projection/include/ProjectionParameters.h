#ifndef PROJ_PARAMS_H
#define PROJ_PARAMS_H

namespace Projection
{
	class PointValue
	{
		public:
			PointValue(double x_, double y_, double value_);
			double x, y, value;
	};
	
	void RotatePointValue(PointValue *pointValue, double theta);
	
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
  double Project_x(PointValue s, PointValue d, double x, double px_width);
  double Project_y(PointValue s, PointValue d, double y, double px_width);
}

#endif  // PROJ_PARAMS_H

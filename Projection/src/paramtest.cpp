#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include "ProjectionParameters.h"

void printPoint(double point[], std::string prefix="", std::string suffix="")
{
  std::cout << prefix << "(" << point[0] << "," << point[1] << ")" << suffix;
}

double projectPoint(double src[2], double pt[2], double y)
{
  double m = (src[1] - pt[1]) / (src[0] - pt[0]);
  double b = pt[1] - m*pt[0];
  return m*y + b;
}

void projectInterval(double src[2], double pt1[2], double pt2[2], double y, double interval[2])
{
  interval[0] = projectPoint(src, pt1, y);
  interval[1] = projectPoint(src, pt2, y);
  std::swap(interval[0], interval[1]);
}

bool intervalsIntersect(double interval1[2], double interval2[2])
{
	return interval1[0] < interval2[1] && interval2[0] < interval1[1];
}

int main(int argc, const char* argv[])
{
  bool outputAll = false;
  int sysSize = (argc <= 1) ? 2 : std::atoi(argv[1]);
  
  auto params = Projection::ProjectionParameters(50.0, 55.0, sysSize, sysSize, sysSize, 10.0, 360.0);
  Projection::PrintProjectionParameters(params);
  
  double det_len = params.detector_length/sysSize;
  double det_begin = 0.5*params.detector_length;
  double src[] = {params.scanning_radius, 0.0};
  double col_begin = params.phantom_radius;
  double px_width = 2.0*params.phantom_radius/params.num_pixels;
  printPoint(src, "src ", "\n\n");
  
  for(int d=0; d<params.num_detectors; d++)
  {
    double det1[] = {-params.scanning_radius, det_begin - d*det_len};
    double det2[] = {-params.scanning_radius, det_begin - (d+1)*det_len};
    std::cout << "det: " << d << "\n";
    printPoint(det1, "d1: ", "\n");
    printPoint(det2, "d2: ", "\n");
    
    double det_proj_interval[2];
    projectInterval(src, det1, det2, 0, det_proj_interval);
    printPoint(det_proj_interval, "det proj interval: ", "\n\n");
    
    for(int c=0; c<params.num_pixels; c++)
    {
      double px_x = col_begin - (c+0.5)*px_width; 
      std::cout << "col " << c << ": " << px_x << "\n";
      
      double cos_correction1 = src[0]/std::sqrt((det_proj_interval[0] - src[1])*(det_proj_interval[0] - src[1]) + src[0]*src[0]);
      double cos_correction2 = src[0]/std::sqrt((det_proj_interval[1] - src[1])*(det_proj_interval[1] - src[1]) + src[0]*src[0]);
      double cos_correction = std::abs(0.5*(cos_correction1 + cos_correction2));
      std::cout << "cos " << cos_correction1 << " " << cos_correction2 << " " << cos_correction << "\n";
//      double px_loc = (col_begin - detp1[1])/px_width;
//      std::cout << "px w " << px_width << "\n";

      for(int r=0; r<params.num_pixels; r++)
      {
        double px1[] = {px_x, col_begin - (r)*px_width};
        double px2[] = {px_x, col_begin - (r+1)*px_width};
		    double px_proj_interval[2];
				projectInterval(src, px1, px2, 0, px_proj_interval);
				
				if(intervalsIntersect(px_proj_interval, det_proj_interval))
				{
					double det_width = det_proj_interval[1] - det_proj_interval[0];
					double det_px_overlap = std::min(det_proj_interval[1], px_proj_interval[1]) - std::max(det_proj_interval[0], px_proj_interval[0]); 
    			printPoint(px_proj_interval, "px proj interval: ", " in\n");
    			std::cout << "weight: " << det_px_overlap / (det_width * cos_correction) << "\n";
    		}
    		else
    		{
    			printPoint(px_proj_interval, "px proj interval: ", " out\n");
    		}
    		
      }
      std::cout << "\n";
    }
  }
	
  return 0;
}

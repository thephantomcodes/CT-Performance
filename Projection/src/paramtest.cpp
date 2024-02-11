#include <iostream>
#include <vector>
#include <cmath>
#include "ProjectionParameters.h"
#include "Image.h"
#include "Scanner.h"

void printPoint(double point[], std::string prefix="", std::string suffix="")
{
  std::cout << prefix << "(" << point[0] << "," << point[1] << ")" << suffix;
}

void projectPoint(double src[2], double det[2], double col, double point[2])
{
  double m = (src[1] - det[1]) / (src[0] - det[0]);
  double b = det[1] - m*det[0];
  std::cout << "y = "<< m << "x + " << b << "\n";
  point[0] = col;
  point[1] = m*point[0] + b;
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
  //for(int d=0; d<1; d++)
  {
    double det1[] = {-params.scanning_radius, det_begin - d*det_len};
    double det2[] = {-params.scanning_radius, det_begin - (d+1)*det_len};
    std::cout << d << "\n";
    printPoint(det1, "d1: ", "\n");
    printPoint(det2, "d1: ", "\n");
    
    for(int c=0; c<params.num_pixels; c++)
//    for(int c=0; c<2; c++)
    {
      double col_x = col_begin - (c+0.5)*px_width;
      double point1[2];
      double point2[2];
      projectPoint(src, det1, col_x, point1);
      projectPoint(src, det2, col_x, point2);
      std::cout << "col " << c << ": " << col_x << "\n";
      printPoint(point1, "p1: ", "\n");
      printPoint(point2, "p2: ", "\n");
      double cos_correction1 = src[0]/std::sqrt((point1[1] - src[1])*(point1[1] - src[1]) + src[0]*src[0]);
      double cos_correction2 = src[0]/std::sqrt((point2[1] - src[1])*(point2[1] - src[1]) + src[0]*src[0]);
//      double corrected_y = point1[1]/std::abs(0.5*(cos_correction1 + cos_correction2));
      double px_loc = (col_begin - point1[1])/px_width;
      std::cout << "px " << px_loc << "\n";
//      std::cout << "px w " << px_width << "\n";
      std::cout << "\n";
    }
  }
  
	
  return 0;
}

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>
#include "ProjectionParameters.h"

void printPoint(double point[], std::string prefix="", std::string suffix="")
{
  //std::cout << prefix << "(" << point[0] << "," << point[1] << ")" << suffix;
}

double projectPoint(double src[2], double pt[2], double x)
{
  if(src[0] == pt[0]) return pt[0];
  
  double m = (src[1] - pt[1]) / (src[0] - pt[0]);
  double b = pt[1] - m*pt[0];
  return m*x + b;
}

void projectInterval(double src[2], double pt1[2], double pt2[2], double x, double interval[2])
{
  interval[0] = projectPoint(src, pt1, x);
  interval[1] = projectPoint(src, pt2, x);
  if(interval[0] > interval[1])
    std::swap(interval[0], interval[1]);
}

bool intervalsIntersect(double interval1[2], double interval2[2])
{
  return interval1[0] < interval2[1] && interval2[0] < interval1[1];
}

double dot(std::vector<double> a, std::vector<double> b)
{
  if(a.size() != b.size())
  {
    std::cout << "dot: size mismatch " << a.size() << " != " << b.size() << "\n";
    return 0.0;
  }
  
  double dot_product = 0.0;
  auto itb = b.begin();
  for(auto ita=a.begin(); ita != a.end(); ita++)
  	dot_product += *ita * *itb++;
  return dot_product;
}

void rotatePoint(double point[2], double theta)
{
  theta = 3.14159*theta/180.0;
  double point_copy[2] = {point[0], point[1]};
  double vec[2] = {cos(theta) , -sin(theta)};
  point[0] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
  vec[0] = sin(theta);
  vec[1] = cos(theta);
  point[1] = point_copy[0]*vec[0] + point_copy[1]*vec[1];
}

void writePpmHeader(std::string ofname, int width, int height)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  ofs << "P2\n" << width << ' ' << height << "\n255\n";
  ofs.close();
}

void writePpmData(std::string ofname, std::vector<double> data, int width, double max_val)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary | std::ios_base::app);
  for (int w = 0; w < width; ++w)
  {
    ofs << (int)(255.0*data[w]/max_val) << '\n';
  }
  ofs.close();
}

void readFile(std::string fname, std::vector<double>& vec, int size)
{
  std::fstream fs;
  fs.open(fname, std::fstream::in | std::fstream::binary);
  if (!fs.is_open())
  {
    std::cerr << "Can't find input file " << fname << "\n";
    exit(-1);
  }
  
  double buffer;
  for(int i=0; i<size; i++)
  {
    fs.read(reinterpret_cast<char*>(&buffer), 8);
    vec[i] = buffer;
  }
  fs.close();
}

int main(int argc, const char* argv[])
{
  bool outputAll = false;
  int sysSize = (argc <= 1) ? 2 : std::atoi(argv[1]);
  double fov = (argc <= 2) ? 180.0 : (double)std::atof(argv[2]);
  double phase = (argc <= 3) ? 0.0 : (double)std::atof(argv[3]);
  std::string in_file_prefix = "input/unit_disc_";
  std::string out_file_prefix = "output/sino_unit_disc_";
  char input_img = (argc <= 4) ? 'u' : *argv[4];
  if(input_img == 'p')
  {
    in_file_prefix = "input/phantom_";
    out_file_prefix = "output/sino_phantom_";
  }
  
  
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  auto params = Projection::ProjectionParameters(50.0, 55.0, sysSize, sysSize, sysSize, 10.0, fov);
  
  
  double det_len = params.detector_length/sysSize;
  double det_begin = 0.5*params.detector_length;
  double col_begin = params.phantom_radius;
  double px_width = 2.0*params.phantom_radius/params.num_pixels;
  double rotation_delta = params.field_of_view/params.num_views;
  
  // Dim = number of pixels X detectors * views.
  int P = params.num_pixels*params.num_pixels;
  int D = params.num_detectors*params.num_views;
  
  std::vector<double> img(P);
  readFile(in_file_prefix + std::to_string(params.num_pixels) + ".dat", img, P);
//  for(int q=0; q<params.num_pixels; q++)
//  {
//    for(int r=0; r<params.num_pixels; r++)
//      std::cout << img[q*params.num_pixels + r] << " ";
//    std::cout << '\n';
//  }
  std::vector<double> sinogram(D);
  
//  writePpmHeader("A.ppm", P, D);
  
  
  for(int v=0; v<params.num_views; v++)
  {
  	std::chrono::time_point<std::chrono::system_clock> vstart, vend;
  	vstart = std::chrono::system_clock::now();
    double theta = std::fmod(v*rotation_delta + phase, 360.0);
    bool swap_indices = false;
    bool swap_cols = false;
    if((theta > 45.0 && theta <= 135.0) || (theta > 225.0 && theta <= 315.0))
    {
      theta -= 90.0;
      swap_indices = true;
    }

    double src[] = {params.scanning_radius, 0.0};
    rotatePoint(src, theta);
    
    for(int d=0; d<params.num_detectors; d++)
    {
      double det1[] = {-params.scanning_radius, det_begin - d*det_len};
      double det2[] = {-params.scanning_radius, det_begin - (d+1)*det_len};
      //std::cout << "det: " << d << "\n";
      rotatePoint(det1, theta);
      rotatePoint(det2, theta);
      printPoint(src, "src: ", "\n");
      printPoint(det1, "d1: ", "\n");
      printPoint(det2, "d2: ", "\n");
      
      double det_proj_interval[2];
      projectInterval(src, det1, det2, 0, det_proj_interval);
      printPoint(det_proj_interval, "det proj interval: ", "\n\n");
      
      for(int c=0; c<params.num_pixels; c++)
      {
        double px_x = (c+0.5)*px_width - col_begin;
        if(swap_cols)
        	px_x = col_begin - (c+0.5)*px_width;
        //std::cout << "col " << c << ": " << px_x << "\n";
      
        double cos_correction1 = src[0]/std::sqrt((det_proj_interval[0] - src[1])*(det_proj_interval[0] - src[1]) + src[0]*src[0]);
        double cos_correction2 = src[0]/std::sqrt((det_proj_interval[1] - src[1])*(det_proj_interval[1] - src[1]) + src[0]*src[0]);
        double cos_correction = std::abs(0.5*(cos_correction1 + cos_correction2));
//        std::cout << "cos " << cos_correction1 << " " << cos_correction2 << " " << cos_correction << "\n";
				
				int sinogram_index = 0;
        for(int r=0; r<params.num_pixels; r++)
        {
          double px1[] = {px_x, (r)*px_width - col_begin};
          double px2[] = {px_x, (r+1)*px_width - col_begin};
          double px_proj_interval[2];
          projectInterval(src, px1, px2, 0, px_proj_interval);
          
          if(intervalsIntersect(px_proj_interval, det_proj_interval))
          {
            double det_width = det_proj_interval[1] - det_proj_interval[0];
            double det_px_overlap = std::min(det_proj_interval[1], px_proj_interval[1]) - std::max(det_proj_interval[0], px_proj_interval[0]); 
            printPoint(px_proj_interval, "px proj interval: ", " in\n");
            double weight = det_px_overlap / (det_width * cos_correction);
            sinogram_index = (v+1)*params.num_views - d-1;
            int img_index = (r+1)*params.num_pixels - c-1;
            if(swap_indices)
              img_index = (c)*params.num_pixels + r;
            sinogram[sinogram_index] += weight * img[img_index];
          }
        }
      }
//      writePpmData("A.ppm", A, P);

    }
    vend = std::chrono::system_clock::now(); 
    std::chrono::duration<double> velapsed_seconds = vend - vstart;
    std::cout << "velapsed time: " << v << " " << velapsed_seconds.count() << "s\n";
  }
  end = std::chrono::system_clock::now(); 
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  std::cout << std::setprecision(2) << std::fixed;
  
  writePpmHeader(out_file_prefix + std::to_string(params.num_pixels) + ".ppm", params.num_detectors, params.num_views);
  auto it = std::max_element(sinogram.begin(), sinogram.end());
  std::cout << "max val sino " << *it << '\n';
  writePpmData(out_file_prefix + std::to_string(params.num_pixels) + ".ppm", sinogram, D, *it);
  
//  for(auto row : A)
//  {
//    for(auto entry : row)
//      std::cout << entry << " ";
//    std::cout << '\n';
//  }
//  int i=0;
//  for(auto s : sinogram)
//  {
//    std::cout << s << " ";
//    if(!(++i % 32))
//    	std::cout << '\n';
//  }
//  std::cout << '\n';
  
  return 0;
}

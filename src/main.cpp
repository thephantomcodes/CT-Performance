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
#include "Scanner.h"

void printSums(CT::Scanner& params)
{
  double row_sum_total = 0.0;
  double col_sum_total = 0.0;

  std::cout << "Row Sums Size: " << params.row_sums.size() << "\n";
  for(auto row_sum : params.row_sums)
  {
    // std::cout << row_sum << " ";
    row_sum_total += row_sum;
  }
  std::cout << "Row Sum Total: " << row_sum_total << " " << std::endl;
  
  std::cout << "Col Sums Size: " << params.col_sums.size() << "\n";
  for(auto  col_sum : params.col_sums)
  {
    col_sum_total += col_sum;
  }
  std::cout << "Col Sum Total: " << col_sum_total << " " << std::endl;
  
  std::cout << "Grand Total: " << params.grand_total << std::endl;
  std::cout << "Row Diff: " << row_sum_total - params.grand_total << std::endl;
  std::cout << "Col Diff: " << col_sum_total - params.grand_total << std::endl;
}

void writeWeightData(std::string ofname, CT::Scanner& params)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  ofs << params.grand_total;
  for(auto row_sum : params.row_sums) ofs << row_sum;
  for(auto col_sum : params.col_sums) ofs << col_sum;
  ofs.close();
}

void printPoint(double point[], std::string prefix="", std::string suffix="")
{
  std::cout << prefix << "(" << point[0] << "," << point[1] << ")" << suffix;
}

void writePpmHeader(std::string ofname, int width, int height)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  ofs << "P2\n" << width << ' ' << height << "\n255\n";
  ofs.close();
}

void writePpmData(std::string ofname, std::vector<double> data, int width, double max_val, double min_val)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary | std::ios_base::app);
  for (int w = 0; w < width; ++w)
  {
    ofs << (int)(255.0*(data[w] - min_val)/(max_val - min_val)) << '\n';
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
  int sysSize = (argc <= 1) ? 128 : std::atoi(argv[1]);
  double fov = (argc <= 2) ? 360.0 : (double)std::atof(argv[2]);
  char input_img = (argc <= 3) ? 'u' : *argv[3];
  std::string in_file_prefix = "input/unit_disc_";
  std::string out_file_prefix = "output/sino_unit_disc_";
  std::string img_out_file_prefix = "output/img_unit_disc_";
  std::string sart_weight_prefix = "sart_weights/sart_weight_";
  if(input_img == 'p')
  {
    in_file_prefix = "input/phantom_";
    out_file_prefix = "output/sino_phantom_";
    img_out_file_prefix = "output/img_phantom_";
  }
  
  auto params = CT::Scanner(50.0, 55.0, sysSize, sysSize, sysSize, 10.0, fov, 0.0);
  int total_pixels = params.num_pixels*params.num_pixels;
  int total_detectors = params.num_detectors*params.num_views;
  std::vector<double> img(total_pixels);
  std::vector<double> sinogram(total_detectors);
  readFile(in_file_prefix + std::to_string(params.num_pixels) + ".dat", img, total_pixels);
  sart_weight_prefix
    .append(std::to_string(params.num_pixels))
    .append("_")
    .append(std::to_string((int)params.field_of_view))
    .append(".dat");

#ifdef GEN_SART_WEIGHTS
  params.project(&img, &sinogram, 0, params.num_views, CT::ProjectionDirection::Forward);
  printSums(params);
  writeWeightData(sart_weight_prefix, params);
  return 0;
#endif
  
////////////////////////
// Forward Projection
////////////////////////
  
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  std::thread t1(&CT::Scanner::project, params, &img, &sinogram, 0, params.num_views/4, CT::ProjectionDirection::Forward);
  std::thread t2(&CT::Scanner::project, params, &img, &sinogram, params.num_views/4, params.num_views/2, CT::ProjectionDirection::Forward);
  std::thread t3(&CT::Scanner::project, params, &img, &sinogram, params.num_views/2, 3*params.num_views/4, CT::ProjectionDirection::Forward);
  std::thread t4(&CT::Scanner::project, params, &img, &sinogram, 3*params.num_views/4, params.num_views, CT::ProjectionDirection::Forward);
  
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  
  end = std::chrono::system_clock::now(); 
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
  writePpmHeader(out_file_prefix + std::to_string(params.num_pixels) + ".ppm", params.num_detectors, params.num_views);
  double sino_max = *std::max_element(sinogram.begin(), sinogram.end());
  writePpmData(out_file_prefix + std::to_string(params.num_pixels) + ".ppm", sinogram, total_detectors, sino_max, 0.0);
  
////////////////////////
// Ramp Filtering
////////////////////////

  // start = std::chrono::system_clock::now();
  // params.rampFilter(params.num_detectors, sinogram.data());
  // end = std::chrono::system_clock::now(); 
  // elapsed_seconds = end - start;
  // std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
////////////////////////
// Back Projection
////////////////////////

  start = std::chrono::system_clock::now();
  
  std::thread u1(&CT::Scanner::project, params, &img, &sinogram, 0, params.num_views/4, CT::ProjectionDirection::Backward);
  std::thread u2(&CT::Scanner::project, params, &img, &sinogram, params.num_views/4, params.num_views/2, CT::ProjectionDirection::Backward);
  std::thread u3(&CT::Scanner::project, params, &img, &sinogram, params.num_views/2, 3*params.num_views/4, CT::ProjectionDirection::Backward);
  std::thread u4(&CT::Scanner::project, params, &img, &sinogram, 3*params.num_views/4, params.num_views, CT::ProjectionDirection::Backward);
  
  u1.join();
  u2.join();
  u3.join();
  u4.join();
  
  end = std::chrono::system_clock::now(); 
  elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
  writePpmHeader(img_out_file_prefix + std::to_string(params.num_pixels) + ".ppm", params.num_detectors, params.num_views);
  double img_max = *std::max_element(img.begin(), img.end());
  double img_min = *std::min_element(img.begin(), img.end());
  writePpmData(img_out_file_prefix + std::to_string(params.num_pixels) + ".ppm", img, total_pixels, img_max, img_min);
  
  return 0;
}
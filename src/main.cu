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
#include "../include/Scanner.h"


void writeWeightData(std::string ofname, Scanner& scanner)
{
  // std::fstream ofs;
  // ofs.open(ofname, std::fstream::out | std::fstream::binary);
  // for(auto row_sum : scanner.row_sums) ofs.write(reinterpret_cast<char*>(&row_sum), sizeof(double));
  // for(auto col_sum : scanner.col_sums) ofs.write(reinterpret_cast<char*>(&col_sum), sizeof(double));
  // ofs.close();
}

void writeVector(std::string ofname, std::vector<double>& vec)
{
  // std::fstream ofs;
  // ofs.open(ofname, std::fstream::out | std::fstream::binary);
  // for(auto element : vec) ofs.write(reinterpret_cast<char*>(&element), sizeof(double));
  // ofs.close();
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

void writePpmData(std::string ofname, double *data, int width, double max_val, double min_val)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary | std::ios_base::app);
  for (int w = 0; w < width; ++w)
  {
    ofs << (int)(255.0*(data[w] - min_val)/(max_val - min_val)) << '\n';
  }
  ofs.close();
}

void readFile(std::string fname, double *vec, int size)
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
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    vec[i] = buffer;
  }
  fs.close();
}

void readWeightData(std::string fname, Scanner &scanner, double relax_param)
{
  // const double weight_cap = 1.0;
  // double buffer;
  // std::fstream fs;
  // fs.open(fname, std::fstream::in | std::fstream::binary);

  // if (!fs.is_open())
  // {
  //   std::cerr << "Can't find input file " << fname << "\n";
  //   exit(-1);
  // }

  // for(int i=0; i<scanner.num_pixels*scanner.num_pixels; i++)
  // {
  //   fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
  //   scanner.row_sums[i] = 1.0 / (buffer + 0.000001);
  // }

  // for(int i=0; i<scanner.num_views*scanner.num_detectors; i++)
  // {
  //   fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
  //   scanner.col_sums[i] = relax_param / (buffer + 0.000001);
  // }
  // fs.close();
}

void PrintProjectionParameters(Scanner scanner)
{
  std::cout << "scanning_radius: " << scanner.scanning_radius << '\n'
    << "detector_length: " << scanner.detector_length << '\n'
    << "num_pixels: " << scanner.num_pixels << "x" << scanner.num_pixels << '\n'
    << "num_views: " << scanner.num_views << '\n'
    << "num_detectors: " << scanner.num_detectors << '\n'
    << "phantom_radius: " << scanner.phantom_radius << '\n'
    << "field_of_view: " << scanner.field_of_view << std::endl;
}

int main(int argc, const char* argv[])
{
  int sysSize = (argc <= 1) ? 128 : std::atoi(argv[1]);
  double fov = (argc <= 2) ? 360.0 : (double)std::atof(argv[2]);
  char input_img = (argc <= 3) ? 'u' : *argv[3];
  char operation = (argc <= 4) ? 'p' : *argv[4];
  int sart_iter = (argc <= 5) ? 5 : std::atoi(argv[5]);
  double relax_param = (argc <= 6) ? 1.0 : (double)std::atof(argv[6]);

  std::string in_file_prefix = "input/unit_disc_";
  std::string out_file_prefix = "output/sino_unit_disc_";
  std::string img_out_file_prefix = "output/img_unit_disc_";
  std::string sart_out_file_prefix = "sart_output/sart_unit_disc_";
  std::string sart_weight_prefix = "sart_weights/sart_weight_";

  if(input_img == 'p')
  {
    in_file_prefix = "input/phantom_";
    out_file_prefix = "output/sino_phantom_";
    img_out_file_prefix = "output/img_phantom_";
    sart_out_file_prefix = "sart_output/sart_phantom_";
  }
  
  Scanner scanner = { 50.0, 40.0, sysSize, sysSize, sysSize, 10.0, fov };
  PrintProjectionParameters(scanner);
  int total_pixels = scanner.num_pixels*scanner.num_pixels;
  int total_detectors = scanner.num_detectors*scanner.num_views;
  
  double *img;
  double *sinogram;
  cudaMallocManaged(&img, total_pixels*sizeof(double));
  cudaMallocManaged(&sinogram, total_detectors*sizeof(double));

  readFile(in_file_prefix + std::to_string(scanner.num_pixels) + ".dat", img, total_pixels);
  sart_weight_prefix
    .append(std::to_string(scanner.num_pixels))
    .append("_")
    .append(std::to_string((int)scanner.field_of_view))
    .append("_")
    .append(std::to_string((int)scanner.detector_length))
    .append(".dat");

  // printSums(scanner);
  
////////////////////////
// Forward Projection
////////////////////////

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  start = std::chrono::system_clock::now();
  // project(scanner, img.data(), sinogram.data(), 0, scanner.num_views, ProjectionDirection::Forward);
  project<<<1, 1>>>(scanner, img, sinogram, 0, scanner.num_views, ProjectionDirection::Forward);
  cudaDeviceSynchronize();
  end = std::chrono::system_clock::now(); 
  elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
  writePpmHeader(out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", scanner.num_detectors, scanner.num_views);
  double sino_max = *std::max_element(sinogram, sinogram + total_detectors);
  writePpmData(out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", sinogram, total_detectors, sino_max, 0.0);

////////////////////////
// Ramp Filtering
////////////////////////

//   std::chrono::time_point<std::chrono::system_clock> start, end;
//   std::chrono::duration<double> elapsed_seconds;

//   if(operation == 'f')
//   {
//     start = std::chrono::system_clock::now();
//     scanner.rampFilter(scanner.num_detectors, sinogram.data());
//     end = std::chrono::system_clock::now(); 
//     elapsed_seconds = end - start;
//     std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
//   }
  
// ////////////////////////
// // Back Projection
// ////////////////////////

//   if((operation == 'f') || (operation == 'b'))
//   {
//     project(scanner, &img, &sinogram, ProjectionDirection::Backward);
    
//     writePpmHeader(img_out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", scanner.num_detectors, scanner.num_views);
//     double img_max = *std::max_element(img.begin(), img.end());
//     double img_min = *std::min_element(img.begin(), img.end());
//     writePpmData(img_out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", img, total_pixels, img_max, img_min);
//   }

// ////////////////////////
// // SART
// ////////////////////////
  
//   if(operation == 's')
//   {
//     readWeightData(sart_weight_prefix, scanner, relax_param);

//     std::vector<double> sinogram_error, img_error;
//     sinogram_error.resize(total_detectors);
//     img_error.resize(total_pixels);
//     std::string sart_out_fname;

//     std::fill(img.begin(), img.end(), 0.0); 
//     // project(scanner, &img, &sinogram, ProjectionDirection::Backward);

//     for(int i=0; i<sart_iter; i++)
//     {
//       std::cout << "\nSART Iter: " << i << "\n";
//       //compute error
//       project(scanner, &img, &sinogram_error, ProjectionDirection::Forward);
//       std::transform(sinogram_error.begin(), sinogram_error.end(), sinogram.begin(), sinogram_error.begin(), std::minus<double>());

//       //apply weights and project error
//       std::transform(sinogram_error.begin(), sinogram_error.end(), scanner.row_sums.begin(), sinogram_error.begin(), std::multiplies<double>()); 
//       project(scanner, &img_error, &sinogram_error, ProjectionDirection::Backward);
//       std::transform(img_error.begin(), img_error.end(), scanner.col_sums.begin(), img_error.begin(), std::multiplies<double>());

//       //update img
//       std::transform(img.begin(), img.end(), img_error.begin(), img.begin(), std::minus<double>());

//       sart_out_fname = sart_out_file_prefix;
//       sart_out_fname.append(std::to_string(scanner.num_pixels))
//         .append("_")
//         .append(std::to_string((int)scanner.field_of_view))
//         .append("_")
//         .append(std::to_string(i))
//         .append(".ppm");

//       writePpmHeader(sart_out_fname, scanner.num_pixels, scanner.num_pixels);
//       double _max = *std::max_element(img.begin(), img.end());
//       double _min = *std::min_element(img.begin(), img.end());
//       std::cout<< "img "  << _min << " - " << _max << "\n";
//       writePpmData(sart_out_fname, img, total_pixels, _max, _min);
//     }
  // }

  return 0;
}
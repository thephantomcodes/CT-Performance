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

void printSums(CT::Scanner& scanner)
{
  double row_sum_total = 0.0;
  double col_sum_total = 0.0;

  std::cout << "Row Sums Size: " << scanner.row_sums.size() << "\n";
  for(auto row_sum : scanner.row_sums)
  {
    // std::cout << row_sum << " ";
    row_sum_total += row_sum;
  }
  std::cout << "Row Sum Total: " << row_sum_total << " " << std::endl;
  
  std::cout << "Col Sums Size: " << scanner.col_sums.size() << "\n";
  for(auto  col_sum : scanner.col_sums)
  {
    col_sum_total += col_sum;
  }
  std::cout << "Col Sum Total: " << col_sum_total << " " << std::endl;
  
  std::cout << "Grand Total: " << scanner.grand_total << std::endl;
  std::cout << "Row Diff: " << row_sum_total - scanner.grand_total << std::endl;
  std::cout << "Col Diff: " << col_sum_total - scanner.grand_total << std::endl;
}

void writeWeightData(std::string ofname, CT::Scanner& scanner)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  ofs.write(reinterpret_cast<char*>(&scanner.grand_total), sizeof(double));
  for(auto row_sum : scanner.row_sums) ofs.write(reinterpret_cast<char*>(&row_sum), sizeof(double));;
  for(auto col_sum : scanner.col_sums) ofs.write(reinterpret_cast<char*>(&col_sum), sizeof(double));;
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
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    vec[i] = buffer;
  }
  fs.close();
}

void readSumFile(std::string fname, CT::Scanner &scanner)
{
  std::fstream fs;
  fs.open(fname, std::fstream::in | std::fstream::binary);
  if (!fs.is_open())
  {
    std::cerr << "Can't find input file " << fname << "\n";
    exit(-1);
  }
  
  double buffer;
  fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
  scanner.grand_total = buffer;

  for(int i=0; i<scanner.num_pixels*scanner.num_pixels; i++)
  {
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    scanner.row_sums[i] = buffer;
  }

  for(int i=0; i<scanner.num_views*scanner.num_detectors; i++)
  {
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    scanner.col_sums[i] = buffer;
  }
  fs.close();
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
  std::string sart_weight_prefix = "sart_weights/sart_weight_";
  if(input_img == 'p')
  {
    in_file_prefix = "input/phantom_";
    out_file_prefix = "output/sino_phantom_";
    img_out_file_prefix = "output/img_phantom_";
  }
  
  auto scanner = CT::Scanner(50.0, 55.0, sysSize, sysSize, sysSize, 10.0, fov, 0.0);
  int total_pixels = scanner.num_pixels*scanner.num_pixels;
  int total_detectors = scanner.num_detectors*scanner.num_views;
  std::vector<double> img(total_pixels);
  std::vector<double> sinogram(total_detectors);
  readFile(in_file_prefix + std::to_string(scanner.num_pixels) + ".dat", img, total_pixels);
  sart_weight_prefix
    .append(std::to_string(scanner.num_pixels))
    .append("_")
    .append(std::to_string((int)scanner.field_of_view))
    .append(".dat");

#ifdef GEN_SART_WEIGHTS
  scanner.project(&img, &sinogram, 0, scanner.num_views, CT::ProjectionDirection::Forward);
  printSums(scanner);
  writeWeightData(sart_weight_prefix, scanner);
  std::cout << "GT " << scanner.grand_total << " RS " << scanner.row_sums[0] << " CS " << scanner.col_sums[0] << "\n";
  return 0;
#endif

  readSumFile(sart_weight_prefix, scanner);
  std::cout << "GT " << scanner.grand_total << " RS " << scanner.row_sums[0] << " CS " << scanner.col_sums[0] << "\n";
  
////////////////////////
// Forward Projection
////////////////////////
  
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  std::thread t1(&CT::Scanner::project, scanner, &img, &sinogram, 0, scanner.num_views/4, CT::ProjectionDirection::Forward);
  std::thread t2(&CT::Scanner::project, scanner, &img, &sinogram, scanner.num_views/4, scanner.num_views/2, CT::ProjectionDirection::Forward);
  std::thread t3(&CT::Scanner::project, scanner, &img, &sinogram, scanner.num_views/2, 3*scanner.num_views/4, CT::ProjectionDirection::Forward);
  std::thread t4(&CT::Scanner::project, scanner, &img, &sinogram, 3*scanner.num_views/4, scanner.num_views, CT::ProjectionDirection::Forward);
  
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  
  end = std::chrono::system_clock::now(); 
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  
  writePpmHeader(out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", scanner.num_detectors, scanner.num_views);
  double sino_max = *std::max_element(sinogram.begin(), sinogram.end());
  writePpmData(out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", sinogram, total_detectors, sino_max, 0.0);
  
////////////////////////
// Ramp Filtering
////////////////////////

  if(operation == 'f')
  {
    start = std::chrono::system_clock::now();
    scanner.rampFilter(scanner.num_detectors, sinogram.data());
    end = std::chrono::system_clock::now(); 
    elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  }
  
////////////////////////
// Back Projection
////////////////////////

  if(operation == 'f' || operation == 'b')
  {
    start = std::chrono::system_clock::now();
    
    std::thread u1(&CT::Scanner::project, scanner, &img, &sinogram, 0, scanner.num_views/4, CT::ProjectionDirection::Backward);
    std::thread u2(&CT::Scanner::project, scanner, &img, &sinogram, scanner.num_views/4, scanner.num_views/2, CT::ProjectionDirection::Backward);
    std::thread u3(&CT::Scanner::project, scanner, &img, &sinogram, scanner.num_views/2, 3*scanner.num_views/4, CT::ProjectionDirection::Backward);
    std::thread u4(&CT::Scanner::project, scanner, &img, &sinogram, 3*scanner.num_views/4, scanner.num_views, CT::ProjectionDirection::Backward);
    
    u1.join();
    u2.join();
    u3.join();
    u4.join();
    
    end = std::chrono::system_clock::now(); 
    elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    
    writePpmHeader(img_out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", scanner.num_detectors, scanner.num_views);
    double img_max = *std::max_element(img.begin(), img.end());
    double img_min = *std::min_element(img.begin(), img.end());
    writePpmData(img_out_file_prefix + std::to_string(scanner.num_pixels) + ".ppm", img, total_pixels, img_max, img_min);
  }
  
  return 0;
}
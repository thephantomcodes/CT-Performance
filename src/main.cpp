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
#include <random>
#include "Scanner.h"

void printSums(CT::Scanner& scanner)
{
  double row_sum_total = 0.0;
  double col_sum_total = 0.0;

  std::cout << "Row Sums Size: " << scanner.row_sums.size() << "\n";
  for(auto row_sum : scanner.row_sums) row_sum_total += row_sum;
  std::cout << "Row Sum Total: " << row_sum_total << " " << std::endl;
  
  std::cout << "Col Sums Size: " << scanner.col_sums.size() << "\n";
  for(auto  col_sum : scanner.col_sums) col_sum_total += col_sum;
  std::cout << "Col Sum Total: " << col_sum_total << " " << std::endl;
  
  std::cout << "Diff: " << row_sum_total - col_sum_total << std::endl;
}

void writeWeightData(std::string ofname, CT::Scanner& scanner)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  for(auto row_sum : scanner.row_sums) ofs.write(reinterpret_cast<char*>(&row_sum), sizeof(double));
  for(auto col_sum : scanner.col_sums) ofs.write(reinterpret_cast<char*>(&col_sum), sizeof(double));
  ofs.close();
}

void writeVector(std::string ofname, std::vector<double>& vec)
{
  std::fstream ofs;
  ofs.open(ofname, std::fstream::out | std::fstream::binary);
  for(auto element : vec) ofs.write(reinterpret_cast<char*>(&element), sizeof(double));
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

void writeFile(std::string fname, std::vector<double>& vec, int size)
{
  std::ofstream fs;
  fs.open(fname, std::fstream::out | std::fstream::binary);
  if (!fs.is_open())
  {
    std::cerr << "Can't find output file " << fname << "\n";
    exit(-1);
  }
  
  for(int i=0; i<size; i++)
  {
    // fs << vec[i];
    fs.write (reinterpret_cast<char*>(&vec[i]),sizeof(double));
  }
  fs.close();
}

void readWeightData(std::string fname, CT::Scanner &scanner, double relax_param)
{
  const double weight_cap = 1.0;
  double buffer;
  std::fstream fs;
  fs.open(fname, std::fstream::in | std::fstream::binary);

  if (!fs.is_open())
  {
    std::cerr << "Can't find input file " << fname << "\n";
    exit(-1);
  }

  for(int i=0; i<scanner.num_pixels*scanner.num_pixels; i++)
  {
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    scanner.row_sums[i] = 1.0 / (buffer + 0.000001);
  }

  for(int i=0; i<scanner.num_views*scanner.num_detectors; i++)
  {
    fs.read(reinterpret_cast<char*>(&buffer), sizeof(double));
    scanner.col_sums[i] = relax_param / (buffer + 0.000001);
  }
  fs.close();
}

void project(CT::Scanner& scanner, std::vector<double> *img, std::vector<double> *sinogram, CT::ProjectionDirection projectionDirection, int thread_count)
{
  std::thread *t = new std::thread[thread_count];

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  
  std::cout << "threads: " << thread_count << " views: " << scanner.num_views << '\n';
  for(int i=0; i<thread_count; i++)
  {
    t[i] = std::thread(&CT::Scanner::project, scanner, img, sinogram, i*scanner.num_views/thread_count, (i+1)*scanner.num_views/thread_count, projectionDirection);
    std::cout << "thread: " << i << " start: " << (i)*scanner.num_views/thread_count << " end: " << (i+1)*scanner.num_views/thread_count << '\n';
  }

  for(int i=0; i<thread_count; i++)
  {
    t[i].join();
  }
  
  end = std::chrono::system_clock::now(); 
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

  delete[] t;
}

int main(int argc, const char* argv[])
{
  int num_pixels = (argc <= 1) ? 256 : std::atoi(argv[1]);
  int num_views = (argc <= 2) ? 256 : std::atoi(argv[2]);
  int num_detectors = (argc <= 3) ? 256 : std::atoi(argv[3]);
  double fov = (argc <= 4) ? 360.0 : (double)std::atof(argv[4]);
  std::string input_img = (argc <= 5) ? "phantom_256" : argv[5];
  std::string operation = (argc <= 6) ? "p" : argv[6];
  int thread_count = (argc <= 7) ? 1 : std::atoi(argv[7]);
  int sart_iter = (argc <= 8) ? 5 : std::atoi(argv[8]);
  double relax_param = (argc <= 9) ? 1.0 : (double)std::atof(argv[9]);

  std::string in_file_prefix = "input/" + input_img;
  std::string out_file_prefix = "output/sino_" + input_img;
  std::string img_out_file_prefix = "output/img_" + input_img;
  std::string sart_out_file_prefix = "sart_output/sart_" + input_img;
  std::string sart_weight_prefix = "sart_weights/sart_weight_" + input_img;
  
  auto scanner = CT::Scanner(50.0, 40.0, num_pixels, num_views, num_detectors, 10.0, fov, 0.0);
  scanner.PrintProjectionParameters();
  int total_pixels = scanner.num_pixels*scanner.num_pixels;
  int total_detectors = scanner.num_detectors*scanner.num_views;
  std::vector<double> img(total_pixels);
  std::vector<double> sinogram(total_detectors);
  readFile(in_file_prefix + ".dat", img, total_pixels);
  sart_weight_prefix
    .append(std::to_string(scanner.num_pixels))
    .append("_")
    .append(std::to_string((int)scanner.field_of_view))
    .append("_")
    .append(std::to_string((int)scanner.detector_length))
    .append(".dat");

#ifdef GEN_SART_WEIGHTS
  scanner.project(&img, &sinogram, 0, scanner.num_views, CT::ProjectionDirection::Forward, thread_count);
  printSums(scanner);
  writeWeightData(sart_weight_prefix, scanner);
  return 0;
#endif

  // printSums(scanner);

////////////////////////
// Gaussian Noise
////////////////////////

  if(operation == "g")
  {
    std::cout << "Adding Gaussian noise.\n";
    //use fov as std_dev
    double std_dev = fov;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution d{0.0, std_dev};
    for(int i=0; i<total_pixels; i++)
    {
      img[i] += d(gen);
    }
    std::string guass_out_prefix = in_file_prefix + "_gauss_" + std::to_string(std_dev);
    writeFile(guass_out_prefix + ".dat", img, total_pixels);

    writePpmHeader(guass_out_prefix + ".ppm", scanner.num_pixels, scanner.num_pixels);
    double img_max = *std::max_element(img.begin(), img.end());
    double img_min = *std::min_element(img.begin(), img.end());
    writePpmData(guass_out_prefix + ".ppm", img, total_pixels, img_max, img_min);
    return 0;
  }
  
////////////////////////
// Forward Projection
////////////////////////
  
  std::cout << "Forward projection\n";
  project(scanner, &img, &sinogram, CT::ProjectionDirection::Forward, thread_count);
  
  std::string fname_fp_out = out_file_prefix + "_" + std::to_string(scanner.num_views) + "_" + std::to_string(scanner.num_detectors) + ".ppm";
  writePpmHeader(fname_fp_out, scanner.num_detectors, scanner.num_views);
  double sino_max = *std::max_element(sinogram.begin(), sinogram.end());
  double sino_min = *std::min_element(sinogram.begin(), sinogram.end());
  writePpmData(fname_fp_out, sinogram, total_detectors, sino_max, sino_min);
  
////////////////////////
// Ramp Filtering
////////////////////////

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<double> elapsed_seconds;

  if(operation == "f")
  {
    std::cout << "Ramp Filtering\n";
    start = std::chrono::system_clock::now();
    scanner.rampFilterHilbert(sinogram.data());
    end = std::chrono::system_clock::now(); 
    elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    std::string fname_filt_out = fname_fp_out + "_filt.ppm";
    writePpmHeader(fname_filt_out, scanner.num_detectors, scanner.num_views);
    sino_max = *std::max_element(sinogram.begin(), sinogram.end());
    sino_min = *std::min_element(sinogram.begin(), sinogram.end());
    writePpmData(fname_filt_out, sinogram, total_detectors, sino_max, sino_min);
  }
  
////////////////////////
// Back Projection
////////////////////////

  if((operation == "f") || (operation == "b"))
  {
    std::cout << "Back projection\n";
    project(scanner, &img, &sinogram, CT::ProjectionDirection::Backward, thread_count);
    
    std::string fname_bp_out = img_out_file_prefix + "_" + std::to_string(scanner.num_views) + "_" + std::to_string(scanner.num_detectors) + ".ppm";
    writePpmHeader(fname_bp_out, scanner.num_pixels, scanner.num_pixels);
    // std::cout << "Back projection header done\n";
    double img_max = *std::max_element(img.begin(), img.end());
    double img_min = *std::min_element(img.begin(), img.end());
    writePpmData(fname_bp_out, img, total_pixels, img_max, img_min);
  }

////////////////////////
// SART
////////////////////////
  
  if(operation == "s")
  {
    readWeightData(sart_weight_prefix, scanner, relax_param);

    std::vector<double> sinogram_error, img_error;
    sinogram_error.resize(total_detectors);
    img_error.resize(total_pixels);
    std::string sart_out_fname;

    std::fill(img.begin(), img.end(), 0.0); 
    // project(scanner, &img, &sinogram, CT::ProjectionDirection::Backward);

    for(int i=0; i<sart_iter; i++)
    {
      std::cout << "\nSART Iter: " << i << "\n";
      //compute error
      project(scanner, &img, &sinogram_error, CT::ProjectionDirection::Forward, thread_count);
      std::transform(sinogram_error.begin(), sinogram_error.end(), sinogram.begin(), sinogram_error.begin(), std::minus<double>());

      //apply weights and project error
      std::transform(sinogram_error.begin(), sinogram_error.end(), scanner.row_sums.begin(), sinogram_error.begin(), std::multiplies<double>()); 
      project(scanner, &img_error, &sinogram_error, CT::ProjectionDirection::Backward, thread_count);
      std::transform(img_error.begin(), img_error.end(), scanner.col_sums.begin(), img_error.begin(), std::multiplies<double>());

      //update img
      std::transform(img.begin(), img.end(), img_error.begin(), img.begin(), std::minus<double>());

      sart_out_fname = sart_out_file_prefix;
      sart_out_fname.append(std::to_string(scanner.num_views))
        .append("_")
        .append(std::to_string(scanner.num_detectors))
        .append("_")
        .append(std::to_string((int)scanner.field_of_view))
        .append("_")
        .append(std::to_string(i))
        .append(".ppm");

      writePpmHeader(sart_out_fname, scanner.num_pixels, scanner.num_pixels);
      double _max = *std::max_element(img.begin(), img.end());
      double _min = *std::min_element(img.begin(), img.end());
      std::cout<< "img "  << _min << " - " << _max << "\n";
      writePpmData(sart_out_fname, img, total_pixels, _max, _min);
    }
  }

  return 0;
}
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O3 -Wall -DNDEBUG -march=native -o projection_O3_native.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O0 -Wall -DNDEBUG -o projection_O0.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -lbenchmark -O0 -Wall -DNDEBUG -DINSTR_RDTSC -o projection_O0_rdtsc.o
#gcc -march=broadwell -Q --help=target -v

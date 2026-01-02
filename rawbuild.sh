g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O0 -Wall -DNDEBUG -o projection_O0.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O0 -Wall -DNDEBUG -DINSTR_RDTSC -o projection_O0_rdtsc.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O1 -Wall -DNDEBUG -o projection_O1.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O2 -Wall -DNDEBUG -o projection_O2.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O3 -Wall -DNDEBUG -o projection_O3.o
g++ src/main.cpp src/Scanner.cpp -I include/ -lfftw3 -lpthread -O3 -Wall -DNDEBUG -march=native -o projection_O3_native.o

all:
	qimcifa qimcifa_cl

qimcifa:
	g++ -std=c++11 src/qimcifa.cpp -o qimcifa -l pthread

qimcifa_cl:
	g++ -std=c++11 src/qimcifa_gpu.cpp -o qimcifa_cl -l pthread

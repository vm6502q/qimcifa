all: qimcifa qimcifa_cl

qimcifa:
	g++ -std=c++11 src/qimcifa.cpp -o qimcifa -l pthread

qimcifa_cl:
	g++ -std=c++11 -Iinclude src/common/oclengine.cpp src/qimcifa_cl.cpp -o qimcifa_cl -lOpenCL -lpthread

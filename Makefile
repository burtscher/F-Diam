CXX = g++
CXXFLAGS = -O3 -fopenmp
TARGET = fdiam
SRC = pdiam11.cpp

all: $(TARGET) setup

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $< -o $@

setup:
	bash setup.sh

clean:
	rm -f $(TARGET)

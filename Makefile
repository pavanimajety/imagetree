CXX = g++ 
CXXFLAGS = -g -std=c++11  `pkg-config --cflags --libs opencv`
LDFLAGS = -L/usr/local/lib -I /usr/local/include
TARGET = main
SOURCE = $(TARGET).cpp Imagetree.cpp Imagetree_lib.cpp

TARGET:
		$(CXX) -o  $(TARGET).exe $(SOURCE) $(CXXFLAGS) $(LDFLAGS)
		./$(TARGET).exe > tester.txt
preprocessor:
	     	$(CXX) -E $(SOURCE) $(CXXFLAGS) $(LDFLAGS) > preproc.txt
		

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET).exe tester.txt

CXX = g++ 
CXXFLAGS = -g -std=c++11 -O5 `pkg-config --cflags --libs opencv`
LDFLAGS = -L/usr/local/lib -I /usr/local/include
TARGET = main
SOURCE = $(TARGET).cpp Imagetree.cpp Imagetreel.cpp

TARGET:
		$(CXX) -o  $(TARGET) $(SOURCE) $(CXXFLAGS) $(LDFLAGS)
		./$(TARGET) > test.txt
preprocessor:
	     	$(CXX) -E $(SOURCE) $(CXXFLAGS) $(LDFLAGS) > preproc.txt
		

.PHONY: clean
clean:
	rm -f *.o *~ $(TARGET)

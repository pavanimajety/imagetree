/****************************************************************
 * This is a header file to create a quadtree out of an image
 * Author: Pavani Majety (Can be used for academic purposes)
 * **************************************************************/

#ifndef IMAGETREE_H
#define IMAGETREE_H
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "opencv2/core/types.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
using namespace std;
#define NUMBER_OF_QUADRANTS 4

/**************************************************
 * Abstract class for quadtrees
 * can store pixels - pixel is a struct that stores the Red,
 * Green and Blue channel values. The quadnode and quad leaf
 * are derived from the abstract class for quadtrees.
 **************************************************/

enum directions_t { NorthWest = 0, NorthEast, SouthEast, SouthWest };

class Point_xy {
 public:
  Point_xy(int x_, int y_) : x(x_), y(y_) {}
  int x;
  int y;
};

class Pixel {
 public:
  Pixel(int R_, int G_, int B_, int x_, int y_)
      : R(R_), G(G_), B(B_), x(x_), y(y_) {}
  Pixel(int R_, int G_, int B_) : R(R_), G(G_), B(B_), x(0), y(0) {}
  Pixel() : R(0), G(0), B(0), x(0), y(0) {}
  int x;
  int y;
  int R;
  int G;
  int B;
  
};

void operator<<(ostream &out, const Pixel  &P);

struct PixelDepthQuad {
  Pixel P;
  int depth;
  int quad;
  bool isLeaf;
  PixelDepthQuad* parentPixelDepth;
};
// typedef Pixel T;
struct Print_tree_parameters {
  int depth;
  int quad;
  int x;
  int y;
  int rows_image;
  int cols_image;
};
class Imagetree {
 public:
  ~Imagetree() {}
  Imagetree() {}
  virtual bool isLeaf() = 0;

  virtual bool isNode() { return !isLeaf(); }

  virtual int numberOfLeaves() = 0;

  virtual int numberOfNodes() = 0;

  virtual int numberOfSubTrees() = 0;

  virtual Pixel& value() = 0;

  virtual Imagetree*& son(directions_t dir) = 0;

  virtual void killAllSons() = 0;

  virtual void print_tree(vector<PixelDepthQuad*>& pixVector,
                          PixelDepthQuad* parent,
                          Print_tree_parameters params) = 0;

  void reconstructImageFromTree(vector<PixelDepthQuad*>& pixVector,
                                Mat* reconImage);

  virtual int heightOfTree() = 0;
};
class QuadLeaf : public Imagetree {
  Pixel pixel_val;

 public:
  // QuadLeaf specific functions

  QuadLeaf();
  QuadLeaf(Pixel pixel_in);
  // QuadLeaf(const QuadLeaf & ql_in);

  // overloaded functions
  bool isLeaf();
  int numberOfNodes();
  int numberOfLeaves();
  int numberOfSubTrees();
  Pixel& value();

  Imagetree* const& son(directions_t dir) const;
  Imagetree*& son(directions_t dir);
  void killAllSons();
  void print_tree(vector<PixelDepthQuad*>& pixVector, PixelDepthQuad* parent,
                  Print_tree_parameters params);
  void reconstructImageFromTree(vector<PixelDepthQuad*>& pixVector,
                                Mat* reconImage);

  int heightOfTree() { return 1; };
};
/*********************************************
 * A branching node of a quadtree, whose leaves contain
 * a value of type T
 * *********************************************/

class QuadNode : public Imagetree {
  // A quadnode can have either a node or a leaf as its
  // child nodes
 private:
  Imagetree* sons[NUMBER_OF_QUADRANTS];
  Pixel average_pixel;

 public:
  // Constructors are defined below
  QuadNode();
  QuadNode(Imagetree* inputs[]);
  QuadNode(Imagetree* in1, Imagetree* in2, Imagetree* in3, Imagetree* in4);
  ~QuadNode();

  // Overloaded functions from Imagetree
  bool isLeaf();
  int numberOfLeaves();
  int numberOfNodes();
  int numberOfSubTrees();
  Pixel& value();

  Imagetree* const& son(directions_t dir) const;

  Imagetree*& son(directions_t dir);

  bool killSon(directions_t dir);
  void killAllSons();

  int heightOfTree();

  void print_tree(vector<PixelDepthQuad*>& pixVector, PixelDepthQuad* parent,
                  Print_tree_parameters params);

  void reconstructImageFromTree(vector<PixelDepthQuad*>& pixVector,
                                Mat* reconImage);
};

class ComparePixVector {
 public:
  bool operator()(const PixelDepthQuad PDQ1, const PixelDepthQuad PDQ2) {
    // if(PDQ1.depth == PDQ2.depth){
    //   if((PDQ1.isLeaf)  && (!PDQ2.isLeaf)) return true;
    //   else if((PDQ2.isLeaf)  && (!PDQ1.isLeaf)) return false;
    //   else return false;
    // }
    // else
    return PDQ1.depth < PDQ2.depth;
  }
};

#endif
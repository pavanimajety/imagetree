#include "Imagetree.h"
#include <vector>
// A quadnode can have either a node or a leaf as its
// child nodes

/**************************************************************************************************************/
//    QUADNODE
//    QUADNODE
//    QUADNODE
/**************************************************************************************************************/

// Constructors are defined below
QuadNode::QuadNode() : average_pixel(Pixel(0, 0, 0)) {
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) sons[i] = nullptr;
}

QuadNode::QuadNode(Imagetree* inputs[]) : average_pixel(Pixel(0, 0, 0)) {
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) sons[i] = inputs[i];
}

QuadNode::QuadNode(Imagetree* in1, Imagetree* in2, Imagetree* in3,
                   Imagetree* in4)
    : average_pixel(Pixel(0, 0, 0)) {
  sons[0] = in1;
  sons[1] = in2;
  sons[2] = in3;
  sons[3] = in4;
}

QuadNode::~QuadNode() {
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i)
    if (sons[i]) delete sons[i];
}

/********************************************************************************************************************/
bool QuadNode::isLeaf() { return false; }

int QuadNode::numberOfLeaves() {
  int n = 0;
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) {
    if (sons[i] != nullptr) n += sons[i]->numberOfLeaves();
  }
  return n;
}

int QuadNode::numberOfNodes() {
  int n = 0;
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) {
    if (sons[i] != nullptr) n += sons[i]->numberOfNodes();
  }

  return n;
}

int QuadNode::numberOfSubTrees() {
  int n = 1;
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) {
    n += sons[i]->numberOfSubTrees();
  }
  return n;
}

Pixel& QuadNode::value() { return average_pixel; }
/*************************************************************************************************************/
Imagetree* const& QuadNode::son(directions_t dir) const { return sons[dir]; }

Imagetree*& QuadNode::son(directions_t dir) { return sons[dir]; }

bool QuadNode::killSon(directions_t dir) {
  if (sons[dir]->isLeaf()) delete sons[dir];
  sons[dir] = nullptr;

  return sons[dir] == nullptr;
}

void QuadNode::killAllSons() {
  for (int i = 0; i < NUMBER_OF_QUADRANTS; ++i) {
    sons[i]->killAllSons();
    delete sons[i];
  }
}

void QuadNode::print_tree(vector<PixelDepthQuad*>& pixVector,
                          PixelDepthQuad* parent,
                          Print_tree_parameters params) {
#if DEBUG
  std::cout << "Node: Pixel: [ " << average_pixel.R << "," << average_pixel.G
            << "," << average_pixel.B << "  ],\tdepth: " << params.depth
            << "\tx: " << params.x << "\ty: " << params.y << endl;
#endif
  PixelDepthQuad* ptr =
      new PixelDepthQuad{Pixel(average_pixel.R, average_pixel.G,
                               average_pixel.B, params.x, params.y),
                         params.depth, params.quad, 0, parent};

  pixVector.push_back(ptr);
  int sizeV = pixVector.size();
  PixelDepthQuad* newParent = pixVector[sizeV - 1];
  int differ = pow(2, params.depth + 1);
  int width = (params.cols_image / differ);
  int height = (params.rows_image / differ);

  Print_tree_parameters quad_params[4]({params, params, params, params});
  quad_params[0].depth = params.depth + 1;
  quad_params[0].quad = 0;
  quad_params[0].x = params.x;
  quad_params[0].y = params.y;
  sons[0]->print_tree(pixVector, newParent, quad_params[0]);

  quad_params[1].depth = params.depth + 1;
  quad_params[1].quad = 1;
  quad_params[1].x = params.x;
  quad_params[1].y = params.y + width;
  sons[1]->print_tree(pixVector, newParent, quad_params[1]);

  quad_params[2].depth = params.depth + 1;
  quad_params[2].quad = 2;
  quad_params[2].x = params.x + height;
  quad_params[2].y = params.y + width;
  sons[2]->print_tree(pixVector, newParent, quad_params[2]);

  quad_params[3].depth = params.depth + 1;
  quad_params[3].quad = 3;
  quad_params[3].x = params.x + height;
  quad_params[3].y = params.y;
  sons[3]->print_tree(pixVector, newParent, quad_params[3]);
}

int QuadNode::heightOfTree() {
  if (this == nullptr) {
    return 0;
  } else {
    return 1 + max(sons[0]->heightOfTree(),
                   max(sons[1]->heightOfTree(),
                       max(sons[2]->heightOfTree(), sons[3]->heightOfTree())));
  }
}


void Imagetree::reconstructImageFromTree(vector<PixelDepthQuad*>& pixVector,
                                         Mat* reconImage) {
  for (int i = 0; i < pixVector.size(); ++i) {
    if (pixVector[i]->isLeaf) {
      int factor = pow(2, pixVector[i]->depth);
// #ifdef DEBUG
//       // cout<<"Factor: "<<factor<<" Depth: "<<pixVector[i]->depth<<endl;
// #endif
      for (int row = pixVector[i]->P.x;
           row <
           pixVector[i]->P.x + (reconImage->rows / pow(2, pixVector[i]->depth));
           ++row) {
        for (int col = pixVector[i]->P.y;
             col < pixVector[i]->P.y +
                       (reconImage->cols / pow(2, pixVector[i]->depth));
             ++col) {
          reconImage->at<Vec3b>(row, col) =
              Vec3b(pixVector[i]->P.R, pixVector[i]->P.G, pixVector[i]->P.B);
        }
      }
    }
  }
}

// static fill_pixel(Mat* reconImage,)

// static QuadNode

/**************************************************************************************************************/
//    QUADLEAF
//    QUADLEAF
//    QUADLEAF
/**************************************************************************************************************/

QuadLeaf::QuadLeaf() : pixel_val(Pixel(0, 0, 0)) {}

QuadLeaf::QuadLeaf(Pixel pixel_in) { pixel_val = pixel_in; }

/**************************************************************************************************************/
bool QuadLeaf::isLeaf() { return true; }

int QuadLeaf::numberOfNodes() { return 0; }

int QuadLeaf::numberOfLeaves() { return 1; }

int QuadLeaf::numberOfSubTrees() { return 1; }

Pixel& QuadLeaf::value() { return pixel_val; }
/**************************************************************************************************************/
// const Pixel&  QuadLeaf::value() const{ return pixel_val; }

Imagetree* const& QuadLeaf::son(directions_t dir) const {
  throw std::domain_error("Not a QuadNode, doesn't have sons");
}

Imagetree*& QuadLeaf::son(directions_t dir) {
  throw std::domain_error("Not a QuadNode, doesn't have sons");
}

void QuadLeaf::killAllSons() {}
/*****************************************/

void QuadLeaf::print_tree(vector<PixelDepthQuad*>& pixVector,
                          PixelDepthQuad* parent,
                          Print_tree_parameters params) {
#if DEBUG
  std::cout << "Leaf: Pixel: [ " << pixel_val.R << "," << pixel_val.G << ","
            << pixel_val.B << "  ],\tdepth: " << params.depth
            << "\tx: " << params.x << "\ty: " << params.y << std::endl;
#endif
  PixelDepthQuad* ptr = new PixelDepthQuad{
      Pixel(pixel_val.R, pixel_val.G, pixel_val.B, params.x, params.y),
      params.depth, params.quad, 1, parent};
  pixVector.push_back(ptr);
  int sizeV = pixVector.size();
  PixelDepthQuad* newParent = pixVector[sizeV - 1];
}
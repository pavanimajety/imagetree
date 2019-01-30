#ifndef IMAGETREEL_H
#define IMAGETREEL_H

#include <iostream>
#include <string>
#include "opencv2/core/types.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "Imagetree.h"
#include <cmath>

// #Defining the color of the object in the frame
// # Color  = Blue
// # define the list of boundaries
// # Red:    ([17, 15, 100], [50, 56, 200])
// # Blue:   ([86, 31, 4], [220, 88, 50])
// # yellow: ([25, 146, 190], [62, 174, 250])
// # Gray:   ([103, 86, 65], [145, 133, 128])
// # Green:  ([29,86,6],[64,255,255])

Imagetree * calculateQuadtreeBasedOnMask(Mat *mask, Mat *input_image,
                                         Imagetree *itree);
Imagetree * insert(Mat *mask, Mat *image, Imagetree *itree, int depth);
void print(vector<double> &v);
Vec3b calculate_average_pixel(cv::Mat *image);
void clear_the_pixVector(vector<PixelDepthQuad *> &pixVectors);
double calculate_RMSE(Mat *Original, Mat * reconstructed);

//Requires: An image that is opened and is present in the memory
//Modifies: The pixVector
//Effects:  The pixVector stores the values of the nodes of the quadtree representation of the image
//          The PixelDepthQuad represents the characteristics of each node -
//          R,G,B Values of each node, the x,y values of each pixel, the depth, the quadrant number
void construct_RMSE_LinearQuadtree(vector<PixelDepthQuad *> &pixVectors, Mat * imOriginal);
double calculate_average_pixel_RMSE(Mat *imOriginal, Pixel &p_avg);
double calculate_average_pixel_RMSE_vec3b(Mat *imOriginal, Vec3b &p_avg);
int calculate_max_depth_of_quadtree(Mat * img);
void reconstructImageFromTree(vector<PixelDepthQuad*>& pixVector,
                                         Mat* reconImage);
#endif
#include <iostream>
#include <string>
#include <memory>
#include "Imagetree.h"
#include "opencv2/core/types.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <cstdlib>
#include <pthread.h>
#include "Imagetreel.h"


using namespace std;
using namespace cv;
/******************************************************************
 *  VERSION 2.0
 * AUTHOR: PAVANI MAJETY
 *  
 * CHANGELOG:   1. USES RMSE TO CALCULATE THE QUADTREE
 *              2. USES A LINEAR QUADTREE   
 *              3. USES MULTITHREADING
 *
 * *****************************************************************/
#define CAM 1
#define DEBUG 0
#define DEBUGCAM 0
#define PRINTTREE 1
#define IMAGEE 0

// #Defining the color of the object in the frame
// # Color  = Blue
// # define the list of boundaries
// # Red:    ([17, 15, 100], [50, 56, 200])
// # Blue:   ([86, 31, 4], [220, 88, 50])
// # yellow: ([25, 146, 190], [62, 174, 250])
// # Gray:   ([103, 86, 65], [145, 133, 128])
// # Green:  ([29,86,6],[64,255,255])



int main(int argc, char const *argv[]) {
  Imagetree *imL, *imN;


// 
  Mat imOriginal = imread("Images/facebook.png", CV_LOAD_IMAGE_COLOR);
  // Rect roi(8, 8, imOriginal_.rows - 16, imOriginal_.cols - 16);
  // imOriginal_ = imOriginal_(roi);
  Mat hsvImg;
  Mat threshImg;
  Mat *mask = &threshImg;
  Mat &img = hsvImg;
  if (!imOriginal.data)  // Check for
  {
    cout << "Could not open or find the image" << std::endl;
    return -1;
  } else {
    cout << "opened image" << endl;
  }


  // Setting the threshold for white color
  Scalar lowRange(0, 0, 200);
  Scalar highRange(200, 0, 255);

  Mat imOriginal_ = Mat(512, 512, CV_8UC3); //required size for quadtree decomposition - add the required observations.
/*****STEP 1: Preprocessing of the image -- DO NOT CALCULATE THE MASK FOR THIS VERSION   ***********************/
/* This version uses the RMSE error for calculation */

  resize(imOriginal, imOriginal, imOriginal_.size(), 0, 0, INTER_AREA);
  cout << imOriginal.rows << " " << imOriginal.cols << endl;

  cvtColor(imOriginal, hsvImg, CV_BGR2HSV);
  imshow("Original Image", imOriginal);
  imwrite("Images/original_image.png", imOriginal);

  char charCheckForEscKey = cv::waitKey(1);
//   //thresholding the image
//   inRange(hsvImg, lowRange, highRange, threshImg);
//   imwrite("Images/Mask.png", threshImg);  inRange(hsvImg, lowRange, highRange, threshImg);

/**********************CHECK WHAT HAPPENS WITHOUT GAUSSIAN BLUR ***********************/
/* Adds the coefficients with decimals, and blurs out the edges 
*/

  // GaussianBlur(threshImg, threshImg, Size(3, 3), 0);
  // dilate(threshImg, threshImg, 0);
  // erode(threshImg, threshImg, 0);
  // threshold(threshImg, threshImg, 0, 255, THRESH_BINARY);


/****STEP 2:  Construction of the Linear Quadtree *************/

  // is this required ?
  std::unique_ptr<Imagetree> itree(new QuadNode(0, 0, 0, 0));
  int depth =  calculate_max_depth_of_quadtree(&imOriginal);
  vector<PixelDepthQuad*> pixVectors;//((1+4*depth));
  construct_RMSE_LinearQuadtree(pixVectors, &imOriginal);


/****STEP 3:  Reconstruction of the image *************/
  Mat reconstrutedImage = Mat::zeros(threshImg.rows, threshImg.cols, CV_8UC3);

 // itree->reconstructImageFromTree(pixVectors, &reconstrutedImage);
//   imshow("Reconstructed Image", reconstrutedImage);
//   cout << "Root Mean Squared Error: "<<calculate_RMSE_ERROR(&imOriginal,&reconstrutedImage)<<endl;
//   imwrite("Images/recostructed_img.png", reconstrutedImage);
//   waitKey(1);


//   itree->killAllSons();

  return 0;
}
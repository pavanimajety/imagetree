#include <iostream>
#include <string>
#include "Imagetree.h"
#include "Imagetreel.h"
#include "opencv2/core/types.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace cv;

#define CAM 1
#define DEBUG 0
#define DEBUGCAM 0
#define PRINTTREE 1
#define IMAGEE 0


int main(int argc, char const *argv[]) {
  Imagetree *imL, *imN, *itree;

#if CAM
#if !IMAGEE
  cv::VideoCapture capWebCam;
  capWebCam.open(-1);

  if (capWebCam.isOpened() == false) {
    cout << "error: the webcam is not configured" << endl;
    return 0;
  } else {
    cout << "opened webcam" << endl;
  }

  Mat imOriginal;
  Mat hsvImg;
  Mat threshImg;
  Mat *mask = &threshImg;
  Mat &img = hsvImg;
  char charCheckForEscKey = 0;
#endif  //! IMAGE
#endif  // CAM

#if IMAGEE
  Mat imOriginal = imread("original_image.png", CV_LOAD_IMAGE_COLOR);
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
#endif

  // Setting the threshold for white color
  Scalar lowRange(0, 0, 200);
  Scalar highRange(200, 0, 255);
  // Scalar lowRange(17, 15, 100);
  // Scalar highRange(50, 56, 200);

  Mat imOriginal_ = Mat(512, 512, CV_8UC3);
#if CAM
  while (capWebCam.read(imOriginal)) {
  resize(imOriginal, imOriginal, imOriginal_.size(), 0, 0, INTER_AREA);
  cout << imOriginal.rows << " " << imOriginal.cols << endl;

  cvtColor(imOriginal, hsvImg, CV_BGR2HSV);
  imshow("Original Image", imOriginal);
  imwrite("Images/original_image.png", imOriginal);

  char charCheckForEscKey = cv::waitKey(1);

  inRange(hsvImg, lowRange, highRange, threshImg);
  imwrite("Images/Mask.png", threshImg);
  // GaussianBlur(threshImg, threshImg, Size(3, 3), 0);
  // dilate(threshImg, threshImg, 0);
  // erode(threshImg, threshImg, 0);
  // threshold(threshImg, threshImg, 0, 255, THRESH_BINARY);

#endif
#if DEBUGCAM
  imshow("Original Image", imOriginal);
  charCheckForEscKey = cv::waitKey(1);

  imshow("Mask", threshImg);
  imwrite("Images/Mask.png", threshImg);
  waitKey();
#endif

#if !CAM
  itree = new QuadNode(0, 0, 0, 0);
  Mat C = Mat::zeros(16, 16, CV_64FC1);
  Mat *mask = &C;
  mask->at<double>(5, 5) = 255.0d;
  mask->at<double>(5, 6) = 255.0d;
  mask->at<double>(6, 5) = 255.0d;
  mask->at<double>(6, 6) = 255.0d;
  Mat img = Mat::zeros(16, 16, CV_8UC3);
  img.at<Vec3b>(5, 5) = Vec3b(0, 0, 255);  // BGR Format
  img.at<Vec3b>(5, 6) = Vec3b(0, 255, 0);
  img.at<Vec3b>(6, 5) = Vec3b(255, 0, 0);
  img.at<Vec3b>(6, 6) = Vec3b(30, 0, 30);

  imshow("Fake Image", img);
  waitKey();
  imshow("Fake Mask", C);
  waitKey();

  cout << "C:\n" << C << endl;

  itree = calculateQuadtreeBasedOnMask(mask, &img, itree);

#if PRINTTREE
  Mat reconstrutedImage = Mat::zeros(mask->rows, mask->cols, CV_8UC3);
  int levels = itree->heightOfTree();
  cout << "Depth of tree: " << levels << endl;
  vector<PixelDepthQuad *> pixVectors;

  Print_tree_parameters params{0, 0, 0, 0, mask->rows, mask->cols};
  itree->print_tree(pixVectors, 0, params);
  // sort(pixVectors.begin(),pixVectors.end(),comparePixVector());
#if DEBUG
  for (int i = 1; i < pixVectors.size(); ++i) {
    cout << "i: " << i << "\t[ " << pixVectors[i]->P.B << ","
         << pixVectors[i]->P.G << "," << pixVectors[i]->P.R
         << " ]\tQ:" << pixVectors[i]->quad << "\tD:" << pixVectors[i]->depth
         << "\tisLeaf: " << pixVectors[i]->isLeaf
         << "\t Parent: " << pixVectors[i]->parentPixelDepth
         << "\tMe:" << pixVectors[i] << endl;
  }
#endif  // DEBUG
  itree->reconstructImageFromTree(pixVectors, &reconstrutedImage);
  // cout<<reconstrutedImage;
  // cvtColor(img, cv2.COLOR_BGR2RGB)
  imshow("Reconstructed Image", reconstrutedImage);
  waitKey(1);
  // clear_the_pixVector(pixVectors);
#endif  // PRINTREE
#endif  // CAM
// #endif
#if CAM
  itree = new QuadNode(0, 0, 0, 0);

  calculateQuadtreeBasedOnMask(&threshImg, &imOriginal, itree);

#if PRINTTREE
  Mat reconstrutedImage = Mat::zeros(threshImg.rows, threshImg.cols, CV_8UC3);
  // int levels = itree->heightOfTree();
  // cout << "New Frame" << endl;
  // cout << "Depth of tree: " << levels << endl;
  vector<PixelDepthQuad *> pixVectors;
  Print_tree_parameters params{0, 0, 0, 0, mask->rows, mask->cols};
  itree->print_tree(pixVectors, 0, params);

#if DEBUG
  for (int i = 1; i < pixVectors.size(); ++i) {
    cout << "i: " << i << "\t[ " << pixVectors[i]->P.B << ","
         << pixVectors[i]->P.G << "," << pixVectors[i]->P.R
         << " ]\tQ:" << pixVectors[i]->quad << "\tD:" << pixVectors[i]->depth
         << "\tisLeaf: " << pixVectors[i]->isLeaf
         << "\t Parent: " << pixVectors[i]->parentPixelDepth
         << "\tMe:" << pixVectors[i] << endl;
  }
#endif  // DEBUG
#endif  // PRINTTREE
  itree->reconstructImageFromTree(pixVectors, &reconstrutedImage);
  // cout<<reconstrutedImage;
  // cvtColor(reconstrutedImage,reconstrutedImage,COLOR_HSV2BGR);
  imshow("Reconstructed Image", reconstrutedImage);
  cout << "Root Mean Squared Error: "<<calculate_RMSE_ERROR(&imOriginal,&reconstrutedImage)<<endl;
  imwrite("Images/recostructed_img.png", reconstrutedImage);
  waitKey(1);
  }
// 
  capWebCam.release();
#endif  // CAM

  itree->killAllSons();
  delete itree;

  return 0;
}

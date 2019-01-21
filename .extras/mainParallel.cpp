#include <iostream>
#include <string>
#include "Imagetree.h"
#include "opencv2/core/types.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace std;
using namespace cv;
/******************************************************************
 * Modifications required ---
 * Move insert to Imagetree.h
 * change the vector pushbacks to be modified during inserts
 *  ------------- this way, extra computation during printing tree will
 *  ------------- be avoided
 *
 *
 * *****************************************************************/
#define CAM 1
#define DEBUG 0
#define DEBUGCAM 0
#define PRINTTREE 1
#define IMAGEE 1
// #Defining the color of the object in the frame
// # Color  = Blue
// # define the list of boundaries
// # Red:    ([17, 15, 100], [50, 56, 200])
// # Blue:   ([86, 31, 4], [220, 88, 50])
// # yellow: ([25, 146, 190], [62, 174, 250])
// # Gray:   ([103, 86, 65], [145, 133, 128])
// # Green:  ([29,86,6],[64,255,255])

Imagetree *calculateQuadtreeBasedOnMask(Mat *mask, Mat *input_image,
                                        Imagetree *itree);
Imagetree *insert(Mat *mask, Mat *image, Imagetree *itree, int depth);
void print(vector<double> &v);
Vec3b calculate_average_pixel(cv::Mat *image);
void clear_the_pixVector(vector<PixelDepthQuad *> &pixVectors);

int main(int argc, char const *argv[]) {
  Imagetree *imL, *imN, *itree;

#if CAM
  cv::VideoCapture capWebCam;
  capWebCam.open(-1);

  if (capWebCam.isOpened() == false) {
    cout << "error: the webcam is not configured" << endl;
    return 0;
  } else {
    cout << "opened webcam" << endl;
  }
#endif

  Mat imOriginal;
  Mat hsvImg;
  Mat threshImg;
  Mat *mask = &threshImg;
  Mat &img = hsvImg;
  char charCheckForEscKey = 0;
  // #ifdef IMAGEE
  //   Mat imOriginal__;  // = imread("trial_cp.png", CV_LOAD_IMAGE_COLOR);
  //   // Rect roi(8, 8, imOriginal_.rows - 16, imOriginal_.cols - 16);
  //   // imOriginal_ = imOriginal_(roi);

  //   if (!imOriginal__.data)  // Check for
  //   {
  //     cout << "Could not open or find the image" << std::endl;
  //     return -1;
  //   } else {
  //     cout << "opened image" << endl;
  //   }
  // #endif

  // Setting the threshold for white color
  Scalar lowRange(0, 0, 0);
  Scalar highRange(0, 0, 255);
  // Scalar lowRange(17, 15, 100);
  // Scalar highRange(50, 56, 200);

  Mat imOriginal_ = Mat(512, 512, CV_8UC3);

  while (capWebCam.read(imOriginal)) {
    resize(imOriginal, imOriginal, imOriginal_.size(), 0, 0, INTER_AREA);
    cout << imOriginal.rows << " " << imOriginal.cols << endl;

    cvtColor(imOriginal, hsvImg, CV_BGR2HSV);
    imshow("Original Image", imOriginal);
    charCheckForEscKey = cv::waitKey(50);

    inRange(hsvImg, lowRange, highRange, threshImg);
    GaussianBlur(threshImg, threshImg, Size(3, 3), 0);
    dilate(threshImg, threshImg, 0);
    erode(threshImg, threshImg, 0);
    threshold(threshImg, threshImg, 125, 255, THRESH_BINARY);

// #if DEBUGCAM
//     imshow("Original Image", imOriginal);
//     charCheckForEscKey = cv::waitKey(10);

//     imshow("Mask", threshImg);
//     imwrite("mask.png", threshImg);
//     waitKey(10);
// #endif

// #if !CAM
//     itree = new QuadNode(0, 0, 0, 0);
//     Mat C = Mat::zeros(16, 16, CV_64FC1);
//     Mat *mask = &C;
//     mask->at<double>(5, 5) = 255.0d;
//     mask->at<double>(5, 6) = 255.0d;
//     mask->at<double>(6, 5) = 255.0d;
//     mask->at<double>(6, 6) = 255.0d;
//     Mat img = Mat::zeros(16, 16, CV_8UC3);
//     img.at<Vec3b>(5, 5) = Vec3b(0, 0, 255);  // BGR Format
//     img.at<Vec3b>(5, 6) = Vec3b(0, 255, 0);
//     img.at<Vec3b>(6, 5) = Vec3b(255, 0, 0);
//     img.at<Vec3b>(6, 6) = Vec3b(30, 0, 30);

//     imshow("Fake Image", img);
//     waitKey(10);
//     imshow("Fake Mask", C);
//     waitKey(10);

//     cout << "C:\n" << C << endl;

//     itree = calculateQuadtreeBasedOnMask(mask, &img, itree);

// #if PRINTTREE
//     Mat reconstrutedImage = Mat::zeros(mask->rows, mask->cols, CV_8UC3);
//     int levels = itree->heightOfTree();
//     cout << "Depth of tree: " << levels << endl;
//     vector<PixelDepthQuad *> pixVectors;

//     Print_tree_parameters params{0, 0, 0, 0, mask->rows, mask->cols};
//     itree->print_tree(pixVectors, 0, params);
//     // sort(pixVectors.begin(),pixVectors.end(),comparePixVector());
// #if DEBUG
//     for (int i = 1; i < pixVectors.size(); ++i) {
//       cout << "i: " << i << "\t[ " << pixVectors[i]->P.B << ","
//            << pixVectors[i]->P.G << "," << pixVectors[i]->P.R
//            << " ]\tQ:" << pixVectors[i]->quad << "\tD:" <<
//            pixVectors[i]->depth
//            << "\tisLeaf: " << pixVectors[i]->isLeaf
//            << "\t Parent: " << pixVectors[i]->parentPixelDepth
//            << "\tMe:" << pixVectors[i] << endl;
//     }
// #endif
//     itree->reconstructImageFromTree(pixVectors, &reconstrutedImage);
//     // cout<<reconstrutedImage;
//     // cvtColor(img, cv2.COLOR_BGR2RGB)
//     imshow("Reconstructed Image", reconstrutedImage);
//     waitKey(10);
//     // clear_the_pixVector(pixVectors);
// #endif

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
    // for (int i = 1; i < pixVectors.size(); ++i) {
    //   cout << "i: " << i << "\t[ " << pixVectors[i]->P.B << ","
    //        << pixVectors[i]->P.G << "," << pixVectors[i]->P.R
    //        << " ]\tQ:" << pixVectors[i]->quad << "\tD:" << pixVectors[i]->depth
    //        << "\tisLeaf: " << pixVectors[i]->isLeaf
    //        << "\t Parent: " << pixVectors[i]->parentPixelDepth
    //        << "\tMe:" << pixVectors[i] << endl;
    // }
#endif  // DEBUG
#endif  // PRINTTREE
    itree->reconstructImageFromTree(pixVectors, &reconstrutedImage);
    // cout<<reconstrutedImage;
    // cvtColor(reconstrutedImage,reconstrutedImage,COLOR_HSV2BGR);
    imshow("Reconstructed Image", reconstrutedImage);
    imwrite("recostructed_img.png", reconstrutedImage);
    waitKey(10);
  }

  capWebCam.release();
#endif //CAM

  itree->killAllSons();
  delete itree;

  return 0;
}

Imagetree *calculateQuadtreeBasedOnMask(Mat *mask, Mat *image,
                                        Imagetree *itree) {
  Mat *temp_mask = mask;
  int depth = 0;

  return insert(mask, image, itree, depth);
}

Imagetree *insert(Mat *mask, Mat *image, Imagetree *itree, int depth) {
  vector<double> sum(4);
  for (int i = 0; i < (mask->rows / 2); ++i) {
    for (int j = 0; j < (mask->cols / 2); ++j) {
      if (mask->at<double>(i, j)) ++sum[0];
    }
  }

  for (int i = 0; i < (mask->rows / 2); ++i) {
    for (int j = (mask->cols / 2); j < (mask->cols); ++j) {
      if (mask->at<double>(i, j)) ++sum[1];
    }
  }

  for (int i = mask->rows / 2; i < (mask->rows); ++i) {
    for (int j = mask->cols / 2; j < (mask->cols); ++j) {
      if (mask->at<double>(i, j)) ++sum[2];
    }
  }
  for (int i = mask->rows / 2; i < (mask->rows); ++i) {
    for (int j = 0; j < (mask->cols / 2); ++j) {
      if (mask->at<double>(i, j)) ++sum[3];
    }
  }

#if DEBUG
  print(sum);
#endif

  cv::Rect NW(0, 0, mask->cols / 2, mask->rows / 2);
  cv::Rect NE(mask->cols / 2, 0, mask->cols / 2, mask->rows / 2);
  cv::Rect SE(mask->rows / 2, mask->cols / 2, mask->rows / 2, mask->cols / 2);
  cv::Rect SW(0, mask->rows / 2, mask->cols / 2, mask->rows / 2);

  cv::Mat ImageNW(*image, NW);
  cv::Mat MaskNW(*mask, NW);
  cv::Mat ImageNE(*image, NE);
  cv::Mat MaskNE(*mask, NE);
  cv::Mat ImageSE(*image, SE);
  cv::Mat MaskSE(*mask, SE);
  cv::Mat ImageSW(*image, SW);
  cv::Mat MaskSW(*mask, SW);
#if DEBUGCAM
  imshow("NW", ImageNW);
  waitKey();
  imshow("NE", ImageNE);
  waitKey();
  imshow("SE", ImageSE);
  waitKey();
  imshow("SW", ImageSW);
  waitKey();
#endif

  if (sum[0] > 1) {
    // imshow("NW", ImageNW);
    // auto charCheckForEscKey = cv::waitKey(1);
    auto tmp_ptr = new QuadNode(0, 0, 0, 0);
    auto avg_pixel = calculate_average_pixel(&ImageNW);
    tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
    itree->son(NorthWest) = tmp_ptr;
    insert(&MaskNW, &ImageNW, itree->son(NorthWest), depth + 1);
  } else if (sum[0] == 1) {
    if ((ImageNW.rows > 1) || (ImageNW.cols > 1)) {
      auto tmp_ptr = new QuadNode(0, 0, 0, 0);
      auto avg_pixel = calculate_average_pixel(&ImageNW);
      tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
      itree->son(NorthWest) = tmp_ptr;
      insert(&MaskNW, &ImageNW, itree->son(NorthWest), depth + 1);
    } else {
      auto curr_pixel = ImageNW.at<Vec3b>(0, 0);
      auto tmp_ptr =
          new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
      itree->son(NorthWest) = tmp_ptr;
#if DEBUG
      cout << "Current Depth in NW is: " << depth << endl;
      cout << "Current pixel value: " << curr_pixel << "Size: " << ImageNW.rows
           << " " << ImageNW.cols << endl;
#endif
    }
  } else {
#if DEBUG
    cout << "Reached Null in NW at depth" << depth << endl;
#endif
    auto curr_pixel = calculate_average_pixel(&ImageNW);

    auto tmp_ptr =
        new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
    itree->son(NorthWest) = tmp_ptr;
    // return;
  }

  if (sum[1] > 1) {
    // imshow("NE", ImageNE);
    // auto charCheckForEscKey = cv::waitKey(1);
    auto tmp_ptr = new QuadNode(0, 0, 0, 0);
    auto avg_pixel = calculate_average_pixel(&ImageNE);
    tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
    itree->son(NorthEast) = tmp_ptr;
    insert(&MaskNE, &ImageNE, itree->son(NorthEast), depth + 1);
  } else if (sum[1] == 1) {
    if ((ImageNE.rows > 1) || (ImageNE.cols > 1)) {
      auto tmp_ptr = new QuadNode(0, 0, 0, 0);
      auto avg_pixel = calculate_average_pixel(&ImageNE);
      tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
      itree->son(NorthEast) = tmp_ptr;
      insert(&MaskNE, &ImageNE, itree->son(NorthEast), depth + 1);
    } else {
      auto curr_pixel = ImageNE.at<Vec3b>(0, 0);
      auto tmp_ptr =
          new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
      itree->son(NorthEast) = tmp_ptr;
#if DEBUG
      cout << "Current Depth in NE is: " << depth << endl;
      cout << "Current pixel value: " << curr_pixel << "Size: " << ImageNE.rows
           << " " << ImageNE.cols << endl;

#endif
    }
  } else {
    auto curr_pixel = calculate_average_pixel(&ImageNE);
#if DEBUG
    cout << "Reached Null in NE at depth" << depth << endl;
    cout << "Current Pixel: " << curr_pixel << endl;

#endif

    auto tmp_ptr =
        new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
    itree->son(NorthEast) = tmp_ptr;
  }
  if (sum[2] > 1) {
    // imshw("SE", ImageSE);
    // auto charCheckForEscKey = cv::waitKey(1);
    auto tmp_ptr = new QuadNode(0, 0, 0, 0);
    auto avg_pixel = calculate_average_pixel(&ImageSE);
    tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
    itree->son(SouthEast) = tmp_ptr;
    insert(&MaskSE, &ImageSE, itree->son(SouthEast), depth + 1);
  } else if (sum[2] == 1) {
    if ((ImageSE.rows > 1) || (ImageSE.cols > 1)) {
      auto tmp_ptr = new QuadNode(0, 0, 0, 0);
      auto avg_pixel = calculate_average_pixel(&ImageSE);
      tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
      itree->son(SouthEast) = tmp_ptr;
      insert(&MaskSE, &ImageSE, itree->son(SouthEast), depth + 1);
    } else {
      auto curr_pixel = ImageSE.at<Vec3b>(0, 0);
      auto tmp_ptr =
          new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
      itree->son(SouthEast) = tmp_ptr;
#if DEBUG
      cout << "Current Depth in SE is" << depth << endl;
      cout << "Current pixel value" << curr_pixel << "Size: " << ImageSE.rows
           << " " << ImageSE.cols << endl;
#endif
    }
  } else {
    auto curr_pixel = calculate_average_pixel(&ImageSE);
#if DEBUG
    cout << "Reached Null in SE at depth" << depth << endl;
    cout << "Current Pixel: " << curr_pixel << endl;
#endif
    auto tmp_ptr =
        new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
    itree->son(SouthEast) = tmp_ptr;
  }
  if ((sum[3] > 1)) {
    // imshow("SW", ImageSW);
    // auto charCheckForEscKey = cv::waitKey(1);
    auto tmp_ptr = new QuadNode(0, 0, 0, 0);
    auto avg_pixel = calculate_average_pixel(&ImageSW);
    tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
    itree->son(SouthWest) = tmp_ptr;
    insert(&MaskSW, &ImageSW, itree->son(SouthWest), depth + 1);
  } else if (sum[3] == 1) {
    if ((ImageSW.rows > 1) || (ImageSW.cols > 1)) {
      auto tmp_ptr = new QuadNode(0, 0, 0, 0);
      auto avg_pixel = calculate_average_pixel(&ImageSW);
      tmp_ptr->value() = Pixel(avg_pixel[0], avg_pixel[1], avg_pixel[2]);
      itree->son(SouthWest) = tmp_ptr;
      insert(&MaskSW, &ImageSW, itree->son(SouthWest), depth + 1);
    } else {
      auto curr_pixel = ImageSW.at<Vec3b>(0, 0);
      auto tmp_ptr =
          new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
      itree->son(SouthWest) = tmp_ptr;

#if DEBUG
      cout << "Current Depth in SW is" << depth << endl;
      cout << "Current pixel value" << curr_pixel << "Size: " << ImageSW.rows
           << " " << ImageSW.cols << endl;
#endif
    }
  } else {
    auto curr_pixel = calculate_average_pixel(&ImageSW);

    auto tmp_ptr =
        new QuadLeaf(Pixel(curr_pixel[0], curr_pixel[1], curr_pixel[2]));
#if DEBUG
    cout << "Reached Null in SW at depth" << depth << endl;
    cout << "Current Pixel: " << curr_pixel << endl;

#endif
    itree->son(SouthWest) = tmp_ptr;
  }

  return itree;
}

void print(vector<double> &v) {
  for (int i = 0; i < v.size(); ++i) {
    cout << v[i] << " ";
  }
  cout << endl;
}
// #if DEBUGCAM
// int x = 0;
// #endif
Vec3b calculate_average_pixel(cv::Mat *image) {
  vector<double> sum(3);
  sum[0] = 0;
  sum[1] = 0;
  sum[2] = 0;
  // #if DEBUGCAM
  //   imshow(to_string(++x), *image);
  //   auto charCheckForEscKey = cv::waitKey(1000);
  // #endif
  int total_values = 0;

  for (int i = 0; i < (image->rows); ++i) {
    for (int j = 0; j < (image->cols); ++j) {
      auto pixel = image->at<Vec3b>(i, j);
      sum[0] += pixel[0];
      sum[1] += pixel[1];
      sum[2] += pixel[2];
      ++total_values;
    }
  }
  sum[0] = sum[0] / total_values;
  sum[1] = sum[1] / total_values;
  sum[2] = sum[2] / total_values;

  Vec3b final_sum;
  final_sum[0] = sum[0];
  final_sum[1] = sum[1];
  final_sum[2] = sum[2];
  return final_sum;
}
void clear_the_pixVector(vector<PixelDepthQuad *> &pixVectors) {
  for (int i = 0; i < pixVectors.size(); ++i) {
    delete pixVectors[i];
  }
}
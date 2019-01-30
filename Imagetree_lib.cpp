
#include "Imagetree_lib.h"
#include "Imagetree.h"

Imagetree *calculateQuadtreeBasedOnMask(Mat *mask, Mat *image,
                                        Imagetree *itree) {
  Mat *temp_mask = mask;
  int depth = 0;

  return insert(mask, image, itree, depth);
}

Imagetree *insert(Mat *mask, Mat *image, Imagetree *itree, int depth) {
  int WIDTH = mask->cols;
  int HEIGHT = mask->rows;
  int HWIDTH = mask->cols / 2;
  int HHEIGHT = mask->rows / 2;
  vector<double> sum = {0, 0, 0, 0};
#if DEBUG
  cout << "Quad WIDTH: " << HWIDTH << " Quad Height: " << HHEIGHT << endl;
#endif
  int x = 0;
  int y = 0;
  cv::Rect NW(0, 0, HWIDTH, HHEIGHT);
  for (int i = x; i < x + HHEIGHT; ++i) {
    for (int j = y; j < y + HWIDTH; ++j) {
      if (mask->at<uchar>(i, j) == 255) ++sum[0];
    }
  }
#if DEBUG
  cout << "Sum [0]:" << sum[0] << " x: " << x << " y: " << y << endl;
#endif

  /*****************************************************************/
  x = 0;
  y = HWIDTH;
  cv::Rect NE(mask->cols / 2, 0, mask->cols / 2, mask->rows / 2);

  for (int i = x; i < x + HHEIGHT; ++i) {
    for (int j = y; j < y + HWIDTH; ++j) {
      if (mask->at<uchar>(i, j) == 255) ++sum[1];
    }
  }
#if DEBUG
  cout << "Sum [1]:" << sum[1] << " x: " << x << " y: " << y << endl;
#endif
  x = HHEIGHT;
  y = HWIDTH;
  cv::Rect SE(mask->rows / 2, mask->cols / 2, mask->rows / 2, mask->cols / 2);
  for (int i = x; i < x + HHEIGHT; ++i) {
    for (int j = y; j < y + HWIDTH; ++j) {
      if (mask->at<uchar>(i, j) == 255) ++sum[2];
    }
  }
#if DEBUG
  cout << "Sum [2]:" << sum[2] << " x: " << x << " y: " << y << endl;
#endif
  x = HHEIGHT;
  y = 0;
  Rect SW(0, mask->rows / 2, mask->cols / 2, mask->rows / 2);

  for (int i = x; i < x + HHEIGHT; ++i) {
    for (int j = y; j < y + HWIDTH; ++j) {
      if (mask->at<uchar>(i, j) == 255) ++sum[3];
    }
  }
#if DEBUG
  cout << "Sum [3]:" << sum[3] << " x: " << x << " y: " << y << endl;
#endif

#if DEBUG
  print(sum);
#endif

  cv::Mat ImageNW(*image, NW);
  cv::Mat MaskNW(*mask, NW);
  cv::Mat ImageNE(*image, NE);
  cv::Mat MaskNE(*mask, NE);
  cv::Mat ImageSE(*image, SE);
  cv::Mat MaskSE(*mask, SE);
  cv::Mat ImageSW(*image, SW);
  cv::Mat MaskSW(*mask, SW);
#if DEBUGCAM
  // imshow("NW", ImageNW);
  // waitKey();
  // imshow("NE", ImageNE);
  // waitKey();
  // imshow("SE", ImageSE);
  // waitKey();
  // imshow("SW", ImageSW);
  // waitKey();
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

double calculate_RMSE(Mat *original, Mat *reconstructed) {
  double error = 0;
  // #if DEBUGCAM
  //   imshow(to_string(++x), *image);
  //   auto charCheckForEscKey = cv::waitKey(1000);
  // #endif
  int total_values = 0;

  // calculate the
  for (int i = 0; i < (reconstructed->rows); ++i) {
    for (int j = 0; j < (reconstructed->cols); ++j) {
      auto pixel = reconstructed->at<Vec3b>(i, j) - original->at<Vec3b>(i, j);
      error += pixel[0] * pixel[0] + pixel[1] * pixel[1] + pixel[2] * pixel[2];
    }
  }
  error = error / (reconstructed->rows * reconstructed->cols);

  return error;
}

double calculate_average_pixel_RMSE(Mat *original, Pixel &p_avg) {
  double error = 0;
  // #if DEBUGCAM
  //   imshow(to_string(++x), *image);
  //   auto charCheckForEscKey = cv::waitKey(1000);
  // #endif
  int total_values = 0;

  // calculate the
  for (int i = 0; i < (original->rows); ++i) {
    for (int j = 0; j < (original->cols); ++j) {
      vector<int> pixel = {p_avg.R - original->at<Vec3b>(i, j)[0],
                           p_avg.G - original->at<Vec3b>(i, j)[1],
                           p_avg.B - -original->at<Vec3b>(i, j)[2]};
      error += pixel[0] * pixel[0] + pixel[1] * pixel[1] + pixel[2] * pixel[2];
    }
  }
  error = error / (original->rows * original->cols);

  return error;
}

double calculate_average_pixel_RMSE_vec3b(Mat *original, Vec3b &p_avg) {
  double error = 0;
  // #if DEBUGCAM
  //   imshow(to_string(++x), *image);
  //   auto charCheckForEscKey = cv::waitKey(1000);
  // #endif
  int total_values = 0;

  // calculate the
  for (int i = 0; i < (original->rows); ++i) {
    for (int j = 0; j < (original->cols); ++j) {
      auto pixel = p_avg - original->at<Vec3b>(i, j);
      error += pixel[0] * pixel[0] + pixel[1] * pixel[1] + pixel[2] * pixel[2];
    }
  }
  error = error / (original->rows * original->cols);

  return error;
}
void construct_RMSE_LinearQuadtree(vector<PixelDepthQuad *> &pixVectors,
                                   Mat *imOriginal) {
  // get the number of levels
  int DEPTH = calculate_max_depth_of_quadtree(imOriginal);

  // calculate the average pixel value
  // Add the first PixelDepthQuad for representing the entire image
  auto pix_whole = calculate_average_pixel(imOriginal);
  // storing the first value of the pixVector storing the average pixel of the
  // entire value
  PixelDepthQuad p0{Pixel{pix_whole[0], pix_whole[1], pix_whole[2], 0, 0}, 0, 0,
                    0, 0};
  pixVectors.push_back(&p0);
  // calculate the RMSE
  auto rmse = calculate_average_pixel_RMSE(imOriginal, pixVectors[0]->P);
  cout << rmse << " is the rmse error" << endl;
  // if the RMSE is greater than the required keep subdividing
  int WIDTH = imOriginal->cols;
  int HEIGHT = imOriginal->rows;
  int HWIDTH = imOriginal->cols / 2;
  int HHEIGHT = imOriginal->rows / 2;
  int depth = 1;
  vector<Mat> vec_mat;
  vec_mat.push_back(*imOriginal);
  int matrix_index = 0;
  while (rmse > 0.01) {
    // divide into four quadrants, calculate the average pixel and subdivide
    // further if necessary.
    vector<Vec3b> avgPixQ(4);

    Rect NW(0, 0, HWIDTH, HHEIGHT);
    Rect NE(0, HWIDTH, HWIDTH, HHEIGHT);
    Rect SE(HHEIGHT, 0, HWIDTH, HHEIGHT);
    Rect SW(HHEIGHT, HWIDTH, HWIDTH, HHEIGHT);
    Mat quadImages[4];
    auto my_current_index = matrix_index;
    quadImages[0] = Mat(vec_mat[matrix_index], NW);
    quadImages[1] = Mat(vec_mat[matrix_index], NE);
    quadImages[2] = Mat(vec_mat[matrix_index], SE);
    quadImages[3] = Mat(vec_mat[matrix_index++], SW);

    int x = 0, y = 0;
    rmse = 0;
    PixelDepthQuad *temp;
    for (int i = 0; i < 4; ++i) {
      vec_mat.push_back(quadImages[i]);
      avgPixQ[i] = calculate_average_pixel(&quadImages[i]);
      // add the four quadrants to the pixvector table
      temp = new PixelDepthQuad{Pixel{avgPixQ[i][0], avgPixQ[i][1],
                                      avgPixQ[i][2], x * HWIDTH, y * HHEIGHT},
                                depth, i, 0, pixVectors[my_current_index]};
      pixVectors.push_back(temp);
      int ind = pixVectors.size() - 1;
      rmse += calculate_average_pixel_RMSE_vec3b(&quadImages[i], avgPixQ[i]);
      cout << "[" << pixVectors[ind]->P.B << "," << pixVectors[ind]->P.G << ","
           << pixVectors[ind]->P.R << "]\tQ:" << pixVectors[ind]->quad
           << "\tI:" << pixVectors[ind]->depth << "\tRMSE: " << rmse
           << "\t Parent: " << pixVectors[i]->parentPixelDepth
           << "\tMe:" << pixVectors[i] << endl;
      if (i % 2)
        ++y;
      else
        ++x;
    }
    HWIDTH /= 2;
    HHEIGHT /= 2;

    depth = 1 + matrix_index / 4;
    cout << "depth: " << depth << " rmse: " << rmse
         << " vector size: " << vec_mat.size() << "\n\n"
         << endl;
  }
}

int calculate_max_depth_of_quadtree(Mat *img) {
  auto height = img->rows;
  int depth = 0;

  while (height) {
    ++depth;
    height /= 2;
  }
  return depth;
}

void reconstructImageFromTree(vector<PixelDepthQuad *> &pixVector,
                              Mat *reconImage) {
  for (int i = 1; i < pixVector.size(); ++i) {
    // if (pixVector[i]->isLeaf) {
    int factor = pow(2, pixVector[i]->depth);
    // #ifdef DEBUG
    cout << "Factor: " << factor << " Depth: " << pixVector[i]->depth << endl;

    // #endif
    for (int row = pixVector[i]->P.x;
         row < pixVector[i]->P.x + (reconImage->rows / pow(2, pixVector[i]->depth));
         ++row) {
      for (int col = pixVector[i]->P.y;
           col < pixVector[i]->P.y + (reconImage->cols / pow(2, pixVector[i]->depth));
           ++col) {
        reconImage->at<Vec3b>(row, col) =
            Vec3b(pixVector[i]->P.R, pixVector[i]->P.G, pixVector[i]->P.B);
      }
    }
    //}
  }
}

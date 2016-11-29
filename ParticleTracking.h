/*
 * ParticleTracking.h
 *
 *  Created on: Nov 26, 2016
 *      Author: lengelen
 */

#ifndef PARTICLETRACKING_H_
#define PARTICLETRACKING_H_

#include "objectTypes.h"

/// Simple functions to compare variables etc

inline bool compare(float a, float b); // Compare-function of two float variables to check for equality (lower than threshold)
inline bool operator==(const Point2f& pt1, const Point2f& pt2); // Definition of operator to compare two points (2D), necessary for removedupes
inline void removedupes(std::vector<Point2f> & vec); // Removes doubles in list of 2D points
inline bool point_comparator( cv::Point3f a, cv::Point3f b); // Comparison function of Point3f

/// Functions to detect or process feature points

vector<Point2f> create_points(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, Mat cameraMatrix, Mat distCoeffs, bool showCorners); // Do all operations on image to obtain a sorted list of corner points
vector<Point2f> sortPattern(vector<Point2f> corners, Size boardSize);
vector<Corner> findCorrespondances(vector<Point2f> candidates, vector<Corner> prevCorners);
vector<Point2f> undistortCorners(vector<Point2f> distortedCorners, Mat cameraMatrix, Mat distCoeffs);
vector<Point2f> detectPotentialCorners(Mat img, float ResponseThreshold, float minDistance, int detectionRadius);
vector<Point2f> removeDoubles(vector<Point2f> points, Mat response, float mindist);

/// Functions to compute and process corner strength

inline float computeResponse(Mat points);
inline Mat getPoints5(Mat image, size_t x, size_t y);
inline Mat getPoints10(Mat image, size_t x, size_t y);
inline Mat corner_detect5(const size_t h, const size_t w,  Mat image);
inline Mat corner_detect10(const size_t h, const size_t w,  Mat image);

/// Master functions to detect and update corners in images

vector<Corner> createCornerList(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size PatternSize, Mat cameraMatrix, Mat distCoeffs, bool showCorners);
vector<Corner> updateCornerlist(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Mat cameraMatrix, Mat distCoeffs, bool showCorners, vector<Corner> prevCorners);
vector<Corner> readFeaturesImage(string name, CameraParams cam, string OutputName, Settings s, vector<Corner> prevCorners);
vector<Corner> readFeaturesFirstImage(string name, CameraParams cam, string OutputName, Settings s);

#endif /* PARTICLETRACKING_H_ */

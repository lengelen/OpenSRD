/*
	OpenSRD 1.0.0	
	Copyright (C) 2017  Lukas Engelen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
vector<Point2f> create_points_OpenCV(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, Mat cameraMatrix, Mat distCoeffs, bool showCorners); // Do all operations on image to obtain a sorted list of corner points
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

#endif /* PARTICLETRACKING_H_ */

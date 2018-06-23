#ifndef MORPHING_H
#define MORPHING_H
#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui.hpp>


cv::Mat PointVec2HomogeneousMat(const std::vector<cv::Point2f>& pts);


// Morph points
void MorphPoints(const std::vector<cv::Point2f>& srcPts1, const std::vector<cv::Point2f>& srcPts2, std::vector<cv::Point2f>& dstPts, float s = 0.5);

void GetTriangleVertices(const cv::Subdiv2D& sub_div, const std::vector<cv::Point2f>& points, std::vector<cv::Vec3i>& triangle_vertices);


void TransTrianglerPoints(const std::vector<cv::Vec3i>& triangle_vertices, const std::vector<cv::Point2f>& points,std::vector<std::vector<cv::Point2f>>& triangler_pts);

void PaintTriangles(cv::Mat& img, const std::vector<std::vector<cv::Point2f>>& triangles);
///// for debug /////
void DrawTriangles(cv::Mat& img, const std::vector<std::vector<cv::Point2f>>& triangles);
//////////////////////

void SolveHomography(const std::vector<cv::Point2f>& src_pts1, const std::vector<cv::Point2f>& src_pts2, cv::Mat& H);


void SolveHomography(const std::vector<std::vector<cv::Point2f>>& src_pts1,	const std::vector<std::vector<cv::Point2f>>& src_pts2, std::vector<cv::Mat>& Hmats);


// Morph homography matrix
void MorphHomography(const cv::Mat& Hom, cv::Mat& MorphHom1, cv::Mat& MorphHom2, float blend_ratio);



// Morph homography matrix
void MorphHomography(const std::vector<cv::Mat>& Homs, std::vector<cv::Mat>& MorphHoms1, std::vector<cv::Mat>& MorphHoms2, float blend_ratio);


// create a map for cv::remap()
void CreateMap(const cv::Mat& TriangleMap, const std::vector<cv::Mat>& HomMatrices, cv::Mat& map_x, cv::Mat& map_y);


//! Image Morphing
/*!
\param[in] src_img1 Input image 1
\param[in] src_points1 Points on the image 1
\param[in] src_img2 Input image 2
\param[in] src_points2 Points on the image 2, which must be corresponded to src_point1
\param[out] dst_img Morphed output image
\param[out] dst_points Morphed points on the output image
\param[in] shape_ratio blending ratio (0.0 - 1.0) of shape between image 1 and 2.  If it is 0.0, output shape is same as src_img1.
\param[in] color_ratio blending ratio (0.0 - 1.0) of color between image 1 and 2.  If it is 0.0, output color is same as src_img1. If it is negative, it is set to shape_ratio.
*/
void ImageMorphing(const cv::Mat& src_img1, const std::vector<cv::Point2f>& src_points1,
	const cv::Mat& src_img2, const std::vector<cv::Point2f>& src_points2,
	cv::Mat& dst_img, std::vector<cv::Point2f>& dst_points,
	float shape_ratio = 0.5, float color_ratio = -1);

#endif
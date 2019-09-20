#ifndef TOOL_H
#define TOOL_H

#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <ctime>
#include <io.h>
#include <list>
#include <direct.h>

using namespace cv;
using namespace std;

namespace Tiling_tiles{

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

	extern vector<pair<string, Scalar>> colorbar;

	class Line_Seg{
	public:
		Line_Seg();
		Line_Seg(Point2f a, Point2f b)
		{
			start = a;
			end = b;
		}
		Point2f start;
		Point2f end;
	};

	//all kinds tools
	//draw tool
	void MyLine(Mat img, Point2f start, Point2f end, string color1);
	Mat draw_polygen(string win_name, vector<Point2f> contour_s);
	void draw_poly(Mat &drawing_, vector<Point2f> contour_s, Point2f center, int color = 0);
	void draw_allplane(Mat &drawing_, vector<Point2f> contour_, vector<int> vec_, double scale = 1, int type = 0);
	void draw_result(Mat &drawing_, vector<Point2f> contour_, vector<int> vec_, double scale = 1, int type = 0, Point2f shift = Point2f(0, 0));

	void draw_two(Mat &drawing_, vector<Point2f> &contour_1, vector<int> vec_1, vector<Point2f> &contour_2, vector<int> vec_2, double scale = 1, int type = 0);//draw connection area
	//math tool
	Point2f center_p(vector<Point2f> contour_);
	double contour_length(vector<Point2f> contour);
	double length_two_point2f(Point2f &u, Point2f &v);
	double length_two_point_tar(vector<double> &p1, vector<double> &p2);
	void move_con(vector<Point2f> &con, Point2f sh);


	int cur_char_length(char a, char b);
	double cur_length_two_p(double cur1, double cur2);
	vector<Point2f> sampling(vector<Point2f> &contour_, int points_num);
	vector<Point2f> sampling_ave(vector<Point2f> &contour_, int points_num, vector<int> &contour_sam_index);
	vector<double> curvature_com_k(vector<Point2f> &contour_sam);
	vector<double> curvature_com(const vector<Point2f> &contour_sam); //记录cos值
	vector<int> most_convex_p(vector<Point2f> contour_, vector<double> cont_c, int max_cur_num);
	vector<int> feature_points(vector<Point2f> contour_, double dmin, double dmax, double angle_cos);//amax 165; dmin 0.015Ls
	vector<int> simple_with_feature(vector<Point2f> contour_,int totalnum);


	//translate tool
	string int2string(int number);
	string double2string(double number);

	//flipping
	vector<Point2f> flip_only_coord(vector<Point2f> cont_s, int flag = 0);

	//cross points
	int line_intersection(Line_Seg line1, Line_Seg line2, Point2f &cross_p);
	vector<Point2f> line_polygon(Line_Seg line1, vector<Point2f> contour);
	vector<Point2f> poly_poly(vector<Point2f> contour, vector<Point2f> contour_);
	bool self_intersect(vector<Point2f> &contour_, int &first, int &second);

	//vector and cosin/ sin
	Point2f unit_vec(Point2f vec);
	Point2f vertical_vec(Point2f vec);
	double cos_3edges(double l1, double l2, double l3);
	double cos_two_vector(Point2f &v0, Point2f &v1);
	double sin_two_vector(Point2f &v0, Point2f &v1);
	double tar_sin_2vector(Point2f &v0, Point2f &v1);
	vector<double> recover_consin(const vector<double> &former);

	void sort_comb(vector<double> vect, vector<int> &index_num);
	//void sort_bub(vector<int> &target);
	template<typename T>
	void sort_bub(vector<T> &target)  //从小到大
	{
		int i, j;
		T temp;
		for (i = 0; i < target.size() - 1; i++)
			for (j = 0; j < target.size() - 1 - i; j++)
				if (target[j] > target[j + 1])
				{
					temp = target[j];
					target[j] = target[j + 1];
					target[j + 1] = temp;

				}
	}
	//file cout
	void fileout(string filepath, vector<Point> contour_);
	void write_obj(string filepath, vector<Point2f> contour, double height);

	//bounding box
	void bbx_center_point(vector<vector<Point2f>> all_point, vector<Point2f> &five_p);
	vector<Point2f> b_box(vector<Point2f> contour);//返回的点是从左上方逆时针
	vector<Point2f> b_box_int(vector<Point> contour);//返回的点是从左上方逆时针

	//openglwindow

	void OpenWindow(int w, int h, vector<Point2f> contour1, vector<Point2f> contour2);


	// Morph points
	void MorphPoints(const std::vector<cv::Point2f>& srcPts1, const std::vector<cv::Point2f>& srcPts2, std::vector<cv::Point2f>& dstPts, float s);
	//void merge_close_p(vector<Point2f> &contour_);



}

#endif
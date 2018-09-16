#ifndef COMPARE_H
#define COMPARE_H

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include "tilingOpt.h"

using namespace std;
using namespace cv;

struct Event {
	int id;
	double x, y; //x = len

	bool operator<(const Event& rhs) const {
		return x < rhs.x;
	}
};
//compare
bool is_isolated(const Mat& M, int r, int c);
vector<Point2d> import_image(const char* filename);
vector<double> edge_length(const vector<Point2d>& poly);
vector<double> rot_angle(const vector<Point2d>& poly);
vector<Event> get_events(int id, const vector<double>& rot_angle, const vector<double>& edge_len, int pos);
double calculate_dist(const vector<Event>& func);
void output_events(const vector<Event>& events, const char* filename);
void compare_shapes(vector<Point2d> shape1, vector<Point2d> shape2);
double to_deg(double rad);



#endif
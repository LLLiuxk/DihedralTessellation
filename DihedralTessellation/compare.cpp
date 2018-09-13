#include "compare.h"

bool is_isolated(const Mat& M, int r, int c) {
	if (M.at<uchar>(r, c) == 0)
		return false;
	if (r > 0 && M.at<uchar>(r - 1, c) != 0)
		return false;
	if (c > 0 && M.at<uchar>(r, c - 1) != 0)
		return false;
	if (r + 1 < M.rows && M.at<uchar>(r + 1, c) != 0)
		return false;
	if (c + 1 > M.cols && M.at<uchar>(r, c + 1) != 0)
		return false;
	return true;
}

vector<Point2d> import_image(const char* filename) {
	Mat src, contour_out;
	vector<vector<Point>> contours;

	src = imread(filename, 0);
	src = src < 128;
	/*
	for (int r = 0; r < src.rows; ++r) {
	for (int c = 0; c < src.cols; ++c) {
	if (is_isolated(src, r, c))
	src.at<uchar>(r, c) = 0;
	}
	}
	*/
	findContours(src, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

	int max_id = 0;
	for (int i = 1; i < (int)contours.size(); ++i) {
		if (contours[i].size() > contours[max_id].size())
			max_id = i;
	}

	vector<Point2d> res;
	for (int i = 1; i < (int)contours[max_id].size(); ++i) {
		res.push_back(contours[max_id][i]);
	}
	// printf("%zd\t%s\n", res.size(), filename);
	return res;
}

vector<double> edge_length(const vector<Point2d>& poly) {
	int sz = poly.size();
	vector<double> res;
	double totlen = 0;
	for (int i = 0; i < sz; ++i) {
		int j = (i + 1) % sz;
		Vec2d v = poly[j] - poly[i];
		double vlen = norm(v);
		totlen += vlen;
		res.push_back(vlen);
	}
	for (int i = 0; i < sz; ++i) {
		res[i] /= totlen;
	}
	return res;
}

vector<double> rot_angle(const vector<Point2d>& poly) {
	int sz = poly.size();
	vector<double> res;
	for (int i = 0; i < sz; ++i) {
		int j = (i + 1) % sz;
		int k = (i + sz - 1) % sz;
		Vec3d v1 = Point3d(poly[i]) - Point3d(poly[k]);
		Vec3d v2 = Point3d(poly[j]) - Point3d(poly[i]);
		double v1len = norm(v1);
		double v2len = norm(v2);
		double cos_val = v1.dot(v2) / (v1len * v2len);
		double sin_val = (v1.cross(v2))[2] / (v1len * v2len);
		double angle = acos(cos_val);
		if (sin_val < 0)
			angle = -angle;
		res.push_back(angle);
	}
	return res;
}

vector<Event> get_events(int id, const vector<double>& rot_angle, const vector<double>& edge_len, int pos) {
	int sz = rot_angle.size();
	vector<Event> res;
	double angle = 0;
	double len = 0;
	int i = pos;
	do {
		angle += rot_angle[i];
		res.push_back({ id, len, angle });
		len += edge_len[i];
		i = (i + 1) % sz;
	} while (i != pos);
	return res;
}

double calculate_dist(const vector<Event>& func) {
	double prev_x = 0;
	double prev_y[2] = { -1, -1 };
	double sum = 0;
	for (const Event& e : func) {
		//printf("\tid = %d, x = %f, y = %f\n", e.id, e.x, e.y);
		if (prev_y[e.id - 1] >= 0) {
			double dx = e.x - prev_x;
			double dy = prev_y[0] - prev_y[1];
			//printf("dx = %f, dy = %f\n", dx, dy);
			sum += dx * dy * dy;
		}
		//printf("id = %d, x = %f, y = %f\n", e.id, e.x, e.y);
		prev_x = e.x;
		prev_y[e.id - 1] = e.y;
		//printf("(%f, %f)\n", prev_y[0], prev_y[1]);
	}
	double dx = 1 - prev_x;
	double dy = prev_y[0] - prev_y[1];
	sum += dx * dy * dy;
	double res = sqrt(sum);
	//printf("res = %f\n", res);
	return res;
}

void output_events(const vector<Event>& events, const char* filename) {
	FILE *fp = fopen(filename, "w");
	for (const Event& e : events) {
		fprintf(fp, "%f %f\n", e.x, e.y);
	}
	fclose(fp);
}

void compare_shapes(vector<Point2d> shape1, vector<Point2d> shape2) {
	vector<double> rot_angle1 = rot_angle(shape1);
	if (std::accumulate(rot_angle1.begin(), rot_angle1.end(), 0.0) < 0) {
		reverse(shape1.begin(), shape1.end());
		rot_angle1 = rot_angle(shape1);
	}
	vector<double> edge_len1 = edge_length(shape1);

	vector<double> rot_angle2 = rot_angle(shape2);
	if (std::accumulate(rot_angle2.begin(), rot_angle2.end(), 0.0) < 0) {
		reverse(shape2.begin(), shape2.end());
		rot_angle2 = rot_angle(shape2);
	}
	vector<double> edge_len2 = edge_length(shape2);

	/*
	puts("shape1");
	for (const Point2d& p : shape1) {
	printf("(%.0f, %.0f)\n", p.x, p.y);
	}
	puts("shape2");
	double sum_len = 0;
	double sum_angle = 0;
	for (int i = 0; i < (int)shape2.size(); ++i) {
	const Point2d& p = shape2[i];
	double len = edge_len2[i];
	double angle = rot_angle2[i];
	sum_len += len;
	sum_angle += angle;
	printf("(%.0f, %.0f)\t%f\t%f\t%f\t%f\n", p.x, p.y, len, sum_len, angle, sum_angle);
	}
	exit(0);
	*/

	double dist = HUGE_VAL;
	for (int i = 0; i < (int)shape1.size(); ++i) {
		double cur_dist = HUGE_VAL;
		for (int j = 0; j < (int)shape2.size(); ++j) {
			vector<Event> events1 = get_events(1, rot_angle1, edge_len1, i);
			vector<Event> events2 = get_events(2, rot_angle2, edge_len2, j);
			vector<Event> func(events1.size() + events2.size());
			merge(events1.begin(), events1.end(), events2.begin(), events2.end(), func.begin());
			/*
			if (i == 0 && j == 0) {
			for (const Event& e : func) {
			printf("(%d, %f, %f)\n", e.id, e.x, e.y);
			}
			// exit(0);
			}
			*/
			double ret = calculate_dist(func);
			cur_dist = min(cur_dist, ret);
			// printf("(%d, %d)\t%.6f\n", i, j, ret);
		}
		printf("%d:\t%.6f\n", i, cur_dist);
		dist = min(dist, cur_dist);
	}

	printf("dist = %.6f\n", dist);
}

double to_deg(double rad) {
	return rad / PI * 180;
}

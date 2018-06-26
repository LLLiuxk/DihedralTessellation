#ifndef TILINGOPT_H
#define TILINGOPT_H


#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "Morphing.h"


using namespace cv;
using namespace std;

namespace Tiling_tiles{

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

	

	class Prototile{
	public:
		Prototile();

		void loadTileData(string tile_data);
		void contour_sam_cur();
		void cur_normalize();

		//math tool
		double cos_two_vector(Point &v0, Point &v1);
		void sort_cos(vector<double> &vect, vector<int> &index_num);
		vector<double> curvature_com_k(vector<Point2f> &contour_sam);
		//vector<double> curvature_com(vector<Point2f> &contour_sam);

		//io polygon
		void imgtocout();
		void readTxt(vector<Point2f> &con_point);

		

		string contourname;
		vector<Point2f> contour;
		vector<vector<Point2f>> contour_sample;
		vector<vector<Point2f>> contour_sample_inver;
		vector<vector<double>> contour_curva;
		vector<vector<double>> contour_curva_inver;
		char cur_string[6][600];
		//vector<vector<Point2f>> contour_saliency;  //œ‘÷¯–‘ø…—°
		//vector<vector<Point2f>> contour_saliency_inver;

		double c_length;
		Point2f center_point;

	};


	class Tiling_opt{
	public:
		Tiling_opt();
		~Tiling_opt();
		void com_score(string imagename1, string imagename2);
		double com_scale_factor();
		double com_each_pair(vector<Point2f> &first_interval, vector<Point2f> &second_interval, int &flag);


		void divide_intervals(int intervals, vector<Point2f> &contour, vector<vector<Point2f>> &proto_interval, int delay);


		double com_optimal_score(vector<vector<Point2f>> &proto_interval_1, vector<vector<char>> &proto_first_char, 
			vector<vector<Point2f>> &proto_interval_2, vector<vector<char>> &proto_second_char, vector<pair<int, int>> &order_type);
		
		double com_tra_sim(vector<Point2f> &first_interval, vector<char>&first_char, vector<Point2f> &second_interval, vector<char>&second_char, int &flag);

		//simulation
		double com_score_manual(string imagename1, string imagename2);
	

		//math
		double warpAff_tra(vector<Point2f> &input_, vector<Point2f> &output_);
		double warpAff_tra_sec(vector<Point2f> &input_, vector<Point2f> &output_, vector<char>&second_char, vector<char>&second_char_out);
		double warpAff_tra_ref_y(vector<Point2f> &input_, vector<Point2f> &output_);

		double re_warp_Aff(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double re_warp_Aff_sec(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double re_warp_Aff_ref(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);

		double warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_,Point2f start, Point2f end);

		double Aff_place(vector<Point2f> &input1, vector<Point2f> &input2, vector<vector<Point2f>> &prototwo);
			 
		double DTW(vector<Point2f> &first_arr, vector<Point2f> &second_arr);
		void printPath(double d[][50], double dp[][50], int i, int j, vector<Point2f> &first_arr, vector<Point2f> &second_arr, vector<pair<int, int>>& path);
		double quadr_mismatch(vector<Point2f> &first_arr, vector<Point2f> &second_arr, vector<char> &first_char, vector<char> &second_char);

	private:
		vector<pair<int, int>> dp_path;
		int dp[600][600];
		int dp_inver[600][600];
		Prototile *prototile_first;
		Prototile *prototile_second;
	
	};


	//draw tool
	void MyLine(Mat img, Point2f start, Point2f end, string color1);

	//math tool
	double contour_length(vector<Point2f> contour);
	double length_two_point2f(Point2f &u, Point2f &v);
	int cur_char_length(char a, char b);
	double cur_length_two_p(double cur1, double cur2, double zeta);
	Point2f unit_vec(Point2f vec);

}

#endif
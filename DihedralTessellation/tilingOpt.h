#ifndef TILINGOPT_H
#define TILINGOPT_H


#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include "Morphing.h"


using namespace cv;
using namespace std;

namespace Tiling_tiles{

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
	
	typedef struct candPat_angle_index_error
	{
		int number;
		bool isFilp;
		int angle;
		int index;
		double mismatch;
	}CandPat;

	class Prototile{
	public:
		Prototile();

		void loadTileData(string tile_data);
		void contour_sam_cur();
		vector<int> convex_p(int max_cur_num);                 //求轮廓上值最大的10个不临近的凸点
		vector<int> partition_points(string imaname);  //求得用做划分的点
		
		void cur_normalize();

		//flipping
		vector<Point2f> flip_contour(vector<Point2f> cont_s, int flag);

		//math tool
		vector<double> curvature_com_k(vector<Point2f> &contour_sam);
		vector<double> curvature_com(vector<Point2f> &contour_sam); //记录cos值
		
		//io polygon
		void imgtocout(string tile_image);
		vector<Point2f> readTxt();
		void loadPoints(vector<Point2f> con_point);

		

		string contourname;
		vector<Point2f> contour;
		vector<double> cconvex;
		vector<vector<Point2f>> contour_sample;
		vector<vector<Point2f>> contour_sample_flip;
		vector<vector<double>> contour_curva;
		vector<vector<double>> contour_curva_flip;
		char cur_string[6][600];
		//vector<vector<Point2f>> contour_saliency;  //显著性可选
		//vector<vector<Point2f>> contour_saliency_inver;

		double c_length;
		Point2f center_point;

	};


	class Tiling_opt{
	public:
		Tiling_opt();
		~Tiling_opt();
		void com_score(string imagename1, string imagename2);

		void points_dividing(string imaname);
		bool one_situ_div(vector<int> results, vector<Point2f> contour_s, vector<Point2f> &return_B, vector<int> &return_p);

		//collision
		bool coll_detection(Point2f shifting1, Point2f shifting2, vector<Point2f> &contour_s);
		bool collision_pixel(vector<Point2f> dis_p, vector<Point2f> contour_s);
		
		//load dataset
		void load_dataset();

		//shapes comparing and candidate contour choosing
		vector<CandPat> compare_shapes(vector<Point2f> inner_c, vector<int> &mid_index);
		CandPat min_mismatch(vector<Point2f> inner, vector<Point2f> cand,  vector<double> inner_c, vector<double> cand_c, int theone, bool isFilp);
		double quadr_mismatch(vector<Point2f> first_arr, vector<Point2f> second_arr, vector<double> first_c, vector<double> second_c);
		vector<int> search_align_p(Point2f cent, Point2f end, vector<Point2f> cand_temp);
		
		//from CandPat to contour
		vector<Point2f> CandP2Contour(CandPat candp);
		


		double com_each_pair(vector<Point2f> &first_interval, vector<Point2f> &second_interval, int &flag);
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

		double Aff_place(vector<Point2f> &input1, vector<Point2f> &input2, vector<vector<Point2f>> &prototwo, vector<Point2f> &protoAff, int flag);
			 
		double DTW(vector<Point2f> &first_arr, vector<Point2f> &second_arr);
		void printPath(double d[][50], double dp[][50], int i, int j, vector<Point2f> &first_arr, vector<Point2f> &second_arr, vector<pair<int, int>>& path);

	//private:
		vector<pair<int, int>> dp_path;
		int dp[600][600];
		int dp_inver[600][600];
		Prototile *prototile_first;
		Prototile *prototile_mid;
		Prototile *prototile_second;
		vector<vector<Point2f>> contour_dataset;
	
	};

	//all kinds tools
	//draw tool
	void MyLine(Mat img, Point2f start, Point2f end, string color1);
	Mat draw_polygen(string win_name, vector<Point2f> contour_s);
	//math tool
	Point2f center_p(vector<Point2f> contour_);
	double contour_length(vector<Point2f> contour);
	double length_two_point2f(Point2f &u, Point2f &v);
	int cur_char_length(char a, char b);
	double cur_length_two_p(double cur1, double cur2, double zeta);
	void sort_comb(vector<double> vect, vector<int> &index_num);
	
	//cross points
	int line_intersection(Point2f start1, Point2f end1, Point2f start2, Point2f end2, Point2f &cross_p);
	
	//vector and cosin
	Point2f unit_vec(Point2f vec);
	double cos_two_vector(Point2f &v0, Point2f &v1);

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

	//bounding box
	void bbx_center_point(vector<vector<Point2f>> all_point, vector<Point2f> &five_p);
	vector<Point2f> b_box(vector<Point2f> contour);//返回的点是从左上方逆时针

	//morphing
	vector<Point2f> morphing_2_patterns(vector<Point2f> contour1, vector<Point2f> contour2);

	//skeleton
	Mat thinImage(const Mat & src, const int maxIterations);
	void filterOver(Mat thinSrc);
	vector<cv::Point2f> getPoints(const Mat &thinSrc, unsigned int raudis, unsigned int thresholdMax, unsigned int thresholdMin);
	vector<Point2f> get_Skeleton(string imaname, vector<Point2f> &skeleton);
	


}

#endif
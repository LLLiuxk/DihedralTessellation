#ifndef TILINGOPT_H
#define TILINGOPT_H


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
#include "gl/glut.h"
//#include "Morphing.h"


using namespace cv;
using namespace std;

namespace Tiling_tiles{

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

	extern vector<pair<string, Scalar>> colorbar;

	typedef struct candPat_angle_index_error
	{
		int number;
		bool isFilp;
		int angle;
		int index;
		double mismatch;
	}CandPat;

	typedef struct innerPat_pare
	{
		int type; //0:trans,1:rota,2:flip(13),3:flip(24)
		vector<int> in_interval;
		vector<Point2f> in_contour;
	}inPat;

	typedef struct jointPat_four
	{
		int type; //0:trans,1:rota,2:flip(13),3:flip(24)
		vector<int> interval;
		vector<vector<Point2f>> four_contour;
	}jointPat;

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

	class Prototile{
	public:
		Prototile();
		Prototile(string rname,string tpath);
		void Pro_clear();
		
		void getpath();
		void loadTileData(string tile_data);
		void contour_sam_cur();	
		vector<vector<double>> compute_TAR(vector<Point2f> &contour_,double &shape_complexity, double frac = 0.25);

		vector<int> cand_tiling_v(int max_cur_num);       //求轮廓上值最大的10个不临近的凸点
		vector<int> partition_points(string imaname);  //求得用做划分的点

		//void cur_normalize();

		//flipping
		vector<Point2f> flip_contour(vector<Point2f> cont_s, int flag = 0);

		//io polygon
		void imgtocout(string tile_image, int raw=0);
		vector<Point2f> readTxt();
		void loadPoints(vector<Point2f> con_point);



		string contourname;
		string dataroot;
		string txtpath;
		vector<Point2f> contour;
		//vector<double> cconvex;
		vector<vector<Point2f>> contour_sample;
		vector<vector<Point2f>> contour_sample_flip;
		//vector<vector<double>> contour_curva;
		//vector<vector<double>> contour_curva_flip;
		//char cur_string[6][600];
		//vector<vector<Point2f>> contour_saliency;  //显著性可选
		//vector<vector<Point2f>> contour_saliency_inver;

		double c_length;
		Point2f center_point;
	};


	class Tiling_opt{
	public:
		Tiling_opt();
		~Tiling_opt();
		void Tiling_clear();
		
		void points_dividing(string imaname);
		void tiliing_generation(string imaname);

		//collision
		
		bool coll_detec_bbx(vector<Point2f> contour1, vector<Point2f> contour2, int threshold);
		bool vertex_angle(vector<Point2f> angle1, vector<Point2f> angle2);

		//load dataset
		void load_dataset();
		void com_all_TARs(int num_c);
		void check_Repetitive_pattern();
		
		//three placement rules
		vector<vector<int>> find_rota_tilingV(vector<Point2f> cont, vector<int> mark_13);
		bool translation_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
		bool rotation_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
		bool flipping_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname,int type);
		vector<Point2f> extract_contour(vector<Point2f> contour_, vector<int> mark_p, vector<int> &midmark_p, vector<vector<Point2f>> &four_place, int type);

		//shapes comparing and candidate contour choosing
		vector<pair<int, bool>> compare_choose_TAR(vector<Point2f> inner_c); //得到选择出的pattern的序号和是否翻转的标志
		vector<pair<int, bool>> quick_choose_TAR(vector<Point2f> inner_c); //得到选择出的pattern的序号和是否翻转的标志

		double tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width = 4);//点对应匹配的筛选框宽度
		void print_TAR_Path(double d[][202], double dp[][202], int i, int j, vector<pair<int, int>>& path);
		vector<int> search_align_p(Point2f cent, Point2f end, vector<Point2f> cand_temp);

		//relocation
		vector<int> joint_relocate(vector<Point2f> contour_, vector<int> joint_index, int num_c);
		
		//morphing
		//提高cos值权重+重采样
		vector<Point2f> morphing_2_patterns(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<double> &concur1, vector<double> &concur2, vector<int> &mid_inter, float shape_ratio);
		//morphing by tar
		vector<Point2f> morphing_tar(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<int> &mid_inter, vector<pair<int, int>> &path, int shift);
		void contour_fine_tuning(vector<Point2f> &contour_, int first,int second );

		//迭代morph	
		vector<Point2f> morphing_patterns_iter(vector<Point2f> contour1, vector<Point2f> contour2, vector<double> concur1, vector<double> concur2, float shape_ratio);//, vector<int> mid_inter, float shape_ratio);

		double evalua_deformation(vector<Point2f> contour1, vector<Point2f> contour2);
		
		//simulation 
		vector<Point2f> simulation_mid(string imaname, int inner_one, int cand_one);
		jointPat simulation_tar(string imaname, int inner_one, int cand_one);

		//compute joint
		void pattern_joint(jointPat pattern);

		//math


		//old
		void com_score(string imagename1, string imagename2);
		vector<CandPat> compare_shapes(vector<Point2f> inner_c, int num_c);
		CandPat min_mismatch(vector<Point2f> inner, vector<Point2f> cand, vector<double> inner_c, vector<double> cand_c, int theone, bool isFilp);
		double quadr_mismatch(vector<Point2f> first_arr, vector<Point2f> second_arr, vector<double> first_c, vector<double> second_c, vector<pair<int, int>>& path, double zeta = 1.0);//zeta 是曲率值权重与距离值权重的倍数
		bool one_situ_div(vector<int> results, vector<Point2f> contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
		//from CandPat to contour
		vector<Point2f> CandP2Contour(CandPat candp, int num);
		bool coll_detection(Point2f shifting1, Point2f shifting2, vector<Point2f> &contour_s);
		bool collision_pixel(vector<Point2f> dis_p, vector<Point2f> contour_s);
		double warpAff_tra(vector<Point2f> &input_, vector<Point2f> &output_);
		double warpAff_tra_sec(vector<Point2f> &input_, vector<Point2f> &output_, vector<char>&second_char, vector<char>&second_char_out);
		double warpAff_tra_ref_y(vector<Point2f> &input_, vector<Point2f> &output_);
		double re_warp_Aff(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double re_warp_Aff_sec(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double re_warp_Aff_ref(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end);
		double Aff_place(vector<Point2f> &input1, vector<Point2f> &input2, vector<vector<Point2f>> &prototwo, vector<Point2f> &protoAff, int flag);
		double DTW(vector<Point2f> &first_arr, vector<Point2f> &second_arr);
		void printPath(double d[][202], double d_c[][202], double dp[][202], int i, int j, vector<pair<int, int>>& path);
		double com_each_pair(vector<Point2f> &first_interval, vector<Point2f> &second_interval, int &flag);
		double com_optimal_score(vector<vector<Point2f>> &proto_interval_1, vector<vector<char>> &proto_first_char,
			vector<vector<Point2f>> &proto_interval_2, vector<vector<char>> &proto_second_char, vector<pair<int, int>> &order_type);
		double com_tra_sim(vector<Point2f> &first_interval, vector<char>&first_char, vector<Point2f> &second_interval, vector<char>&second_char, int &flag);
		double com_score_manual(string imagename1, string imagename2);
		//private:
		//vector<pair<int, int>> dp_path;
		//int dp[600][600];
		//int dp_inver[600][600];
		double dis[202][202];//两组点之间的坐标差异
		double dis_cur[202][202];//两组点之间的曲率差异
		double distance[202][202];
		int step[202][202];//记录总的步数
		int all_types;
		Prototile *prototile_first;
		Prototile *prototile_mid;
		Prototile *prototile_second;
		Prototile *prototile_tem;
		vector<vector<Point2f>> contour_dataset;
		vector<vector<vector<double>>> all_con_tars;
		vector<vector<vector<double>>> all_con_tars_flip;
		vector<vector<vector<double>>> all_fea_tars;
		vector<vector<vector<double>>> all_fea_tars_flip;
		vector<double> all_shape_complexity;
	};

	//all kinds tools
	//draw tool
	void MyLine(Mat img, Point2f start, Point2f end, string color1);
	Mat draw_polygen(string win_name, vector<Point2f> contour_s);
	void draw_poly(Mat &drawing_, vector<Point2f> contour_s, Point2f center,int color=0);
	void draw_allplane(Mat &drawing_, vector<Point2f> contour_, vector<int> vec_, double scale = 1,int type = 0);
	void draw_result(Mat &drawing_, vector<Point2f> contour_, vector<int> vec_, double scale = 1, int type = 0, Point2f shift = Point2f(0, 0));
	
	void draw_two(Mat &drawing_, vector<Point2f> &contour_1, vector<int> vec_1, vector<Point2f> &contour_2, vector<int> vec_2, double scale = 1, int type = 0);//draw connection area
	//math tool
	Point2f center_p(vector<Point2f> contour_);
	double contour_length(vector<Point2f> contour);
	double length_two_point2f(Point2f &u, Point2f &v);
	double length_two_point_tar(vector<double> &p1,vector<double> &p2);
	void move_con(vector<Point2f> &con, Point2f sh);


	int cur_char_length(char a, char b);
	double cur_length_two_p(double cur1, double cur2);
	vector<Point2f> sampling(vector<Point2f> &contour_, int points_num);
	vector<double> curvature_com_k(vector<Point2f> &contour_sam);
	vector<double> curvature_com(const vector<Point2f> &contour_sam); //记录cos值
	vector<int> most_convex_p(vector<Point2f> contour_, vector<double> cont_c, int max_cur_num);
	vector<int> feature_points(vector<Point2f> contour_, double dmin, double dmax, double angle_cos);//amax 165; dmin 0.015Ls
	//translate tool
	string int2string(int number);
	string double2string(double number);

	//flipping
	vector<Point2f> flip_only_coord(vector<Point2f> cont_s, int flag = 0);

	//cross points
	int line_intersection(Line_Seg line1, Line_Seg line2, Point2f &cross_p);
	vector<Point2f> line_polygon(Line_Seg line1, vector<Point2f> contour);
 	bool self_intersect(vector<Point2f> &contour_, int &first, int &second);

	//vector and cosin/ sin
	Point2f unit_vec(Point2f vec);
	Point2f vertical_vec(Point2f vec);
	double cos_3edges(double l1,double l2,double l3);
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
	//bounding box
	void bbx_center_point(vector<vector<Point2f>> all_point, vector<Point2f> &five_p);
	vector<Point2f> b_box(vector<Point2f> contour);//返回的点是从左上方逆时针
	
	//openglwindow

	void OpenWindow(int w, int h, vector<Point2f> contour1, vector<Point2f> contour2);

	
	// Morph points
	void MorphPoints(const std::vector<cv::Point2f>& srcPts1, const std::vector<cv::Point2f>& srcPts2, std::vector<cv::Point2f>& dstPts, float s);
	//void merge_close_p(vector<Point2f> &contour_);
	



}

#endif
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
#include "Tool.h" 
//#include "Morphing.h"


using namespace cv;
using namespace std;

namespace Tiling_tiles{

#define Feature_Min 1
#define Feature_Max 3

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
		vector<int> feature_index;
		vector<int> in_interval;
		vector<Point2f> in_contour;
	}inPat;

	typedef struct jointPat_four
	{
		int type; //0:trans,1:rota,2:flip(13),3:flip(24)
		vector<int> interval;
		vector<vector<Point2f>> four_contour;
	}jointPat;

	typedef struct Point_feature
	{
		Point2f point;
		int type; //0:普通点 1:候选点 2:特征点
	}Point_f;

	class Prototile{
	public:
		Prototile();
		Prototile(string rname,string tpath);
		void Pro_clear();
		
		void setpath();
		void setname(string con_name);
		void loadTileData(string tile_data);
		void contour_sam_cur(int show_mat = 0);	
		vector<vector<double>> compute_TAR(vector<Point2f> &contour_,double &shape_complexity, double frac = 0.25); //frac*num是计算的点数目

		//vector<int> cand_tiling_v(int max_cur_num);       //求轮廓上值最大的10个不临近的凸点
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
		vector<Point2f> contour_r;
		vector<Point_f> contour_f;
		//vector<double> cconvex;
		vector<vector<Point2f>> contour_sample;
		vector<vector<Point2f>> contour_sample_flip;
		vector<int> feature_p;
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
		

		//three tiling rules
		int Rotation_rule(vector<int> part_points_index, vector<Point2f> &contour_s, string rootname);
		int Tanslation_rule(vector<int> part_points_index, vector<Point_f> &contour_s, string rootname);
		int Flipping_rule(vector<int> part_points_index, vector<Point2f> &contour_s, string rootname);

		vector<vector<Point2f>> find_rota_tilingV(vector<Point2f> &cont, vector<int> mark_13, vector<pair<Point2f, int>> &all_insert_points);
		bool translation_placement(vector<int> &feature_p, vector<int> results, vector<Point_f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
		bool rotation_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
		bool flipping_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname,int type);
		vector<Point2f> extract_contour_f(vector<Point_f> contour_, vector<int> mark_p, vector<int> &feature_p, vector<int> &midmark_p, vector<vector<Point_f>> &four_place, int type);
		vector<Point2f> extract_contour(vector<Point2f> contour_, vector<int> mark_p, vector<int> &midmark_p, vector<vector<Point2f>> &four_place, int type);
		//match candidate patterns
		void match_candidate(int Tiling_index);

		//shapes comparing and candidate contour choosing
		vector<pair<int, bool>> compare_choose_TAR(vector<Point2f> inner_c); //得到选择出的pattern的序号和是否翻转的标志
		vector<pair<int, bool>> quick_choose_TAR(vector<Point2f> inner_c); //得到选择出的pattern的序号和是否翻转的标志

		double tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width = 4);//点对应匹配的筛选框宽度
		void print_TAR_Path(double d[][202], double dp[][202], int i, int j, vector<pair<int, int>>& path);
		vector<int> search_align_p(Point2f cent, Point2f end, vector<Point2f> cand_temp);
		//relocation
		vector<int> joint_relocate(vector<Point2f> contour_, vector<int> joint_index, int num_c);

		//morphing
		vector<int> morphed_results(vector<Point2f> &morphed_A, int Candidate_index, int Tiling_index);
		void contour_fine_tuning(vector<Point2f> &contour_, int first, int second);
		
		//morphing
		//提高cos值权重+重采样
		vector<Point2f> morphing_2_patterns(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<double> &concur1, vector<double> &concur2, vector<int> &mid_inter, float shape_ratio);
		//morphing by tar
		vector<Point2f> morphing_tar(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<int> &mid_inter, vector<pair<int, int>> &path, int shift);

		//迭代morph	
		vector<Point2f> morphing_patterns_iter(vector<Point2f> contour1, vector<Point2f> contour2, vector<double> concur1, vector<double> concur2, float shape_ratio);//, vector<int> mid_inter, float shape_ratio);

		double evalua_deformation(vector<Point2f> contour1, vector<Point2f> contour2);
		
		//simulation 
		vector<Point2f> simulation_mid(string imaname, int inner_one, int cand_one);
		jointPat simulation_tar(string imaname, int inner_one, int cand_one);

		//compute joint
		void pattern_joint(jointPat pattern);
		vector<Point2f> construct_joint(jointPat pattern, int &mid);
		//math
		vector<Point2f> p2f2p_f(vector<Point_f> origin);

		//from CandPat to contour
		vector<Point2f> CandP2Contour(CandPat candp, int num);
		//bool coll_detection(Point2f shifting1, Point2f shifting2, vector<Point2f> &contour_s);
		//bool collision_pixel(vector<Point2f> dis_p, vector<Point2f> contour_s);

		//old
		//void com_score(string imagename1, string imagename2);
		//vector<CandPat> compare_shapes(vector<Point2f> inner_c, int num_c);
		//CandPat min_mismatch(vector<Point2f> inner, vector<Point2f> cand, vector<double> inner_c, vector<double> cand_c, int theone, bool isFilp);
		//double quadr_mismatch(vector<Point2f> first_arr, vector<Point2f> second_arr, vector<double> first_c, vector<double> second_c, vector<pair<int, int>>& path, double zeta = 1.0);//zeta 是曲率值权重与距离值权重的倍数
		//bool one_situ_div(vector<int> results, vector<Point2f> contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname);
				
		/*double warpAff_tra(vector<Point2f> &input_, vector<Point2f> &output_);
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
		double com_score_manual(string imagename1, string imagename2);*/

		//private:
		//vector<pair<int, int>> dp_path;
		//int dp[600][600];
		//int dp_inver[600][600];
		double dis[202][202];//两组点之间的坐标差异
		double dis_cur[202][202];//两组点之间的曲率差异
		double distance[202][202];
		int step[202][202];//记录总的步数
		int all_types;
		int sampling_num;
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
		vector<inPat> all_inner_conts;

		vector<vector<Point2f>> candidate_contours;
		vector<vector<pair<int, int>>> cand_paths;
		vector<int> mid_inter; //每个tiling placement对应一个mid_inter
		vector<int> mid_inter_morphed; //变形后的mid_inter


	};

}

#endif
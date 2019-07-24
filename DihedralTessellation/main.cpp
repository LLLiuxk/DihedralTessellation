#include "tilingOpt.h"
#include  <stdio.h>
#include  <stdlib.h>

//#include "stdafx.h"

using namespace Tiling_tiles;

int main(int argc, char** argv)
{
	clock_t start, midtime,finish;
	start = clock();
	
	int width = 200;

	Tiling_tiles::Tiling_opt *tiling_opt;
	tiling_opt = new Tiling_tiles::Tiling_opt();
	Tiling_tiles::Prototile *prototile_first;
	prototile_first = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_second;
	prototile_second = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_third;
	prototile_third = new Tiling_tiles::Prototile();
	//////prototile_first->imgtocout(imagename1);
	int f = 0;
	if (f == 12)
	{
		string fff = "hahahah.txt";
		cout << fff.size()<<endl;
		string t = fff.substr(fff.size() - 3, fff.size() - 1);
		cout << t;
	}
	if (f ==0) //已有dataset
	{

		tiling_opt->tiliing_generation("488");
		//vector<int> p_p_index = prototile_first->partition_points("test");
		//Mat drawing5 = Mat(3*width, 3*width, CV_8UC3, Scalar(255, 255, 255));
		//vector<Point2f> conr = prototile_first->contour_sample[1];
		//int contoursize = conr.size();
		////cout << contoursize << endl;
		//for (int j = 0; j < contoursize; j++)
		//{
		//	circle(drawing5, conr[j], 2, Scalar(0, 0, 0), -1);

		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//imshow("sample", drawing5);
		//int iii[20] = {48,75,80,124,220,218,212,228,248,251,280,312,317};//{484,483,482,482,480,479,478,474,472,467,451}; // { 0, 36, 278, 293, 298, 312, 357, 454, 452, 443, 429, 364 };//{458,451,439,439,426,357,359,59,70,75,102,128,169,175,223,312};//{205,210,215,219,230,233,180};//
		//for (int i = 5; i < 486; i++)
		//{
		//	tiling_opt->Tiling_clear();
		////	tiling_opt->tiliing_generation(int2string(iii[i]));
		////	
		//    tiling_opt->tiliing_generation(int2string(i));
		//	system("cls");
		//}

		
		// ________________feature_points___________________________

		//Mat drawing2 = Mat(600, 600, CV_8UC1, Scalar(255));
		//prototile_first->loadTileData("0");
		//vector<Point2f> a = prototile_first->contour_sample[2];
		//Point2f cent = center_p(a);	
		//draw_poly(drawing2, a, cent);
		//circle(drawing2, Point2f(300,300), 2, Scalar(255), -1);
		//circle(drawing2, cent, 2, Scalar(255), -1);
		//imshow("as", drawing2);
		//vector<Point2f> a = prototile_first->contour_sample[prototile_first->contour_sample.size()-1];
		//vector<double> b = curvature_com(a);
		//cout << a.size() << endl;
		//vector<int> c = most_convex_p(a, b, 20);
		//cout << "most num:  "<<c.size() << endl;
		//Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//for (int i = 0; i < a.size(); i++)
		//{
		//	circle(drawing4, a[i], 2, Scalar(0, 0, 0), -1);
		//}
		//for (int j = 0; j < c.size(); j++)
		//{
		//	circle(drawing4, a[c[j]], 4, Scalar(0, 0, 255), -1);
		//}
		//imshow("most", drawing4);
		//vector<int> d = feature_points(a, 1, 3, cos(PI * 160 / 180)); //dmin ,dmax 指点的数目,200个点每个点之间距离约为0.005
		//cout <<"feature:  "<< d.size() << endl;
		//Mat drawing3 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//for (int i = 0; i < a.size(); i++)
		//{
		//	circle(drawing3, a[i], 2, Scalar(0, 0, 0), -1);
		//}
		//for (int j = 0; j < d.size(); j++)
		//{
		//	circle(drawing3, a[d[j]], 4, Scalar(255, 0, 0), -1);
		//}
		//imshow("fea", drawing3);
		
		//--------------------test morphing------------

		//tiling_opt->load_dataset();
		//prototile_first->loadTileData("test20");
		////imgtocout();
		//vector<Point2f> contour_inner = prototile_first->contour_sample[1];
		//vector<double> coninner_cur = curvature_com(contour_inner); //prototile_first->contour_curva[1];
		//vector<CandPat> candida_contours;
		//candida_contours = tiling_opt->compare_shapes(prototile_first->contour, 1);
		//CandPat tem = candida_contours[0];
		////prototile_second->Pro_clear();
		////prototile_second->loadPoints(tiling_opt->contour_dataset[tem.number]);
		//vector<Point2f> contour_cand = tiling_opt->CandP2Contour(tem, 1);
		//vector<double> concand_cur = curvature_com(contour_cand);
		//vector<Point2f> inter_;
		//vector<int> mid_inter;
		//mid_inter.push_back(0);
		//mid_inter.push_back(60);
		//mid_inter.push_back(120);
		//mid_inter.push_back(175);
		//inter_ = tiling_opt->morphing_2_patterns(contour_inner, contour_cand, coninner_cur, concand_cur, mid_inter, 0.5);
		////inter_ = tiling_opt->morphing_patterns_iter(contour_inner, contour_cand, coninner_cur, concand_cur, 0.25);
		//cout <<"inter_: "<< inter_.size() << endl;
		////MorphPoints(contour_inner, contour_cand, inter_, 0.5);
		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
	 //   Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
	 //   draw_poly(ttt, inter_, Point2f(400, 400));
		//draw_poly(drawing_pro1, contour_cand, Point2f(400, 400));
		//draw_poly(tt, contour_inner, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//for (int i = 0; i < 4; i++)
		//{
		//	circle(tt, contour_inner[mid_inter[i]], 5, Scalar(255, 0, 0), -1);
		//}
		////
		//imshow("only point", ttt);
	 //   imshow("aaab", tt);
		//imshow("aaa", drawing_pro1);

		//-------------测试min_mismatch函数-----------

		//prototile_first->loadTileData("test1");
		//////Mat img1;
		////double leng = prototile_first->c_length;
		////vector<Point2f> a = prototile_first->contour;
		//vector<Point2f> b = prototile_first->contour_sample[1];
		//vector<double> b_c = prototile_first->contour_curva[1];
		//vector<Point2f> a = prototile_first->contour_sample_flip[1];
		//vector<double> a_c = prototile_first->contour_curva_flip[1];
		//
		////for (int i = 0; i < 130; i++)
		////{
		////	cout << "b: " <<b[i]<<" b_c[i]: "<< b_c[i] << "   a: " <<a[i]<<"  a_c: "<< a_c[a_c.size() - 1 - i] << endl;
		////}
		//prototile_second->loadTileData("686");
		//vector<Point2f> con2 = prototile_second->contour_sample[1];
		//vector<double> con_c = prototile_second->contour_curva[1];
		////if (b.size() == b_c.size()) cout << "b.size:" << b.size() << endl;
		////if (con2.size() == con_c.size()) cout << "b.size:" << con2.size() << endl;
	 //   CandPat hahah = tiling_opt->min_mismatch(a, con2, a_c, con_c,22, true);
		//CandPat hahah1 = tiling_opt->min_mismatch(b,con2, b_c, con_c,22,false);
		//cout <<hahah.angle<<" : "<< hahah.mismatch << endl;
		//cout << hahah1.angle << " no: " << hahah1.mismatch << endl;

        //------------test quadr_mismatch-----------------

		//prototile_first->loadTileData("685");
		//////Mat img1;
		////double leng = prototile_first->c_length;
		////vector<Point2f> a = prototile_first->contour;
		//vector<Point2f> b = prototile_first->contour_sample[1];
		//vector<double> b_c = prototile_first->contour_curva[1];
		//vector<Point2f> a = prototile_first->contour_sample_flip[1];
		//vector<double> a_c = prototile_first->contour_curva_flip[1];
		////for (int i = 0; i < 130; i++)
		////{
		////	cout << "b: " <<b[i]<<" b_c[i]: "<< b_c[i] << "   a: " <<a[i]<<"  a_c: "<< a_c[a_c.size() - 1 - i] << endl;
		////}
		//prototile_second->loadTileData("19");
		//vector<Point2f> con2 = prototile_second->contour_sample[1];
		//vector<double> con_c = prototile_second->contour_curva[1];
		////if (b.size() == b_c.size()) cout << "b.size:" << b.size() << endl;
		////if (con2.size() == con_c.size()) cout << "b.size:" << con2.size() << endl;
		//vector<pair<int, int>> dppath;
		////int arr[800][800] = {0};
		//double hahah = tiling_opt->quadr_mismatch(b, con2, b_c, con_c,dppath);
		//
		//cout << hahah <<" : "  << endl;

        //------------test getRotationMatrix2D-----------------

		//Mat rot_mat = getRotationMatrix2D(center_p(con2), hahah.angle, 1);
		//vector<Point2f> d;
		//transform(a, d, rot_mat);
		//int col = 800;
		//int row = 800;
		//Mat drawing_pro = Mat(col, row, CV_8UC3, Scalar(255, 255, 255));
		//int n = d.size();
		////cout << "n: " << n << endl;
		//Point rook_points[1][2000];
		//for (int t = 0; t < n; t++)
		//{
		//	rook_points[0][t] = d[t];
		//}
		//const Point* ppt[1] = { rook_points[0] };
		//int npt[] = { n };
		//fillPoly(drawing_pro,
		//	ppt,
		//	npt,
		//	1,
		//	Scalar(0, 0, 0) //黑色
		//	//Scalar(255, 255, 255) //白色
		//	);
		//circle(drawing_pro, d[0], 4, Scalar(255), 3);
		//circle(drawing_pro, d[hahah.index], 4, Scalar(255, 0, 255), 3);
		//Point2f cenpo = center_p(d);
		//circle(drawing_pro, cenpo, 4, Scalar(0, 255, 255), -1);
		//imshow("4: ", drawing_pro);
  //  	if (hahah.mismatch < hahah1.mismatch)
		//{
		//	cout << "flip" << endl;
		//	cout << "angle: " << hahah.angle << "  index: " << hahah.index << "  mismatch: " << hahah.mismatch << endl;
		//}
		//else
		//{
		//	cout << "right" << endl;
		//	cout << "angle: " << hahah1.angle << "  index: " << hahah1.index << "  mismatch: " << hahah1.mismatch << endl;
		//}
		//------------------------

        //------------------测试compare_shapes函数---------------

		//tiling_opt->load_dataset();
		//midtime = clock();
		//cout << endl << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		//prototile_first->contourname = "19";
		//prototile_first->contour = prototile_first->readTxt();
		////cout << prototile_first->contour[5] << endl;
		//
		////prototile_first->contour_sam_cur();
		//vector<CandPat> candida_contours;
		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(tt, prototile_first->contour, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//imshow("aaa", tt); 
		//candida_contours = tiling_opt->compare_shapes(prototile_first->contour, 1);
        ////-------------------------------------------
			
	}
	else if (f == 1)   //test draw_allplane
	{
		//int i = 2;
		//while (i-- != 0)
		//{
			//prototile_first->imgtocout("777");
			//prototile_first->loadTileData("600");
			Mat show = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
			//Mat show1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
			//Mat show2 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		    vector<Point2f> a;
			a.push_back(Point2f(0, 0));
			a.push_back(Point2f(50, 20));
			a.push_back(Point2f(55, 40));
			a.push_back(Point2f(50, 60));
			a.push_back(Point2f(0, 40));
			a.push_back(Point2f(5, 20));
		    /*a.push_back(Point2f(0, 40));
		    a.push_back(Point2f(30, 0));
		    a.push_back(Point2f(80, 40));
		    a.push_back(Point2f(30, 70));*/
			//vector<vector<Point2f>> b;
			vector<int> tt;
			//vector<int> ttt;
			//vector<int> tttt;
			////vector<Point2f> ccc = prototile_first->contour_sample[3];
			tt.push_back(0);
			tt.push_back(1);
			tt.push_back(3);
			tt.push_back(4);
			//vector<Point2f> c = tiling_opt->extract_contour(a, tt, ttt, b, 3);
			//draw_allplane(show, a, tt, 1, 1);
			//draw_poly(show, a, center_p(a),1);
			//draw_allplane(show1, c, ttt, 1, 3);
			draw_result(show, a, tt, 1, 0);
			imshow("hahah",show);
			//imshow("aahahah", show1);
			//vector<Point2f> d = tiling_opt->extract_contour(c, ttt, tttt, b, 3);
			//draw_allplane(show2, d, tttt, 1, 3);
			//imshow("final", show2);
		    //--------------------test line_intersection 函数--------------------------------------
			//vector<Point2f> a;
			//a.push_back(Point2f(178,60));
			//a.push_back(Point2f(175.068,61));
			//a.push_back(Point2f(33, 121.33));
			//a.push_back(Point2f(33, 121.33));
			//Point2f c;
			///*a.push_back(Point2f(119, 93.1645));
			//a.push_back(Point2f(119, 96.97));
			//a.push_back(Point2f(119, 100.775));
			//a.push_back(Point2f(119, 104.581));*/
			////Point2f c;
			//Line_Seg line1(a[2], a[0]);
			//Line_Seg line2(a[3], a[1]);
			//cout << line1.start << " " << line1.end << endl << line2.start << " " << line2.end << endl;
			//if (line_intersection(line1, line2, c)==1) cout << c << endl;
			//else cout << "line=0";
			//-------------------------------------------------------------------------------------

			//for (int i = 0; i < ccc.size(); i++)
			//{
			//	circle(show,ccc[i],1,Scalar(255),-1);
			//}
			////show = draw_polygen("aaa", ccc);
			//cout << ccc[14] << "  " << ccc[15] << endl
			//	<< ccc[16] << " " << ccc[17] << endl;
			////circle(show, ccc[14], 1, Scalar(0, 0, 255), -1);
			//circle(show, ccc[15], 2, Scalar(0, 0, 255), -1);
			//circle(show, ccc[16], 2, Scalar(0, 255, 0), -1);
			////circle(show, ccc[17], 1, Scalar(0, 255, 0), -1);

			////MyLine(show, ccc[14], ccc[15], "red");
			////MyLine(show, ccc[16], ccc[17], "green");
			////circle(show, center_p(prototile_first->contour),4,Scalar(255,0,0),-1);
			//imshow("777", show);
		//}
		//prototile_first->txtpath = "D:\\VisualStudioProjects\\images\\txt\\";
		//prototile_first->dataroot = "D:\\VisualStudioProjects\\images\\scr\\";
		//prototile_first->imgtocout("13",1);

	}	
	else if (f==2){  //批量读图
		//int iii[] = {714,763,736};// {610, 614, 618, 636, 637, 638};
		int t =  0;
		//int i = 308;
		//for (int i = 485; i < 486; i++)
		{
			//string image = int2string(iii[i]);
			string image = "488";//int2string(i);
			////cout << image << endl;
			/*string image1 = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + int2string(i) + ".png"; //"C:\\Users\\liuxk\\Desktop\\datasetnew\\zym\\datasetnew\\" + image + ".png";
			Mat src = imread(image1);
			if (src.empty())
			{
				cout << "empty!" << endl;
				continue;
			}
			else
			{
			    cout << t << endl;
				string image2 = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + int2string(t) + ".png";
				imwrite(image2, src);
				t++;
			}*/
			prototile_first->txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";
			prototile_first->dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";
			prototile_first->imgtocout(image);


			//像素点操作
			//threshold(src, src, 128, 255, cv::THRESH_BINARY);
			//imshow("???:", src);
			//int rows = src.rows;
			//int cols = src.cols;
			////////int count = 0;
			//for (int i = 0; i <rows; i++)
			//	for (int j = 0; j < cols; j++)
			//	{
			//		//cout << (int)src.at<uchar>(i, j) << " ";
			////		//if (src.at<uchar>(i, j) == 0 || src.at<uchar>(i, j) == 255) count++;
			//		if (src.at<uchar>(i, j) == 0)
			//		{
			////			count++;
			//			src.at<uchar>(i, j) = 255;
			//		}
			//		else if (src.at<uchar>(i, j) == 255) src.at<uchar>(i, j) = 0;
			//	}
			//imshow("!!!:", src);
			//////cout << count << endl;


		    //string image2 = int2string(i+534);
			//imwrite("D:\\VisualStudioProjects\\images\\new\\" + image2 + ".png",src);
			//imwrite("D:\\VisualStudioProjects\\p16.png", src);

			//prototile_first->Pro_clear();
			//prototile_first->txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";
			//prototile_first->dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";

			////prototile_first->txtpath =  "C:\\Users\\liuxk\\Desktop\\shape\\txt\\";
			////prototile_first->dataroot = "C:\\Users\\liuxk\\Desktop\\shape\\new\\";
			//prototile_first->imgtocout(image,1);

			//////
			//////string image = int2string(i);
			//////cout << image << endl;
			//image = "3 (" + image + ")";

			//if (prototile_first->contour.size() < 300)
			//cout << image + ".png may be error" << endl;
		}
	}
	else if (f == 3) // 将txt文件保存为黑色图像
	{		
		for (int i = 750; i < 751; i++)
		{
			string image = int2string(i);
			prototile_first->Pro_clear();
			prototile_first->txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";
			prototile_first->dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";
			//prototile_first->txtpath = "C:\\Users\\liuxk\\Desktop\\shape\\txt\\";
			//prototile_first->dataroot = "C:\\Users\\liuxk\\Desktop\\shape\\new\\";
			prototile_first->contourname = image;
			prototile_first->contour = prototile_first->readTxt();
			Mat drawing_ = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
			draw_poly(drawing_, prototile_first->contour, Point2f(300, 300));

			string image2 = int2string(i + 605);
			imwrite("D:\\VisualStudioProjects\\DihedralTessellation\\datasetnew\\" + image + ".png", drawing_);
		}
	}
	else if (f == 4)  //simulation of one tiling
	{
		/*Point2f s_axis(10,10);
		Point2f line1(0, -5);
		double cos_ax_1_3 = cos_two_vector(s_axis, line1);
		Point2f line3 = s_axis;
		if (cos_ax_1_3 < 0) line3 = -s_axis;
		line3 = line3*(abs(line1.x*line3.x + line1.y*line3.y) / (line3.x*line3.x + line3.y*line3.y));
		cout << line3 << endl;*/
		//for (int i = 0; i < 60; i++)
		//{
			tiling_opt->Tiling_clear();
		//	//vector<Point2f> sim_mid = tiling_opt->simulation_mid("19", 20, 0);
			vector<Point2f> sim_mid = tiling_opt->simulation_tar("489", 176, 3);
		//}
		//-----------------test compute_TAR----------------------

		/*vector<Point2f> a;
		a.push_back(Point2f(100,10));
		a.push_back(Point2f(150, 10));
		a.push_back(Point2f(150, 60));
		a.push_back(Point2f(200, 60));
		a.push_back(Point2f(210, 50));
		a.push_back(Point2f(200, 80));
		a.push_back(Point2f(100, 120));
		a.push_back(Point2f(100, 110));
		a.push_back(Point2f(50, 110));
		a.push_back(Point2f(50,60));
		a.push_back(Point2f(100, 60));
		Mat drawing_ = Mat(300, 300, CV_8UC3, Scalar(255,255,255));
		draw_poly(drawing_,a,Point2f(150,150));*/
		//double shape_com;
		//vector<vector<double>> tar_all = prototile_first->compute_TAR(a, shape_com);
		//double b = sin_2vector_convexc(Point2f(5,0), Point2f(0,5));
		//cout << b << endl;
		//for (int i = 0; i < tar_all.size(); i++)
		//{
		//	cout << tar_all[i][0] << "  ";
		//}
		//imshow("aaaa", drawing_);	
		//prototile_second->loadTileData("157");
		//double shape_com;
		//double shape_com_flip;
		//vector<vector<double>> tar_all = prototile_second->compute_TAR(prototile_second->contour_sample[1], shape_com);//(num_c+1)*100 points
		//vector<vector<double>> tar_all_flip = prototile_second->compute_TAR(prototile_second->contour_sample_flip[1], shape_com_flip);
		//cout << "shape_com: " << shape_com << "  shape_com_flip: " << shape_com_flip << endl;

		//-------------------test tar_mismatch函数--------------------

		//tiling_opt->load_dataset();
		//tiling_opt->com_all_TARs(1);
		//prototile_first->loadTileData("test12");
		//vector<Point2f> a = prototile_first->contour;// prototile_first->contour_sample[1];
		//Mat drawing_ = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//int asize = a.size();
		//cout << asize << endl;
		//for (int i = 0; i < asize; i++)
		//{
		//	circle(drawing_, a[i], 2, Scalar(0, 255, 0), -1);
		//}
		//imshow("asa", drawing_);
		////a.pop_back();
		//double sca;
		//vector<vector<double>> tar_all = prototile_first->compute_TAR(a, sca);
		//vector<int> cand_points_index = most_convex_p(a, curvature_com(a), 30);
		//cout << "feature :" << cand_points_index.size() << endl;
		//vector<vector<double>> tar_fea;
		//for (int j = 0; j < cand_points_index.size(); j++)
		//{
		//	tar_fea.push_back(tar_all[cand_points_index[j]]);
		//}
		///*prototile_second->loadTileData("23");
		//vector<Point2f> b = prototile_second->contour_sample[1];
		//vector<Point2f> c = prototile_second->contour_sample_flip[1];
		//double scb;
		//vector<vector<double>> tar_all1 = prototile_second->compute_TAR(b, scb);
		//vector<vector<double>> tar_all2 = prototile_second->compute_TAR(c, scb);*/
		//int num = 122;
		//prototile_second->loadPoints(tiling_opt->contour_dataset[num]);
		//vector<Point2f> b = prototile_second->contour_sample[1];
		//vector<Point2f> c = prototile_second->contour_sample_flip[1];
		//double scb = tiling_opt->all_shape_complexity[num];
		//cout << "sca: " << sca << "   scb:  " << scb << endl;
		//vector<vector<double>> tar_all1 = tiling_opt->all_con_tars[num];//prototile_second->compute_TAR(b, scb);
		//vector<vector<double>> tar_all2 = tiling_opt->all_con_tars_flip[num];//prototile_second->compute_TAR(c, scb);
		//vector<vector<double>> tar_fea2;
		//vector<vector<double>> tar_fea_flip2;
		//cand_points_index = most_convex_p(b, curvature_com(b), 30);// feature_points(b, 1, 3, cos(PI * 160 / 180));
		//cout << "feature :" << cand_points_index.size() << endl;
		//for (int j = 0; j < cand_points_index.size(); j++)
		//{
		//	tar_fea2.push_back(tar_all1[cand_points_index[j]]);
		//}
		//cand_points_index = most_convex_p(c, curvature_com(c), 30); //feature_points(c, 1, 3, cos(PI * 160 / 180));
		//cout << "feature :" << cand_points_index.size() << endl;
		//for (int j = 0; j < cand_points_index.size(); j++)
		//{
		//	tar_fea_flip2.push_back(tar_all2[cand_points_index[j]]);
		//}
		//vector<pair<int, int>> path;
		////double re = length_two_point_tar(tar_all[0], tar_all1[0]);
		//int shift = 0;
		//double re = tiling_opt->tar_mismatch(tar_all, tar_all1, path, shift);//点对应匹配的筛选框宽度
		//cout << "result: " << re << "  shift:" << shift << endl;
	 //   double re2 = tiling_opt->tar_mismatch(tar_all, tar_all2, path, shift);
		//cout << "result2: " << re2 << "  shift:" << shift << endl;
		//double re3 = tiling_opt->tar_mismatch(tar_fea, tar_fea2, path, shift);
		//cout << "result3: " << re3 << "  shift:" << shift << endl;
		//double re4 = tiling_opt->tar_mismatch(tar_fea, tar_fea_flip2, path, shift);
		//cout << "result4: " << re4 << "  shift:" << shift << endl;
		//re = re / (1 + sca + scb);
		//re2 = re2 / (1 + sca + scb);
		//re3 = re3 / (1 + sca + scb);
		//re4 = re4 / (1 + sca + scb);
		//
		//cout << "result: " <<re<<"  result2: "<< re2 << "  shift:" << shift << endl;
		//cout << "result3: " << re3 << "  result4: " << re4 << "  shift:" << shift << endl;
		//cout << path.size() << endl;
		/*for (int i = 0; i < path.size(); i++)
		{
		cout << path[i].first << "--"<<path[i].second << endl;
		}*/

		//---------------  test compare_choose_TAR-------------------

		/*tiling_opt->load_dataset();
		tiling_opt->com_all_TARs(1);
		midtime = clock();
		
		prototile_first->loadTileData("test50");
		Mat draw = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(draw, prototile_first->contour, Point2f(400, 400));
		imshow("asd", draw);
		vector<pair<int, bool>> all_total = tiling_opt->compare_choose_TAR(prototile_first->contour);
		cout << endl << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;*/

		//--------------------test morphing_tar by all------------
		//start = clock();
		//tiling_opt->load_dataset();
		//tiling_opt->com_all_TARs(1);
		////cout << "tiling_opt" << tiling_opt->contour_dataset.size()<<endl
		////	<< "tars: " << tiling_opt->all_con_tars.size()<<endl;
		//prototile_first->loadTileData("test50");
		//vector<Point2f> contour_inner = prototile_first->contour_sample[1];
		////cout << "prototile_first->contour" << prototile_first->contour.size() << "  contour_inner: " << contour_inner.size()<<endl;
		////midtime = clock();
		////cout << endl << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		//double sc_inner = 0;
		//vector<vector<double>> inner_tar = prototile_first->compute_TAR(contour_inner, sc_inner);
		//cout << "sc_inner: " << sc_inner << endl;
		//vector<pair<int, bool>> cand = tiling_opt->compare_choose_TAR(prototile_first->contour);
		////vector<pair<int, bool>> cand = tiling_opt->quick_choose_TAR(prototile_first->contour);
		//cout << "candsize: " << cand.size() << endl;
		//midtime = clock();
		//cout << endl << "All time consumption: " << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;

		//int t = 0;
		//prototile_second->loadPoints(tiling_opt->contour_dataset[cand[t].first]);
		//vector<pair<int, int>> path;
		//vector<vector<double>> tar_sec;
		//vector<Point2f> contour_cand;
		//int shift = 0;
		//if (cand[t].second)
		//{
		//	tar_sec = tiling_opt->all_con_tars_flip[cand[t].first];
		//	contour_cand = prototile_second->contour_sample_flip[1];
		//}
		//else
		//{
		//	tar_sec = tiling_opt->all_con_tars[cand[t].first];
		//	contour_cand = prototile_second->contour_sample[1];
		//}
		//cout << "first: " << cand[t].first << " second: " << cand[t].second << endl;
		//int width = 6;
		//double re = tiling_opt->tar_mismatch(inner_tar, tar_sec, path, shift, width);
		//finish = clock();
		//cout << endl << "once consumption: " << (double)(finish - midtime) / CLOCKS_PER_SEC << " s " << endl;
		//vector<int> mid_inter;
		//mid_inter.push_back(0);
		//mid_inter.push_back(60);
		//mid_inter.push_back(120);
		//mid_inter.push_back(175);
		//vector<Point2f> mor_result = tiling_opt->morphing_tar(contour_inner, contour_cand, mid_inter, path, shift);
		//start = clock();
		//cout << endl << "once morphing consumption: " << (double)(start - finish) / CLOCKS_PER_SEC << " s " << endl;
		////imgtocout();
		//cout << "inter_: " << mor_result.size() << endl;
		////MorphPoints(contour_inner, contour_cand, inter_, 0.5);
		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(ttt, mor_result, Point2f(400, 400));
		//draw_poly(drawing_pro1, contour_cand, Point2f(400, 400));
		//draw_poly(tt, contour_inner, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//for (int i = 0; i < 4; i++)
		//{
		//	circle(tt, contour_inner[mid_inter[i]], 5, Scalar(255, 0, 0), -1);
		//}
		////
		//imshow("result", ttt);
		//imshow("inner", tt);
		//imshow("cand", drawing_pro1);

		//--------------------test morphing_tar by one------------

  //      start = clock();
		//prototile_first->loadTileData("test50");
		//vector<Point2f> contour_inner = prototile_first->contour_sample[1];
		//double sc_inner = 0;
		//vector<vector<double>> inner_tar = prototile_first->compute_TAR(contour_inner, sc_inner);
		//prototile_second->loadTileData("485");
		//vector<Point2f> contour_cand = prototile_second->contour_sample[1];
		//double sc_cand = 0;
		//vector<vector<double>> cand_tar = prototile_second->compute_TAR(contour_cand, sc_cand);
		//vector<Point2f> contour_cand_f = prototile_second->contour_sample_flip[1];
		//double sc_cand_f = 0;
		//vector<vector<double>> cand_tar_f = prototile_second->compute_TAR(contour_cand_f, sc_cand_f);
		//vector<pair<int, int>> path;
		//int shift;	
		//int width = 4;

		//double re = tiling_opt->tar_mismatch(inner_tar, cand_tar, path, shift, width);
		//
		//double re2 = tiling_opt->tar_mismatch(inner_tar, cand_tar_f, path, shift, width);
		//if (re < re2)
		//{
		//	re2 = tiling_opt->tar_mismatch(inner_tar, cand_tar, path, shift, width);
		//	contour_cand_f = prototile_second->contour_sample[1];
		//}
		//cout << "re2: " << re2 <<"shift: "<<shift<< endl;
		//cout << "load over!" << endl;
		//vector<int> mid_inter;
		//mid_inter.push_back(0);
		//mid_inter.push_back(35);
		//mid_inter.push_back(105);
		//mid_inter.push_back(170);
		//vector<Point2f> mor_result = tiling_opt->morphing_tar(contour_inner, contour_cand_f, mid_inter, path, shift);
		//mid_inter.swap(vector<int>());
		//mid_inter.push_back(0);
		//mid_inter.push_back(35);
		//mid_inter.push_back(105);
		//mid_inter.push_back(170);
		////imgtocout();
		//finish = clock();
		//cout << endl << "once morphing consumption: " << (double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
		//cout << "inter_: " << mor_result.size() << endl;
		////MorphPoints(contour_inner, contour_cand, inter_, 0.5);
		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		////Mat drawing_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(ttt, mor_result, Point2f(400, 400));
		//draw_poly(drawing_pro1, contour_cand_f, Point2f(400, 400));
		//draw_poly(tt, contour_inner, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//for (int i = 0; i < 4; i++)
		//{
		//	circle(tt, contour_inner[mid_inter[i]], 5, Scalar(255, 0, 0), -1);
		//}
		////
		//vector<int> return_p;
		//vector<vector<Point2f>> four_;
		//vector<Point2f> morphed_A = tiling_opt->extract_contour(mor_result, mid_inter, return_p, four_, 0);
		//Mat drawing_ex = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		//Point2f sss = Point2f(400, 400) - center_p(morphed_A);
		//for (int t = 0; t < morphed_A.size(); t++)
		//{
		//	circle(drawing_ex, morphed_A[t] + sss, 2, Scalar(0, 0, 255), -1);
		//}
		////MyLine(drawing_ex, morphed_A[128] + sss, morphed_A[129] + sss, "blue");
		////MyLine(drawing_ex, morphed_A[131] + sss, morphed_A[132] + sss, "green");
		////circle(drawing_ex, morphed_A[128] + sss, 2, Scalar(0, 255, 0), -1);
		////circle(drawing_ex, morphed_A[129] + sss, 2, Scalar(0, 255, 0), -1);
		////circle(drawing_ex, morphed_A[131] + sss, 2, Scalar(255, 0, 0), -1);
		////circle(drawing_ex, morphed_A[132] + sss, 2, Scalar(255, 0, 0), -1);
		//draw_poly(drawing_ex, morphed_A, Point2f(1200, 400));
		//
		//int first = 0;
		//int second = 0;
		////if (self_intersect(morphed_A,first,second)) cout << "self_intersect" << endl;
		//while (self_intersect(morphed_A, first, second))
		//{
		//	cout << "self_intersect" << endl;
		//	//tiling_opt->contour_fine_tuning(morphed_A, return_p, first, second);
		//}
		//for (int t = 0; t < morphed_A.size(); t++)
		//{
		//	circle(drawing_ex, morphed_A[t] + sss, 2, Scalar(0, 255, 0), -1);
		//}
		//Mat final_re = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(final_re, morphed_A, Point2f(400, 400));
		//imshow("final_re: ", final_re);
		//cout << "zahuishia: "<<morphed_A.size() << endl;
		//imshow("morphed_A: ", drawing_ex);
		//imshow("morphing result", ttt);
		//imshow("contour_inner", tt);
		//imshow("cand_pattern", drawing_pro1);
        //-----------------------------------------------------------------

		//for (int i = 0; i < 202; i++)
		//{
		//	prototile_first->~Prototile();
		//	string image = int2string(i);
		//	int thresh = 100;
		//	int max_thresh = 255;
		//	String imageName("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + image + ".png");
		//	cout << image << endl;
		//	Mat src = imread(imageName);
		//	imshow("src", src);
		//	Mat src_g;
		//	cvtColor(src, src_g, COLOR_BGR2GRAY);
		//	blur(src_g, src_g, Size(3, 3));
		//	Mat canny_output;
		//	//考虑到可能有多个轮廓
		//	vector<vector<Point> > contours;
		//	vector<Vec4i> hierarchy;

		//	//用candy法由灰度图求出掩码图
		//	Canny(src_g, canny_output, thresh, thresh * 2, 3);
		//	//imshow("canny_output", canny_output);
		//	//由掩码图求出有序轮廓点
		//	findContours(canny_output, contours, hierarchy, CV_RETR_EXTERNAL, CHAIN_APPROX_SIMPLE, Point(0, 0));
		//	cout << "contours num:" << contours.size() << endl;
		//	vector<Point2f> sampling_;
		//	for (int t = 0; t < contours[0].size(); t++)
		//	{
		//		sampling_.push_back(contours[0][t]);
		//	}
		//	vector<double> sam_cur = curvature_com(sampling_);
		//	int cur_p_num = 30;
		//	vector<int> max_order;
		//	max_order = most_convex_p(sampling_, sam_cur, cur_p_num);
		//	int contoursize = sampling_.size();
		//	double c_length = arcLength(sampling_,true);
		//	vector<int> all_order = max_order;
		//	double dist = 1000;
		//	int index = 0;
		//	for (int i = 0; i < contoursize; i = i + 10)
		//	{
		//		int flag = 0;
		//		for (vector<int>::iterator it = max_order.begin(); it != max_order.end(); it++)
		//		{
		//			double leng = length_two_point2f(sampling_[i], sampling_[*it]);
		//			if (leng < 0.005*c_length)
		//			{
		//				flag = 1;
		//				break;
		//			}
		//		}

		//		if (flag == 0)
		//		{
		//			all_order.push_back(i);
		//		}
		//	}
		//	cout << "all_order.size:" << all_order.size() << endl;
		//	sort_bub(all_order);
		//	vector<Point2f> res;
		//	for (int j = 0; j < all_order.size(); j++)
		//	{
		//		res.push_back(sampling_[all_order[j]]);
		//	}
		//	Mat tt = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//	draw_poly(tt, res, Point2f(300, 300));
		//	imshow("??", tt);
		//	imwrite("D:\\VisualStudioProjects\\DihedralTessellation\\other\\" + image + "s.png", tt);
		//}
		//
	}
	else if (f == 5)
	{
		tiling_opt->load_dataset();
		tiling_opt->check_Repetitive_pattern();

	}
	else if (f == 6)
	{
		//------test 角点度数检测----------------------
		/*vector<Point2f> a;
		vector<Point2f> b;
		a.push_back(Point2f(2, 2));
		a.push_back(Point2f(5, 2));
		a.push_back(Point2f(2, 5));
		b.push_back(Point2f(2, 2));
		b.push_back(Point2f(2, 5));
		b.push_back(Point2f(5, 2.5));

		if (tiling_opt->vertex_angle(a, b)) cout << "pengzhuang";*/
		vector<Point2f> con_point;
		//读取一个存有轮廓点的文件，格式对应上一步计算轮廓点保存的文件
		string filepath = "D:\\VisualStudioProjects\\DihedralTessellation\\simulation\\307_morA_14.txt";
		ifstream in(filepath);
		if (!in.is_open())
		{
			cout << filepath << endl;
			cout << "Error opening file" << endl;
		}
		//挨个处理每个字符
		//cout << "Opening file!!!" << endl;
		vector<char> each_point;
		int aa = 0;
		int bb = 0;
		int nn = 0;
		char cc;
		char buf[200];
		//for ()
		in.getline(buf, 200);
		//cout << "num: " << buf << endl;
		while (!in.eof())
		{
			aa = 0;
			bb = 0;
			nn = 0;
			int f = 0;
			in.getline(buf, 200);
			cc = buf[nn++];
			while ((cc >= '0' && cc <= '9') || cc == ',' || cc == ' ')
			{
				f = 1;
				if ((cc >= '0' && cc <= '9'))
				{
					each_point.push_back(cc);
				}
				if (cc == ',')
				{
					for (int i = 0; i < each_point.size(); i++)
					{
						aa = aa * 10 + (each_point[i] - 48);
					}
					each_point.swap(vector<char>());
				}
				cc = buf[nn++];
			}
			if (f)
			{
				for (int i = 0; i < each_point.size(); i++)
				{
					bb = bb * 10 + (each_point[i] - 48);
				}
				each_point.swap(vector<char>());
				con_point.push_back(Point(aa, bb));
			}
		}
		in.close();
		Mat final_re = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(final_re, con_point, Point2f(400, 400));
		imshow("dads", final_re);

			

	}
	else if (f == 7)
	{
		//------------------------test translation, rotation and reflection----------------
		string imaname = "307";
		//prototile_first->imgtocout(imaname);
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		vector<Point2f> a = prototile_first->contour;
		vector<int> mark;
		//mark.push_back(0);
		
		mark.push_back(132);
		mark.push_back(234);
		mark.push_back(310);
		mark.push_back(440);
		Mat draw_p = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//Mat draw2 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		
		for (int i = 0; i < a.size(); i++)
		{
			circle(draw_p, a[i], 1, Scalar(0, 0, 0), -1);
		}
		for (int jj = 0; jj < p_p_index.size(); jj++)
		{
			circle(draw_p, a[p_p_index[jj]], 4, Scalar(0, 0, 255), -1);
		}
		for (int jj = 0; jj < 4; jj++)
		{
			circle(draw_p, a[mark[jj]], 6, Scalar(0, 255, 0), -1);
		}
		//Point2f t = a[479] - a[480];
		//Point2f tt = a[481] - a[480];
		//cout << t << " " << tt << "  " << cos_two_vector(t, tt) <<"  "<<sin_two_vector(t,tt)<< endl;;
		vector<int> midmark;
		vector<Point2f> return_B;
		vector<vector<Point2f>> four_place;
		//Mat draw1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		return_B = tiling_opt->extract_contour(a, mark, midmark, four_place,0);
		Mat drawing_ttt = Mat(2000, 2000, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift = Point2f(800, 800) - center_p(return_B);
		for (int i = 0; i < 4; i++)
		{
			draw_poly(drawing_ttt, four_place[i], center_p(four_place[i]) + shift);
		}
		draw_poly(drawing_ttt, return_B, center_p(return_B) + shift, 2);
		for (int i = 0; i < return_B.size(); i++)
			MyLine(drawing_ttt, return_B[i] + shift, return_B[(i + 1) % return_B.size()] + shift, "green2");
		imshow("extracted contour", drawing_ttt);
		imwrite("D:\\ex_contour.png", drawing_ttt);
		imshow("sample points", draw_p);
		imwrite("D:\\samplepoints.png", draw_p);
		//Mat rot_mat = getRotationMatrix2D(center_p(a), 50, 1);
		//transform(a, return_B, rot_mat);
		//if (tiling_opt->flipping_placement(mark, a, return_B, midmark, draw1,1)) cout << " pengzhuang " << endl;
		////if (tiling_opt->translation_placement(mark, a, return_B, midmark, draw1)) cout << "pngzhuang" << endl;
		////if (tiling_opt->rotation_placement(mark, a, return_B, midmark, draw1)) cout << " pengzhuang " << endl;
		//else
		//{
		//	cout << "return_B" << return_B.size()<< endl;
		//	draw_poly(draw2, return_B, center_p(return_B));
		//	//MyLine(draw1, return_B[12], return_B[144], "green");
		//	//imshow("hahah", draw1);
		//	//imwrite("D:\\XXXX.PNG", draw1);
		//}
		/*for (int i = 0; i < a.size(); i++)
			circle(draw2, a[i], 2, Scalar(0, 0, 255), -1);
		circle(draw2, a[12], 6, Scalar(0, 255, 0), -1);
		circle(draw2, return_B[12], 6, Scalar(0, 255, 0), -1);
		for (int j = 1; j < 4; j++)
		{
			circle(draw2, a[mark[j]], 4, Scalar(0, 255, 0), -1);
			circle(draw2, return_B[mark[j]], 4, Scalar(0, 255, 0), -1);
		}
		MyLine(draw2, a[12], a[144], "green");
		draw_poly(draw1, return_B, center_p(return_B));*/
		//imshow("hahah", draw1);
		//imshow("asdad", draw2);
		//---------------------------------------------------
        
	
		////简单的检测
		//cv::Mat ima = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname);
		//cv::Mat imb = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname1);
		//cv::Mat imc = Mat::zeros(imb.rows,imb.cols,imb.type());
		//vector<cv::Point2f> ima_p;
		//vector<cv::Point2f> imb_p;
		//vector<cv::Point2f> imc_p;
		//MyLine(ima, Point2f(100, 100), Point2f(150, 50),"red");
		//MyLine(ima, Point2f(150, 50), Point2f(200, 50), "red");
		//MyLine(ima, Point2f(200, 50), Point2f(250, 100), "red");


		//MyLine(imb, Point2f(100, 200), Point2f(150, 250), "blue");
		//MyLine(imb, Point2f(150, 250), Point2f(200, 250), "blue");
		//MyLine(imb, Point2f(200, 250), Point2f(250, 200), "blue");
		//cv::Point2f a = cv::Point2f(100,100);
		//cv::Point2f b = cv::Point2f(150, 100);
		//cv::Point2f c = cv::Point2f(200, 100);
		//cv::Point2f d = cv::Point2f(250, 100);
		//cv::Point2f aa = cv::Point2f(100, 200);
		//cv::Point2f bb = cv::Point2f(150, 200);
		//cv::Point2f cc = cv::Point2f(200, 200);
		//cv::Point2f dd = cv::Point2f(250, 200);
		////cv::Point2f e = cv::Point2f(130, 150);
		//ima_p.push_back(a);
		//ima_p.push_back(b);
		//ima_p.push_back(c);
		//ima_p.push_back(d);
		////ima_p.push_back(e);

		//imb_p.push_back(aa);
		//imb_p.push_back(bb);
		//imb_p.push_back(cc);
		//imb_p.push_back(dd);
		////imb_p.push_back(c);
		////imshow("hahahah",ima);
		//ImageMorphing(ima, ima_p, imb, imb_p, imc, imc_p,0.5,0.5);
		//imshow("ima",ima);
		//imshow("imb",imb);
		//imshow("imc",imc);


		//collision
		/*vector<Point2f> a;
		vector<Point2f> b;
		a.push_back(Point2f(0.5, 0.9));
		a.push_back(Point2f(1.2, 1.9));
		a.push_back(Point2f(2.5, 2.9));
		a.push_back(Point2f(3.5, 3.9));
		a.push_back(Point2f(5.5, 4.9));

		b.push_back(Point2f(0.2, 0.7));
		b.push_back(Point2f(1.2, 1.7));
		b.push_back(Point2f(2.2, 2.7));
		b.push_back(Point2f(3.2, 3.7));
		b.push_back(Point2f(5.2, 4.7));
		if (tiling_opt->collision_pixel(Point2f(10, 10), Point2f(0, 0), a, b))
		cout << "pengzhuang" << endl;
		else cout << "no";*/
	}
	else if (f == 8)
	{
		//show collision
		vector<int> max_order = prototile_first->partition_points("312");
		vector<Point2f> contours = prototile_first->contour;	
		int csize = contours.size();
		int results[4] = {0,150,261,451};
		Point2f line1 = contours[results[2]] - contours[results[0]];
		Point2f line2 = contours[results[3]] - contours[results[1]];
		
		vector<Point2f> one_loca;
		//translation
		for (int j = 0; j < csize; j++)
		{
			one_loca.push_back(contours[j] + line2);
		}
		//rotation
		/*Point2f rota_cent = contours[results[2]];
		Mat rot_mat = getRotationMatrix2D(rota_cent, 180, 1);
		transform(contours, one_loca, rot_mat);*/

		Mat drwa = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat drwa1 = Mat(1200, 1200, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i < contours.size(); i++)
		{
			circle(drwa, contours[i], 1, Scalar(0, 0, 0), -1);
		}
		/*for (int i = 0; i < 4; i++)
		{
			circle(drwa, contours[results[i]], 6, Scalar(0, 255, 0), -1);
		}*/
		
		for (int i = 0; i < max_order.size(); i++)
		{
			circle(drwa, contours[max_order[i]], 4, Scalar(0, 0, 255), -1);
		}
		circle(drwa, contours[results[1]], 6, Scalar(0, 255, 0), -1);
		circle(drwa, contours[results[3]], 6, Scalar(0, 255, 0), -1);
		imshow("contour", drwa);
		Point2f cent = center_p(contours);
		Point2f onecent = center_p(one_loca);
		Point2f shift = Point2f(600, 600) - 0.5*cent - 0.5*onecent;
		for (int i = 0; i < csize; i++)
		{
			contours[i] += shift;
			one_loca[i] += shift;
		}
		cent = cent + shift;// center_p(contours);
		onecent = onecent + shift;// center_p(one_loca);
		vector<Point2f> bbox1 = b_box(contours);
		vector<Point2f> bbox2 = b_box(one_loca);
		draw_poly(drwa1, contours, cent);
		draw_poly(drwa1, one_loca, onecent);
		for (int i = 0; i < csize; i++)
		{
			//circle(drwa1, contours[i], 2, Scalar(255, 255, 255), -1);
			MyLine(drwa1, contours[i], contours[(i + 1) % csize], "grey");
		}
		
		for (int i = 0; i < csize; i++)
		{
			//circle(drwa1, one_loca[i], 2, Scalar(255, 255, 255), -1);
			MyLine(drwa1, one_loca[i], one_loca[(i + 1) % csize], "grey");
		}
		circle(drwa1, one_loca[results[1]], 6, Scalar(0, 255, 0), -1);
		circle(drwa1, one_loca[results[3]], 6, Scalar(0, 255, 0), -1);
		circle(drwa1, contours[results[1]], 6, Scalar(0, 255, 0), -1);
		circle(drwa1, contours[results[3]], 6, Scalar(0, 255, 0), -1);
		//circle(drwa1, rota_cent + shift, 3, Scalar(0, 255, 0), -1);
		
	
		for (int i = 0; i < bbox1.size(); i++)
		{
			MyLine(drwa1, bbox1[i], bbox1[(i + 1) % bbox1.size()], "red");
			MyLine(drwa1, bbox2[i], bbox2[(i + 1) % bbox2.size()], "red");
		}
		imshow("draw3", drwa1);
		imwrite("D:\\pic.png", drwa1);

	}
	else if (f==9)
	{

		//tiling_opt->Tiling_clear();
		//vector<Point2f> sim_mid = tiling_opt->simulation_tar("302", 96, 0);

		Mat src, erosion_dst;
		int erosion_elem = 1;
		int erosion_size = 1;
		int const max_elem = 2;
		int const max_kernel_size = 21;
		
		string name = "D:\\model.png";//"D:\\result.png";
		//string name = "D:\\print.png";
		src = imread(name, IMREAD_COLOR);
		if (src.empty())
		{
			return -1;
		}
		
		namedWindow("Erosion Demo", WINDOW_AUTOSIZE);
		int erosion_type = 0;
		if (erosion_elem == 0){ erosion_type = MORPH_RECT; }
		else if (erosion_elem == 1){ erosion_type = MORPH_CROSS; }
		else if (erosion_elem == 2) { erosion_type = MORPH_ELLIPSE; }
		Mat element = getStructuringElement(erosion_type,
			Size(2 * erosion_size + 1, 2 * erosion_size + 1),
			Point(erosion_size, erosion_size));
		erode(src, erosion_dst, element);
		imwrite("D:\\Erosion Demo.png", erosion_dst);
		int row = erosion_dst.rows;
		int col = erosion_dst.cols;
		Mat flip1 = Mat(row, col, CV_8UC3, Scalar(0, 0, 0));
		Mat flip2 = Mat(row, col, CV_8UC3, Scalar(0, 0, 0));
		cout << row << " " << col << endl;
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
			{
				//src.at<Vec3b>(i, j) = Vec3b(255, 255, 255) - src.at<Vec3b>(i, j);
				if ((int)src.at<Vec3b>(i, j)[0] == 255 && (int)src.at<Vec3b>(i, j)[1] == 255 && (int)src.at<Vec3b>(i, j)[2] == 255)
				{
					//cout << "youle" << " ";
					flip1.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
				}
				if ((int)src.at<Vec3b>(i, j)[0] == 0 && (int)src.at<Vec3b>(i, j)[1] == 0 && (int)src.at<Vec3b>(i, j)[2] == 0)
				{
					//cout << "youle" << " ";
					flip2.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
				}
				//cout << (int)erosion_dst.at<Vec3b>(i, j)[0] << ',';
				if ((int)erosion_dst.at<Vec3b>(i, j)[0] == 255 && (int)erosion_dst.at<Vec3b>(i, j)[1] == 255 && (int)erosion_dst.at<Vec3b>(i, j)[2] == 255)
				{
					//cout << "youle" << " ";
					src.at<Vec3b>(i, j) = Vec3b(0, 0, 0);
				}
			}
		//erode(src, src, element);
		imwrite("D:\\printre.png", src);
		int erosion_size1 = 0;
		int erosion_size2 = 0;
		Mat element1 = getStructuringElement(erosion_type,
			Size(2 * erosion_size1 + 1, 2 * erosion_size1 + 1),
			Point(erosion_size1, erosion_size1));
		Mat element2 = getStructuringElement(erosion_type,
			Size(2 * erosion_size2 + 1, 2 * erosion_size2 + 1),
			Point(erosion_size2, erosion_size2));
		erode(flip1, flip1, element1);
		erode(flip2, flip2, element2);
		/*for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
			{
				flip1.at<Vec3b>(i, j) = Vec3b(255, 255, 255) - flip1.at<Vec3b>(i, j);
				flip2.at<Vec3b>(i, j) = Vec3b(255, 255, 255) - flip2.at<Vec3b>(i, j);
			}*/
		imwrite("D:\\flip1.png", flip1);
		imwrite("D:\\flip2.png", flip2);

	}
	else if (f == 10)
	{
		Mat src;
		string name = "D:\\input.png";
		src = imread(name, IMREAD_COLOR);
		if (src.empty())
		{
			return -1;
		}
		Mat dst = src;
		cvtColor(src, src, COLOR_BGR2GRAY);
	
		threshold(src, src, 128, 255, cv::THRESH_BINARY); //255 white
		//imshow("zala", src);
		//imshow("nine", dst);
		Vec3b color_p[20];
		color_p[0] = Vec3b(194, 194, 194);//gray
		color_p[1] = Vec3b(251, 228, 169); //blue1
		color_p[2] = Vec3b(251, 204, 176); //blue2
		color_p[3] = Vec3b(130, 174, 89); //green1
		color_p[4] = Vec3b(222, 250, 167); //green2
		color_p[5] = Vec3b(127, 110, 174); //red1
		color_p[6] = Vec3b(70, 124, 217); //orange1
		color_p[7] = Vec3b(251, 181, 105); //blue3
		color_p[8] = Vec3b(237, 171, 245); //pink
		color_p[9] = Vec3b(3, 142, 249); //orange2
		color_p[10] = Vec3b(175, 211, 249); //lightorange2
		color_p[11] = Vec3b(244, 211, 247); //lightpink
		//color_p[12] = Vec3b(237, 171, 245); //pink
		//color_p[13] = Vec3b(237, 171, 245); //pink

		for (int i = 0; i < src.rows; i++)
			for (int j = 0; j < src.cols; j++)
			{
				//cout << (int)src.at<uchar>(i, j) << " ";
				if ((int)src.at<uchar>(i, j) == 0) dst.at<Vec3b>(i, j) = color_p[6];
				if ((int)src.at<uchar>(i, j) == 255) dst.at<Vec3b>(i, j) = color_p[7];
			}
		imwrite("D:\\dst.png", dst);
	}

	else if (f == 11)
	{
		glutInit(&argc, argv);
		vector<int> p_p_index = prototile_first->partition_points("test14");
		vector<Point2f> conr = prototile_first->contour_sample[1];
		vector<int> p_p_index1 = prototile_second->partition_points("test");
		vector<Point2f> conr1 = prototile_second->contour_sample[1];
		int contoursize = conr.size();

		GLsizei wh = 800, ww = 1200;
		OpenWindow(ww, wh, conr, conr1);
		//glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
		//glutInitWindowSize(ww, wh);
		////glutInitWindowPosition(100, 120);
		//glutCreateWindow("test");
		////glutReshapeFunc(myReshape);

		//glutDisplayFunc(display);

		////	glutDisplayFunc(Display);
		////glEnable(GL_DEPTH_TEST);
		//glutMouseFunc(myMouse);
		//initial();
		//glutMainLoop();
		
	}
    finish = clock();
	cout << endl<< "All time consumption: "<<(double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
	waitKey(0);
	return 0;
}


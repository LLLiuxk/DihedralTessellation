#include "tilingOpt.h"
#include  <stdio.h>
#include  <stdlib.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_2.h>
//#include "stdafx.h"

using namespace Tiling_tiles;

int main(int argc, char** argv)
{
	clock_t start, midtime, finish;
	start = clock();

	Tiling_tiles::Tiling_opt *tiling_opt;
	tiling_opt = new Tiling_tiles::Tiling_opt();
	Tiling_tiles::Prototile *prototile_first;
	prototile_first = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_second;
	prototile_second = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_third;
	prototile_third = new Tiling_tiles::Prototile();
	//////prototile_first->imgtocout(imagename1);
	int f = 1;
	//0:result  1:simulation  2:批量读图  3:feature points  4:compute_TAR  5:min_minsmatch  6:extract_contour  7:compare and choose
	//8:morphing  9:draw  10:math  11:check  12:thickness  13:color  14:windows 15:evalua_deformation
	//17:2Dtriangle  18:求差 19：三角化  20:contour_dilate
	if (f == 111) //test
	{
		//Point2f a(600, 200);
		//vector<vector<Point2f>> contours = extract_contours("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\502.png");
		//Mat drawing_pro = Mat(20000, 20000, CV_8UC3, Scalar(255, 255, 255));
		//int outc_num = contours.size();
		//cout << " outc_num: " << outc_num << endl;
		//for (int j = 0; j < outc_num; j++)
		//{
		//	int out1size = contours[j].size();
		//	cout << " out1size: " << out1size << endl;
		//	for (int t = 0; t < out1size; t++)
		//	{
		//		//cout << " contours: " << contours[j][t] << endl;
		//		//circle(drawing_pro, contours[j][t], 10, Scalar(0, 255, 0), -1);
		//		MyLine(drawing_pro, contours[j][t], contours[j][(t + 1) % out1size], "red");
		//	}
		//}
		//imwrite("D:\\VisualStudioProjects\\DihedralTessellation\\halftone22.png", drawing_pro);
		////true表示点到轮廓的距离
		vector<Point2f> c1;
		c1.push_back(Point2f(200,200));
		c1.push_back(Point2f(400, 200));
		c1.push_back(Point2f(400, 400));
		//c1.push_back(Point2f(300, 400));
		vector<Point2f> c2;
		c2.push_back(Point2f(35, 35));
		c2.push_back(Point2f(70, 40));
		c2.push_back(Point2f(80, 80));
		for (int i = 0; i < c2.size(); i++)
		{
			c2[i] += Point2f(100, 100);
		}
		//c2.push_back(Point2f(500, 800));
		Mat M_seg1 = getAffineTransform(c1, c2);
		cout << "M_seg1: " << M_seg1 << endl;
		c1.push_back(Point2f(300, 400));
		cv::transform(c1, c2, M_seg1);
		for (int i = 0; i < c2.size(); i++)
		{
			cout << "c2: " << c2[i] << endl;
		}
		//cv::transform(seg_sam1, seg_sam1, M_seg1);
		//c1 = contours[0];
		//cout << "c1_size: " << c1.size() << endl;
		//double a0 = pointPolygonTest(c1, a, true);
		////false表示计算点与轮廓的位置关系-1表示外部，0在轮廓上，1在轮廓内
		//double b0 = pointPolygonTest(c1, a, false);
		//cout << a0 << "  :  " << b0 << endl;
		//Mat drawing5 = draw_polygen("hahahah", c1);
		//circle(drawing5,a,6,Scalar(0,255,0),-1);
		//imshow("rrr", drawing5);

		//--------------sampling_seg------------------
		/*vector<Point2f> gg = sampling_seg(ttt,10);
		int gs = gg.size();
		cout << gs << endl;
		for (int g = 0; g < gs; g++)
		{
			cout << gg[g] << endl;
		}*/
		

		/*Point2f cen1 = Point2f(10, 25);
		double angle = 90;
		Mat rot_mat(2, 3, CV_32FC1);
		rot_mat = getRotationMatrix2D(cen1, angle, 1);
		transform(ttt, ttt, rot_mat);
		cout << ttt[0] << "  " << ttt[1] << endl;*/

		//prototile_first->setname("swan2");
		//vector<Point2f> new_c = prototile_first->readTxt();
		//string file = "D:\\swan2.obj";
		//write_obj(file, new_c, 40);

		/*prototile_second->setname("gezi_m");
		vector<Point2f> new_c2 = prototile_second->readTxt();
		string file2 = "D:\\show\\4\\gezi_m.obj";
		
		write_obj(file2, new_c2, 20);*/

		//string fff = "hahahah.txt";
		//cout << fff.size()<<endl;
		//string t = fff.substr(fff.size() - 3, fff.size() - 1);
		//cout << t;
		//读一张图像数据
		//Mat drawing5 = Mat(1800, 900, CV_8UC3, Scalar(255, 255, 255));
		//vector<Point2f> p = { Point2f(10, 10), Point2f(20, 10), Point2f(20, 15), Point2f(25, 15), Point2f(20, 20), Point2f(10, 20) };
		//Line_Seg t(Point2f(20,0),Point2f(20,25));
		//line_polygon(t,p);
		//draw_poly(drawing5,new_c,Point2f(300,900));

		//vector<int> p_p_index = prototile_first->partition_points("test14");
		//Mat drawing5 = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//vector<Point2f> conr = prototile_first->contour_sample[1];
		//int contoursize = conr.size();
		////cout << contoursize << endl;
		//for (int j = 0; j < contoursize; j++)
		//{
		//	circle(drawing5, conr[j], 2, Scalar(0, 0, 0), -1);
		//	//MyLine(drawing5, conr[j], conr[(j + 1) % contoursize], "pink");
		//}
		//imshow("sample", drawing5);
	}
	if (f == 0) //已有dataset, 计算结果
	{

		tiling_opt->tiliing_generation("293");

		//批量计算
		//int iii[20] = {48,75,80,124,220,218,212,228,248,251,280,312,317};//{484,483,482,482,480,479,478,474,472,467,451}; 
		//for (int i = 5; i < 486; i++)
		//{
		//	tiling_opt->Tiling_clear();
		////	tiling_opt->tiliing_generation(int2string(iii[i]));
		////	
		//	system("cls");
		//}
	}
	if (f == 1)  //simulation
	{
		tiling_opt->Tiling_clear();
		int inner_one = 1;
		int cand_one = 1;
		jointPat sim_mid = tiling_opt->simulation_tar("293", inner_one, cand_one);  //

		//int mid = 0;
		//vector<Point2f> new_c = tiling_opt->construct_joint(sim_mid,mid);
		///*string file = "D:\\show\\model1.obj";
		//write_obj(file, new_c,30);
		//string file2 = "D:\\show\\cylinder.obj";
		//vector<Point2f> cylinder;
		//int r = 10;
		//for (int i = 0; i < 12; i++)
		//{
		//	cylinder.push_back(Point2f(r*cos(2 * PI / 12 * i), r*sin(2 * PI / 12 * i)));
		//}
		//double length = length_two_point2f(new_c[mid], new_c[0]);
		//write_obj(file2, cylinder, length);*/
		////draw four_
		//vector<Point2f> boxbox= b_box(new_c);
		//int raw = abs(boxbox[2].y - boxbox[0].y) + 600;
		//int col = abs(boxbox[2].x - boxbox[0].x) + 600;
		//Mat drawing_four = Mat( raw,col, CV_8UC3, Scalar(255, 255, 255));

		//draw_poly(drawing_four, new_c, Point2f(col / 2, raw / 2));
		//Point2f shift = Point2f(col/2, raw/2) - center_p(new_c);
		////circle(drawing_four, new_c[0] + shift, 3, Scalar(0, 0, 255), -1);
		////circle(drawing_four, new_c[mid] + shift, 3, Scalar(0, 0, 255), -1);
		//imshow("aaa", drawing_four);
		//imwrite("D:\\show\\4\\111.png", drawing_four);

		//Point2f shh = 0.25 * (center_p(sim_mid.four_contour[0]) + center_p(sim_mid.four_contour[1]) + center_p(sim_mid.four_contour[2]) + center_p(sim_mid.four_contour[3]));
		//shh = Point2f(800, 800) - shh;
		//cout << "mid:" << sim_mid.interval[0] << " " << sim_mid.interval[1] << " " << sim_mid.interval[2] << " " << sim_mid.interval[3] << endl;
		//cout << "type: " << sim_mid.type << endl;
		//for (int four_i = 0; four_i < 4; four_i++)
		//{
		//	for (int t = 0; t < 4; t++)
		//	{
		//		circle(drawing_four, sim_mid.four_contour[four_i][sim_mid.interval[t]] + shh, 3, Scalar(0, 0, 250), -1);
		//	}
		//	//draw_poly(drawing_four, sim_mid.four_contour[four_i], center_p(sim_mid.four_contour[four_i]) + shh);
		//	string filepathname = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\joint\\" +
		//		int2string(inner_one) + "_" + int2string(cand_one) + "_four_" + int2string(four_i) + ".txt";
		//	vector<Point> con;
		//	for (int it = 0; it < sim_mid.four_contour[four_i].size()/3; it++)
		//	{
		//		circle(drawing_four, sim_mid.four_contour[four_i][it] + shh, 2, Scalar(0), -1);
		//		//con.push_back((Point)sim_mid.four_contour[four_i][it]);
		//	}
		//	for (int it = sim_mid.four_contour[four_i].size() / 3; it < sim_mid.four_contour[four_i].size() /2; it++)
		//	{
		//		circle(drawing_four, sim_mid.four_contour[four_i][it] + shh, 2, Scalar(255,0,0), -1);
		//		//con.push_back((Point)sim_mid.four_contour[four_i][it]);
		//	}
		//	for (int it = sim_mid.four_contour[four_i].size() / 2; it < sim_mid.four_contour[four_i].size(); it++)
		//	{
		//		circle(drawing_four, sim_mid.four_contour[four_i][it] + shh, 2, Scalar(0, 255, 0), -1);
		//		//con.push_back((Point)sim_mid.four_contour[four_i][it]);
		//	}
		//	//fileout(filepathname, con);
		//}
		//imshow("four", drawing_four);


		
	}
	if (f == 2) //批量读图
	{
		string image = "502";//int2string(i);

		prototile_first->txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";
		prototile_first->dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";
		prototile_first->imgtocout(image); //default  raw=0
		
		//prototile_first->setname(image);
		//string filepath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\aaa.txt";
		//prototile_first->contour=prototile_first->readTxt();
		//vector<Point2f> ta = sampling_ave(prototile_first->contour,500);
		//vector<Point> aa;
		//for (int g = 0; g < ta.size(); g++)
		//	aa.push_back((Point)ta[g]);
		//fileout(filepath, aa);

		//读取txt文件保存为新的png
		//prototile_first->contourname = image;
		//prototile_first->contour = prototile_first->readTxt();
		//Mat drawing_ = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(drawing_, prototile_first->contour, Point2f(300, 300));
		//imwrite("D:\\VisualStudioProjects\\DihedralTessellation\\datasetnew\\" + image + ".png", drawing_);

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
	}
	if (f == 3)
	{
		// ________________feature_points___________________________

		Mat drawing2 = Mat(600, 600, CV_8UC1, Scalar(255));
		prototile_first->loadTileData("6");
		vector<Point2f> a = prototile_first->contour_sample[2];
		Point2f cent = center_p(a);	
		draw_poly(drawing2, a, cent);
		circle(drawing2, Point2f(300,300), 2, Scalar(255), -1);
		circle(drawing2, cent, 2, Scalar(255), -1);
		imshow("as", drawing2);
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

	}
	if (f == 4)
	{
		//-----------------test compute_TAR----------------------

		vector<Point2f> a;
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
		Mat drawing_ = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(drawing_, a, Point2f(150, 150));

		//double shape_com;
		//vector<vector<double>> tar_all = prototile_first->compute_TAR(a, shape_com,0.5);
		//cout << "shape_com" << shape_com<<endl;
		////double b = sin_2vector_convexc(Point2f(5,0), Point2f(0,5));
		////cout << b << endl;
		//for (int i = 0; i < tar_all.size(); i++)
		//{
		//	for (int j = 0; j < tar_all[i].size();j++)
		//		cout << tar_all[i][j] << "  ";
		//	cout << endl;
		//}
		imshow("aaaa", drawing_);	
		prototile_second->loadTileData("157");
		double shape_com;
		double shape_com_flip;
		vector<vector<double>> tar_all = prototile_second->compute_TAR(prototile_second->contour_sample[1], shape_com);//(num_c+1)*100 points
		vector<vector<double>> tar_all_flip = prototile_second->compute_TAR(prototile_second->contour_sample_flip[1], shape_com_flip);
		cout << "shape_com: " << shape_com << "  shape_com_flip: " << shape_com_flip << endl;


	}
	if (f == 5)
	{
		//-------------------test tar_mismatch函数--------------------

		tiling_opt->load_dataset();
		tiling_opt->com_all_TARs(1);
		prototile_first->loadTileData("15");
		vector<Point2f> a = prototile_first->contour;// prototile_first->contour_sample[1];
		vector<Point2f> aaa = prototile_first->contour_sample_flip[1];
		Mat drawing_ = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));
		int asize = a.size();
		cout << asize << endl;
		for (int i = 0; i < asize; i++)
		{
			circle(drawing_, a[i], 2, Scalar(0, 255, 0), -1);
		}
		imshow("asa", drawing_);
		//a.pop_back();
		double sca;
		vector<vector<double>> tar_all = prototile_first->compute_TAR(a, sca);
		//double scaaa;
		//vector<vector<double>> tar_allaa = prototile_first->compute_TAR(aaa, scaaa);
		/*vector<int> cand_points_index = most_convex_p(a, curvature_com(a), 30);
		cout << "feature :" << cand_points_index.size() << endl;
		vector<vector<double>> tar_fea;
		for (int j = 0; j < cand_points_index.size(); j++)
		{
			tar_fea.push_back(tar_all[cand_points_index[j]]);
		}*/
		/*prototile_second->loadTileData("23");
		vector<Point2f> b = prototile_second->contour_sample[1];
		vector<Point2f> c = prototile_second->contour_sample_flip[1];
		double scb;
		vector<vector<double>> tar_all1 = prototile_second->compute_TAR(b, scb);
		vector<vector<double>> tar_all2 = prototile_second->compute_TAR(c, scb);*/
		int num = 122;
		prototile_second->loadPoints(tiling_opt->contour_dataset[num]);
		vector<Point2f> b = prototile_second->contour_sample[1];
		vector<Point2f> c = prototile_second->contour_sample_flip[1];
		double scb = tiling_opt->all_shape_complexity[num];
		cout << "sca: " << sca <<  "   scb:  " << scb << endl;
		vector<vector<double>> tar_all1 = tiling_opt->all_con_tars[num];//prototile_second->compute_TAR(b, scb);
		vector<vector<double>> tar_all2 = tiling_opt->all_con_tars_flip[num];//prototile_second->compute_TAR(c, scb);
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
		vector<pair<int, int>> path;
		//double re = length_two_point_tar(tar_all[0], tar_all1[0]);
		int shift = 0;
		double re = tiling_opt->tar_mismatch(tar_all, tar_all1, path, shift);//点对应匹配的筛选框宽度
		cout << "result: " << re << "  shift:" << shift << endl;
		double re2 = tiling_opt->tar_mismatch(tar_all, tar_all2, path, shift);
		cout << "result2: " << re2 << "  shift:" << shift << endl;
		/*double re3 = tiling_opt->tar_mismatch(tar_fea, tar_fea2, path, shift);
		cout << "result3: " << re3 << "  shift:" << shift << endl;
		double re4 = tiling_opt->tar_mismatch(tar_fea, tar_fea_flip2, path, shift);
		cout << "result4: " << re4 << "  shift:" << shift << endl;*/
		re = re / (1 + sca + scb);
		re2 = re2 / (1 + sca + scb);
		//re3 = re3 / (1 + sca + scb);
		//re4 = re4 / (1 + sca + scb);
		
		cout << "result: " <<re<<"  result2: "<< re2 << "  shift:" << shift << endl;
		//cout << "result3: " << re3 << "  result4: " << re4 << "  shift:" << shift << endl;
		cout << path.size() << endl;
		for (int i = 0; i < path.size(); i++)
		{
		  cout << path[i].first << "--"<<path[i].second << endl;
		}

	}
	if (f == 6)
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
		return_B = tiling_opt->extract_contour(a, mark, midmark, four_place, 0);
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


	}
	if (f == 7)
	{
		//---------------  test compare_choose_TAR-------------------

		tiling_opt->load_dataset();
		tiling_opt->com_all_TARs(1);
		midtime = clock();

		prototile_first->loadTileData("test50");
		Mat draw = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(draw, prototile_first->contour, Point2f(400, 400));
		imshow("asd", draw);
		vector<pair<int, bool>> all_total = tiling_opt->compare_choose_TAR(prototile_first->contour);
		cout << endl << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		//-------------------------------------------

	}
	if (f == 8)
	{
		//--------------------test morphing_tar from all------------
		start = clock();
		tiling_opt->load_dataset();
		tiling_opt->com_all_TARs(1);
		prototile_first->loadTileData("test50");
		vector<Point2f> contour_inner = prototile_first->contour_sample[1];
		//cout << "prototile_first->contour" << prototile_first->contour.size() << "  contour_inner: " << contour_inner.size()<<endl;
		//midtime = clock();
		//cout << endl << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		double sc_inner = 0;
		vector<vector<double>> inner_tar = prototile_first->compute_TAR(contour_inner, sc_inner);
		cout << "sc_inner: " << sc_inner << endl;
		vector<pair<int, bool>> cand = tiling_opt->compare_choose_TAR(prototile_first->contour);
		//vector<pair<int, bool>> cand = tiling_opt->quick_choose_TAR(prototile_first->contour);
		cout << "candsize: " << cand.size() << endl;
		midtime = clock();
		cout << endl << "All time consumption: " << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;

		int t = 0;
		prototile_second->loadPoints(tiling_opt->contour_dataset[cand[t].first]);
		vector<pair<int, int>> path;
		vector<vector<double>> tar_sec;
		vector<Point2f> contour_cand;
		int shift = 0;
		if (cand[t].second)
		{
			tar_sec = tiling_opt->all_con_tars_flip[cand[t].first];
			contour_cand = prototile_second->contour_sample_flip[1];
		}
		else
		{
			tar_sec = tiling_opt->all_con_tars[cand[t].first];
			contour_cand = prototile_second->contour_sample[1];
		}
		cout << "first: " << cand[t].first << " second: " << cand[t].second << endl;
		int width = 6;
		double re = tiling_opt->tar_mismatch(inner_tar, tar_sec, path, shift, width);
		finish = clock();
		cout << endl << "once consumption: " << (double)(finish - midtime) / CLOCKS_PER_SEC << " s " << endl;
		vector<int> mid_inter;
		mid_inter.push_back(0);
		mid_inter.push_back(60);
		mid_inter.push_back(120);
		mid_inter.push_back(175);
		vector<Point2f> mor_result = tiling_opt->morphing_tar(contour_inner, contour_cand, mid_inter, path, shift);
		start = clock();
		cout << endl << "once morphing consumption: " << (double)(start - finish) / CLOCKS_PER_SEC << " s " << endl;
		//imgtocout();
		cout << "inter_: " << mor_result.size() << endl;
		//MorphPoints(contour_inner, contour_cand, inter_, 0.5);
		Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(ttt, mor_result, Point2f(400, 400));
		draw_poly(drawing_pro1, contour_cand, Point2f(400, 400));
		draw_poly(tt, contour_inner, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		for (int i = 0; i < 4; i++)
		{
			circle(tt, contour_inner[mid_inter[i]], 5, Scalar(255, 0, 0), -1);
		}
		//
		imshow("result", ttt);
		imshow("inner", tt);
		imshow("cand", drawing_pro1);

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

	}
	if (f == 9)   //draw tool
	{
		//--------------------test draw_allplane----------------------------
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
		imshow("hahah", show);
		//imshow("aahahah", show1);
		//vector<Point2f> d = tiling_opt->extract_contour(c, ttt, tttt, b, 3);
		//draw_allplane(show2, d, tttt, 1, 3);
		//imshow("final", show2);


	}
	if (f == 10) //math tool
	{
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
	}
	if (f == 11) //check
	{
		tiling_opt->load_dataset();
		tiling_opt->check_Repetitive_pattern();
		//show collision
		vector<int> max_order = prototile_first->partition_points("312");
		vector<Point2f> contours = prototile_first->contour;
		int csize = contours.size();
		int results[4] = { 0, 150, 261, 451 };
		Point2f line1 = contours[results[2]] - contours[results[0]];
		Point2f line2 = contours[results[3]] - contours[results[1]];


	}
	if (f == 12)  //增加连接区域厚度
	{
		//tiling_opt->Tiling_clear();
		//vector<Point2f> sim_mid = tiling_opt->simulation_tar("302", 96, 0);
		//腐蚀白色区域
		Mat src, erosion_dst;
		int erosion_elem = 2;
		int erosion_size = 6;
		int const max_elem = 2;
		int const max_kernel_size = 21;

		string name = "L:\\mmm.png";//"D:\\result.png";
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
			cv::Point(erosion_size, erosion_size));
		erode(src, erosion_dst, element);
		imwrite("L:\\Erosion Demo.png", erosion_dst);

		int row = erosion_dst.rows;
		int col = erosion_dst.cols;
		Mat flip1 = Mat(row, col, CV_8UC3, Scalar(0, 0, 0));
		Mat flip2 = Mat(row, col, CV_8UC3, Scalar(0, 0, 0));
		cout << row << " " << col << endl;
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
			{
				//cout << src.at<Vec3b>(i, j);// = Vec3b(255, 255, 255) - src.at<Vec3b>(i, j);
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
		imwrite("L:\\printre.png", src);

		int erosion_size1 = 3;
		int erosion_size2 = 3;
		Mat element1 = getStructuringElement(erosion_type,
			Size(2 * erosion_size1 + 1, 2 * erosion_size1 + 1),
			cv::Point(erosion_size1, erosion_size1));
		Mat element2 = getStructuringElement(erosion_type,
			Size(2 * erosion_size2 + 1, 2 * erosion_size2 + 1),
			cv::Point(erosion_size2, erosion_size2));
		erode(flip1, flip1, element1);
		erode(flip2, flip2, element2);
		/*for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
		{
		flip1.at<Vec3b>(i, j) = Vec3b(255, 255, 255) - flip1.at<Vec3b>(i, j);
		flip2.at<Vec3b>(i, j) = Vec3b(255, 255, 255) - flip2.at<Vec3b>(i, j);
		}*/
		imwrite("L:\\flip1.png", flip1);
		imwrite("L:\\flip2.png", flip2);

	}

	if (f == 13)  //换成不同配色
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

		threshold(src, src, 155, 255, cv::THRESH_BINARY); //255 white
		//imshow("zala", src);
		//imshow("nine", dst);
		Vec3b color1, color2;
		//Scalar t = colorbar[6].second;
		//color1 = Vec3b(t.val[0], t.val[1], t.val[2]);
		//t = colorbar[7].second;
		//color2 = Vec3b(t.val[0], t.val[1], t.val[2]);

		color1 = Vec3b(130, 149, 225);
		color2 = Vec3b(215, 223, 233);

		for (int i = 0; i < src.rows; i++)
			for (int j = 0; j < src.cols; j++)
			{
				//cout << (int)src.at<uchar>(i, j) << " ";
				if ((int)src.at<uchar>(i, j) == 0) dst.at<Vec3b>(i, j) = color1;
				if ((int)src.at<uchar>(i, j) == 255) dst.at<Vec3b>(i, j) = color2;
			}
		imwrite("D:\\dst.png", dst);
	}
	if (f == 14)
	{
		glutInit(&argc, argv);
		vector<int> p_p_index = prototile_first->partition_points("test14");
		vector<Point2f> conr = prototile_first->contour_sample[1];
		vector<int> p_p_index1 = prototile_second->partition_points("test");
		vector<Point2f> conr1 = prototile_second->contour_sample[1];
		int contoursize = conr.size();

		//GLsizei wh = 800, ww = 1200;
		//OpenWindow(ww, wh, conr, conr1);
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
	if (f == 15)
	{
		vector<Point2f> ww;
		ww.push_back(Point2f(10, 10));
		ww.push_back(Point2f(10, 110));
		ww.push_back(Point2f(110, 110));
		ww.push_back(Point2f(110, 10));
		prototile_first->loadTileData("161");
		//prototile_second->loadTileData("308");
		//tiling_opt->evalua_deformation(prototile_first->contour, prototile_second->contour);
		compute_TAR_new(prototile_first->contour_sample[1]);//prototile_first->contour_sample[1]
	}
	if (f == 16)
	{
		//vector<Point2f> sim_mid = tiling_opt->simulation_tar("302", 96, 0);
		//腐蚀白色区域
		Mat src, dilation_dst;
		int dilation_elem = 2;
		int dilation_size = 6;
		//int const max_elem = 2;
		//int const max_kernel_size = 21;

		string name = "L:\\bmm.png";//"D:\\result.png";
		//string name = "D:\\print.png";
		src = imread(name, IMREAD_COLOR);
		if (src.empty())
		{
			return -1;
		}
		int dilation_type = 0;
		namedWindow("Erosion Demo", WINDOW_AUTOSIZE);
		if (dilation_elem == 0){ dilation_type = MORPH_RECT; }
		else if (dilation_elem == 1){ dilation_type = MORPH_CROSS; }
		else if (dilation_elem == 2) { dilation_type = MORPH_ELLIPSE; }
		Mat element = getStructuringElement(dilation_type,
			Size(2 * dilation_size + 1, 2 * dilation_size + 1),
			cv::Point(dilation_size, dilation_size));
		dilate(src, dilation_dst, element);
		imwrite("L:\\Erosion Demo2.png", dilation_dst);

		//int row = dilation_dst.rows;
		//int col = dilation_dst.cols;
	
	}
	if (f == 17)
	{
		string pathname = "D:/2dtriangle.txt";
		vector<Point2f> ttt;
		vector<vector<int>> ttt4;
		read_2dtriangle(pathname, ttt, ttt4);

		cout << ttt4[0][0] << "  " << ttt4[0][1] << "  " << ttt4[0][2] << endl;
		Mat drawing_pro = Mat(2000, 2000, CV_8UC3, Scalar(255, 255, 255));
		vector<Point2f> bbx1 = b_box(ttt);
		double l_max = max(abs(bbx1[0].y - bbx1[1].y), abs(bbx1[1].x - bbx1[2].x));
		double sca = 1800 / l_max;
		Point2f shift = Point2f(10, 10) - sca*bbx1[1];
		for (int i = 0; i < ttt.size(); i++)
		{
			//cout << ttt[i] << endl;
			ttt[i] = ttt[i] * sca + shift;
		}
		for (int i = 0; i < ttt4.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				MyLine(drawing_pro, ttt[ttt4[i][j]], ttt[ttt4[i][(j + 1) % 3]], "black");
			}
		}
		imshow("2dtriangle", drawing_pro);
		imwrite("D://2dtriangle.png", drawing_pro);
	}
	if (f == 18)
	{
		Polygon_set_2 S;
		Polygon_2 ttt;
		ttt.push_back(Point_2(0, 0));
		ttt.push_back(Point_2(1, 0));
		ttt.push_back(Point_2(0, 2));
		Polygon_2 t;
		t.push_back(Point_2(2, 0));
		t.push_back(Point_2(4, 0));
		t.push_back(Point_2(2, 2));

		std::ostringstream os;
		os << ttt[0];
		std::string str = "100,100";//= os.str();
		cout << str<< endl;
		double aa = stod(str);
		cout << aa << endl;
		S.insert(ttt);
		S.insert(t);
		Polygon_2 rect2;
		rect2.push_back(Point_2(1, 0));
		rect2.push_back(Point_2(3, 0));
		rect2.push_back(Point_2(3, 2));
		rect2.push_back(Point_2(1, 2));
		S = A_difference_B(S, rect2);


		std::list<Polygon_with_holes_2> res;
		std::list<Polygon_with_holes_2>::const_iterator it;
		S.polygons_with_holes(std::back_inserter(res));
		for (it = res.begin(); it != res.end(); ++it) {
			std::cout << "--> ";
			print_polygon_with_holes(*it);
		}
	}
	if (f == 19)
	{
		//std::vector<Point> points = { Point(0, 0), Point(1, 0), Point(0, 1) };
		//Triangulation T;
		//T.insert(points.begin(), points.end());
		//std::cout << "Triangulation_2::Finite_vertices_iterator is like a  Triangulation_2::Vertex_handle\n";
		//for (Finite_vertices_iterator it = T.finite_vertices_begin();
		//	it != T.finite_vertices_end();
		//	++it){
		//	std::cout << it->point() << std::endl;
		//}
		//std::cout << "Triangulation_2::Finite_vertex_handles::iterator dereferences to Triangulation_2::Vertex_handle\n";
		//Finite_vertex_handles::iterator b, e;
		//std::tie(b, e) = T.finite_vertex_handles();
		//for (; b != e; ++b){
		//	Vertex_handle vh = *b; // you must dereference the iterator to get a handle
		//	std::cout << vh->point() << std::endl;
		//}

		//std::cout << "and you can use a C++11 for loop\n";
		//for (Vertex_handle vh : T.finite_vertex_handles()){
		//	std::cout << vh->point() << std::endl;
		//}
	}
	if (f == 20)
	{
		Polygon2 poly;
		poly.push_back(Point2(100, 100));
		poly.push_back(Point2(105, 100));
		poly.push_back(Point2(195, 100));
		poly.push_back(Point2(200, 100));
		poly.push_back(Point2(200, 105));
		poly.push_back(Point2(200, 195));
		poly.push_back(Point2(200, 200));
		poly.push_back(Point2(100, 200));
		
		if (poly.is_simple()) cout << "poly.is_simple();" << endl;
		ifstream in("D://swan_after.txt");
		if (!in) cout << "file error!" << endl;
		
		//in >> poly;

		//cout <<  poly<<endl;
		Polygon2 off_c = offset_poly(-20, poly);//正值为腐蚀
		cout << "off_c.size!" << off_c.size() << endl;
		vector<Point2f> vec_c = Polygon2vector(poly);
		if (contour_is_simple(vec_c)) cout << "contour_is_simple" << endl;
		cout << vec_c.size()<< endl;
		vector<Point2f> vec_c2 = contour_dilate(vec_c,17);//Polygon2vector(off_c);
		cout <<vec_c2.size() << endl;
		Mat drawing_pro = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(drawing_pro, vec_c, Point2f(400,400));
		draw_poly(drawing_pro, vec_c2, Point2f(1200, 400));
		imwrite("D://asd.png", drawing_pro);
		

		//Polygon2 out = offset_poly(15, poly);
		//int p_size = out.size();
		//cout << "out_contour" << endl;
		//for (int i = 0; i < p_size; i++)
		//{
		//	cout << out[i] << endl;
		//}
		//vector<Point2f> tile_erode = contour_erode(all_tiles[m], 8);
		////draw_poly(drawing_tesse, tile_dilate, center_p(tile_dilate));
		////draw_poly(drawing_tesse, tile_erode, center_p(tile_erode), 1);
		
	}
	finish = clock();
	cout << endl << "All time consumption: " << (double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
	waitKey(0);
	return 0;
}


#include "tilingOpt.h"
#include<ctime>
using namespace Tiling_tiles;

int main(int argc, char** argv)
{
	clock_t start, finish;
	start = clock();
	
	string imagename1 = "22"; 
	string imagename2 = "1";
	//string imagename1 = "Boat";
	//string imagename2 = "fish8";
	//string imagename1 = "fish5";
	string imagename3 = "19";
	//string txtname = "D:/images/111.png";
	//string txtname1 = "D:/images/fish3.png";

	//vector<Point2f> ske_points;
	//ske_points = get_Skeleon(imagename1);

	Tiling_tiles::Tiling_opt *tiling_opt;
	tiling_opt = new Tiling_tiles::Tiling_opt();
	Tiling_tiles::Prototile *prototile_first;
	prototile_first = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_second;
	prototile_second = new Tiling_tiles::Prototile();
	Tiling_tiles::Prototile *prototile_third;
	prototile_third = new Tiling_tiles::Prototile();
	//////prototile_first->imgtocout(imagename1);
	int f = 3;
	if (f == 0) //已有dataset
	{

		for (int i = 404; i < 409; i++)
		{ 
			tiling_opt->Tiling_clear();
			string name_ = int2string(i);
			cout << name_ << endl;
			tiling_opt->points_dividing(name_);
			system("cls");
		}
		
		//test morphing
		//--------------------
		//tiling_opt->load_dataset();
		//prototile_first->loadTileData("test");
	
		//vector<Point2f> contour_inner = prototile_first->contour_sample[1];
		//vector<CandPat> candida_contours;
		//candida_contours = tiling_opt->compare_shapes(prototile_first->contour, 1);
		//CandPat tem = candida_contours[0];
		//prototile_second->~Prototile();
		//prototile_second->loadPoints(tiling_opt->contour_dataset[tem.number]);
		//vector<Point2f> contour_cand = tiling_opt->CandP2Contour(tem, 1);
		//vector<Point2f> inter_;
		//vector<int> mid_inter;
		//inter_ = tiling_opt->morphing_2_patterns(contour_inner, contour_cand, mid_inter, 0.5);
		////MorphPoints(contour_inner, contour_cand, inter_, 0.5);


		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//Mat drawing_ = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(ttt, inter_, Point2f(400, 400));
		//draw_poly(drawing_pro1, contour_cand, Point2f(400, 400));
		//draw_poly(tt, contour_inner, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//
		//imshow("only point", ttt);
		//imshow("aaab", tt);
		//imshow("aaa", drawing_pro1);




		//vector<Point2f> inter_mid;
		//ImageMorphing(tt, contour_inner, drawing_pro1, contour_cand, drawing_dst, inter_mid, 0.5);
		//imshow("dst", drawing_dst);
		//draw_poly(drawing_, inter_mid, Point2f(400, 400));
		//imshow("dst_p", drawing_);
		//vector<Point2f> inter_mid = morphing_2_patterns(contour_inner, contour_cand, 0.5);
		
		/*
		Mat drawing_src1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_src2 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_src3 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_dst = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		prototile_first->loadTileData("test");
		prototile_second->loadTileData("128");
		vector<Point2f> src1_points;
		vector<Point2f> src2_points;
		vector<Point2f> output_final;

		vector<Point2f> a;
		vector<Point2f> b;
		src1_points = prototile_first->contour_sample[0];
		src2_points = prototile_second->contour_sample[0];
		src1_points.pop_back();
		int size = src1_points.size();

		for (int i = 100; i < 150; i++)
		{
			a.push_back(src1_points[i]);
			b.push_back(src2_points[i]);
		}
		Mat rot_mat;
		rot_mat = getRotationMatrix2D(a[10], 45, 1);
		transform(a, a, rot_mat);
		int asize = a.size();
		//for (int i = 0; i < size; i++)
		//{
		//	src1_points[i].x = src1_points[i].x +100;
		//	src1_points[i].y = src1_points[i].y +100;
		//}
		//src1_points.push_back(Point2f(100, 100));
		//src1_points.push_back(Point2f(200, 100));
		//src1_points.push_back(Point2f(200, 200));
		//src1_points.push_back(Point2f(150, 300));
		//src1_points.push_back(Point2f(100, 200));

		//src2_points.push_back(Point2f(300, 100));
		//src2_points.push_back(Point2f(400, 100));
		//src2_points.push_back(Point2f(400, 200));
		//src2_points.push_back(Point2f(350, 220));
		//src2_points.push_back(Point2f(300, 200));
		//cout << src1_points.size() << "   " << src2_points.size() << endl;
		for (int i = 0; i < asize-1; i++)
		{
			MyLine(drawing_src1, a[i], a[i + 1], "red");
			MyLine(drawing_src2, b[i], b[i + 1], "green");
			MyLine(drawing_src3, a[i], a[i + 1], "red");
			MyLine(drawing_src3, b[i], b[i + 1], "green");
			
		}
		//drawing_src1 = draw_polygen("drawing_src1", a);
		//drawing_src2 = draw_polygen("drawing_src2", b);
		ImageMorphing(drawing_src1, a, drawing_src2, b, drawing_dst, output_final, 0.5);
		imshow("out1: ", drawing_dst);
		for (int i = 0; i < output_final.size()-1; i++)
		{
			cout << output_final[i] << endl;
			MyLine(drawing_dst, output_final[i], output_final[i + 1], "black");
			MyLine(drawing_src3, output_final[i], output_final[i + 1], "black");
		}
		//drawing_dst = draw_polygen("out: ", output_final);
		imshow("1: ", drawing_src1);
		imshow("2: ", drawing_src2);
		imshow("3: ", drawing_src3);
		imshow("out: ", drawing_dst);
		//------------------------

		//------------------------
		//测试min_mismatch函数
		//prototile_first->loadTileData(imagename3);
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
		//prototile_second->loadTileData(imagename1);
		//vector<Point2f> con2 = prototile_second->contour_sample[1];
		//vector<double> con_c = prototile_second->contour_curva[1];
		////if (b.size() == b_c.size()) cout << "b.size:" << b.size() << endl;
		////if (con2.size() == con_c.size()) cout << "b.size:" << con2.size() << endl;
	 //   CandPat hahah = tiling_opt->min_mismatch(con2, a, con_c, a_c);
		//CandPat hahah1 = tiling_opt->min_mismatch(con2, b, con_c, b_c);

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
    	//if (hahah.mismatch < hahah1.mismatch)
		//{
		//	cout << "flip" << endl;
		//	cout << "angle: " << hahah.angle << "  index: " << hahah.index << "  mismatch: " << hahah.mismatch << endl;
		//}
		//else
		//{
		//	cout << "right" << endl;
		//	cout << "angle: " << hahah1.angle << "  index: " << hahah1.index << "  mismatch: " << hahah1.mismatch << endl;
		//}*/
		//------------------------

		//
		//
		//cout << center_p(b) << endl;
		//for (int i = 0; i < b.size(); i++)
		//{
		//	b[i].x = b[i].x*0.8;
		//	b[i].y = b[i].y*0.8;
		//}
		//cout << center_p(b) << endl;
		//double len = arcLength(b,true);
		//cout << "leng: " << leng
		//	<< endl << "len: " << len << endl;
		//vector<Point2f> b;
		//a.push_back(Point2f(100, 100));
		//a.push_back(Point2f(200, 100));
		//a.push_back(Point2f(200, 200));
		//a.push_back(Point2f(180, 200));
		//a.push_back(Point2f(100, 200));
		
		//a.push_back(Point2f(100, 200));
			
		//draw_polygen("111", a);


		//b.push_back(Point2f(300, 300));
		//b.push_back(Point2f(400, 300));
		//b.push_back(Point2f(400, 400));
		//b.push_back(Point2f(300, 400));
		//
		//draw_polygen("hahahah",b);
		//Point2f cen = center_p(b);
		//vector<Point2f> c;
		//Mat rot_mat = getRotationMatrix2D(cen, 45, 1);
		//transform(b, c, rot_mat);
		//		
		//vector<int> haha=tiling_opt->search_align_p(cen, b[0], c);
		//cout << haha[0]<<endl<<c[haha[0]] << endl;

		//warpAffine(img, img1, rot_mat, Size(800,800));
		//draw_polygen("111", a);
		//Moments mu = moments(a);
		//Point2f cen = Point2f(mu.m10 / mu.m00, mu.m01 / mu.m00);
		//Mat rot_mat(2, 3, CV_32FC1);
		////Point center = Point(300, 100);
		//double angle = 90.0;
		//double scale = 1;
		//cout << center_p(a) << endl;
		//rot_mat = getRotationMatrix2D(center_p(a), angle, scale);
		//cv::transform(a, b, rot_mat);
		
		//cout << center_p << endl;
		
		
		//
		//cout << cen << endl;
		//draw_polygen("222", a);
		//imshow("show: ", img1);
		//tiling_opt->points_dividing(imagename3);
		/*vector<Point2f> po;
		vector<Point2f> po1;
		po.push_back(Point2f(650, 650));
		po.push_back(Point2f(800, 650));
		po.push_back(Point2f(900, 750));
		po.push_back(Point2f(750, 750));
		po1.push_back(Point2f(650, 650));
		po1.push_back(Point2f(900, 750));
		po1.push_back(Point2f(750, 750));
		po1.push_back(Point2f(800, 650));
		
		vector<double>  ter(4, 0);
		vector<double>  ter1(4, 0);
		//tiling_opt->min_mismatch(po, po1, ter, ter1);*/
		//draw_polygen("win", po);
		//String imageName("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\22.png"); // by default
		//Mat src = imread(imageName, 0);
		//Mat src1 = src(Range(100,600), Range(0, 507));
		////threshold(src, src, 128, 1, cv::THRESH_BINARY);
		//cout << src.cols << "  " << src.rows << endl;
		//for (int i = 400; i <src.cols; i++)
		//	for (int j = 0; j < src.rows; j++)
		//	{
		//		src.at<uchar>(j, i)=(int)src.at<uchar>(j, i) + 1;
		//		//cout << (int)src.at<uchar>(i, j)<< endl;;
		//	}
		//
		////threshold(src, src, 1.5, 255, cv::THRESH_BINARY);
		////namedWindow("result", 1);
		////threshold(src, src, 128, 255, cv::THRESH_BINARY);
		//imshow("result", src1);	
		//imshow("win", src);

        //--------------------------------
        //测试compare_shapes函数
		//tiling_opt->load_dataset();
		//prototile_first->contourname = "test";
		//prototile_first->contour = prototile_first->readTxt();
		////cout << prototile_first->contour[5] << endl;
		//
		////prototile_first->contour_sam_cur();
		//vector<CandPat> candida_contours;
		//Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));;
		//draw_poly(tt, prototile_first->contour, Point2f(400, 400));//draw_polygen("hhhh", prototile_first->contour);
		//imshow("aaa", tt); 
		//candida_contours = tiling_opt->compare_shapes(prototile_first->contour, 2);*/
        ////-------------------------------------------


		//cout << order.first << "  " << order.second << endl;
		/*vector<Point2f> con1 = prototile_first->contour_sample[0];
		prototile_second->loadTileData(imagename2);
		vector<Point2f> con2 = prototile_second->contour_sample[0];
		prototile_third->loadTileData(imagename1);
		vector<Point2f> con3 = prototile_third->contour_sample[0];*/
		

		//prototile_first->loadTileData(imagename3);
		//vector<Point2f> con1 = prototile_first->contour_sample[0];
		//prototile_second->loadTileData(imagename2);
		//vector<Point2f> con2 = prototile_second->contour_sample[0];
		//prototile_third->loadTileData(imagename1);
		//vector<Point2f> con3 = prototile_third->contour_sample[0];

		
	}
	else if (f == 1)
	{
		//int i = 2;
		//while (i-- != 0)
		//{
			prototile_first->imgtocout("13");
			prototile_first->~Prototile();
			cout << prototile_first->contourname << endl;
			if (prototile_first->contour.empty())
			{
				cout << prototile_first->contour.size() << endl;
				cout << "empty" << endl;
			}
			else cout << "no empty" << prototile_first->contour.size() << endl;
		//}
		//prototile_first->txtpath = "D:\\VisualStudioProjects\\images\\txt\\";
		//prototile_first->dataroot = "D:\\VisualStudioProjects\\images\\scr\\";
		//prototile_first->imgtocout("13",1);

	}	
	else if (f==2){  //批量读图
		//int iii[] = {63};
		//for (int i = 1; i < 74; i++)
		//{
			//string image = int2string(i);
			//cout << image << endl;
			string image1 = "D:\\VisualStudioProjects\\p16.png";
			//string image1 = "D:\\VisualStudioProjects\\images\\scr\\new\\" + image+".png";

			Mat src = imread(image1, IMREAD_GRAYSCALE);


			threshold(src, src, 128, 255, cv::THRESH_BINARY);
			//imshow("???:", src);
			int rows = src.rows;
			int cols = src.cols;
			//////int count = 0;
			for (int i = 0; i <rows; i++)
				for (int j = 0; j < cols; j++)
				{
					//cout << (int)src.at<uchar>(i, j) << " ";
			//		//if (src.at<uchar>(i, j) == 0 || src.at<uchar>(i, j) == 255) count++;
					if (src.at<uchar>(i, j) == 0)
					{
			//			count++;
						src.at<uchar>(i, j) = 255;
					}
					else if (src.at<uchar>(i, j) == 255) src.at<uchar>(i, j) = 0;
				}
			imshow("!!!:", src);
			//////cout << count << endl;


		    //string image2 = int2string(i+534);
			//imwrite("D:\\VisualStudioProjects\\images\\new\\" + image2 + ".png",src);
			imwrite("D:\\VisualStudioProjects\\p16.png", src);
			//prototile_first->Pro_clear();
			//prototile_first->txtpath = "D:\\VisualStudioProjects\\images\\txt\\";
			//prototile_first->dataroot = "D:\\VisualStudioProjects\\images\\scr\\new\\";
			//////
			//////string image = int2string(i);
			//////cout << image << endl;
			//image = "3 (" + image + ")";
			//prototile_first->imgtocout(image);


			//prototile_first->contourname = image;
			//prototile_first->contour = prototile_first->readTxt();
			

			//if (prototile_first->contour.size() < 300)
			//cout << image + ".png may be error" << endl;
		//}
	
	}
	else if (f == 3)
	{
		tiling_opt->Tiling_clear();
		vector<Point2f> sim_mid = tiling_opt->simulation_mid(imagename3,12,0);



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
	//Tiling_tiles::Prototile *prototile_second;
	//prototile_second = new Tiling_tiles::Prototile();
	//prototile_first->loadTileData(imagename2);

	//Tiling_tiles::Tiling_opt *tiling;
	//tiling = new Tiling_tiles::Tiling_opt();
	//////tiling->com_cur_string(imagename1, imagename2);
	//////tiling->com_score(imagename1, imagename1);
	//tiling->com_score_manual(imagename1, imagename2);

	//检测各种仿射变换
	/*vector<Point2f> a;
	vector<Point2f> b;
	a.push_back(Point2f(1, 1));
	a.push_back(Point2f(2, 1));
	a.push_back(Point2f(2, 2));
	a.push_back(Point2f(3, 2));
	a.push_back(Point2f(3, 1));
	a.push_back(Point2f(4, 1));
	a.push_back(Point2f(5, 0));
	a.push_back(Point2f(6, 1));
	a.push_back(Point2f(7, 1));

	b.push_back(Point2f(1, 1));
	b.push_back(Point2f(2, 1));
	b.push_back(Point2f(3, 0));
	b.push_back(Point2f(4, 1));
	b.push_back(Point2f(5, 1));
	b.push_back(Point2f(5, 2));
	b.push_back(Point2f(6, 2));
	b.push_back(Point2f(6, 1));
	b.push_back(Point2f(7, 1));

	vector<vector<Point2f>> prototwo;
	tiling->Aff_place(a, b, prototwo);*/
	//Point2f s(5, 1);
	//Point2f e(1, 5);
	////double ab = tiling->re_warp_Aff(a, b, s, e);
	//double ab = tiling->warpAff_sca(a, b, s, e);
	//for (int i = 0; i < a.size(); i++)
	//{
	//	cout << "output: " << b[i] << endl;

	//}

	//warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)

	//tiling->DTW(a, b);

	
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
    finish = clock();
	cout << endl<< (finish - start) / CLOCKS_PER_SEC << " s " << endl;
	waitKey(0);
	return 0;
}


#include "tilingOpt.h"

using namespace Tiling_tiles;

int main(int argc, char** argv)
{

	string imagename1 = "15"; 
	string imagename2 = "19";
	//string imagename1 = "Boat";
	//string imagename2 = "fish8";
	//string imagename1 = "fish5";
	string imagename3 = "22";
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
	int f = 0;
	if (f == 0) //已有dataset
	{
		prototile_first->contourname = imagename3;
		prototile_first->contour = prototile_first->readTxt();
		prototile_first->contour_sam_cur();
		//Mat img = imread("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\22.png", 0);
		//Mat img1;
		double leng = prototile_first->c_length;
		vector<Point2f> a = prototile_first->contour;
		vector<Point2f> b = prototile_first->contour_sample[1];
		cout << center_p(b) << endl;
		for (int i = 0; i < b.size(); i++)
		{
			b[i].x = b[i].x*0.8;
			b[i].y = b[i].y*0.8;
		}
		cout << center_p(b) << endl;
		//double len = arcLength(b,true);
		//cout << "leng: " << leng
		//	<< endl << "len: " << len << endl;
		//vector<Point2f> b;
		//a.push_back(Point2f(100, 100));
		//a.push_back(Point2f(200, 100));
		//a.push_back(Point2f(200, 200));
		//a.push_back(Point2f(180, 200));
		//a.push_back(Point2f(100, 200));
		
		/*a.push_back(Point2f(100, 200));
			
		draw_polygen("111", a);


		b.push_back(Point2f(1, 1));
		b.push_back(Point2f(2, 1));
		b.push_back(Point2f(3, 0));
		b.push_back(Point2f(4, 1));
		b.push_back(Point2f(5, 1));
		b.push_back(Point2f(5, 2));
		b.push_back(Point2f(6, 2));

		*/
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


		/*tiling_opt->load_dataset();
		prototile_first->contourname = imagename3;
		prototile_first->contour = prototile_first->readTxt();
		prototile_first->contour_sam_cur();
		vector<Point2f> test1;
		vector<double> test11;
		for (int i = 0; i < 100; i++)
		{
			test1.push_back(prototile_first->contour_sample[0][i]);
			test11.push_back(prototile_first->contour_curva[0][i]);
		}
		vector<int> order = tiling_opt->compare_shapes(prototile_first->contour);

		prototile_second->loadPoints(tiling_opt->contour_dataset[order[1]]);
		vector<Point2f> test2;
		vector<double> test22;
		for (int i = 0; i < 90; i++)
		{
			test2.push_back(prototile_second->contour_sample[0][i]);
			test22.push_back(prototile_second->contour_curva[0][i]);
		}
		//cout << "test1:" << test11.size()
		//	<< "test2:" << test22.size() << endl;
		tiling_opt->min_mismatch(test1, test2, test11, test22);
		//cout << re << endl;

		*/
		/*vector<double> sss;
		vector<int> index_s;
		prototile_first->loadTileData(imagename3);
		vector<Point2f> contour_mid = prototile_first->contour_sample[0];
		int total_num = tiling_opt->contour_dataset.size();
		for (int can_num = 0; can_num < total_num; can_num++)
		{
			index_s.push_back(can_num);
			prototile_second->~Prototile();
			prototile_second->loadPoints(tiling_opt->contour_dataset[can_num]);
			vector<Point2f> contour_second = prototile_second->contour_sample[0];
			int method_ = 3;
			double score;
			score = matchShapes(contour_mid, contour_second, method_, 0);
			sss.push_back(score);
		}
		sort_comb(sss,index_s);
		for (int i = index_s.size()-1; i > 180; i--)
		{
			cout << index_s[i] << " " << sss[index_s[i]] << endl;
		}*/
		//cout << order.first << "  " << order.second << endl;
		/*vector<Point2f> con1 = prototile_first->contour_sample[0];
		prototile_second->loadTileData(imagename2);
		vector<Point2f> con2 = prototile_second->contour_sample[0];
		prototile_third->loadTileData(imagename1);
		vector<Point2f> con3 = prototile_third->contour_sample[0];*/
		
		/*double score[3][3];
		cout << endl << imagename3 << " v " << imagename2 <<": "<< endl;
		for (int i = 1; i < 4; i++)
		{
			score[0][i-1] = matchShapes(con1, con2, i, 0);
			cout << score[0][i - 1] << " " ;
		}
		cout << endl<< imagename1 << " v " << imagename2 << ": " << endl;
		for (int i = 1; i < 4; i++)
		{
			score[1][i - 1] = matchShapes(con3, con2, i, 0);
			cout << score[1][i - 1] << " ";
		}
		cout << endl << imagename3<<" v " << imagename1 << ": " << endl;
		for (int i = 1; i < 4; i++)
		{
			score[2][i - 1] = matchShapes(con1, con3, i, 0);
			cout << score[2][i - 1] << " ";
		}*/

		//prototile_first->loadTileData(imagename3);
		//vector<Point2f> con1 = prototile_first->contour_sample[0];
		//prototile_second->loadTileData(imagename2);
		//vector<Point2f> con2 = prototile_second->contour_sample[0];
		//prototile_third->loadTileData(imagename1);
		//vector<Point2f> con3 = prototile_third->contour_sample[0];

		
	}
	else if (f == 1)
	{
		prototile_first->imgtocout(imagename3);
	}	
	else if (f==2){  //批量读图
		for (int i = 201; i < 205; i++)
		{
			prototile_first->~Prototile();
			char ch[4];
			//for (int j = 0; j < 2; j++)
			//cout << ch[j] << endl;
			if (i < 100)
			{
				if (i / 10 == 0)
				{
					ch[0] = i % 10 + 48;
					ch[1] = '\0';
				}
				else
				{
					ch[0] = i / 10 + 48;
					ch[1] = i % 10 + 48;
					ch[2] = '\0';
				}
			}
			else
			{
				ch[0] = i / 100 + 48;
				ch[1] = ( i % 100 ) / 10+ 48;
				ch[2] = i % 10 + 48;
				ch[3] = '\0';

			}
			

			string image = ch;
			cout << image << endl;
			prototile_first->imgtocout(image);
			//if (prototile_first->contour.size() < 300)
			//cout << image + ".png may be error" << endl;
		}
	
	}
	else if (f == 3)
	{
		//Point2f xy;
		//if (get_line_intersection(Point2f(0, 0),Point2f( 10, 10),Point2f(0, 0), Point2f(6, 6), xy) == 0)
		//	cout << "no" << endl;
		//cout << xy.x << " " << xy.y<<endl;
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

	waitKey(0);
	getchar();
	return 0;
}


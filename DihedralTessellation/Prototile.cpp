#include "tilingOpt.h"

using namespace std;

namespace Tiling_tiles{
	Prototile::Prototile(){
		c_length = 0;
	};

	void Prototile::loadTileData(string tile_data)
	{
		//vector<double> cos_vec_contour;
		//vector<int> index_cos;
		contourname = tile_data;
		imgtocout();
	    readTxt(contour);
		//c_length = contour_length(contour);
		////采样并求曲率
		//contour_sam_cur();
		//cur_normalize();
	}

	void Prototile::imgtocout()
	{
		Mat src;
		Mat src_gray;
		int thresh = 100;
		int max_thresh = 255;

		//read image
		String imageName("D:/dataset/" + contourname + ".png"); // by default
		src = imread(imageName, IMREAD_COLOR);

		if (src.empty())
		{
			cerr << "No image supplied ..." << endl;
			return;
		}
		cvtColor(src, src_gray, COLOR_BGR2GRAY);
		blur(src_gray, src_gray, Size(3, 3));

		Mat canny_output;
		//考虑到可能有多个轮廓
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;

		//用candy法由灰度图求出掩码图
		Canny(src_gray, canny_output, thresh, thresh * 2, 3);

		//由掩码图求出有序轮廓点
		findContours(canny_output, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0));
		cout << "contours num:" << contours.size() << endl;
		Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);
		//Mat drawing2 = Mat::zeros(canny_output.size(), CV_8UC3);

		////output the oral points
		//Mat drawing3 = Mat::zeros(800,800, CV_8UC3);
		//for (int j = 0; j < contours[1].size(); j++)
		//{
		//	circle(drawing3, contours[1][j], 1, Scalar(255, 0, 0), -1);
		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//imshow("oral contour: ", drawing3);

		//output two 
		//namedWindow("contour one", WINDOW_AUTOSIZE);

		//分别处理多个轮廓，此时按两个算，第一个每两个点采样一个画出来
		int n = contours[0].size() / 2;
		Point rook_points[1][1000];
		for (int t = 0; t < n; t++)
		{
			rook_points[0][t] = contours[0][2 * t];
		}
		const Point* ppt[1] = { rook_points[0] };
		int npt[] = { n };
		polylines(drawing,
			ppt,
			npt,
			1,
			true,
			Scalar(255, 255, 255)
			);
		//imshow("contour" + contourname + " one", drawing);

		//存储数据时还是原样存储,现在是顺时针
		ofstream out("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (out.is_open())
		{
			out << contours[0].size()+1 << endl;//contours[0].size()
			for (int j = contours[0].size()-1; j >= 0; j--)
				out << contours[0][j].x << "," << contours[0][j].y <<endl;
			out << contours[0][contours[0].size() - 1].x << "," << contours[0][contours[0].size() - 1].y << endl;  //首尾连起来
		}
		cout << "contours[0].size(): " << contours[0].size()<< endl;
		out.close();

	}


	void Prototile::readTxt(vector<Point2f> &con_point)
	{
		//读取一个存有轮廓点的文件，格式对应上一步计算轮廓点保存的文件
		ifstream in("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (!in.is_open())
		{
			cout << "Error opening file";
			return;
		}
		//挨个处理每个字符
		cout << "Error opening file  ok!!!";
		vector<char> each_point;
		int aa = 0;
		int bb = 0;
		int nn = 0;
		char cc;
		char buf[200];
		//for ()
		in.getline(buf, 200);
		cout << "num: " << buf << endl;
		while (!in.eof())
		{
			aa = 0;
			bb = 0;
			nn = 0;
			int f = 0;
			in.getline(buf, 200);
			cc = buf[nn++];
			while ((cc >= '0' && cc <= '9')|| cc == ','|| cc == ' ')
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
		//while (in >> cc)
		//{
		//	if ((cc == '[') || (cc == ' '))
		//		continue;
		//	if ((cc >= '0' && cc <= '9'))
		//	{
		//		each_point.push_back(cc);
		//		continue;
		//	}
		//	if (cc == ',')
		//	{
		//		for (int i = 0; i < each_point.size(); i++)
		//		{
		//			aa = aa * 10 + (each_point[i] - 48);
		//		}
		//		//cout << "aa:  " << aa << endl;
		//		each_point.swap(vector<char>());
		//		continue;
		//	}

		//	if (cc == ']')
		//	{
		//		for (int i = 0; i < each_point.size(); i++)
		//		{
		//			bb = bb * 10 + (each_point[i] - 48);
		//		}
		//		//cout << "bb:  " << bb << endl;
		//		each_point.swap(vector<char>());
		//		con_point.push_back(Point(aa, bb));
		//		aa = 0;
		//		bb = 0;
		//	}
		//else continue;
		Mat drwa = Mat::zeros(800, 800, CV_8UC3);
		int i = 0;
		for (; i < con_point.size() / 4; i++)
		{
			circle(drwa, con_point[i], 1, Scalar(255, 0, 0), -1);
		}
		for (; i < con_point.size() / 2; i++)
		{
			circle(drwa, con_point[i], 1, Scalar(0, 255, 0), -1);
		}
		for (; i < con_point.size(); i++)
		{
			circle(drwa, con_point[i], 1, Scalar(0, 0, 255), -1);
		}
		imshow("contour" + contourname, drwa);
		in.close();

	}


	void Prototile::contour_sam_cur()
	{
		int iter_num = 0;
		double Lambda = 0;
		int sam_num = 0;
		// center point
		for (int j = 0; j < contour.size(); j++)
		{
			center_point.x += contour[j].x;
			center_point.y += contour[j].y;
		}
		center_point.x = center_point.x / contour.size();
		center_point.y = center_point.y / contour.size();

		//sampling and computing curvature
		for (int i = 1; i < 6; i++)
		{

			Lambda = 0;
			sam_num = i * 100;
			vector<Point2f> contour_sam;
			vector<Point2f> contour_sam_inver;
			Point2f sample;

			if (contour.size() < sam_num)
			{
				cout << "Lost the " << sam_num << " and more samples!" << endl;
				//break;
			}
			Lambda = c_length / sam_num;
			contour_sam.push_back(contour[0]);
			sample = contour[0];
			for (int t = 1; t < contour.size(); t++)
			{
				double length_ = length_two_point2f(sample, contour[t]);
				if (length_ > Lambda)
				{
					Point2f vec = unit_vec(contour[t] - sample);
					sample = sample + Lambda * vec;
					contour_sam.push_back(sample);
					t = t - 1;
				}
				else if (t < contour.size() - 1)
				{
					while ((length_ + length_two_point2f(contour[t], contour[t + 1])) < Lambda)
					{
						length_ = length_ + length_two_point2f(contour[t], contour[t + 1]);
						t++;
						if (t >= (contour.size() - 1)) break;
					}
					if (t >= (contour.size() - 1)) break;
					Point2f vec = unit_vec(contour[t + 1] - contour[t]);
					sample = contour[t] + (Lambda - length_) * vec;
					contour_sam.push_back(sample);
				}
			}
			if (contour_sam[0] == contour_sam[contour_sam.size() - 1]) contour_sam.pop_back();
			contour_sample.push_back(contour_sam);
			contour_curva.push_back(curvature_com_k(contour_sam));
			// invertion
			for (int i = contour_sam.size() - 1; i >= 0; i--)
			{
				contour_sam_inver.push_back(contour_sam[i]);
			}
			contour_sample_inver.push_back(contour_sam_inver);
			contour_curva_inver.push_back(curvature_com_k(contour_sam_inver));


			////_________________________show the result

			////Mat drawing4 = Mat::zeros(800, 800, CV_8UC3);
			//Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

			//for (int j = 0; j < contour_sam.size(); j++)
			//{
			//	circle(drawing4, contour_sam[j], 1, Scalar(0, 0, 0), -1);
			//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
			//}
			////circle(drawing4, contour_sam[contour_sam.size() - 2], 1, Scalar(255, 255, 0), -1);
			////circle(drawing4, contour_sam[contour_sam.size() - 1], 1, Scalar(255, 255, 0), -1);
			////namedWindow("simple contour", WINDOW_AUTOSIZE);

			/////*	int n = contour_sam.size();
			////	cout << "n: " << n << endl;
			////	Point rook_points[1][600];
			////	for (int t = 0; t < n; t++)
			////	{
			////		rook_points[0][t] = contour_sam[t];
			////	}
			////	const Point* ppt[1] = { rook_points[0] };
			////	int npt[] = { n };
			////	polylines(drawing4,
			////		ppt,
			////		npt,
			////		1,
			////		true,
			////		Scalar(255, 255, 255)
			////		);*/
			//string name = "the ";
			//name = name + char(i + 48) + " simple contour ";
			//imshow(name, drawing4);
			////________________________ show over
		}

		//求一点两个方向的曲率的均值
		for (int i = 0; i < 5; i++)
		{
			for (int j = 0; j < contour_curva[i].size(); j++)
			{
				contour_curva[i][j] = (contour_curva[i][j] + contour_curva_inver[i][contour_curva[i].size() - 1 - j]) / 2;
			}
			//cout << "cur: " << contour_curva[0][i] << "  cur_inver: " << contour_curva_inver[0][i] << endl;
		}
	}

	void Prototile::cur_normalize()
	{
		for (int n = 0; n < 5; n++)
		{
			double min = 100;
			double max = 0;
			for (int i = 0; i < contour_curva[n].size(); i++)
			{
				//cout << contour_curva[n][i] << " ";
				if (contour_curva[n][i]>100) continue;
				if (contour_curva[n][i] < min) min = contour_curva[n][i];
				if (contour_curva[n][i] > max) max = contour_curva[n][i];
			}

			///////////
			double step = (max - min) / 52;
			double mid = (max + min) / 2;
			//cout << "step: " << step << "\n max: " << max << "\n mid: " << mid << endl;

			//Z~A0a~z
			cout << endl;
			for (int i = 0; i < contour_curva[n].size(); i++)
			{
				if (contour_curva[n][i]>max) cur_string[n][i] = 'z';
				if (contour_curva[n][i] == mid) cur_string[n][i] = '0';
				else if (contour_curva[n][i] < mid)
				{
					int type = (mid - contour_curva[n][i]) / step - 1;
					char c = 'A' + type;
					cur_string[n][i] = c;
				}
				else if (contour_curva[n][i] > mid)
				{
					int type = (contour_curva[n][i] - mid) / step - 1;
					char c = 'a' + type;
					cur_string[n][i] = c;
				}
			}
			//cout << endl;
			//for (int i = 0; i < contour_curva[4].size(); i++)
			//{
			//	cout <<cur_string[4][i] ;
			//}
		}
		//  不需要再逆向求曲率字符串
		//for (int n = 0; n < 5; n++)
		//{
		//	double min = 1000;
		//	double max = 0;
		//	for (int i = 0; i < contour_curva_inver[n].size(); i++)
		//	{
		//		//cout << contour_curva_inver[n][i] << " ";
		//		if (contour_curva_inver[n][i] < min) min = contour_curva_inver[n][i];
		//		if (contour_curva_inver[n][i] > max) max = contour_curva_inver[n][i];
		//	}

		//	double step = (max - min) / 52;
		//	double mid = (max + min) / 2;
		//	cout << "step: " << step << "\n max: " << max << "\n mid: " << mid << endl;

		//	//Z~A0a~z
		//	cout << endl;
		//	for (int i = 0; i < contour_curva_inver[n].size(); i++)
		//	{
		//		if (contour_curva_inver[n][i] == mid) cur_string_inver[n][i] = '0';
		//		else if (contour_curva_inver[n][i] < mid)
		//		{
		//			int type = (mid - contour_curva_inver[n][i]) / step - 1;
		//			char c = 'A' + type;
		//			cur_string_inver[n][i] = c;
		//		}
		//		else if (contour_curva_inver[n][i] > mid)
		//		{
		//			int type = (contour_curva_inver[n][i] - mid) / step - 1;
		//			char c = 'a' + type;
		//			cur_string_inver[n][i] = c;
		//		}
		//	}
		//}


	}


}
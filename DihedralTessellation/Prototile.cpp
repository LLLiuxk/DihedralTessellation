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
		//imgtocout();
		contour = readTxt();
		c_length = contour_length(contour);
		//采样并求曲率
		contour_sam_cur();
		
		//cur_normalize();
	}

	void Prototile::loadPoints(vector<Point2f> con_point)
	{
		contour.swap(con_point);
		c_length = contour_length(contour);
		contour_sam_cur();
	}

	void Prototile::imgtocout(string tile_image)
	{
		Mat src;
		Mat src_gray;
		Mat src_2;
		int thresh = 100;
		int max_thresh = 255;

		//read image
		String imageName("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + tile_image + ".png"); // by default
		src = imread(imageName, IMREAD_COLOR);
		if (src.empty())
		{
			cerr << "No image supplied ..." << endl;
			return;
		}
		/*这里的处理过程：
		1.将任意一张图像转化成灰度图 
		2.模糊，利于下一步进行筛选过滤 
		3.转化二值图 
		4.模糊方便提取轮廓*/
		int raw = 0;
		if (raw == 1)
		{
		    cvtColor(src, src_gray, COLOR_BGR2GRAY);
			blur(src_gray, src_gray, Size(3, 3));
			//未经处理的图需要四步，之前处理好的用四步会处理过头
			threshold(src_gray, src_gray, 128, 255, cv::THRESH_BINARY);
			if (imread("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + tile_image + "2.png", IMREAD_COLOR).empty())
				imwrite("D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + tile_image + "2.png", src_gray);
			else cout << tile_image + "2.png has already exist!!!" << endl;
			blur(src_gray, src_gray, Size(3, 3));
			imshow("src_gray_blur", src_gray);
		}
		else
		{
			cvtColor(src, src_gray, COLOR_BGR2GRAY);
			blur(src_gray, src_gray, Size(3, 3));
			imshow("src_gray_blur", src_gray);
		}
		Mat canny_output;
		//考虑到可能有多个轮廓
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;

		//用candy法由灰度图求出掩码图
		Canny(src_gray, canny_output, thresh, thresh * 2, 3);
		imshow("canny_output", canny_output);
		//由掩码图求出有序轮廓点
		findContours(canny_output, contours, hierarchy, CV_RETR_EXTERNAL, CHAIN_APPROX_SIMPLE, Point(0, 0));
		cout << "contours num:" << contours.size() << endl;
		
		
		Mat drwa = Mat::zeros(800, 800, CV_8UC3);
		int i = 0;
		for (; i < contours[0].size() / 4; i++)
		{
			circle(drwa, contours[0][i], 1, Scalar(255, 0, 0), -1);
		}
		for (; i < contours[0].size() / 2; i++)
		{
			circle(drwa, contours[0][i], 1, Scalar(0, 255, 0), -1);
		}
		for (; i < contours[0].size(); i++)
		{
			circle(drwa, contours[0][i], 1, Scalar(0, 0, 255), -1);
		}
		imshow("contour" + tile_image, drwa);

		//output two 
		//namedWindow("contour one", WINDOW_AUTOSIZE);

		//分别处理多个轮廓，此时按两个算，第一个每两个点采样一个画出来
		/*Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);
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
			);*/
		//imshow("contour" + contourname + " one", drawing);

		//顺时针存储
		/*ofstream out("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (out.is_open())
		{
			out << contours[0].size()+1 << endl;//contours[0].size()
			for (int j = contours[0].size()-1; j >= 0; j--)
				out << contours[0][j].x << "," << contours[0][j].y <<endl;
			out << contours[0][contours[0].size() - 1].x << "," << contours[0][contours[0].size() - 1].y << endl;  //首尾连起来
		}
		cout << "contours[0].size(): " << contours[0].size()<< endl;
		out.close();*/

		//逆时针存储
		ofstream out("D:\\VisualStudioProjects\\contours\\" + tile_image + ".txt");
		if (out.is_open())
		{
			out << contours[0].size() + 1 << endl;//contours[0].size()
			for (int j = 0; j < contours[0].size(); j++)
				out << contours[0][j].x << "," << contours[0][j].y << endl;
			out << contours[0][0].x << "," << contours[0][0].y << endl;  //首尾连起来
		}
		cout << "contours[0].size(): " << contours[0].size() << endl;
		out.close();

	}


	vector<Point2f> Prototile::readTxt()
	{
		vector<Point2f> con_point;
		//读取一个存有轮廓点的文件，格式对应上一步计算轮廓点保存的文件
		ifstream in("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (!in.is_open())
		{
			cout << "Error opening file" << endl;
			return con_point;
		}
		//挨个处理每个字符
		cout << "Opening file!!!" << endl;
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
		/*Mat drwa = Mat::zeros(800, 800, CV_8UC3);
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
		imshow("contour" + contourname, drwa);*/
		in.close();
		return con_point;

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

		for (int i = 3; i < 7; i++)  //确定采样点数，此处为600点
		{			
			Lambda = 0;
			sam_num = i * 100;
			vector<Point2f> contour_sam;
			vector<Point2f> contour_sam_inver;
			Point2f sample;

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
			contour_curva.push_back(curvature_com(contour_sam));
			
			
			//// invertion
			//for (int i = contour_sam.size() - 1; i >= 0; i--)
			//{
			//	contour_sam_inver.push_back(contour_sam[i]);
			//}
			//contour_sample_inver.push_back(contour_sam_inver);
			//contour_curva_inver.push_back(curvature_com_k(contour_sam_inver));

			////求一点两个方向的曲率的均值
			//for (int j = 0; j < contour_curva[i-1].size(); j++)
			//{
			//	contour_curva[i-1][j] = (contour_curva[i-1][j] + contour_curva_inver[i-1][contour_curva[i-1].size() - 1 - j]) / 2;
			//}
			////cout << "cur: " << contour_curva[0][i] << "  cur_inver: " << contour_curva_inver[0][i] << endl;


			//_________________________show the result

			Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

			for (int j = 0; j < contour_sam.size(); j++)
			{
				circle(drawing4, contour_sam[j], 1, Scalar(0, 0, 0), -1);
				//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
			}
			string name = "the ";
			name = name + char(i + 48) + " sample contour ";
			imshow(name, drawing4);
			//________________________ show over
		}

		
	}

	vector<int> Prototile::convex_p(int max_cur_num)
	{
		//排序，找最大的max_cur_num个凹凸点
		contour.swap(vector<Point2f>());
		int sam_num = contour_sample.size();	
		contour = contour_sample[sam_num - 2];
		cconvex = contour_curva[sam_num - 2];
		cout << "contour_sample num: " << contour_sample[sam_num - 2].size() << endl;
		vector<int> index_num;
		int contoursize = contour.size();
		for (int i = 0; i < contoursize; i++)
		{
			index_num.push_back(i);

		}
		sort_cos(cconvex, index_num);

		vector<int> cand_points_index;
		int t = 1;
		cand_points_index.push_back(index_num[0]);
		cout << "length: " << c_length << endl;
		
		for (int i = 1; i < contoursize; i++)
		{
			if (t >= max_cur_num) break;
			else
			{ 
				int flag = 0;
				for (vector<int>::iterator it = cand_points_index.begin(); it != cand_points_index.end(); it++)
				{
					double leng = length_two_point2f(contour[index_num[i]], contour[*it]);				
					if (leng < 0.01*c_length)
					{
						flag = 1;
						break;
					}
				}

				if (flag == 0)
				{
					cand_points_index.push_back(index_num[i]);
					t++;
				}

			}
		}
		
		cout << "cand_points_index: " << cand_points_index.size()<<endl;
		//// show convex points
		//Mat drawing5 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

		//for (int j = 0; j < contour.size(); j++)
		//{
		//	circle(drawing5, contour[j], 1, Scalar(0, 0, 0), -1);

		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//for (int j = 0; j < 20; j++)
		//{
		//	circle(drawing5, contour[index_num[j]], 4, Scalar(0, 0, 255), -1);
		//}		
		//imshow("convex points: ", drawing5);
		
		return cand_points_index;
	}

	vector<int> Prototile::partition_points(string imaname)
	{
		
		vector<Point2f> ske_points;
		vector<Point2f> skeleton_points;
		ske_points = get_Skeleton(imaname, skeleton_points);
		cout << ske_points.size() << endl;

		int cur_p_num = 10;
		vector<int> max_order;
		imgtocout(imaname);
		loadTileData(imaname);
		max_order = convex_p(cur_p_num );

		int contoursize = contour.size();
		vector<int> cct;
		for (int i = 0; i < ske_points.size(); i++)
		{
			double dist = 1000;
			int index = 0;
			for (int j = 0; j < max_order.size(); j++)
			{
				if (dist > length_two_point2f(ske_points[i], contour[max_order[j]]))
				{
					dist = length_two_point2f(ske_points[i], contour[max_order[j]]);
					index = max_order[j];
				}
			}
			cct.push_back(index);
		}
		sort_bub(cct);
		for (int t = 0; t < cct.size() - 1; t++)
		{
			max_order.push_back((cct[t] + cct[t + 1]) / 2);
		}
		int mid = ((cct[0] + contour.size() - cct[cct.size() - 1]) / 2 + cct[cct.size() - 1]) % contour.size();
		max_order.push_back(mid);

		cct.swap(vector<int>());
		cct.push_back(max_order[0]);
		for (int i = 1; i < max_order.size(); i++)
		{
			int flag = 0;
			for (vector<int>::iterator it = cct.begin(); it != cct.end(); it++)
			{
				double leng = length_two_point2f(contour[max_order[i]], contour[*it]);
				if (leng < 0.01*c_length)
				{
					flag = 1;
					break;
				}
			}

			if (flag == 0)
			{
				cct.push_back(max_order[i]);
			}		
		}
		cct.swap(max_order);
		cout << max_order.size() << endl;
		sort_bub(max_order);
		for (int t = 0; t < max_order.size(); t++)
		{
			cout << max_order[t] << " ";
		}
		// show convex points
		Mat drawing5 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

		for (int j = 0; j < contoursize; j++)
		{
			circle(drawing5, contour[j], 1, Scalar(0, 0, 0), -1);

			//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		}
		for (int j = 0; j < max_order.size(); j++)
		{
			circle(drawing5, contour[max_order[j]], 4, Scalar(0, 0, 255), -1);
		}
		
		for (int j = 0; j < skeleton_points.size(); j++)
		{
			circle(drawing5, skeleton_points[j], 1, Scalar(128, 128, 128), -1);
		}
		for (int j = 0; j < ske_points.size(); j++)
		{
			circle(drawing5, ske_points[j], 4, Scalar(255, 0, 0), -1);
		}
		imshow("convex points: ", drawing5);

		return max_order;
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
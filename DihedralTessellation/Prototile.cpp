#include "tilingOpt.h"

using namespace std;

namespace Tiling_tiles{

	Prototile::Prototile(){
		c_length = 0;
		dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";
		txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";
	};

	Prototile::Prototile(string rname, string fpath){
		c_length = 0;
		dataroot = rname;
		txtpath = fpath;
	};
	
	void Prototile::Pro_clear()
	{
		contourname.clear();
		dataroot.clear();
		txtpath.clear();
		contour.swap(vector<Point2f>());
		cconvex.swap(vector<double>());
		contour_sample.swap(vector<vector<Point2f>>());
		contour_sample_flip.swap(vector<vector<Point2f>>());
		//contour_curva.swap(vector<vector<double>>());
		//contour_curva_flip.swap(vector<vector<double>>());
		c_length = 0;
		center_point = Point2f(0,0);
	}

	void Prototile::getpath()
	{
		if (dataroot.empty()) dataroot = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\";
		if (txtpath.empty()) txtpath = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\";

	}

	void Prototile::loadTileData(string tile_data)
	{
		getpath();
		//vector<double> cos_vec_contour;
		//vector<int> index_cos;		
		contourname = tile_data;
		//imgtocout();
		contour = readTxt();	
		//采样并求曲率
		contour_sam_cur();
		
		//cur_normalize();
	}

	void Prototile::loadPoints(vector<Point2f> con_point)
	{
		contour.swap(con_point);
		contour_sam_cur();
	}

	void Prototile::imgtocout(string tile_image, int raw)
	{
		getpath();
		Mat src;
		Mat src_gray;
		Mat src_2;
		int thresh = 100;
		int max_thresh = 255;

		//read image
		String imageName = dataroot + tile_image + ".png"; // by default
		cout << imageName << endl;
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
		if (raw == 1)
		{
		    cvtColor(src, src_gray, COLOR_BGR2GRAY);
			//blur(src_gray, src_gray, Size(3, 3));
			GaussianBlur(src_gray, src_gray, Size(3,3), 10, 10);
			//未经处理的图需要四步，之前处理好的用四步会处理过头
			threshold(src_gray, src_gray, 128, 255, cv::THRESH_BINARY); //255 white
			//imwrite(dataroot +"new\\"+ tile_image + ".png", src_gray);
			imwrite(dataroot + tile_image + ".png", src_gray);

			/*if (imread(dataroot + tile_image + "2.png", IMREAD_COLOR).empty())
				imwrite(dataroot + tile_image + "2.png", src_gray);
			else
			{
				cout << tile_image + "2.png has already exist!!!" << endl;
				imwrite(dataroot + tile_image + "0000.png", src_gray);
			}*/
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
		string filepath = txtpath + tile_image + ".txt";
		fileout(filepath, contours[0]);

	}


	vector<Point2f> Prototile::readTxt()
	{
		vector<Point2f> con_point;
		//读取一个存有轮廓点的文件，格式对应上一步计算轮廓点保存的文件
		string filepath = txtpath + contourname + ".txt";
		ifstream in(filepath);
		if (!in.is_open())
		{
			cout << filepath << endl;
			cout << "Error opening file" << endl;
			return con_point;
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
		double Lambda = 0;
		int sam_num = 0;
		c_length = contour_length(contour);
		// center point
		center_point = center_p(contour);

		//sampling and computing curvature
		for (int i = 1; i < 6; i++)  //确定采样点数，此处为500点
		{			
			vector<Point2f> contour_sam;
			vector<Point2f> contour_sam_flip;
			contour_sam = sampling(contour,i); //点数为 i*100
			
			contour_sample.push_back(contour_sam);
			//contour_curva.push_back(curvature_com(contour_sam));

			contour_sam_flip = flip_contour(contour_sam,0); //0是水平翻转
			contour_sample_flip.push_back(contour_sam_flip);
			//contour_curva_flip.push_back(curvature_com(contour_sam_flip));
			//_________________________show the result

			//Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

			//for (int j = 0; j < contour_sam.size(); j++)
			//{
			//	circle(drawing4, contour_sam[j], 1, Scalar(0, 0, 0), -1);
			//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
			//}
			//string name = "the ";
			//name = name + char(i + 48) + " sample contour ";
			//imshow(name, drawing4);
			//________________________ show over
		}		
	}

	vector<Point2f> Prototile::flip_contour(vector<Point2f> cont_s, int flag)
	{
		Point2f ccen = center_p(cont_s);
		int cont_size = cont_s.size();

		//flag==0 水平翻转
		if (flag == 0)
		{
			for (int i = 0; i < cont_size; i++)
			{
				cont_s[i].x = 2 * ccen.x - cont_s[i].x;
			}
			for (int i = 0; i < cont_size / 2; i++)
			{
				Point2f mid = cont_s[i];
				cont_s[i] = cont_s[cont_size - 1 - i];
				cont_s[cont_size - 1 - i] = mid;
			}
		}
		else if (flag == 1)
		{
			for (int i = 0; i < cont_size; i++)
			{
				cont_s[i].y = 2 * ccen.y - cont_s[i].y;
			}
			for (int i = 0; i < cont_size / 2; i++)
			{
				Point2f mid = cont_s[i];
				cont_s[i] = cont_s[cont_size - 1 - i];
				cont_s[cont_size - 1 - i] = mid;
			}
		}
		return cont_s;
	}

	vector<vector<double>> Prototile::compute_TAR(vector<Point2f> &contour_)
	{
		vector<vector<double>> all_tar;
		int consize = contour_.size();
		int tar_num = consize / 2 - 1;
		cout << "consize: " << consize << " tar_num: " << tar_num << endl;
		vector<double> maxtar(tar_num, 0);
		for (int i = 0; i < consize; i++)
		{
			vector<double> one_p_tar;
			for (int j = 0; j < tar_num; j++)
			{
				Point2f vpsubts_vp = contour_[(i - j - 1 + consize) % consize] - contour_[i];
				Point2f vpplusts_vp = contour_[(i + j + 1) % consize] - contour_[i];
				double tar = 0.5 * tar_sin_2vector(vpplusts_vp,vpsubts_vp);
				//cout << vpsubts_vp << "  " << vpplusts_vp << endl;
				one_p_tar.push_back(tar);
				if (abs(tar) > maxtar[j]) maxtar[j] = abs(tar);
			}
			all_tar.push_back(one_p_tar);
		}
		for (int i = 0; i < consize; i++)
		{
			for (int j = 0; j < tar_num; j++)
			{
				all_tar[i][j] = all_tar[i][j] / maxtar[j];
			}
		}
		cout << all_tar[0].size() << endl;
		return all_tar;
	}

	vector<int> Prototile::convex_p(int max_cur_num)
	{
		//排序，找最大的max_cur_num个凹凸点
		contour.swap(vector<Point2f>());
		int sam_num = contour_sample.size();	
		contour = contour_sample[sam_num - 1];
		cconvex = curvature_com(contour);
		cout << "contour_sample num: " << contour_sample[sam_num - 1].size() << endl;
		vector<int> cand_points_index;
		//cand_points_index = most_convex_p(contour, cconvex, max_cur_num);
		cand_points_index = feature_points(contour, 1, 3, cos(PI * 160 / 180));
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
		int cur_p_num = 20;   //cur_p_num 个不相邻的最大cos值点
		int margin = 12;      //margin个点的采样间隔
		double ratio = 0.012; //筛选间隔与周长之比
		vector<int> max_order;
		imgtocout(imaname);
		loadTileData(imaname);
		max_order = convex_p(cur_p_num );
		int contoursize = contour.size();
		vector<int> all_order = max_order;

		for (int i = 0; i < contoursize; i = i + margin)
		{
			int flag = 0;
			for (vector<int>::iterator it = max_order.begin(); it != max_order.end(); it++)
			{
				double leng = length_two_point2f(contour[i], contour[*it]);
				if (leng < ratio*c_length)
				{
					flag = 1;
					break;
				}
			}

			if (flag == 0)
			{
				all_order.push_back(i);
			}
		}		
		cout << "all_order.size:"<<all_order.size() << endl;
		sort_bub(all_order);
		/*for (int t = 0; t < all_order.size(); t++)
		{
			cout << all_order[t] << " ";
		}*/
		// show convex points
		Mat drawing5 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

		for (int j = 0; j < contoursize; j++)
		{
			circle(drawing5, contour[j], 1, Scalar(0, 0, 0), -1);

			//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		}
		for (int j = 0; j < max_order.size(); j++)
		{
			circle(drawing5, contour[max_order[j]], 6, Scalar(0, 0, 255), -1);
		}
		
		for (int j = 0; j < all_order.size(); j++)
		{
			circle(drawing5, contour[all_order[j]], 3, Scalar(128, 128, 128), -1);
		}
		
		imshow("convex points: ", drawing5);

		return all_order;
	}

	//void Prototile::cur_normalize()
	//{
	//	for (int n = 0; n < 5; n++)
	//	{
	//		double min = 100;
	//		double max = 0;
	//		for (int i = 0; i < contour_curva[n].size(); i++)
	//		{
	//			//cout << contour_curva[n][i] << " ";
	//			if (contour_curva[n][i]>100) continue;
	//			if (contour_curva[n][i] < min) min = contour_curva[n][i];
	//			if (contour_curva[n][i] > max) max = contour_curva[n][i];
	//		}

	//		///////////
	//		double step = (max - min) / 52;
	//		double mid = (max + min) / 2;
	//		//cout << "step: " << step << "\n max: " << max << "\n mid: " << mid << endl;

	//		//Z~A0a~z
	//		cout << endl;
	//		for (int i = 0; i < contour_curva[n].size(); i++)
	//		{
	//			if (contour_curva[n][i]>max) cur_string[n][i] = 'z';
	//			if (contour_curva[n][i] == mid) cur_string[n][i] = '0';
	//			else if (contour_curva[n][i] < mid)
	//			{
	//				int type = (mid - contour_curva[n][i]) / step - 1;
	//				char c = 'A' + type;
	//				cur_string[n][i] = c;
	//			}
	//			else if (contour_curva[n][i] > mid)
	//			{
	//				int type = (contour_curva[n][i] - mid) / step - 1;
	//				char c = 'a' + type;
	//				cur_string[n][i] = c;
	//			}
	//		}
	//		//cout << endl;
	//		//for (int i = 0; i < contour_curva[4].size(); i++)
	//		//{
	//		//	cout <<cur_string[4][i] ;
	//		//}
	//	}
	//}


}
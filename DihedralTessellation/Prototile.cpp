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
	    readTxt(contour);
		c_length = contour_length(contour);
		//������������
		contour_sam_cur();
		

		//cur_normalize();
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
		/*����Ĵ�����̣�
		1.������һ��ͼ��ת���ɻҶ�ͼ 
		2.ģ����������һ������ɸѡ���� 
		3.ת����ֵͼ 
		4.ģ��������ȡ����*/
		int raw = 0;
		if (raw == 1)
		{
		    cvtColor(src, src_gray, COLOR_BGR2GRAY);
			blur(src_gray, src_gray, Size(3, 3));
			//δ�������ͼ��Ҫ�Ĳ���֮ǰ����õ����Ĳ��ᴦ���ͷ
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
		//���ǵ������ж������
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;

		//��candy���ɻҶ�ͼ�������ͼ
		Canny(src_gray, canny_output, thresh, thresh * 2, 3);
		imshow("canny_output", canny_output);
		//������ͼ�������������
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

		//�ֱ�������������ʱ�������㣬��һ��ÿ���������һ��������
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

		//˳ʱ��洢
		/*ofstream out("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (out.is_open())
		{
			out << contours[0].size()+1 << endl;//contours[0].size()
			for (int j = contours[0].size()-1; j >= 0; j--)
				out << contours[0][j].x << "," << contours[0][j].y <<endl;
			out << contours[0][contours[0].size() - 1].x << "," << contours[0][contours[0].size() - 1].y << endl;  //��β������
		}
		cout << "contours[0].size(): " << contours[0].size()<< endl;
		out.close();*/

		//��ʱ��洢
		ofstream out("D:\\VisualStudioProjects\\contours\\" + tile_image + ".txt");
		if (out.is_open())
		{
			out << contours[0].size() + 1 << endl;//contours[0].size()
			for (int j = 0; j < contours[0].size(); j++)
				out << contours[0][j].x << "," << contours[0][j].y << endl;
			out << contours[0][0].x << "," << contours[0][0].y << endl;  //��β������
		}
		cout << "contours[0].size(): " << contours[0].size() << endl;
		out.close();

	}


	void Prototile::readTxt(vector<Point2f> &con_point)
	{
		//��ȡһ��������������ļ�����ʽ��Ӧ��һ�����������㱣����ļ�
		ifstream in("D:\\VisualStudioProjects\\contours\\" + contourname + ".txt");
		if (!in.is_open())
		{
			cout << "Error opening file" << endl;
			return;
		}
		//��������ÿ���ַ�
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
		for (int i = 1; i < 5; i++)  //ȷ�������������˴�Ϊ800��
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
			contour_curva.push_back(curvature_com_k(contour_sam));
			
			
			// invertion
			for (int i = contour_sam.size() - 1; i >= 0; i--)
			{
				contour_sam_inver.push_back(contour_sam[i]);
			}
			contour_sample_inver.push_back(contour_sam_inver);
			contour_curva_inver.push_back(curvature_com_k(contour_sam_inver));

			//��һ��������������ʵľ�ֵ
			for (int j = 0; j < contour_curva[i-1].size(); j++)
			{
				contour_curva[i-1][j] = (contour_curva[i-1][j] + contour_curva_inver[i-1][contour_curva[i-1].size() - 1 - j]) / 2;
			}
			//cout << "cur: " << contour_curva[0][i] << "  cur_inver: " << contour_curva_inver[0][i] << endl;


			//_________________________show the result

			//Mat drawing4 = Mat::zeros(800, 800, CV_8UC3);
			Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
			double max_cur = 0.0;
			for (int j = 0; j < contour_sam.size(); j++)
			{
				if (contour_curva[i-1][j]>0.1)
				{
					max_cur = max_cur + 1;
					circle(drawing4, contour_sam[j], 1, Scalar(0, 0, 255), -1);
				}
				else 
					circle(drawing4, contour_sam[j], 1, Scalar(0, 0, 0), -1);

				//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
			}
			cout << "contour_curva["<<i-1<<"].max= " << max_cur << endl;
			//circle(drawing4, contour_sam[contour_sam.size() - 2], 1, Scalar(255, 255, 0), -1);
			//circle(drawing4, contour_sam[contour_sam.size() - 1], 1, Scalar(255, 255, 0), -1);
			//namedWindow("simple contour", WINDOW_AUTOSIZE);

			///*	int n = contour_sam.size();
			//	cout << "n: " << n << endl;
			//	Point rook_points[1][600];
			//	for (int t = 0; t < n; t++)
			//	{
			//		rook_points[0][t] = contour_sam[t];
			//	}
			//	const Point* ppt[1] = { rook_points[0] };
			//	int npt[] = { n };
			//	polylines(drawing4,
			//		ppt,
			//		npt,
			//		1,
			//		true,
			//		Scalar(255, 255, 255)
			//		);*/
			string name = "the ";
			name = name + char(i + 48) + " simple contour ";
			imshow(name, drawing4);
			//________________________ show over
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
		//  ����Ҫ�������������ַ���
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
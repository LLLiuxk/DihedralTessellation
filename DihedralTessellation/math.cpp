#include<stdio.h>
#include<math.h>
#include "tilingOpt.h"

namespace Tiling_tiles{

	void MyLine(Mat img, Point2f start, Point2f end, string color1)
	{
		Scalar color;
		if (color1.compare("red") == 0) color = Scalar(0, 0, 255);
		else if (color1.compare("blue") == 0) color = Scalar(255, 0, 0);
		else if (color1.compare("green") == 0) color = Scalar(0, 255, 0);
		else if (color1.compare("cyan") == 0) color = Scalar(255, 255, 0);
		else if (color1.compare("grey") == 0) color = Scalar(190, 190, 190);
		else if (color1.compare("yellow") == 0) color = Scalar(0, 255, 255);
		else if (color1.compare("purple") == 0) color = Scalar(60, 32, 240);
		else if (color1.compare("white") == 0) color = Scalar(255, 255, 255);
		else if (color1.compare("black") == 0) color = Scalar(0, 0, 0);
		else color = Scalar(0, 165, 255);
		int thickness = 2;
		int lineType = 8;
		line(img,
			start,
			end,
			color,
			thickness,
			lineType);
	}
	Mat draw_polygen(string win_name, vector<Point2f> contour_s)
	{
		int col = 800;
		int row = 800;
		Mat drawing_pro = Mat(row, col, CV_8UC3, Scalar(255, 255, 255));
		int n = contour_s.size();
		//cout << "n: " << n << endl;
		Point rook_points[1][2000];
		for (int t = 0; t < n; t++)
		{
			rook_points[0][t] = contour_s[t];
		}
		const Point* ppt[1] = { rook_points[0] };
		int npt[] = { n };
		fillPoly(drawing_pro,
			ppt,
			npt,
			1,
			Scalar(0, 0, 0) //黑色
			//Scalar(255, 255, 255) //白色
			);
		//imshow(win_name, drawing_pro);
		return drawing_pro;
	}
	void draw_poly(Mat &drawing_, vector<Point2f> contour_s, Point2f center)
	{
		Point2f shift = center - center_p(contour_s);
		int n = contour_s.size();
		//cout << "n: " << n << endl;
		Point rook_points[1][2000];
		for (int t = 0; t < n; t++)
		{
			rook_points[0][t] = contour_s[t] + shift;
		}
		const Point* ppt[1] = { rook_points[0] };
		int npt[] = { n };
		fillPoly(drawing_,
			ppt,
			npt,
			1,
			Scalar(0, 0, 0) //黑色
			//Scalar(255, 255, 255) //白色
			);
		//circle(drawing_, contour_s[0] + shift, 4, Scalar(255), 3);
		//circle(drawing_, center, 4, Scalar(0, 255, 255), -1);
	}

	void draw_allplane(Mat &drawing_, vector<Point2f> contour_, vector<int> vec_, int type, double scale)
	{
		int drawrow = drawing_.rows;
		int drawcol = drawing_.cols;
		int border = 300 * scale;
		int con_size = contour_.size();
		if (scale != 1)
		{
			for (int i = 0; i < con_size; i++)
			{
				contour_[i] = scale*contour_[i];
			}
		}
		Point2f center = center_p(contour_);
		vector<Point2f> allcent_in_plane;
		//vector<vector<Point2f>> four_place;
		//vector<Point2f> one_loca;
		if (type == 0) //translation
		{
			stack<Point2f> cent_stack;
			Point2f step[6];
			step[0] = contour_[vec_[2]] - contour_[vec_[0]];
			step[1] = contour_[vec_[3]] - contour_[vec_[1]];
			step[2] = step[0] + step[1];
			step[3] = -step[0];
			step[4] = -step[1];
			step[5] = -step[2];
			Point2f cenp = Point2f(drawrow / 2, drawcol/2);
			draw_poly(drawing_, contour_, cenp);
			cent_stack.push(cenp);
			while (!cent_stack.empty())
			{
				Point2f top_p = cent_stack.top();
				cent_stack.pop();
				for (int i = 0; i < 6; i++)
				{
					cenp = top_p + step[i];
					if ((cenp.x<drawcol + border) && (cenp.x>-border) && (cenp.y<drawrow + border) && (cenp.y>-border))
					{
						int flag = 0;
						for (vector<Point2f>::iterator it = allcent_in_plane.begin(); it != allcent_in_plane.end(); it++)
						{
							double leng = length_two_point2f(cenp, *it);
							if (leng < 10)
							{
								flag = 1;
								break;
							}
						}
						if (flag == 0)
						{
							cent_stack.push(cenp);
							allcent_in_plane.push_back(cenp);
							draw_poly(drawing_, contour_, cenp);
						}

						
					}
				}
			}
			
		}
		else if (type == 1) //rotation
		{
			cout << "no";
		}
	}

	Point2f center_p(vector<Point2f> contour_)
	{
		//利用轮廓的矩
		Moments mu = moments(contour_);
		return Point2f(mu.m10 / mu.m00, mu.m01 / mu.m00);
	}

	double contour_length(vector<Point2f> contour)
	{
		double length = 0;
		int i = 0;
		for (; i < contour.size() - 1; i++)
		{
			length += length_two_point2f(contour[i], contour[i + 1]);
		}
		length += length_two_point2f(contour[i], contour[0]);
		return length;
	}

	double length_two_point2f(Point2f &u, Point2f &v)
	{
		return sqrt((u.x - v.x)*(u.x - v.x) + (u.y - v.y)*(u.y - v.y));
	}

	//double area_poly(vector<Point2f> &cont)
	//{
	//	int csize = cont.size();
	//	if (csize < 3) return 0;
	//	double sum = 0;
	//	for (int i = 0; i < csize; i++)
	//		sum += cont[i].x * cont[(i + 1) % csize].y - cont[i].y * cont[(i + 1) % csize].x;
	//	return fabs(sum / 2.0);
	//}
	

	void sort_comb(vector<double> vect, vector<int> &index_num) //将下标和数值联合排序，只保留下标的排序,从大到小
	{
		int i, j;
		double temp;
		int num;
		for (i = 0; i < vect.size() - 1; i++)
			for (j = 0; j < vect.size() - 1 - i; j++)
				if (vect[j] < vect[j + 1])
				{
					temp = vect[j];
					vect[j] = vect[j + 1];
					vect[j + 1] = temp;
					num = index_num[j];
					index_num[j] = index_num[j + 1];
					index_num[j + 1] = num;
				}
	}
	int line_intersection(Line_Seg line1, Line_Seg line2, Point2f &cross_p)
	{
		Point2f start1 = line1.start;
		Point2f end1 = line1.end;
		Point2f start2 = line2.start;
		Point2f end2 = line2.end;
		Point2f s10 = end1 - start1;
		Point2f s32 = end2 - start2;
		Point2f s02 = start1 - start2;
		float s_numer, t_numer, denom, t;

		denom = s10.x * s32.y - s32.x * s10.y;
		s_numer = s10.x * s02.y - s10.y * s02.x;
		t_numer = s32.x * s02.y - s32.y * s02.x;

		if (denom == 0)//平行或共线
		{
			if (s_numer == 0)//Collinear,返回离end1最近的点
			{
				double dis1 = sqrt((start2.x - end1.x)*(start2.x - end1.x) + (start2.y - end1.y)*(start2.y - end1.y));
				double dis2 = sqrt((end2.x - end1.x)*(end2.x - end1.x) + (end2.y - end1.y)*(end2.y - end1.y));
				if (dis1 > dis2)
				{
					cross_p = end2;
				}
				else{
					cross_p = start2;
				}
				return 2;
			}
			else return 0; // parallel
			cout << "denom == 0" << endl;
		}
		bool denomPositive = denom > 0;

		if ((s_numer < 0) == denomPositive)//参数是大于等于0且小于等于1的，分子分母必须同号且分子小于等于分母
			return 0; // No collision


		if ((t_numer < 0) == denomPositive)
			return 0; // No collision

		if (fabs(s_numer) > fabs(denom) || fabs(t_numer) > fabs(denom))
			return 0; // No collision
		// Collision detected
		t = t_numer / denom;

		//cout << "t:" << t << endl;
		cross_p.x = start1.x + t * s10.x;
		cross_p.y = start1.y + t * s10.y;
		return 1;
	}

	bool seif_intersect(vector<Point2f> &contour_)
	{
		int sizec = contour_.size();
		for (int i = 0; i < sizec-2; i++)
		{
			if (i == 0)
			{
				for (int j = 2; j < sizec-1; j++)
				{
					Point2f crossp;
					Line_Seg line1(contour_[i], contour_[i + 1]);
					Line_Seg line2(contour_[j], contour_[j + 1]);
					if (line_intersection(line1, line2, crossp)==1)
					{
						//cout << i << "  " << j << endl;
						return true;
					}
				}
			}
			else{
				for (int j = i + 2; j < sizec; j++)
				{
					Point2f crossp;
					Line_Seg line1(contour_[i], contour_[i + 1]);
					Line_Seg line2(contour_[j], contour_[(j + 1) % sizec]);
					if (line_intersection(line1, line2, crossp)==1)
					{
						//cout << i << "  " << j << endl;
						return true;
					}
				}		
			}	
		}
		return false;
	}

	Point2f unit_vec(Point2f vec)
	{
		double fenmu = sqrt(vec.x*vec.x + vec.y*vec.y);
		Point2f unit = Point2f(vec.x / fenmu, vec.y / fenmu);
		return unit;
	}

	double cos_two_vector(Point2f &v0, Point2f &v1)
	{
		return unit_vec(v0).x*unit_vec(v1).x + unit_vec(v0).y*unit_vec(v1).y;
	}

	double sin_2vector_convexc(Point2f &v0, Point2f &v1)
	{
		return unit_vec(v0).x*unit_vec(v1).y - unit_vec(v0).y*unit_vec(v1).x;
	}

	double cur_length_two_p(double cur1, double cur2)
	{
		double mis = 0;
		if ((cur1 < 0 && cur2 < 0) || (cur1>0 && cur2>0))
		{
			mis = (cur1 - cur2)*(cur1 - cur2);
		}
		else 
		{
			double factor = cur1 + cur2 + 2;//
			if (factor > 2.45) factor = 2.45;
			mis = factor*factor;
		}		
		return mis;
	}

	int cur_char_length(char a, char b)
	{
		if (a > 'A'&&a < 'Z')
		{
			if (b == '0')
				return a - 64;
			if (b > 'A'&& b < 'Z')
				return abs(a - b);
			if (b > 'a'&& b < 'z')
				return a + b - 160;
		}
		if (a > 'a'&& a < 'z')
		{
			if (b == '0')
				return a - 96;
			if (b > 'A'&& b < 'Z')
				return a + b - 160;
			if (b > 'a'&& b < 'z')
				return abs(a - b);
		}
		if (a == '0')
		{
			if (b == '0')
				return 0;
			if (b > 'A'&& b < 'Z')
				return b - 64;
			if (b > 'a'&& b < 'z')
				return abs(a - b);
		}

	}

	void fileout(string filepath, vector<Point> contour_)
	{
		ofstream out(filepath);
		if (out.is_open())
		{
			out << contour_.size() + 1 << endl;//contours[0].size()
			for (int j = 0; j < contour_.size(); j++)
				out << contour_[j].x << "," << contour_[j].y << endl;
			out << contour_[0].x << "," << contour_[0].y << endl;  //首尾连起来
		}
		cout << "contours[0].size(): " << contour_.size() << endl;
		out.close();
	}

	void bbx_center_point(vector<vector<Point2f>> all_point, vector<Point2f> &five_p)
	{
		five_p.swap(vector<Point2f>());
		vector<Point2f> contour;
		for (int i = 0; i < all_point.size(); i++)
		{
			for (int j = 0; j < all_point[i].size(); j++)
			{
				contour.push_back(all_point[i][j]);
			}
		}
		double bbx_max_x = -10000;
		double bbx_max_y = -10000;
		double bbx_min_x = 10000;
		double bbx_min_y = 10000;
		double center_x = 0;
		double center_y = 0;

		for (int i = 0; i < contour.size(); i++)
		{
			center_x += contour[i].x;
			center_y += contour[i].y;
			if (contour[i].x < bbx_min_x) bbx_min_x = contour[i].x;
			if (contour[i].x > bbx_max_x) bbx_max_x = contour[i].x;
			if (contour[i].y < bbx_min_y) bbx_min_y = contour[i].y;
			if (contour[i].y > bbx_max_y) bbx_max_y = contour[i].y;
		}
		center_x = center_x / contour.size();
		center_y = center_y / contour.size();
		five_p.push_back(Point2f(center_x, center_y));
		five_p.push_back(Point2f(bbx_min_x, bbx_max_y));
		five_p.push_back(Point2f(bbx_min_x, bbx_min_y));
		five_p.push_back(Point2f(bbx_max_x, bbx_min_y));
		five_p.push_back(Point2f(bbx_max_x, bbx_max_y));

	}

	vector<Point2f> b_box(vector<Point2f> contour)
	{
		vector<Point2f> four_cor;
		double bbx_max_x = -10000;
		double bbx_max_y = -10000;
		double bbx_min_x = 10000;
		double bbx_min_y = 10000;
		for (int i = 0; i < contour.size(); i++)
		{
			if (contour[i].x < bbx_min_x) bbx_min_x = contour[i].x;
			if (contour[i].x > bbx_max_x) bbx_max_x = contour[i].x;
			if (contour[i].y < bbx_min_y) bbx_min_y = contour[i].y;
			if (contour[i].y > bbx_max_y) bbx_max_y = contour[i].y;
		}
		four_cor.push_back(Point2f(bbx_min_x, bbx_max_y));
		four_cor.push_back(Point2f(bbx_min_x, bbx_min_y));
		four_cor.push_back(Point2f(bbx_max_x, bbx_min_y));
		four_cor.push_back(Point2f(bbx_max_x, bbx_max_y));
		return four_cor;
	}


	vector<double> curvature_com_k(vector<Point2f> &contour_sam) //k等于角度的差值比弧长的差值（近似为三角的斜边）
	{
		vector<double> eachOfcurvature;
		int c_s = contour_sam.size();
		for (int i = 0; i < c_s; i++)
		{
			double angle1_tan = (contour_sam[(i + 1) % c_s].y - contour_sam[i].y) / (contour_sam[(i + 1) % c_s].x - contour_sam[i].x);
			double angle2_tan = (contour_sam[(i + 2) % c_s].y - contour_sam[(i + 1) % c_s].y) / (contour_sam[(i + 2) % c_s].x - contour_sam[(i + 1) % c_s].x);
			double angle_incre = abs(atan(angle2_tan) - atan(angle1_tan));
			double curvature = angle_incre / sqrt((contour_sam[(i + 1) % c_s].x - contour_sam[i].x)*(contour_sam[(i + 1) % c_s].x - contour_sam[i].x) + (contour_sam[(i + 1) % c_s].y - contour_sam[i].y)*(contour_sam[(i + 1) % c_s].y - contour_sam[i].y));
			eachOfcurvature.push_back(curvature);
		}
		return eachOfcurvature;
	}

	vector<double> curvature_com(const vector<Point2f> &contour_sam)
	{
		vector<double> eachOfcurvature;
		int c_s = contour_sam.size();
		//sin_2vector_convexc>0 为一凸点,cos值+1.1; <0为一凹点,cos值-1.1
		//使用1.1是为了防止在无限接近0时来回加减出现近似误差
		//此处需要注意，opencv图像与正常坐标y轴相反
		//提取轮廓点在图上看是逆时针，实际在正常坐标系为顺时针，因此应用顺时针来计算凹凸
		for (int i = 0; i < c_s; i++)
		{
			double curvature = cos_two_vector(contour_sam[(i + c_s - 1) % c_s] - contour_sam[i], contour_sam[(i + 1) % c_s] - contour_sam[i]);
			if (sin_2vector_convexc(contour_sam[(i + c_s - 1) % c_s] - contour_sam[i], contour_sam[(i + 1) % c_s] - contour_sam[i]) > 0)
				eachOfcurvature.push_back(curvature + 1.1);
			else eachOfcurvature.push_back(curvature - 1.1);
		}
		return eachOfcurvature;
	}

	vector<double> recover_consin(const vector<double> &former)
	{
		vector<double> real_cos;
		int size_f = former.size();
		for (int i = 0; i < size_f; i++)
		{
			if (former[i] <= 0) real_cos.push_back(former[i] + 1.1);
			else real_cos.push_back(former[i] - 1.1);
		}
		return real_cos;
	}


	vector<int> most_convex_p(vector<Point2f> contour_, vector<double> cont_c, int max_cur_num)
	{
		vector<int> index_num;
		int contoursize = contour_.size();
		for (int i = 0; i < contoursize; i++)
		{
			index_num.push_back(i);
		}
		vector<double> cont_c1 = recover_consin(cont_c);
		//sort_comb(cont_c1, index_num);
		sort_comb(cont_c, index_num);
		//for (int i = 0; i < 50; i++)
		//	cout << cont_c[i] << endl;

		vector<int> cand_points_index;
		int t = 1;
		cand_points_index.push_back(index_num[0]);
		double length = contour_length(contour_);

		for (int i = 1; i < contoursize; i++)
		{
			if (t >= max_cur_num) break;
			else
			{
				int flag = 0;
				for (vector<int>::iterator it = cand_points_index.begin(); it != cand_points_index.end(); it++)
				{
					double leng = length_two_point2f(contour_[index_num[i]], contour_[*it]);
					if (leng < 0.008*length)
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
		if (cand_points_index.size()<max_cur_num)
			cout << "cand_points_index: " << cand_points_index.size() <<"  <max_cur_num" <<endl;
		return cand_points_index;
	}

	vector<Point2f> extract_contour(vector<Point2f> contour_, vector<int> mark_p, vector<int> &midmark_p)
	{
		Point2f line1 = contour_[mark_p[2]] - contour_[mark_p[0]];
		Point2f line2 = contour_[mark_p[3]] - contour_[mark_p[1]];
		//Point2f cente = center_p(inter_mid);
		int csize = contour_.size();
		vector<vector<Point2f>> four_place;
		vector<Point2f> one_loca;
		// 提取围成的轮廓，目前为止只考虑正向摆放，不考虑旋转和翻转
		four_place.push_back(contour_);
		for (int i = 0; i < csize; i++)
		{
			one_loca.push_back(contour_[i] + line1);
		}
		four_place.push_back(one_loca);
		one_loca.swap(vector<Point2f>());
		for (int i = 0; i < csize; i++)
		{
			one_loca.push_back(contour_[i] + line2);
		}
		four_place.push_back(one_loca);
		one_loca.swap(vector<Point2f>());
		for (int i = 0; i < csize; i++)
		{
			one_loca.push_back(contour_[i] + line1 + line2);
		}
		four_place.push_back(one_loca);
		int total_num = 0;
		midmark_p.push_back(0);
		vector<Point2f> morphed_B;
		for (int t = mark_p[3]; t > mark_p[2]; t--)
		{
			total_num++;
			morphed_B.push_back(four_place[0][t]);
		}
		midmark_p.push_back(total_num);
		for (int t = mark_p[0]; t >= 0; t--)
		{
			total_num++;
			morphed_B.push_back(four_place[1][t]);
		}
		for (int t = csize - 1; t > mark_p[3]; t--)
		{
			total_num++;
			morphed_B.push_back(four_place[1][t]);
		}
		midmark_p.push_back(total_num);
		for (int t = mark_p[1]; t > mark_p[0]; t--)
		{
			total_num++;
			morphed_B.push_back(four_place[3][t]);
		}
		midmark_p.push_back(total_num);
		for (int t = mark_p[2]; t > mark_p[1]; t--)
		{
			morphed_B.push_back(four_place[2][t]);
		}
		return morphed_B;
	}

	double Tiling_opt::warpAff_tra(vector<Point2f> &input_, vector<Point2f> &output_)
	{

		Point2f srcTri[3];
		Point2f dstTri[3];

		srcTri[0] = input_[0];
		srcTri[1] = input_[input_.size() - 1];
		dstTri[0] = Point2f(0, 0);//Point2f(50, 50); 
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		dstTri[1] = Point2f(length_two_point2f(srcTri[0], srcTri[1]), 0); //Point2f(length_two_point2f(srcTri[0], srcTri[1]) + 50, 50);

		int half = input_.size() / 2 - 1;
		double dis_line = 0;
		double dis_x = 0;
		double line_k = 0;
		double line_b = 0;
		double line_x = 0;
		if (abs(srcTri[1].x - srcTri[0].x) <= 0.01) //line 为x=srcTri[1].x
		{
			line_x = srcTri[1].x;
			dis_line = input_[half].x - line_x;
			while ((dis_line == 0) && (half < (input_.size() - 1))) //dis_line==0 说明共线
			{
				half++;
				dis_line = input_[half].x - line_x;
			}
			if (dis_line == 0)
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			dis_x = abs(input_[0].y - input_[half].y);
			srcTri[2] = input_[half];
			if (dis_line > 0)
			{
				if (srcTri[0].y > srcTri[1].y) dstTri[2] = Point2f(dis_x, dis_line);
				else dstTri[2] = Point2f(dis_x, -dis_line);
			}
			else if (dis_line < 0)
			{
				if (srcTri[0].y > srcTri[1].y) dstTri[2] = Point2f(dis_x, dis_line);
				else dstTri[2] = Point2f(dis_x, -dis_line);
			}
		}
		else{                                //line 为y=kx+b
			line_k = (srcTri[1].y - srcTri[0].y) / (srcTri[1].x - srcTri[0].x);
			line_b = srcTri[1].y - line_k*srcTri[1].x;
			dis_line = (line_k*input_[half].x - input_[half].y + line_b) / sqrt(line_k*line_k + 1);
			while ((dis_line == 0) && (half < (input_.size() - 1))) //dis_line==0 说明共线
			{
				half++;
				dis_line = (line_k*input_[half].x - input_[half].y + line_b) / sqrt(line_k*line_k + 1);
			}
			if (dis_line == 0)
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			dis_x = sqrt(((input_[half].x - srcTri[0].x)*(input_[half].x - srcTri[0].x) + (input_[half].y - srcTri[0].y)*(input_[half].y - srcTri[0].y)) - dis_line*dis_line);
			srcTri[2] = input_[half];
			if (dis_line > 0)
			{
				if (srcTri[0].x < srcTri[1].x) dstTri[2] = Point2f(dis_x, -dis_line);
				else dstTri[2] = Point2f(dis_x, dis_line);
			}
			else if (dis_line < 0)
			{
				if (srcTri[0].x < srcTri[1].x) dstTri[2] = Point2f(dis_x, -dis_line);
				else dstTri[2] = Point2f(dis_x, dis_line);
			}
		}


		Mat warp_mat(2, 3, CV_32FC1);

		/*for (int i = 0; i < 3; i++)
		{
		cout << "srcTri["<<i<<"]: " << srcTri[i] << endl;
		cout << "dstTri[" << i << "]: " << dstTri[i] << endl;
		}*/
		warp_mat = getAffineTransform(srcTri, dstTri);

		//cout << warp_mat << endl;

		cv::transform(input_, output_, warp_mat);//transform
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		//for (int i = 0; i < input_.size(); i++)
		//{
		//	cout << "src[" << i << "]:" << input_[i] << endl;
		//	cout << "dst[" << i << "]:" << output_[i] << endl;

		//}
		return length_two_point2f(srcTri[0], srcTri[1]) / 2;
	}

	double Tiling_opt::warpAff_tra_sec(vector<Point2f> &input_, vector<Point2f> &output_, vector<char>&second_char, vector<char>&second_char_out)
	{
		double shift = warpAff_tra(input_, output_);
		int num = output_.size();
		for (int i = 0; i < output_.size(); i++)
		{
			output_[i].x = output_[output_.size() - 1].x - output_[i].x;
			output_[i].y = -output_[i].y;
		}
		for (int j = 0; j < output_.size() / 2; j++)
		{
			Point2f temp = output_[j];
			output_[j] = output_[num - j - 1];
			output_[num - j - 1] = temp;
		}
		for (int t = second_char.size() - 1; t >= 0; t--)
		{
			second_char_out.push_back(second_char[t]);
		}
		return shift;
	}

	double Tiling_opt::warpAff_tra_ref_y(vector<Point2f> &input_, vector<Point2f> &output_)
	{
		double shift = warpAff_tra(input_, output_);
		for (int i = 0; i < output_.size(); i++)
		{
			output_[i].y = -output_[i].y;
		}
		return shift;
	}

	double Tiling_opt::re_warp_Aff(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		Point2f srcTri[3];
		Point2f dstTri[3];
		srcTri[0] = input_[0];
		srcTri[1] = input_[input_.size() - 1];
		if (abs(srcTri[0].y) <= 0.01 && abs(srcTri[0].x) <= 0.01 && abs(srcTri[1].y) <= 0.01) cout << "ok________________" << endl;
		else {
			cout << "gg" << endl;
			return 0;
		}
		dstTri[0] = start;//Point2f(50, 50); 
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		dstTri[1] = end; //Point2f(length_two_point2f(srcTri[0], srcTri[1]) + 50, 50);

		int half = input_.size() / 2 - 1;
		double dis_line = input_[half].y;
		double line_k = 0;
		double line_b = 0;

		while ((abs(dis_line) <= 0.01) && (half < (input_.size() - 1))) //dis_line==0 说明共线
		{
			half++;
			dis_line = input_[half].y;
		}
		if (dis_line == 0)
		{
			cout << "Error: no enough Noncollinear Point!" << endl;
			return 0;
		}
		srcTri[2] = input_[half];

		//dis_line正负不影响算式
		if (abs(start.x - end.x) <= 0.001)  //line 为x=start.x
		{
			if (start.y>end.y)
			{
				dstTri[2].x = start.x + dis_line;
				dstTri[2].y = start.y - input_[half].x;//此处默认input_为0点对齐
			}
			else
			{
				dstTri[2].x = start.x - dis_line;
				dstTri[2].y = start.y + input_[half].x;//此处默认input_为0点对齐
			}
		}
		else   //line为y=kx+b,点到直线的距离为|kx-y+b|/sart(k*k+1)
		{
			line_k = (end.y - start.y) / (end.x - start.x);
			line_b = end.y - end.x*line_k;
			Point2f vec = unit_vec(Point2f(1, line_k));
			Point2f mid;
			if (start.x < end.x)
			{
				mid = start + input_[half].x*vec;
			}
			else
			{
				mid = start - input_[half].x*vec;
			}
			Point2f v = unit_vec(Point2f(1, -1 / line_k));
			if (((start.x - end.x) *line_k)<0)  //异号，s.x<e.x且k>0 或者s.x>e.x且k<0  ,p=mid-d*v
			{
				dstTri[2] = mid - dis_line * v;
			}
			else
			{
				dstTri[2] = mid + dis_line * v;
			}
		}

		Mat warp_mat(2, 3, CV_32FC1);

		/*for (int i = 0; i < 3; i++)
		{
		cout << "srcTri["<<i<<"]: " << srcTri[i] << endl;
		cout << "dstTri[" << i << "]: " << dstTri[i] << endl;
		}*/
		warp_mat = getAffineTransform(srcTri, dstTri);

		//cout << warp_mat << endl;

		cv::transform(input_, output_, warp_mat);//transform
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		//for (int i = 0; i < input_.size(); i++)
		//{
		//	cout << "src[" << i << "]:" << input_[i] << endl;
		//	cout << "dst[" << i << "]:" << output_[i] << endl;

		//}
		return length_two_point2f(dstTri[0], dstTri[1]) / 2;
	}

	double Tiling_opt::re_warp_Aff_sec(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		int num = input_.size();
		for (int i = 0; i < input_.size(); i++)
		{
			input_[i].x = input_[input_.size() - 1].x - input_[i].x;
			input_[i].y = -input_[i].y;
		}
		for (int j = 0; j < input_.size() / 2; j++)
		{
			Point2f temp = input_[j];
			input_[j] = input_[num - j - 1];
			input_[num - j - 1] = temp;
		}
		//for (int t = second_char.size() - 1; t >= 0; t--)
		//{
		//	second_char_out.push_back(second_char[t]);
		//}
		double shift = re_warp_Aff(input_, output_, start, end);
		return shift;
	}

	double Tiling_opt::re_warp_Aff_ref(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		for (int i = 0; i < input_.size(); i++)
		{
			input_[i].y = -input_[i].y;
		}
		double shift = re_warp_Aff(input_, output_, start, end);
		return shift;
	}

	double Tiling_opt::warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		if (input_[0] == start) cout << "ok!!----" << endl;
		double scale = length_two_point2f(start, end) / length_two_point2f(input_[0], input_[input_.size() - 1]);
		Point2f srcTri[3];
		Point2f dstTri[3];
		srcTri[0] = input_[0];
		srcTri[1] = input_[input_.size() - 1];
		dstTri[0] = start;//Point2f(50, 50); 
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		dstTri[1] = end; //Point2f(length_two_point2f(srcTri[0], srcTri[1]) + 50, 50);

		//先求dis_line和dis_x,在这里加上一个flag判断第三点在线的哪一侧
		double dis_line = 0;
		double line_k = 0;
		double line_b = 0;
		double dis_x;
		int half = input_.size() / 2 - 1;
		bool left = true;
		if (abs(srcTri[0].x - srcTri[1].x) < 0.001)  //line 为x=start.x
		{
			dis_line = input_[half].x - srcTri[0].x;
			while ((abs(dis_line) <= 0.01) && (half < (input_.size() - 1))) //dis_line==0 说明共线
			{
				half++;
				dis_line = input_[half].y;
			}
			if (dis_line == 0)
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			dis_x = sqrt(length_two_point2f(input_[half], srcTri[0])*length_two_point2f(input_[half], srcTri[0]) - dis_line*dis_line);
			if (dis_line > 0)
			{
				if (srcTri[0].y < srcTri[1].y)
					left = false;
			}
			else
			{
				if (srcTri[0].y > srcTri[1].y)
					left = false;
			}
			dis_line = abs(dis_line);
		}
		else  //line为y=kx+b,点到直线的距离为|kx-y+b|/sart(k*k+1)
		{
			line_k = (srcTri[1].y - srcTri[0].y) / (srcTri[1].x - srcTri[0].x);
			line_b = srcTri[0].y - srcTri[0].x*line_k;
			dis_line = (line_k*input_[half].x - input_[half].y + line_b) / sqrt(line_k*line_k + 1); //正数说明在下方
			while ((dis_line == 0) && (half < (input_.size() - 1))) //dis_line==0 说明共线
			{
				half++;
				dis_line = (line_k*input_[half].x - input_[half].y + line_b) / sqrt(line_k*line_k + 1);
			}
			if (dis_line == 0)
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			dis_x = sqrt(length_two_point2f(input_[half], srcTri[0])*length_two_point2f(input_[half], srcTri[0]) - dis_line*dis_line);
			if (dis_line > 0)
			{
				if (srcTri[0].x < srcTri[1].x)
					left = false;
			}
			else
			{
				if (srcTri[0].x > srcTri[1].x)
					left = false;
			}
			dis_line = abs(dis_line);
		}
		srcTri[2] = input_[half];
		dis_x = dis_x*scale;
		dis_line = dis_line*scale;
		//确定第三点的目标位置
		if (abs(start.x - end.x) <= 0.001)  //line 为x=start.x
		{
			if (start.y>end.y)
			{
				if (left)
				{
					dstTri[2].x = start.x + dis_line;
					dstTri[2].y = start.y - dis_x;
				}
				else
				{
					dstTri[2].x = start.x - dis_line;
					dstTri[2].y = start.y - dis_x;
				}
			}
			else
			{
				if (left)
				{
					dstTri[2].x = start.x - dis_line;
					dstTri[2].y = start.y + dis_x;
				}
				else
				{
					dstTri[2].x = start.x + dis_line;
					dstTri[2].y = start.y + dis_x;
				}
			}
		}
		else   //line为y=kx+b,点到直线的距离为|kx-y+b|/sart(k*k+1)
		{
			line_k = (end.y - start.y) / (end.x - start.x);
			Point2f vec = unit_vec(Point2f(1, line_k));
			Point2f mid;
			if (start.x < end.x)
			{
				mid = start + dis_x*vec;
			}
			else
			{
				mid = start - dis_x*vec;
			}
			Point2f v = unit_vec(Point2f(1, -1 / line_k));
			if (!left) dis_line = -dis_line;

			if (((start.x - end.x) *line_k)<0)  //异号，s.x<e.x且k>0 或者s.x>e.x且k<0  ,p=mid-d*v
			{
				dstTri[2] = mid - dis_line * v;
			}
			else
			{
				dstTri[2] = mid + dis_line * v;
			}
		}

		Mat warp_mat(2, 3, CV_32FC1);

		/*for (int i = 0; i < 3; i++)
		{
		cout << "srcTri["<<i<<"]: " << srcTri[i] << endl;
		cout << "dstTri[" << i << "]: " << dstTri[i] << endl;
		}*/
		warp_mat = getAffineTransform(srcTri, dstTri);

		//cout << warp_mat << endl;

		cv::transform(input_, output_, warp_mat);//transform
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		//for (int i = 0; i < input_.size(); i++)
		//{
		//	cout << "src[" << i << "]:" << input_[i] << endl;
		//	cout << "dst[" << i << "]:" << output_[i] << endl;

		//}
		return length_two_point2f(dstTri[0], dstTri[1]) / 2;
	}

	double Tiling_opt::Aff_place(vector<Point2f> &input1, vector<Point2f> &input2, vector<vector<Point2f>> &prototwo, vector<Point2f> &protoAff, int flag = 0)
	{
		Point2f srcTri[3];
		Point2f dstTri[3];
		vector<Point2f> output;
		if (flag == 0)  //flag=0,sec
		{
			int half = input2.size() / 2 - 1;
			srcTri[0] = input2[0];
			srcTri[1] = input2[input2.size() - 1];

			dstTri[0] = input1[input1.size() - 1];//Point2f(50, 50); 
			dstTri[1] = input1[0];

			while ((unit_vec(input2[half] - input2[0]) == unit_vec(input2[input2.size() - 1] - input2[half])) && half<input2.size() - 1)
			{
				half++;
			}
			if (half == (input2.size() - 1))
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			srcTri[2] = input2[half];
			dstTri[2] = input1[input1.size() - 1 - half];
		}
		else   //flag=1,ref
		{
			int half = input2.size() / 2 - 1;
			srcTri[0] = input2[0];
			srcTri[1] = input2[input2.size() - 1];

			dstTri[0] = input1[0];
			dstTri[1] = input1[input1.size() - 1];

			while ((unit_vec(input2[half] - input2[0]) == unit_vec(input2[input2.size() - 1] - input2[half])) && half<input2.size() - 1)
			{
				half++;
			}
			if (half == (input2.size() - 1))
			{
				cout << "Error: no enough Noncollinear Point!" << endl;
				return 0;
			}
			srcTri[2] = input2[half];
			dstTri[2] = input1[half];
		}
		Mat warp_mat(2, 3, CV_32FC1);

		warp_mat = getAffineTransform(srcTri, dstTri);

		////cout << warp_mat << endl;
		////检测output与input2是否一致
		//cv::transform(input2, output, warp_mat);//transform
		//Mat drawing1 = Mat::zeros(800, 800, CV_8UC3);
		//for (int j = 0; j < input1.size(); j++)
		//{
		//	circle(drawing1, input1[j], 1, Scalar(255, 0, 0), -1);
		//}
		//Mat drawing2 = Mat::zeros(800, 800, CV_8UC3);
		//for (int j = 0; j < input2.size(); j++)
		//{
		//	circle(drawing2, input2[j], 1, Scalar(0, 255, 0), -1);
		//}
		//Mat drawing3 = Mat::zeros(800, 800, CV_8UC3);
		//for (int j = 0; j < output.size(); j++)
		//{
		//	circle(drawing3, output[j], 1, Scalar(0, 0, 255), -1);
		//}
		//imshow("1:", drawing1);
		//imshow("2:", drawing2);
		//imshow("out:", drawing3);


		//将整个prototwo与input1映齐
		Mat drawing2 = Mat::zeros(1600, 1600, CV_8UC3);
		Mat drawing3 = Mat::zeros(1600, 1600, CV_8UC3);
		for (int i = 0; i < 4; i++)
		{
			vector<Point2f> out_;
			cv::transform(prototwo[i], out_, warp_mat);
			if (!protoAff.empty())
			{
				Point2f pop = protoAff[protoAff.size() - 1];
				if (pop == out_[0])
				{
					protoAff.pop_back();
				}
			}
			for (int j = 0; j < out_.size(); j++)
			{
				protoAff.push_back(out_[j]);
			}

			//for (int j = 0; j < prototwo[i].size() - 1; j++)
			//{
			//	cout << "proto[i]:" << prototwo[i][j] << endl;
			//	MyLine(drawing3, prototwo[i][j], prototwo[i][j + 1], "green");
			//	//circle(drawing2, protoAff[i][j], 1, Scalar(0, 255, 0), -1);
			//}
			//for (int j = 0; j < protoAff[i].size()-1; j++)
			//{
			//	cout << "protoAff[i]:" << protoAff[i][j] << endl;
			//	MyLine(drawing2, protoAff[i][j], protoAff[i][j + 1], "blue");
			//	//circle(drawing2, protoAff[i][j], 1, Scalar(0, 255, 0), -1);
			//}
			//cv::transform(prototwo[i], protoAff[i], warp_mat);
		}
		//imshow("out11:", drawing2);
		//imshow("out111:", drawing3);

		return 0;
	}

	double Tiling_opt::DTW(vector<Point2f> &first_arr, vector<Point2f> &second_arr)
	{
		//ofstream out("D:\\VisualStudioProjects\\manual_record\\dtw.txt");
		int first_num = first_arr.size();
		int second_num = second_arr.size();

		double dis[50][50];
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis[i][j] = abs(length_two_point2f(first_arr[i], second_arr[j]));
			}
		}
		/*for (int i = 0; i < first_num; i++)
		{
		cout << "dis[" << i << "]: ";
		for (int j = 0; j < second_num; j++)
		{
		if (j % 10 == 0)
		cout << endl;
		cout <<  dis[i][j]<<"  " ;

		}
		cout << endl;
		}*/
		//out.close();
		double distance[50][50];
		int step[50][50];//记录总的步数
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				distance[i][j] = 0;
				step[i][j] = 0;
			}
		}
		distance[0][0] = (double)2 * dis[0][0];
		step[0][0] = 0;

		for (int i = 1; i < first_num; i++)
		{
			distance[i][0] = distance[i - 1][0] + dis[i][0];
			step[i][0] = step[i - 1][0] + 1;
		}
		for (int i = 1; i < second_num; i++)
		{
			distance[0][i] = distance[0][i - 1] + dis[0][i];
			step[0][i] = step[0][i - 1] + 1;
		}

		for (int j = 1; j < second_num; j++)
		{

			for (int i = 1; i < first_num; i++)//(int i = istart; i <= imax; i++)
			{
				double g1 = distance[i - 1][j] + dis[i][j];
				double g2 = distance[i - 1][j - 1] + 2 * dis[i][j];
				double g3 = distance[i][j - 1] + dis[i][j];
				if (g1 < g2)
				{
					if (g1 < g3)
					{
						distance[i][j] = g1;
						step[i][j] = step[i - 1][j] + 1;
					}
					else
					{
						distance[i][j] = g3;
						step[i][j] = step[i][j - 1] + 1;
					}
				}
				else
				{
					if (g2 < g3)
					{
						distance[i][j] = g2;
						step[i][j] = step[i - 1][j - 1];
					}
					else
					{
						distance[i][j] = g3;
						step[i][j] = step[i][j - 1] + 1;
					}
				}
			}

		}

		/*for (int i = 0; i < first_num; i++)
		{
		cout << "distance[" << i << "]: ";
		for (int j = 0; j < second_num; j++)
		{
		cout << distance[i][j] << " ";

		}
		cout << endl;
		}	*/
		/*for (int i = 0; i < first_num; i++)
		{
		cout << "step[" << i << "]: ";
		for (int j = 0; j < second_num; j++)
		{
		cout << step[i][j] << " ";

		}
		cout << endl;
		}*/
		//dp_path.swap(vector<pair<int, int>>());
		//printPath(dis, distance, first_num - 1, second_num - 1, first_arr, second_arr, dp_path);
		// show the dp_path
		/*for (int i = 0; i < dp_path.size(); i++)
		{
		cout << first_arr[dp_path[i].first] << " - " << second_arr[dp_path[i].second] << endl;
		}*/
		return distance[first_num - 1][second_num - 1];
	}

		

	
	// Morph points
	void MorphPoints(const vector<Point2f>& srcPts1, const vector<Point2f>& srcPts2, vector<Point2f>& dstPts, float s)
	{
		assert(srcPts1.size() == srcPts2.size());

		int num_pts = srcPts1.size();

		dstPts.resize(num_pts);
		for (int i = 0; i<num_pts; i++){
			dstPts[i].x = (1.0 - s) * srcPts1[i].x + s * srcPts2[i].x;
			dstPts[i].y = (1.0 - s) * srcPts1[i].y + s * srcPts2[i].y;
		}
	}


	string int2string(int number)
	{
		char ch[8];
		int  divisor = 1000;
		int index = 0;
		if (number == 0)
		{
			ch[index] = '0';
			ch[index + 1] = '\0';
			return ch;
		}
		while (divisor != 0)
		{
			if (number / divisor == 0)
			{
				divisor = divisor / 10;
				if (index != 0)
				{
					ch[index] = '0';
					index++;
				}

			}
			else
			{
				ch[index] = number / divisor + '0';
				index++;
				number = number % divisor;
				divisor = divisor / 10;
			}
		}
		ch[index] = '\0';
		return ch;
	}

	vector<Point2f> sampling(vector<Point2f> &contour_, int points_num)
	{
		double length = contour_length(contour_);
		int sam_num = points_num * 100;
		double Lambda = length / sam_num;

		vector<Point2f> contour_sam;
		Point2f sample;

		contour_sam.push_back(contour_[0]);
		sample = contour_[0];
		for (int t = 1; t < contour_.size(); t++)
		{
			double length_ = length_two_point2f(sample, contour_[t]);
			if (length_ > Lambda)
			{
				Point2f vec = unit_vec(contour_[t] - sample);
				sample = sample + Lambda * vec;
				contour_sam.push_back(sample);
				t = t - 1;
			}
			else if (t < contour_.size() - 1)
			{
				while ((length_ + length_two_point2f(contour_[t], contour_[t + 1])) < Lambda)
				{
					length_ = length_ + length_two_point2f(contour_[t], contour_[t + 1]);
					t++;
					if (t >= (contour_.size() - 1)) break;
				}
				if (t >= (contour_.size() - 1)) break;
				Point2f vec = unit_vec(contour_[t + 1] - contour_[t]);
				sample = contour_[t] + (Lambda - length_) * vec;
				contour_sam.push_back(sample);
			}
		}
		if (length_two_point2f(contour_sam[0], contour_sam[contour_sam.size() - 1])<1) contour_sam.pop_back();

		return contour_sam;
	}
}
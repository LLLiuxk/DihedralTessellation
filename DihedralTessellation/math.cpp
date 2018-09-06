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

	void sort_bub(vector<int> &target)  //从小到大
	{
		int i, j;
		int temp;
		for (i = 0; i < target.size() - 1; i++)
			for (j = 0; j < target.size() - 1 - i; j++)
				if (target[j] > target[j + 1])
				{
					temp = target[j];
					target[j] = target[j + 1];
					target[j + 1] = temp;
					
				}
	}

	
	double cur_length_two_p(double cur1, double cur2, double zeta)
	{
		double mis = (cur1 - cur2)*(cur1 - cur2)*zeta;
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

	vector<double> Prototile::curvature_com_k(vector<Point2f> &contour_sam) //k等于角度的差值比弧长的差值（近似为三角的斜边）
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

	vector<double> Prototile::curvature_com(vector<Point2f> &contour_sam) 
	{
		vector<double> eachOfcurvature;
		int c_s = contour_sam.size();
		double curvature = cos_two_vector(contour_sam[c_s - 1] - contour_sam[0], contour_sam[1] - contour_sam[0]);
		eachOfcurvature.push_back(curvature);
		int i = 1;
		for (; i < c_s - 1; i++)
		{
			curvature = cos_two_vector(contour_sam[i - 1] - contour_sam[i], contour_sam[i + 1] - contour_sam[i]);
			eachOfcurvature.push_back(curvature);
		}
		curvature = cos_two_vector(contour_sam[i - 1] - contour_sam[i], contour_sam[0] - contour_sam[i]);
		eachOfcurvature.push_back(curvature);
		return eachOfcurvature;
	}

	void Prototile::sort_cos(vector<double> vect, vector<int> &index_num) //只保留下标的排序,从小到大
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
		for (int t = second_char.size()-1; t >= 0; t--)
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
		double shift = re_warp_Aff(input_, output_,start,end);
		return shift;
	}

	double Tiling_opt::re_warp_Aff_ref(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		for (int i = 0; i < input_.size(); i++)
		{
			input_[i].y = -input_[i].y;
		}
		double shift = re_warp_Aff(input_, output_,start,end);
		return shift;
	}

	double Tiling_opt::warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
	{
		if (input_[0] == start) cout << "ok!!----" << endl;
		double scale = length_two_point2f(start, end) / length_two_point2f(input_[0], input_[input_.size()-1]);
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
		dp_path.swap(vector<pair<int, int>>());
		printPath(dis, distance, first_num - 1, second_num - 1, first_arr, second_arr, dp_path);
		// show the dp_path
		/*for (int i = 0; i < dp_path.size(); i++)
		{
			cout << first_arr[dp_path[i].first] << " - " << second_arr[dp_path[i].second] << endl;
		}*/
		return distance[first_num - 1][second_num - 1];
	}

	void Tiling_opt::printPath(double d[][50], double dp[][50], int i, int j, vector<Point2f> &first_arr, vector<Point2f> &second_arr, vector<pair<int, int>>& path)
	{

		if (i == 0 && j == 0) {
			//cout << first_arr[i] << " - " << second_arr[j] << endl;
			path.push_back(make_pair(i, j));
			return;
			//cout make_pair(i,j);
		}

		if (abs(dp[i][j] - (dp[i - 1][j - 1] + 2 * d[i][j])) < 0.1){
			printPath(d, dp, i - 1, j - 1, first_arr, second_arr, path);

		}
		else if (abs(dp[i][j] - (dp[i][j - 1] + d[i][j])) < 0.1) {
			printPath(d, dp, i, j - 1, first_arr, second_arr, path);

		}
		else {
			printPath(d, dp, i - 1, j, first_arr, second_arr, path);
		}
		path.push_back(make_pair(i, j));
		//cout << first_arr[i] << " - " << second_arr[j] << endl;

	}

	double Tiling_opt::quadr_mismatch(vector<Point2f> &first_arr, vector<Point2f> &second_arr, vector<char> &first_char, vector<char> &second_char)
	{
		//ofstream out("D:\\VisualStudioProjects\\manual_record\\dtw.txt");
		int first_num = first_arr.size();
		int second_num = second_arr.size();

		//cout << "\n        first.size: " << first_num << "  -----    chararr_size" << first_char.size() << endl;

		double dis[50][50];//两组点之间的坐标差异
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis[i][j] = length_two_point2f(first_arr[i], second_arr[j]);
			}
		}
		int dis_cur[50][50];//两组点之间的曲率差异
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis_cur[i][j] = cur_char_length(first_char[i], second_char[j]);
				//dis_cur[i][j] = dis_cur[i][j] * dis_cur[i][j];
			}
		}

		double distance[100][100];
		int step[100][100];//记录总的步数
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				distance[i][j] = 0;
				step[i][j] = 0;
			}
		}
		distance[0][0] = (double)2 * (dis[0][0] + dis_cur[0][0]);
		step[0][0] = 0;

		for (int i = 1; i < first_num; i++)
		{
			distance[i][0] = distance[i - 1][0] + (dis[i][0] + dis_cur[i][0]);
			step[i][0] = step[i - 1][0] + 1;
		}
		for (int i = 1; i < second_num; i++)
		{
			distance[0][i] = distance[0][i - 1] + (dis[0][i] + dis_cur[0][i]);
			step[0][i] = step[0][i - 1] + 1;
		}

		for (int j = 1; j < second_num; j++)
		{

			for (int i = 1; i < first_num; i++)//(int i = istart; i <= imax; i++)
			{
				double g1 = distance[i - 1][j] + dis[i][j] + dis_cur[i][j];// +kappa;
				double g2 = distance[i - 1][j - 1] + 2 * (dis[i][j] + dis_cur[i][j]);
				double g3 = distance[i][j - 1] + dis[i][j] + dis_cur[i][j];// +kappa;
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
		return distance[first_num - 1][second_num - 1];
	}


	
}
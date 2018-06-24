#include "tilingOpt.h"


namespace Tiling_tiles{

	Tiling_opt::Tiling_opt()
	{
		prototile_first = new Prototile();
		prototile_second = new Prototile();
		//memset(dp, 0, sizeof(dp));
		//memset(dp_inver, 0, sizeof(dp_inver));

	}

	Tiling_opt::~Tiling_opt()
	{
		delete[] prototile_first;
		delete[] prototile_second;
	}

	double Tiling_opt::com_scale_factor()
	{
		double scale_factor; //   first/second
		cout << "length1: " << prototile_first->c_length << "length2: " << prototile_second->c_length << endl;
		scale_factor = prototile_first->c_length / prototile_second->c_length;
		cout << "scale factor: " << scale_factor << endl;
		return scale_factor;
	}

	//每一个轮廓分成四段
	void Tiling_opt::com_score(string imagename1, string imagename2)
	{
		//vector<Each_interval> initial_pair;
		//Each_interval candidate[50];

		//prototile_first->loadTileData(imagename1);
		//prototile_second->loadTileData(imagename2);

		//int sam_num = 0;

		//vector<Point2f> contour1 = prototile_first->contour_sample[sam_num];
		//vector<Point2f> contour2 = prototile_second->contour_sample[sam_num];
		//vector<double> contour1_cur = prototile_first->contour_curva[sam_num];
		//vector<double> contour2_cur = prototile_second->contour_curva[sam_num];
		//char * contour1_str = prototile_first->cur_string[sam_num];
		//char * contour2_str = prototile_second->cur_string[sam_num];
		////对齐
		//double factor = com_scale_factor();
		//for (int i = 0; i <contour2.size(); i++)
		//{
		//	contour2[i].x = contour2[i].x * factor;
		//	contour2[i].y = contour2[i].y * factor;
		//}

		//int interval_num1_first = contour1.size() / 10;


		////int intervals_second = 4; //暂定分为四段
		////int all_num_second = 0;
		//int interval_num1_second = contour2.size() / 10;///prototile_second->contour_sample[0].size() / intervals_second;

		//for (int n = 0; n < contour1.size(); n++)
		//{
		//	for (int m = 0; m < contour2.size(); m++)
		//	{
		//		int candi_length = 32;
		//		vector<Point2f> first_;
		//		vector<Point2f> second_;
		//		//vector<double> first_cur;
		//		//vector<double> second_cur;
		//		vector<char> first_char;
		//		vector<char> second_char;


		//		int flag = 1; //0 代表sec,正序无翻折,1 代表ref
		//		for (int t = 0; t < interval_num1_first; t++)
		//		{
		//			first_.push_back(contour1[(n + t) % contour1.size()]);
		//			//first_cur.push_back(contour1_cur[(n + t) % contour1_cur.size()]);
		//			first_char.push_back(contour1_str[(n + t) % contour1.size()]);
		//		}
		//		for (int s = 0; s < interval_num1_second; s++)
		//		{
		//			second_.push_back(contour2[(m + s) % contour2.size()]);
		//			//second_cur.push_back(contour2_cur[(m + s) % contour2_cur.size()]);
		//			second_char.push_back(contour1_str[(m + s) % contour2.size()]);
		//		}

		//		double score = 0;
		//		//double score2 = 0;
		//		vector<Point2f> first_after_aff;
		//		vector<Point2f> second_after_aff;
		//		vector<Point2f> second_after_aff_ref;
		//		//用仿射变化将两段轮廓对齐，消除初始位置的影响

		//		double re_first = warpAff_tra(first_, first_after_aff);
		//		double re_second = warpAff_tra_sec(second_, second_after_aff);
		//		double re_second_ref = warpAff_tra_ref_y(second_, second_after_aff_ref);
		//		double length = (first_.size() + second_.size()) / 2;
		//		//score = DTW(first_after_aff, second_after_aff_ref);
		//		score = quadr_mismatch(first_after_aff, second_after_aff, first_char, second_char);
		//		Each_interval result_sec(n, interval_num1_first, m, interval_num1_second, score, 0);
		//		//initial_pair.push_back(result_sec);
		//		for (int t = 0; t < candi_length; t++)
		//		{
		//			if (score<candidate[t].score)
		//			{
		//				for (int s = candi_length - 1; s > t; s--)
		//				{
		//					candidate[s] = candidate[s - 1];
		//				}
		//				candidate[t] = result_sec;
		//				break;
		//			}
		//		}
		//		score = quadr_mismatch(first_after_aff, second_after_aff, first_char, second_char);
		//		Each_interval result_ref(n, interval_num1_first, m, interval_num1_second, score, 0);
		//		//initial_pair.push_back(result_ref);

		//		for (int t = 0; t < candi_length; t++)
		//		{
		//			if (score<candidate[t].score)
		//			{
		//				for (int s = candi_length - 1; s > t; s--)
		//				{
		//					candidate[s] = candidate[s - 1];
		//				}
		//				candidate[t] = result_ref;
		//				break;
		//			}
		//		}
		//		//double score = com_each_pair(first_, second_,flag);
		//		//Each_interval result(n, interval_num1_first, m, interval_num1_second, score, flag);

		//	}

		//}
		////排序并取前32对
		////sort_Inter(initial_pair);

		//Mat drawing5 = Mat::zeros(800, 800, CV_8UC3);
		//for (int j = 0; j < contour1.size(); j++)
		//{
		//	circle(drawing5, contour1[j], 1, Scalar(255, 0, 0), -1);
		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//for (int i = 0; i < 32; i++)
		//{
		//	cout << candidate[i].index_x << " ";
		//	circle(drawing5, contour1[candidate[i].index_x], 1, Scalar(255, 255, 0), -1);
		//}
		//imshow("candidate: ", drawing5);


		//Mat drawing6 = Mat::zeros(800, 800, CV_8UC3);
		//for (int j = 0; j < contour2.size(); j++)
		//{
		//	circle(drawing6, contour2[j], 1, Scalar(255, 0, 0), -1);
		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//for (int i = 0; i < 32; i++)
		//{
		//	cout << candidate[i].index_x << " ";
		//	circle(drawing6, contour2[candidate[i].index_y], 1, Scalar(255, 255, 0), -1);
		//}
		//imshow("candidate2: ", drawing6);
		////cout << "min: " << min << endl
		////	<< "n: " << min_x << endl
		////	<< "m: " << min_y << endl;
		////for (int i = 0; i < 4; i++ )
		////{
		////	cout << i << " -- " << all_types[min_x][min_y].second[i].first << " type: " << all_types[min_x][min_y].second[i].second<< endl;
		////}
		////

		//////把结果显示出来
		////proto_interval_first.swap(vector<vector<Point2f>>());
		////proto_interval_second.swap(vector<vector<Point2f>>());
		////divide_intervals(4, prototile_first->contour_sample[0], proto_interval_first, min_x);
		////divide_intervals(4, prototile_second->contour_sample[0], proto_interval_second, min_y);
		////
		////vector<Point2f> output_first;
		////vector<Point2f> output_second;

		////
		////for (int j = 0; j < 4; j++)
		////{
		////	//all_types[max_x][max_y].second[i];
		////	double re_first = warpAff_tra(proto_interval_first[j], output_first);
		////	double re_second = 0;
		////	if (all_types[min_x][min_y].second[j].second == 0)
		////	{
		////		re_second = warpAff_tra_sec(proto_interval_second[all_types[min_x][min_y].second[j].first], output_second);
		////	}
		////	else if (all_types[min_x][min_y].second[j].second == 1)
		////	{
		////		re_second = warpAff_tra_ref_y(proto_interval_second[all_types[min_x][min_y].second[j].first], output_second);
		////	}
		////	
		////	int size = 800;//展示图片的尺寸

		////	Mat drawing = Mat::zeros(size, size, CV_8UC3);
		////	for (int i = 0; i < output_first.size() - 1; i++)
		////	{
		////		//MyLine(drawing, first_interval[i], first_interval[i + 1], 1);
		////		MyLine(drawing, output_first[i] + Point2f(size / 2 - re_first, size / 2), output_first[i + 1] + Point2f(size / 2 - re_first, size / 2), "green");
		////	}

		////	for (int i = 0; i < output_second.size() - 1; i++)
		////	{
		////		//MyLine(drawing, first_interval[i], first_interval[i + 1], 1);
		////		MyLine(drawing, output_second[i] + Point2f(size / 2 - re_second, size / 2), output_second[i + 1] + Point2f(size / 2 - re_second, size / 2), "cyan");
		////	}
		////	//char i;
		////	string name = "the ";
		////	name = name + char(j + 48) + " pair: ";
		////	imshow(name, drawing);
		////}

		////

		////
		////int row = 800;
		////int col = 1600;
		////Mat drawing4 = Mat::zeros(row, col, CV_8UC3);
		//////Mat drawing5 = Mat::zeros(800, 800, CV_8UC3);
		//////namedWindow("simple111 contour", WINDOW_AUTOSIZE);

		//////用shift将两张图对齐到图像的左右中点
		////Point2f shift1 = Point2f(prototile_first->center_point.x - col/4, prototile_first->center_point.y - row/2);
		////Point2f shift2 = prototile_second->center_point - Point2f(3 * col / 4, row / 2);

		////for (int i = 0; i < 4; i++)
		////{
		////	//cout << "proto_interval_first[" << i << "]  num:" << proto_interval_first[i].size() << endl;
		////	string color;
		////	if (i == 0) color = "red";
		////	else if (i == 1) color = "blue";
		////	else if (i == 2) color = "green";
		////	else if (i == 3) color = "cyan";
		////	for (int j = 0; j < proto_interval_first[i].size() - 1; j++)
		////	{
		////		MyLine(drawing4, proto_interval_first[i][j] - shift1, proto_interval_first[i][j + 1] - shift1, color);

		////	}			
		////	for (int j = 0; j < proto_interval_second[all_types[min_x][min_y].second[i].first].size() - 1; j++)
		////	{
		////		MyLine(drawing4, proto_interval_second[all_types[min_x][min_y].second[i].first][j] - shift2, proto_interval_second[all_types[min_x][min_y].second[i].first][j + 1] - shift2, color);

		////	}
		////}
		////imshow("mapping two contours", drawing4);  

		////imshow("second contour: " + imagename2 + " intervals", drawing5);	
		////namedWindow("simple111 contour", WINDOW_AUTOSIZE);

	}
	double Tiling_opt::com_each_pair(vector<Point2f> &first_interval, vector<Point2f> &second_interval, int &flag)
	{
		double score = 0;
		double score2 = 0;
		vector<Point2f> first_after_aff;
		vector<Point2f> second_after_aff;
		vector<Point2f> second_after_aff_ref;
		//用仿射变化将两段轮廓对齐，消除初始位置的影响

		double re_first = warpAff_tra(first_interval, first_after_aff);
		//double re_second = warpAff_tra_sec(second_interval, second_after_aff);
		double re_second_ref = warpAff_tra_ref_y(second_interval, second_after_aff_ref);
		double length = (contour_length(first_interval) + contour_length(second_interval)) / 2;

		score = DTW(first_after_aff, second_after_aff);
		//score2 = DTW(first_after_aff, second_after_aff_ref);
		//cout << "sec score: " << score << endl << "ref score: " << score2 << endl;
		if (score < score2)
		{
			flag = 0;
			return score / length;
		}
		else
		{
			flag = 1;
			return score2 / length;
		}

	}

	void Tiling_opt::divide_intervals(int intervals, vector<Point2f> &contour, vector<vector<Point2f>> &proto_interval, int delay) //此处还是平均分的方法，需要做的： 优化分段的方法
	{
		int interval_num = contour.size() / intervals; //每段的点数
		int all_num = 0;
		for (int i = 0; i < intervals - 1; i++)
		{
			vector<Point2f> each_inter;
			for (int j = 0; j < interval_num; j++)
			{
				each_inter.push_back(contour[all_num + delay]);
				all_num++;
			}
			each_inter.push_back(contour[all_num + delay]);
			proto_interval.push_back(each_inter);
		}
		vector<Point2f> each_inter;
		for (; all_num < contour.size(); all_num++)
		{
			each_inter.push_back(contour[(all_num + delay) % (contour.size())]);
		}
		each_inter.push_back(contour[(all_num + delay) % (contour.size())]);
		proto_interval.push_back(each_inter);

		//cout << "proto_interval  num:" << proto_interval.size() << endl;
	}


	double Tiling_opt::com_optimal_score(vector<vector<Point2f>> &proto_interval_1, vector<vector<char>> &proto_first_char,
		vector<vector<Point2f>> &proto_interval_2, vector<vector<char>> &proto_second_char, vector<pair<int, int>> &order_type)
	{
		vector<int> order;
		pair<double, int> score_order[5][5];

		if (order.empty())
		{
			for (int i = 0; i < 4; i++)
				order.push_back(-1);
		}
		else {
			for (int i = 0; i < order.size(); i++)
				order[i] = -1;

		}

		for (int i = 0; i < proto_interval_1.size(); i++)
		{
			for (int j = 0; j < proto_interval_2.size(); j++)
			{
				pair<double, int> each_score;
				int flag = 0; //0 代表sec,1 代表ref
				double score = com_tra_sim(proto_interval_1[i], proto_first_char[i], proto_interval_2[j], proto_second_char[j], flag);
				each_score = make_pair(score, flag);
				score_order[i][j] = each_score;
			}
		}

		// 全部列举
		double total_score = 100000;
		for (int i = 0; i < 4; i++)
		{
			double total_sco = 0;
			total_sco += score_order[0][i].first;
			double num1 = total_sco;
			for (int j = 0; j < 4; j++)
			{
				total_sco = num1;
				if (j == i) continue;
				total_sco += score_order[1][j].first;
				double num2 = total_sco;
				for (int t = 0; t < 4; t++)
				{
					total_sco = num2;
					if (t == i || t == j) continue;
					total_sco += score_order[2][t].first;
					for (int n = 0; n < 4; n++)
					{
						if (n == i || n == j || n == t) continue;
						total_sco += score_order[3][n].first;
						if (total_sco < total_score)
						{
							total_score = total_sco;
							order[0] = i;
							order[1] = j;
							order[2] = t;
							order[3] = n;
						}
					}
				}
			}
		}

		for (int i = 0; i < 4; i++)
		{
			order_type.push_back(make_pair(order[i], score_order[i][order[i]].second));
		}

		// 将数组结果排序输出
		vector<pair<pair<int, int>, pair<double, int>>> score_;
		for (int i = 0; i <4; i++) //
		{
			for (int j = 0; j < 4; j++)
			{
				score_.push_back(make_pair(make_pair(i, j), make_pair(score_order[i][j].first, score_order[i][j].second)));
			}
		}
		pair<pair<int, int>, pair<double, int>> temp;
		for (int i = 0; i < score_.size() - 1; i++)
		{
			for (int j = 0; j < score_.size() - 1 - i; j++)
			{
				if (score_[j].second.first > score_[j + 1].second.first)
				{
					temp = score_[j];
					score_[j] = score_[j + 1];
					score_[j + 1] = temp;
				}
			}
		}
		for (int i = 0; i < score_.size(); i++) //show the order result
		{
			cout << score_[i].first.first << "---" << score_[i].first.second << ": " << score_[i].second.first << "  type: " << score_[i].second.second << endl;
		}
		//ofstream out("D:\\VisualStudioProjects\\manual_record\\manual_record.txt", ios::app);
		for (int i = 0; i <4; i++) //show the order result
		{
			for (int j = 0; j < 4; j++)
			{
				cout << i << "---" << j << " : " << score_order[i][j].first << "  type: " << score_order[i][j].second << endl;
			}
		}
		cout << "min: " << total_score << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << i << " -- " << order[i] << endl;
		}
		//out.close();
		cout << "total_score: " << total_score << endl;


		return total_score;
			
	}


	double Tiling_opt::com_tra_sim(vector<Point2f> &first_interval, vector<char>&first_char, vector<Point2f> &second_interval, vector<char>&second_char, int &flag) //compute trajectory similarity by DTW(Dynamic time warping)计算一对的相似得分
	{
		double score = 0;
		double score2 = 0;
		
		vector<Point2f> first_after_aff;
		vector<Point2f> second_after_aff;
		vector<Point2f> second_after_aff_ref;

		vector<char> second_char_out;
		//用仿射变化将两段轮廓对齐，消除初始位置的影响，其中warpAff_tra和warpAff_tra_ref_y都不影响对应曲率字符串的顺序

		double re_first = warpAff_tra(first_interval, first_after_aff);
		double re_second = warpAff_tra_sec(second_interval, second_after_aff, second_char, second_char_out);
		double re_second_ref = warpAff_tra_ref_y(second_interval, second_after_aff_ref);
		//int col = 800;//展示图片的尺寸
		//Mat drawing = Mat::zeros(col, col, CV_8UC3);
		//for (int i = 0; i < first_interval.size() - 1; i++)
		//{
		//	MyLine(drawing, first_interval[i], first_interval[i + 1], "red");
		//	MyLine(drawing, first_after_aff[i] + Point2f(col / 2 - re_first, col / 2), first_after_aff[i + 1] + Point2f(col / 2 - re_first, col / 2), "green");
		//}
		//for (int i = 0; i < second_interval.size() - 1; i++)
		//{
		//	MyLine(drawing, second_interval[i], second_interval[i + 1], "red");
		//	MyLine(drawing, second_after_aff[i] + Point2f(col / 2 - re_second, col / 2), second_after_aff[i + 1] + Point2f(col / 2 - re_second, col / 2), "lightblue");
		//}

		//imshow( "aff_sec: ", drawing);
		//
		//Mat drawing1 = Mat::zeros(col, col, CV_8UC3);
		//for (int i = 0; i < first_interval.size() - 1; i++)
		//{
		//	MyLine(drawing1, first_interval[i], first_interval[i + 1], "red");
		//	MyLine(drawing1, first_after_aff[i] + Point2f(col / 2 - re_first, col / 2), first_after_aff[i + 1] + Point2f(col / 2 - re_first, col / 2), "green");
		//}
		//for (int i = 0; i < second_interval.size() - 1; i++)
		//{
		//	MyLine(drawing1, second_interval[i], second_interval[i + 1], "red");
		//	MyLine(drawing1, second_after_aff_ref[i] + Point2f(col / 2 - re_second_ref, col / 2), second_after_aff_ref[i + 1] + Point2f(col / 2 - re_second_ref, col / 2), "lightblue");
		//}

		//imshow("aff_ref: ", drawing1);
		double length = (contour_length(first_interval) + contour_length(second_interval)) / 2;
		score = quadr_mismatch(first_after_aff, second_after_aff, first_char, second_char_out);
		score2 = quadr_mismatch(first_after_aff, second_after_aff, first_char, second_char);
		//score = DTW(first_after_aff, second_after_aff);
		//score2 = DTW(first_after_aff, second_after_aff_ref);
		//cout << "sec score: " << score << endl << "ref score: " << score2 << endl;
		if (score < score2)
		{
			flag = 0;
			return score;// / length;
		}
		else
		{
			flag = 1;
			return score2;// / length;
		}
	}


	double Tiling_opt::com_score_manual(string imagename1, string imagename2)
	{
		ofstream out("D:\\VisualStudioProjects\\manual_record\\manual_record.txt");
		out << "0 代表sec,1 代表ref" << endl;
		prototile_first->loadTileData(imagename1);
		prototile_second->loadTileData(imagename2);

		vector<Point2f> contour1 = prototile_first->contour_sample[0];
		vector<Point2f> contour2 = prototile_second->contour_sample[0];

		//vector<double> contour1_cur = prototile_first->contour_curva[0];
		//vector<double> contour2_cur = prototile_second->contour_curva[0];
		

		double factor = com_scale_factor();
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i].x = contour2[i].x * factor;
			contour2[i].y = contour2[i].y * factor;
		}
		/*double scale = 0.5;
		for (int i = 0; i < contour1.size(); i++)
		{
		contour1[i].x = contour1[i].x * scale;
		contour1[i].y = contour1[i].y * scale;
		}
		for (int i = 0; i < contour2.size(); i++)
		{
		contour2[i].x = contour2[i].x * scale;
		contour2[i].y = contour2[i].y * scale;
		}*/
		int parem[2][4];
		parem[0][0] = 21;
		parem[0][1] = 15;
		parem[0][2] = 21;
		parem[0][3] = 0;
		int all_num = 0;
		vector<vector<Point2f>> proto_interval_first;
		vector<vector<char>> proto_first_char;
		vector<Point2f> each_inter;
		vector<char> each_cur_char;
		for (int j = 0; j < parem[0][0]; j++)
		{
			each_inter.push_back(contour1[all_num + parem[0][3]]);
			each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
			all_num++;
		}
		each_inter.push_back(contour1[all_num + parem[0][3]]);
		each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
		proto_interval_first.push_back(each_inter);
		proto_first_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());

		for (int j = 0; j < parem[0][1]; j++)
		{
			each_inter.push_back(contour1[all_num + parem[0][3]]);
			each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
			all_num++;
		}
		each_inter.push_back(contour1[all_num + parem[0][3]]);
		each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
		proto_interval_first.push_back(each_inter);
		proto_first_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());


		for (int j = 0; j < parem[0][2]; j++)
		{
			each_inter.push_back(contour1[all_num + parem[0][3]]);
			each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
			all_num++;
		}
		each_inter.push_back(contour1[all_num + parem[0][3]]);
		each_cur_char.push_back(prototile_first->cur_string[0][all_num + parem[0][3]]);
		proto_interval_first.push_back(each_inter);
		proto_first_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());


		for (; all_num < contour1.size(); all_num++)
		{
			each_inter.push_back(contour1[(all_num + parem[0][3]) % (contour1.size())]);
			each_cur_char.push_back(prototile_first->cur_string[0][(all_num + parem[0][3]) % (contour1.size())]);
		}
		each_inter.push_back(contour1[(all_num + parem[0][3]) % (contour1.size())]);
		each_cur_char.push_back(prototile_first->cur_string[0][(all_num + parem[0][3]) % (contour1.size())]);
		proto_interval_first.push_back(each_inter);
		proto_first_char.push_back(each_cur_char);

		if (out.is_open())
		{
			out << "first image " << imagename1 << "  second image " << imagename2 << endl;
			out << "first num1: " << parem[0][0] << " num2: " << parem[0][1] << " num3: " << parem[0][2] <<
				" all num: " << all_num << " delay: " << parem[0][3] << endl;
		}
		parem[1][0] = 20;
		parem[1][1] = 15;
		parem[1][2] = 24;
		parem[1][3] = 9;
		all_num = 0;
		vector<vector<Point2f>> proto_interval_second;
		vector<vector<char>> proto_second_char;
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());
		for (int j = 0; j < parem[1][0]; j++)
		{
			each_inter.push_back(contour2[all_num + parem[1][3]]);
			each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
			all_num++;
		}
		each_inter.push_back(contour2[all_num + parem[1][3]]);
		each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
		proto_interval_second.push_back(each_inter);
		proto_second_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());

		for (int j = 0; j < parem[1][1]; j++)
		{
			each_inter.push_back(contour2[all_num + parem[1][3]]);
			each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
			all_num++;
		}
		each_inter.push_back(contour2[all_num + parem[1][3]]);
		each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
		proto_interval_second.push_back(each_inter);
		proto_second_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());

		for (int j = 0; j < parem[1][2]; j++)
		{
			each_inter.push_back(contour2[all_num + parem[1][3]]);
			each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
			all_num++;
		}
		each_inter.push_back(contour2[all_num + parem[1][3]]);
		each_cur_char.push_back(prototile_second->cur_string[0][all_num + parem[1][3]]);
		proto_interval_second.push_back(each_inter);
		proto_second_char.push_back(each_cur_char);
		each_inter.swap(vector<Point2f>());
		each_cur_char.swap(vector<char>());

		for (; all_num < contour2.size(); all_num++)
		{
			each_inter.push_back(contour2[(all_num + parem[1][3]) % (contour2.size())]);
			each_cur_char.push_back(prototile_second->cur_string[0][(all_num + parem[1][3]) % (contour2.size())]);
		}
		each_inter.push_back(contour2[(all_num + parem[1][3]) % (contour2.size())]);
		each_cur_char.push_back(prototile_second->cur_string[0][(all_num + parem[1][3]) % (contour2.size())]);
		proto_interval_second.push_back(each_inter);
		proto_second_char.push_back(each_cur_char);

		//int flag1;
		//double score = com_tra_sim(proto_interval_first[2], proto_interval_second[1], flag1);
		//cout << score << endl;

		//检查字符串与曲率个数是否对应
		//for (int i = 0; i < 4; i++)
		//{
		//	cout << "proto_interval_first: " << proto_interval_first[i].size() << "   ------   proto_first_char: " << proto_first_char[i].size() << endl;
		//	cout << "proto_interval_second: " << proto_interval_second[i].size() << "   ------   proto_second_char: " << proto_second_char[i].size() << endl;
		//}

		if (out.is_open())
		{
			out << "second num1: " << parem[1][0] << " num2: " << parem[1][1] << " num3: " << parem[1][2] <<
				" all num: " << all_num << " delay: " << parem[1][3] << endl;

		}
		out.close();
		vector<pair<int, int>> order;

		double score = com_optimal_score(proto_interval_first, proto_first_char,proto_interval_second, proto_second_char, order); //返回一种排列顺序以及分数并依次记录，order[i].first：轮廓1中第i段对应的2中的段数，order[i].second:匹配类型是反转还是正常
		
		//求得最优结果后重新计算得到排列方式，这里仍有一个问题：上一步com_optimal里的仍是通过仿射变化一段某一端点对齐的结果，下面展示的时候则是把对应边中点对齐
		vector<Point2f> output_first;
		vector<Point2f> output_second;
		vector<Point2f> output_final;
		vector<char> sec_char_out;
		//vector<Point2f> output_second_ref;
		
		
		//在这里进行变形已经还原，可以将第一层for换成while加终止条件
		for (int i = 0; i < 4; i++)
		{
			//all_types[max_x][max_y].second[i];
			Point2f extre[2][2];
			extre[0][0] = proto_interval_first[i][0];
			extre[0][1] = proto_interval_first[i][proto_interval_first[i].size()-1];
			extre[1][0] = proto_interval_second[order[i].first][0];
			extre[1][1] = proto_interval_second[order[i].first][proto_interval_second[order[i].first].size()-1];
			Point2f sae[2][2];

			double re_first = warpAff_tra(proto_interval_first[i], output_first);
			double re_second = 0;
			if (order[i].second == 0)
			{
				re_second = warpAff_tra_sec(proto_interval_second[order[i].first], output_second, proto_second_char[order[i].first], sec_char_out);
			}
			else if (order[i].second == 1)
			{
				re_second = warpAff_tra_ref_y(proto_interval_second[order[i].first], output_second);
			}

			//double re_second_ref = warpAff_tra_ref_y(proto_interval_second[order[i]], output_second_ref);

			int size = 800;//展示图片的尺寸

			Point2f shift1 = Point2f(size / 2 - re_first, size / 2);
			Point2f shift2 = Point2f(size / 2 - re_second, size / 2);
			//Mat drawing = Mat::zeros(size, size, CV_8UC3);  //黑色背景
			Mat drawing_src1 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));  // 白色背景
			Mat drawing_src2 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));
			Mat drawing_src3 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));
			Mat drawing_dst = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));
			for (int t = 0; t < output_first.size(); t++)
			{
				output_first[t] = output_first[t] + shift1;
			}
			for (int t = 0; t < output_second.size(); t++)
			{
				output_second[t] = output_second[t] + shift2;
			}
			for (int t = 0; t < output_first.size() - 1; t++)
			{
				MyLine(drawing_src1, output_first[t], output_first[t + 1], "green");
				MyLine(drawing_src3, output_first[t], output_first[t + 1], "green");  //总的展示
			}
			for (int t = 0; t < output_second.size() - 1; t++)
			{
				//	MyLine(drawing, first_interval[i], first_interval[i + 1], 1);
				MyLine(drawing_src2, output_second[t], output_second[t + 1], "orange");
				MyLine(drawing_src3, output_second[t], output_second[t + 1], "orange");  //总的展示
			}
			
			//重新计算
			DTW(output_first, output_second);
			
			for (int i = 0; i < dp_path.size(); i++)
			{
				MyLine(drawing_src3, output_first[dp_path[i].first], output_second[dp_path[i].second], "grey"); 
			}

			//---------------------以下为morphing阶段
			// 首先找到一一对应的点序列
			vector<Point2f> src1_points;
			vector<Point2f> src2_points;
			src1_points.push_back(output_first[0]);
			src2_points.push_back(output_second[0]);
			if (output_first.size() < output_second.size())
			{	
				int f = 1;
				int j = 0;
				for (int i = 1; i < output_first.size() - 1;)
				{
			
					double min = 10000;
					while (f < dp_path.size())
					{

						if ((dp_path[f].first < i) || dp_path[f].second == j) //检查下一条路径的另一端点是否已被使用
						{
							f++;
							if (dp_path[f].first > i) i++;
							continue;
						}
						if (dp_path[f].first == i)
						{
							if (length_two_point2f(output_first[dp_path[f].first], output_second[dp_path[f].second]) < min)
							{
								min = length_two_point2f(output_first[dp_path[f].first], output_second[dp_path[f].second]);
								j = dp_path[f].second;
							}
							f++;
						}
						if (dp_path[f].first > i)
						{
							src1_points.push_back(output_first[i]);
							src2_points.push_back(output_second[j]);
							i++;
							j++;
							break;
						}

					}
				}
			}
			else
			{			
				int f = 1;
				int j = 0;
				for (int i = 1; i < output_second.size() - 1; )
				{
					double min = 10000;

					while (f < dp_path.size())
					{

						if ((dp_path[f].second < i) || dp_path[f].first < j) 
						{
							
							f++;
							if (dp_path[f].second > i) i++;
							continue;
						}
						if (dp_path[f].second == i)
						{
							if (length_two_point2f(output_first[dp_path[f].first], output_second[dp_path[f].second]) < min)
							{
								min = length_two_point2f(output_first[dp_path[f].first], output_second[dp_path[f].second]);
								j = dp_path[f].first;
							}
							f++;
						}
						if (dp_path[f].second > i)
						{
							src1_points.push_back(output_first[j]);
							src2_points.push_back(output_second[i]);
							i++;
							j++;
							break;
						}

					}
				}
			}
			src1_points.push_back(output_first[dp_path[dp_path.size() - 1].first]);
			src2_points.push_back(output_second[dp_path[dp_path.size() - 1].second]);
			
			cout << endl;
			if (src1_points.size() == src2_points.size())
			{
				for (int i = 0; i < src1_points.size(); i++)
				{
					circle(drawing_src3, src1_points[i], 1, Scalar(255, 0, 0), -1);
					circle(drawing_src3, src2_points[i], 1, Scalar(0, 0, 255), -1);
					MyLine(drawing_src3, src1_points[i], src2_points[i], "orange");
					cout << src1_points[i] << " - " << src2_points[i] << endl;
				}
			}
			else cout << "num error!" << endl;

			ImageMorphing(drawing_src1, src1_points, drawing_src2, src2_points, drawing_dst, output_final, 0.5);
			//char i;
			double length_final = length_two_point2f(output_final[0], output_final[output_final.size()-1]);
			for (int i = 0; i < output_final.size()-1; i++)
			{
				MyLine(drawing_src3, output_final[i], output_final[i+1], "green");
			}
			//imshow("drawing_src1", drawing_src1);
			//imshow("drawing_src2", drawing_src2);
			//imshow("drawing_dst", drawing_dst);

			string name = "the ";
			name = name + char(i + 48) + " pair: ";
			imshow(name, drawing_src3);

			//在这里将out_final映射回原来的位置，根据两个端点的坐标
			Point2f mid1;
			mid1.x = (extre[0][0].x + extre[0][1].x) / 2;
			mid1.y = (extre[0][0].y + extre[0][1].y) / 2;
			sae[0][0] = mid1 + (length_final / 2) * unit_vec(extre[0][0] - mid1);
			sae[0][1] = mid1 + (length_final / 2) * unit_vec(extre[0][1] - mid1);
			Point2f mid2;
			mid2.x = (extre[1][0].x + extre[1][1].x) / 2;
			mid2.y = (extre[1][0].y + extre[1][1].y) / 2;
			sae[1][0] = mid2 + (length_final / 2) * unit_vec(extre[1][0] - mid2);
			sae[1][1] = mid2 + (length_final / 2) * unit_vec(extre[1][1] - mid2);
			
			vector<Point2f> out1;
			re_warp_Aff(output_final, out1, sae[0][0], sae[0][1]);
			proto_interval_first[i].swap(out1);
			
			//找到哪一条边的哪一个点是是改变后的sae点

		}

		int row = 800;
		int col = 1600;
		//Mat drawing4 = Mat::zeros(row, col, CV_8UC3); //黑色背景
		Mat drawing4 = Mat(row, col, CV_8UC3, Scalar(255, 255, 255)); //白色背景
		//Mat drawing5 = Mat::zeros(800, 800, CV_8UC3);
		//namedWindow("simple111 contour", WINDOW_AUTOSIZE);

		//用shift将两张图对齐到图像的左右中点
		Point2f shift1 = Point2f(prototile_first->center_point.x - col / 4, prototile_first->center_point.y - row / 2);

		Point2f shift2 = prototile_second->center_point - Point2f(3 * col / 4, row / 2);

		/*for (int s = 0; s < 4; s++)
		{
		Point2f srcTri[3];
		Point2f dstTri[3];
		srcTri[0] = proto_interval_first[s][0];
		srcTri[1] = proto_interval_first[s][proto_interval_first.size() - 1];
		double midd = (srcTri[1].x + srcTri[0].x) / 2;
		Point2f midp = Point2f((srcTri[1].x + srcTri[0].x) / 2, (srcTri[1].y + srcTri[0].y) / 2;);
		srcTri[2] = Point2f(midd,)
		dstTri[0] = Point2f(0, 0);//Point2f(50, 50);
		//cout << "dstTri[0]: " << dstTri[0] << endl;
		dstTri[1] = Point2f(length_two_point2f(srcTri[0], srcTri[1]), 0);
		for (int i = 0; i < 4; i++)
		{
		//cout << "proto_interval_first[" << i << "]  num:" << proto_interval_first[i].size() << endl;
		string color;
		if (i == 0) color = "red";
		else if (i == 1) color = "blue";
		else if (i == 2) color = "green";
		else if (i == 3) color = "cyan";

		if (s == 0)
		{
		for (int j = 0; j < proto_interval_first[i].size() - 1; j++)
		{
		MyLine(drawing4, proto_interval_first[i][j] - shift1, proto_interval_first[i][j + 1] - shift1, color);

		}
		MyLine(drawing4, proto_interval_first[i][0] - shift1, proto_interval_first[i][proto_interval_first[i].size() - 1] - shift1, "orange");
		}
		for (int j = 0; j < proto_interval_second[order[i].first].size() - 1; j++)
		{
		MyLine(drawing4, proto_interval_second[order[i].first][j] - shift2, proto_interval_second[order[i].first][j + 1] - shift2, color);

		}
		MyLine(drawing4, proto_interval_second[i][0] - shift2, proto_interval_second[i][proto_interval_second[i].size() - 1] - shift2, "orange");
		}
		}*/

		for (int i = 0; i < 4; i++)
		{
			//cout << "proto_interval_first[" << i << "]  num:" << proto_interval_first[i].size() << endl;
			string color;
			if (i == 0) color = "red";
			else if (i == 1) color = "blue";
			else if (i == 2) color = "green";
			else if (i == 3) color = "cyan";
			Point2f srcTri[3];
			Point2f dstTri[3];
			for (int j = 0; j < proto_interval_first[i].size() - 1; j++)
			{
				MyLine(drawing4, proto_interval_first[i][j] - shift1, proto_interval_first[i][j + 1] - shift1, color);

			}
			MyLine(drawing4, proto_interval_first[i][0] - shift1, proto_interval_first[i][proto_interval_first[i].size() - 1] - shift1, "orange");
			for (int j = 0; j < proto_interval_second[order[i].first].size() - 1; j++)
			{
				MyLine(drawing4, proto_interval_second[order[i].first][j] - shift2, proto_interval_second[order[i].first][j + 1] - shift2, color);

			}
			MyLine(drawing4, proto_interval_second[i][0] - shift2, proto_interval_second[i][proto_interval_second[i].size() - 1] - shift2, "orange");
		}
		imshow("first contour: " + imagename1 + " intervals", drawing4);


		return 0;
	}

}
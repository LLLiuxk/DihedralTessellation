/* 1. 求得轮廓，采样，计算
2. 求得凹凸性进行划分
3. 将图案进行组装，从dataset中寻找与间隙相似度高的图案
4. 变形和decorate
*/
#include "tilingOpt.h"

//static int N = 800;
namespace Tiling_tiles{

	Tiling_opt::Tiling_opt()
	{
		prototile_first = new Prototile();
		prototile_mid = new Prototile();
		prototile_second = new Prototile();
		prototile_tem = new Prototile();
		all_types = 500;
		sampling_num = 3;
		//memset(dp, 0, sizeof(dp));
		//memset(dp_inver, 0, sizeof(dp_inver));

	}

	Tiling_opt::~Tiling_opt()
	{
	}

	void Tiling_opt::Tiling_clear()
	{
		memset(dis, 0, sizeof(dis));
		memset(dis_cur, 0, sizeof(dis_cur));
		memset(distance, 0, sizeof(distance));
		memset(step, 0, sizeof(step));
		prototile_first->Pro_clear();
		prototile_mid->Pro_clear();
		prototile_second->Pro_clear();
		prototile_tem->Pro_clear();
		contour_dataset.swap(vector<vector<Point2f>>());
	}

	void Tiling_opt::load_dataset()
	{
		if (!contour_dataset.empty()) contour_dataset.swap(vector<vector<Point2f>>());
		for (int i = 0; i <= all_types; i++)
		{
			prototile_second->Pro_clear();
			prototile_second->setpath();
			string image = int2string(i);
			vector<Point2f> data_;
			prototile_second->setname(image);

			//imgtocout();
			data_ = prototile_second->readTxt();
			//cout << image << endl;	
			if (!data_.empty())
				contour_dataset.push_back(data_);
		}
		cout << "load over!" << endl;
	}

	void Tiling_opt::com_all_TARs(int num_c)
	{
		if (!all_con_tars.empty()) 
		{
			cout << "the TARs have been exist" << endl;
			return;
		}
		int all_num = contour_dataset.size();
		cout << "all_num?" << all_num << endl;
		for (int i = 0; i < all_num; i++)
		{
			prototile_second->Pro_clear();
			prototile_second->loadPoints(contour_dataset[i]);
			double shape_com;
			vector<vector<double>> tar_all = prototile_second->compute_TAR(prototile_second->contour_sample[num_c], shape_com);//(num_c+1)*100 points
			vector<vector<double>> tar_all_flip = prototile_second->compute_TAR(prototile_second->contour_sample_flip[num_c], shape_com);
			////将特征点的tar值记录下来用作快速筛选算法
			//vector<vector<double>> tar_fea;
			//vector<vector<double>> tar_fea_flip;
			//vector<int> cand_points_index = feature_points(prototile_second->contour_sample[num_c], 1, 3, cos(PI * 160 / 180));
			//for (int j = 0; j < cand_points_index.size(); j++)
			//{
			//	tar_fea.push_back(tar_all[cand_points_index[j]]);
			//}
			//all_fea_tars.push_back(tar_fea);
			//cand_points_index = feature_points(prototile_second->contour_sample_flip[num_c], 1, 3, cos(PI * 160 / 180));
			//for (int j = 0; j < cand_points_index.size(); j++)
			//{
			//	tar_fea_flip.push_back(tar_all_flip[cand_points_index[j]]);
			//}
			//all_fea_tars_flip.push_back(tar_fea_flip);

			vector<vector<double>> tar_fea;
			vector<vector<double>> tar_fea_flip;
			vector<int> cand_points_index = most_convex_p(prototile_second->contour_sample[num_c], curvature_com(prototile_second->contour_sample[num_c]),30);
			
			for (int j = 0; j < cand_points_index.size(); j++)
			{
				tar_fea.push_back(tar_all[cand_points_index[j]]);
			}
			all_fea_tars.push_back(tar_fea);
			cand_points_index = most_convex_p(prototile_second->contour_sample_flip[num_c], curvature_com(prototile_second->contour_sample_flip[num_c]), 30);
			for (int j = 0; j < cand_points_index.size(); j++)
			{
				tar_fea_flip.push_back(tar_all_flip[cand_points_index[j]]);
			}
			all_fea_tars_flip.push_back(tar_fea_flip);
			
			all_con_tars.push_back(tar_all);
			all_con_tars_flip.push_back(tar_all_flip);
			all_shape_complexity.push_back(shape_com);

			//double shape_com1;
			//vector<vector<double>> tar_all1 = prototile_second->compute_TAR(prototile_second->contour_sample[0], shape_com1);//(num_c+1)*100 points
			//vector<vector<double>> tar_all_flip1 = prototile_second->compute_TAR(prototile_second->contour_sample_flip[0], shape_com1);
		}
		cout << "all TARs of contours have been computed" << endl;
	}

	void Tiling_opt::check_Repetitive_pattern()
	{
		for (int i = 0; i <= all_types; i++)
		{
			prototile_second->Pro_clear();
			prototile_second->loadPoints(contour_dataset[i]);
			vector<Point2f> contour_sec = prototile_second->contour_sample[2];
			double score = 10000;
			int samey = 0;
			for (int j = 0; j < all_types; j++)
			{
				if (j == i) continue;
				prototile_tem->Pro_clear();
				prototile_tem->loadPoints(contour_dataset[j]);
				vector<Point2f> contour_tem = prototile_tem->contour_sample[2];
				double score1 = matchShapes(contour_sec, contour_tem, 1, 0);
				if (score1 < score)
				{
					score = score1;
					samey = j;
				}
			}
			if (score < 0.02) cout << i << " Probably similar with " << samey << " : " << score << endl;

		}
	}


	void Tiling_opt::tiliing_generation(string imaname)
	{
		clock_t start, midtime, finish;
		start = clock();
		bool check_self_intersect = true;
		int num_c = 1;//选择(num_c+1)*100个点
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		cout << "p_p_index: " << p_p_index.size() << endl;
		vector<Point2f> contour_ = prototile_first->contour;
		vector<Point2f> cont_orig = prototile_first->contour_sample[1];
		vector<Point2f> cont_rota = prototile_first->contour_sample[4]; //旋转暂时用500个点
		int contsize = contour_.size();
		Point2f cent_cont = center_p(contour_);
		//vector<Point2f> contour_ = prototile_first->contour;
		load_dataset();
		com_all_TARs(num_c);
		string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result300\\" + prototile_first->contourname;
		//string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result_tuning\\" + prototile_first->contourname;
		const char *na = rootname.c_str();
		if (_access(na, 0) != -1)

		{
			printf("The  file/dir had been Exisit ");
			return;
		}
		mkdir(na);

		int trans = Tanslation_rule(p_p_index, contour_, rootname);
		int rotas = Rotation_rule(p_p_index, contour_, rootname);
		int flips = Flipping_rule(p_p_index, contour_, rootname);
		int count = trans + rotas + flips;
		cout << "succeed count: " << count << " trans: " << trans << " rotat: " << rotas << " flips: " << flips << endl;
		//midtime = clock();
		//cout << "Time consumption: " << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		if (count == 0)
		{
			cout << "no right placement" << endl;
			return;
		}
	}

	int Tiling_opt::Tanslation_rule(vector<int> part_points_index, vector<Point2f> &contour_s, string rootname)
	{
		int trans = 0;
		int ppindex = part_points_index.size();
		int margin = contour_s.size() / 20;
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(contour_s);

		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				if (abs(part_points_index[j] - part_points_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(part_points_index[m] - part_points_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						//cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
						if (abs(part_points_index[n] - part_points_index[m]) < margin) continue;
						vector<Point2f> inner_contour;
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						vector<int> result_index; //1.translation：以下所有该类摆放都以1-3,2-4为轴摆放
						result_index.push_back(part_points_index[i]);
						result_index.push_back(part_points_index[j]);
						result_index.push_back(part_points_index[m]);
						result_index.push_back(part_points_index[n]);
						Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (!translation_placement(result_index, contour_s, inner_contour, mid_interval, drawing1))
						{
							cout << ++trans << " Translation succeed" << endl;
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 0;
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(400, 400) - cent_cont;
							for (int jj = 0; jj < contsize; jj++)
							{
								circle(drawing1, contour_s[jj] + shift2, 1, Scalar(0, 0, 0), -1);

								//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
							}
							for (int jj = 0; jj < ppindex; jj++)
							{
								circle(drawing1, contour_s[part_points_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
							}
							for (int jj = 0; jj < 4; jj++)
							{
								circle(drawing1, contour_s[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
							}

							string filename = rootname + "\\" + int2string(trans - 1) + "transPlacingResult.png";
							imwrite(filename, drawing1);
							//Mat draw = drawing1(Rect(800, 0, 800, 800));
							//all_tiling_Mat.push_back(draw);
						}
					}
				}
			}
		}
		return trans;
	}

	int Tiling_opt::Rotation_rule(vector<int> part_points_index, vector<Point2f> &contour_s, string rootname)
	{
		int rotas = 0;
		int ppindex = part_points_index.size();
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(contour_s);
		vector<pair<Point2f, int>> insert_points;
		vector<vector<Point2f>> all_result_points;
		vector<vector<int>> all_result_index; //通过重新定位确定点的下标序列
		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				vector<Point2f> inner_contour;
				vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
				vector<int> mark_13;
				mark_13.push_back(part_points_index[i]);
				mark_13.push_back(part_points_index[j]);

				vector<vector<Point2f>> all_result_p = find_rota_tilingV(contour_s, mark_13, insert_points);
				int allresultsize = all_result_p.size();
				for (int num = 0; num < allresultsize; num++)
				{
					all_result_points.push_back(all_result_p[num]);
				}
			}
		}
		//将新点插入轮廓点集，先对insert_points进行排序，然后从后往前排，因为insert里记录的也是序列，如果从前往后就会出错

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


		for (int num = 0; num < allresultsize; num++)
		{
			Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
			if (!rotation_placement(all_result_index[num], contour_s, inner_contour, mid_interval, drawing1))
			{
				cout << ++rotas << " Rotation succeed" << endl;
				inPat one_situation;
				one_situation.in_contour = inner_contour;
				one_situation.in_interval = mid_interval;
				one_situation.type = 1;
				all_inner_conts.push_back(one_situation);
				inner_contour.swap(vector<Point2f>());
				mid_interval.swap(vector<int>());

				Point2f shift2 = Point2f(400, 400) - cent_cont;
				for (int jj = 0; jj < contsize; jj++)
				{
					circle(drawing1, contour_s[jj] + shift2, 1, Scalar(0, 0, 0), -1);

					//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
				}
				for (int jj = 0; jj < ppindex; jj++)
				{
					circle(drawing1, contour_s[part_points_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
				}
				for (int jj = 0; jj < 4; jj++)
				{
					circle(drawing1, contour_s[all_result_index[num][jj]] + shift2, 8, Scalar(0, 255, 0), -1);
				}
				string filename = rootname + "\\" + int2string(rotas - 1) + "rota_PlacingResult.png";
				imwrite(filename, drawing1);
				//Mat draw = drawing1(Rect(800, 0, 800, 800));
				//all_tiling_Mat.push_back(draw);
			}
		}
		return rotas;

	}

	vector<vector<Point2f>> Tiling_opt::find_rota_tilingV(vector<Point2f>&cont, vector<int> mark_13, vector<pair<Point2f, int>> &all_insert_points)
	{
		vector<vector<Point2f>> all_result_points;
		int contsize = cont.size();
		double line_length = arcLength(cont, 1) / 2;
		double margin = 2 * line_length / contsize;
		Point2f center = 0.5 * cont[mark_13[0]] + 0.5 * cont[mark_13[1]];
		Point2f unitv = unit_vec(cont[mark_13[1]] - cont[mark_13[0]]);
		Point2f st = center - line_length*unitv;
		Point2f end = center + line_length*unitv;
		vector<Point2f> line1;
		line1.push_back(st);
		line1.push_back(end);
		vector<Point2f> line24;
		Mat rot_mat(2, 3, CV_32FC1);

		for (int angle = 1; angle < 360; angle++)
		{
			rot_mat = getRotationMatrix2D(center, angle, 1);
			transform(line1, line24, rot_mat);
			vector<Point2f> left_sec;
			vector<int> left_sec_index; 
			vector<Point2f> right_sec; 
			vector<int> right_sec_index;
			//首先计算两侧所有的交点，然后逐一配对，检查长度
			for (int i = mark_13[0]; i < mark_13[1]; i++)
			{
				Point2f sec_p;
				int line_sec_p = line_intersection(Line_Seg(line24[0], line24[1]), Line_Seg(cont[i], cont[i + 1]), sec_p);
				if (line_sec_p == 1)
				{
					left_sec.push_back(sec_p);
					left_sec_index.push_back(i);
				}
			}
			for (int i = mark_13[1]; i < mark_13[0] + contsize; i++)
			{
				Point2f sec_p;
				int line_sec_p = line_intersection(Line_Seg(line24[0], line24[1]), Line_Seg(cont[i % contsize], cont[(i + 1) % contsize]), sec_p);
				if (line_sec_p == 1)
				{
					right_sec.push_back(sec_p);
					right_sec_index.push_back(i % contsize);
				}
			}
			//接下来需要分别判断交点属于原有的点还是需要插值的点
			//首先判断长度是否合适
			for (int j = 0; j < left_sec.size(); j++)
				for (int t = 0; t < right_sec.size(); t++)
				{
					double leng1 = length_two_point2f(left_sec[j], center);
					double leng2 = length_two_point2f(right_sec[t], center);
					if (abs(leng1 - leng2) < 0.01)  //长度符合要求，进行下一步处理
					{
						vector<Point2f> result_points;   //记录四个点的坐标，存储到all_result_points中
						result_points.push_back(cont[mark_13[0]]);
						//判断左边点是否是原有的点
						if (length_two_point2f(left_sec[j], cont[left_sec_index[j]]) < 0.01)//可以看做是原有的点
						{
							result_points.push_back(cont[left_sec_index[j]]);
						}
						else if (length_two_point2f(left_sec[j], cont[left_sec_index[j] + 1]) < 0.01)//可以看做是原有的点
						{
							result_points.push_back(cont[left_sec_index[j] + 1]);
						}
						else  //说明是插入的新点，存到all_insert_points
						{
							result_points.push_back(left_sec[j]);
							all_insert_points.push_back(make_pair(left_sec[j], left_sec_index[j]));
						}
						result_points.push_back(cont[mark_13[1]]);
						//判断右边点是否是原有的点
						if (length_two_point2f(right_sec[t], cont[right_sec_index[t]]) < 0.01) //可以看做是原有的点
						{
							result_points.push_back(cont[right_sec_index[t]]);
						}
						else if (length_two_point2f(right_sec[t], cont[right_sec_index[t] + 1]) < 0.01) //可以看做是原有的点
						{
							result_points.push_back(cont[right_sec_index[t] + 1]);
						}
						else
						{
							result_points.push_back(right_sec[t]);
							all_insert_points.push_back(make_pair(right_sec[t], right_sec_index[t]));
						}
						all_result_points.push_back(result_points);

					}
				}
			//计算在mark参数下所有的平行四边形，将左边全部返回

		}
		cout << "rotation all_result_points num:" << all_result_points.size() << endl;
		return all_result_points;
	}


	int Tiling_opt::Flipping_rule(vector<int> part_points_index, vector<Point2f> &contour_s, string rootname)
	{
		int flips = 0;
		int ppindex = part_points_index.size();
		int margin = contour_s.size() / 20;
		int contsize = contour_s.size();
		Point2f cent_cont = center_p(contour_s);

		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				if (abs(part_points_index[j] - part_points_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(part_points_index[m] - part_points_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						//cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
						if (abs(part_points_index[n] - part_points_index[m]) < margin) continue;
						vector<Point2f> inner_contour;
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						vector<int> result_index; //1.translation：以下所有该类摆放都以1-3,2-4为轴摆放
						result_index.push_back(part_points_index[i]);
						result_index.push_back(part_points_index[j]);
						result_index.push_back(part_points_index[m]);
						result_index.push_back(part_points_index[n]);
						Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						Mat drawing2 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));

						if (!flipping_placement(result_index, contour_s, inner_contour, mid_interval, drawing1, 0))
						{

							cout << ++flips << " Flipping succeed" << endl;
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 2;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());

							Point2f shift2 = Point2f(400, 400) - cent_cont;
							for (int jj = 0; jj < contsize; jj++)
							{
								circle(drawing1, contour_s[jj] + shift2, 1, Scalar(0, 0, 0), -1);

								//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
							}
							for (int jj = 0; jj < ppindex; jj++)
							{
								circle(drawing1, contour_s[part_points_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
							}
							for (int jj = 0; jj < 4; jj++)
							{
								circle(drawing1, contour_s[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
							}

							string filename = rootname + "\\" + int2string(flips - 1) + "flip(1-3)PlacingResult.png";
							imwrite(filename, drawing1);
							//Mat draw = drawing1(Rect(800, 0, 800, 800));
							//all_tiling_Mat.push_back(draw);
						}
						if (!flipping_placement(result_index, contour_s, inner_contour, mid_interval, drawing2, 1))
						{
							cout << ++flips << " Flipping succeed" << endl;
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 3;
							all_inner_conts.push_back(one_situation);

							Point2f shift2 = Point2f(400, 400) - cent_cont;
							for (int jj = 0; jj < contsize; jj++)
							{
								circle(drawing2, contour_s[jj] + shift2, 1, Scalar(0, 0, 0), -1);
								//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
							}
							for (int jj = 0; jj < ppindex; jj++)
							{
								circle(drawing2, contour_s[part_points_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
							}
							for (int jj = 0; jj < 4; jj++)
							{
								circle(drawing2, contour_s[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
							}

							string filename = rootname + "\\" + int2string(flips - 1) + "flip(2-4)PlacingResult.png";
							imwrite(filename, drawing2);
							//Mat draw = drawing2(Rect(800, 0, 800, 800));
							//all_tiling_Mat.push_back(draw);
						}
					}
				}
			}
		}
		return flips;
	}

	void Tiling_opt::match_candidate(int Tiling_index)
	{
		if (all_inner_conts.empty())
		{
			cout << "No correct middle patterns!" << endl;
			return;
		}
		int num_c = 1;//采样点数
		vector<pair<int, bool>> candidate_patterns;
		candidate_patterns = compare_choose_TAR(all_inner_conts[Tiling_index].in_contour);

		mid_inter = joint_relocate(all_inner_conts[Tiling_index].in_contour, all_inner_conts[Tiling_index].in_interval, num_c);

		vector<Point2f> contour_inner = prototile_mid->contour_sample[1];
		double sc_inner = 0;
		vector<vector<double>> inner_tar = prototile_mid->compute_TAR(contour_inner, sc_inner);
		candidate_contours.swap(vector<vector<Point2f>>());
		cand_paths.swap(vector<vector<pair<int, int>>>());

		for (int j = 0; j < 8; j++) //只要候选图案里的前8个
		{
			//将所有的结果保存下来
			prototile_second->loadPoints(contour_dataset[candidate_patterns[j].first]);

			vector<Point2f> contour_cand;// = prototile_second->contour_sample[1];
			vector<vector<double>> cand_tar;
			if (candidate_patterns[j].second)
			{
				//cout << "it is flip" << endl;
				cand_tar = all_con_tars_flip[candidate_patterns[j].first];
				contour_cand = prototile_second->contour_sample_flip[1];
			}
			else
			{
				//cout << "it is not flip" << endl;
				cand_tar = all_con_tars[candidate_patterns[j].first];
				contour_cand = prototile_second->contour_sample[1];
			}
			vector<pair<int, int>> path;
			int shift = 0;
			int width = 6;
			double re = tar_mismatch(inner_tar, cand_tar, path, shift, width);

			vector<Point2f> contour_mid;
			int c2size = contour_cand.size();
			for (int i = shift; i < shift + c2size; i++)
			{
				contour_mid.push_back(contour_cand[i % c2size]);
			}

			cout << contour_inner.size() << "   " << c2size << "  c3: " << contour_mid.size() << endl;
			contour_cand = contour_mid;
			double scale = arcLength(contour_inner, true) / arcLength(contour_cand, true);
			//cout << "scale: " << scale << endl;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] * scale;
			}
			Point2f cen1 = center_p(contour_inner);
			Point2f cen2 = center_p(contour_cand);
			Point2f shift2 = cen1 - cen2;
			for (int i = 0; i < contour_cand.size(); i++)
			{
				contour_cand[i] = contour_cand[i] + shift2;
			}
			Point2f v0 = contour_inner[0] - cen1;
			Point2f v1 = contour_cand[0] - cen1;
			double angle = acos(cos_two_vector(v0, v1)) / PI * 180;
			if (sin_two_vector(v0, v1) < 0) angle = -angle;
			Mat rot_mat(2, 3, CV_32FC1);
			rot_mat = getRotationMatrix2D(cen1, angle, 1);
			transform(contour_cand, contour_mid, rot_mat);
			contour_cand = contour_mid;

			candidate_contours.push_back(contour_cand);
			cand_paths.push_back(path);
		}
	}

	vector<int> Tiling_opt::morphed_results(vector<Point2f> &morphed_A, int Candidate_index, int Tiling_index)
	{
		cout << "cand_paths num: " << cand_paths.size() << endl;
		vector<Point2f> final_pettern;
		vector<int> mid_inter_new;
		vector<pair<int, int>> path = cand_paths[Candidate_index];

		int psize = path.size();
		cout << "path size: " << psize << endl;
		vector<Point2f> contour1 = prototile_mid->contour_sample[1];
		vector<Point2f> contour2 = prototile_second->contour_sample[1];
		int t = 0;
		for (int i = 0; i < psize; i++)
		{
			int first = path[i].first;
			int sec = path[i].second;
			Point2f fin = 0.5 * contour1[first] + 0.5 * contour2[sec];;
			if (t<4 && first == mid_inter[t])
			{
				//cout << t << " " << first << "  " << i << endl;
				t++;
				mid_inter_new.push_back(i);
			}
			final_pettern.push_back(fin);
		}
		mid_inter_new = joint_relocate(final_pettern, mid_inter_new, 1);
		final_pettern = sampling(final_pettern, 2);

		if (mid_inter_new.size() != 4)
		{
			cout << "lack of tiling vertices!" << endl;
			return mid_inter_new;
		}

		int first = 0;
		int second = 0;

		vector<int> return_p;
		vector<vector<Point2f>> four_place;
		morphed_A = extract_contour(final_pettern, mid_inter_new, return_p, four_place, all_inner_conts[Tiling_index].type);   //

		int times = 0;
		while (self_intersect(morphed_A, first, second) && (times < 3))
		{
			++times;
			cout << "morphed_A self_intersect is repairing" << endl;
			contour_fine_tuning(morphed_A, first, second);
		}
		//cout << "times" << times << endl;
		four_place.swap(vector<vector<Point2f>>());
		mid_inter_new.swap(vector<int>());
		final_pettern = extract_contour(morphed_A, return_p, mid_inter_new, four_place, all_inner_conts[Tiling_index].type);
		cout << "lack of mid_inter_new: " << mid_inter_new.size() << endl;
		prototile_tem->loadPoints(final_pettern);
		return mid_inter_new;
	}

	/*void Tiling_opt::tiliing_generation(string imaname)
	{
		clock_t start, midtime, finish;
		start = clock();
		bool check_self_intersect = true;
		int num_c = 1;//选择(num_c+1)*100个点
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		cout << "p_p_index: " << p_p_index.size() << endl;
		vector<Point2f> contour_ = prototile_first->contour;
		vector<Point2f> cont_orig = prototile_first->contour_sample[1];
		vector<Point2f> cont_rota = prototile_first->contour_sample[4]; //旋转暂时用500个点
		int contsize = contour_.size();   //这里现在应该是200个
		Point2f cent_cont = center_p(contour_);
		//vector<Point2f> contour_ = prototile_first->contour;
		load_dataset();
		com_all_TARs(num_c);
		string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result300\\" + prototile_first->contourname;
		//string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result_tuning\\" + prototile_first->contourname;
		const char *na = rootname.c_str();
		if (_access(na, 0) != -1)

		{
			printf("The  file/dir had been Exisit ");
			return;
		}
		mkdir(na);
		int ppindex = p_p_index.size();
		int margin = prototile_first->contour.size() / 10;
		cout << "margin: " << margin << endl;
		int count = 0;
		int trans = 0;
		int rotas = 0;
		int flips = 0;
		vector<inPat> all_inner_conts;
		//vector<vector<Point2f>> inner_conts;
		//vector<vector<int>> all_situation_index;
		//vector<vector<int>> mid_interval_index;

		vector<vector<int>> all_result_index;
		for (int i = 0; i < ppindex; i++)    //计算旋转情况
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				 = find_rota_tilingV(cont_rota, mark_13);
			}
		}

		for (int i = 0; i < ppindex; i++)    //计算旋转情况
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				vector<Point2f> inner_contour;
				vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
				vector<int> mark_13;
				mark_13.push_back(p_p_index[i]);
				mark_13.push_back(p_p_index[j]);
				vector<vector<int>> all_result_index = find_rota_tilingV(cont_rota, mark_13);
				int allresultsize = all_result_index.size();
				for (int num = 0; num < allresultsize; num++)
				{
					Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
					if (!rotation_placement(all_result_index[num], contour_, inner_contour, mid_interval, drawing1))
					{
						++rotas;
						cout << ++count << " succeed" << endl;
						inPat one_situation;
						one_situation.in_contour = inner_contour;
						one_situation.in_interval = mid_interval;
						one_situation.type = 1;
						all_inner_conts.push_back(one_situation);
						inner_contour.swap(vector<Point2f>());
						mid_interval.swap(vector<int>());

						Point2f shift2 = Point2f(400, 400) - cent_cont;
						for (int jj = 0; jj < contsize; jj++)
						{
							circle(drawing1, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

							//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
						}
						for (int jj = 0; jj < ppindex; jj++)
						{
							circle(drawing1, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
						}
						for (int jj = 0; jj < 4; jj++)
						{
							circle(drawing1, contour_[all_result_index[num][jj]] + shift2, 8, Scalar(0, 255, 0), -1);
						}
						string filename = rootname + "\\" + int2string(count - 1) + "rota_PlacingResult.png";
						imwrite(filename, drawing1);
					}
				}
			}
		}

		//计算旋转和翻转的情况
		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				vector<Point2f> inner_contour;
				vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
				//vector<int> mark_13;
				//mark_13.push_back(p_p_index[i]);
				//mark_13.push_back(p_p_index[j]);
				//vector<vector<int>> all_result_index = find_rota_tilingV(cont_rota, mark_13);
				//int allresultsize = all_result_index.size();
				//for (int num = 0; num < allresultsize; num++)
				//{
				//	Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
				//	if (!rotation_placement(all_result_index[num], contour_, inner_contour, mid_interval, drawing1))
				//	{
				//		++rotas;
				//		cout << ++count << " succeed" << endl;
				//		inPat one_situation;
				//		one_situation.in_contour = inner_contour;
				//		one_situation.in_interval = mid_interval;
				//		one_situation.type = 1;
				//		all_inner_conts.push_back(one_situation);
				//		inner_contour.swap(vector<Point2f>());
				//		mid_interval.swap(vector<int>());
				//		Point2f shift2 = Point2f(400, 400) - cent_cont;
				//		for (int jj = 0; jj < contsize; jj++)
				//		{
				//			circle(drawing1, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);
				//			//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
				//		}
				//		for (int jj = 0; jj < ppindex; jj++)
				//		{
				//			circle(drawing1, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
				//		}
				//		for (int jj = 0; jj < 4; jj++)
				//		{
				//			circle(drawing1, contour_[all_result_index[num][jj]] + shift2, 8, Scalar(0, 255, 0), -1);
				//		}
				//		string filename = rootname + "\\" + int2string(count - 1) + "rota_PlacingResult.png";
				//		imwrite(filename, drawing1);
				//	}
				//}
				if (abs(p_p_index[j] - p_p_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(p_p_index[m] - p_p_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						//cout << i<<" "<<j<<" "<<m<<" "<<n << endl;
						if (abs(p_p_index[n] - p_p_index[m]) < margin) continue;                                                     
						vector<Mat> all_mat;
						vector<string> all_png;
						vector<int> result_index;
						result_index.push_back(p_p_index[i]);                                                                       //总共讨论一下三种摆放规律，
						result_index.push_back(p_p_index[j]);                                                                       //1.translation：以下所有该类摆放都以1-3,2-4为轴摆放 
						result_index.push_back(p_p_index[m]);                                                                       //2.rotation：以下所有该类摆放都按照1-2-3-4的顺序摆放
						result_index.push_back(p_p_index[n]);                                                                       //3:flipping:一下所有该类摆放都依次以1-3,2-4为轴旋转摆放			
						Mat drawing2 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						Mat drawing3 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						Mat drawing4 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));

						//其中flip又分为两种，分别沿1-3和2-4
						if (!translation_placement(result_index, contour_, inner_contour, mid_interval, drawing2))
						{
							cout << ++count << " succeed" << endl;
							++trans;
							all_mat.push_back(drawing2);
							all_png.push_back("trans");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 0;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
						}
						if (!flipping_placement(result_index, contour_, inner_contour, mid_interval, drawing3, 0))
						{
							++flips;
							cout << ++count << " succeed" << endl;
							all_mat.push_back(drawing3);
							all_png.push_back("flip(1-3)");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 2;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
						}
						if (!flipping_placement(result_index, contour_, inner_contour, mid_interval, drawing4, 1))
						{
							++flips;
							cout << ++count << " succeed" << endl;
							all_mat.push_back(drawing4);
							all_png.push_back("flip(2-4)");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 3;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
						}
						for (int num = 0; num < all_mat.size(); num++)
						{
							//show the marked points
							Point2f shift2 = Point2f(400, 400) - center_p(contour_);
							for (int jj = 0; jj < contsize; jj++)
							{
								circle(all_mat[num], contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

								//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
							}
							for (int jj = 0; jj < ppindex; jj++)
							{
								circle(all_mat[num], contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
							}
							for (int jj = 0; jj < 4; jj++)
							{
								circle(all_mat[num], contour_[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
							}

							string filename = rootname + "\\" + int2string(count - all_mat.size() + num) + all_png[num] + "PlacingResult.png";
							imwrite(filename, all_mat[num]);
						}						
					}
				}
				
			}
		}
		cout << "succeed count: " << count <<" trans: "<<trans<<" rotat: "<<rotas<<" flips: "<<flips<< endl;
		//midtime = clock();
		//cout << "Time consumption: " << (double)(midtime - start) / CLOCKS_PER_SEC << " s " << endl;
		if (count == 0)
		{
			cout << "no right placement" << endl;
			return;
		}

		for (int i = 0; i < all_inner_conts.size(); i++) //inner_conts.size()
		{
			cout << "count: " << i << "/" << count - 1 << endl;
			vector<pair<int, bool>> all_total = compare_choose_TAR(all_inner_conts[i].in_contour);
			//cout << "candida_contours" << candida_contours.size()<< endl;
			//string conut_name = rootname + "\\placement " + int2string(i);
			vector<int> mid_inter_ = joint_relocate(all_inner_conts[i].in_contour, all_inner_conts[i].in_interval, num_c);
			
			//if (abs(mid_inter[1] - mid_inter[0]) < 2 || abs(mid_inter[2] - mid_inter[1]) < 2 || abs(mid_inter[3] - mid_inter[2]) < 2 || abs(mid_inter[0] - mid_inter[3]) < 2)
			//{
			//	cout << "two tiling vertices are too closed!" << endl;
			//	continue;
			//}
			prototile_mid->Pro_clear();
			prototile_mid->loadPoints(all_inner_conts[i].in_contour);
			vector<Point2f> contour_inner = prototile_mid->contour_sample[1];
			double sc_inner = 0;
			vector<vector<double>> inner_tar = prototile_mid->compute_TAR(contour_inner, sc_inner);

			for (int j = 0; j <5; j++) //candida_contours.size()6                        
			{
				//将所有的结果保存下来
				vector<int> mid_inter = mid_inter_;
				prototile_second->Pro_clear();
				prototile_second->loadPoints(contour_dataset[all_total[j].first]);
				vector<Point2f> contour_cand;// = prototile_second->contour_sample[1];
				vector<vector<double>> cand_tar;
				if (all_total[j].second)
				{
					//cout << "it is flip" << endl;
					cand_tar = all_con_tars_flip[all_total[j].first];
					contour_cand = prototile_second->contour_sample_flip[1];
				}
				else
				{
					//cout << "it is not flip" << endl;
					cand_tar = all_con_tars[all_total[j].first];
					contour_cand = prototile_second->contour_sample[1];
				}
				vector<pair<int, int>> path;
				int shift = 0;
				int width = 6;
				double re = tar_mismatch(inner_tar, cand_tar, path, shift, width);
				vector<Point2f> mor_result = morphing_tar(contour_inner, contour_cand, mid_inter, path, shift);


				if (mid_inter.size() != 4)
				{
					cout << "lack of tiling vertices!" << endl;
					continue;
				}
				int first = 0;
				int second = 0;
				//if (self_intersect(mor_result, first, second)) continue;
				Mat drawing_pro = Mat(800, 2400, CV_8UC3, Scalar(255, 255, 255));
				draw_poly(drawing_pro, contour_inner, Point2f(400, 400));
				draw_poly(drawing_pro, contour_cand, Point2f(1200, 400));
				draw_poly(drawing_pro, mor_result, Point2f(2000, 400));

				Mat drawing_mid = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
				Mat drawing_mA = Mat(1600, 2400, CV_8UC3, Scalar(255, 255, 255));
				draw_allplane(drawing_mid, mor_result, mid_inter, 0.4, all_inner_conts[i].type);

				vector<int> return_p;
				vector<vector<Point2f>> four_place;
				vector<Point2f> morphed_A = extract_contour(mor_result, mid_inter, return_p, four_place, all_inner_conts[i].type);   //
				//if (self_intersect(morphed_A, first, second))
				//{
				//	cout << "morphed_A self intersected before morphing" << endl;
				//	continue;
				//}
				int times = 0;
				if (check_self_intersect)
				{
					while (self_intersect(morphed_A, first, second))
					{
						++times;
						cout << "morphed_A self_intersect is repairing" << endl;
						contour_fine_tuning(morphed_A, first, second);
					}
					if (times > 3) 
					{
						cout << "Too many tuning!" << endl;
						continue;
					}
					//cout << "times" << times << endl;
					four_place.swap(vector<vector<Point2f>>());
					mor_result = extract_contour(morphed_A, return_p, mid_inter, four_place, all_inner_conts[i].type);
					if (self_intersect(mor_result, first, second))
					{
						cout << "mor_result self intersected" << endl;
						continue;
					}
				}
				//if (self_intersect(morphed_A, first, second)) continue;
				double score_fir_r = evalua_deformation(morphed_A, cont_orig);
				double score_sec_r = evalua_deformation(mor_result, contour_cand);
				double re_score = sqrt((score_fir_r*score_fir_r + score_sec_r*score_sec_r)/2);
				cout << "score_fir_r: " << score_fir_r << "score_sec_r: " << score_sec_r << endl;

				draw_allplane(drawing_mA, morphed_A, return_p, 0.4, all_inner_conts[i].type);
				string filename = rootname + "\\";
				string file2 = filename + int2string(i) + "_" + int2string(j) + "result_" + double2string(re_score) + ".png";
				filename = filename + int2string(i) + "_Candidate_" + int2string(j) + ".png";
				imwrite(filename, drawing_pro);
				imwrite(file2, drawing_mA);
				imshow("tiling_result_of_mid: ", drawing_mid);
			}
		}
		finish = clock();
		std::cout << "Time consumption of " << imaname<<" : " << (double)(finish - start) / CLOCKS_PER_SEC << " s " << endl;
		
	}

	*/
	jointPat Tiling_opt::simulation_tar(string imaname, int inner_one, int cand_one)
	{
		bool check_self_intersect = true;
		int num_c = 1;//选择(num_c+1)*100个点
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		vector<Point2f> contour_ = prototile_first->contour;
		vector<Point2f> con_ori = prototile_first->contour_sample[1];
		int contsize = contour_.size();
		Point2f cent_cont = center_p(contour_);
		//vector<Point2f> contour_ = prototile_first->contour;
		load_dataset();
		com_all_TARs(num_c);
		string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result\\simula\\" + prototile_first->contourname;
		const char *na = rootname.c_str();
		mkdir(na);
		int ppindex = p_p_index.size();
		int margin = prototile_first->contour.size() / 10;
		cout << "margin: " << margin << endl;
		int count = 0;
		int trans = 0;
		int rotas = 0;
		int flips = 0;
		vector<inPat> all_inner_conts;

		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				vector<Point2f> inner_contour;
				vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
				vector<int> mark_13;
				mark_13.push_back(p_p_index[i]);
				mark_13.push_back(p_p_index[j]);
				vector<vector<int>> all_result_index = find_rota_tilingV(contour_, mark_13);
				int allresultsize = all_result_index.size();
				for (int num = 0; num < allresultsize; num++)
				{
					Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
					if (!rotation_placement(all_result_index[num], contour_, inner_contour, mid_interval, drawing1))
					{
						++rotas;
						cout << ++count << " succeed" << endl;
						inPat one_situation;
						one_situation.in_contour = inner_contour;
						one_situation.in_interval = mid_interval;
						one_situation.type = 1;
						all_inner_conts.push_back(one_situation);
						inner_contour.swap(vector<Point2f>());
						mid_interval.swap(vector<int>());
						if ((count - inner_one) == 1)
						{
							Point2f shift2 = Point2f(400, 400) - cent_cont;
							for (int jj = 0; jj < contsize; jj++)
							{
								circle(drawing1, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

								//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
							}
							for (int jj = 0; jj < ppindex; jj++)
							{
								circle(drawing1, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
							}
							for (int jj = 0; jj < 4; jj++)
							{
								circle(drawing1, contour_[all_result_index[num][jj]] + shift2, 8, Scalar(0, 255, 0), -1);
							}
							string filename = rootname + "\\" + int2string(count - 1) + "rota_PlacingResult_sim.png";
							imwrite(filename, drawing1);
							cout << "Point set: " << endl;
							for (int i = 0; i < 4; i++)
								cout << all_result_index[num][i] << "  ";
							cout << endl;
							break;
						}		
					}
				}
				if (abs(p_p_index[j] - p_p_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(p_p_index[m] - p_p_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						if (abs(p_p_index[n] - p_p_index[m]) < margin) continue;
						vector<Mat> all_mat;
						vector<string> all_png;
						vector<int> result_index;
						result_index.push_back(p_p_index[i]);                                                                       //总共讨论一下三种摆放规律，
						result_index.push_back(p_p_index[j]);                                                                       //1.translation：以下所有该类摆放都以1-3,2-4为轴摆放 
						result_index.push_back(p_p_index[m]);                                                                       //2.rotation：以下所有该类摆放都按照1-2-3-4的顺序摆放
						result_index.push_back(p_p_index[n]);                                                                       //3:flipping:一下所有该类摆放都依次以1-3,2-4为轴旋转摆放			
						Mat drawing2 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						Mat drawing3 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						Mat drawing4 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));

						//其中flip又分为两种，分别沿1-3和2-4
						if (!translation_placement(result_index, contour_, inner_contour, mid_interval, drawing2))
						{
							cout << ++count << " succeed" << endl;
							++trans;
							all_mat.push_back(drawing2);
							all_png.push_back("trans");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 0;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
							if ((count - inner_one) == 1)
							{
								Point2f shift2 = Point2f(400, 400) - cent_cont;
								for (int jj = 0; jj < contsize; jj++)
								{
									circle(drawing2, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

									//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
								}
								for (int jj = 0; jj < ppindex; jj++)
								{
									circle(drawing2, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
								}
								for (int jj = 0; jj < 4; jj++)
								{
									circle(drawing2, contour_[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
								}
								string filename = rootname + "\\" + int2string(count - 1)  + "transPlacingResult_sim.png";
								imwrite(filename, drawing2);
								cout << "Point set: " << endl;
								cout << p_p_index[i] << " " << p_p_index[j] << " " << p_p_index[m] << " " << p_p_index[n] << endl;
								break;
							}
						}

						if (!flipping_placement(result_index, contour_, inner_contour, mid_interval, drawing3, 0))
						{
							++flips;
							//cout << p_p_index[i] << "  " << p_p_index[j] << "  " << p_p_index[m] << "  " << p_p_index[n] << endl;
							cout << ++count << " succeed" << endl;
							all_mat.push_back(drawing3);
							all_png.push_back("flip(1-3)");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 2;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
							if ((count - inner_one) == 1)
							{
								Point2f shift2 = Point2f(400, 400) - cent_cont;
								for (int jj = 0; jj < contsize; jj++)
								{
									circle(drawing3, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

									//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
								}
								for (int jj = 0; jj < ppindex; jj++)
								{
									circle(drawing3, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
								}
								for (int jj = 0; jj < 4; jj++)
								{
									circle(drawing3, contour_[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
								}
								string filename = rootname + "\\" + int2string(count - 1) + "flip(1-3)PlacingResult_sim.png";
								imwrite(filename, drawing3);
								cout << "Point set: " << endl;
								cout << p_p_index[i] << " " << p_p_index[j] << " " << p_p_index[m] << " " << p_p_index[n] << endl;
								break;
							}
						}
						if (!flipping_placement(result_index, contour_, inner_contour, mid_interval, drawing4, 1))
						{
							++flips;
							//cout << p_p_index[i] << "  " << p_p_index[j] << "  " << p_p_index[m] << "  " << p_p_index[n] << endl;
							cout << ++count << " succeed" << endl;
							all_mat.push_back(drawing4);
							all_png.push_back("flip(2-4)");
							inPat one_situation;
							one_situation.in_contour = inner_contour;
							one_situation.in_interval = mid_interval;
							one_situation.type = 3;
							all_inner_conts.push_back(one_situation);
							inner_contour.swap(vector<Point2f>());
							mid_interval.swap(vector<int>());
							if ((count - inner_one) == 1)
							{
								Point2f shift2 = Point2f(400, 400) - cent_cont;
								for (int jj = 0; jj < contsize; jj++)
								{
									circle(drawing4, contour_[jj] + shift2, 1, Scalar(0, 0, 0), -1);

									//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
								}
								for (int jj = 0; jj < ppindex; jj++)
								{
									circle(drawing4, contour_[p_p_index[jj]] + shift2, 4, Scalar(0, 0, 255), -1);
								}
								for (int jj = 0; jj < 4; jj++)
								{
									circle(drawing4, contour_[result_index[jj]] + shift2, 8, Scalar(0, 255, 0), -1);
								}
								string filename = rootname + "\\" + int2string(count - 1) + "flip(2-4)PlacingResult_sim.png";
								imwrite(filename, drawing4);
								cout << "Point set: " << endl;
								cout << p_p_index[i] << " " << p_p_index[j] << " " << p_p_index[m] << " " << p_p_index[n] << endl;
								break;
							}
						}
					}
					if ((count - inner_one) == 1) break;
				}
				if ((count - inner_one) == 1) break;
			}
			if ((count - inner_one) == 1) break;
		}
		cout << "succeed count: " << count << endl;
		if (count == 0)
		{
			cout << "no right placement" << endl;
			exit(0);
		}
		
		string filepathname = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test" + int2string(inner_one)+".txt";
		vector<Point> con;
		for (int i = 0; i < all_inner_conts[inner_one].in_contour.size(); i++)
		{
			con.push_back((Point)all_inner_conts[inner_one].in_contour[i]);
		}
		fileout(filepathname, con);
		
		//Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(ttt, inner_conts[inner_one], Point2f(400, 400));
		//imwrite("D:\\ttt.png", ttt);
		cout << all_inner_conts.size() << "  :  " << inner_one << endl;
		/*vector<Point2f> inner_c = all_inner_conts[inner_one].in_contour;
		cout << "inner_c: " << inner_c.size() << endl;
		Mat drawing_ = Mat(600, 600, CV_8UC3, Scalar(255, 255, 255));

		for (int i = 0; i < inner_c.size(); i++)
		{
			circle(drawing_, inner_c[i], 2, Scalar(0, 255, 0), -1);
		}
		imwrite("D:\\asa.png", drawing_);
		imshow("asa", drawing_);
		return vector<Point2f>();*/
		vector<pair<int, bool>> all_total = compare_choose_TAR(all_inner_conts[inner_one].in_contour);
		for (int t = 0; t < 5; t++)
		{
			cout << all_total[t].first << endl;
		}
		//cout << "all_tital: " << all_total.size() << endl;
		cout << "all_inner_conts[inner_one].in_interval[i]: " << endl;
		for (int i = 0; i < 4; i++)
			cout << all_inner_conts[inner_one].in_interval[i] << "  ";
		cout << endl;
		vector<int> mid_inter = joint_relocate(all_inner_conts[inner_one].in_contour, all_inner_conts[inner_one].in_interval, num_c);
		cout << "mid_inter??: " << mid_inter.size() << endl;
		for (int i = 0; i < 4; i++)
			cout << mid_inter[i] << "  ";
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(all_inner_conts[inner_one].in_contour);
		vector<Point2f> contour_inner = prototile_mid->contour_sample[1];
		double sc_inner = 0;
		vector<vector<double>> inner_tar = prototile_mid->compute_TAR(contour_inner, sc_inner);

		prototile_second->Pro_clear();
		prototile_second->loadPoints(contour_dataset[all_total[cand_one].first]);
		vector<Point2f> contour_cand;// = prototile_second->contour_sample[1];
		vector<vector<double>> cand_tar;
		if (all_total[cand_one].second)
		{
			//cout << "it is flip" << endl;
			cand_tar = all_con_tars_flip[all_total[cand_one].first];
			contour_cand = prototile_second->contour_sample_flip[1];
		}
		else
		{
			//cout << "it is ont flip" << endl;
			cand_tar = all_con_tars[all_total[cand_one].first];
			contour_cand = prototile_second->contour_sample[1];
		}
		vector<pair<int, int>> path;
		int shift = 0;
		int width = 6;
		double re = tar_mismatch(inner_tar, cand_tar, path, shift, width);
		vector<Point2f> mor_result = morphing_tar(contour_inner, contour_cand, mid_inter, path, shift);
		if (mid_inter.size() != 4)
		{
			cout << "lack of tiling vertices!" << endl;
			exit(0);
		}
		int first = 0;
		int second = 0;
		/*while (self_intersect(mor_result, first, second))
		{
			cout << "mor_result self_intersect1" << endl;
			contour_fine_tuning(mor_result, mid_inter, first, second);
		}*/
		Point2f cente = center_p(mor_result);
		vector<int> return_p;
		vector<vector<Point2f>> four_;
		vector<Point2f> morphed_A = extract_contour(mor_result, mid_inter, return_p, four_, all_inner_conts[inner_one].type);
		
		int times = 0;
		if (check_self_intersect)
		{
			while (self_intersect(morphed_A, first, second))
			{
				++times;
				cout << " morphed_A.size: " << morphed_A.size() << endl;
				contour_fine_tuning(morphed_A, first, second);
			}
			if (times > 3)
			{
				cout << "Too many tuning!" << endl;
				exit(0);
			}
			cout << "times" << times<<endl;
			four_.swap(vector<vector<Point2f>>());
			mor_result = extract_contour(morphed_A, return_p, mid_inter, four_, all_inner_conts[inner_one].type);

			//if (self_intersect(mor_result, first, second)) return vector<Point2f>();
		}
		
		double score_fir_r = evalua_deformation(morphed_A, con_ori);
		double score_sec_r = evalua_deformation(mor_result, contour_cand);
		cout << "score_fir_r: " << score_fir_r << "score_sec_r: " << score_sec_r << endl;


		string filep = "D:\\VisualStudioProjects\\DihedralTessellation\\simulation\\";
		filep = filep + imaname + "_morA_" + int2string(inner_one) + ".txt";
		vector<Point> con1;
		for (int i = 0; i < morphed_A.size(); i++)
		{
			con1.push_back((Point)morphed_A[i]);
		}
		fileout(filep, con1);

		//将所有的结果保存下来
		Mat drawing_pro = Mat(800, 2400, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_mid = Mat(1200, 2400, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_mA = Mat(3000, 3000, CV_8UC3, Scalar(255, 255, 255));
		
		draw_poly(drawing_pro, contour_inner, Point2f(400, 400));
		draw_poly(drawing_pro, contour_cand, Point2f(1200, 400));
		draw_poly(drawing_pro, mor_result, Point2f(2000, 400));
		////show the result		
		//draw_allplane(drawing_mid, mor_result, mid_inter, 0.4, all_inner_conts[inner_one].type);
		//draw_result(drawing_mid, mor_result, mid_inter, 0.4, all_inner_conts[inner_one].type);
		//draw_two(drawing_mid, morphed_A, return_p, mor_result, mid_inter, 0.4, all_inner_conts[inner_one].type);
		Point2f shif1 = Point2f(400, 600) - center_p(mor_result);
		Point2f shif2 = Point2f(1200, 600) - center_p(morphed_A);
		Point2f shif3 = Point2f(2000, 600) - center_p(con_ori);
		draw_poly(drawing_mid, mor_result, Point2f(400, 600), 9);
		draw_poly(drawing_mid, morphed_A, Point2f(1200, 600), 10);
		draw_poly(drawing_mid, con_ori, Point2f(2000, 600), 13);
		for (int i = 0; i < mor_result.size(); i++)
		{
			//circle(drawing_mid, mor_result[i] + shif1, 3, Scalar(0, 0, 230), -1);
			MyLine(drawing_mid, mor_result[i] + shif1, mor_result[(i + 1) % mor_result.size()] + shif1, "black");
		}
		for (int i = 0; i < morphed_A.size(); i++)
		{
			MyLine(drawing_mid, morphed_A[i] + shif2, morphed_A[(i + 1) % morphed_A.size()] + shif2, "black");
			//circle(drawing_mid, morphed_A[i] + shif2, 3, Scalar(250, 100, 100), -1);
		}
		for (int i = 0; i < con_ori.size(); i++)
		{
			MyLine(drawing_mid, con_ori[i] + shif3, con_ori[(i + 1) % con_ori.size()] + shif3, "black");
			//circle(drawing_mid, con_ori[i] + shif3, 3, Scalar(230, 0, 0), -1);
		}
		
		imwrite("D:\\1.png", drawing_mid);
		//将该proto1以及相邻四个proto2展示出来
		
		Mat drawing_ex = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_ttt = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
		for (int t = 0; t < four_[0].size(); t++)
		{
			circle(drawing_ttt, four_[0][t], 2, Scalar(0, 0, 255), -1);
		}
		for (int t = 0; t < four_[1].size(); t++)
		{
			circle(drawing_ttt, four_[1][t], 2, Scalar(128, 0, 128), -1);
		}
		for (int t = 0; t < four_[2].size(); t++)
		{
			circle(drawing_ttt, four_[2][t], 2, Scalar(0, 255, 0), -1);
		}
		for (int t = 0; t < four_[3].size(); t++)
		{
			circle(drawing_ttt, four_[3][t], 2, Scalar(255, 128, 0), -1);
		}
		for (int t = 0; t < return_p.size(); t++)
		{
			circle(drawing_ttt, morphed_A[return_p[t]], 12, Scalar(0, 0, 0), -1);
		}
		imshow("all", drawing_ttt);
		cout << "zahuishia!" << endl;
		/*Point2f sss = Point2f(400, 400) - center_p(morphed_A);
		for (int t = 0; t < morphed_A.size(); t++)
		{
			circle(drawing_ex, morphed_A[t] + sss, 2, Scalar(0, 0, 255), -1);
		}
		draw_poly(drawing_ex, morphed_A,Point2f(1200,400));
		imshow("morphed_A: ", drawing_ex);*/
		

		//if (self_intersect(morphed_A, first, second)) cout << "self_intersect" << endl;
		
		draw_allplane(drawing_mA, morphed_A, return_p, 0.4, all_inner_conts[inner_one].type);
		//draw_result(drawing_mA, morphed_A, return_p, 0.4, all_inner_conts[inner_one].type);
		imwrite("D:\\model.png", drawing_mA);

		//Mat drawing_re = Mat(3600, 3600, CV_8UC3, Scalar(180, 180, 180));
		//draw_two(drawing_re, morphed_A, return_p, mor_result, mid_inter, 0.5, all_inner_conts[inner_one].type);
		//imwrite("D:\\result.png", drawing_re);

		string filename = rootname + "\\Candidate_" + int2string(cand_one) + ".png";
		string file2 = rootname + "\\Cand_" + int2string(cand_one) + "tiling_result.png";

		imwrite(filename, drawing_pro);
		imwrite(file2, drawing_mA);
		imshow("result_mid_show: ", drawing_mid);
		jointPat four_pattern;
		four_pattern.four_contour = four_;
		four_pattern.interval = return_p;
		four_pattern.type = all_inner_conts[inner_one].type;
		return four_pattern;
		//return contour_inner;
	}

	
	
	bool Tiling_opt::translation_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname)
	{
		int csize = contour_s.size();
		Point2f line1 = contour_s[results[2]] - contour_s[results[0]];
		Point2f line2 = contour_s[results[3]] - contour_s[results[1]];
		vector<Point2f> shifting;
		shifting.push_back(line1);
		shifting.push_back(line2);
		shifting.push_back(line1+line2);
		// 提取平移围成的轮廓
		vector<vector<Point2f>> four_place;
		vector<Point2f> one_loca;		
		four_place.push_back(contour_s);
		for (int i = 0; i < 3; i++)
		{
			Point2f move = shifting[i];
			for (int j = 0; j < csize; j++)
			{
				one_loca.push_back(contour_s[j] + move);
			}
			four_place.push_back(one_loca);
			one_loca.swap(vector<Point2f>());
		}	

		vector<Point2f> a[2];
		vector<Point2f> b[2];
		for (int i = 0; i < 4; i++)
		{                                                                                                     //translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
			if (i == 0)                                                                                       
			{																								  //              ①     ③            共有①②③④四个图案，
				a[0].push_back(four_place[0][results[2]]);                                                    //               /\    /\			   检测①中的顶点3与②中顶点1，	
				a[1].push_back(four_place[0][results[3]]);                                                    //              /1 \  /1 \           以及①中顶点4与③中顶点2的角度之和
				b[0].push_back(four_place[1][results[0]]);                                                    //             /    \/    \         
				b[1].push_back(four_place[2][results[1]]);                                                    //            2\   4/\2   /4
			}                                                                                                 //              \  /  \  /
			else                                                                                              //              3\/    \/3
			{                                                                                                 //               /\    /\     
				a[0].push_back(four_place[0][(results[2] - i + csize) % csize]);                              //              /1 \  /1 \ 
				a[1].push_back(four_place[0][(results[3] - i + csize) % csize]);                              //             /    \/    \          为保证在一定范围内可靠，
				a[0].push_back(four_place[0][(results[2] + i) % csize]);                                      //            2\   4/\2  4/          比较时计算∠Pt-iPtPt+i的角度，
				a[1].push_back(four_place[0][(results[3] + i) % csize]);                                      //              \  /  \  /           i从1取到3
				b[0].push_back(four_place[1][(results[0] - i + csize) % csize]);                              //              3\/    \/3
				b[1].push_back(four_place[2][(results[1] - i + csize) % csize]);                              //               ②    ④
				b[0].push_back(four_place[1][(results[0] + i) % csize]);                                      
				b[1].push_back(four_place[2][(results[1] + i) % csize]);
			}                                                                                                  
		}                                                                                         

		if (vertex_angle(a[0], b[0]) || vertex_angle(a[1], b[1])) return true;  //到这里用时0.009s
		if (coll_detec_bbx(four_place[0], four_place[1], 10) || coll_detec_bbx(four_place[0], four_place[2], 10) || coll_detec_bbx(four_place[0], four_place[3], 0) || coll_detec_bbx(four_place[1], four_place[2], 0))
			return true; 
		return_B = extract_contour(contour_s, results, return_p, four_place,0);
		//visual presentation
		Point2f shift1 = Point2f(1200, 400) - (0.4*center_p(contour_s) + 0.2*line1 + 0.2*line2);
		for (int i = 0; i < 4; i++)
		{
			vector<Point2f> one_;
			for (int j = 0; j < four_place[i].size(); j++)
			{
				Point2f p = four_place[i][j] * 0.4 + shift1;
				one_.push_back(p);
			}
			draw_poly(countname, one_, center_p(one_));
		}
		
		return false;
	}

	bool Tiling_opt::rotation_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname)
	{
		
		Line_Seg line1(contour_s[results[2]], contour_s[results[0]]);
		Line_Seg line2(contour_s[results[3]], contour_s[results[1]]);
		//cout << results[0] << " " << results[1] << " " << results[2] << " " << results[3] << endl;
		//cout << line1.start << " " << line1.end << endl << line2.start << " " << line2.end << endl;
		Point2f cross_p;
		if (line_intersection(line1, line2, cross_p) != 1) return true;   //rotation 里
		if (abs(length_two_point2f(contour_s[results[2]], cross_p) - length_two_point2f(contour_s[results[0]], cross_p)) > 5) return true;
		if (abs(length_two_point2f(contour_s[results[3]], cross_p) - length_two_point2f(contour_s[results[1]], cross_p)) > 5) return true;
		int csize = contour_s.size();
		// 提取旋转围成的轮廓
		vector<vector<Point2f>> four_place;
		vector<Point2f> one_loca = contour_s;
		four_place.push_back(one_loca);
		
		for (int i = 0; i < 3; i++)
		{
			Point2f rota_cent = one_loca[results[i]];
			Mat rot_mat = getRotationMatrix2D(rota_cent, 180, 1);
			transform(one_loca, one_loca, rot_mat);
			four_place.push_back(one_loca);                                                                       //rotation的提取规律为1的1-4，4的4-3,3的3-2,2的2-1
		}
		vector<Point2f> a[4];                                                                                     //               ②    ③            共有①②③④四个图案，
		for (int i = 0; i < 4; i++)                                                                               //               /\    /\			   对于旋转的情况来说，接触点的角为同一个角，
		{                                                                                                         //              /3 \  /1 \           因此只需检测①和②中的顶点3，以此类推
			if (i == 0)                                                                                           //             /    \/    \ 
			{                                                                                                     //            4\   2/\2   /4
				for (int j = 0; j < 4; j++)                                                                       //              \  /  \  /
				{																								  //              1\/    \/3
					a[j].push_back(four_place[0][results[j]]);                                                    //               /\    /\						                                                                      
				}	                                                                                              //              /1 \  /3 \  
		}		                                                                                                  //             /    \/    \          为保证在一定范围内可靠，
			else                                                                                                  //            2\   4/\4  2/          比较时计算∠Pt-iPtPt+i的角度，
			{	                                                                                                  //              \  /  \  /           i从1取到3
			    for (int j = 0; j < 4; j++)                                                                       //              3\/    \/1          
			    {                                                                                                 //               ①    ④
					a[j].push_back(four_place[0][(results[j] - i + csize) % csize]);
					a[j].push_back(four_place[0][(results[j] + i) % csize]);                                                                                            
				}	                              
			}		                              
		}

		if (vertex_angle(a[0], a[0]) || vertex_angle(a[1], a[1]) || vertex_angle(a[2], a[2]) || vertex_angle(a[3], a[3])) return true; //到这里用时0.009s

		for (int i = 0; i < 4; i++)
		{
			for (int j = i + 1; j < 4; j++)
			{
				if (coll_detec_bbx(four_place[i], four_place[j], 10))
					return true;
			}
		}
		return_B = extract_contour(contour_s, results, return_p, four_place,1);
		//visual presentation
		vector<Point2f> all_cent;
		Point2f allcenter = Point2f(0, 0);
		for (int i = 0; i < 4; i++)
		{
			Point2f t = center_p(four_place[i]);
			all_cent.push_back(t);
			allcenter += t;
		}
		Point2f shift = Point2f(1200, 400) - 0.1*allcenter;
		for (int i = 0; i < 4; i++)
		{
			vector<Point2f> one_;
			for (int j = 0; j < four_place[i].size(); j++)
			{
				Point2f p = four_place[i][j] * 0.4 + shift;
				one_.push_back(p);
			}
			draw_poly(countname, one_, all_cent[i] * 0.4 + shift);
		}
		return false;
	}
	
	bool Tiling_opt::flipping_placement(vector<int> results, vector<Point2f> &contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname, int type)
	{
		//0:以1-3对应的中线为轴翻转 1:以2-4为轴 
		//不一定需要以1-3为轴，只要1,3两点到对称轴的距离相同
		//因此第一步是根据1,3两点确定翻转轴
		int csize = contour_s.size();
		Point2f cent = center_p(contour_s);
		Point2f line1 = contour_s[results[2]] - contour_s[results[0]];
		Point2f line2 = contour_s[results[3]] - contour_s[results[1]];	
		//if (abs(cos_two_vector(line1, line2))>0.06) return true;
		//cout << line1 << endl << line2 << endl << cos_two_vector(line1, line2) << endl;
		vector<vector<Point2f>> four_place;
		four_place.push_back(contour_s);
		vector<Point2f> one_loca;
		vector<Point2f> one_loca_;
		vector<Point2f> one_loca_f;
		// 提取翻转围成的轮廓
		if (type == 0)
		{
			//实际的对称轴应该是过1，3中点且与2，4垂直的直线
			Point2f line1_cent = 0.5 * contour_s[results[2]] + 0.5 * contour_s[results[0]];
			Point2f s_axis(line2.y, -line2.x);
			Line_Seg symmetry_axis(line1_cent + 5 * s_axis, line1_cent - 5 * s_axis);
			//计算翻转后的轮廓，因为直接任意轴对称翻转比较耗时，因此先以x轴为对称轴翻转，然后旋转相应角度
			double rota_angle = acos(cos_two_vector(Point2f(0, 1), s_axis)) / PI * 180;
			if (sin_two_vector(Point2f(0, 1), s_axis) > 0) rota_angle = -rota_angle;
			//cout << rota_angle << endl;
			Point2f cross;
			if (line_intersection(Line_Seg(symmetry_axis), Line_Seg(cent - 5 * line2, cent + 5 * line2), cross) != 1) return true;  //没有交点说明出错了
			Point2f cent_shift = 2 * (cross - cent);
			one_loca = flip_only_coord(contour_s);
			Mat rot_mat = getRotationMatrix2D(cent, 2 * rota_angle, 1);
			transform(one_loca, one_loca, rot_mat);
			//将以上所得翻转的轮廓沿相对应的轴平移，计算该位置
			double cos_ax_1_3 = cos_two_vector(s_axis, line1);
			Point2f line3 = s_axis;
			if (cos_ax_1_3 < 0) line3 = -s_axis;
			//这里line3代表对称轴上与line1成锐角夹角方向的向量，长度为line1在对称轴上的投影
			line3 = line3*(abs(line1.x*line3.x + line1.y*line3.y) / (line3.x*line3.x + line3.y*line3.y)); 
			Point2f shift = cent_shift + line3;
			for (int i = 0; i < csize; i++)
			{
				one_loca[i] += shift;
				one_loca_.push_back(contour_s[i] + line2);
				one_loca_f.push_back(one_loca[i] + line2);
			}
			four_place.push_back(one_loca);
			four_place.push_back(one_loca_);
			four_place.push_back(one_loca_f);
		}
		else if (type == 1)
		{
			//实际的对称轴应该是过2，4中点且与1，3垂直的直线
			Point2f line2_cent = 0.5 * contour_s[results[3]] + 0.5 * contour_s[results[1]];
			Point2f s_axis(line1.y, -line1.x);
			Line_Seg symmetry_axis(line2_cent + 5 * s_axis, line2_cent - 5 * s_axis);

			double rota_angle = acos(cos_two_vector(Point2f(0, 1), s_axis)) / PI * 180;
			if (sin_two_vector(Point2f(0, 1), s_axis) > 0) rota_angle = -rota_angle;
			//cout << rota_angle << endl;
			Point2f cross;
			if (line_intersection(Line_Seg(symmetry_axis), Line_Seg(cent - 5 * line1, cent + 5 * line1), cross) != 1) return true;  //没有交点说明出错了
			Point2f cent_shift = 2 * (cross - cent);
			one_loca = flip_only_coord(contour_s);
			Mat rot_mat = getRotationMatrix2D(cent, 2 * rota_angle, 1);
			transform(one_loca, one_loca, rot_mat);

			double cos_ax_2_4 = cos_two_vector(s_axis, line2);
			Point2f line3 = s_axis;
			if (cos_ax_2_4 < 0) line3 = -s_axis;
			line3 = line3*(abs(line2.x*line3.x + line2.y*line3.y) / (line3.x*line3.x + line3.y*line3.y)); //这里line3代表对称轴上与line2成锐角夹角方向的向量，长度为line2在对称轴上的投影
			Point2f shift = cent_shift + line3;
			for (int i = 0; i < csize; i++)
			{
				one_loca[i] += shift;
				one_loca_.push_back(contour_s[i] + line1);
				one_loca_f.push_back(one_loca[i] + line1);
			}
			four_place.push_back(one_loca_);
			four_place.push_back(one_loca);			
			four_place.push_back(one_loca_f);
		}
		//Mat draw1 = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
		////for (int i = 1; i < 3; i++)
		//for (int j = 0; j < four_place[0].size(); j++)
		//{
		//	circle(draw1, four_place[0][j], 2, Scalar(255, 0, 0), -1);
		//}
		//imshow("dadas",draw1);
		//cout << four_place.size() << endl;
		vector<Point2f> a[2];                                                                                 // flippling(1-3)的提取规律为1的4-3, 2的1-2, 4的4-1, 3的3-2
		vector<Point2f> b[2];                                                                                 // flippling(2-4)的提取规律为1的4-3, 2的1-4, 4的2-3, 3的1-2
		for (int i = 0; i < 4; i++)
		{
			if (i == 0)                                                                                       //     沿1-3轴翻转        沿2-4轴翻转
			{																								  //       ①    ③          ①    ③          共有①②③④四个图案，
				a[0].push_back(four_place[0][results[2]]);                                                    //       /\    /\	         /\    /\		   检测①中的顶点3与②中顶点1，	
				a[1].push_back(four_place[0][results[3]]);                                                    //      /1 \  /1 \        /1 \  /3 \         以及①中顶点4与③中顶点2的角度之和
				b[0].push_back(four_place[1][results[0]]);                                                    //     /    \/    \      /    \/    \         
				b[1].push_back(four_place[2][results[1]]);                                                    //    2\   4/\2   /4    2\   4/\2   /4
			}                                                                                                 //      \  /  \  /        \  /  \  /
			else                                                                                              //      3\/    \/3        3\/    \/1
			{                                                                                                 //       /\    /\          /\    /\     
				a[0].push_back(four_place[0][(results[2] - i + csize) % csize]);                              //      /1 \  /1 \        /1 \  /3 \ 
				a[1].push_back(four_place[0][(results[3] - i + csize) % csize]);                              //     /    \/    \      /    \/    \        为保证在一定范围内可靠，
				a[0].push_back(four_place[0][(results[2] + i) % csize]);                                      //    4\   2/\4  2/     2\   4/\2  4/          比较时计算∠Pt-iPtPt+i的角度，
				a[1].push_back(four_place[0][(results[3] + i) % csize]);                                      //      \  /  \  /        \  /  \  /           i从1取到3
				b[0].push_back(four_place[1][(results[0] - i + csize) % csize]);                              //      3\/   3\/         3\/   1\/
				b[1].push_back(four_place[2][(results[1] - i + csize) % csize]);                              //       ②    ④          ②    ④
				b[0].push_back(four_place[1][(results[0] + i) % csize]);
				b[1].push_back(four_place[2][(results[1] + i) % csize]);
			} 
		}
		if (vertex_angle(a[0], b[0]) || vertex_angle(a[1], b[1])) return true;  //到这里用时0.009s
		if (coll_detec_bbx(four_place[0], four_place[1], 10) || coll_detec_bbx(four_place[0], four_place[2], 10) || coll_detec_bbx(four_place[0], four_place[3], 0) || coll_detec_bbx(four_place[1], four_place[2], 0))
			return true;
		if (type == 0) return_B = extract_contour(contour_s, results, return_p, four_place, 2);
		if (type == 1) return_B = extract_contour(contour_s, results, return_p, four_place, 3);
		////visual presentation
		Point2f shift1 = Point2f(1200, 400) - (0.4*center_p(contour_s) + 0.2*line1 + 0.2*line2);
		for (int i = 0; i < 4; i++)
		{
			vector<Point2f> one_;
			for (int j = 0; j < four_place[i].size(); j++)
			{
				Point2f p = four_place[i][j] * 0.4 + shift1;
				one_.push_back(p);
			}
			draw_poly(countname, one_, center_p(one_));
		}
		return false;
	}

	vector<Point2f> Tiling_opt::extract_contour(vector<Point2f> contour_, vector<int> mark_p, vector<int> &midmark_p, vector<vector<Point2f>> &four_place, int type)
	{   
		//  0:translation  1:rotation  2:flipping(1-3)   3:flipping(2-4)
		int csize = contour_.size();
		Point2f cent = center_p(contour_);
		vector<Point2f> morphed_B;
		Point2f line1 = contour_[mark_p[2]] - contour_[mark_p[0]];
		Point2f line2 = contour_[mark_p[3]] - contour_[mark_p[1]];
		if (type == 0)   // 提取translation围成的轮廓
		{
			if (four_place.empty())
			{
				vector<Point2f> shifting;
				shifting.push_back(line1);
				shifting.push_back(line2);
				shifting.push_back(line1 + line2);
				vector<Point2f> one_loca;
				four_place.push_back(contour_);
				for (int i = 0; i < 3; i++)
				{
					Point2f move = shifting[i];
					for (int j = 0; j < csize; j++)
					{
						one_loca.push_back(contour_[j] + move);
					}
					four_place.push_back(one_loca);
					one_loca.swap(vector<Point2f>());
				}
			}
			//cout << four_place.size() << endl;
			/*Mat drawing_ttt = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
			for (int i = 0; i < 4; i++)
				for (int t = 0; t < four_place[i].size(); t++)
				{
					circle(drawing_ttt, four_place[i][t], 2, Scalar(0, 0, 255), -1);
				}*/
			int total_num = 0;
			midmark_p.push_back(0);
			for (int t = mark_p[3]; t > mark_p[2]; t--)                   //translation的提取规律为1的4-3,2的1-4,4的2-1,3的3-2
			{
				total_num++;
				morphed_B.push_back(four_place[0][t]);
				//circle(drawing_ttt, four_place[0][t], 2, Scalar(0, 255, 0), -1);
			}
			midmark_p.push_back(total_num);
			for (int t = mark_p[0] + csize; t > mark_p[3]; t--)
			{
				total_num++;
				morphed_B.push_back(four_place[1][t % csize]);
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
		}
		else if (type == 1)   // 提取rotation围成的轮廓
		{
			if (four_place.empty())
			{
				vector<Point2f> one_loca = contour_;
				four_place.push_back(one_loca);
				for (int i = 0; i < 3; i++)
				{
					Point2f rota_cent = one_loca[mark_p[i]];
					Mat rot_mat = getRotationMatrix2D(rota_cent, 180, 1);
					cv::transform(one_loca, one_loca, rot_mat);
					four_place.push_back(one_loca);
				}
			}
			//Mat drawing_ttt = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
			//for (int i = 0; i < 4; i++)
			//{
			//	//circle(drawing_ttt, four_place[0][mark_p[i]] + Point2f(300, 0), 15, Scalar(0, 124, 0), -1);
			//	for (int t = 0; t < four_place[i].size(); t++)
			//	{
			//		circle(drawing_ttt, four_place[i][t] + Point2f(300, 0), 2, Scalar(0, 0, 255), -1);
			//	}
			//}
			//circle(drawing_ttt, four_place[0][mark_p[0]] + Point2f(300, 0), 15, Scalar(255, 124, 0), -1);
			int total_num = 0;
			midmark_p.push_back(0);
			if (mark_p[0] < mark_p[3])
			{
				for (int t = mark_p[0] + csize; t > mark_p[3]; t--)                         //rotation的提取规律为1的1-4，4的4-3,3的3-2,2的2-1
				{
					total_num++;
					morphed_B.push_back(four_place[0][t % csize]);
					//circle(drawing_ttt, four_place[0][t % csize] + Point2f(300, 0), 2, Scalar(0, 255, 0), -1);
				}
				midmark_p.push_back(total_num);
				for (int t = mark_p[3]; t > mark_p[2]; t--)
				{
					total_num++;
					morphed_B.push_back(four_place[3][t]);
				}
				midmark_p.push_back(total_num);
			}
			else
			{
				for (int t = mark_p[0]; t > mark_p[3]; t--)                         //rotation的提取规律为1的1-4，4的4-3,3的3-2,2的2-1，1234顺序翻转
				{
					total_num++;
					morphed_B.push_back(four_place[0][t]);
				}
				midmark_p.push_back(total_num);
				for (int t = mark_p[3] + csize; t > mark_p[2]; t--)
				{
					total_num++;
					morphed_B.push_back(four_place[3][t % csize]);
				}
				midmark_p.push_back(total_num);
			}
			for (int t = mark_p[2]; t > mark_p[1]; t--)
			{
				total_num++;
				morphed_B.push_back(four_place[2][t]);
			}
			midmark_p.push_back(total_num);
			for (int t = mark_p[1]; t > mark_p[0]; t--)
			{
				morphed_B.push_back(four_place[1][t]);
			}	
			//imwrite("D:\\rotation.png ", drawing_ttt);
		}
		else if (type == 2)
		{
			if (four_place.empty())
			{
				four_place.push_back(contour_);
				vector<Point2f> one_loca;
				vector<Point2f> one_loca_;
				vector<Point2f> one_loca_f;
				Point2f line1_cent = 0.5 * contour_[mark_p[2]] + 0.5 * contour_[mark_p[0]];
				Point2f s_axis(line2.y, -line2.x);
				Line_Seg symmetry_axis(line1_cent + 5 * s_axis, line1_cent - 5 * s_axis);
				//计算翻转后的轮廓，因为直接任意轴对称翻转比较耗时，因此先以x轴为对称轴翻转，然后旋转相应角度
				double rota_angle = acos(cos_two_vector(Point2f(0, 1), s_axis)) / PI * 180;
				if (sin_two_vector(Point2f(0, 1), s_axis) > 0) rota_angle = -rota_angle;
				//cout << rota_angle << endl;
				Point2f cross;
				if (line_intersection(Line_Seg(symmetry_axis), Line_Seg(cent - 5 * line2, cent + 5 * line2), cross) != 1) cout << "Error!" << endl;  //没有交点说明出错了
				Point2f cent_shift = 2 * (cross - cent);
				one_loca = flip_only_coord(contour_);
				Mat rot_mat = getRotationMatrix2D(cent, 2 * rota_angle, 1);
				transform(one_loca, one_loca, rot_mat);
				//将以上所得翻转的轮廓沿相对应的轴平移，计算该位置
				double cos_ax_1_3 = cos_two_vector(s_axis, line1);
				Point2f line3 = s_axis;
				if (cos_ax_1_3 < 0) line3 = -s_axis;
				line3 = line3*(abs(line1.x*line3.x + line1.y*line3.y) / (line3.x*line3.x + line3.y*line3.y));
				Point2f shift = cent_shift + line3;
				for (int i = 0; i < csize; i++)
				{
					one_loca[i] += shift;
					one_loca_.push_back(contour_[i] + line2);
					one_loca_f.push_back(one_loca[i] + line2);
				}
				four_place.push_back(one_loca);
				four_place.push_back(one_loca_);
				four_place.push_back(one_loca_f);
			}
			int total_num = 0;
			midmark_p.push_back(0);                                                        //flippling(1-3)的提取规律为1的4-3,2的1-2,4的4-1,3的3-2		
			for (int t = mark_p[3]; t > mark_p[2]; t--)
			{
				total_num++;
				morphed_B.push_back(four_place[0][t]);
			}
			midmark_p.push_back(total_num);
			for (int t = mark_p[0]; t < mark_p[1]; t++)
			{
				total_num++;
				morphed_B.push_back(four_place[1][t]);
			}
			midmark_p.push_back(total_num);
			for (int t = mark_p[3]; t <mark_p[0]+csize; t++)
			{
				total_num++;
				morphed_B.push_back(four_place[3][t%csize]);
			}
			midmark_p.push_back(total_num);		
			for (int t = mark_p[2]; t > mark_p[1]; t--)
			{
				morphed_B.push_back(four_place[2][t]);
			}
				
		}
		else if (type == 3)
		{
			if (four_place.empty())
			{
				four_place.push_back(contour_);
				vector<Point2f> one_loca;
				vector<Point2f> one_loca_;
				vector<Point2f> one_loca_f;
				Point2f line2_cent = 0.5 * contour_[mark_p[3]] + 0.5 * contour_[mark_p[1]];
				Point2f s_axis(line1.y, -line1.x);
				Line_Seg symmetry_axis(line2_cent + 5 * s_axis, line2_cent - 5 * s_axis);

				double rota_angle = acos(cos_two_vector(Point2f(0, 1), s_axis)) / PI * 180;
				if (sin_two_vector(Point2f(0, 1), s_axis) > 0) rota_angle = -rota_angle;
				//cout << rota_angle << endl;
				Point2f cross;
				if (line_intersection(Line_Seg(symmetry_axis), Line_Seg(cent - 5 * line1, cent + 5 * line1), cross) != 1) cout<<"Error!"<<endl;  //没有交点说明出错了
				Point2f cent_shift = 2 * (cross - cent);
				one_loca = flip_only_coord(contour_);
				Mat rot_mat = getRotationMatrix2D(cent, 2 * rota_angle, 1);
				transform(one_loca, one_loca, rot_mat);

				double cos_ax_2_4 = cos_two_vector(s_axis, line2);
				Point2f line3 = s_axis;
				if (cos_ax_2_4 < 0) line3 = -s_axis;
				line3 = line3*(abs(line2.x*line3.x + line2.y*line3.y) / (line3.x*line3.x + line3.y*line3.y));
				Point2f shift = cent_shift + line3;
				for (int i = 0; i < csize; i++)
				{
					one_loca[i] += shift;
					one_loca_.push_back(contour_[i] + line1);
					one_loca_f.push_back(one_loca[i] + line1);
				}
				four_place.push_back(one_loca_);
				four_place.push_back(one_loca);				
				four_place.push_back(one_loca_f);
			}
			int total_num = 0;
			midmark_p.push_back(0);                                                      //flippling(2-4)的提取规律为1的4-3,2的1-4,4的2-3,3的1-2
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
			for (int t = mark_p[1]; t < mark_p[2]; t++)
			{
				total_num++;
				morphed_B.push_back(four_place[3][t]);
			}
			midmark_p.push_back(total_num);
			for (int t = mark_p[0]; t < mark_p[1]; t++)
			{
				morphed_B.push_back(four_place[2][t]);
			}
		}
		/*Mat drawing_ttt = Mat(1600, 1600, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i < 4; i++)
		{
			draw_poly(drawing_ttt, four_place[i], center_p(four_place[i]));
		}
		for (int i = 0; i < morphed_B.size(); i++)
			MyLine(drawing_ttt, morphed_B[i], morphed_B[(i + 1) % morphed_B.size()], "green2");
		imshow("extracted contour",drawing_ttt);*/

		return morphed_B;
	}

	bool Tiling_opt::coll_detec_bbx(vector<Point2f> contour1, vector<Point2f> contour2,int threshold)
	{
		int csize = contour1.size();
		vector<Point2f> bbox1 = b_box(contour1);
		Point2f shift = Point2f(0, 0) - bbox1[1];
		for (int i = 0; i < csize; i++)
		{
			contour1[i] += shift;
			contour2[i] += shift;
		}
		Point2f cen1 = center_p(contour1);
		Point2f cen2 = center_p(contour2);
		bbox1 = b_box(contour1);
		vector<Point2f> bbox2 = b_box(contour2);
		/*Mat drawing3 = Mat(1200, 1200, CV_8UC1, Scalar(255));
		draw_poly(drawing3, contour1, cen1);		
		draw_poly(drawing3, contour2, cen2);
		for (int i = 0; i < bbox1.size(); i++)
		{
			MyLine(drawing3, bbox1[i], bbox1[(i + 1) % bbox1.size()], "");
			MyLine(drawing3, bbox2[i], bbox2[(i + 1) % bbox2.size()], "");
		}
		imshow("draw3", drawing3);*/
		double margin = max(length_two_point2f(bbox1[0], bbox1[2]), length_two_point2f(bbox2[0], bbox2[2]));
		//cout << "margin "<<margin<<endl;
		vector<Point2f> search_boundary;
		vector<int> num_in1;
		vector<int> num_in2;
		for (int i = 0; i < 4; i++)
		{
			if (bbox1[i].x <= bbox2[3].x && bbox1[i].x >= bbox2[1].x)
				if (bbox1[i].y <= bbox2[3].y && bbox1[i].y >= bbox2[1].y)
				{
					num_in1.push_back(i);
				}
		}
		for (int i = 0; i < 4; i++)
		{
			if (bbox2[i].x <= bbox1[3].x && bbox2[i].x >= bbox1[1].x)
				if (bbox2[i].y <= bbox1[3].y && bbox2[i].y >= bbox1[1].y)
				{
					num_in2.push_back(i);
				}
		}
		//cout << "num_in1.size() and num_in2.size()" << num_in1.size() << "   " << num_in2.size() << endl;
		if (num_in1.size() == 0 && num_in2.size() == 0)
		{
			if (length_two_point2f(cen1, cen2) < margin)
			{
				search_boundary.push_back(Point2f(0, 0));
				search_boundary.push_back(Point2f(599, 599));
			}
			else return false;
		}
		else if (num_in1.size() == 1 && num_in2.size() == 1)
		{
			double bbx_max_x = max(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x);
			double bbx_max_y = max(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y);
			double bbx_min_x = min(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x);
			double bbx_min_y = min(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y);
			if ((abs(bbx_max_x - bbx_min_x) < 1) || (abs(bbx_max_y - bbx_min_y) < 1)) return false;
			search_boundary.push_back(Point2f(bbx_min_x, bbx_min_y));
			search_boundary.push_back(Point2f(bbx_max_x, bbx_max_y));
		}
		else if (num_in1.size() == 2 && num_in2.size() == 2)
		{
			double bbx_max_x = max(max(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x), bbox1[num_in1[1]].x);
			double bbx_max_y = max(max(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y), bbox1[num_in1[1]].y);
			double bbx_min_x = min(min(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x), bbox1[num_in1[1]].x);
			double bbx_min_y = min(min(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y), bbox1[num_in1[1]].x);
			if ((abs(bbx_max_x - bbx_min_x) < 1) || (abs(bbx_max_y - bbx_min_y) < 1)) return false;
			search_boundary.push_back(Point2f(bbx_min_x, bbx_min_y));
			search_boundary.push_back(Point2f(bbx_max_x, bbx_max_y));
		}
		else if (num_in1.size() != num_in2.size())
		{
			search_boundary.push_back(Point2f(0, 0));
			search_boundary.push_back(Point2f(599, 599));
		}
		//将第一个图画到画布上，计算搜索区域的点

		//MyLine(drawing3, search_boundary[0], search_boundary[1], "black");

		Mat drawing1 = Mat(600, 600, CV_8UC1, Scalar(255));
		Mat drawing2 = Mat(600, 600, CV_8UC1, Scalar(255));
		draw_poly(drawing1, contour1, cen1);
		draw_poly(drawing2, contour2, cen2);

		//imshow("draw1111", drawing1);
		//imshow("draw2222", drawing2);

		//imwrite("D:\\tu.png", drawing3);
		//cout << "1: " << search_boundary[0].x << "  2: " << search_boundary[1].x << "  3:" << search_boundary[0].y << "   4:" << search_boundary[1].y << endl;
		if (search_boundary[0].x < 0 || search_boundary[1].x>600 || search_boundary[0].y < 0 || search_boundary[1].y>600) cout << "search_boundary beyond!" << endl;
		int count = 0;
		for (int i = search_boundary[0].y; i <= search_boundary[1].y + 0.5; i++)
			for (int j = search_boundary[0].x; j <= search_boundary[1].x + 0.5; j++)
			{
				if ((int)drawing1.at<uchar>(i, j) == 0 && (int)drawing2.at<uchar>(i, j) == 0)
				   count++;
			}
		//cout << count << endl;
		if (count > threshold) return true;
		else return false;

	}

	bool Tiling_opt::vertex_angle(vector<Point2f> angle1, vector<Point2f> angle2)
	{
		//  >360为true,<360为false
		int size = (angle1.size() - 1) / 2;
		for (int i = 0; i < size; i++)
		{
			double ang1 = 0;
			double ang2 = 0;
			double angle_cos = cos_two_vector((angle1[2 * i + 1] - angle1[0]), (angle1[2 * i + 2] - angle1[0]));
			double angle_sin = sin_two_vector((angle1[2 * i + 1] - angle1[0]), (angle1[2 * i + 2] - angle1[0]));
			if (angle_sin > 0) ang1 = acos(angle_cos) / PI * 180;
			else ang1 = 360 - acos(angle_cos) / PI * 180;

			angle_cos = cos_two_vector((angle2[2 * i + 1] - angle2[0]), (angle2[2 * i + 2] - angle2[0]));
			angle_sin = sin_two_vector((angle2[2 * i + 1] - angle2[0]), (angle2[2 * i + 2] - angle2[0]));
			if (angle_sin > 0) ang2 = acos(angle_cos) / PI * 180;
			else ang2 = 360 - acos(angle_cos) / PI * 180;
			if (ang1 + ang2 > 361)
			{
				//cout << ang1 <<"  "<<ang2<< endl;
				return true;
			}
		}
		
		return false;
	}

	vector<pair<int,bool>> Tiling_opt::compare_choose_TAR(vector<Point2f> inner_c)
	{
		int method = 1;
		int match_width = 4;
		vector<pair<int, bool>> all_total;

		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_c);
		vector<Point2f> contour_mid = prototile_mid->contour_sample[1];

		if (method == 1) all_total = quick_choose_TAR(contour_mid);
		else if (method == 2)
		{
			for (int i = 0; i < all_types; i++) 
				all_total.push_back(make_pair(i, true));
		}		
		vector<pair<int, bool>> all_final;
		vector<double> all_result;
		int total_num = all_total.size();
		//cout << "all total" << total_num << endl;
	
		double shape_com_mid;
		vector<vector<double>> tar_mid = prototile_mid->compute_TAR(contour_mid, shape_com_mid);
		//cout << "contour_mid: " << contour_mid.size() << "  tar_mid: " << tar_mid.size() << endl;
		for (int can_num = 0; can_num < total_num; can_num++)
		{
			int index = all_total[can_num].first;
			//cout << index << "  : ";
			vector<vector<double>> tar_sec = all_con_tars[index];
			vector<vector<double>> tar_sec_f = all_con_tars_flip[index];
			vector<pair<int, int>> path;
			int shift = 0;
			double re = tar_mismatch(tar_mid, tar_sec, path, shift, match_width);
			double re2 = tar_mismatch(tar_mid, tar_sec_f, path, shift, match_width);
			re = re / (1 + shape_com_mid + all_shape_complexity[index]);
			re2 = re2 / (1 + shape_com_mid + all_shape_complexity[index]);

			if (re < re2)
			{
				all_result.push_back(re);
				all_final.push_back(make_pair(index, false));
				//cout << re << endl;
			}
			else
			{
				all_result.push_back(re2);
				all_final.push_back(make_pair(index, true));
				//cout << re2 << endl;
			}
		}
		double temp;
		pair<int,bool> tempp;
		int all_size = all_result.size();
		for (int i = 0; i < all_size - 1; i++)
			for (int j = 0; j < all_size - 1 - i; j++)
				if (all_result[j] < all_result[j + 1])
				{
					temp = all_result[j];
					all_result[j] = all_result[j + 1];
					all_result[j + 1] = temp;
					tempp = all_final[j];
					all_final[j] = all_final[j + 1];
					all_final[j + 1] = tempp;
				}
		//cout << "the fianl order: " << endl;
		vector<pair<int, bool>> all_total_mid;
		for (int t = total_num - 1; t > total_num - 150; t--)
		{
			all_total_mid.push_back(all_final[t]);
			cout << "order: " << all_final[t].first << "  flip: " << all_final[t].second << " value: " << all_result[t] << " complxeity: " << all_shape_complexity[all_final[t].first] << endl;
		}
		all_final.swap(all_total_mid);
			
		return all_final;
	}


	vector<pair<int, bool>> Tiling_opt::quick_choose_TAR(vector<Point2f> inner_c) //得到选择出的pattern的序号和是否翻转的标志
	{
		int match_width = 4;
		vector<pair<int, bool>> all_total; //未翻转:false,翻转:true
		vector<double> all_result;
		int total_num = contour_dataset.size();
		vector<Point2f> contour_mid = inner_c;
		double shape_com_mid;
		vector<vector<double>> tar_mid = prototile_mid->compute_TAR(contour_mid, shape_com_mid);
		vector<int> cand_points_index = most_convex_p(contour_mid, curvature_com(contour_mid), 30); 
		//cout << "feature :" << cand_points_index.size()  << "total_num:  " << total_num<<endl;
		vector<vector<double>> tar_fea;
		for (int j = 0; j < cand_points_index.size(); j++)
		{
			tar_fea.push_back(tar_mid[cand_points_index[j]]);
		}

		//cout << "contour_mid: " << contour_mid.size() << "  tar_mid: " << tar_mid.size() << endl;
		for (int can_num = 0; can_num < total_num; can_num++)
		{
			//vector<vector<double>> tar_sec = all_con_tars[can_num];
			//vector<vector<double>> tar_sec_f = all_con_tars_flip[can_num];
			vector<vector<double>> tar_sec = all_fea_tars[can_num];
			vector<vector<double>> tar_sec_f = all_fea_tars_flip[can_num];
			//cout << "can_num: " << can_num <<  " tar_sec :" << tar_sec.size() << "   tar_sec_f: " << tar_sec_f.size() << endl;
			vector<pair<int, int>> path;
			int shift = 0;
			//double re = tar_mismatch(tar_mid, tar_sec, path, shift, match_width);
			//double re2 = tar_mismatch(tar_mid, tar_sec_f, path, shift, match_width);
			double re = tar_mismatch(tar_fea, tar_sec, path, shift, match_width);
			double re2 = tar_mismatch(tar_fea, tar_sec_f, path, shift, match_width);
			re = re / (1 + shape_com_mid + all_shape_complexity[can_num]);
			re2 = re2 / (1 + shape_com_mid + all_shape_complexity[can_num]);

			if (re < re2)
			{
				all_result.push_back(re);
				all_total.push_back(make_pair(can_num, false));
			}
			else
			{
				all_result.push_back(re2);
				all_total.push_back(make_pair(can_num, true));
			}
		}

		cout << "all_result.size :  " << all_result.size() << "all_total: " << all_total.size() << endl;
		double temp;
		pair<int, bool> tempp;
		int all_size = all_result.size();
		for (int i = 0; i < all_size - 1; i++)
			for (int j = 0; j < all_size - 1 - i; j++)
				if (all_result[j] < all_result[j + 1])
				{
					temp = all_result[j];
					all_result[j] = all_result[j + 1];
					all_result[j + 1] = temp;
					tempp = all_total[j];
					all_total[j] = all_total[j + 1];
					all_total[j + 1] = tempp;
				}
		vector<pair<int, bool>> all_total_mid;
		for (int t = all_total.size() - 1; t > total_num - 350; t--)
		{
			all_total_mid.push_back(all_total[t]);
			//cout << "order: " << all_total[t].first << "  flip: " << all_total[t].second << " value: " << all_result[t] << " complxeity: " << all_shape_complexity[all_total[t].first] << endl;
		}
		int midsize = all_total_mid.size();
		cout <<" midsize: "<< midsize << endl;
		all_total.swap(all_total_mid);
		return all_total;
	}

	vector<Point2f> Tiling_opt::morphing_tar(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<int> &mid_inter, vector<pair<int, int>> &path, int shift) //必须先用shift转换之后才能使path里的对应关系合法
	{
		vector<Point2f> final_pettern;
		vector<int> mid_inter_new;
		vector<Point2f> contour_mid;
		int c2size = contour2.size();
		for (int i = shift; i < shift + c2size; i++)
		{
			contour_mid.push_back(contour2[i % c2size]);
		}

		cout << contour1.size() << "   " << c2size << "  c3: " << contour_mid .size()<< endl;
		contour2 = contour_mid;
		double scale = arcLength(contour1, true) / arcLength(contour2, true);
		//cout << "scale: " << scale << endl;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] * scale;
		}
		Point2f cen1 = center_p(contour1);
		Point2f cen2 = center_p(contour2);
		Point2f shift2 = cen1 - cen2;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] + shift2;
		}

		Point2f v0 = contour1[0] - cen1;
		Point2f v1 = contour2[0] - cen1;
		double angle = acos(cos_two_vector(v0, v1))/PI*180;
		if (sin_two_vector(v0, v1) < 0) angle = -angle;
		Mat rot_mat(2, 3, CV_32FC1);
		rot_mat = getRotationMatrix2D(cen1, angle, 1);
		transform(contour2, contour_mid, rot_mat);
		contour2 = contour_mid;
		/*cout << "mid_inter_inner: " << mid_inter.size() << endl;
		for (int i = 0; i < 4; i++)
			cout << mid_inter[i] << "  ";*/
		int psize = path.size();
		int t = 0;
		for (int i = 0; i < psize; i++)
		{
			int first = path[i].first;
			int sec = path[i].second;
			Point2f fin = 0.5 * contour1[first] + 0.5 * contour2[sec];;
			if (t<4 && first == mid_inter[t])
			{
				//cout << t << " " << first << "  " << i << endl;
				t++;
				mid_inter_new.push_back(i);
			} 
			final_pettern.push_back(fin);

			/*Point2f fin;
			if (t<4 && first == mid_inter[t])
			{
				t++;
				mid_inter_new.push_back(i);
				fin = 0.95 * contour1[first] + 0.05 * contour2[sec];
			}
			else fin = 0.5 * contour1[first] + 0.5 * contour2[sec];
			final_pettern.push_back(fin);*/
			
		}
		Point2f shiftting = Point2f(400, 400) - cen1;
		Point2f shiftting2 = Point2f(1200, 400) - cen1;
		Mat ta = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat tttt = Mat(800, 3200, CV_8UC3, Scalar(255, 255, 255));
		//draw_poly(tttt, contour1, Point2f(400, 400), 3);
		//draw_poly(tttt, contour2, Point2f(1200, 400), 4);
		//draw_poly(tttt, final_pettern, Point2f(2000, 400), 9);
		for (int t = 0; t < contour1.size(); t++)
		{
			circle(tt, contour1[t] + shiftting, 3, Scalar(100, 230, 0), -1);
			circle(ttt, contour1[t] + shiftting, 3, Scalar(100, 230, 0), -1);
			circle(ta, contour1[t] + shiftting, 3, Scalar(100, 230, 0), -1);
			MyLine(tttt, contour1[t] + shiftting, contour1[(t + 1) % contour1.size()] + shiftting,"black");
		}
		for (int t = 0; t < contour2.size(); t++)
		{
			circle(tt, contour2[t] + shiftting, 3, Scalar(0, 255, 0), -1);
			circle(ttt, contour2[t] + shiftting, 3, Scalar(0, 255, 0), -1);
			circle(ta, contour2[t] + shiftting2, 3, Scalar(0, 255, 0), -1);
			MyLine(tttt, contour2[t] + shiftting2, contour2[(t + 1) % contour2.size()] + shiftting2, "black");
		}
		for (int t = 0; t < psize; t++)
		{
			//circle(ttt, src1_points[t], 4, Scalar(255, 120, 120), -1);
			//circle(ta, final_pettern[t] + shiftting, 3, Scalar(255, 120, 120), -1);
			MyLine(tt, contour1[path[t].first] + shiftting, contour2[path[t].second] + shiftting, "grey");
		}
		/*for (int t = 0; t < 20; t++)
		{
			MyLine(tt, contour1[path[t].first], contour2[path[t].second], "grey");
		}*/
		mid_inter_new = joint_relocate(final_pettern, mid_inter_new, 1);
		mid_inter = mid_inter_new;
		final_pettern = sampling(final_pettern, 2);
		Point2f shiftting3 = Point2f(2000, 400) - center_p(final_pettern);
		for (int i = 0; i < final_pettern.size(); i++)
		{
			//circle(tttt, final_pettern[i] + shiftting, 3, Scalar(0, 230, 130), -1);
			MyLine(tttt, final_pettern[i] + shiftting3, final_pettern[(i + 1) % final_pettern.size()] + shiftting3, "black");
			circle(ttt, final_pettern[i] + shiftting, 3, Scalar(0, 230, 130), -1);
		}
		//for (int i = 0; i < path.size(); i++)
		//{
			//cout << path[i].first << "   " << path[i].second << endl;
			//MyLine(ta, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
			//MyLine(tt, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
		//}
		//imshow("??:", ta);
		//imshow("???:", tt);
		//imshow("????:", ttt);
		
		imwrite("D:\\mor_all.png", tttt);
		return final_pettern;
	}

	void Tiling_opt::contour_fine_tuning(vector<Point2f> &contour_, int first, int second)
	{
		vector<Point2f> mid_con;
		int csize = contour_.size();
		vector<Point2f> pre;
		vector<Point2f> post;
		cout << csize<<" "<< first << " " << second << endl;
		if (second - first < 4)
		{
			Point2f tem = contour_[first + 1];
			contour_[first + 1] = contour_[second];
			contour_[second] = tem;
		}
		else 
		{
			int mid = (first + second) / 2;
			for (int t = first; t < mid; t++)
			{
				pre.push_back(contour_[t]);
			}
			for (int t = mid; t <= second + 1; t++)
			{
				post.push_back(contour_[t]);
			}
			for (int t = 0; t < first; t++)
			{
				mid_con.push_back(contour_[t]);
			}
			if (mid_con.empty()) mid_con.push_back(0.5*contour_[csize - 1] + 0.5*post[post.size() - 2]);
			else mid_con.push_back(0.5*mid_con.back() + 0.5*post[post.size() - 2]);
			for (int t = post.size() - 2; t >= 0; t--)
			{
				mid_con.push_back(post[t]);
			}
			for (int t = pre.size() - 1; t > 0; t--)
			{
				mid_con.push_back(pre[t]);
			}
			mid_con.push_back(0.5*mid_con.back() + 0.5*contour_[(second + 2) % csize]);
			for (int t = second + 2; t < csize; t++)
			{
				mid_con.push_back(contour_[t]);
			}
			contour_.swap(mid_con);
		}		
	}

	vector<int> Tiling_opt::joint_relocate(vector<Point2f> contour_, vector<int> joint_index, int num_c) //将原始轮廓上的划分点对应到采样后的轮廓上
	{
		vector<Point2f> contour_mid = sampling(contour_, num_c + 1);//选择最少的点进行比较
		int con_size = contour_mid.size();
		vector<int> mid_in;
		for (int i = 0; i <joint_index.size(); i++)
		{
			int index_ = 0;
			double dist = 10000;
			for (int t = 0; t <con_size; t++)
			{
				if (length_two_point2f(contour_[joint_index[i]], contour_mid[t]) < dist)
				{
					dist = length_two_point2f(contour_[joint_index[i]], contour_mid[t]);
					index_ = t;
				}
			}
			mid_in.push_back(index_);
		}
		return mid_in;
	}
		

	double Tiling_opt::tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width) //点对应匹配的筛选框宽度
	{
		int first_num = first_arr.size();
		int second_num = second_arr.size(); //first 作为y轴 ,second为x轴
		//cout <<  first_num << " - " << second_num << endl;
		//width = 30;
		if (first_num != second_num)
		{
			cout <<"The sampling points of two contours are not equal: "<< first_num<<" - "<<second_num<< endl;
			//return 0;
		}
		if (first_arr[0].size() != second_arr[0].size())
		{
			cout << "The tar num of each point is not equal: " << first_arr[0].size() << " - " << second_arr[0].size()<< endl;
			return 0;
		}
		double min_mis = 10000;
		vector<pair<int, int>> path_min;
		//double distance[202][202];
		//int step[202][202];//记录总的步数
		//cout << "first: " << first_num << "second: " << second_num << endl;
		for (int shift = 0; shift < second_num; shift++) //将first固定，分别对齐second的起点
		{
			//int ccc = 0;
			//int ddd = 0;
			for (int i = 0; i < first_num; i++)
			{
				for (int j = 0; j < second_num; j++)
				{
					//ddd++;
					dis[i][j] = 0;
					if (max(0, i - width) <= j && j <= min(second_num - 1, i + width))
					{
						//ccc++;
						distance[i][j] = 0;
					}
					else distance[i][j] = 100000;

				}
			}
			//cout << "ccc :  " << ccc <<"ddd: "<<ddd<< endl;
			//distance[0][0]记录的是对齐的第一个点
			dis[0][0] = length_two_point_tar(first_arr[0], second_arr[shift]);//
			distance[0][0] = dis[0][0];

			for (int i = 1; i < first_num; i++)
			{
				if (distance[i][0] == 100000) continue;
				dis[i][0] = length_two_point_tar(first_arr[i], second_arr[shift]);
				distance[i][0] = distance[i - 1][0] + dis[i][0];
			}
			for (int i = 1; i < second_num; i++)
			{
				if (distance[0][i] == 100000) continue;
				dis[0][i] = length_two_point_tar(first_arr[0], second_arr[(i + shift)%second_num]);
				distance[0][i] = distance[0][i - 1] + dis[0][i];
			}
			
			for (int i = 1; i < first_num; i++)
			{
				for (int j = 1; j < second_num; j++)
					//(int i = istart; i <= imax; i++)
				{
					if (distance[i][j] == 100000) continue;
					dis[i][j] = length_two_point_tar(first_arr[i], second_arr[(j + shift) % second_num]);
					double g1 = distance[i - 1][j] + dis[i][j];
					double g2 = distance[i - 1][j - 1] + dis[i][j];
					double g3 = distance[i][j - 1] + dis[i][j];
					if (g1 < g2)
					{
						if (g1 < g3) distance[i][j] = g1;
						else distance[i][j] = g3;
					}
					else
					{
						if (g2 < g3) distance[i][j] = g2;
						else distance[i][j] = g3;
					}
				}
			}
			//for (int i = 0; i < first_num; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < second_num; j++)
			//		//(int i = istart; i <= imax; i++)
			//	{
			//		if (distance[i][j] == 100000) cout <<  "  0   ";
			//		else cout << distance[i][j] << "  ";
			//	}
			//}
			//cout << "shift " << shift << "  distance[first_num - 1][second_num - 1]:  " << distance[first_num - 1][second_num - 1] << endl;
			if (distance[first_num - 1][second_num - 1] < min_mis)
			{
				path_min.swap(vector<pair<int, int>>());
				min_mis = distance[first_num - 1][second_num - 1];
				print_TAR_Path(dis, distance, first_num - 1, second_num - 1, path_min);
				sec_shift = shift;
			}
		}
		path = path_min;
		return min_mis;
	}

	void Tiling_opt::print_TAR_Path(double d[][202], double dp[][202], int i, int j, vector<pair<int, int>>& path)
	{
		if (i == 0 && j == 0) {
			//cout << first_arr[i] << " - " << second_arr[j] << endl;
			path.push_back(make_pair(i, j));
			return;
		}

		if (abs(dp[i][j] - (dp[i - 1][j - 1] + dis[i][j])) < 0.001){
			print_TAR_Path(d, dp, i - 1, j - 1, path);

		}
		else if (abs(dp[i][j] - (dp[i][j - 1] + dis[i][j])) < 0.001) {
			print_TAR_Path(d, dp, i, j - 1, path);

		}
		else {
			print_TAR_Path(d, dp, i - 1, j, path);
		}
		path.push_back(make_pair(i, j));
	}

	double Tiling_opt::evalua_deformation(vector<Point2f> contour1, vector<Point2f> contour2)
	{
		//储存顺序为 mid sec result
		
		double total_score = 0;
		prototile_second->Pro_clear();
		prototile_second->loadPoints(contour1);
		double shape_com1 = 0;
		vector<vector<double>> tar_all = prototile_second->compute_TAR(prototile_second->contour_sample[1], shape_com1);//(num_c+1)*100 points
		prototile_tem->Pro_clear();
		prototile_tem->loadPoints(contour2);
		double shape_com2;
		vector<vector<double>> tar_all1 = prototile_tem->compute_TAR(prototile_tem->contour_sample[1], shape_com2);
		vector<vector<double>> tar_all2 = prototile_tem->compute_TAR(prototile_tem->contour_sample_flip[1], shape_com2);
		vector<pair<int, int>> path;
		vector<pair<int, int>> path2;
		int shift = 0;
		int shift2 = 0;
		int Lab = 0;
		vector<Point2f> contour_1 = prototile_second->contour_sample[1];
		vector<Point2f> contour_2 = prototile_tem->contour_sample[1];
		//计算出来最匹配的路径	
		double re = tar_mismatch(tar_all, tar_all1, path, shift);
		double re2 = tar_mismatch(tar_all, tar_all2, path2, shift2);
		if (re > re2)
		{
			re = re2;
			shift = shift2;
			path = path2;
			contour_2 = prototile_tem->contour_sample_flip[1];
		}
		Lab = path.size();
		
		//tar 对比得分
		//total_score = re / (1 + shape_com1 + shape_com2);
		total_score = 1 - (re / (1 + shape_com1 + shape_com2)) / (2 * Lab);
		//将contour_2进行变换，缩放及旋转
		vector<Point2f> contour_mid;
		int c2size = contour_2.size();
		for (int i = shift; i < shift + c2size; i++)
		{
			contour_mid.push_back(contour_2[i % c2size]);
		}
		//cout << contour_1.size() << "   " << c2size << "  c3: " << contour_mid.size() << endl;
		contour_2 = contour_mid;
		double scale = arcLength(contour_1, true) / arcLength(contour_2, true);

		for (int i = 0; i < c2size; i++)
		{
			contour_2[i] = contour_2[i] * scale;
		}
		cout << "shift:  " << shift << "  Lab: " << Lab<<"   scale : " << scale << endl;
		
		Point2f cen1 = center_p(contour_1);
		Point2f cen2 = center_p(contour_2);
		Point2f shift_p = cen1 - cen2;
		for (int i = 0; i < c2size; i++)
		{
			contour_2[i] = contour_2[i] + shift_p;
		}
		Point2f v0 = contour_1[0] - cen1;
		Point2f v1 = contour_2[0] - cen1;
		double angle = acos(cos_two_vector(v0, v1)) / PI * 180;
		if (sin_two_vector(v0, v1) < 0) angle = -angle;
		Mat rot_mat(2, 3, CV_32FC1);
		rot_mat = getRotationMatrix2D(cen1, angle, 1);
		transform(contour_2, contour_mid, rot_mat);
		contour_2 = contour_mid;

		//使用像素级的方法计算两个形状的面积差
		/*vector<Point2f> bbx1 = b_box(contour_1);
		vector<Point2f> bbx2 = b_box(contour_2);
		cout << bbx1[0] << "  " << bbx1[2];
		cout << bbx2[0] << "  " << bbx2[2];
		int raw = (abs(bbx1[0].y - bbx1[1].y) > abs(bbx2[0].y - bbx2[1].y)) ? abs(bbx1[0].y - bbx1[1].y) : abs(bbx2[0].y - bbx2[1].y);
		int col = (abs(bbx1[1].x - bbx1[2].x) > abs(bbx2[1].x - bbx2[2].x)) ? abs(bbx1[1].x - bbx1[2].x) : abs(bbx2[1].x - bbx2[2].x);
		cout << raw << "  " << col;*/
		int raw = 1000;
		int col = 1000;
		Mat drawing_1 = Mat(raw, col, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_2 = Mat(raw, col, CV_8UC3, Scalar(255, 255, 255));
		Mat drawing_3 = Mat(raw, col, CV_8UC3, Scalar(255, 255, 255));
		draw_poly(drawing_1, contour_1, Point2f(col / 2, raw / 2), 0);
		draw_poly(drawing_2, contour_2, Point2f(col / 2, raw / 2), 0);
		imshow("poly1", drawing_1);
		imshow("poly2", drawing_2);
		cvtColor(drawing_1, drawing_1, COLOR_BGR2GRAY);
		threshold(drawing_1, drawing_1, 128, 1, cv::THRESH_BINARY);
		cvtColor(drawing_2, drawing_2, COLOR_BGR2GRAY);
		threshold(drawing_2, drawing_2, 128, 1, cv::THRESH_BINARY);
		int poly1 = 0;
		int poly2 = 0;
		int poly_ = 0;
		for (int i = 0; i < col; i++)
		{
			for (int j = 0; j < raw; j++)
			{
				if ((int)drawing_1.at<uchar>(i, j) ==0 && (int)drawing_2.at<uchar>(i, j) ==0)
				{
					poly_++;
					//drawing_3.at<uchar>
				}
				else if ((int)drawing_1.at<uchar>(i, j) == 0)
				{
					poly1++;
				}
				else if ((int)drawing_2.at<uchar>(i, j) ==0)
				{
					poly2++;
				}
			}
		}
		double area_score = 1 - (poly1 + poly2) / poly_;

		cout << "dianshu:  " << poly1 << "    " << poly2 << "    " << poly_ << endl;
		cout << "total_score: " << total_score << " area_score" << area_score << endl;

		
		//score_mid_r = contourArea(contour[2]) / contourArea(contour[0]);
		//score_sec_r = contourArea(contour[1]) / contourArea(contour[0]);
		return total_score;
	}

	vector<Point2f> Tiling_opt::simulation_mid(string imaname, int inner_one, int cand_one)
	{
		int num_c = 1;//选择(num_c+1)*100个点
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		vector<Point2f> contour_ = prototile_first->contour;
		//vector<Point2f> contour_ = prototile_first->contour;
		load_dataset();
		string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result\\simula\\" + prototile_first->contourname;
		const char *na = rootname.c_str();
		mkdir(na);
		int ppindex = p_p_index.size();
		int margin = prototile_first->contour.size() / 10;
		cout << "margin: " << margin << endl;
		int count = 0;

		vector<vector<Point2f>> inner_conts;
		vector<vector<int>> all_situation_index;
		vector<vector<int>> mid_interval_index;
		vector<int> result_index;
		Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				if (abs(p_p_index[j] - p_p_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(p_p_index[m] - p_p_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						if (abs(p_p_index[n] - p_p_index[m]) < margin) continue;
						vector<Point2f> inner_contour;
						result_index.swap(vector<int>());
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						result_index.push_back(p_p_index[i]);
						result_index.push_back(p_p_index[j]);
						result_index.push_back(p_p_index[m]);
						result_index.push_back(p_p_index[n]);
						//cout << "  i: " << p_p_index[i]
						//	<< "   j: " << p_p_index[j]
						//	<< "   m: " << p_p_index[m]
						//	<< "   n: " << p_p_index[n] << endl;	
						////one_situ_div(result_index, contour_, inner_contour);					
						//string image = int2string(count);		
						Mat drawing_1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						if (!one_situ_div(result_index, contour_, inner_contour, mid_interval, drawing_1))
						{
							//cout << "-------------collision-------------" << endl;
							continue;
						}

						count++;
						cout << endl << count << " succeed";
						cout << "    inner_contour.size : " << inner_contour.size() << endl;

						inner_conts.push_back(inner_contour);
						all_situation_index.push_back(result_index);
						mid_interval_index.push_back(mid_interval);
						drawing1 = drawing_1;
						if ((count - inner_one) == 1) break;

					}
					if ((count - inner_one) == 1) break;
				}
				if ((count - inner_one) == 1) break;
			}
			if ((count - inner_one) == 1) break;
		}
		cout << "succeed count: " << count << endl;
		// search the right image
		/*
		1.将周长调整为一致
		2.搜索最相近的图案以及角度
		3.变形
		*/
		if (count == 0)
		{
			cout << "no right placement" << endl;
			return vector<Point2f>();
		}
		Point2f shift2 = Point2f(400, 400) - center_p(contour_);
		for (int j = 0; j < contour_.size(); j++)
		{
			circle(drawing1, contour_[j] + shift2, 1, Scalar(0, 0, 0), -1);

			//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		}
		for (int j = 0; j < p_p_index.size(); j++)
		{
			circle(drawing1, contour_[p_p_index[j]] + shift2, 4, Scalar(0, 0, 255), -1);
		}
		circle(drawing1, contour_[result_index[0]] + shift2, 8, Scalar(0, 255, 0), -1);
		circle(drawing1, contour_[result_index[1]] + shift2, 8, Scalar(0, 255, 0), -1);
		circle(drawing1, contour_[result_index[2]] + shift2, 8, Scalar(0, 255, 0), -1);
		circle(drawing1, contour_[result_index[3]] + shift2, 8, Scalar(0, 255, 0), -1);
		string conut_name = rootname + "\\PlacingResult_" + int2string(inner_one) + ".png";
		imwrite(conut_name, drawing1);

		/*cout << "inner_one: " << inner_one  << endl;
		string filepathname = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test" + int2string(inner_one)+".txt";
		vector<Point> con;
		for (int i = 0; i < inner_conts[inner_one].size(); i++)
		con.push_back((Point)inner_conts[inner_one][i]);
		fileout(filepathname, con);*/


		vector<CandPat> candida_contours;

		//ofstream out("D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test.txt");
		//if (out.is_open())
		//{
		//	out << inner_conts[i].size() + 1 << endl;//contours[0].size()
		//	for (int j = 0; j < inner_conts[i].size(); j++)
		//		out << (int)inner_conts[i][j].x << "," << (int)inner_conts[i][j].y << endl;
		//	out << (int)inner_conts[i][0].x << "," << (int)inner_conts[i][0].y << endl;  //首尾连起来
		//}
		//cout << "contours[0].size(): " << i << endl;
		//out.close();
		cout << inner_conts.size() << "  :  " << inner_one << endl;
		candida_contours = compare_shapes(inner_conts[inner_one], num_c);
		//cout << "candida_contours" << candida_contours.size()<< endl;

		vector<int> mid_inter = joint_relocate(inner_conts[inner_one], mid_interval_index[inner_one], num_c);
		//将所有的结果保存下来
		Mat drawing_pro = Mat(800, 2400, CV_8UC3, Scalar(255, 255, 255));
		CandPat tem = candida_contours[cand_one];
		prototile_second->Pro_clear();
		prototile_second->loadPoints(contour_dataset[tem.number]);
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_conts[inner_one]);

		vector<Point2f> contour_inner = prototile_mid->contour_sample[num_c]; //选择最少的点进行比较
		vector<double> contour_inner_c = curvature_com(contour_inner);// prototile_mid->contour_curva[0];
		vector<Point2f> contour_cand = CandP2Contour(tem, num_c);
		vector<double> contour_cand_c = curvature_com(contour_cand);

		//inner and cand morph into the final pettern	
		//!!!!!!在下一步中将morph改为分段!!!此时没有考虑中心平移和缩放倍数，因为之前是一一对应进行的变形
		int num = 1;
		float ratio = 2;
		while (num-- != 0)
		{
			ratio = ratio + 1;
			vector<Point2f> inter_mid = morphing_2_patterns(contour_inner, contour_cand, contour_inner_c, contour_cand_c, mid_inter, ratio / 10);
			//show the result
			draw_poly(drawing_pro, contour_inner, Point2f(400, 400));
			draw_poly(drawing_pro, contour_cand, Point2f(1200, 400));
			draw_poly(drawing_pro, inter_mid, Point2f(2000, 400));

			Mat drawing_mid = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
			Mat drawing_mA = Mat(1600, 2400, CV_8UC3, Scalar(255, 255, 255));

			draw_allplane(drawing_mid, inter_mid, mid_inter, 0.4);

			//将该proto1以及相邻四个proto2展示出来
			Point2f cente = center_p(inter_mid);
			vector<int> return_p;
			vector<vector<Point2f>> four_;
			vector<Point2f> morphed_A = extract_contour(inter_mid, mid_inter, return_p, four_,0);
			//evalua_deformation();
			draw_allplane(drawing_mA, morphed_A, return_p, 0.4);

			string filename = rootname + "\\0.";
			string file2 = filename + char(ratio + 48) + "_Cand_" + int2string(cand_one) + "tiling.png";
			filename = filename + char(ratio + 48) + "_Candidate_" + int2string(cand_one) + ".png";

			imwrite(filename, drawing_pro);
			imwrite(file2, drawing_mA);
			imshow("result_mid_show: ", drawing_mid);

		}
		return contour_inner;
	}


	void Tiling_opt::pattern_joint(jointPat pattern)
	{
		double pattern_width = length_two_point2f(pattern.four_contour[0][pattern.interval[0]], pattern.four_contour[0][pattern.interval[2]]);
		double pattern_length = length_two_point2f(pattern.four_contour[0][pattern.interval[1]], pattern.four_contour[0][pattern.interval[3]]);
		//实际中把较长的一边设为5cm
		double W = 0.2 * pattern_width; //对应的实际长度为0.2
		double L = 0.2 * pattern_length;

		double safe_thickness_l = 0.2 * L;
		double safe_thickness_w = 0.2 * W;
		double safe_margin = 0.1 * L;

		Point2f dir13 = unit_vec(pattern.four_contour[0][pattern.interval[2]] - pattern.four_contour[0][pattern.interval[0]]);//第一个图案的1,3向量
		Point2f dir24 = unit_vec(pattern.four_contour[0][pattern.interval[3]] - pattern.four_contour[0][pattern.interval[1]]);//第一个图案的2,4向量

		//  0:translation  1:rotation  2:flipping(1-3)   3:flipping(2-4)		
		if (pattern.type == 0)
		{
		}
		else if (pattern.type == 1)
		{
		}
		else if (pattern.type == 2)
		{
		}
		else if (pattern.type == 3)  // 3时，关节位于1,3pattern的2,4交点
		{
			vector<Point2f> insert;
			Point2f intersect_p = pattern.four_contour[0][pattern.interval[3]];
			Point2f endp1 = intersect_p + 0.5 * W * dir13;
			Point2f endp2 = intersect_p - 0.5 * W * dir13;
			Point2f mid1 = intersect_p + 0.25 * W * dir13;
			Point2f mid2 = intersect_p - 0.25 * W * dir13;
			Line_Seg cut_line_mid = Line_Seg(mid1, mid2);
			Line_Seg cut_line = Line_Seg(endp1, endp2);

			//找到两个正向向量
			Point2f v_vec = vertical_vec(dir13);
			if (cos_two_vector(v_vec, dir24) < 0) v_vec = -v_vec;

			vector<Point2f> all_intersect = line_polygon(cut_line, pattern.four_contour[0]);
			int intersize = all_intersect.size();
			if (intersize == 1)
			{
				vector<Point2f> all1 = line_polygon(Line_Seg(mid1 + safe_thickness_l*v_vec, mid1 - safe_thickness_l*v_vec), pattern.four_contour[0]);

				vector<Point2f> all2 = line_polygon(Line_Seg(mid2 + safe_thickness_l*v_vec, mid2 - safe_thickness_l*v_vec), pattern.four_contour[0]);
				if (all1.empty() && all2.empty())
				{

				}
			}


		}
	}

	vector<Point2f> Tiling_opt::construct_joint(jointPat pattern, int &mid)
	{
		int n = 3;
		int margin = 0;
		vector<Point2f> origin = pattern.four_contour[2];
		int size_o = origin.size();
		int start = pattern.interval[0];
		int end = pattern.interval[2];
		Point2f vec13 = origin[end] - origin[start]; //3-1

		vector<Point2f> new_contour;
		vector<Point2f> new_contour2;
		vector<Point2f> shift_contour = origin;
		new_contour.push_back(origin[start]);
		for (int m = 1; m < margin; m++)
		{
			new_contour.push_back(origin[start + m]);
			new_contour2.push_back(origin[(start + size_o - m) % size_o]);
		}
		mid = margin - 1;
		for (int i = 0; i < n; i++)
		{
			if (i > 0)
			{
				shift_contour.swap(vector<Point2f>());
				Point2f shift = i*vec13;
				for (int d = 0; d < size_o; d++)
				{
					shift_contour.push_back(origin[d] + shift);
				}
			}
			for (int j = start + margin; j <= end - margin; j++)
			{
				new_contour.push_back(shift_contour[j]);
				mid++;
			}
			for (int j = size_o + start - margin; j >= end + margin; j--)
			{
				new_contour2.push_back(shift_contour[j%size_o]);
			}
		}
		for (int t = 0; t < 2 * margin - 1; t++)
		{
			new_contour.push_back(shift_contour[end - 2 + t]);
		}
		mid = mid + margin;
		for (int i = new_contour2.size(); i > 0; i--)
		{
			new_contour.push_back(new_contour2[i - 1]);
		}
		return new_contour;
	}


	vector<CandPat> Tiling_opt::compare_shapes(vector<Point2f> inner_c, int num_c)
	{
		vector<int> order_total;
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_c);
		vector<Point2f> contour_mid = prototile_mid->contour_sample[num_c]; //选择最少的点进行比较
		cout << "contour_mid.size()" << contour_mid.size() << endl;
		vector<double> contour_mid_c = curvature_com(contour_mid);

		int total_num = contour_dataset.size();
		cout << "total_num:" << total_num << endl;
		vector<vector<double>> score_3types(3);
		vector<vector<int>> index_s(3);
		//int internum = 0;
		for (int can_num = 0; can_num < total_num; can_num++)
		{
			prototile_second->Pro_clear();
			prototile_second->loadPoints(contour_dataset[can_num]);
			vector<Point2f> contour_second = prototile_second->contour_sample[num_c];
			for (int method_ = 1; method_ <= 3; method_++)
			{
				double score;
				score = matchShapes(contour_mid, contour_second, method_, 0);
				score_3types[method_ - 1].push_back(score);
				index_s[method_ - 1].push_back(can_num);
			}

		}

		for (int i = 0; i < 3; i++)
		{
			sort_comb(score_3types[i], index_s[i]);
			//cout
			/*cout << "3zhongpingfen" << endl;
			for (int t = index_s[i].size() - 1; t > total_num-20; t--)
			{
			cout << index_s[i][t] << " " << score_3types[i][index_s[i][t]] << endl;
			}	*/
			//添加到order数组里
			for (int t = index_s[i].size() - 1; t > total_num - 15; t--)
			{
				int f = 0;
				for (int j = 0; j < order_total.size(); j++)
				{
					if (index_s[i][t] == order_total[j])
					{
						f = 1;
						break;
					}
				}
				if (f == 0) order_total.push_back(index_s[i][t]);
			}
		}
		vector<CandPat> final_order;
		vector<CandPat> score_total;
		int ordersize = order_total.size();
		//cout 合并后的候选项
		cout << "ordersize: " << ordersize << endl;
		for (int i = 0; i < ordersize; i++)
		{
			cout << "order_total: " << order_total[i] << endl;
		}

		for (int t = 0; t < ordersize; t++)  //展示所有order里的候选图案的对比结果,筛选出用min_mismatch值最小的结果
		{
			prototile_second->Pro_clear();
			prototile_second->loadPoints(contour_dataset[order_total[t]]);
			vector<Point2f> contour_second = prototile_second->contour_sample[num_c];
			vector<double> contour_second_c = curvature_com(contour_second); //prototile_second->contour_curva[num_c];
			vector<Point2f> contour_sec_f = prototile_second->contour_sample_flip[num_c];
			vector<double> contour_sec_c_f = curvature_com(contour_sec_f); //prototile_second->contour_curva_flip[num_c];
			//cout << "t: " <<t<< endl;
			CandPat right_mis = min_mismatch(contour_mid, contour_second, contour_mid_c, contour_second_c, order_total[t], false);
			CandPat flip_mis = min_mismatch(contour_mid, contour_sec_f, contour_mid_c, contour_sec_c_f, order_total[t], true);
			if (right_mis.mismatch < flip_mis.mismatch)
			{
				score_total.push_back(right_mis);
			}
			else
			{
				score_total.push_back(flip_mis);
			}
			//cout << "youzale???" << endl;
		}

		CandPat temp;
		for (int i = 0; i < score_total.size() - 1; i++)
		{
			for (int j = 0; j < score_total.size() - 1 - i; j++)
				if (score_total[j].mismatch > score_total[j + 1].mismatch)
				{
					temp = score_total[j];
					score_total[j] = score_total[j + 1];
					score_total[j + 1] = temp;

				}
		}
		//cout << "mismatch order" << endl;
		//for (int i = 0; i < score_total.size(); i++)
		//{
		//	cout << score_total[i].number << "  " << score_total[i].mismatch << "  " << endl;
		//}
		cout << "the top ten" << endl;
		for (int i = 0; i < 10; i++)
		{
			final_order.push_back(score_total[i]);
			cout << final_order[i].number << "   " << final_order[i].mismatch << endl;
		}

		return final_order;
	}

	CandPat Tiling_opt::min_mismatch(vector<Point2f> inner, vector<Point2f> cand, vector<double> inner_c, vector<double> cand_c, int theone, bool isFilp)
	{
		//将两个轮廓的周长质心对齐
		// 这里要保证inner和cand的size差不多，因为会出现截尾现象，如果差距太大会造成较大误差
		//cout << "1: " << cand.size() << "  2: " << cand_f.size()<<endl;
		//cout << "c: " << center_p(cand) << "  d: " << center_p(cand_f) << endl;
		double scale = arcLength(inner, true) / arcLength(cand, true);
		//cout << "scale: " << scale << endl;
		for (int i = 0; i < cand.size(); i++)
		{
			cand[i].x = cand[i].x * scale;
			cand[i].y = cand[i].y * scale;
		}
		Point2f Ccen = center_p(inner);
		Point2f shift = Ccen - center_p(cand);
		for (int i = 0; i < cand.size(); i++)
		{
			cand[i] = cand[i] + shift;
		}
		//对cand进行360度的旋转
		int total_num = inner.size() < cand.size() ? inner.size() : cand.size();
		int each_num = 201;
		vector<Point2f> cand_tem;
		Mat rot_mat(2, 3, CV_32FC1);
		vector<Point2f> test1;
		vector<double> test11;
		vector<Point2f> test2;
		vector<double> test22;

		double min_mismatch = 100000;
		int min_angle = 0;
		int min_index = 0;

		for (int angle = 0; angle < 360; angle = angle + 5)
		{
			if (angle != 0)
			{
				rot_mat = getRotationMatrix2D(Ccen, angle, 1);
				transform(cand, cand_tem, rot_mat);
			}
			else cand_tem = cand;
			vector<int> cand_n;
			cand_n = search_align_p(Ccen, inner[0], cand_tem);
			int cand_n_s = cand_n.size();
			int cand_tem_size = cand_tem.size();

			//inner质点和第一个点的连接线与cand交点的个数
			//cout << "num:_can" << cand_n_s << endl;
			//for (int i = 0; i < cand_n_s; i++)
			//{
			//	cout << cand_n[i] << endl;
			//}
			//
			for (int t = 0; t < cand_n_s; t++)
			{
				int num_now = 0;
				double accumu_mis = 0;
				while (num_now < total_num)
				{
					test1.swap(vector<Point2f>());
					test11.swap(vector<double>());
					test2.swap(vector<Point2f>());
					test22.swap(vector<double>());
					//cout << "num_now: " << num_now << endl;

					if (total_num - num_now > each_num)
					{
						int mid_now = num_now + each_num;
						for (int i = num_now; i < mid_now; ++i)
						{
							++num_now;
							test1.push_back(inner[i]);
							test11.push_back(inner_c[i]);
							//cout << "cand_n[t]: " << i<< endl;
							test2.push_back(cand_tem[(i + cand_n[t]) % cand_tem_size]);
							test22.push_back(cand_c[(i + cand_n[t]) % cand_tem_size]);
						}
					}
					else
					{
						for (int i = num_now; i < total_num; i++)
						{
							++num_now;
							test1.push_back(inner[i]);
							test11.push_back(inner_c[i]);
							test2.push_back(cand_tem[(i + cand_n[t]) % cand_tem_size]);
							test22.push_back(cand_c[(i + cand_n[t]) % cand_tem_size]);
						}
					}
					//cout << "test1:" << test11.size()
					//	<< "test2:" << test22.size() << endl;
					vector<pair<int, int>> dppath;
					accumu_mis += quadr_mismatch(test1, test2, test11, test22, dppath); //因为quadr函数里用的数组是100x100，所以需要截取输入
				}
				//cout << "angle  "<<angle<<" mismatch: "<<accumu_mis << endl;
				if (accumu_mis < min_mismatch)
				{
					min_angle = angle;
					min_mismatch = accumu_mis;
					min_index = cand_n[t];
				}

			}
		}
		//cout << min_mismatch << endl;
		CandPat min_pat = { theone, isFilp, min_angle, min_index, min_mismatch };
		return min_pat;
	}

	double Tiling_opt::quadr_mismatch(vector<Point2f> first_arr, vector<Point2f> second_arr, vector<double> first_c, vector<double> second_c, vector<pair<int, int>>& path, double zeta) //每次只能比较100个点以内
	{
		//ofstream out("D:\\VisualStudioProjects\\manual_record\\dtw.txt");
		int first_num = first_arr.size();
		int second_num = second_arr.size();

		//cout << "\n        first.size: " << first_c[50] << "  -----    chararr_size" << second_c[50] << endl;

		//double dis[202][202];//两组点之间的坐标差异
		double max_dis = 0;
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis[i][j] = length_two_point2f(first_arr[i], second_arr[j]);
				if (dis[i][j] > max_dis) max_dis = dis[i][j];
			}
		}
		//坐标差值归一化
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis[i][j] = dis[i][j] / max_dis;
				//cout << "dis[i][j]: " << dis[i][j] << endl;
			}
		}
		//for (int i = 0; i < 10; i++)
		//{
		//	for (int j = 0; j < 10; j++)
		//	{
		//		cout << "dis[i][j]: " << dis[i][j] << endl;
		//	}
		//}

		//double dis_cur[202][202];//两组点之间的曲率差异
		double max_dis_cur = 0;
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis_cur[i][j] = cur_length_two_p(first_c[i], second_c[j]);

				//cout << "dis_cur[i][j]: " << dis_cur[i][j] << endl;
				if (dis_cur[i][j]>max_dis_cur)
				{
					//cout << "hahah:" << endl;
					max_dis_cur = dis_cur[i][j];


				}
			}
		}
		//曲率差值归一化
		//cout << "max_dis_cur: " << max_dis_cur << endl;
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis_cur[i][j] = zeta * dis_cur[i][j] / max_dis_cur;
				//cout << "dis_cur[i][j]: " << dis_cur[i][j] << endl;
			}
		}
		//double distance[202][202];
		//int step[202][202];//记录总的步数
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				distance[i][j] = 0;
				//step[i][j] = 0;
			}
		}
		distance[0][0] = (double)2 * (dis[0][0] + dis_cur[0][0]);
		//step[0][0] = 0;

		for (int i = 1; i < first_num; i++)
		{
			distance[i][0] = distance[i - 1][0] + (dis[i][0] + dis_cur[i][0]);
			//step[i][0] = step[i - 1][0] + 1;
		}
		for (int i = 1; i < second_num; i++)
		{
			distance[0][i] = distance[0][i - 1] + (dis[0][i] + dis_cur[0][i]);
			//step[0][i] = step[0][i - 1] + 1;
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
						//step[i][j] = step[i - 1][j] + 1;
					}
					else
					{
						distance[i][j] = g3;
						//step[i][j] = step[i][j - 1] + 1;
					}
				}
				else
				{
					if (g2 < g3)
					{
						distance[i][j] = g2;
						//step[i][j] = step[i - 1][j - 1];
					}
					else
					{
						distance[i][j] = g3;
						//step[i][j] = step[i][j - 1] + 1;
					}
				}
			}

		}
		printPath(dis, dis_cur, distance, first_num - 1, second_num - 1, path);
		// show the dp_path
		//for (int i = 0; i < dp_path.size(); i++)
		//{
		//cout << first_arr[dp_path[i].first] << " - " << second_arr[dp_path[i].second] << endl;
		//}
		return distance[first_num - 1][second_num - 1];
	}

	void Tiling_opt::printPath(double d[][202], double d_c[][202], double dp[][202], int i, int j, vector<pair<int, int>>& path)
	{

		if (i == 0 && j == 0) {
			//cout << first_arr[i] << " - " << second_arr[j] << endl;
			path.push_back(make_pair(i, j));
			return;
			//cout make_pair(i,j);
		}

		if (abs(dp[i][j] - (dp[i - 1][j - 1] + 2 * (d[i][j] + d_c[i][j]))) < 0.1){
			printPath(d, d_c, dp, i - 1, j - 1, path);

		}
		else if (abs(dp[i][j] - (dp[i][j - 1] + d[i][j] + d_c[i][j])) < 0.1) {
			printPath(d, d_c, dp, i, j - 1, path);

		}
		else {
			printPath(d, d_c, dp, i - 1, j, path);
		}
		path.push_back(make_pair(i, j));
		//cout << first_arr[i] << " - " << second_arr[j] << endl;

	}

	vector<int> Tiling_opt::search_align_p(Point2f cent, Point2f end, vector<Point2f> cand_temp)
	{
		//保证线段足够长来求交点，将线段长度放大3倍
		Point2f dis = end - cent;
		Point2f endf = Point2f(cent.x + 2 * dis.x, cent.y + 2 * dis.y);
		int ctsize = cand_temp.size();
		vector<int> cand_index;
		for (int i = 0; i < ctsize; i++)
		{
			Point2f crosP;
			if (line_intersection(Line_Seg(endf, cent), Line_Seg(cand_temp[i], cand_temp[(i + 1) % ctsize]), crosP) == 0) continue;
			else 
			{
				if (length_two_point2f(end, cand_temp[i]) < length_two_point2f(end, cand_temp[(i + 1) % ctsize]))
				{
					int flag = 0;
					for (int j = 0; j < cand_index.size(); j++)
					{
						if (cand_index[j] == i) flag = 1;
					}
					if (flag == 0) cand_index.push_back(i);
					else break;
				}
				else
				{
					int flag = 0;
					for (int j = 0; j < cand_index.size(); j++)
					{
						if (cand_index[j] == (i + 1) % ctsize) flag = 1;
					}
					if (flag == 0) cand_index.push_back((i + 1) % ctsize);
					else break;
				}
			}
			
		}
		//返回对齐点，可能只有一个，可能有多个,删掉重复点

		return cand_index;
		/*double leng = 10000;
		int min_in = 0;
		for (int t = 0; t < cand_index.size(); t++)
		{
		if (length_two_point2f(cand_temp[cand_index[t]],end) < leng)
		{
		leng = length_two_point2f(cand_temp[cand_index[t]], end);
		min_in = t;
		}

		}*/
		//return min_in;

	}

	vector<Point2f> Tiling_opt::CandP2Contour(CandPat candp, int num)
	{
		prototile_second->Pro_clear();
		prototile_second->loadPoints(contour_dataset[candp.number]);
		vector<Point2f> contour_cand;
		vector<Point2f> cand_tem;
		if (candp.isFilp) contour_cand = prototile_second->contour_sample_flip[num];
		else contour_cand = prototile_second->contour_sample[num];
		Point2f Ccen = center_p(contour_cand);
		Mat rot_mat;
		rot_mat = getRotationMatrix2D(Ccen, candp.angle, 1);
		transform(contour_cand, contour_cand, rot_mat);
		int sizec = contour_cand.size();
		for (int i = 0; i < sizec; i++)
		{
			cand_tem.push_back(contour_cand[(i + candp.index) % sizec]);
		}
		return cand_tem;
	}

	vector<Point2f> Tiling_opt::morphing_patterns_iter(vector<Point2f> contour1, vector<Point2f> contour2, vector<double> concur1, vector<double> concur2, float shape_ratio)//, vector<int> mid_inter, float shape_ratio)
	{
		vector<Point2f> final_pettern;
		cout << contour1.size() << "   " << contour2.size() << endl;
		double scale = arcLength(contour1, true) / arcLength(contour2, true);
		//cout << "scale: " << scale << endl;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] * scale;
		}
		Point2f cen1 = center_p(contour1);
		Point2f cen2 = center_p(contour2);
		Point2f shift2 = cen1 - cen2;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] + shift2;
		}

		int count = 100;
		while (count > 5)
		{
			count = 0;
			vector<pair<int, int>> dppath;
			double accumu_mis = quadr_mismatch(contour1, contour2, concur1, concur2, dppath, 10);
			vector<double> concur11 = recover_consin(concur1);
			vector<double> concur22 = recover_consin(concur2);
			int psize = dppath.size();
			int first = -1;
			int sec = -1;
			for (int i = 0; i < psize; i++)
			{
				if (first == dppath[i].first || sec == dppath[i].second) continue;
				else{
					first = dppath[i].first;
					sec = dppath[i].second;
					double ratios = (concur11[first] + 2) / ((concur11[first] + 2) + (concur22[sec] + 2));

					if (length_two_point2f(contour2[sec], contour1[first])> 0.2)
					{
						count++;
						Point2f vec = contour2[sec] - contour1[first];
						contour1[first] += shape_ratio * (1 - ratios)*vec;
						contour2[sec] += shape_ratio * ratios * (-vec);
					}
				}
				//Point2f fin = ratios * contour1[first] + (1.0 - ratios) * contour2[sec];
				//final_pettern.push_back(fin);
			}
			concur1 = curvature_com(contour1);
			concur2 = curvature_com(contour2);
		}

		Mat ta = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		for (int t = 0; t < contour1.size(); t++)
		{
			circle(tt, contour1[t], 3, Scalar(255, 0, 0), -1);
			circle(ttt, contour1[t], 3, Scalar(255, 0, 0), -1);
			circle(ta, contour1[t], 3, Scalar(255, 0, 0), -1);
		}
		for (int t = 0; t < contour2.size(); t++)
		{
			circle(tt, contour2[t], 3, Scalar(0, 255, 0), -1);
			circle(ttt, contour2[t], 3, Scalar(0, 255, 0), -1);
			circle(ta, contour2[t], 3, Scalar(0, 255, 0), -1);
		}
		cout << contour1.size() << "   " << contour2.size() << endl;
		//for (int t = 0; t < psize; t++)
		//{
		//	//circle(ttt, src1_points[t], 4, Scalar(255, 120, 120), -1);
		//	circle(ta, final_pettern[t], 3, Scalar(255, 120, 120), -1);
		//	//circle(ttt, src2_points[t], 4, Scalar(120, 200, 120), -1);
		//	//circle(tt, src2_points[t], 6, Scalar(120, 200, 120), -1);
		//	//MyLine(ttt, src1_points[t], src2_points[t], "grey");
		//	MyLine(tt, contour1[dppath[t].first], contour2[dppath[t].second], "grey");
		//}
		final_pettern = sampling(contour1,2);// sampling(final_pettern, 2);
		for (int i = 0; i < final_pettern.size(); i++)
		{
			circle(ttt, final_pettern[i], 3, Scalar(0, 0, 255), -1);
		}
		//for (int i = 0; i < qmpath.size(); i++)
		//{
		//	cout << qmpath[i].first << "   " << qmpath[i].second << endl;
		//	MyLine(ta, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
		//	//MyLine(tt, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
		//}
		imshow("??:", ta);
		imshow("???:", tt);
		imshow("????:", ttt);
		return final_pettern;
	}
	
	void Tiling_opt::points_dividing(string imaname) //
	{
		vector<int> p_p_index = prototile_first->partition_points(imaname);
		vector<Point2f> contour_ = prototile_first->contour;
		//vector<Point2f> contour_ = prototile_first->contour;
		load_dataset();
		string rootname = "D:\\VisualStudioProjects\\DihedralTessellation\\result\\" + prototile_first->contourname;
		const char *na = rootname.c_str();
		mkdir(na);
		int ppindex = p_p_index.size();
		int margin = prototile_first->contour.size() / 10;
		cout << "margin: " << margin << endl;
		int count = 0;

		vector<vector<Point2f>> inner_conts;
		vector<vector<int>> all_situation_index;
		vector<vector<int>> mid_interval_index;
		for (int i = 0; i < ppindex; i++)
		{
			for (int j = i + 1; j < ppindex; j++)
			{
				if (abs(p_p_index[j] - p_p_index[i]) < margin) continue;
				//cout << "i: " << p_p_index[i] << "   j: " << p_p_index[j% ppindex] << endl;
				for (int m = j + 1; m < ppindex; m++)
				{
					if (abs(p_p_index[m] - p_p_index[j]) < margin) continue;
					for (int n = m + 1; n < ppindex; n++)
					{
						if (abs(p_p_index[n] - p_p_index[m]) < margin) continue;
						vector<Point2f> inner_contour;
						vector<int> result_index;
						vector<int> mid_interval; //mid_interval存储的是组合成的inner 的分段连接点
						result_index.push_back(p_p_index[i]);
						result_index.push_back(p_p_index[j]);
						result_index.push_back(p_p_index[m]);
						result_index.push_back(p_p_index[n]);
						//cout << "  i: " << p_p_index[i]
						//	<< "   j: " << p_p_index[j]
						//	<< "   m: " << p_p_index[m]
						//	<< "   n: " << p_p_index[n] << endl;	
						////one_situ_div(result_index, contour_, inner_contour);					
						//string image = int2string(count);
						/*string conut_name = rootname + "\\placement " + int2string(count);
						na = conut_name.c_str();
						mkdir(na);*/
						Mat drawing1 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255));
						if (!one_situ_div(result_index, contour_, inner_contour, mid_interval, drawing1))
						{
							//cout << "-------------collision-------------" << endl;
							continue;
						}
						count++;
						//show the marked points
						Point2f shift2 = Point2f(400, 400) - center_p(contour_);
						for (int j = 0; j < contour_.size(); j++)
						{
							circle(drawing1, contour_[j] + shift2, 1, Scalar(0, 0, 0), -1);

							//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
						}
						for (int j = 0; j < p_p_index.size(); j++)
						{
							circle(drawing1, contour_[p_p_index[j]] + shift2, 4, Scalar(0, 0, 255), -1);
						}
						circle(drawing1, contour_[result_index[0]] + shift2, 8, Scalar(0, 255, 0), -1);
						circle(drawing1, contour_[result_index[1]] + shift2, 8, Scalar(0, 255, 0), -1);
						circle(drawing1, contour_[result_index[2]] + shift2, 8, Scalar(0, 255, 0), -1);
						circle(drawing1, contour_[result_index[3]] + shift2, 8, Scalar(0, 255, 0), -1);
						//imshow("result_mid", drawing1);
						string filename = rootname + "\\" + int2string(count - 1) + "PlacingResult.png";
						imwrite(filename, drawing1);


						cout << endl << count << " succeed";
						cout << "    inner_contour.size : " << inner_contour.size() << endl;

						inner_conts.push_back(inner_contour);
						all_situation_index.push_back(result_index);
						mid_interval_index.push_back(mid_interval);
						//if (count == 2) break;

					}
					//if (count == 2) break;
				}
				//if (count == 2) break;
			}
			//if (count == 2) break;
		}
		//string conut_name = rootname + "\\placement " + int2string(count);
		//rmdir(conut_name.c_str());
		cout << "succeed count: " << count << endl;
		// search the right image
		/*
		1.将周长调整为一致
		2.搜索最相近的图案以及角度
		3.变形
		*/
		if (count == 0)
		{
			cout << "no right placement" << endl;
			return;
		}

		for (int i = 0; i < inner_conts.size(); i++) //inner_conts.size()
		{
			cout << "count: " << i << "/" << count - 1 << endl;
			int num_c = 1;//选择(num_c+1)*100个点
			vector<CandPat> candida_contours;

			/*cout << "inner_one: " << inner_one  << endl;
			string filepathname = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test" + int2string(i)+".txt";
			vector<Point> con;
			for (int i = 0; i < inner_conts[inner_one].size(); i++)
			con.push_back((Point)inner_conts[inner_one][i]);
			fileout(filepathname, con);*/


			//compare_choose_TAR(prototile_first->contour)
			candida_contours = compare_shapes(inner_conts[i], num_c);
			//cout << "candida_contours" << candida_contours.size()<< endl;
			//string conut_name = rootname + "\\placement " + int2string(i);
			vector<int> mid_inter = joint_relocate(inner_conts[i], mid_interval_index[i], num_c);
			for (int j = 0; j <candida_contours.size(); j++) //candida_contours.size()
			{
				//将所有的结果保存下来
				Mat drawing_pro = Mat(800, 3200, CV_8UC3, Scalar(255, 255, 255));
				CandPat tem = candida_contours[j];
				prototile_second->Pro_clear();
				prototile_second->loadPoints(contour_dataset[tem.number]);
				prototile_mid->Pro_clear();
				prototile_mid->loadPoints(inner_conts[i]);

				vector<Point2f> contour_inner = prototile_mid->contour_sample[num_c]; //选择最少的点进行比较
				vector<double> contour_inner_c = curvature_com(contour_inner);// prototile_mid->contour_curva[0];
				vector<Point2f> contour_cand = CandP2Contour(tem, num_c);
				vector<double> contour_cand_c = curvature_com(contour_cand);
				//inner and cand morph into the final pettern	
				//!!!!!!在下一步中将morph改为分段!!!此时没有考虑中心平移和缩放倍数，因为之前是一一对应进行的变形
				int num = 1;
				float ratio = 4;
				while (num-- != 0)
				{
					vector<int> mid_inter_new = mid_inter;
					ratio = ratio + 1;
					vector<Point2f> inter_mid = morphing_2_patterns(contour_inner, contour_cand, contour_inner_c, contour_cand_c, mid_inter_new, ratio / 10);

					//show the result
					draw_poly(drawing_pro, contour_inner, Point2f(400, 400));
					draw_poly(drawing_pro, contour_cand, Point2f(1200, 400));
					draw_poly(drawing_pro, inter_mid, Point2f(2000, 400));
					//draw_polygen("cand pattern: ", contour_cand);

					vector<int> return_p;
					vector<vector<Point2f>> four_;
					vector<Point2f> morphed_B = extract_contour(inter_mid, mid_inter_new, return_p, four_, 0);
					//evalua_deformation();
					Mat drawing_pro2 = Mat(1600, 2400, CV_8UC3, Scalar(255, 255, 255));
					draw_allplane(drawing_pro2, morphed_B, return_p, 0.4);

					string filename = rootname + "\\";
					string file2 = filename + int2string(i) + "_Candidate_" + int2string(j) + "tiling_result.png";
					filename = filename + int2string(i) + "_Candidate_" + int2string(j) + ".png";
					imwrite(filename, drawing_pro);
					imwrite(file2, drawing_pro2);
					imshow("tiling_result: ", drawing_pro2);
				}

			}
		}
	}

	vector<Point2f> Tiling_opt::morphing_2_patterns(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<double> &concur1, vector<double> &concur2, vector<int> &mid_inter, float shape_ratio)
	{
		//传统morphing是由start和end求得一系列intermediate状态，这里是通过c1和c2获得最终的变形结果

		vector<Point2f> final_pettern;
		vector<int> mid_inter_new;
		/*int difference = contour1.size() - contour2.size();
		difference = abs(difference);
		if (difference != 0)
		{
			if (contour1.size() < contour2.size())
			{
				for (int i = 0; i < difference; i++)
				{
					contour2.pop_back();
				}
			}
			else if (contour1.size() > contour2.size())
			{
				for (int i = 0; i < difference; i++)
				{
					contour1.pop_back();
				}
			}
		}*/
		cout << contour1.size() << "   " << contour2.size() << endl;
		double scale = arcLength(contour1, true) / arcLength(contour2, true);
		//cout << "scale: " << scale << endl;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] * scale;
		}
		Point2f cen1 = center_p(contour1);
		Point2f cen2 = center_p(contour2);
		Point2f shift2 = cen1 - cen2;
		for (int i = 0; i < contour2.size(); i++)
		{
			contour2[i] = contour2[i] + shift2;
		}

		vector<pair<int, int>> dppath;

		double accumu_mis = quadr_mismatch(contour1, contour2, concur1, concur2, dppath,10);
		vector<double> concur11 = recover_consin(concur1);
		vector<double> concur22 = recover_consin(concur2);
		/*for (int i = 0; i < mid_inter.size(); i++)
		{
			concur11[mid_inter[i]] += 10;
		}*/
		int psize = dppath.size();
		int t = 0;
		for (int i = 0; i < psize; i++)
		{
			int first = dppath[i].first;
			int sec = dppath[i].second; 
			//cout << "first: " << first << " sec: " << sec  << endl;
			double ratios = (concur11[first] + 2) / ((concur11[first] + 2) + (concur22[sec] + 2));
			Point2f fin = ratios * contour1[first] + (1.0 - ratios) * contour2[sec];
			final_pettern.push_back(fin);
			if (t<4 && first == mid_inter[t])
			{
				t++;
				mid_inter_new.push_back(i);
			}
		}
		//cout << "mid_inter_new: " << mid_inter_new.size() << endl;
		Mat ta = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat tt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Mat ttt = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		for (int t = 0; t < contour1.size(); t++)
		{
			circle(tt, contour1[t], 3, Scalar(255, 0, 0), -1);
			circle(ttt, contour1[t], 3, Scalar(255, 0, 0), -1);
			circle(ta, contour1[t], 3, Scalar(255, 0, 0), -1);
		}
		for (int t = 0; t < contour2.size(); t++)
		{
			circle(tt, contour2[t], 3, Scalar(0, 255, 0), -1);
			circle(ttt, contour2[t], 3, Scalar(0, 255, 0), -1);
			circle(ta, contour2[t], 3, Scalar(0, 255, 0), -1);
		}
		cout << contour1.size() << "   " << contour2.size() << endl;
		for (int t = 0; t < psize; t++)
		{
			//circle(ttt, src1_points[t], 4, Scalar(255, 120, 120), -1);
			circle(ta, final_pettern[t], 3, Scalar(255, 120, 120), -1);
			//circle(ttt, src2_points[t], 4, Scalar(120, 200, 120), -1);
			//circle(tt, src2_points[t], 6, Scalar(120, 200, 120), -1);
			//MyLine(ttt, src1_points[t], src2_points[t], "grey");
			MyLine(tt, contour1[dppath[t].first], contour2[dppath[t].second], "grey");
		}
		mid_inter_new = joint_relocate(final_pettern, mid_inter_new,1);
		mid_inter = mid_inter_new;
		final_pettern = sampling(final_pettern,2);
		for (int i = 0; i < final_pettern.size(); i++)
		{
			circle(ttt, final_pettern[i], 3, Scalar(0, 0, 255), -1);
		}
		//for (int i = 0; i < qmpath.size(); i++)
		//{
		//	cout << qmpath[i].first << "   " << qmpath[i].second << endl;
		//	MyLine(ta, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
		//	//MyLine(tt, cont1[qmpath[i].first], cont2[qmpath[i].second], "grey");
		//}
		imshow("??:", ta);
		imshow("???:",tt);
		imshow("????:", ttt);
		////vector<int> mid_in = joint_relocate(contour1, mid_inter, 1);
		//vector<Point2f> cont1 = sampling(contour1, 1);
		//vector<Point2f> cont2 = sampling(contour2, 1);
		//vector<double> cont1_c = curvature_com(cont1);
		//vector<double> cont2_c = curvature_com(cont2);

		//vector<int> cand_points_index = most_convex_p(cont1, cont1_c,20);
		//sort_bub(cand_points_index);
		//vector<pair<int, int>> qmpath;
		//quadr_mismatch(cont1, cont2, cont1_c, cont2_c, qmpath);


		//// 首先找到一一对应的点序列
		//vector<Point2f> src1_points;
		//vector<Point2f> src2_points;
		//src1_points.push_back(cont1[0]);
		//src2_points.push_back(cont2[0]);

		////for (int i = 0; i < dp_path.size(); i++)
		////{
		////	cout << dp_path[i].first << " -- " << dp_path[i].second << endl;
		////}
		//if (cont1.size() < cont2.size())
		//{
		//	int f = 1;
		//	int j = 0;
		//	for (int i = 1; i < cont1.size() - 1;)
		//	{

		//		double min = 10000;
		//		while (f < qmpath.size())
		//		{

		//			if ((qmpath[f].first < i) || qmpath[f].second < j) //检查下一条路径的另一端点是否已被使用
		//			{
		//				f++;
		//				if (qmpath[f].first > i) i++;
		//				continue;
		//			}
		//			if (qmpath[f].first == i)
		//			{
		//				if (length_two_point2f(cont1[qmpath[f].first], cont2[qmpath[f].second]) < min)
		//				{
		//					min = length_two_point2f(cont1[qmpath[f].first], cont2[qmpath[f].second]);
		//					j = qmpath[f].second;
		//				}
		//				f++;
		//			}
		//			if (qmpath[f].first > i)
		//			{
		//				src1_points.push_back(cont1[i]);
		//				src2_points.push_back(cont2[j]);
		//				i++;
		//				j++;
		//				break;
		//			}

		//		}
		//	}
		//}
		//else
		//{
		//	int f = 1;
		//	int j = 0;
		//	for (int i = 1; i < cont2.size() - 1;)
		//	{
		//		double min = 10000;

		//		while (f < qmpath.size())
		//		{

		//			if ((qmpath[f].second < i) || qmpath[f].first < j)
		//			{

		//				f++;
		//				if (qmpath[f].second > i) i++;
		//				continue;
		//			}
		//			if (qmpath[f].second == i)
		//			{
		//				if (length_two_point2f(cont1[qmpath[f].first], cont2[qmpath[f].second]) < min)
		//				{
		//					min = length_two_point2f(cont1[qmpath[f].first], cont2[qmpath[f].second]);
		//					j = qmpath[f].first;
		//				}
		//				f++;
		//			}
		//			if (qmpath[f].second > i)
		//			{
		//				src1_points.push_back(cont1[j]);
		//				src2_points.push_back(cont2[i]);
		//				i++;
		//				j++;
		//				break;
		//			}

		//		}
		//	}
		//}
		//src1_points.push_back(cont1[qmpath[qmpath.size() - 1].first]);
		//src2_points.push_back(cont2[qmpath[qmpath.size() - 1].second]);

		

		//MorphPoints(contour1, contour2, final_pettern, shape_ratio);


		//MorphPoints(src1_points, src2_points, final_pettern, shape_ratio);
		//cout << "contour1: "<<contour1.size() << "  contour2: " << contour2.size() << "  final_pettern: " << final_pettern.size() << endl;

		return final_pettern;
	}



	bool Tiling_opt::one_situ_div(vector<int> results, vector<Point2f> contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname) //检测一种划分情况的结果
	{
		Point2f line1 = contour_s[results[2]] - contour_s[results[0]];
		Point2f line2 = contour_s[results[3]] - contour_s[results[1]];



		vector<Point2f> shift_4;
		shift_4.push_back(Point2f(0, 0));
		shift_4.push_back(line1);
		shift_4.push_back(line2);
		shift_4.push_back(line1 + line2);

		//先在这里进行接触点角度的统计，大于360视为会发生碰撞


		vector<Point2f> bbox_s = b_box(contour_s);
		return_B.swap(vector<Point2f>());

		int fpsize = shift_4.size();
		//cout <<"four_place.size(): "<< endl;

		for (int i = 0; i < fpsize; i++)
		{
			for (int j = i + 1; j < fpsize; j++)
			{
				if (coll_detection(shift_4[i], shift_4[j], contour_s))
				{
					return false;
				}
				//else {
				//	cout << "i: " << i << "j: " << j << endl;
				//}
			}

		}
		vector<vector<Point2f>> four_place;
		vector<Point2f> one_loca;
		// 提取围成的轮廓，目前为止只考虑正向摆放，不考虑旋转和翻转
		four_place.push_back(contour_s);

		for (int i = 0; i < contour_s.size(); i++)
		{
			one_loca.push_back(contour_s[i] + line1);
		}
		four_place.push_back(one_loca);
		one_loca.swap(vector<Point2f>());

		for (int i = 0; i < contour_s.size(); i++)
		{
			one_loca.push_back(contour_s[i] + line2);
		}
		four_place.push_back(one_loca);
		one_loca.swap(vector<Point2f>());

		for (int i = 0; i < contour_s.size(); i++)
		{
			one_loca.push_back(contour_s[i] + line1 + line2);
		}
		four_place.push_back(one_loca);

		//将该proto1以及相邻四个proto2展示出来
		//Mat drawing_pro = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));

		//Point2f shift1 = Point2f(1200 - 0.2 * (2*center_p(contour_s).x + line1.x + line2.x), 400 - 0.2 * (2*center_p(contour_s).y + line1.y + line2.y));
		Point2f shift1 = Point2f(1200, 400) - (0.4*center_p(contour_s) + 0.2*line1 + 0.2*line2);
		// show bbox
		//for (int i = 0; i < 4; i++)
		//{
		//	for (int j = 0; j < 4; j++)
		//	{
		//		MyLine(drawing_pro, f_f_cor[i][j] * 0.4 + shift1, f_f_cor[i][(j + 1) % 4] * 0.4 + shift1, "red");
		//	}
		//}

		//show four proto
		//vector<vector<Point2f>> four_place_;

		for (int i = 0; i < 4; i++)
		{
			vector<Point2f> one_;
			//MyLine(drawing_pro, four_cor[i]*0.4+shift1, four_cor[(i+1)%4]*0.4+shift1, "red");
			//prototwoAff_place.swap(vector<Point2f>());
			for (int j = 0; j < four_place[i].size(); j++)
			{
				one_.push_back(four_place[i][j] * 0.4 + shift1);
				//MyLine(countname, four_place[i][j] * 0.4 + shift1, four_place[i][j + 1] * 0.4 + shift1, "green");
			}
			//four_place_.push_back(one_);
			draw_poly(countname, one_, center_p(one_));
		}
		int total_num = 0;
		return_p.push_back(0);
		for (int t = results[3] - 1; t > results[2] + 1; t--)
		{
			total_num++;
			return_B.push_back(four_place[0][t]);
		}
		return_p.push_back(total_num);
		for (int t = results[0] - 1; t > 0; t--)
		{
			total_num++;
			return_B.push_back(four_place[1][t]);
		}
		for (int t = four_place[1].size() - 1; t > results[3] + 1; t--)
		{
			total_num++;
			return_B.push_back(four_place[1][t]);
		}
		return_p.push_back(total_num);
		for (int t = results[1] - 1; t > results[0] + 1; t--)
		{
			total_num++;
			return_B.push_back(four_place[3][t]);
		}
		return_p.push_back(total_num);
		for (int t = results[2] - 1; t > results[1] + 1; t--)
		{
			return_B.push_back(four_place[2][t]);
		}
		cout << "num: " << return_B.size();
		Point2f cent_p = center_p(return_B);
		Point2f shift3 = Point2f(400, 400) - cent_p;
		for (int z = 0; z < return_B.size(); z++)
		{
			return_B[z] = return_B[z] + shift3;
		}
		//show middle pattern
		//Mat drawing7 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));				
		//for (int z = 0; z < return_B.size(); z++)
		//{
		//	circle(drawing7, return_B[z], 1, Scalar(0, 0, 0), -1);
		//	//MyLine(drawing4, prototile_first->contour_sample[sam_num][j] - shift1, prototile_first->contour_sample[sam_num][j + 1] - shift1, "red");
		//}
		//imshow("b:", drawing7);
		//draw_polygen("B:re_show", return_B);
		return true;
	}

	bool Tiling_opt::coll_detection(Point2f shifting1, Point2f shifting2, vector<Point2f> &contour_s)
	{
		//首先通过包围盒求得粗糙的相交区域，然后通过像素相交求是否产生碰撞
		vector<Point2f> bbox_s = b_box(contour_s);
		vector<Point2f> bbox1;
		vector<Point2f> bbox2;
		for (int i = 0; i < bbox_s.size(); i++)
		{
			bbox1.push_back(bbox_s[i] + shifting1);
			bbox2.push_back(bbox_s[i] + shifting2);
		}
		vector<int> num_in1;
		vector<int> num_in2;
		for (int i = 0; i < 4; i++)
		{
			if (bbox1[i].x <= bbox2[3].x && bbox1[i].x >= bbox2[1].x)
				if (bbox1[i].y <= bbox2[3].y && bbox1[i].y >= bbox2[1].y)
				{
					num_in1.push_back(i);
				}
		}
		for (int i = 0; i < 4; i++)
		{
			if (bbox2[i].x <= bbox1[3].x && bbox2[i].x >= bbox1[1].x)
				if (bbox2[i].y <= bbox1[3].y && bbox2[i].y >= bbox1[1].y)
				{
					num_in2.push_back(i);
				}
		}
		int dis = num_in1.size() - num_in2.size();
		//cout << "num_in1.size() and num_in2.size()" << num_in1.size() << "   " << num_in2.size() << endl;
		if (num_in1.size() == 0 && num_in2.size() == 0)
		{
			return false;
		}
		else if (num_in1.size() == 1 && num_in2.size() == 1)
		{
			double bbx_max_x = max(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x);
			double bbx_max_y = max(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y);
			double bbx_min_x = min(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x);
			double bbx_min_y = min(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y);
			Point2f max_p = Point2f(bbx_max_x + 0.5, bbx_max_y + 0.5);
			Point2f min_p = Point2f(bbx_min_x - 0.5, bbx_min_y - 0.5);
			//cout << max_p << endl << min_p << endl;
			vector<Point2f> dis_p;
			dis_p.push_back(min_p - bbox1[0]);
			dis_p.push_back(max_p - bbox1[0]);
			dis_p.push_back(min_p - bbox2[0]);
			dis_p.push_back(max_p - bbox2[0]);
			/*dis_p.push_back(Point2f(min_p.x - (int)bbox1[0].x, min_p.y - (int)bbox1[0].y));
			dis_p.push_back(Point2f(max_p.x - (int)bbox1[0].x, max_p.y - (int)bbox1[0].y));
			dis_p.push_back(Point2f(min_p.x - (int)bbox2[0].x, min_p.y - (int)bbox2[0].y));
			dis_p.push_back(Point2f(max_p.x - (int)bbox2[0].x, max_p.y - (int)bbox2[0].y));*/
			return collision_pixel(dis_p, contour_s);
		}
		else if (num_in1.size() == 2 && num_in2.size() == 2)
		{
			double bbx_max_x = max(max(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x), bbox1[num_in1[1]].x);
			double bbx_max_y = max(max(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y), bbox1[num_in1[1]].y);
			double bbx_min_x = min(min(bbox1[num_in1[0]].x, bbox2[num_in2[0]].x), bbox1[num_in1[1]].x);
			double bbx_min_y = min(min(bbox1[num_in1[0]].y, bbox2[num_in2[0]].y), bbox1[num_in1[1]].x);
			Point2f max_p = Point2f(bbx_max_x + 0.5, bbx_max_y + 0.5);
			Point2f min_p = Point2f(bbx_min_x - 0.5, bbx_min_y - 0.5);
			vector<Point2f> dis_p;
			dis_p.push_back(min_p - bbox1[0]);
			dis_p.push_back(max_p - bbox1[0]);
			dis_p.push_back(min_p - bbox2[0]);
			dis_p.push_back(max_p - bbox2[0]);
			return collision_pixel(dis_p, contour_s);
		}
		else if (abs(dis) == 2)
		{
			return true;

		}
		else if (abs(dis) == 4)
		{
			double bbx_max_x;
			double bbx_max_y;
			double bbx_min_x;
			double bbx_min_y;
			if (num_in1.size() == 4)
			{
				bbx_max_x = max(max(bbox1[num_in1[0]].x, bbox1[num_in1[1]].x), bbox1[num_in1[2]].x);
				bbx_max_y = max(max(bbox1[num_in1[0]].y, bbox1[num_in1[1]].y), bbox1[num_in1[2]].y);
				bbx_min_x = min(min(bbox1[num_in1[0]].x, bbox1[num_in1[1]].x), bbox1[num_in1[2]].x);
				bbx_min_y = min(min(bbox1[num_in1[0]].y, bbox1[num_in1[1]].y), bbox1[num_in1[2]].y);
			}
			else if (num_in2.size() == 4)
			{
				bbx_max_x = max(max(bbox2[num_in2[0]].x, bbox2[num_in2[1]].x), bbox2[num_in2[2]].x);
				bbx_max_y = max(max(bbox2[num_in2[0]].y, bbox2[num_in2[1]].y), bbox2[num_in2[2]].y);
				bbx_min_x = min(min(bbox2[num_in2[0]].x, bbox2[num_in2[1]].x), bbox2[num_in2[2]].x);
				bbx_min_y = min(min(bbox2[num_in2[0]].y, bbox2[num_in2[1]].y), bbox2[num_in2[2]].y);
			}
			Point2f max_p = Point2f(bbx_max_x + 0.5, bbx_max_y + 0.5);
			Point2f min_p = Point2f(bbx_min_x - 0.5, bbx_min_y - 0.5);
			vector<Point2f> dis_p;
			dis_p.push_back(min_p - bbox1[0]);
			dis_p.push_back(max_p - bbox1[0]);
			dis_p.push_back(min_p - bbox2[0]);
			dis_p.push_back(max_p - bbox2[0]);
			return collision_pixel(dis_p, contour_s);
		}
		return true;
	}

	bool Tiling_opt::collision_pixel(vector<Point2f> dis_p, vector<Point2f> contour_s)
	{
		int coll_num = 10;

		vector<Point2f> bbox_s = b_box(contour_s);
		Point2f min1 = bbox_s[0] + dis_p[0];
		Point2f max1 = bbox_s[0] + dis_p[1];
		Mat drawing_mid = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		int n = contour_s.size();
		Point rook_points[1][800];
		for (int t = 0; t < n; t++)
		{
			rook_points[0][t] = contour_s[t];
		}
		const Point* ppt[1] = { rook_points[0] };
		int npt[] = { n };
		fillPoly(drawing_mid,
			ppt,
			npt,
			1,
			Scalar(0, 0, 0) //黑色
			);
		imshow("coli 1", drawing_mid);
		//cout << "zahuishia " << endl;
		if (min1.y<0 || max1.y<0 || min1.x<0 || max1.x<0) return true;
		//cout << "min1.y " << min1.y << "  max1.y " << max1.y << " min1.x " << min1.x << " max1.x" << max1.x<<endl;
		Mat draw1 = drawing_mid(Range(min1.y, max1.y), Range(min1.x, max1.x));
		cvtColor(draw1, draw1, COLOR_BGR2GRAY);
		threshold(draw1, draw1, 128, 1, cv::THRESH_BINARY);

		Point2f min2 = bbox_s[0] + dis_p[2];
		Point2f max2 = bbox_s[0] + dis_p[3];
		if (min2.y<0 || max2.y<0 || min2.x<0 || max2.x<0) return true;
		//cout << "min2.y " << min2.y << "  max2.y " << max2.y << " min2.x " << min2.x << " max2.x" << max2.x<<endl;
		Mat draw2 = drawing_mid(Range(min2.y, max2.y), Range(min2.x, max2.x));
		cvtColor(draw2, draw2, COLOR_BGR2GRAY);
		threshold(draw2, draw2, 128, 1, cv::THRESH_BINARY);

		int count = 0;
		int rows = draw1.rows < draw2.rows ? draw1.rows : draw2.rows;
		int cols = draw1.cols < draw2.cols ? draw1.cols : draw2.cols;
		//cout << "rows: " << rows << "cols: " << cols << endl;

		//rows = draw2.rows;
		//cols = draw2.cols;
		//cout << "rows: " << rows << "cols: " << cols << endl;
		for (int i = 0; i <rows; i++)
			for (int j = 0; j < cols; j++)
			{
				draw1.at<uchar>(i, j) = (int)draw1.at<uchar>(i, j) + (int)draw2.at<uchar>(i, j);
				if ((int)draw1.at<uchar>(i, j) == 0) count++;
			}
		//cout << "col num:  " << count << endl;

		if (count > coll_num) return true;
		else
		{
			threshold(draw1, draw1, 0.5, 255, cv::THRESH_BINARY);
			imshow("coli result", draw1);
			return false;
		}


		/*Mat drawing4 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
		Point2f shift1 = Point2f(400, 400) - contour1[0]*0.4;
		for (int i = 0; i < contour2.size(); i++)
		{
		circle(drawing4, contour1[i]*0.4+shift1, 1, Scalar(0, 0, 0), -1);
		circle(drawing4, contour2[i]*0.4+shift1, 1, Scalar(0, 0, 0), -1);
		}
		MyLine(drawing4, Point2f(min_p.x, max_p.y)*0.4 + shift1, Point2f(min_p.x, min_p.y)*0.4 + shift1, "red");
		MyLine(drawing4, Point2f(min_p.x, min_p.y)*0.4 + shift1, Point2f(max_p.x, min_p.y)*0.4 + shift1, "red");
		MyLine(drawing4, Point2f(max_p.x, min_p.y)*0.4 + shift1, Point2f(max_p.x, max_p.y)*0.4 + shift1, "red");
		MyLine(drawing4, Point2f(max_p.x, max_p.y)*0.4 + shift1, Point2f(min_p.x, max_p.y)*0.4 + shift1, "red");
		imshow("cand_fram",drawing4);*/

		//

	}

}
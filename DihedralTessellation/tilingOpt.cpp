/* 1. �������������������
2. ��ð�͹�Խ��л���
3. ��ͼ��������װ����dataset��Ѱ�����϶���ƶȸߵ�ͼ��
4. ���κ�decorate
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
		all_types = 695;
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
			prototile_second->getpath();
			string image = int2string(i);
			vector<Point2f> data_;
			prototile_second->contourname = image;

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
		for (int i = 0; i <= all_types; i++)
		{
			prototile_second->Pro_clear();
			prototile_second->loadPoints(contour_dataset[i]);
			double shape_com;
			vector<vector<double>> tar_all = prototile_second->compute_TAR(prototile_second->contour_sample[num_c], shape_com);//(num_c+1)*100 points
			all_con_tars.push_back(tar_all);
			all_shape_complexity.push_back(shape_com);
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
						vector<int> mid_interval; //mid_interval�洢������ϳɵ�inner �ķֶ����ӵ�
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
						Mat drawing1 = Mat(800, 1600, CV_8UC3, CV_8UC3);
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
						string filename = rootname + "\\" + int2string(count-1)+ "PlacingResult.png";
						imwrite(filename, drawing1);

						
						cout << endl << count<< " succeed" ;
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
		1.���ܳ�����Ϊһ��
		2.�����������ͼ���Լ��Ƕ�
		3.����
		*/
		if (count == 0)
		{
			cout << "no right placement" << endl;
			return;
		}

		for (int i = 0; i < inner_conts.size(); i++) //inner_conts.size()
		{
			cout << "count: " << i<<"/"<<count-1 << endl;
			int num_c = 1;//ѡ��(num_c+1)*100����
			vector<CandPat> candida_contours;

			//ofstream out("D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test.txt");
			//if (out.is_open())
			//{
			//	out << inner_conts[i].size() + 1 << endl;//contours[0].size()
			//	for (int j = 0; j < inner_conts[i].size(); j++)
			//		out << (int)inner_conts[i][j].x << "," << (int)inner_conts[i][j].y << endl;
			//	out << (int)inner_conts[i][0].x << "," << (int)inner_conts[i][0].y << endl;  //��β������
			//}
			//cout << "contours[0].size(): " << i << endl;
			//out.close();

			candida_contours = compare_shapes(inner_conts[i], num_c);
			//cout << "candida_contours" << candida_contours.size()<< endl;
			//string conut_name = rootname + "\\placement " + int2string(i);
			vector<int> mid_inter = joint_relocate(inner_conts[i], mid_interval_index[i], num_c);
			for (int j = 0; j <candida_contours.size(); j++) //candida_contours.size()
			{
				//�����еĽ����������
				Mat drawing_pro = Mat(800, 3200, CV_8UC3, Scalar(255, 255, 255));
				CandPat tem = candida_contours[j];
				prototile_second->Pro_clear();
				prototile_second->loadPoints(contour_dataset[tem.number]);
				prototile_mid->Pro_clear();
				prototile_mid->loadPoints(inner_conts[i]);

				vector<Point2f> contour_inner = prototile_mid->contour_sample[num_c]; //ѡ�����ٵĵ���бȽ�
				vector<double> contour_inner_c = curvature_com(contour_inner);// prototile_mid->contour_curva[0];
				vector<Point2f> contour_cand = CandP2Contour(tem, num_c);
				vector<double> contour_cand_c = curvature_com(contour_cand);
				//inner and cand morph into the final pettern	
				//!!!!!!����һ���н�morph��Ϊ�ֶ�!!!��ʱû�п�������ƽ�ƺ����ű�������Ϊ֮ǰ��һһ��Ӧ���еı���
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
					
					//����proto1�Լ������ĸ�proto2չʾ����
					Point2f line1 = inter_mid[mid_inter_new[2]] - inter_mid[mid_inter_new[0]];
					Point2f line2 = inter_mid[mid_inter_new[3]] - inter_mid[mid_inter_new[1]];
					Point2f cente = center_p(inter_mid);
					vector<vector<Point2f>> four_place;
					vector<Point2f> one_loca;
					// ��ȡΧ�ɵ�������ĿǰΪֹֻ��������ڷţ���������ת�ͷ�ת
					four_place.push_back(inter_mid);
					for (int i = 0; i < inter_mid.size(); i++)
					{
						one_loca.push_back(inter_mid[i] + line1);
					}
					four_place.push_back(one_loca);
					one_loca.swap(vector<Point2f>());
					for (int i = 0; i < inter_mid.size(); i++)
					{
						one_loca.push_back(inter_mid[i] + line2);
					}
					four_place.push_back(one_loca);
					one_loca.swap(vector<Point2f>());
					for (int i = 0; i < inter_mid.size(); i++)
					{
						one_loca.push_back(inter_mid[i] + line1 + line2);
					}
					four_place.push_back(one_loca);
					int total_num = 0;
					vector<int> return_p;// 
					return_p.push_back(0);
					vector<Point2f> morphed_B;
					for (int t = mid_inter_new[3] - 1; t > mid_inter_new[2] + 1; t--)
					{
						total_num++;
						morphed_B.push_back(four_place[0][t]);
					}
					return_p.push_back(total_num);
					for (int t = mid_inter_new[0] - 1; t > 0; t--)
					{
						total_num++;
						morphed_B.push_back(four_place[1][t]);
					}
					for (int t = four_place[1].size() - 1; t > mid_inter_new[3] + 1; t--)
					{
						total_num++;
						morphed_B.push_back(four_place[1][t]);
					}
					return_p.push_back(total_num);
					for (int t = mid_inter_new[1] - 1; t > mid_inter_new[0] + 1; t--)
					{
						total_num++;
						morphed_B.push_back(four_place[3][t]);
					}
					return_p.push_back(total_num);
					for (int t = mid_inter_new[2] - 1; t > mid_inter_new[1] + 1; t--)
					{
						morphed_B.push_back(four_place[2][t]);
					}
					Mat drawing_pro2 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
					Point2f line3 = morphed_B[return_p[2]] - morphed_B[return_p[0]];
					Point2f line4 = morphed_B[return_p[3]] - morphed_B[return_p[1]];
					vector<vector<Point2f>> four_place1;
					vector<Point2f> one_loca1;
					four_place1.push_back(morphed_B);
					//circle(drawing_pro2, morphed_B[return_p[2]], 4, Scalar(255, 0, 255), -1);
					//circle(drawing_pro2, morphed_B[return_p[0]], 4, Scalar(0, 255, 0), -1);
					//draw_poly(drawing_pro2, morphed_B, center_p(morphed_B));
					for (int i = 0; i < morphed_B.size(); i++)
					{
						one_loca1.push_back(morphed_B[i] + line3);
					}
					//draw_poly(drawing_pro2, one_loca1, center_p(one_loca1));
					four_place1.push_back(one_loca1);
					one_loca1.swap(vector<Point2f>());
					for (int i = 0; i < morphed_B.size(); i++)
					{
						one_loca1.push_back(morphed_B[i] + line4);
					}
					four_place1.push_back(one_loca1);
					one_loca1.swap(vector<Point2f>());
					for (int i = 0; i < morphed_B.size(); i++)
					{
						one_loca1.push_back(morphed_B[i] + line3 + line4);
					}
					four_place1.push_back(one_loca1);


					Mat drawing_pro1 = Mat(800, 800, CV_8UC3, Scalar(255, 255, 255));
					Point2f shift1 = Point2f(300 - cente.x*0.4, 400 - cente.y*0.4);
					Point2f shift3 = Point2f(2800, 400) - (0.4*center_p(morphed_B) + 0.2*line3 + 0.2*line4);


					//show four proto
					for (int i = 0; i < 4; i++)
					{

						vector<Point2f> ttt;
						for (int t = 0; t <four_place1[i].size(); t++)
						{
							ttt.push_back(four_place1[i][t] * 0.4 + shift3);
						}
						draw_poly(drawing_pro, ttt, center_p(ttt));
						//MyLine(drawing_pro, four_cor[i]*0.4+shift1, four_cor[(i+1)%4]*0.4+shift1, "red");
						//prototwoAff_place.swap(vector<Point2f>());
						for (int j = 0; j < four_place[i].size() - 1; j++)
						{
							MyLine(drawing_pro1, four_place[i][j] * 0.4 + shift1, four_place[i][j + 1] * 0.4 + shift1, "green");
						}
					}
					string filename = rootname + "\\";
					filename = filename + int2string(i) + "_Candidate_" + int2string(j) + ".png";
					imwrite(filename, drawing_pro);

					imshow("result_mid_show: ", drawing_pro1);
					imshow("result_mid: ", drawing_pro2);
				}

			}
		}


	}


	vector<Point2f> Tiling_opt::simulation_mid(string imaname, int inner_one, int cand_one)
	{
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
						vector<int> mid_interval; //mid_interval�洢������ϳɵ�inner �ķֶ����ӵ�
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
		1.���ܳ�����Ϊһ��
		2.�����������ͼ���Լ��Ƕ�
		3.����
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
		string conut_name = rootname + "\\PlacingResult_"+ int2string(inner_one)+".png";
		imwrite(conut_name, drawing1);
		
		cout << "inner_one: " << inner_one  << endl;
		string filepathname = "D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test" + int2string(inner_one)+".txt";
		vector<Point> con;
		for (int i = 0; i < inner_conts[inner_one].size(); i++)
			con.push_back((Point)inner_conts[inner_one][i]);
		fileout(filepathname, con);

		int num_c = 1;//ѡ��(num_c+1)*100����
		vector<CandPat> candida_contours;

		//ofstream out("D:\\VisualStudioProjects\\DihedralTessellation\\contours\\test.txt");
		//if (out.is_open())
		//{
		//	out << inner_conts[i].size() + 1 << endl;//contours[0].size()
		//	for (int j = 0; j < inner_conts[i].size(); j++)
		//		out << (int)inner_conts[i][j].x << "," << (int)inner_conts[i][j].y << endl;
		//	out << (int)inner_conts[i][0].x << "," << (int)inner_conts[i][0].y << endl;  //��β������
		//}
		//cout << "contours[0].size(): " << i << endl;
		//out.close();
		cout << inner_conts.size() << "  :  " << inner_one << endl;
		candida_contours = compare_shapes(inner_conts[inner_one], num_c);
		//cout << "candida_contours" << candida_contours.size()<< endl;

		vector<int> mid_inter = joint_relocate(inner_conts[inner_one], mid_interval_index[inner_one], num_c);
		//�����еĽ����������
		Mat drawing_pro = Mat(800, 2400, CV_8UC3, Scalar(255, 255, 255));
		CandPat tem = candida_contours[cand_one];
		prototile_second->Pro_clear();
		prototile_second->loadPoints(contour_dataset[tem.number]);
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_conts[inner_one]);

		vector<Point2f> contour_inner = prototile_mid->contour_sample[num_c]; //ѡ�����ٵĵ���бȽ�
		vector<double> contour_inner_c = curvature_com(contour_inner);// prototile_mid->contour_curva[0];
		vector<Point2f> contour_cand = CandP2Contour(tem, num_c);
		vector<double> contour_cand_c = curvature_com(contour_cand);

		//inner and cand morph into the final pettern	
		//!!!!!!����һ���н�morph��Ϊ�ֶ�!!!��ʱû�п�������ƽ�ƺ����ű�������Ϊ֮ǰ��һһ��Ӧ���еı���
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

			draw_allplane(drawing_mid, inter_mid, mid_inter, 0, 0.4);

			//����proto1�Լ������ĸ�proto2չʾ����
			Point2f cente = center_p(inter_mid);
			vector<int> return_p;
			vector<Point2f> morphed_A = extract_contour(inter_mid, mid_inter, return_p);
			//evalua_deformation();
			draw_allplane(drawing_mA, morphed_A, return_p,0,0.4);

			string filename = rootname + "\\0.";
			string file2 = filename + char(ratio + 48) + "_Cand_" + int2string(cand_one) + "tiling.png";
			filename = filename + char(ratio + 48) + "_Candidate_" + int2string(cand_one) + ".png";
			
			imwrite(filename, drawing_pro);
			imwrite(file2, drawing_mA);
			imshow("result_mid_show: ", drawing_mid);
			
		}
		return contour_inner;
	}



	bool Tiling_opt::one_situ_div(vector<int> results, vector<Point2f> contour_s, vector<Point2f> &return_B, vector<int> &return_p, Mat &countname) //���һ�ֻ�������Ľ��
	{
		Point2f line1 = contour_s[results[2]] - contour_s[results[0]];
		Point2f line2 = contour_s[results[3]] - contour_s[results[1]];

		vector<Point2f> shift_4;
		shift_4.push_back(Point2f(0, 0));
		shift_4.push_back(line1);
		shift_4.push_back(line2);
		shift_4.push_back(line1 + line2);
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
		// ��ȡΧ�ɵ�������ĿǰΪֹֻ��������ڷţ���������ת�ͷ�ת
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

		//����proto1�Լ������ĸ�proto2չʾ����
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
		//����ͨ����Χ����ôֲڵ��ཻ����Ȼ��ͨ�������ཻ���Ƿ������ײ
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
			Scalar(0, 0, 0) //��ɫ
			);
		imshow("coli 1", drawing_mid);
		//cout << "zahuishia " << endl;
		if (min1.y<0||  max1.y<0 || min1.x<0|| max1.x<0) return true;
		//cout << "min1.y " << min1.y << "  max1.y " << max1.y << " min1.x " << min1.x << " max1.x" << max1.x<<endl;
		Mat draw1 = drawing_mid(Range(min1.y, max1.y), Range(min1.x, max1.x));
		cvtColor(draw1, draw1, COLOR_BGR2GRAY);
		threshold(draw1, draw1, 128, 1, cv::THRESH_BINARY);
		
		Point2f min2 = bbox_s[0] + dis_p[2];
		Point2f max2 = bbox_s[0] + dis_p[3];
		if (min2.y<0|| max2.y<0 || min2.x<0 || max2.x<0) return true;
		//cout << "min2.y " << min2.y << "  max2.y " << max2.y << " min2.x " << min2.x << " max2.x" << max2.x<<endl;
		Mat draw2 = drawing_mid(Range(min2.y, max2.y), Range(min2.x, max2.x));
		cvtColor(draw2, draw2, COLOR_BGR2GRAY);
		threshold(draw2, draw2, 128, 1, cv::THRESH_BINARY);

		int count = 0;
		int rows = draw1.rows < draw2.rows ? draw1.rows : draw2.rows;
		int cols = draw1.cols < draw2.cols ? draw1.cols :  draw2.cols;
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

	vector<int> Tiling_opt::compare_choose_TAR(vector<Point2f> inner_c)
	{
		vector<int> all_total;
		vector<double> all_result;
		int total_num = contour_dataset.size();
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_c);
		vector<Point2f> contour_mid = prototile_mid->contour_sample[1];
		double shape_com_mid;
		vector<vector<double>> tar_mid = prototile_mid->compute_TAR(contour_mid, shape_com_mid);
		for (int can_num = 0; can_num < total_num; can_num++)
		{
			vector<vector<double>> tar_sec = all_con_tars[can_num];
			vector<pair<int, int>> path;
			int shift = 0;
			double re = tar_mismatch(tar_mid, tar_sec, path, shift,2);
			re = re / (1 + shape_com_mid + all_shape_complexity[can_num]);
			all_result.push_back(re);
			all_total.push_back(can_num);
		}
		sort_comb(all_result, all_total);
		for (int t = all_total.size() - 1; t > total_num - 30; t--)
		{
			cout << "order: " << all_total[t] << endl;
		}
		return all_total;
	}

	vector<CandPat> Tiling_opt::compare_shapes(vector<Point2f> inner_c, int num_c)
	{
		vector<int> order_total;
		prototile_mid->Pro_clear();
		prototile_mid->loadPoints(inner_c);
		vector<Point2f> contour_mid = prototile_mid->contour_sample[num_c]; //ѡ�����ٵĵ���бȽ�
		cout << "contour_mid.size()"<<contour_mid.size() << endl;
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
			//��ӵ�order������
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
		//cout �ϲ���ĺ�ѡ��
		cout << "ordersize: " << ordersize << endl;
		for (int i = 0; i < ordersize; i++)
		{
	    	cout << "order_total: " << order_total[i] << endl;
		}

		for (int t = 0; t < ordersize; t++)  //չʾ����order��ĺ�ѡͼ���ĶԱȽ��,ɸѡ����min_mismatchֵ��С�Ľ��
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

	//��360����Ѱ��ʹ�����С�ĽǶȼ���ʱ�����
	CandPat Tiling_opt::min_mismatch(vector<Point2f> inner, vector<Point2f> cand, vector<double> inner_c, vector<double> cand_c, int theone, bool isFilp)
	{
		//�������������ܳ����Ķ���
		// ����Ҫ��֤inner��cand��size��࣬��Ϊ����ֽ�β����������̫�����ɽϴ����
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
		//��cand����360�ȵ���ת
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

			//inner�ʵ�͵�һ�������������cand����ĸ���
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
					accumu_mis += quadr_mismatch(test1, test2, test11, test22, dppath); //��Ϊquadr�������õ�������100x100��������Ҫ��ȡ����
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

	double Tiling_opt::tar_mismatch(vector<vector<double>> first_arr, vector<vector<double>> second_arr, vector<pair<int, int>>& path, int &sec_shift, int width) //���Ӧƥ���ɸѡ����
	{
		int first_num = first_arr.size();
		int second_num = second_arr.size(); //first ��Ϊy�� ,secondΪx��
		if (first_num != second_num)
		{
			cout << "The sampling points of two contours are not equal " << endl;
			return 0;
		}
		int each_pnum = first_num / 2 - 1;
		double min_mis = 10000;
		vector<pair<int, int>> path_min;
		//double distance[202][202];
		//int step[202][202];//��¼�ܵĲ���
		for (int shift = 0; shift < second_num; shift++) //��first�̶����ֱ����second�����
		{
			for (int i = 0; i < first_num; i++)
			{
				for (int j = 0; j < second_num; j++)
				{
					dis[i][j] = 0;
					if (max(0, i - width) <= j && j <= min(first_num - 1, i + width))
					{
						distance[i][j] = 0;
					}
					else distance[i][j] = 100000;

				}
			}
			//distance[0][0]��¼���Ƕ���ĵ�һ����
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

		if (abs(dp[i][j] - (dp[i - 1][j - 1] + dis[i][j])) < 0.1){
			print_TAR_Path(d, dp, i - 1, j - 1, path);

		}
		else if (abs(dp[i][j] - (dp[i][j - 1] + dis[i][j])) < 0.1) {
			print_TAR_Path(d, dp, i, j - 1, path);

		}
		else {
			print_TAR_Path(d, dp, i - 1, j, path);
		}
		path.push_back(make_pair(i, j));
	}

	double Tiling_opt::quadr_mismatch(vector<Point2f> first_arr, vector<Point2f> second_arr, vector<double> first_c, vector<double> second_c, vector<pair<int, int>>& path, double zeta) //ÿ��ֻ�ܱȽ�100��������
	{
		//ofstream out("D:\\VisualStudioProjects\\manual_record\\dtw.txt");
		int first_num = first_arr.size();
		int second_num = second_arr.size();

		//cout << "\n        first.size: " << first_c[50] << "  -----    chararr_size" << second_c[50] << endl;

		//double dis[202][202];//�����֮����������
		double max_dis = 0;
		for (int i = 0; i < first_num; i++)
		{
			for (int j = 0; j < second_num; j++)
			{
				dis[i][j] = length_two_point2f(first_arr[i], second_arr[j]);
				if (dis[i][j] > max_dis) max_dis = dis[i][j];
			}
		}
		//�����ֵ��һ��
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

		//double dis_cur[202][202];//�����֮������ʲ���
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
		//���ʲ�ֵ��һ��
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
		//int step[202][202];//��¼�ܵĲ���
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
		//��֤�߶��㹻�����󽻵㣬���߶γ��ȷŴ�3��
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
		//���ض���㣬����ֻ��һ���������ж��,ɾ���ظ���

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

	vector<int> Tiling_opt::joint_relocate(vector<Point2f> contour_, vector<int> joint_index, int num_c) //��ԭʼ�����ϵĻ��ֵ��Ӧ���������������
	{
		vector<Point2f> contour_mid = sampling(contour_, num_c+1);//ѡ�����ٵĵ���бȽ�
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
	
	vector<Point2f> Tiling_opt::morphing_2_patterns(vector<Point2f> &contour1, vector<Point2f> &contour2, vector<double> &concur1, vector<double> &concur2, vector<int> &mid_inter, float shape_ratio)
	{
		//��ͳmorphing����start��end���һϵ��intermediate״̬��������ͨ��c1��c2������յı��ν��

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


		//// �����ҵ�һһ��Ӧ�ĵ�����
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

		//			if ((qmpath[f].first < i) || qmpath[f].second < j) //�����һ��·������һ�˵��Ƿ��ѱ�ʹ��
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


		/*	//---------------------����Ϊmorphing�׶�
		// �����ҵ�һһ��Ӧ�ĵ�����
		vector<Point2f> src1_points;
		vector<Point2f> src2_points;
		src1_points.push_back(output_first[0]);
		src2_points.push_back(output_second[0]);

		//for (int i = 0; i < dp_path.size(); i++)
		//{
		//	cout << dp_path[i].first << " -- " << dp_path[i].second << endl;
		//}
		if (output_first.size() < output_second.size())
		{
		int f = 1;
		int j = 0;
		for (int i = 1; i < output_first.size() - 1;)
		{

		double min = 10000;
		while (f < dp_path.size())
		{

		if ((dp_path[f].first < i) || dp_path[f].second < j) //�����һ��·������һ�˵��Ƿ��ѱ�ʹ��
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
		for (int i = 1; i < output_second.size() - 1;)
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

		if (src1_points.size() == src2_points.size())
		{
		for (int i = 0; i < src1_points.size(); i++)
		{
		circle(drawing_src3, src1_points[i], 1, Scalar(255, 0, 0), -1);
		circle(drawing_src3, src2_points[i], 1, Scalar(0, 0, 255), -1);
		MyLine(drawing_src3, src1_points[i], src2_points[i], "orange");
		//cout << src1_points[i] << " - " << src2_points[i] << endl;
		}
		}
		else cout << "num error!" << endl;
		//����ϵ���ı�仯�̶ȣ���Ϊ�Ż��ı���
		ImageMorphing(drawing_src1, src1_points, drawing_src2, src2_points, drawing_dst, output_final, 0.5);
		//char i;

		double length_final = length_two_point2f(output_final[0], output_final[output_final.size() - 1]);
		for (int i = 0; i < output_final.size() - 1; i++)
		{
		MyLine(drawing_src3, output_final[i], output_final[i + 1], "red");
		}

		//imshow("drawing_src1", drawing_src1);
		//imshow("drawing_src2", drawing_src2);
		//imshow("drawing_dst", drawing_dst);

		string name = "the ";
		name = name + char(i + 48) + " pair: ";
		imshow(name, drawing_src3);

		*/

		return final_pettern;
	}

	double evalua_deformation(vector<vector<Point2f>> contour, vector<vector<double>> curvature)
	{
		//����˳��Ϊ mid sec result
		double total_score = 0;
		double score_mid_r = 0;
		double score_sec_r = 0;
		score_mid_r = contourArea(contour[2]) / contourArea(contour[0]);
		score_sec_r = contourArea(contour[1]) / contourArea(contour[0]);
		return total_score;
	}





	/*
	//ÿһ�������ֳ��Ķ�
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
	////����
	//double factor = com_scale_factor();
	//for (int i = 0; i <contour2.size(); i++)
	//{
	//	contour2[i].x = contour2[i].x * factor;
	//	contour2[i].y = contour2[i].y * factor;
	//}

	//int interval_num1_first = contour1.size() / 10;


	////int intervals_second = 4; //�ݶ���Ϊ�Ķ�
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


	//		int flag = 1; //0 ����sec,�����޷���,1 ����ref
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
	//		//�÷���仯�������������룬������ʼλ�õ�Ӱ��

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
	////����ȡǰ32��
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

	//////�ѽ����ʾ����
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
	////	int size = 800;//չʾͼƬ�ĳߴ�

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

	//////��shift������ͼ���뵽ͼ��������е�
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
	//�÷���仯�������������룬������ʼλ�õ�Ӱ��

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
	int flag = 0; //0 ����sec,1 ����ref
	double score = com_tra_sim(proto_interval_1[i], proto_first_char[i], proto_interval_2[j], proto_second_char[j], flag);
	each_score = make_pair(score, flag);
	score_order[i][j] = each_score;
	}
	}

	// ȫ���о�
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

	// ���������������
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


	double Tiling_opt::com_tra_sim(vector<Point2f> &first_interval, vector<char>&first_char, vector<Point2f> &second_interval, vector<char>&second_char, int &flag) //compute trajectory similarity by DTW(Dynamic time warping)����һ�Ե����Ƶ÷�
	{
	double score = 0;
	double score2 = 0;

	vector<Point2f> first_after_aff;
	vector<Point2f> second_after_aff;
	vector<Point2f> second_after_aff_ref;

	vector<char> second_char_out;
	//�÷���仯�������������룬������ʼλ�õ�Ӱ�죬����warpAff_tra��warpAff_tra_ref_y����Ӱ���Ӧ�����ַ�����˳��

	double re_first = warpAff_tra(first_interval, first_after_aff);
	double re_second = warpAff_tra_sec(second_interval, second_after_aff, second_char, second_char_out);
	double re_second_ref = warpAff_tra_ref_y(second_interval, second_after_aff_ref);
	//int col = 800;//չʾͼƬ�ĳߴ�
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
	out << "0 ����sec,1 ����ref" << endl;
	prototile_first->loadTileData(imagename1);
	prototile_second->loadTileData(imagename2);

	//ѡ��100��Ĳ�������
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
	//cout << "contour1_length(): " << contour_length(contour1) << endl << "contour2_length(): " << contour_length(contour2) << endl;
	//����������һ�����ţ������������ʾ�׶�
	//double scale = 0.5;
	//for (int i = 0; i < contour1.size(); i++)
	//{
	//contour1[i].x = contour1[i].x * scale;
	//contour1[i].y = contour1[i].y * scale;
	//}
	//for (int i = 0; i < contour2.size(); i++)
	//{
	//contour2[i].x = contour2[i].x * scale;
	//contour2[i].y = contour2[i].y * scale;
	//}
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

	//����ַ��������ʸ����Ƿ��Ӧ
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

	double score = com_optimal_score(proto_interval_first, proto_first_char, proto_interval_second, proto_second_char, order); //����һ������˳���Լ����������μ�¼��order[i].first������1�е�i�ζ�Ӧ��2�еĶ�����order[i].second:ƥ�������Ƿ�ת��������


	//Mat drawing41 = Mat(800, 1600, CV_8UC3, Scalar(255, 255, 255)); //��ɫ����
	////��shift������ͼ���뵽ͼ��������е�
	//Point2f shift11 = Point2f(prototile_first->center_point.x - 400, prototile_first->center_point.y - 400);
	//Point2f shift21 = prototile_second->center_point - Point2f(1 * 400 + 700, 400);
	//for (int i = 0; i < 4; i++)
	//{
	//	//cout << "proto_interval_first[" << i << "]  num:" << proto_interval_first[i].size() << endl;
	//	string color;
	//	if (i == 0) color = "red";
	//	else if (i == 1) color = "blue";
	//	else if (i == 2) color = "green";
	//	else if (i == 3) color = "cyan";

	//	for (int j = 0; j < proto_interval_first[i].size() - 1; j++)
	//	{
	//		MyLine(drawing41, proto_interval_first[i][j] - shift11, proto_interval_first[i][j + 1] - shift11, color);

	//	}
	//	MyLine(drawing41, proto_interval_first[i][0] - shift11, proto_interval_first[i][proto_interval_first[i].size() - 1] - shift11, "orange");
	//	for (int j = 0; j < proto_interval_second[order[i].first].size() - 1; j++)
	//	{
	//		MyLine(drawing41, proto_interval_second[order[i].first][j] - shift21, proto_interval_second[order[i].first][j + 1] - shift21, color);

	//	}
	//	MyLine(drawing41, proto_interval_second[i][0] - shift21, proto_interval_second[i][proto_interval_second[i].size() - 1] - shift21, "orange");
	//}
	//imshow("final111 " + imagename1 + " and  " + imagename2, drawing41);

	//������Ž�������¼���õ����з�ʽ����������һ�����⣺��һ��com_optimal�������ͨ������仯һ��ĳһ�˵����Ľ��������չʾ��ʱ�����ǰѶ�Ӧ���е����
	vector<Point2f> output_first;
	vector<Point2f> output_second;
	vector<Point2f> output_final;
	vector<char> sec_char_out;
	//vector<Point2f> output_second_ref;

	//��������б��εĻ�ԭ�����Խ���һ��for����while����ֹ����
	//�任��ȥ��һ����û�ж���0��
	//dtw��dp����һ����û��
	for (int t = 0; t < 4; t++)
	{
	int i = t % 4;
	//all_types[max_x][max_y].second[i];
	Point2f extre[2][2];
	extre[0][0] = proto_interval_first[i][0];
	extre[0][1] = proto_interval_first[i][proto_interval_first[i].size()-1];
	extre[1][0] = proto_interval_second[order[i].first][0];
	extre[1][1] = proto_interval_second[order[i].first][proto_interval_second[order[i].first].size()-1];
	Point2f sae[2][2];
	//for (int j = 0; j < proto_interval_first[i].size(); j++)
	//{
	//cout << "proto_interval_first[i]: " << proto_interval_first[i][j]  << endl;
	//}
	//for (int j = 0; j < proto_interval_first[i].size(); j++)
	//{
	//cout << "proto_interval_second[i]: " << proto_interval_second[order[i].first][j] << endl;

	//}

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

	int size = 800;//չʾͼƬ�ĳߴ�

	Point2f shift1 = Point2f(size / 2 - re_first, size / 2);
	Point2f shift2 = Point2f(size / 2 - re_second, size / 2);
	//Mat drawing = Mat::zeros(size, size, CV_8UC3);  //��ɫ����
	Mat drawing_src1 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));  // ��ɫ����
	Mat drawing_src2 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));
	Mat drawing_src3 = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));
	Mat drawing_dst = Mat(size, size, CV_8UC3, Scalar(255, 255, 255));

	// ��һ����output���뵽�м�λ��
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
	MyLine(drawing_src3, output_first[t], output_first[t + 1], "green");  //�ܵ�չʾ
	}
	for (int t = 0; t < output_second.size() - 1; t++)
	{
	//	MyLine(drawing, first_interval[i], first_interval[i + 1], 1);
	MyLine(drawing_src2, output_second[t], output_second[t + 1], "orange");
	MyLine(drawing_src3, output_second[t], output_second[t + 1], "orange");  //�ܵ�չʾ
	}

	cout << "dp_path: " << dp_path.size() << endl;
	//���¼���
	DTW(output_first, output_second); // �����Ƿ�����ַ����Ƚ�

	for (int i = 0; i < dp_path.size(); i++)
	{
	MyLine(drawing_src3, output_first[dp_path[i].first], output_second[dp_path[i].second], "grey");
	}

	//---------------------����Ϊmorphing�׶�
	// �����ҵ�һһ��Ӧ�ĵ�����
	vector<Point2f> src1_points;
	vector<Point2f> src2_points;
	src1_points.push_back(output_first[0]);
	src2_points.push_back(output_second[0]);

	//for (int i = 0; i < dp_path.size(); i++)
	//{
	//	cout << dp_path[i].first << " -- " << dp_path[i].second << endl;
	//}
	if (output_first.size() < output_second.size())
	{
	int f = 1;
	int j = 0;
	for (int i = 1; i < output_first.size() - 1;)
	{

	double min = 10000;
	while (f < dp_path.size())
	{

	if ((dp_path[f].first < i) || dp_path[f].second < j) //�����һ��·������һ�˵��Ƿ��ѱ�ʹ��
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
	for (int i = 1; i < output_second.size() - 1;)
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

	if (src1_points.size() == src2_points.size())
	{
	for (int i = 0; i < src1_points.size(); i++)
	{
	circle(drawing_src3, src1_points[i], 1, Scalar(255, 0, 0), -1);
	circle(drawing_src3, src2_points[i], 1, Scalar(0, 0, 255), -1);
	MyLine(drawing_src3, src1_points[i], src2_points[i], "orange");
	//cout << src1_points[i] << " - " << src2_points[i] << endl;
	}
	}
	else cout << "num error!" << endl;
	//����ϵ���ı�仯�̶ȣ���Ϊ�Ż��ı���
	ImageMorphing(drawing_src1, src1_points, drawing_src2, src2_points, drawing_dst, output_final, 0.5);
	//char i;

	double length_final = length_two_point2f(output_final[0], output_final[output_final.size()-1]);
	for (int i = 0; i < output_final.size()-1; i++)
	{
	MyLine(drawing_src3, output_final[i], output_final[i+1], "red");
	}

	//imshow("drawing_src1", drawing_src1);
	//imshow("drawing_src2", drawing_src2);
	//imshow("drawing_dst", drawing_dst);

	string name = "the ";
	name = name + char(i + 48) + " pair: ";
	imshow(name, drawing_src3);


	//�����ｫout_finalӳ���ԭ����λ�ã����������˵������
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
	//�˴���output_final���뵽0��
	Point2f shift3 = Point2f(length_final / 2 - size / 2, - size / 2);
	//for (int i = 0; i < output_final.size(); i++)
	//{
	//	cout << "output_final_origin[i]: " << output_final[i] << endl;
	//}
	for (int i = 0; i < output_final.size(); i++)
	{
	output_final[i] = output_final[i] + shift3;
	}
	//for (int i = 0; i < output_final.size(); i++)
	//{
	//	cout << "output_final[i]: " << output_final[i] << endl;
	//}
	re_warp_Aff(output_final, out1, sae[0][0], sae[0][1]);
	//���ת��ǰ���ı仯
	//cout << "out1[0]: " << out1[0] << "  out1[1]: " << out1[out1.size()-1] << endl;
	//cout << "out1[0]: " << proto_interval_first[i][0] << "  out1[1]: " << proto_interval_first[i][proto_interval_first[i].size()-1] << endl;
	proto_interval_first[i].swap(out1);
	if (order[i].second == 0)  //��ͬ��ֵ���ò�ͬ�ĺ���
	{
	vector<Point2f> out2;
	re_warp_Aff_sec(output_final, out2, sae[1][0], sae[1][1]);
	proto_interval_second[order[i].first].swap(out2);
	}
	else if (order[i].second == 1)
	{
	vector<Point2f> out2;
	re_warp_Aff_ref(output_final, out2, sae[1][0], sae[1][1]);
	proto_interval_second[order[i].first].swap(out2);
	}

	//�����κ�ı߻�ԭ��ԭ��λ�ã����������߽���΢��
	Point2f start;
	Point2f end;
	//cout << "out1[0]: " << proto_interval_first[1][0] << "  out1[1]: " << proto_interval_first[1][proto_interval_first[1].size() - 1] << endl;

	if (i == 0)
	{
	vector<Point2f> output_;
	warpAff_sca(proto_interval_first[3], output_, proto_interval_first[3][0], proto_interval_first[i][0]);
	proto_interval_first[3].swap(output_);
	output_.swap(vector<Point2f>());
	warpAff_sca(proto_interval_first[1], output_, proto_interval_first[i][proto_interval_first[i].size() - 1], proto_interval_first[1][proto_interval_first[1].size() - 1]);
	proto_interval_first[1].swap(output_);
	}
	else
	{
	vector<Point2f> output_;
	warpAff_sca(proto_interval_first[i-1], output_, proto_interval_first[i-1][0], proto_interval_first[i][0]);
	proto_interval_first[i-1].swap(output_);
	output_.swap(vector<Point2f>());
	warpAff_sca(proto_interval_first[(i + 1) % 4], output_, proto_interval_first[i][proto_interval_first[i].size() - 1], proto_interval_first[(i + 1) % 4][proto_interval_first[(i + 1) % 4].size() - 1]);
	proto_interval_first[(i + 1) % 4].swap(output_);
	}
	//��i=4ʱ�˴�������
	if (order[i].first == 0)
	{
	vector<Point2f> output_;
	cout << "hahahah" << endl;
	warpAff_sca(proto_interval_second[3], output_, proto_interval_second[3][0], proto_interval_second[order[i].first][0]);
	proto_interval_second[3].swap(output_);
	output_.swap(vector<Point2f>());
	cout << "hahahah" << endl;
	warpAff_sca(proto_interval_second[1], output_, proto_interval_second[order[i].first][proto_interval_second[order[i].first].size() - 1], proto_interval_second[1][proto_interval_second[1].size() - 1]);
	proto_interval_second[1].swap(output_);
	}
	else
	{
	vector<Point2f> output_;
	warpAff_sca(proto_interval_second[order[i].first - 1], output_, proto_interval_second[order[i].first - 1][0], proto_interval_second[order[i].first][0]);
	proto_interval_second[order[i].first - 1].swap(output_);
	output_.swap(vector<Point2f>());
	warpAff_sca(proto_interval_second[(order[i].first + 1) % 4], output_, proto_interval_second[order[i].first][proto_interval_second[order[i].first].size() - 1], proto_interval_second[(order[i].first + 1) % 4][proto_interval_second[(order[i].first + 1) % 4].size() - 1]);
	proto_interval_second[(order[i].first + 1) % 4].swap(output_);
	}
	//cout << "out1_[0]: " << proto_interval_first[1][0] << "  out1[1]: " << proto_interval_first[1][proto_interval_first[1].size() - 1] << endl;

	//�ҵ���һ���ߵ���һ�������Ǹı���sae��
	}

	int row = 800;
	int col = 1600;
	//Mat drawing4 = Mat::zeros(row, col, CV_8UC3); //��ɫ����
	Mat drawing4 = Mat(row, col, CV_8UC3, Scalar(255, 255, 255)); //��ɫ����
	//��shift������ͼ���뵽ͼ��������е�
	Point2f shift1 = Point2f(prototile_first->center_point.x - col / 4, prototile_first->center_point.y - row / 2);
	Point2f shift2 = prototile_second->center_point - Point2f(1 * col / 4+700, row / 2);
	for (int i = 0; i < 4; i++)
	{
	//cout << "proto_interval_first[" << i << "]  num:" << proto_interval_first[i].size() << endl;
	string color;
	if (i == 0) color = "red";
	else if (i == 1) color = "blue";
	else if (i == 2) color = "green";
	else if (i == 3) color = "cyan";

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
	imshow("final " + imagename1 + " and  " + imagename2, drawing4);

	//�������ڷţ��õ����յ�����ͼ
	int num1 = 800;
	int num2 = 800;
	//Mat drawing5 = Mat(num1, num2, CV_8UC3, Scalar(255, 255, 255)); //��ɫ����
	Mat drawing5 = Mat(num1, num2, CV_8UC3, Scalar(0, 0, 0)); //��ɫ����
	double scale_ = 0.4;

	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < proto_interval_first[i].size(); j++)
	{
	proto_interval_first[i][j] = proto_interval_first[i][j] * scale_;
	}
	for (int j = 0; j < proto_interval_second[i].size(); j++)
	{
	proto_interval_second[i][j] = proto_interval_second[i][j] * scale_;
	}
	}

	//�˴���Ҫ��һ���ж�proto1��ν������У��Ƿ��з�ת�����Ը������յõ���1��2�Ķ�Ӧ��ϵ�õ���������ʱĬ��Ϊ˳������
	Point2f hori_axis_length = proto_interval_first[2][0] - proto_interval_first[0][0];
	Point2f vert_axis_length = proto_interval_first[3][0] - proto_interval_first[1][0];
	int fff_sign[30][30] = {0};

	Point2f sign[4];
	sign[0] = hori_axis_length;
	sign[1] = -hori_axis_length;
	sign[2] = vert_axis_length;
	sign[3] = -vert_axis_length;

	vector<Point2f> prototwoAff_place;
	vector<Point2f> protoFirst_bbx;//����5���㣬center��up,down,left and right
	pair<int, int> protoFirst;
	vector<vector<Point2f>> bbx_all;
	vector<pair<int, int>> sign_all;
	vector<vector<Point2f>> proto_first;
	bbx_center_point(proto_interval_first, protoFirst_bbx);
	bbx_all.push_back(protoFirst_bbx);
	int n1 = 14;
	int n2 = 14;
	fff_sign[n1][n2] = 1;
	protoFirst = make_pair(n1, n2);
	sign_all.push_back(protoFirst);
	int type = 0;//����proto1�ĸ�������
	int g = 0;
	while (!bbx_all.empty())
	{
	g++;
	cout << "g: " << g << endl;
	vector<Point2f> first_bbx;
	pair<int, int> first_sign;
	if (type == 0)
	{
	first_bbx = bbx_all[bbx_all.size() - 1];
	bbx_all.pop_back();
	first_sign = sign_all[sign_all.size() - 1];
	sign_all.pop_back();
	cout << "first_bbx: " << first_bbx.size() << endl;
	cout << "bbx_all:  " << bbx_all.size() << endl;
	}

	proto_first.swap(vector<vector<Point2f>>());
	Point2f diff = first_bbx[0] - protoFirst_bbx[0];
	cout << "diff: " << diff << endl;
	//�õ���proto1������
	for (int i = 0; i < 4; i++)
	{
	vector<Point2f> first_inter;
	for (int j = 0; j < proto_interval_first[i].size(); j++)
	{
	first_inter.push_back(proto_interval_first[i][j]+diff);
	}
	proto_first.push_back(first_inter);
	}

	//����proto1�������ĸ�proto1�����ջ
	for (int i = 0; i < 4; i++)
	{
	pair<int, int> sign_;
	if (i == 0) sign_ = make_pair(first_sign.first + 1, first_sign.second);
	if (i == 1) sign_ = make_pair(first_sign.first - 1, first_sign.second);
	if (i == 2) sign_ = make_pair(first_sign.first, first_sign.second + 1);
	if (i == 3) sign_ = make_pair(first_sign.first, first_sign.second - 1);
	if (fff_sign[sign_.first][sign_.second] == 0)
	{
	vector<Point2f> firstbbx;
	int fff = 0;
	for (int t = 0; t < 5; t++)
	{
	Point2f sigh_re = first_bbx[t] + sign[i];
	firstbbx.push_back(sigh_re);
	if (firstbbx[t].x < num1 && firstbbx[t].x>0 && firstbbx[t].y>0 && firstbbx[t].y<num2)
	{
	fff = 1;
	}
	}
	if (fff == 1)
	{
	bbx_all.push_back(firstbbx);
	sign_all.push_back(sign_);
	fff_sign[sign_.first][sign_.second] = 1;
	}
	}

	cout << "bbx_all:  " << bbx_all.size() << endl;
	}
	//����proto1�Լ������ĸ�proto2չʾ����
	vector<Point2f> pro_st;
	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < proto_first[i].size() - 1; j++)
	{
	pro_st.push_back(proto_first[i][j]);
	}
	}
	int n = pro_st.size();
	cout << "n: " << n << endl;
	Point rook_points[1][800];
	for (int t = 0; t < n; t++)
	{
	rook_points[0][t] = pro_st[t];
	}
	const Point* ppt[1] = { rook_points[0] };
	int npt[] = { n };
	fillPoly(drawing5,
	ppt,
	npt,
	1,
	//Scalar(0, 0, 0) //��ɫ
	Scalar(255, 255, 255) //��ɫ
	);
	for (int i = 0; i < 4; i++)
	{
	prototwoAff_place.swap(vector<Point2f>());
	for (int j = 0; j < proto_first[i].size() - 1; j++)
	{
	MyLine(drawing5, proto_first[i][j], proto_first[i][j + 1], "black");
	}
	}

	//for (int i = 0; i < 4; i++)
	//{
	//	prototwoAff_place.swap(vector<Point2f>());
	//	for (int j = 0; j < proto_first[i].size() - 1; j++)
	//	{
	//		//pro_st.push_back(proto_first[i][j]);
	//		MyLine(drawing5, proto_first[i][j], proto_first[i][j + 1], "blue");
	//	}
	//	Aff_place(proto_first[i], proto_interval_second[order[i].first], proto_interval_second, prototwoAff_place, order[i].second);
	//	int n = prototwoAff_place.size();
	//	cout << "n: " << n << endl;
	//	Point rook_points[1][800];
	//	for (int t = 0; t < n; t++)
	//	{
	//		rook_points[0][t] = prototwoAff_place[t];
	//	}
	//	const Point* ppt[1] = { rook_points[0] };
	//	int npt[] = { n };
	//	fillPoly(drawing5,
	//		ppt,
	//		npt,
	//		1,
	//		Scalar(0, 0, 0) //��ɫ
	//		//Scalar(255, 255, 255) //��ɫ
	//		);
	//}
	//if (g == 2)break;
	}
	imshow("result", drawing5);
	imwrite("D:\\images\\111\\111.png", drawing5);


	return 0;
	}
	*/
}
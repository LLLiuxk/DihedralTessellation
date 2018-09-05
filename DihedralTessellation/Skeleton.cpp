#include "tilingOpt.h"

using namespace std;

namespace Tiling_tiles{
	/**
	* @brief ������ͼ�����ϸ��,������
	* @param srcΪ����ͼ��,��cvThreshold�����������8λ�Ҷ�ͼ���ʽ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
	* @param maxIterations���Ƶ���������������������ƣ�Ĭ��Ϊ-1���������Ƶ���������ֱ��������ս��
	* @return Ϊ��srcϸ��������ͼ��,��ʽ��src��ʽ��ͬ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
	*/
	Mat thinImage(const cv::Mat & src, const int maxIterations = -1)
	{
		assert(src.type() == CV_8UC1);
		Mat dst;
		int width = src.cols;
		int height = src.rows;
		src.copyTo(dst);
		int count = 0;  //��¼��������  
		while (true)
		{
			count++;
			if (maxIterations != -1 && count > maxIterations) //���ƴ������ҵ�����������  
				break;
			vector<uchar *> mFlag; //���ڱ����Ҫɾ���ĵ�  
			//�Ե���  
			for (int i = 0; i < height; ++i)
			{
				uchar * p = dst.ptr<uchar>(i);
				for (int j = 0; j < width; ++j)
				{
					//��������ĸ����������б��  
					//  p9 p2 p3  
					//  p8 p1 p4  
					//  p7 p6 p5  
					uchar p1 = p[j];
					if (p1 != 1) continue;
					uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
					uchar p8 = (j == 0) ? 0 : *(p + j - 1);
					uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);
					uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);
					uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);
					uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);
					uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);
					uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);

					if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
					{
						int ap = 0;
						if (p2 == 0 && p3 == 1) ++ap;
						if (p3 == 0 && p4 == 1) ++ap;
						if (p4 == 0 && p5 == 1) ++ap;
						if (p5 == 0 && p6 == 1) ++ap;
						if (p6 == 0 && p7 == 1) ++ap;
						if (p7 == 0 && p8 == 1) ++ap;
						if (p8 == 0 && p9 == 1) ++ap;
						if (p9 == 0 && p2 == 1) ++ap;

						if (ap == 1 && p2 * p4 * p8 == 0 && p2 * p6 * p8 == 0)
						{
							//���  
							mFlag.push_back(p + j);
						}
					}
				}
			}
			//����ǵĵ�ɾ��  
			for (vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
			{
				**i = 0;
			}
			//ֱ��û�е����㣬�㷨����  
			if (mFlag.empty())
			{
				break;
			}
			else
			{
				mFlag.clear();//��mFlag���  
			}
			
			//�Ե���  
			for (int i = 0; i < height; ++i)
			{
				uchar * p = dst.ptr<uchar>(i);
				for (int j = 0; j < width; ++j)
				{
					//��������ĸ����������б��  
					//  p9 p2 p3  
					//  p8 p1 p4  
					//  p7 p6 p5  
					uchar p1 = p[j];
					if (p1 != 1) continue;
					uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
					uchar p8 = (j == 0) ? 0 : *(p + j - 1);
					uchar p2 = (i == 0) ? 0 : *(p - dst.step + j);
					uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - dst.step + j + 1);
					uchar p9 = (i == 0 || j == 0) ? 0 : *(p - dst.step + j - 1);
					uchar p6 = (i == height - 1) ? 0 : *(p + dst.step + j);
					uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + dst.step + j + 1);
					uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + dst.step + j - 1);
					if ((p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) >= 2 && (p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) <= 6)
					{
						int ap = 0;
						if (p2 == 0 && p3 == 1) ++ap;
						if (p3 == 0 && p4 == 1) ++ap;
						if (p4 == 0 && p5 == 1) ++ap;
						if (p5 == 0 && p6 == 1) ++ap;
						if (p6 == 0 && p7 == 1) ++ap;
						if (p7 == 0 && p8 == 1) ++ap;
						if (p8 == 0 && p9 == 1) ++ap;
						if (p9 == 0 && p2 == 1) ++ap;

						if (ap == 1 && p2 * p4 * p6 == 0 && p4 * p6 * p8 == 0)
						{
							//���  
							mFlag.push_back(p + j);
						}
					}
				}
			}

			//����ǵĵ�ɾ��  
			for (vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
			{
				**i = 0;
			}

			//ֱ��û�е����㣬�㷨����  
			if (mFlag.empty())
			{
				break;
			}
			else
			{
				mFlag.clear();//��mFlag���  
			}

		}
		return dst;
	}

	/**
	* @brief �Թ�����ͼ���ݽ��й��ˣ�ʵ��������֮�����ٸ�һ���հ�����
	* @param thinSrcΪ����Ĺ�����ͼ��,8λ�Ҷ�ͼ���ʽ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
	*/
	void filterOver(Mat thinSrc)
	{
		assert(thinSrc.type() == CV_8UC1);
		int width = thinSrc.cols;
		int height = thinSrc.rows;
		for (int i = 0; i < height; ++i)
		{
			uchar * p = thinSrc.ptr<uchar>(i);
			for (int j = 0; j < width; ++j)
			{
				// ʵ��������֮�����ٸ�һ������
				//  p9 p2 p3  
				//  p8 p1 p4  
				//  p7 p6 p5  
				uchar p1 = p[j];
				if (p1 != 1) continue;
				uchar p4 = (j == width - 1) ? 0 : *(p + j + 1);
				uchar p8 = (j == 0) ? 0 : *(p + j - 1);
				uchar p2 = (i == 0) ? 0 : *(p - thinSrc.step + j);
				uchar p3 = (i == 0 || j == width - 1) ? 0 : *(p - thinSrc.step + j + 1);
				uchar p9 = (i == 0 || j == 0) ? 0 : *(p - thinSrc.step + j - 1);
				uchar p6 = (i == height - 1) ? 0 : *(p + thinSrc.step + j);
				uchar p5 = (i == height - 1 || j == width - 1) ? 0 : *(p + thinSrc.step + j + 1);
				uchar p7 = (i == height - 1 || j == 0) ? 0 : *(p + thinSrc.step + j - 1);
				if (p2 + p3 + p8 + p9 >= 1)
				{
					p[j] = 0;
				}
			}
		}
	}

	/**
	* @brief �ӹ��˺�Ĺ�����ͼ����Ѱ�Ҷ˵�ͽ����
	* @param thinSrcΪ����Ĺ��˺������ͼ��,8λ�Ҷ�ͼ���ʽ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
	* @param raudis����뾶���Ե�ǰ���ص�λԲ�ģ���Բ��Χ���жϵ��Ƿ�Ϊ�˵�򽻲��
	* @param thresholdMax�������ֵ���������ֵΪ�����
	* @param thresholdMin�˵���ֵ��С�����ֵΪ�˵�
	* @return Ϊ��srcϸ��������ͼ��,��ʽ��src��ʽ��ͬ��Ԫ����ֻ��0��1,1������Ԫ�أ�0����Ϊ�հ�
	*/
	vector<Point2f> getPoints(const Mat &thinSrc, unsigned int raudis = 4, unsigned int thresholdMax = 6, unsigned int thresholdMin = 4)
	{
		assert(thinSrc.type() == CV_8UC1);
		int width = thinSrc.cols;
		int height = thinSrc.rows;
		Mat tmp;
		thinSrc.copyTo(tmp);
		vector<Point2f> points;
		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				if (*(tmp.data + tmp.step * i + j) == 0)
				{
					continue;
				}
				int count = 0;
				for (int k = i - raudis; k < i + raudis + 1; k++)
				{
					for (int l = j - raudis; l < j + raudis + 1; l++)
					{
						if (k < 0 || l < 0 || k>height - 1 || l>width - 1)
						{
							continue;

						}
						else if (*(tmp.data + tmp.step * k + l) == 1)
						{
							count++;
						}
					}
				}

				if (count > thresholdMax || count<thresholdMin)
				{
					Point2f Point2f(j, i);
					points.push_back(Point2f);
				}
			}
		}
		return points;
	}


	vector<Point2f> get_Skeleton(string imaname, vector<Point2f> &skeleton)
	{
		Mat src;
		//��ȡͼ��  
		string name = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + imaname + ".png";

		src = imread(name, cv::IMREAD_GRAYSCALE);

		if (src.empty())
		{
			cout << "��ȡ�ļ�ʧ�ܣ�" << std::endl;
			exit;
		}

		//��ԭͼ��ת��Ϊ��ֵͼ��  
		threshold(src, src, 128, 1, cv::THRESH_BINARY);
		//cout << (int)src.at<uchar>(155, 155) << endl;
		for (int i = 0; i < src.rows; i++)
		{
			for (int j = 0; j < src.cols; j++)
			{
				//cout << (int)src.at<uchar>(i, j);
				if (src.at<uchar>(i, j) == 1) src.at<uchar>(i, j) = (uchar)0;
				else if (src.at<uchar>(i, j) == 0) src.at<uchar>(i, j) = (uchar)1;
			}
		}
		//ͼ��ϸ����������  
		Mat dst = thinImage(src);
		
		//����ϸ�����ͼ��
		//filterOver(dst);
		
		//���Ҷ˵�ͽ����  
		vector<Point2f> points = getPoints(dst, 6, 34, 12);

		//��ֵͼת���ɻҶ�ͼ���������ҵ��ĵ�
		dst = dst * 255;
		src = src * 255;

		for (int i = 0; i < dst.rows; i++)
		{
			for (int j = 0; j < dst.cols; j++)
			{
				//cout << (int)src.at<uchar>(i, j);
				if (dst.at<uchar>(i, j) == 255) skeleton.push_back(Point2f(j,i));
			}
		}
		/*for (int i = 1; i < contoursize; i++)
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
		}*/

		vector<Point2f>::iterator it = points.begin();
		for (; it != points.end(); it++)
		{
			circle(dst, *it, 4, 255, -1);
		}
		//imwrite("dst.jpg", dst);

		//��ʾͼ��  
		cv::namedWindow("src1", CV_WINDOW_AUTOSIZE);
		cv::namedWindow("dst1", CV_WINDOW_AUTOSIZE);
		cv::imshow("src1", src);
		cv::imshow("dst1", dst);
		return points;
	}

}
#include "tilingOpt.h"

using namespace std;

namespace Tiling_tiles{
	/**
	* @brief 对输入图像进行细化,骨骼化
	* @param src为输入图像,用cvThreshold函数处理过的8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白
	* @param maxIterations限制迭代次数，如果不进行限制，默认为-1，代表不限制迭代次数，直到获得最终结果
	* @return 为对src细化后的输出图像,格式与src格式相同，元素中只有0与1,1代表有元素，0代表为空白
	*/
	Mat thinImage(const cv::Mat & src, const int maxIterations = -1)
	{
		assert(src.type() == CV_8UC1);
		Mat dst;
		int width = src.cols;
		int height = src.rows;
		src.copyTo(dst);
		int count = 0;  //记录迭代次数  
		while (true)
		{
			count++;
			if (maxIterations != -1 && count > maxIterations) //限制次数并且迭代次数到达  
				break;
			vector<uchar *> mFlag; //用于标记需要删除的点  
			//对点标记  
			for (int i = 0; i < height; ++i)
			{
				uchar * p = dst.ptr<uchar>(i);
				for (int j = 0; j < width; ++j)
				{
					//如果满足四个条件，进行标记  
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
							//标记  
							mFlag.push_back(p + j);
						}
					}
				}
			}
			//将标记的点删除  
			for (vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
			{
				**i = 0;
			}
			//直到没有点满足，算法结束  
			if (mFlag.empty())
			{
				break;
			}
			else
			{
				mFlag.clear();//将mFlag清空  
			}
			
			//对点标记  
			for (int i = 0; i < height; ++i)
			{
				uchar * p = dst.ptr<uchar>(i);
				for (int j = 0; j < width; ++j)
				{
					//如果满足四个条件，进行标记  
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
							//标记  
							mFlag.push_back(p + j);
						}
					}
				}
			}

			//将标记的点删除  
			for (vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
			{
				**i = 0;
			}

			//直到没有点满足，算法结束  
			if (mFlag.empty())
			{
				break;
			}
			else
			{
				mFlag.clear();//将mFlag清空  
			}

		}
		return dst;
	}

	/**
	* @brief 对骨骼化图数据进行过滤，实现两个点之间至少隔一个空白像素
	* @param thinSrc为输入的骨骼化图像,8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白
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
				// 实现两个点之间至少隔一个像素
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
	* @brief 从过滤后的骨骼化图像中寻找端点和交叉点
	* @param thinSrc为输入的过滤后骨骼化图像,8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白
	* @param raudis卷积半径，以当前像素点位圆心，在圆范围内判断点是否为端点或交叉点
	* @param thresholdMax交叉点阈值，大于这个值为交叉点
	* @param thresholdMin端点阈值，小于这个值为端点
	* @return 为对src细化后的输出图像,格式与src格式相同，元素中只有0与1,1代表有元素，0代表为空白
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
		//获取图像  
		string name = "D:\\VisualStudioProjects\\DihedralTessellation\\dataset\\" + imaname + ".png";

		src = imread(name, cv::IMREAD_GRAYSCALE);

		if (src.empty())
		{
			cout << "读取文件失败！" << std::endl;
			exit;
		}

		//将原图像转换为二值图像  
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
		//图像细化，骨骼化  
		Mat dst = thinImage(src);
		
		//过滤细化后的图像
		//filterOver(dst);
		
		//查找端点和交叉点  
		vector<Point2f> points = getPoints(dst, 6, 34, 12);

		//二值图转化成灰度图，并绘制找到的点
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

		//显示图像  
		cv::namedWindow("src1", CV_WINDOW_AUTOSIZE);
		cv::namedWindow("dst1", CV_WINDOW_AUTOSIZE);
		cv::imshow("src1", src);
		cv::imshow("dst1", dst);
		return points;
	}

}
//#include "tilingOpt.h"
//
//using namespace Tiling_tiles;
//int main(int argc, char** argv)
//{
//
//	//string imagename1 = "bird4"; 
//	//string imagename2 = "fish3";
//	//string imagename1 = "Boat";
//	//string imagename2 = "fish8";
//	string imagename1 = "fish5";
//	string imagename2 = "bird5";
//	//string txtname = "D:/images/111.png";
//	//string txtname1 = "D:/images/fish3.png";
//
//	//Tiling_tiles::Prototile *prototile_first;
//	//prototile_first = new Tiling_tiles::Prototile();
//	//////prototile_first->imgtocout(imagename1);
//	//prototile_first->loadTileData(imagename1);
//	//Tiling_tiles::Prototile *prototile_second;
//	//prototile_second = new Tiling_tiles::Prototile();
//	//prototile_first->loadTileData(imagename2);
//
//	Tiling_tiles::Tiling_opt *tiling;
//	tiling = new Tiling_tiles::Tiling_opt();
//	////tiling->com_cur_string(imagename1, imagename2);
//	////tiling->com_score(imagename1, imagename1);
//	tiling->com_score_manual(imagename1, imagename2);
//
//	//检测各种仿射变换
//	/*vector<Point2f> a;
//	vector<Point2f> b;
//	a.push_back(Point2f(1, 1));
//	a.push_back(Point2f(2, 1));
//	a.push_back(Point2f(2, 2));
//	a.push_back(Point2f(3, 2));
//	a.push_back(Point2f(3, 1));
//	a.push_back(Point2f(4, 1));
//	a.push_back(Point2f(5, 0));
//	a.push_back(Point2f(6, 1));
//	a.push_back(Point2f(7, 1));
//
//	b.push_back(Point2f(1, 1));
//	b.push_back(Point2f(2, 1));
//	b.push_back(Point2f(3, 0));
//	b.push_back(Point2f(4, 1));
//	b.push_back(Point2f(5, 1));
//	b.push_back(Point2f(5, 2));
//	b.push_back(Point2f(6, 2));
//	b.push_back(Point2f(6, 1));
//	b.push_back(Point2f(7, 1));
//
//	vector<vector<Point2f>> prototwo;
//	tiling->Aff_place(a, b, prototwo);*/
//	//Point2f s(5, 1);
//	//Point2f e(1, 5);
//	////double ab = tiling->re_warp_Aff(a, b, s, e);
//	//double ab = tiling->warpAff_sca(a, b, s, e);
//	//for (int i = 0; i < a.size(); i++)
//	//{
//	//	cout << "output: " << b[i] << endl;
//
//	//}
//
//	//warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)
//
//	//tiling->DTW(a, b);
//
//	
//	////简单的检测
//	//cv::Mat ima = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname);
//	//cv::Mat imb = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname1);
//	//cv::Mat imc = Mat::zeros(imb.rows,imb.cols,imb.type());
//	//vector<cv::Point2f> ima_p;
//	//vector<cv::Point2f> imb_p;
//	//vector<cv::Point2f> imc_p;
//
//	//MyLine(ima, Point2f(100, 100), Point2f(150, 50),"red");
//	//MyLine(ima, Point2f(150, 50), Point2f(200, 50), "red");
//	//MyLine(ima, Point2f(200, 50), Point2f(250, 100), "red");
//
//
//	//MyLine(imb, Point2f(100, 200), Point2f(150, 250), "blue");
//	//MyLine(imb, Point2f(150, 250), Point2f(200, 250), "blue");
//	//MyLine(imb, Point2f(200, 250), Point2f(250, 200), "blue");
//	//cv::Point2f a = cv::Point2f(100,100);
//	//cv::Point2f b = cv::Point2f(150, 100);
//	//cv::Point2f c = cv::Point2f(200, 100);
//	//cv::Point2f d = cv::Point2f(250, 100);
//	//cv::Point2f aa = cv::Point2f(100, 200);
//	//cv::Point2f bb = cv::Point2f(150, 200);
//	//cv::Point2f cc = cv::Point2f(200, 200);
//	//cv::Point2f dd = cv::Point2f(250, 200);
//	////cv::Point2f e = cv::Point2f(130, 150);
//	//ima_p.push_back(a);
//	//ima_p.push_back(b);
//	//ima_p.push_back(c);
//	//ima_p.push_back(d);
//	////ima_p.push_back(e);
//
//	//imb_p.push_back(aa);
//	//imb_p.push_back(bb);
//	//imb_p.push_back(cc);
//	//imb_p.push_back(dd);
//	////imb_p.push_back(c);
//	////imshow("hahahah",ima);
//	//ImageMorphing(ima, ima_p, imb, imb_p, imc, imc_p,0.5,0.5);
//	//imshow("ima",ima);
//	//imshow("imb",imb);
//	//imshow("imc",imc);
//
//
//	waitKey(0);
//	getchar();
//	return 0;
//}
//

#include <opencv2/opencv.hpp>  
#include <opencv2/core/core.hpp>  
#include <iostream>  
#include <vector>  
using namespace cv;
using namespace std;

/**
* @brief 对输入图像进行细化,骨骼化
* @param src为输入图像,用cvThreshold函数处理过的8位灰度图像格式，元素中只有0与1,1代表有元素，0代表为空白
* @param maxIterations限制迭代次数，如果不进行限制，默认为-1，代表不限制迭代次数，直到获得最终结果
* @return 为对src细化后的输出图像,格式与src格式相同，元素中只有0与1,1代表有元素，0代表为空白
*/
cv::Mat thinImage(const cv::Mat & src, const int maxIterations = -1)
{
	assert(src.type() == CV_8UC1);
	cv::Mat dst;
	int width = src.cols;
	int height = src.rows;
	src.copyTo(dst);
	int count = 0;  //记录迭代次数  
	while (true)
	{
		count++;
		if (maxIterations != -1 && count > maxIterations) //限制次数并且迭代次数到达  
			break;
		std::vector<uchar *> mFlag; //用于标记需要删除的点  
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
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
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

					if (ap == 1 && p2 * p4 * p8 == 0 && p2 * p6 * p8 == 0)
					{
						//标记  
						mFlag.push_back(p + j);
					}
				}
			}
		}

		//将标记的点删除  
		for (std::vector<uchar *>::iterator i = mFlag.begin(); i != mFlag.end(); ++i)
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
void filterOver(cv::Mat thinSrc)
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
std::vector<cv::Point> getPoints(const cv::Mat &thinSrc, unsigned int raudis = 4, unsigned int thresholdMax = 6, unsigned int thresholdMin = 4)
{
	assert(thinSrc.type() == CV_8UC1);
	int width = thinSrc.cols;
	int height = thinSrc.rows;
	cv::Mat tmp;
	thinSrc.copyTo(tmp);
	std::vector<cv::Point> points;
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
				Point point(j, i);
				points.push_back(point);
			}
		}
	}
	return points;
}


int main(int argc, char*argv[])
{
	cv::Mat src;
	//获取图像  
	if (argc != 2)
	{
		src = cv::imread("D:\\images\\222.png", cv::IMREAD_GRAYSCALE);
	}
	else
	{
		src = cv::imread(argv[1], cv::IMREAD_GRAYSCALE);
	}
	if (src.empty())
	{
		std::cout << "读取文件失败！" << std::endl;
		return -1;
	}
	
	//将原图像转换为二值图像  
	cv::threshold(src, src, 128, 1, cv::THRESH_BINARY);
	//cout << (int)src.at<uchar>(155, 155) << endl;
	/*for (int i = 0; i < src.cols; i++)
	{
		for (int j = 0; j < src.rows; j++)
		{
			if (src.at<uchar>(i, j) == 1)src.at<uchar>(i, j) = 0;
			else if (src.at<uchar>(i, j) == 0) src.at<uchar>(i, j) = 1;
		}
	}*/
	//图像细化，骨骼化  
	cv::Mat dst = thinImage(src);
	//过滤细化后的图像
	filterOver(dst);
	//查找端点和交叉点  
	std::vector<cv::Point> points = getPoints(dst, 6, 9, 6);
	//二值图转化成灰度图，并绘制找到的点
	dst = dst * 255;
	src = src * 255;
	vector<cv::Point>::iterator it = points.begin();
	for (; it != points.end(); it++)
	{
		circle(dst, *it, 4, 255, 1);
	}
	imwrite("dst.jpg", dst);
	//显示图像  
	cv::namedWindow("src1", CV_WINDOW_AUTOSIZE);
	cv::namedWindow("dst1", CV_WINDOW_AUTOSIZE);
	cv::imshow("src1", src);
	cv::imshow("dst1", dst);
	cv::waitKey(0);
}

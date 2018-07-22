#include "tilingOpt.h"

using namespace Tiling_tiles;
int main(int argc, char** argv)
{

	string imagename1 = "bird4"; 
	string imagename2 = "fish3";
	//string imagename1 = "Boat";
	//string imagename2 = "fish8";
	//string imagename1 = "fish5";
	//string imagename2 = "bird5";
	//string txtname = "D:/images/111.png";
	//string txtname1 = "D:/images/fish3.png";

	//vector<Point2f> ske_points;
	//ske_points = get_Skeleton(imagename1);


	Tiling_tiles::Prototile *prototile_first;
	prototile_first = new Tiling_tiles::Prototile();
	//////prototile_first->imgtocout(imagename1);
	prototile_first->loadTileData(imagename2);
	//Tiling_tiles::Prototile *prototile_second;
	//prototile_second = new Tiling_tiles::Prototile();
	//prototile_first->loadTileData(imagename2);

	//Tiling_tiles::Tiling_opt *tiling;
	//tiling = new Tiling_tiles::Tiling_opt();
	//////tiling->com_cur_string(imagename1, imagename2);
	//////tiling->com_score(imagename1, imagename1);
	//tiling->com_score_manual(imagename1, imagename2);

	//¼ì²â¸÷ÖÖ·ÂÉä±ä»»
	/*vector<Point2f> a;
	vector<Point2f> b;
	a.push_back(Point2f(1, 1));
	a.push_back(Point2f(2, 1));
	a.push_back(Point2f(2, 2));
	a.push_back(Point2f(3, 2));
	a.push_back(Point2f(3, 1));
	a.push_back(Point2f(4, 1));
	a.push_back(Point2f(5, 0));
	a.push_back(Point2f(6, 1));
	a.push_back(Point2f(7, 1));

	b.push_back(Point2f(1, 1));
	b.push_back(Point2f(2, 1));
	b.push_back(Point2f(3, 0));
	b.push_back(Point2f(4, 1));
	b.push_back(Point2f(5, 1));
	b.push_back(Point2f(5, 2));
	b.push_back(Point2f(6, 2));
	b.push_back(Point2f(6, 1));
	b.push_back(Point2f(7, 1));

	vector<vector<Point2f>> prototwo;
	tiling->Aff_place(a, b, prototwo);*/
	//Point2f s(5, 1);
	//Point2f e(1, 5);
	////double ab = tiling->re_warp_Aff(a, b, s, e);
	//double ab = tiling->warpAff_sca(a, b, s, e);
	//for (int i = 0; i < a.size(); i++)
	//{
	//	cout << "output: " << b[i] << endl;

	//}

	//warpAff_sca(vector<Point2f> &input_, vector<Point2f> &output_, Point2f start, Point2f end)

	//tiling->DTW(a, b);

	
	////¼òµ¥µÄ¼ì²â
	//cv::Mat ima = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname);
	//cv::Mat imb = Mat(300, 300, CV_8UC3, Scalar(255, 255, 255));//imread(txtname1);
	//cv::Mat imc = Mat::zeros(imb.rows,imb.cols,imb.type());
	//vector<cv::Point2f> ima_p;
	//vector<cv::Point2f> imb_p;
	//vector<cv::Point2f> imc_p;

	//MyLine(ima, Point2f(100, 100), Point2f(150, 50),"red");
	//MyLine(ima, Point2f(150, 50), Point2f(200, 50), "red");
	//MyLine(ima, Point2f(200, 50), Point2f(250, 100), "red");


	//MyLine(imb, Point2f(100, 200), Point2f(150, 250), "blue");
	//MyLine(imb, Point2f(150, 250), Point2f(200, 250), "blue");
	//MyLine(imb, Point2f(200, 250), Point2f(250, 200), "blue");
	//cv::Point2f a = cv::Point2f(100,100);
	//cv::Point2f b = cv::Point2f(150, 100);
	//cv::Point2f c = cv::Point2f(200, 100);
	//cv::Point2f d = cv::Point2f(250, 100);
	//cv::Point2f aa = cv::Point2f(100, 200);
	//cv::Point2f bb = cv::Point2f(150, 200);
	//cv::Point2f cc = cv::Point2f(200, 200);
	//cv::Point2f dd = cv::Point2f(250, 200);
	////cv::Point2f e = cv::Point2f(130, 150);
	//ima_p.push_back(a);
	//ima_p.push_back(b);
	//ima_p.push_back(c);
	//ima_p.push_back(d);
	////ima_p.push_back(e);

	//imb_p.push_back(aa);
	//imb_p.push_back(bb);
	//imb_p.push_back(cc);
	//imb_p.push_back(dd);
	////imb_p.push_back(c);
	////imshow("hahahah",ima);
	//ImageMorphing(ima, ima_p, imb, imb_p, imc, imc_p,0.5,0.5);
	//imshow("ima",ima);
	//imshow("imb",imb);
	//imshow("imc",imc);


	waitKey(0);
	getchar();
	return 0;
}


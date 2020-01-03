#include "cgal_tool.h"

Polygon_set_2 A_difference_B(Polygon_set_2 A, Polygon_2 B)
{
	A.difference(B);
	return A;
	// Construct the two initial polygons and the clipping rectangle.
	Polygon_2 ttt;
	ttt.push_back(Point_2(-2, 0));
	ttt.push_back(Point_2(2, 0));
	ttt.push_back(Point_2(-2, 4));

	Polygon_2 rect;
	rect.push_back(Point_2(0, 0));
	rect.push_back(Point_2(3, 0));
	rect.push_back(Point_2(3, 2));
	rect.push_back(Point_2(0, 2));

	Polygon_2 rect2;
	rect2.push_back(Point_2(-3, 1));
	rect2.push_back(Point_2(3, 1));
	rect2.push_back(Point_2(3, 2));
	rect2.push_back(Point_2(-2, 2));

	// Perform a sequence of operations.
	Polygon_set_2 S;
	S.insert(ttt);
	//S.join (Q);                   // Compute the union of P and Q.
	//S.complement();               // Compute the complement.
	//S.intersection (rect);        // Intersect with the clipping rectangle.
	S.difference(rect2);
	// Print the result.
	//std::list<Polygon_with_holes_2> res;
	//std::list<Polygon_with_holes_2>::const_iterator it;

	std::cout << "The result contains " << S.number_of_polygons_with_holes()
		<< " components:" << std::endl;

	/*S.polygons_with_holes(std::back_inserter(res));
	for (it = res.begin(); it != res.end(); ++it) {
		std::cout << "--> ";
		print_polygon_with_holes(*it);
	}*/

	//return S;
}


vector<Point2f> Polygon_2vector(Polygon_2 poly)
{
	vector<Point2f> contour;
	int p_size = poly.size();
	for (int i = 0; i < p_size; i++)
	{
		contour.push_back(To_cvp(poly[i]));
	}

	return contour;
}

vector<Point2f> Polygon2vector(Polygon2 poly)
{
	vector<Point2f> contour;
	int p_size = poly.size();
	for (int i = 0; i < p_size; i++)
	{
		contour.push_back(To_cvp(poly[i]));
	}
	return contour;
}

Polygon_2 vectorPolygon_2(vector<Point2f> contour)
{
	Polygon_2 poly;
	int c_size = contour.size();
	for (int i = 0; i < c_size; i++)
	{
		poly.push_back(Point_2(contour[i].x, contour[i].y));
	}
	return poly;
}

Polygon2 vectorPolygon2(vector<Point2f> contour)
{
	Polygon2 poly;
	int c_size = contour.size();
	for (int i = 0; i < c_size; i++)
	{
		poly.push_back(Point2(contour[i].x, contour[i].y));
	}
	return poly;
}

Polygon2 offset_poly(double Offset, Polygon2 poly) //ÕýÖµÎª¸¯Ê´
{
	SsPtr ss = CGAL::create_interior_straight_skeleton_2(poly);
	PolygonPtrVector offset_polygons = CGAL::create_offset_polygons_2<Polygon2>(Offset, *ss);
	
	std::cout << "Polygon list with " << offset_polygons.size() << " polygons" << std::endl;
	
	PolygonPtrVector::const_iterator pi = offset_polygons.begin();
	Polygon2 outp = **pi;
	for (; pi != offset_polygons.end(); ++pi)
	{
		Polygon2 pt= **pi;
		std::cout << "Polygon with " << pt.size() << " vertices" << std::endl;
		for (Polygon2::Vertex_const_iterator vi = pt.vertices_begin(); vi != pt.vertices_end(); ++vi)
		{
			std::cout << "(" << (*vi).x() << "," << (*vi).y() << ")" << std::endl;
		}
	}
	return outp;
}


Point2f To_cvp(Point2 p)
{
	std::ostringstream osx,osy;
	osx << p.x();
	osy << p.y();
	return Point2f(stod(osx.str()), stod(osy.str()));
}

Point2f To_cvp(Point_2 p)
{
	std::ostringstream osx, osy;
	osx << p.x();
	osy << p.y();
	return Point2f(stod(osx.str()), stod(osy.str()));
}
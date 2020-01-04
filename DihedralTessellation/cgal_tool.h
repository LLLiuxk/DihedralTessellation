#ifndef cgal_tool_h
#define cgal_tool_h

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Triangulation_2.h>
//boolean
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


//skeleton and offset
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<boost/shared_ptr.hpp>
#include<CGAL/create_offset_polygons_2.h>

#include <list>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Polygon_set_2<Kernel>                       Polygon_set_2;


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                   Point2;
typedef CGAL::Polygon_2<K>           Polygon2;
typedef CGAL::Straight_skeleton_2<K> Ss;
typedef boost::shared_ptr<Polygon2> PolygonPtr;
typedef boost::shared_ptr<Ss> SsPtr;
typedef std::vector<PolygonPtr> PolygonPtrVector;

//typedef CGAL::Triangulation_2<K>      Triangulation;
//typedef Triangulation::Vertex_handle  Vertex_handle;
//typedef Triangulation::Point          Point;
//typedef Triangulation::Finite_vertex_handles    Finite_vertex_handles;
// The following types are different
// Its value type is Triangulation_2::Vertex
//typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
// Its value type is Triangulation_2::Vertex_handle
//typedef Finite_vertex_handles::iterator         Finite_vertex_handles_iterator;

#include "print_utils.h"


Polygon_set_2  A_difference_B(Polygon_set_2 A, Polygon_2 B);

Polygon2 offset_poly(double Offset, Polygon2 poly);//ÕýÖµÎª¸¯Ê´

#endif
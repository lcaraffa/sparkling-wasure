/*
 * This file is part of the watertight surface reconstruction code https://github.com/lcaraffa/spark-ddt
 * Copyright (c) 2024 Caraffa Laurent, Mathieu Brédif.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef DDT_CGAL_TRAITS_3_HPP
#define DDT_CGAL_TRAITS_3_HPP


#include <limits>
#include <functional>
#include <math.h>
#include <unordered_set>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include "cgal_traits_base.hpp"
#include "Facet_const_iterator_3.hpp"


#include "../io/number_parser.hpp"
namespace ddt
{


struct Cgal_traits_3_Raw
{
    enum { D = 3 };
    typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
    typedef CGAL::Delaunay_triangulation_3<K>                 Delaunay_triangulation;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_handle_raw;
    typedef typename Delaunay_triangulation::Cell_handle                            Cell_handle_raw;

    typedef typename Delaunay_triangulation::Vertex_iterator                    Vertex_const_iterator;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_const_handle;
    typedef typename Delaunay_triangulation::Vertex_iterator                          Vertex_iterator;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_handle;

    typedef typename Delaunay_triangulation::Cell_iterator                 Cell_const_iterator;
    typedef typename Delaunay_triangulation::Cell_handle                   Cell_const_handle;
    typedef typename Delaunay_triangulation::Cell_iterator                       Cell_iterator;
    typedef typename Delaunay_triangulation::Cell_handle                         Cell_handle;

    typedef typename Delaunay_triangulation::Facet_iterator                 Facet_iterator;
    typedef Facet_iterator                   Facet_const_handle;
    typedef  K::Vector_3 Vector;

    typedef typename Delaunay_triangulation::Point                                    Point;
    typedef std::pair<Cell_const_handle, int>                      Facet;

    typedef typename CGAL::Unique_hash_map<Vertex_handle_raw, uint> v_hmap_uint;


    Delaunay_triangulation triangulation(int dimension) const
    {
        return Delaunay_triangulation();
    }

    inline size_t number_of_cells(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_cells();
    }
    inline size_t number_of_vertices(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_vertices();
    }

    inline Cell_const_iterator cells_begin(const Delaunay_triangulation& dt) const
    {
        return dt.cells_begin();
    }
    inline Cell_const_iterator cells_end(const Delaunay_triangulation& dt) const
    {
        return dt.cells_end();
    }

    inline Facet_iterator facets_begin(Delaunay_triangulation& dt)
    {
        return dt.facets_begin();
    }
    inline Facet_iterator facets_end( Delaunay_triangulation& dt)
    {
        return dt.facets_end();
    }

    inline Vertex_const_iterator vertices_begin(const Delaunay_triangulation& dt) const
    {
        return dt.vertices_begin();
    }
    inline Vertex_const_iterator vertices_end(const Delaunay_triangulation& dt) const
    {
        return dt.vertices_end();
    }

    inline Vertex_iterator vertices_begin(Delaunay_triangulation& dt) const
    {
        return dt.vertices_begin();
    }
    inline Vertex_iterator vertices_end(Delaunay_triangulation& dt) const
    {
        return dt.vertices_end();
    }

    inline int index_of_covertex(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        return f->second;
    }

    inline Cell_const_handle full_cell(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        return f->first;
    }

    inline Cell_const_iterator neighbor(const Delaunay_triangulation& dt, Cell_const_iterator c, int i) const
    {
        return c->neighbor(i);
    }

    inline Point circumcenter(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        return  CGAL::circumcenter(c->vertex(0)->point(),
                                   c->vertex(1)->point(),
                                   c->vertex(2)->point(),
                                   c->vertex(3)->point());
    }


    template<typename PP>
    double squared_dist(const PP & p1, const PP & p2, int dim) const
    {
        double acc = 0;
        for(int d = 0; d < dim; d++)
            acc += (p1[d] - p2[d])*(p1[d] - p2[d]);
        return acc;
    }


    template<int D>
    bool is_inside(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c) const
    {
        typename Delaunay_triangulation::Geom_traits::FT res = 0, delta;
        auto point  = c->vertex(0)->point();
        auto center = circumcenter(dt, c);
        double dist = sqrt(squared_dist(point,center,D));
        for(int d = 0; d < D; d++)
        {
            if(dist > fabs(center[d] - bbox.max(d)) || dist > fabs(center[d] - bbox.min(d)))
                return false;
        }
        return true;
    }




};



template<typename DV, typename DC,typename G = bool>
struct Cgal_traits_3
{
    enum { D = 3 };

    typedef DV                                                     Data_V;
    typedef DC                                                     Data_C;
    typedef typename DV::Id                                             Id;
    typedef typename DV::Flag                                        Flag_V;
    typedef typename DC::Flag                                        Flag_C;


    typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
    typedef CGAL::Triangulation_vertex_base_with_info_3<Data_V, K>   Vb;

    typedef CGAL::Triangulation_cell_base_with_info_3<Data_C,K>            Cb;
    typedef CGAL::Triangulation_data_structure_3<Vb,Cb>               TDS;
    //#endif
    typedef typename K::Point_3                                    Point;
    typedef typename K::Sphere_3                                    Sphere;
    typedef typename K::Plane_3                                    Plane;
    typedef std::pair<Point,Id>                                    Point_id;
    typedef std::tuple<Point,Id,Id>                                Point_id_id;
    typedef  K::Vector_3 Vector;

    typedef typename TDS::Vertex_iterator                          Vertex_const_iterator;
    typedef typename TDS::Vertex_handle                            Vertex_const_handle;
    typedef typename TDS::Vertex_iterator                          Vertex_iterator;
    typedef typename TDS::Vertex_handle                            Vertex_handle;

    typedef typename TDS::Cell_iterator                            Cell_const_iterator;
    typedef typename TDS::Cell_handle                              Cell_const_handle;
    typedef typename TDS::Cell_iterator                            Cell_iterator;
    typedef typename TDS::Cell_handle                              Cell_handle;

    typedef std::pair<Cell_const_handle, int>                      Facet;
    typedef Facet_const_iterator_3<TDS>                            Facet_const_iterator;
    typedef Facet_const_iterator                                   Facet_const_handle;
    typedef Facet_const_iterator                                   Facet_iterator;
    typedef Facet_const_iterator                                   Facet_handle;


    typedef CGAL::Delaunay_triangulation_3<K, TDS>                 Delaunay_triangulation;
    typedef CGAL::Random_points_in_sphere_3<Point>                 Random_points_in_ball;

    typedef  CGAL::Unique_hash_map<Vertex_const_handle,bool> v_hmap_bool;
    typedef  CGAL::Unique_hash_map<Vertex_handle,uint> v_hmap_uint;



    template<typename Iterator> static Point make_point(Iterator it)
    {
        Iterator x = it;
        Iterator y = ++it;
        Iterator z = ++it;
        return Point(*x,*y,*z);
    }


    Delaunay_triangulation triangulation(int dimension) const
    {
        assert(dimension == D);
        return Delaunay_triangulation();
    }

    inline Id    id  (Vertex_const_handle v) const
    {
        return v->info().id;
    }
    inline Flag_V& flag(Vertex_const_handle v) const
    {
        return v->info().flag;
    }

    inline Id    gid  (Vertex_const_handle v) const
    {
        return v->info().gid;
    }

    inline Id    gid  (Cell_const_handle v) const
    {
        return v->info().gid;
    }


    inline Flag_C& flag(Cell_const_handle c) const
    {
        return c->info().flag;
    }

    inline Data_C& data(Cell_const_handle c)  const
    {
        return c->info();
    }


    inline Data_V& data(Vertex_const_handle v)  const
    {
        return v->info();
    }

    inline int current_dimension(const Delaunay_triangulation& dt) const
    {
        return dt.dimension();
    }
    inline int maximal_dimension(const Delaunay_triangulation& dt) const
    {
        return D;
    }
    inline size_t number_of_cells(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_cells();
    }
    inline size_t number_of_vertices(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_vertices();
    }
    inline Vertex_const_handle vertex(const Delaunay_triangulation& dt, Cell_const_handle c, int i) const
    {
        return c->vertex(i);
    }
    inline Vertex_const_iterator vertices_begin(const Delaunay_triangulation& dt) const
    {
        return dt.all_vertices_begin();
    }
    inline Vertex_const_iterator vertices_end(const Delaunay_triangulation& dt) const
    {
        return dt.all_vertices_end();
    }
    inline Vertex_iterator vertices_begin(Delaunay_triangulation& dt) const
    {
        return dt.all_vertices_begin();
    }
    inline Vertex_iterator vertices_end(Delaunay_triangulation& dt) const
    {
        return dt.all_vertices_end();
    }
    inline Facet_const_iterator facets_begin(const Delaunay_triangulation& dt) const
    {
        return Facet_const_iterator(dt.tds());
    }
    inline Facet_const_iterator facets_end(const Delaunay_triangulation& dt) const
    {
        return Facet_const_iterator();
    }
    inline Cell_const_iterator cells_begin(const Delaunay_triangulation& dt) const
    {
        return dt.all_cells_begin();
    }
    inline Cell_const_iterator cells_end(const Delaunay_triangulation& dt) const
    {
        return dt.all_cells_end();
    }

    inline Vertex_handle infinite_vertex(const Delaunay_triangulation& dt) const
    {
        return dt.infinite_vertex();
    }

    inline void clear(Delaunay_triangulation& dt) const
    {
        return dt.clear();
    }


    void get_list_vertices(const Cell_const_handle fch,std::list<Vertex_const_handle> & lp) const
    {
        for(int dd = 0; dd < D+1; dd++)
        {
            lp.push_back(fch->vertex(dd));
        }
    }


    std::vector<double> get_cell_barycenter(Cell_const_handle ch) const
    {
        std::vector<double>  coords(D,0);
        for(int dd = 0; dd < D+1; dd++)
        {
            for(uint d = 0; d < D; d++)
            {
                coords[d] += (ch->vertex(dd)->point())[d];
            }
        }
        for(uint d = 0; d < D; d++)
            coords[d] /= ((double)D+1);
        return coords;
    }


    std::vector<double> get_cell_barycenter_const(Cell_const_handle ch) const
    {
        std::vector<double>  coords(D,0);
        for(int dd = 0; dd < D+1; dd++)
        {
            for(uint d = 0; d < D; d++)
            {
                coords[d] += (ch->vertex(dd)->point())[d];
            }
        }
        for(uint d = 0; d < D; d++)
            coords[d] /= ((double)D+1);
        return coords;
    }


    Cell_const_handle locate_cell_point(const Delaunay_triangulation& dt,Point & pp) const
    {
        Cell_handle cc = dt.locate(pp);
        return cc;
    }


    Cell_const_handle locate_cell(const Delaunay_triangulation& dt,Cell_const_handle c) const
    {
        std::vector<double> coords = get_cell_barycenter_const(c);
        auto pp = Point(coords[0],coords[1],coords[2]);
        return locate_cell_point(dt,pp);
    }

    template<class It> Vertex_handle insert_splay_one(Delaunay_triangulation& dt, Id id, It it, Vertex_handle hint,
            std::map<Id, std::vector<Vertex_const_handle>>& out) const
    {
        Vertex_handle v = hint;
        Id pid = it->second;
        std::unordered_set<Vertex_handle> vertices;
        std::unordered_set<int> ids;
        if (dt.dimension() < 3)
        {
            v = dt.insert(it->first);
            dt.finite_adjacent_vertices(v, std::inserter(vertices, vertices.end()));
            for(auto vit = vertices.begin(); vit != vertices.end(); ++vit)
                ids.insert((*vit)->info().id);
        }
        else
        {
            typename Delaunay_triangulation::Locate_type lt;
            int li, lj;
            Cell_handle c = dt.locate(it->first, lt, li, lj, hint);
            if (lt == Delaunay_triangulation::VERTEX)
                return v;
            std::vector<Cell_handle> cells;
            cells.reserve(32);
            std::vector<Facet> facets;
            facets.reserve(32);
            dt.find_conflicts(it->first, c,
                              std::back_inserter(facets), // Conflict facets
                              std::back_inserter(cells)); // Conflict cells
            for(auto fit = facets.begin(); fit != facets.end(); ++fit)
            {
                bool foreign = true;
                for(int i=0; i<=D; ++i)
                {
                    if (i==fit->second) continue;
                    Vertex_handle vi = fit->first->vertex(i);
                    if (dt.is_infinite(vi)) continue;
                    Id vid = vi->info().id;
                    ids.insert(vid);
                    foreign = foreign && vid != id;
                }
                if (pid != id && !foreign)
                    for(int i=0; i<=D; ++i)
                    {
                        if (i==fit->second) continue;
                        Vertex_handle vi = fit->first->vertex(i);
                        if (dt.is_infinite(vi)) continue;
                        vertices.insert(vi);
                    }
            }
            if (pid != id && ids.find(id) == ids.end())
                return v; // future star of "it" is entirely foreign, skipping
            v = dt.insert_in_hole(it->first, cells.begin(), cells.end(), facets.begin()->first, facets.begin()->second);
        }
        for(auto vit = vertices.begin(); vit != vertices.end(); ++vit)
            if(!dt.is_infinite(*vit) && pid != (*vit)->info().id)
                out[pid].push_back(*vit);
        // send it to all ids (if != id)
        for(auto iit = ids.begin(); iit != ids.end(); ++iit)
            if(*iit != id)
                out[*iit].push_back(v);
        v->info().id = pid;
        return v;
    }

    template<class It> int insert_splay(Delaunay_triangulation& dt, Id id, It begin, It end, std::map<Id, std::vector<Vertex_const_handle>>& out) const
    {
        int res = 0;
        Vertex_handle hint;
        for(It it = begin; it != end ; ++it)
        {
            Vertex_handle v = insert_splay_one(dt, id, it, hint, out);
            if (v != Vertex_handle())
            {
                hint = v;
                ++res;
            }
        }
        return res;
    }

    Vertex_handle insert_simplified_one(Delaunay_triangulation& dt, Id id, Id pid, const Point& p, Vertex_handle hint) const
    {
        Vertex_handle v;
        if (id == pid || dt.dimension() < D)
        {
            v = dt.insert(p);
        }
        else
        {
            // Locate the point
            typename Delaunay_triangulation::Locate_type lt;
            int li, lj;
            Cell_handle c = dt.locate(p, lt, li, lj, hint);
            if (lt == Delaunay_triangulation::VERTEX)
            {
                return v; // Point already exists
            }
            // Get the cells that conflict with p in a vector V,
            // and a facet on the boundary of this hole in f.
            std::vector<Cell_handle> cells;
            cells.reserve(32);
            std::vector<Facet> facets;
            facets.reserve(32);
            dt.find_conflicts(p, c,
                              std::back_inserter(facets), // Conflict facets
                              std::back_inserter(cells)); // Conflict cells
            bool foreign = true;
            for(auto fit = facets.begin(); fit != facets.end(); ++fit)
            {
                for(int i=0; i<=D; ++i)
                {
                    if (i==fit->second) continue;
                    Vertex_handle vi = fit->first->vertex(i);
                    if (!dt.is_infinite(vi) && id == vi->info().id) { foreign = false; break; }
                }
            }
            if (foreign)
            {
                return v; // future star of "it" is entirely foreign, skipping
            }
            v = dt.insert_in_hole(p, cells.begin(), cells.end(), facets.begin()->first, facets.begin()->second);
        }
        v->info().id = pid;
        return v;
    }


    template<class It> inline void insert(Delaunay_triangulation& dt, It begin, It end) const
    {
        dt.insert(begin, end);
    }

    template<class It> inline void remove(Delaunay_triangulation& dt, It begin, It end) const
    {
        dt.remove_cluster(begin, end);
    }

    template<class It> inline It get_local_convex_hull(const Delaunay_triangulation& dt, Id id, It out) const
    {
        std::vector<Vertex_const_handle> hull;
        dt.adjacent_vertices(dt.infinite_vertex(), std::back_inserter(hull)); // get convex hull points
        for(auto vh : hull) if( vh->info().id == id) *out++ = vh; // filter local points
        return out;
    }

    inline bool is_visible(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c, int i) const
    {
        const Point& p0 = c->vertex(i            )->point();
        const Point& p1 = c->vertex((i+1) % (D+1))->point();
        const Point& p2 = c->vertex((i+2) % (D+1))->point();
        int res = (i%2==0) ? CGAL::ON_NEGATIVE_SIDE : CGAL::ON_POSITIVE_SIDE;
        return
            dt.orientation(p0, p1, p2, Point(bbox.min(0), bbox.min(1), bbox.min(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.max(0), bbox.min(1), bbox.min(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.min(0), bbox.max(1), bbox.min(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.max(0), bbox.max(1), bbox.min(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.min(0), bbox.min(1), bbox.max(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.max(0), bbox.min(1), bbox.max(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.min(0), bbox.max(1), bbox.max(2))) == res &&
            dt.orientation(p0, p1, p2, Point(bbox.max(0), bbox.max(1), bbox.max(2))) == res ;
    }

    bool is_touching(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c) const
    {
        typename Delaunay_triangulation::Geom_traits::FT res = 0, delta;
        auto point  = c->vertex(0)->point();
        auto center = circumcenter(dt, c);
        for(int d = 0; d < D; d++)
        {
            res -= CGAL::square(center[d] - point[d]);
            delta = center[d] - bbox.max(d);
            if(CGAL::is_positive(delta))
                res += CGAL::square(delta);
            else
            {
                delta = center[d] - bbox.min(d);
                if(CGAL::is_negative(delta))
                    res += CGAL::square(delta);
            }
        }
        return CGAL::is_negative(res);
    }

    inline Point circumcenter(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        return  CGAL::circumcenter(c->vertex(0)->point(),
                                   c->vertex(1)->point(),
                                   c->vertex(2)->point(),
                                   c->vertex(3)->point());
    }

    inline bool vertex_is_infinite(const Delaunay_triangulation& dt, Vertex_const_handle v) const
    {
        return dt.is_infinite(v);
    }

    inline bool facet_is_infinite(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        for(int i = 0; i<=D; ++i)
            if(i!=f->second && dt.is_infinite(f->first->vertex(i)))
                return true;
        return false;
    }

    inline bool cell_is_infinite(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        for(int i = 0; i<=D; ++i)
            if(dt.is_infinite(c->vertex(i)))
                return true;
        return false;
    }

    inline const Point& point(const Delaunay_triangulation& dt, Vertex_const_handle v) const
    {
        return v->point();
    }

    inline double coord(const Delaunay_triangulation& dt, const Point& p, int i) const
    {
        return CGAL::to_double(p[i]);
    }

    inline double scalar_product(const Point& p, const Point& q) const
    {
        return CGAL::to_double(p.x()*q.x()+p.y()*q.y()+p.z()*q.z());
    }

    template<typename V>
    bool are_vertices_equal(const Delaunay_triangulation& t1, V v1, const Delaunay_triangulation& t2, V v2) const
    {
        bool inf1 = vertex_is_infinite(t1, v1);
        bool inf2 = vertex_is_infinite(t2, v2);
        return (inf1 || inf2) ? (inf1 == inf2) : v1->point() == v2->point();
    }

    template<typename C>
    bool are_cells_equal(const Delaunay_triangulation& t1, C c1, const Delaunay_triangulation& t2, C c2) const
    {
        for(int i1=0; i1<=D; ++i1)
        {
            Vertex_handle v1 = c1->vertex(i1);
            bool is_equal = false;
            for(int i2=0; i2<=D; ++i2)
            {
                Vertex_handle v2 = c2->vertex(i2);
                if(are_vertices_equal(t1, v1, t2, v2))
                {
                    is_equal = true;
                    break;
                }
            }
            if(!is_equal)
                return false;
        }
        return true;
    }

    inline int index_of_covertex(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        return f->second;
    }

    inline Cell_const_handle full_cell(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        return f->first;
    }

    inline int mirror_index(const Delaunay_triangulation& dt, Cell_const_handle c, int i) const
    {
        return dt.mirror_index(c, i);
    }

    inline Cell_const_iterator neighbor(const Delaunay_triangulation& dt, Cell_const_iterator c, int i) const
    {
        return c->neighbor(i);
    }

    Facet_const_iterator facet(const Delaunay_triangulation& dt, Cell_const_iterator c, int i) const
    {
        Facet f(c, i);
        return Facet_const_iterator(dt.tds(), f);
    }

    std::istream& read_cgal(std::istream&  ifile, Delaunay_triangulation& dt, bool do_data = true, bool is_ascii = true)
    {
        return ifile;
    }

    std::ostream & write_cgal(std::ostream & ofile, const Delaunay_triangulation& dt, bool do_data = true, bool is_ascii = true) const
    {
        return ofile;
    }


    void incident_cells(const Delaunay_triangulation& dt, Vertex_const_handle v, std::vector<Cell_const_handle> & cells) const
    {
        cells.reserve(32);
        dt.incident_cells(v, std::back_inserter(cells));
    }

    template<typename Unary_op>
    void finite_adjacent_vertices(const Delaunay_triangulation& dt, Vertex_const_handle v, Unary_op op) const
    {
        std::vector<Vertex_const_handle> vertices;
        vertices.reserve(32);
        dt.finite_adjacent_vertices(v, std::back_inserter(vertices));
        for(auto c : vertices) op(c);
    }


    template<typename PP>
    double squared_dist(const PP & p1, const PP & p2, int dim) const
    {
        double acc = 0;
        for(int d = 0; d < dim; d++)
            acc += (p1[d] - p2[d])*(p1[d] - p2[d]);
        return acc;
    }


    template<int D>
    bool is_inside(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c) const
    {
        typename Delaunay_triangulation::Geom_traits::FT res = 0, delta;
        auto point  = c->vertex(0)->point();
        auto center = circumcenter(dt, c);
        double dist = sqrt(squared_dist(point,center,D));
        for(int d = 0; d < D; d++)
        {
            if(dist > fabs(center[d] - bbox.max(d)) || dist > fabs(center[d] - bbox.min(d)))
                return false;
        }
        return true;
    }



    template<typename DDT_DATA>
    void export_tri_to_data(Delaunay_triangulation& tri,DDT_DATA & data)
    {
        typedef typename DDT_DATA::Data_ply Data_ply;
        int D = data.D;
        data.dmap[data.xyz_name] = Data_ply(data.xyz_name,"vertex",D,D,data.get_float_type());
        data.dmap[data.simplex_name] = Data_ply(data.simplex_name,"face",D+1,D+1,data.get_int_type());
        data.dmap[data.nb_name] = Data_ply(data.nb_name,"face",D+1,D+1,data.get_int_type());
        data.dmap[data.vid_name] = Data_ply(data.vid_name,"vertex",1,1,data.get_int_type());
        data.dmap[data.flag_vertex_name] = Data_ply(data.flag_vertex_name,"vertex",1,1,data.get_int_type());
        data.dmap[data.flag_simplex_name] = Data_ply(data.flag_simplex_name,"face",1,1,data.get_int_type());
        std::vector<double> v_xyz;
        std::vector<Id> v_simplex,v_nb;
        std::vector<Id> v_vid;
        //    std::vector<Id> v_cid;
        std::vector<Flag_V> v_flagv;
        std::vector<Flag_C> v_flags;
        uint n = tri.number_of_vertices();
        std::map<Vertex_const_handle, uint> vertex_map;
        vertex_map[tri.infinite_vertex()] = 0;
        v_vid.push_back(0);
        v_flagv.push_back(0);
        for(int d = 0; d < D; d++)
        {
            v_xyz.push_back(0.0);
        }
        uint i = 1;
        for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
        {
            if(tri.is_infinite(vit))
            {
                continue;
            }
            int ii = i;
            for(int d = 0; d < D; d++)
            {
                double pv = vit->point()[d];
                v_xyz.push_back(pv);
            }
            v_vid.push_back(vit->info().id);
            v_flagv.push_back(vit->info().flag);
            vertex_map[vit] = ii;
            i++;
        }
        // write the number of cells
        n = tri.number_of_cells();
        v_simplex.resize(n*(D+1));
        v_flags.resize(n);
        v_nb.resize(n*(D+1));
        // write the cells
        std::map<Cell_const_handle, uint> cell_map;
        i = 0;
        int max_id = 0;
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            if(it->info().gid > max_id)
                max_id = it->info().gid;
        }
        max_id++;
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            int ii =  it->info().gid;
            if(ii == -1)
                ii = max_id++;
            cell_map[it] = ii;
            for(int d = 0; d < D+1; d++)
            {
                int vertex_id = vertex_map[it->vertex(d)] ;
                v_simplex[ii*(D+1)+d] = vertex_id;
                // write the id
            }
            for(int j = 0; j < D+1; j++)
            {
                int nb_id = cell_map[it->neighbor(j)];
                v_nb[ii*(D+1) +j] = nb_id;
            }
            v_flags[ii] = it->info().flag;
            ++i;
        }
        data.dmap[data.xyz_name].fill_full_uint8_vect(v_xyz);
        data.dmap[data.simplex_name].fill_full_uint8_vect(v_simplex);
        data.dmap[data.nb_name].fill_full_uint8_vect(v_nb);
        data.dmap[data.vid_name].fill_full_uint8_vect(v_vid);
        data.dmap[data.flag_vertex_name].fill_full_uint8_vect(v_flagv);
        data.dmap[data.flag_simplex_name].fill_full_uint8_vect(v_flags);
    }

    template<typename DDT_DATA>
    void build_tri_from_data(Delaunay_triangulation& tri,DDT_DATA & data,int tid)
    {
        std::vector<double> v_xyz;
        std::vector<int> v_simplex,v_nb,v_vid,v_cid,v_flagv,v_flags;
        data.dmap[data.xyz_name].extract_full_input(v_xyz,false);
        data.dmap[data.simplex_name].extract_full_input(v_simplex,false);
        data.dmap[data.nb_name].extract_full_input(v_nb,false);
        bool do_data = false;
        if(data.dmap.find(data.flag_vertex_name) != data.dmap.end() &&
                data.dmap.find(data.flag_simplex_name) != data.dmap.end() &&
                data.dmap.find(data.vid_name) != data.dmap.end() &&
                data.dmap.find(data.cid_name) != data.dmap.end() )
        {
            do_data = true;
            data.dmap[data.flag_vertex_name].extract_full_input(v_flagv,false);
            data.dmap[data.flag_simplex_name].extract_full_input(v_flags,false);
            data.dmap[data.vid_name].extract_full_input(v_vid,false);
            data.dmap[data.cid_name].extract_full_input(v_cid,false);
        }
        int D = 3;
        tri.clear();
        tri.tds().set_dimension(D);
        auto cit = tri.all_cells_begin();
        tri.tds().delete_cell(cit);
        // 4) read the number of vertices
        uint num_v = data.dmap[data.xyz_name].get_nbe_input();
        std::vector<Vertex_handle> vertex_map(num_v);
        vertex_map[0] = tri.infinite_vertex();
        // 5) read and create the vertices
        for(uint i = 1; i < num_v; ++i)
        {
            int ii = i;
            std::vector<double> coords_v(D);
            for(uint d = 0; d < D; d++)
            {
                coords_v[d] = v_xyz[ii*D +d];
            }
            Point p(coords_v[0],coords_v[1],coords_v[2]);
            vertex_map[ii] = tri.tds().create_vertex(p);
            if(do_data)
            {
                vertex_map[ii]->info().id = v_vid[ii];
                vertex_map[ii]->info().flag = v_flagv[ii];
            }
            else
            {
                vertex_map[ii]->info().id = tid;
                vertex_map[ii]->info().flag = 1;
            }
        }
        // 6) read the number of faces
        uint num_c = data.dmap[data.simplex_name].get_nbe_input();///(D+1);
        // 7) read and create the faces
        std::vector<Cell_handle> cell_map(num_c);
        uint ik;
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            Cell_handle ch = tri.tds().create_face();
            for (uint d = 0; d < D+1; d++)
            {
                ik = v_simplex[i*(D+1)+d];
                ch->set_vertex(d, vertex_map[ik]);
                vertex_map[ik]->set_cell(ch);
            }
            cell_map[ii] = ch;
            if(do_data)
            {
                ch->info().flag = v_flags[ii];
                ch->info().id = ii;
            }
            else
            {
                ch->info().flag = 1;
                ch->info().id = ii;
            }
        }
        // 8) read and construct neighbourhood relationships for faces
        for(uint j = 0; j < num_c; ++j)
        {
            Cell_handle ch  = cell_map[j];
            for(uint d = 0; d < D+1; d++)
            {
                ik = v_nb[j*(D+1)+d];
                ch->set_neighbor(d, cell_map[ik]);
            }
        }
        tri.set_infinite_vertex(vertex_map[0]);
        CGAL_triangulation_assertion(tri.is_valid());
    }



    std::istream& deserialize_b64_cgal(Delaunay_triangulation& tri,std::istream&  ifile)
    {
        std::vector<double> v_double;
        std::vector<int> v_int;
        std::vector<int> v_char;
        char cc;
        uint ldim;
        int num_c,num_v;
        ifile >> ldim;
        ifile >> num_v >> num_c;
        int D = ldim;
        tri.clear();
        tri.tds().set_dimension(D);
        auto cit = tri.all_cells_begin();
        tri.tds().delete_cell(cit);
        // 5) read and create the vertices
        deserialize_b64_vect<double>(v_double,ifile);
        std::vector<Vertex_handle> vertex_map(num_v+1);
        vertex_map[0] = tri.infinite_vertex();
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            std::vector<double> coords_v(D);
            for(uint d = 0; d < D; d++)
            {
                coords_v[d] = v_double[ii*D +d];
            }
            Point p(coords_v[0],coords_v[1],coords_v[2]);
            vertex_map[ii] = tri.tds().create_vertex(p);
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->info().id = v_int[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->info().gid = v_int[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->info().flag = v_int[ii];
        }
        uint ik;
        deserialize_b64_vect(v_int,ifile);
        num_c = v_int.size()/(D+1);
        std::vector<Cell_handle> cell_map(num_c);
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            Cell_handle ch = tri.tds().create_face();
            for (uint d = 0; d < D+1; d++)
            {
                ik = v_int[i*(D+1)+d];
                ch->set_vertex(d, vertex_map[ik]);
                vertex_map[ik]->set_cell(ch);
            }
            cell_map[ii] = ch;
        }
        deserialize_b64_vect(v_char,ifile);
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            cell_map[ii]->info().flag =  v_char[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            cell_map[ii]->info().gid =  v_int[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint j = 0; j < num_c; ++j)
        {
            Cell_handle ch  = cell_map[j];
            for(uint d = 0; d < D+1; d++)
            {
                ik = v_int[j*(D+1)+d];
                ch->set_neighbor(d, cell_map[ik]);
            }
        }
        assert(tri.is_valid());
        return ifile;
    }


    std::ostream & serialize_b64_cgal( const Delaunay_triangulation& tri,std::ostream & ofile) const
    {
        int D = 3;
        int num_c = tri.number_of_cells();
        uint num_v = tri.number_of_vertices();
        ofile << D << " " << num_v << " " << num_c << " ";
        std::vector<double> v_double;
        std::vector<int> v_int;
        std::vector<int> v_char;
        CGAL::Unique_hash_map<Vertex_const_handle, uint> vertex_map;
        vertex_map[tri.infinite_vertex()] = 0;
        uint i = 1;
        for(int d = 0; d < D; d++)
            v_double.push_back(d*100);
        for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
        {
            if(tri.is_infinite(vit))
            {
                continue;
            }
            int ii = i;
            for(int d = 0; d < D; d++)
            {
                double pv = vit->point()[d];
                v_double.push_back(pv);
            }
            vertex_map[vit] = ii;
            i++;
        }
        serialize_b64_vect(v_double,ofile);
        v_int.push_back(0);
        for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
        {
            if(tri.is_infinite(vit))
                continue;
            v_int.push_back(vit->info().id);
        }
        serialize_b64_vect(v_int,ofile);
        v_int.push_back(0);
        for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
        {
            if(tri.is_infinite(vit))
                continue;
            v_int.push_back(vit->info().gid);
        }
        serialize_b64_vect(v_int,ofile);
        v_int.push_back(0);
        for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
        {
            if(tri.is_infinite(vit))
                continue;
            v_int.push_back(vit->info().flag);
        }
        serialize_b64_vect(v_int,ofile);
        CGAL::Unique_hash_map<Cell_const_handle, uint> cell_map;
        i = 0;
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            int ii = i;
            cell_map[it] = ii;
            ++i;
            for(int d = 0; d < D+1; d++)
            {
                int vertex_id = vertex_map[it->vertex(d)] ;
                v_int.push_back(vertex_id);
            }
        }
        serialize_b64_vect(v_int,ofile);
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            v_char.push_back(it->info().flag);
        }
        serialize_b64_vect(v_char,ofile);
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            v_int.push_back(it->info().gid);
        }
        serialize_b64_vect(v_int,ofile);
        for(auto it = tri.cells_begin(); it != tri.cells_end(); ++it)
        {
            for(int j = 0; j < D+1; j++)
            {
                int nb_id = cell_map[it->neighbor(j)];
                v_int.push_back(nb_id);
            }
        }
        serialize_b64_vect(v_int,ofile);
        return ofile;
    }





};

}

#endif // DDT_CGAL_TRAITS_3_HPP

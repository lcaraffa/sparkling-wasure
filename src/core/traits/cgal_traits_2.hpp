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
#ifndef DDT_CGAL_TRAITS_2_HPP
#define DDT_CGAL_TRAITS_2_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/property_map.h>

#include "cgal_traits_base.hpp"
#include "Facet_const_iterator_2.hpp"

namespace ddt
{

template<typename I, typename F, typename G = bool>
struct Cgal_traits_2
{
    enum { D = 2 };
    typedef I                                                      Id;
    typedef F                                                      Flag_V;
    typedef G                                                      Flag_C;
    typedef ddt::Data<Id, Flag_V>                                  Data_V;
    typedef ddt::Data<Id, Flag_C>                                  Data_C;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
    typedef CGAL::Triangulation_vertex_base_with_info_2<Data_V, K> Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<Data_C, K>   Cb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Cb>            TDS;
    typedef typename K::Point_2                                    Point;
    typedef std::pair<Point,Id>                                    Point_id;
    typedef std::tuple<Point,Id,Id>                                Point_id_id;

    typedef typename TDS::Vertex_iterator                          Vertex_const_iterator;
    typedef typename TDS::Vertex_handle                            Vertex_const_handle;
    typedef typename TDS::Vertex_iterator                          Vertex_iterator;
    typedef typename TDS::Vertex_handle                            Vertex_handle;

    typedef typename TDS::Face_iterator                            Cell_const_iterator;
    typedef typename TDS::Face_handle                              Cell_const_handle;
    typedef typename TDS::Face_iterator                            Cell_iterator;
    typedef typename TDS::Face_handle                              Cell_handle;

    typedef std::pair<Cell_const_handle, int>                      Facet;
    typedef Facet_const_iterator_2<TDS>                            Facet_const_iterator;
    typedef Facet_const_iterator                                   Facet_const_handle;
    typedef Facet_const_iterator                                   Facet_iterator;
    typedef Facet_const_iterator                                   Facet_handle;

    typedef CGAL::Delaunay_triangulation_2<K, TDS>                 Delaunay_triangulation;
    typedef CGAL::Random_points_in_disc_2<Point>                   Random_points_in_ball;

    template<typename Iterator> static Point make_point(Iterator it)
    {
        Iterator x = it;
        Iterator y = ++it;
        return Point(*x,*y);
    }

    struct Random_points_in_box : CGAL::Random_points_in_square_2<Point>
    {
        Random_points_in_box(int d, double g, CGAL::Random& rnd = CGAL::get_default_random()) : CGAL::Random_points_in_square_2<Point>(g, rnd)
        {
            CGAL_assertion(d==2);
        }
        Random_points_in_box(double g, CGAL::Random& rnd = CGAL::get_default_random()) : CGAL::Random_points_in_square_2<Point>(g, rnd) {}
    };

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
    inline Flag_C& flag(Cell_const_handle c) const
    {
        return c->info().flag;
    }

    inline const Data_C& data(Cell_const_handle c)  const
    {
        return c->info();
    }

    inline const Data_V& data(Vertex_const_handle v)  const
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
        return dt.number_of_faces();
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
        return dt.all_faces_begin();
    }
    inline Cell_const_iterator cells_end(const Delaunay_triangulation& dt) const
    {
        return dt.all_faces_end();
    }

    inline Vertex_handle infinite_vertex(const Delaunay_triangulation& dt) const
    {
        return dt.infinite_vertex();
    }

    inline void clear(Delaunay_triangulation& dt) const
    {
        return dt.clear();
    }

    template<class It, class Out> inline int insert_splay(Delaunay_triangulation& dt, Id id, It begin, It end, Out out) const
    {
        assert(false);
        dt.insert(begin, end);
        return 0; //FIXME
    }

    template<typename PP>
    double squared_dist(const PP & p1, const PP & p2, int dim) const
    {
        double acc = 0;
        for(int d = 0; d < dim; d++)
            acc += (p1[d] - p2[d])*(p1[d] - p2[d]);
        return acc;
    }


    Vertex_handle insert_simplified_one(Delaunay_triangulation& dt, Id id, const Point& p, Id pid, Cell_handle hint, int& exist_count, int& foreign_count) const
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
            int li;
            Cell_handle c = dt.locate(p, lt, li, hint);
            if (lt == Delaunay_triangulation::VERTEX)
            {
                ++exist_count;
                return v; // Point already exists
            }
            std::vector<typename Delaunay_triangulation::Edge> facets;
            facets.reserve(32);
            dt.get_boundary_of_conflicts(p, std::back_inserter(facets), c);
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
                ++foreign_count;
                return v; // future star of "it" is entirely foreign, skipping
            }
            v = dt.insert(p, c);
        }
        v->info().id = pid;
        return v;
    }

    template<class It, class Out> void insert_simplified(Delaunay_triangulation& dt, Id id, It begin, It end, Out out) const
    {
        std::vector<std::size_t> indices;
        std::vector<Point> points;
        std::vector<Id> infos;
        std::size_t index=0;
        for (It it=begin; it!=end; ++it)
        {
            points.push_back( it->first  );
            infos.push_back ( it->second );
            indices.push_back(index++);
        }
        typedef typename CGAL::Pointer_property_map<Point>::type Pmap;
        typedef CGAL::Spatial_sort_traits_adapter_2<K,Pmap> Search_traits;
        spatial_sort(indices.begin(),
                     indices.end(),
                     Search_traits(make_property_map(points), dt.geom_traits()));
        int exist_count = 0;
        int foreign_count = 0;
        int inserted_count = 0;
        int vertex_count = dt.number_of_vertices();
        Cell_handle hint;
        for(auto it = indices.begin(); it != indices.end() ; ++it)
        {
            Vertex_handle v = insert_simplified_one(dt, id, points[*it], infos[*it], hint, exist_count, foreign_count);
            if (v != Vertex_handle())
            {
                hint = v->face();
                *out++ = v;
                ++inserted_count;
            }
        }
    }
    template<class It> inline void insert(Delaunay_triangulation& dt, It begin, It end) const
    {
        Vertex_handle hint;
        for(auto it=begin; it!=end; ++it)
        {
            {
                hint = dt.insert(it->first);
            }
            hint->info() = it->second;
        }
    }

    template<class It> inline void remove(Delaunay_triangulation& dt, It begin, It end) const
    {
        for(It it=begin; it!=end; ++it) dt.remove(*it);
    }

    template<class It> It get_local_convex_hull(const Delaunay_triangulation& dt, Id id, It out) const
    {
        auto vc = dt.incident_vertices(dt.infinite_vertex()),
             done(vc);
        if(vc != 0) do
            {
                if( vc->info().id == id ) *out++ = vc;
            }
            while (++vc != done);
        return out;
    }

    inline bool is_visible(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c, int i) const
    {
        return true;
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
        return dt.dual(c);
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

    inline bool vertex_is_infinite(const Delaunay_triangulation& dt, Vertex_const_handle v) const
    {
        return dt.is_infinite(v);
    }



    inline bool facet_is_infinite(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        for(int i = 0; i<=D; ++i)
            if(dt.is_infinite(vertex(dt, f->first, i)))
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
        return CGAL::to_double(p.x()*q.x()+p.y()*q.y());
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
        return c->neighbor(i)->index(c);
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

    template<typename Unary_op>
    void incident_cells(const Delaunay_triangulation& dt, Vertex_const_handle v, Unary_op op) const
    {
        auto c = dt.incident_faces(v), done(c);
        if (c != 0) do { op(c); }
            while(++c != done);
    }

    template<typename Unary_op>
    void finite_adjacent_vertices(const Delaunay_triangulation& dt, Vertex_const_handle v, Unary_op op) const
    {
        auto c = dt.incident_vertices(v), done(c);
        if (c != 0) do { if (!dt.is_infinite(c)) op(c); }
            while(++c != done);
    }



    template<typename DDT_DATA>
    void export_tri_to_data(Delaunay_triangulation& tri,DDT_DATA & data,bool do_init_id,bool add_data)
    {
        typedef typename DDT_DATA::Data_ply Data_ply;
        int D = data.D;
        data.dmap[data.xyz_name] = Data_ply(data.xyz_name,"vertex",D,D,data.get_float_type());
        data.dmap[data.simplex_name] = Data_ply(data.simplex_name,"face",D+1,D+1,data.get_int_type());
        data.dmap[data.nb_name] = Data_ply(data.nb_name,"face",D+1,D+1,data.get_int_type());
        if(add_data)
        {
            data.dmap[data.vid_name] = Data_ply(data.vid_name,"vertex",1,1,data.get_int_type());
            data.dmap[data.cid_name] = Data_ply(data.cid_name,"face",1,1,data.get_int_type());
            data.dmap[data.flag_vertex_name] = Data_ply(data.flag_vertex_name,"vertex",1,1,data.get_int_type());
            data.dmap[data.flag_simplex_name] = Data_ply(data.flag_simplex_name,"face",1,1,data.get_int_type());
        }
        std::vector<double> v_xyz;
        std::vector<int> v_simplex,v_nb,v_vid,v_cid,v_flagv,v_flags;
        std::map<Vertex_const_handle, uint> vertex_map;
        vertex_map[tri.infinite_vertex()] = 0;
        v_vid.push_back(0);
        v_flagv.push_back(0);
        for(int d = 0; d < D; d++)
        {
            v_xyz.push_back(0.0);
        }
        uint i = 1;
        for(auto vit = tri.all_vertices_begin(); vit != tri.all_vertices_end(); ++vit)
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
            if(add_data)
            {
                v_vid.push_back(vit->info().id);
                v_flagv.push_back(vit->info().flag);
            }
            vertex_map[vit] = ii;
            i++;
        }
        std::map<Cell_const_handle, uint> cell_map;
        i = 0;
        for(auto it = tri.all_faces_begin(); it != tri.all_faces_end(); ++it)
        {
            int ii = i;
            if(!do_init_id)
                int ii =  it->info().id;
            cell_map[it] = ii;
            ++i;
            for(int d = 0; d < D+1; d++)
            {
                int vertex_id = vertex_map[it->vertex(d)] ;
                v_simplex.push_back(vertex_id);
            }
            if(add_data)
            {
                v_flags.push_back(it->info().flag);
                v_cid.push_back(ii);
            }
        }
        for(auto it = tri.all_faces_begin(); it != tri.all_faces_end(); ++it)
        {
            for(int j = 0; j < D+1; j++)
            {
                int nb_id = cell_map[it->neighbor(j)];
                v_nb.push_back(nb_id);
            }
        }
        data.dmap[data.xyz_name].fill_full_output(v_xyz);
        data.dmap[data.simplex_name].fill_full_output(v_simplex);
        data.dmap[data.nb_name].fill_full_output(v_nb);
        if(add_data)
        {
            data.dmap[data.vid_name].fill_full_output(v_vid);
            data.dmap[data.cid_name].fill_full_output(v_cid);
            data.dmap[data.flag_vertex_name].fill_full_output(v_flagv);
            data.dmap[data.flag_simplex_name].fill_full_output(v_flags);
            data.dmap[data.vid_name].do_exist = true;
            data.dmap[data.cid_name].do_exist = true;
            data.dmap[data.flag_vertex_name].do_exist = true;
            data.dmap[data.flag_simplex_name].do_exist = true;
        }
        data.dmap[data.xyz_name].do_exist = true;
        data.dmap[data.simplex_name].do_exist = true;
        data.dmap[data.nb_name].do_exist = true;
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
        int D = 2;
        tri.clear();
        tri.tds().set_dimension(D);
        auto cit = tri.all_faces_begin();
        Cell_handle inf_ch = cit;
        tri.tds().delete_face(cit);
        uint num_v = data.dmap[data.xyz_name].get_nbe_input();
        std::vector<Vertex_handle> vertex_map(num_v);
        vertex_map[0] = tri.infinite_vertex();
        for(uint i = 1; i < num_v; ++i)
        {
            int ii = i;
            std::vector<double> coords_v(D);
            for(uint d = 0; d < D; d++)
            {
                coords_v[d] = v_xyz[ii*D +d];
            }
            Point p(coords_v[0],coords_v[1]);
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
        uint num_c = data.dmap[data.simplex_name].get_nbe_input();
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
                vertex_map[ik]->set_face(ch);
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

};

}

#endif // DDT_CGAL_TRAITS_2_HPP

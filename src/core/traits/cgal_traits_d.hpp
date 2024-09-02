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
#ifndef DDT_CGAL_TRAITS_D_HPP
#define DDT_CGAL_TRAITS_D_HPP

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Triangulation_ds_vertex.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>

#include "cgal_traits_base.hpp"
#include "Facet_const_iterator_d.hpp"
#include <functional>
#include <unordered_set>
#include <limits>
#include <math.h>
#include "../io/number_parser.hpp"
typedef std::numeric_limits< double > dbl;

namespace ddt
{


template< typename Dim = CGAL::Dynamic_dimension_tag>
struct Cgal_traits_raw_d
{
    typedef Dim                                                    Dim_tag;
    typedef CGAL::Epick_d<Dim_tag>                                 K;
    typedef CGAL::Delaunay_triangulation<K>                 Delaunay_triangulation;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_handle_raw;
    typedef typename Delaunay_triangulation::Full_cell_handle                            Cell_handle_raw;

    typedef typename Delaunay_triangulation::Vertex_iterator                    Vertex_const_iterator;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_const_handle;
    typedef typename Delaunay_triangulation::Vertex_iterator                          Vertex_iterator;
    typedef typename Delaunay_triangulation::Vertex_handle                            Vertex_handle;

    typedef typename Delaunay_triangulation::Full_cell_const_iterator                 Cell_const_iterator;
    typedef typename Delaunay_triangulation::Full_cell_const_handle                   Cell_const_handle;
    typedef typename Delaunay_triangulation::Full_cell_iterator                       Cell_iterator;
    typedef typename Delaunay_triangulation::Full_cell_handle                         Cell_handle;
    typedef typename Delaunay_triangulation::Full_cell::Vertex_handle_iterator Vertex_h_iterator;


    typedef typename Delaunay_triangulation::Facet_iterator                 Facet_iterator;
    typedef Facet_iterator                   Facet_const_handle;

    typedef typename K::Point_d                                    Point;
    typedef std::pair<Cell_const_handle, int>                      Facet;

    typedef typename CGAL::Unique_hash_map<Vertex_handle_raw, uint> v_hmap_uint;


    inline int current_dimension(const Delaunay_triangulation& dt) const
    {
        return dt.current_dimension();
    }

    Delaunay_triangulation triangulation(int dimension) const
    {
        return Delaunay_triangulation(dimension);
    }

    inline size_t number_of_cells(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_full_cells();
    }
    inline size_t number_of_vertices(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_vertices();
    }

    inline Cell_const_iterator cells_begin(const Delaunay_triangulation& dt) const
    {
        return dt.full_cells_begin();
    }
    inline Cell_const_iterator cells_end(const Delaunay_triangulation& dt) const
    {
        return dt.full_cells_end();
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
        return f->index_of_covertex();
    }

    inline Cell_const_handle full_cell(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        return f->full_cell();
    }

    inline Cell_const_iterator neighbor(const Delaunay_triangulation& dt, Cell_const_iterator c, int i) const
    {
        return c->neighbor(i);
    }

    int
    Cell2lp(const Cell_const_handle & ch,   std::vector<Point> & lp) const
    {
        for(auto cit = ch->vertices_begin();
                cit != ch->vertices_end();
                ++cit)
        {
            lp.push_back((*cit)->point());
        }
        return 0;
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
            if(dist >= fabs(center[d] - bbox.max(d)) || dist >= fabs(center[d] - bbox.min(d)))
                return false;
        }
        return true;
    }

    Point get_center( Cell_const_handle  ch) const
    {
        int dim = ch->maximal_dimension();
        std::vector<Point> lpf;
        Cell2lp(ch,lpf);
        CGAL::Sphere_d<K> sph1(dim,lpf.begin(),lpf.end());
        return sph1.center();
    }

    inline Point circumcenter(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        return get_center(c);
    }
};

template<typename DV, typename DC, typename Dim = CGAL::Dynamic_dimension_tag>
struct Cgal_traits_d
{
    typedef Dim                                                    Dim_tag;
    typedef DV                                                     Data_V;
    typedef DC                                                     Data_C;

    typedef typename DV::Id                                          Id;
    typedef typename DV::Flag                                        Flag_V;
    typedef typename DC::Flag                                        Flag_C;



    typedef CGAL::Epick_d<Dim_tag>                                 K;
    typedef CGAL::Triangulation_vertex<K,Data_V>                   Vb;
    typedef CGAL::Triangulation_full_cell<K,Data_C>                Cb;
    typedef CGAL::Triangulation_data_structure<Dim_tag,Vb,Cb>      TDS;
    typedef typename K::Point_d                                    Point;
    typedef std::pair<Point,Id>                                    Point_id;
    typedef std::tuple<Point,Id,Id>                                Point_id_id;

    typedef typename TDS::Face                                     Face;
    typedef typename TDS::Vertex_const_iterator                          Vertex_const_iterator;
    typedef typename TDS::Vertex_const_handle                            Vertex_const_handle;
    typedef typename TDS::Vertex_iterator                          Vertex_iterator;
    typedef typename TDS::Vertex_handle                            Vertex_handle;


    typedef typename TDS::Full_cell_const_iterator                 Cell_const_iterator;
    typedef typename TDS::Full_cell_const_handle                   Cell_const_handle;
    typedef typename TDS::Full_cell_iterator                       Cell_iterator;
    typedef typename TDS::Full_cell_handle                         Cell_handle;

    typedef std::pair<Cell_const_handle, int>                      Facet;
    typedef Facet_const_iterator_d<TDS>                            Facet_const_iterator;
    typedef Facet_const_iterator                                   Facet_const_handle;
    typedef Facet_const_iterator                                   Facet_iterator;
    typedef Facet_const_iterator                                   Facet_handle;

    typedef CGAL::Delaunay_triangulation<K,TDS>                    Delaunay_triangulation;
    typedef CGAL::Random_points_in_cube_d<Point>                   Random_points_in_box;

    typedef  CGAL::Unique_hash_map<Vertex_const_handle,bool> v_hmap_bool;
    typedef  CGAL::Unique_hash_map<Vertex_handle,uint> v_hmap_uint;

    typedef typename Delaunay_triangulation::Full_cell::Vertex_handle_iterator Vertex_h_iterator;

    Delaunay_triangulation triangulation(int dimension) const
    {
        return Delaunay_triangulation(dimension);
    }

    inline Id    id  (Vertex_const_handle v) const
    {
        return v->data().id;
    }
    inline Id    gid  (Vertex_const_handle v) const
    {
        return v->data().gid;
    }

    inline Id    id  (Cell_const_handle v) const
    {
        return v->data().id;
    }
    inline Id    gid  (Cell_const_handle v) const
    {
        return v->data().gid;
    }



    inline Flag_V& flag(Vertex_const_handle v) const
    {
        return v->data().flag;
    }
    inline Flag_C& flag(Cell_const_handle c) const
    {
        return c->data().flag;
    }

    inline const Data_C& data(Cell_const_handle c)  const
    {
        return c->data();
    }

    inline const Data_V& data(Vertex_const_handle v)  const
    {
        return v->data();
    }

    inline int current_dimension(const Delaunay_triangulation& dt) const
    {
        return dt.current_dimension();
    }
    inline int maximal_dimension(const Delaunay_triangulation& dt) const
    {
        return dt.maximal_dimension();
    }
    inline size_t number_of_cells(const Delaunay_triangulation& dt) const
    {
        return dt.number_of_full_cells();
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
        return dt.full_cells_begin();
    }
    inline Cell_const_iterator cells_end(const Delaunay_triangulation& dt) const
    {
        return dt.full_cells_end();
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
        insert(dt, begin, end);
        return 0; //FIXME
    }

    template<class It, class Out> void insert_simplified(Delaunay_triangulation& dt, Id id, It begin, It end, Out out) const
    {
        Vertex_handle hint;
        for(auto it=begin; it!=end; ++it)
        {
            if (hint != Vertex_handle())
            {
                hint = dt.insert(it->first, hint);
            }
            else
            {
                hint = dt.insert(it->first);
            }
            hint->data() = it->second;
            *out++ = hint;
            //FIXME : simplify
        }
    }
    template<class It> inline void insert(Delaunay_triangulation& dt, It begin, It end) const
    {
        Vertex_handle hint;
        for(auto it=begin; it!=end; ++it)
        {
            if (hint != Vertex_handle())
            {
                hint = dt.insert(it->first, hint);
            }
            else
            {
                hint = dt.insert(it->first);
            }
            hint->data() = it->second;
        }
    }

    template<class It> inline void remove(Delaunay_triangulation& dt, It begin, It end) const
    {
        dt.remove(begin, end);
    }

    template<class It> inline It get_local_convex_hull(const Delaunay_triangulation& dt, Id id, It out) const
    {
        std::vector<Point> chull;
        std::vector<Face> edges;
        std::back_insert_iterator<std::vector<Face>> out_e(edges);
        dt.incident_faces(dt.infinite_vertex(), 1, out_e);
        for(auto it : edges)
        {
            for(int df = 0; df < 2; df++)
            {
                Vertex_handle vh = it.vertex(df);
                if( ! dt.is_infinite(vh) )
                {
                    *out++ = vh;
                }
            }
        }
        return out;
    }



    void get_list_vertices(Cell_const_handle fch,std::list<Vertex_const_handle> & lp)
    {
        for(auto vht = fch->vertices_begin() ;
                vht != fch->vertices_end() ;
                ++vht)
        {
            Vertex_handle v = *vht;
            lp.push_back(v);
        }
    }

    template<int D>
    inline bool is_visible(const Delaunay_triangulation& dt, const Bbox<D>& bbox, Cell_const_handle c, int i) const
    {
        return true;
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
            if(dist >= fabs(center[d] - bbox.max(d)) || dist >= fabs(center[d] - bbox.min(d)))
                return false;
        }
        return true;
    }


    template<int D>
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


    int
    Cell2lp(const Cell_const_handle & ch,   std::vector<Point> & lp) const
    {
        for(auto cit = ch->vertices_begin();
                cit != ch->vertices_end();
                ++cit)
        {
            lp.push_back((*cit)->point());
        }
        return 0;
    }


    Point get_center( Cell_const_handle  ch) const
    {
        int dim = ch->maximal_dimension();
        std::vector<Point> lpf;
        Cell2lp(ch,lpf);
        CGAL::Sphere_d<K> sph1(dim,lpf.begin(),lpf.end());
        return sph1.center();
    }


    inline Point circumcenter(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        return get_center(c);
    }

    inline bool vertex_is_infinite(const Delaunay_triangulation& dt, Vertex_const_handle v) const
    {
        return dt.is_infinite(v);
    }

    inline bool facet_is_infinite(const Delaunay_triangulation& dt, Facet_const_handle f) const
    {
        for(int i = 0; i<=dt.current_dimension(); ++i)
            if(i!=f->second && dt.is_infinite(vertex(dt, f->first, i)))
                return true;
        return false;
    }

    inline bool cell_is_infinite(const Delaunay_triangulation& dt, Cell_const_handle c) const
    {
        for(int i = 0; i<=dt.current_dimension(); ++i)
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
        std::cerr << "scalar_product not implemented in cgal_traits_d" << std::endl;
        assert(false);
        return 0;
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
        for(auto v1 = (c1)->vertices_begin(); v1 != (c1)->vertices_end(); ++v1 )
        {
            bool is_equal = false;
            for(auto v2 = (c2)->vertices_begin(); v2 != (c2)->vertices_end(); ++v2 )
            {
                if(are_vertices_equal(t1, *v1, t2, *v2))
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
        return c->mirror_index(i);
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




    std::vector<double> get_cell_barycenter(Cell_const_handle ch)
    {
        int D = ch->maximal_dimension();
        std::vector<double>  coords(D,0.0);
        for(auto vht = ch->vertices_begin() ;
                vht != ch->vertices_end() ;
                ++vht)
        {
            Vertex_handle v = *vht;
            for(uint d = 0; d < D; d++)
            {
                coords[d] += (v->point())[d];
            }
        }
        for(uint d = 0; d < D; d++)
            coords[d] /= ((double)D+1);
        return  coords;
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
        tri.set_current_dimension(D);
        auto cit = tri.full_cells_begin();
        Cell_handle inf_ch = cit;
        tri.tds().delete_full_cell(inf_ch);
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
            Point p(D,coords_v.begin(),coords_v.end());
            vertex_map[ii] = tri.new_vertex(p);
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->data().id = v_int[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->data().gid = v_int[ii];
        }
        deserialize_b64_vect(v_int,ifile);
        for(uint i = 1; i <= num_v; ++i)
        {
            int ii = i;
            vertex_map[ii]->data().flag = v_int[ii];
        }
        uint ik;
        deserialize_b64_vect(v_int,ifile);
        num_c = v_int.size()/(D+1);
        std::vector<Cell_handle> cell_map(num_c);
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            Cell_handle ch = tri.new_full_cell();
            for (uint d = 0; d < D+1; d++)
            {
                ik = v_int[i*(D+1)+d];
                ch->set_vertex(d, vertex_map[ik]);
                vertex_map[ik]->set_full_cell(ch);
            }
            cell_map[ii] = ch;
        }
        for(uint i = 0; i < num_c; ++i)
        {
            int ii = i;
            cell_map[ii]->data() =  0;
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
        for(uint j = 0; j < num_c; ++j)
        {
            Cell_handle s  = cell_map[j];
            for( uint j = 0; j <= D; ++j )
            {
                if( -1 != s->mirror_index(j) )
                    continue;
                Cell_handle n = s->neighbor(j);
                int k = 0;
                Cell_handle nn = n->neighbor(k);
                while( s != nn )
                    nn = n->neighbor(++k);
                s->set_mirror_index(j,k);
                n->set_mirror_index(k,j);
            }
        }
        assert(tri.is_valid());
        return ifile;
    }


    std::ostream & serialize_b64_cgal( const Delaunay_triangulation& tri,std::ostream & ofile) const
    {
        int D = tri.current_dimension();
        int num_c = tri.number_of_full_cells();
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
        if(true)
        {
            v_int.push_back(0);
            for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
            {
                if(tri.is_infinite(vit))
                    continue;
                v_int.push_back(vit->data().id);
            }
            serialize_b64_vect(v_int,ofile);
            v_int.push_back(0);
            for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
            {
                if(tri.is_infinite(vit))
                    continue;
                v_int.push_back(vit->data().gid);
            }
            serialize_b64_vect(v_int,ofile);
            v_int.push_back(0);
            for(auto vit = tri.vertices_begin(); vit != tri.vertices_end(); ++vit)
            {
                if(tri.is_infinite(vit))
                    continue;
                v_int.push_back(vit->data().flag);
            }
            serialize_b64_vect(v_int,ofile);
        }
        CGAL::Unique_hash_map<Cell_const_handle, uint> cell_map;
        i = 0;
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
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
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
        {
            v_char.push_back(it->data().flag);
        }
        serialize_b64_vect(v_char,ofile);
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
        {
            v_int.push_back(it->data().gid);
        }
        serialize_b64_vect(v_int,ofile);
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
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
            v_vid.push_back(vit->data().id);
            v_flagv.push_back(vit->data().flag);
            vertex_map[vit] = ii;
            i++;
        }
        n = tri.number_of_full_cells();
        v_simplex.resize(n*(D+1));
        v_flags.resize(n);
        v_nb.resize(n*(D+1));
        std::map<Cell_const_handle, uint> cell_map;
        i = 0;
        int max_id = 0;
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
        {
            if(it->data().gid > max_id)
                max_id = it->data().gid;
        }
        max_id++;
        for(auto it = tri.full_cells_begin(); it != tri.full_cells_end(); ++it)
        {
            int ii =  it->data().gid;
            if(ii == -1)
                ii = max_id++;
            cell_map[it] = ii;
            for(int d = 0; d < D+1; d++)
            {
                int vertex_id = vertex_map[it->vertex(d)] ;
                v_simplex[ii*(D+1)+d] = vertex_id;
            }
            for(int j = 0; j < D+1; j++)
            {
                int nb_id = cell_map[it->neighbor(j)];
                v_nb[ii*(D+1) +j] = nb_id;
            }
            v_flags[ii] = it->data().flag;
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
    void build_tri_from_data(Delaunay_triangulation& tri,DDT_DATA & data,bool do_clean_data,int tid)
    {
        uint num_v = data.dmap[data.xyz_name].get_nbe_input();
        uint num_c = data.dmap[data.simplex_name].get_nbe_input();///(D+1);
        int D = data.D;
        std::vector<double> v_xyz;
        std::vector<int> v_simplex,v_nb,v_vid,v_cid,v_flagv,v_flags;
        data.dmap[data.xyz_name].extract_full_input(v_xyz,do_clean_data);
        data.dmap[data.simplex_name].extract_full_input(v_simplex,do_clean_data);
        data.dmap[data.nb_name].extract_full_input(v_nb,do_clean_data);
        data.dmap[data.flag_vertex_name].extract_full_input(v_flagv,do_clean_data);
        data.dmap[data.flag_simplex_name].extract_full_input(v_flags,do_clean_data);
        data.dmap[data.vid_name].extract_full_input(v_vid,do_clean_data);
        tri.set_current_dimension(D);
        auto cit = tri.full_cells_begin();
        Cell_handle inf_ch = cit;
        tri.tds().delete_full_cell(inf_ch);
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
            Point p(D,coords_v.begin(),coords_v.end());
            vertex_map[ii] = tri.new_vertex(p);
            vertex_map[ii]->data().id = v_vid[ii];
            vertex_map[ii]->data().flag = v_flagv[ii];
        }
        std::vector<Cell_handle> cell_map(num_c);
        uint ik;
        for(uint ii = 0; ii < num_c; ++ii)
        {
            Cell_handle ch = tri.new_full_cell();
            for (uint d = 0; d < D+1; d++)
            {
                ik = v_simplex[ii*(D+1)+d];
                ch->set_vertex(d, vertex_map[ik]);
                vertex_map[ik]->set_full_cell(ch);
            }
            cell_map[ii] = ch;
            ch->data().flag = v_flags[ii];
            ch->data().id = ii;
        }
        for(uint jj = 0; jj < num_c; ++jj)
        {
            Cell_handle ch  = cell_map[jj];
            for(uint d = 0; d < D+1; d++)
            {
                ik = v_nb[jj*(D+1)+d];
                ch->set_neighbor(d, cell_map[ik]);
            }
        }
        // compute the mirror indices
        for(uint j = 0; j < num_c; ++j)
        {
            Cell_handle s  = cell_map[j];
            for( uint j = 0; j <= D; ++j )
            {
                if( -1 != s->mirror_index(j) )
                    continue;
                Cell_handle n = s->neighbor(j);
                int k = 0;
                Cell_handle nn = n->neighbor(k);
                while( s != nn )
                    nn = n->neighbor(++k);
                s->set_mirror_index(j,k);
                n->set_mirror_index(k,j);
            }
        }
    }



    template<typename Unary_op>
    void incident_cells(const Delaunay_triangulation& dt, Vertex_const_handle v, Unary_op op) const
    {
        std::vector<Cell_const_handle> cells;
        cells.reserve(32);
        dt.incident_full_cells(v, std::back_inserter(cells));
        for(auto c : cells) op(c);
    }

    template<typename Unary_op>
    void finite_adjacent_vertices(const Delaunay_triangulation& dt, Vertex_const_handle v, Unary_op op) const
    {
        std::unordered_set<Vertex_handle> vertices;
        uint ldim = dt.maximal_dimension();
        incident_cells(dt, v, [&vertices,ldim](Cell_const_handle c)
        {
            for(uint d = 0; d < ldim+1; d++)
                vertices.insert(c->vertex(d));
        });
        for(auto w : vertices) if(v != w && !dt.is_infinite(w)) op(w);
    }



};



template<unsigned int N, typename DataV, typename DataC>
struct Cgal_traits : public Cgal_traits_d<DataV,DataC,CGAL::Dimension_tag<N>>
{
    enum { D = N };
    typedef typename Cgal_traits_d<DataV,DataC,CGAL::Dimension_tag<N>>::Point Point;
    template<typename Iterator> Point make_point(Iterator it) const
    {
        return Point(D, it, it+D);
    }
};

template<unsigned int N>
struct Cgal_traits_raw : public Cgal_traits_raw_d<CGAL::Dimension_tag<N>>
{
    enum { D = N };
};
}


#endif // DDT_CGAL_TRAITS_D_HPP

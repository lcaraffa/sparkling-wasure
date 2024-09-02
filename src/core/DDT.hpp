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
#ifndef DDT_DDT_HPP
#define DDT_DDT_HPP

#include <string>
#include <unordered_map>
#include <type_traits>

#include "iterator/Vertex_const_iterator.hpp"
#include "iterator/Facet_const_iterator.hpp"
#include "iterator/Cell_const_iterator.hpp"
#include "iterator/Map_iterators.hpp"

#include "tile.hpp"

namespace ddt
{

template<typename _Traits, typename _Tile = ddt::Tile<_Traits>>
class DDT
{
public:
    typedef _Traits                                  Traits;
    typedef _Tile                                    Tile;

    typedef typename Traits::Point                   Point;
    typedef typename Traits::Id                      Id;
    typedef typename Traits::Delaunay_triangulation  DT;
    typedef typename Traits::Vertex_handle           Tile_vertex_handle;
    typedef typename Traits::Vertex_iterator         Tile_vertex_iterator;
    typedef typename Traits::Vertex_const_handle     Tile_vertex_const_handle;
    typedef typename Traits::Vertex_const_iterator   Tile_vertex_const_iterator;
    typedef typename Traits::Cell_handle             Tile_cell_handle;
    typedef typename Traits::Cell_const_handle       Tile_cell_const_handle;
    typedef typename Traits::Cell_const_iterator     Tile_cell_const_iterator;
    typedef typename Traits::Facet_handle            Tile_facet_handle;
    typedef typename Traits::Facet_const_handle      Tile_facet_const_handle;
    typedef typename Traits::Facet_const_iterator    Tile_facet_const_iterator;

    typedef std::pair<Tile_cell_const_handle,Id>     Tile_cell_const_handle_and_id;
    typedef std::pair<Tile_vertex_const_handle,Id>   Tile_vertex_const_handle_and_id;

    typedef ddt::Vertex_const_iterator<DDT>          Vertex_const_iterator;
    typedef ddt::Facet_const_iterator <DDT>          Facet_const_iterator;
    typedef ddt::Cell_const_iterator  <DDT>          Cell_const_iterator;

    typedef std::map<Id, Tile>                                              Tile_container;
    typedef Mapped_const_iterator<typename Tile_container::const_iterator>  Tile_const_iterator ;
    typedef Mapped_iterator<typename Tile_container::iterator>              Tile_iterator ;
    typedef Key_const_iterator<typename Tile_container::const_iterator>     Tile_id_const_iterator ;

    enum { D = Traits::D };

    inline int maximal_dimension() const
    {
        return D;
    }

    DDT() :
        tiles(),
        number_of_vertices_(0),
        number_of_facets_  (0),
        number_of_cells_   (0)
    {
    }

    DDT(const DDT& ddt) :
        tiles(ddt.tiles),
        number_of_vertices_(ddt.number_of_vertices_),
        number_of_facets_  (ddt.number_of_facets_  ),
        number_of_cells_   (ddt.number_of_cells_   )
    {
    }

    inline size_t number_of_cells   () const { return number_of_cells_;    }
    inline size_t number_of_vertices() const { return number_of_vertices_; }
    inline size_t number_of_facets  () const { return number_of_facets_;   }
    inline size_t number_of_tiles   () const { return tiles.size();   }

    Vertex_const_iterator vertices_begin() const { return Vertex_const_iterator(tiles_begin(), tiles_end()); }
    Vertex_const_iterator vertices_end  () const { return Vertex_const_iterator(tiles_begin(), tiles_end(), tiles_end()); }

    Cell_const_iterator cells_begin() const { return Cell_const_iterator(tiles_begin(), tiles_end()); }
    Cell_const_iterator cells_end  () const { return Cell_const_iterator(tiles_begin(), tiles_end(), tiles_end()); }

    Facet_const_iterator facets_begin() const { return Facet_const_iterator(tiles_begin(), tiles_end()); }
    Facet_const_iterator facets_end  () const { return Facet_const_iterator(tiles_begin(), tiles_end(), tiles_end()); }


    Tile_vertex_const_iterator tile_vertices_begin(Id tid) const { return tiles[tid].vertices_begin();}
    Tile_vertex_const_iterator tile_vertices_end  (Id tid) const { return tiles[tid].vertices_end();}

    Tile_cell_const_iterator tile_cells_begin(Id tid) const { return tiles[tid].cells_begin() ;}
    Tile_cell_const_iterator tile_cells_end  (Id tid) const { return tiles[tid].cells_end() ;}

    Tile_facet_const_iterator tile_facets_begin(Id tid) const { return tiles[tid].facets_begin();}
    Tile_facet_const_iterator tile_facets_end  (Id tid) const { return tiles[tid].facets_end();}


    Tile_id_const_iterator tile_ids_begin() const { return tiles.begin(); }
    Tile_id_const_iterator tile_ids_end  () const { return tiles.end  (); }

    Tile_const_iterator tiles_begin  () const { return tiles.begin (); }
    Tile_const_iterator tiles_end    () const { return tiles.end   (); }
    Tile_const_iterator get_const_tile(Id id) const { return tiles.find(id); }

    Tile_iterator tiles_begin  () { return tiles.begin (); }
    Tile_iterator tiles_end    () { return tiles.end   (); }
    Tile_iterator get_tile(Id id) { return tiles.find(id); }

    int vertex_id(Vertex_const_iterator v) const
    {
        if (v->is_infinite()) return -1;
        return std::distance(vertices_begin(), v->main());
    }

    int cell_id(Cell_const_iterator c) const
    {
        return std::distance(cells_begin(), c->main());
    }

    template<class F, class Scheduler> void for_each(Scheduler& s, F&& f = F())       { return s.for_each(tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> void for_each(Scheduler& s, F&& f = F()) const { return s.for_each(tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile      &)>::type for_each_rec (Scheduler& s, F&& f = F())       { return s.for_each_rec (tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile const&)>::type for_each_rec (Scheduler& s, F&& f = F()) const { return s.for_each_rec (tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile      &)>::type for_each_rec2(Scheduler& s, F&& f = F())       { return s.for_each_rec2(tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile const&)>::type for_each_rec2(Scheduler& s, F&& f = F()) const { return s.for_each_rec2(tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile      &)>::type transform_sum(Scheduler& s, F&& f = F())       { return s.transform_sum(tiles_begin(), tiles_end(), f); }
    template<class F, class Scheduler> typename std::result_of<F&&(Tile const&)>::type transform_sum(Scheduler& s, F&& f = F()) const { return s.transform_sum(tiles_begin(), tiles_end(), f); }

    template<class Scheduler> int insert_splay(Scheduler& s, bool do_simplify=true) { return transform_sum(s, s.insert_splay(do_simplify)); }
    template<class Scheduler> int insert(Scheduler& s, bool do_simplify=true) { return transform_sum(s, s.insert(do_simplify)); }
    template<class F, class Scheduler> int send_all(Scheduler& s, F f=F()) const { return transform_sum(s, s.send_all_func(tile_ids_begin(), tile_ids_end(), f)); }
    template<class F, class Scheduler> int insert_simplified(Scheduler& s, bool do_simplify=true, F f=F()) { return transform_sum(s, s.insert_simplified(f, do_simplify)); }
    template<class F, class Scheduler> int splay(Scheduler& s, bool do_simplify=true, F f=F())    { return transform_sum(s, s.splay_func(f, do_simplify)); }
    template<class F, class Scheduler> int splay_rec(Scheduler& s, bool do_simplify=true, F f=F())    { return for_each_rec(s, s.splay_func(f, do_simplify)); }
    template<class F, class Scheduler> int splay_rec2(Scheduler& s, bool do_simplify=true, F f=F())    { return for_each_rec2(s, s.splay_func(f, do_simplify)); }

    void init(int id)
    {
        tiles.emplace(id, id);
    }


    bool tile_is_loaded(int id)
    {
        return tiles.find(id) != tiles.end();
    }

    void clear(int NT)
    {
        tiles.clear();
        for(int id=0; id<NT; ++id)
        {
            tiles.emplace(id, id);
        }
    }
    template<class Scheduler, class Iterator, class Partitioner>
    int send_points(Scheduler& s, Iterator it, int count, Partitioner& part)
    {
        int res = count;
        for(; count; --count, ++it)
        {
            Point p {*it};
            int id = part(p);
            if (tiles.find(id) == tiles.end())
                tiles.emplace(id, id);
            s.send(p,id);
        }
        return res;
    }

    void get_adjacency_graph(std::unordered_multimap<Id,Id>& edges) const
    {
        for(Tile_const_iterator tile = tiles_begin(); tile != tiles_end(); ++tile)
        {
            std::set<Id> out_edges;
            tile->get_adjacency_graph_edges(out_edges);
            Id source = tile->id();
            for(Id target : out_edges)
                edges.emplace(source, target);
        }
    }

    bool is_adjacency_graph_symmetric() const
    {
        std::unordered_multimap<Id,Id> edges;
        std::unordered_multimap<Id,Id> reversed;
        get_adjacency_graph(edges);
        for(auto& edge : edges)
            reversed.emplace(edge.second, edge.first);
        return edges == reversed;
    }

    void get_ring(Cell_const_iterator c, int deg, std::set<Cell_const_iterator>& cset) const
    {
        std::set<Cell_const_iterator> seeds;
        c = c->main();
        cset.insert(c);
        seeds.insert(c);
        for(int i=0; i<deg; ++i)
        {
            std::set<Cell_const_iterator> next;
            next_ring(seeds, next);
            cset.insert(next.begin(), next.end());
            seeds.swap(next);
        }
    }

    void next_ring(const std::set<Cell_const_iterator>& seeds, std::set<Cell_const_iterator>& next) const
    {
        for(auto seed : seeds)
        {
            for(int d = 0; d <= D; d++)
            {
                auto c = seed->neighbor(d)->main();
                if(seeds.find(c) == seeds.end())
                    next.insert(c);
            }
        }
    }

    template<class Scheduler> void finalize(Scheduler &s)
    {
        finalize_vertices(s);
        finalize_facets(s);
        finalize_cells(s);
    }

    template<class Scheduler> void finalize_vertices(Scheduler &s) { number_of_vertices_ = transform_sum(s, std::mem_fn(&Tile::finalize_vertices)); }
    template<class Scheduler> void finalize_facets(Scheduler &s)   { number_of_facets_ = transform_sum(s, std::mem_fn(&Tile::finalize_facets)); }
    template<class Scheduler> void finalize_cells(Scheduler &s)    { number_of_cells_ = transform_sum(s, std::mem_fn(&Tile::finalize_cells)); }

    bool tile_is_valid(const Tile& tile) const
    {
        if (!tile.is_valid())
        {
            return false;
        }
        for(auto v = tile.vertices_begin(); v != tile.vertices_end(); ++v)
        {
            if(tile.vertex_is_infinite(v)) continue;
            Id tid = tile.id(v);
            if(tid == tile.id()) continue;
            auto t = get_tile(tid);
            if(t->locate_vertex(tile, v) == t->vertices_end())
            {
                return false;
            }
        }
        for(auto f = tile.facets_begin(); f != tile.facets_end(); ++f)
        {
            if(!tile.facet_is_mixed(f)) continue;
            std::set<Id> tids;
            for(int d = 0; d <= tile.current_dimension(); ++d)
            {
                if(d==tile.index_of_covertex(f)) continue;
                auto c = tile.full_cell(f);
                auto v = tile.vertex(c, d);
                if(tile.vertex_is_infinite(v)) continue;
                Id tid = tile.id(v);
                if(tid == tile.id()) continue;
                tids.insert(tid);
            }
            for(auto tid : tids)
            {
                auto t = get_tile(tid);
                if(t->locate_facet(tile, f) == t->facets_end())
                {
                    return false;
                }
            }
        }
        for(auto c = tile.cells_begin(); c != tile.cells_end(); ++c)
        {
            if(!tile.cell_is_mixed(c)) continue;
            std::set<Id> tids;
            for(int d = 0; d <= tile.current_dimension(); ++d)
            {
                auto v = tile.vertex(c, d);
                if(tile.vertex_is_infinite(v)) continue;
                Id tid = tile.id(v);
                if(tid == tile.id()) continue;
                tids.insert(tid);
            }
            for(auto tid : tids)
            {
                auto t = get_tile(tid);
                if(t->locate_cell(tile, c) == t->cells_end())
                {
                    std::cerr << "locate_cell failed" << std::endl;
                    return false;
                }
            }
        }
        return true;
    }



    void init_local_id()
    {
        typedef typename Traits::Data_C                    Data_C;
        typedef typename Traits::Data_V                    Data_V;
        int nextid = 0;
        int D = Traits::D;
        int acc = 0;
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit)
        {
            const Data_V & vd = vit->vertex_data();
            Data_V & vd_quickndirty = const_cast<Data_V &>(vd);
            Id main_id = vit->main_id();
            vd_quickndirty.gid =  (acc++);
        }
        acc = 0;
        for(auto iit = cells_begin(); iit != cells_end(); ++iit)
        {
            const Data_C & cd = iit->cell_data();
            Data_C & cd_quickndirty = const_cast<Data_C &>(cd);
            Id main_id = iit->main_id();
            cd_quickndirty.gid =  (acc++);
        }
    }


    template<class Scheduler> bool is_valid(Scheduler &s) const
    {
        bool valid = true;
        size_t invalid_tiles = transform_sum(s, [this](const Tile& t) { return size_t(!this->tile_is_valid(t)); });
        if (invalid_tiles) { std::cerr << invalid_tiles << " invalid tile(s)" << std::endl; valid = false; }
        if (number_of_vertices_ != transform_sum(s, std::mem_fn(&Tile::number_of_main_vertices))) { std::cerr << "incorrect number_of_vertices" << std::endl; valid = false; }
        if (number_of_facets_ != transform_sum(s, std::mem_fn(&Tile::number_of_main_facets))) { std::cerr << "incorrect number_of_facets" << std::endl; valid = false; }
        if (number_of_cells_ != transform_sum(s, std::mem_fn(&Tile::number_of_main_cells))) { std::cerr << "incorrect number_of_cells" << std::endl; valid = false; }
        return valid;
    }

    template<class Scheduler> bool sanity_check(Scheduler &s) const
    {
        typedef std::map<std::vector<Id>,size_t> count_map;
        typedef std::unordered_map<Id,count_map> count_map_map;
        std::unordered_map<Id,count_map_map> counts1, counts2;
        std::for_each(tile_ids_begin(), tile_ids_end(), [&counts1,&counts2](Id id) { counts1[id]; counts2[id]; });
        s.for_each(tiles_begin(), tiles_end(), [&counts1](const Tile& t) { t.sanity_check_send(counts1.at(t.id())); });
        s.for_each(tile_ids_begin(), tile_ids_end(), [&counts1,&counts2](Id id)
        {
            for(auto it = counts1.begin(); it != counts1.end(); ++it)
            {
                auto it1 = it->second.find(id);
                if ( it1 != it->second.end() ) counts2[id][it->first].swap(it1->second);
            }
        });
        return s.all_of(tiles_begin(), tiles_end(), [&counts2](const Tile& t) { return t.sanity_check_recv(counts2.at(t.id())); });
    }






private:

    Tile_container tiles;
    size_t number_of_vertices_;
    size_t number_of_facets_;
    size_t number_of_cells_;
};

}

#endif // DDT_DDT_HPP

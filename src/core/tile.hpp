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
#ifndef DDT_TILE_HPP
#define DDT_TILE_HPP

#include <vector>
#include <unordered_set>
#include <set>
#include <map>
#include <list>
#include "bbox.hpp"
#include "ddt_exeptions.hpp"


namespace ddt
{


template<class Map, typename K, typename T>
const T& at_with_default(const Map& map, K k, const T& t)
{
    auto it = map.find(k);
    return it == map.end() ? t : it->second;
}


template<class T>
class Tile
{
public:
    typedef T                                        Traits;
    typedef typename Traits::Id                      Id;
    typedef typename Traits::Point                   Point;
    typedef typename Traits::Flag_V                   Flag_V;
    typedef typename Traits::Flag_C                   Flag_C;
    typedef typename Traits::Data_C                   Data_C;
    typedef typename Traits::Data_C                   Data_V;
    typedef typename Traits::Delaunay_triangulation  DT;
    typedef typename Traits::Vertex_handle           Vertex_handle;
    typedef typename Traits::Vertex_iterator         Vertex_iterator;
    typedef typename Traits::Vertex_const_handle     Vertex_const_handle;
    typedef typename Traits::Vertex_const_iterator   Vertex_const_iterator;
    typedef std::pair<Vertex_const_handle,Id>        Vertex_const_handle_and_id;
    typedef typename Traits::Point_id                      Point_id;
    typedef typename Traits::Point_id_id                   Point_id_id;
    typedef typename Traits::Cell_handle             Cell_handle;
    typedef typename Traits::Cell_const_handle       Cell_const_handle;
    typedef typename Traits::Cell_const_iterator     Cell_const_iterator;
    typedef typename Traits::Facet_const_iterator    Facet_const_iterator;
    typedef typename Traits::Facet_const_iterator    Facet_const_handle;
    typedef std::pair<Cell_const_handle,Id>          Cell_const_handle_and_id;
    enum { D = Traits::D };

    Tile(int id, int dimension = D)
        : id_(id),
          dt_(traits_.triangulation(dimension)),
          number_of_main_vertices_(0),
          number_of_main_facets_(0),
          number_of_main_cells_(0)
    {
        tile_ids =  std::vector<int>(3,0);
    }

    template<class C> Cell_const_handle deref(C c) const { return *c; }
    Cell_const_handle deref(Cell_const_handle c) const { return c; }

    inline DT& tri() { return dt_; }
    inline DT& triangulation() { return dt_; }
    inline const DT& triangulation() const { return dt_; }
    inline const Traits& traits() const { return traits_; }

    inline Id id() const { return id_; }
    inline void set_id(Id i) { id_ = i; }

    inline int maximal_dimension() const { return traits_.maximal_dimension(dt_); }
    inline int current_dimension() const { return traits_.current_dimension(dt_); }




    inline Cell_const_iterator cells_begin() const { return traits_.cells_begin(dt_); }
    inline Cell_const_iterator cells_end  () const { return traits_.cells_end  (dt_); }

    inline Vertex_const_iterator vertices_begin() const { return traits_.vertices_begin(dt_); }
    inline Vertex_const_iterator vertices_end  () const { return traits_.vertices_end  (dt_); }

    inline Vertex_iterator vertices_begin() { return traits_.vertices_begin(dt_); }
    inline Vertex_iterator vertices_end  () { return traits_.vertices_end  (dt_); }

    inline Facet_const_iterator  facets_begin()  const { return traits_.facets_begin(dt_); }
    inline Facet_const_iterator  facets_end  ()  const { return traits_.facets_end  (dt_); }

    inline size_t number_of_vertices() const { return traits_.number_of_vertices(dt_); }
    inline size_t number_of_cells   () const { return traits_.number_of_cells   (dt_); }

    inline size_t number_of_main_vertices() const { return number_of_main_vertices_; }
    inline size_t number_of_main_facets  () const { return number_of_main_facets_;   }
    inline size_t number_of_main_cells   () const { return number_of_main_cells_;    }

    inline bool flag(Vertex_const_handle v,int f) const
    {
        return flagv(v) & static_cast<int>(1 << f);
    }

    inline void flag(Vertex_const_handle v,int f,bool val) const
    {
        Flag_V & fvi = flagv(v);
        val ? (fvi |= static_cast<int>(1 << f)) : (fvi &= ~static_cast<int>(1 << f));
    }

    inline bool flag(Cell_const_handle c,int f) const
    {
        return flagc(c) & static_cast<int>(1 << f);
    }

    inline bool flag(Flag_C  fci,int f) const
    {
        return fci & static_cast<int>(1 << f);
    }


    inline void get_list_vertices(const Cell_const_handle c,std::list<Vertex_const_handle> & ll) const  {return traits_.get_list_vertices  (c,ll);}

    inline void flag(Cell_const_handle c,int f,bool val) const
    {
        Flag_C & fci = flagc(c);
        val ? (fci |= static_cast<int>(1 << f)) : (fci &= ~static_cast<int>(1 << f));
    }


    inline Id   id  (Vertex_const_handle v) const { assert(!vertex_is_infinite(v)); return traits_.id  (v); }
    inline Id   id  (Cell_const_handle c) const { assert(!vertex_is_infinite(c)); return traits_.id  (c); }
    inline Id   gid  (Vertex_const_handle v) const
    {
        assert(!vertex_is_infinite(v));
        return tile_ids[0] + traits_.gid  (v);
    }
    inline Id   gid  (Cell_const_handle c) const { assert(!cell_is_infinite(c)); return tile_ids[2] + traits_.gid  (c); }
    inline Id   lid  (Cell_const_handle c) const { assert(!cell_is_infinite(c)); return  traits_.gid  (c); }
    inline Flag_V& flagv(Vertex_const_handle v) const { assert(!vertex_is_infinite(v)); return traits_.flag(v); }
    inline Flag_C& flagc(Cell_const_handle c) const {  return traits_.flag(c); }
    inline const Data_C& datac(Cell_const_handle c) const  {  return traits_.data(c); }
    inline const Data_V& datav(Vertex_const_handle v) const  {  return traits_.data(v); }





    Id cell_main_id(Cell_const_handle c) const
    {
        Id cid = -1;
        int D = current_dimension();
        for(int i=0; i<=D; ++i)
        {
            auto v = vertex(c, i);
            if(vertex_is_infinite(v)) continue;
            Id vid = id(v);
            if (cid==-1 || vid < cid) cid = vid;
        }
        return cid;
    }

    inline void clear() { traits_.clear(dt_); }
    template<class It> inline void insert(It begin, It end) { traits_.insert(dt_, begin, end); }
    template<class It> inline int insert_splay(It begin, It end, std::map<Id, std::vector<Vertex_const_handle>>& out) { return traits_.insert_splay(dt_, id_, begin, end, out); }
    template<class It, class Out> inline void insert_simplified(It begin, It end, Out out) { traits_.insert_simplified(dt_, id_, begin, end, out); }
    template<class It> inline void remove(It begin, It end) { traits_.remove(dt_, begin, end); }

    inline Vertex_handle infinite_vertex() const { return traits_.infinite_vertex(dt_); }
    inline const Point& point(Vertex_const_handle v) const { return traits_.point(dt_, v); }
    inline double coord(const Point& p, int i) const { return traits_.coord(dt_, p, i); }
    inline double scalar_product(const Point& p, const Point&q) const { return traits_.scalar_product(p,q); }

    inline bool vertex_is_infinite(Vertex_const_handle v) const { return traits_.vertex_is_infinite(dt_, v); }
    inline bool facet_is_infinite (Facet_const_handle  f) const { return traits_.facet_is_infinite (dt_, f); }
    inline bool cell_is_infinite  (Cell_const_handle   c) const { return traits_.cell_is_infinite  (dt_, c); }

    inline std::vector<double > get_cell_barycenter(const Cell_const_handle  c) const { return traits_.get_cell_barycenter  ( c); }

    // Facet functions
    inline int index_of_covertex(Facet_const_handle f) const { return traits_.index_of_covertex(dt_, f); }
    inline Cell_const_handle full_cell(Facet_const_handle f) const { return traits_.full_cell(dt_, f); }

    // Cell functions
    inline Vertex_const_handle vertex(Cell_const_handle c, int i) const { return traits_.vertex(dt_, c, i); }
    inline Facet_const_iterator facet(Cell_const_handle c, int i) const { return traits_.facet(dt_, c, i); }
    inline int mirror_index(Cell_const_handle c, int i) const { return traits_.mirror_index(dt_, c, i); }
    inline Cell_const_handle neighbor(Cell_const_handle c, int i) const { return traits_.neighbor(dt_, c, i); }
    inline bool cell_is_visible(const Bbox<D>& bbox, Cell_const_handle c, int i) const { return traits_.is_visible(dt_, bbox, c, i); }
    inline bool cell_is_touching(const Bbox<D>& bbox, Cell_const_handle c) const { return traits_.is_touching(dt_, bbox, c); }

    inline const Point circumcenter(Cell_const_handle c) const { return traits_.circumcenter(dt_,c); }


    inline void incident_cells(Vertex_const_handle v, std::vector<Cell_const_handle> & cells  ) const { traits_.incident_cells(dt_, v, cells); }
    template<typename Unary_op> inline void finite_adjacent_vertices(Vertex_const_handle v, Unary_op op) const { traits_.finite_adjacent_vertices(dt_, v, op); }


    // Local => all vertex ids are equal to the tile id
    // Mixed  => some vertex ids are equal to the tile id
    // Foreign => no vertex ids are equal to the tile id
    inline bool vertex_is_local(Vertex_const_handle v) const { assert(!vertex_is_infinite(v)); return id(v) == id(); }
    inline bool vertex_is_local_in_tile(Vertex_const_handle v, int tid) const { assert(!vertex_is_infinite(v)); return id(v) == tid; }
    inline bool vertex_is_foreign(Vertex_const_handle v) const { return !vertex_is_local(v); }


    template<typename F>
    int facet_local_score(F f) const
    {
        int icv = index_of_covertex(f);
        auto c = full_cell(f);
        int local = 0;
        for(int i=0; i<=current_dimension(); ++i)
        {
            if (i == icv) continue;
            auto v = vertex(c,i);
            local += vertex_is_local(v);
        }
        return local;
    }

    template<typename F>
    bool facet_is_local(F f) const
    {
        int icv = index_of_covertex(f);
        auto c = full_cell(f);
        for(int i=0; i<=current_dimension(); ++i)
        {
            if (i == icv) continue;
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_foreign(v) ) return false;
        }
        return true;
    }

    template<typename F>
    bool facet_is_mixed(F f) const
    {
        int icv = index_of_covertex(f);
        auto c = full_cell(f);
        bool local_found = false;
        bool foreign_found = false;
        for(int i=0; i <= current_dimension(); ++i)
        {
            if ( i == icv ) continue;
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_local(v) )
            {
                if (foreign_found) return true;
                local_found = true;
            }
            else
            {
                if (local_found) return true;
                foreign_found = true;
            }
        }
        return false;
    }

    template<typename F>
    bool facet_is_foreign(F f) const
    {
        int icv = index_of_covertex(f);
        auto c = full_cell(f);
        for(int i=0; i<=current_dimension(); ++i)
        {
            if ( i == icv ) continue;
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_local(v) ) return false;
        }
        return true;
    }

    template<typename C>
    int cell_local_score(C c) const
    {
        int local = 0;
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            local += vertex_is_local(v);
        }
        return local;
    }


    template<typename C>
    bool cell_is_local(C c) const
    {
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_foreign(v) ) return false;
        }
        return true;
    }

    template<typename C>
    bool cell_is_mixed(C c) const
    {
        bool local_found = false;
        bool foreign_found = false;
        for(int i=0; i <= current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_local(v) )
            {
                if (foreign_found) return true;
                local_found = true;
            }
            else
            {
                if (local_found) return true;
                foreign_found = true;
            }
        }
        return false;
    }

    template<typename C>
    bool cell_is_foreign(C c) const
    {
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_local(v) ) return false;
        }
        return true;
    }


    template<typename C>
    bool cell_is_foreign_in_tile(C c, int tid) const
    {
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if ( vertex_is_infinite(v) ) continue;
            if ( vertex_is_local_in_tile(v,tid) ) return false;
        }
        return true;
    }

    // Main
    template<typename V>
    bool vertex_is_main(V v) const
    {
        // TODO: define somehow the main infinite vertex
        return !vertex_is_infinite(v) && vertex_is_local(v) ;
    }





    // Check if the facet is main
    template<typename F>
    bool facet_is_main(F f) const
    {
        // return true;
        int icv = index_of_covertex(f);
        auto c1 = full_cell(f);
        auto c2 = neighbor(c1,icv);
        std::vector<int> lid;
        for(int i=0; i<=current_dimension(); ++i)
        {
            if (i == icv) continue;
            auto v = vertex(c1,i);
            if (vertex_is_infinite(v)) continue;
            Id vid = id(v);
            lid.push_back(vid);
        }
        std::sort(lid.begin(),lid.end());
        for(int i=0; i<lid.size(); ++i)
        {
            int vid = lid[i];
            if(!cell_is_foreign_in_tile(c1,vid) && !cell_is_foreign_in_tile(c2,vid))
                return vid == id();
        }
        return false;
    }




    template<typename C>
    bool cell_is_main(C c) const
    {
        bool foreign = true;
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if (vertex_is_infinite(v))
            {
                continue;
            }
            Id vid = id(v);
            if ( vid < id() )
                return false;
            else if (vid == id())
                foreign = false;
        }
        return !foreign;
    }


    // Check if the facet is main
    template<typename F>
    bool facet_has_id(F f,Id ii) const
    {
        int icv = index_of_covertex(f);
        auto c = full_cell(f);
        bool foreign = true;
        for(int i=0; i<=current_dimension(); ++i)
        {
            if (i == icv) continue;
            auto v = vertex(c,i);
            if (vertex_is_infinite(v)) continue;
            Id vid = id(v);
            if(vid == ii)
                return true;
        }
        return false;
    }


    template<typename C>
    bool has_id(C c,Id ii) const
    {
        bool foreign = true;
        for(int i=0; i<=current_dimension(); ++i)
        {
            auto v = vertex(c,i);
            if (vertex_is_infinite(v))
            {
                continue;
            }
            Id vid = id(v);
            if ( vid == ii )
                return true;
        }
        return false;
    }


    // Active
    template<typename C>
    bool cell_is_active(const Bbox<D>& bbox, C c) const
    {
        for(int i=0; i<=current_dimension(); ++i)
            if (vertex_is_infinite(vertex(c,i)))
                return cell_is_visible(bbox, c, i);
        return cell_is_touching(bbox, c);
    }


    // Active
    template<typename C>
    bool cell_is_inside(Bbox<D> &  cur_bb, C c) const
    {
        return traits_.is_inside(dt_,cur_bb, c);
    }

    // Active
    template<typename C>
    bool cell_is_inside(C c) const
    {
        return traits_.is_inside(dt_,bbox(id_), c);
    }



    int simplify()
    {
        return simplify_range2(vertices_begin(), vertices_end(), cells_begin(), cells_end());
    }

    template<class VertexIterator, class CellIterator>
    int simplify_range2(VertexIterator vbegin, VertexIterator vend, CellIterator cbegin, CellIterator cend)
    {
        // initialize flags to 1
        for(auto vit = vbegin; vit != vend; ++vit)
            if(!vertex_is_infinite(vit))
                flag(vit,0,1);
        // set flags of vertices incident to non-foreign cells to 0
        for(auto cit = cbegin; cit != cend; ++cit)
        {
            if(cell_is_foreign(cit)) continue;
            for(int i=0; i<=current_dimension(); ++i)
            {
                Vertex_const_handle v = vertex(cit, i);
                if(!vertex_is_infinite(v))
                {
                    flag(v,0,0);
                }
            }
        }
        // gather vertices that are to be removed
        std::vector<Vertex_handle> todo;
        for(auto vit = vbegin; vit != vend; ++vit)
            if(!vertex_is_infinite(vit) && flag(vit,0))
                todo.push_back(vit);
        // remove these vertices
        remove(todo.begin(), todo.end());
        return todo.size();
    }

    template<class Iterator>
    void simplify_range(Iterator begin, Iterator end, std::unordered_set<Vertex_handle>& removed)
    {
        // get vertices adjacent to new inserted foreign vertices (in the range)
        std::unordered_set<Vertex_handle> neighbors;
        for(auto it = begin; it != end; ++it)
            if (vertex_is_foreign(*it))
                finite_adjacent_vertices(*it, [&neighbors](Vertex_handle v) { neighbors.insert(v); });
        // remove if they have a foreign star
        std::vector<Vertex_handle> vertices;
        vertices.reserve(32);
        for(auto it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            if (vertex_is_local(*it)) continue;
            vertices.clear();
            finite_adjacent_vertices(*it, [&vertices](Vertex_handle v) { vertices.push_back(v); });
            auto vit = vertices.begin();
            for(; vit != vertices.end(); ++vit)
                if (vertex_is_local(*vit))
                    break;
            if (vit == vertices.end())
                removed.insert(*it);
        }
        remove(removed.begin(), removed.end());
    }

    template<typename Iterator>
    void get_neighbors_range(Iterator begin, Iterator end, std::map<Id, std::unordered_set<Vertex_const_handle>>& outbox) const
    {
        for(auto cit = begin; cit != end; ++cit)
            for(int i=0; i<=current_dimension(); ++i)
            {
                Cell_const_handle c = deref(cit);
                Vertex_const_handle v = vertex(c, i);
                if(vertex_is_infinite(v)) continue;
                Id idv = id(v);
                if(idv != id()) // v is foreign
                    for(int j=0; j<=current_dimension(); ++j)
                    {
                        Vertex_const_handle w = vertex(c, j);
                        if(!vertex_is_infinite(w) && id(w) != idv)  // implies i!=j
                        {
                            outbox[idv].insert(w);
                        }
                    }
            }
    }

    int send_one(
        std::map<Id, std::vector<Point_id>>& inbox,
        const std::vector<Vertex_const_handle_and_id>& outbox
    )
    {
        int count = 0;
        for(auto& vi : outbox)
            count += send_vertex(inbox, vi.first, id(), vi.second);
        return count;
    }

    template<class Iterator>
    int send_one(
        std::map<Id, std::vector<Point_id>>& inbox,
        const std::map<Id, std::vector<Vertex_const_handle>>& outbox
    )
    {
        int count = 0;
        for(auto& pair : outbox)
            for(auto& vi : pair.second)
                count += send_vertex(inbox, vi, id(), pair.first);
        return count;
    }

    template<class Iterator, class Out>
    int send_one(Id target, Iterator begin, Iterator end, Out out)
    {
        int count = 0;
        for(Iterator it = begin; it != end; ++it)
        {
            Vertex_const_handle vh = *it;
            if(target != id() && !vertex_is_infinite(vh))
            {
                Id vid = id(vh);
                if(target != vid && sent_[target].insert(vh).second)
                {
                    *out++ = std::make_pair(point(vh),vid);
                    ++count;
                }
            }
        }
        return count;
    }

    template<typename Id_iterator>
    int send_all(
        std::map<Id, std::vector<Point_id>>& inbox,
        const std::vector<Vertex_const_handle>& outbox,
        Id_iterator begin,
        Id_iterator end
    ) const
    {
        int count = 0;
        for(auto v : outbox)
            for(Id_iterator target = begin; target != end; ++target)
                count += send_vertex(inbox, v, id(), *target);
        return count;
    }


    inline bool send_vertex(
        std::map<Id, std::vector<Point_id>>& inbox, Vertex_const_handle vh, Id source, Id target) const
    {
        if(target == source || vertex_is_infinite(vh))
            return false;
        Id vid = id(vh);
        if(target==vid || !sent_[target].insert(vh).second)
            return false;
        inbox[target].emplace_back(point(vh),vid);
        return true;
    }

    int insert(const std::vector<Point_id>& inbox, bool do_simplify = true)
    {
        if(inbox.empty()) return 0;
        for(auto& v : inbox)
        {
            bbox_[v.second] += v.first;
        }
        insert(inbox.begin(), inbox.end());
        int s = 0;
        if(do_simplify)
            s = simplify();
        return inbox.size() - s;
    }


    // Algo
    int insert_points_id_id(std::vector<Point_id_id> & vpis, int tid, bool do_insert_local = false)
    {
        std::vector<Point_id> v_pai;
        int nbi1 = 0;
        bool do_simplify=true;
        for(auto & pp : vpis)
        {
            Id idp =  std::get<1>(pp);
            Id id_source =  std::get<2>(pp);
            if( (id_source==tid && !do_insert_local) )
            {
                continue;
            }
            else
            {
                v_pai.emplace_back(std::get<0>(pp), idp);
            }
        }
        if(v_pai.size()> 0)
        {
            nbi1 += insert(v_pai,do_simplify);
        }
        return nbi1;
    }


    int insert_simplified(const std::vector<Point_id>& inbox, std::vector<Vertex_handle>& inserted)
    {
        if(inbox.empty()) return 0;
        for(auto& v : inbox)
        {
            bbox_[v.second] += v.first;
        }
        insert_simplified(inbox.begin(), inbox.end(), std::back_inserter(inserted));
        return inserted.size();
    }

    int insert_splay(const std::vector<Point_id>& inbox,
                     std::map<Id, std::vector<Vertex_const_handle>>& outbox, bool do_simplify = true)
    {
        if(inbox.empty()) return 0;
        for(auto& v : inbox)
        {
            bbox_[v.second] += v.first;
        }
        int s = insert_splay(inbox.begin(), inbox.end(), outbox);
        if(do_simplify)
            s -= simplify();
        return s;
    }

    void get_mixed_cells(std::vector<Cell_const_handle_and_id>& out) const
    {
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
            if(cell_is_mixed(cit))
                for(int i=0; i<=current_dimension(); ++i)
                    if(!vertex_is_infinite(vertex(cit,i)))
                        out.emplace_back(cit, id(vertex(cit,i)));
    }

    void get_adjacency_graph_edges(std::set<Id>& out_edges) const
    {
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
            if(cell_is_mixed(cit))
                for(int i=0; i<=current_dimension(); ++i)
                {
                    Vertex_const_handle v = vertex(cit,i);
                    if(!vertex_is_infinite(v) && id(v) != id())
                        out_edges.insert(id(v));
                }
    }

    void get_merge_graph_edges(std::set<Id>& out_edges, std::vector<Cell_const_handle>& finalized) const
    {
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
        {
            if(cell_is_foreign(cit)) continue;
            bool active = false;
            for(auto pair : bbox_)
                if((active || out_edges.find(pair.first)==out_edges.end()) && cell_is_active(pair.second, cit) )
                {
                    active = true;
                    out_edges.insert(pair.first);
                }
            if(!active)
                finalized.push_back(cit);
        }
    }


    bool cell_is_finalized(Cell_const_handle c) const
    {
        if(cell_is_foreign(c)) return false;
        // TODO: acceleration data structure !!!
        for(auto pair : bbox_)
            if(pair.first != id() && cell_is_active(pair.second, c)) return false;
        return true;
    }

    const std::map<Id, Bbox<D>>& bbox() const { return bbox_; }
    std::map<Id, Bbox<D>>& bbox()     { return bbox_; }
    const Bbox<D>& bbox(Id ii) const { return bbox_.at(ii); }

    Vertex_const_handle locate_vertex(const Tile& t, Vertex_const_handle v) const
    {
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit )
        {
            if(traits_.are_vertices_equal(t.dt_, v, dt_, vit))
                return vit;
        }
        assert(false);
        return vertices_end();
    }

    Facet_const_handle locate_facet(const Tile& t, Facet_const_handle f) const
    {
        assert(!t.cell_is_foreign(t.full_cell(f)));
        Cell_const_handle c = full_cell(f);
        Cell_const_handle d = locate_cell(t, c);
        Vertex_const_handle v = vertex(c, index_of_covertex(f));
        for(int i=0; i<=current_dimension(); ++i)
        {
            if(traits_.are_vertices_equal(t.dt_, v, dt_, vertex(d, i)))
                return facet(d, i);
        }
        assert(false);
        return facets_end();
    }

    Cell_const_handle locate_cell_point(const Tile& t, Point & pp ) const
    {
        auto cc = traits_.locate_cell_point(dt_,pp);
        return cc;
    }

    // locate the cell by searching the barycenter with the local function on the oposite delaunay triangulation
    // It it does not work, loop on all the cell
    Cell_const_handle locate_cell(const Tile& t, Cell_const_handle c) const
    {
        assert(!t.cell_is_foreign(c));
        auto cc = traits_.locate_cell(dt_,c);
        if(traits_.are_cells_equal(t.dt_, c, dt_, cc))
            return cc;
        else
            return locate_cell_slow(t,c);
    }

    // slow localisation on the
    Cell_const_handle locate_cell_slow(const Tile& t, Cell_const_handle c) const
    {
        assert(!t.cell_is_foreign(c));
        for(auto cit = cells_begin(); cit != cells_end(); ++cit )
        {
            if(traits_.are_cells_equal(t.dt_, c, dt_, cit))
                return cit;
        }
        std::cerr << "[Warning] locate cell not found" << std::endl;
        throw DDT_exeption("ex_localte_cell_not_found");
        assert(false);
        return cells_end();
    }

    size_t finalize_vertices()
    {
        number_of_main_vertices_ = 0;
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit)
            if(vertex_is_main(vit))
                ++number_of_main_vertices_;
        return number_of_main_vertices_;
    }

    size_t finalize_facets()
    {
        number_of_main_facets_ = 0;
        for(auto fit = facets_begin(); fit != facets_end(); ++fit )
            if(facet_is_main(fit))
                ++number_of_main_facets_;
        return number_of_main_facets_;
    }

    size_t finalize_cells()
    {
        number_of_main_cells_ = 0;
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
            if(cell_is_main(cit))
                ++number_of_main_cells_;
        return number_of_main_cells_;
    }

    void finalize()
    {
        finalize_vertices();
        finalize_facets();
        finalize_cells();
    }

    bool is_valid() const
    {
        return dt_.is_valid();
    }


    void init_local_id_tile()
    {
        int D = Traits::D;
        int acc = 0;
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit)
        {
            if(!vertex_is_main(vit))
                continue;
            const Data_V & vd = datav(vit);
            Data_V & vd_quickndirty = const_cast<Data_V &>(vd);
            vd_quickndirty.gid =  (acc++);
        }
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit)
        {
            if(vertex_is_main(vit))
                continue;
            const Data_V & vd = datav(vit);
            Data_V & vd_quickndirty = const_cast<Data_V &>(vd);
            vd_quickndirty.gid =  (acc++);
        }
        acc = 0;
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
        {
            if(!cell_is_main(cit)) continue;
            const Data_C & cd = datac(cit);
            Data_C & cd_quickndirty = const_cast<Data_C &>(cd);
            cd_quickndirty.gid =  (acc++);
        }
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
        {
            if(cell_is_main(cit)) continue;
            const Data_C & cd = datac(cit);
            Data_C & cd_quickndirty = const_cast<Data_C &>(cd);
            cd_quickndirty.gid =  (acc++);
        }
    }

    void update_local_flag() const
    {
        int D = Traits::D;
        int acc = 0;
        for(auto vit = vertices_begin(); vit != vertices_end(); ++vit)
        {
            if(!vertex_is_main(vit))
                flagv(vit) = vertex_is_local(vit);
        }
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
        {
            if(!cell_is_main(cit)) continue;
            int local = 0;
            for(int i=0; i<=D+1; ++i)
            {
                auto v = cit->vertex(i % (D+1));
                if(i>0)
                {
                    local += vertex_is_local(v);
                }
            }
            flagc(cit) = local;
        }
    }


    template<class Map>
    void sanity_check_send(Map& map) const
    {
        for(auto cit = cells_begin(); cit != cells_end(); ++cit)
        {
            std::vector<Id> ids;
            bool foreign = true;
            for(int i=0; i<=current_dimension(); ++i)
            {
                Vertex_const_handle v = vertex(cit, i);
                if (vertex_is_infinite(v)) continue;
                Id vid = id(v);
                ids.push_back(vid);
                if (vid == id()) foreign = false;
            }
            if (foreign) continue;
            std::sort(ids.begin(), ids.end());
            for(auto i: ids) ++map[i][ids];
        }
    }

    template<class Map>
    bool sanity_check_recv(const Map& map) const
    {
        Id id1 = id();
        std::map<std::vector<Id>,size_t> empty_map;
        // self-assessed counts of mixed and local cells of this tile
        auto& c1 = at_with_default(map, id1, empty_map);
        for(const auto& pair: map)
        {
            Id id2 = pair.first;
            if(id1==id2) continue;
            auto& c2 = pair.second;
            // check that maps have the same counts on keys that reference both id1 and id2
            auto it1 = c1.begin();
            auto it2 = c2.begin();
            while(it1 != c1.end() && it2 != c2.end())
            {
                if (!std::binary_search(it1->first.begin(), it1->first.end(), id2)) { ++it1; continue; }
                if (!std::binary_search(it2->first.begin(), it2->first.end(), id1)) { ++it2; continue; }
                if (!(it1->first == it2->first && it1->second == it2->second))
                    return false;
                ++it1;
                ++it2;
            }
            // check left overs
            while(it1 != c1.end()) if (std::binary_search(it1->first.begin(), it1->first.end(), id2)) return false;
                else ++it1;
            while(it2 != c2.end()) if (std::binary_search(it2->first.begin(), it2->first.end(), id1)) return false;
                else ++it2;
        }
        return true;
    }

    std::map<Id, std::set<Point>> points_sent_;
    std::vector<int> tile_ids;
private:
    Traits traits_;
    Id id_;
    DT dt_;

    size_t number_of_main_vertices_;
    size_t number_of_main_facets_;
    size_t number_of_main_cells_;


    std::map<Id, Bbox<D>> bbox_;
    mutable std::map<Id, std::unordered_set<Vertex_const_handle>> sent_;
};

}

#endif // DDT_TILE_HPP

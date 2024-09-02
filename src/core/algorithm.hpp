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
#ifndef DDT_ALGORITHM_HPP
#define DDT_ALGORITHM_HPP

#include <vector>
#include <unordered_set>
#include <cmath>
#include "bbox.hpp"

namespace ddt
{

float sqrt3(const float x)
{
    union
    {
        int i;
        float x;
    } u;
    u.x = x;
    u.i = (1<<29) + (u.i >> 1) - (1<<22);
    return u.x;
}


double doubleRand()
{
    return double(std::rand()) / (double(RAND_MAX) + 1.0);
}



template <typename DTC,typename CHR,typename TR>
bool is_inside_bbox(DTC & tri,CHR cit,  ddt::Bbox<TR::D> & tri_bbox, TR ttr)
{
    if(tri.is_infinite(cit))
    {
        return false;
    }
    else
    {
        auto point  = cit->vertex(0)->point();
        auto center = ttr.circumcenter(tri,cit);
        double dist = 0;
        int D  = TR::D;
        for(int d = 0; d < D; d++)
            dist += (point[d] - center[d])*(point[d] - center[d]);
        dist = sqrt3(dist);
        for(int d = 0; d < D; d++)
        {
            if(dist  >= fabs(center[d] - tri_bbox.max(d)) || dist >= fabs(center[d] - tri_bbox.min(d)))
                return false;
        }
        return true;
    }
}


// filter the cell, if cells are inside or outside the bounding box.
template <typename T>
struct filter_cell
{
    filter_cell(ddt::Bbox<T::D> bb) : tri_bbox(bb) {}
    template <typename TTr,typename DTC,typename CHR>
    bool do_keep(DTC & tri,CHR cit,TTr & traits)
    {
        return traits.is_inside(tri,tri_bbox,cit);
    }
    ddt::Bbox<T::D> tri_bbox;
};

// the same filter as above
template <typename T>
struct filter_cell_ddt
{
    filter_cell_ddt(ddt::Bbox<T::D> bb, int ii) : tri_bbox(bb),tid(ii) {}
    template <typename TTr,typename DTC,typename CHR>
    bool do_keep(DTC & tri,CHR cit,TTr & traits)
    {
        int local_score = 0;
        bool is_main = true;
        if(tri.is_infinite(cit))
            return false;
        for(int d = 0; d < (T::D) +1; d++)
        {
            int pid = traits.id(cit->vertex(d));
            if(pid < tid)
                is_main = false;
            if(pid == tid)
                local_score++;
        }
        return (!traits.is_inside(tri,tri_bbox,cit) && local_score != 0 && is_main);
    }
    ddt::Bbox<T::D> tri_bbox;
    int tid;
};



template<typename DT, typename OutputIterator,typename TR>
OutputIterator get_bbox_points_raw(DT& tri, OutputIterator out,TR traits)
{
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Point Point;
    const int D = TR::D;
    Vertex_handle v[2*D];
    auto vit = tri.vertices_begin();
    // first local point
    for(; vit != tri.vertices_end(); ++vit)
    {
        if (!tri.is_infinite(vit) )
        {
            for(int i=0; i<2*D; ++i) v[i] = vit;
            break;
        }
    }
    if(vit == tri.vertices_end()) return out; // no local points
    // other local points
    for(; vit != tri.vertices_end(); ++vit)
    {
        if (!tri.is_infinite(vit))
        {
            const Point& p = vit->point();
            for(int i=0; i<D; ++i)
            {
                if(p[i] < v[i  ]->point()[i]) v[i  ] = vit;
                if(p[i] > v[i+D]->point()[i]) v[i+D] = vit;
            }
        }
    }
    // remove duplicates (O(D^2) implementation, should we bother ?)
    for(int i=0; i<2*D; ++i)
    {
        int j = 0;
        for(; j<i; ++j) if(v[j]==v[i]) break;
        if(i==j) *out++ = v[i];
    }
    return out;
}


struct get_bbox_points
{
    template<typename Tile, typename OutputIterator>
    OutputIterator operator()(const Tile& tile, OutputIterator out) const
    {
        typedef typename Tile::Vertex_const_handle Vertex_const_handle;
        typedef typename Tile::Point Point;
        const int D = tile.maximal_dimension();
        Vertex_const_handle v[2*D];
        auto vit = tile.vertices_begin();
        // first local point
        for(; vit != tile.vertices_end(); ++vit)
        {
            if (!tile.vertex_is_infinite(vit) && tile.vertex_is_local(vit))
            {
                for(int i=0; i<2*D; ++i) v[i] = vit;
                break;
            }
        }
        if(vit == tile.vertices_end()) return out; // no local points
        // other local points
        for(; vit != tile.vertices_end(); ++vit)
        {
            if (!tile.vertex_is_infinite(vit) && tile.vertex_is_local(vit))
            {
                const Point& p = tile.point(vit);
                for(int i=0; i<D; ++i)
                {
                    if(p[i] < tile.point(v[i  ])[i]) v[i  ] = vit;
                    if(p[i] > tile.point(v[i+D])[i]) v[i+D] = vit;
                }
            }
        }
        // remove duplicates (O(D^2) implementation, should we bother ?)
        for(int i=0; i<2*D; ++i)
        {
            int j = 0;
            for(; j<i; ++j) if(v[j]==v[i]) break;
            if(i==j) *out++ = v[i];
        }
        return out;
    }
};



struct get_neighbors
{
    template<typename Tile, typename OutputIterator>
    inline OutputIterator operator()(const Tile& tile, OutputIterator out) const
    {
        typedef typename Tile::Vertex_const_handle Vertex_const_handle;
        typedef typename Tile::Id Id;
        std::map<Id, std::unordered_set<Vertex_const_handle>> outbox;
        tile.get_neighbors_range(tile.cells_begin(), tile.cells_end(), outbox);
        for(auto&& pair : outbox)
            for(auto vh : pair.second)
                *out++ = std::make_pair(vh, pair.first);
        return out;
    }
};



}

#endif // DDT_ALGORITHM_HPP

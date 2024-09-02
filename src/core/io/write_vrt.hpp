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
#ifndef DDT_WRITE_VRT_HPP
#define DDT_WRITE_VRT_HPP

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "conf_header/conf.hpp"
#include "algorithm.hpp"

namespace ddt
{


void add_qgis_style(const std::string& filename,const std::string& stylename);
void write_vrt_header_vert (std::ofstream& csv, const std::string& filename);
void write_vrt_header_facet(std::ofstream& csv, const std::string& filename);
void write_vrt_header_cell (std::ofstream& csv, const std::string& filename);
void write_vrt_header_tin  (std::ofstream& csv, const std::string& filename);
void write_vrt_header_bbox (std::ofstream& csv, const std::string& filename);


template<typename Tile>
void write_csv_vert(const Tile& tile, std::ostream& csv, bool main_only=false)
{
    int D = tile.maximal_dimension();
    for(auto vit = tile.vertices_begin(); vit != tile.vertices_end(); ++vit)
    {
        if(tile.vertex_is_infinite(vit)) continue;
        if(main_only && !tile.vertex_is_main(vit)) continue;
        csv << "POINT( ";
        for(int d=0; d<D; ++d)
            csv << tile.point(vit)[d] << " ";
        csv << ")," << int(tile.id()) << "," << int(tile.id(vit)) << "\n";
    }
}

template<typename Tile>
void write_csv_facet(const Tile& tile, std::ostream& csv, bool main_only=false)
{
    int D = tile.maximal_dimension();
    for(auto fit = tile.facets_begin(); fit != tile.facets_end(); ++fit)
    {
        if(tile.facet_is_infinite(fit)) continue;
        if(main_only && !tile.facet_is_main(fit)) continue;
        auto cit = tile.full_cell(fit);
        int idx = tile.index_of_covertex(fit);
        csv << "\"LINESTRING(";
        int local = 0;
        int j = 0;
        for(int i=0; i<=D; ++i)
        {
            if(i == idx) continue;
            auto v = tile.vertex(cit,i);
            local += tile.vertex_is_local(v);
            for(int d=0; d<D; ++d)
                csv << tile.point(v)[d] << " ";
            if (++j < D) csv << ",";
        }
        csv << ")\"," << int(tile.id()) << "," << local << "\n";
    }
}

template<typename Tile>
void write_csv_cell(const Tile& tile, std::ostream& csv, bool main_only=false)
{
    int D = tile.maximal_dimension();
    for(auto cit = tile.cells_begin(); cit != tile.cells_end(); ++cit)
    {
        if(tile.cell_is_infinite(cit)) continue;
        if(main_only && !tile.cell_is_main(cit))continue;
        csv << "\"POLYGON((";
        int local = 0;
        for(int i=0; i<=D; ++i)
        {
            const typename Tile::Vertex_const_handle v = tile.vertex(cit,i);
            local += tile.vertex_is_local(v);
            for(int d=0; d<D; ++d)
                csv << tile.point(v)[d] << " ";
            csv << ",";
        }
        for(int d=0; d<D; ++d) // repeat first to close the polygon
            csv << tile.point(tile.vertex(cit, 0))[d] << " ";
        csv << "))\"," << int(tile.id()) << "," << int(local) << "," << int(tile.cell_is_main(cit)) << "\n";
    }
}

template<typename Tile>
void write_tile_vrt_vertices(const Tile& tile, const std::string& vrt_name, bool main_only=false)
{
    std::ofstream csv;
    write_vrt_header_vert(csv,vrt_name);
    add_qgis_style(vrt_name,std::string("tri1.qml"));
    write_csv_vert(tile, csv,main_only);
}

template<typename Tile>
void write_tile_vrt_cells(const Tile& tile, const std::string& vrt_name, bool main_only=false)
{
    std::ofstream csv;
    write_vrt_header_cell(csv,vrt_name);
    add_qgis_style(vrt_name,std::string("tri1.qml"));
    write_csv_cell(tile, csv,main_only);
}

template<typename Tile>
void write_csv_tin(const Tile& tile, std::ostream& csv, bool main_only=false)
{
    int D = tile.maximal_dimension();
    csv << "\"TIN (";
    bool first = true;
    for(auto cit = tile.cells_begin(); cit != tile.cells_end(); ++cit)
    {
        if(tile.cell_is_infinite(cit)) continue;
        if(main_only && !tile.cell_is_main(cit))continue;
        if(!first) csv << ", ";
        first = false;
        csv << "((";
        typename Tile::Vertex_const_handle v;
        for(int i=0; i<=D; ++i)
        {
            v = tile.vertex(cit,i);
            for(int d=0; d<D; ++d)
                csv << tile.point(v)[d] << " ";
            csv /*<< int(tile.id(v))*/ << ",";
        }
        v = tile.vertex(cit,0);
        for(int d=0; d<D; ++d) // repeat first to close the polygon
            csv << tile.point(v)[d] << " ";
        csv /*<< int(tile.id(v))*/ << "))";
    }
    csv << ")\"," << int(tile.id()) << "\n";
}

template<typename Tile>
void write_csv_bbox(const Tile& tile, std::ostream& csv)
{
    for(auto& pair : tile.bbox())
    {
        auto bboxid = pair.first;
        auto bbox   = pair.second;
        csv << "\"POLYGON((";
        csv << bbox.min(0) << " "<< bbox.min(1) << ", ";
        csv << bbox.max(0) << " "<< bbox.min(1) << ", ";
        csv << bbox.max(0) << " "<< bbox.max(1) << ", ";
        csv << bbox.min(0) << " "<< bbox.max(1) << ", ";
        csv << bbox.min(0) << " "<< bbox.min(1);
        csv << "))\"," << int(tile.id()) << "," << int(bboxid) << "\n";
    }
}

template<typename Tile>
void write_csv_bbox_vert(const Tile& tile, std::ostream& csv)
{
    int D = tile.maximal_dimension();
    std::vector<typename Tile::Vertex_const_handle> points;
    get_bbox_points()(tile, std::back_inserter(points));
    for(auto it : points)
    {
        csv << "POINT( ";
        for(int d=0; d<D; ++d)
            csv << tile.point(it)[d] << " ";
        csv << ")," << int(tile.id()) << "," << int(tile.id(it)) << "\n";
    }
}


template<typename Iterator>
void write_vrt_vert_range(Iterator begin, Iterator end, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_vert(csv, filename);
    typedef typename Iterator::value_type Vertex_const_iterator;
    typedef typename Vertex_const_iterator::Traits Traits;
    int D = Traits::D;
    for(auto vit = begin; vit != end; ++vit)
    {
        if(vit->is_infinite()) continue;
        csv << "POINT( ";
        for(int d=0; d<D; ++d)
            csv << vit->point()[d] << " ";
        csv << ")," << int(vit->tile()->id()) << "," << int(vit->main_id()) << "\n";
    }
}

template<typename Iterator>
void write_vrt_facet_range(Iterator begin, Iterator end, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_facet(csv, filename);
    typedef typename Iterator::value_type Facet_const_iterator;
    typedef typename Facet_const_iterator::Traits Traits;
    int D = Traits::D;
    for(auto fit = begin; fit != end; ++fit)
    {
        if(fit->is_infinite()) continue;
        auto cit = fit->full_cell();
        int idx = fit->index_of_covertex();
        csv << "\"LINESTRING(";
        int local = 0;
        int j = 0;
        for(int i=0; i<=D; ++i)
        {
            if(i == idx) continue;
            auto v = cit->vertex(i);
            local += v->is_local();
            for(int d=0; d<D; ++d)
                csv << v->point()[d] << " ";
            if (++j < D) csv << ",";
        }
        csv << ")\"," << int(fit->tile()->id()) << "," << local << "\n";
    }
}

template<typename Iterator>
void write_vrt_cell_range(Iterator begin, Iterator end, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_cell(csv, filename);
    typedef typename Iterator::value_type Cell_const_iterator;
    typedef typename Cell_const_iterator::Traits Traits;
    std::map<Cell_const_iterator, int> cmap;
    int nextid = 0;
    int D = Traits::D;
    for(auto iit = begin; iit != end; ++iit)
    {
        if(iit->is_infinite()) continue;
        csv << "\"POLYGON((";
        int local = 0;
        for(int i=0; i<=D+1; ++i) // repeat first to close the polygon
        {
            auto v = iit->vertex(i % (D+1));
            if(i>0)
            {
                csv << ",";
                local += v->is_local();
            }
            auto p = v->point();
            for(int d=0; d<D; ++d) csv << p[d] << " ";
        }
        csv << "))\"," << int(iit->tile()->id()) << "," << int(local) << "," << int(iit->main_id());
        auto n0 = iit->neighbor(0)->main();
        auto n1 = iit->neighbor(1)->main();
        auto n2 = iit->neighbor(2)->main();
        if(!cmap.count(*iit)) cmap[*iit] = nextid++;
        if(!cmap.count(n0)) cmap[n0] = nextid++;
        if(!cmap.count(n1)) cmap[n1] = nextid++;
        if(!cmap.count(n2)) cmap[n2] = nextid++;
        csv << "," << cmap[*iit] << "," << cmap[n0] << "," << cmap[n1] << "," << cmap[n2];
        csv << "," << 0 << "," << 0 << "," << 0 << "," << 0;
        csv << "\n";
    }
}

// VRT+CSV writers (DDT)

template<typename DDT>
void write_vrt_vert(const DDT& ddt, const std::string& filename)
{
    write_vrt_vert_range(ddt.vertices_begin(), ddt.vertices_end(), filename);
}

template<typename DDT>
void write_vrt_facet(const DDT& ddt, const std::string& filename)
{
    write_vrt_facet_range(ddt.facets_begin(), ddt.facets_end(), filename);
}

template<typename DDT>
void write_vrt_cell(const DDT& ddt, const std::string& filename)
{
    write_vrt_cell_range(ddt.cells_begin(), ddt.cells_end(), filename);
}

template<typename DDT>
void write_vrt_tin(const DDT& tri, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_tin(csv, filename);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_csv_tin(*tile, csv);
}

template<typename DDT>
void write_vrt_bbox(const DDT& tri, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_bbox(csv, filename);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_csv_bbox(*tile, csv);
}

template<typename DDT>
void write_vrt_bbox_vert(const DDT& tri, const std::string& filename)
{
    std::ofstream csv;
    write_vrt_header_vert(csv, filename);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_csv_bbox_vert(*tile, csv);
}


template<typename DDT>
void write_vrt_verts(const DDT& tri, const std::string& dirname)
{
    if(boost::filesystem::exists(dirname))
    {
        for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        {
            std::string filename(dirname + std::string("/tile_verts") + std::to_string(tile->id()) + ".vrt");
            std::ofstream csv;
            write_vrt_header_vert(csv, filename);
            write_csv_vert(*tile, csv, true);
        }
    }
    else
    {
        std::cerr << "[ERROR]" << dirname << " does not exist, create it before writing" << std::endl;
    }
}

template<typename DDT>
void write_vrt_facets(const DDT& tri, const std::string& dirname)
{
    if(boost::filesystem::exists(dirname))
    {
        for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        {
            std::string filename(dirname + std::string("/tile_facets") + std::to_string(tile->id()) + ".vrt");
            std::ofstream csv;
            write_vrt_header_facet(csv, filename);
            write_csv_facet(*tile, csv, true);
        }
    }
    else
    {
        std::cerr << "[ERROR]" << dirname << " does not exist, create it before writing" << std::endl;
    }
}

template<typename DDT>
void write_vrt_cells(const DDT& tri, const std::string& dirname)
{
    if(boost::filesystem::exists(dirname))
    {
        for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        {
            std::string filename(dirname + std::string("/tile_cells") + std::to_string(tile->id()) + ".vrt");
            std::ofstream csv;
            write_vrt_header_cell(csv, filename);
            write_csv_cell(*tile, csv, true);
        }
    }
    else
    {
        std::cerr << "[ERROR]" << dirname << " does not exist, create it before writing" << std::endl;
    }
}

}

#endif // DDT_WRITE_VRT_HPP

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
#ifndef DDT_WRITE_PLY_HPP
#define DDT_WRITE_PLY_HPP

#include <string>
#include <fstream>
#include <map>
#include <CGAL/Random.h>
#include <CGAL/Unique_hash_map.h>

#include "number_parser.hpp"
namespace ddt
{

void write_ply_header_begin(std::ostream& out);
void write_ply_header_end(std::ostream& out);



template<typename DDT>
void write_ply_element_cell(const DDT& tri, std::ostream& out)
{
    int nc = tri.number_of_cells();
    out << "element face " << nc << std::endl;
    out << "property list uint8 int vertex_indices" << std::endl;
    out << "property uint8 tile" << std::endl;
    out << "property uint8 local" << std::endl;
}

template<typename DDT>
void write_ply_element_vert(const DDT& tri, std::ostream& out)
{
    int D = tri.maximal_dimension();
    int nv = tri.number_of_vertices();
    out << "element vertex " << nv << std::endl;
    out << "property float32 x" << std::endl;
    out << "property float32 y" << std::endl;
    if(D>2) out << "property float32 z" << std::endl;
    out << "property uint8 tile" << std::endl;
    out << "property uint8 id" << std::endl;
}

template<typename Tile>
void write_ply_property_cell(const Tile& tile, std::ostream& out)
{
    unsigned char N = (unsigned char)(tile.maximal_dimension()+1);
    unsigned char tid = tile.id();
    std::map<typename Tile::Vertex_const_handle, int> dict;
    int id = 0;
    for(auto it = tile.vertices_begin(); it != tile.vertices_end(); ++it)
        if(!tile.vertex_is_infinite(it)) dict[it] = id++;
    for(auto cit = tile.cells_begin(); cit != tile.cells_end(); ++cit)
    {
        if(tile.cell_is_infinite(cit)) continue;
        unsigned char local = 0;
        out.write(reinterpret_cast<char *>(&N), sizeof(N));
        for(int i=0; i<N; ++i)
        {
            typename Tile::Vertex_const_handle v = tile.vertex(cit, i);
            int id = dict[v];
            out.write(reinterpret_cast<char *>(&id), sizeof(id));
            local += tile.id(v) == tid;
        }
        out.write(reinterpret_cast<char *>(&tid), sizeof(tid));
        out.write(reinterpret_cast<char *>(&local), sizeof(local));
    }
}

template<typename Tile>
void write_ply_property_vert(const Tile& tile, std::ostream& out)
{
    int D = tile.maximal_dimension();
    unsigned char tid = tile.id();
    for(auto it = tile.vertices_begin(); it != tile.vertices_end(); ++it)
    {
        if(tile.vertex_is_infinite(it)) continue;
        unsigned char id = tile.id(it);
        for(int d=0; d<D; ++d)
        {
            float coord = float(tile.coord(tile.point(it),d));
            out.write(reinterpret_cast<char *>(&coord), sizeof(coord));
        }
        out.write(reinterpret_cast<char *>(&tid), sizeof(tid));
        out.write(reinterpret_cast<char *>(&id), sizeof(id));
    }
}

template<typename DDT>
void write_ply_cell(const DDT& tri, const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    write_ply_header_begin(out);
    write_ply_element_cell(tri, out);
    write_ply_header_end(out);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_ply_property_cell(*tile, out);
    out.close();
}

template<typename DDT>
void write_ply_vert(const DDT& tri, const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    write_ply_header_begin(out);
    write_ply_element_vert(tri, out);
    write_ply_header_end(out);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_ply_property_vert(*tile, out);
    out.close();
}

template<typename DDT>
void write_ply(const DDT& tri, const std::string& filename)
{
    std::ofstream out(filename, std::ios::binary);
    write_ply_header_begin(out);
    write_ply_element_vert(tri, out);
    write_ply_element_cell(tri, out);
    write_ply_header_end(out);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_ply_property_vert(*tile, out);
    for(auto tile = tri.tiles_begin(); tile != tri.tiles_end(); ++tile)
        write_ply_property_cell(*tile, out);
    out.close();
}


// Create a ply from the cgal structure
template <typename TTr,typename DTC, typename FTC,typename Id>
std::ostream & cgal2ply_split(std::ostream & ofile,DTC & tri, FTC &filter, int nbc_finalized,std::string dump_mode,Id tid)
{
    TTr traits;
    typedef typename TTr::Vertex_handle                            Vertex_handle_raw;
    int D = TTr::D;
    char buffer[kBufferSize];
    double_conversion::StringBuilder builder(buffer, kBufferSize);
    double_conversion::DoubleToStringConverter dc(flags_deser, "Infinity", "NaN", 'e', 0, 0, 0, 0);
    int full_bufflen;
    char * buffer_char;
    int nbb;
    int pos = 0;
    int acc_pose = 0;
    char cc;
    bool do_simplex = false;
    int NB_DIGIT_OUT_PLY  = 3;
    CGAL::Unique_hash_map<Vertex_handle_raw, uint> vertex_map;
    int nb_vert = tri.number_of_vertices();
    int nb_cell = nbc_finalized;
    int nb_cell2 = 0;
    int nb_ply = 0;
    int max_bnc = 1000000;
    if(dump_mode == "TRIANGLE_SOUP")
    {
        nb_cell = nbc_finalized*4;
        max_bnc = 4000000;
    }
    std::string buffer_header("ply;");
    buffer_header.append("format ascii 1.0;");
    buffer_header.append("comment tid_" + std::to_string(tid) + "_0 ;");
    buffer_header.append("element vertex " + std::to_string(nb_vert) + ";");
    buffer_header.append("property double x;");
    buffer_header.append("property double y;");
    buffer_header.append("property double z;");
    buffer_header.append("element face " + std::to_string(nb_cell) + ";");
    buffer_header.append("property list uchar int vertex_index ;");
    buffer_header.append("end_header;");
    ofile << buffer_header ;
    std::stringstream sstr_d;
    full_bufflen = kBufferSize*nb_vert*D;
    buffer_char  = new char[full_bufflen];
    pos = 0;
    int ii=0;
    for(auto vit = traits.vertices_begin(tri); vit != traits.vertices_end(tri); ++vit)
    {
        if(tri.is_infinite(vit))
        {
            continue;
        }
        for(int d = 0; d < D; d++)
        {
            double dd = vit->point()[d];
            builder.Reset();
            dc.ToFixed(dd,NB_DIGIT_OUT_PLY,&builder);
            int pp = builder.position();
            memcpy( buffer_char + pos, buffer, pp );
            buffer_char[pos+pp] = ' ';
            pos += (pp+1);
        }
        if(dump_mode == "TRIANGLE_SOUP")
            buffer_char[pos-1] = ';';
        vertex_map[vit] = ii++;
    }
    ofile.write(buffer_char,pos);
    ofile << " ";
    delete[] buffer_char;
    if(nb_cell < max_bnc)
        full_bufflen = nb_cell*kBufferSize*(D+1);
    else
        full_bufflen = max_bnc*kBufferSize*(D+1);
    buffer_char = new char[full_bufflen];
    pos = 0;
    int acc = 0;
    if(dump_mode == "SIMPLEX_SOUP")
    {
        for(auto cit = traits.cells_begin(tri); cit != traits.cells_end(tri); ++cit)
        {
            if(!filter.do_keep(tri,cit,traits))
                continue;
            buffer_char[pos++] = '1'+D;
            buffer_char[pos++] = ' ';
            for(int d = 0; d < D+1; d++)
            {
                Id vid = vertex_map[cit->vertex(d)] ;
                pos += u32toa_countlut(vid,buffer_char + pos);
            }
            acc++;
            if(acc > max_bnc)
            {
                ofile.write(buffer_char,pos);
                ofile << std::endl;
                ofile << "ply tid_" << std::to_string(tid) << "_" << std::to_string(++nb_ply) << " ;";
                std::fill(buffer_char, buffer_char + full_bufflen, ' ');
                pos = 0;
                acc = 0;
            }
        }
    }
    else
    {
        for(auto cit = traits.cells_begin(tri); cit != traits.cells_end(tri); ++cit)
        {
            if(!filter.do_keep(tri,cit,traits))
                continue;
            for(int ss = 0; ss < D+1; ss++)
            {
                buffer_char[pos++] = '0'+D;
                buffer_char[pos++] = ' ';
                for(int d = 0; d < D+1; d++)
                {
                    if(d == ss)
                        continue;
                    Id vid = vertex_map[cit->vertex(d)] ;
                    pos += u32toa_countlut(vid,buffer_char + pos);
                }
                buffer_char[pos-1] = ';';
                acc++;
                if(acc > max_bnc)
                {
                    ofile.write(buffer_char,pos);
                    ofile << std::endl;
                    ofile << "ply tid_" << std::to_string(tid) << "_" << std::to_string(++nb_ply) << " ;";
                    std::fill(buffer_char, buffer_char + full_bufflen, ' ');
                    pos = 0;
                    acc = 0;
                }
            }
        }
    }
    ofile.write(buffer_char,pos);
    delete [] buffer_char;
    return ofile;
}


}

#endif // DDT_WRITE_PLY_HPP

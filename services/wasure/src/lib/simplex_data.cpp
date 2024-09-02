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
#include <string>

#include "simplex_data.hpp"
#include <CGAL/IO/io.h>
using namespace iqlib;

cell_data::cell_data() : bag_of_flags<std::uint64_t>(), idx(-1), tile_idx(-1), state(-1) { }

cell_data::cell_data(const cell_data& dt)
{
    copy(dt);
}

void cell_data::copy(const cell_data& dt)
{
    flags = dt.flags;
    idx = dt.idx;
    tile_idx = dt.tile_idx;
    state = dt.state;
}

cell_data& cell_data::operator =(const cell_data& dt)
{
    copy(dt);
    return *this;
}

bool cell_data::is_shared()
{
    return state == 1;
}

void cell_data::set_main(int val)
{
    bag_of_flags::flag(MAIN_FLAG,val);
}

bool cell_data::is_main()
{
    return bag_of_flags::flag(MAIN_FLAG);
}


void cell_data::write(std::ostream& os,bool only_iq,bool is_ascii) const
{
    if(is_ascii)
    {
        os << flags << " "  << idx << " " << tile_idx << " " << state << " ";
    }
    else
    {
        os.write((char*)(&(flags)), sizeof(flags));
        os.write((char*)(&(idx)), sizeof(idx));
        os.write((char*)(&(tile_idx)), sizeof(tile_idx));
        os.write((char*)(&(state)), sizeof(state));
    }
}

void cell_data::read(std::istream& is,bool only_iq, bool is_ascii)
{
    if(is_ascii)
    {
        is >> flags >> idx >> tile_idx >> state;
    }
    else
    {
        is.read((char*)(&(flags)), sizeof(flags));
        is.read((char*)(&(idx)), sizeof(idx));
        is.read((char*)(&(tile_idx)), sizeof(tile_idx));
        is.read((char*)(&(state)), sizeof(state));
    }
}


std::ostream& iqlib::operator<<(std::ostream& os, const cell_data& dt)
{
    os << std::to_string(dt.flags);
    return os;
}

std::istream& iqlib::operator>>(std::istream& is, cell_data& dt)
{
    std::string read{""};
    is >> read;
    dt.flags = std::stoull(read);
    return is;
}

bool iqlib::operator==(const cell_data& left, const cell_data& right)
{
    return left.flags == right.flags;
}

bool iqlib::operator!=(const cell_data& left, const cell_data& right)
{
    return !(left == right);
}

vertex_data::vertex_data() : bag_of_flags<std::uint64_t>(), idx(-1),tile_idx(-1), state(-1) { }

vertex_data::vertex_data(const vertex_data& dt) : bag_of_flags (dt.flags), idx(-1),tile_idx(-1), state(-1) { }


void vertex_data::copy(const vertex_data& dt)
{
    flags = dt.flags;
    idx = dt.idx;
    tile_idx = dt.idx;
    state = dt.state;
}

vertex_data& vertex_data::operator =(const vertex_data& dt)
{
    copy(dt);
    return *this;
}

bool vertex_data::is_shared()
{
    return state == 1;
}

void vertex_data::set_main(int val)
{
    bag_of_flags::flag(MAIN_FLAG,val);
}

bool vertex_data::is_main()
{
    return bag_of_flags::flag(MAIN_FLAG);
}


void vertex_data::write(std::ostream& os,bool only_iq,bool is_ascii) const
{
    if(is_ascii)
    {
        os << flags << " " << idx << " " << tile_idx << " " << state << " ";
    }
    else
    {
        os.write((char*)(&(flags)), sizeof(flags));
        os.write((char*)(&(idx)), sizeof(idx));
        os.write((char*)(&(tile_idx)), sizeof(tile_idx));
        os.write((char*)(&(state)), sizeof(state));
    }
}

void vertex_data::read(std::istream& is,bool only_iq,bool is_ascii)
{
    if(is_ascii)
    {
        is >> flags >> idx >> tile_idx >> state;
    }
    else
    {
        is.read((char*)(&(flags)), sizeof(flags));
        is.read((char*)(&(idx)), sizeof(idx));
        is.read((char*)(&(tile_idx)), sizeof(tile_idx));
        is.read((char*)(&(state)), sizeof(state));
    }
}

std::ostream& iqlib::operator<<(std::ostream& os, const vertex_data& dt)
{
    os << std::to_string(dt.flags);
    return os;
}

std::istream& iqlib::operator>>(std::istream& is, vertex_data& dt)
{
    std::string read{""};
    is >> read;
    dt.flags = std::stoull(read);
    return is;
}



bool iqlib::operator==(const vertex_data& left, const vertex_data& right)
{
    return left.flags == right.flags;
}

bool iqlib::operator!=(const vertex_data& left, const vertex_data& right)
{
    return !(left == right);
}



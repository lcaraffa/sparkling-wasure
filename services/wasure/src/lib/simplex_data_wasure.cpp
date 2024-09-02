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
#include <CGAL/IO/io.h>
#include "simplex_data_wasure.hpp"

using namespace iqlib;

#define NBD (4)

cell_data_wasure::cell_data_wasure() : cell_data::cell_data(), vpe(0), vpo(0),vpu(1),acc(0),lab(-1),dat(std::vector<typedat>(NBD))
{
    init_dat();
}

cell_data_wasure::cell_data_wasure(const cell_data_wasure& dt) : dat(std::vector<typedat>(NBD))
{
    copy(dt);
}

void cell_data_wasure::copy(const cell_data_wasure& dt)
{
    cell_data::copy(dt);
    lab = dt.lab;
    vpe = dt.vpe;
    vpo = dt.vpo;
    vpu = dt.vpu;
    acc = dt.acc;
    for(int i = 0; i < dat.size(); i++)
    {
        dat[i] = dt.dat[i];
    }
}

cell_data_wasure& cell_data_wasure::operator =(const cell_data_wasure& dt)
{
    copy(dt);
    return *this;
}

void cell_data_wasure::init_dat()
{
    for(int i = 0; i < dat.size(); i++)
    {
        dat[i] = 0;
    }
}


void cell_data_wasure::resize(int nb_samples)
{
    vpe.resize(nb_samples,0);
    vpo.resize(nb_samples,0);
    vpu.resize(nb_samples,1);
    pts.resize(nb_samples);
}

void cell_data_wasure::write(std::ostream& os,bool only_iq,bool is_ascii) const
{
    cell_data::write(os,only_iq, is_ascii);
    if(!only_iq)
    {
        if(is_ascii)
        {
            os << lab << " ";
            for(int i = 0; i < dat.size(); i++)
            {
                os << dat[i] << " ";
            }
            os << acc << " ";
        }
        else
        {
            os.write((char*)(&(lab)), sizeof(lab));
            for(int i = 0; i < dat.size(); i++)
            {
                os.write(reinterpret_cast<const char *>(&dat[i]), sizeof(dat[i]));
            }
            os.write((char*)(&(acc)), sizeof(acc));
        }
    }
}

void cell_data_wasure::read(std::istream& is,bool only_iq,bool is_ascii)
{
    cell_data::read(is,only_iq,is_ascii);
    if(!only_iq)
    {
        if(is_ascii)
        {
            is >> lab;
            for(int i = 0; i < dat.size(); i++)
                is >> dat[i];
            is >> acc;
        }
        else
        {
            is.read((char*)(&(lab)), sizeof(lab));
            for(int i = 0; i < dat.size(); i++)
                is.read(reinterpret_cast<char *>(&dat[i]), sizeof(dat[i]));
            is.read((char*)(&(acc)), sizeof(acc));
        }
    }
}

std::ostream& iqlib::operator<<(std::ostream& os, const cell_data_wasure& dt)
{
    return os;
}

std::istream& iqlib::operator>>(std::istream& is, cell_data_wasure& dt)
{
    std::string read{""};
    is >> read;
    return is;
}


vertex_data_wasure::vertex_data_wasure() : vertex_data(), tacc(0), tweig(0), acc(0) { }

void vertex_data_wasure::copy(const vertex_data_wasure& dt)
{
    vertex_data::copy(dt);
}

vertex_data_wasure::vertex_data_wasure(const vertex_data_wasure& dt)
{
    copy(dt);
}




void vertex_data_wasure::write(std::ostream& os,bool only_iq,bool is_ascii ) const
{
    vertex_data::write(os,only_iq,is_ascii);
}

void vertex_data_wasure::read(std::istream& is,bool only_iq,bool is_ascii )
{
    vertex_data::read(is,only_iq,is_ascii);
}

std::ostream& iqlib::operator<<(std::ostream& os, const vertex_data_wasure& dt)
{
    os << std::to_string(dt.flags);
    return os;
}

std::istream& iqlib::operator>>(std::istream& is, vertex_data_wasure& dt)
{
    std::string read{""};
    is >> read;
    dt.flags = std::stoull(read);
    return is;
}

bool iqlib::operator==(const vertex_data_wasure& left, const vertex_data_wasure& right)
{
    return left.flags == right.flags;
}

bool iqlib::operator!=(const vertex_data_wasure& left, const vertex_data_wasure& right)
{
    return !(left == right);
}




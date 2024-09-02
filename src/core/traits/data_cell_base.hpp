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
#ifndef DATA_CELL_BASE_HPP
#define DATA_CELL_BASE_HPP

#include "../bbox.hpp"

namespace ddt
{

template<typename I>
struct Data_cell_base
{
    typedef I Id;

    Data_cell_base(Id i=0) : id(i) {}

    void write(std::ostream& os,bool is_ascii) const
    {
        if(is_ascii)
        {
            os << id << " ";
        }
        else
        {
            os.write((char*)(&(id)), sizeof(id));
        }
    }
    void read(std::istream& is, bool is_ascii)
    {
        if(is_ascii)
        {
            is >> id;
        }
        else
        {
            is.read((char*)(&(id)), sizeof(id));
        }
    }


    Id id;
};


}

#endif // DATA_CELL_BASE_HPP

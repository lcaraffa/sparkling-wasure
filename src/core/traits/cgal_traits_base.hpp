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
#ifndef DDT_CGAL_TRAITS_BASE_HPP
#define DDT_CGAL_TRAITS_BASE_HPP

#include "../bbox.hpp"

namespace ddt
{

template<typename I, typename F>
struct Data
{
    typedef I Id;
    typedef F Flag;

    Data(Id i=0,Id gid=-1, Flag f=0) : id(i),gid(i), flag(f) {}

    void write(std::ostream& os,bool is_ascii) const
    {
        if(is_ascii)
        {
            os << id << " "  << gid << " " <<  flag << " ";
        }
        else
        {
            os.write((char*)(&(id)), sizeof(id));
            os.write((char*)(&(gid)), sizeof(id));
            os.write((char*)(&(flag)), sizeof(flag));
        }
    }
    void read(std::istream& is, bool is_ascii)
    {
        if(is_ascii)
        {
            is >> id >> gid >> flag;
        }
        else
        {
            is.read((char*)(&(id)), sizeof(id));
            is.read((char*)(&(gid)), sizeof(id));
            is.read((char*)(&(flag)), sizeof(flag));
        }
    }


    Id id,gid;
    mutable Flag flag;
};



template<typename I, typename F>
struct DataF
{
    typedef I Id;
    typedef F Flag;

    DataF() :  flag(0) {}
    DataF(Flag f) :  flag(f) {}
    mutable Flag flag;
};


template<typename T1, typename T2>
struct Pair2nd : public std::unary_function<const T1&, std::pair<T1,T2>>
{
    T2 _2;
    Pair2nd(const T2& x) : _2(x) {}
    inline std::pair<T1,T2> operator()(const T1& _1) const
    {
        return std::make_pair(_1, _2);
    }
};

}

#endif // DDT_CGAL_TRAITS_BASE_HPP

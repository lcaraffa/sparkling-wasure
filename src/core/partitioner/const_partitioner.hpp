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
#ifndef DDT_CONST_PARTITIONER_HPP
#define DDT_CONST_PARTITIONER_HPP

namespace ddt
{

template<typename Traits>
class const_partitioner
{
public:
    typedef typename Traits::Point Point;
    typedef typename Traits::Id    Id;

    const_partitioner(Id id) : id_(id) {}

    inline Id operator()(const Point& p) const
    {
        return id_;
    }

private:
    Id id_;
};

}

#endif // DDT_CONST_PARTITIONER_HPP

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
#ifndef DDT_GRID_PARTITIONER_HPP
#define DDT_GRID_PARTITIONER_HPP

#include "../bbox.hpp"

namespace ddt
{

template<typename Traits>
class grid_partitioner
{
public:
    typedef typename Traits::Point Point;
    typedef typename Traits::Id    Id;
    enum { D = Traits::D };

    template<typename Iterator>
    grid_partitioner(const Bbox<D>& bbox, Iterator it, Iterator end)
    {
        M = 1;
        int n = 1;
        for(int i=0; i<D; ++i)
        {
            if (it!=end) n = *it++;
            N[i] = n;
            M *= n;
            inv_step[i] = n/(bbox.max(i)-bbox.min(i));
            origin[i] = bbox.min(i);
        }
    }

    grid_partitioner(const Bbox<D>& bbox, int n)
    {
        M = 1;
        for(int i=0; i<D; ++i)
        {
            N[i] = n;
            M *= n;
            inv_step[i] = n/(bbox.max(i)-bbox.min(i));
            origin[i] = bbox.min(i);
        }
    }

    Id operator()(const Point& p) const
    {
        int id = 0;
        for(int i=0; i<D; ++i)
            id = id*N[i] + (int((p[i]-origin[i])*inv_step[i]) % N[i]);
        return id;
    }
    const int *begin() const { return N; }
    const int *end() const { return N+D; }
    int size() const { return M; }

private:
    int N[D], M;
    double inv_step[D];
    double origin[D];
};

}

#endif // DDT_GRID_PARTITIONER_HPP

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
#ifndef DDT_RANDOM_PARTITIONER_HPP
#define DDT_RANDOM_PARTITIONER_HPP

#include <chrono>
#include <random>

namespace ddt
{

template<typename Traits, typename Generator = std::default_random_engine>
class random_partitioner
{
public:
    typedef typename Traits::Point Point;
    typedef typename Traits::Id    Id;

    random_partitioner(Id a, Id b, unsigned int seed = 0) : distribution(a,b), generator(seed)
    {
        if(seed == 0)
        {
            seed = std::chrono::system_clock::now().time_since_epoch().count();
            generator.seed(seed);
        }
        std::cout << "random_partitioner seed: " << seed << std::endl;
    }

    random_partitioner(Id a, Id b, const Generator& g) : distribution(a,b), generator(g)
    {
    }

    inline Id operator()(const Point& p)
    {
        return distribution(generator);
    }

private:
    std::uniform_int_distribution<Id> distribution;
    Generator generator;
};

}

#endif // DDT_RANDOM_PARTITIONER_HPP

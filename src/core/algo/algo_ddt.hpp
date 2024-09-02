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
#ifndef algo_ddt_algo_ddt_HPP
#define algo_ddt_algo_ddt_HPP

#include <string>
#include <unordered_map>

#include "io/logging.hpp"



namespace ddt
{

template<typename DDT, typename Scheduler>
class algo_ddt
{
public:

    typedef typename DDT::Traits                                        Traits;
    typedef typename DDT::Tile                   Tile;



    algo_ddt(DDT & ddt_, Scheduler & sch_) : ddt(ddt_),sch(sch_)
    {
        log(0, sch.number_of_threads(), " thread(s)\n");
    }



    int for_each(const std::string& step, const std::string& type, const std::function<int(Tile&, bool)>& func)
    {
        log.step(step,type);
        return sch.for_each(ddt.tiles_begin(), ddt.tiles_end(), func);
    }
    int splay()
    {
        log.step("Splay", "Rcv   ");
        return sch.for_each(ddt.tiles_begin(), ddt.tiles_end(), sch.splay_func());
    }
    int send_all_bbox_points()
    {
        log.step("Send", "Loc+BB");
        return sch.for_each(ddt.tiles_begin(), ddt.tiles_end(), sch.send_all_func(&Tile::get_bbox_points, ddt.tile_ids_begin(), ddt.tile_ids_end()));
    }
    int send_all_local_convex_hull()
    {
        log.step("Send", "CHull ");
        return sch.for_each(ddt.tiles_begin(), ddt.tiles_end(), sch.send_all_func(&Tile::get_local_convex_hull, ddt.tile_ids_begin(), ddt.tile_ids_end()));
    }
    int splay_rec_neighbors()
    {
        log.step("Splay", "Star  ");
        return sch.for_each_rec(ddt.tiles_begin(), ddt.tiles_end(), sch.splay_one_func(&Tile::get_neighbors));
    }
    int splay_rec2_local_neighbors()
    {
        log.step("Splay", "Star  ");
        return sch.for_each_rec2(ddt.tiles_begin(), ddt.tiles_end(), sch.splay_one_func(&Tile::get_local_neighbors));
    }
    int splay_rec2_neighbors()
    {
        log.step("Splay", "Star  ");
        return sch.for_each_rec2(ddt.tiles_begin(), ddt.tiles_end(), sch.splay_one_func(&Tile::get_neighbors));
    }
    int splay_one_active_points()
    {
        log.step("Splay", "Active");
        return sch.for_each(ddt.tiles_begin(), ddt.tiles_end(), sch.splay_one_func(&Tile::get_active_points));
    }



private:


    Scheduler sch;
    logging log;
    DDT ddt;

};

}

#endif // algo_ddt_DDT_HPP

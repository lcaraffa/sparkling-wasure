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

#ifndef DDT_TYPEDEFS_H
#define DDT_TYPEDEFS_H


typedef long int Id;
typedef int FlagV;
typedef int FlagC;



#include <CGAL/Random.h>

#include <stdio.h>      /* printf */
#include <math.h>
#include <io/stream_api.hpp>
#include "io/write_stream.hpp"
#include "io/write_vrt.hpp"
#include "io/read_stream.hpp"

#include "tile.hpp"
#include "traits/traits.hpp"
#include "scheduler/scheduler.hpp"
#include "partitioner/grid_partitioner.hpp"
#include "partitioner/const_partitioner.hpp"
#include "DDT.hpp"
#include "io/logging_stream.hpp"

#include "ddt_data.hpp"
#include "ddt_spark_utils.hpp"
#include "io_ddt_stream.hpp"


typedef ddt::Data<Id,FlagV>                                  Data_V;
typedef ddt::Data<Id,FlagC>                                     Data_C;
typedef ddt::Traits<Data_V,Data_C> Traits;

typedef ddt::Tile<Traits> Tile;
typedef ddt::Scheduler<Tile> Scheduler;
typedef ddt::DDT<Traits> DDT;
typedef ddt::grid_partitioner<Traits> Grid_partitioner;


//typedef Traits::Random_points_in_box Random_points;
typedef typename DDT::Tile_const_iterator  Tile_const_iterator ;
typedef typename DDT::Tile_iterator  Tile_iterator ;

typedef typename Traits::Point       Point;
//typedef typename Traits::Id          Id;
typedef typename Traits::Point_id    Point_id;
typedef typename Traits::Point_id_id Point_id_id;


typedef ddt::Traits_raw Traits_raw;

typedef typename ddt::Traits_raw::Delaunay_triangulation       DT_raw;
typedef typename Traits_raw::Vertex_handle_raw        Vertex_handle_raw;
typedef typename Traits_raw::Cell_handle_raw        Cell_handle_raw;


typedef typename Tile::Vertex_const_handle_and_id Vertex_const_handle_and_id;
typedef typename Tile::Vertex_const_handle Vertex_const_handle;


typedef typename Tile::Point_id Point_and_id;
typedef std::tuple<Point,Id,Id>                  Point_id_source;


#endif


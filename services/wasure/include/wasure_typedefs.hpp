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
#ifndef WASURE_TYPEDEFS_H
#define WASURE_TYPEDEFS_H

#include <CGAL/Random.h>

#include <stdio.h>      /* printf */
#include <math.h>

#include "tile.hpp"
#include "traits/traits.hpp"
#include "traits/data_cell_base.hpp"
#include "scheduler/scheduler.hpp"
#include "partitioner/grid_partitioner.hpp"
#include "partitioner/const_partitioner.hpp"
#include "simplex_data_wasure.hpp"

#include "DDT.hpp"



typedef int Id;
typedef int FlagV;
typedef int FlagC;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef ddt::Data<Id,FlagV>                                  Data_V;
typedef ddt::Data<Id,FlagC>                                     Data_C;
typedef ddt::Traits<Data_V,Data_C> Traits;
typedef ddt::Traits_raw Traits_raw;


typedef ddt::Tile<Traits> Tile;
typedef ddt::Scheduler<Tile> Scheduler;
typedef ddt::DDT<Traits> DTW;


typedef ddt::grid_partitioner<Traits> Grid_partitioner;



typedef typename DTW::Tile_const_iterator  Tile_const_iterator ;
typedef typename DTW::Tile_cell_const_handle Tile_cell_const_handle;
typedef typename DTW::Tile_iterator  Tile_iterator ;

typedef typename Traits::Delaunay_triangulation DT;
typedef typename ddt::Traits_raw::Delaunay_triangulation       DT_raw;
typedef typename Traits::Point       Point;
typedef typename Traits::Sphere       Sphere;
typedef typename Traits::Plane       Plane;
typedef typename Traits::Vector       Vector;
typedef typename Traits::Id          Id;
typedef typename Traits::Point_id    Point_id;
typedef typename Traits::Point_id_id Point_id_id;


typedef typename Traits::Point                                    Point;
typedef typename Traits::Point                                    Facet;
typedef typename Traits::Cell_handle                                    Cell_handle;
typedef typename Traits::Vertex_handle                                    Vertex_handle;


typedef typename DTW::Facet_const_iterator  Facet_const_iterator;
typedef typename DTW::Vertex_const_iterator  Vertex_const_iterator;
typedef typename DTW::Tile_cell_const_handle              Tile_cell_const_handle;

typedef typename DTW::Cell_const_iterator                 Cell_const_iterator;
typedef typename DTW::Facet_const_iterator                Facet_const_iterator;




typedef typename Tile::Vertex_const_handle_and_id Vertex_const_handle_and_id;
typedef typename Tile::Vertex_const_handle Vertex_const_handle;




typedef typename Tile::Point_id Point_id;
typedef typename Tile::Point_id_id      Point_id_source;



#endif

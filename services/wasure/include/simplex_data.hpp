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
#ifndef SIMPLEX_DATA_HPP
#define SIMPLEX_DATA_HPP

#include "flags.hpp"
#include "typedefs_type.hpp"
#include <iostream>
#include <cstdint>

#define MAIN_FLAG (0)
namespace iqlib
{

/**
   Cell data class
   2D space.

   @author Laurent Caraffa
*/
class cell_data : public bag_of_flags<std::uint64_t>
{
public :
    int idx,tile_idx;
    int state;
    cell_data();                                // default conclassor
    cell_data(const cell_data& dt);             // copy conclassor
    cell_data& operator=(const cell_data& dt);  // assignment operator
    void copy (const cell_data& dt);


    bool is_shared();
    void set_main(int val);
    bool is_main();

    void write(std::ostream& os,bool only_iq,bool is_ascii) const;
    void read(std::istream& is,bool only_iq,bool is_ascii);

    friend std::ostream& operator<<(std::ostream& os, const cell_data& dt);
    friend std::istream& operator>>(std::istream& is, cell_data& dt);
    friend bool operator==(const cell_data& left, const cell_data& right);
    friend bool operator!=(const cell_data& left, const cell_data& right);
};

std::ostream& operator<<(std::ostream& os, const cell_data& dt);
std::istream& operator>>(std::istream& is, cell_data& dt);
bool operator==(const cell_data& left, const cell_data& right);
bool operator!=(const cell_data& left, const cell_data& right);

/**
   Vertex data class
   2D space.

   @author Laurent Caraffa
*/
class vertex_data : public bag_of_flags<std::uint64_t>
{
public :
    int idx,tile_idx;
    int state;
    vertex_data();
    vertex_data(const vertex_data& dt);
    vertex_data& operator=(const vertex_data& dt);
    void copy(const vertex_data& dt);

    bool is_shared();
    void set_main(int val);
    bool is_main();

    void write(std::ostream& os,bool only_iq,bool is_ascii ) const;
    void read(std::istream& is,bool only_iq,bool is_ascii );

    friend std::ostream& operator<<(std::ostream& os, const vertex_data& dt);
    friend std::istream& operator>>(std::istream& is, vertex_data& dt);
    friend bool operator==(const vertex_data& left, const vertex_data& right);
    friend bool operator!=(const vertex_data& left, const vertex_data& right);
};

std::ostream& operator<<(std::ostream& os, const vertex_data& dt);
std::istream& operator>>(std::istream& is, vertex_data& dt);
bool operator==(const vertex_data& left, const vertex_data& right);
bool operator!=(const vertex_data& left, const vertex_data& right);



}      // namespace iqlib
#endif // SIMPLEX_DATA_HPP

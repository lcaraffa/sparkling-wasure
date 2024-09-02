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
#ifndef SIMPLEX_DATA_WASURE_HPP
#define SIMPLEX_DATA_WASURE_HPP

#include "./simplex_data.hpp"
#include "./flags.hpp"
#include <vector>
#include <iostream>


namespace iqlib
{


class cell_data_wasure : public cell_data
{
public :

    // typelab lab;
    std::vector<typedat> dat;
    std::vector<std::vector<typedat> > pts;
    std::vector<typedat> vpe,vpo,vpu;
    int acc,seg;
    int opt_idx;
    int lab;

    void init_dat();
    void resize(int nb_samples);
    int get_lab() {return lab;};

    cell_data_wasure();                                // default conclassor
    cell_data_wasure(const cell_data_wasure& dt);             // copy conclassor
    cell_data_wasure& operator=(const cell_data_wasure& dt);  // assignment operator
    void copy(const cell_data_wasure& dt);


    void write(std::ostream& os,bool only_iq,bool is_ascii ) const;
    void read(std::istream& is,bool only_iq,bool is_ascii );

    friend std::ostream& operator<<(std::ostream& os, const cell_data_wasure& dt);
    friend std::istream& operator>>(std::istream& is, cell_data_wasure& dt);
    friend bool operator==(const cell_data_wasure& left, const cell_data_wasure& right);
    friend bool operator!=(const cell_data_wasure& left, const cell_data_wasure& right);

};

template<typename Traits>
std::ostream& operator<<(std::ostream& os, const cell_data_wasure& dt);
template<typename Traits>
std::istream& operator>>(std::istream& is, cell_data_wasure& dt);
template<typename Traits>
bool operator==(const cell_data_wasure& left, const cell_data_wasure& right);
template<typename Traits>
bool operator!=(const cell_data_wasure& left, const cell_data_wasure& right);

class vertex_data_wasure : public vertex_data
{
public :

    std::vector<double> tacc;
    std::vector<double> tweig;
    int acc;
    vertex_data_wasure();
    vertex_data_wasure(const vertex_data_wasure& dt);
    vertex_data_wasure& operator=(const vertex_data_wasure& dt);
    void copy(const vertex_data_wasure& dt);


    void write(std::ostream& os,bool only_iq,bool is_ascii) const;
    void read(std::istream& is,bool only_iq,bool is_ascii);

    friend std::ostream& operator<<(std::ostream& os, const vertex_data_wasure& dt);
    friend std::istream& operator>>(std::istream& is, vertex_data_wasure& dt);
    friend bool operator==(const vertex_data_wasure& left, const vertex_data_wasure& right);
    friend bool operator!=(const vertex_data_wasure& left, const vertex_data_wasure& right);

};

std::ostream& operator<<(std::ostream& os, const vertex_data_wasure& dt);
std::istream& operator>>(std::istream& is, vertex_data_wasure& dt);
bool operator==(const vertex_data_wasure& left, const vertex_data_wasure& right);
bool operator!=(const vertex_data_wasure& left, const vertex_data_wasure& right);



}      // namespace iqlib
#endif // SIMPLEX_DATA_WASURE_HPP

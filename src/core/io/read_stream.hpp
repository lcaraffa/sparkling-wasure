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
#ifndef DDT_READ_STREAM_HPP
#define DDT_READ_STREAM_HPP

#include <unordered_map>
#include <map>
#include <set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "io/stream_api.hpp"
#include <algorithm>
#include <cctype>
#include <locale>

namespace ddt
{




template <typename Traits>
std::istream & read_point_set_serialized(std::vector<typename Traits::Point> & lp, std::istream & ifile, const Traits& traits)
{
    char buffer[kBufferSize];
    char cc;
    uint D;
    int nbp,pos,nb_char,processed;
    double ddv;
    ifile >> D >> nbp ;
    ifile.get(cc);
    double coords[Traits::D];
    std::vector<double> input_v;
    deserialize_b64_vect(input_v,ifile);
    for(int n = 0; n< input_v.size()/(D); n++)
    {
        for(int d = 0; d < D; d++)
            coords[d] = input_v[n*D+d];
        lp.push_back(traits.make_point(coords));
    }
    input_v.clear();
    return ifile;
}



inline double read_double(std::istream & ifile)
{
    double dd;
    if(IS_BINARY)
        ifile.read(reinterpret_cast<char *>(&dd), sizeof(dd));
    else
        ifile >> dd;
    return dd;
}


template <typename Traits>
std::istream & read_points_id_source_serialized(std::vector<typename Traits::Point_id_id> & lp,
        std::istream & ifile, const Traits& traits)
{
    int D = Traits::D;
    double coords[Traits::D];
    std::vector<double> input_v;
    deserialize_b64_vect(input_v,ifile);
    for(int n = 0; n< input_v.size()/(D+2); n++)
    {
        int id1 = input_v[n*(D+2)];
        int id2 = input_v[n*(D+2)+1];
        for(int d = 0; d < D; d++)
            coords[d] = input_v[n*(D+2)+2+d];
        lp.push_back(std::make_tuple(traits.make_point(coords),id1,id2));
    }
    return ifile;
}




template<typename Tile, typename Id>
std::istream& read_json_stream(Tile & tile,std::istream&  ifile)
{
    boost::property_tree::ptree root_node;
    boost::property_tree::read_xml(ifile, root_node);
    auto & bbox = tile->bbox();
    for (auto its : root_node.get_child("bbox"))
    {
        std::string input = its.first;
        trim(input);
        int iid;
        try
        {
            if (input.empty())
            {
                throw std::invalid_argument("Input string is empty");
            }
            iid  = std::stoi(input);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "[read json stream] Invalid argument: " << e.what() << " for input: " << input << std::endl;
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "Out of range: " << e.what() << " for input: " << input << std::endl;
        }
        Id id = iid;
        std::stringstream ss (its.second.data());
        ss >> bbox[id];
    }
    return ifile;
}



template <typename Traits>
typename Traits::Point read_point(std::istream & ifile, const Traits& traits)
{
    double coords[Traits::D];
    for(int d = 0 ; d < Traits::D; d++)
        ifile >> coords[d];
    return traits.make_point(coords);
}


template <typename Traits>
typename Traits::Point_id_id read_point_id_source(std::istream & ifile, const Traits& traits)
{
    double coords[Traits::D];
    typename Traits::Id id1,id2;
    ifile >> id1;
    ifile >> id2;
    for(int d = 0 ; d < Traits::D; d++)
        ifile >> coords[d];
    return  std::make_tuple(traits.make_point(coords),id1,id2);
}

template <typename Traits>
std::istream & read_points_stream(std::vector<typename Traits::Point> & lp, std::istream & ifile, const Traits& traits)
{
    int nbp = -1;
    uint D;
    ifile >> D >> nbp;
    assert(D == Traits::D);
    if(nbp == -1)
        return ifile;
    for(int i = 0; i < nbp; i++)
    {
        lp.push_back(read_point(ifile, traits));
    }
    return ifile;
}

template <typename Traits>
std::istream & read_points_stream(std::set<typename Traits::Point> & lp, std::istream & ifile, const Traits& traits)
{
    int nbp = -1;
    uint D;
    ifile >> D >> nbp;
    assert(D == Traits::D);
    if(nbp == -1)
        return ifile;
    for(int i = 0; i < nbp; i++)
    {
        lp.insert(read_point(ifile, traits));
    }
    return ifile;
}



template <typename Traits>
std::istream & read_points_id_source_stream(std::vector<typename Traits::Point_id_id> & lp,
        std::istream & ifile, const Traits& traits)
{
    int nbp = -1;
    uint D;
    ifile >> D >> nbp;
    assert(D == Traits::D);
    if(nbp == -1)
        return ifile;
    for(int i = 0; i < nbp; i++)
    {
        lp.push_back(read_point_id_source(ifile, traits));
    }
    return ifile;
}

template <typename Traits>
std::istream & read_map_stream(std::map<typename Traits::Id, std::set<typename Traits::Point>> & mp, std::istream & ifile, const Traits& traits)
{
    int nb_elem;
    ifile >> nb_elem;
    for(int i = 0; i < nb_elem; i++)
    {
        typename Traits::Id id;
        ifile >> id;
        std::set<typename Traits::Point> sp;
        read_points_stream(sp, ifile, traits);
        mp[id] = sp;
    }
    return ifile;
}


template<typename DDTT>
int read_tile_stream(DDTT & ddt, std::istream & ifile, typename DDTT::Id tid, bool do_data = true, bool is_ascii = true)
{
    ddt.init(tid);
    auto tile  = ddt.get_tile(tid);
    tile->read_cgal(ifile,do_data,is_ascii);
    read_map_stream(tile->points_sent_,ifile,tile->traits());
    read_json_stream<typename DDTT::Tile_iterator, typename DDTT::Id>(tile,ifile);
    tile->set_id(tid);
    tile->finalize();
    return 0;
}

template<typename DDTT>
int read_full_stream(DDTT & ddt, std::istream & ifile, int nb_dat, bool do_data = true, bool is_ascii = true)
{
    for(int i = 0; i < nb_dat; i++)
    {
        stream_data_header hpi;
        hpi.parse_header(ifile);
        if(hpi.get_lab() == "t")
        {
            read_tile_stream(ddt, hpi.get_input_stream(), hpi.get_id(0),do_data,is_ascii);
        }
    }
    return 0;
}

template<typename DDTT>
void read_stream(DDTT & ddt, std::string ss)
{
    ddt::stream_app_header sah;
    std::istringstream iss(ss);
    sah.parse_header(iss);
    if(sah.is_void())
        return;
    read_full_stream(ddt, iss, sah.get_nb_dat());
}

}

#endif

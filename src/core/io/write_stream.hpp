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
#ifndef DDT_WRITE_STREAM_HPP
#define DDT_WRITE_STREAM_HPP

#include <unordered_map>
#include <map>
#include <set>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "number_parser.hpp"
#include "io/stream_api.hpp"


namespace ddt
{


template <typename Point>
std::ostream & write_point_set_serialized(const std::vector<Point> & lv, std::ostream & ofile, uint D)
{
    int num_v = lv.size();
    ofile << D << " " << num_v << " ";
    std::vector<double> outputv;
    for(auto pp : lv)
    {
        for(int d = 0 ; d < D; d++)
            outputv.emplace_back(pp[d]);
    }
    serialize_b64_vect(outputv,ofile);
    return ofile;
}



inline  std::ostream & write_double(double dd,std::ostream & ofile)
{
    if(IS_BINARY)
        ofile.write(reinterpret_cast<char *>(&dd), sizeof(dd));
    else
        ofile << dd << " ";
    return ofile;
}





template<typename Tile>
std::ostream & write_json_stream(const Tile & tile,std::ostream & ofile)
{
    boost::property_tree::ptree root_node;
    boost::property_tree::ptree bbox_node;
    root_node.put("id", tile->id());
    root_node.put("nbmc", tile->number_of_main_cells());
    root_node.put("nbmv", tile->number_of_main_vertices());
    root_node.put("nbmf", tile->number_of_main_facets());
    auto & bbox = tile->bbox();
    for(auto iter = bbox.begin(); iter != bbox.end(); ++iter)
    {
        std::stringstream ss;
        ss << iter->second;
        bbox_node.put(std::to_string(iter->first),ss.str());
    }
    root_node.add_child("bbox", bbox_node);
    std::stringstream ss;
    boost::property_tree::write_xml(ss, root_node);
    std::string str = ss.str();
    str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
    ofile << str;
    return ofile;
}

template <typename Point>
std::ostream & write_point(Point & pp,std::ostream & ofile,int D)
{
    for(int d = 0 ; d < D; d++)
        write_double(pp[d],ofile);
    return ofile;
}


template <typename Point_id_source,typename Point>
std::ostream & write_point_id_source(Point_id_source & pp,std::ostream & ofile,int D)
{
    ofile << std::get<1>(pp) << " " << std::get<2>(pp) << " ";
    for(int d = 0 ; d < D; d++)
        write_double(std::get<0>(pp)[d],ofile);
    return ofile;
}
template <typename Vertex_const_handle>
int write_vch_stream(const std::vector<Vertex_const_handle> & lp, std::ostream & ofile, uint D)
{
    ofile << D << " " << lp.size() << " ";
    for(auto pp : lp)
    {
        write_point(pp->point(),ofile,D);
    }
    return 0;
}
template <typename set_pts>
int write_points_stream(const set_pts & lp, std::ostream & ofile, uint D)
{
    ofile << D << " " << lp.size() << " ";
    for(auto pp : lp)
    {
        write_point(pp,ofile,D);
    }
    return 0;
}


template <typename Point_id_source,typename Point>
int write_points_id_source_serialized(const std::vector<Point_id_source> & lp, std::ostream & ofile, uint D)
{
    std::vector<double> outputv;
    for(auto pp : lp)
    {
        outputv.emplace_back(std::get<1>(pp));
        outputv.emplace_back(std::get<2>(pp));
        for(int d = 0 ; d < D; d++)
            outputv.emplace_back(std::get<0>(pp)[d]);
    }
    serialize_b64_vect(outputv,ofile);
    return 0;
}


template <typename Point_id_source,typename Point>
int write_points_id_source_stream(const std::vector<Point_id_source> & lp, std::ostream & ofile, uint D)
{
    std::vector<double> outputv;
    for(auto pp : lp)
    {
        outputv.emplace_back(std::get<1>(pp));
        outputv.emplace_back(std::get<2>(pp));
        for(int d = 0 ; d < D; d++)
            outputv.emplace_back(std::get<0>(pp)[d]);
    }
    serialize_b64_vect(outputv,ofile);
    return 0;
}



template <typename Point>
std::ostream & write_vertex_set_stream(const std::set<Point> & lv, std::ostream & ofile, uint D)
{
    ofile << D << " " << lv.size() << " ";
    for(auto vv : lv)
    {
        write_point(vv->point(),ofile,D);
    }
    return ofile;
}


template<typename Id,typename Point>
std::ostream & write_map_stream( const std::map<Id, std::set<Point>> & mp, std::ostream & ofile, int D)
{
    ofile << mp.size() << " ";
    for ( auto it = mp.begin(); it != mp.end(); it++ )
    {
        ofile << it->first << " ";
        write_points_stream<std::set<Point>>(it->second,ofile,D);
    }
    return ofile;
}


template<typename DDT>
std::ostream & write_tile_stream(const DDT& ddt, std::ostream & ofile, int tid)
{
    auto tile  = ddt.get_tile(tid);
    tile->write_cgal(ofile);
    write_map_stream(tile->points_sent_,ofile,tile->current_dimension());
    write_json_stream(tile,ofile);
    return ofile; // FIXME ?
}



template<typename DDT>
int write_full_stream(const DDT& ddt, std::istream & ifile, int nb_dat)
{
    std::vector<int> lnodes;
    for(int i = 0; i < nb_dat; i++)
    {
        stream_data_header hpi;
        hpi.parse_header(ifile);
        if(hpi.get_lab() == "t")
        {
            write_tile_stream(ddt, hpi.get_input_stream(), hpi.get_id(0));
            lnodes.push_back(hpi.get_id(0));
        }
    }
    return 0;
}

template<typename DDT>
std::string write_stream(const DDT& ddt, std::string output_dir, int tid)
{
    stream_data_header oqh("t","f",tid);
    oqh.set_file_name(output_dir + "/tri" + std::to_string(tid) + ".bin");
    std::stringstream ss;
    oqh.write_header(ss);
    write_tile_stream(ddt, oqh.get_output_stream(), tid);
    oqh.finalize();
    return ss.str();
}





}

#endif // DDT_WRITE_STREAM_HPP

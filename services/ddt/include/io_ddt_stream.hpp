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


#ifndef IO_DDT_STREAM_HPP
#define IO_DDT_STREAM_HPP

#include "ddt_data.hpp"
#include "io/write_stream.hpp"
#include "io/read_stream.hpp"
#include "io/logging_stream.hpp"


namespace ddt
{



template<typename DDTT,typename Scheduler >
int read_ddt_full_stream(DDTT & ddt, std::istream & ifile, int nb_dat,ddt::logging_stream & log)
{
    for(int i = 0; i < nb_dat; i++)
    {
        stream_data_header hpi;
        hpi.parse_header(ifile);
        bool do_serialize = false;
        bool do_clean_data = true;
        if(hpi.get_lab() == "t")
        {
            read_ddt_stream(ddt,hpi.get_input_stream(), hpi.get_id(0),do_serialize,do_clean_data,log);
        }
    }
    Scheduler sch(1);
    ddt.finalize(sch);
    return 0; // FIXME ?
}



template<typename DDT,typename WD>
int read_ddt_stream(DDT & ddt, WD & wd,  std::istream & ifile, typename DDT::Id tid,bool do_serialize, bool do_clean_data, ddt::logging_stream & log)
{
    int dd;
    ifile >> dd;
    if(dd == 1)
        wd.read_serialized_stream(ifile);
    read_ddt_stream_core(ddt,ifile, tid,do_serialize,  do_clean_data, log);
    return 0;
}

template<typename DDT>
int read_ddt_stream(DDT & ddt,  std::istream & ifile, typename DDT::Id tid,bool do_serialize, bool do_clean_data, ddt::logging_stream & log)
{
    int dd;
    ifile >> dd;
    if(dd != 0)
        std::cerr <<  "DD should be 0!!";
    read_ddt_stream_core(ddt,ifile, tid,do_serialize,  do_clean_data, log);
    return 0;
}


template<typename DDT>
int read_ddt_stream_core(DDT & ddt,  std::istream & ifile, typename DDT::Id tid,bool do_serialize, bool do_clean_data, ddt::logging_stream & log)
{
    typename DDT::Traits ttraits;
    ddt.init(tid);
    auto  tile  = ddt.get_tile(tid);
    typename DDT::Traits::Delaunay_triangulation & ttri = tile->tri();
    ttraits.deserialize_b64_cgal(ttri,ifile);
    deserialize_b64_vect(tile->tile_ids,ifile);
    std::string input;
    std::getline(ifile, input);
    std::stringstream ifile2(input);
    ddt::read_map_stream(tile->points_sent_,ifile2,tile->traits());
    ddt::read_json_stream<typename DDT::Tile_iterator, typename DDT::Id>(tile,ifile2);
    tile->set_id(tid);
    tile->finalize();
    return 0; // FIXME ?
}

template<typename DDT,typename WD>
std::ostream & write_ddt_stream(const DDT& ddt,  WD& wd, std::ostream & ofile, int tid,bool do_serialize,ddt::logging_stream & log)
{
    ofile << " 1 ";
    wd.write_serialized_stream(ofile);
    write_ddt_stream_core(ddt, ofile, tid,do_serialize, log);
    return ofile;
}


template<typename DDT>
std::ostream & write_ddt_stream(const DDT& ddt, std::ostream & ofile, int tid,bool do_serialize,ddt::logging_stream & log)
{
    ofile << " 0 ";
    write_ddt_stream_core(ddt, ofile, tid,do_serialize, log);
    return ofile;
}

template<typename DDT>
std::ostream & write_ddt_stream_core(const DDT& ddt, std::ostream & ofile, int tid,bool do_serialize,ddt::logging_stream & log)
{
    //std::cout << "write tile [id:" << tid << "]" << std::endl;
    typename DDT::Traits ttraits;
    auto tile  = ddt.get_const_tile(tid);
    tile->update_local_flag();
    auto ttri = tile->triangulation();
    ttraits.serialize_b64_cgal(ttri,ofile);
    serialize_b64_vect(tile->tile_ids,ofile);
    //tile->write_cgal(ofile);
    ddt::write_map_stream(tile->points_sent_,ofile,tile->current_dimension());
    ddt::write_json_stream(tile,ofile);
    return ofile; // FIXME ?
}

}
#endif

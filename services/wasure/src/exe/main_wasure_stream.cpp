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
#include <iostream>
#include <string>
#include <typeinfo>
#include <fstream>
#include <filesystem>

#include <libxml/tree.h>
#include <libxml/xmlwriter.h>

#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_xyz_points.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/read_las_points.h>
#include <CGAL/IO/write_las_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/jet_estimate_normals.h>

#include "io/stream_api.hpp"
#include "io/write_stream.hpp"
#include "io/write_vrt.hpp"
#include "io/read_stream.hpp"
#include "io/logging_stream.hpp"

#include "wasure_typedefs.hpp"
#include "write_geojson_wasure.hpp"
#include "wasure_data.hpp"
#include "wasure_algo.hpp"
#include "tbmrf_reco.hpp"
#include "io_ddt_stream.hpp"
#include "graph_cut.hpp"
#include "input_params.hpp"

#include "ddt_spark_utils.hpp"
#include <scanline_orient_normals.hpp>

typedef std::map<Id,wasure_data<Traits> > D_MAP;
typedef std::map<Id,std::list<wasure_data<Traits>> > D_LMAP;
typedef std::tuple<Id,double,double,double>                                SharedData;
typedef std::tuple<Id,Point,Point,Point,Point,double,double,double>           SharedDataDst;
typedef Id                                EdgeData;

int write_id_double_serialized(const std::map<Id,SharedData>  & lp, std::ostream & ofile,bool do_print = false)
{
    do_print = false;
    std::vector<double> outputv;
    for(auto pp : lp)
    {
        outputv.emplace_back(pp.first);
        outputv.emplace_back(std::get<0>(pp.second));
        outputv.emplace_back(std::get<1>(pp.second));
        outputv.emplace_back(std::get<2>(pp.second));
        outputv.emplace_back(std::get<3>(pp.second));
    }
    serialize_b64_vect(outputv,ofile);
    return 0;
}

int write_id_dst_serialized(const std::map<Id,SharedDataDst>  & lp, std::ostream & ofile,bool do_print = false)
{
    std::vector<double> outputv;
    for(auto pp : lp)
    {
        outputv.emplace_back(pp.first);
        outputv.emplace_back(std::get<0>(pp.second));
        for(int i = 0; i < 3; i++)
            outputv.emplace_back(std::get<1>(pp.second)[i]);
        for(int i = 0; i < 3; i++)
            outputv.emplace_back(std::get<2>(pp.second)[i]);
        for(int i = 0; i < 3; i++)
            outputv.emplace_back(std::get<3>(pp.second)[i]);
        for(int i = 0; i < 3; i++)
            outputv.emplace_back(std::get<4>(pp.second)[i]);
        outputv.emplace_back(std::get<5>(pp.second));
        outputv.emplace_back(std::get<6>(pp.second));
        outputv.emplace_back(std::get<7>(pp.second));
    }
    serialize_b64_vect(outputv,ofile);
    return 0;
}



std::istream & read_id_dst_serialized(std::map<Id,SharedDataDst> & lp, std::istream & ifile, bool do_print = false)
{
    Traits  traits;
    do_print = false;
    std::vector<double> input_v;
    deserialize_b64_vect(input_v,ifile);
    int nbe = 17;
    for(int n = 0; n < input_v.size()/nbe; n++)
    {
        Id id1 = input_v[n*nbe];
        double coords1[Traits::D];
        double coords2[Traits::D];
        double coords3[Traits::D];
        double coords4[Traits::D];
        for(int i = 0; i < 3; i++)
        {
            coords1[i] = input_v[n*nbe+2+i];
            coords2[i] = input_v[n*nbe+5+i];
            coords3[i] = input_v[n*nbe+8+i];
            coords4[i] = input_v[n*nbe+11+i];
        }
        auto pp1 = traits.make_point(coords1);
        auto pp2 = traits.make_point(coords2);
        auto pp3 = traits.make_point(coords3);
        auto pp4 = traits.make_point(coords4);
        lp[id1] = std::make_tuple(input_v[n*nbe+1],pp1,pp2,pp3,pp4,input_v[n*nbe+14],input_v[n*nbe+15],input_v[n*nbe+17]);
    }
    return ifile;
}



std::istream & read_id_double_serialized(std::map<Id,SharedData> & lp, std::istream & ifile, bool do_print = false)
{
    std::vector<double> input_v;
    deserialize_b64_vect(input_v,ifile);
    int nbe = 5;
    for(int n = 0; n< input_v.size()/nbe; n++)
    {
        Id id1 = input_v[n*nbe];
        lp[id1] = std::make_tuple(input_v[n*nbe+1],input_v[n*nbe+2],input_v[n*nbe+3],input_v[n*nbe+4]);
    }
    return ifile;
}


int write_edges_serialized(const std::map<Id,EdgeData>  & lp, std::ostream & ofile)
{
    std::vector<double> outputv;
    for(auto pp : lp)
    {
        outputv.emplace_back(pp.first);
        outputv.emplace_back(pp.second);
    }
    serialize_b64_vect(outputv,ofile);
    return 0;
}



std::istream & read_edges_serialized(std::map<Id,EdgeData> & lp, std::istream & ifile)
{
    std::vector<double> input_v;
    deserialize_b64_vect(input_v,ifile);
    int nbe = 2;
    for(int n = 0; n< input_v.size()/nbe; n++)
    {
        Id id1 = input_v[n*nbe];
        lp[id1] = input_v[n*nbe+1];
    }
    return ifile;
}



int simplify(Id tid,wasure_params & params,int nb_dat)
{
    std::cout.setstate(std::ios_base::failbit);
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        if(hpi.get_lab() == "p")
        {
            ddt_data<Traits> w_datas;
            w_datas.read_ply_stream(hpi.get_input_stream(),hpi.get_nl_char());
            hpi.finalize();
            std::cout.clear();
            Id id = hpi.get_id(0);
            ddt::stream_data_header oqh("p","s",id);
            std::string filename(params.output_dir + "/" + params.slabel +"_id_"+ std::to_string(tid) + "_" + std::to_string(id));
            if(params.dump_ply)
                oqh.write_into_file(filename,".ply");
            oqh.write_header(std::cout);
            w_datas.write_ply_stream(oqh.get_output_stream(),oqh.get_nl_char());
            oqh.finalize();
            std::cout << std::endl;
        }
    }
    return 0;
}

template <typename T>
auto normalize(T const& V)
{
    auto const slen = V.squared_length();
    auto const d = CGAL::approximate_sqrt(slen);
    return V / d;
}


int dim_splitted(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    wasure_algo w_algo;
    int D = Traits::D;
    Traits  traits;
    D_LMAP w_datas_map;
    bool do_splitted =  false;
    wasure_data<Traits>  w_datas_full;
    std::vector<Point> p_simp_full;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        log.step("read");
        if(hpi.get_lab() == "z" )
        {
            w_datas_map[hid].push_back(wasure_data<Traits>());
            if(hpi.is_file())
            {
                w_datas_map[hid].back().read_ply_stream(hpi.get_input_stream(),hpi.get_nl_char());
                w_datas_map[hid].back().shpt2uint8();
            }
            else
            {
                w_datas_map[hid].back().read_serialized_stream(hpi.get_input_stream());
            }
        }
        wasure_data<Traits> & w_datas = w_datas_map[hid].back();
        hpi.finalize();
        w_datas.dmap[w_datas.xyz_name].extract_full_uint8_vect(w_datas.format_points,false);
        w_datas.dmap[w_datas.center_name].extract_full_uint8_vect(w_datas.format_centers,false);
        w_datas.dmap[w_datas.normal_name].extract_full_uint8_vect(w_datas.format_normals,false);
        w_datas.extract_flags(w_datas.format_flags,false);
        std::vector<Point> p_simp;
        log.step("compute_dim");
        if(w_datas.format_centers.size() == 0)
        {
            double coords[Traits::D];
            for(auto pp : w_datas.format_points)
            {
                for(int d = 0; d < D; d++)
                {
                    if(d < D-1)
                        coords[d] = pp[d];
                    else
                        coords[d] = pp[d] + 30;
                }
                w_datas.format_centers.push_back(traits.make_point(coords));
            }
        }
        if(do_splitted)
        {
            w_algo.compute_dim(w_datas.format_points,
                               w_datas.format_egv,
                               w_datas.format_sigs,log);
            w_algo.flip_dim_ori(w_datas.format_points,
                                w_datas.format_egv,
                                w_datas.format_centers);
            w_datas_full.format_points.insert(w_datas_full.format_points.end(),w_datas.format_points.begin(),w_datas.format_points.end());
            w_datas_full.format_egv.insert(w_datas_full.format_egv.end(),w_datas.format_egv.begin(),w_datas.format_egv.end());
            w_datas_full.format_sigs.insert(w_datas_full.format_sigs.end(),w_datas.format_sigs.begin(),w_datas.format_sigs.end());
            p_simp_full.insert(p_simp_full.end(),p_simp.begin(),p_simp.end());
        }
        else
        {
            if(i == 0)
            {
                w_datas_full = w_datas;
            }
            else
            {
                w_datas_full.insert(w_datas);
            }
        }
    }
    if(!do_splitted)
    {
        w_datas_full.dmap[w_datas_full.xyz_name].extract_full_uint8_vect(w_datas_full.format_points,false);
        w_datas_full.dmap[w_datas_full.center_name].extract_full_uint8_vect(w_datas_full.format_centers,false);
        w_datas_full.dmap[w_datas_full.normal_name].extract_full_uint8_vect(w_datas_full.format_normals,false);
        w_algo.compute_dim(w_datas_full.format_points,
                           w_datas_full.format_egv,
                           w_datas_full.format_sigs,log);
        for(auto pp : w_datas_full.format_normals)
        {
            break;
            CGAL::Vector_3<typename Traits::K> v1{pp[0],pp[1],pp[2]};
            CGAL::Vector_3<typename Traits::K> u = {-v1[1], v1[0], 0.0};
            auto v2 = normalize(CGAL::cross_product(v1,u));
        }
        if(w_datas_full.format_centers.size() == 0)
        {
            double coords[Traits::D];
            for(auto pp : w_datas_full.format_points)
            {
                for(int d = 0; d < D; d++)
                {
                    if(d < D-1)
                        coords[d] = pp[d];
                    else
                        coords[d] = pp[d] + 30;
                }
                w_datas_full.format_centers.push_back(traits.make_point(coords));
            }
        }
        w_algo.flip_dim_ori(w_datas_full.format_points,
                            w_datas_full.format_egv,
                            w_datas_full.format_centers);
    }
    if(params.pscale < 1)
    {

        w_algo.tessel_adapt(w_datas_full.format_points,
                            p_simp_full,
                            w_datas_full.format_egv,
                            w_datas_full.format_sigs,
                            w_datas_full.format_normals,
                            NUM_ITER_TESSEL,params.pscale,D,tid
                           );
    }
    else
    {
        p_simp_full.insert(p_simp_full.end(),w_datas_full.format_points.begin(),w_datas_full.format_points.end());
    }
    log.step("dump");
    if(do_splitted)
    {
        for ( auto it = w_datas_map.begin(); it != w_datas_map.end(); it++ )
        {
            int acc = 0;
            for(auto & w_datas : w_datas_map[it->first])
            {
                w_datas.fill_egv(w_datas.format_egv);
                w_datas.fill_sigs(w_datas.format_sigs);
                std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(it->first) + "_" + std::to_string(acc++) +  "_dim");
                std::cout.clear();
                log.step("write");
                ddt::stream_data_header oth("z","s",tid);
                oth.write_header(std::cout);
                if(params.dump_ply)
                    w_datas.write_ply_stream(oth.get_output_stream(),oth.get_nl_char(),true);
                else
                    w_datas.write_serialized_stream(oth.get_output_stream());
                oth.finalize();
                std::cout << std::endl;
            }
        }
    }
    else
    {
        w_datas_full.fill_egv(w_datas_full.format_egv);
        w_datas_full.fill_sigs(w_datas_full.format_sigs);
        w_datas_full.dmap[w_datas_full.xyz_name].fill_full_uint8_vect(w_datas_full.format_points);
        std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(tid) +  "_dim");
        std::cout.clear();
        log.step("write");
        ddt::stream_data_header oth("z","s",tid);
        if(params.dump_ply)
            oth.write_into_file(ply_name,".ply");
        oth.write_header(std::cout);
        if(params.dump_ply)
            w_datas_full.write_ply_stream(oth.get_output_stream(),oth.get_nl_char(),true);
        else
            w_datas_full.write_serialized_stream(oth.get_output_stream());
        oth.finalize();
        std::cout << std::endl;
    }
    bool do_debug = false;
    if(p_simp_full.size() > 0)
    {
        ddt::stream_data_header oxh("x","z",tid);
        ddt_data<Traits> datas_out;
        datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(p_simp_full);
        std::string ply_name(params.output_dir +  "/simp_id_" + std::to_string(tid) + "_simp");
        if(!do_debug)
        {
            if(params.dump_ply)
                oxh.write_into_file(ply_name,".ply");
            oxh.write_header(std::cout);
            if(params.dump_ply)
                datas_out.write_ply_stream(oxh.get_output_stream(),oxh.get_nl_char());
            else
                datas_out.write_serialized_stream(oxh.get_output_stream());
        }
        else
        {
            datas_out.write_ply_stream(oxh.get_output_stream(),oxh.get_nl_char());
        }
        oxh.finalize();
        std::cout << std::endl;
    }
    return 0;
}




int dst(const Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_LMAP w_datas_pts;
    D_MAP w_datas_tri;
    std::map<Id,ddt_data<Traits> > d_datas_tri;
    auto w_data_full = wasure_data<Traits>();
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri, hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
        }
        if(hpi.get_lab() == "s")
        {
            std::vector<int> vv(3);
            for(int d = 0; d < 3; d++)
            {
                hpi.get_input_stream() >> vv[d];
            }
            w_datas_tri[hid].tile_ids = vv;
        }
        if(hpi.get_lab() == "z")
        {
            w_datas_pts[hid].push_back(wasure_data<Traits>());
            if(hpi.is_file())
            {
                w_datas_pts[hid].back().read_ply_stream(hpi.get_input_stream(),hpi.get_nl_char());
                w_datas_pts[hid].back().shpt2uint8();
            }
            else
            {
                w_datas_pts[hid].back().read_serialized_stream(hpi.get_input_stream());
            }
        }
        hpi.finalize();
    }
    log.step("preprocess");
    for ( auto it = w_datas_pts.begin(); it != w_datas_pts.end(); it++ )
    {
        for(auto & wpt : w_datas_pts[it->first])
        {
            wpt.dmap[wpt.xyz_name].extract_full_uint8_vect(wpt.format_points,false);
            w_data_full.format_points.insert(w_data_full.format_points.end(),wpt.format_points.begin(),wpt.format_points.end());
            wpt.dmap[wpt.center_name].extract_full_uint8_vect(wpt.format_centers,false);
            wpt.dmap[wpt.normal_name].extract_full_uint8_vect(wpt.format_normals,false);
            w_data_full.format_centers.insert(w_data_full.format_centers.end(),wpt.format_centers.begin(),wpt.format_centers.end());
            w_data_full.format_normals.insert(w_data_full.format_normals.end(),wpt.format_normals.begin(),wpt.format_normals.end());
            wpt.extract_sigs(w_data_full.format_sigs,false);
            wpt.extract_egv(w_data_full.format_egv,false);
        }
    }
    Tile_iterator  tile1  = tri.get_tile(tid);
    int nbs = tile1->number_of_cells();
    if(w_data_full.format_flags.size() == 0)
    {
        for(int ss = 0; ss < nbs ; ss++)
        {
            w_data_full.format_flags.push_back(0);
        }
    }
    std::vector<std::vector<double>>  & format_dst = w_datas_tri[tid].format_dst; ;
    if(format_dst.size() == 0)
    {
        for(int ss = 0; ss < nbs ; ss++)
        {
            format_dst.push_back(std::vector<double>({0.0,0.0,1.0}));
        }
    }
    log.step("compute");
    DT & tri_tile  = tri.get_tile(tid)->triangulation();
    w_algo.compute_dst_with_center(tri,w_datas_tri[tid],w_data_full,params,tid);
    log.step("finalize");
    w_datas_tri[tid].fill_dst(w_datas_tri[tid].format_dst);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    return 0;
}


int dst_conflict(const Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_LMAP w_datas_pts;
    D_MAP w_datas_tri;
    std::map<Id,ddt_data<Traits> > d_datas_tri;
    auto w_data_full = wasure_data<Traits>();
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
        }
        if(hpi.get_lab() == "s")
        {
            std::vector<int> vv(3);
            for(int d = 0; d < 3; d++)
            {
                hpi.get_input_stream() >> vv[d];
            }
            w_datas_tri[hid].tile_ids = vv;
        }
        if(hpi.get_lab() == "z")
        {
            w_datas_pts[hid].push_back(wasure_data<Traits>());
            if(hpi.is_file())
            {
                w_datas_pts[hid].back().read_ply_stream(hpi.get_input_stream(),hpi.get_nl_char());
                w_datas_pts[hid].back().shpt2uint8();
            }
            else
            {
                w_datas_pts[hid].back().read_serialized_stream(hpi.get_input_stream());
            }
        }
        hpi.finalize();
    }
    log.step("preprocess");
    for ( auto it = w_datas_pts.begin(); it != w_datas_pts.end(); it++ )
    {
        for(auto & wpt : w_datas_pts[it->first])
        {
            wpt.dmap[w_data_full.xyz_name].extract_full_uint8_vect(w_data_full.format_points,false);
            wpt.extract_egv(wpt.format_egv,false);
            wpt.extract_sigs(wpt.format_sigs,false);
            wpt.dmap[w_data_full.center_name].extract_full_uint8_vect(w_data_full.format_centers,false);
        }
    }
    std::vector<std::vector<double>>  & format_dst = w_datas_tri[tid].format_dst; ;
    int nbs = w_datas_tri[tid].nb_simplex_uint8_vect();
    if(format_dst.size() == 0)
    {
        for(int ss = 0; ss < nbs ; ss++)
        {
            format_dst.push_back(std::vector<double>({0.0,0.0,1.0}));
        }
    }
    log.step("compute");
    DT & tri_tile  = tri.get_tile(tid)->triangulation();
    std::list<wasure_data<Traits>> l_tri;
    for ( auto it = w_datas_pts.begin(); it != w_datas_pts.end(); it++ )
    {
        for(auto & wpt : w_datas_pts[it->first])
        {
            l_tri.push_back(wasure_data<Traits>());
            auto & wd_tri = l_tri.back();
            std::vector<std::vector<double>>  & format_dst_pts = wd_tri.format_dst; ;
            if(format_dst_pts.size() == 0)
            {
                for(int ss = 0; ss < nbs ; ss++)
                {
                    format_dst_pts.push_back(std::vector<double>({0.0,0.0,1.0}));
                }
            }
            w_algo.compute_dst_with_center(tri,wd_tri,wpt,params,tid);
            for(int ss = 0; ss < nbs ; ss++)
            {
                double & vpe = format_dst_pts[ss][0];
                double & vpo = format_dst_pts[ss][1];
                double & vpu = format_dst_pts[ss][2];
                if(vpe > vpo)
                {
                    vpe = vpe-vpo;
                    vpo = 0;
                    if(vpe > 0.95)
                        vpe = 0.95;
                    vpu =1-vpe;
                }
                else
                {
                    vpo = vpo-vpe;
                    vpe = 0;
                    if(vpo > 0.95)
                        vpo = 0.95;
                    vpu = 1-vpo;
                }
            }
        }
    }
    int acc = 0;
    for(auto & wd_tri : l_tri)
    {
        std::cout << "dst_acc:" << acc << std::endl;
        std::vector<std::vector<double>>  & format_dst_pts = wd_tri.format_dst;
        for(int ss = 0; ss < nbs ; ss++)
        {
            double  vpe = format_dst[ss][0];
            double  vpo = format_dst[ss][1];
            double  vpu = format_dst[ss][2];
            w_algo.ds_score(vpe,vpo,vpu,
                            format_dst_pts[ss][0],format_dst_pts[ss][1],format_dst_pts[ss][2],
                            format_dst[ss][0],format_dst[ss][1],format_dst[ss][2]);
            regularize(format_dst[ss][0],format_dst[ss][1],format_dst[ss][2]);
            acc++;
        }
    }
    log.step("finalize");
    w_datas_tri[tid].fill_dst(w_datas_tri[tid].format_dst);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    return 0;
}



int dst_good(const Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_LMAP w_datas_pts;
    D_MAP w_datas_tri;
    std::map<Id,ddt_data<Traits> > d_datas_tri;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
        }
        if(hpi.get_lab() == "s")
        {
            std::vector<int> vv(3);
            for(int d = 0; d < 3; d++)
            {
                hpi.get_input_stream() >> vv[d];
            }
            w_datas_tri[hid].tile_ids = vv;
        }
        if(hpi.get_lab() == "z")
        {
            w_datas_pts[hid].push_back(wasure_data<Traits>());
            if(hpi.is_file())
            {
                w_datas_pts[hid].back().read_ply_stream(hpi.get_input_stream(),hpi.get_nl_char());
                w_datas_pts[hid].back().shpt2uint8();
            }
            else
            {
                w_datas_pts[hid].back().read_serialized_stream(hpi.get_input_stream());
            }
        }
        hpi.finalize();
    }
    log.step("preprocess");
    for ( auto it = w_datas_pts.begin(); it != w_datas_pts.end(); it++ )
    {
        for(auto & wpt : w_datas_pts[it->first])
        {
            wpt.dmap[wpt.xyz_name].extract_full_uint8_vect(wpt.format_points,false);
            wpt.extract_egv(wpt.format_egv,false);
            wpt.extract_sigs(wpt.format_sigs,false);
            wpt.dmap[wpt.center_name].extract_full_uint8_vect(wpt.format_centers,false);
        }
    }
    std::vector<std::vector<double>>  & format_dst = w_datas_tri[tid].format_dst; ;
    if(format_dst.size() == 0)
    {
        int nbs = w_datas_tri[tid].nb_simplex_uint8_vect();
        for(int ss = 0; ss < nbs ; ss++)
        {
            format_dst.push_back(std::vector<double>({0.0,0.0,1.0}));
        }
    }
    log.step("compute");
    DT & tri_tile  = tri.get_tile(tid)->triangulation();
    for(auto wpt = w_datas_pts[tid].begin() ; wpt != w_datas_pts[tid].end() ; wpt++)
    {
        w_algo.compute_dst_with_center(tri,w_datas_tri[tid],*wpt,params,tid);
    }
    log.step("finalize");
    w_datas_tri[tid].fill_dst(w_datas_tri[tid].format_dst);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    return 0;
}




int regularize_slave_focal(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    Tile_iterator  tile_k  = tri.get_tile(tid);
    for(auto cit = tile_k->cells_begin();
            cit != tile_k->cells_end();
            cit++)
    {
        if(!tile_k->cell_is_mixed(cit) || tile_k->cell_is_main(cit))
            continue;
        Id main_tid = tile_k->cell_main_id(cit);
        Tile_iterator main_tile = tri.get_tile(main_tid);
        Id lid1 = tile_k->lid(cit);
        auto main_cell = main_tile->locate_cell(*tile_k,cit);
        Id lid2 = main_tile->lid(main_cell);
        w_datas_tri[tid].replace_attribute(w_datas_tri[main_tid],lid1,lid2);
    }
    Id tid_k = tid;
    Tile_const_iterator  tilec_k  = tri.get_const_tile(tid);
    std::map<Id,std::map<Id,SharedData> > shared_data_map;
    // Loop over each shared cell to extracts id relation
    // lid_l , lid_l <-> lid_k
    for( auto cit_k = tilec_k->cells_begin();
            cit_k != tilec_k->cells_end(); ++cit_k )
    {
        Cell_const_iterator fch = Cell_const_iterator(tilec_k,tilec_k, tilec_k, cit_k);
        if(!tile_k->cell_is_mixed(cit_k) || tile_k->cell_is_infinite(cit_k))
            continue;
        Id lid_k = tile_k->lid(cit_k);
        int D = tile_k->current_dimension();
        // Loop on shared tiles
        std::unordered_set<Id> idSet ;
        for(int i=0; i<=D; ++i)
        {
            // Select only id != tid_k and new ids
            Id tid_l = tile_k->id(tile_k->vertex(cit_k,i));
            if ((idSet.find(tid_l) != idSet.end()) || tid_l == tid_k )
            {
                continue;
            }
            idSet.insert(tid_l);
            if(tid_l == tid_k)
                continue;
            Tile_iterator tile_l = tri.get_tile(tid_l);
            auto cit_l = tile_l->locate_cell(*tile_k,cit_k);
            Id lid_l = tile_l->lid(cit_l);
            if(shared_data_map.find(tid_l) == shared_data_map.end())
                shared_data_map[tid_l] = std::map<Id,SharedData>();
            // The current data structure => local_id of the shared tet, lag and tau
            shared_data_map[tid_l][lid_k] = std::make_tuple(lid_l,0,1,0);
        }
    }
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    for(auto ee : shared_data_map)
    {
        Id tid2 = ee.first;
        if(tid2 == tid)
            continue;
        ddt::stream_data_header hto("e","z",std::vector<int> {tid,tid2});
        std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid) + "_nid" + std::to_string(tid2));
        hto.write_header(std::cout);
        write_id_double_serialized(ee.second,hto.get_output_stream());
        hto.finalize();
        std::cout << std::endl;
    }
    return 0;
}



int regularize_slave_extract(Id tid_1,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    std::map<Id,std::map<Id,SharedDataDst> > edges_dst_map;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    Tile_const_iterator  tilec_1  = tri.get_const_tile(tid_1);
    std::vector<std::vector<double>> & v_dst = w_datas_tri[tid_1].format_dst;
    // Loop over each shared cell to extracts id relation
    // lid_2 , lid_2 <-> lid_1
    for( auto cit_1 = tilec_1->cells_begin();
            cit_1 != tilec_1->cells_end(); ++cit_1 )
    {
        Cell_const_iterator fch = Cell_const_iterator(tilec_1,tilec_1, tilec_1, cit_1);
        if(!tilec_1->cell_is_mixed(cit_1) || tilec_1->cell_is_infinite(cit_1))
            continue;
        Id lid_1 = tilec_1->lid(cit_1);
        int D = tilec_1->current_dimension();
        // Loop on shared tiles
        std::unordered_set<Id> idSet ;
        for(int i=0; i<=D; ++i)
        {
            // Select only id != tid_1 and new ids
            Id tid_2 = tilec_1->id(tilec_1->vertex(cit_1,i));
            if ((idSet.find(tid_2) != idSet.end()) || tid_2 == tid_1 )
            {
                continue;
            }
            idSet.insert(tid_2);
            if(edges_dst_map.find(tid_2) == edges_dst_map.end())
                edges_dst_map[tid_2] = std::map<Id,SharedDataDst>();
            const Point& a = tilec_1->vertex(cit_1,0)->point();
            const Point& b = tilec_1->vertex(cit_1,1)->point();
            const Point& c = tilec_1->vertex(cit_1,2)->point();
            const Point& d = tilec_1->vertex(cit_1,3)->point();
            // The current data structure => local_id of the shared tet, lag and tau
            edges_dst_map[tid_2][lid_1] = std::make_tuple(lid_1,a,b,c,d,v_dst[lid_1][0],v_dst[lid_1][1],v_dst[lid_1][2]);
        }
    }
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid_1);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid_1], oth.get_output_stream(),tid_1,false,log);
    oth.finalize();
    std::cout << std::endl;
    for(auto ee : edges_dst_map)
    {
        Id tid_2 = ee.first;
        ddt::stream_data_header hto("f","z",std::vector<int> {tid_1,tid_2});
        std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1) + "_nid" + std::to_string(tid_2));
        hto.write_header(std::cout);
        write_id_dst_serialized(ee.second,hto.get_output_stream(),tid_1 ==3);
        hto.finalize();
        std::cout << std::endl;
    }
    return 0;
}




int regularize_slave_insert(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    int D = Traits::D;
    D_MAP w_datas_tri;
    Traits  traits;
    std::map<Id,std::map<Id,SharedDataDst> > edges_dst_map;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,true);
        }
        if(hpi.get_lab() == "f")
        {
            Id eid1 = hpi.get_id(0);
            edges_dst_map[eid1] = std::map<Id,SharedDataDst>();
            std::map<Id,SharedDataDst>  & lp = edges_dst_map[eid1];
            read_id_dst_serialized(lp, hpi.get_input_stream(),tid == 3);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    Id tid_k = tid;
    Tile_const_iterator  tile_k  = tri.get_const_tile(tid);
    std::map<Id,std::map<Id,SharedData> > shared_data_map;
    std::vector<std::vector<double>>  & format_dst = w_datas_tri[tid].format_dst; ;
    for(auto ee : edges_dst_map)
    {
        Id tid_l = ee.first;
        for(auto ee_map : ee.second)
        {
            auto edm_l = ee_map.second;
            Id lid_l = std::get<0>(edm_l);
            std::vector<double> dstv{std::get<5>(edm_l),std::get<6>(edm_l),std::get<7>(edm_l)};
            std::vector<Point> vpp
            {
                std::get<1>(edm_l),
                std::get<2>(edm_l),
                std::get<3>(edm_l),
                std::get<4>(edm_l)
            };
            auto pp1 = std::get<1>(edm_l);
            auto pp2 = std::get<2>(edm_l);
            auto pp3 = std::get<3>(edm_l);
            auto pp4 = std::get<4>(edm_l);
            std::vector<double> bary
            {
                (pp1[0]+pp2[0]+pp3[0]+pp4[0])/4,
                (pp1[1]+pp2[1]+pp3[1]+pp4[1])/4,
                (pp1[2]+pp2[2]+pp3[2]+pp4[2])/4
            };
            auto pp_bary = traits.make_point(bary.begin());
            auto main_cell = tile_k->locate_cell_point(*tile_k,pp_bary);
            bool is_cell_equal = true;
            for(int d1 = 0; d1 <= D; d1++)
            {
                bool is_pts_equal = false;
                auto pd1 = vpp[d1];
                for(int d2 = 0; d2 <= D; d2++)
                {
                    auto pd2 = main_cell->vertex(d2)->point();
                    is_pts_equal = (is_pts_equal || pd1 == pd2);
                }
                is_cell_equal = is_cell_equal && is_pts_equal;
            }
            Id cmid = tile_k->cell_main_id(main_cell);
            Id lid_k = tile_k->lid(main_cell);
            if(shared_data_map.find(tid_l) == shared_data_map.end())
                shared_data_map[tid_l] = std::map<Id,SharedData>();
            shared_data_map[tid_l][lid_k] = std::make_tuple(lid_l,0,1,0);
            if(cmid == tid_l)
            {
                for(int i = 0 ; i < 3; i++)
                    format_dst[lid_k][i] = dstv[i];
            }
        }
    }
    w_datas_tri[tid].fill_dst(w_datas_tri[tid].format_dst);
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    // Dum edges
    for(auto ee : shared_data_map)
    {
        Id tid2 = ee.first;
        if(tid2 == tid)
            continue;
        ddt::stream_data_header hto("e","z",std::vector<int> {tid,tid2});
        std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid) + "_nid" + std::to_string(tid2));
        hto.write_header(std::cout);
        write_id_double_serialized(ee.second,hto.get_output_stream());
        hto.finalize();
        std::cout << std::endl;
    }
    return 0;
}

int extract_surface(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::vector<Facet_const_iterator> lft;
    std::vector<bool> lbool;
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    int D = Traits::D;
    D_MAP w_datas_tri;
    ddt_data<Traits> datas_out;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t" || hpi.get_lab() == "u")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    log.step("compute");
    int mode = params.mode;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        try
        {
            if(fit->main_id() != tid || fit->is_infinite())
                continue;
            Cell_const_iterator tmp_fch = fit.full_cell();
            int tmp_idx = fit.index_of_covertex();
            Cell_const_iterator tmp_fchn = tmp_fch->neighbor(tmp_idx);
            if(!tri.tile_is_loaded(tmp_fch->main_id()) ||
                    !tri.tile_is_loaded(tmp_fchn->main_id()))
            {
                std::cerr << "ERROR tile not loaded" << std::endl;
                continue;
            }
            bool is_on_convex = false;
            if(tmp_fch->is_infinite() ||  tmp_fchn->is_infinite() )
                is_on_convex = true;
            Cell_const_iterator fch = tmp_fch->main();
            int id_cov = fit.index_of_covertex();
            Cell_const_iterator fchn = tmp_fchn->main();
            int cccid = fch->lid();
            int cccidn = fchn->lid();
            int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
            int chnlab = w_datas_tri[fchn->tile()->id()].format_labs[cccidn];
            if(
                (ch1lab != chnlab)  || (mode == 0 && (is_on_convex && (ch1lab == 0 || chnlab == 0)))
            )
            {
                lft.push_back(*fit);
                const Point& a = fch->vertex((id_cov+1)&3)->point();
                const Point& b = fch->vertex((id_cov+2)&3)->point();
                const Point& c = fch->vertex((id_cov+3)&3)->point();
                const Point& d = fch->vertex((id_cov)&3)->point();
                datas_out.bbox += a;
                datas_out.bbox += b;
                datas_out.bbox += c;
                bool bl =
                    (CGAL::orientation(a,b,c,d) == 1 && chnlab == 0) ||
                    (CGAL::orientation(a,b,c,d) == -1 && chnlab == 1);
                lbool.push_back(!bl);
            }
        }
        catch (ddt::DDT_exeption& e)
        {
            std::cerr << "!! WARNING !!!" << std::endl;
            std::cerr << "Exception catched : " << e.what() << std::endl;
            continue;
        }
    }
    log.step("write");
    std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(tid) + "_surface");
    std::cout.clear();
    ddt::stream_data_header oth("p","z",tid);
    if(D == 2)
    {
        oth.write_into_file(ply_name,".geojson");
        oth.write_header(std::cout);
    }
    switch (D)
    {
    case 2 :
    {
        wasure::dump_2d_surface_geojson<DTW>(lft,oth.get_output_stream());
        break;
    }
    case 3 :
    {
        std::vector<Point>  format_points;
        std::vector<int> v_simplex;
        std::map<Vertex_const_iterator, uint> vertex_map;
        int acc = 0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            for(int i = 0; i < D+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    if(vertex_map.find(v) == vertex_map.end())
                    {
                        vertex_map[v] = acc++;
                        format_points.push_back(v->point());
                    }
                }
            }
        }
        acc=0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            Id ida = (id_cov+1)&3;
            Id idb = (id_cov+2)&3;
            Id idc = (id_cov+3)&3;
            v_simplex.push_back(vertex_map[fch->vertex(ida)]);
            if(!lbool[acc])
            {
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
            }
            else
            {
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
            }
            acc++;
        }
        datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",D,D,tinyply::Type::INT32);
        datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
        datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
        datas_out.write_ply_stream(oth.get_output_stream(),oth.get_nl_char(),false,false,true);
        break;
    }
    default :
    {
        return 1;
        break;
    }
    }
    oth.finalize();
    std::cout << std::endl;
    return 0;
}


int extract_surface_new(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::vector<Facet_const_iterator> lft;
    std::vector<bool> lbool;
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    int D = Traits::D;
    D_MAP w_datas_tri;
    ddt_data<Traits> datas_out;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t" || hpi.get_lab() == "u")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    log.step("compute");
    int mode = params.mode;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        try
        {
            if(fit->is_infinite())
                continue;
            Cell_const_iterator tmp_fch = fit.full_cell();
            int tmp_idx = fit.index_of_covertex();
            Cell_const_iterator tmp_fchn = tmp_fch->neighbor(tmp_idx);
            bool is_on_convex = false;
            if(tmp_fch->is_infinite() ||  tmp_fchn->is_infinite() )
            {
                is_on_convex = true;
                if(mode == 1)
                    continue;
            }
            Cell_const_iterator fch = tmp_fch;
            int id_cov = fit.index_of_covertex();
            Cell_const_iterator fchn = tmp_fchn;
            int cccid = fch->lid();
            int cccidn = fchn->lid();
            int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
            int chnlab = w_datas_tri[fchn->tile()->id()].format_labs[cccidn];
            if(
                (ch1lab != chnlab)  || (mode == 0 && (is_on_convex && (ch1lab == 0 || chnlab == 0)))
            )
            {
                const Point& a = fch->vertex((id_cov+1)&3)->point();
                const Point& b = fch->vertex((id_cov+2)&3)->point();
                const Point& c = fch->vertex((id_cov+3)&3)->point();
                const Point& d = fch->vertex((id_cov)&3)->point();
                bool skip  = false;
                double dist1 = CGAL::approximate_sqrt(CGAL::squared_distance(fch->vertex((id_cov+1)&3)->point(),fch->vertex((id_cov+2)&3)->point()));
                double dist2 = CGAL::approximate_sqrt(CGAL::squared_distance(fch->vertex((id_cov+2)&3)->point(),fch->vertex((id_cov+3)&3)->point()));
                double dist3 = CGAL::approximate_sqrt(CGAL::squared_distance(fch->vertex((id_cov+1)&3)->point(),fch->vertex((id_cov+3)&3)->point()));
                if(dist1 > MAX_TRIANGLE_SIZE ||
                        dist2 > MAX_TRIANGLE_SIZE ||
                        dist3 > MAX_TRIANGLE_SIZE )
                {
                    skip = true;
                    continue;
                }
                lft.push_back(*fit);
                datas_out.bbox += a;
                datas_out.bbox += b;
                datas_out.bbox += c;
                bool bl =
                    (CGAL::orientation(a,b,c,d) == 1 && chnlab == 0) ||
                    (CGAL::orientation(a,b,c,d) == -1 && chnlab == 1);
                lbool.push_back(!bl);
            }
        }
        catch (ddt::DDT_exeption& e)
        {
            std::cerr << "!! WARNING !!!" << std::endl;
            std::cerr << "Exception catched : " << e.what() << std::endl;
            continue;
        }
    }
    log.step("write");
    std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(tid) + "_surface");
    std::cout.clear();
    ddt::stream_data_header oth("p","z",tid);
    if(D == 2)
    {
        oth.write_into_file(ply_name,".geojson");
        oth.write_header(std::cout);
    }
    switch (D)
    {
    case 2 :
    {
        wasure::dump_2d_surface_geojson<DTW>(lft,oth.get_output_stream());
        break;
    }
    case 3 :
    {
        std::vector<Point>  format_points;
        std::vector<int> v_simplex;
        std::map<Vertex_const_iterator, uint> vertex_map;
        int acc = 0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            for(int i = 0; i < D+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    if(vertex_map.find(v) == vertex_map.end())
                    {
                        vertex_map[v] = acc++;
                        format_points.push_back(v->point());
                    }
                }
            }
        }
        acc=0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            Id ida = (id_cov+1)&3;
            Id idb = (id_cov+2)&3;
            Id idc = (id_cov+3)&3;
            v_simplex.push_back(vertex_map[fch->vertex(ida)]);
            if(!lbool[acc])
            {
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
            }
            else
            {
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
            }
            acc++;
        }
        datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",D,D,tinyply::Type::INT32);
        datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
        datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
        datas_out.write_ply_stream(oth.get_output_stream(),oth.get_nl_char(),false,false,true);
        break;
    }
    default :
    {
        return 1;
        break;
    }
    }
    oth.finalize();
    std::cout << std::endl;
    return 0;
}





int extract_surface_area(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::vector<Facet_const_iterator> lft;
    std::vector<bool> lbool;
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    int D = Traits::D;
    D_MAP w_datas_tri;
    log.step("read");
    std::vector<Id> lid;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        lid.push_back(hid);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    log.step("compute");
    int mode = params.mode;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        if(fit->is_infinite() || (fit->main_id() != lid[0] &&  fit->main_id() != lid[1]))
            continue;
        Cell_const_iterator tmp_fch = fit.full_cell();
        int tmp_idx = fit.index_of_covertex();
        Cell_const_iterator tmp_fchn = tmp_fch->neighbor(tmp_idx);
        bool do_debug = false;
        for(int i=0; i<=D; ++i)
        {
            auto pp = tmp_fch->vertex(i)->point();
            if(std::fabs(pp[0] - (-11.15)) < 0.02 &&
                    std::fabs(pp[1] - 0.26    )  < 0.02  &&
                    std::fabs(pp[2] - 0.37    ) < 0.02)
            {
                do_debug = true;
            }
        }
        try
        {
            if(params.area_processed == 2  )
            {
                bool is_edge =  false;
                if((tmp_fch->has_id(lid[0]) && tmp_fchn->has_id(lid[1])) ||
                        (tmp_fch->has_id(lid[1]) && tmp_fchn->has_id(lid[0]))
                  )
                {
                    is_edge = true;
                }
                if(!is_edge )
                    continue;
            }
            if(params.area_processed == 1)
            {
                if(tmp_fch->is_mixed() || tmp_fchn->is_mixed())
                    continue;
            }
            bool do_keep_local = false;
            Cell_const_iterator fch = tmp_fch;
            int id_cov = fit.index_of_covertex();
            Cell_const_iterator fchn = tmp_fchn;
            int cccid = fch->lid();
            int cccidn = fchn->lid();
            int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
            int chnlab = w_datas_tri[fchn->tile()->id()].format_labs[cccidn];
            if(
                (ch1lab != chnlab)
            )
            {
                lft.push_back(*fit);
                const Point& a = fch->vertex((id_cov+1)&3)->point();
                const Point& b = fch->vertex((id_cov+2)&3)->point();
                const Point& c = fch->vertex((id_cov+3)&3)->point();
                const Point& d = fch->vertex((id_cov)&3)->point();
                bool bl =
                    (CGAL::orientation(a,b,c,d) == 1 && chnlab == 0) ||
                    (CGAL::orientation(a,b,c,d) == -1 && chnlab == 1);
                lbool.push_back(!bl);
            }
        }
        catch (ddt::DDT_exeption& e)
        {
            if(do_debug)
            {
                std::cerr << "!! WARNING !!!" << std::endl;
                std::cerr << "params.area_processed : " << params.area_processed << std::endl;
                std::cerr << "Exception catched : " << e.what() << std::endl;
            }
            continue;
        }
    }
    log.step("write");
    std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(tid) + "_surface");
    std::cout.clear();
    ddt::stream_data_header oth("p","z",tid);
    if(D == 2)
    {
        oth.write_into_file(ply_name,".geojson");
        oth.write_header(std::cout);
    }
    if(lft.size() == 0)
        return 0;
    switch (D)
    {
    case 2 :
    {
        wasure::dump_2d_surface_geojson<DTW>(lft,oth.get_output_stream());
        break;
    }
    case 3 :
    {
        std::vector<Point>  format_points;
        std::vector<int> v_simplex;
        std::map<Vertex_const_iterator, uint> vertex_map;
        ddt_data<Traits> datas_out;
        int acc = 0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            for(int i = 0; i < D+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    if(vertex_map.find(v) == vertex_map.end())
                    {
                        vertex_map[v] = acc++;
                        format_points.push_back(v->point());
                    }
                }
            }
        }
        acc=0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            Id ida = (id_cov+1)&3;
            Id idb = (id_cov+2)&3;
            Id idc = (id_cov+3)&3;
            const Point& a = fch->vertex((id_cov+1)&3)->point();
            const Point& b = fch->vertex((id_cov+2)&3)->point();
            const Point& c = fch->vertex((id_cov+3)&3)->point();
            const Point& d = fch->vertex((id_cov)&3)->point();
            datas_out.bbox += a;
            datas_out.bbox += b;
            datas_out.bbox += c;
            v_simplex.push_back(vertex_map[fch->vertex(ida)]);
            if(!lbool[acc])
            {
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
            }
            else
            {
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
            }
            acc++;
        }
        datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",D,D,tinyply::Type::INT32);
        datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
        datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
        datas_out.write_ply_stream(oth.get_output_stream(),oth.get_nl_char(),false,false,true);
        break;
    }
    default :
    {
        return 1;
        break;
    }
    }
    oth.finalize();
    std::cout << std::endl;
    return 0;
}



int gc_on_stream(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.clear();
    gc_on_stream(std::cin,std::cout);
    return 1;
}



int extract_graph(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    log.step("read");
    std::map<int,std::vector<int>> tile_ids;;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
            std::vector<int>  & format_labs = w_datas_tri[hid].format_labs ;
            if(format_labs.size() == 0)
            {
                int nbs = w_datas_tri[hid].format_dst.size();
                for(int ss = 0; ss < nbs ; ss++)
                {
                    format_labs.push_back(0);
                }
            }
        }
        if(hpi.get_lab() == "s")
        {
            std::vector<int> vv(3);
            for(int d = 0; d < 3; d++)
            {
                hpi.get_input_stream() >> vv[d];
            }
            tile_ids[hid] = vv;
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    log.step("compute");
    tbmrf_reco<DTW,D_MAP> mrf(params.nb_labs,&tri,&w_datas_tri);
    mrf.lambda = params.lambda;
    mrf.set_mode(params.mode);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid),osh("s","s",tid);;
    int nbc = 0;
    nbc = mrf.extract_stream_graph_v2(1,tri,w_datas_tri,tile_ids,oth.get_output_stream(),tid,params.graph_type,params.area_processed,params.coef_mult);
    oth.finalize();
    std::cout << std::endl;
    if(params.graph_type == 0)
    {
        osh.write_header(std::cout);
        osh.get_output_stream() << 0 << " " ;
        osh.get_output_stream() << 0 << " " ;
        osh.get_output_stream() << nbc << " " ;
        osh.finalize();
        std::cout << std::endl;
    }
    return 0;
}






int fill_graph(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    int D = Traits::D;
    std::map<int,std::vector<int>> tile_ids;
    std::map<int,int> labs_map;
    bool tri_parsed = false;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
            w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
            tri_parsed = true;
        }
        if(hpi.get_lab() == "s")
        {
            std::vector<int> vv(3);
            for(int d = 0; d < 3; d++)
            {
                hpi.get_input_stream() >> vv[d];
            }
            tile_ids[hid] = vv;
        }
        if(hpi.get_lab() == "l")
        {
            int id,lab_val,nb_elems;
            hpi.get_input_stream() >> nb_elems;
            for(int ss = 0; ss < nb_elems ; ss++)
            {
                hpi.get_input_stream() >> id;
                hpi.get_input_stream() >> lab_val;
                labs_map[id] = lab_val;
            }
        }
        hpi.finalize();
    }
    tri.finalize(sch);
    std::map<int,int> g2lmap;
    for( auto cit = tri.cells_begin();
            cit != tri.cells_end(); ++cit )
    {
        int lid = cit->lid();
        int gid = cit->gid();
        g2lmap[gid] = lid;
    }
    int nbs = tri.number_of_cells();
    std::vector<int>  & format_labs = w_datas_tri[tid].format_labs ;
    format_labs.resize(nbs);
    for ( auto it = labs_map.begin(); it != labs_map.end(); it++ )
    {
        format_labs[g2lmap[it->first]] = it->second;
    }
    log.step("compute");
    w_datas_tri[tid].fill_labs(w_datas_tri[tid].format_labs);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    return 0;
}



int seg_lagrange(Id tid_1,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    double tau = params.tau;
    bool is_first = true;
    log.step("read");
    std::map<Id,std::map<Id,SharedData> > shared_data_map;
    std::map<Id,std::map<Id,SharedData> > edges_data_map;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
            std::vector<int>  & format_labs = w_datas_tri[hid].format_labs ;
            if(w_datas_tri[hid].dmap[w_datas_tri[hid].labseg_name].do_exist)
            {
                is_first = false;
                w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
            }
            else
            {
                int nbs = w_datas_tri[hid].format_dst.size();
                for(int ss = 0; ss < nbs ; ss++)
                {
                    format_labs.push_back(0);
                }
            }
        }
        if(hpi.get_lab() == "e")
        {
            Id eid2 = hpi.get_id(1);
            shared_data_map[eid2] = std::map<Id,SharedData>();
            std::map<Id,SharedData>  & lp = shared_data_map[eid2];
            read_id_double_serialized(lp, hpi.get_input_stream(),tid_1 == 3);
        }
        if(hpi.get_lab() == "f")
        {
            Id eid1 = hpi.get_id(0);
            edges_data_map[eid1] = std::map<Id,SharedData>();
            std::map<Id,SharedData>  & lp = edges_data_map[eid1];
            read_id_double_serialized(lp, hpi.get_input_stream(),tid_1 == 3);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    Tile_const_iterator tile_1 = tri.get_const_tile(tid_1);
    int D = tile_1->current_dimension();
    Tile_iterator main_tile = tri.get_tile(tid_1);
    tbmrf_reco<DTW,D_MAP> mrf(params.nb_labs,&tri,&w_datas_tri);
    mrf.lambda = params.lambda;
    mrf.set_mode(params.mode);
    int acc_tot = 0;
    int acc_diff = 0;
    for( auto cit_1 = tile_1->cells_begin();
            cit_1 != tile_1->cells_end(); ++cit_1 )
    {
        Cell_const_iterator fch = Cell_const_iterator(tile_1,tile_1, tile_1, cit_1);
        if(!tile_1->cell_is_mixed(cit_1) || tile_1->cell_is_infinite(cit_1))
            continue;
        Id lid_1 = tile_1->lid(cit_1);
        int lab_1 = 0;
        if(!is_first)
        {
            lab_1 = w_datas_tri[tid_1].format_labs[lid_1];
        }
        int D = tile_1->current_dimension();
        std::unordered_set<Id> idSet ;
        for(int i=0; i<=D; ++i)
        {
            Id tid_2 = tile_1->id(tile_1->vertex(cit_1,i));
            if (idSet.find(tid_2) != idSet.end())
            {
                continue;
            }
            idSet.insert(tid_2);
            if(tid_2 == tid_1)
                continue;
            double new_lagrange = 0;
            double cur_tau = mrf.get_volume_reco(fch) * params.lambda * params.coef_mult;
            Id lid_2 = std::get<0>(shared_data_map[tid_2][lid_1]);
            double v_diff = 0;
            if(!is_first)
            {
                new_lagrange = std::get<1>(shared_data_map[tid_2][lid_1]);
                cur_tau = std::get<2>(shared_data_map[tid_2][lid_1]);
                double old_lagrange = new_lagrange;
                double old_diff = std::get<3>(shared_data_map[tid_2][lid_1]);
                if(edges_data_map[tid_2].find(lid_2) == edges_data_map[tid_2].end())
                    std::cerr << "!! Warning !! edges_data_map[tid_2].find(lid_2) == edges_data_map[tid_2].end()" << std::endl;
                int lab_2 = std::get<1>(edges_data_map[tid_2][lid_2]);
                if(lab_2 != lab_1)
                    acc_diff++;
                acc_tot++;
                if(tid_1 < tid_2)
                    v_diff = lab_1-lab_2;
                else
                    v_diff = lab_2-lab_1;
                new_lagrange  = new_lagrange + cur_tau*(v_diff);
                if(old_diff != v_diff)
                    cur_tau /=2;
                if(params.do_finalize)
                {
                    Id mid = fch->main_id();
                    if(tid_2 == mid)
                        w_datas_tri[tid_1].format_labs[lid_1] = lab_2;
                }
            }
            else
            {
                if(shared_data_map.find(tid_2) == shared_data_map.end())
                    shared_data_map[tid_2] = std::map<Id,SharedData>();
            }
            shared_data_map[tid_2][lid_1] = std::make_tuple(lid_2,new_lagrange,cur_tau,v_diff);
        }
    }
    log.step("compute");
    if(!params.do_finalize)
        mrf.opt_gc_lagrange(1,tri,w_datas_tri,shared_data_map,tid_1,params.use_weight);
    log.step("write");
    std::cout.clear();
    // ==== Build new edges =====
    if(!params.do_finalize)
    {
        edges_data_map.clear();
        for( auto cit_1 = tile_1->cells_begin();
                cit_1 != tile_1->cells_end(); ++cit_1 )
        {
            Cell_const_iterator fch = Cell_const_iterator(tile_1,tile_1, tile_1, cit_1);
            // lagrangian only on mixed cell
            if(!tile_1->cell_is_mixed(cit_1) || tile_1->cell_is_infinite(cit_1))
                continue;
            Id lid_1 = tile_1->lid(cit_1);
            int lab_1 = w_datas_tri[tid_1].format_labs[lid_1];
            int D = tile_1->current_dimension();
            std::unordered_set<Id> idSet ;
            for(int i=0; i<=D; ++i)
            {
                // Select only id != tid_1
                Id tid_2 = tile_1->id(tile_1->vertex(cit_1,i));
                Id lid_2 = std::get<0>(shared_data_map[tid_2][lid_1]);
                if (idSet.find(tid_2) != idSet.end())
                {
                    continue;
                }
                idSet.insert(tid_2);
                if(tid_2 == tid_1)
                    continue;
                if(edges_data_map.find(tid_2) == edges_data_map.end())
                    edges_data_map[tid_2] = std::map<Id,SharedData>();
                //  send to the neighbor tile the value of the tets at LID_1
                edges_data_map[tid_2][lid_1] = std::make_tuple(lid_1,((double)lab_1),0,0);
            }
        }
    }
    else
    {
        for( auto cit_1 = tile_1->cells_begin();
                cit_1 != tile_1->cells_end(); ++cit_1 )
        {
            Cell_const_iterator fch = Cell_const_iterator(tile_1,tile_1, tile_1, cit_1);
            Id mid = fch->main_id();
            // lagrangian only on mixed cell
            if(!tile_1->cell_is_mixed(cit_1) || tile_1->cell_is_infinite(cit_1))
                continue;
            Id lid_1 = tile_1->lid(cit_1);
            int lab_1 = w_datas_tri[tid_1].format_labs[lid_1];
            int D = tile_1->current_dimension();
            std::unordered_set<Id> idSet ;
            for(int i=0; i<=D; ++i)
            {
                // Select only id != tid_1
                Id tid_2 = tile_1->id(tile_1->vertex(cit_1,i));
                Id lid_2 = std::get<0>(shared_data_map[tid_2][lid_1]);
                int lab_2 = std::get<1>(edges_data_map[tid_2][lid_2]);
                if(tid_2 == tid_1)
                    continue;
                if (idSet.find(tid_2) != idSet.end())
                {
                    continue;
                }
                idSet.insert(tid_2);
                if(tid_2 == mid)
                    w_datas_tri[tid_1].format_labs[lid_1] = lab_2;
            }
        }
    }
    log.step("finalize");
    w_datas_tri[tid_1].fill_labs(w_datas_tri[tid_1].format_labs);
    // ============== Dump results =================
    ddt::stream_data_header oth("t","z",tid_1);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid_1], oth.get_output_stream(),tid_1,false,log);
    oth.finalize();
    std::cout << std::endl;
    bool do_dump_u = false;
    if(params.do_finalize && do_dump_u)
    {
        D_MAP w_datas_tri_fz;
        for(auto ee : edges_data_map)
        {
            Id tid2 = ee.first;
            if(tid2 > tid_1)
            {
                ddt::stream_data_header oth("u","z",std::vector<int> {tid_1,tid2});
                std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1));
                oth.write_header(std::cout);
                ddt::write_ddt_stream(tri, w_datas_tri[tid_1], oth.get_output_stream(),tid_1,false,log);
                oth.finalize();
                std::cout << std::endl;
            }
        }
    }
    for(auto ee : shared_data_map)
    {
        Id tid2 = ee.first;
        if(tid2 == tid_1)
            continue;
        ddt::stream_data_header hto("e","z",std::vector<int> {tid_1,tid2});
        std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1) + "_nid" + std::to_string(tid2));
        hto.write_header(std::cout);
        write_id_double_serialized(ee.second,hto.get_output_stream());
        hto.finalize();
        std::cout << std::endl;
    }
    for(auto ee : edges_data_map)
    {
        Id tid2 = ee.first;
        ddt::stream_data_header hto("f","z",std::vector<int> {tid_1,tid2});
        std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid_1) + "_nid" + std::to_string(tid2));
        hto.write_header(std::cout);
        write_id_double_serialized(ee.second,hto.get_output_stream());
        hto.finalize();
        std::cout << std::endl;
    }
    ddt::stream_data_header sth("s","z",tid_1);
    sth.write_header(std::cout);
    std::cout << "[diff:tot] " << acc_diff << " " << acc_tot ;
    sth.finalize();
    std::cout << std::endl;
    return 0;
}



int seg(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
            std::vector<int>  & format_labs = w_datas_tri[hid].format_labs ;
            if(format_labs.size() == 0)
            {
                int nbs = w_datas_tri[hid].format_dst.size();
                for(int ss = 0; ss < nbs ; ss++)
                {
                    format_labs.push_back(0);
                }
            }
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    // ===== Init the id of each cell
    log.step("compute");
    tbmrf_reco<DTW,D_MAP> mrf(params.nb_labs,&tri,&w_datas_tri);
    mrf.lambda = params.lambda;
    mrf.set_mode(params.mode);
    mrf.opt_gc(1,tri,w_datas_tri);
    log.step("finalize");
    w_datas_tri[tid].fill_labs(w_datas_tri[tid].format_labs);
    log.step("write");
    std::cout.clear();
    ddt::stream_data_header oth("t","z",tid);
    std::string filename(params.output_dir + "/" + params.slabel + "_id" + std::to_string(tid));
    oth.write_header(std::cout);
    ddt::write_ddt_stream(tri, w_datas_tri[tid], oth.get_output_stream(),tid,false,log);
    oth.finalize();
    std::cout << std::endl;
    return 0;
}




int seg_global_extract(Id tid,wasure_params & params,int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    D_MAP w_datas_tri_ori;
    D_MAP w_datas_tri;
    std::map<Id,std::vector<int>> v_map;
    log.step("read");
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_dst(w_datas_tri[hid].format_dst,false);
            std::vector<int>  & format_labs = w_datas_tri[hid].format_labs ;
            if(w_datas_tri[hid].dmap[w_datas_tri[hid].labseg_name].do_exist)
            {
                w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
            }
            else
            {
                int nbs = w_datas_tri[hid].format_dst.size();
                for(int ss = 0; ss < nbs ; ss++)
                {
                    format_labs.push_back(0);
                }
            }
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    tbmrf_reco<DTW,D_MAP> mrf_ori(params.nb_labs,&tri,&w_datas_tri);
    mrf_ori.lambda = params.lambda;
    mrf_ori.set_mode(params.mode);
    double energy_largrange = mrf_ori.get_energy(tri,w_datas_tri);
    for (auto it = w_datas_tri.begin(); it != w_datas_tri.end(); it++)
    {
        int hid = it->first;
        std::vector<int>  & format_labs = w_datas_tri[hid].format_labs ;
        std::vector<int> format_labs_new(format_labs);
        v_map[hid] = format_labs_new;
        std::fill(format_labs.begin(), format_labs.end(), 0);
    }
    log.step("compute");
    tbmrf_reco<DTW,D_MAP> mrf(params.nb_labs,&tri,&w_datas_tri);
    mrf.lambda = params.lambda;
    mrf.set_mode(params.mode);
    mrf.opt_gc(1,tri,w_datas_tri);
    double energy_full = mrf.get_energy(tri,w_datas_tri);
    int acc_tot = 0;
    int acc_diff = 0;
    double wacc_tot = 0;
    double wacc_diff = 0;
    for( auto cit = tri.cells_begin();
            cit != tri.cells_end(); ++cit )
    {
        if(cit->is_foreign())
        {
            continue;
        }
        int cccid = cit->lid();
        Cell_const_iterator fch = *cit;
        double cur_tau =  mrf.get_volume_reco(fch);
        int lcurr = w_datas_tri[fch->tile()->id()].format_labs[cccid];
        int lnew = v_map[fch->tile()->id()][cccid];
        if(lcurr != lnew)
        {
            acc_diff++;
            wacc_diff+=cur_tau;
        }
        acc_tot++;
        wacc_tot+=cur_tau;
    }
    // =============
    // SURFACE EXTRACTION
    std::vector<Facet_const_iterator> lft;
    std::vector<bool> lbool;
    std::cout.setstate(std::ios_base::failbit);
    int D = Traits::D;
    int mode = params.mode;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        try
        {
            if(fit->is_infinite())
                continue;
            Cell_const_iterator tmp_fch = fit.full_cell();
            int tmp_idx = fit.index_of_covertex();
            Cell_const_iterator tmp_fchn = tmp_fch->neighbor(tmp_idx);
            if(!tri.tile_is_loaded(tmp_fch->main_id()) ||
                    !tri.tile_is_loaded(tmp_fchn->main_id()))
            {
                std::cerr << "ERROR tile not loaded" << std::endl;
                continue;
            }
            bool is_on_convex = false;
            if(tmp_fch->is_infinite() ||  tmp_fchn->is_infinite() )
                is_on_convex = true;
            Cell_const_iterator fch = tmp_fch->main();
            int id_cov = fit.index_of_covertex();
            Cell_const_iterator fchn = tmp_fchn->main();
            int cccid = fch->lid();
            int cccidn = fchn->lid();
            int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
            int chnlab = w_datas_tri[fchn->tile()->id()].format_labs[cccidn];
            if(
                (ch1lab != chnlab)  || (mode == 0 && (is_on_convex && (ch1lab == 0 || chnlab == 0)))
            )
            {
                lft.push_back(*fit);
                const Point& a = fch->vertex((id_cov+1)&3)->point();
                const Point& b = fch->vertex((id_cov+2)&3)->point();
                const Point& c = fch->vertex((id_cov+3)&3)->point();
                const Point& d = fch->vertex((id_cov)&3)->point();
                bool bl =
                    (CGAL::orientation(a,b,c,d) == 1 && chnlab == 0) ||
                    (CGAL::orientation(a,b,c,d) == -1 && chnlab == 1);
                lbool.push_back(!bl);
            }
        }
        catch (ddt::DDT_exeption& e)
        {
            std::cerr << "!! WARNING !!!" << std::endl;
            std::cerr << "Exception catched : " << e.what() << std::endl;
            continue;
        }
    }
    std::string ply_name(params.output_dir +  "/" + params.slabel + "_id_surface");
    std::cout.clear();
    ddt::stream_data_header sth("s","z",tid);
    sth.write_header(std::cout);
    std::cout << "[diff:tot:e_ori:e_new] " << wacc_diff << " " << wacc_tot << " " << energy_largrange << " " << energy_full;
    sth.finalize();
    std::cout << std::endl;
    ddt::stream_data_header oth("p","z",0);
    if(D == 2)
    {
        oth.write_into_file(ply_name,".geojson");
        oth.write_header(std::cout);
    }
    switch (D)
    {
    case 2 :
    {
        wasure::dump_2d_surface_geojson<DTW>(lft,oth.get_output_stream());
        break;
    }
    case 3 :
    {
        std::vector<Point>  format_points;
        std::vector<int> v_simplex;
        std::map<Vertex_const_iterator, uint> vertex_map;
        ddt_data<Traits> datas_out;
        int acc = 0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            for(int i = 0; i < D+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    if(vertex_map.find(v) == vertex_map.end())
                    {
                        vertex_map[v] = acc++;
                        format_points.push_back(v->point());
                    }
                }
            }
        }
        acc=0;
        for(auto fit = lft.begin(); fit != lft.end(); ++fit)
        {
            Cell_const_iterator fch = fit->full_cell();
            int id_cov = fit->index_of_covertex();
            Id ida = (id_cov+1)&3;
            Id idb = (id_cov+2)&3;
            Id idc = (id_cov+3)&3;
            v_simplex.push_back(vertex_map[fch->vertex(ida)]);
            if(!lbool[acc])
            {
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
            }
            else
            {
                v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                v_simplex.push_back(vertex_map[fch->vertex(idb)]);
            }
            acc++;
        }
        datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",D,D,tinyply::Type::INT32);
        datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
        datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
        datas_out.write_ply_stream(oth.get_output_stream(),oth.get_nl_char());
        break;
    }
    default :             // Note the colon, not a semicolon
    {
        return 1;
        break;
    }
    }
    oth.finalize();
    std::cout << std::endl;
    return 0;
}




int tri2geojson(Id tid,wasure_params & params, int nb_dat,ddt::logging_stream & log)
{
    std::cout.setstate(std::ios_base::failbit);
    DTW tri;
    Scheduler sch(1);
    wasure_algo w_algo;
    int D = Traits::D;
    D_MAP w_datas_tri;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        if(hpi.get_lab() == "t")
        {
            w_datas_tri[hid] = wasure_data<Traits>();
            bool do_clean_data = false;
            bool do_serialize = false;
            ddt::read_ddt_stream(tri,w_datas_tri[hid], hpi.get_input_stream(),hid,do_serialize,do_clean_data,log);
            w_datas_tri[hid].extract_labs(w_datas_tri[hid].format_labs,false);
        }
        tri.finalize(sch);
        hpi.finalize();
    }
    std::string json_name(params.output_dir +  "/" + params.slabel + "_" + std::to_string(tid) + "_tri");
    std::cout.clear();
    ddt::stream_data_header oth("j","h",tid);
    oth.write_into_file(json_name,".geojson");
    oth.write_header(std::cout);
    wasure::write_geojson_tri_wasure(tri,w_datas_tri,oth.get_output_stream());
    oth.finalize();
    std::cout << std::endl;
    return 0;
}





int ply2geojson(Id tile_id,wasure_params & params,int nb_dat)
{
    std::cout.setstate(std::ios_base::failbit);
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        wasure_data<Traits> w_datas;
        hpi.parse_header(std::cin);
        if(hpi.get_lab() == "z")
        {
            w_datas.read_serialized_stream(hpi.get_input_stream());
        }
        hpi.finalize();
        std::string json_name(params.output_dir +  "/" + params.slabel + "_" + std::to_string(tile_id) + "_geo");
        std::cout.clear();
        Id id = hpi.get_id(0);
        ddt::stream_data_header oqh_1("p","s",id),oqh_2("p","s",id),oqh_3("p","s",id);
        std::string filename(params.output_dir + "/" + params.slabel +"_id_"+ std::to_string(tile_id) + "_" + std::to_string(id));
        oqh_1.write_into_file(filename,"_pts.geojson");
        oqh_1.write_header(std::cout);
        oqh_2.write_into_file(filename,"_spx.geojson");
        oqh_2.write_header(std::cout);
        oqh_3.write_into_file(filename,"_nrm.geojson");
        oqh_3.write_header(std::cout);
        w_datas.write_geojson_norms(oqh_3.get_output_stream());
        oqh_1.finalize();
        oqh_2.finalize();
        oqh_3.finalize();
        ddt::add_qgis_style(oqh_2.get_file_name(), std::string("cell_style_flag.qml"));
        std::cout << std::endl;
    }
    return 0;
}


int hello(Id tile_id,wasure_params & params,int nb_dat)
{
    std::cout.setstate(std::ios_base::failbit);
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        wasure_data<Traits> w_datas;
        hpi.parse_header(std::cin);
        hpi.finalize();
        std::cout.clear();
        ddt::stream_data_header oth("j","f",tile_id);
        oth.write_into_file("hello",".json");
        oth.write_header(std::cout);
        oth.finalize();
        std::cout << std::endl;
    }
    return 0;
}


using Kernel = CGAL::Simple_cartesian<double>;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Point_with_info = std::tuple<Point_3, Vector_3, float,unsigned short>;
using Point_with_info_2 = std::tuple<Point_3, Vector_3, Vector_3>;
using Point_map = CGAL::Nth_of_tuple_property_map<0, Point_with_info>;
using Normal_map = CGAL::Nth_of_tuple_property_map<1, Point_with_info>;
using Scan_angle_map = CGAL::Nth_of_tuple_property_map<2, Point_with_info>;
using Scanline_id_map = CGAL::Nth_of_tuple_property_map<3, Point_with_info>;
using Point_map_2 = CGAL::Nth_of_tuple_property_map<0, Point_with_info_2>;
using Normal_map_2 = CGAL::Nth_of_tuple_property_map<1, Point_with_info_2>;
using Ori_map = CGAL::Nth_of_tuple_property_map<2, Point_with_info_2>;
using Conf_map = CGAL::Nth_of_tuple_property_map<3, Point_with_info_2>;


void dump_2(const char* filename, const std::vector<Point_with_info_2>& points)
{
    std::ofstream ofile (filename, std::ios::binary);
    CGAL::IO::set_binary_mode(ofile);
    CGAL::IO::write_PLY_with_properties
    (ofile, points,
     CGAL::make_ply_point_writer (Point_map_2()),
     std::make_tuple(Normal_map_2(),
                     CGAL::IO::PLY_property<double>("nx"),
                     CGAL::IO::PLY_property<double>("ny"),
                     CGAL::IO::PLY_property<double>("nz")),
     std::make_tuple(Ori_map(),
                     CGAL::IO::PLY_property<double>("x_origin"),
                     CGAL::IO::PLY_property<double>("y_origin"),
                     CGAL::IO::PLY_property<double>("z_origin")));
    return;
}


std::string extractDirectoryPath(const std::string& filePath)
{
    size_t found = filePath.find_last_of("/\\");
    if (found != std::string::npos)
    {
        return filePath.substr(0, found);
    }
    return "";
}


int compute_bbox(Id tid,wasure_params & params, int nb_dat)
{
    std::cout.setstate(std::ios_base::failbit);
    Traits  traits;
    wasure_algo w_algo;
    int D = Traits::D;
    std::map<Id,wasure_data<Traits> > datas_map;
    std::map<Id,std::string > fname_map;
    std::map<Id,std::string > fname_map2;
    ddt::Bbox<Traits::D> full_bbox;
    int full_nbp = 0;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        std::string fname = hpi.get_file_name();
        auto & ifile = hpi.get_input_stream();
        std::vector<Point_with_info> points;
        std::vector<Point_with_info_2> points_2;
        if(hpi.get_lab() == "p")
        {
            if (!ifile ||
                    !CGAL::IO::read_LAS_with_properties
                    (ifile, std::back_inserter (points),
                     CGAL::IO::make_las_point_reader (Point_map()),
                     std::make_pair (Scan_angle_map(),
                                     CGAL::IO::LAS_property::Scan_angle()),
                     std::make_pair (Scanline_id_map(),
                                     CGAL::IO::LAS_property::Point_source_ID())))
            {
                std::cerr << "Can't read " << fname << std::endl;
                return EXIT_FAILURE;
            }
            hpi.finalize();
        }
        std::cout.clear();
        CGAL::Bbox_3 bbox;
        double alpha = 200.0;
        for (auto pp : points)
        {
            auto vx = std::get<0>(pp)[0];
            auto vy = std::get<0>(pp)[1];
            auto vz = std::get<0>(pp)[2];
            bbox = bbox +  CGAL::Bbox_3(vx,vy,vz,
                                        vx,vy,vz);
        }
        ddt::stream_data_header sth("s","z",0);
        sth.write_header(std::cout);
        std::cout <<  bbox.xmin() << " " <<  bbox.xmax() <<  " " << bbox.ymin() << " " << bbox.ymax()<<  " " << bbox.zmin() << " " << bbox.zmax() << " " << points.size();
        sth.finalize();
        std::cout << std::endl;
    }
    return 0;
}


int preprocess(Id tid,wasure_params & params, int nb_dat)
{
    std::cout.setstate(std::ios_base::failbit);
    Traits  traits;
    wasure_algo w_algo;
    int D = Traits::D;
    std::map<Id,wasure_data<Traits> > datas_map;
    std::map<Id,std::string > fname_map;
    std::map<Id,std::string > fname_map2;
    ddt::Bbox<Traits::D> bbox_ori;
    std::stringstream ss;
    ss << params.bbox_string;
    ss >> bbox_ori;
    ddt::Bbox<Traits::D> full_bbox;
    int full_nbp = 0;
    for(int i = 0; i < nb_dat; i++)
    {
        ddt::stream_data_header hpi;
        hpi.parse_header(std::cin);
        Id hid = hpi.get_id(0);
        std::string fname = hpi.get_file_name();
        auto & ifile = hpi.get_input_stream();
        std::vector<Point_with_info> points;
        std::vector<Point_with_info_2> points_2;
        if(hpi.get_lab() == "p")
        {
            if (!ifile ||
                    !CGAL::IO::read_LAS_with_properties
                    (ifile, std::back_inserter (points),
                     CGAL::IO::make_las_point_reader (Point_map()),
                     std::make_pair (Scan_angle_map(),
                                     CGAL::IO::LAS_property::Scan_angle()),
                     std::make_pair (Scanline_id_map(),
                                     CGAL::IO::LAS_property::Point_source_ID())))
            {
                std::cerr << "Can't read " << fname << std::endl;
                return EXIT_FAILURE;
            }
            hpi.finalize();
        }
        std::cout.clear();
        // CGAL::jet_estimate_normals<CGAL::Parallel_if_available_tag>
        // (points, 100,
        //  CGAL::parameters::point_map (Point_map()).
        //  normal_map (Normal_map()));
        CGAL::scanline_orient_normals
        (points,
         CGAL::parameters::point_map (Point_map()).
         normal_map (Normal_map()).
         scan_angle_map (Scan_angle_map()).
         scanline_id_map (Scanline_id_map()));
        CGAL::Bbox_3 bbox;
        double alpha = 200.0;
        int acc=0;
        int max_ppf = 2000000;
        for (auto pp : points)
        {
            auto vx = std::get<0>(pp)[0];
            auto vy = std::get<0>(pp)[1];
            auto vz = std::get<0>(pp)[2];
            auto lx = std::get<1>(pp)[0];
            auto ly = std::get<1>(pp)[1];
            auto lz = std::get<1>(pp)[2];
            auto ID = std::get<2>(pp);
            Kernel::Vector_3 ori(
                std::get<0>(pp)[0] + alpha*lx - bbox_ori.min(0),
                std::get<0>(pp)[1] + alpha*ly - bbox_ori.min(1),
                std::get<0>(pp)[2] + alpha*lz - bbox_ori.min(2)
            );
            Kernel::Point_3 pp_new(
                std::get<0>(pp)[0] - bbox_ori.min(0),
                std::get<0>(pp)[1] - bbox_ori.min(1),
                std::get<0>(pp)[2] - bbox_ori.min(2)
            );
            bbox = bbox +  CGAL::Bbox_3(pp_new[0],pp_new[1],pp_new[2],
                                        pp_new[0],pp_new[1],pp_new[2]);
            double conf = 0.9;
            points_2.push_back(std::make_tuple(pp_new,std::get<1>(pp),ori));
            acc++;
            if(points_2.size() > max_ppf or acc == points.size() -1)
            {
                std::string oname(params.output_dir + "/" + params.slabel +"_id_"+ std::to_string(hid) + "_" + std::to_string(acc) + ".ply");
                dump_2(oname.c_str(), points_2);
                ddt::stream_data_header sth("s","z",0);
                sth.write_header(std::cout);
                std::cout <<  bbox.xmin() << " " <<  bbox.xmax() <<  " " << bbox.ymin() << " " << bbox.ymax()<<  " " << bbox.zmin() << " " << bbox.zmax() << " " << points.size();
                sth.finalize();
                std::cout << std::endl;
                bbox = CGAL::Bbox_3();
                points_2.clear();
            }
        }
    }
    return 0;
}



int main(int argc, char **argv)
{
    std::cout.setstate(std::ios_base::failbit);
    wasure_params params;
    params.parse(argc,argv);
    int rv = 0;
    int acc = 0;
    bool do_dump_log = true;
    while(true)
    {
        ddt::stream_app_header sah;
        Id tile_id = -1;
        int nb_dat = -1;
        if(!params.skip_app_header)
        {
            sah.parse_header(std::cin);
            if(sah.is_void())
                return 0;
            tile_id = ((Id)sah.tile_id);
            nb_dat = sah.get_nb_dat();
        }
        srand(time(NULL));
        ddt::logging_stream log(std::to_string(tile_id) + "_" + params.algo_step, params.log_level);
        if(acc == 0)
        {
            std::cerr << " ======================================================= " << std::endl;
            std::cerr << "     [MAIN_DDT_STREAM_LOG] stream  : " << params.slabel << std::endl;
            std::cerr << "     [MAIN_DDT_STREAM_LOG] tile_id : " << tile_id <<  std::endl;
            std::cerr << "     [MAIN_DDT_STREAM_LOG] step    : " << params.algo_step << std::endl;
            std::cerr << "     [MAIN_DDT_STREAM_LOG] acc     : " << acc++ << std::endl;
            std::cerr << "     [MAIN_DDT_STREAM_LOG] nb dat  : " << nb_dat << std::endl;
            std::cerr << " ======================================================= " << std::endl;
        }
        try
        {
            if(params.algo_step == std::string("hello"))
            {
                rv = hello(tile_id,params,nb_dat);
            }
            else if(params.algo_step == std::string("dim"))
            {
                rv = dim_splitted(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("simplify"))
            {
                rv = simplify(tile_id,params,nb_dat);
            }
            else if(params.algo_step == std::string("regularize_slave_focal"))
            {
                rv = regularize_slave_focal(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("preprocess"))
            {
                rv = preprocess(tile_id,params,nb_dat);
            }
            else if(params.algo_step == std::string("compute_bbox"))
            {
                rv = compute_bbox(tile_id,params,nb_dat);
            }
            else if(params.algo_step == std::string("regularize_slave_extract"))
            {
                rv = regularize_slave_extract(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("regularize_slave_insert"))
            {
                rv = regularize_slave_insert(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("seg"))
            {
                rv = seg(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("seg_global_extract"))
            {
                rv = seg_global_extract(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("seg_lagrange_raw"))
            {
                params.use_weight = false;
                rv = seg_lagrange(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("seg_lagrange_weight"))
            {
                params.use_weight = true;
                rv = seg_lagrange(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("extract_surface"))
            {
                if(params.area_processed == 0)
                    rv = extract_surface_new(tile_id,params,nb_dat,log);
                else
                    rv = extract_surface_area(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("extract_graph"))
            {
                rv = extract_graph(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("fill_graph"))
            {
                rv = fill_graph(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("dst"))
            {
                rv = dst(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("gc_on_stream"))
            {
                rv = gc_on_stream(tile_id,params,nb_dat,log);
                return 0;
            }
            else if(params.algo_step == std::string("tri2geojson"))
            {
                rv = tri2geojson(tile_id,params,nb_dat,log);
            }
            else if(params.algo_step == std::string("ply2geojson"))
            {
                if(Traits::D == 2)
                {
                    rv = ply2geojson(tile_id,params,nb_dat);
                }
                else
                {
                    rv = 0;
                }
            }
            else
            {
                std::cerr << "==== exe : no params ===" << std::endl;
                return 1;
            }
            if(rv != 0) return rv;
            if(do_dump_log)
            {
                log.dump_log(std::cerr);
            }
        }
        catch (std::exception& e)
        {
            std::cerr << "Exception catched : " << e.what() << std::endl;
            std::cerr << "tile_id               : " << tile_id << std::endl;
        }
    }
    std::cerr << "[MAIN_DDT_STREAM_LOG] end exe " << std::endl;
    return rv;
}






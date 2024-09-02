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
#ifndef WASURE_DATA
#define WASURE_DATA

#include <map>
#include "tinyply.h"
#include "ddt_data.hpp"


using namespace tinyply;

typedef typename ddt_data<Traits>::Data_ply Data_ply;
template<typename Traits>
class wasure_data : public ddt_data<Traits>
{
public :

    typedef ddt_data<Traits> ddtd;

    wasure_data() : ddt_data<Traits>(), tile_ids(3)
    {
        init_wasure_name();
        init_wasure_map();
    }

    wasure_data(std::map<std::vector<std::string>, Data_ply > & init_dmap)
    {
        init_wasure_name();
        for ( const auto &ee : init_dmap )
        {
            ddt_data<Traits>::dmap[ee.first] =  Data_ply(ee.first,ee.second.part,Traits::D,ee.second.get_vsize(),ee.second.type);
        }
    }



    void init_wasure_map()
    {
        int D = Traits::D;
        ddt_data<Traits>::dmap[center_name] = Data_ply(center_name,"vertex",D,D,DATA_FLOAT_TYPE);
        ddt_data<Traits>::dmap[normal_name] = Data_ply(normal_name,"vertex",D,D,DATA_FLOAT_TYPE);
        ddt_data<Traits>::dmap[egv_name] = Data_ply(egv_name,"vertex",D*D,D*D,DATA_FLOAT_TYPE);
        ddt_data<Traits>::dmap[flags_name] = Data_ply(flags_name,"vertex",1,1,tinyply::Type::INT32);
        ddt_data<Traits>::dmap[dst_name] = Data_ply(dst_name,"face",3,3,DATA_FLOAT_TYPE);
        ddt_data<Traits>::dmap[labseg_name] = Data_ply(labseg_name,"face",1,1,tinyply::Type::INT32);
        ddt_data<Traits>::dmap[gid_name] = Data_ply(gid_name,"face",1,1,tinyply::Type::INT32);
        ddt_data<Traits>::dmap[sig_name] = Data_ply(sig_name,"vertex",D,D,DATA_FLOAT_TYPE);
    }


    void init_wasure_name()
    {
        int D = Traits::D;
        sig_name = ddt_data<Traits>::subvect({"sigma1","sigma2","sigma3","sigma4"},D);
        dst_name = {"vpe","vpo","vpu"};
        labseg_name = {"seg"};
        gid_name = {"gid"};
        glob_scale_name = {"glob_scale"};
        flags_name = {"flags"};
        center_name = ddt_data<Traits>::subvect({"x_origin","y_origin","z_origin","t_origin"},D);
        normal_name = ddt_data<Traits>::subvect({"nx","ny","nz","nt"},D);
        switch(D)
        {
        case 1 :
            egv_name = {"eigenVector1x"}; // prints "1"
            break;       // and exits the switch
        case 2 :
            egv_name = {"eigenVector1x","eigenVector1y","eigenVector2x","eigenVector2y"};
            break;
        case 3 :
            egv_name = {"eigenVector1x","eigenVector1y","eigenVector1z",
                        "eigenVector2x","eigenVector2y","eigenVector2z",
                        "eigenVector3x","eigenVector3y","eigenVector3z"
                       };
            break;
        }
    }


    void write_geojson_norms(std::ostream & ofs, bool is_first = true)
    {
        std::vector<Point> raw_points;
        std::vector<Point> raw_centers;
        std::vector<double> raw_egv;
        std::vector<double> raw_sigs;
        ddt_data<Traits>::dmap[ddt_data<Traits>::xyz_name].extract_full_uint8_vect(raw_points,false);
        ddt_data<Traits>::dmap[center_name].extract_full_uint8_vect(raw_centers,false);
        ddt_data<Traits>::dmap[sig_name].extract_raw_uint8_vect(raw_sigs,false);
        ddt_data<Traits>::dmap[egv_name].extract_raw_uint8_vect(raw_egv,false);
        bool do_ray = false;
        int D = Traits::D;
        std::vector<std::string> lab_color = {"\"red\"","\"green\"","\"blue\""};
        ofs << "{" << std::endl;
        ofs << "\"type\": \"FeatureCollection\"," << std::endl;
        ofs << "\"features\": [" << std::endl;
        for(int id = 0; id < ddtd::nb_pts_uint8_vect(); id++)
        {
            int id_pts = id*D;
            int id_sigs = id*D;
            int id_egv = id*D*D;
            if(!is_first)
                ofs << "," << std::endl;
            is_first=false;
            // Center
            if(raw_centers.size() > 0)
            {
                ofs << "{" << std::endl;
                ofs << "\"type\": \"Feature\"," << std::endl;
                ofs << "\"geometry\": {" << std::endl;
                ofs << "\"type\": \"Point\"," << std::endl;
                ofs << "\"coordinates\": [";
                for(int d=0; d<D-1; ++d)
                    ofs << raw_centers[id ][d] << ",";
                ofs << raw_centers[id][D-1] << "]" << std::endl;
                ofs << "}," << std::endl;
                ofs << "\"properties\": {" << std::endl;
                ofs << "\"prop1\": { \"this\": \"that\" }" << std::endl;
                ofs << "}" << std::endl;
                ofs << "}" << std::endl;
                // Ray
                if(do_ray)
                {
                    ofs << "," << std::endl;
                    ofs << "{" << std::endl;
                    ofs << "\"type\": \"Feature\"," << std::endl;
                    ofs << "\"geometry\": {" << std::endl;
                    ofs << "\"type\": \"LineString\"," << std::endl;
                    ofs << "\"coordinates\": [[";
                    for(int d=0; d<D-1; ++d)
                        ofs << raw_centers[id ][d] << ",";
                    ofs << raw_centers[id][D-1] << "]," << std::endl;
                    ofs << "[";
                    for(int d=0; d<D-1; ++d)
                        ofs << raw_points[id][d] << ",";
                    ofs << raw_points[id][D-1] << "]]" << std::endl;
                    ofs << "}," << std::endl;
                    ofs << "\"properties\": {" << std::endl;
                    if(raw_centers.size() > 0)
                    {
                        ofs << "\"type\":\"ray\"," << std::endl;
                    }
                    ofs << "\"prop1\": { \"this\": \"that\" }" << std::endl;
                    ofs << "}" << std::endl;
                    ofs << "}" << std::endl;
                }
            }
            if(raw_egv.size() > 0)
            {
                ofs << "," << std::endl;
                ofs << "{" << std::endl;
                ofs << "\"type\": \"Feature\"," << std::endl;
                ofs << "\"geometry\": {" << std::endl;
                ofs << "\"type\": \"LineString\"," << std::endl;
                ofs << "\"coordinates\": [[";
                for(int d=0; d<D-1; ++d)
                    ofs << raw_points[id ][d] << ",";
                ofs << raw_points[id][D-1] << "]," << std::endl;
                ofs << "[";
                for(int d=0; d<D-1; ++d)
                    ofs << (raw_points[id ][d] + (raw_sigs[id_sigs+ D-1]/3.0)*raw_egv[id_egv+D+d])  << ",";
                ofs << (raw_points[id][D-1] + (raw_sigs[id_sigs+ D-1]/3.0)*raw_egv[id_egv+D+D-1]) << "]]" << std::endl;
                ofs << "}," << std::endl;
                ofs << "\"properties\": {" << std::endl;
                if(raw_egv.size() > 0)
                {
                    ofs << "\"type\":\"norms\"," << std::endl;
                }
                ofs << "\"prop1\": { \"this\": \"that\" }" << std::endl;
                ofs << "}" << std::endl;
                ofs << "}" << std::endl;
            }
        }
        ofs << "]" << std::endl;
        ofs << "}" << std::endl;
    }






    void extract_egv( std::vector<std::vector<Point> >  & format_egv,bool do_clear = true)
    {
        int D = Traits::D;
        std::vector<double> raw_egv;
        ddt_data<Traits>::dmap[egv_name].extract_raw_uint8_vect(raw_egv,false);
        double coords[D];
        std::vector<Point> act_vect;
        for(int i = 0 ; i < raw_egv.size(); i++)
        {
            coords[i%D] = raw_egv[i];
            if(i%D == D-1)
            {
                act_vect.push_back(ddt_data<Traits>::traits.make_point(coords));
            }
            if(i % ((int)D*D) == ((int)D*D)-1)
            {
                format_egv.push_back(act_vect);
                act_vect.clear();
            }
        }
        if(do_clear)
        {
            raw_egv.clear();
        }
    }

    void fill_egv(std::vector<std::vector<Point> > & format_egv, bool do_clear = true)
    {
        int D = Traits::D;
        std::vector<double> raw_egv;
        for(int i = 0 ; i < format_egv.size(); i++)
        {
            for(int d1 = 0 ; d1 < D; d1++)
            {
                for(int d2 = 0 ; d2 < D; d2++)
                {
                    raw_egv.push_back(format_egv[i][d1][d2]);
                }
            }
        }
        ddt_data<Traits>::dmap[egv_name].fill_full_uint8_vect(raw_egv);
        ddt_data<Traits>::dmap[egv_name].do_exist = true;
        if(do_clear)
            format_egv.clear();
        raw_egv.clear();
    }

    void extract_sigs( std::vector<std::vector<double>>  & format_sigs,bool do_clear = true)
    {
        int D = Traits::D;
        std::vector<double> raw_sigs;
        ddt_data<Traits>::dmap[sig_name].extract_raw_uint8_vect(raw_sigs,false);
        double coords[Traits::D];
        std::vector<double> act_vsig;
        for(int i = 0 ; i < raw_sigs.size(); i++)
        {
            act_vsig.push_back(raw_sigs[i]);
            if(i%D == D-1)
            {
                format_sigs.push_back(act_vsig);
                act_vsig.clear();
            }
        }
        if(do_clear)
            raw_sigs.clear();
    }

    void fill_sigs(std::vector<std::vector<double>>  & format_sigs, bool do_clear = true)
    {
        int D = Traits::D;
        std::vector<double> raw_sigs;
        for(int i = 0 ; i < format_sigs.size(); i++)
        {
            for(int d = 0 ; d < D; d++)
                raw_sigs.push_back(format_sigs[i][d]);
        }
        ddt_data<Traits>::dmap[sig_name].fill_full_uint8_vect(raw_sigs);
        if(do_clear)
        {
            format_sigs.clear();
        }
        raw_sigs.clear();
    }




    void extract_dst( std::vector<std::vector<double>>  & format_dst,bool do_clear = true)
    {
        std::vector<double> raw_dst;
        ddt_data<Traits>::dmap[dst_name].extract_raw_uint8_vect(raw_dst,false);
        std::vector<double> act_vdst;
        for(int i = 0 ; i < raw_dst.size(); i++)
        {
            act_vdst.push_back(raw_dst[i]);
            if(i%3 == 2)
            {
                format_dst.push_back(act_vdst);
                act_vdst.clear();
            }
        }
        if(do_clear)
            raw_dst.clear();
    }



    void fill_dst(std::vector<std::vector<double>>  & format_dst, bool do_clear = true)
    {
        ddt_data<Traits>::dmap[dst_name] = Data_ply(dst_name,"face",3,3,DATA_FLOAT_TYPE);
        std::vector<double> raw_dst;
        for(int i = 0 ; i < format_dst.size(); i++)
        {
            for(int d = 0 ; d < 3; d++)
            {
                raw_dst.push_back(format_dst[i][d]);
            }
        }
        ddt_data<Traits>::dmap[dst_name].fill_full_uint8_vect(raw_dst);
        ddt_data<Traits>::dmap[dst_name].do_exist = true;
        if(do_clear)
        {
            format_dst.clear();
        }
        raw_dst.clear();
    }



    void fill_labs(std::vector<int> & format_labs,bool do_clear = true)
    {
        ddt_data<Traits>::dmap[labseg_name] = Data_ply(labseg_name,"face",1,1,tinyply::Type::INT32);
        std::vector<int> raw_labs;
        for(int i = 0 ; i < format_labs.size(); i++)
        {
            raw_labs.push_back(format_labs[i]);
        }
        ddt_data<Traits>::dmap[labseg_name].fill_full_uint8_vect(raw_labs);
        ddt_data<Traits>::dmap[labseg_name].do_exist = true;
        if(do_clear)
        {
            format_labs.clear();
        }
        raw_labs.clear();
    }

    void extract_labs(std::vector<int> & format_labs,bool do_clear = true)
    {
        int D = Traits::D;
        std::vector<int> raw_labs;
        ddt_data<Traits>::dmap[labseg_name].extract_raw_uint8_vect(raw_labs,false);
        uint num_s = ddt_data<Traits>::dmap[labseg_name].get_nbe_shpt_vect();
        for(int i = 0 ; i < raw_labs.size(); i++)
        {
            format_labs.push_back(raw_labs[i]);
        }
        if(do_clear)
            raw_labs.clear();
    }




    void extract_flags(std::vector<int> & format_flags,bool do_clear = true)
    {
        int D = Traits::D;
        auto ddt_dat = ddt_data<Traits>::dmap[flags_name];
        if(tinyply::Type::INT32 == ddt_dat.get_type() )
        {
            ddt_dat.extract_full_uint8_vect(format_flags,do_clear);
        }
        else
        {
            std::vector<char> raw_flags;
            ddt_data<Traits>::dmap[flags_name].extract_raw_uint8_vect(raw_flags,false);
            uint num_s = ddt_data<Traits>::dmap[flags_name].get_nbe_shpt_vect();
            for(int i = 0 ; i < raw_flags.size(); i++)
            {
                format_flags.push_back((int)raw_flags[i]);
            }
            if(do_clear)
                raw_flags.clear();
        }
    }






    std::vector<Point>  format_points;
    std::vector<Point>  format_centers;
    std::vector<Point>  format_normals;
    std::vector<std::vector<Point> >  format_egv;
    std::vector<std::vector<double>>  format_sigs ;
    std::vector<std::vector<double>>  format_dst ;
    std::vector<int>  format_labs ;
    std::vector<int>  format_flags ;
    std::vector<double>  format_glob_scale ;
    std::vector<int> tile_ids;

    std::vector<std::string> egv_name,sig_name,labseg_name,gid_name,center_name,normal_name,dst_name,glob_scale_name,flags_name;

};

#endif

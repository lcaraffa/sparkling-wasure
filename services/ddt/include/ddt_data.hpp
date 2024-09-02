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

#ifndef DATA_H
#define DATA_H

#include <map>
#include <algorithm>

#include "tinyply.h"

#include "io/number_parser.hpp"
#include "io/logging_stream.hpp"
#include "io/base64.hpp"


#define DATA_FLOAT_TYPE (tinyply::Type::FLOAT64)
#define NB_DIGIT_OUT (5)




template<typename Traits>
class ddt_data
{
public :
    typedef typename Traits::Point                                    Point;
    typedef typename Traits::Delaunay_triangulation                                    Dt;
    typedef typename Traits::Vertex_const_handle                                    Vertex_const_handle;
    typedef typename Traits::Cell_const_handle                                    Cell_const_handle;
    typedef typename Traits::Cell_handle                                    Cell_handle;
    typedef typename Traits::Vertex_handle                                    Vertex_handle;


    class Data_ply
    {
    public:

        Data_ply(std::vector<std::string> vv, std::string pp, int dim, int vs, tinyply::Type tt) : vname(vv),part(pp),D(dim),vsize(vs),do_exist(false),type(tt) {};
        Data_ply() : do_exist(false) {};

        std::vector<std::string> get_name()
        {
            return vname;
        }

        std::string get_name(int ii, bool do_suffix = false)
        {
            if(vname.size() == 1 )
                return vname[0];
            return vname[0]  + "_" + std::to_string(ii);
        }


        tinyply::Type get_type()
        {
            return type;
        }

        bool is_init() const
        {
            return(uint8_vect.size() > 0);
        }


        int get_nbe_uint8_vect() const
        {
            return uint8_vect.size()/((tinyply::PropertyTable[type].stride)*vsize);
        }


        int get_nbe_shpt_vect() const
        {
            if(shpt_vect != nullptr)
                return get_nbe_shpt_vect_shp();
            return 0;
        }


        int get_nbe() const
        {
            if(shpt_vect != nullptr)
                return get_nbe_shpt_vect();
            if(uint8_vect.size() != 0)
                return get_nbe_uint8_vect();
            return 0;
        }

        int get_nbe_shpt_vect_shp() const
        {
            if(shpt_vect == nullptr)
                return 0;
            else
                return shpt_vect->buffer.size_bytes()/((tinyply::PropertyTable[type].stride)*vsize);
        }



        bool has_label(std::string v1) const
        {
            for(auto v2 : vname)
            {
                if( v1 == v2)
                    return true;
            }
            return false;
        }

        void print_elems(std::ostream & ss) const
        {
            int acc=0;
            std::cout << "==============" << std::endl;
            for(auto ee : vname)
            {
                ss << ee << "\t";
            }
            ss << std::endl;
        }

        int get_vsize() const
        {
            return vsize;
        }



        int get_vnbb() const
        {
            return vsize*tinyply::PropertyTable[type].stride;
        }

        int size_bytes()
        {
            if(shpt_vect != nullptr)
                return shpt_vect->buffer.size_bytes();
            else
                return uint8_vect.size();
        }

        void shpt_vect2uint8_vect(bool do_clean = true)
        {
            uint8_vect.resize(size_bytes());
            std::memcpy(uint8_vect.data(), shpt_vect->buffer.get(),size_bytes());
            if(do_clean)
                shpt_vect.reset();
        }

        Point extract_pts(int id)
        {
            return traits.make_point( reinterpret_cast< double * >(shpt_vect->buffer.get())+id*D);
        }



        template<typename DT>
        Point extract_vect(std::vector<DT> & formated_data,int id)
        {
            int vnbb =  get_vnbb();
            formated_data.resize(vsize);
            std::memcpy(formated_data.data(), (shpt_vect->buffer.get())+id*D,vnbb);
        }

        template<typename DT>
        void  extract_value( int id, DT & vv, int i=0) const
        {
            int vnbb =  get_vnbb();
            int szd = get_nbe_shpt_vect_shp()*vsize;
            if(szd > 0)
            {
                vv =  reinterpret_cast<DT &>(shpt_vect->buffer.get()[id*vnbb+i*tinyply::PropertyTable[type].stride]);
            }
            else
            {
                std::memcpy(&vv,&uint8_vect[id*vnbb+i*tinyply::PropertyTable[type].stride],tinyply::PropertyTable[type].stride);
            }
        }


        template<typename DT>
        void extract_full_shpt_vect(std::vector<DT> & formated_data, bool do_clean = true)
        {
            int szd = get_nbe_shpt_vect_shp()*vsize;
            if(szd > 0)
            {
                formated_data.resize(szd);
                std::memcpy(formated_data.data(), shpt_vect->buffer.get(),size_bytes());
                if(do_clean)
                    shpt_vect.reset();
                return;
            }
            szd = get_nbe_uint8_vect()*vsize;
            if(szd > 0)
            {
                formated_data.resize(szd);
                std::memcpy(formated_data.data(), &uint8_vect[0],size_bytes());
                if(do_clean)
                    uint8_vect.clear();
            }
            return;
        }


        template<typename DT>
        void fill_full_uint8_vect(std::vector<DT> & formated_data, bool do_clean = true)
        {
            int szb = sizeof(DT)*formated_data.size();
            uint8_vect.resize(szb);
            std::memcpy(uint8_vect.data(), formated_data.data(),szb);
            if(do_clean)
                formated_data.clear();
            do_exist = true;
        }



        template<typename DT>
        void extract_full_uint8_vect(std::vector<DT> & formated_data, bool do_clean = true)
        {
            int nbe = get_nbe_uint8_vect();
            int szb = sizeof(uint8_t)*uint8_vect.size();
            formated_data.resize(nbe);
            std::memcpy(formated_data.data(),uint8_vect.data(),szb);
            if(do_clean)
            {
                uint8_vect.clear();
                do_exist = false;
            }
        }


        template<typename DT>
        void extract_raw_uint8_vect(std::vector<DT> & formated_data, bool do_clean = true)
        {
            int nbe = get_nbe_uint8_vect()*vsize;
            int szb = sizeof(uint8_t)*uint8_vect.size();
            formated_data.resize(nbe);
            std::memcpy(formated_data.data(),uint8_vect.data(),szb);
            if(do_clean)
            {
                uint8_vect.clear();
                do_exist = false;
            }
        }





        void clean_shpt_vect()
        {
            shpt_vect.reset();
        }

        void clean_uint8_vect()
        {
            uint8_vect.clear();
        }

        void set_exist(bool bb)
        {
            do_exist = bb;
        }

        int get_dim() const
        {
            return D;
        }

        bool do_exist;
        std::string part;
        std::vector<std::string>  vname;
        std::shared_ptr<tinyply::PlyData> shpt_vect;
        std::vector<uint8_t> uint8_vect;
        tinyply::Type type;
    protected :
        int D,vsize;
    };


    std::vector<std::string> subvect(std::vector<std::string> vname, int dd)
    {
        return std::vector<std::string>(vname.begin(),vname.begin()+dd);
    }




    void init_map()
    {
        int D = Traits::D;
        dmap[xyz_name] = Data_ply(xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
        dmap[conf_vertex_name] = Data_ply(conf_vertex_name,"vertex",1,1,DATA_FLOAT_TYPE);
        dmap[simplex_name] = Data_ply(simplex_name,"face",D+1,D+1,tinyply::Type::INT32);
        dmap[nb_name] = Data_ply(nb_name,"face",D+1,D+1,tinyply::Type::INT32);
        dmap[center_name] = Data_ply(center_name,"vertex",D,D,DATA_FLOAT_TYPE);
        dmap[normal_name] = Data_ply(normal_name,"vertex",D,D,DATA_FLOAT_TYPE);
    }



    void init_name()
    {
        int D = Traits::D;
        xyz_name = subvect({"x","y","z","t"},D);
        simplex_name = {"vertex_index"};
        vid_name = {"vid"};
        vtileid_name = {"vtid"};
        ctileid_name = {"ctid"};
        cid_name = {"cid"};
        flag_vertex_name = {"flag_v"};
        conf_vertex_name = {"conf"};
        flag_simplex_name = {"flag_s"};
        nb_name = {"nb_indices"};
        gid_name = {"gid"};
        center_name = subvect({"x_origin","y_origin","z_origin","t_origin"},D);
        normal_name = subvect({"nx","ny","nz","nt"},D);
    }


    ddt_data()
    {
        init_name();
        init_map();
    }




    ddt_data(std::map<std::vector<std::string>, Data_ply > & init_dmap)
    {
        init_name();
        for ( const auto &ee : init_dmap )
        {
            dmap[ee.first] =  Data_ply(ee.first,ee.second.part,D,ee.second.get_vsize(),ee.second.type);
        }
    }


    void clear ()
    {
        for ( const auto &ee : dmap )
        {
            std::vector<uint8_t>().swap(dmap[ee.first].uint8_vect);
        }
    }



    void write_geojson_header(std::ostream & ss)
    {
        ss << "{" << std::endl;
        ss << "\"type\": \"FeatureCollection\"," << std::endl;
        ss << "\"features\": [" << std::endl;
    }

    void write_geojson_footer(std::ostream & ss)
    {
        ss << "]" << std::endl;
        ss << "}" << std::endl;
    }



    void write_geojson_tri(std::ostream & ofs_pts,std::ostream & ofs_spx, bool is_full = true)
    {
        ofs_spx << std::fixed << std::setprecision(5);
        int D = Traits::D;
        std::vector<std::string> lab_color = {"\"red\"","\"green\"","\"blue\""};
        bool is_first = is_full;
        if(is_full)
        {
            write_geojson_header(ofs_pts);
            write_geojson_header(ofs_spx);
        }
        std::vector<double> v_xyz;
        std::vector<int> v_simplex;
        dmap[xyz_name].extract_full_shpt_vect(v_xyz,false);
        dmap[simplex_name].extract_full_shpt_vect(v_simplex,false);
        int nb_pts = v_xyz.size()/D;
        for(int id = 0; id < nb_pts; id++)
        {
            int id_pts = id*D;
            if(!is_first)
                ofs_pts << "," << std::endl;
            is_first=false;
            ofs_pts << "{" << std::endl;
            ofs_pts << "\"type\": \"Feature\"," << std::endl;
            ofs_pts << "\"geometry\": {" << std::endl;
            ofs_pts << "\"type\": \"Point\"," << std::endl;
            ofs_pts << "\"coordinates\": [";
            for(int d=0; d<D-1; ++d)
                ofs_pts << v_xyz[id_pts +d] << ",";
            ofs_pts << v_xyz[id_pts + D-1] << "]" << std::endl;
            ofs_pts << "}," << std::endl;
            ofs_pts << "\"properties\": {" << std::endl;
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part == "vertex" && ee.second.do_exist)
                {
                    for(int nn = 0 ; nn < dmap[ee.first].get_vsize(); nn++)
                    {
                        if(dmap[ee.first].type == tinyply::Type::INT32)
                        {
                            int vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_pts << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                        else  if(dmap[ee.first].type == tinyply::Type::UINT32)
                        {
                            uint vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_pts << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                        else  if(dmap[ee.first].type == DATA_FLOAT_TYPE)
                        {
                            double vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_pts << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                    }
                }
            }
            ofs_pts << "\"datatype\":\"point\"," << std::endl;
            ofs_pts << "\"prop1\": { \"this\": \"that\" }" << std::endl;
            ofs_pts << "}" << std::endl;
            ofs_pts << "}" << std::endl;
        }
        is_first=true;
        uint num_c = dmap[simplex_name].get_nbe();///(D+1);
        for(int id = 0; id < num_c; id++)
        {
            int local = 0;
            bool is_infinite = false;
            for(int i=0; i<=D+1; ++i)
            {
                if(v_simplex[id*(D+1)+(i%(D+1))] == 0)
                    is_infinite = true;
            }
            if(is_infinite)
                continue;
            if(!is_first)
                ofs_spx << "," << std::endl;
            is_first = false;
            ofs_spx << "{" << std::endl;
            ofs_spx << "\"type\": \"Feature\"," << std::endl;
            ofs_spx << "\"geometry\": {" << std::endl;
            ofs_spx << "\"type\": \"Polygon\"," << std::endl;
            ofs_spx << "\"coordinates\": [" << std::endl;
            ofs_spx << "[[";
            for(int i=0; i<=D+1; ++i) // repeat first to close the polygon
            {
                if(i>0)
                {
                    ofs_spx << "],[";
                }
                int id_pp = v_simplex[id*(D+1)+(i%(D+1))];
                for(int d=0; d<D-1; ++d) ofs_spx << v_xyz[id_pp*D + d] << ",";
                ofs_spx << v_xyz[id_pp*D + D-1];
            }
            ofs_spx << "]]";
            ofs_spx << "]";
            ofs_spx << "}," << std::endl;
            ofs_spx << "\"properties\": {" << std::endl;
            switch(local)
            {
            case 0 :
                ofs_spx << "\"fill\":\"red\"," << std::endl;
                break;
            case 1 :
                ofs_spx << "\"fill\":\"green\"," << std::endl;
                break;
            case 2 :
                ofs_spx << "\"fill\":\"blue\"," << std::endl;
                break;
            }
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part == "face" && ee.second.do_exist)
                {
                    for(int nn = 0 ; nn < dmap[ee.first].get_vsize(); nn++)
                    {
                        if(dmap[ee.first].type == tinyply::Type::INT32)
                        {
                            int vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_spx << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                        else  if(dmap[ee.first].type == tinyply::Type::UINT32)
                        {
                            uint vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_spx << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                        else  if(dmap[ee.first].type == DATA_FLOAT_TYPE)
                        {
                            double vv;
                            dmap[ee.first].extract_value(id,vv,nn);
                            ofs_spx << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                        }
                    }
                }
            }
            ofs_spx << "\"stroke-width\":\"2\"," <<  std::endl;
            ofs_spx << "\"prop1\": { \"this\": \"that\" }" << std::endl;
            ofs_spx << "}" << std::endl;
            ofs_spx << "}" << std::endl;
        }
        if(is_full)
        {
            write_geojson_footer(ofs_pts);
            write_geojson_footer(ofs_spx);
        }
    }


    void shpt2uint8()
    {
        for ( const auto &ee : dmap )
        {
            if(ee.second.do_exist)
            {
                if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                        dmap[ee.first].get_nbe_shpt_vect() != 0)
                {
                    dmap[ee.first].shpt_vect2uint8_vect();
                }
            }
        }
    }

    void write_serialized_stream( std::ostream & ss)
    {
        int nn = 0;
        for ( const auto &ee : dmap )
        {
            if(ee.second.do_exist)
            {
                nn++;
                if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                        dmap[ee.first].get_nbe_shpt_vect() != 0)
                {
                    dmap[ee.first].shpt_vect2uint8_vect();
                }
            }
        }
        ss << nn << " ";
        for ( const auto &ee : dmap )
        {
            if(ee.second.do_exist)
            {
                auto vv = dmap[ee.first].uint8_vect;
                ss << dmap[ee.first].vname.size() << " ";
                for(auto nn : dmap[ee.first].vname)
                {
                    ss << nn << " ";
                }
                ss << ee.second.part << " ";
                ss << ee.second.get_vsize() << " ";
                ss << ((int) ee.second.type) << " ";
                serialize_b64_vect(vv,ss);
                ss << " ";
            }
        }
    }

    void read_serialized_stream(std::istream & ss)
    {
        int nbe;
        ss >> nbe;
        for(int i = 0 ; i < nbe; i++)
        {
            std::vector<std::string> data_name;
            std::string tt_name("vertex");
            int dim,vs,dn_size;
            ss >> dn_size;
            std::string nnn;
            for(int i = 0; i < dn_size; i++)
            {
                ss >> nnn;
                data_name.push_back(nnn);
            }
            ss >> tt_name;
            ss >> vs;
            int ttti;
            ss >> ttti;
            dmap[data_name] = Data_ply(data_name,tt_name,dim,vs,static_cast<tinyply::Type>(ttti));
            dmap[data_name].set_exist(true);
            deserialize_b64_vect(dmap[data_name].uint8_vect,ss);
        }
    }

    void write_ply_stream( std::ostream & ss,char nl_char = '\n',bool is_binary = true,bool do_elem_newline = false,bool do_export_bbox = true)
    {
        try
        {
            tinyply::PlyFile file_out;
            std::vector<std::string> & comments_string = file_out.get_comments();
            if(do_export_bbox)
            {
                std::stringstream ss;
                ss << "bbox " << bbox;
                comments_string.push_back(ss.str());
            }
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part == "vertex")
                {
                    if(ee.second.do_exist)
                    {
                        if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                                dmap[ee.first].get_nbe_shpt_vect() != 0)
                        {
                            dmap[ee.first].shpt_vect2uint8_vect();
                        }
                        int nbe = dmap[ee.first].get_nbe_uint8_vect();
                        uint8_t * vv = dmap[ee.first].uint8_vect.data();
                        if(nbe > 0)
                        {
                            file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                               dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::INVALID, 0);
                        }
                    }
                }
            }
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part != "vertex")
                {
                    if(ee.second.do_exist)
                    {
                        if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                                dmap[ee.first].get_nbe_shpt_vect() != 0)
                        {
                            dmap[ee.first].shpt_vect2uint8_vect();
                        }
                        int nbe = dmap[ee.first].get_nbe_uint8_vect();
                        uint8_t * vv = dmap[ee.first].uint8_vect.data();
                        if(nbe > 0)
                        {
                            if(dmap[ee.first].part == "face")
                            {
                                if(ee.first[0] == simplex_name[0] || ee.first[0] == nb_name[0])
                                    file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                                       dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::UINT8, dmap[ee.first].get_vsize());
                                else
                                    file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                                       dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::INVALID, 0);
                            }
                        }
                    }
                }
            }
            file_out.write(ss,is_binary,nl_char,do_elem_newline);
        }
        catch (const std::exception & e)
        {
            std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
        }
    }



    void write_dataset_stream( std::ostream & ss,char nl_char,int tid)
    {
        int nbp = nb_pts();
        std::vector<int> raw_ids_vertex(nbp,tid);
        std::vector<int> raw_ids_simplex(nbp,tid);
        dmap[vtileid_name] = Data_ply(vtileid_name,"vertex",1,1,tinyply::Type::INT32);
        dmap[ctileid_name] = Data_ply(ctileid_name,"face",1,1,tinyply::Type::INT32);
        dmap[vtileid_name].fill_full_uint8_vect(raw_ids_vertex);
        dmap[ctileid_name].fill_full_uint8_vect(raw_ids_simplex);
        dmap[vtileid_name].do_exist = true;
        dmap[ctileid_name].do_exist = true;
        try
        {
            tinyply::PlyFile file_out;
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part == "vertex")
                {
                    if(ee.second.do_exist)
                    {
                        if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                                dmap[ee.first].get_nbe_shpt_vect() != 0)
                        {
                            dmap[ee.first].shpt_vect2uint8_vect();
                        }
                        int nbe = dmap[ee.first].get_nbe_uint8_vect();
                        uint8_t * vv = dmap[ee.first].uint8_vect.data();
                        if(nbe > 0)
                        {
                            file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                               dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::INVALID, 0);
                        }
                    }
                }
            }
            for ( const auto &ee : dmap )
            {
                if(dmap[ee.first].part != "vertex")
                {
                    if(ee.second.do_exist)
                    {
                        if(dmap[ee.first].get_nbe_uint8_vect() == 0 &&
                                dmap[ee.first].get_nbe_shpt_vect() != 0)
                        {
                            dmap[ee.first].shpt_vect2uint8_vect();
                        }
                        int nbe = dmap[ee.first].get_nbe_uint8_vect();
                        uint8_t * vv = dmap[ee.first].uint8_vect.data();
                        if(nbe > 0)
                        {
                            if(dmap[ee.first].part == "face")
                            {
                                if(ee.first[0] == simplex_name[0] || ee.first[0] == nb_name[0])
                                    file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                                       dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::UINT8, dmap[ee.first].get_vsize());
                                else
                                    file_out.add_properties_to_element(dmap[ee.first].part, dmap[ee.first].get_name(),
                                                                       dmap[ee.first].type, nbe, reinterpret_cast<uint8_t*>(vv), tinyply::Type::INVALID, 0);
                            }
                        }
                    }
                }
            }
            file_out.write(ss,false,nl_char,true);
        }
        catch (const std::exception & e)
        {
            std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
        }
    }

    void read_ply_stream(std::istream & ss,char nl_char = '\n')
    {
        try
        {
            tinyply::PlyFile file;
            file.parse_header(ss,nl_char);
            for (auto e : file.get_elements())
            {
                for (auto p : e.properties)
                {
                    std::vector<std::string> pname({p.name});
                    bool do_exist = false;
                    for ( const auto &ee : dmap )
                    {
                        if(ee.second.has_label(p.name))
                        {
                            do_exist = true;
                            dmap[ee.first].do_exist = true;
                            dmap[ee.first].type = p.propertyType;
                            break;
                        }
                    }
                    if(!do_exist)
                    {
                        dmap[pname] = Data_ply(pname,e.name,D,1,p.propertyType);
                        dmap[pname].do_exist = true;
                    }
                }
            }
            for ( const auto &ee : dmap )
            {
                if(ee.second.do_exist)
                {
                    try { dmap[ee.first].shpt_vect = file.request_properties_from_element(ee.second.part, dmap[ee.first].get_name() ); }
                    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
                }
            }
            file.read(ss);
        }
        catch (const std::exception & e)
        {
            std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
        }
    }






    void write_b64_generic_stream( std::ostream & ss,ddt::logging_stream & log)
    {
        int acc = 0;
        for ( const auto &ee : dmap )
        {
            if(ee.second.do_exist)
            {
                acc++;
            }
        }
        ss << acc << " ";
        for ( const auto &ee : dmap )
        {
            if(ee.second.do_exist)
            {
                std::vector<std::string> cname = dmap[ee.first].get_name();
                std::string part = dmap[ee.first].part;
                int dim = dmap[ee.first].get_dim();
                int vsize = dmap[ee.first].get_vsize();
                int tt =  static_cast<int>(dmap[ee.first].type);
                ss << cname.size() << " ";
                for(auto ccn : cname)
                    ss << ccn << " ";
                ss <<  " " << part << " " << dim << " " << vsize << " " << tt << " ";
                write_dataply_stream_fast(dmap[cname],ss,log);
            }
        }
    }

    void read_b64_generic_stream(std::istream & ifile)
    {
        int acc;
        ifile >> acc;
        for(int i = 0; i < acc; i++)
        {
            std::vector<std::string> cname;
            int name_size;
            std::string part;
            int dim,vsize,tt;
            ifile >> name_size;
            for(int j = 0; j < name_size; j++)
            {
                std::string lname;
                ifile >> lname;
                cname.emplace_back(lname);
            }
            ifile >> part >> dim >> vsize >> tt;
            dmap[cname] = Data_ply(cname,part,dim,vsize,static_cast<tinyply::Type>(tt));
            read_b64_data_fast(dmap[cname],ifile);
        }
    };



    tinyply::Type get_int_type()
    {
        return tinyply::Type::INT32;
    }

    tinyply::Type get_float_type()
    {
        return DATA_FLOAT_TYPE;
    }



    std::shared_ptr<tinyply::PlyData> & get_raw_points_ref()
    {
        return dmap[xyz_name].shpt_vect;
    }


    void copy_attribute(ddt_data & wd, int id, std::string ll)
    {
        for ( const auto &ee : wd.dmap )
        {
            if(dmap[ee.first].has_label(ll))
            {
                if(ee.second.do_exist)
                {
                    int vnbb =  ee.second.get_vnbb();
                    for(int i = 0 ; i < vnbb; i++)
                    {
                        dmap[ee.first].uint8_vect.push_back(wd.dmap[ee.first].uint8_vect[id*vnbb+i]);
                    }
                    dmap[ee.first].do_exist = true;
                }
            }
        }
    }





    void insert(ddt_data & wd)
    {
        for ( const auto &ee : wd.dmap )
        {
            if(ee.second.do_exist)
            {
                dmap[ee.first].uint8_vect.insert(dmap[ee.first].uint8_vect.end(),wd.dmap[ee.first].uint8_vect.begin(),wd.dmap[ee.first].uint8_vect.end());
                dmap[ee.first].do_exist = true;
            }
        }
    }




    void copy_point(ddt_data & wd, int id)
    {
        for ( const auto &ee : wd.dmap )
        {
            if(ee.second.do_exist)
            {
                int vnbb =  ee.second.get_vnbb();
                for(int i = 0 ; i < vnbb; i++)
                {
                    dmap[ee.first].uint8_vect.push_back(wd.dmap[ee.first].uint8_vect[id*vnbb+i]);
                }
                dmap[ee.first].do_exist = true;
            }
        }
    }

    void replace_attribute(ddt_data & wd, int id1,int id2)
    {
        for ( const auto &ee : wd.dmap )
        {
            if(ee.second.do_exist)
            {
                int vnbb =  ee.second.get_vnbb();
                for(int i = 0 ; i < vnbb; i++)
                {
                    dmap[ee.first].uint8_vect[id1*vnbb+i] = wd.dmap[ee.first].uint8_vect[id2*vnbb+i];
                }
            }
        }
    }



    Point get_pts(int id,std::vector<std::string> & name)
    {
        std::shared_ptr<tinyply::PlyData> & rp = dmap[name].shpt_vect;
        return traits.make_point( reinterpret_cast< double * >(rp->buffer.get())+id*D);
    }


    Point get_pts(int id)
    {
        std::shared_ptr<tinyply::PlyData> & rp = get_raw_points_ref();
        return traits.make_point( reinterpret_cast< double * >(rp->buffer.get())+id*D);
    }

    int nb_pts_shpt_vect()
    {
        return dmap[xyz_name].get_nbe_shpt_vect();
    }

    int nb_pts_uint8_vect()
    {
        return dmap[xyz_name].get_nbe_uint8_vect();
    }
    int nb_pts()
    {
        return std::max(nb_pts_shpt_vect(),nb_pts_uint8_vect());
    }


    template<typename DT>
    int extract_full_shpt_vect(std::vector<std::string> & name,std::vector<DT> & formated_data, bool do_clean = true)
    {
        dmap[name].extract_full_shpt_vect(formated_data,do_clean);
    }

    template<typename DT>
    int fill_full_uint8_vect(std::vector<std::string> & name,std::vector<DT> & formated_data, bool do_clean = true)
    {
        dmap[name].fill_full_uint8_vect(formated_data,do_clean);
    }

    int nb_simplex_shpt_vect()
    {
        return dmap[simplex_name].get_nbe_shpt_vect();
    }

    int nb_simplex_uint8_vect()
    {
        return dmap[simplex_name].get_nbe_uint8_vect();
    }

    void extract_ptsvect(std::vector<std::string> & name,std::vector<Point> & formated_data, bool do_clean = true)
    {
        int nbv = dmap[name].get_nbe_shpt_vect();
        if(nbv > 0)
        {
            for(int nn = 0; nn < nbv ; nn++)
                formated_data.push_back(get_pts(nn,name));
            if(do_clean)
                dmap[name].clean_shpt_vect();
        }
        nbv = dmap[name].get_nbe_uint8_vect();
        if(nbv > 0)
        {
            dmap[name].extract_full_uint8_vect(formated_data,true);
        }
    }


    Traits traits;
    int D = Traits::D;
    std::string stream_lab;
    std::map<std::vector<std::string>, Data_ply > dmap;
    ddt::Bbox<Traits::D> bbox;
    std::vector<std::string> xyz_name,
        vtileid_name,
        ctileid_name,
        vid_name,
        cid_name,
        flag_vertex_name,
        conf_vertex_name,
        flag_simplex_name,
        gid_name,
        normal_name,
        simplex_name,nb_name,center_name;
};

#endif

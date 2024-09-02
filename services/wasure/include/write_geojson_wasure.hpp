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

#ifndef WRITE_GEOJSON_WASURE_HPP
#define WRITE_GEOJSON_WASURE_HPP


#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "conf_header/conf.hpp"
#include "wasure_data.hpp"

// https://en.wikipedia.org/wiki/GeoJSON
namespace wasure
{


template <typename DTW>
void plot_vertex_pts(typename DTW::Vertex_const_iterator  v,std::ostream & ofs)
{
    using Point = typename DTW::Point;
    ofs.setf( std::ios::fixed, std::ios::floatfield );
    ofs.precision(10);
    Point pt = v->point();
    for(uint i = 0; i < 2; i++)
    {
        ofs << pt[i];
        if(i < 1)
        {
            ofs  << "," ;
        }
    }
}


template <typename DTW>
int dump_2d_surface_geojson(std::vector<typename DTW::Facet_const_iterator> & lft, std::ostream & ofs)
{
    using Vertex_const_iterator = typename DTW::Vertex_const_iterator;
    using Point = typename DTW::Point;
    using Facet = typename DTW::Facet_const_iterator;
    bool is_first = true,is_first2 = true;
    ofs << "{" << std::endl;
    ofs << "\"type\": \"FeatureCollection\"," << std::endl;
    ofs << "\"features\": [" << std::endl;
    for(auto fit = lft.begin();  fit != lft.end(); ++fit)
    {
        Facet ft = *fit;
        std::vector<Point> pts;
        if(!is_first)
            ofs << ",";
        ofs << "{" << std::endl;
        ofs << "\"type\": \"Feature\"," << std::endl;
        ofs << "\"geometry\": {" << std::endl;
        ofs << "\"type\": \"Polygon\"," << std::endl;
        ofs << "\"coordinates\": [" << std::endl;
        ofs << "[[";
        for(int i = 0; i < 3; ++i)
        {
            if(i != ft.index_of_covertex())
            {
                Vertex_const_iterator v = ft.full_cell()->vertex(i);
                if(!is_first2)
                    ofs << "],[";
                plot_vertex_pts<DTW>(v,ofs);
                is_first2 = false;
            }
        }
        ofs << "]]";
        ofs << "]";
        ofs << "}" << std::endl;
        ofs << "}" << std::endl;
        is_first2 = true;
        is_first = false;
    }
    ofs << "]" << std::endl;
    ofs << "}" << std::endl;
    return 0;
}





template<typename Iterator>
void write_geojson_facet_range_wasure(Iterator begin, Iterator end, std::ostream & ofs, bool is_first = true)
{
    typedef typename Iterator::value_type Facet_const_iterator;
    typedef typename Facet_const_iterator::Traits Traits;
    int D = Traits::D;
    for(auto fit = begin; fit != end; ++fit)
    {
        if(!is_first)
            ofs << "," << std::endl;
        auto cit = fit->full_cell();
        int idx = fit->index_of_covertex();
        ofs << "{" << std::endl;
        ofs << "\"type\": \"Feature\"," << std::endl;
        ofs << "\"geometry\": {" << std::endl;
        ofs << "\"type\": \"LineString\"," << std::endl;
        ofs << "\"coordinates\": [[";
        int local = 0;
        int j = 0;
        for(int i=0; i<=D; ++i)
        {
            if(i == idx) continue;
            auto v = cit->vertex(i);
            local += v->is_local();
            for(int d=0; d<D; ++d)
            {
                ofs << v->point()[d];
                if(d < 1)
                    ofs << ",";
            }
            if (++j < D) ofs << "],[";
        }
        ofs << "]]";
        ofs << "}," << std::endl;
        ofs << "\"properties\": {" << std::endl;
        ofs << "\"is_local\": " << int(fit->is_local()) << "," << std::endl;
        ofs << "\"local_score\": " << int(fit->local_score()) << "," << std::endl;
        ofs << "\"is_infinite\": " << int(fit->is_infinite()) << "," << std::endl;
        ofs << "\"is_main\": " << int(fit->is_main()) << "," << std::endl;
        ofs << "\"main_id\": " <<  int(fit->main_id())  << "," << std::endl;
        ofs << "\"prop1\": { \"this\": \"that\" }" << std::endl;
        ofs << "}" << std::endl;
        ofs << "}" << std::endl;
    }
}

template<typename Iterator>
void write_geojson_vert_range_wasure(Iterator begin, Iterator end, std::ostream & ofs, bool is_first = true)
{
    typedef typename Iterator::value_type Vertex_const_iterator;
    typedef typename Vertex_const_iterator::Traits Traits;
    int D = Traits::D;
    for(auto vit = begin; vit != end; ++vit)
    {
        if(vit->is_infinite()) continue;
        if(!is_first)
            ofs << "," << std::endl;
        is_first=false;
        ofs << "{" << std::endl;
        ofs << "\"type\": \"Feature\"," << std::endl;
        ofs << "\"geometry\": {" << std::endl;
        ofs << "\"type\": \"Point\"," << std::endl;
        ofs << "\"coordinates\": [";
        for(int d=0; d<D-1; ++d)
            ofs << vit->point()[d] << ",";
        ofs << vit->point()[D-1] << "]" << std::endl;;
        ofs << "}," << std::endl;
        ofs << "\"properties\": {" << std::endl;
        ofs << "\"fill\":" << (vit->is_local() ? "\"red\"" : "\"blue\"") <<  "," << std::endl;
        ofs << "\"tid\": " << int(vit->tile()->id()) <<  "," << std::endl;
        ofs << "\"main_id\": " <<  int(vit->main_id())  << std::endl;
        ofs << "}" << std::endl;
        ofs << "}" << std::endl;
    }
}






template<typename Iterator>
void write_geojson_cell_range_wasure(Iterator begin, Iterator end, std::ostream & ofs, std::map<Id,wasure_data<Traits> > & w_datas_tri,bool is_first = true)
{
    typedef typename Iterator::value_type Cell_const_iterator;
    std::map<Cell_const_iterator, int> cmap;
    int nextid = 0;
    int D = Traits::D;
    for(auto iit = begin; iit != end; ++iit)
    {
        if(!is_first)
            ofs << "," << std::endl;
        is_first=false;
        ofs << "{" << std::endl;
        ofs << "\"type\": \"Feature\"," << std::endl;
        ofs << "\"geometry\": {" << std::endl;
        ofs << "\"type\": \"Polygon\"," << std::endl;
        ofs << "\"coordinates\": [" << std::endl;
        int local = 0;
        ofs << "[[";
        for(int i=0; i<=D+1; ++i)
        {
            auto v = iit->vertex(i % (D+1));
            if(i>0)
            {
                ofs << "],[";
                local += v->is_local();
            }
            auto p = v->point();
            for(int d=0; d<D-1; ++d) ofs << p[d] << ",";
            ofs << p[D-1];
        }
        int tid = int(iit->tile()->id());
        int id = iit->cell_data().id;
        ofs << "]]";
        ofs << "]";
        ofs << "}," << std::endl;
        ofs << "\"properties\": {" << std::endl;
        ofs << "\"local_score\": " << int(iit->local_score()) << "," << std::endl;
        ofs << "\"is_local\": " << int(iit->is_local()) << "," << std::endl;
        ofs << "\"is_infinite\": " << int(iit->is_infinite()) << "," << std::endl;
        ofs << "\"main_id\": " <<  int(iit->main_id())  << "," << std::endl;
        ofs << "\"tid\": " << tid  << "," << std::endl;
        ofs << "\"data_id\": " << int(iit->cell_data().id)  << "," << std::endl;
        if(!cmap.count(*iit)) cmap[*iit] = nextid++;
        ofs << "\"opt_id\": " << cmap[*iit] << "," << std::endl;
        for(int i = 0 ; i < D+1; i++)
        {
            int iid = -1;
            try
            {
                auto nb0 = iit->neighbor(i);
                auto n0 = nb0->main();
                if(!cmap.count(n0)) cmap[n0] = nextid++;
                iid = cmap[n0];
            }
            catch (...)
            {
            }
            ofs << "\"neigbhor " << i << "\": " << iid << "," << std::endl;
        }
        auto dmap = w_datas_tri[tid].dmap;
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
                        ofs << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                    }
                    else  if(dmap[ee.first].type == tinyply::Type::UINT32)
                    {
                        uint vv;
                        dmap[ee.first].extract_value(id,vv,nn);
                        ofs << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                    }
                    else  if(dmap[ee.first].type == DATA_FLOAT_TYPE)
                    {
                        double vv;
                        dmap[ee.first].extract_value(id,vv,nn);
                        ofs << "\"" << dmap[ee.first].get_name(nn,true) << "\":" << vv << "," << std::endl;
                    }
                }
            }
        }
        ofs << "\"prop1\": { \"this\": \"that\" }" << std::endl;
        ofs << "}" << std::endl;
        ofs << "}" << std::endl;
    }
}

template<typename DDT>
void write_geojson_tri_wasure(const DDT& ddt,   std::map<Id,wasure_data<Traits> > & w_datas_tri, std::ostream & ofs)
{
    ofs << "{" << std::endl;
    ofs << "\"type\": \"FeatureCollection\"," << std::endl;
    ofs << "\"features\": [" << std::endl;
    write_geojson_vert_range_wasure(ddt.vertices_begin(), ddt.vertices_end(), ofs,true);
    write_geojson_cell_range_wasure(ddt.cells_begin(), ddt.cells_end(), ofs,w_datas_tri,false);
    write_geojson_facet_range_wasure(ddt.facets_begin(), ddt.facets_end(), ofs,false);
    ofs << "]" << std::endl;
    ofs << "}" << std::endl;
}


template<typename DDT>
void write_geojson_tri_wasure(const DDT& ddt, std::ostream & ofs,wasure_data<typename DDT::Traits> ddtd)
{
    ofs << "{" << std::endl;
    ofs << "\"type\": \"FeatureCollection\"," << std::endl;
    ofs << "\"features\": [" << std::endl;
    write_geojson_vert_range_wasure(ddt.vertices_begin(), ddt.vertices_end(), ofs,true);
    write_geojson_cell_range_wasure(ddt.cells_begin(), ddt.cells_end(), ofs,ddtd,false);
    write_geojson_facet_range_wasure(ddt.facets_begin(), ddt.facets_end(), ofs,false);
    ofs << "]" << std::endl;
    ofs << "}" << std::endl;
}




}



#endif // WRITE_GEOJSON_WASURE_HPP

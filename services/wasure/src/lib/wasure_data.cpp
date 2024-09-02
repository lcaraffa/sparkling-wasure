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

#include "wasure_data.hpp"
#include "utils/utils.hpp"
#include "io/iqlib_io.hpp"
#include "tinyply.h"
#include "spark_api.hpp"
#include "typedefs_type.hpp"
#include <CGAL/Epick_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/predicates_d.h>

int wasure_data::dump_surface(std::vector<Facet> & lft, int nblabs, std::string fileName)
{
    int dim = 3;
    std::map<typename DTW::Vertex_handle, uint> vertex_map;
    std::vector<Point> point_list;
    int acc = 0;
    Vertex_iterator fvit;
    for(auto fit = lft.begin();  fit != lft.end(); ++fit)
    {
        Facet ft = *fit;
        Full_cell_handle fch = ft.full_cell();
        int idx = ft.index_of_covertex();
        Full_cell_handle fchn = fch->neighbor(idx);
        for(int i = 0; i < dim+1; ++i)
        {
            if(i != ft.index_of_covertex())
            {
                Vertex_handle v = ft.full_cell()->vertex(i);
                if(vertex_map.find(v) == vertex_map.end())
                {
                    vertex_map[v] = acc++;
                    point_list.push_back(v->point());
                }
            }
        }
    }
    int nbv = vertex_map.size();
    int nbf = lft.size();
    std::ofstream fo;
    fo.open (fileName.c_str());
    std::cout << "\t Writing " << fileName << "..." << std::endl;
    fo << "ply" << std::endl;
    fo << "format ascii 1.0" << std::endl;
    fo << "comment VCGLIB generated" << std::endl;
    fo << "element vertex " << nbv << std::endl;
    fo << "property float x" << std::endl;
    fo << "property float y" << std::endl;
    fo << "property float z" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "element face " << nbf << std::endl;
    fo << "property list uchar int vertex_indices" << std::endl;
    fo << "property uchar red" << std::endl;
    fo << "property uchar green" << std::endl;
    fo << "property uchar blue" << std::endl;
    fo << "end_header" << std::endl;
    fo << std::setprecision(12);
    std::cout << "wasure_data:dumping point...." << std::endl;
    for(auto const &pi : point_list)
    {
        fo << pi[0] << " " <<  pi[1]  << " " <<  pi[2]  << " " <<  255  <<  " " << 255 << " " << 255 << std::endl;
    }
    std::cout << "wasure_data:dumping facets...." << std::endl;
    for(auto fit = lft.begin();  fit != lft.end(); ++fit)
    {
        Facet ft = *fit;
        int acc = 0;
        Full_cell_handle fch = ft.full_cell();
        int idx = ft.index_of_covertex();
        Full_cell_handle fchn = fch->neighbor(idx);
        int lab = fch->data().lab;
        int nlab = fchn->data().lab;
        double plab = ((double)lab/((double)(nblabs-1)));
        double pnlab = ((double)nlab/((double)(nblabs-1)));
        std::vector<int> lid;
        std::vector<Point> lp;
        for(int i = 0 ; i < dim +1; i++)
        {
            if(i != ft.index_of_covertex())
            {
                lid.push_back(vertex_map[ft.full_cell()->vertex(i)]);
                lp.push_back(ft.full_cell()->vertex(i)->point());
            }
        }
        lid.push_back(ft.index_of_covertex());
        lp.push_back(ft.full_cell()->vertex(ft.index_of_covertex())->point());
        assert(lid.size() == 4 && lp.size() == 4);
        int cr = 200;
        int cg = fabs(plab - pnlab)*255;
        int cb = fabs(plab - pnlab)*255;
        bool is_inf = false;
        for(auto ppp : lp)
        {
            if(ppp.dimension() < dim)
            {
                std::cout << "WARNING, infinit point found" << std::endl;
                is_inf = true;
                break;
            }
        }
        int o1 = !is_inf ? CGAL::orientation(lp.begin(),lp.end()) :  0 ;
        bool bl =  ((o1 == -1 && pnlab >= 0.5) || ( o1 == 1 && pnlab < 0.5));
        if(bl)
            fo << "3 " << lid[0] << " " << lid[2] << " " << lid[1] << " " << cr << " " << cg << " " << cb << std::endl;
        else
            fo << "3 " << lid[0] << " " << lid[1] << " " << lid[2] << " " << cr << " " << cg << " " << cb << std::endl;
    }
    std::cout << "dumping ok" << std::endl;
    fo.close();
    return 0;
}


int wasure_data::ech_data(double ech)
{
    source_label = ech_vector(source_label,ech);
    points = ech_vector(points,ech);
    centers_pts = ech_vector(centers_pts,ech);
    norms_pts = ech_vector(norms_pts,ech);
    glob_scale = ech_vector(glob_scale,ech);
    dims_norms = ech_vector(dims_norms,ech);
    dims_scales = ech_vector(dims_scales,ech);
}



int wasure_data::ech_data(std::vector<bool> & bv)
{
    source_label = ech_vector(source_label,bv);
    points = ech_vector(points,bv);
    centers_pts = ech_vector(centers_pts,bv);
    norms_pts = ech_vector(norms_pts,bv);
    glob_scale = ech_vector(glob_scale,bv);
    dims_norms = ech_vector(dims_norms,bv);
    dims_scales = ech_vector(dims_scales,bv);
}



void dump_vrt_norms(std::string name_out_vrt,std::string name_out_csv)
{
    std::ofstream f;
    f.open(name_out_vrt.c_str());
    if(!f)
    {
        std::cerr << "ERROR :" << name_out_vrt << " cannot be open" << std::endl;
        return;
    }
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << remove_path(remove_extension(name_out_csv)) <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << remove_path(name_out_csv) << "</SrcDataSource>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbLineString</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"idx\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
}


void wasure_data::dump2qgis(std::string namefile)
{
    std::string name_out_norms_csv = remove_extension(namefile) + "_norms.csv";
    std::string name_out_norms_vrt = remove_extension(namefile) + "_norms.vrt";
    dump_vrt_norms(name_out_norms_vrt,name_out_norms_csv);
    uint D = 2;
    std::ofstream fo;
    fo.open (name_out_norms_csv.c_str());
    fo << "idx,geom"<< std::endl;
    fo.precision(10);
    for(int i = 0 ; i < points.size(); i++)
    {
        fo << "1,\"MULTILINESTRING((";
        fo << points[i][0] << " " << points[i][1] << ",";
        fo << points[i][0] + dims_norms[i][1][0]*dims_scales[i][1]/3.0 << " " << points[i][1] + dims_norms[i][1][1]*dims_scales[i][1]/3.0 ;
        fo << "))\""  << std::endl;
    }
}


inline  std::ostream & wasure_data::dump_double(double dd,std::ostream & ofile)
{
    ofile << dd << " ";
}

inline  std::ostream & wasure_data::dump_int(int ii,std::ostream & ofile)
{
    ofile << ii << " ";
}



inline double wasure_data::read_double(std::istream & ifile)
{
    double dd;
    ifile >> dd;
    return dd;
}


inline int wasure_data::read_int(std::istream & ifile)
{
    int ii;
    ifile >> ii;
    return ii;
}



std::ostream & wasure_data::dump_point(Point & pp,std::ostream & ofile)
{
    for(int d = 0 ; d < D; d++)
        dump_double(pp[d],ofile);
}



Point wasure_data::read_point(std::istream & ifile)
{
    std::vector<double> coords(D);
    for(int d = 0 ; d < D; d++)
        ifile >> coords[d];
    return Point(coords);
}



std::ostream & wasure_data::dump_stream(std::ostream & ofile)
{
    ofile << std::setprecision(NB_DIGIT_OUT);
    ofile << D << " " << points.size() << " ";
    for(auto pp : source_label)
    {
        dump_int(pp,ofile);
    }
    for(auto pp : points)
    {
        dump_point(pp,ofile);
    }
    for(auto pp : centers_pts)
    {
        dump_point(pp,ofile);
    }
    for(auto dd : glob_scale)
    {
        dump_double(dd,ofile);
    }
    for(auto lpts : dims_norms)
    {
        for(auto pp : lpts)
        {
            dump_point(pp,ofile);
        }
    }
    for(auto lpts : dims_scales)
    {
        for(auto dd : lpts)
        {
            dump_double(dd,ofile);
        }
    }
}


std::ostream & wasure_data::dump_xyz(std::ostream & ofile)
{
    ofile << std::setprecision(NB_DIGIT_OUT);
    for(int i = 0 ; i < points.size() ; i++)
    {
        dump_int(source_label[i],ofile);
        Point & p1 = points[i];
        dump_point(p1,ofile);
        if(centers_pts.size() > 0)
        {
            Point & p2 = centers_pts[i];
            dump_point(p2,ofile);
        }
        if(dims_norms.size() > 0)
        {
            for(int dn = 0; dn < dims_norms[i].size(); dn++)
            {
                for(int d = 0 ; d < D; d++)
                {
                    Point & p = dims_norms[i][dn];
                    dump_point(p,ofile);
                }
            }
        }
        if(dims_scales.size() > 0)
        {
            for(int ds = 0; ds < dims_scales[i].size(); ds++)
            {
                dump_double(dims_scales[i][ds],ofile);
            }
        }
        if(glob_scale.size() > 0)
        {
            double & gs = glob_scale[i];
            dump_double(gs,ofile);
        }
        ofile << "\n";
    }
}



std::istream & wasure_data::read_seq_stream(std::istream & ifile)
{
    int nbp;
    ifile >> nbp;
    for(int i = 0; i < nbp; i ++)
        read_stream(ifile);
    return ifile;
}


std::istream & wasure_data::read_iq_stream(std::istream & ifile, int & tile_id)
{
    stream_data_header iqh;
    iqh.parse_header(std::cin);
    tile_id = iqh.get_id(0);
    read_stream(std::cin);
}


std::ostream & wasure_data::write_iq_stream(std::ostream & ofile, int tile_id)
{
    stream_data_header iqh("p","s",tile_id);
    iqh.write_header(std::cout);
    dump_stream(std::cout);
}

std::istream & wasure_data::read_stream(std::istream & ifile)
{
    int nbp = -1,nbc;
    ifile >> D >> nbp;
    if(nbp == -1)
        return ifile;
    for(int i = 0; i < nbp; i++)
    {
        source_label.push_back(read_int(ifile));
    }
    for(int i = 0; i < nbp; i++)
    {
        points.push_back(read_point(ifile));
    }
    for(int i = 0; i < nbp; i++)
    {
        centers_pts.push_back(read_point(ifile));
    }
    for(int i = 0; i < nbp; i++)
    {
        glob_scale.push_back(read_double(ifile));
    }
    for(int i = 0; i < nbp; i++)
    {
        std::vector<Point> ll;
        for(int d = 0; d < D; d++)
        {
            ll.push_back(read_point(ifile));
        }
        dims_norms.push_back(ll);
    }
    for(int i = 0; i < nbp; i++)
    {
        std::vector<double> ll;
        for(int d = 0; d < D; d++)
        {
            ll.push_back(read_double(ifile));
        }
        dims_scales.push_back(ll);
    }
}


std::istream & wasure_data::read_stream_raw(std::istream & ifile)
{
    int lab = -1;
    while(ifile >> lab)
    {
        if(lab == -1)
            return ifile;
        source_label.push_back(lab);
        points.push_back(read_point(ifile));
        centers_pts.push_back(read_point(ifile));
        glob_scale.push_back(read_double(ifile));
        std::vector<Point> llp;
        for(int d = 0; d < D; d++)
        {
            llp.push_back(read_point(ifile));
        }
        dims_norms.push_back(llp);
        std::vector<double> lld;
        for(int d = 0; d < D; d++)
        {
            lld.push_back(read_double(ifile));
        }
        dims_scales.push_back(lld);
    }
}

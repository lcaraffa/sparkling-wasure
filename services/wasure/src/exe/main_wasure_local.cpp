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

#include <unordered_map>

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

typedef std::map<Id,wasure_data<Traits> > D_MAP;
typedef std::map<Id,std::list<wasure_data<Traits>> > D_LMAP;
typedef typename Traits::Vertex_const_handle     Vertex_const_handle;
typedef typename Traits::Cell_const_handle     Cell_const_handle;


void seg_tri(Cell_const_handle & cc, int id_seg,int lab_seg,
             Vertex_const_handle & vv,
             std::unordered_map<int,int> & id_map,Tile_iterator & tile, std::vector<int> & format_labs)
{
    int cid = tile->lid(cc);
    int ch1lab = format_labs[cid];
    auto id_seg_cc = id_map.find(cid);
    if(ch1lab != lab_seg || id_seg_cc != id_map.end() )
        return;
    id_map[cid] = id_seg;
    for(int i = 0; i < 4; i++)
    {
        auto nn = tile->neighbor(cc,i);
        int vid = cc->index(vv);
        if(i == vid)
            continue;
        seg_tri(nn,id_seg,lab_seg,vv,id_map,tile,format_labs);
    }
    return;
}

std::string time_in_HH_MM_SS_MMM()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << std::put_time(&bt, "%d-%m-%Y-%H-%M-%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}

void init_local_ids( DTW & tri1)
{
    int acc = 0;
    for(auto iit = tri1.cells_begin(); iit != tri1.cells_end(); ++iit)
    {
        const Data_C & cd = iit->cell_data();
        Data_C & cd_quickndirty = const_cast<Data_C &>(cd);
        cd_quickndirty.id = acc;
        cd_quickndirty.gid = acc++;
    }
}


int main(int argc, char **argv)
{
    wasure_params params;
    params.parse(argc,argv);
    Scheduler sch(1);
    Traits traits;
    int D = Traits::D;
    Id tid = 0;
    ddt::logging_stream log(std::to_string(tid) + "_" + params.algo_step, params.log_level);
    wasure_data<Traits> wdp;
    wasure_data<Traits> wds;
    std::ifstream ifile;
    if(!boost::filesystem::exists(params.filename))
    {
        std::cerr << params.filename << " does not exist"  << std::endl;
        return 1;
    }
    else
    {
        std::cout << "Start processing " << params.filename <<  " ..." << std::endl;
    }
    ifile.open(params.filename);
    wdp.read_ply_stream(ifile);
    ifile.close();
    wdp.shpt2uint8();
    int count = wdp.nb_pts_uint8_vect();
    wdp.dmap[wdp.xyz_name].extract_full_uint8_vect(wdp.format_points);
    wdp.dmap[wdp.center_name].extract_full_uint8_vect(wdp.format_centers);
    wdp.dmap[wdp.normal_name].extract_full_uint8_vect(wdp.format_normals);
    wdp.format_flags.clear();
    wdp.format_flags.resize(count,0);
    std::string filename_cen(params.output_dir + "/centrs.xyz");
    std::string filename_dim(params.output_dir + "/dim.ply");
    std::string filename_tes(params.output_dir + "/tessel.ply");
    std::string filename_dst(params.output_dir + "/dst.ply");
    std::string filename_tri(params.output_dir + "/tri.stream");
    std::ofstream ofile;
    ofile.open(filename_cen);
    for(auto cc : wdp.format_centers)
    {
        ofile << cc << std::endl;
    }
    ofile.close();
    wasure_algo w_algo;
    if(boost::filesystem::exists(filename_dim) &&
            boost::filesystem::exists(filename_tes)
      )
    {
        ifile.open(filename_dim);
        wdp.read_ply_stream(ifile);
        wdp.shpt2uint8();
        ifile.close();
        ifile.open(filename_tes);
        wds.read_ply_stream(ifile);
        wds.shpt2uint8();
        ifile.close();
        wdp.dmap[wdp.xyz_name].extract_full_uint8_vect(wdp.format_points);
        wdp.extract_sigs(wdp.format_sigs);
        wdp.extract_egv(wdp.format_egv);
        wds.dmap[wdp.xyz_name].extract_full_uint8_vect(wds.format_points);
    }
    else
    {
        std::cout << "compute dim ..." << std::endl;
        w_algo.compute_dim(wdp.format_points,
                           wdp.format_egv,
                           wdp.format_sigs,log);
        std::cout << "flip dim ..." << std::endl;
        w_algo.flip_dim_ori(wdp.format_points,
                            wdp.format_egv,
                            wdp.format_centers);
        std::cout << "start tessel ..." << std::endl;
        w_algo.tessel_adapt(wdp.format_points,
                            wds.format_points,
                            wdp.format_egv,
                            wdp.format_sigs,
                            wdp.format_normals,
                            NUM_ITER_TESSEL,params.pscale,D,tid
                           );
        wdp.fill_egv(wdp.format_egv,false);
        wdp.fill_sigs(wdp.format_sigs,false);
        wdp.dmap[wdp.xyz_name].fill_full_uint8_vect(wdp.format_points,false);
        wds.dmap[wds.xyz_name].fill_full_uint8_vect(wds.format_points,false);
        ofile.open(filename_dim);
        wdp.write_ply_stream(ofile,'\n',true);
        ofile.close();
        ofile.open(filename_tes);
        wds.write_ply_stream(ofile,'\n',true);
        ofile.close();
    }
    DTW tri1;
    D_MAP w_datas_tri;
    w_datas_tri[tid] = wasure_data<Traits>();
    tri1.init(tid);
    int nbs;
    if(boost::filesystem::exists(filename_tri))
    {
        ifile.open(filename_tri);
        bool do_clean_data = false;
        bool do_serialize = false;
        read_ddt_stream(tri1,w_datas_tri[tid],ifile,tid,do_serialize,do_clean_data,log);
        w_datas_tri[tid].extract_dst(w_datas_tri[tid].format_dst,false);
        auto tile = tri1.get_tile(tid);
        nbs = tile->number_of_cells();
        init_local_ids(tri1);
        tri1.finalize(sch);
    }
    else
    {
        std::vector<Point_id>  vp;
        Tile_iterator tci = tri1.get_tile(tid);
        for(auto pp : wds.format_points)
        {
            vp.emplace_back(std::make_pair(pp,tid));
        }
        int nbi1 = tci->insert(vp,false);
        tri1.finalize(sch);
        DT & tri_tile  = tri1.get_tile(tid)->triangulation();
        auto tile = tri1.get_tile(tid);
        std::vector<std::vector<double>>  & format_dst = w_datas_tri[tid].format_dst; ;
        nbs = tile->number_of_cells();
        // Init each simplex at "unknown"
        // 0 0 1 => 0% in, 0% out, 100% unknown
        if(format_dst.size() == 0)
        {
            for(int ss = 0; ss < nbs ; ss++)
            {
                format_dst.push_back(std::vector<double>({0.0,0.0,1.0}));
            }
        }
        // ===== Init the id of each cell
        init_local_ids(tri1);
        std::cout << "compute dst ..." << std::endl;
        w_algo.compute_dst_with_center(tri1,w_datas_tri[tid],wdp,params,tid);
        w_datas_tri[tid].fill_dst(w_datas_tri[tid].format_dst,false);
        std::ofstream ofile;
        ofile.open(filename_tri);
        ddt::write_ddt_stream(tri1, w_datas_tri[tid], ofile,tid,false,log);
        ofile.close();
    }
    // ===== Segmentation =====
    std::vector<int>  & format_labs = w_datas_tri[tid].format_labs ;
    format_labs.resize(nbs);
    tbmrf_reco<DTW,D_MAP> mrf(params.nb_labs,&tri1,&w_datas_tri);
    std::vector<double> lambda_list({0.01,0.05,0.07,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.9,1,1.5,2,4,8,10,20});
    std::vector<int> opt_mode({0});
    // Mode 0 => outdoor scene
    // Mode 1 => indoor scene
    mrf.set_mode(0);
    std::cout << "regularize ..." << std::endl;
    for(auto opt_id : opt_mode)
    {
        for(auto ll : lambda_list)
        {
            mrf.lambda = ll;
            for(int ii = 0; ii < nbs; ii++)
                format_labs[ii] = 0;
            if(opt_id == 0)
            {
                mrf.opt_gc(1,tri1,w_datas_tri);
            }
            else
                mrf.opt_belief(1,tri1,w_datas_tri);
            w_datas_tri[tid].fill_labs(w_datas_tri[tid].format_labs);
            if(D == 2)
            {
                DT & tri_tile  = tri1.get_tile(tid)->triangulation();
                traits.export_tri_to_data(tri_tile,w_datas_tri[tid]);
                ddt::stream_data_header oqh_1("p","s",tid),oqh_2("p","s",tid);
                std::string filename(params.output_dir +  "/" + params.slabel + "_id_" + std::to_string(tid) + "_seg");
                oqh_1.write_into_file(filename,"_pts.geojson");
                oqh_1.write_header(std::cout);
                oqh_2.write_into_file(filename,"_spx.geojson");
                oqh_2.write_header(std::cout);
                w_datas_tri[tid].write_geojson_tri(oqh_1.get_output_stream(),oqh_2.get_output_stream());
                oqh_1.finalize();
                oqh_2.finalize();
                ddt::add_qgis_style(oqh_2.get_file_name(),"tri_seg.qml");
                std::cout << std::endl;
            }
            std::vector<Facet_const_iterator> lft;
            std::vector<bool> lbool;
            mrf.extract_surface(tid,lft,w_datas_tri);
            std::string string_name = time_in_HH_MM_SS_MMM();
            std::string ply_name(params.output_dir +  "/" + params.slabel + "_" + string_name + "_ll_" + std::to_string(ll) + "_ww_" + std::to_string(params.ray_weight) + "_ss_" + std::to_string(params.rat_ray_sample) +  "_surface_" + std::to_string(opt_id));
            ddt::stream_data_header oth("p","f",tid);
            if(D == 2)
                oth.write_into_file(ply_name,".geojson");
            else
                oth.write_into_file(ply_name,".ply");
            oth.write_header(std::cout);
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
                std::vector<int> vid;
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
                for(auto fit = lft.begin(); fit != lft.end(); ++fit)
                {
                    Cell_const_iterator fch = fit->full_cell();
                    int id_cov = fit->index_of_covertex();
                    Cell_const_iterator fchn = fch->neighbor(id_cov);
                    int id_covn = fch->mirror_index(id_cov);
                    int cccid = fch->lid();
                    int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
                    if(fch->is_infinite())
                    {
                        fch = fchn;
                        id_cov = id_covn;
                    }
                    const Point& a = fch->vertex((id_cov+1)&3)->point();
                    const Point& b = fch->vertex((id_cov+2)&3)->point();
                    const Point& c = fch->vertex((id_cov+3)&3)->point();
                    const Point& d = fch->vertex((id_cov)&3)->point();
                    bool bl =
                        (CGAL::orientation(a,b,c,d) == 1 && ch1lab == 0) ||
                        (CGAL::orientation(a,b,c,d) == -1 && ch1lab == 1);
                    Id ida = (id_cov+1)&3;
                    Id idb = (id_cov+2)&3;
                    Id idc = (id_cov+3)&3;
                    vid.push_back(fch->lid());
                    v_simplex.push_back(vertex_map[fch->vertex(ida)]);
                    if(!bl)
                    {
                        v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                        v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                    }
                    else
                    {
                        v_simplex.push_back(vertex_map[fch->vertex(idc)]);
                        v_simplex.push_back(vertex_map[fch->vertex(idb)]);
                    }
                }
                datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",D,D,DATA_FLOAT_TYPE);
                datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",D,D,tinyply::Type::INT32);
                datas_out.dmap[std::vector<std::string> {"id"}] = ddt_data<Traits>::Data_ply(std::vector<std::string> {"id"},"face",1,1,tinyply::Type::INT32);
                datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
                datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
                datas_out.dmap[std::vector<std::string> {"id"}].fill_full_uint8_vect(vid);
                datas_out.write_ply_stream(oth.get_output_stream(),'\n',true);
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
        }
    }
    return 0;
}

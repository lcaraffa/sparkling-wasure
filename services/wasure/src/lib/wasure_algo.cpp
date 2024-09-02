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
#include "wasure_algo.hpp"

#include <algorithm>
#include <random>

#include "ANN/ANN.h"
#include "wasure_maths.hpp"
#include "input_params.hpp"
#include "wasure_typedefs.hpp"

template <typename T>
auto normalize(T const& V)
{
    auto const slen = V.squared_length();
    auto const d = CGAL::approximate_sqrt(slen);
    return V / d;
}


void wasure_algo::tessel_adapt(std::vector<Point> & points,std::vector<Point> & vps,std::vector<std::vector<Point>> & norms,std::vector<std::vector<double>> & scales,std::vector<Point> & v_los, int maxit, double target_err, int D, int tid)
{
    std::vector<int> lidx(points.size());
    std::iota(lidx.begin(), lidx.end(), 0);
    int nb_inserted = 0;
    random_shuffle(std::begin(lidx), std::end(lidx));
    typedef DT_raw::Locate_type    Locate_type;
    typedef DT_raw::Point          Point;
    DT_raw  tri = traits_raw.triangulation(D) ;
    for(int i = 0; i < points.size(); i++)
    {
        int pidx = lidx[i];
        Point p1 = points[pidx];
        std::vector<Point> & pts_norm = norms[pidx];
        std::vector<double> & pts_scale = scales[pidx];
        Locate_type lt;
        int li, lj;
        auto loc = tri.locate(p1,lt, li, lj);
        if(lt != DT_raw::CELL)
        {
            tri.insert(p1,loc);
            nb_inserted++;
            continue;
        }
        bool do_insert = true;
        for(uint dd = 0; dd < D+1; dd++)
        {
            Point & pii = loc->vertex(dd)->point();
            std::vector<double> pii_coefs = compute_base_coef<Point>(p1,pii,pts_norm,D);
            int is_close_enough = 0;
            for(int d = 0; d < D; d++)
            {
                if(fabs(pii_coefs[d]) < pts_scale[d]*target_err)
                {
                    is_close_enough++;
                    break;
                }
            }
            if(is_close_enough > 0  )
            {
                do_insert = false;
                break;
            }
        }
        if(do_insert)
        {
            tri.insert(p1,loc);
            nb_inserted++;
        }
    }
    tessel(tri,points,vps,norms,scales,v_los,maxit,tid);
}


int dump_tessel(DT_raw  & tri, int it, int max_it,int tid, double * bbox_min,double * bbox_max)
{
    Traits_raw traits_raw;
    int nbp_face = 10;;
    typedef typename DT_raw::Vertex_handle                            Vertex_const_handle;
    std::vector<std::vector<double>> l_cir;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        auto tmp_fch = fit->first;
        int tmp_idx = fit->second;
        auto tmp_fchn = tmp_fch->neighbor(tmp_idx);
        auto bb1 = traits_raw.circumcenter(tri,tmp_fch);
        auto bb2 = traits_raw.circumcenter(tri,tmp_fchn);
        for(int i = 0; i < nbp_face; i++)
        {
            std::vector<double> pp;
            bool do_insert = true;
            for(int j = 0; j < 3; j++)
            {
                pp.push_back((bb1[j] + (bb2[j]-bb1[j])*(i/((double)nbp_face))));
            }
            for(int j = 0; j < 3; j++)
            {
                if(pp[j] < bbox_min[j] || pp[j] > bbox_max[j])
                {
                    do_insert = false;
                }
            }
            if(do_insert)
                l_cir.push_back(pp);
        }
    }
    l_cir.clear();
    std::ofstream myfile;
    std::string filename("/tmp/tessel_" + std::to_string(it) + "_" + std::to_string(tid) + ".ply");
    myfile.open (filename);
    myfile << "ply" <<  std::endl;
    myfile << "format ascii 1.0" << std::endl;
    myfile << "element vertex "  << tri.number_of_vertices() + l_cir.size()  << std::endl;
    myfile << "property float x" << std::endl;
    myfile << "property float y" << std::endl;
    myfile << "property float z" << std::endl;
    myfile << "property uchar red " << std::endl;
    myfile << "property uchar green" << std::endl;
    myfile << "property uchar blue" << std::endl;
    myfile << "element face " << tri.number_of_finite_facets() << std::endl;
    myfile << "property list uchar int vertex_index " << std::endl;
    myfile << "end_header                " << std::endl;
    double ccol = ((double)it)/((double)max_it);
    CGAL::Unique_hash_map<Vertex_const_handle, int> vertex_map;
    int acc = 0;
    for(auto vv = traits_raw.vertices_begin(tri); vv != traits_raw.vertices_end(tri) ; ++vv)
    {
        if(tri.is_infinite(vv))
            continue;
        myfile << vv->point() << " " << ((int)(255*ccol)) << " " << ((int)(255*ccol)) << " " << ((int)(255*ccol)) << std::endl;
        vertex_map[vv] = acc++;
    }
    for(auto pp : l_cir)
        myfile << pp[0] << " " << pp[1] << " " << pp[2] << " 0 255 0" << std::endl;;
    for(auto fit = tri.facets_begin();  fit != tri.facets_end(); ++fit)
    {
        if(tri.is_infinite(*fit))
            continue;
        auto fch = fit->first;
        int id_cov = fit->second;
        myfile << "3 ";
        for(int i = 0; i < 4; ++i)
        {
            if(i != id_cov)
            {
                auto v = fch->vertex(i);
                myfile << vertex_map[v] << " ";
            }
        }
        myfile << std::endl;
    }
    myfile.close();
    return 0;
}




int dump_vector_pts(std::vector<Point> vp, int it, int tid)
{
    std::ofstream myfile;
    std::string filename("/tmp/extra_" + std::to_string(it) + "_" + std::to_string(tid) + ".ply");
    myfile.open (filename);
    myfile << "ply" <<  std::endl;
    myfile << "format ascii 1.0" << std::endl;
    myfile << "element vertex "  << vp.size() << std::endl;
    myfile << "property float x" << std::endl;
    myfile << "property float y" << std::endl;
    myfile << "property float z" << std::endl;
    myfile << "end_header                " << std::endl;
    for(auto vv : vp)
    {
        myfile << vv << " " << std::endl;
    }
    myfile.close();
    return 0;
}


int
wasure_algo::tessel(DT_raw  & tri,
                    std::vector<Point> & points,  std::vector<Point> & vps,
                    std::vector<std::vector<Point> > & norms, std::vector<std::vector<double>> & scales, std::vector<Point> & v_los,int max_it, Id tid)
{
    typedef typename DT_raw::Vertex_handle                            Vertex_const_handle;
    double bbox_min[Traits::D];
    double bbox_max[Traits::D];
    for(int d = 0; d < D; d++)
    {
        bbox_min[d] = 100000;
        bbox_max[d] = -100000;
    }
    for(auto p1 : points)
    {
        for(int d = 0; d < D; d++)
        {
            if(p1[d] < bbox_min[d])
                bbox_min[d] = p1[d];
            if(p1[d] > bbox_max[d])
                bbox_max[d] = p1[d];
        }
    }
    vps.clear();
    std::vector<Point> extra_pts;
    for(int it = 0; it < max_it; it++)
    {
        CGAL::Unique_hash_map<Vertex_const_handle, std::vector<double>> vertex_map;
        CGAL::Unique_hash_map<Vertex_const_handle, std::vector<double>> norm_map;
        CGAL::Unique_hash_map<Vertex_const_handle, std::vector<double>> scale_map;
        CGAL::Unique_hash_map<Vertex_const_handle, std::vector<double>> los_map;
        for(auto vv = traits_raw.vertices_begin(tri); vv != traits_raw.vertices_end(tri) ; ++vv)
        {
            if(tri.is_infinite(vv))
                continue;
            vertex_map[vv] = std::vector<double>(D+1,0);
            norm_map[vv] = std::vector<double>(D,0);
            scale_map[vv] = std::vector<double>(D,0);
            los_map[vv] = std::vector<double>(D,0);
        }
        for(int ii = 0; ii < points.size(); ii++)
        {
            Point pp = points[ii];
            double ss = 1;
            auto vv = tri.nearest_vertex(pp);
            for(int d = 0; d < D; d++)
            {
                vertex_map[vv][d] += ss*pp[d];
                norm_map[vv][d] += ss*norms[ii][D-1][d];
                scale_map[vv][d] += ss*scales[ii][d];
                if (v_los.size() > 0)
                    los_map[vv][d] += ss*v_los[ii][d];
            }
            vertex_map[vv][D]+=ss;
        }
        for(auto vv = traits_raw.vertices_begin(tri); vv != traits_raw.vertices_end(tri) ; ++vv)
        {
            if(tri.is_infinite(vv))
                continue;
            auto vpo = vv->point();
            auto vpn = vertex_map[vv];
            auto vn = norm_map[vv];
            auto vs = scale_map[vv];
            if(vpn[D] == 0)
                continue;
            for(int d = 0; d < D; d++)
            {
                vpn[d]=vpn[d]/vpn[D];
                vn[d]=vn[d]/vpn[D];
                vs[d]=vs[d]/vpn[D];
            }
            auto pp = Point(vpo[0] +(vpn[0]-vpo[0]),
                            vpo[1] +(vpn[1]-vpo[1]),
                            vpo[2] +(vpn[2]-vpo[2]));
            bool do_insert = true;
            for(int d = 0; d < D; d++)
                if(pp[d] > bbox_max[d] || pp[d] < bbox_min[d])
                    do_insert = false;
            if(do_insert)
            {
                tri.move(vv,pp);
            }
            double sp = 1;
            if(v_los.size() > 0)
            {
                auto vlos = los_map[vv];
                CGAL::Vector_3<typename Traits::K> v1{vn[0],vn[1],vn[2]};
                CGAL::Vector_3<typename Traits::K> v2{vlos[0],vlos[1],vlos[2]};
                sp = (normalize(v1)*normalize(v2));
            }
            double scll = std::max({vs[0], vs[1], vs[2]})/3;
            if(((double) rand() / (RAND_MAX)) > 0.8 && it == max_it-1 && sp > 0.7)
            {
                if(v_los.size() > 0)
                {
                    auto vlos = los_map[vv];
                    auto pp2 = Point(vpn[0]-vlos[0]*scll,vpn[1]-vlos[1]*scll,vpn[2]-vlos[2]*scll);
                    extra_pts.push_back(pp2);
                }
                else
                {
                    auto pp2 = Point(vpn[0]-vn[0]*scll,vpn[1]-vn[1]*scll,vpn[2]-vn[2]*scll);
                    extra_pts.push_back(pp2);
                }
            }
        }
        // if(it == max_it-1)
        // {
        //dump_tessel(tri,it,max_it,tid,bbox_min,bbox_max);
        //dump_vector_pts(extra_pts,it,tid);
        //}
    }
    for(auto vv = traits_raw.vertices_begin(tri); vv != traits_raw.vertices_end(tri) ; ++vv)
    {
        bool do_insert = true;
        for(int d = 0; d < D; d++)
            if(vv->point()[d] > bbox_max[d] || vv->point()[d] < bbox_min[d])
                do_insert = false;
        if(do_insert)
            vps.push_back(vv->point());
    }
    for(auto pp : extra_pts)
    {
        bool do_insert = true;
        for(int d = 0; d < D; d++)
            if(pp[d] > bbox_max[d] || pp[d] < bbox_min[d])
                do_insert = false;
        if(do_insert)
            vps.push_back(pp);
    }
    return 0;
}




int
wasure_algo::simplify(std::vector<Point> & points, std::vector<bool> & do_keep, double dist )
{
    int D = Traits::D;
    int nbp = points.size();
    double eps = 0;
    int K_T = DIM_SIZE_NB;
    if(K_T > points.size() -1)
        K_T = points.size() - 1;
    ANNpointArray	dataPts;
    ANNpoint	queryPt;
    ANNidxArray	nnIdx;
    ANNdistArray	dists;
    ANNkd_tree*	kdTree;
    queryPt = annAllocPt(D);
    dataPts = annAllocPts(nbp, D);
    nnIdx = new ANNidx[K_T];
    dists = new ANNdist[K_T];
    int npts =0;
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        for(int d = 0; d < D; d++)
        {
            dataPts[npts][d] = p1[d];
        }
        do_keep[npts] = false;
        npts++;
    }
    kdTree = new ANNkd_tree(dataPts,npts,D);
    int cur_id = 0;
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        for(int d = 0; d < D; d++)
            queryPt[d] = p1[d];
        kdTree->annkSearch(queryPt,K_T, nnIdx,dists,eps);
        for(int k = 1; k < K_T; k++)
        {
            double kdist = dists[k];
            if(kdist > dist)
            {
                do_keep[cur_id] = true;
                break;
            }
            if(do_keep[nnIdx[k]])
            {
                break;
            }
        }
        cur_id++;
    }
    delete kdTree;
    delete [] nnIdx;
    delete [] dists;
    return 0;
}



// =======================  DIM ================================
void
wasure_algo::compute_svd(int K_T, const  ANNidxArray & nnIdx, const std::vector<Point> & points,std::vector<Point> & pts_norms, std::vector<double> &coords_scale)
{
    int D = Traits::D;
    std::vector<Point> loc_pts;
    for(int i = 0; i < K_T; i++)
    {
        int idx = nnIdx[i];
        loc_pts.push_back(points[idx]);
    }
    int N = loc_pts.size();
    Eigen::MatrixXf mat(N,D);
    int acc = 0;
    for(std::vector<Point>::iterator pit = loc_pts.begin(); pit != loc_pts.end(); ++pit)
    {
        Point p = *pit;
        for(int d = 0 ; d < D ; d++)
            mat(acc,d) = p[d];
        acc++;
    }
    Eigen::MatrixXf m = mat.rowwise() - mat.colwise().mean();
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXf sv = svd.singularValues() ;
    Eigen::MatrixXf ev =  svd.matrixV();
    for(int d1 = 0; d1 < D; d1++)
    {
        double coords_norm[Traits::D];
        for(int d2 = 0; d2 < D; d2++)
        {
            coords_norm[d2] = ev(d2,d1);
        }
        pts_norms.push_back(traits.make_point(coords_norm));
    }
    double v_min = 0.0000001;
    for(int d = 0; d < D; d++)
    {
        coords_scale[d] = sv(d);
        if(coords_scale[d] <= v_min || std::isnan(coords_scale[d]))
            coords_scale[d] = v_min;
    }
}



int
wasure_algo::compute_dim(std::vector<Point> & points, std::vector<std::vector<Point> > & norms, std::vector<std::vector<double>> & scales,ddt::logging_stream & loggin)
{
    int D = Traits::D;
    auto is_nan = [](Point x)
    {
        int D = Traits::D;
        for(int d = 0; d < D; d++)
        {
            if(x[d] != x[d])
            {
                return true;
            }
        }
        return false;
    };
    auto new_end = std::remove_if(points.begin(), points.end(), is_nan);
    points.erase(new_end, points.end());
    int nbp = points.size();
    double eps = 0;
    int K_T = 150;
    int step = 10;
    if(K_T > points.size() -1)
    {
        K_T = points.size() -1;
    }
    if(K_T/step < 1)
        step = K_T;
    int K_T_min = 10;
    ANNpointArray	dataPts;
    ANNpoint	queryPt;
    ANNidxArray	nnIdx;
    ANNdistArray	dists;
    ANNkd_tree*	kdTree;
    queryPt = annAllocPt(D);
    dataPts = annAllocPts(nbp, D);
    nnIdx = new ANNidx[K_T];
    dists = new ANNdist[K_T];
    int npts =0;
    loggin.step("build_ann");
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        // last-minute quick n diry hack for open source release discoverd during stress testing the algo
        for(int d = 0; d < D; d++)
            dataPts[npts][d] = p1[d];
        npts++;
    }
    kdTree = new ANNkd_tree(dataPts,npts,D);
    loggin.step("compute_dim_on_pts");
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        std::vector<Point> pts_norms;
        std::vector<double> coords_scale(D);
        for(int d1 = 0; d1 < D; d1++)
        {
            double coords_norm[Traits::D];
            for(int d2 = 0; d2 < D; d2++)
            {
                if(d1 != d2)
                    coords_norm[d2] = 0;
                else
                    coords_norm[d2] = 1;
            }
            pts_norms.push_back(traits.make_point(coords_norm));
            coords_scale[d1] = 1;
        }
        double entropy = 1000000000;
        for(int d = 0; d < D; d++)
            queryPt[d] = p1[d];
        kdTree->annkSearch(queryPt,K_T, nnIdx,dists,eps);
        for(int k = K_T_min; k < K_T; k++)
        {
            if(dists[k] < 0.02 && k < 30)
                continue;
            std::vector<Point> cur_pts_norms;
            std::vector<double> cur_coords_scale(D);
            compute_svd(k, nnIdx, points,cur_pts_norms,cur_coords_scale);
            double cst = 0;
            std::vector<double> dim(D);
            for(int d = 0 ; d < D; d++)
            {
                cst += cur_coords_scale[d];
            }
            for(int d = 0; d < D; d++)
            {
                dim[d] = (d == D-1) ? cur_coords_scale[d]/cst : (cur_coords_scale[d] - cur_coords_scale[d+1])/cst;
            }
            double cur_entropy = 0;
            for(int d = 0; d < D; d++)
            {
                cur_entropy += -dim[d]*log(dim[d]);
            }
            if(cur_entropy < entropy)
            {
                entropy = cur_entropy;
                pts_norms = cur_pts_norms;
                coords_scale = cur_coords_scale;
            }
        }
        if(pts_norms.size() < D )
        {
            std::cerr << "Error: " << "dim:" << D << " kt:" << K_T << " entropy:" << entropy  << std::endl;
            std::vector<Point> loc_pts;
            for(int i = 0; i < K_T; i++)
            {
                int idx = nnIdx[i];
                loc_pts.push_back(points[idx]);
            }
            bool are_equal = true;
            for(int i = 1; i < loc_pts.size(); i++)
                if(loc_pts[i] != loc_pts[i-1])
                {
                    are_equal = false;
                    break;
                }
            if(are_equal)
            {
                std::cerr << "Point " << loc_pts[0] << " duplicated " << K_T << " time, svd not defined!" << std::endl;
                std::cerr << "remove duplicted point to continue" << std::endl;
                return 1;
            }
        }
        else
        {
            norms.push_back(pts_norms);
            scales.push_back(coords_scale);
        }
    }
    delete kdTree;
    delete [] nnIdx;
    delete [] dists;
    return 0;
}


int
wasure_algo::compute_dim_with_simp(std::vector<Point> & points, std::vector<std::vector<Point> > & norms, std::vector<std::vector<double>> & scales,std::vector<Point> & simp,double pscale)
{
    int D = Traits::D;
    std::vector<double> entropy_vect;
    int nbp = points.size();
    double eps = 0;
    int K_T = 150;
    if(K_T > points.size() -1)
        K_T = points.size() - 1;
    ANNpointArray	dataPts;
    ANNpoint	queryPt;
    ANNidxArray	nnIdx;
    ANNdistArray	dists;
    ANNkd_tree*	kdTree;
    queryPt = annAllocPt(D);
    dataPts = annAllocPts(nbp, D);
    nnIdx = new ANNidx[K_T];
    dists = new ANNdist[K_T];
    double bbox_min[Traits::D];
    double bbox_max[Traits::D];
    for(int d = 0; d < D; d++)
    {
        bbox_min[d] = 100000;
        bbox_max[d] = -100000;
    }
    int npts =0;
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        for(int d = 0; d < D; d++)
        {
            dataPts[npts][d] = p1[d];
            if(p1[d] < bbox_min[d])
                bbox_min[d] = p1[d];
            if(p1[d] > bbox_max[d])
                bbox_max[d] = p1[d];
        }
        npts++;
    }
    kdTree = new ANNkd_tree(dataPts,npts,D);
    for(std::vector<Point>::iterator pit = points.begin(); pit != points.end(); ++pit)
    {
        Point p1 = *pit;
        std::vector<Point> pts_norms;
        std::vector<double> coords_scale(D);
        double entropy = 1000000000;
        for(int d = 0; d < D; d++)
            queryPt[d] = p1[d];
        kdTree->annkSearch(queryPt,K_T, nnIdx,dists,eps);
        for(int k = 3; k < K_T; k++)
        {
            std::vector<Point> cur_pts_norms;
            std::vector<double> cur_coords_scale(D);
            compute_svd(k, nnIdx, points,cur_pts_norms,cur_coords_scale);
            double cst = 0;
            std::vector<double> dim(D);
            for(int d = 0 ; d < D; d++)
            {
                cst += cur_coords_scale[d];
            }
            for(int d = 0; d < D; d++)
            {
                dim[d] = (d == D-1) ? cur_coords_scale[d]/cst : (cur_coords_scale[d] - cur_coords_scale[d+1])/cst;
            }
            double cur_entropy = 0;
            for(int d = 0; d < D; d++)
            {
                cur_entropy += -dim[d]*log(dim[d]);
            }
            if(cur_entropy < entropy)
            {
                entropy = cur_entropy;
                pts_norms = cur_pts_norms;
                coords_scale = cur_coords_scale;
            }
        }
        entropy_vect.push_back(entropy);
        if(pts_norms.size() < D)
        {
            std::cerr << "Error :" << "dim:" << D << " kt:" << K_T << " entropy:" << entropy  << std::endl;
            std::vector<Point> loc_pts;
            for(int i = 0; i < K_T; i++)
            {
                int idx = nnIdx[i];
                loc_pts.push_back(points[idx]);
            }
            bool are_equal = true;
            for(int i = 1; i < loc_pts.size(); i++)
                if(loc_pts[i] != loc_pts[i-1])
                {
                    are_equal = false;
                    break;
                }
            if(are_equal)
            {
                std::cerr << "Point " << loc_pts[0] << " duplicated " << K_T << " time, svd not defined!" << std::endl;
                std::cerr << "remove duplicted point to continue" << std::endl;
                return 1;
            }
        }
        else
        {
            norms.push_back(pts_norms);
            scales.push_back(coords_scale);
        }
    }
    delete kdTree;
    delete [] nnIdx;
    delete [] dists;
    return 0;
}


void wasure_algo::flip_dim_ori( std::vector<Point> & points, std::vector<std::vector<Point> > & norms, std::vector<Point> &  ori)
{
    int nbp = points.size();
    int D = Traits::D;
    int nb_flip = 0;
    for(int n = 0; n < nbp; n++)
    {
        Point pi  = ori[n];
        Point p1 = points[n];
        std::vector<double> pts_coefs = compute_base_coef<Point>(p1,pi,norms[n],D);
        if(pts_coefs[D-1] < 0)
        {
            double coords_norm[Traits::D];
            for(int d = 0; d < D; d++)
            {
                coords_norm[d] = -norms[n][D-1][d];
            }
            norms[n][D-1] = traits.make_point(coords_norm);
            nb_flip++;
        }
        else
        {
        }
    }
}

void wasure_algo::flip_dim( std::vector<Point> & points, std::vector<std::vector<Point> > & norms, Point p1)
{
    int nbp = points.size();
    int D = Traits::D;
    for(int n = 0; n < nbp; n++)
    {
        Point pi  = points[n];
        std::vector<double> pts_coefs = compute_base_coef<Point>(p1,pi,norms[n],D);
        if(pts_coefs[D-1] < 0)
        {
            double coords_norm[Traits::D];
            for(int d = 0; d < D; d++)
            {
                coords_norm[d] = -norms[n][D-1][d];
            }
            norms[n][D-1] = traits.make_point(coords_norm);
        }
        else
        {
        }
    }
}





std::vector<double>  wasure_algo::Pick_2d(const Point & v0,const Point & v1,const Point & v2)
{
    int D=2;
    std::vector<double> coords(D);
    double r1 = (((double) rand() / (RAND_MAX)));
    double r2 = (((double) rand() / (RAND_MAX)));
    for(int d = 0; d < D; d++)
        coords[d] = (1 - sqrt(r1)) * v0[d] + (sqrt(r1) * (1 - r2)) * v1[d] + (sqrt(r1) * r2) * v2[d];
    return coords;
}

std::vector<double>  wasure_algo::Pick_3d(const Point & v0,const Point & v1,const Point & v2,const Point & v3)
{
    int D = Traits::D;
    double s = (((double) rand() / (RAND_MAX)));
    double t = (((double) rand() / (RAND_MAX)));
    double u = (((double) rand() / (RAND_MAX)));
    if(s+t>1.0)
    {
        s = 1.0 - s;
        t = 1.0 - t;
    }
    if(t+u>1.0)
    {
        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    }
    else if(s+t+u>1.0)
    {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }
    double a=1-s-t-u; // a,s,t,u are the barycentric coordinates of the random point.
    std::vector<double> coords(D);
    coords[0] = (v0[0]*a + v1[0]*s + v2[0]*t + v3[0]*u);
    coords[1] = (v0[1]*a + v1[1]*s + v2[1]*t + v3[1]*u);
    coords[2] = (v0[2]*a + v1[2]*s + v2[2]*t + v3[2]*u);
    return coords;
}




void wasure_algo::init_sample(DT & tri,int nb_samples, int dim)
{
    // void
}


void wasure_algo::finalize_sample(DT & tri, int nb_samples)
{
    // void
}

double get_min_scale(std::vector<double> &v)
{
    size_t n = v.size() / 10;
    return v[n];
}

double get_median_scale(std::vector<double> &v)
{
    size_t n = v.size()*3 / 4;
    return v[n];
}





void wasure_algo::ds_score(double v_e1,double v_o1,double v_u1,double v_e2,double v_o2,double v_u2,double & v_e3,double & v_o3,double & v_u3)
{
    double vK = v_o1*v_e2 + v_e1*v_o2;
    v_e3 = (v_e1*v_e2 + v_e1*v_u2 + v_u1*v_e2)/(1-vK);
    v_o3 = (v_o1*v_o2 + v_o1*v_u2 + v_u1*v_o2)/(1-vK);
    v_u3 = (v_u1*v_u2)/(1-vK);
}






void wasure_algo::compute_dst_mass_norm(std::vector<double> coefs, std::vector<double> scales, double coef_conf, double pdfs_e,double pdfs_o, double & v_e1, double & v_o1, double & v_u1)
{
    int D = scales.size();
    double c3 = coefs[D-1];
    double nscale = scales[D-1];
    if(nscale <= 0) nscale = 0.000001;
    if(c3 > 0)
    {
        v_e1 =  1-0.5*(exp(-fabs(c3)/nscale));
        v_o1 = 0.5*(exp(-fabs(c3)/nscale));
    }
    else if (c3 < 0)
    {
        v_e1 = 0.5*(exp(-fabs(c3)/nscale));
        v_o1 = 1-0.5*(exp(-fabs(c3)/nscale));
    }
    else
    {
        v_e1 = v_o1 = 0.5;
    }
    for(int d = 0; d < D-1; d++)
    {
        if(scales[d] <= 0)
        {
            v_e1 = v_o1 = 0;
        }
        else
        {
            v_e1 = v_e1*score_pdf(coefs[d],scales[d]);
            v_o1 = v_o1*score_pdf(coefs[d],scales[d]);
        }
    }
    v_e1 = v_e1*exp(-(fabs(c3)/(pdfs_e))*(fabs(c3)/(pdfs_e)));
    v_o1 = v_o1*exp(-(fabs(c3)/(pdfs_o))*(fabs(c3)/(pdfs_o)));
    v_e1 = v_e1*coef_conf;
    v_o1 = v_o1*coef_conf;
    v_u1 = 1-v_e1-v_o1;
    regularize(v_e1,v_o1,v_u1);
}





double
get_conf_volume(const std::vector<double> & pts_scales, int dim)
{
    double c1 = MIN(pts_scales[dim-1]/pts_scales[0],1);
    return 1-c1;
}

void
wasure_algo::get_params_surface_dst(const std::vector<double> & pts_scales,double glob_scale,double min_scale, double & pdf_smooth,double & coef_conf, int D)
{
    double data_scale = *std::max_element(pts_scales.begin(),pts_scales.end());
    if(glob_scale > 0)
    {
        pdf_smooth = glob_scale*3;
    }
    else
    {
        pdf_smooth = data_scale*GLOB_SMOOTH;
    }
    coef_conf = 1;
    if(min_scale > 0)
        coef_conf = coef_conf*MIN(MAX(min_scale/data_scale,0.000001),1);
}


void
wasure_algo::get_params_conflict_dst(const std::vector<double> & pts_scales,double glob_scale,double min_scale, double & pdf_smooth,double & coef_conf, int D)
{
    double data_scale = *std::max_element(pts_scales.begin(),pts_scales.end());
    if(glob_scale > 0)
    {
        pdf_smooth = glob_scale/3;
    }
    else
    {
        pdf_smooth = data_scale/3;
    }
    coef_conf = 1;
}




void
wasure_algo::compute_dst_tri(DTW & tri, wasure_data<Traits>  & datas_tri, wasure_data<Traits>  & datas_pts, wasure_params & params)
{
    std::vector<Point> & points_dst =  datas_pts.format_points;
    std::vector<Point> & v_los =  datas_pts.format_normals;
    std::vector<std::vector<Point>> & norms = datas_pts.format_egv;
    std::vector<std::vector<double>> & scales = datas_pts.format_sigs;
    std::vector<int> & format_flags = datas_pts.format_flags;
    std::vector<std::vector<double>> & v_dst = datas_tri.format_dst;
    int D = datas_pts.D;
    int nbp = points_dst.size();
    int K_T = 30;
    if(K_T >= nbp)
        K_T = nbp-1;
    double eps = 0;
    ANNpointArray	dataPts;
    ANNpoint	queryPt;
    ANNidxArray	nnIdx;
    ANNdistArray	dists;
    ANNkd_tree*	kdTree;
    queryPt = annAllocPt(D);
    dataPts = annAllocPts(nbp, D);
    nnIdx = new ANNidx[K_T];
    dists = new ANNdist[K_T];
    int npts =0;
    for(std::vector<Point>::iterator pit = points_dst.begin(); pit != points_dst.end(); ++pit)
    {
        Point p1 = *pit;
        for(int d = 0; d < D; d++)
            dataPts[npts][d] = p1[d];
        npts++;
    }
    kdTree = new ANNkd_tree(dataPts,npts,D);
    std::vector<double> v_scale(scales.size());
    for(int i = 0; i < scales.size() ; i++)
    {
        v_scale[i] = scales[i][D-1];
    }
    std::sort(v_scale.begin(), v_scale.end());
    if(params.dst_scale < 0)
    {
        params.min_scale = -1;
    }
    else
    {
        params.min_scale = params.dst_scale;
    }
    bool do_debug = false;
    for( auto cit = tri.cells_begin();
            cit != tri.cells_end(); ++cit )
    {
        if( cit->is_infinite() )
            continue;
        int cid = cit->lid();
        double  vpe = v_dst[cid][0];
        double  vpo = v_dst[cid][1];
        double  vpu = v_dst[cid][2];
        double pe_acc=0;
        double po_acc=0;
        double pu_acc=0;
        for(int x = 0; x < params.nb_samples; x++)
        {
            std::vector<double>  C =  (x == 0) ? cit->barycenter() : Pick(cit->full_cell(),D);
            Point  PtSample = traits.make_point(C.begin());
            for(int d = 0; d < D; d++)
            {
                queryPt[d] = PtSample[d];
            }
            std::ofstream myfile;
            kdTree->annkSearch(queryPt,K_T, nnIdx,dists,eps);
            vpe = vpo = 0;
            vpu = 1;
            for(int k = 0; k < K_T; k++)
            {
                double pe1,po1,pu1,pe2,po2,pu2;
                int idx = nnIdx[k];
                Point Pt3d = points_dst[idx];
                std::vector<Point> pts_norms = norms[idx];
                std::vector<double> pts_scales = scales[idx];
                if(((int)format_flags[idx]) < 0)
                    continue;
                double pdf_smooth = -1;
                double coef_conf = -1;
                double gbl_scale = -1;
                get_params_surface_dst(pts_scales,gbl_scale,params.min_scale,pdf_smooth,coef_conf,D);
                std::vector<double> pts_coefs = compute_base_coef(Pt3d,PtSample,pts_norms,D);
                if(params.dst_scale > 0)
                {
                    if(((int)format_flags[idx]) > 0)
                    {
                        coef_conf = 0.2*coef_conf;
                    }
                }
                double sp = 1;
                if(v_los.size() > 0)
                {
                    Point & pts_los = v_los[idx];
                    CGAL::Vector_3<typename Traits::K> v1{pts_los[0],pts_los[1],pts_los[2]};
                    CGAL::Vector_3<typename Traits::K> v2{pts_norms[D-1][0],pts_norms[D-1][1],pts_norms[D-1][2]};
                    //sp = (normalize(v1)*normalize(v2));
                }
                coef_conf = coef_conf*sp;
                if(coef_conf != coef_conf)
                    coef_conf = 0;
                compute_dst_mass_norm(pts_coefs,pts_scales,coef_conf,pdf_smooth,pdf_smooth,pe2, po2,pu2);
                pe1=vpe;
                po1=vpo;
                pu1=vpu;
                ds_score(pe1,po1,pu1,pe2,po2,pu2,pe1,po1,pu1);
                regularize(pe1,po1,pu1);
                vpe = pe1;
                vpo = po1;
                vpu = pu1;
            }
            pe_acc += vpe;
            po_acc += vpo;
            pu_acc += vpu;
            if(do_debug)
                myfile.close();
        }
        // NAN check
        if(vpe == vpe &&
                vpo == vpo &&
                vpu == vpu)
        {
            v_dst[cid][0] = pe_acc/params.nb_samples;
            v_dst[cid][1] = po_acc/params.nb_samples;
            v_dst[cid][2] = pu_acc/params.nb_samples ;
        }
        else
        {
            std::cerr << "NAN ERROR DST" << std::endl;
        }
    }
}


double
wasure_algo::compute_angle_rad(const Point & a, const Point & b, const Point & c, int dim)
{
    std::vector<double> v1(dim);
    std::vector<double> v2(dim);
    for(int d = 0; d < dim; d++)
    {
        v1[d] = c[d]-a[d];
        v2[d] = c[d]-b[d];
    }
    double num=0,accv1=0,accv2=0;
    for(int d = 0; d < dim; d++)
    {
        num += v1[d]*v2[d];
        accv1 += v1[d]*v1[d];
        accv2 += v2[d]*v2[d];
    }
    double cosine = num/sqrt(accv1)/sqrt(accv2);
    return std::acos(cosine);
}


void
wasure_algo::compute_dst_mass_beam(std::vector<double> & coefs, std::vector<double> & scales,double angle,double angle_scale, double coef_conf, double & v_e1, double & v_o1, double & v_u1)
{
    int D = scales.size();
    double c3 = coefs[D-1];
    double nscale = scales[D-1];
    if(nscale <= 0) nscale = 0.000001;
    if(c3 >= 0)
    {
        v_e1 =  1-(exp(-fabs(c3)/(nscale*2)));
        v_o1 = 0;
    }
    else if (c3 < 0)
    {
        v_e1 = 0;
        v_o1 = 0;
    }
    v_e1 = v_e1*coef_conf;
    v_o1 = v_o1*coef_conf;
    regularize(v_e1,v_o1,v_u1);
}



void
wasure_algo::sample_cell(Cell_handle & ch,Point &  Pt3d, Point & PtCenter, wasure_data<Traits>  & datas_tri, wasure_data<Traits>  & datas_pts, wasure_params & params, int rid, int dim)
{
    Id cid = traits.gid(ch);
    std::vector<Point> & pts_norm = datas_pts.format_egv[rid];
    std::vector<double> & pts_scale = datas_pts.format_sigs[rid];
    std::vector<std::vector<double>> & v_dst = datas_tri.format_dst;
    Traits  traits;
    for(int x = 0; x < params.nb_samples; x++)
    {
        std::vector<double>  C = (x == 0) ? traits.get_cell_barycenter(ch) : Pick(ch,D);
        Point  PtSample = traits.make_point(C.begin());
        double pe1,po1,pu1,pe2,po2,pu2;
        pe1 = v_dst[cid][0];
        po1 = v_dst[cid][1];
        pu1 = v_dst[cid][2];
        std::vector<double> pts_coefs = compute_base_coef<Point>(Pt3d,PtSample,pts_norm,dim);
        double angle = compute_angle_rad(PtSample,Pt3d,PtCenter,dim);
        double pdf_smooth = -1;
        double coef_conf = -1;
        double gbl_scale = -1;
        get_params_surface_dst(pts_scale,gbl_scale,params.min_scale,pdf_smooth,coef_conf,dim);
        compute_dst_mass_beam(pts_coefs,pts_scale,angle,ANGLE_SCALE,coef_conf,pe2,po2,pu2);
        ds_score(pe1,po1,pu1,pe2,po2,pu2,pe1,po1,pu1);
        regularize(pe1,po1,pu1);
        if(pe1 == pe1 &&
                po1 == po1 &&
                pu1 == pu1)
        {
            v_dst[cid][0] = pe1;
            v_dst[cid][1] = po1 ;
            v_dst[cid][2] = pu1 ;
        }
    }
}




void
wasure_algo::sample_cell_raw(Cell_handle & ch,Point &  Pt3d, Point & PtCenter, wasure_data<Traits>  & datas_tri, wasure_data<Traits>  & datas_pts, wasure_params & params, int rid, int dim)
{
    Id cid = traits.gid(ch);
    Traits  traits;
    std::vector<std::vector<double>> & v_dst = datas_tri.format_dst;
    std::vector<Point> & pts_norm = datas_pts.format_egv[rid];
    std::vector<double> & pts_scale = datas_pts.format_sigs[rid];
    for(int x = 0; x < params.nb_samples; x++)
    {
        std::vector<double>  C = (x == 0) ? traits.get_cell_barycenter(ch) : Pick(ch,D);
        Point  PtSample = traits.make_point(C.begin());
        double pe1,po1,pu1,pe2,po2,pu2;
        std::vector<double> pts_coefs = compute_base_coef<Point>(Pt3d,PtSample,pts_norm,dim);
        pe1 = v_dst[cid][0];
        po1 = v_dst[cid][1];
        pu1 = v_dst[cid][2];
        double angle=0;
        double pdf_smooth = -1;
        double coef_conf = -1;
        double gbl_scale = -1;
        get_params_surface_dst(pts_scale,gbl_scale,params.min_scale,pdf_smooth,coef_conf,dim);
        compute_dst_mass_beam(pts_coefs,pts_scale,angle,ANGLE_SCALE,coef_conf,pe2,po2,pu2);
        pe1 = pe1*params.ray_weight;
        po1 = 0;
        pu1 = 1-pe1;
        ds_score(pe1,po1,pu1,pe2,po2,pu2,pe1,po1,pu1);
        regularize(pe1,po1,pu1);
        if(pe1 == pe1 &&
                po1 == po1 &&
                pu1 == pu1)
        {
            v_dst[cid][0] = pe1;
            v_dst[cid][1] = po1 ;
            v_dst[cid][2] = pu1 ;
        }
    }
}




Cell_handle wasure_algo::walk_locate(DT & tri,
                                     Point & Pt3d,  Point & Ptcenter, Point & Ptcenter_mir,
                                     wasure_data<Traits>  & datas_tri,
                                     wasure_data<Traits>  & datas_pts,
                                     wasure_params & params,
                                     int idr,
                                     Cell_handle & start_cell
                                    )
{
#if defined(DDT_CGAL_TRAITS_D)
    typedef typename DT::Face Face;
    DT::Locate_type  loc_type;
    Face face(tri.maximal_dimension());
    DT::Facet  facet;
    Cell_handle start = tri.locate(Ptcenter,start_cell);
    start_cell = start;
    DT::Geom_traits geom_traits;
    CGAL::Random                      rng_;
    Traits::K::Orientation_d orientation_pred = geom_traits.orientation_d_object();
    int cur_dim = tri.current_dimension();
    std::vector<CGAL::Oriented_side>  orientations_(cur_dim+1);
    if( cur_dim == -1 )
    {
        loc_type = DT::OUTSIDE_AFFINE_HULL;
        return Cell_handle();
    }
    else if( cur_dim == 0 )
    {
        Vertex_handle vit = tri.infinite_full_cell()->neighbor(0)->vertex(0);
        if( CGAL::EQUAL != geom_traits.compare_lexicographically_d_object()(Ptcenter_mir, vit->point()) )
        {
            loc_type = DT::OUTSIDE_AFFINE_HULL;
            return Cell_handle();
        }
        else
        {
            loc_type = DT::ON_VERTEX;
            face.set_full_cell(vit->full_cell());
            face.set_index(0, 0);
            return vit->full_cell();
        }
    }
    Cell_handle s;
    // if we don't know where to start, we start from any bounded full_cell
    if( Cell_handle() == start )
    {
        // THE HACK THAT NOBODY SHOULD DO... BUT DIFFICULT TO WORK AROUND
        // THIS... TODO: WORK AROUND IT
        Cell_handle inf_c = tri.infinite_full_cell();
        int inf_v_index = inf_c->index(tri.infinite_vertex());
        s = inf_c->neighbor(inf_v_index);
    }
    else
    {
        s = start;
        if( tri.is_infinite(s) )
        {
            int inf_v_index = s->index(tri.infinite_vertex());
            s = s->neighbor(inf_v_index);
        }
    }
    Cell_handle previous = Cell_handle();
    bool full_cell_not_found = true;
    while(full_cell_not_found) // we walk until we locate the query point |p|
    {
#ifdef CGAL_TRIANGULATION_STATISTICS
        ++walk_size_;
#endif
        // For the remembering stochastic walk, we need to start trying
        // with a random index:
        int j, i = rng_.get_int(0, cur_dim);
        // we check |p| against all the full_cell's hyperplanes in turn
        for(j = 0; j <= cur_dim; ++j, i = (i + 1) % (cur_dim + 1) )
        {
            Cell_handle next = s->neighbor(i);
            if( previous == next )
            {
                // no need to compute the orientation, we already know it
                orientations_[i] = CGAL::POSITIVE;
                continue; // go to next full_cell's facet
            }
            CGAL::Substitute_point_in_vertex_iterator<
            DT::Full_cell::Vertex_handle_const_iterator>
            spivi(s->vertex(i), &Ptcenter_mir);
            orientations_[i] = orientation_pred(boost::make_transform_iterator(s->vertices_begin(), spivi), boost::make_transform_iterator(s->vertices_begin() + cur_dim + 1, spivi));
            if( orientations_[i] != CGAL::NEGATIVE )
            {
                continue;
            }
            // At this point, we know that we have to jump to the |next|
            // full_cell because orientation_[i] == NEGATIVE
            previous = s;
            s = next;
            if( tri.is_infinite(next) )
            {
                // we have arrived OUTSIDE the convex hull of the triangulation,
                // so we stop the search
                full_cell_not_found = false;
                loc_type = DT::OUTSIDE_CONVEX_HULL;
                face.set_full_cell(s);
            }
            else
            {
                sample_cell_raw(s,Pt3d,Ptcenter,datas_tri,datas_pts,params,idr, D);
                for(int ii = 0; ii <= cur_dim; ii++)
                {
                    Cell_handle s_nbr = s->neighbor(ii);
                    if(!tri.is_infinite(s_nbr))
                        sample_cell(s_nbr,Pt3d,Ptcenter,datas_tri,datas_pts,params,idr, D);
                }
            }
            break;
        } // end of the 'for' loop
        if( ( cur_dim + 1 ) == j ) // we found the full_cell containing |p|
            full_cell_not_found = false;
    }
    // Here, we know in which full_cell |p| is in.
    // We now check more precisely where |p| landed:
    // vertex, facet, face or full_cell.
    if( ! tri.is_infinite(s) )
    {
        face.set_full_cell(s);
        int num(0);
        int verts(0);
        for(int i = 0; i < cur_dim; ++i)
        {
            if( orientations_[i] == CGAL::COPLANAR )
            {
                ++num;
                facet = DT::Facet(s, i);
            }
            else
                face.set_index(verts++, i);
        }
        //-- We could put the if{}else{} below in the loop above, but then we would
        // need to test if (verts < cur_dim) many times... we do it only once
        // here:
        if( orientations_[cur_dim] == CGAL::COPLANAR )
        {
            ++num;
            facet = DT::Facet(s, cur_dim);
        }
        else if( verts < cur_dim )
            face.set_index(verts, cur_dim);
        //-- end of remark above //
        if( 0 == num )
        {
            loc_type = DT::IN_FULL_CELL;
            face.clear();
        }
        else if( cur_dim == num )
            loc_type = DT::ON_VERTEX;
        else if( 1 == num )
            loc_type = DT::IN_FACET;
        else
            loc_type = DT::IN_FACE;
    }
    return s;
#endif
#if defined(DDT_CGAL_TRAITS_3)
    int n_of_turns = 10000;
    auto seg =  Traits::K::Segment_3(Ptcenter, Pt3d);
    if(tri.dimension() < 3)
        return start_cell;
    Cell_handle start = tri.locate(Ptcenter,start_cell);
    start_cell = start;
    // Make sure we continue from here with a finite cell.
    if(start == Cell_handle())
        start = tri.infinite_cell();
    int ind_inf;
    if(start->has_vertex(tri.infinite_vertex(), ind_inf))
        start = start->neighbor(ind_inf);
    // CGAL_triangulation_precondition(start != Cell_handle());
    //     CGAL_triangulation_precondition(! start->has_vertex(infinite));
    // We implement the remembering visibility walk.
    // In this phase, no need to be stochastic
    // Remembers the previous cell to avoid useless orientation tests.
    Cell_handle previous = Cell_handle();
    Cell_handle c = start;
    // Now treat the cell c.
try_next_cell:
    n_of_turns--;
    // We know that the 4 vertices of c are positively oriented.
    // So, in order to test if p is seen outside from one of c's facets,
    // we just replace the corresponding point by p in the orientation
    // test.  We do this using the array below.
    const Point* pts[4] = { &(c->vertex(0)->point()),
                            &(c->vertex(1)->point()),
                            &(c->vertex(2)->point()),
                            &(c->vertex(3)->point())
                          };
    // (non-stochastic) visibility walk
    Traits::K::Tetrahedron_3 tet(*pts[0],*pts[1],*pts[2],*pts[3]);
    if(CGAL::do_intersect(seg, tet))
    {
        sample_cell_raw(c,Pt3d,Ptcenter,datas_tri,datas_pts,params,idr, D);
    }
    for(int i=0; i != 4; ++i)
    {
        Cell_handle next = c->neighbor(i);
        if(previous == next) continue;
        // We temporarily put p at i's place in pts.
        const Point* backup = pts[i];
        pts[i] = &Ptcenter_mir;
        if(tri.inexact_orientation(*pts[0], *pts[1], *pts[2], *pts[3]) != CGAL::NEGATIVE)
        {
            pts[i] = backup;
            continue;
        }
        if(next->has_vertex(tri.infinite_vertex()))
        {
            // We are outside the convex hull.
            return start_cell;
        }
        previous = c;
        c = next;
        if(n_of_turns) goto try_next_cell;
    }
#endif
    return start_cell;
}


void wasure_algo::center_dst(DTW & ttri, wasure_data<Traits>  & datas_tri,std::vector<Point> & center_pts, Id tid)
{
    auto & tri = ttri.get_tile(tid)->tri();
    auto  tile = ttri.get_tile(tid);
    std::vector<std::vector<double>> & v_dst = datas_tri.format_dst;
    Cell_handle cloc_start =  Cell_handle();
    for(auto pit : center_pts)
    {
        Cell_handle cloc = tri.locate(pit,cloc_start);
        cloc_start = cloc;
        // Quick and dity hack
        int cid = tile->lid(cloc);
        double pe1,po1,pu1,pe2,po2,pu2;
        pe1 = v_dst[cid][0];
        po1 = v_dst[cid][1];
        pu1 = v_dst[cid][2];
        pe2 = 0.05;
        po2 = 0;
        pu2 = 0.95;
        ds_score(pe1,po1,pu1,pe2,po2,pu2,pe1,po1,pu1);
        regularize(pe1,po1,pu1);
        v_dst[cid][0] = pe1;
        v_dst[cid][1] = po1;
        v_dst[cid][2] = pu1;
    }
}



void wasure_algo::compute_dst_ray(DT & tri, wasure_data<Traits>  & datas_tri,wasure_data<Traits>  & datas_pts, wasure_params & params)
{
    Traits traits;
    std::vector<Point> & points = datas_pts.format_points;
    std::vector<Point> & centers = datas_pts.format_centers;
    std::vector<int> & format_flags = datas_pts.format_flags;
    double rat_ray_sample = params.rat_ray_sample;
    if(rat_ray_sample == 0)
        return;
    int D = Traits::D;
    int nb_pts = points.size();
    int acc = 0;
    Cell_handle  start_walk = Cell_handle();
    for(int n = 1 ; n < nb_pts; n++)
    {
        if(acc++ % ((int)(1.0/(rat_ray_sample)))  == 0)
        {
            Point & Pt3d = points[n];
            Point & Ptcenter = centers[n];
            double coords[Traits::D];
            for(int d = 0; d < D; d++)
                coords[d] = Pt3d[d] -(Ptcenter[d] - Pt3d[d]);
            Point Ptcenter_mir = traits.make_point(coords);
            if(((int)format_flags[n]) ==0)
            {
                walk_locate(tri,Pt3d,Ptcenter,Pt3d,datas_tri,datas_pts,params,n,start_walk);
            }
            else if(((int)format_flags[n]) ==1)
            {
                walk_locate(tri,Pt3d,Ptcenter,Ptcenter_mir,datas_tri,datas_pts,params,n,start_walk);
            }
        }
    }
}


void wasure_algo::compute_dst_with_center(DTW & tri, wasure_data<Traits>  & datas_tri, wasure_data<Traits>  & datas_pts, wasure_params & params, Id tid)
{
    compute_dst_tri(tri,datas_tri,datas_pts,params);
    DT & tri_tile  = tri.get_tile(tid)->triangulation();
    compute_dst_ray(tri_tile,datas_tri,datas_pts,params);
    center_dst(tri,datas_tri,datas_pts.format_centers,tid);
}



void wasure_algo::extract_surface(DTW & tri, std::map<Id,wasure_data<Traits> >  & w_datas_tri)
{
    int dim = Traits::D;
    std::vector<Point>  format_points;
    std::vector<int> v_simplex;
    std::map<Vertex_const_iterator, uint> vertex_map;
    int acc = 0;
    for(auto fit = tri.facets_begin(); fit != tri.facets_end(); fit++)
    {
        Cell_const_iterator fch = fit.full_cell();
        int id_cov = fit.index_of_covertex();
        Cell_const_iterator fchn = fch->neighbor(id_cov);
        int cccid = fch->cell_data().id;
        int cccidn = fchn->cell_data().id;
        if( w_datas_tri.find(fch->main_id()) == w_datas_tri.end() ||
                w_datas_tri.find(fchn->main_id()) == w_datas_tri.end()
          )
        {
            continue;
        }
        int ch1lab = w_datas_tri[fch->main_id()].format_labs[cccid];
        int chnlab = w_datas_tri[fchn->main_id()].format_labs[cccidn];
        if(ch1lab != chnlab)
        {
            for(int i = 0; i < dim+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    if(vertex_map.find(v) == vertex_map.end() && v->is_main() )
                    {
                        vertex_map[v] = acc++;
                        format_points.push_back(v->point());
                    }
                }
            }
        }
    }
    for(auto fit = tri.facets_begin(); fit != tri.facets_end(); fit++)
    {
        Cell_const_iterator fch = fit.full_cell();
        int id_cov = fit.index_of_covertex();
        Cell_const_iterator fchn = fch->neighbor(id_cov);
        int cccid = fch->cell_data().id;
        int cccidn = fchn->cell_data().id;
        if( w_datas_tri.find(fch->main_id()) == w_datas_tri.end() ||
                w_datas_tri.find(fchn->main_id()) == w_datas_tri.end()
          )
        {
            continue;
        }
        bool do_skip = false;
        int ch1lab = w_datas_tri[fch->main_id()].format_labs[cccid];
        int chnlab = w_datas_tri[fchn->main_id()].format_labs[cccidn];
        if(ch1lab != chnlab && !do_skip)
        {
            for(int i = 0; i < dim+1; ++i)
            {
                if(i != id_cov)
                {
                    Vertex_const_iterator v = fch->vertex(i);
                    v_simplex.push_back(vertex_map[v]);
                }
            }
        }
    }
    ddt_data<Traits> datas_out;
    datas_out.dmap[datas_out.xyz_name] = ddt_data<Traits>::Data_ply(datas_out.xyz_name,"vertex",dim,dim,DATA_FLOAT_TYPE);
    datas_out.dmap[datas_out.simplex_name] = ddt_data<Traits>::Data_ply(datas_out.simplex_name,"face",dim,dim,tinyply::Type::INT32);
    datas_out.dmap[datas_out.xyz_name].fill_full_uint8_vect(format_points);
    datas_out.dmap[datas_out.simplex_name].fill_full_uint8_vect(v_simplex);
}

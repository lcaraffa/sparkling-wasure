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
#ifndef TBMRFRECO_H
#define TBMRFRECO_H


#include "tbmrf.hpp"

template<typename DTW,typename D_MAP>
class tbmrf_reco : public tbmrf<DTW,D_MAP>
{
public :

    typedef typename DTW::Cell_const_iterator                 Cell_const_iterator;
    typedef typename DTW::Traits                              Traits;
    typedef typename Traits::Point                            Point;
    typedef typename Traits::Vertex_const_handle              Vertex_const_handle;
    typedef typename Traits::Cell_const_handle                Cell_const_handle;
    typedef typename Traits::Cell_handle                      Cell_handle;
    typedef typename Traits::Vertex_handle                    Vertex_handle;


    tbmrf_reco(int nblabs, DTW * t,D_MAP * dm) : mode(-1),pLabsIn(nblabs),pLabsOut(nblabs),pLabsUnk(nblabs),tbmrf<DTW,D_MAP>(nblabs,t,dm)
    {
        for(int i = 0; i < nblabs; i++)
        {
            double prob = ((double)i)/(double(nblabs-1));
            tbmrf<DTW,D_MAP>::labs[i] = i;
            tbmrf<DTW,D_MAP>::pLabs[i] = prob;
            pLabsIn[i] = (prob >= 0.5) ? (prob-0.5)*2 : 0;
            pLabsOut[i] = (prob < 0.5) ? (0.5-prob)*2 : 0;
            pLabsUnk[i] = 1-fabs(prob-0.5)*2;
        }
    }

    void set_mode(int mm)
    {
        mode = mm;
    }



    double get_volume_reco(Cell_const_iterator fch)
    {
        return log(1+tbmrf<DTW,D_MAP>::get_volume(fch));
    }

    double get_score_linear(Cell_const_iterator fch,int label,D_MAP & data_map)
    {
        double volume = 1;
        int D = Traits::D;
        try
        {
            if(fch->is_infinite())
            {
                if(mode != 1)
                    return 0;
                for(int d = 0; d < D+1; d++)
                {
                    if(fch->neighbor(d)->is_mixed())
                        std::cerr << "mixed detected" << std::endl;
                }
                //int id_cov;
                auto chn = fch;//->full_cell();
                //int mir = fch->mirror_index(id_cov);
                int ii = 0;
                const Point& a = chn->vertex((ii+1)&3)->point();
                const Point& b = chn->vertex((ii+2)&3)->point();
                const Point& c = chn->vertex((ii+3)&3)->point();
                const Point& d = chn->vertex(ii)->point();
                CGAL::Vector_3<typename Traits::K> v1 = a - c;
                CGAL::Vector_3<typename Traits::K> v2 = b - c;
                auto vn = CGAL::cross_product(v1,v2);
                double coef_proj = compute_coef_proj(c,d,vn,D);
                auto bary = chn->barycenter();
                if((coef_proj > 0 && vn[D-1] < 0) ||
                        (coef_proj < 0 && vn[D-1] > 0) )
                {
                    if(pLabsOut[label] > 0.5)
                        return 100000000;
                    else
                        return 0;
                }
                else
                {
                    if(pLabsOut[label] > 0.5)
                        return 0;
                    else
                        return 100000000;
                }
            }
        }
        catch (ddt::DDT_exeption& e)
        {
            std::cerr << "!! WARNING !!!" << std::endl;
            std::cerr << "Exception catched : " << e.what() << std::endl;
            return 1000000;
        }
        double nbe = 1;
        double coef = volume/nbe;
        int cell_id = fch->lid();
        int tile_id = fch->tile()->id();
        double PIn = data_map[tile_id].format_dst[cell_id][0];
        double POut = data_map[tile_id].format_dst[cell_id][1];
        double PUnk = data_map[tile_id].format_dst[cell_id][2];
        double scoreCurr = fabs(pLabsIn[label] - PIn) + fabs(pLabsOut[label] - POut) + fabs(pLabsUnk[label] - PUnk);
        return coef*scoreCurr;
    }


    int extract_surface(int tid, std::vector<Facet_const_iterator> & lft, D_MAP & w_datas_tri)
    {
        for(auto fit = this->tri->facets_begin();  fit != this->tri->facets_end(); ++fit)
        {
            try
            {
                if(fit->main_id() != tid || fit->is_infinite())
                    continue;
                Cell_const_iterator tmp_fch = fit.full_cell();
                int tmp_idx = fit.index_of_covertex();
                Cell_const_iterator tmp_fchn = tmp_fch->neighbor(tmp_idx);
                if(!this->tri->tile_is_loaded(tmp_fch->main_id()) ||
                        !this->tri->tile_is_loaded(tmp_fchn->main_id()))
                    continue;
                bool is_on_convex = false;
                if(tmp_fch->is_infinite() ||  tmp_fchn->is_infinite() )
                    is_on_convex = true;
                Cell_const_iterator fch = tmp_fch->main();
                Cell_const_iterator fchn = tmp_fchn->main();
                int cccid = fch->lid();
                int cccidn = fchn->lid();
                int ch1lab = w_datas_tri[fch->tile()->id()].format_labs[cccid];
                int chnlab = w_datas_tri[fchn->tile()->id()].format_labs[cccidn];
                if(
                    (ch1lab != chnlab) ||  (mode == 0 && (is_on_convex && (ch1lab == 0 || chnlab == 0))))
                {
                    lft.push_back(*fit);
                }
            }
            catch (ddt::DDT_exeption& e)
            {
                std::cerr << "!! WARNING !!!" << std::endl;
                std::cerr << "Exception catched : " << e.what() << std::endl;
                continue;
            }
        }
        return 0;
    }

    int mode;
    std::vector<double>  pLabsIn;
    std::vector<double>  pLabsOut;
    std::vector<double> pLabsUnk;

};

#endif

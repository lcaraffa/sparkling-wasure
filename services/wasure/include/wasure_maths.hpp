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
#ifndef WASURE_MATHS_HPP
#define WASURE_MATHS_HPP

#include "wasure_typedefs.hpp"
#include <Eigen/SVD>
#include <Eigen/Core>

inline double mEstimator(double a,double t)
{
    return ((1.0/a)*(pow((1+t),a) - 1.0));
};

// First order derivative of the Mestimator
//  y = (1+t).^(a-1);
// Input  a : power
//        t : Squared residual on data
// Output   : 1er order Error
inline double dmEstimator(double a,double t)
{
    return pow((1+t),(a-1.0));
};


inline int factorial(int x)
{
    return (x == 1 ? x : x * factorial(x - 1));
}


inline double  m_pdf(double x,double a)
{
    return exp(-0.5*mEstimator(x,a));
}


double score_pdf(double a,double scale);


template<typename T_PTS>
double
n_volume(std::list<T_PTS> & lp, int D)
{
    T_PTS p0 = lp.front();
    int nbp = lp.size();
    lp.pop_front();
    int acc=0;
    Eigen::MatrixXd mat(nbp-1,D);
    for(typename std::list<T_PTS>::iterator pit = lp.begin(); pit != lp.end(); ++pit)
    {
        T_PTS pn = *pit;
        for(int d = 0 ; d < D ; d++)
            mat(d,acc) = pn[d] - p0[d];
        acc++;
    }
    return fabs(mat.determinant()/factorial(D));
}


template<typename T_PTS>
T_PTS
compute_pts_proj(const T_PTS & A, const T_PTS & C, const T_PTS & v1, int D)
{
    double lambda = compute_coef_proj(A,C,v1,D);
    T_PTS p(D);
    for(int d = 0; d < D; d++)
    {
        p[d] = A[d] + lambda*v1[d];
    }
    return p;
}


template<typename T_PTS,typename T_VECT>
double
compute_coef_proj(const T_PTS & A, const T_PTS & C, const T_VECT & v1, int D)
{
    double acc1 =0;
    double acc2 =0;
    for(int d = 0; d < D; d++)
    {
        acc1 += (C[d]-A[d])*v1[d];
        acc2 += v1[d]*v1[d];
    }
    return acc1/acc2;
}



template<typename T_PTS>
std::vector<double>
compute_base_coef(const T_PTS & A, const T_PTS & C, const std::vector<T_PTS> & norms, int D)
{
    std::vector<double> coefs(D);
    for(int d = 0; d < D; d++)
    {
        coefs[d] = compute_coef_proj(A,C,norms[d],D);
    }
    return coefs;
}


template <typename Point>
double squared_dist(Point & p1, Point & p2, int D)
{
    double acc = 0;
    for(int d = 0; d < D; d++)
        acc += (p1[d] - p2[d])*(p1[d] - p2[d]);
    return acc;
}



void regularize(double & a, double & b, double & c);


template <typename Point>
std::vector<double>  get_barycenter(std::list<Point> & lp, int D)
{
    std::vector<double> coords(D);
    for(int d = 0; d < D; d++)
        coords[d] = 0;
    for(typename std::list<Point>::iterator pit = lp.begin(); pit != lp.end(); ++pit)
    {
        Point pn = *pit;
        for(int d = 0 ; d < D ; d++)
            coords[d] += pn[d];
    }
    for(int d = 0; d < D; d++)
        coords[d] /= lp.size();
    return coords;
}

template <typename Point,typename Traits>
Point get_barycenter_pts(std::list<Point> & lp, int D)
{
    Traits traits;
    double coords[Traits::D];
    for(int d = 0; d < D; d++)
        coords[d] = 0;
    for(typename std::list<Point>::iterator pit = lp.begin(); pit != lp.end(); ++pit)
    {
        Point pn = *pit;
        for(int d = 0 ; d < D ; d++)
            coords[d] += pn[d];
    }
    for(int d = 0; d < D; d++)
        coords[d] /= lp.size();
    return traits.make_point(coords);
}







template <typename Point,typename Traits>
double n_surface(typename std::list<Point> & lp, int D)
{
    Traits traits;
    Point p0 = get_barycenter_pts<Point,Traits>(lp,D);
    int acc=0;
    Eigen::MatrixXd mat(D,D);
    for(auto pit = lp.begin(); pit != lp.end(); ++pit)
    {
        Point pn = *pit;
        for(int d = 0 ; d < D ; d++)
            mat(acc,d) = pn[d] - p0[d];
        acc++;
    }
    Eigen::MatrixXd m = mat.rowwise() - mat.colwise().mean();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd sv = svd.singularValues() ;
    Eigen::MatrixXd ev =  svd.matrixV();
    std::vector<Point> nb;
    double coords[Traits::D];
    std::vector<double> coords_v(D);
    for(int i = 0; i < D; i++)
    {
        for(int d = 0 ; d < D ; d++)
            coords[d] = ev(d,i);
        nb.push_back(traits.make_point(coords));
    }
    std::list<Point> nbl;
    Point ori = p0;
    for(auto pit = lp.begin();
            pit != lp.end();
            ++pit)
    {
        Point p = *pit;
        coords_v = compute_base_coef(ori,p,nb,D);
        coords_v.pop_back();
        for(int d = 0; d < D; d++)
            coords[d] = coords_v[d];
        nbl.push_back(traits.make_point(coords));
    }
    return n_volume(nbl,D-1);
}












#endif

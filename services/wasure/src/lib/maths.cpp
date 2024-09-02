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
#include "wasure_maths.hpp"

double score_pdf(double a,double scale)
{
    return exp(-(a/scale)*(a/scale));
}

void regularize2(double & a)
{
    if(a > 1) a = 1;
    if(a < 0) a = 0;
}

void regularize(double & a, double & b, double & c)
{
    regularize2(a);
    regularize2(b);
    c = 1-b-a;
    if(c < 0)
    {
        c = 0;
    }
    double norm = a+b+c;
    a = a/norm;
    b = b/norm;
    c = c/norm;
}









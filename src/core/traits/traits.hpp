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
#ifndef DDT_TRAITS_HPP
#define DDT_TRAITS_HPP

#if defined(DDT_CGAL_TRAITS_2)

#include "cgal_traits_2.hpp"
namespace ddt
{
template <typename DataV,typename DataC> using Traits = ddt::Cgal_traits_2<DataV,DataC>;
}

#elif defined(DDT_CGAL_TRAITS_3)

#include "cgal_traits_3.hpp"
namespace ddt
{
template <typename DataV,typename DataC> using Traits = ddt::Cgal_traits_3<DataV,DataC>;
using Traits_raw = ddt::Cgal_traits_3_Raw;
}

#elif defined(DDT_CGAL_TRAITS_D)

#include "cgal_traits_d.hpp"
namespace ddt
{
template <typename DataV,typename DataC> using Traits = ddt::Cgal_traits<DDT_CGAL_TRAITS_D,DataV,DataC>;
using Traits_raw = ddt::Cgal_traits_raw<DDT_CGAL_TRAITS_D>;
}

#endif

#endif // DDT_TRAITS_HPP

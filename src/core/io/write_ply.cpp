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
#include "write_ply.hpp"

namespace ddt
{

void write_ply_header_begin(std::ostream& out)
{
    out << "ply" << std::endl;
    out << "format binary_little_endian 1.0" << std::endl;
    out << "comment creator: Mathieu Bredif" << std::endl;
}

void write_ply_header_end(std::ostream& out)
{
    out << "end_header" << std::endl;
}

}


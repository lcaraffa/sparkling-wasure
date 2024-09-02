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
#include "logging_stream.hpp"
#include <iostream>

namespace ddt
{

logging_stream::logging_stream( const std::string& s, int l) : level(l), time(), overall(s)
{
    time = last = start = std::chrono::system_clock::now();
}

void logging_stream::dump_log(std::ostream & oos)
{
    step("[total]");
    time = std::chrono::system_clock::now();
    operator()(0, std::chrono::duration<float>(time-start).count(), "");
    last = time;
    oos << sstr.str();
}

void logging_stream::step(const std::string& s)
{
    time = std::chrono::system_clock::now();
    if(last!=start)
    {
        operator()(2, "");
        operator()(0, std::chrono::duration<float>(time-last).count(), ";");
    }
    last = time;
    operator()(0, overall + "_"+ s, ":");
}

void logging_stream::do_log()
{
    return;
}

}

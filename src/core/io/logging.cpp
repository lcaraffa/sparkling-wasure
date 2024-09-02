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
#include "logging.hpp"
#include <iostream>

namespace ddt
{

logging::logging(const std::string& s, int l) : level(l), time(), overall(s)
{
    time = last = start = std::chrono::system_clock::now();
}

logging::~logging()
{
    step(overall);
    time = std::chrono::system_clock::now();
    operator()(0, std::chrono::duration<float>(time-start).count(), "\n");
    last = time;
    std::cerr << std::endl;
}

void logging::step(const std::string& s) const
{
    time = std::chrono::system_clock::now();
    if(last!=start)
    {
        operator()(2, "\n");
        operator()(0, std::chrono::duration<float>(time-last).count(), "\n");
    }
    last = time;
    operator()(0, s, "\t");
}

void logging::do_log() const
{
    std::cerr << std::flush;
}

}

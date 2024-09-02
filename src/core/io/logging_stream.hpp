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
#ifndef DDT_LOGGING_STREAM_HPP
#define DDT_LOGGING_STREAM_HPP

#include <chrono>
#include <iostream>
#include <sstream>

namespace ddt
{

class logging_stream
{
public:
    logging_stream(const std::string& s, int l);

    template<typename...Args>
    void operator()(int l, Args&&... args)
    {
        if(level>=l) do_log(args...);
    }
    void step(const std::string& s) ;
    void dump_log(std::ostream & oos) ;
    int level;
private:

    void do_log() ;
    template<typename Arg, typename...Args>
    void do_log(Arg&& arg, Args&&... args)
    {
        sstr << arg;
        do_log(args...);
    }

    mutable std::chrono::time_point<std::chrono::system_clock> time, last, start;
    std::stringstream sstr;
    std::string overall;
};

}

#endif // DDT_LOGGING_STREAM_HPP

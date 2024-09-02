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

#ifndef DDT_SPARK_UTILS_H
#define DDT_SPARK_UTILS_H


void parse_string(const std::string& str, std::vector<std::string>& elements, char delimeter)
{
    std::stringstream ss(str);
    std::string element;
    while(std::getline(ss, element, delimeter))
    {
        elements.push_back(element);
    }
}


std::string double_2_string(const double dbl)
{
    std::ostringstream strs;
    strs << dbl;
    std::string str = strs.str();
    std::replace(str.begin(),str.end(),'.',',');
    return strs.str();
}


std::string remove_path(const std::string & filename)
{
    return filename.substr(filename.find_last_of("\\/")+1,filename.size()-1);
}

std::string remove_extension(const std::string& filename)
{
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}


std::string get_bname(const std::string & filename)
{
    return remove_path(remove_extension(filename));
}


std::string extract_last_digit(const std::string & str)
{
    std::size_t last_index = str.find_last_not_of("0123456789");
    std::string result = str.substr(last_index + 1);
    return result;
}


int get_tri_idx(const std::string & str)
{
    return std::stoi(extract_last_digit(remove_extension(str)));
}

std::string switch_extension(const std::string& filename,const std::string& ext)
{
    std::stringstream ss2;
    std::string bname =  remove_extension(filename);
    ss2 << bname << "." << ext;
    return std::string(ss2.str());
}

#endif




















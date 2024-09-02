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
#ifndef IQ_HEADER_H
#define IQ_HEADER_H

#include <iomanip>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <utility>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <chrono>
#include <thread>
#include <sys/stat.h>

#include <random>
#include <algorithm>
#include <ctime>
#include <chrono>

#include <limits>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "io/logging_stream.hpp"
#include "conf_header/conf.hpp"

#define MAXBUFLEN (10000000)
#define SPARK_BUF_SIZE (65536)
#define IS_BINARY false

typedef std::numeric_limits< double > dbl;

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


// trim from start (in place)
static inline void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
    {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
    {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s)
{
    ltrim(s);
    rtrim(s);
}



#define OUTPUT_DIR ("/mnt/shared_spark/out")
namespace ddt
{

class stream_data_header
{
public :

    stream_data_header() : fis(NULL),fos(NULL),ifile(NULL),ofile(NULL),bool_dump(true),log(NULL)
    {
    }

    stream_data_header(std::string ll,std::string tt,std::vector<int> v) : stream_data_header()
    {
        lab = ll;
        type = tt;
        lidx = v;
    };

    stream_data_header(std::string ll,std::string tt,int id)  : stream_data_header()
    {
        lab = ll;
        type = tt;
        lidx = std::vector<int>(1,id);
    };


    std::string random_string( size_t length );

    void finalize();

    std::string  get_lab()
    {
        return lab;
    }

    void  set_lab(std::string ll )
    {
        lab = ll;
    }


    std::string get_ext()
    {
        return filename.substr(filename.find_last_of(".") + 1);
    }


    int get_id(int i)
    {
        return lidx[i];
    }

    int get_id()
    {
        return lidx[0];
    }

    bool is_stream();
    bool is_serialized();
    void serialize(bool);

    bool is_file() ;

    bool is_hdfs() ;

    bool fexists(const std::string& filename)
    {
        std::ifstream iff(filename.c_str());
        return (bool)iff;
    }


    char get_nl_char()
    {
        if(is_file())
            return '\n';
        return ';';
    }

    std::string get_file_name()
    {
        return filename;
    }

    void write_into_file(std::string root_dir,std::string ext, bool rand_ext = false);
    void set_file_name(std::string fname);


    std::istream & parse_header(std::istream & ifs, bool is_binary = true);

    std::ostream & write_header(std::ostream & ost, bool is_binary = true);

    void set_logger(  ddt::logging_stream * ll)
    {
        log = ll;
    }

    std::istream & get_input_stream()
    {
        if(is_file()) // File
            return *fis;
        else  if(is_hdfs()) // Hdfs
            return *ifile;
        else if(is_stream())
            return *ifile;
        else return std::cin;
    }

    std::ostream & get_output_stream()
    {
        if(is_file())
        {
            fos->precision(dbl::max_digits10);
            return *fos;
        }
        else if(is_hdfs())
        {
            ofile->precision(dbl::max_digits10);
            return *ofile;
        }
        else
        {
            std::cout.precision(dbl::max_digits10);
            return std::cout;
        }
    }


private:
    std::string lab,type,filename;
    std::vector<int> lidx;
    std::ifstream * fis;
    std::ofstream * fos;
    std::stringstream * ifile;
    std::ostringstream * ofile;
    bool bool_dump;
    ddt::logging_stream * log;
};



class stream_app_header
{
public :

    stream_app_header() : tile_id(-1),nbd(-1)
    {
    }

    bool is_void()
    {
        return tile_id == -1;
    }
    int get_nb_dat()
    {
        return nbd;
    }
    std::istream & parse_header(std::istream & ist)
    {
        std::string fip;
        trim(fip);
        ist >> fip;
        if(fip.empty())
        {
            return ist;
        }
        try
        {
            if (fip.empty())
            {
                throw std::invalid_argument("Input string is empty");
            }
            if(is_number(fip))
            {
                tile_id = std::stoi(fip);
            }
            else
            {
                ist >> fip;
                tile_id = std::stoi(fip);
            }
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "[parse header] Invalid argument: " << e.what() << " for input: " << fip << std::endl;
        }
        catch (const std::out_of_range& e)
        {
            std::cerr << "Out of range: " << e.what() << " for input: " << fip << std::endl;
        }
        ist >> nbd;
        return ist;
    }
    int tile_id,nbd;
};

void regularize_sspath(std::string & ss)
{
    boost::replace_all(ss,"//","/");
}

int intRand(const int & min, const int & max)
{
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}


std::string time_in_HH_MM_SS_MMM()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << std::put_time(&bt, "%d-%m-%Y-%H-%M-%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    return oss.str();
}






std::string stream_data_header::random_string( size_t length )
{
    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d_%m_%Y_%H_%M_%S");
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[intRand(0,max_index )];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}

void stream_data_header::write_into_file(std::string root_name,std::string ext, bool rand_ext)
{
    std::string curname;
    if(rand_ext)
        curname = root_name  + "_" + random_string(10)  + ext;
    else
        curname = root_name   + ext;
    set_file_name(curname);
}



void stream_data_header::set_file_name(std::string fname)
{
    filename = fname;
    regularize_sspath(filename);
    type = "f";
}


void stream_data_header::finalize()
{
    if(is_file())
    {
        if(fis != NULL)
            fis->close();
        if(fos != NULL)
        {
            fos->flush();
            fos->close();
            delete fos;
        }
    }
    else
    {
        if(ifile != NULL)
        {
            delete ifile;
        }
    }
}
bool stream_data_header::is_stream()
{
    return (type.compare("s") == 0);
}

bool stream_data_header::is_serialized()
{
    return (type.compare("z") == 0);
}

void stream_data_header::serialize(bool bb)
{
    if(bb)
        type = "z";
}

bool stream_data_header::is_file()
{
    return (type.compare("f") == 0);
}

bool stream_data_header::is_hdfs()
{
    return (type.compare("h") == 0);
}





std::istream & stream_data_header::parse_header(std::istream & ist, bool is_binary)
{
    if(log != NULL) log->step("[read]parse_header");
    ist >> lab;
    int size_c;
    ist >> size_c;
    int idx;
    for(int i = 0; i < size_c; i++)
    {
        ist >> idx;
        lidx.push_back(idx);
    }
    if(lab.empty() || lidx.size() == 0)
    {
        std::cerr << "[ERROR] error during header parsing" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    ist >> type;
    if(is_file())
    {
        ist >> filename;
        fis = new std::ifstream();
        int acc_op = 0;
        bool is_open = false;
        while (acc_op++ < 100)
        {
            if(is_binary)
            {
                fis->open(filename,std::ios::binary);
            }
            else
            {
                fis->open(filename);
            }
            if (fis->is_open())
            {
                is_open = true;
                break;
            }
            else
            {
                std::cerr << "[ERROR] read header ERROR during opening: " << filename << " [" << acc_op << "]" << std::endl;
                std::cerr << filename << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(2000));
            }
        }
        if(!is_open)
        {
            std::exit(EXIT_FAILURE);
        }
    }
    else if (is_hdfs())
    {
        std::cerr << "[error] not suported anymore" << std::endl;
    }
    else if (is_stream())
    {
        if(log != NULL) log->step("[read]cin2sstream");
        std::string input;
        std::getline(std::cin, input);
        ifile = new std::stringstream(input);
    }
    else
    {
        std::cout << " ";
        if(ifile != NULL)
        {
            delete ifile;
        }
    }
    return ist;
}
std::ostream & stream_data_header::write_header(std::ostream & ost, bool is_binary)
{
    if(lab.empty() || lidx.size() == 0)
    {
        std::cerr << "[ERROR] error during header writing" << std::endl;
        std::exit (EXIT_FAILURE);
    }
    ost << lab << " "  << lidx.size() << " ";
    for(auto ll : lidx)
        ost << ll << " ";
    ost <<  type << " ";
    if(is_file())
    {
        ost << filename << " ";
        boost::filesystem::path path(filename);
        fos = new std::ofstream();
        int acc_op = 0;
        bool is_open = false;
        while (acc_op++ < 20)
        {
            if(is_binary)
            {
                fos->open(filename,std::ios::binary);
            }
            else
            {
                fos->open(filename, std::ios::out);
            }
            if(fos->is_open())
            {
                is_open = true;
                break;
            }
            else
            {
                std::cerr << "[ERROR] ERROR durring writing header: " << filename << " " << std::endl;
                std::cerr << filename << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            }
        }
        if(!is_open)
        {
            std::exit(EXIT_FAILURE);
        }
    }
    else if (is_hdfs())
    {
        ost << filename << " ";
        ofile = new std::ostringstream();
    }
    return ost;
}
}

#endif




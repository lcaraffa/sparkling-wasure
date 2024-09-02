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
#ifndef NUMBER_PARSER_HPP
#define NUMBER_PARSER_HPP

#include "count.hpp"
#include "string-to-double.h"
#include "double-to-string.h"
#include "base64.hpp"



#define kBufferSize (10)
#define NB_DIGIT_OUT (5)


const  int flags_ser = double_conversion::DoubleToStringConverter::UNIQUE_ZERO |
                       double_conversion::DoubleToStringConverter::EMIT_POSITIVE_EXPONENT_SIGN;

const int flags_deser = double_conversion::DoubleToStringConverter::UNIQUE_ZERO |
                        double_conversion::DoubleToStringConverter::EMIT_POSITIVE_EXPONENT_SIGN;

const char gDigitsLut[200] =
{
    '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
    '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
    '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
    '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
    '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
    '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
    '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
    '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
    '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
    '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9'
};


const char cdigits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };


template<typename TT>
void  deserialize_b64_vect(std::vector<TT> & vv,std::istream & ifile)
{
    int bufflen_in,nb_dat;
    long unsigned int bufflen_out;
    char cc;
    ifile >> nb_dat >> bufflen_in;
    ifile.get(cc);
    unsigned char * buffer_char_in = new unsigned char[bufflen_in];
    unsigned char * buffer_char_out;
    ifile.read(((char *)buffer_char_in),bufflen_in);
    buffer_char_out = fast_base64_decode(buffer_char_in,bufflen_in,&bufflen_out);
    vv.clear();
    vv.resize(nb_dat);
    memcpy(&vv[0], buffer_char_out, bufflen_out);
    delete [] buffer_char_out;
    delete [] buffer_char_in;
}


template<typename TT>
void serialize_b64_vect(std::vector<TT> & vv,std::ostream & ss)
{
    int buff_size = 65536;
    int nbdat = vv.size();
    int bufflen_in = sizeof(TT)*nbdat;
    unsigned long int bufflen_out;
    unsigned char * buff_in = new unsigned char[bufflen_in];
    memcpy(buff_in, &vv[0], bufflen_in);
    vv.clear();
    unsigned char * buff_out = fast_base64_encode(buff_in, bufflen_in,&bufflen_out);
    ss << nbdat << " " <<  bufflen_out << " ";
    int acc = 0;
    ss.write(reinterpret_cast<char*>(buff_out+acc), bufflen_out - acc);
    ss << " ";
    delete [] buff_in;
    delete [] buff_out;
}


template<typename TT>
void serialize_b64_vect(const std::vector<TT> & vv,std::ostream & ss)
{
    int buff_size = 65536;
    int nbdat = vv.size();
    int bufflen_in = sizeof(TT)*nbdat;
    unsigned long int bufflen_out;
    unsigned char * buff_in = new unsigned char[bufflen_in];
    memcpy(buff_in, &vv[0], bufflen_in);
    unsigned char * buff_out = fast_base64_encode(buff_in, bufflen_in,&bufflen_out);
    ss << nbdat << " " <<  bufflen_out << " ";
    int acc = 0;
    ss.write(reinterpret_cast<char*>(buff_out+acc), bufflen_out - acc);
    ss << " ";
    delete [] buff_in;
    delete [] buff_out;
}


int u32toa_countlut(uint32_t value, char* buffer);
int count_nbd(const char * str);
inline int u32toa_count(uint32_t value, char* buffer);
int i32toa_count(int32_t value, char* buffer);
void u32toa_lut(int value, char* buffer) ;


#endif

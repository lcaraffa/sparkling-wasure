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
#include "number_parser.hpp"



int u32toa_countlut(uint32_t value, char* buffer)
{
    unsigned digit = CountDecimalDigit32(value);
    buffer += digit;
    *buffer = ' ';
    while (value >= 100)
    {
        const unsigned i = (value % 100) << 1;
        value /= 100;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }
    if (value < 10)
    {
        *--buffer = char(value) + '0';
    }
    else
    {
        const unsigned i = value << 1;
        *--buffer = gDigitsLut[i + 1];
        *--buffer = gDigitsLut[i];
    }
    buffer[digit++] = ' ';
    return digit++;
}

int count_nbd(const char * str)
{
    int val = 0;
    while( *str != ' ' )
    {
        str++;
        val++;
    }
    return val;
}


inline int u32toa_count(uint32_t value, char* buffer)
{
    unsigned digit = CountDecimalDigit32(value);
    buffer += digit;
    *buffer = ' ';
    do
    {
        *--buffer = char(value % 10) + '0';
        value /= 10;
    }
    while (value > 0);
    buffer[digit++] = ' ';
    return digit++;
}

int i32toa_count(int32_t value, char* buffer)
{
    uint32_t u = static_cast<uint32_t>(value);
    if (value < 0)
    {
        *buffer++ = '-';
        u = ~u + 1;
    }
    return u32toa_count(u, buffer);
}


void u32toa_lut(int value, char* buffer)
{
    char temp[10];
    char* p = temp;
    while (value >= 100)
    {
        const unsigned i = (value % 100) << 1;
        value /= 100;
        *p++ = gDigitsLut[i + 1];
        *p++ = gDigitsLut[i];
    }
    if (value < 10)
        *p++ = char(value) + '0';
    else
    {
        const unsigned i = value << 1;
        *p++ = gDigitsLut[i + 1];
        *p++ = gDigitsLut[i];
    }
    do
    {
        *buffer++ = *--p;
    }
    while (p != temp);
    *buffer = '\0';
}




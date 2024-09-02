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

#ifndef ALGO_PARAMS
#define ALGO_PARAMS

#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <list>
#include <algorithm>

enum Estep
{
    ERR = -1,
    INI = 0,
    DIM = 1,
    TES = 2,
    TRI = 3,
    DST = 4,
    SEG = 5,
    REG = 6,
    EXT = 7,
};




class algo_params
{
public :
    algo_params() : lambda(1), terr(0.1), rat_ray_sample(0.1), max_it(10), nb_labs(6), dim(2), nb_threads(1),rat_extra_pts(0.2),nb_samples(1),
        center_type(1),min_scale(0.1),verbose_flag(0),extra_datas(0),tile_id(-1), estep(ERR),ech_input(1),area_processed(0),
        input_tri_dir(std::string("")), output_tri_dir(std::string("")),  input_wasure_dir(std::string("")), output_wasure_dir(std::string("")),
        bbox_string(std::string("")), tiling_name(std::string("")) {};

    double lambda,terr,rat_ray_sample,rat_extra_pts,min_scale,ech_input;
    int nb_samples,max_it,nb_labs,dim,nb_threads,area_processed;
    int center_type,tile_id;
    int verbose_flag,extra_datas;
    Estep estep;
    std::list<std::string> list_name;

    std::string bbox_string,tiling_name,input_tri_dir,output_tri_dir,input_wasure_dir,output_wasure_dir;
    std::ostream& operator<<(std::ostream& os)
    {
        os << "lambda         : " << lambda << std::endl;
        os << "terr           : " << terr << std::endl;
        os << "rat_ray_sample : " << rat_ray_sample << std::endl;
        os << "max_it         : " << max_it << std::endl;
        os << "nb_labs        : " << nb_labs << std::endl;
        os << "dim            : " << nb_labs << std::endl;
        return os;
    }


    // ---
    // --- help function ----
    void help_wasure()
    {
        std::cerr << "--------------------------------------------------" << std::endl
                  << "INPUTS" << std::endl
                  << "--input_dir <path> : intput_data_dir " << std::endl
                  << "--output_dir <path> : output_data_dir " << std::endl
                  << "--graph <filenames.xml> : the merging graph " << std::endl
                  << "--tiling <filenames.xml> : the tiling xml " << std::endl
                  << "--dim <int> : dimensionality " << std::endl
                  << "[ --lambda <value> ] weight of the global prior " << std::endl
                  << "[ --help ]  : print this message" << std::endl
                  << "--------------------------------------------------" << std::endl ;
    }



    int parse(int argc, char **argv)
    {
        int cc;
        int errflg = 5;
        static struct option long_options[] =
        {
            // /* These options set a flag. */
            {"verbose", no_argument,       &verbose_flag, 1},
            //{"files", required_argument,     0, 'i'},
            {"input_wasure_dir", required_argument,     0, 'i'},
            {"output_wasure_dir", required_argument,     0, 'o'},
            {"input_tri_dir", required_argument,     0, 'u'},
            {"output_tri_dir", required_argument,     0, 'v'},
            {"bbox", required_argument,     0, 'g'},
            {"tiling", required_argument,     0, 't'},
            {"tile_id", required_argument,     0, 'x'},
            {"step",  required_argument, 0, 's'},
            {"center_type",  required_argument, 0, 'c'},
            {"area_processed",  required_argument, 0, 'q'},


            {"ech_input",  required_argument, 0, 'j'},
            {"lambda",  required_argument, 0, 'l'},
            {"terr",  required_argument, 0, 'e'},
            {"nb_samples",  required_argument, 0, 'f'},
            {"rat_ray_sample",  required_argument, 0, 'k'},
            {"tessel_target_err",  required_argument, 0, 'b'},
            {"rat_extra_pts",  required_argument, 0, 'p'},
            {"max_it",  required_argument, 0, 'm'},

            {"nb_labs",  required_argument, 0, 'n'},
            {"dim",  required_argument, 0, 'd'},
            {"nb_threads",  required_argument, 0, 'a'},
            {"extra_datas",  required_argument, 0, 'z'},

            {"help",  no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        while ((cc = getopt_long(argc, argv, "i:o:g:t:x:s:l:e:j:f:k:m:q:n:d:a:b:u:v:w:h",long_options,&option_index)) != -1)
        {
            switch (cc)
            {
            case 'i':
                input_wasure_dir = std::string(optarg);
                errflg--;
                break;
            case 'o':
                output_wasure_dir = std::string(optarg);
                errflg--;
                break;
            case 'u':
                input_tri_dir = std::string(optarg);
                errflg--;
                break;
            case 'v':
                output_tri_dir = std::string(optarg);
                errflg--;
                break;
            case 'g':
                bbox_string = std::string(optarg);
                std::replace( bbox_string.begin(), bbox_string.end(), ':', ' '); // replace all 'x' to 'y'
                std::replace( bbox_string.begin(), bbox_string.end(), 'x', ' ');
                errflg--;
                break;
            case 't':
                tiling_name = std::string(optarg);
                errflg--;
                break;
            case 's':
                estep = static_cast<Estep>(atoi(optarg));
                errflg--;
                break;
            case 'c':
                center_type = (atoi(optarg));
                break;
            case 'q':
                area_processed = atoi(optarg);
                break;
            case 'l':
                lambda = atof(optarg);
                break;
            case 'j':
                ech_input = atof(optarg);
                break;
            case 'e':
                terr = atof(optarg);
                break;
            case 'f':
                nb_samples = atoi(optarg);
                break;
            case 'k':
                rat_ray_sample = atof(optarg);
                break;
            case 'p':
                rat_extra_pts = atof(optarg);
                break;
            case 'm':
                max_it = atof(optarg);
                break;
            case 'n':
                nb_labs = atoi(optarg);
                break;
            case 'd':
                dim = atoi(optarg);;
                break;
            case 'x':
                tile_id = atoi(optarg);;
                break;
            case 'a':
                nb_threads = atoi(optarg);;
                break;
            case 'b':
                terr = atof(optarg);;
                break;
            case 'z':
                extra_datas = atoi(optarg);;
                break;
            case 'h':
                char * optc = (argv[optind]);
                std::cout << "value:" << (char)cc << "  :  " << optc  << std::endl;
                help_wasure();
                return 0;
            }
        }
        if(errflg != 0)
        {
            //help();
            // return 1;
        }
    }

};





#endif  // ALGO_PARAMS

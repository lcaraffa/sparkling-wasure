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
#ifndef WASURE_PARAMS
#define WASURE_PARAMS

#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <list>
#include <random>
#include <iostream>

class wasure_params
{
public :
    wasure_params() : verbose_flag(0),nbp(0),log_level(2),nbt_side(-1),dump_ply(false),area_processed(0),coef_mult(1),dst_scale(-1),
        input_dir(std::string("")),output_dir(std::string("")),  algo_step(std::string("")),slabel(std::string("")),mode(-1),rat_ray_sample(0),pscale(1),nb_samples(1),lambda(1),nb_labs(2),graph_type(0),tau(0.5),skip_app_header(false),ray_weight(0)
    {
    };
    int verbose_flag,seed,nbp,log_level,id_padding,graph_type,area_processed;
    bool show_ghost,skip_app_header,dump_ply,process_only_shared;
    double lambda,terr,rat_ray_sample,rat_extra_pts,min_scale,pscale,coef_mult,tau,dst_scale,ray_weight;
    int nb_samples,max_it,nb_labs,nb_threads,mode;
    int center_type,tile_id,nbt_side;
    bool use_weight = true;
    bool dump_debug  = false;
    bool do_finalize = false;
    std::string bbox_string,bba_ori_string,input_dir,output_dir,algo_step,slabel,filename;
    std::ostream& operator<<(std::ostream& os)
    {
        std::default_random_engine er((unsigned int)time(0));
        seed = er();
        os << "input_dir : " << input_dir << std::endl;
        os << "algo_step      : " << algo_step << std::endl;
        return os;
    }

    void help_param()
    {
        std::cerr << "--------------------------------------------------" << std::endl
                  << "INPUTS" << std::endl
                  << "[TODO]" << std::endl
                  << "[ --help ]  : print this message" << std::endl
                  << "--------------------------------------------------" << std::endl ;
    }


    int parse(int argc, char **argv)
    {
        int cc;
        int errflg = 0;
        static struct option long_options[] =
        {
            // /* These options set a flag. */
            {"seed",  required_argument, 0, 'a'},
            {"bbox",  required_argument, 0, 'b'},
            {"lambda",  required_argument, 0, 'c'},
            {"dim",  required_argument, 0, 'd'},
            {"ray_weight", required_argument,0, 'e'},
            {"pscale", required_argument,0, 'f'},
            {"step",  required_argument, 0, 's'},
            {"nbp",  required_argument, 0, 'n'},
            {"filename",  required_argument, 0, 'o'},
            {"bbox_ori",  required_argument, 0, 'i'},
            {"input_dir",  required_argument, 0, 'r'},
            {"mode",  required_argument, 0, 'm'},
            {"output_dir",  required_argument, 0, 'w'},
            {"label",  required_argument, 0, 'l'},
            {"coef_mult",  required_argument, 0, 'k'},
            {"rat_ray_sample",  required_argument, 0, 'u'},
            {"dst_scale", required_argument, 0, 't'},
            {"nb_samples",  required_argument, 0, 'j'},
            {"area_processed",  required_argument, 0, 'p'},
            {"nbt_side", required_argument,0, 'g'},
            {"dump_debug",  no_argument, 0, 'q'},
            {"dump_ply", no_argument,0, 'x'},
            {"do_finalize", no_argument,0, 'z'},
            {"verbose", no_argument, &verbose_flag, 1},
            {"help",  no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        while ((cc = getopt_long(argc, argv, "s:a:c:k:n:i:d:u:m:e:o:f:j:l:b:p:r:w:t:g:qxzh",long_options,&option_index)) != -1)
        {
            switch (cc)
            {
            case 's':
                algo_step = std::string(optarg);
                break;
            case 'r':
                input_dir = std::string(optarg);
                break;
            case 'o':
                filename = std::string(optarg);
                break;
            case 'w':
                output_dir = std::string(optarg);
                break;
            case 'm':
                mode = atoi(optarg);
                break;
            case 'l':
                slabel = std::string(optarg);
                break;
            case 'q':
                dump_debug = true;
                break;
            case 'b':
                bbox_string = std::string(optarg);
                std::replace( bbox_string.begin(), bbox_string.end(), ':', ' ');
                std::replace( bbox_string.begin(), bbox_string.end(), 'x', ' ');
                break;
            case 'n':
                nbp = atoi(optarg);
                break;
            case 'c':
                lambda = atof(optarg);
                break;
            case 'k':
                coef_mult = atof(optarg);
                break;
            case 'e':
                ray_weight = atof(optarg);
                break;
            case 'u':
                rat_ray_sample = atof(optarg);
                break;
            case 'd':
                break;
            case 'p':
                area_processed = atoi(optarg);
                break;
            case 'f':
                pscale = atof(optarg);
                break;
            case 'a':
                seed = atoi(optarg);
                break;
            case 'i':
                bba_ori_string = std::string(optarg);
                std::replace( bba_ori_string.begin(), bba_ori_string.end(), ':', ' ');
                std::replace( bba_ori_string.begin(), bba_ori_string.end(), 'x', ' ');
                break;
            case 't':
                dst_scale = atof(optarg);
                break;
            case 'j':
                nb_samples = atoi(optarg);
                break;
            case 'g':
                nbt_side = atoi(optarg);
                break;
            case 'x':
                dump_ply=true;
                break;
            case 'z':
                do_finalize=true;
                break;
            case 'h':
                char * optc = (argv[optind]);
                std::cout << "value:" << (char)cc << "  :  " << optc  << std::endl;
                help_param();
                return 0;
            }
        }
        if(errflg != 0)
        {
            help_param();
            return 1;
        }
        return 0;
    }

};





#endif  // WASURE_PARAMS


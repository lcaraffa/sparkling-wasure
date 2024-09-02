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

#ifndef DDT_STREAM_PARAMS
#define DDT_STREAM_PARAMS

#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <list>
#include <random>
#include <iostream>

class algo_params
{
public :
    algo_params() : nbt_side(1),dump_ply(false),verbose_flag(0),nbp(0),log_level(2),ech_input(1),min_ppt(0),plot_lvl(0),max_ppt(200000),show_ghost(false),finalize_tri(false),extract_edg_nbrs(false),extract_tri_crown(false),io_mode(0),area_processed(0),dump_mode("NONE"),do_send_empty_edges(false),do_simple_output(false),input_dir(std::string("")),output_dir(std::string("")),  algo_step(std::string("")),
        style(std::string("tri1.qml")),slabel(std::string("")) {};
    double ech_input;
    int nbt_side,verbose_flag,io_mode,seed,nbp,log_level,min_ppt,max_ppt,plot_lvl,area_processed;
    bool show_ghost,do_simple_output,extract_edg_nbrs,extract_tri_crown,do_send_empty_edges,finalize_tri,dump_ply;
    std::string bbox_string,input_dir,output_dir,algo_step,slabel,style,dump_mode;
    std::ostream& operator<<(std::ostream& os)
    {
        std::default_random_engine er((unsigned int)time(0));
        seed = er();
        os << "input_dir : " << input_dir << std::endl;
        os << "algo_step      : " << algo_step << std::endl;
        return os;
    }

    // ---
    // --- help function ----
    void help_param()
    {
        std::cerr << "--------------------------------------------------" << std::endl
                  << "INPUTS" << std::endl
                  << "--algo_step <string> : input_data_dir " << std::endl
                  << "--dim <int> : output_data_dir " << std::endl
                  << "[ --help ]  : print this message" << std::endl
                  << "--------------------------------------------------" << std::endl ;
    }


    int parse(int argc, char **argv)
    {
        int cc;
        int errflg = 1;
        static struct option long_options[] =
        {
            // /* These options set a flag. */
            {"step",  required_argument, 0, 's'},
            {"nbp",  required_argument, 0, 'n'},
            {"seed",  required_argument, 0, 'a'},
            {"input_dir",  required_argument, 0, 'r'},
            {"ech_input",  required_argument, 0, 'j'},
            {"output_dir",  required_argument, 0, 'w'},
            {"label",  required_argument, 0, 'l'},
            {"bbox",  required_argument, 0, 'b'},
            {"area_processed",  required_argument, 0, 'q'},
            {"nbt_side",  required_argument, 0, 't'},
            {"min_ppt",  required_argument, 0, 'p'},
            {"max_ppt",  required_argument, 0, 'k'},
            {"style",  required_argument, 0, 'u'},
            {"dump_mode", required_argument,0, 'z'},
            {"plot_lvl", required_argument,0, 'y'},
            {"io_mode", required_argument,0, 'x'},
            {"show_ghost", no_argument,0, 'g'},
            {"send_empty_edges", no_argument,0, 'i'},
            {"dump_ply",  no_argument, 0, 'd'},
            {"extract_edg_nbrs", no_argument,0, 'e'},
            {"extract_tri_crown", no_argument,0, 'c'},
            {"finalize_tri", no_argument,0, 'm'},
            {"do_simple_output", no_argument,0, 'f'},
            {"verbose", no_argument, &verbose_flag, 1},
            {"help",  no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        while ((cc = getopt_long(argc, argv, "s:a:n:u:t:l:b:z:j:p:q:k:y:r:w:m:x:gidecfh",long_options,&option_index)) != -1)
        {
            switch (cc)
            {
            case 's':
                algo_step = std::string(optarg);
                errflg--;
                break;
            case 'r':
                input_dir = std::string(optarg);
                break;
            case 'w':
                output_dir = std::string(optarg);
                break;
            case 'u':
                style = std::string(optarg);
                break;
            case 'l':
                slabel = std::string(optarg);
                break;
            case 'b':
                bbox_string = std::string(optarg);
                std::replace( bbox_string.begin(), bbox_string.end(), ':', ' '); // replace all 'x' to 'y'
                std::replace( bbox_string.begin(), bbox_string.end(), 'x', ' ');
                break;
            case 'n':
                nbp = atoi(optarg);
                break;
            case 'p':
                min_ppt = atoi(optarg);
                break;
            case 'd':
                dump_ply = true;
                break;
            case 'k':
                max_ppt = atoi(optarg);
                break;
            case 'j':
                ech_input = atof(optarg);
                break;
            case 'a':
                seed = atoi(optarg);
                break;
            case 'y':
                plot_lvl = atoi(optarg);
                break;
            case 'q':
                area_processed = atoi(optarg);
                break;
            case 't':
                nbt_side = atoi(optarg);
                break;
            case 'g':
                show_ghost=true;
                break;
            case 'm':
                finalize_tri=true;
                break;
            case 'i':
                do_send_empty_edges=true;
                break;
            case 'e':
                extract_edg_nbrs=true;
                break;
            case 'z':
                dump_mode= std::string(optarg);
                break;
            case 'c':
                extract_tri_crown=true;
                break;
            case 'f':
                do_simple_output=true;
                break;
            case 'x':
                io_mode = atoi(optarg);
                break;
            case 'h':
                break;
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





#endif  // DDT_STREAM_PARAMS

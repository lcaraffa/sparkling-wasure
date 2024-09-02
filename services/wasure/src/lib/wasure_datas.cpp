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

#include "wasure_data.hpp"

void wasure_data::parse_points_ply(std::string inputName)
{
    std::cout << "\t [loading file lbb]open " << inputName << "..." << std::endl;
    // Tinyply can and will throw exceptions at you!
    try
    {
        // Read the file and create a std::istringstream suitable
        // for the lib -- tinyply does not perform any file i/o.
        std::ifstream ss(inputName, std::ios::binary);
        // Parse the ASCII header fields
        PlyFile file(ss);
        for (auto e : file.get_elements())
        {
            std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
            for (auto p : e.properties)
            {
                std::cout << "\tproperty - " << p.name << " (" << PropertyTable[p.propertyType].str << ")" << std::endl;
            }
        }
        std::cout << std::endl;
        // Define containers to hold the extracted data. The type must match
        // the property type given in the header. Tinyply will interally allocate the
        // the appropriate amount of memory.
        std::vector<float> verts;
        std::vector<float> norms;
        std::vector<float> origins;
        std::vector<float> dims_n;
        std::vector<float> dims_s;
        uint32_t vertexCount, normalCount, originCount,dims_nCount,dims_sCount;
        vertexCount = normalCount = originCount = dims_nCount = dims_sCount;
        // The count returns the number of instances of the property group. The vectors
        // above will be resized into a multiple of the property group size as
        // they are "flattened"... i.e. verts = {x, y, z, x, y, z, ...}
        vertexCount = file.request_properties_from_element("vertex", {"x", "y", "z"}, verts);
        normalCount = file.request_properties_from_element("vertex", {"nx", "ny", "nz"}, norms);
        originCount = file.request_properties_from_element("vertex", {"x_origin", "y_origin", "z_origin"}, origins);
        dims_nCount = file.request_properties_from_element("vertex", {"eigenVector1x", "eigenVector1y", "eigenVector1z","eigenVector2x", "eigenVector2y", "eigenVector2z","eigenVector3x", "eigenVector3y", "eigenVector3z"}, dims_n);
        dims_sCount = file.request_properties_from_element("vertex", {"sigma1", "sigma2", "sigma3"}, dims_s);
        // Now populate the vectors...
        //timepoint before = now();
        file.read(ss);
        //timepoint after = now();
        // Good place to put a breakpoint!
        //std::cout << "Parsing took " << difference_micros(before, after) << "μs: " << std::endl;
        std::cout << "\tRead " << verts.size() << " total vertices (" << vertexCount << " properties)." << std::endl;
        std::cout << "\tRead " << norms.size() << " total normals (" << normalCount << " properties)." << std::endl;
        std::cout << "\tRead " << origins.size() << " total vertex origins (" << originCount << " properties)." << std::endl;
        std::cout << "\tRead " << dims_s.size() << " total vertex dims_s (" << dims_sCount << " properties)." << std::endl;
        std::cout << "\tRead " << dims_n.size() << " total vertex dims_n (" << dims_nCount << " properties)." << std::endl;
        int acc = 0;
        std::vector<double> coords_v(3);
        for(int i = 0 ; i < verts.size(); i++)
        {
            coords_v[acc%3] = verts[i];
            if(++acc %3 == 0)
            {
                Point p(coords_v[0],coords_v[1],coords_v[2]);
                points.push_back(p);
            }
        }
        acc = 0;
        for(int i = 0 ; i < norms.size(); i++)
        {
            coords_v[acc%3] = norms[i];
            if(++acc %3 == 0)
            {
                Point p(coords_v[0],coords_v[1],coords_v[2]);
                norms_pts.push_back(p);
            }
        }
        acc=0;
        for(int i = 0 ; i < origins.size(); i++)
        {
            coords_v[acc%3] = origins[i];
            if(++acc %3 == 0)
            {
                Point p(coords_v[0],coords_v[1],coords_v[2]);
                centers_pts.push_back(p);
            }
        }
        std::vector<Point> act_vect;
        //dims_norms.push_back(
        acc=0;
        for(int i = 0 ; i < dims_n.size(); i++)
        {
            coords_v[acc%3] = dims_n[i];
            if(++acc %3 == 0)
            {
                Point p(coords_v[0],coords_v[1],coords_v[2]);
                act_vect.push_back(p);
            }
            if(acc %9 == 0)
            {
                dims_norms.push_back(act_vect);
                act_vect.clear();
            }
        }
        std::vector<double> act_scale;
        acc = 0;
        for(int i = 0 ; i < dims_s.size(); i++)
        {
            act_scale.push_back(dims_s[i]);
            if(++acc %3 == 0)
            {
                dims_scales.push_back(act_scale);
                act_scale.clear();
            }
        }
        std::cout << "checkv : " << points.size() << " " << centers_pts.size() << " " << dims_norms.size() << " " << dims_scales.size() << " " << std::endl;
    }
    catch (const std::exception & e)
    {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }
}

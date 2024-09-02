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
#include "write_vrt.hpp"

namespace ddt
{

void add_qgis_style(const std::string& filename,const std::string& stylename)
{
    boost::filesystem::path path(filename);
    boost::filesystem::path path_style_target = path.replace_extension(".qml");
    boost::filesystem::path path_style_source(std::string(STYLE_SOURCE_DIR) + stylename);
    boost::filesystem::copy_file(path_style_source, path_style_target, boost::filesystem::copy_option::overwrite_if_exists);
}

void write_vrt_header_vert(std::ofstream& csv, const std::string& filename)
{
    boost::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::ofstream f(filename);
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << stem <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << stem << ".csv</SrcDataSource>" << std::endl;
    f <<    "<SrcLayer>" << stem <<  "</SrcLayer>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbPoint</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"tid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"id\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
    add_qgis_style(filename, std::string("vert1.qml"));
    csv.open(path.replace_extension("csv").string());
    csv << "geom,tid,id" << std::endl;
}

void write_vrt_header_facet(std::ofstream& csv, const std::string& filename)
{
    boost::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::ofstream f(filename);
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << stem <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << stem << ".csv</SrcDataSource>" << std::endl;
    f <<    "<SrcLayer>" << stem <<  "</SrcLayer>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbLineString</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"tid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"local\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
    csv.open(path.replace_extension("csv").string());
    csv << "geom,tid,local" << std::endl;
}

void write_vrt_header_cell(std::ofstream& csv, const std::string& filename)
{
    boost::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::ofstream f(filename);
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << stem <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << stem << ".csv</SrcDataSource>" << std::endl;
    f <<    "<SrcLayer>" << stem <<  "</SrcLayer>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbPolygon</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"tid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"local\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"main\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"cid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"cid0\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"cid1\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"cid2\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
    add_qgis_style(filename,std::string("tri1.qml"));
    csv.open(path.replace_extension("csv").string());
    csv << "geom,tid,local,main,cid,cid0,cid1,cid2" << std::endl;
}

void write_vrt_header_tin(std::ofstream& csv, const std::string& filename)
{
    boost::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::ofstream f(filename);
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << stem <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << stem << ".csv</SrcDataSource>" << std::endl;
    f <<    "<SrcLayer>" << stem <<  "</SrcLayer>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbTIN</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"tid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"id\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
    csv.open(path.replace_extension("csv").string());
    csv << "geom,tid,id" << std::endl;
}

void write_vrt_header_bbox(std::ofstream& csv, const std::string& filename)
{
    boost::filesystem::path path(filename);
    std::string stem = path.stem().string();
    std::ofstream f(filename);
    f <<"<OGRVRTDataSource>" << std::endl;
    f <<  "<OGRVRTLayer name=\"" << stem <<  "\">" << std::endl;
    f <<    "<SrcDataSource relativeToVRT=\"1\">" << stem << ".csv</SrcDataSource>" << std::endl;
    f <<    "<SrcLayer>" << stem <<  "</SrcLayer>" << std::endl;
    f <<    "<LayerSRS>IGNF:LAMB93</LayerSRS> " << std::endl;
    f <<    "<GeometryType>wkbPolygon</GeometryType> " << std::endl;
    f <<    "<GeometryField encoding=\"WKT\" field=\"geom\"/> " << std::endl;
    f <<    "<Field name=\"tid\" type=\"Integer\"/>" << std::endl;
    f <<    "<Field name=\"bboxid\" type=\"Integer\"/>" << std::endl;
    f <<  "</OGRVRTLayer>" << std::endl;
    f <<"</OGRVRTDataSource>" << std::endl;
    add_qgis_style(filename, std::string("bbox1.qml"));
    csv.open(path.replace_extension("csv").string());
    csv << "geom,tid,bboxid" << std::endl;
}

}

#include <numeric>
#include <cmath>
#include <algorithm>
#include <duckdb.hpp>
#include "mining_sim.h"

namespace astay
{
    Datatwin::Datatwin()
    {
    }

    void Datatwin::setDBfilepath(const std::string &filepath)
    {
        db_filepath_ = filepath;
    }

    void Datatwin::setModelBlockfilepath(const std::string &filepath)
    {
        mb_filepath_ = filepath;
    }

    void Datatwin::setTopographyfilepath(const std::string &filepath)
    {
        tp_filepath_ = filepath;
    }

    void Datatwin::setSVGDir(const std::string &filepath)
    {
        svg_dir_ = filepath;
    }

    void Datatwin::finalizeAndStoreCut(MiningCut &cut, int level, int day)
    {
        int id = cuts_.size();
        cut.setLevel(level);
        cut.setExtractionDay(day);
        cut.setId(id);
        cut.correct();
        cut.computeBox();
        cuts_.push_back(cut);
        // Assign the cut id to its corresponding level.
        cut_id_by_level_[cut.level()].push_back(cut.id());
    }

    void Datatwin::finalizeAndStorePolyline(Polyline &polyline, int level)
    {
        int id = polylines_.size();
        polyline.setLevel(level);
        polyline.setId(id);
        polyline.computeCentroid();
        polylines_.push_back(polyline);
        polyline_id_by_level_[level].push_back(id);
    }

    void Datatwin::finalizeAndStoreGeoPoly(MiningGeoPoly &geopoly, int level, std::string &name)
    {
        int id = geopolygons_.size();
        geopoly.setName(name);
        geopoly.setLevel(level);
        geopoly.setId(id);
        geopoly.correct();
        geopoly.computeBox();
        geopolygons_.push_back(geopoly);
    }

    void Datatwin::ReadMiningCuts(bool is_verbose)
    {
        // Open the database in read-only mode.
        SQLite::Database db(db_filepath_, SQLite::OPEN_READONLY);

        // Vectors to store the extracted data.
        std::vector<float_t> x_cuts, y_cuts, z_cuts;
        std::vector<int> vertex_numbers, period_days;

        // SQL query to retrieve mining cuts data.
        std::string query_cuts = "SELECT x,y,z,vertice,period FROM tb_PlanDiario ORDER BY period ASC, vertice ASC";
        SQLite::Statement query(db, query_cuts);

        // Execute the query and populate vectors with the results.
        while (query.executeStep())
        {
            auto x = query.getColumn(0).getDouble();
            auto y = query.getColumn(1).getDouble();
            auto z = query.getColumn(2).getInt();
            auto vertice = query.getColumn(3).getInt();
            auto period = query.getColumn(4).getString();
            auto day = std::stoi(period.substr(period.find('_') + 1));

            x_cuts.push_back(x);
            y_cuts.push_back(y);
            z_cuts.push_back(z);
            vertex_numbers.push_back(vertice);
            period_days.push_back(day);
        }

        // Check if enough data is retrieved.
        if (vertex_numbers.size() < 3)
        {
            if (is_verbose)
                std::cout << "[ERROR] ReadMiningCuts | Doesn't have any lines." << std::endl;
            return;
        }
        if (is_verbose)
            std::cout << "[DEBUG] ReadMiningCuts | Has " << vertex_numbers.size() << " lines." << std::endl;

        // Calculate the offsets for x and y coordinates.
        xOffset_ = std::accumulate(x_cuts.begin(), x_cuts.end(), 0) / x_cuts.size();
        yOffset_ = std::accumulate(y_cuts.begin(), y_cuts.end(), 0) / y_cuts.size();
        for (auto &x : x_cuts)
            x -= xOffset_;
        for (auto &y : y_cuts)
            y -= yOffset_;

        // Process and store the mining cuts.
        int p_vertex_n = vertex_numbers[0] - 1;
        MiningCut cut;
        for (size_t i = 0; i < vertex_numbers.size(); i++)
        {
            int c_vertex_n = vertex_numbers[i];
            // Check for new cut based on vertex number.
            if (c_vertex_n != p_vertex_n + 1)
            {
                finalizeAndStoreCut(cut, z_cuts[i - 1], period_days[i - 1]);
                cut.clear();
            }
            // Add point to the current cut.
            point_t p{x_cuts[i], y_cuts[i]};
            cut.add(p);
            p_vertex_n = c_vertex_n;
        }

        // Finalize and store the last mining cut.
        finalizeAndStoreCut(cut, z_cuts.back(), period_days.back());
        cut.clear();

        // Output debug information if verbose mode is enabled.
        if (is_verbose)
        {
            std::cout << std::setprecision(2) << std::fixed;
            std::cout << "Offset: (" << xOffset_ << ", " << yOffset_ << ")" << std::endl;
            std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            for (const auto &cut : cuts_)
            {
                std::cout << "- " << cut << std::endl;
            }
            std::cout << "------------------------------------------------" << std::endl;
        }
    }

    void Datatwin::ReadMiningBlocks(bool is_verbose)
    {
        // Open the database in read-only mode.
        SQLite::Database sqlite(db_filepath_, SQLite::OPEN_READONLY);

        // Vectors to store the extracted data.
        float_t azimuth, degrees, radians;
        int x_size, y_size, z_size;

        // SQL query to retrieve mining cuts data.
        std::string query_setup = "SELECT x_size,y_size,z_size,azimuth FROM tb_SetupBlock WHERE Type = 'ModelBlocks'";
        SQLite::Statement query(sqlite, query_setup);

        // Execute the query and populate vectors with the results.
        while (query.executeStep())
        {
            x_size = query.getColumn(0).getInt();
            y_size = query.getColumn(1).getInt();
            z_size = query.getColumn(2).getInt();
            azimuth = query.getColumn(3).getDouble();
            // Extract radians from azimuth
            degrees = fmod(360 - azimuth, 360);
            radians = degrees * (M_PI / 180);
        }
        mb_height_ = z_size;
        // std::cout << x_size << "," << y_size << "," << mb_height_ << std::endl; // [DEBUG]
        // std::cout << azimuth << "," << degrees << "," << radians << std::endl;  // [DEBUG]

        // Set boundaries before read data from MiningBlocks
        auto bounding_area = computeBoundingArea(1.2); // margin
        std::string x_lower = std::to_string(xOffset_ + bg::get<bg::min_corner, 0>(bounding_area));
        std::string x_upper = std::to_string(xOffset_ + bg::get<bg::max_corner, 0>(bounding_area));
        std::string y_lower = std::to_string(yOffset_ + bg::get<bg::min_corner, 1>(bounding_area));
        std::string y_upper = std::to_string(yOffset_ + bg::get<bg::max_corner, 1>(bounding_area));

        // DuckDB object to connect with the parquet file.
        duckdb::DuckDB db(nullptr);
        duckdb::Connection con(db);

        // Add all the levels for the query
        std::string z_levels = "(";
        for (const auto &set : cut_id_by_level_)
        {
            auto level = set.first;
            z_levels += std::to_string(level + mb_height_ / 2) + ",";
        }
        z_levels += ")";

        // SQL query to retrieve mining cuts data.
        std::string query_filter = "ELEVATION IN " + z_levels + " AND EASTING>" + x_lower + " AND EASTING<" + x_upper + " AND NORTHING>" + y_lower + " AND NORTHING<" + y_upper;
        std::string query_mb = "SELECT * FROM read_parquet('" + mb_filepath_ + "') WHERE " + query_filter;
        auto result = con.Query(query_mb);
        auto count = result->RowCount();
        // std::cout << query_mb << " | #rows : " << count << std::endl; // [DEBUG]

        // Iterate over the chunks and store information
        std::vector<float_t> x_pts, y_pts;
        std::vector<int> z_pts, otype_v;
        std::vector<float_t> topo_v, dens_v, au_v, ag_v, cu_v, tcm_v, s_v;

        auto chunk = result->Fetch();
        while (chunk != nullptr)
        {
            // Extract data
            auto x_vec = chunk->data[0];
            auto y_vec = chunk->data[1];
            auto z_vec = chunk->data[2];
            auto topo_vec = chunk->data[3];
            auto otype_vec = chunk->data[4];
            auto dens_vec = chunk->data[5];
            auto au_vec = chunk->data[6];
            auto ag_vec = chunk->data[7];
            auto cu_vec = chunk->data[8];
            auto tcm_vec = chunk->data[9];
            auto s_vec = chunk->data[10];
            for (int i = 0; i < chunk->size(); i++)
            {
                auto x = x_vec.GetValue(i).GetValue<double>();
                auto y = y_vec.GetValue(i).GetValue<double>();
                auto z = z_vec.GetValue(i).GetValue<int>();
                auto topo = topo_vec.GetValue(i).GetValue<double>();
                auto otype = otype_vec.GetValue(i).GetValue<int>();
                auto dens = dens_vec.GetValue(i).GetValue<double>();
                auto au = au_vec.GetValue(i).GetValue<double>();
                auto ag = ag_vec.GetValue(i).GetValue<double>();
                auto cu = cu_vec.GetValue(i).GetValue<double>();
                auto tcm = tcm_vec.GetValue(i).GetValue<double>();
                auto s = s_vec.GetValue(i).GetValue<double>();
                // Store information on vectors
                x_pts.push_back(x - xOffset_);
                y_pts.push_back(y - yOffset_);
                z_pts.push_back(z);
                topo_v.push_back(topo);
                otype_v.push_back(otype);
                dens_v.push_back(dens);
                au_v.push_back(au);
                ag_v.push_back(ag);
                cu_v.push_back(cu);
                tcm_v.push_back(tcm);
                s_v.push_back(s);
            }
            chunk = result->Fetch();
        }

        // Check if enough data is retrieved.
        if (x_pts.size() < 3)
        {
            if (is_verbose)
                std::cout << "[ERROR] ReadModelBlocks | Doesn't have any blocks." << std::endl;
            return;
        }
        if (is_verbose)
            std::cout << "[DEBUG] ReadModelBlocks | Has " << x_pts.size() << " points." << std::endl;

        // Process and store the MiningBlocks.
        float_t delta_x = x_size / 2 * cos(radians);
        float_t delta_y = y_size / 2 * sin(radians);
        // std::cout << "delta_x: " << delta_x << ", delta_y: " << delta_y << std::endl; // [DEBUG]
        // std::cout << "x_pts: " << x_pts[0] << ", y_pts: " << y_pts[0] << std::endl; // [DEBUG]
        for (int i = 0; i < x_pts.size(); i++)
        {
            MiningBlock block;
            // Populate the mining block.
            point_t vertex_a{x_pts[i] + delta_x, y_pts[i] + delta_y};
            point_t vertex_b{x_pts[i] + delta_x, y_pts[i] - delta_y};
            point_t vertex_c{x_pts[i] - delta_x, y_pts[i] + delta_y};
            point_t vertex_d{x_pts[i] - delta_x, y_pts[i] - delta_y};
            block.add(vertex_a);
            block.add(vertex_b);
            block.add(vertex_c);
            block.add(vertex_d);
            block.correct();

            // Add particular parameters of blocks
            block.setType(ag_v[i]);
            block.setTopo(ag_v[i]);
            block.setDensity(ag_v[i]);
            block.setAu(ag_v[i]);
            block.setAg(ag_v[i]);
            block.setCu(ag_v[i]);
            block.setTcm(ag_v[i]);
            block.setS(ag_v[i]);

            // Add both level and unique id to each block
            int id = blocks_.size();
            int level = z_pts[i] - mb_height_ / 2;
            block.setLevel(level);
            block.setId(id);

            // Asociate the mining block with the cuts that intersects.
            for (int cut_id : cut_id_by_level_[level])
            {
                if (cuts_[cut_id].containsPolygon(block.geometry()))
                {
                    block.setStatus(true); // false by default
                    block_id_by_cut_id_[cut_id].push_back(id);
                }
            }
            blocks_.push_back(block);
        }

        if (is_verbose)
        {
            std::cout << std::setprecision(2) << std::fixed;
            std::cout << "Number of Model Blocks: " << blocks_.size() << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            for (const auto &block_in_cut : block_id_by_cut_id_)
            {
                auto cut_id = block_in_cut.first;
                auto level = cuts_[cut_id].level();
                auto num_blocks = block_in_cut.second.size();
                std::cout << "- Mining Cut - ID (" << cut_id << ") :\tLevel (" << level << ") | # Mining Model Blocks (" << num_blocks << ")" << std::endl;
            }
            std::cout << "------------------------------------------------" << std::endl;
        }
    }

    void Datatwin::ReadTopography(bool is_verbose)
    {
        auto bounding_area = computeBoundingArea(1.5, is_verbose); // margin
        std::string x_lower = std::to_string(xOffset_ + bg::get<bg::min_corner, 0>(bounding_area));
        std::string x_upper = std::to_string(xOffset_ + bg::get<bg::max_corner, 0>(bounding_area));
        std::string y_lower = std::to_string(yOffset_ + bg::get<bg::min_corner, 1>(bounding_area));
        std::string y_upper = std::to_string(yOffset_ + bg::get<bg::max_corner, 1>(bounding_area));

        // DuckDB object to connect with the parquet file.
        duckdb::DuckDB db(nullptr);
        duckdb::Connection con(db);

        // Add all the levels for the query
        std::set<int> levels;
        for (const auto &set : cut_id_by_level_)
        {
            auto level = set.first;
            // Insert current level and next level to the set
            levels.insert(level);
            levels.insert(level + mb_height_);
        }

        std::string z_levels = "(";
        for (const auto &level : levels)
        {
            z_levels += std::to_string(level) + ",";
        }
        z_levels += ")";

        // Design the query to get the polycurve points.
        std::string query_filter = "z IN " + z_levels + " AND x>" + x_lower + " AND x<" + x_upper + " AND y>" + y_lower + " AND y<" + y_upper;
        std::string query = "SELECT x,y,z,vertice,objeto FROM read_parquet('" + tp_filepath_ + "') WHERE " + query_filter + " ORDER BY z, objeto, vertice ASC";

        // Execute query and get the points.
        auto result = con.Query(query);
        auto count = result->RowCount();

        // Iterate over the chunks and store information
        std::vector<float_t> x_pts, y_pts;
        std::vector<int> z_pts, vertex_numbers, polyline_ids;

        auto chunk = result->Fetch();
        while (chunk != nullptr)
        {
            // Extract data
            auto x_vec = chunk->data[0];
            auto y_vec = chunk->data[1];
            auto z_vec = chunk->data[2];
            auto vertex_vec = chunk->data[3];
            auto id_vec = chunk->data[4];
            for (int i = 0; i < chunk->size(); i++)
            {
                auto x = x_vec.GetValue(i).GetValue<double>();
                auto y = y_vec.GetValue(i).GetValue<double>();
                auto z = z_vec.GetValue(i).GetValue<int>();
                auto vertex = vertex_vec.GetValue(i).GetValue<int>();
                auto id = id_vec.GetValue(i).GetValue<int>();
                x_pts.push_back(x - xOffset_);
                y_pts.push_back(y - yOffset_);
                z_pts.push_back(z);
                vertex_numbers.push_back(vertex);
                polyline_ids.push_back(id);
            }
            chunk = result->Fetch();
        }

        // Check if enough data is retrieved.
        if (vertex_numbers.size() < 2)
        {
            if (is_verbose)
                std::cout << "[ERROR] ReadTopography | Doesn't have any segments." << std::endl;
            return;
        }
        if (is_verbose)
            std::cout << "[DEBUG] ReadTopography | Has " << vertex_numbers.size() << " points." << std::endl;

        // Process and store the polylines.
        int p_vertex_n = vertex_numbers[0] - 1;
        Polyline polyline;
        for (size_t i = 0; i < vertex_numbers.size(); i++)
        {
            int c_vertex_n = vertex_numbers[i];
            // Check for new cut based on vertex number.
            if (c_vertex_n != p_vertex_n + 1)
            {
                finalizeAndStorePolyline(polyline, z_pts[i - 1]);
                polyline.clear();
            }
            // Add point to the current cut.
            point_t p{x_pts[i], y_pts[i]};
            polyline.add(p);
            p_vertex_n = c_vertex_n;
        }

        // Finalize and store the last mining cut.
        finalizeAndStorePolyline(polyline, z_pts.back());
        polyline.clear();

        // Print the query if verbose mode is enabled
        if (is_verbose)
        {
            std::cout << "Query: " << query << std::endl;
            std::cout << std::setprecision(2) << std::fixed;
            std::cout << "Number of Polylines: " << polylines_.size() << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            for (const auto &polyline : polylines_)
            {
                std::cout << "- " << polyline << std::endl;
            }
            std::cout << "------------------------------------------------" << std::endl;
        }
    }

    void Datatwin::ReadMiningGeoPolygons(bool is_verbose)
    {
        // Open the database in read-only mode.
        SQLite::Database db(db_filepath_, SQLite::OPEN_READONLY);
        for (const auto &set : cut_id_by_level_)
        {
            // Extract all the geopoly polygons that belongs to each level
            int level = set.first;

            std::vector<float_t> x_geo, y_geo;
            std::vector<int> z_geo, vertex_numbers;
            std::vector<std::string> geopoly_name;

            std::string query_geo = "SELECT x,y,z,vertice,typename FROM tb_Polygon WHERE Z=" + std::to_string(level) + " ORDER BY typename ASC, vertice ASC";
            SQLite::Statement query(db, query_geo);
            while (query.executeStep())
            {
                auto x = query.getColumn(0).getDouble();
                auto y = query.getColumn(1).getDouble();
                auto z = query.getColumn(2).getInt();
                auto vertice = query.getColumn(3).getInt();
                auto type_name = query.getColumn(4).getString();
                // Populate vectors
                x_geo.push_back(x - xOffset_);
                y_geo.push_back(y - yOffset_);
                z_geo.push_back(z);
                vertex_numbers.push_back(vertice);
                geopoly_name.push_back(type_name);
            }

            // Check if enough data is retrieved.
            if (vertex_numbers.size() < 3)
            {
                if (is_verbose)
                    std::cout << "[ERROR] Read GeoPolygons | Level " << level << " doesn't have any lines." << std::endl;
                continue;
            }
            if (is_verbose)
                std::cout << "[DEBUG] Read GeoPolygons | Level " << level << " have " << vertex_numbers.size() << " lines." << std::endl;

            // init extraction
            int p_vertex_n = vertex_numbers[0] - 1;
            MiningGeoPoly geopoly;
            for (size_t i = 0; i < vertex_numbers.size(); i++)
            {
                int c_vertex_n = vertex_numbers[i];
                if (c_vertex_n != p_vertex_n + 1) // if true, reading new geopoly
                {
                    finalizeAndStoreGeoPoly(geopoly, z_geo[i - 1], geopoly_name[i - 1]);
                    for (int cut_id : set.second)
                    {
                        if (cuts_[cut_id].containsPolygon(geopoly.geometry()))
                        {
                            geopoly.setStatus(true); // false by default
                            geopoly_id_by_cut_id_[cut_id].push_back(geopoly.id());
                        }
                    }
                    geopoly.clear();
                }
                point_t p{x_geo[i], y_geo[i]};
                geopoly.add(p);
                p_vertex_n = c_vertex_n;
            }

            // Finalize and store the last mining geopolygon.
            finalizeAndStoreGeoPoly(geopoly, z_geo.back(), geopoly_name.back());
            for (int cut_id : set.second)
            {
                if (cuts_[cut_id].containsPolygon(geopoly.geometry()))
                {
                    geopoly.setStatus(true);
                    geopoly_id_by_cut_id_[cuts_[cut_id].id()].push_back(geopoly.id());
                }
            }
            geopoly.clear();
        }

        if (is_verbose)
        {
            std::cout << std::setprecision(2) << std::fixed;
            std::cout << "Offset: (" << xOffset_ << ", " << yOffset_ << ")" << std::endl;
            std::cout << "Number of GeoPolygons: " << geopolygons_.size() << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            for (const auto &geopoly_in_cut : geopoly_id_by_cut_id_)
            {
                auto cut_id = geopoly_in_cut.first;
                auto level = cuts_[cut_id].level();
                auto num_polygons = geopoly_in_cut.second.size();
                std::cout << "- Mining Cut - ID (" << cut_id << ") :\tLevel (" << level << ") | # Mining GeoPolygons (" << num_polygons << ")" << std::endl;
                std::cout << "\t > Geopolygons IDs | ";
                for (const auto &geopoly_id : geopoly_in_cut.second)
                {
                    std::cout << geopoly_id << " | ";
                }
                std::cout << std::endl;
            }
            std::cout << "------------------------------------------------" << std::endl;
        }
    }

    void Datatwin::IntersectMiningCutsWithPolylines(bool is_verbose)
    {
        for (auto &set : cut_id_by_level_)
        {
            auto level = set.first;
            for (const auto &polyline_id : polyline_id_by_level_[level])
            {
                Polyline &polyline = polylines_[polyline_id];

                // Calculate the nearest centroid for the next level polyline.
                float_t min_dist{DBL_MAX};
                point_t next_lvl_centroid;
                for (const auto &next_polyline_id : polyline_id_by_level_[level + mb_height_])
                {
                    Polyline &next_polyline = polylines_[next_polyline_id];
                    auto dist = bg::distance(polyline.centroid(), next_polyline.centroid());
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        next_lvl_centroid = next_polyline.centroid();
                    }
                }
                // [DEBUG]
                // std::cout << "near centroid: " << next_lvl_centroid.x() << "," << next_lvl_centroid.y() << std::endl;

                // Modify structure of cuts that intersect with polyline
                std::vector<int> erase_cut_id_vec;
                std::vector<int> added_cut_id_vec;
                for (auto cut_id : set.second)
                {
                    MiningCut &cut = cuts_[cut_id];

                    // Points of intersection between the polyline and the polygon.
                    std::vector<point_t> inter_pts;
                    boost::geometry::intersection(polyline.linestring(), cut.geometry(), inter_pts);

                    // Continue to next cut if the polyline doesn't cross the polygon.
                    int n_cross = inter_pts.size() / 2;
                    if (n_cross == 0)
                        continue;

                    // Store the direction sense of each line's subsegment that cross the polygon.
                    MiningCut new_cut;
                    std::vector<bool> cross_sense(n_cross);
                    std::vector<linestring_t> cross_linestring(2 * n_cross);
                    for (int i = 0; i < n_cross; i++)
                    {
                        // Compute direction of sense for each cross.
                        bool sense = isLeftOfLine(inter_pts[2 * i], inter_pts[2 * i + 1], next_lvl_centroid);
                        cross_sense[i] = sense;

                        // Compute all the points inside each cross.
                        linestring_t linestring = linestringFromPolyline(polyline, inter_pts[2 * i], inter_pts[2 * i + 1]);
                        cross_linestring[2 * i] = linestring;
                        bg::reverse(linestring);
                        cross_linestring[2 * i + 1] = linestring;
                        // [DEBUG]
                        // std::cout << std::setprecision(2) << std::fixed;
                        // std::cout << "Cross " << i << ": next level centroid is left? - " << cross_sense[i] << std::endl;
                    }

                    // Check each cut vertex and add any vertex that is on the correct sense from the segment.
                    polygon_t poly = cut.geometry();
                    point_t prev_point;
                    int prev_flag{0};
                    for (auto it = boost::begin(bg::exterior_ring(poly)); it != boost::end(bg::exterior_ring(poly)); ++it)
                    {
                        // Check status (valid or not) from each segment perspective.
                        int valid_ctr = 0;
                        for (int i = 0; i < n_cross; i++)
                        {
                            bool sense = isLeftOfLine(inter_pts[2 * i], inter_pts[2 * i + 1], *it);
                            valid_ctr += int(sense == cross_sense[i]);
                        }

                        // If a vertex is set as valid from all the segments perspective, add the vertex.
                        if (valid_ctr == n_cross)
                        {
                            // std::cout << "[VALID] \t: (" << bg::get<0>(*it) << "," << bg::get<1>(*it) << ")" << std::endl;
                            // First add the line point.
                            if (prev_flag < 0) // I2V
                            {
                                float_t min_dist{FLT_MAX};
                                point_t cross_point;
                                for (auto &point : inter_pts)
                                {
                                    float_t dist = abs(distanceFromLine(prev_point, *it, point));
                                    if (dist < min_dist)
                                    {
                                        min_dist = dist;
                                        cross_point = point;
                                    }
                                }
                                new_cut.add(cross_point);
                            }
                            new_cut.add(*it);
                            // Update previous values.
                            prev_point = *it;
                            prev_flag = 1;
                        }
                        else
                        {
                            // std::cout << "[INVALID] \t: (" << bg::get<0>(*it) << "," << bg::get<1>(*it) << ")" << std::endl;
                            // Add the segment cross-point.
                            if (prev_flag > 0) // V2I
                            {
                                float_t min_dist{FLT_MAX};
                                point_t cross_point;
                                int cross_index;
                                for (int i = 0; i < inter_pts.size(); i++)
                                {
                                    float_t dist = abs(distanceFromLine(prev_point, *it, inter_pts[i]));
                                    if (dist < min_dist)
                                    {
                                        min_dist = dist;
                                        cross_index = i;
                                        cross_point = inter_pts[i];
                                    }
                                }
                                new_cut.add(cross_point);

                                // Add all the points inside the cross linestring.
                                for (auto point : cross_linestring[cross_index])
                                {
                                    new_cut.add(point);
                                }
                            }
                            // Update previous values.
                            prev_point = *it;
                            prev_flag = -1;
                        }
                    }

                    // Add the new cut to the cuts_ vector.
                    int id = cuts_.size();
                    new_cut.setLevel(cut.level());
                    new_cut.setExtractionDay(cut.extractionDay());
                    new_cut.setId(id);
                    new_cut.correct();
                    new_cut.computeBox();
                    cuts_.push_back(new_cut);

                    // Add the new_cut id to the later added cuts list.
                    added_cut_id_vec.push_back(new_cut.id());

                    // Add the cut id to the later erased cuts list.
                    erase_cut_id_vec.push_back(cut.id());

                    // [DEBUG]
                    // for (const auto &point : inter_pts)
                    // {
                    //     std::cout << std::setprecision(2) << std::fixed;
                    //     std::cout << "Mining Cut [" << cut.id() << "] : (" << point.x() << "," << point.y() << ")" << std::endl;
                    // }
                    // std::cout << "- Mining Cut: \t" << cut << " | " << bg::dsv(cut.geometry()) << std::endl;
                    // std::cout << "- New Mining Cut: \t" << new_cut << " | " << bg::dsv(new_cut.geometry()) << std::endl;
                }

                // Delete all the cuts that are on the list.
                for (auto erased_cut_id : erase_cut_id_vec)
                {
                    // std::cout << "[DEBUG] Erased Mining Cut [" << erased_cut_id << "]" << std::endl;
                    auto it = set.second.begin();
                    while (it != set.second.end())
                    {
                        if (*it == erased_cut_id)
                        {
                            it = set.second.erase(it);
                            break;
                        }
                        it++;
                    }
                }

                // Add all the cuts that are on the list.
                for (auto added_cut_id : added_cut_id_vec)
                {
                    // std::cout << "[DEBUG] Added Mining Cut [" << added_cut_id << "]" << std::endl;
                    set.second.push_back(added_cut_id);
                }

                if (is_verbose)
                {
                    std::cout << "Polyline [" << polyline.id() << "]" << std::endl;
                    std::cout << "# Added Cuts (" << added_cut_id_vec.size() << ") | "
                              << "# Erased Cuts (" << erase_cut_id_vec.size() << ")" << std::endl;
                    for (const auto &cut_id : set.second)
                    {
                        std::cout << std::setprecision(2) << std::fixed;
                        std::cout << "- " << cuts_[cut_id] << std::endl;
                    }
                }
            }
        }
        // if (is_verbose)
        // {
        //     std::cout << std::setprecision(2) << std::fixed;
        //     std::cout << "Number of Polylines: " << polylines_.size() << std::endl;
        //     std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
        // }
    }

    void Datatwin::Draw(bool is_verbose)
    {
        std::string cuts_dir = svg_dir_ + "mining_cuts/";
        std::string blocks_dir = svg_dir_ + "mining_blocks/";
        std::string geopoly_dir = svg_dir_ + "mining_geopolygons/";

        // Draw Cuts and include plan points
        for (const auto &set : cut_id_by_level_)
        {
            auto bounding_area = computeBoundingArea(1.5); // margin
            auto x_min = bg::get<bg::min_corner, 0>(bounding_area);
            auto y_min = bg::get<bg::min_corner, 1>(bounding_area);
            auto x_max = bg::get<bg::max_corner, 0>(bounding_area);
            auto y_max = bg::get<bg::max_corner, 1>(bounding_area);
            auto level = set.first;

            // Calculate the bounding centroid and the svg window (height & width)
            point_t centroid{(x_max + x_min) / 2,
                             (y_max + y_min) / 2};
            auto svg_w = int(x_max - x_min);
            auto svg_h = int(y_max - y_min);

            // Set svg filepath and window
            std::string svg_file = cuts_dir + std::to_string(level) + ".svg";
            std::string width_height = "width=\"" + std::to_string(svg_w) + "\" height=\"" + std::to_string(svg_h) + "\"";

            // Create the svg file and draw the Mining-Cuts
            Drawer drawer(svg_file, svg_w, svg_h, width_height);
            if (is_verbose)
            {
                std::cout << "[../mining_cuts/" << level << ".svg] | centroid ("
                          << centroid.x() << ", " << centroid.y() << ") | box [("
                          << x_min << ", " << y_min << "), ("
                          << x_max << ", " << y_max << ")] | size ("
                          << svg_w << ", " << svg_h << ")" << std::endl;
            }
            for (const auto &cut_id : set.second)
            {
                MiningCut &cut = cuts_[cut_id];
                MiningCut new_cut = ApplyOffsetToPolygon(cut, centroid);
                point_t cut_centroid, new_cut_centroid;
                bg::centroid(cut.geometry(), cut_centroid);
                bg::centroid(new_cut.geometry(), new_cut_centroid);
                drawer.AddCut(new_cut, cut.extractionDay());
                drawer.AddNumber(cut.extractionDay(), new_cut_centroid);
                if (is_verbose)
                {
                    std::cout << "- " << cut
                              << " | area (" << bg::area(cut.geometry())
                              << ") | perimeter (" << bg::perimeter(cut.geometry())
                              << ") | centroid ("
                              << cut_centroid.x() << ","
                              << cut_centroid.y() << ") \t| affter offset ("
                              << new_cut_centroid.x() << ","
                              << new_cut_centroid.y() << ")" << std::endl;
                }
            }
            for (const auto &polyline_id : polyline_id_by_level_[level])
            {
                Polyline &polyline = polylines_[polyline_id];
                Polyline new_polyline = ApplyOffsetToPolyline(polyline, centroid);
                drawer.AddPolyline(new_polyline);
                if (is_verbose)
                {
                    std::cout << "- " << polyline << std::endl;
                }
            }
            for (const auto &polyline_id : polyline_id_by_level_[level + mb_height_])
            {
                Polyline &polyline = polylines_[polyline_id];
                Polyline new_polyline = ApplyOffsetToPolyline(polyline, centroid);
                drawer.AddPolyline(new_polyline);
                if (is_verbose)
                {
                    std::cout << "- " << polyline << std::endl;
                }
            }
        }
    }

    /**
     * Constructs a linestring from a polyline based on proximity to two points.
     *
     * This function processes a given polyline, creating a linestring that starts
     * from the point in the polyline closest to either of two specified points
     * (point_a or point_b) and continues until it reaches the second point.
     *
     * Parameters:
     *   polyline - The polyline to be processed.
     *   point_a - The first reference point.
     *   point_b - The second reference point.
     *
     * Returns:
     *   linestring_t - The constructed linestring from the polyline.
     */
    Datatwin::linestring_t Datatwin::linestringFromPolyline(Polyline &polyline, point_t point_a, point_t point_b)
    {
        linestring_t out_linestring;
        point_t prev_point;
        point_t other_point{point_b};
        bool invert_order{false};
        int ctr{0}, num_found{0};
        for (auto point : polyline.linestring())
        {
            // Initialization for the first point.
            if (ctr == 0)
            {
                prev_point = point;
                ctr++;
                continue;
            }

            // Process points to find the closest to point_a or point_b.
            if (num_found == 0)
            {
                auto dist_to_a = abs(distanceFromLine(prev_point, point, point_a));
                auto dist_to_b = abs(distanceFromLine(prev_point, point, point_b));
                auto min_dist = std::min(dist_to_a, dist_to_b);
                if (min_dist < 1.0E-05)
                {
                    num_found++;
                    if (dist_to_b < dist_to_a)
                    {
                        other_point = point_a;
                        invert_order = true;
                    }
                }
                // std::cout << "[INVALID] : (" << point.x() << "," << point.y() << ")" << std::endl;
            }
            else if (num_found == 1)
            {
                // Append points to the linestring until reaching the second point.
                auto dist = abs(distanceFromLine(prev_point, point, other_point));
                if (dist < 1.0E-05)
                    num_found++;
                bg::append(out_linestring, point);
                // std::cout << "[VALID] : (" << point.x() << "," << point.y() << ")" << std::endl;
            }
            else
            {
                // std::cout << "[INVALID] : (" << point.x() << "," << point.y() << ")" << std::endl;
                break;
            }
            prev_point = point;
            ctr++;
        }

        // Invert the order of the linestring to respect the order point_a -> linestring -> point_b.
        if (invert_order)
        {
            // std::cout << "[INVERT]" << std::endl;
            bg::reverse(out_linestring);
        }

        return out_linestring;
    }

    /**
     * Calculates the signed perpendicular distance from a point to a line.
     *
     * This function computes the signed distance from a point to a line defined
     * by two points (line_a and line_b). The distance is positive if the point
     * is on the left side of the line and negative if it's on the right side.
     *
     * Parameters:
     *   line_a - The first point defining the line.
     *   line_b - The second point defining the line.
     *   point_c - The point from which the distance to the line is calculated.
     *
     * Returns:
     *   float_t - The signed perpendicular distance from the point to the line.
     */
    Datatwin::float_t Datatwin::distanceFromLine(point_t &line_a, point_t &line_b, point_t &point_c)
    {
        point_t line_vec{line_b.x() - line_a.x(), line_b.y() - line_a.y()};
        point_t point_vec{point_c.x() - line_a.x(), point_c.y() - line_a.y()};
        float_t dist = line_vec.x() * point_vec.y() -
                       line_vec.y() * point_vec.x();
        return dist;
    }

    /**
     * Determines if a point is to the left of a line.
     *
     * By using the signed distance from the line (calculated in DistanceFromLine),
     * this function checks if a given point (point_c) is on the left side of a
     * line defined by two points (line_a and line_b).
     *
     * Parameters:
     *   line_a - The first point defining the line.
     *   line_b - The second point defining the line.
     *   point_c - The point to check its position relative to the line.
     *
     * Returns:
     *   bool - True if point_c is to the left of the line, False otherwise.
     */
    bool Datatwin::isLeftOfLine(point_t &line_a, point_t &line_b, point_t &point_c)
    {
        float_t dist = distanceFromLine(line_a, line_b, point_c);
        return (dist > 0);
    }

    /**
     * Compute the bounding box encompassing all cuts with a margin.
     *
     * This function calculates the outer area considering all cuts in each level,
     * expanding the area by a specified margin.
     *
     * Parameters:
     *   margin: The factor by which to expand the bounding box of the cuts.
     *           A margin of 1.0 means no expansion, 2.0 doubles the size, etc.
     *   is_verbose: If true, the function prints detailed output for debugging.
     *
     * Returns:
     *   A box_t representing the computed area.
     */
    Datatwin::box_t Datatwin::computeBoundingArea(float margen, bool is_verbose)
    {
        Datatwin::box_t bounding_area;
        float_t x_min{DBL_MAX}, y_min{DBL_MAX};
        float_t x_max{-DBL_MAX}, y_max{-DBL_MAX};

        // Iterate through each cut to find the extreme coordinates
        for (const auto &cut : cuts_)
        {
            auto cut_box = cut.box();
            float_t x_min_cut{bg::get<bg::min_corner, 0>(cut_box)};
            float_t y_min_cut{bg::get<bg::min_corner, 1>(cut_box)};
            float_t x_max_cut{bg::get<bg::max_corner, 0>(cut_box)};
            float_t y_max_cut{bg::get<bg::max_corner, 1>(cut_box)};

            // Update the overall min and max coordinates
            x_min = std::min(x_min_cut, x_min);
            y_min = std::min(y_min_cut, y_min);
            x_max = std::max(x_max_cut, x_max);
            y_max = std::max(y_max_cut, y_max);
        }

        // Calculate the deltas for expanding the box based on the margin
        float_t delta_x = (x_max - x_min) * (margen - 1) / 2;
        float_t delta_y = (y_max - y_min) * (margen - 1) / 2;

        // Set the expanded coordinates of the bounding area
        bg::set<bg::min_corner, 0>(bounding_area, x_min - delta_x);
        bg::set<bg::min_corner, 1>(bounding_area, y_min - delta_y);
        bg::set<bg::max_corner, 0>(bounding_area, x_max + delta_x);
        bg::set<bg::max_corner, 1>(bounding_area, y_max + delta_y);

        // Print the bounding area if verbose mode is enabled
        if (is_verbose)
        {
            std::cout << "Bounding Area: " << bg::dsv(bounding_area) << std::endl;
        }

        return bounding_area;
    }

    /*
        // -- read


        void Datatwin::ReadBlocks(bool is_verbose)
        {
            ModelBlock model_block_(modelblocks_db_);
            model_block_.set_offset(xOffset_, yOffset_);
            model_block_.read_setup();
            std::set<int> blocks_id;
            for (auto &cut : cuts_)
            {
                std::vector<MiningBlock> temp_blocks = model_block_.get_blocks(cut);
                for (MiningBlock block : temp_blocks)
                {
                    int id = block.id();
                    if (blocks_id.find(id) != blocks_id.end())
                        continue;
                    else
                    {
                        blocks_id.insert(id);
                        blocks_by_level_[block.level()].push_back(block);
                    }
                }
            }
        }

        void Datatwin::ReadPlan(bool is_verbose)
        {
            // Clear
            plan_.clear();
            // Extract from DB
            SQLite::Database db(db_file_path_, SQLite::OPEN_READONLY);
            std::string query_plan = "SELECT Source_,X_Source,Y_Source,Z_Source from Vista_Material WHERE Project_Id=" + std::to_string(project_id_);
            SQLite::Statement query(db, query_plan);
            while (query.executeStep())
            {
                std::string source = query.getColumn(0).getString();
                std::string loader = source.substr(0, source.find('_'));
                float_t x = query.getColumn(1).getDouble() - xOffset_;
                float_t y = query.getColumn(2).getDouble() - yOffset_;
                float_t z = query.getColumn(3).getDouble();
                // Populate vectors
                plan_.add_row(x, y, z, loader);
            }
            if (is_verbose)
            {
                std::cout << "Number of rows: " << plan_.size() << std::endl;
                std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
                std::cout << plan_ << std::endl;
                std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            }
        }

        // -- compute
        void Datatwin::SortCuts(bool is_verbose)
        {
            // Clear vectors
            cuts_by_loader_.clear();
            row_loaders_.clear();
            row_cuts_.clear();
            for (auto &cut : cuts_)
            {
                // Extract cut centroid
                point_t centroid;
                bg::centroid(cut.geometry(), centroid);

                int cut_vert = bg::num_points(cut.geometry());

                // Get the loader whose some position is nearest to the cut
                std::string loader;
                float_t min_dist = FLT_MAX;
                for (int i = 0; i < plan_.size(); i++)
                {
                    // Distance between cut centroid and each plan
                    float dist = sqrt(pow(plan_.x()[i] - centroid.x(), 2) +
                                      pow(plan_.y()[i] - centroid.y(), 2));
                    // Store neariest plan loader
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        loader = plan_.loaders()[i];
                    }
                }
                int cut_idx = cuts_by_loader_[loader].size() + 1;
                cut.setId(cut_idx);
                // Populate loader per vertice
                row_loaders_.insert(row_loaders_.end(), cut_vert, loader);
                row_cuts_.insert(row_cuts_.end(), cut_vert, cut_idx);
                // Populate a dict of cuts by loader
                cuts_by_loader_[loader].push_back(cut);
            }
            if (is_verbose)
            {
                std::cout << "Number of Loaders: " << cuts_by_loader_.size() << std::endl;
                std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
                std::cout << "Number of Vertices: " << row_loaders_.size() << std::endl;
                std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
                for (const auto &loader_cuts : cuts_by_loader_)
                {
                    std::cout << loader_cuts.first << ":" << std::endl;
                    for (const auto &cut : loader_cuts.second)
                    {
                        point_t cut_centroid;
                        bg::centroid(cut.geometry(), cut_centroid);
                        std::cout << "- " << cut
                                  << " | centroid ("
                                  << cut_centroid.x() << ","
                                  << cut_centroid.y() << ")" << std::endl;
                    }
                }
                std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            }
        }

        void Datatwin::ValidateDB(bool is_verbose)
        {
            SQLite::Database db(db_file_path_, SQLite::OPEN_READONLY);
            std::string query_direction = "SELECT COUNT(Project_ID) FROM MiningDirection WHERE Project_ID=" + std::to_string(project_id_);
            SQLite::Statement query(db, query_direction);
            int count;
            while (query.executeStep())
            {
                count = query.getColumn(0).getInt();
            }
            if (count > 0)
                is_sorted_ = true;
            else
                is_sorted_ = false;

            if (is_verbose)
                std::cout << "Validation Count: " << count << std::endl;
        }

        void Datatwin::PopulateDB(bool is_verbose)
        {
            if (is_sorted_)
            {
                if (is_verbose)
                    std::cout << "Database already updated" << std::endl;
            }
            else
            {
                SQLite::Database db(db_file_path_, SQLite::OPEN_READWRITE);
                SQLite::Transaction transaction(db);
                // Populate MiningCuts
                int idx = 0;
                for (const auto &loader : row_loaders_)
                {
                    int rowid = row_ids_[idx];
                    int cutid = row_cuts_[idx];
                    // populate 'Loader' and 'Cut_ID'
                    std::string query_loader = "UPDATE MiningCuts SET Loader='" + loader +
                                               "' WHERE rowid=" + std::to_string(rowid);
                    std::string query_cut = "UPDATE MiningCuts SET Cut_ID='" + std::to_string(cutid) +
                                            "' WHERE rowid=" + std::to_string(rowid);
                    db.exec(query_loader);
                    db.exec(query_cut);
                    idx++;
                    if (is_verbose)
                    {
                        std::cout << "Loader query:" << query_loader << std::endl;
                        std::cout << "Cut query:" << query_cut << std::endl;
                    }
                }
                // Populate MiningDirection
                for (const auto &loader_cuts : cuts_by_loader_)
                {
                    std::string shovel = loader_cuts.first;
                    for (const auto &cut : loader_cuts.second)
                    {
                        std::string query_direction = "INSERT INTO MiningDirection (Shovel,Cut_ID,Dir_X,Dir_Y,Project_ID) VALUES ('" +
                                                      shovel + "'," + std::to_string(cut.id()) + ",0.0,0.0," + std::to_string(project_id_) + ")";
                        db.exec(query_direction);
                        if (is_verbose)
                        {
                            std::cout << "Direction query:" << query_direction << std::endl;
                        }
                    }
                }
                transaction.commit();
            }
        }

        void Datatwin::Draw(bool is_verbose)
        {
            std::string cortes_dir = svg_dir_ + "cortes/";
            std::string blocks_dir = svg_dir_ + "blocks/";
            std::string ores_dir = svg_dir_ + "ores/";
            // Draw Cuts | include plan points
            for (const auto &loader_cuts : cuts_by_loader_)
            {
                // Extract loader centroid and outer box
                point_t min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
                point_t loader_centroid;
                // Compute box
                for (const auto &cut : loader_cuts.second)
                {
                    float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                    float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                    auto cut_box = cut.box();
                    float_t x_min_cut{bg::get<bg::min_corner, 0>(cut_box)};
                    float_t y_min_cut{bg::get<bg::min_corner, 1>(cut_box)};
                    float_t x_max_cut{bg::get<bg::max_corner, 0>(cut_box)};
                    float_t y_max_cut{bg::get<bg::max_corner, 1>(cut_box)};
                    if (x_min_cut < x_min)
                        bg::set<0>(min_corner, x_min_cut);
                    if (y_min_cut < y_min)
                        bg::set<1>(min_corner, y_min_cut);
                    if (x_max_cut > x_max)
                        bg::set<0>(max_corner, x_max_cut);
                    if (y_max_cut > y_max)
                        bg::set<1>(max_corner, y_max_cut);
                }
                bg::set<0>(loader_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
                bg::set<1>(loader_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
                int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
                int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
                // Draw after apply offset
                std::string svg_file = cortes_dir + loader_cuts.first + ".svg";
                std::string width_height = "width=\"" + std::to_string(svg_w) +
                                           "\" height=\"" + std::to_string(svg_h) + "\"";
                Drawer drawer(svg_file, svg_w, svg_h, width_height);
                if (is_verbose)
                {
                    std::cout << loader_cuts.first << " | centroid ("
                              << loader_centroid.x() << ", "
                              << loader_centroid.y() << ") | box [("
                              << bg::get<0>(min_corner) << ", "
                              << bg::get<1>(min_corner) << "), ("
                              << bg::get<0>(max_corner) << ", "
                              << bg::get<1>(max_corner) << ")] | size ("
                              << svg_w << ", " << svg_h << ")" << std::endl;
                }
                for (int i = 0; i < plan_.size(); i++)
                {
                    // apply offset
                    point_t plan_point;
                    bg::set<0>(plan_point, (plan_.x()[i] - loader_centroid.x()));
                    bg::set<1>(plan_point, (plan_.y()[i] - loader_centroid.y()));
                    // draw
                    drawer.AddPoint(plan_point);
                }
                for (const auto &cut : loader_cuts.second)
                {
                    point_t cut_moved_centroid;
                    MiningCut cut_moved = ApplyOffsetToPolygon(cut, loader_centroid);
                    bg::centroid(cut_moved.geometry(), cut_moved_centroid);
                    drawer.AddCut(cut_moved, cut.period_days());
                    drawer.AddNumber(cut.id(), cut_moved_centroid);
                    if (is_verbose)
                    {
                        point_t cut_centroid;
                        bg::centroid(cut.geometry(), cut_centroid);
                        std::cout << "- " << cut
                                  << " | centroid ("
                                  << cut_centroid.x() << ","
                                  << cut_centroid.y() << ") | modified ("
                                  << cut_moved_centroid.x() << ","
                                  << cut_moved_centroid.y() << ")" << std::endl;
                    }
                }

            }
            // Draw Blocks
            for (const auto &level_block : blocks_by_level_)
            {
                // Extract loader centroid and outer box
                point_t min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
                point_t level_centroid;
                // Compute box
                for (const auto &block : level_block.second)
                {
                    float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                    float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                    auto block_box = block.geometry();
                    float_t x_min_block{bg::get<bg::min_corner, 0>(block_box)};
                    float_t y_min_block{bg::get<bg::min_corner, 1>(block_box)};
                    float_t x_max_block{bg::get<bg::max_corner, 0>(block_box)};
                    float_t y_max_block{bg::get<bg::max_corner, 1>(block_box)};
                    if (x_min_block < x_min)
                        bg::set<0>(min_corner, x_min_block);
                    if (y_min_block < y_min)
                        bg::set<1>(min_corner, y_min_block);
                    if (x_max_block > x_max)
                        bg::set<0>(max_corner, x_max_block);
                    if (y_max_block > y_max)
                        bg::set<1>(max_corner, y_max_block);
                }
                bg::set<0>(level_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
                bg::set<1>(level_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
                int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
                int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
                // Draw after apply offset
                std::string svg_file = blocks_dir + std::to_string(level_block.first) + ".svg";
                std::string width_height = "width=\"" + std::to_string(svg_w) +
                                           "\" height=\"" + std::to_string(svg_h) + "\"";
                Drawer drawer(svg_file, svg_w, svg_h, width_height);
                if (is_verbose)
                {
                    std::cout << level_block.first << " | centroid ("
                              << level_centroid.x() << ", "
                              << level_centroid.y() << ") | box [("
                              << bg::get<0>(min_corner) << ", "
                              << bg::get<1>(min_corner) << "), ("
                              << bg::get<0>(max_corner) << ", "
                              << bg::get<1>(max_corner) << ")] | size ("
                              << svg_w << ", " << svg_h << ")" << std::endl;
                }
                for (const auto &block : level_block.second)
                {
                    bool paint = true;
                    MiningBlock block_moved = ApplyOffsetToPolygon(block, level_centroid);
                    block_moved.set_cu(block.cu());
                    drawer.AddBlock(block_moved, paint);
                    if (is_verbose)
                    {
                        point_t block_centroid, block_moved_centroid;
                        bg::centroid(block.geometry(), block_centroid);
                        bg::centroid(block_moved.geometry(), block_moved_centroid);
                        std::cout << "- " << block
                                  << " | centroid ("
                                  << block_centroid.x() << ","
                                  << block_centroid.y() << ") | modified ("
                                  << block_moved_centroid.x() << ","
                                  << block_moved_centroid.y() << ")" << std::endl;
                    }
                }
            }
            // Draw Ores
            for (const auto &level_ores : geopolygons_by_level_)
            {
                // Extract level centroid and outer box
                point_t min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
                point_t level_centroid;
                // Compute box
                for (const auto &geopoly : level_ores.second)
                {
                    float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                    float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                    auto ore_box = geopoly.box();
                    float_t x_min_ore{bg::get<bg::min_corner, 0>(ore_box)};
                    float_t y_min_ore{bg::get<bg::min_corner, 1>(ore_box)};
                    float_t x_max_ore{bg::get<bg::max_corner, 0>(ore_box)};
                    float_t y_max_ore{bg::get<bg::max_corner, 1>(ore_box)};
                    if (x_min_ore < x_min)
                        bg::set<0>(min_corner, x_min_ore);
                    if (y_min_ore < y_min)
                        bg::set<1>(min_corner, y_min_ore);
                    if (x_max_ore > x_max)
                        bg::set<0>(max_corner, x_max_ore);
                    if (y_max_ore > y_max)
                        bg::set<1>(max_corner, y_max_ore);
                }
                bg::set<0>(level_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
                bg::set<1>(level_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
                int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
                int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
                // Draw after apply offset
                std::string svg_file = ores_dir + std::to_string(level_ores.first) + ".svg";
                std::string width_height = "width=\"" + std::to_string(svg_w) +
                                           "\" height=\"" + std::to_string(svg_h) + "\"";
                Drawer drawer(svg_file, svg_w, svg_h, width_height);
                if (is_verbose)
                {
                    std::cout << level_ores.first << " | centroid ("
                              << level_centroid.x() << ", "
                              << level_centroid.y() << ") | box [("
                              << bg::get<0>(min_corner) << ", "
                              << bg::get<1>(min_corner) << "), ("
                              << bg::get<0>(max_corner) << ", "
                              << bg::get<1>(max_corner) << ")] | size ("
                              << svg_w << ", " << svg_h << ")" << std::endl;
                }
                for (const auto &geopoly : level_ores.second)
                {
                    MiningGeoPoly ore_moved = ApplyOffsetToPolygon(geopoly, level_centroid);
                    drawer.AddOre(ore_moved, geopoly.status());
                    if (is_verbose)
                    {
                        point_t ore_centroid, ore_moved_centroid;
                        bg::centroid(geopoly.geometry(), ore_centroid);
                        bg::centroid(ore_moved.geometry(), ore_moved_centroid);
                        std::cout << "- " << geopoly
                                  << " | centroid ("
                                  << ore_centroid.x() << ","
                                  << ore_centroid.y() << ") | modified ("
                                  << ore_moved_centroid.x() << ","
                                  << ore_moved_centroid.y() << ")" << std::endl;
                    }
                }
            }
        }
    */
}
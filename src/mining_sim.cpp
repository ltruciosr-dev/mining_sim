#include <numeric>
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

    // Helper function to finalize and store a Mining Cut.
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

    // Helper function to finalize and store a Polyline.
    void Datatwin::finalizeAndStorePolyline(Polyline &polyline, int level)
    {
        int id = polylines_.size();
        polyline.setLevel(level);
        polyline.setId(id);
        polyline.computeCentroid();
        polylines_.push_back(polyline);
        polyline_id_by_level_[level].push_back(id);
    }

    // Helper function to finalize and store a Mining GeoPolygon
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

    // Reads mining cut data from a database and processes it.
    //
    // This function opens a database connection, retrieves mining cut data,
    // and stores the data in internal structures for further processing.
    // It also calculates offsets and processes each mining cut.
    //
    // Parameters:
    //   is_verbose - If true, additional debug information will be printed to standard output.
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
                            geopoly_id_by_cut_id_[cuts_[cut_id].id()].push_back(geopoly.id());
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
        for (const auto &set : cut_id_by_level_)
        {
            auto level = set.first;
            for (const auto &polyline_id : polyline_id_by_level_[level])
            {
                Polyline &polyline = polylines_[polyline_id];

                // Calculate the nearest centroid for the next level polyline.
                float_t min_dist{DBL_MAX};
                point_t near_centroid;
                for (const auto &next_polyline_id : polyline_id_by_level_[level + mb_height_])
                {
                    Polyline &next_polyline = polylines_[next_polyline_id];
                    auto dist = bg::distance(polyline.centroid(), next_polyline.centroid());
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        near_centroid = next_polyline.centroid();
                    }
                }
                // [DEBUG]
                std::cout << "near centroid: " << near_centroid.x() << "," << near_centroid.y() << std::endl;

                // Modify structure of cuts that intersect with polyline
                for (const auto &cut_id : set.second)
                {
                    MiningCut &cut = cuts_[cut_id];
                    std::vector<point_t> inter_pts;
                    boost::geometry::intersection(polyline.linestring(), cut.geometry(), inter_pts);

                    // Check if the polyline cross the mining-cut
                    int n_cross = inter_pts.size() / 2;
                    if (n_cross > 0)
                    {
                        // Get vector of distances
                        std::vector<point_t> line_vectors(n_cross);
                        std::vector<int> line_dirs(n_cross);
                        for (int i = 0; i < n_cross; i++)
                        {
                            point_t line_vec{inter_pts[2 * i + 1].x() - inter_pts[2 * i].x(),
                                             inter_pts[2 * i + 1].y() - inter_pts[2 * i].y()};
                            point_t point_vec{near_centroid.x() - inter_pts[2 * i].x(),
                                              near_centroid.y() - inter_pts[2 * i].y()};
                            float_t dist = line_vec.x() * point_vec.y() -
                                           line_vec.y() * point_vec.x();
                            line_vectors[i] = line_vec;
                            line_dirs[i] = int(dist > 0);
                            // [DEBUG]
                            std::cout << std::setprecision(2) << std::fixed;
                            std::cout << i << ": [" << line_dirs[i] << "] | (" << line_vec.x() << "," << line_vec.y() << ")" << std::endl;
                        }

                        // Check each cut vertex to determine if they are inside or outside the cut
                        polygon_t poly = cut.geometry();
                        for (auto it = boost::begin(bg::exterior_ring(poly)); it != boost::end(bg::exterior_ring(poly)); ++it)
                        {
                            for (int i = 0; i < n_cross; i++)
                            {
                                point_t point_vec{bg::get<0>(*it) - inter_pts[2 * i].x(),
                                                  bg::get<1>(*it) - inter_pts[2 * i].y()};
                                float_t dist = line_vectors[i].x() * point_vec.y() -
                                               line_vectors[i].y() * point_vec.x();
                                // If higher than 0, the geometry vertex is valid.
                                if (dist * line_dirs[i] > 0)
                                {
                                    std::cout << "[VALID] \t: (" << bg::get<0>(*it) << "," << bg::get<1>(*it) << ")" << std::endl;
                                }
                                else
                                {
                                    std::cout << "[INVALID] \t: (" << bg::get<0>(*it) << "," << bg::get<1>(*it) << ")" << std::endl;
                                }
                            }
                            // float_t x = bg::get<0>(*it) - centroid.x();
                            // float_t y = bg::get<1>(*it) - centroid.y();
                            // point_t p{x, y};
                            // out_polygon.add(p);
                        }

                        // [DEBUG]
                        for (const auto &point : inter_pts)
                        {
                            std::cout << std::setprecision(2) << std::fixed;
                            std::cout << "Mining Cut [" << cut.id() << "] : (" << point.x() << "," << point.y() << ")" << std::endl;
                        }
                    }
                }
            }
        }

        if (is_verbose)
        {
            std::cout << std::setprecision(2) << std::fixed;
            std::cout << "Number of Polylines: " << polylines_.size() << std::endl;
            std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
        }
    }

    void Datatwin::Draw(bool is_verbose)
    {
        std::string cuts_dir = svg_dir_ + "mining_cuts/";
        std::string blocks_dir = svg_dir_ + "mining_blocks/";
        std::string geopoly_dir = svg_dir_ + "mining_geopolygons/";

        // Draw Cuts | include plan points
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
                std::cout << svg_file << " | centroid ("
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
                              << " \t| centroid ("
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

    /*
    void Datatwin::ReadMiningBlocks(bool is_verbose)
    {
        duckdb::DuckDB db(nullptr);
        duckdb::Connection con(db);

        // SQL query to retrieve mining cuts data.
        std::string query_cuts = "SELECT x,y,z,vertice,geometry FROM read_parquet('" + mb_filepath_ + "') WHERE z=4060";
    }
    */

    // Compute the bounding box encompassing all cuts with a margin.
    // The function calculates the outer area considering all cuts in each level,
    // expanding the area by a specified margin.
    //
    // Parameters:
    //   margin: The factor by which to expand the bounding box of the cuts.
    //           A margin of 1.0 means no expansion, 2.0 doubles the size, etc.
    //   is_verbose: If true, the function prints detailed output for debugging.
    //
    // Returns:
    //   A box_t representing the computed area.
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
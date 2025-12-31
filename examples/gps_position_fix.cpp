/**
 * gps_position_fix.cpp
 * 
 * A straightforward example showing how to compute a GPS position fix.
 * We'll set up some satellites, generate measurements, and solve for position.
 * 
 * Hamoud, 2025
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "../include/gnss_solver.hpp"

using namespace gnss;

// Quick helper to make a simple GPS constellation
// This isn't real satellite data - just something that works for testing
std::vector<SatelliteData> makeGPSConstellation(double time_sec = 0.0) {
    std::vector<SatelliteData> satellites;
    
    const double GPS_ALT = 20200000.0;  // ~20,200 km up
    const double GPS_R = constants::WGS84_A + GPS_ALT;
    const int SATS_PER_PLANE = 4;
    const int NUM_PLANES = 6;
    
    int prn = 1;
    for (int plane = 0; plane < NUM_PLANES; ++plane) {
        double raan = plane * 60.0 * constants::DEG2RAD;  // space out the planes
        double inc = 55.0 * constants::DEG2RAD;  // GPS orbit angle
        
        for (int slot = 0; slot < SATS_PER_PLANE; ++slot) {
            // Where is this satellite in its orbit?
            double M = slot * 90.0 * constants::DEG2RAD + plane * 15.0 * constants::DEG2RAD;
            double n = std::sqrt(constants::MU_EARTH / (GPS_R * GPS_R * GPS_R));
            M += n * time_sec;  // let it move a bit
            
            double nu = M;  // close enough for circular orbits
            
            // Position in orbit plane
            double x_orb = GPS_R * std::cos(nu);
            double y_orb = GPS_R * std::sin(nu);
            
            // Rotate into Earth-fixed frame
            double cos_raan = std::cos(raan);
            double sin_raan = std::sin(raan);
            double cos_i = std::cos(inc);
            double sin_i = std::sin(inc);
            
            double x = cos_raan * x_orb - sin_raan * cos_i * y_orb;
            double y = sin_raan * x_orb + cos_raan * cos_i * y_orb;
            double z = sin_i * y_orb;
            
            // Satellites have slightly different clock offsets
            double clock_bias = (prn % 7 - 3) * 1e-8;
            
            satellites.emplace_back(prn++, ECEF(x, y, z), clock_bias);
        }
    }
    
    return satellites;
}

int main() {
    std::cout << "=====================================\n";
    std::cout << "   GPS Position Fix Demo\n";
    std::cout << "=====================================\n\n";

    // Let's say we're standing in Austin, Texas
    LLA actual_pos(30.2672 * constants::DEG2RAD,
                   -97.7431 * constants::DEG2RAD,
                   150.0);
    
    ECEF actual_pos_ecef = llaToECEF(actual_pos);
    double actual_clock_bias = 1000.0;  // receiver clock is off by ~3 microseconds
    
    std::cout << "Where we actually are:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << actual_pos.latDeg() << "°N, " << actual_pos.lonDeg() << "°E\n";
    std::cout << std::setprecision(1);
    std::cout << "  " << actual_pos.alt << " m altitude\n";
    std::cout << "  ECEF: [" << actual_pos_ecef.x << ", " 
              << actual_pos_ecef.y << ", " << actual_pos_ecef.z << "]\n\n";

    // Make up a satellite constellation
    auto sats = makeGPSConstellation(0.0);
    std::cout << "Set up " << sats.size() << " GPS satellites\n\n";

    // Simulate measurements with some noise (real GPS is noisy!)
    double noise = 5.0;  // ~5 meters is typical
    auto measurements = simulatePseudoranges(actual_pos_ecef, actual_clock_bias, 
                                             sats, noise, 12345);
    
    std::cout << "Got measurements from " << measurements.size() << " satellites\n";
    std::cout << "Noise level: " << noise << " m\n\n";
    
    std::cout << "Measurements:\n";
    for (const auto& m : measurements) {
        std::cout << "  PRN " << std::setw(2) << m.prn << ": " 
                  << std::setprecision(1) << m.range << " m\n";
    }
    std::cout << "\n";

    // Solve it
    SolverConfig cfg;
    cfg.convergence_threshold = 1e-4;
    cfg.max_iterations = 20;
    cfg.apply_earth_rotation = true;  // account for earth spinning while signal travels
    
    GNSSSolver solver(cfg);
    NavSolution result = solver.solve(sats, measurements);

    std::cout << "Solution:\n";
    std::cout << "---------\n";
    printSolution(result);

    // How far off are we?
    ECEF err = result.position - actual_pos_ecef;
    double err_3d = err.norm();
    
    // Convert error to local coordinates (easier to understand)
    double s_lat = std::sin(actual_pos.lat);
    double c_lat = std::cos(actual_pos.lat);
    double s_lon = std::sin(actual_pos.lon);
    double c_lon = std::cos(actual_pos.lon);
    
    double east = -s_lon * err.x + c_lon * err.y;
    double north = -s_lat * c_lon * err.x - s_lat * s_lon * err.y + c_lat * err.z;
    double up = c_lat * c_lon * err.x + c_lat * s_lon * err.y + s_lat * err.z;
    double horiz = std::sqrt(east*east + north*north);
    
    std::cout << "\nHow accurate is it?\n";
    std::cout << std::setprecision(2);
    std::cout << "  3D error: " << err_3d << " m\n";
    std::cout << "  Horizontal: " << horiz << " m\n";
    std::cout << "  Vertical: " << std::abs(up) << " m\n";
    std::cout << "  (East: " << east << " m, North: " << north << " m)\n";
    
    double clock_err = result.clock_bias - actual_clock_bias;
    std::cout << "  Clock error: " << clock_err << " m\n";

    // Run it multiple times to see how consistent it is
    std::cout << "\n=====================================\n";
    std::cout << "   Running 100 trials...\n";
    std::cout << "=====================================\n\n";
    
    std::vector<double> errors_3d, errors_h, errors_v;
    
    for (int i = 0; i < 100; ++i) {
        auto m = simulatePseudoranges(actual_pos_ecef, actual_clock_bias,
                                      sats, noise, i * 1000);
        NavSolution sol = solver.solve(sats, m);
        
        if (sol.valid) {
            ECEF e = sol.position - actual_pos_ecef;
            
            double ee = -s_lon * e.x + c_lon * e.y;
            double en = -s_lat * c_lon * e.x - s_lat * s_lon * e.y + c_lat * e.z;
            double eu = c_lat * c_lon * e.x + c_lat * s_lon * e.y + s_lat * e.z;
            
            errors_3d.push_back(e.norm());
            errors_h.push_back(std::sqrt(ee*ee + en*en));
            errors_v.push_back(std::abs(eu));
        }
    }
    
    // Stats
    auto avg = [](const std::vector<double>& v) {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };
    auto rms = [](const std::vector<double>& v) {
        double sq = 0;
        for (double x : v) sq += x * x;
        return std::sqrt(sq / v.size());
    };
    auto pct = [](std::vector<double> v, double p) {
        std::sort(v.begin(), v.end());
        return v[static_cast<size_t>(p * v.size())];
    };
    
    std::cout << "Results:\n";
    std::cout << "  3D error: avg=" << avg(errors_3d) << "m, rms=" 
              << rms(errors_3d) << "m, 95%=" << pct(errors_3d, 0.95) << "m\n";
    std::cout << "  Horizontal: avg=" << avg(errors_h) << "m, rms=" 
              << rms(errors_h) << "m, 95%=" << pct(errors_h, 0.95) << "m\n";
    std::cout << "  Vertical: avg=" << avg(errors_v) << "m, rms=" 
              << rms(errors_v) << "m, 95%=" << pct(errors_v, 0.95) << "m\n";
    
    std::cout << "\nBased on DOP (HDOP=" << result.dop.hdop << ", VDOP=" 
              << result.dop.vdop << "), expect:\n";
    std::cout << "  Horizontal: ~" << result.dop.hdop * noise << " m\n";
    std::cout << "  Vertical: ~" << result.dop.vdop * noise << " m\n";

    std::cout << "\n=====================================\n";
    std::cout << "   Done!\n";
    std::cout << "=====================================\n";

    return 0;
}

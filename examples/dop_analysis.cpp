/**
 * dop_analysis.cpp
 * 
 * Playing around with satellite geometry to see how it affects
 * positioning accuracy. DOP = Dilution of Precision.
 * 
 * Hamoud, 2025
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/gnss_solver.hpp"

using namespace gnss;

// Create satellites at specific angles from the receiver
std::vector<SatelliteData> placeSatellites(const ECEF& rx_pos,
                                             const std::vector<std::pair<double, double>>& positions) {
    std::vector<SatelliteData> sats;
    LLA rx_lla = ecefToLLA(rx_pos);
    
    double s_lat = std::sin(rx_lla.lat);
    double c_lat = std::cos(rx_lla.lat);
    double s_lon = std::sin(rx_lla.lon);
    double c_lon = std::cos(rx_lla.lon);
    
    const double RANGE = 22000000.0;  // typical GPS range
    
    int prn = 1;
    for (const auto& [azimuth, elevation] : positions) {
        double az_rad = azimuth * constants::DEG2RAD;
        double el_rad = elevation * constants::DEG2RAD;
        
        // Direction in local ENU coordinates
        double e = std::sin(az_rad) * std::cos(el_rad);
        double n = std::cos(az_rad) * std::cos(el_rad);
        double u = std::sin(el_rad);
        
        // Convert to ECEF
        double dx = -s_lon * e - s_lat * c_lon * n + c_lat * c_lon * u;
        double dy = c_lon * e - s_lat * s_lon * n + c_lat * s_lon * u;
        double dz = c_lat * n + s_lat * u;
        
        ECEF sat_pos(rx_pos.x + RANGE * dx,
                     rx_pos.y + RANGE * dy,
                     rx_pos.z + RANGE * dz);
        
        sats.emplace_back(prn++, sat_pos, 0.0);
    }
    
    return sats;
}

// Calculate DOP for a configuration
DOP getDOP(const ECEF& rx_pos, const std::vector<SatelliteData>& sats) {
    if (sats.size() < 4) {
        return DOP();  // not enough satellites
    }
    
    size_t n = sats.size();
    Matrix H(n, 4);
    
    for (size_t i = 0; i < n; ++i) {
        ECEF diff = sats[i].position - rx_pos;
        double r = diff.norm();
        if (r < 1e-6) return DOP();
        
        H(i, 0) = -diff.x / r;
        H(i, 1) = -diff.y / r;
        H(i, 2) = -diff.z / r;
        H(i, 3) = 1.0;
    }
    
    try {
        LLA rx_lla = ecefToLLA(rx_pos);
        return computeDOP(H, rx_lla);
    } catch (...) {
        DOP bad;
        bad.gdop = bad.pdop = bad.hdop = bad.vdop = bad.tdop = 99.99;
        return bad;
    }
}

int main() {
    std::cout << "=====================================\n";
    std::cout << "   DOP Analysis\n";
    std::cout << "=====================================\n\n";

    // Pick somewhere in the middle of the US
    LLA rx_lla(40.0 * constants::DEG2RAD, -100.0 * constants::DEG2RAD, 0.0);
    ECEF rx_pos = llaToECEF(rx_lla);

    std::cout << std::fixed << std::setprecision(2);

    // Test 1: Different 4-satellite setups
    std::cout << "Test 1: Four satellites, different arrangements\n";
    std::cout << "------------------------------------------------\n\n";

    // Spread out nicely
    std::vector<std::pair<double, double>> spread = {
        {0, 30}, {90, 50}, {180, 35}, {270, 55}
    };
    
    // All bunched together
    std::vector<std::pair<double, double>> bunched = {
        {0, 40}, {20, 45}, {340, 50}, {10, 35}
    };
    
    // All high up (bad for altitude)
    std::vector<std::pair<double, double>> overhead = {
        {0, 70}, {90, 75}, {180, 65}, {270, 70}
    };

    auto s1 = placeSatellites(rx_pos, spread);
    auto s2 = placeSatellites(rx_pos, bunched);
    auto s3 = placeSatellites(rx_pos, overhead);

    DOP d1 = getDOP(rx_pos, s1);
    DOP d2 = getDOP(rx_pos, s2);
    DOP d3 = getDOP(rx_pos, s3);

    std::cout << "Setup          GDOP   PDOP   HDOP   VDOP\n";
    std::cout << "------------   ----   ----   ----   ----\n";
    std::cout << "Spread out     " << std::setw(4) << d1.gdop << "   "
              << std::setw(4) << d1.pdop << "   " << std::setw(4) << d1.hdop 
              << "   " << std::setw(4) << d1.vdop << "\n";
    std::cout << "Bunched        " << std::setw(4) << d2.gdop << "   "
              << std::setw(4) << d2.pdop << "   " << std::setw(4) << d2.hdop 
              << "   " << std::setw(4) << d2.vdop << "\n";
    std::cout << "Overhead       " << std::setw(4) << d3.gdop << "   "
              << std::setw(4) << d3.pdop << "   " << std::setw(4) << d3.hdop 
              << "   " << std::setw(4) << d3.vdop << "\n\n";

    // Test 2: More satellites = better?
    std::cout << "Test 2: Adding more satellites\n";
    std::cout << "------------------------------\n\n";

    std::cout << "Count   GDOP   PDOP   HDOP   VDOP\n";
    std::cout << "-----   ----   ----   ----   ----\n";

    for (int n = 4; n <= 12; n += 2) {
        std::vector<std::pair<double, double>> cfg;
        
        // Spread them around evenly
        for (int i = 0; i < n; ++i) {
            double az = (360.0 * i) / n;
            double el = 30.0 + 20.0 * (i % 3);  // mix of heights
            cfg.push_back({az, el});
        }
        
        auto sats = placeSatellites(rx_pos, cfg);
        DOP dop = getDOP(rx_pos, sats);
        
        std::cout << "  " << std::setw(2) << n << "    "
                  << std::setw(4) << dop.gdop << "   " << std::setw(4) << dop.pdop 
                  << "   " << std::setw(4) << dop.hdop << "   " 
                  << std::setw(4) << dop.vdop << "\n";
    }
    std::cout << "\n";

    // Test 3: Elevation cutoff
    std::cout << "Test 3: Elevation mask (blocking low satellites)\n";
    std::cout << "------------------------------------------------\n\n";

    std::cout << "Cutoff   GDOP   PDOP   HDOP   VDOP\n";
    std::cout << "------   ----   ----   ----   ----\n";

    for (int cutoff = 5; cutoff <= 45; cutoff += 10) {
        std::vector<std::pair<double, double>> cfg;
        
        for (int i = 0; i < 8; ++i) {
            double az = (360.0 * i) / 8;
            double el = cutoff + 5.0 + 40.0 * (i % 3) / 2.0;
            if (el > 85) el = 85;
            cfg.push_back({az, el});
        }
        
        auto sats = placeSatellites(rx_pos, cfg);
        DOP dop = getDOP(rx_pos, sats);
        
        std::cout << " " << std::setw(2) << cutoff << "°    "
                  << std::setw(4) << dop.gdop << "   " << std::setw(4) << dop.pdop 
                  << "   " << std::setw(4) << dop.hdop << "   " 
                  << std::setw(4) << dop.vdop << "\n";
    }
    std::cout << "\n";

    // Test 4: Urban canyon (buildings block satellites)
    std::cout << "Test 4: Urban canyons\n";
    std::cout << "---------------------\n\n";

    // Open sky
    std::vector<std::pair<double, double>> open = {
        {0, 30}, {45, 45}, {90, 35}, {135, 50},
        {180, 40}, {225, 55}, {270, 45}, {315, 35}
    };

    // N-S street (can only see north/south)
    std::vector<std::pair<double, double>> ns_street = {
        {0, 35}, {15, 50}, {165, 45}, {180, 40}, {195, 55}, {345, 35}
    };

    // E-W street (can only see east/west)
    std::vector<std::pair<double, double>> ew_street = {
        {75, 35}, {90, 50}, {105, 45}, {255, 45}, {270, 40}, {285, 55}
    };

    auto open_sats = placeSatellites(rx_pos, open);
    auto ns_sats = placeSatellites(rx_pos, ns_street);
    auto ew_sats = placeSatellites(rx_pos, ew_street);

    DOP open_dop = getDOP(rx_pos, open_sats);
    DOP ns_dop = getDOP(rx_pos, ns_sats);
    DOP ew_dop = getDOP(rx_pos, ew_sats);

    std::cout << "Situation         GDOP   PDOP   HDOP   VDOP\n";
    std::cout << "---------------   ----   ----   ----   ----\n";
    std::cout << "Open sky (8)      " << std::setw(4) << open_dop.gdop << "   "
              << std::setw(4) << open_dop.pdop << "   " << std::setw(4) << open_dop.hdop 
              << "   " << std::setw(4) << open_dop.vdop << "\n";
    std::cout << "N-S street (6)    " << std::setw(4) << ns_dop.gdop << "   "
              << std::setw(4) << ns_dop.pdop << "   " << std::setw(4) << ns_dop.hdop 
              << "   " << std::setw(4) << ns_dop.vdop << "\n";
    std::cout << "E-W street (6)    " << std::setw(4) << ew_dop.gdop << "   "
              << std::setw(4) << ew_dop.pdop << "   " << std::setw(4) << ew_dop.hdop 
              << "   " << std::setw(4) << ew_dop.vdop << "\n\n";

    std::cout << "=====================================\n";
    std::cout << "   What did we learn?\n";
    std::cout << "=====================================\n\n";
    std::cout << "- Spreading satellites out helps more than just having more\n";
    std::cout << "- High-only satellites give terrible vertical accuracy\n";
    std::cout << "- Urban canyons make accuracy worse in certain directions\n";
    std::cout << "- Low satellites help vertical, but watch out for multipath\n";
    std::cout << "- DOP < 2 is great, 2-5 is fine, > 6 is sketchy\n";

    std::cout << "\n=====================================\n";
    std::cout << "   Done!\n";
    std::cout << "=====================================\n";

    return 0;
}

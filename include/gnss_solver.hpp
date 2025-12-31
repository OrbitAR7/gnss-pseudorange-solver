#ifndef GNSS_SOLVER_HPP
#define GNSS_SOLVER_HPP

#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace gnss {

// ============================================================================
// Constants
// ============================================================================

namespace constants {
    constexpr double C = 299792458.0;           // Speed of light [m/s]
    constexpr double PI = 3.14159265358979323846;
    constexpr double DEG2RAD = PI / 180.0;
    constexpr double RAD2DEG = 180.0 / PI;
    
    // WGS84 Ellipsoid
    constexpr double WGS84_A = 6378137.0;           // Semi-major axis [m]
    constexpr double WGS84_F = 1.0 / 298.257223563; // Flattening
    constexpr double WGS84_B = WGS84_A * (1.0 - WGS84_F);  // Semi-minor axis
    constexpr double WGS84_E2 = 2.0 * WGS84_F - WGS84_F * WGS84_F; // First eccentricity squared
    
    // GPS specific
    constexpr double OMEGA_EARTH = 7.2921151467e-5;  // Earth rotation rate [rad/s]
    constexpr double MU_EARTH = 3.986005e14;         // Earth gravitational parameter [m^3/s^2]
}

// ============================================================================
// Linear Algebra Utilities (minimal, no external deps)
// ============================================================================

/**
 * @brief Simple matrix class for least-squares operations
 */
class Matrix {
public:
    std::vector<double> data;
    size_t rows, cols;

    Matrix() : rows(0), cols(0) {}
    Matrix(size_t r, size_t c) : data(r * c, 0.0), rows(r), cols(c) {}
    Matrix(size_t r, size_t c, double val) : data(r * c, val), rows(r), cols(c) {}

    double& operator()(size_t i, size_t j) { return data[i * cols + j]; }
    double operator()(size_t i, size_t j) const { return data[i * cols + j]; }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j)
                result(j, i) = (*this)(i, j);
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows)
            throw std::runtime_error("Matrix dimension mismatch for multiplication");
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < other.cols; ++j)
                for (size_t k = 0; k < cols; ++k)
                    result(i, j) += (*this)(i, k) * other(k, j);
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::runtime_error("Matrix dimension mismatch for subtraction");
        Matrix result(rows, cols);
        for (size_t i = 0; i < data.size(); ++i)
            result.data[i] = data[i] - other.data[i];
        return result;
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::runtime_error("Matrix dimension mismatch for addition");
        Matrix result(rows, cols);
        for (size_t i = 0; i < data.size(); ++i)
            result.data[i] = data[i] + other.data[i];
        return result;
    }

    // 4x4 matrix inverse using adjugate method (for navigation solution)
    Matrix inverse4x4() const {
        if (rows != 4 || cols != 4)
            throw std::runtime_error("inverse4x4 only works for 4x4 matrices");

        const auto& m = data;
        Matrix inv(4, 4);
        auto& i = inv.data;

        i[0] = m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
        i[4] = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
        i[8] = m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15] + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
        i[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14] - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
        
        i[1] = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
        i[5] = m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
        i[9] = -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15] - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
        i[13] = m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14] + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
        
        i[2] = m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15] + m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
        i[6] = -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
        i[10] = m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
        i[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
        
        i[3] = -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
        i[7] = m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
        i[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11] - m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
        i[15] = m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10] + m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

        double det = m[0]*i[0] + m[1]*i[4] + m[2]*i[8] + m[3]*i[12];
        if (std::abs(det) < 1e-15)
            throw std::runtime_error("Matrix is singular, cannot invert");

        for (auto& val : inv.data)
            val /= det;

        return inv;
    }
};

// ============================================================================
// Position Representations
// ============================================================================

/**
 * @brief ECEF (Earth-Centered Earth-Fixed) position
 */
struct ECEF {
    double x, y, z;  // [m]

    ECEF() : x(0), y(0), z(0) {}
    ECEF(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    
    ECEF operator-(const ECEF& other) const {
        return ECEF(x - other.x, y - other.y, z - other.z);
    }
};

/**
 * @brief Geodetic coordinates (latitude, longitude, height)
 */
struct LLA {
    double lat;  // Latitude [rad]
    double lon;  // Longitude [rad]
    double alt;  // Altitude above ellipsoid [m]

    LLA() : lat(0), lon(0), alt(0) {}
    LLA(double lat_, double lon_, double alt_) : lat(lat_), lon(lon_), alt(alt_) {}

    double latDeg() const { return lat * constants::RAD2DEG; }
    double lonDeg() const { return lon * constants::RAD2DEG; }
};

/**
 * @brief Convert ECEF to LLA (geodetic)
 */
inline LLA ecefToLLA(const ECEF& ecef) {
    const double a = constants::WGS84_A;
    const double e2 = constants::WGS84_E2;
    
    double x = ecef.x, y = ecef.y, z = ecef.z;
    double lon = std::atan2(y, x);
    
    double p = std::sqrt(x*x + y*y);
    double lat = std::atan2(z, p * (1.0 - e2));  // Initial estimate
    
    // Iterative refinement (Bowring's method)
    for (int i = 0; i < 10; ++i) {
        double sin_lat = std::sin(lat);
        double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
        lat = std::atan2(z + e2 * N * sin_lat, p);
    }
    
    double sin_lat = std::sin(lat);
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    double alt = p / std::cos(lat) - N;
    
    return LLA(lat, lon, alt);
}

/**
 * @brief Convert LLA to ECEF
 */
inline ECEF llaToECEF(const LLA& lla) {
    const double a = constants::WGS84_A;
    const double e2 = constants::WGS84_E2;
    
    double sin_lat = std::sin(lla.lat);
    double cos_lat = std::cos(lla.lat);
    double sin_lon = std::sin(lla.lon);
    double cos_lon = std::cos(lla.lon);
    
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    
    double x = (N + lla.alt) * cos_lat * cos_lon;
    double y = (N + lla.alt) * cos_lat * sin_lon;
    double z = (N * (1.0 - e2) + lla.alt) * sin_lat;
    
    return ECEF(x, y, z);
}

// ============================================================================
// Satellite and Measurement Data
// ============================================================================

/**
 * @brief Satellite vehicle data
 */
struct SatelliteData {
    int prn;            // PRN number
    ECEF position;      // Satellite ECEF position [m]
    double clock_bias;  // Satellite clock bias [s]
    
    SatelliteData() : prn(0), clock_bias(0) {}
    SatelliteData(int prn_, const ECEF& pos_, double clk_ = 0.0) 
        : prn(prn_), position(pos_), clock_bias(clk_) {}
};

/**
 * @brief Pseudorange measurement
 */
struct Pseudorange {
    int prn;            // Satellite PRN
    double range;       // Measured pseudorange [m]
    double cn0;         // Carrier-to-noise ratio [dB-Hz] (optional)
    double weight;      // Measurement weight (1/variance)
    
    Pseudorange() : prn(0), range(0), cn0(45), weight(1.0) {}
    Pseudorange(int prn_, double range_, double cn0_ = 45.0) 
        : prn(prn_), range(range_), cn0(cn0_), weight(1.0) {}
};

// ============================================================================
// DOP (Dilution of Precision)
// ============================================================================

/**
 * @brief Dilution of Precision values
 */
struct DOP {
    double gdop;  // Geometric DOP
    double pdop;  // Position DOP
    double hdop;  // Horizontal DOP
    double vdop;  // Vertical DOP
    double tdop;  // Time DOP

    DOP() : gdop(0), pdop(0), hdop(0), vdop(0), tdop(0) {}
};

/**
 * @brief Compute DOP from geometry matrix
 */
inline DOP computeDOP(const Matrix& H, const LLA& receiver_lla) {
    DOP dop;
    
    // Q = (H^T * H)^-1
    Matrix HtH = H.transpose() * H;
    Matrix Q = HtH.inverse4x4();
    
    // GDOP = sqrt(trace(Q))
    dop.gdop = std::sqrt(Q(0,0) + Q(1,1) + Q(2,2) + Q(3,3));
    dop.pdop = std::sqrt(Q(0,0) + Q(1,1) + Q(2,2));
    dop.tdop = std::sqrt(Q(3,3));
    
    // Transform to local ENU for HDOP/VDOP
    double sin_lat = std::sin(receiver_lla.lat);
    double cos_lat = std::cos(receiver_lla.lat);
    double sin_lon = std::sin(receiver_lla.lon);
    double cos_lon = std::cos(receiver_lla.lon);
    
    // Rotation matrix ECEF -> ENU
    Matrix R(3, 3);
    R(0,0) = -sin_lon;           R(0,1) = cos_lon;            R(0,2) = 0;
    R(1,0) = -sin_lat*cos_lon;   R(1,1) = -sin_lat*sin_lon;   R(1,2) = cos_lat;
    R(2,0) = cos_lat*cos_lon;    R(2,1) = cos_lat*sin_lon;    R(2,2) = sin_lat;
    
    // Q_pos = 3x3 position submatrix
    Matrix Q_pos(3, 3);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            Q_pos(i, j) = Q(i, j);
    
    // Q_enu = R * Q_pos * R^T
    Matrix Q_enu = R * Q_pos * R.transpose();
    
    dop.hdop = std::sqrt(Q_enu(0,0) + Q_enu(1,1));  // East + North
    dop.vdop = std::sqrt(Q_enu(2,2));               // Up
    
    return dop;
}

// ============================================================================
// Navigation Solution
// ============================================================================

/**
 * @brief Navigation solution result
 */
struct NavSolution {
    ECEF position;          // Estimated receiver position [m]
    double clock_bias;      // Receiver clock bias [m] (c * dt)
    LLA lla;                // Geodetic position
    DOP dop;                // Dilution of precision
    int num_iterations;     // Iterations to converge
    double residual_rms;    // RMS of post-fit residuals [m]
    bool valid;             // Solution validity flag
    std::vector<double> residuals;  // Post-fit residuals per satellite

    NavSolution() : clock_bias(0), num_iterations(0), residual_rms(0), valid(false) {}
};

/**
 * @brief Configuration for the solver
 */
struct SolverConfig {
    double convergence_threshold = 1e-4;  // Position change threshold [m]
    int max_iterations = 20;
    bool apply_earth_rotation = true;     // Sagnac correction
    bool weighted_least_squares = false;  // Use measurement weights
    ECEF initial_position = ECEF(0, 0, 0); // Initial guess (0,0,0 = center of Earth)
};

// ============================================================================
// GNSS Position Solver
// ============================================================================

/**
 * @brief Iterative least-squares GNSS position solver
 * 
 * Solves the navigation equations using Gauss-Newton iteration:
 *   rho_i = ||r_sat_i - r_rx|| + c*dt_rx - c*dt_sat_i + corrections
 */
class GNSSSolver {
public:
    SolverConfig config;

    GNSSSolver() = default;
    explicit GNSSSolver(const SolverConfig& cfg) : config(cfg) {}

    /**
     * @brief Solve for receiver position and clock bias
     * @param satellites Vector of satellite positions/clock data
     * @param pseudoranges Vector of pseudorange measurements
     * @return Navigation solution
     */
    NavSolution solve(const std::vector<SatelliteData>& satellites,
                      const std::vector<Pseudorange>& pseudoranges) const {
        
        NavSolution solution;
        
        // Need at least 4 satellites
        if (satellites.size() < 4 || pseudoranges.size() < 4) {
            solution.valid = false;
            return solution;
        }

        // Match satellites with pseudoranges
        std::vector<std::pair<SatelliteData, Pseudorange>> matched;
        for (const auto& pr : pseudoranges) {
            for (const auto& sat : satellites) {
                if (sat.prn == pr.prn) {
                    matched.push_back({sat, pr});
                    break;
                }
            }
        }

        if (matched.size() < 4) {
            solution.valid = false;
            return solution;
        }

        size_t n = matched.size();

        // State vector: [x, y, z, cdt]
        double x = config.initial_position.x;
        double y = config.initial_position.y;
        double z = config.initial_position.z;
        double cdt = 0.0;  // Receiver clock bias in meters

        // Gauss-Newton iteration
        for (int iter = 0; iter < config.max_iterations; ++iter) {
            Matrix H(n, 4);      // Geometry matrix
            Matrix dz(n, 1);     // Measurement residuals
            Matrix W(n, n);      // Weight matrix

            for (size_t i = 0; i < n; ++i) {
                const auto& sat = matched[i].first;
                const auto& pr = matched[i].second;

                // Satellite position (apply Earth rotation correction if enabled)
                ECEF sat_pos = sat.position;
                if (config.apply_earth_rotation) {
                    double tau = pr.range / constants::C;
                    double omega_tau = constants::OMEGA_EARTH * tau;
                    double cos_ot = std::cos(omega_tau);
                    double sin_ot = std::sin(omega_tau);
                    sat_pos.x = cos_ot * sat.position.x + sin_ot * sat.position.y;
                    sat_pos.y = -sin_ot * sat.position.x + cos_ot * sat.position.y;
                    sat_pos.z = sat.position.z;
                }

                // Geometric range
                double dx = sat_pos.x - x;
                double dy = sat_pos.y - y;
                double dz_pos = sat_pos.z - z;
                double rho = std::sqrt(dx*dx + dy*dy + dz_pos*dz_pos);

                // Predicted pseudorange
                double rho_pred = rho + cdt - constants::C * sat.clock_bias;

                // Measurement residual
                dz(i, 0) = pr.range - rho_pred;

                // Geometry matrix row (partial derivatives)
                H(i, 0) = -dx / rho;  // d(rho)/dx
                H(i, 1) = -dy / rho;  // d(rho)/dy
                H(i, 2) = -dz_pos / rho;  // d(rho)/dz
                H(i, 3) = 1.0;        // d(rho)/d(cdt)

                // Weight (inverse variance)
                W(i, i) = config.weighted_least_squares ? pr.weight : 1.0;
            }

            // Weighted least squares: dx = (H^T W H)^-1 H^T W dz
            Matrix Ht = H.transpose();
            Matrix HtWH = Ht * W * H;
            Matrix HtWdz = Ht * W * dz;
            
            Matrix HtWH_inv = HtWH.inverse4x4();
            Matrix delta = HtWH_inv * HtWdz;

            // Update state
            x += delta(0, 0);
            y += delta(1, 0);
            z += delta(2, 0);
            cdt += delta(3, 0);

            solution.num_iterations = iter + 1;

            // Check convergence
            double pos_change = std::sqrt(delta(0,0)*delta(0,0) + 
                                         delta(1,0)*delta(1,0) + 
                                         delta(2,0)*delta(2,0));
            if (pos_change < config.convergence_threshold) {
                break;
            }
        }

        // Store solution
        solution.position = ECEF(x, y, z);
        solution.clock_bias = cdt;
        solution.lla = ecefToLLA(solution.position);

        // Compute post-fit residuals
        solution.residuals.resize(n);
        double sum_sq = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const auto& sat = matched[i].first;
            const auto& pr = matched[i].second;
            
            ECEF sat_pos = sat.position;
            if (config.apply_earth_rotation) {
                double tau = pr.range / constants::C;
                double omega_tau = constants::OMEGA_EARTH * tau;
                sat_pos.x = std::cos(omega_tau) * sat.position.x + std::sin(omega_tau) * sat.position.y;
                sat_pos.y = -std::sin(omega_tau) * sat.position.x + std::cos(omega_tau) * sat.position.y;
            }
            
            double rho = (sat_pos - solution.position).norm();
            double rho_pred = rho + cdt - constants::C * sat.clock_bias;
            solution.residuals[i] = pr.range - rho_pred;
            sum_sq += solution.residuals[i] * solution.residuals[i];
        }
        solution.residual_rms = std::sqrt(sum_sq / n);

        // Compute DOP
        Matrix H_final(n, 4);
        for (size_t i = 0; i < n; ++i) {
            const auto& sat = matched[i].first;
            double dx = sat.position.x - x;
            double dy = sat.position.y - y;
            double dz_pos = sat.position.z - z;
            double rho = std::sqrt(dx*dx + dy*dy + dz_pos*dz_pos);
            H_final(i, 0) = -dx / rho;
            H_final(i, 1) = -dy / rho;
            H_final(i, 2) = -dz_pos / rho;
            H_final(i, 3) = 1.0;
        }
        solution.dop = computeDOP(H_final, solution.lla);
        solution.valid = true;

        return solution;
    }
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Generate simulated pseudoranges from true position
 */
inline std::vector<Pseudorange> simulatePseudoranges(
    const ECEF& true_position,
    double true_clock_bias,
    const std::vector<SatelliteData>& satellites,
    double noise_std = 1.0,
    unsigned int seed = 42) {
    
    std::vector<Pseudorange> measurements;
    
    // Simple pseudo-random noise (deterministic for reproducibility)
    auto noise = [&seed]() {
        seed = seed * 1103515245 + 12345;
        double u1 = (seed % 10000) / 10000.0;
        seed = seed * 1103515245 + 12345;
        double u2 = (seed % 10000) / 10000.0;
        return std::sqrt(-2.0 * std::log(u1 + 1e-10)) * std::cos(2.0 * constants::PI * u2);
    };
    
    for (const auto& sat : satellites) {
        // Check if satellite is above horizon (elevation > 5°)
        ECEF diff = sat.position - true_position;
        double range = diff.norm();
        
        LLA user_lla = ecefToLLA(true_position);
        double sin_lat = std::sin(user_lla.lat);
        double cos_lat = std::cos(user_lla.lat);
        double sin_lon = std::sin(user_lla.lon);
        double cos_lon = std::cos(user_lla.lon);
        
        // Up vector in ECEF
        double up_x = cos_lat * cos_lon;
        double up_y = cos_lat * sin_lon;
        double up_z = sin_lat;
        
        // Elevation angle
        double elevation = std::asin((diff.x * up_x + diff.y * up_y + diff.z * up_z) / range);
        
        if (elevation > 5.0 * constants::DEG2RAD) {
            double true_range = range + true_clock_bias - constants::C * sat.clock_bias;
            double measured = true_range + noise_std * noise();
            measurements.emplace_back(sat.prn, measured, 45.0);
        }
    }
    
    return measurements;
}

/**
 * @brief Print navigation solution
 */
inline void printSolution(const NavSolution& sol, const std::string& label = "") {
    if (!label.empty()) {
        std::cout << "=== " << label << " ===" << std::endl;
    }
    
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Position (ECEF): [" << sol.position.x << ", " 
              << sol.position.y << ", " << sol.position.z << "] m" << std::endl;
    std::cout << std::setprecision(8);
    std::cout << "Position (LLA):  " << sol.lla.latDeg() << "° N, " 
              << sol.lla.lonDeg() << "° E, " 
              << std::setprecision(2) << sol.lla.alt << " m" << std::endl;
    std::cout << std::setprecision(3);
    std::cout << "Clock bias:      " << sol.clock_bias << " m (" 
              << sol.clock_bias / constants::C * 1e9 << " ns)" << std::endl;
    std::cout << "Iterations:      " << sol.num_iterations << std::endl;
    std::cout << "Residual RMS:    " << sol.residual_rms << " m" << std::endl;
    std::cout << std::setprecision(2);
    std::cout << "DOP: GDOP=" << sol.dop.gdop << " PDOP=" << sol.dop.pdop 
              << " HDOP=" << sol.dop.hdop << " VDOP=" << sol.dop.vdop << std::endl;
}

} // namespace gnss

#endif // GNSS_SOLVER_HPP

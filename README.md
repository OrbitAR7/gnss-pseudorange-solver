# GNSS Pseudorange Solver

A simple C++ library for computing GPS/GNSS positions from pseudorange measurements. It's header-only (just one file to include) and uses iterative least squares to figure out where you are. I built this mainly for learning and tinkering with GPS concepts.

## What it does

- **Position solving** — Uses Gauss-Newton iteration to estimate position and clock bias
- **Weighted least squares** — Can weight measurements if some are better than others
- **DOP values** — Computes GDOP, PDOP, HDOP, VDOP, TDOP
- **Coordinate math** — Converts between ECEF and lat/lon/altitude (WGS84)
- **Earth rotation** — Accounts for Sagnac effect (earth spins while signal travels)
- **Testing tools** — Generate fake measurements for experiments
- **Header-only** — Just include one file, no dependencies

## Quick example

```cpp
#include "gnss_solver.hpp"

using namespace gnss;

int main() {
    // Set up some satellites (normally from ephemeris data)
    std::vector<SatelliteData> satellites;
    satellites.emplace_back(1, ECEF(15600000, 7540000, 20140000), 0.0);
    satellites.emplace_back(2, ECEF(18760000, 2750000, 18610000), 0.0);
    satellites.emplace_back(3, ECEF(17610000, 14630000, 13480000), 0.0);
    satellites.emplace_back(4, ECEF(19170000, 610000, 18390000), 0.0);
    satellites.emplace_back(5, ECEF(14760000, 14270000, 15450000), 0.0);

    // Measurements from receiver
    std::vector<Pseudorange> measurements;
    measurements.emplace_back(1, 20200000.0);
    measurements.emplace_back(2, 21300000.0);
    measurements.emplace_back(3, 22100000.0);
    measurements.emplace_back(4, 21500000.0);
    measurements.emplace_back(5, 20800000.0);

    // Solve it
    GNSSSolver solver;
    NavSolution solution = solver.solve(satellites, measurements);

    if (solution.valid) {
        std::cout << "Position: " << solution.lla.latDeg() << "°N, "
                  << solution.lla.lonDeg() << "°E\n";
        std::cout << "HDOP: " << solution.dop.hdop << "\n";
    }

    return 0;
}
```

## Building

### What you need
- C++17 compiler
- CMake 3.14+ (optional)

### Using CMake

```bash
mkdir build && cd build
cmake ..
make

# Try the examples
./gps_position_fix
./dop_analysis

# Run tests
./test_gnss_solver
```

### Without CMake

```bash
g++ -std=c++17 -O2 -I include examples/gps_position_fix.cpp -o gps_position_fix
```

## Examples

### GPS Position Fix (`examples/gps_position_fix.cpp`)

Shows a complete position fix:
- Makes a fake 24-satellite GPS constellation
- Generates measurements with realistic noise
- Solves for position
- Runs 100 trials to check consistency

### DOP Analysis (`examples/dop_analysis.cpp`)

Plays with satellite geometry to see how it affects accuracy:
- Good vs bad satellite arrangements
- More satellites vs fewer satellites
- Urban canyon effects (buildings blocking satellites)
- Elevation mask tradeoffs

## API

### Main types

**ECEF** - Earth-centered coordinates (meters)
```cpp
struct ECEF {
    double x, y, z;
};
```

**LLA** - Latitude/longitude/altitude
```cpp
struct LLA {
    double lat;  // radians
    double lon;  // radians
    double alt;  // meters
};
```

**SatelliteData** - Info about one satellite
```cpp
struct SatelliteData {
    int prn;            // which satellite
    ECEF position;      // where it is
    double clock_bias;  // its clock offset (seconds)
};
```

**Pseudorange** - One measurement
```cpp
struct Pseudorange {
    int prn;        // which satellite
    double range;   // measured range (meters)
    double weight;  // for weighted least squares
};
```

**NavSolution** - The answer
```cpp
struct NavSolution {
    ECEF position;
    double clock_bias;      // meters
    LLA lla;
    DOP dop;
    int num_iterations;
    double residual_rms;
    bool valid;
};
```

**DOP** - Accuracy indicators
```cpp
struct DOP {
    double gdop;  // geometric
    double pdop;  // position
    double hdop;  // horizontal
    double vdop;  // vertical
    double tdop;  // time
};
```

### Configuration

```cpp
struct SolverConfig {
    double convergence_threshold = 1e-4;  // meters
    int max_iterations = 20;
    bool apply_earth_rotation = true;
    bool weighted_least_squares = false;
    ECEF initial_position;  // starting guess
};
```

### Coordinate conversion

```cpp
LLA ecefToLLA(const ECEF& ecef);
ECEF llaToECEF(const LLA& lla);
```

### Helper functions

```cpp
// Make fake measurements for testing
std::vector<Pseudorange> simulatePseudoranges(
    const ECEF& true_position,
    double true_clock_bias,
    const std::vector<SatelliteData>& satellites,
    double noise_std = 1.0);

// Print solution nicely
void printSolution(const NavSolution& sol, const std::string& label = "");
```

## How it works

### The math

Pseudorange equation:
```
ρ = ||sat_pos - rx_pos|| + c·clock_bias + noise
```

We linearize and solve iteratively:
```
δx = (H'WH)^(-1) H'W δz
```

Where:
- H = geometry matrix (direction cosines)
- W = weight matrix
- δz = measurement residuals

### DOP calculation

Derived from covariance matrix:
```
Q = (H'H)^(-1)
GDOP = √(Q₁₁ + Q₂₂ + Q₃₃ + Q₄₄)
PDOP = √(Q₁₁ + Q₂₂ + Q₃₃)
```

HDOP/VDOP need rotation to local coordinates.

### Earth rotation correction

While the signal travels (~70ms), Earth rotates a bit:
```
angle = ω·τ  where τ = range/c
```

Rotate satellite position by this angle.

## DOP cheat sheet

| DOP | Quality | Notes |
|-----|---------|-------|
| < 1 | Perfect | Rarely happens |
| 1-2 | Great | Good geometry |
| 2-5 | Fine | Most situations |
| 5-10 | Meh | Use with caution |
| > 10 | Bad | Don't trust it |

## What's missing

This is for learning/prototyping. Real GPS receivers also do:

- Tropospheric correction (signal delay through atmosphere)
- Ionospheric correction (charged particles slow signal)
- Precise satellite clocks (interpolation between samples)
- Carrier phase (way more accurate than pseudorange)
- RAIM (detecting bad satellites)
- Multiple constellations (GLONASS, Galileo, BeiDou)

## References

1. Kaplan & Hegarty - *Understanding GPS/GNSS* (the standard textbook)
2. Misra & Enge - *Global Positioning System* (another good one)
3. IS-GPS-200 (official GPS spec)

## License

MIT - do whatever you want with it

## About

Made by Hamoud Alshammari  
[GitHub](https://github.com/OrbitAR7) | [LinkedIn](https://www.linkedin.com/in/hamoud-alshammari-7415b617a/)

---

*If you find this useful or have questions, feel free to open an issue!*

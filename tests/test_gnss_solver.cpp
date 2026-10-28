// Validation tests for gnss_solver.hpp. Returns non-zero if any check fails.
//  1. Noise-free pseudoranges, no Earth rotation: exact recovery of position and clock.
//  2. Sagnac-consistent pseudoranges: solver with Earth-rotation correction recovers truth.
//  3. 1 m white noise, 500 Monte Carlo runs: 3-D RMS error ~ PDOP * sigma.
#include "gnss_solver.hpp"
#include <cstdio>
#include <cmath>
using namespace gnss;
static int failures = 0;
static void check(bool ok, const char* what) { std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what); if (!ok) ++failures; }
int main(){
  LLA rx_lla; rx_lla.lat=24.7*constants::DEG2RAD; rx_lla.lon=46.7*constants::DEG2RAD; rx_lla.alt=600;
  ECEF rx=llaToECEF(rx_lla); const double b=1500.0;   // receiver clock bias [m]
  // 8 satellites at GPS radius spread in az/el above the receiver
  double lat=rx_lla.lat, lon=rx_lla.lon; double R=26560e3;
  double E[3]={-sin(lon),cos(lon),0}, N[3]={-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat)}, U[3]={cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)};
  std::vector<SatelliteData> sats; double az[8]={10,60,110,170,220,280,330,0}, el[8]={15,40,25,70,20,50,30,85};
  for(int i=0;i<8;i++){ double a=az[i]*constants::DEG2RAD,e=el[i]*constants::DEG2RAD; double u[3]; for(int k=0;k<3;k++) u[k]=cos(e)*sin(a)*E[k]+cos(e)*cos(a)*N[k]+sin(e)*U[k];
    double px=rx.x,py=rx.y,pz=rx.z, B=2*(px*u[0]+py*u[1]+pz*u[2]), C=px*px+py*py+pz*pz-R*R; double t=(-B+sqrt(B*B-4*C))/2;
    sats.emplace_back(i+1, ECEF(px+t*u[0],py+t*u[1],pz+t*u[2]), 0.0); }
  // Test 1: no Earth rotation, noise-free
  { std::vector<Pseudorange> m; for(auto&s:sats){ m.emplace_back(s.prn,(s.position-rx).norm()+b);} SolverConfig c; c.apply_earth_rotation=false; GNSSSolver g(c); auto sol=g.solve(sats,m);
    printf("no-rotation, noise-free: |dr|=%.2e m, clock err=%.2e m, iters=%d\n",(sol.position-rx).norm(),sol.clock_bias-b,sol.num_iterations);
    check(sol.valid && (sol.position-rx).norm()<1e-4 && std::fabs(sol.clock_bias-b)<1e-4, "exact recovery without Earth rotation"); }
  // Test 2: Sagnac-consistent truth (satellite positions given in ECEF at transmit time)
  { std::vector<Pseudorange> m; for(auto&s:sats){ double tau=0.07; for(int k=0;k<10;k++){ double ot=constants::OMEGA_EARTH*tau; ECEF p(cos(ot)*s.position.x+sin(ot)*s.position.y,-sin(ot)*s.position.x+cos(ot)*s.position.y,s.position.z); tau=(p-rx).norm()/constants::C;} 
      double ot=constants::OMEGA_EARTH*tau; ECEF p(cos(ot)*s.position.x+sin(ot)*s.position.y,-sin(ot)*s.position.x+cos(ot)*s.position.y,s.position.z); m.emplace_back(s.prn,(p-rx).norm()+b);} 
    GNSSSolver g; auto sol=g.solve(sats,m); printf("with Sagnac, noise-free: |dr|=%.3f m, clock err=%.3f m, PDOP=%.2f\n",(sol.position-rx).norm(),sol.clock_bias-b,sol.dop.pdop);
    check((sol.position-rx).norm()<0.05, "Sagnac-corrected solution within 5 cm");
    SolverConfig c; c.apply_earth_rotation=false; GNSSSolver g2(c); auto s2=g2.solve(sats,m); printf("same data, Sagnac correction OFF: |dr|=%.1f m (size of the effect)\n",(s2.position-rx).norm());
    check((s2.position-rx).norm()>5.0, "ignoring Earth rotation produces a large error"); }
  // Test 3: 1 m noise Monte Carlo, compare RMS error with PDOP*sigma
  { double se=0; int N=500; double pdop=0; for(int t=0;t<N;t++){ auto m=simulatePseudoranges(rx,b,sats,1.0,1000+t); SolverConfig c; c.apply_earth_rotation=false; GNSSSolver g(c); auto sol=g.solve(sats,m); se+=pow((sol.position-rx).norm(),2); pdop=sol.dop.pdop;} 
    printf("1 m noise, %d runs: 3-D RMS error %.2f m vs PDOP*sigma %.2f m\n",N,sqrt(se/N),pdop);
    check(std::fabs(sqrt(se/N)-pdop)/pdop<0.15, "Monte Carlo RMS error within 15% of PDOP*sigma"); }
  std::printf("%s\n", failures ? "SOME TESTS FAILED" : "All tests passed");
  return failures ? 1 : 0;
}

#include "kernel.h"

inline real W(const realvec dr, const real dr2){
  constexpr real H = slen / 2.0;
#if DAM_2D
  constexpr real COEF = 10.0 / 7.0 / pi / pow(H, 2);
#else
  constexpr real COEF = 1.0 / pi / pow(H, 3);
#endif
  const real s = sqrt(dr2) / H;
  real v;
  if (s < 1.0) {
    v = 1.0 - 1.5 * pow(s, 2) + 0.75 * pow(s, 3);
  } else if (s < 2.0) {
    v = 0.25 * pow(2.0 - s, 3);
  }
  return COEF * v;
}

inline realvec gradW(const realvec dr, const real dr2){
  constexpr real H = slen / 2.0;
#if DAM_2D
  constexpr real COEF = 45.0 / 14.0 / pi / pow(H, 4);
#else
  constexpr real COEF = 2.25 / pi / pow(H, 5);
#endif
  const real s = sqrt(dr2) / H;
  realvec v;
  if (s < 1.0) {
    v = (s - 4.0 / 3.0) * dr;
  } else if (s < 2.0) {
    v = - pow(2.0 - s, 2) / 3.0 / s * dr;
  }
  return COEF * v;
}

/* Force Functors */
void CalcDensity::operator () (const EP* const ep_i, const PS::S32 Nip,
    const EP* const ep_j, const PS::S32 Njp, Dens* const dens){
  asm volatile("# dens begins");
  constexpr real slen2 = slen * slen;
  for(PS::S32 i = 0 ; i < Nip ; ++i){
    dens[i].clear();
    for(PS::S32 j = 0 ; j < Njp ; ++j){
      const realvec dr  = ep_i[i].pos - ep_j[j].pos;
      const real    dr2 = dr * dr;
      if (dr2 >= slen2) continue;
      const real W_ij = W(dr, dr2);
      dens[i].dens += ep_j[j].mass * W_ij;
    }
    dens[i].pres = calcPressure(dens[i].dens);
  }
  asm volatile("# dens ends");
}

void CalcHydroForce::operator () (const EP* const ep_i, const PS::S32 Nip,
    const EP* const ep_j, const PS::S32 Njp, Hydro* const hydro){
  asm volatile("# hydro begins");
  constexpr real slen2 = slen * slen;
  real tmp_pd_j[Njp];
  for(PS::S32 j = 0; j < Njp ; ++j){
    tmp_pd_j[j] = ep_j[j].pres / pow(ep_j[j].dens, 2);
  }
  for(PS::S32 i = 0; i < Nip ; ++ i){
    hydro[i].clear();
    /* if (ep_i[i].type != FLUID) continue; */
    const real tmp_pd_i = ep_i[i].pres / pow(ep_i[i].dens, 2);
    for(PS::S32 j = 0; j < Njp ; ++j){
      const realvec dr  = ep_i[i].pos - ep_j[j].pos;
      const real    dr2 = dr * dr;
      if (dr2 >= slen2) continue;
      const realvec gradW_ij = gradW(dr, dr2);
      const realvec dv = ep_i[i].vel - ep_j[j].vel;
      const real vr = dv * dr;
      const real AV = (vr <= 0) ? 0 : - VISC * vr / (dr2 + 0.01 * slen2);
      hydro[i].acc -= ep_j[j].mass * (tmp_pd_i + tmp_pd_j[j] + AV) * gradW_ij;
    }
    hydro[i].acc += gravity;
#if CFL_DT
    hydro[i].f = ep_i[i].mass * sqrt(hydro[i].acc * hydro[i].acc);
#endif
  }
  asm volatile("# hydro ends");
}

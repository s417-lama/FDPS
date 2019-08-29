#pragma once

#define DAM_2D 1
#define LARGE  0
#define CFL_DT 0
#define REUSE  1
#define OUTPUT 0
#define DOUBLE 1

#define DISABLE_STEAL 0

#define MAX_STEP 1000

#if DOUBLE
#define real    PS::F64
#define realvec PS::F64vec
#else
#define real    PS::F32
#define realvec PS::F32vec
#endif

#if DAM_2D
#define PARTICLE_SIMULATOR_TWO_DIMENSION 1
#endif

// Include FDPS header
#include <ps_defs.hpp>
namespace PS = ParticleSimulator;

// Include the standard C++ headers
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <sys/time.h>

enum { FLUID = 1, WALL = 2 };

/* Parameters */
#if OUTPUT
constexpr PS::U32 OUTPUT_INTERVAL = 100;
#endif
constexpr real END_TIME        = 1.5;
#if DAM_2D
#if LARGE
constexpr real l0              = 0.55 / 30 / 8;
#else
constexpr real l0              = 0.55 / 30 / 4;
#endif
#else
constexpr real l0              = 0.55 / 30;
#endif
constexpr real slen            = l0 * 2.1;
#if REUSE
constexpr real skin            = slen * 0.3;
#else
constexpr real skin            = 0.0;
#endif
constexpr real DENS0           = 1000.0;
constexpr real SND             = 31.3;
constexpr real C_B             = DENS0 * SND * SND / 7;
constexpr real ALPHA           = 0.1;
constexpr real VISC            = ALPHA * slen * SND / DENS0;
constexpr real DT              = 0.4 * slen / SND / (1 + 0.6 * ALPHA);
/* const real DT              = 0.0005; */
constexpr real pi              = atan(1.0) * 4.0;
#if DAM_2D
const realvec gravity(0.0, -9.81);
#else
const realvec gravity(0.0, 0.0, -9.81);
#endif

inline PS::U64 get_time_in_usec() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec * 1000 * 1000 + tv.tv_usec;
}

/* Class Definitions */
//** Force Class (Result Class)
class Dens{
  public:
    real dens;
    real pres;
    void clear(){
      dens = 0;
      pres = 0;
    }
};
class Hydro{
  public:
    realvec acc;
#if CFL_DT
    real f;
#endif
    inline void clear(){
      acc = 0;
#if CFL_DT
      f = 0.0;
#endif
    }
};

//** Full Particle Class
struct FP{
  real mass;
  realvec pos;
#if REUSE
  realvec prev_pos;
#endif
  realvec vel;
  realvec acc;
  real dens;
  real pres;
  realvec vel_half;
#if CFL_DT
  real f;
#endif
  PS::S32 type;
  inline void copyFromForce(const Dens& dens){
    this->dens = dens.dens;
    this->pres = dens.pres;
  }
  inline void copyFromForce(const Hydro& force){
    this->acc = force.acc;
#if CFL_DT
    this->f = force.f;
#endif
  }
  inline real getCharge() const{
    return this->mass;
  }
  inline realvec getPos() const{
    return this->pos;
  }
  inline real getRSearch() const{
    return slen + skin;
  }
  inline void setPos(const realvec& pos){
    this->pos = pos;
  }
  inline void writeAscii(FILE* fp) const{
#if DAM_2D
    fprintf(fp, "%f %f %d\n", this->pos.x, this->pos.y, this->type);
#else
    fprintf(fp, "%f %f %f %d\n", this->pos.x, this->pos.y, this->pos.z, this->type);
#endif
  }
  inline void readAscii(FILE* fp){
#if DAM_2D
    fscanf(fp,
#if DOUBLE
          "%lf %lf %d\n",
#else
          "%f %f %d\n",
#endif
          &this->pos.x, &this->pos.y, &this->type);
#else
    fscanf(fp,
#if DOUBLE
           "%lf %lf %lf %d\n",
#else
           "%f %f %f %d\n",
#endif
           &this->pos.x, &this->pos.y, &this->pos.z, &this->type);
#endif
  }
};

//** Essential Particle Class
struct EP{
  realvec pos;
  realvec vel;
  real    mass;
  real    dens;
  real    pres;
  PS::S32 type;
  inline void copyFromFP(const FP& rp){
    this->pos  = rp.pos;
    this->vel  = rp.vel;
    this->mass = rp.mass;
    this->dens = rp.dens;
    this->pres = rp.pres;
    this->type = rp.type;
  }
  inline realvec getPos() const{
    return this->pos;
  }
  inline real getRSearch() const{
    return slen + skin;
  }
  inline void setPos(const realvec& pos){
    this->pos = pos;
  }
};

class CalcDensity{
  public:
    void operator () (const EP* const ep_i, const PS::S32 Nip,
        const EP* const ep_j, const PS::S32 Njp, Dens* const dens);
};

class CalcHydroForce{
  public:
    void operator () (const EP* const ep_i, const PS::S32 Nip,
        const EP* const ep_j, const PS::S32 Njp, Hydro* const hydro);
};

inline real calcPressure(const real dens) {
  return std::max(0.0, C_B * (pow(dens / DENS0, 7) - 1));
}

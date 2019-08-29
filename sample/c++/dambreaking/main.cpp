#include "kernel.h"

// Include FDPS header
#include <particle_simulator.hpp>

void makeOutputDirectory(char * dir_name) {
  struct stat st;
  PS::S32 ret;
  if (PS::Comm::getRank() == 0) {
    if (stat(dir_name, &st) != 0) {
      ret = mkdir(dir_name, 0777);
    } else {
      ret = 0; // the directory named dir_name already exists.
    }
  }
  PS::Comm::broadcast(&ret, 1);
  if (ret == 0) {
    if (PS::Comm::getRank() == 0)
      fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
  } else {
    if (PS::Comm::getRank() == 0)
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
    PS::Abort();
  }
}

void SetupIC(PS::ParticleSystem<FP>& sph_system){
  // Read SPH particles
#if DAM_2D
#if LARGE
  sph_system.readParticleAscii("data/data2d_large.txt");
#else
  sph_system.readParticleAscii("data/data2d.txt");
#endif
#else
  sph_system.readParticleAscii("data/data3d.txt");
#endif
  for(PS::U32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    sph_system[i].dens = DENS0;
    sph_system[i].pres = calcPressure(sph_system[i].dens);
#if DAM_2D
    sph_system[i].mass = DENS0 * pow(l0, 2);
#else
    sph_system[i].mass = DENS0 * pow(l0, 3);
#endif
  }
  std::cout << "# of ptcls is... " << sph_system.getNumberOfParticleLocal() << std::endl;
  // Fin.
  std::cout << "setup..." << std::endl;
}

PS::F64 getTimeStep(const PS::ParticleSystem<FP>& sph_system){
#if CFL_DT
   real fmax = 0.0;
   for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
     fmax = std::max(fmax, sph_system[i].f);
   }
   if (fmax == 0.0) {
     return DT;
   } else {
     return std::min(0.25 * slen / fmax, DT);
   }
#else
   return DT;
#endif
}

void InitialKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
#ifdef PARTICLE_SIMULATOR_TASK_PARALLEL
  mtbb::parallel_for(0, sph_system.getNumberOfParticleLocal(), 1, CUTOFF_PFOR, [&] (PS::S32 a, PS::S32 b) {
    for(PS::S32 i = a ; i < b ; ++ i){
      if (sph_system[i].type == FLUID) {
        sph_system[i].vel_half = sph_system[i].vel + 0.5 * dt * sph_system[i].acc;
      }
    }
  });
#else
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    if (sph_system[i].type == FLUID) {
      sph_system[i].vel_half = sph_system[i].vel + 0.5 * dt * sph_system[i].acc;
    }
  }
#endif
}

#if REUSE
bool FullDrift(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
  bool reuse = true;
  // time becomes t + dt;
#ifdef PARTICLE_SIMULATOR_TASK_PARALLEL
  mtbb::parallel_for(0, sph_system.getNumberOfParticleLocal(), 1, CUTOFF_PFOR, [&] (PS::S32 a, PS::S32 b) {
    for(PS::S32 i = a; i < b; ++ i){
      if (sph_system[i].type == FLUID) {
        sph_system[i].pos += dt * sph_system[i].vel_half;
        // check whether we should reuse the list
        realvec dp = sph_system[i].pos - sph_system[i].prev_pos;
        if (sqrt(dp * dp) >= skin / 2.0) reuse = false;
      }
    }
  });
#else
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    if (sph_system[i].type == FLUID) {
      sph_system[i].pos += dt * sph_system[i].vel_half;
      // check whether we should reuse the list
      realvec dp = sph_system[i].pos - sph_system[i].prev_pos;
      if (sqrt(dp * dp) >= skin / 2.0) reuse = false;
    }
  }
#endif
  return reuse;
}
#else
void FullDrift(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
  // time becomes t + dt;
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    if (sph_system[i].type == FLUID) {
      sph_system[i].pos += dt * sph_system[i].vel_half;
    }
  }
}
#endif

void Predict(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    if (sph_system[i].type == FLUID) {
      sph_system[i].vel += dt * sph_system[i].acc;
    }
  }
}

void FinalKick(PS::ParticleSystem<FP>& sph_system, const PS::F64 dt){
#ifdef PARTICLE_SIMULATOR_TASK_PARALLEL
  mtbb::parallel_for(0, sph_system.getNumberOfParticleLocal(), 1, CUTOFF_PFOR, [&] (PS::S32 a, PS::S32 b) {
    for(PS::S32 i = a ; i < b ; ++ i){
      if (sph_system[i].type == FLUID) {
        sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * sph_system[i].acc;
      }
    }
  });
#else
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    if (sph_system[i].type == FLUID) {
      sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * sph_system[i].acc;
    }
  }
#endif
}

#if REUSE
void setPrevPos(PS::ParticleSystem<FP>& sph_system){
#ifdef PARTICLE_SIMULATOR_TASK_PARALLEL
  mtbb::parallel_for(0, sph_system.getNumberOfParticleLocal(), 1, CUTOFF_PFOR, [&] (PS::S32 a, PS::S32 b) {
    for(PS::S32 i = a; i < b; ++ i){
      sph_system[i].prev_pos = sph_system[i].pos;
    }
  });
#else
#ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
#pragma omp parallel for
#endif
  for(PS::S32 i = 0 ; i < sph_system.getNumberOfParticleLocal() ; ++ i){
    sph_system[i].prev_pos = sph_system[i].pos;
  }
#endif
}
#endif

int main(int argc, char* argv[]){
  // Initialize FDPS
  PS::Initialize(argc, argv);
  // Make a directory
  char dir_name[1024];
  sprintf(dir_name,"./result");
  makeOutputDirectory(dir_name);
  // Display # of MPI processes and threads
  PS::S32 nprocs = PS::Comm::getNumberOfProc();
  PS::S32 nthrds = PS::Comm::getNumberOfThread();
  std::cout << "===========================================" << std::endl
            << " This is a sample program of "               << std::endl
            << " Smoothed Particle Hydrodynamics on FDPS!"   << std::endl
            << " # of processes is " << nprocs               << std::endl
            << " # of thread is    " << nthrds               << std::endl
            << "===========================================" << std::endl;
  // Make an instance of ParticleSystem and initialize it
  PS::ParticleSystem<FP> sph_system;
  sph_system.initialize();
#if PARTICLE_SIMULATOR_TASK_PARALLEL
#if DISABLE_STEAL
  myth_adws_set_stealable(0);
#endif
#endif
  // Define local variables
  // Make an initial condition and initialize the particle system
  SetupIC(sph_system);
  // Make an instance of DomainInfo and initialize it
  PS::DomainInfo dinfo;
  dinfo.initialize();
  // Set the boundary condition
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
  // Perform domain decomposition
  dinfo.decomposeDomainAll(sph_system);
  // Exchange the SPH particles between the (MPI) processes
  sph_system.exchangeParticle(dinfo);
  // Make two tree structures
  // (one is for the density calculation and
  //  another is for the force calculation.)
  PS::TreeForForceShort<Dens, EP, EP>::Gather dens_tree;
  dens_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());

  PS::TreeForForceShort<Hydro, EP, EP>::Symmetry hydr_tree;
  hydr_tree.initialize(3 * sph_system.getNumberOfParticleGlobal());
  // Compute density, pressure, acceleration due to pressure gradient
#if REUSE
  bool reuse = false;
  PS::S32 reuse_count = 0;
#endif
  // Main loop for time integration
  PS::F64 dt;
  PS::S32 step = 0;
  PS::S64 t_all = 0;
  for(PS::F64 time = 0 ; time < END_TIME && step < MAX_STEP ; time += dt, ++ step){
    PS::U64 t1 = get_time_in_usec();

    if (step > 0) {
      // Leap frog: Initial Kick & Full Drift
      InitialKick(sph_system, dt);
#if REUSE
      reuse = FullDrift(sph_system, dt);
#else
      FullDrift(sph_system, dt);
#endif
    }
    // Leap frog: Predict
    // Predict(sph_system, dt);
    // Perform domain decomposition again
    // dinfo.decomposeDomainAll(sph_system);
    // Exchange the SPH particles between the (MPI) processes
    // sph_system.exchangeParticle(dinfo);
    // Compute density, pressure, acceleration due to pressure gradient
#if REUSE
    if (reuse) {
      dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo, false, PS::REUSE_LIST);
      hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, false, PS::REUSE_LIST);
    } else {
      setPrevPos(sph_system);
      dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo, false, PS::MAKE_LIST_FOR_REUSE);
      hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo, false, PS::MAKE_LIST_FOR_REUSE);
      reuse_count++;
    }
#else
    dens_tree.calcForceAllAndWriteBack(CalcDensity()   , sph_system, dinfo);
    hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
#endif
    if (step > 0) {
      // Leap frog: Final Kick
      FinalKick(sph_system, dt);
    }
    // Get a new timestep
    dt = getTimeStep(sph_system);

    PS::U64 t2 = get_time_in_usec();
    t_all += t2 - t1;

    // Output result files
#if OUTPUT
    if(step % OUTPUT_INTERVAL == 0){
      char filename[256];
#if DAM_2D
      sprintf(filename, "result/dambreaking2d.txt.%d", step / OUTPUT_INTERVAL);
#else
      sprintf(filename, "result/dambreaking3d.txt.%d", step / OUTPUT_INTERVAL);
#endif
      sph_system.writeParticleAscii(filename);
      if (PS::Comm::getRank() == 0){
        std::cout << "================================" << std::endl;
        std::cout << "output " << filename << "." << std::endl;
        std::cout << "================================" << std::endl;
      }
#ifdef RECORD_SPLIT
      {
        sprintf(filename, "result/split.txt.%d", step / OUTPUT_INTERVAL);
        std::ofstream fout(filename);
        hydr_tree.dumpSplit(fout);
        fout.close();
      }
#endif
#ifdef RECORD_TIMELINE
      {
        sprintf(filename, "result/timeline.txt.%d", step / OUTPUT_INTERVAL);
        std::ofstream fout(filename);
        hydr_tree.dumpTimeline(fout);
        fout.close();
      }
#endif
    }
#endif
    // Output information to STDOUT
    if (PS::Comm::getRank() == 0){
      std::cout << "time: "    << std::fixed   << std::setprecision(5) << time    << " [s] "
                << "step: "    << std::setw(5) << std::right           << step    << " "
                << "elapsed: " << std::setw(7) << std::right           << t2 - t1 << " [usec] "
                << "reuse: "   << (reuse ? "true" : "false")                      << std::endl;
    }
  }
#if REUSE
  std::cout << "reuse = " << reuse_count << std::endl;
  std::cout << "total time = " << (PS::F64)t_all / 1000 / 1000 << " sec" << std::endl;
#endif
  // Finalize FDPS
  PS::Finalize();
  return 0;
}

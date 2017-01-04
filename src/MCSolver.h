/*
 @file    MCSolver.h
 @brief   MCSolver class, a subclass of solver
 @author  Luke Eure
 @date    May 11 2016
*/

#ifndef MONTE_CARLO_SOLVER_H
#define MONTE_CARLO_SOLVER_H

#include <cmath>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>

#include "Point.h"
#include "Universe.h"
#include "Material.h"
#include "Geometry.h"
#include "LocalCoords.h"
#include "Solver.h"

#include "Tally.h"
#include "Neutron.h"
#include "Fission.h"

class MCSolver : public Solver {

enum tally_names {CROWS, NUM_CROWS, LEAKS, ABSORPTIONS, FISSIONS};
enum fission_bank_names {OLD, NEW};

public:

  MCSolver();
  virtual ~MCSolver();

  void setGeometry(Geometry* geometry);

  void sampleLocation(Neutron* neutron);

  void computeEigenvalue(int n_histories, int num_batches, int num_groups);

  void transportNeutron(std::vector <Tally> &tallies, bool first_round,
      Fission* fission_banks, int neutron_num, int batch,
      int write_neutron = -1, int write_batch = -1, 
      Neutron* input_neutron = NULL);
  void transportNeutronWithTrack(std::vector <Tally> &tallies, bool first_round,
      Fission* fission_banks, int neutron_num, int batch,
      int write_neutron = -1, int write_batch = -1, 
      Neutron* input_neutron = NULL);
  
  Geometry* getGeometry();
  virtual FP_PRECISION getFlux(int fsr_id, int group);
  FP_PRECISION getKeff();

  void saveBadNeutron(Neutron* neutron, int neutron_num, int batch);
  void trackSingleNeutron();

  void initializeLocks();
  void setNumThreads(int num_threads);
  void setFixedSourceByCell(Cell* cell, int group, float source);

  // functions that make MCSolver compatable with Solver
  virtual void computeFSRFissionRates(double* fission_rates, int num_FSRs);
  virtual void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes);
  virtual void setFluxes(FP_PRECISION* in_fluxes, int num_fluxes);
  virtual void flattenFSRFluxes(FP_PRECISION value);
  virtual void computeFSRFissionSources();
  virtual void computeFSRScatterSources();
  virtual void initializeSourceArrays();
  virtual void addSourceToScalarFlux();
  virtual void initialize();
  virtual void initializeFluxArrays();
  virtual void computeFSRSources();
  virtual void zeroTrackFluxes();
  virtual void normalizeFluxes();
  virtual void transportSweep();
  virtual void storeFSRFluxes();
  virtual void computeKeff();
  virtual double computeResidual(residualType res_type);

private:

  double _k_eff;
  double _num_threads;
  Geometry* _geometry;
  Universe* _root_universe;
  FP_PRECISION* _FSR_volumes;
  Material** _FSR_materials;
  Cell* _root_cell;
  int _num_groups;
  int _num_FSRs;
  int _num_materials;
  FP_PRECISION* _scalar_flux;
  FP_PRECISION* _cumulative_scalar_flux;

  omp_lock_t* _leak_lock;
  omp_lock_t* _absorption_lock;
  omp_lock_t* _fission_lock;
  omp_lock_t* _crow_lock;
  std::vector <omp_lock_t*> _flux_locks;

  /** A mapping of fixed sources keyed by the pair (FSR ID, energy group) */
  std::map< std::pair<Cell*, int>, FP_PRECISION > _fix_src_cell_map;

  bool _fixed_sources_exist;

};

#endif

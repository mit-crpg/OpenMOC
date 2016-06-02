/* 
 @file      MCSolver.h
 @brief     MCSolver class, a subclass of solver
 @author    Luke Eure
 @date      May 11 2016
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

    MCSolver(TrackGenerator* track_generator=NULL);
    virtual ~MCSolver();

    void setGeometry(Geometry* geometry, Cell* root_cell);
    void initialize(Lattice* lattice);

    void sampleLocation(Neutron* neutron);

    void computeEigenValue(int n_histories, int num_batches);

    void transportNeutron(std::vector <Tally> &tallies, bool first_round,
            Fission* fission_banks, int neutron_num);
    
    Geometry* getGeometry();
    virtual FP_PRECISION getFlux(int fsr_id, int group);
    FP_PRECISION getKeff();

    // functions that make MCSolver compatable with Solver
    virtual void computeFSRFissionRates(double* fission_rates, int num_FSRs);
    virtual void getFluxes(FP_PRECISION* out_fluxes, int num_fluxes);
    virtual void setFluxes(FP_PRECISION* in_fluxes, int num_fluxes);
    virtual void flattenFSRFluxes(FP_PRECISION value);
    virtual void computeFSRFissionSources();
    virtual void computeFSRScatterSources();
    virtual void initializeSourceArrays();
    virtual void addSourceToScalarFlux();
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
    Geometry* _geometry;
    Universe* _root_universe;
    FP_PRECISION* _FSR_volumes;
    Material** _FSR_materials;
    Cell* _root_cell;
    int _num_groups;
    int _num_FSRs;
    int _num_materials;
    FP_PRECISION* _scalar_flux;

};

#endif

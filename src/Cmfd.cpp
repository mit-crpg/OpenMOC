#include "Cmfd.h"

/**
 * Cmfd constructor
 * @param geom pointer to the geometry
 * @param solve_methor
 * @param relaxation factor for CMFD acceleration
 */
Cmfd::Cmfd(Geometry* geometry, solveType solve_method, double relax_factor) {

  /* Create objects */
  _geometry = geometry;
  _mesh = geometry->getMesh();
  _cells_x = _mesh->getCellsX();
  _cells_y = _mesh->getCellsY();
  _quad = new Quadrature(TABUCHI);
  _timer = new Timer();
  
  /* Boolean and Enum flags to toggle features */
  _optically_thick = false;
  _solve_method = solve_method;
  _relax_factor = relax_factor;
  _assemble_M = true;
  _flux_method = PRIMAL;
  
  /* Global variables used in solving Cmfd problem */
  _l2_norm = 1.0;
  _l2_norm_old = 1.0;
  _keff_conv = 1e-8;
  
  /* General problem parameters */
  _num_groups = _mesh->getNumGroups();
  _num_fsrs = _mesh->getNumFSRs();
  
  /* Create Cmfd matrix objects */
  int petsc_err;
  petsc_err = createAMPhi();
  
  /* If solving diffusion problem, create arrays for FSR parameters */
  if (_solve_method == DIFFUSION){
    _FSR_volumes = (FP_PRECISION*)calloc(_num_fsrs, sizeof(FP_PRECISION));
    _FSR_materials = new Material*[_num_fsrs];
    _FSR_fluxes = (FP_PRECISION*)calloc(_num_fsrs*_num_groups, sizeof(FP_PRECISION));
    initializeFSRs();
  }
}

/**
 * cmfd Destructor clears all memory
 */
Cmfd::~Cmfd() {
}


/**
 * Create the loss matrix (A), fission matrix (A), and flux vector (phi)
 * @param number of columns needed in A matrix
 * @param number of rows needed in A matrix and M matrix
 * @param number of rows needed in phi vector
 */
int Cmfd::createAMPhi(){

  log_printf(DEBUG, "Creating AMPhi...");
  
  int petsc_err = 0;
  PetscInt size1 = _cells_x * _cells_y * _num_groups;
  PetscInt size2 = 4 + _num_groups;
  
  /* Create _A (loss) matrix, _A_array array to pre-store matrix values,
   * and _indices_y_A to store indices of matrix elements */
  petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &_A);
  _A_array = new PetscScalar[size2];
  _indices_y_A = new PetscInt[size2];
  
  /* Create _M (gain) matrix, _M_array array to pre-store matrix values,
   * and _indices_y_M to store indices of matrix elements */
  size2 = size2 - 4;
  petsc_err = MatCreateSeqAIJ(PETSC_COMM_WORLD, size1, size1, size2, PETSC_NULL, &_M);
  _M_array = new PetscScalar[size2];
  _indices_y_M = new PetscInt[size2];
  
  /* Create flux, source, and residual vectors */
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_phi_new);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_phi_old);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_source_old);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_sold);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_snew);
  petsc_err = VecCreateSeq(PETSC_COMM_WORLD, _cells_x*_cells_y*_num_groups, &_res);
  CHKERRQ(petsc_err);
  
  CHKERRQ(petsc_err);
  
  return petsc_err;
}



/* Create cross sections and fluxes for each cmfd cell by
 * energy condensing and volume averaging cross sections from
 * the MOC sweep.
 * @param pointer to an array of fsrs
 */
void Cmfd::computeXS(){

  log_printf(INFO, "computing cmfd cross sections...");

  /* split corner currents to side surfaces */
  if (_solve_method == MOC)
    _mesh->splitCorners();
  
  /* initialize variables for FSR properties*/
  double volume, flux, abs, tot, nu_fis, chi, buckle;
  double* scat;
  Material** materials = _mesh->getMaterials();
  FP_PRECISION* fluxes = _mesh->getFluxes("old_flux");
  
  /* initialize tallies for each parameter */
  double abs_tally, nu_fis_tally, dif_tally, rxn_tally, vol_tally, tot_tally;
  double scat_tally[_num_groups];
  
  /* interator to loop over fsrs in each mesh cell */
  std::vector<int>::iterator iter;
  
  /* create pointers to objects */
  Material* fsr_material;
  Material* cell_material;
  
  /* loop over mesh cells */
  for (int i = 0; i < _cells_x * _cells_y; i++){
    
    /* loop over energy groups */
    for (int e = 0; e < _num_groups; e++) {
      
      /* zero tallies for this group */
      abs_tally = 0;
      nu_fis_tally = 0;
      dif_tally = 0;
      rxn_tally = 0;
      vol_tally = 0;
      tot_tally = 0;
      
      /* zero each group to group scattering tally */
      for (int g = 0; g < _num_groups; g++){
	scat_tally[g] = 0;
      }
      
      /* loop over FSRs in mesh cell */
      for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter){
	
	/* Gets FSR volume, material, and cross sections */
	fsr_material = _FSR_materials[*iter];
	cell_material = materials[i];
	volume = _FSR_volumes[*iter];
	flux = _FSR_fluxes[(*iter)*_num_groups+e];
	abs = fsr_material->getSigmaA()[e];
	tot = fsr_material->getSigmaT()[e];
	scat = fsr_material->getSigmaS();
	buckle = fsr_material->getBuckling()[e];

	/* if material has a diffusion coefficient, use it; otherwise
	 * estimate the diffusion coefficient with 1 / (3 * sigma_t) */
        //if (fsr_material->getDifCoef()[e] > 1e-6){
	//  dif_tally += fsr_material->getDifCoef()[e] * flux * volume;
	//}
	//else
	dif_tally += flux * volume / (3.0 * tot);
	
	/* if material has a chi, use it; otherwise set to 0 */
	if (fsr_material->getChi() != NULL)
	  chi = fsr_material->getChi()[e];
	else
	  chi = 0.0;
	
	/* if material has a nu_sig_f, use it; otherwise set to 0 */
	if (fsr_material->getNuSigmaF() != NULL)
	  nu_fis = fsr_material->getNuSigmaF()[e];
	else
	  nu_fis = 0.0;
	
	/* increment tallies for this group */
	abs_tally += abs * flux * volume;
	tot_tally += tot * flux * volume;
	nu_fis_tally += nu_fis * flux * volume;
	rxn_tally += flux * volume;
	vol_tally += volume;
	for (int g = 0; g < _num_groups; g++)
	  scat_tally[g] += scat[g*_num_groups + e] * flux * volume;
	
	/* choose a chi for this group */
	if (chi >= cell_material->getChi()[e])
	  cell_material->setChiByGroup(chi,e);
      }
      

      /* set the mesh cell properties with the tallies */
      _mesh->setVolume(vol_tally, i);
      cell_material->setSigmaAByGroup(abs_tally / rxn_tally, e);
      cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e);
      cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e);
      cell_material->setDifCoefByGroup(dif_tally / rxn_tally, e);
      fluxes[i*_num_groups+e] = rxn_tally / vol_tally;
      cell_material->setBucklingByGroup(buckle, e);

      log_printf(DEBUG, "cell: %i, group: %i, vol: %f, siga: %f, sigt: %f, nu_sigf: %f, dif_coef: %f, flux: %f", i, e, vol_tally, abs_tally / rxn_tally, 
		 tot_tally / rxn_tally, nu_fis_tally / rxn_tally, dif_tally / rxn_tally, rxn_tally / vol_tally);
      
      for (int g = 0; g < _num_groups; g++)
	cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, e, g);
    }
  }
}


/* Compute the diffusion coefficients (d_dif - straight diffusion coefficient,
 * d_hat - surface diffusion coefficient, and d_tilde - surface diffusion coefficient
 * correction factor) for each mesh while ensuring neutron balance is achieved.
 */
void Cmfd::computeDs(){

  log_printf(INFO, "computing cmfd Ds...");
  
  /* initialize variables */
  double d = 0, d_next = 0, d_hat = 0, d_tilde = 0;
  double current = 0, flux = 0, flux_next = 0, f = 1, f_next = 1;
  double length, length_perpen, next_length_perpen;
  double sense;
  int next_surface;
  int cell, cell_next;
  
  Material** materials = _mesh->getMaterials();
  FP_PRECISION* cell_flux = _mesh->getFluxes("old_flux");
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();
  FP_PRECISION* currents = _mesh->getCurrents();
  
  /* loop over mesh cells in y direction */
  for (int y = 0; y < _cells_y; y++){
    
    /* loop over mesh cells in x direction */
    for (int x = 0; x < _cells_x; x++){
      
      cell = y*_cells_x+x;
      
      /* loop over surfaces in a cell */
      for (int surface = 0; surface < 4; surface++){
	
	/* loop over groups */
	for (int e = 0; e < _num_groups; e++){
	  
	  /* get diffusivity and flux for mesh cell */
	  d = materials[cell]->getDifCoef()[e];
	  flux = cell_flux[cell*_num_groups+e];
	  cell_next = _mesh->getCellNext(cell, surface);
	  
	  /* set sense of the surface */
	  if (surface == 0 || surface == 3)
	    sense = -1.0;
	  else
	    sense = 1.0;
	  
	  /* set the length of this surface and the
	   * perpendicular surface */
	  if (surface == 0 || surface== 2){
	    length = heights[y];
	    length_perpen = widths[x];
	  }
	  else if (surface == 1 || surface == 3){
	    length = widths[x];
	    length_perpen = heights[y];
	  }
	  
	  /* compute the optical thickness correction factor */
	  f = computeDiffCorrect(d, length_perpen);
	  
	  /* if surface is on a boundary, choose appropriate BCs */
	  if (cell_next == -1){
	    
	    /* REFLECTIVE BC */
	    if (_mesh->getBoundary(surface) == REFLECTIVE){
	      
	      /* set d's */
	      current = 0.0;
	      d_hat = 0.0;
	      d_tilde = 0.0;
	    }
	    /* VACUUM BC */
	    else if (_mesh->getBoundary(surface) == VACUUM){
	      current = sense * currents[cell*_num_groups*8 + surface*_num_groups + e];
	      
	      /* set d's */
	      d_hat =  2 * d*f / length_perpen / (1 + 4 * d*f / length_perpen);
	      d_tilde = (sense * d_hat * flux - current / length) / flux;
	    }
	    /* ZERO_FLUX BC */
	    else if (_mesh->getBoundary(surface) == ZERO_FLUX){
	      current = sense * currents[cell*_num_groups*8 + surface*_num_groups + e];
	      
	      /* set d's */
	      d_hat = 2 * d*f / length_perpen;
	      d_tilde = (sense * d_hat * flux - current / length) / flux;
	    }
	  }
	  /* if surface is an interface, use finite differencing */
	  else{
	    
	    /* set properties for cell next to surface */
	    if (surface == 0){
	      next_length_perpen = _mesh->getLengthsX()[cell_next % _cells_x];
	      next_surface = 2;
	    }
	    else if (surface == 1){
	      next_length_perpen = _mesh->getLengthsY()[cell_next / _cells_x];
	      next_surface = 3;
	    }
	    else if (surface == 2){
	      next_length_perpen = _mesh->getLengthsX()[cell_next % _cells_x];
	      next_surface = 0;
	    }
	    else if (surface == 3){
	      next_length_perpen = _mesh->getLengthsY()[cell_next / _cells_x];
	      next_surface = 1;
	    }
	    
	    /* set diffusion coefficient and flux for neighboring cell */
	    d_next = materials[cell_next]->getDifCoef()[e];
	    flux_next = _mesh->getFluxes("old_flux")[cell_next*_num_groups + e];
	    
	    /* get optical thickness correction term for meshCellNext */
	    f_next = computeDiffCorrect(d_next, next_length_perpen);

	    /* compute d_hat */
	    d_hat = 2.0 * d * f * d_next * f_next / (length_perpen * d * f + next_length_perpen * d_next*f_next);
	    
	    /* compute net current */
	    current = sense * currents[cell*_num_groups*8 + surface*_num_groups + e]
	      - sense * currents[cell_next*_num_groups*8 + next_surface*_num_groups + e];
	    
	    /* compute d_tilde */
	    d_tilde = -(sense * d_hat * (flux_next - flux) + current  / length) / (flux_next + flux);
	  }
	  
	  /* if the magnitude of d_tilde is greater than the magnitude of d_hat,
	   * select new values d_tilde and d_hat to ensure the course mesh equations
	   * are guaranteed to be diagonally dominant */
	  if (fabs(d_tilde) > fabs(d_hat)){

	    if (sense == -1){
	      /* d_tilde is positive */
	      if (1 - fabs(d_tilde)/d_tilde < 1e-8){
		d_hat   = - current/(2*flux*length);
		d_tilde = - current/(2*flux*length);
	      }
	      /* if d_tilde is negative */
	      else{
		d_hat   = current/(2*flux_next*length);
		d_tilde = - current/(2*flux_next*length);
	      }
	    }
	    else{
	      /* d_tilde is positive */
	      if (1 - fabs(d_tilde)/d_tilde < 1e-8){
		d_hat   = - current/(2*flux_next*length);
		d_tilde = - current/(2*flux_next*length);
	      }
	      /* if d_tilde is negative */
	      else{
		d_hat   = current/(2*flux*length);
		d_tilde = - current/(2*flux*length);
	      }
	    }
	  }
	  

	  
	  /* set d_hat and d_tilde */
	  d_tilde = materials[cell]->getDifTilde()[surface*_num_groups + e] * (1 - _relax_factor) + _relax_factor * d_tilde;
	  materials[cell]->setDifHatByGroup(d_hat, e, surface);
	  materials[cell]->setDifTildeByGroup(d_tilde, e, surface);

	  log_printf(DEBUG, "cell: %i, group: %i, side: %i, flux: %f, current: %f, d: %f, dhat: %f, dtilde: %f", y*_cells_x + x, e, surface, flux, current, d, d_hat, d_tilde);
	  
	}
      }
    }
  }
}


/*
 * CMFD solver that solves the diffusion problem
 * @param solve methed - either diffusion or cmfd (acceleration)
 * @param iteration number of in MOC solver - used for plotting
 * @return k-effective
 */
double Cmfd::computeKeff(){

  log_printf(INFO, "Running diffusion solver...");

  if (_solve_method == DIFFUSION)
    _timer->startTimer();
  
  int petsc_err = 0;
  int max_outer, iter = 0;
  PetscScalar sumold, sumnew, scale_val, eps;
  PetscReal rtol = 1e-6;
  PetscReal atol = 1e-6;
  Vec sold, snew, res;
  float criteria = 1e-8;
  
  KSP ksp;

  computeXS();
  computeDs();
  
  petsc_err = MatZeroEntries(_A);
  if (_assemble_M)
    petsc_err = MatZeroEntries(_M);
    
  petsc_err = constructMatrices();
  CHKERRQ(petsc_err);

  petsc_err = VecAssemblyBegin(_phi_old);
  petsc_err = VecAssemblyEnd(_phi_old);
  petsc_err = VecAssemblyBegin(_phi_new);
  petsc_err = VecAssemblyEnd(_phi_new);
  petsc_err = VecAssemblyBegin(_source_old);
  petsc_err = VecAssemblyEnd(_source_old);
  petsc_err = VecAssemblyBegin(_sold);
  petsc_err = VecAssemblyEnd(_sold);
  petsc_err = VecAssemblyBegin(_snew);
  petsc_err = VecAssemblyEnd(_snew);
  petsc_err = VecAssemblyBegin(_res);
  petsc_err = VecAssemblyEnd(_res);
  petsc_err = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  petsc_err = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);  

  if (_assemble_M){
    petsc_err = MatAssemblyBegin(_M, MAT_FINAL_ASSEMBLY);
    petsc_err = MatAssemblyEnd(_M, MAT_FINAL_ASSEMBLY);
  }
  CHKERRQ(petsc_err);
  
  petsc_err = KSPCreate(PETSC_COMM_WORLD, &ksp);
  petsc_err = KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, PETSC_DEFAULT);
  petsc_err = KSPSetType(ksp, KSPGMRES);
  petsc_err = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  petsc_err = KSPSetOperators(ksp, _A, _A, SAME_NONZERO_PATTERN);
  petsc_err = KSPSetUp(ksp);
  petsc_err = KSPSetFromOptions(ksp);
  CHKERRQ(petsc_err);  

  petsc_err = VecCopy(_phi_old, _phi_new);
  petsc_err = MatMult(_M, _phi_new, _snew);
  petsc_err = VecSum(_snew, &sumnew);
  petsc_err = MatMult(_A, _phi_new, _sold);
  petsc_err = VecSum(_sold, &sumold);
  CHKERRQ(petsc_err);
  _keff = float(sumnew/sumold);

  petsc_err = MatMult(_M, _phi_new, _sold);
  petsc_err = VecSum(_sold, &sumold);
  scale_val = (_cells_x * _cells_y * _num_groups) / sumold;
  petsc_err = VecScale(_sold, scale_val);
  sumold = _cells_x * _cells_y * _num_groups;
  CHKERRQ(petsc_err);

  /* diffusion solver */
  for (iter = 0; iter < 1000; iter++){

    petsc_err = KSPSolve(ksp, _sold, _phi_new);

    /* compute the new source */
    petsc_err = MatMult(_M, _phi_new, _snew);
    petsc_err = VecSum(_snew, &sumnew);
    CHKERRQ(petsc_err);
    
    /* compute and set keff */
    _keff = sumnew / sumold;
    petsc_err = VecScale(_sold, _keff);
    
    /* compute the L2 norm of source error */
    scale_val = 1e-15;
    petsc_err = VecShift(_snew, scale_val);
    petsc_err = VecPointwiseDivide(_res, _sold, _snew);
    scale_val = -1;
    petsc_err = VecShift(_res, scale_val);
    CHKERRQ(petsc_err);
    petsc_err = VecNorm(_res, NORM_2, &eps);
    
    /* compute error */
    eps = eps / (_cells_x * _cells_y * _num_groups);
    scale_val = (_cells_x * _cells_y * _num_groups) / sumnew;
    petsc_err = VecScale(_snew, scale_val);
    CHKERRQ(petsc_err);
    
    /* set old source to new source */
    petsc_err = VecCopy(_snew, _sold);
    CHKERRQ(petsc_err);
    
    log_printf(INFO, "CMFD iter: %i, keff: %f, error: %f", iter + 1, _keff, eps);

    /* check for convergence */
    if (iter > 5 && eps < criteria){
      break;
    }
  }

  log_printf(INFO, "CMFD KEFF: %f", _keff);

  petsc_err = rescaleFlux();
  petsc_err = setMeshCellFlux();
  updateMOCFlux();  

  petsc_err = VecCopy(_snew, _source_old);
  std::vector<int> fsrs;
  std::vector<int>::iterator fsr_iter;
  
  if (_solve_method == DIFFUSION){
    for (int i = 0; i < _cells_x*_cells_y; i++){
      fsrs = _mesh->getCellFSRs()->at(i);
      for (fsr_iter = fsrs.begin(); fsr_iter != fsrs.begin(); fsr_iter++){
	for (int e = 0; e < _num_groups; e++){
	  _FSR_fluxes[*fsr_iter*_num_groups+e] = _mesh->getFluxes("new_flux")[i*_num_groups + e];
	  _mesh->getFluxes("old_flux")[i*_num_groups + e] = _mesh->getFluxes("new_flux")[i*_num_groups + e];
	}
      }
    }
  }
  

  if (_solve_method == DIFFUSION){
    std::string msg_string;
    log_printf(TITLE, "TIMING REPORT");
    _timer->stopTimer();
    _timer->recordSplit("Total time to solve diffusion eigenvalue problem");
    
    double tot_time = _timer->getSplit("Total time to solve diffusion eigenvalue problem");
    msg_string = "Total time to solve diffusion eigenvalue problem";
    msg_string.resize(53, '.');
    log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);
  }
  
  return _keff;
}


/* rescale flux */
int Cmfd::rescaleFlux(){

  int petsc_err = 0;
  PetscScalar sumnew, sumold, scale_val;


  /* rescale the new and old flux to have an avg source of 1.0 */
  petsc_err = MatMult(_M, _phi_new, _snew);
  petsc_err = VecSum(_snew, &sumnew);
  scale_val = _cells_x*_cells_y*_num_groups / sumnew;
  petsc_err = VecScale(_phi_new, scale_val);
  petsc_err = MatMult(_M, _phi_old, _sold);
  petsc_err = VecSum(_sold, &sumold);
  scale_val = _cells_x*_cells_y*_num_groups / sumold;
  petsc_err = VecScale(_phi_old, scale_val);
  CHKERRQ(petsc_err);
  
  return 0;
}


int Cmfd::setMeshCellFlux(){

  /* initialize variables */

  int petsc_err = 0;
  PetscScalar *old_phi;
  PetscScalar *new_phi;
  petsc_err = VecGetArray(_phi_old, &old_phi);
  petsc_err = VecGetArray(_phi_new, &new_phi);
  CHKERRQ(petsc_err);
  

  for (int i = 0; i < _cells_x*_cells_y; i++){
    for (int e = 0; e < _num_groups; e++){
      _mesh->getFluxes("old_flux")[i*_num_groups + e] = double(old_phi[i*_num_groups + e]);
      _mesh->getFluxes("new_flux")[i*_num_groups + e] = double(new_phi[i*_num_groups + e]);
    }
  }
  

  petsc_err = VecRestoreArray(_phi_old, &old_phi);
  petsc_err = VecRestoreArray(_phi_new, &new_phi);
  CHKERRQ(petsc_err);

  return 0;
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @param A matrix
 * @param M matrix
 * @param old flux vector
 * @param solve methed - either DIFFUSION or CMFD
 * @return petsc error indicator
 */
int Cmfd::constructMatrices(){

  log_printf(INFO,"Constructing AMPhi...");
  

  /* initialized variables */
  int petsc_err = 0;

  PetscInt indice_x;
  PetscScalar value, phi, b_prime = 0;
  int cell;
  
  Material** materials = _mesh->getMaterials();
  FP_PRECISION* old_flux = _mesh->getFluxes("old_flux");
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();
  
  for (int y = 0; y < _cells_y; y++){
    

    for (int x = 0; x < _cells_x; x++){
      
      cell = y*_cells_x + x;
      

      for (int e = 0; e < _num_groups; e++){
	

	for (int f = 0; f < _num_groups; f++){
	  _M_array[f] = 0.0;
	}
	for (int f = 0; f < _num_groups+4; f++){
	  _A_array[f] = 0.0;
	}
	

	indice_x = (y*_cells_x + x)*_num_groups+e;
	if (x != 0)
	  _indices_y_A[0] = (y*_cells_x + x-1)*_num_groups+e;
	else
	  _indices_y_A[0] = -1;
	
	if (y != _cells_y - 1)
	  _indices_y_A[1] = ((y+1)*_cells_x + x)*_num_groups+e;
	else
	  _indices_y_A[1] = -1;
	
	if (x != _cells_x-1)
	  _indices_y_A[_num_groups+2] = (y*_cells_x + x+1)*_num_groups+e;
	else
	  _indices_y_A[_num_groups+2] = -1;
	
	if (y != 0)
	  _indices_y_A[_num_groups+3] = ((y-1)*_cells_x + x)*_num_groups+e;
	else
	  _indices_y_A[_num_groups+3] = -1;
	

	for (int g = 0; g < _num_groups; g++){
	  _indices_y_A[g+2] = (y*_cells_x + x)*_num_groups+g;
	  _indices_y_M[g] = (y*_cells_x + x)*_num_groups+g;
	}
	

	phi = old_flux[cell*_num_groups + e];

	value = materials[cell]->getSigmaA()[e] * _mesh->getVolumes()[cell];
	_A_array[2 + e] += value;

	/* set buckling */
	//value = materials[cell]->getDifCoef()[e] * materials[cell]->getBuckling()[e] * _mesh->getVolumes()[cell];
	//_A_array[2 + e] += value;

	for (int g = 0; g < _num_groups; g++){
	  if (e != g){
	    value = materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell];
	    _A_array[2 + e] += value;
	  }
	}
	

	if (_flux_method == PRIMAL){
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){
	      value = - materials[cell]->getSigmaS()[g*_num_groups + e] * _mesh->getVolumes()[cell];
	      _A_array[2 + g] += value;
	    }
	  }
	}
	else{
	  for (int g = 0; g < _num_groups; g++){
	    if (e != g){
	      value = - materials[cell]->getSigmaS()[e*_num_groups + g] * _mesh->getVolumes()[cell];
	      _A_array[2 + g] += value;
	    }
	  }
	}
	
	/* RIGHT SURFACE */
       
	/* set transport term on diagonal */
	if (_solve_method == MOC)
	  value = (materials[cell]->getDifHat()[2*_num_groups + e]      - materials[cell]->getDifTilde()[2*_num_groups + e]) * heights[cell / _cells_x];
	else if (_solve_method == DIFFUSION)
	  value = materials[cell]->getDifHat()[2*_num_groups + e] * heights[cell / _cells_x];
	
	_A_array[2 + e] += value;
	
	/* set transport term on off diagonal */
	if (x != _cells_x - 1){
	  if (_solve_method == MOC)
	    value = - (materials[cell]->getDifHat()[2*_num_groups + e] + materials[cell]->getDifTilde()[2*_num_groups + e]) * heights[cell / _cells_x];
	  else if (_solve_method == DIFFUSION)
	    value = - materials[cell]->getDifHat()[2*_num_groups + e] * heights[cell / _cells_x];
	
	  _A_array[2 + _num_groups] += value;
	}
	
	/* LEFT SURFACE */
	
	/* set transport term on diagonal */
	if (_solve_method == MOC)
	  value = (materials[cell]->getDifHat()[0*_num_groups + e]      + materials[cell]->getDifTilde()[0*_num_groups + e]) * heights[cell / _cells_x];
	else if (_solve_method == DIFFUSION)
	  value = materials[cell]->getDifHat()[0*_num_groups + e] * heights[cell / _cells_x];
	
	_A_array[2 + e] += value;
	
	/* set transport term on off diagonal */
	if (x != 0){
	  if (_solve_method == MOC)
	    value = - (materials[cell]->getDifHat()[0*_num_groups + e] - materials[cell]->getDifTilde()[0*_num_groups + e]) * heights[cell / _cells_x];
	  else if (_solve_method == DIFFUSION)
	    value = - materials[cell]->getDifHat()[0*_num_groups + e] * heights[cell / _cells_x];
	  
	  _A_array[0] += value;
	}
	
	/* BOTTOM SURFACE */
	
	/* set transport term on diagonal */
	if (_solve_method == MOC)
	  value = (materials[cell]->getDifHat()[1*_num_groups + e]      - materials[cell]->getDifTilde()[1*_num_groups + e]) * widths[cell % _cells_x];
	else if (_solve_method == DIFFUSION)
	  value = materials[cell]->getDifHat()[1*_num_groups + e] * widths[cell % _cells_x];
	
	_A_array[2 + e] += value;
	
	/* set transport term on off diagonal */
	if (y != _cells_y - 1){
	  if (_solve_method == MOC)
	    value = - (materials[cell]->getDifHat()[1*_num_groups + e] + materials[cell]->getDifTilde()[1*_num_groups + e]) * widths[cell % _cells_x];
	  else if (_solve_method == DIFFUSION)
	    value = - materials[cell]->getDifHat()[1*_num_groups + e] * widths[cell % _cells_x];
	  
	  _A_array[1] += value;
	}
	
	/* TOP SURFACE */
	
	/* set transport term on diagonal */
	if (_solve_method == MOC)
	  value = (materials[cell]->getDifHat()[3*_num_groups + e]      + materials[cell]->getDifTilde()[3*_num_groups + e]) * widths[cell % _cells_x];
	else if (_solve_method == DIFFUSION)
	  value = materials[cell]->getDifHat()[3*_num_groups + e] * widths[cell % _cells_x];
	
	_A_array[2 + e] += value;
	
	/* set transport term on off diagonal */
	if (y != 0){
	  if (_solve_method == MOC)
	    value = - (materials[cell]->getDifHat()[3*_num_groups + e] - materials[cell]->getDifTilde()[3*_num_groups + e]) * widths[cell % _cells_x];
	  else if (_solve_method == DIFFUSION)
	    value = - materials[cell]->getDifHat()[3*_num_groups + e] * widths[cell % _cells_x];
	  
	  _A_array[3 + _num_groups] += value;
	}
	

	if (_assemble_M){
	  for (int g = 0; g < _num_groups; g++){
	    
	    int f = g;
	    
	    if (_flux_method == ADJOINT){
	      g = e;
	      e = f;
	    }
	    
	    _M_array[f] += materials[cell]->getChi()[e] * materials[cell]->getNuSigmaF()[g] * _mesh->getVolumes()[cell];
	    
	    if (_flux_method == ADJOINT){
	      e = g;
	      g = f;
	    }
	  }
	}
	

	if (_assemble_M){
	  petsc_err = MatSetValues(_M,1,&indice_x,_num_groups,_indices_y_M,_M_array,INSERT_VALUES);
	  CHKERRQ(petsc_err);
	}

	log_printf(DEBUG, "cell: %i, group: %i, 0: %f, 1: %f, 2: %f, 3: %f, 4: %f, 5: %f", cell, e,
		   _A_array[0], _A_array[1], _A_array[2], _A_array[3], _A_array[4], _A_array[5]); 

	petsc_err = MatSetValues(_A,1,&indice_x,_num_groups+4,_indices_y_A,_A_array,INSERT_VALUES);
	CHKERRQ(petsc_err);
	

	petsc_err = VecSetValue(_phi_old, indice_x, phi, INSERT_VALUES);
	CHKERRQ(petsc_err);
	
      }
    }
  }
  
  log_printf(INFO,"Done constructing AMPhi...");

  
  return petsc_err;
}


/* Update the MOC flux in each FSR
 * @param MOC iteration number
 */
void Cmfd::updateMOCFlux(){

  log_printf(INFO, "Updating MOC flux...");
  
  /* initialize variables */
  std::vector<int>::iterator iter;
  FP_PRECISION* old_flux = _mesh->getFluxes("old_flux");
  FP_PRECISION* new_flux = _mesh->getFluxes("new_flux");
  double old_cell_flux, new_cell_flux;
  
  /* loop over mesh cells */
  for (int i = 0; i < _cells_x*_cells_y; i++){
    
    /* loop over groups */
    for (int e = 0; e < _num_groups; e++){
      
      /* get the old and new meshCell flux */
      old_cell_flux = old_flux[i*_num_groups + e];
      new_cell_flux = new_flux[i*_num_groups + e];
      
      /* loop over FRSs in mesh cell */
      for (iter = _mesh->getCellFSRs()->at(i).begin(); iter != _mesh->getCellFSRs()->at(i).end(); ++iter) {
	
	/* set new flux in FSR */
	_FSR_fluxes[*iter*_num_groups+e] = new_cell_flux / old_cell_flux * _FSR_fluxes[*iter*_num_groups+e];
	
	log_printf(DEBUG, "Updating flux in FSR: %i, cell: %i, group: %i, ratio: %f", *iter ,i, e, new_cell_flux / old_cell_flux);
      }
    }
  }
}


/* Compute diffusion correction factors to correct diffusion coefficients
 * in optically thick mesh cells
 * @param old diffusion coefficient
 * @param height of cell
 * @return correction factor
 */
double Cmfd::computeDiffCorrect(double d, double h){

  if (_optically_thick){
    
    /* initialize variables */
    double alpha, mu, expon;
    double rho, F;
    rho = 0.0;
    
    /* loop over polar angles */
    for (int p = 0; p < 3; p++){
      mu = cos(asin(_quad->getSinTheta(p)));
      expon = exp(- h / (3 * d * mu));
      alpha = (1 + expon) / (1 - expon) - 2 * mu / h;
      rho += mu * _quad->getWeight(p) * alpha;
    }
    
    /* compute correction factor, F */
    F = 1.0 + h * rho / (2 * d);
    
    return F;
  }
  else
    return 1.0;
  
}


/* get pointer to loss matrix, A
 * @return pointer to loss matrix, A
 */
Mat Cmfd::getA(){
  return _A;
}



Mat Cmfd::getM(){
  return _M;
}



Vec Cmfd::getPhiNew(){
  return _phi_new;
}


/* get k_eff
 * @return k_eff
 */
double Cmfd::getKeff(){
  return _keff;
}


/* compute the L2 norm of consecutive fission sources
 * @return L2 norm
 */
int Cmfd::fisSourceNorm(Vec snew){


  _l2_norm = 0.0;
  int petsc_err = 0;
  

  PetscScalar *old_source;
  PetscScalar *new_source;
  petsc_err = VecGetArray(_source_old, &old_source);
  petsc_err = VecGetArray(snew, &new_source);
  CHKERRQ(petsc_err);
  
  for (int i = 0; i < _cells_x*_cells_y; i++){
    for (int e = 0; e < _num_groups; e++){
      if (new_source[i*_num_groups+e] != 0.0)
	_l2_norm += pow(new_source[i*_num_groups+e] / old_source[i*_num_groups+e] - 1.0, 2);
    }
  }
  
  _l2_norm = pow(_l2_norm, 0.5);
  

  petsc_err = VecRestoreArray(_source_old, &old_source);
  petsc_err = VecRestoreArray(snew, &new_source);
  CHKERRQ(petsc_err);
  
  return petsc_err;
}


double Cmfd::getL2Norm(){
  return _l2_norm;
}


void Cmfd::initializeFSRs(){

  log_printf(INFO, "Initialize FSRs...");
  
  /* intialize variables */
  int fsr_id;
  CellBasic* cell;
  Material* material;
  Universe* univ_zero = _geometry->getUniverse(0);
  double* heights = _mesh->getLengthsY();
  double* widths = _mesh->getLengthsX();
  
  for (int i = 0; i < _cells_x * _cells_y; i++){
    
    /* get mesh cell and fsr volume */
    fsr_id = _mesh->getCellFSRs()->at(i).front();
    _FSR_volumes[fsr_id] = heights[i / _cells_x] * widths[i % _cells_x];
    
    /* initialize the fsr fluxes to 1.0 */
    for (int e = 0; e < _num_groups; e++)
      _FSR_fluxes[fsr_id*_num_groups+e] = 1.0;
    
    /* Get the cell corresponding to this FSR from the geometry */
    cell = static_cast<CellBasic*>(_geometry->findCell(univ_zero, fsr_id));
    
    /* Get the cell's material and assign it to the FSR */
    material = _geometry->getMaterial(cell->getMaterial());
    _FSR_materials[fsr_id] = material;
    
    log_printf(DEBUG, "cell %i with FSR id = %d has cell id = %d and material id = %d "
	       "and volume = %f", i, fsr_id, cell->getId(),
	       _FSR_materials[fsr_id]->getUid(), _FSR_volumes[fsr_id]);
    
  }
  
  log_printf(INFO, "Done initializing FSRs");
}


void Cmfd::assembleM(bool assembleM){
  _assemble_M = assembleM;
}


solveType Cmfd::getSolveType(){
  return _solve_method;
}

void Cmfd::setRelaxFactor(double relax_factor){
  _relax_factor = relax_factor;
}

void Cmfd::setFSRMaterials(Material** FSR_materials){
  _FSR_materials = FSR_materials;
}


void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes){
        _FSR_volumes = FSR_volumes;

	std::vector<int>::iterator iter;
	double volume;

	/* set volume of mesh cells */
	for (int y = 0; y < _cells_y; y++){
		for (int x = 0; x < _cells_x; x++){
			volume = 0.0;

			for (iter = _mesh->getCellFSRs()->at(y*_cells_x+x).begin(); iter != _mesh->getCellFSRs()->at(y*_cells_x+x).end(); ++iter)
				volume += _FSR_volumes[*iter];

			_mesh->setVolume(volume, y*_cells_x+x);
			log_printf(DEBUG, "set cell %i volume to: %f", y*_cells_x+x, volume);
		}
	}
}

void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux){
  _FSR_fluxes = scalar_flux;
}

Mesh* Cmfd::getMesh(){
  return _mesh;
}

void Cmfd::toggleFluxType(fluxType flux_method){
  _flux_method = flux_method;
}


void Cmfd::opticallyThick(bool thick){
  _optically_thick = thick;
}


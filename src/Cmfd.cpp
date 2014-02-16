#include "Cmfd.h"

/**
 * @brief Constructor initializes boundaries and variables that describe
 *          the cmfd solver object.
 * @details The construcor initializes the many variables that describe
 *          the cmfd mesh, the solve method, and flux type. If solve
 *          method is diffusion theory, the fsr volumes, materials,
 *          and fluxes are initialized.
 * @param geometry pointer to the geometry
 * @param criteria convergence criteria on keff
 */
Cmfd::Cmfd(Geometry* geometry, double criteria) {

    /* Create objects */
    _geometry = geometry;
    _mesh = geometry->getMesh();
    _cx = _mesh->getCellsX();
    _cy = _mesh->getCellsY();
    _quad = new Quadrature(TABUCHI);
    _timer = new Timer();
    _omega = 1.0;
    
    /* Boolean and Enum flags to toggle features */
    _solve_method = _mesh->getSolveType();
    _flux_type = PRIMAL;
    _eigen_method = POWER;
    
    /* Global variables used in solving Cmfd problem */
    _l2_norm = 1.0;
    _conv_criteria = criteria;
    
    /* General problem parameters */
    _ng = _geometry->getNumEnergyGroups();
    _ncg = _ng;
    _num_fsrs = _mesh->getNumFSRs();
    _group_width = 1;

    _A = NULL;
    _M = NULL;
    _phi_temp = NULL;
    _sold = NULL;
    _snew = NULL;
        
    /* If solving diffusion problem, create arrays for FSR parameters */
    if (_solve_method == DIFFUSION){
	try{
	    _FSR_volumes = (FP_PRECISION*)calloc(_num_fsrs, sizeof(FP_PRECISION));
	    _FSR_materials = new Material*[_num_fsrs];
	    _FSR_fluxes = (FP_PRECISION*)calloc(_num_fsrs*_ng, sizeof(FP_PRECISION));
	}
	catch(std::exception &e){
	    log_printf(ERROR, "Could not allocate memory for the CMFD diffusion "
		       "mesh objects. Backtrace:%s", e.what());
	}
	
	_mesh->initializeSurfaceCurrents();
	initializeFSRs();
    }
}


/**
 * @brief Destructor deletes arrays of A and M row insertion arrays
 */
Cmfd::~Cmfd() {

    /* delete matrix and vector objects */
    if (_M != NULL){
	for (int i = 0; i < _cx*_cy; i++)
	    delete [] _M[i];

	delete [] _M;
    }

    if (_A != NULL){
	for (int i = 0; i < _cx*_cy; i++)
	    delete [] _A[i];

	delete [] _A;
    }

    if (_AM != NULL){
	for (int i = 0; i < _cx*_cy; i++)
	    delete [] _AM[i];

	delete [] _AM;
    }
    
    if (_phi_temp != NULL)
	delete [] _phi_temp;

    if (_sold != NULL)
	delete [] _sold;

    if (_snew != NULL)
	delete [] _snew;
}


/**
 * @brief Create cross sections and fluxes for each cmfd cell by
 *        energy condensing and volume averaging cross sections from
 *        the MOC sweep.
 */
void Cmfd::computeXS(){

    log_printf(INFO, "computing cmfd cross sections...");
    
    /* split corner currents to side surfaces */
    _mesh->splitCorners();
    
    /* initialize variables for FSR properties*/
    double volume, flux, abs, tot, nu_fis, chi, dif_coef;
    FP_PRECISION* scat;
    Material** materials = _mesh->getMaterials();
    double* fluxes = _mesh->getFluxes(PRIMAL);
    
    /* initialize tallies for each parameter */
    double abs_tally, nu_fis_tally, dif_tally, rxn_tally;
    double vol_tally, tot_tally, neut_prod_tally;
    double scat_tally[_ncg];
    double chi_groups[_ncg];
    double chi_tally[_ncg];
    
    /* interator to loop over fsrs in each mesh cell */
    std::vector<int>::iterator iter;
    
    /* create pointers to objects */
    Material* fsr_material;
    Material* cell_material;
    
    /* loop over mesh cells */
    #pragma omp parallel for private(volume, flux, abs, tot, nu_fis, chi, \
    	dif_coef, scat, abs_tally, nu_fis_tally, dif_tally, rxn_tally,	\
      vol_tally, tot_tally, scat_tally, iter, fsr_material, cell_material, \
    	neut_prod_tally, chi_tally, chi_groups)
    for (int i = 0; i < _cx * _cy; i++){
	
      cell_material = materials[i];

	/* loop over cmfd energy groups */
	for (int e = 0; e < _ncg; e++) {
	    
	    /* zero tallies for this group */
	    abs_tally = 0;
	    nu_fis_tally = 0;
	    dif_tally = 0;
	    rxn_tally = 0;
	    vol_tally = 0;
	    tot_tally = 0;
	    neut_prod_tally = 0.0;
    
	    /* zero each group to group scattering tally */
	    for (int g = 0; g < _ncg; g++){
		scat_tally[g] = 0;
		chi_tally[g] = 0.0;
	    }

	    /* loop over FSRs in mesh cell */
	    for (iter = _mesh->getCellFSRs()->at(i).begin(); 
		 iter != _mesh->getCellFSRs()->at(i).end(); ++iter){

	        fsr_material = _FSR_materials[*iter];
		volume = _FSR_volumes[*iter];
		scat = fsr_material->getSigmaS();
		vol_tally += volume;

		/* chi tallies */
		for (int b = 0; b < _ncg; b++){
		    chi_groups[b] = 0.0;
		    for (int h = _group_indices[b]; h < _group_indices[b+1]; h++)
		        chi_groups[b] += fsr_material->getChi()[h];

		    for (int h = 0; h < _ng; h++){
		        chi_tally[b] += chi_groups[b] * fsr_material->getNuSigmaF()[h] * 
			  _FSR_fluxes[(*iter)*_ng+h] * volume;
			neut_prod_tally += chi_groups[b] * 
			  fsr_material->getNuSigmaF()[h] * 
			  _FSR_fluxes[(*iter)*_ng+h] * volume;
		    }
		}

		/* loop over groups within this cmfd group */
		for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){
		  
		    /* Gets FSR volume, material, and cross sections */
		    flux = _FSR_fluxes[(*iter)*_ng+h];
		    abs = fsr_material->getSigmaA()[h];
		    tot = fsr_material->getSigmaT()[h];
		    dif_coef = fsr_material->getDifCoef()[h];
		    nu_fis = fsr_material->getNuSigmaF()[h];
		    
		    /* if material has a diffusion coefficient, use it; otherwise
		     * estimate the diffusion coefficient with 1 / (3 * sigma_t) */
		    if (fsr_material->getDifCoef()[h] > 1e-8){
		        dif_tally += fsr_material->getDifCoef()[h] * flux * volume;
		    }
		    else
		        dif_tally += flux * volume / (3.0 * tot);
		    
		    /* increment tallies for this group */
		    abs_tally += abs * flux * volume;
		    tot_tally += tot * flux * volume;
		    nu_fis_tally += nu_fis * flux * volume;
		    rxn_tally += flux * volume;
		    
		    /* scattering tallies */
		    for (int g = 0; g < _ng; g++){
		        scat_tally[std::min(g / _group_width, _ncg-1)] += 
			  scat[g*_ng+h] * flux * volume;
		    }
		}
	    }
	    
	    /* set the mesh cell properties with the tallies */
	    _mesh->setVolume(vol_tally, i);
	    cell_material->setSigmaAByGroup(abs_tally / rxn_tally, e);
	    cell_material->setSigmaTByGroup(tot_tally / rxn_tally, e);
	    cell_material->setNuSigmaFByGroup(nu_fis_tally / rxn_tally, e);
	    cell_material->setDifCoefByGroup(dif_tally / rxn_tally, e);
	    fluxes[i*_ncg+e] = rxn_tally / vol_tally;

	    /* set chi */
	    if (neut_prod_tally != 0.0)
	      cell_material->setChiByGroup(chi_tally[e] / neut_prod_tally,e);	    
	    else
	      cell_material->setChiByGroup(0.0,e);	    

	    log_printf(DEBUG, "cell: %i, group: %i, vol: %f, siga: %f, sigt: %f,"
		       " nu_sigf: %f, dif_coef: %f, flux: %f, chi: %f", i, e, vol_tally, 
		       abs_tally / rxn_tally, tot_tally / rxn_tally, 
		       nu_fis_tally / rxn_tally, dif_tally / rxn_tally, 
		       rxn_tally / vol_tally, chi_tally[e] / (neut_prod_tally+1e-12));
	    
	    for (int g = 0; g < _ncg; g++){
		cell_material->setSigmaSByGroup(scat_tally[g] / rxn_tally, g, e);
		log_printf(DEBUG, "scattering from %i to %i: %f", e, g, 
			   scat_tally[g] / rxn_tally);
	    }
	}
    }
}


/**
 * @brief Compute the diffusion coefficients (d_dif - straight 
 *        diffusion coefficient, d_hat - surface diffusion coefficient, 
 *        and d_tilde - surface diffusion coefficient correction factor)
 *        for each mesh while ensuring neutron balance is achieved.
 */
void Cmfd::computeDs(){

    log_printf(INFO, "computing cmfd Ds...");
    
    /* initialize variables */
    double d, d_next, d_hat, d_tilde;
    double current, flux, flux_next, f, f_next;
    double length, length_perpen, next_length_perpen;
    double sense;
    int next_surface;
    int cell, cell_next;
    
    Material** materials = _mesh->getMaterials();
    double* cell_flux = _mesh->getFluxes(PRIMAL);
    double* lengths_y = _mesh->getLengthsY();
    double* lengths_x = _mesh->getLengthsX();
    double* currents = _mesh->getCurrents();
    
    /* loop over mesh cells in y direction */
    #pragma omp parallel for private(d, d_next, d_hat, d_tilde, current, flux, \
    	flux_next, f, f_next, length, length_perpen, next_length_perpen, \
    	sense, next_surface, cell, cell_next)
    for (int y = 0; y < _cy; y++){
	
	/* loop over mesh cells in x direction */
	for (int x = 0; x < _cx; x++){
	    
	    cell = y*_cx+x;
	    
	    /* loop over surfaces in a cell */
	    for (int surface = 0; surface < 4; surface++){
		
		/* loop over groups */
		for (int e = 0; e < _ncg; e++){
		    
		    /* get diffusivity and flux for mesh cell */
		    d = materials[cell]->getDifCoef()[e];
		    flux = cell_flux[cell*_ncg+e];
		    cell_next = _mesh->getCellNext(cell, surface);
		    
		    /* set sense of the surface */
		    if (surface == 0 || surface == 3)
			sense = -1.0;
		    else
			sense = 1.0;
		    
		    /* set the length of this surface and the
		     * perpendicular surface */
		    if (surface == 0 || surface== 2){
			length = lengths_y[y];
			length_perpen = lengths_x[x];
		    }
		    else if (surface == 1 || surface == 3){
			length = lengths_x[x];
			length_perpen = lengths_y[y];
		    }
		    
		    /* compute the optical thickness correction factor */
		    f = computeDiffCorrect(d, length_perpen);
		    
		    /* if surface is on a boundary, choose appropriate BCs */
		    if (cell_next == -1){
			
			current = sense * currents[cell*_ncg*8 + surface*_ncg + e];

			/* REFLECTIVE BC */
			if (_mesh->getBoundary(surface) == REFLECTIVE){ 
			    
			    /* set d's */ 
			    d_hat = 0.0;
			    d_tilde = 0.0;
			}
			/* VACUUM BC */
			else if (_mesh->getBoundary(surface) == VACUUM){	      
			    
			    /* set d's */
			    d_hat =  2 * d*f / length_perpen / (1 + 4 * d*f / 
								length_perpen);
			    
			    if (_solve_method == MOC)
				d_tilde = (sense * d_hat * flux - current / 
					   length) / flux;
			    else
				d_tilde = 0.0;
			}
			/* ZERO_FLUX BC */
			else if (_mesh->getBoundary(surface) == ZERO_FLUX){
			    
			    /* set d's */
			    d_hat = 2 * d*f / length_perpen;
			    d_tilde = 0.0;
			}
		    }
		    /* if surface is an interface, use finite differencing */
		    else{
			
			/* set properties for cell next to surface */
			if (surface == 0){
			    next_length_perpen = lengths_x[cell_next % _cx];
			    next_surface = 2;
			}
			else if (surface == 1){
			    next_length_perpen = lengths_y[cell_next / _cx];
			    next_surface = 3;
			}
			else if (surface == 2){
			    next_length_perpen = lengths_x[cell_next % _cx];
			    next_surface = 0;
			}
			else if (surface == 3){
			    next_length_perpen = lengths_y[cell_next / _cx];
			    next_surface = 1;
			}
			
			/* set diffusion coefficient and flux for neighboring
			 * cell */
			d_next = materials[cell_next]->getDifCoef()[e];
			flux_next = cell_flux[cell_next*_ncg + e];
			
			/* get optical thickness correction term for 
			 * meshCellNext */
			f_next = computeDiffCorrect(d_next, next_length_perpen);
			
			/* compute d_hat */
			d_hat = 2.0 * d * f * d_next * f_next / (length_perpen 
			        * d * f + next_length_perpen * d_next*f_next);
			
			/* compute net current */
			current = sense * currents[cell*_ncg*8 + surface*_ncg + e]
			    - sense * currents[cell_next*_ncg*8 + 
					       next_surface*_ncg + e];
			
			/* compute d_tilde */
			if (_solve_method == MOC)
			    d_tilde = -(sense * d_hat * (flux_next - flux) + 
			    current  / length) / (flux_next + flux);
			else
			    d_tilde = 0.0;
			
			/* if the magnitude of d_tilde is greater than the 
			 * magnitude of d_hat, select new values d_tilde 
			 * and d_hat to ensure the course mesh equations
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
		    }  
		    
		    /* perform underrelaxation on d_tilde */
		    d_tilde = materials[cell]->getDifTilde()[surface*_ncg + e] * 
			(1 - _mesh->getRelaxFactor()) + _mesh->getRelaxFactor()
			* d_tilde;
		    
		    /* set d_hat and d_tilde */
		    materials[cell]->setDifHatByGroup(d_hat, e, surface);
		    materials[cell]->setDifTildeByGroup(d_tilde, e, surface);
		    
		    log_printf(DEBUG, "cell: %i, group: %i, side: %i, flux: %f,"
			       " current: %f, d: %f, dhat: %f, dtilde: %f", 
			       y*_cx + x, e, surface, flux, current, d, d_hat, 
			       d_tilde);
		}
	    }
	}
    }
}


/*
 * @brief CMFD solver that solves the diffusion problem
 * @return k-effective
 */
double Cmfd::computeKeff(){

    log_printf(INFO, "Running diffusion solver...");

    /* Create matrix and vector objects */
    if (_A == NULL){
        try{
	    _AM = NULL;
	    _M = new double*[_cx*_cy];
	    _A = new double*[_cx*_cy];
	    _phi_temp = new double[_cx*_cy*_ncg];
	    _sold = new double[_cx*_cy*_ncg];
	    _snew = new double[_cx*_cy*_ncg];
	    _mesh->setNumGroups(_ncg);
	    createGroupStructure();
	    _mesh->initializeFlux();

	    
	    if (_solve_method == MOC)
	        _mesh->initializeMaterialsMOC();

	    for (int i = 0; i < _cx*_cy; i++){
	        _M[i] = new double[_ncg*_ncg];
		_A[i] = new double[_ncg*(_ncg+4)];
	    }
	}
	catch(std::exception &e){
	    log_printf(ERROR, "Could not allocate memory for the CMFD mesh objects. "
		       "Backtrace:%s", e.what());
	}
    }   
    
    /* if solving diffusion problem, initialize timer */
    if (_solve_method == DIFFUSION)
	_timer->startTimer();
    
    /* initialize variables */
    double sumnew, sumold, val, norm, scale_val;
    double flux_conv = 1e-8;
    int row;
    double* phi_old = _mesh->getFluxes(PRIMAL);
    double* phi_new = _mesh->getFluxes(PRIMAL_UPDATE);

    /* compute the cross sections and surface 
     * diffusion coefficients */
    if (_solve_method == MOC)
	computeXS();
    
    computeDs();

    /* construct matrices */
    constructMatrices();
    
    if (_eigen_method == POWER){

	/* compute and normalize the initial source */
	matMultM(_M, phi_old, _sold);
	sumold = vecSum(_sold);
	scale_val = (_cx * _cy * _ncg) / sumold;
	vecScale(_sold, scale_val);
	vecCopy(phi_old, phi_new);
	vecScale(phi_new, scale_val);
	sumold = _cx * _cy * _ncg;
	
	/* power iteration diffusion solver */
	for (int iter = 0; iter < 25000; iter++){
	    
	    /* Solve phi = A^-1 * old_source */
	    linearSolve(_A, phi_new, _sold, flux_conv);
	    
	    /* compute the new source */
	    matMultM(_M, phi_new, _snew);
	    sumnew = vecSum(_snew);
	    
	    /* compute and set keff */
	    _k_eff = sumnew / sumold;
	    
	    /* scale the old source by keff */
	    vecScale(_sold, _k_eff);
	    
	    /* compute the L2 norm of source error */
	    norm = 0.0;
	    for (int i = 0; i < _cx*_cy*_ncg; i++){
		if (_snew[i] != 0.0)
		    norm += pow((_snew[i] - _sold[i]) / _snew[i], 2);
	    }
	    
	    norm = sqrt(norm / (_cx*_cy*_ncg));
	    
	    scale_val = (_cx * _cy * _ncg) / sumnew;
	    vecScale(_snew, scale_val);
	    vecCopy(_snew, _sold);
	    
	    log_printf(INFO, "GS POWER iter: %i, keff: %f, error: %f", iter, 
		       _k_eff, norm);
	    
	    /* check for convergence */
	    if (norm < _conv_criteria)
		break;
	}
    }	
    else{

	/* allocate memory for AM matrix */
	if (_AM == NULL){
	    log_printf(INFO, "Allocating memory for AM");
	    try{
		_AM = new double*[_cx*_cy];

		for (int i = 0; i < _cx*_cy; i++){
		    _AM[i] = new double[_ncg*(_ncg+4)];
		}
	    }
	    catch(std::exception &e){
		log_printf(ERROR, "Could not allocate memory for the _AM matrix."
			   " Backtrace:%s", e.what());
	    }
	    log_printf(INFO, "Done allocating memory for AM");
	}

	int max_iter = 100;
	double shift, offset;
	int iter = 0;
	norm = 1.0;
	
	/*copy old flux to new */
	vecCopy(phi_old, phi_new);

	/* compute the initial k_eff */
	_k_eff = rayleighQuotient(phi_new, _snew, _phi_temp);

	/* compute and normalize the initial source */
	matMultM(_M, phi_new, _sold);
	sumold = vecSum(_sold);
	scale_val = (_cx * _cy * _ncg) / sumold;
	vecScale(_sold, scale_val);
	vecScale(phi_new, scale_val);
	sumold = _cx * _cy * _ncg;

	shift = 1.5;

	while (norm > 1e-5){

	    /* reconstruct _AM */
	    matSubtract(_AM, _A, 1.0/shift, _M);

	    /* solve inverse system */
	    linearSolve(_AM, phi_new, _sold, flux_conv);
	    
	    /* compute new flux */
	    vecScale(phi_new, vecMax(phi_new));
	    matMultM(_M, phi_new, _snew);
	    sumnew = vecSum(_snew);
	    vecScale(_snew, (_cx*_cy*_ncg) / sumnew);
	    vecScale(phi_new, (_cx*_cy*_ncg) / sumnew);

	    /* compute new eigenvalue */
	    _k_eff = rayleighQuotient(phi_new, _snew, _phi_temp);

	    /* compute the L2 norm of source error */
	    norm = 0.0;
	    for (int i = 0; i < _cx*_cy*_ncg; i++){
		if (_snew[i] != 0.0)
		    norm += pow((_snew[i] - _sold[i]) / _snew[i], 2);
	    }
	    
	    norm = pow(norm, 0.5);
	    norm = norm / (_cx*_cy*_ncg);
	    vecCopy(_snew, _sold);

	    iter++;

	    log_printf(INFO, "iter: %i, k_eff: %f, norm: %f", iter, 
		       _k_eff, norm);
	}

	offset = 0.05;
	log_printf(INFO, "offset set to: %f", offset);

	while (norm > _conv_criteria){
	    
	    /* reconstruct _AM */
	    matSubtract(_AM, _A, 1.0/(_k_eff + offset), _M);

	    /* solve inverse system */
	    linearSolve(_AM, phi_new, _sold, flux_conv);
	    
	    /* compute new flux */
	    vecScale(phi_new, vecMax(phi_new));
	    matMultM(_M, phi_new, _snew);
	    sumnew = vecSum(_snew);
	    vecScale(_snew, (_cx*_cy*_ncg) / sumnew);
	    vecScale(phi_new, (_cx*_cy*_ncg) / sumnew);

	    /* compute new eigenvalue */
	    _k_eff = rayleighQuotient(phi_new, _snew, _phi_temp);
	    
	    /* compute the L2 norm of source error */
	    norm = 0.0;
	    for (int i = 0; i < _cx*_cy*_ncg; i++){
		if (_snew[i] != 0.0)
		    norm += pow((_snew[i] - _sold[i]) / _snew[i], 2);
	    }
	    
	    norm = pow(norm, 0.5);
	    norm = norm / (_cx*_cy*_ncg);
	    vecCopy(_snew, _sold);

	    iter++;

	    log_printf(INFO, "iter: %i, k_eff: %f, norm: %f", iter, 
		       _k_eff, norm);
	}
    }
    
    /* rescale the old and new flux */
    rescaleFlux();
    
    /* update the MOC flux */
    if (_solve_method == MOC)
	updateMOCFlux();  
    
    if (_flux_type == ADJOINT)
	vecCopy(phi_new, _mesh->getFluxes(ADJOINT));
    
    /* If solving diffusion problem, print timing results */
    if (_solve_method == DIFFUSION){
	std::string msg_string;
	log_printf(TITLE, "TIMING REPORT");
	_timer->stopTimer();
	_timer->recordSplit("Total time to solve diffusion eigenvalue problem");
	
	double tot_time = _timer->getSplit("Total time to solve diffusion "
					   "eigenvalue problem");
	msg_string = "Total time to solve diffusion eigenvalue problem";
	msg_string.resize(53, '.');
	log_printf(RESULT, "%s%1.4E sec", msg_string.c_str(), tot_time);
    }
    
    return _k_eff;
}


/**
 * @brief solve the linear system Ax=b using Gauss Seidel
 *        with SOR
 * @param mat pointer to A matrix
 * @param vec_x pointer to x vector
 * @param vec_b pointer to b vector
 * @param conv flux convergence criteria
 */
void Cmfd::linearSolve(double** mat, double* vec_x, double* vec_b, double conv,
    int max_iter){

    /* initialize variables */
    double norm = 1e10;
    int row, cell;
    double val;
    int iter = 0;
    
    while (norm > conv){

	/* pass new flux to old flux */
	vecCopy(vec_x, _phi_temp);

	/* iteration over red cells */
	#pragma omp parallel for private(row, val, cell)
	for (int y = 0; y < _cy; y++){
	    for (int x = y % 2; x < _cx; x += 2){
		cell = y*_cx+x;
		for (int g = 0; g < _ncg; g++){
		    row = (y*_cx+x)*_ncg + g;
		    val = 0.0;
		    
		    /* previous flux term */
		    val += (1.0 - _omega) * vec_x[row];
		    
		    /* source term */
		    val += _omega * vec_b[row] / mat[cell][g*(_ncg+4)+g+2];
		    
		    
		    /* left surface */
		    if (x != 0)
			val -= _omega * vec_x[row - _ncg] * mat[cell][g*(_ncg+4)] /
			    mat[cell][g*(_ncg+4)+g+2];
		    
		    /* bottom surface */
		    if (y != _cy - 1)
			val -= _omega * vec_x[row + _cx * _ncg] * 
			    mat[cell][g*(_ncg+4)+1] / mat[cell][g*(_ncg+4)+g+2];
		    
		    /* group to group */
		    for (int e = 0; e < _ncg; e++){
			if (e != g)
			    val -= _omega * vec_x[(y*_cx+x)*_ncg+e] * 
				mat[cell][g*(_ncg+4)+2+e] / mat[cell][g*(_ncg+4)+g+2];
		    }
		    
		    /* right surface */
		    if (x != _cx - 1)
			val -= _omega * vec_x[row + _ncg] * 
			    mat[cell][g*(_ncg+4)+_ncg+2] / mat[cell][g*(_ncg+4)+g+2];
		     
		    /* top surface */
		    if (y != 0)
			val -= _omega * vec_x[row - _ncg*_cx] * 
			    mat[cell][g*(_ncg+4)+_ncg+3] / mat[cell][g*(_ncg+4)+g+2];
		    
		    vec_x[row] = val;
		}
	    }
	}

	/* iteration over black cells */
        #pragma omp parallel for private(row, val, cell)
	for (int y = 0; y < _cy; y++){
	    for (int x = 1 - y % 2; x < _cx; x += 2){
		cell = y*_cx+x;
		for (int g = 0; g < _ncg; g++){
		    row = (y*_cx+x)*_ncg + g;
		    val = 0.0;
		    
		    /* previous flux term */
		    val += (1.0 - _omega) * vec_x[row];
		    
		    /* source term */
		    val += _omega * vec_b[row] / mat[cell][g*(_ncg+4)+g+2];
		    
		    
		    /* left surface */
		    if (x != 0)
			val -= _omega * vec_x[row - _ncg] * 
			    mat[cell][g*(_ncg+4)] / mat[cell][g*(_ncg+4)+g+2];
		    
		    /* bottom surface */
		    if (y != _cy - 1)
			val -= _omega * vec_x[row + _cx * _ncg] * 
			    mat[cell][g*(_ncg+4)+1] / mat[cell][g*(_ncg+4)+g+2];
		    
		    /* group to group */
		    for (int e = 0; e < _ncg; e++){
			if (e != g)
			    val -= _omega * vec_x[(y*_cx+x)*_ncg+e] * 
				mat[cell][g*(_ncg+4)+2+e] / mat[cell][g*(_ncg+4)+g+2];
		    }
		    
		    /* right surface */
		    if (x != _cx - 1)
			val -= _omega * vec_x[row + _ncg] * 
			    mat[cell][g*(_ncg+4)+_ncg+2] / mat[cell][g*(_ncg+4)+g+2];
		    
		    /* top surface */
		    if (y != 0)
			val -= _omega * vec_x[row - _ncg*_cx] * 
			    mat[cell][g*(_ncg+4)+_ncg+3] / mat[cell][g*(_ncg+4)+g+2];
		    
		    vec_x[row] = val;
		}
	    }
	}
	
	norm = 0.0;
	for (int i = 0; i < _cx*_cy*_ncg; i++){
	    if (vec_x[i] != 0.0)
		norm += pow((vec_x[i] - _phi_temp[i]) / vec_x[i], 2);	
        }
	
	norm = pow(norm, 0.5) / (_cx*_cy*_ncg);

	iter++;

	log_printf(DEBUG, "GS iter: %i, norm: %f", iter, norm);

	if (iter >= max_iter)
	    break;
    }
    
    log_printf(DEBUG, "linear solver iterations: %i", iter);
}


/**
 * @brief rescale the initial and converged flux arrays
 * @return petsc_err petsc error flag
 */
void Cmfd::rescaleFlux(){
    
    double sumnew, sumold, scale_val;
    double* phi_old = _mesh->getFluxes(PRIMAL);
    double* phi_new = _mesh->getFluxes(PRIMAL_UPDATE);
    
    /* rescale the new and old flux to have an avg source of 1.0 */
    matMultM(_M, phi_new, _snew);
    sumnew = vecSum(_snew);
    scale_val = _cx*_cy*_ncg / sumnew;
    vecScale(phi_new, scale_val);
    matMultM(_M, phi_old, _sold);
    sumold = vecSum(_sold);
    scale_val = _cx*_cy*_ncg / sumold;
    vecScale(phi_old, scale_val);
}


/**
 * @brief scale vector
 * @param vec vector to be scaled
 * @param scale_val value to scale vector
 */
void Cmfd::vecScale(double* vec, double scale_val){

    #pragma omp parallel for
    for (int i = 0; i < _cx*_cy*_ncg; i++)
	vec[i] *= scale_val;

}


/**
 * @brief set vector to value
 * @param vec vector to be set
 * @param val value to set vector to
 */
void Cmfd::vecSet(double* vec, double val){

    #pragma omp parallel for
    for (int i = 0; i < _cx*_cy*_ncg; i++)
	vec[i] = val;

}


/**
 * @brief normalize vector to have avg source of 1.0
 * @param mat source matrix
 * @param vec vector to be normalized
 */
void Cmfd::vecNormal(double** mat, double* vec){
 
    double source, scale_val;
    matMultM(mat, vec, _phi_temp);
    source = vecSum(_phi_temp);
    scale_val = (_cx*_cy*_ncg) / source;
    vecScale(vec, scale_val);
    
}


/**
 * @brief multiply matrix by vector, y = M *x
 * @param mat source matrix
 * @param vec_x x vector
 * @param vec_y y vector
 */
void Cmfd::matMultM(double** mat, double* vec_x, double* vec_y){

    vecSet(vec_y, 0.0);
    
    for (int i = 0; i < _cx*_cy; i++){
	for (int g = 0; g < _ncg; g++){
	    for (int e = 0; e < _ncg; e++){
		vec_y[i*_ncg+g] += mat[i][g*_ncg+e] * vec_x[i*_ncg+e];
	    }
	}
    }
}

/**
 * @brief sum vector
 * @param vec vector to be summed
 * @return sum sum of vector
 */
double Cmfd::vecSum(double* vec){

    double sum = 0.0;
    
    for (int i = 0; i < _cx*_cy*_ncg; i++)
	sum += vec[i];
    
    return sum;
}


/**
 * @brief copy vector
 * @param vec_from vector to be copied
 * @return vec_to vector to receive copy
 */
void Cmfd::vecCopy(double* vec_from, double* vec_to){
    
    #pragma omp parallel for
    for (int i = 0; i < _cx*_cy*_ncg; i++)
	vec_to[i] = vec_from[i];
}


/**
 * @brief zero matrix
 * @param mat matrix to be zeroed
 * @param width width of matrix row
 */
void Cmfd::matZero(double** mat, int width){
    
    #pragma omp parallel for
    for (int i = 0; i < _cx*_cy; i++){
	for (int g = 0; g < _ncg*width; g++)
	    mat[i][g] = 0.0;
    }
}


/* Fill in the values in the A matrix, M matrix, and phi_old vector
 * @return petsc_err petsc error flag
 */
void Cmfd::constructMatrices(){
    
    log_printf(INFO,"Constructing matrices...");
    
    /* initialize variables */
    double value, volume;
    int cell, row;
    Material* material;
    
    /* get arrays */
    double* heights = _mesh->getLengthsY();
    double* widths = _mesh->getLengthsX();

    /* zero _A and _M matrices */
    matZero(_M, _ncg);
    matZero(_A, _ncg+4);
    
    /* loop over cells */
    #pragma omp parallel for private(value, volume, cell, row, material)
    for (int y = 0; y < _cy; y++){
	for (int x = 0; x < _cx; x++){
	    
	    cell = y*_cx + x;
	    material = _mesh->getMaterials()[cell];
	    volume = _mesh->getVolumes()[cell];

	    /* loop over groups */
	    for (int e = 0; e < _ncg; e++){
		
		row = cell*_ncg + e;
		
		/* absorption term */
		value = material->getSigmaA()[e] * volume;
		_A[cell][e*(_ncg+4)+e+2] += value;
		
		/* out (1st) and int (2nd) scattering */
		if (_flux_type == PRIMAL){
		    for (int g = 0; g < _ncg; g++){
			if (e != g){
			    value = material->getSigmaS()[g*_ncg + e] * volume;
			    _A[cell][e*(_ncg+4)+e+2] += value;
			    value = - material->getSigmaS()[e*_ncg + g] * volume;
			    _A[cell][e*(_ncg+4)+g+2] += value;
			}
		    }
		}
		else{
		    for (int g = 0; g < _ncg; g++){
			if (e != g){
			    value = material->getSigmaS()[e*_ncg + g] * volume;
			    _A[cell][e*(_ncg+4)+e+2] += value;
			    value = - material->getSigmaS()[g*_ncg + e] * volume;
			    _A[cell][e*(_ncg+4)+g+2] += value;
			}
		    }
		}
		
		/* RIGHT SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[2*_ncg + e] 
			 - material->getDifTilde()[2*_ncg + e]) 
		    * heights[cell / _cx];
		
		_A[cell][e*(_ncg+4)+e+2] += value;
		
		/* set transport term on off diagonal */
		if (x != _cx - 1){
		    value = - (material->getDifHat()[2*_ncg + e] 
			       + material->getDifTilde()[2*_ncg + e]) 
			* heights[cell / _cx];
		    
		    _A[cell][e*(_ncg+4)+_ncg+2] += value;
		}
		
		/* LEFT SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[0*_ncg + e] 
			 + material->getDifTilde()[0*_ncg + e]) 
		    * heights[cell / _cx];
		
		_A[cell][e*(_ncg+4)+e+2] += value;
		
		/* set transport term on off diagonal */
		if (x != 0){
		    value = - (material->getDifHat()[0*_ncg + e] 
			       - material->getDifTilde()[0*_ncg + e]) 
			* heights[cell / _cx];
		    
		    _A[cell][e*(_ncg+4)] += value;
		}
		
		/* BOTTOM SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[1*_ncg + e] 
			 - material->getDifTilde()[1*_ncg + e]) 
		    * widths[cell % _cx];
		
		_A[cell][e*(_ncg+4)+e+2] += value;
		
		/* set transport term on off diagonal */
		if (y != _cy - 1){
		    value = - (material->getDifHat()[1*_ncg + e] 
			       + material->getDifTilde()[1*_ncg + e]) 
			* widths[cell % _cx];
		    
		    _A[cell][e*(_ncg+4)+1] += value;
		}
		
		/* TOP SURFACE */
		
		/* set transport term on diagonal */
		value = (material->getDifHat()[3*_ncg + e] 
			 + material->getDifTilde()[3*_ncg + e]) 
		    * widths[cell % _cx];
		
		_A[cell][e*(_ncg+4)+e+2] += value;
		
		/* set transport term on off diagonal */
		if (y != 0){
		    value = - (material->getDifHat()[3*_ncg + e] 
		     - material->getDifTilde()[3*_ncg + e]) 
			* widths[cell % _cx];
		    
		    _A[cell][e*(_ncg+4)+_ncg+3] += value;
		}
		
		/* source term */
		for (int g = 0; g < _ncg; g++){	 
		    value = material->getChi()[e] * material->getNuSigmaF()[g] 
			* volume;
		    
		    if (_flux_type == PRIMAL)
			_M[cell][e*_ncg+g] += value;
		    else
			_M[cell][g*_ncg+e] += value;
		}
		
		log_printf(DEBUG, "cel: %i, vol; %f", cell, 
			   _mesh->getVolumes()[cell]);
		
		for (int i = 0; i < _ncg+4; i++)
		    log_printf(DEBUG, "i: %i, A value: %f", i, 
			       _A[cell][e*(_ncg+4)+i]);
		
		for (int i = 0; i < _ncg; i++)
		    log_printf(DEBUG, "i: %i, M value: %f", i, 
			       _M[cell][e*(_ncg+4)+i]);
		
	    }
	}
    }
    
    
    log_printf(INFO,"Done constructing matrices...");
}


/**
 * @brief Update the MOC flux in each FSR
 */
void Cmfd::updateMOCFlux(){
    
    log_printf(INFO, "Updating MOC flux...");
    
    /* initialize variables */
    std::vector<int>::iterator iter;
    double* old_flux = _mesh->getFluxes(PRIMAL);
    double* new_flux = _mesh->getFluxes(PRIMAL_UPDATE);
    double old_cell_flux, new_cell_flux;
    
    /* loop over mesh cells */
    #pragma omp parallel for private(iter, old_cell_flux, new_cell_flux)
    for (int i = 0; i < _cx*_cy; i++){
	
	/* loop over cmfd groups */
	for (int e = 0; e < _ncg; e++){
	    
	    /* get the old and new meshCell flux */
	    old_cell_flux = old_flux[i*_ncg + e];
	    new_cell_flux = new_flux[i*_ncg + e];
	    
	    for (int h = _group_indices[e]; h < _group_indices[e+1]; h++){

	        /* loop over FRSs in mesh cell */
	        for (iter = _mesh->getCellFSRs()->at(i).begin(); 
		     iter != _mesh->getCellFSRs()->at(i).end(); ++iter) {
		
		    /* set new flux in FSR */
		    _FSR_fluxes[*iter*_ng+h] = new_cell_flux / 
		        old_cell_flux * _FSR_fluxes[*iter*_ng+h];
		
		    log_printf(DEBUG, "Updating flux in FSR: %i, cell: %i, group: "
			       "%i, ratio: %f", *iter ,i, h, 
			       new_cell_flux / old_cell_flux);
		}
	    }
	}
    }
}


/**
 * @brief Compute diffusion correction factors to correct 
 * diffusion coefficients in optically thick mesh cells
 * @param d old diffusion coefficient
 * @param h height of cell
 * @return f correction factor
 */
double Cmfd::computeDiffCorrect(double d, double h){

    if (_mesh->getOpticallyThick() && _solve_method == MOC){
	
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


/**
 * @brief get k_eff
 * @return _k_eff k_eff
 */
double Cmfd::getKeff(){
    return _k_eff;
}


/**
 * @brief initialize the fsrs 
 */
void Cmfd::initializeFSRs(){

    log_printf(INFO, "Initialize FSRs...");
    
    /* intialize variables */
    int fsr_id;
    CellBasic* cell;
    Material* material;
    Universe* univ_zero = _geometry->getUniverse(0);
    double* heights = _mesh->getLengthsY();
    double* widths = _mesh->getLengthsX();
    
    for (int i = 0; i < _cx * _cy; i++){
	
	/* get mesh cell and fsr volume */
	fsr_id = _mesh->getCellFSRs()->at(i).front();
	_FSR_volumes[fsr_id] = heights[i / _cx] * widths[i % _cx];
	
	/* initialize the fsr fluxes to 1.0 */
	for (int e = 0; e < _ng; e++)
	    _FSR_fluxes[fsr_id*_ng+e] = 1.0;
	
	/* Get the cell corresponding to this FSR from the geometry */
	cell = _geometry->findCellContainingFSR(fsr_id);
	
	/* Get the cell's material and assign it to the FSR */
	material = _geometry->getMaterial(cell->getMaterial());
	_FSR_materials[fsr_id] = material;
	
	log_printf(DEBUG, "cell %i with FSR id = %d has cell id = %d and material id = %d "
		   "and volume = %f", i, fsr_id, cell->getId(),
		   _FSR_materials[fsr_id]->getUid(), _FSR_volumes[fsr_id]);
	
    }
    
    log_printf(INFO, "Done initializing FSRs");
}


/**
 * @brief Set the fsr materials array pointer
 * @param FSR_materials pointer to fsr materials array
 */
void Cmfd::setFSRMaterials(Material** FSR_materials){
    _FSR_materials = FSR_materials;
}


/**
 * @brief Set the fsr volumes by summing the volumes of 
 *        the fsrs contained in each cell
 * @param FSR_volumes array of fsr volumes
 */
void Cmfd::setFSRVolumes(FP_PRECISION* FSR_volumes){
    _FSR_volumes = FSR_volumes;
    
    std::vector<int>::iterator iter;
    double volume;
    
    /* set volume of mesh cells */
    for (int y = 0; y < _cy; y++){
	for (int x = 0; x < _cx; x++){
	    volume = 0.0;
	    
	    for (iter = _mesh->getCellFSRs()->at(y*_cx+x).begin(); iter != _mesh->getCellFSRs()->at(y*_cx+x).end(); ++iter)
		volume += _FSR_volumes[*iter];
	    
	    _mesh->setVolume(volume, y*_cx+x);
	    log_printf(DEBUG, "set cell %i volume to: %f", y*_cx+x, volume);
	}
    }
}


/**
 * @brief Set pointer to fsr flux array
 * @param scalar_flux pointer to fsr flux array
 */
void Cmfd::setFSRFluxes(FP_PRECISION* scalar_flux){
    _FSR_fluxes = scalar_flux;
}


/**
 * @brief Get pointer to the mesh 
 * @return _mesh pointer to mesh
 */
Mesh* Cmfd::getMesh(){
    return _mesh;
}


/**
 * @brief set the flux type
 * @param flux_type char string representing enum for flux type
 */
void Cmfd::setFluxType(const char* flux_type){
    
    if (strcmp("PRIMAL", flux_type) == 0)
	_flux_type = PRIMAL;
    else if (strcmp("ADJOINT", flux_type) == 0)
	_flux_type = ADJOINT;
    else
	log_printf(ERROR, "Could not recognize flux type: "
		   " the options are PRIMAL and ADJOINT");
}


/**
 * @brief set the eigenvalue method
 * @param flux_type char string representing enum for eigen method
 */
void Cmfd::setEigenMethod(const char* eigen_method){
    
    if (strcmp("POWER", eigen_method) == 0)
	_eigen_method = POWER;
    else if (strcmp("WIELANDT", eigen_method) == 0)
	_eigen_method = WIELANDT;
    else
	log_printf(ERROR, "Could not recognize eigen method: "
		   " the options are POWER and WIELANDT");
}


/**
 * @brief dump a vector to screen
 * @param vec vector to be dumped
 * @param length length of vector
 */
void Cmfd::dumpVec(double* vec, int length){

    log_printf(NORMAL, "dumping vector...");
    
    for (int i = 0; i < length; i++)
	log_printf(NORMAL, "cell: %i, value: %f", i, vec[i]);
    
    log_printf(NORMAL, "done dumping vector...");
    
}


/**
 * @brief set the SOR factor
 */
void Cmfd::setOmega(double omega){
    _omega = omega;
}


double Cmfd::rayleighQuotient(double* x, double* snew, double* sold){

    double numer = 0.0;
    double denom = 0.0;

    matMultA(_A, x, sold);
    matMultM(_M, x, snew);

    for (int i = 0; i < _cx*_cy*_ncg; i++){
	numer += x[i]*snew[i];
	denom += x[i]*sold[i];
    }

    return numer/denom;
}


/**
 * @brief multiply matrix by vector, y = M *x
 * @param mat source matrix
 * @param vec_x x vector
 * @param vec_y y vector
 */
void Cmfd::matMultA(double** mat, double* vec_x, double* vec_y){

    vecSet(vec_y, 0.0);
    int row, cell;
    
    for (int y = 0; y < _cy; y++){
        for (int x = 0; x < _cx; x++){
            cell = y*_cx+x;
            for (int g = 0; g < _ncg; g++){
                row = cell*_ncg + g;
                
		if (x != 0)
		    vec_y[row] += mat[cell][g*(_ncg+4)] * 
			vec_x[(cell-1)*_ncg+g];

		if (y != _cy - 1)
		    vec_y[row] += mat[cell][g*(_ncg+4)+1] * 
			vec_x[(cell+_cx)*_ncg+g];

		if (x != _cx - 1)
		    vec_y[row] += mat[cell][g*(_ncg+4)+_ncg+2] * 
			vec_x[(cell+1)*_ncg+g];

		if (y != 0)
		    vec_y[row] += mat[cell][g*(_ncg+4)+_ncg+3] * 
			vec_x[(cell-_cx)*_ncg+g];

		for (int e = 0; e < _ncg; e++)
                    vec_y[row] += mat[cell][g*(_ncg+4)+2+e] * 
			vec_x[cell*_ncg+e];				    
	    }
	}
    }
}


void Cmfd::matSubtract(double** AM, double** A, double omega, double** M){

    /* copy A to AM */
    for (int i = 0; i < _cx*_cy; i++){
	for (int g = 0; g < _ncg*(_ncg+4); g++)
	    AM[i][g] = A[i][g];
    }

    for (int i = 0; i < _cx*_cy; i++){
	for (int e = 0; e < _ncg; e++){
	    for (int g = 0; g < _ncg; g++)
		AM[i][g*(_ncg+4)+e+2] -= omega*M[i][g*_ncg+e];
	}
    }
}


double Cmfd::vecMax(double* vec){

    double max = vec[0];

    for (int i = 0; i < _cx*_cy*_ncg; i++)
	max = std::max(max, vec[i]);

    return max;
}


int Cmfd::getNumCmfdGroups(){
  return _ncg;
}


int Cmfd::getCmfdGroupWidth(){
  return _group_width;
}


void Cmfd::setNumCmfdGroups(int num_cmfd_groups){
  _ncg = num_cmfd_groups;
  _mesh->setNumGroups(_ncg);
}


void Cmfd::createGroupStructure(){
  
  _group_width = _ng / _ncg;

  _group_indices = new int[_ncg+1];

  for (int i = 0; i < _ncg; i++){
    _group_indices[i] = i*_group_width;
    log_printf(INFO, "group indices %i: %i", i, _group_indices[i]);
  }

  _group_indices[_ncg] = _ng;
  log_printf(INFO, "group indices %i: %i", _ncg, _group_indices[_ncg]);
  log_printf(INFO, "group width: %i", _group_width);
}

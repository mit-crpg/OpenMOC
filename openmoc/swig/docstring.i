
// File: index.xml

// File: classCell.xml


%feature("docstring") Cell "

Represents a Cell inside of a Universe.  

C++ includes: src/Cell.h
";

%feature("docstring") Cell::incrementVolume "
incrementVolume(double volume)  

Increment the volume/area of the Cell by some amount.  

This routine is called by the TrackGenerator during track generation and segmentation.  

Parameters
----------
* volume :  
    the amount to increment the current volume by  
";

%feature("docstring") Cell::getRotationMatrix "
getRotationMatrix() -> double *  

Return pointer to array for the rotation matrix.  

Returns
-------
a pointer to an array of rotation angles  
";

%feature("docstring") Cell::retrieveTranslation "
retrieveTranslation(double *translations, int num_axes)  

Fills an array with the translations along x, y and z.  

This class method is intended to be called by the OpenMOC Python OpenMC compatiblity
module. Although this method appears to require two arguments, in reality it only requires
one due to SWIG and would be called from within Python as follows:  


Parameters
----------
* translation :  
    an array of translations of length 3 for x, y and z  
* num_axes :  
    the number of axes (this must always be 3)  
";

%feature("docstring") Cell::getUid "
getUid() const  -> int  

Return the Cell's unique ID.  

Returns
-------
the Cell's unique ID  
";

%feature("docstring") Cell::getMaxY "
getMaxY() -> double  

Return the maximum reachable y-coordinate in the Cell.  

Returns
-------
the maximum y-coordinate  
";

%feature("docstring") Cell::getMaxX "
getMaxX() -> double  

Return the maximum reachable x-coordinate in the Cell.  

Returns
-------
the maximum x-coordinate  
";

%feature("docstring") Cell::getName "
getName() const  -> char *  

Return the user-defined name of the Cell.  

Returns
-------
the Cell name  
";

%feature("docstring") Cell::setNumSectors "
setNumSectors(int num_sectors)  

Set the Cell's number of sectors.  

Parameters
----------
* num_sectors :  
    the number of sectors in this Cell  
";

%feature("docstring") Cell::getMinY "
getMinY() -> double  

Return the minimum reachable y-coordinate in the Cell.  

Returns
-------
the minimum y-coordinate  
";

%feature("docstring") Cell::getMinX "
getMinX() -> double  

Return the minimum reachable x-coordinate in the Cell.  

Returns
-------
the minimum x-coordinate  
";

%feature("docstring") Cell::getMinZ "
getMinZ() -> double  

Return the minimum reachable z-coordinate in the Cell.  

Returns
-------
the minimum z-coordinate  
";

%feature("docstring") Cell::getNumSurfaces "
getNumSurfaces() const  -> int  

Return the number of Surfaces in the Cell.  

Returns
-------
the number of Surfaces  
";

%feature("docstring") Cell::setFill "
setFill(Material *fill)  
setFill(Universe *fill)  

Overloaded function
-------------------
* setFill(Material *fill)  
    
    Sets the Material filling this Cell.  

    Parameters:  
    * fill :  
        the Material filling this Cell  

* setFill(Universe *fill)  
    
    Sets the Universe filling this Cell.  

    Parameters:  
    * fill :  
        the Universe filling this Cell  
";

%feature("docstring") Cell::Cell "
Cell(int id=0, const char *name=\"\")  

Constructor sets the unique and user-specifed IDs for this Cell.  

Parameters
----------
* id :  
    the user-specified optional Cell ID  
* name :  
    the user-specified optional Cell name  
";

%feature("docstring") Cell::minSurfaceDist "
minSurfaceDist(LocalCoords *coords) -> double  

Computes the minimum distance to a Surface from a point with a given trajectory at a
certain angle stored in a LocalCoords object.  

If the trajectory will not intersect any of the Surfaces in the Cell returns INFINITY.  

Parameters
----------
* coords :  
    a pointer to a localcoords  
";

%feature("docstring") Cell::clone "
clone() -> Cell *  

Create a duplicate of the Cell.  

Returns
-------
a pointer to the clone  
";

%feature("docstring") Cell::getVolume "
getVolume() -> double  

Return the aggregate volume/area of all instances of this Cell.  

The volume/area of the Cell is computed from track segments which overlap this Cell during
track generation.  

Returns
-------
the volume/area of the Cell  
";

%feature("docstring") Cell::getMaxZ "
getMaxZ() -> double  

Return the maximum reachable z-coordinate in the Cell.  

Returns
-------
the maximum z-coordinate  
";

%feature("docstring") Cell::setNumRings "
setNumRings(int num_rings)  

Set the Cell's number of rings.  

Parameters
----------
* num_rings :  
    the number of rings in this Cell  
";

%feature("docstring") Cell::addSurface "
addSurface(int halfspace, Surface *surface)  

Insert a Surface into this Cell's container of bounding Surfaces.  

Parameters
----------
* halfspace :  
    the Surface halfspace (+/-1)  
* surface :  
    a pointer to the Surface  
";

%feature("docstring") Cell::getAllCells "
getAllCells() -> std::map< int, Cell * >  

Returns the std::map of Cell IDs and Cell pointers within any nested Universes filling
this Cell.  

Returns
-------
std::map of Cell IDs and pointers  
";

%feature("docstring") Cell::getOldestAncestor "
getOldestAncestor() -> Cell *  

Get the oldest ancestor Cell for this Cell.  

This method traverses the linked list of parent Cells to find the one at the root node.
The oldest ancestor Cell is likely the one created by the user at runtime, while
intermediate ancestors were created during radial and angular spatial discretization.  

Returns
-------
this Cell's oldest ancestor Cell  
";

%feature("docstring") Cell::~Cell "
~Cell()  

Destructor clears vector of Surface pointers bounding the Cell.  
";

%feature("docstring") Cell::getAllUniverses "
getAllUniverses() -> std::map< int, Universe * >  

Returns the std::map of all nested Universe IDs and Universe pointers filling this Cell.  

Returns
-------
std::map of Universe IDs and pointers  
";

%feature("docstring") Cell::setTranslation "
setTranslation(double *translation, int num_axes)  

Set the Cell's translation along the x, y and z axes.  

This method is a helper function to allow OpenMOC users to assign the Cell's translations
in Python. A user must initialize a length 3 NumPy array as input to this function. This
function then stores the data values in the NumPy array in the Cell's translation array.
An example of how this function might be called in Python is as follows:  


Parameters
----------
* translation :  
    the array of translations  
* num_axes :  
    the number of axes (this must always be 3)  
";

%feature("docstring") Cell::retrieveRotation "
retrieveRotation(double *rotations, int num_axes, std::string units=\"degrees\")  

Fills an array with the rotation angles for x, y and z.  

This class method is intended to be called by the OpenMOC Python OpenMC compatiblity
module. Although this method appears to require two arguments, in reality it only requires
one due to SWIG and would be called from within Python as follows:  


Parameters
----------
* rotation :  
    an array of rotation angles of length 3 for x, y and z  
* num_axes :  
    the number of axes (this must always be 3)  
* units :  
    the angular units in \"radians\" or \"degrees\" (default)  
";

%feature("docstring") Cell::addNeighborCell "
addNeighborCell(Cell *cell)  

Add a neighboring Cell to this Cell's collection of neighbors.  

Parameters
----------
* cell :  
    a pointer to the neighboring Cell  
";

%feature("docstring") Cell::setVolume "
setVolume(double volume)  

Set the volume/area of the Cell.  

Parameters
----------
* volume :  
    the volume/area of the Cell  
";

%feature("docstring") Cell::getType "
getType() const  -> cellType  

Return the Cell type (FILL or MATERIAL).  

Returns
-------
the Cell type  
";

%feature("docstring") Cell::getMinXBoundaryType "
getMinXBoundaryType() -> boundaryType  

Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at the minimum reachable
x-coordinate in the Cell.  

Returns
-------
the boundary condition at the minimum x-coordinate  
";

%feature("docstring") Cell::getMaxXBoundaryType "
getMaxXBoundaryType() -> boundaryType  

Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at the maximum reachable
x-coordinate in the Cell.  

Returns
-------
the boundary condition at the maximum x-coordinate  
";

%feature("docstring") Cell::getNumSectors "
getNumSectors() -> int  

Return the number of sectors in the Cell.  

Returns
-------
the number of sectors  
";

%feature("docstring") Cell::buildNeighbors "
buildNeighbors()  

Build a collection of neighboring Cells for optimized ray tracing.  
";

%feature("docstring") Cell::isTranslated "
isTranslated() -> bool  

Return a boolean indicating whether the Cell has been translated.  

Returns
-------
whether the Cell has been translated  
";

%feature("docstring") Cell::getTranslation "
getTranslation() -> double *  

Return pointer to array for the translations along x, y and z.  

Returns
-------
a pointer to an array of translations  
";

%feature("docstring") Cell::getTheta "
getTheta(std::string units=\"degrees\") -> double  

Get the rotation angle about the y-axis in degrees.  

Parameters
----------
* units :  
    the angular units in \"radians\" or \"degrees\" (default)  

Returns
-------
the rotation angle about the y-axis  
";

%feature("docstring") Cell::getPhi "
getPhi(std::string units=\"degrees\") -> double  

Get the rotation angle about the x-axis in degrees.  

Parameters
----------
* units :  
    the angular units in \"radians\" or \"degrees\" (default)  

Returns
-------
the rotation angle about the x-axis  
";

%feature("docstring") Cell::subdivideCell "
subdivideCell(double max_radius)  

Subdivides a Cell into rings and sectors aligned with the z-axis.  

This method uses the Cell's clone method to produce a vector of this Cell's subdivided
ring and sector Cells.  

Parameters
----------
* max_radius :  
    the maximum allowable radius used in the subdivisions  

Returns
-------
a vector of Cell pointers to the new subdivided Cells  

A container of all Cell clones created for rings and sectors  
";

%feature("docstring") Cell::getNumRings "
getNumRings() -> int  

Return the number of rings in the Cell.  

Returns
-------
the number of rings  
";

%feature("docstring") Cell::getMaxYBoundaryType "
getMaxYBoundaryType() -> boundaryType  

Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at the maximum reachable
y-coordinate in the Cell.  

Returns
-------
the boundary condition at the maximum y-coordinate  
";

%feature("docstring") Cell::containsPoint "
containsPoint(Point *point) -> bool  

Determines whether a Point is contained inside a Cell.  

Queries each Surface inside the Cell to determine if the Point is on the same side of the
Surface. This point is only inside the Cell if it is on the same side of every Surface in
the Cell.  

Parameters
----------
* point :  
    a pointer to a Point  
";

%feature("docstring") Cell::setName "
setName(const char *name)  

Sets the name of the Cell.  

Parameters
----------
* name :  
    the Cell name string  
";

%feature("docstring") Cell::removeSurface "
removeSurface(Surface *surface)  

Removes a Surface from this Cell's container of bounding Surfaces.  

Parameters
----------
* surface :  
    a pointer to the Surface to remove  
";

%feature("docstring") Cell::getPsi "
getPsi(std::string units=\"degrees\") -> double  

Get the rotation angle about the z-axis in degrees.  

Parameters
----------
* units :  
    the angular units in \"radians\" or \"degrees\" (default)  

Returns
-------
the rotation angle about the z-axis  
";

%feature("docstring") Cell::getNumInstances "
getNumInstances() -> int  

Return the number of instances of this Cell in the Geometry.  

The number of instances of this Cell in the Geometry is determined during track
generation.  

Returns
-------
the number of cell instances  
";

%feature("docstring") Cell::toString "
toString() -> std::string  

Convert this Cell's attributes to a string format.  

Returns
-------
a character array of this Cell's attributes  

Add string data for the Surfaces in this Cell  
";

%feature("docstring") Cell::getNeighbors "
getNeighbors() const  -> std::vector< Cell * >  

Return the std::vector of neighbor Cells to this Cell.  

Returns
-------
std::vector of neighbor Cell pointers  
";

%feature("docstring") Cell::setNumInstances "
setNumInstances(int num_instances)  

Set the number of instances of this Cell.  

Parameters
----------
* num_instances :  
    the number of instances of this Cell in the Geometry  
";

%feature("docstring") Cell::getSurfaces "
getSurfaces() const  -> std::map< int, surface_halfspace * >  

Return the std::map of Surface pointers and halfspaces (+/-1) for all surfaces bounding
the Cell.  

Returns
-------
std::map of Surface pointers and halfspaces  
";

%feature("docstring") Cell::getFillUniverse "
getFillUniverse() -> Universe *  

Return a pointer to the Material filling this Cell.  

Returns
-------
the Material fill pointer  
";

%feature("docstring") Cell::setRotation "
setRotation(double *rotation, int num_axes, std::string units=\"degrees\")  

Set the Cell's rotation angles about the x, y and z axes.  

This method is a helper function to allow OpenMOC users to assign the Cell's rotation
angles in Python. A user must initialize a length 3 NumPy array as input to this function.
This function then stores the data values in the NumPy array in the Cell's rotation array.
An example of how this function might be called in Python is as follows:  


Parameters
----------
* rotation :  
    the array of rotation angles  
* num_axes :  
    the number of axes (this must always be 3)  
* units :  
    the angular units in \"radians\" or \"degrees\" (default)  
";

%feature("docstring") Cell::isRotated "
isRotated() -> bool  

Return a boolean indicating whether the Cell has been rotated.  

Returns
-------
whether the Cell has been rotated  
";

%feature("docstring") Cell::setParent "
setParent(Cell *parent)  

Assign a parent Cell to this Cell.  

This is used by Cell cloning when applied for radial and angular discretization.  

Parameters
----------
* parent :  
    a pointer to the parent Cell  
";

%feature("docstring") Cell::incrementNumInstances "
incrementNumInstances()  

Increment the number of instances of this Cell.  

This routine is called by the TrackGenerator during track generation and segmentation.  
";

%feature("docstring") Cell::printString "
printString()  

Prints a string representation of all of the Cell's attributes to the console.  
";

%feature("docstring") Cell::isFissionable "
isFissionable() -> bool  

Returns true if this Cell is filled with a fissionable Material.  

If the Cell is filled by a Material, this method will simply query the filling Material.
If the Cell is filled by a Universe, this method will consider any Materials filling those
Cells contained by the filling Universe. This method should not be called prior to the
calling of the Geometry::computeFissionability() method.  

Returns
-------
true if contains a fissionable Material  
";

%feature("docstring") Cell::hasParent "
hasParent() -> bool  

Return true if the Cell has a parent and false otherwise.  

Returns
-------
whether the Cell has a parent Cell  
";

%feature("docstring") Cell::getMinYBoundaryType "
getMinYBoundaryType() -> boundaryType  

Return the boundary condition (REFLECTIVE, VACUUM, or INTERFACE) at the minimum reachable
y-coordinate in the Cell.  

Returns
-------
the boundary condition at the minimum y-coordinate  
";

%feature("docstring") Cell::getId "
getId() const  -> int  

Return the Cell's user-specified ID.  

Returns
-------
the Cell's user-specified ID  
";

%feature("docstring") Cell::getParent "
getParent() -> Cell *  

Return this Cell's parent Cell.  

If no parent Cell has been assigned from Cell cloning, then NULL is returned.  

Returns
-------
a pointer to the parent Cell  
";

%feature("docstring") Cell::getFillMaterial "
getFillMaterial() -> Material *  

Return a pointer to the Material filling this Cell.  

Returns
-------
the Material fill pointer  
";

%feature("docstring") Cell::containsCoords "
containsCoords(LocalCoords *coords) -> bool  

Determines whether a Point is contained inside a Cell.  

Queries each Surface inside the Cell to determine if the Point is on the same side of the
Surface. This Point is only inside the Cell if it is on the same side of every Surface in
the Cell.  

Parameters
----------
* coords :  
    a pointer to a localcoord  
";

// File: classCmfd.xml


%feature("docstring") Cmfd "

A class for Coarse Mesh Finite Difference (CMFD) acceleration.  

C++ includes: src/Cmfd.h
";

%feature("docstring") Cmfd::initializeGroupMap "
initializeGroupMap()  

Initialize and set array that links the MOC energy groups to the CMFD energy groups.  

This method initializes the _group_indices_map, which is a 1D array of length
_num_moc_groups that maps the MOC energy groups to CMFD energy groups. The indices into
_group_indices_map are the MOC energy groups and the values are the CMFD energy groups.  
";

%feature("docstring") Cmfd::getNumY "
getNumY() -> int  

Get the number of Mesh cells in a column.  

Returns
-------
The number of Mesh cells in a column  
";

%feature("docstring") Cmfd::getNumX "
getNumX() -> int  

Get the number of Mesh cells in a row.  

Returns
-------
The number of Mesh cells in a row  
";

%feature("docstring") Cmfd::~Cmfd "
~Cmfd()  

Destructor.  
";

%feature("docstring") Cmfd::setNumFSRs "
setNumFSRs(int num_fsrs)  

Set the number of FSRs.  

Parameters
----------
* num_fsrs :  
    The number of FSRs  
";

%feature("docstring") Cmfd::setFSRVolumes "
setFSRVolumes(FP_PRECISION *FSR_volumes)  

Set the pointer to the array of FSR_volumes.  

Parameters
----------
* FSR_volumes :  
    Array of FSR volumes  
";

%feature("docstring") Cmfd::getNumCmfdGroups "
getNumCmfdGroups() -> int  

Get the number of coarse CMFD energy groups.  

Returns
-------
The number of CMFD energy groups  
";

%feature("docstring") Cmfd::getCmfdGroup "
getCmfdGroup(int group) -> int  

Get the CMFD group given an MOC group.  

Parameters
----------
* group :  
    The MOC energy group  

Returns
-------
The CMFD energy group  
";

%feature("docstring") Cmfd::setLatticeStructure "
setLatticeStructure(int num_x, int num_y)  

The structure of the Lattice to be used as the CMFD mesh.  

Parameters
----------
* num_x :  
    The number of cells in the x direction.  
* num_y :  
    The number of cells in the y direction.  
";

%feature("docstring") Cmfd::getNumMOCGroups "
getNumMOCGroups() -> int  

Get the number of MOC energy groups.  

Returns
-------
The number of MOC energy groups  
";

%feature("docstring") Cmfd::Cmfd "
Cmfd()  

Constructor initializes boundaries and variables that describe the Cmfd object.  

The construcor initializes the many variables that describe the CMFD mesh and are used to
solve the nonlinear diffusion acceleration problem.  
";

%feature("docstring") Cmfd::setGeometry "
setGeometry(Geometry *geometry)  

Set a pointer to the Geometry.  

Parameters
----------
* goemetry :  
    A pointer to a Geometry object.  
";

%feature("docstring") Cmfd::zeroCurrents "
zeroCurrents()  

Zero the surface currents for each mesh cell and energy group.  
";

%feature("docstring") Cmfd::setGroupStructure "
setGroupStructure(std::vector< std::vector< int > > group_indices)  

Set a coarse energy group structure for CMFD.  

CMFD does not necessarily need to have the same energy group structure as the MOC problem.
This function can be used to set a sparse energy group structure to speed up the CMFD
solve. An example of how this may be called from Python to use a coarse 2-group CMFD
structure atop a fine 7-group MOC structure is illustrated below:  


Parameters
----------
* group_indices :  
    A nested vector of MOC-to-CMFD group mapping  
";

%feature("docstring") Cmfd::setWidthY "
setWidthY(double width)  

Set Mesh width in the y-direction.  

Parameters
----------
* width :  
    Physical width of Mesh in the y-direction  
";

%feature("docstring") Cmfd::setWidthX "
setWidthX(double width)  

Set Mesh width in the x-direction.  

Parameters
----------
* width :  
    Physical width of Mesh in the x-direction  
";

%feature("docstring") Cmfd::setFluxUpdateOn "
setFluxUpdateOn(bool flux_update_on)  

Set flag indicating whether to update the MOC flux.  

Parameters
----------
* flux_update_on :  
    Boolean saying whether to update MOC flux.  
";

%feature("docstring") Cmfd::getCellFSRs "
getCellFSRs() -> std::vector< std::vector< int > > *  

Return a pointer to the vector of vectors that contains the FSRs that lie in each cell.  

Returns
-------
Vector of vectors containing FSR IDs in each cell.  
";

%feature("docstring") Cmfd::setFSRMaterials "
setFSRMaterials(Material **FSR_materials)  

Set the FSR materials array pointer.  

Parameters
----------
* FSR_materials :  
    Pointer to FSR_materials array  
";

%feature("docstring") Cmfd::initializeLattice "
initializeLattice(Point *offset)  

Initialize the CMFD lattice.  
";

%feature("docstring") Cmfd::convertFSRIdToCmfdCell "
convertFSRIdToCmfdCell(int fsr_id) -> int  

Return the CMFD cell ID that an FSR lies in.  

Note that a CMFD cell is not an actual Cell object; rather, a CMFD cell is just a way of
describing each of the rectangular regions that make up a CMFD lattice. CMFD cells are
numbered with 0 in the lower left corner and monotonically increasing from left to right
and from bottom to top. For example, the indices for a 4 x 4 lattice are: 12 13 14 15 8 9
10 11 4 5 6 7 0 1 2 3  

Parameters
----------
* fsr_id :  
    The FSR ID.  

Returns
-------
The CMFD cell ID. Return -1 if cell is not found.  
";

%feature("docstring") Cmfd::setCellFSRs "
setCellFSRs(std::vector< std::vector< int > > *cell_fsrs)  

Set the vector of vectors that contains the FSRs that lie in each cell.  

Parameters
----------
* cell_fsrs :  
    Vector of vectors containing FSR IDs in each cell.  
";

%feature("docstring") Cmfd::tallyCurrent "
tallyCurrent(segment *curr_segment, FP_PRECISION *track_flux, int azim_index, bool fwd)  

Tallies the current contribution from this segment across the the appropriate CMFD mesh
cell surface.  

Parameters
----------
* curr_segment :  
    The current Track segment  
* track_flux :  
    The outgoing angular flux for this segment  
* azim_index :  
    Azimuthal angle index of the current Track  
* fwd :  
    Boolean indicating direction of integration along segment  
";

%feature("docstring") Cmfd::setCentroidUpdateOn "
setCentroidUpdateOn(bool centroid_update_on)  

Set flag indicating whether to use FSR centroids to update the MOC flux.  

Parameters
----------
* centroid_update_on :  
    Flag saying whether to use centroids to update MOC flux.  
";

%feature("docstring") Cmfd::findCmfdSurface "
findCmfdSurface(int cell_id, LocalCoords *coords) -> int  

Find the cmfd surface that a LocalCoords object lies on.  

If the coords is not on a surface, -1 is returned. Otherwise, the surface ID is returned.  

Parameters
----------
* cell_id :  
    The CMFD cell ID that the local coords is in.  
* coords :  
    The coords being evaluated.  

Returns
-------
The surface ID.  
";

%feature("docstring") Cmfd::getBoundary "
getBoundary(int side) -> int  

Get the boundaryType for one side of the CMFD mesh.  

Parameters
----------
* side :  
    The CMFD mesh surface ID.  

Returns
-------
The boundaryType for the surface.  
";

%feature("docstring") Cmfd::setNumMOCGroups "
setNumMOCGroups(int num_moc_groups)  

Set the number of MOC energy groups.  

Parameters
----------
* num_groups :  
    Number of MOC energy groups  
";

%feature("docstring") Cmfd::setSORRelaxationFactor "
setSORRelaxationFactor(FP_PRECISION SOR_factor)  

Set the successive over-relaxation factor for the linear solve within the diffusion
eigenvalue solve.  

Parameters
----------
* SOR_factor :  
    Over-relaxation factor  
";

%feature("docstring") Cmfd::initialize "
initialize()  

Initialize the Matrix and Vector objects, k-nearest stencils, the CMFD cell currents and
MOC materials.  
";

%feature("docstring") Cmfd::setKNearest "
setKNearest(int k_nearest)  

Set a number of k-nearest neighbor cells to use in updating the FSR flux.  

Parameters
----------
* k_nearest :  
    The number of nearest neighbor CMFD cells.  
";

%feature("docstring") Cmfd::setFSRFluxes "
setFSRFluxes(FP_PRECISION *scalar_flux)  

Set pointer to FSR flux array.  

Parameters
----------
* scalar_flux :  
    Pointer to FSR flux array  
";

%feature("docstring") Cmfd::setNumY "
setNumY(int num_y)  

Set the number of Mesh cells in a column.  

Parameters
----------
* num_y :  
    Number of Mesh cells in a column  
";

%feature("docstring") Cmfd::setNumX "
setNumX(int num_x)  

Set the number of Mesh cells in a row.  

Parameters
----------
* num_x :  
    Number of Mesh cells in a row  
";

%feature("docstring") Cmfd::getLattice "
getLattice() -> Lattice *  

Returns the Lattice object used as the CMFD mesh.  

Returns
-------
A pointer to a Lattice object.  
";

%feature("docstring") Cmfd::addFSRToCell "
addFSRToCell(int cell_id, int fsr_id)  

Add an FSR ID to a vector that contains all the FSR IDs contained within a CMFD mesh cell.  

Parameters
----------
* cell_id :  
    The CMFD cell ID.  
* fsr_id :  
    The FSR ID.  
";

%feature("docstring") Cmfd::findCmfdCell "
findCmfdCell(LocalCoords *coords) -> int  

Find the CMFD cell that a LocalCoords object is in.  

Parameters
----------
* coords :  
    The coords being evaluated.  

Returns
-------
The CMFD cell ID.  
";

%feature("docstring") Cmfd::updateBoundaryFlux "
updateBoundaryFlux(Track **tracks, FP_PRECISION *boundary_flux, int num_tracks)  

Update the MOC boundary fluxes.  

The MOC boundary fluxes are updated using the P0 approximation. With this approximation,
the boundary fluxes are updated using the ratio of new to old flux for the cell that the
outgoing flux from the track enters.  

Parameters
----------
* tracks :  
    2D array of Tracks  
* boundary_flux :  
    Array of boundary fluxes  

Returns
-------
The number of Tracks  
";

%feature("docstring") Cmfd::setBoundary "
setBoundary(int side, boundaryType boundary)  

Set the CMFD boundary type for a given surface.  

The CMFD boundary is assumed to be rectangular with the surfaces identified by constants
in the constants.h file.  

Parameters
----------
* side :  
    The CMFD surface UID.  
* boundary :  
    The boundaryType of the surface.  
";

%feature("docstring") Cmfd::setQuadrature "
setQuadrature(Quadrature *quadrature)  

Sets the Quadrature object in use by the MOC Solver.  

Parameters
----------
* quadrature :  
    A Quadrature object pointer from the Solver  
";

%feature("docstring") Cmfd::setSourceConvergenceThreshold "
setSourceConvergenceThreshold(FP_PRECISION source_thresh)  

Sets the threshold for CMFD source convergence (>0)  

Parameters
----------
* the :  
    threshold for source convergence  
";

%feature("docstring") Cmfd::isFluxUpdateOn "
isFluxUpdateOn() -> bool  

Get flag indicating whether to update the MOC flux.  

Returns
-------
Boolean saying whether to update MOC flux.  
";

%feature("docstring") Cmfd::computeKeff "
computeKeff(int moc_iteration) -> FP_PRECISION  

Solve the nonlinear diffusion acceleration problem to accelerate the convergence of the
MOC problem.  

This method uses the information from the last MOC transport sweep and solves a simplified
nonlinear diffusion problem. The diffusion problem is tightly converged and the solution
is used to update the the solution of the MOC problem.  

Parameters
----------
* moc_iteration :  
    MOC iteration number  

Returns
-------
The dominant eigenvalue of the nonlinear diffusion problem  
";

%feature("docstring") Cmfd::initializeCellMap "
initializeCellMap()  

Initializes the vector of vectors that links CMFD cells with FSRs.  

This method is called by the geometry once the CMFD mesh has been initialized by the
geometry. This method allocates a vector for each CMFD cell that is used to store the FSR
ids contained within that cell.  
";

%feature("docstring") Cmfd::isCentroidUpdateOn "
isCentroidUpdateOn() -> bool  

Get flag indicating whether to use FSR centroids to update the MOC flux.  

Returns
-------
Flag saying whether to use centroids to update MOC flux.  
";

%feature("docstring") Cmfd::getNumCells "
getNumCells() -> int  

Get the number of CMFD cells.  

Returns
-------
The number of CMFD cells  
";

// File: classCPUSolver.xml


%feature("docstring") CPUSolver "

This a subclass of the Solver class for multi-core CPUs using OpenMP multi-threading.  

C++ includes: src/CPUSolver.h
";

%feature("docstring") CPUSolver::initializeFSRs "
initializeFSRs()  

Initializes the FSR volumes and Materials array.  

This method gets an array of OpenMP mutual exclusion locks for each FSR for use in the
transport sweep algorithm.  
";

%feature("docstring") CPUSolver::computeResidual "
computeResidual(residualType res_type) -> double  

Computes the residual between source/flux iterations.  

Parameters
----------
* res_type :  
    the type of residuals to compute (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)  

Returns
-------
the average residual in each FSR  
";

%feature("docstring") CPUSolver::initializeSourceArrays "
initializeSourceArrays()  

Allocates memory for FSR source arrays.  

Deletes memory for old source arrays if they were allocated for a previous simulation.  
";

%feature("docstring") CPUSolver::computeFSRFissionRates "
computeFSRFissionRates(double *fission_rates, int num_FSRs)  

Computes the volume-integrated, energy-integrated nu-fission rate in each FSR and stores
them in an array indexed by FSR ID.  

This is a helper method for SWIG to allow users to retrieve FSR nu-fission rates as a
NumPy array. An example of how this method can be called from Python is as follows:  


Parameters
----------
* fission_rates :  
    an array to store the nu-fission rates (implicitly passed in as a NumPy array from
    Python)  
* num_FSRs :  
    the number of FSRs passed in from Python  
";

%feature("docstring") CPUSolver::computeFSRSources "
computeFSRSources()  

Computes the total source (fission, scattering, fixed) in each FSR.  

This method computes the total source in each FSR based on this iteration's current
approximation to the scalar flux.  
";

%feature("docstring") CPUSolver::storeFSRFluxes "
storeFSRFluxes()  

Stores the FSR scalar fluxes in the old scalar flux array.  
";

%feature("docstring") CPUSolver::computeFSRScatterSources "
computeFSRScatterSources()  

Computes the total scattering source in each FSR.  

This method is a helper routine for the openmoc.krylov submodule.  
";

%feature("docstring") CPUSolver::initializeFluxArrays "
initializeFluxArrays()  

Allocates memory for Track boundary angular and FSR scalar fluxes.  

Deletes memory for old flux arrays if they were allocated for a previous simulation.  
";

%feature("docstring") CPUSolver::flattenFSRFluxes "
flattenFSRFluxes(FP_PRECISION value)  

Set the scalar flux for each FSR and energy group to some value.  

Parameters
----------
* value :  
    the value to assign to each FSR scalar flux  
";

%feature("docstring") CPUSolver::initializeFixedSources "
initializeFixedSources()  

Populates array of fixed sources assigned by FSR.  
";

%feature("docstring") CPUSolver::addSourceToScalarFlux "
addSourceToScalarFlux()  

Add the source term contribution in the transport equation to the FSR scalar flux.  
";

%feature("docstring") CPUSolver::transportSweep "
transportSweep()  

This method performs one transport sweep of all azimuthal angles, Tracks, Track segments,
polar angles and energy groups.  

The method integrates the flux along each Track and updates the boundary fluxes for the
corresponding output Track, while updating the scalar flux in each flat source region.  
";

%feature("docstring") CPUSolver::setFluxes "
setFluxes(FP_PRECISION *in_fluxes, int num_fluxes)  

Set the flux array for use in transport sweep source calculations.  This is a helper
method for the checkpoint restart capabilities, as well as the IRAMSolver in the
openmoc.krylov submodule. This routine may be used as follows from within Python:  

  

         NOTE: This routine stores a pointer to the fluxes for the Solver
         to use during transport sweeps and other calculations. Hence, the
         flux array pointer is shared between NumPy and the Solver.  

Parameters
----------
* in_fluxes :  
    an array with the fluxes to use  
* num_fluxes :  
    the number of flux values (# groups x # FSRs)  
";

%feature("docstring") CPUSolver::computeFSRFissionSources "
computeFSRFissionSources()  

Computes the total fission source in each FSR.  

This method is a helper routine for the openmoc.krylov submodule.  
";

%feature("docstring") CPUSolver::getFluxes "
getFluxes(FP_PRECISION *out_fluxes, int num_fluxes)  

Fills an array with the scalar fluxes.  

This class method is a helper routine called by the OpenMOC Python \"openmoc.krylov\"
module for Krylov subspace methods. Although this method appears to require two arguments,
in reality it only requires one due to SWIG and would be called from within Python as
follows:  


Parameters
----------
* fluxes :  
    an array of FSR scalar fluxes in each energy group  
* num_fluxes :  
    the total number of FSR flux values  
";

%feature("docstring") CPUSolver::zeroTrackFluxes "
zeroTrackFluxes()  

Zero each Track's boundary fluxes for each energy group and polar angle in the \"forward\"
and \"reverse\" directions.  
";

%feature("docstring") CPUSolver::CPUSolver "
CPUSolver(TrackGenerator *track_generator=NULL)  

Constructor initializes array pointers for Tracks and Materials.  

The constructor retrieves the number of energy groups and FSRs and azimuthal angles from
the Geometry and TrackGenerator if passed in as parameters by the user. The constructor
initalizes the number of OpenMP threads to a default of 1.  

Parameters
----------
* track_generator :  
    an optional pointer to the TrackGenerator  
";

%feature("docstring") CPUSolver::normalizeFluxes "
normalizeFluxes()  

Normalizes all FSR scalar fluxes and Track boundary angular fluxes to the total fission
source (times $ \\nu $).  
";

%feature("docstring") CPUSolver::computeKeff "
computeKeff()  

Compute $ k_{eff} $ from successive fission sources.  
";

%feature("docstring") CPUSolver::getNumThreads "
getNumThreads() -> int  

Returns the number of shared memory OpenMP threads in use.  

Returns
-------
the number of threads  
";

%feature("docstring") CPUSolver::setNumThreads "
setNumThreads(int num_threads)  

Sets the number of shared memory OpenMP threads to use (>0).  

Parameters
----------
* num_threads :  
    the number of threads  
";

// File: structdev__material.xml


%feature("docstring") dev_material "

A Material's nuclear data to be stored on a GPU.  

Attributes
----------
* _id : int  
    A user-defined ID for each Material created  

* _sigma_t : FP_PRECISION *  
    An array of the total cross-sections for each energy group  

* _sigma_f : FP_PRECISION *  
    A 2D array of the scattering cross-section matrix. The first index is row number and
    second index is column number  

* _nu_sigma_f : FP_PRECISION *  
    An array of the fission cross-sections multiplied by nu $ \\nu $ for each energy group  

* _chi : FP_PRECISION *  
    An array of the chi $ \\chi $ values for each energy group  

* _fiss_matrix : FP_PRECISION *  
    A 2D array of the fission matrix from/into each group  

* _sigma_s : FP_PRECISION *  
    A 2D array of the scattering cross-section matrix from/into each group  

C++ includes: DeviceMaterial.h
";

%feature("docstring") dev_material::dev_material "
dev_material()  

Constructor for a dev_material struct on a GPU.  
";

%feature("docstring") dev_material::~dev_material "
~dev_material()  

Destructor releases data for all Material's cross-sections on GPU.  
";

// File: structdev__segment.xml


%feature("docstring") dev_segment "

A dev_segment represents a line segment within a single flat source region along a track.  

The dev_segment is intended for use on the GPU.  

Attributes
----------
* _length : FP_PRECISION  
    The length of the segment (cm)  

* _material_index : int  
    An index into the _materials array that contains Material pointers  

* _region_uid : int  
    The ID for flat source region in which this segment resides  

C++ includes: DeviceTrack.h
";

// File: structdev__track.xml


%feature("docstring") dev_track "

A dev_track represents a characteristic line across the geometry.  

A dev_track has particular starting and ending points on the boundaries of the geometry
and an azimuthal angle. The dev_track is intended for use on the GPU.  

Attributes
----------
* _uid : int  
    A monotonically increasing unique ID for each Track created  

* _azim_angle_index : int  
    The azimuthal angle index into the global 2D ragged array of Tracks  

* _segments : dev_segment *  
    A vector of segments making up this track  

* _num_segments : int  
    The number of segments making up this Track  

* _track_in : int  
    Index of the next Track when traveling along this Track in the \"forward\" direction.  

* _track_out : int  
    Index of the next Track when traveling along this Track in the \"reverse\" direction.  

* _next_in : bool  
    A boolean to indicate whether to give the flux to the \"forward\" (false) or
    \"reverse\" (true) direction of the next Track going in the \"forward\" direction.  

* _next_out : bool  
    A boolean to indicate whether to give the flux to the \"forward\" (false) or
    \"reverse\" (true) direction of the next Track going in the \"reverse\" direction.  

* _transfer_flux_in : bool  
    A boolean to indicate whether the outgoing angular flux along this Track's \"forward\"
    direction should be transferred to the outgoing Track.  

* _transfer_flux_out : bool  
    A boolean to indicate whether the outgoing angular flux along this Track's \"reverse\"
    direction should be transferred to the incoming Track.  

C++ includes: DeviceTrack.h
";

// File: classEqualAnglePolarQuad.xml


%feature("docstring") EqualAnglePolarQuad "

Equal angle polar quadrature.  

C++ includes: src/Quadrature.h
";

%feature("docstring") EqualAnglePolarQuad::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles  
";

%feature("docstring") EqualAnglePolarQuad::initialize "
initialize()  

Routine to initialize the polar quadrature.  

This routine generates the sine thetas and weights.  
";

%feature("docstring") EqualAnglePolarQuad::EqualAnglePolarQuad "
EqualAnglePolarQuad()  

Dummy constructor calls the parent constructor.  
";

%feature("docstring") EqualAnglePolarQuad::precomputeWeights "
precomputeWeights(bool solve_3D)  

Calculates total weights for every azimuthal/polar combination based on the equal angle
polar quadrature.  

Parameters
----------
* solve_3D :  
    Boolean indicating whether this is a 3D quadrature  
";

// File: classEqualWeightPolarQuad.xml


%feature("docstring") EqualWeightPolarQuad "

Equal weight polar quadrature.  

C++ includes: src/Quadrature.h
";

%feature("docstring") EqualWeightPolarQuad::EqualWeightPolarQuad "
EqualWeightPolarQuad()  

Dummy constructor calls the parent constructor.  
";

%feature("docstring") EqualWeightPolarQuad::initialize "
initialize()  

Routine to initialize the polar quadrature.  

This routine generates the sine thetas and weights.  
";

%feature("docstring") EqualWeightPolarQuad::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles  
";

%feature("docstring") EqualWeightPolarQuad::precomputeWeights "
precomputeWeights(bool solve_3D)  

Calculates total weights for every azimuthal/polar combination based on the equal weight
polar quadrature.  

Parameters
----------
* solve_3D :  
    Boolean indicating whether this is a 3D quadrature  
";

// File: classExpEvaluator.xml


%feature("docstring") ExpEvaluator "

This is a class for evaluating exponentials.  

The ExpEvaluator includes different algorithms to evaluate exponentials with varying
degrees of accuracy and speed. This is a helper class for the Solver and its subclasses
and it not intended to be initialized as a standalone object.  

C++ includes: src/ExpEvaluator.h
";

%feature("docstring") ExpEvaluator::getTableSize "
getTableSize() -> int  

Get the number of entries in the exponential interpolation table.  

Parameters
----------
* entries :  
    in the interpolation table  
";

%feature("docstring") ExpEvaluator::computeExponential "
computeExponential(FP_PRECISION tau, int polar) -> FP_PRECISION  

Computes the exponential term for a optical length and polar angle.  

This method computes $ 1 - exp(-\\tau/sin(\\theta_p)) $ for some optical path length and
polar angle. This method uses either a linear interpolation table (default) or the
exponential intrinsic exp(...) function.  

Parameters
----------
* tau :  
    the optical path length (e.g., sigma_t times length)  
* polar :  
    the polar angle index  

Returns
-------
the evaluated exponential  
";

%feature("docstring") ExpEvaluator::getTableSpacing "
getTableSpacing() -> FP_PRECISION  

Returns the exponential table spacing.  

Returns
-------
exponential table spacing  
";

%feature("docstring") ExpEvaluator::setExpPrecision "
setExpPrecision(FP_PRECISION exp_precision)  

Sets the maximum acceptable approximation error for exponentials.  

This routine only affects the construction of the linear interpolation table for
exponentials, if in use. By default, a value of 1E-5 is used for the table, as recommended
by the analysis of Yamamoto in his 2004 paper on the subject.  

Parameters
----------
* exp_precision :  
    the maximum exponential approximation error  
";

%feature("docstring") ExpEvaluator::setMaxOpticalLength "
setMaxOpticalLength(FP_PRECISION max_optical_length)  

Sets the maximum optical length covered in the exponential interpolation table.  

Parameters
----------
* max_optical_length :  
    the maximum optical length  
";

%feature("docstring") ExpEvaluator::useIntrinsic "
useIntrinsic()  

Use the exponential intrinsic exp(...) to compute exponentials.  
";

%feature("docstring") ExpEvaluator::isUsingInterpolation "
isUsingInterpolation() -> bool  

Returns true if using linear interpolation to compute exponentials.  

Returns
-------
true if so, false otherwise  
";

%feature("docstring") ExpEvaluator::getExpTable "
getExpTable() -> FP_PRECISION *  

Returns a pointer to the exponential interpolation table.  

Returns
-------
pointer to the exponential interpolation table  
";

%feature("docstring") ExpEvaluator::useInterpolation "
useInterpolation()  

Use linear interpolation to compute exponentials.  
";

%feature("docstring") ExpEvaluator::getExpPrecision "
getExpPrecision() -> FP_PRECISION  

Gets the maximum acceptable approximation error for exponentials.  

Returns
-------
the maximum exponential approximation error  
";

%feature("docstring") ExpEvaluator::setQuadrature "
setQuadrature(Quadrature *quadrature)  

Set the Quadrature to use when computing exponentials.  

Parameters
----------
* quadrature :  
    a Quadrature object pointer  
";

%feature("docstring") ExpEvaluator::getMaxOpticalLength "
getMaxOpticalLength() -> FP_PRECISION  

Gets the maximum optical length covered with the exponential interpolation table.  

Returns
-------
max_optical_length the maximum optical length  
";

%feature("docstring") ExpEvaluator::ExpEvaluator "
ExpEvaluator()  

Constructor initializes array pointers to NULL.  

The constructor sets the interpolation scheme as the default for computing exponentials.  
";

%feature("docstring") ExpEvaluator::initialize "
initialize()  

If using linear interpolation, builds the table for each polar angle.  

Parameters
----------
* tolerance :  
    the minimum acceptable interpolation accuracy  
";

%feature("docstring") ExpEvaluator::~ExpEvaluator "
~ExpEvaluator()  

Destructor deletes table for linear interpolation of exponentials.  
";

// File: classFixedHashMap.xml


%feature("docstring") FixedHashMap "

A fixed-size hash map supporting insertion and lookup operations.  

The FixedHashMap class supports insertion and lookup operations but not deletion as
deletion is not needed in the OpenMOC application. This hash table uses chaining for
collisions and does not incorporate concurrency objects except for tracking the number of
entries in the table for which an atomic increment is used. This hash table is not thread
safe but is used as a building block for the ParallelHashMap class. This table guarantees
O(1) insertions and lookups on average.  

C++ includes: src/ParallelHashMap.h
";

%feature("docstring") FixedHashMap::clear "
clear()  

Clears all key/value pairs form the hash table.  
";

%feature("docstring") FixedHashMap::keys "
keys() -> K *  

Returns an array of the keys in the fixed-size table.  

All buckets are scanned in order to form a list of all keys present in the table and then
the list is returned. WARNING: The user is responsible for freeing the allocated memory
once the array is no longer needed.  

Returns
-------
an array of keys in the map whose length is the number of key/value pairs in the table.  
";

%feature("docstring") FixedHashMap::bucket_count "
bucket_count() -> size_t  

Returns the number of buckets in the fixed-size table.  

Returns
-------
number of buckets in the map  
";

%feature("docstring") FixedHashMap::values "
values() -> V *  

Returns an array of the values in the fixed-size table.  

All buckets are scanned in order to form a list of all values present in the table and
then the list is returned. WARNING: The user is responsible for freeing the allocated
memory once the array is no longer needed.  

Returns
-------
an array of values in the map whose length is the number of key/value pairs in the table.  
";

%feature("docstring") FixedHashMap::at "
at(K key) -> V &  

Determine the value associated with a given key in the fixed-size table.  

The linked list in the bucket associated with the key is searched and once the key is
found, the corresponding value is returned. An exception is thrown if the key is not
present in the map.  

Parameters
----------
* key :  
    key whose corresponding value is desired  

Returns
-------
value associated with the given key  
";

%feature("docstring") FixedHashMap::size "
size() -> size_t  

Returns the number of key/value pairs in the fixed-size table.  

Returns
-------
number of key/value pairs in the map  
";

%feature("docstring") FixedHashMap::insert_and_get_count "
insert_and_get_count(K key, V value) -> int  

Inserts a key/value pair into the fixed-size table and returns the order number with which
it was inserted.  

The specified key value pair is inserted into the fixed-size table. If the key already
exists in the table, the pair is not inserted and the function returns -1.  

Parameters
----------
* key :  
    key of the key/value pair to be inserted  
* value :  
    value of the key/value pair to be inserted  

Returns
-------
order number in which key/value pair was inserted, -1 is returned if key was already
present in map.  
";

%feature("docstring") FixedHashMap::insert "
insert(K key, V value)  

Inserts a key/value pair into the fixed-size table.  

The specified key value pair is inserted into the fixed-size table. If the key already
exists in the table, the pair is not inserted and the function returns.  

Parameters
----------
* key :  
    key of the key/value pair to be inserted  
* value :  
    value of the key/value pair to be inserted  
";

%feature("docstring") FixedHashMap::~FixedHashMap "
~FixedHashMap()  

Destructor deletes all nodes in the linked lists associated with each bucket in the fixed-
size table and their pointers.  
";

%feature("docstring") FixedHashMap::print_buckets "
print_buckets()  

Prints the contents of each bucket to the screen.  

All buckets are scanned and the contents of the buckets are printed, which are pointers to
linked lists. If the pointer is NULL suggesting that the linked list is empty, NULL is
printed to the screen.  
";

%feature("docstring") FixedHashMap::contains "
contains(K key) -> bool  

Determine whether the fixed-size table contains a given key.  

The linked list in the bucket associated with the key is searched to determine whether the
key is present.  

Parameters
----------
* key :  
    key to be searched  

Returns
-------
boolean value referring to whether the key is contained in the map  
";

%feature("docstring") FixedHashMap::FixedHashMap "
FixedHashMap(size_t M=64)  

Constructor initializes fixed-size table of buckets filled with empty linked lists.  

The constructor initializes a fixed-size hash map with the size as an input parameter. If
no size is given the default size (64) is used. Buckets are filled with empty linked lists
presented as NULL pointers.  

Parameters
----------
* M :  
    size of fixed hash map  
";

// File: structfsr__data.xml


%feature("docstring") fsr_data "

A fsr_data struct represents an FSR with a unique FSR ID and a characteristic point that
lies within the FSR that can be used to recompute the hierarchical LocalCoords linked
list.  

Attributes
----------
* _fsr_id : int  
    The FSR ID  

* _cmfd_cell : int  
    The CMFD Cell  

* _mat_id : int  
    The Material ID  

* _point : Point *  
    Characteristic point in Root Universe that lies in FSR  

* _centroid : Point *  
    Global numerical centroid in Root Universe  

C++ includes: Geometry.h
";

%feature("docstring") fsr_data::~fsr_data "
~fsr_data()  

Destructor for fsr_data  
";

%feature("docstring") fsr_data::fsr_data "
fsr_data()  

Constructor for FSR data initializes centroids and points to NULL  
";

// File: classGeometry.xml


%feature("docstring") Geometry "

The master class containing references to all geometry-related objects - Surfaces, Cells,
Universes and Lattices - and Materials.  

The primary purpose for the geometry is to serve as a collection of all geometry-related
objects, as well as for ray tracing of characteristic tracks across the Geometry and
computing FSR-to-cell offset maps.  

C++ includes: src/Geometry.h
";

%feature("docstring") Geometry::getAllUniverses "
getAllUniverses() -> std::map< int, Universe * >  

Return a std::map container of Universe IDs (keys) with Unierses pointers (values).  

Returns
-------
a std::map of Universes indexed by Universe ID in the geometry  
";

%feature("docstring") Geometry::initializeFSRVectors "
initializeFSRVectors()  

Initialize key and material ID vectors for lookup by FSR ID  This function initializes and
sets reverse lookup vectors by FSR ID. This is called after the FSRs have all been
identified and allocated during segmentation. This function must be called after
Geometry::segmentize() has completed. It should not be called if tracks are loaded from a
file.  
";

%feature("docstring") Geometry::findCellContainingFSR "
findCellContainingFSR(int fsr_id) -> Cell *  

Finds the Cell containing a given fsr ID.  

Parameters
----------
* fsr_id :  
    an FSR ID.  
";

%feature("docstring") Geometry::getAllMaterialCells "
getAllMaterialCells() -> std::map< int, Cell * >  

Return a std::map container of Cell IDs (keys) with Cells pointers (values).  

Returns
-------
a std::map of Cells indexed by Cell ID in the geometry  
";

%feature("docstring") Geometry::getMaxZ "
getMaxZ() -> double  

Return the maximum z-coordinate contained by the Geometry.  

Returns
-------
the maximum z-coordinate (cm)  
";

%feature("docstring") Geometry::getMaxY "
getMaxY() -> double  

Return the maximum y-coordinate contained by the Geometry.  

Returns
-------
the maximum y-coordinate (cm)  
";

%feature("docstring") Geometry::getMaxX "
getMaxX() -> double  

Return the maximum x-coordinate contained by the Geometry.  

Returns
-------
the maximum x-coordinate (cm)  
";

%feature("docstring") Geometry::initializeFSRs "
initializeFSRs(bool neighbor_cells=false)  

Compute the number of flat source regions in the Geometry and initialize CMFD.  

This method is intended to be called by the user before initiating source iteration. This
method first subdivides all Cells by calling the Geometry::subdivideCells() method. Then
it initializes the CMFD object. neighbor_cells whether to use neighbor cell optimizations  
";

%feature("docstring") Geometry::segmentize "
segmentize(Track *track)  

This method performs ray tracing to create Track segments within each flat source region
in the Geometry.  

This method starts at the beginning of a Track and finds successive intersection points
with FSRs as the Track crosses through the Geometry and creates segment structs and adds
them to the Track.  

Parameters
----------
* track :  
    a pointer to a track to segmentize  
";

%feature("docstring") Geometry::getFSRsToKeys "
getFSRsToKeys() -> std::vector< std::string > &  

Returns the vector that maps FSR IDs to FSR key hashes.  

Returns
-------
_FSR_keys_map map of FSR keys to FSR IDs  
";

%feature("docstring") Geometry::setCmfd "
setCmfd(Cmfd *cmfd)  

Sets the pointer to a CMFD object used for acceleration.  

Parameters
----------
* cmfd :  
    a pointer to the CMFD object  
";

%feature("docstring") Geometry::getFSRId "
getFSRId(LocalCoords *coords) -> int  

Return the ID of the flat source region that a given LocalCoords object resides within.  

Parameters
----------
* coords :  
    a LocalCoords object pointer  

Returns
-------
the FSR ID for a given LocalCoords object  
";

%feature("docstring") Geometry::getFSRCentroid "
getFSRCentroid(int fsr_id) -> Point *  

Return the centroid for a given FSR ID.  

Parameters
----------
* fsr_id :  
    the FSR ID  

Returns
-------
the FSR's centroid  
";

%feature("docstring") Geometry::getFSRKeysMap "
getFSRKeysMap() -> ParallelHashMap< std::string, fsr_data * > &  

Returns a pointer to the map that maps FSR keys to FSR IDs.  

Returns
-------
pointer to _FSR_keys_map map of FSR keys to FSR IDs  
";

%feature("docstring") Geometry::getMinY "
getMinY() -> double  

Return the minimum y-coordinate contained by the Geometry.  

Returns
-------
the minimum y-coordinate (cm)  
";

%feature("docstring") Geometry::getMinX "
getMinX() -> double  

Return the minimum x-coordinate contained by the Geometry.  

Returns
-------
the minimum x-coordinate (cm)  
";

%feature("docstring") Geometry::getMinZ "
getMinZ() -> double  

Return the minimum z-coordinate contained by the Geometry.  

Returns
-------
the minimum z-coordinate (cm)  
";

%feature("docstring") Geometry::getMaxYBoundaryType "
getMaxYBoundaryType() -> boundaryType  

Returns the boundary conditions (REFLECTIVE or VACUUM) at the maximum y-coordinate in the
Geometry.  

Returns
-------
the boundary conditions for the maximum y-coordinate in the Geometry  
";

%feature("docstring") Geometry::getAllCells "
getAllCells() -> std::map< int, Cell * >  

Return a std::map container of Cell IDs (keys) with Cells pointers (values).  

Returns
-------
a std::map of Cells indexed by Cell ID in the geometry  
";

%feature("docstring") Geometry::getNumCells "
getNumCells() -> int  

Returns the number of Cells in the Geometry.  

Returns
-------
the number of Cells  
";

%feature("docstring") Geometry::getFSRKey "
getFSRKey(LocalCoords *coords) -> std::string  

Generate a string FSR \"key\" that identifies an FSR by its unique hierarchical
lattice/universe/cell structure.  

Since not all FSRs will reside on the absolute lowest universe level and Cells might
overlap other cells, it is important to have a method for uniquely identifying FSRs. This
method creates a unique FSR key by constructing a structured string that describes the
hierarchy of lattices/universes/cells.  

Parameters
----------
* coords :  
    a LocalCoords object pointer  

Returns
-------
the FSR key  
";

%feature("docstring") Geometry::subdivideCells "
subdivideCells()  

Subdivides all Cells in the Geometry into rings and angular sectors aligned with the
z-axis.  

This method is called by the Geometry::initializeFSRs() method but may also be called by
the user in Python if needed:  

  
";

%feature("docstring") Geometry::getRootUniverse "
getRootUniverse() -> Universe *  

Returns the Universe at the root node in the CSG tree.  

Returns
-------
the root Universe  
";

%feature("docstring") Geometry::getSpatialDataOnGrid "
getSpatialDataOnGrid(std::vector< double > grid_x, std::vector< double > grid_y, double
    zcoord, const char *domain_type=\"material\") -> std::vector< int >  

Get the material, cell or FSR IDs on a 2D spatial grid.  

This is a helper method for the openmoc.plotter module. This method may also be called by
the user in Python if needed. A user must initialize NumPy arrays with the x and y grid
coordinates input to this function. This function then fills a NumPy array with the domain
IDs for each coordinate. An example of how this function might be called in Python is as
follows:  


Parameters
----------
* grid_x :  
    a NumPy array or list of the x-coordinates  
* num_x :  
    the number of x-coordinates in the grid  
* grid_y :  
    a NumPy array or list of the y-coordinates  
* num_y :  
    the number of y-coordinates in the grid  
* zcoord :  
    the z-coordinate to use to find the domain IDs  
* domain_type :  
    the type of domain ('fsr', 'material', 'cell')  

Returns
-------
a NumPy array or list of the domain IDs  
";

%feature("docstring") Geometry::setRootUniverse "
setRootUniverse(Universe *root_universe)  

Sets the root Universe for the CSG tree.  

Parameters
----------
* root_universe :  
    the root Universe of the CSG tree.  
";

%feature("docstring") Geometry::computeFissionability "
computeFissionability(Universe *univ=NULL)  

Determines the fissionability of each Universe within this Geometry.  

A Universe is determined fissionable if it contains a Cell filled by a Material with a
non-zero fission cross-section. Note that this method recurses through all Universes at
each level in the nested Universe hierarchy. Users should only call this method without a
parameter (the default) from Python as follows to ensure that the recursion starts from
the uppermost Universe level:  


Parameters
----------
* univ :  
    the Universe of interest (default is NULL)  
";

%feature("docstring") Geometry::getNumEnergyGroups "
getNumEnergyGroups() -> int  

Returns the number of energy groups for each Material's nuclear data.  

Returns
-------
the number of energy groups  
";

%feature("docstring") Geometry::Geometry "
Geometry()  

Constructor initializes an empty Geometry.  
";

%feature("docstring") Geometry::findCellContainingCoords "
findCellContainingCoords(LocalCoords *coords) -> Cell *  

Find the Cell that this LocalCoords object is in at the lowest level of the nested
Universe hierarchy.  

This method assumes that the LocalCoords has been initialized with coordinates and a
Universe ID. The method will recursively find the Cell on the lowest level of the nested
Universe hierarchy by building a linked list of LocalCoords from the LocalCoord passed in
as an argument down to the lowest level Cell found. In the process it will set the
coordinates at each level of the hierarchy for each LocalCoord in the linked list for the
Lattice or Universe that it is in. If the LocalCoords is outside the bounds of the
Geometry or on the boundaries this method will return NULL; otherwise it will return a
pointer to the Cell that is found by the recursive Geometry::findCell(...) method.  

Parameters
----------
* coords :  
    pointer to a LocalCoords object  

Returns
-------
returns a pointer to a Cell if found, NULL if no Cell found  
";

%feature("docstring") Geometry::getWidthZ "
getWidthZ() -> double  

Returns the total width in the z-direction of the Geometry in cm.  

Returns
-------
the total width of the Geometry in the z-direction (cm)  
";

%feature("docstring") Geometry::getWidthY "
getWidthY() -> double  

Returns the total width in the y-direction of the Geometry in cm.  

Returns
-------
the total width of the Geometry in the y-direction (cm)  
";

%feature("docstring") Geometry::getWidthX "
getWidthX() -> double  

Returns the total width in the x-direction of the Geometry in cm.  

Returns
-------
the total width of the Geometry in the x-direction (cm)  
";

%feature("docstring") Geometry::findFSRMaterial "
findFSRMaterial(int fsr_id) -> Material *  

Find the Material for a flat source region ID.  

Parameters
----------
* fsr_id :  
    a FSR id  

Returns
-------
a pointer to the Material that this FSR is in  
";

%feature("docstring") Geometry::printString "
printString()  

Prints a string representation of all of the Geometry's attributes to the console.  

This method calls the printString() method for all Materials, Surfaces, Cell, Universes
and Lattices contained by the Geometry.  
";

%feature("docstring") Geometry::getNumFSRs "
getNumFSRs() -> int  

Returns the number of flat source regions in the Geometry.  

Returns
-------
number of FSRs  
";

%feature("docstring") Geometry::findFSRId "
findFSRId(LocalCoords *coords) -> int  

Find and return the ID of the flat source region that a given LocalCoords object resides
within.  

Parameters
----------
* coords :  
    a LocalCoords object pointer  

Returns
-------
the FSR ID for a given LocalCoords object  
";

%feature("docstring") Geometry::getCmfd "
getCmfd() -> Cmfd *  

Returns a pointer to the CMFD object.  

Returns
-------
A pointer to the CMFD object  
";

%feature("docstring") Geometry::withinBounds "
withinBounds(LocalCoords *coords) -> bool  

Determins whether a point is within the bounding box of the geometry.  

Parameters
----------
* coords :  
    a populated LocalCoords linked list  

Returns
-------
boolean indicating whether the coords is within the geometry  
";

%feature("docstring") Geometry::getMaxXBoundaryType "
getMaxXBoundaryType() -> boundaryType  

Returns the boundary conditions (REFLECTIVE or VACUUM) at the maximum x-coordinate in the
Geometry.  

Returns
-------
the boundary conditions for the maximum z-coordinate in the Geometry  
";

%feature("docstring") Geometry::getMinYBoundaryType "
getMinYBoundaryType() -> boundaryType  

Returns the boundary conditions (REFLECTIVE or VACUUM) at the minimum y-coordinate in the
Geometry.  

Returns
-------
the boundary conditions for the minimum y-coordinate in the Geometry  
";

%feature("docstring") Geometry::getAllMaterials "
getAllMaterials() -> std::map< int, Material * >  

Return a std::map container of Material IDs (keys) with Materials pointers (values).  

Returns
-------
a std::map of Materials indexed by Material ID in the geometry  
";

%feature("docstring") Geometry::getAllSurfaces "
getAllSurfaces() -> std::map< int, Surface * >  

Return a std::map container of Surface IDs (keys) with Surfaces pointers (values).  

Returns
-------
a std::map of Surfaces indexed by Surface ID in the geometry  
";

%feature("docstring") Geometry::getMinXBoundaryType "
getMinXBoundaryType() -> boundaryType  

Returns the boundary conditions (REFLECTIVE or VACUUM) at the minimum x-coordinate in the
Geometry.  

Returns
-------
the boundary conditions for the minimum x-coordinate in the Geometry  
";

%feature("docstring") Geometry::getFSRPoint "
getFSRPoint(int fsr_id) -> Point *  

Return the characteristic point for a given FSR ID.  

Parameters
----------
* fsr_id :  
    the FSR ID  

Returns
-------
the FSR's characteristic point  
";

%feature("docstring") Geometry::initializeCmfd "
initializeCmfd()  

This is a method that initializes the CMFD Lattice and sets CMFD parameters.  
";

%feature("docstring") Geometry::toString "
toString() -> std::string  

Converts this Geometry's attributes to a character array.  

This method calls the toString() method for all Surfaces, Cells, Universes and Lattices
contained by the Geometry. Since this routine provides the metadata used by the
TrackGenerator to discriminate between geometries when exporting / importing binary track
files.  

Returns
-------
a character array of this Geometry's class attributes  

Add string data for all Cells  

Add string data for the Surfaces in this Cell  

Add string data for all Universes  
";

%feature("docstring") Geometry::getNumMaterials "
getNumMaterials() -> int  

Returns the number of Materials in the Geometry.  

Returns
-------
the number of Materials  
";

%feature("docstring") Geometry::setFSRCentroid "
setFSRCentroid(int fsr, Point *centroid)  

Sets the centroid for an FSR.  

The _FSR_keys_map stores a hash of a std::string representing the Lattice/Cell/Universe
hierarchy for a unique region and the associated FSR data. _centroid is a point that
represents the numerical centroid of an FSR computed using all segments contained in the
FSR. This method is used by the TrackGenerator to set the centroid after segments have
been created. It is important to note that this method is a helper function for the
TrackGenerator and should not be explicitly called by the user.  

Parameters
----------
* fsr :  
    a FSR ID  
* centroid :  
    a Point representing the FSR centroid  
";

%feature("docstring") Geometry::~Geometry "
~Geometry()  

Destructor clears FSR to Cells and Materials maps.  
";

// File: classGLPolarQuad.xml


%feature("docstring") GLPolarQuad "

Gauss-Legendre's polar quadrature.  

C++ includes: src/Quadrature.h
";

%feature("docstring") GLPolarQuad::GLPolarQuad "
GLPolarQuad()  

Dummy constructor calls the parent constructor.  
";

%feature("docstring") GLPolarQuad::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles (maximum 12)  
";

%feature("docstring") GLPolarQuad::initialize "
initialize()  

Routine to initialize the polar quadrature.  

This routine uses the tabulated values for the Gauss-Legendre polar angle quadrature,
including the sine thetas and weights.  
";

%feature("docstring") GLPolarQuad::precomputeWeights "
precomputeWeights(bool solve_3D)  

Calculates total weights for every azimuthal/polar combination based on the Gauss-Legendre
polar quadrature.  

Parameters
----------
* solve_3D :  
    Boolean indicating whether this is a 3D quadrature  
";

// File: classGPUSolver.xml


%feature("docstring") GPUSolver "

This a subclass of the Solver class for NVIDIA Graphics Processing Units (GPUs).  

The source code for this class includes C++ coupled with compute intensive CUDA kernels
for execution on the GPU.  

C++ includes: openmoc/src/dev/gpu/GPUSolver.h
";

%feature("docstring") GPUSolver::initializeFSRs "
initializeFSRs()  

Initializes the FSR volumes and dev_materials array on the GPU.  

This method assigns each FSR a unique, monotonically increasing ID, sets the Material for
each FSR, and assigns a volume based on the cumulative length of all of the segments
inside the FSR.  
";

%feature("docstring") GPUSolver::zeroTrackFluxes "
zeroTrackFluxes()  

Zero each Track's boundary fluxes for each energy group and polar angle in the \"forward\"
and \"reverse\" directions.  
";

%feature("docstring") GPUSolver::normalizeFluxes "
normalizeFluxes()  

Normalizes all FSR scalar fluxes and Track boundary angular fluxes to the total fission
source (times $ \\nu $).  

Create Thrust vector of fission sources in each FSR  
";

%feature("docstring") GPUSolver::getNumThreadsPerBlock "
getNumThreadsPerBlock() -> int  

Returns the number of threads per block to execute on the GPU.  

Returns
-------
the number of threads per block  
";

%feature("docstring") GPUSolver::computeResidual "
computeResidual(residualType res_type) -> double  

Computes the residual between source/flux iterations.  

Parameters
----------
* res_type :  
    the type of residuals to compute (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)  

Returns
-------
the average residual in each flat source region  
";

%feature("docstring") GPUSolver::setNumThreadsPerBlock "
setNumThreadsPerBlock(int num_threads)  

Sets the number of threads per block (>0) for CUDA kernels.  

Parameters
----------
* num_threads :  
    the number of threads per block  
";

%feature("docstring") GPUSolver::~GPUSolver "
~GPUSolver()  

Solver destructor frees all memory on the device, including arrays for the FSR scalar
fluxes and sources and Track boundary fluxes.  
";

%feature("docstring") GPUSolver::setNumThreadBlocks "
setNumThreadBlocks(int num_blocks)  

Sets the number of thread blocks (>0) for CUDA kernels.  

Parameters
----------
* num_blocks :  
    the number of thread blocks  
";

%feature("docstring") GPUSolver::getNumThreadBlocks "
getNumThreadBlocks() -> int  

Returns the number of thread blocks to execute on the GPU.  

Returns
-------
the number of thread blocks  
";

%feature("docstring") GPUSolver::addSourceToScalarFlux "
addSourceToScalarFlux()  

Add the source term contribution in the transport equation to the FSR scalar flux.  
";

%feature("docstring") GPUSolver::initializeTracks "
initializeTracks()  

Allocates memory for all Tracks on the GPU.  
";

%feature("docstring") GPUSolver::storeFSRFluxes "
storeFSRFluxes()  

Stores the FSR scalar fluxes in the old scalar flux array.  
";

%feature("docstring") GPUSolver::getFSRSource "
getFSRSource(int fsr_id, int group) -> FP_PRECISION  

Returns the source for some energy group for a flat source region.  

This is a helper routine used by the openmoc.process module.  

Parameters
----------
* fsr_id :  
    the ID for the FSR of interest  
* group :  
    the energy group of interest  

Returns
-------
the flat source region source  
";

%feature("docstring") GPUSolver::initializeMaterials "
initializeMaterials(solverMode mode=ADJOINT)  

Allocates all Materials data on the GPU.  

This method loops over the materials in the host_materials map. Since CUDA does not
support std::map data types on the device, the materials map must be converted to an array
and a map created that maps a material ID to an indice in the new materials array. In
initializeTracks, this map is used to convert the Material ID associated with every
segment to an index in the materials array.  

Parameters
----------
* mode :  
    the solution type (FORWARD or ADJOINT)  
";

%feature("docstring") GPUSolver::GPUSolver "
GPUSolver(TrackGenerator *track_generator=NULL)  

Constructor initializes arrays for dev_tracks and dev_materials..  

The constructor initalizes the number of CUDA threads and thread blocks each to a default
of 64.  

Parameters
----------
* track_generator :  
    an optional pointer to the TrackjGenerator  
";

%feature("docstring") GPUSolver::getFlux "
getFlux(int fsr_id, int group) -> FP_PRECISION  

Returns the scalar flux for some FSR and energy group.  

Parameters
----------
* fsr_id :  
    the ID for the FSR of interest  
* group :  
    the energy group of interest  

Returns
-------
the FSR scalar flux  
";

%feature("docstring") GPUSolver::computeFSRFissionSources "
computeFSRFissionSources()  

Computes the fission source in each FSR.  

This method computes the fission source in each FSR based on this iteration's current
approximation to the scalar flux.  
";

%feature("docstring") GPUSolver::computeKeff "
computeKeff()  

Compute $ k_{eff} $ from successive fission sources.  

This method computes the current approximation to the multiplication factor on this
iteration as follows: $ k_{eff} = \\frac{\\displaystyle\\sum_{i \\in I}
\\displaystyle\\sum_{g \\in G} \\nu \\Sigma^F_g \\Phi V_{i}} {\\displaystyle\\sum_{i \\in
I} \\displaystyle\\sum_{g \\in G} (\\Sigma^T_g \\Phi V_{i} - \\Sigma^S_g \\Phi V_{i} -
L_{i,g})} $  
";

%feature("docstring") GPUSolver::transportSweep "
transportSweep()  

This method performs one transport sweep of all azimuthal angles, Tracks, Track segments,
polar angles and energy groups.  

The method integrates the flux along each Track and updates the boundary fluxes for the
corresponding output Track, while updating the scalar flux in each flat source region.  
";

%feature("docstring") GPUSolver::computeFSRSources "
computeFSRSources()  

Computes the total source (fission, scattering, fixed) in each FSR.  

This method computes the total source in each FSR based on this iteration's current
approximation to the scalar flux.  
";

%feature("docstring") GPUSolver::initializeFixedSources "
initializeFixedSources()  

Populates array of fixed sources assigned by FSR.  
";

%feature("docstring") GPUSolver::initializeExpEvaluator "
initializeExpEvaluator()  

does nothing.
";

%feature("docstring") GPUSolver::computeFSRScatterSources "
computeFSRScatterSources()  

Computes the scatter source in each FSR.  

This method computes the scatter source in each FSR based on this iteration's current
approximation to the scalar flux.  
";

%feature("docstring") GPUSolver::initializeFluxArrays "
initializeFluxArrays()  

Allocates memory for Track boundary angular and FSR scalar fluxes.  

Deletes memory for old flux vectors if they were allocated for a previous simulation.  
";

%feature("docstring") GPUSolver::setGeometry "
setGeometry(Geometry *geometry)  

Sets the Geometry for the Solver.  

This is a private setter method for the Solver and is not intended to be called by the
user.  

Parameters
----------
* geometry :  
    a pointer to a Geometry object  
";

%feature("docstring") GPUSolver::getFluxes "
getFluxes(FP_PRECISION *out_fluxes, int num_fluxes)  

Fills an array with the scalar fluxes on the GPU.  

This class method is a helper routine called by the OpenMOC Python \"openmoc.krylov\"
module for Krylov subspace methods. Although this method appears to require two arguments,
in reality it only requires one due to SWIG and would be called from within Python as
follows:  


Parameters
----------
* fluxes :  
    an array of FSR scalar fluxes in each energy group  
* num_fluxes :  
    the total number of FSR flux values  
";

%feature("docstring") GPUSolver::setFluxes "
setFluxes(FP_PRECISION *in_fluxes, int num_fluxes)  

Set the flux array for use in transport sweep source calculations.  This is a helper
method for the checkpoint restart capabilities, as well as the IRAMSolver in the
openmoc.krylov submodule. This routine may be used as follows from within Python:  

  

         NOTE: This routine stores a pointer to the fluxes for the Solver
         to use during transport sweeps and other calculations. Hence, the
         flux array pointer is shared between NumPy and the Solver.  

Parameters
----------
* in_fluxes :  
    an array with the fluxes to use  
* num_fluxes :  
    the number of flux values (# groups x # FSRs)  
";

%feature("docstring") GPUSolver::flattenFSRFluxes "
flattenFSRFluxes(FP_PRECISION value)  

Set the scalar flux for each FSR and energy group to some value.  

Parameters
----------
* value :  
    the value to assign to each FSR scalar flux  
";

%feature("docstring") GPUSolver::computeFSRFissionRates "
computeFSRFissionRates(double *fission_rates, int num_FSRs)  

Computes the volume-averaged, energy-integrated nu-fission rate in each FSR and stores
them in an array indexed by FSR ID.  

This is a helper method for SWIG to allow users to retrieve FSR nu-fission rates as a
NumPy array. An example of how this method can be called from Python is as follows:  


Parameters
----------
* fission_rates :  
    an array to store the nu-fission rates (implicitly passed in as a NumPy array from
    Python)  
* num_FSRs :  
    the number of FSRs passed in from Python  
";

%feature("docstring") GPUSolver::initializeSourceArrays "
initializeSourceArrays()  

Allocates memory for FSR source vectors on the GPU.  

Deletes memory for old source vectors if they were allocated for a previous simulation.  
";

%feature("docstring") GPUSolver::setTrackGenerator "
setTrackGenerator(TrackGenerator *track_generator)  

Sets the Solver's TrackGenerator with characteristic Tracks.  

The TrackGenerator must already have generated Tracks and have used ray tracing to
segmentize them across the Geometry. This should be initated in Python prior to assigning
the TrackGenerator to the Solver:  


Parameters
----------
* track_generator :  
    a pointer to a TrackGenerator object  
";

// File: structisinf__test.xml


%feature("docstring") isinf_test "

A struct used to check if a value on the GPU is equal to INF.  

This is used as a predicate in Thrust routines.  
";

// File: structisnan__test.xml


%feature("docstring") isnan_test "

A struct used to check if a value on the GPU is equal to NaN.  

This is used as a predicate in Thrust routines.  
";

// File: classLattice.xml


%feature("docstring") Lattice "

Represents a repeating 3D Lattice of Universes.  

C++ includes: src/Universe.h
";

%feature("docstring") Lattice::withinBounds "
withinBounds(Point *point) -> bool  

Checks if a Point is within the bounds of a Lattice.  

Parameters
----------
* point :  
    a pointer to the Point of interest  

Returns
-------
true if the Point is in the bounds, false if not  
";

%feature("docstring") Lattice::printString "
printString()  

Prints a string representation of all of the Lattice's attributes to the console.  
";

%feature("docstring") Lattice::subdivideCells "
subdivideCells(double max_radius=INFINITY)  

Subdivides all of the Material-filled Cells within this Lattice into rings and angular
sectors aligned with the z-axis.  

Parameters
----------
* max_radius :  
    the maximum allowable radius used in the subdivisions  
";

%feature("docstring") Lattice::getUniverses "
getUniverses() -> std::vector< std::vector< std::vector< std::pair< int, Universe * > > >
    > *  

Return a 3D vector of the Universes in the Lattice.  

Returns
-------
3D vector of Universes  
";

%feature("docstring") Lattice::getMinZ "
getMinZ() -> double  

Returns the minimum reachable z-coordinate in the Lattice.  

Returns
-------
the minimum reachable z-coordinate  
";

%feature("docstring") Lattice::getMinX "
getMinX() -> double  

Returns the minimum reachable x-coordinate in the Lattice.  

Returns
-------
the minimum reachable x-coordinate  
";

%feature("docstring") Lattice::getMinY "
getMinY() -> double  

Returns the minimum reachable y-coordinate in the Lattice.  

Returns
-------
the minimum reachable y-coordinate  
";

%feature("docstring") Lattice::minSurfaceDist "
minSurfaceDist(LocalCoords *coords) -> double  

Finds the distance to the nearest surface.  

Knowing that a Lattice must be cartesian, this function computes the distance to the
nearest boundary between lattice cells in the direction of the track.  

Parameters
----------
* coords :  
    a pointer to a localcoords object  

Returns
-------
the distance to the nearest Lattice cell boundary  
";

%feature("docstring") Lattice::findCell "
findCell(LocalCoords *coords) -> Cell *  

Finds the Cell within this Lattice that a LocalCoords is in.  

This method first find the Lattice cell, then searches the Universe inside that Lattice
cell. If LocalCoords is outside the bounds of the Lattice, this method will return NULL.  

Parameters
----------
* coords :  
    the LocalCoords of interest  

Returns
-------
a pointer to the Cell this LocalCoord is in or NULL  
";

%feature("docstring") Lattice::getLatZ "
getLatZ(Point *point) -> int  

Finds the Lattice cell z index that a point lies in.  

Parameters
----------
* point :  
    a pointer to a point being evaluated.  

Returns
-------
the Lattice cell z index.  
";

%feature("docstring") Lattice::getLatY "
getLatY(Point *point) -> int  

Finds the Lattice cell y index that a point lies in.  

Parameters
----------
* point :  
    a pointer to a point being evaluated.  

Returns
-------
the Lattice cell y index.  
";

%feature("docstring") Lattice::setUniverses "
setUniverses(int num_z, int num_y, int num_x, Universe **universes)  

Sets the array of Universe pointers filling each Lattice cell.  

This is a helper method for SWIG to allow users to assign Universes to a Lattice using a
3D Python list (list of lists of lists). An example how this method can be called from
Python is as follows:  


Parameters
----------
* num_z :  
    the number of Lattice cells along z  
* num_y :  
    the number of Lattice cells along y  
* num_x :  
    the number of Lattice cells along x  
* universes :  
    the array of Universes for each Lattice cell  
";

%feature("docstring") Lattice::removeUniverse "
removeUniverse(Universe *universe)  

Removes all references to a Universe from the Lattice.  

Parameters
----------
* universe :  
    the Universe to remove  
";

%feature("docstring") Lattice::getUniverse "
getUniverse(int lat_x, int lat_y, int lat_z) const  -> Universe *  

Returns a pointer to the Universe within a specific Lattice cell.  

Parameters
----------
* lat_x :  
    the x index to the Lattice cell  
* lat_y :  
    the y index to the Lattice cell  
* lat_z :  
    the z index to the Lattice cell  

Returns
-------
pointer to a Universe filling the Lattice cell  
";

%feature("docstring") Lattice::buildNeighbors "
buildNeighbors()  

Builds collections of neighboring Cells for all Cells in each Universe in the Lattice for
optimized ray tracing.  
";

%feature("docstring") Lattice::Lattice "
Lattice(const int id=-1, const char *name=\"\")  

Constructor sets the user-specified and unique IDs for this Lattice.  

Parameters
----------
* id :  
    the user-specified optional Lattice (Universe) ID  
* name :  
    the user-specified optional Lattice (Universe) name  
";

%feature("docstring") Lattice::setNumZ "
setNumZ(int num_z)  

Set the number of Lattice cells along the z-axis.  

Parameters
----------
* num_z :  
    the number of Lattice cells along z  
";

%feature("docstring") Lattice::setNumX "
setNumX(int num_x)  

Set the number of Lattice cells along the x-axis.  

Parameters
----------
* num_x :  
    the number of Lattice cells along x  
";

%feature("docstring") Lattice::setNumY "
setNumY(int num_y)  

Set the number of Lattice cells along the y-axis.  

Parameters
----------
* num_y :  
    the number of Lattice cells along y  
";

%feature("docstring") Lattice::getLatX "
getLatX(Point *point) -> int  

Finds the Lattice cell x index that a point lies in.  

Parameters
----------
* point :  
    a pointer to a point being evaluated.  

Returns
-------
the Lattice cell x index.  
";

%feature("docstring") Lattice::setOffset "
setOffset(double x, double y, double z)  

Set the offset in global coordinates for this Lattice.  

A lattice is assumed to be a rectilinear grid with the center/origin of the grid located
in the center of the Lattice's parent universe. The offset represents the offset of the
lattice center/origin with respect to the center of the parent universe. Therefore an
offset of (-1,2) would move the center/origin of the lattice to the left 1 cm and up 2 cm.  

Parameters
----------
* x :  
    the offset in the x direction  
* y :  
    the offset in the y direction  
* z :  
    the offset in the z direction  
";

%feature("docstring") Lattice::getWidthX "
getWidthX() const  -> double  

Return the width of the Lattice along the x-axis.  

Returns
-------
the width of the Lattice cells along x  
";

%feature("docstring") Lattice::getWidthY "
getWidthY() const  -> double  

Return the width of the Lattice along the y-axis.  

Returns
-------
the width of the Lattice cells along y  
";

%feature("docstring") Lattice::getWidthZ "
getWidthZ() const  -> double  

Return the width of the Lattice along the z-axis.  

Returns
-------
the width of the Lattice cells along z  
";

%feature("docstring") Lattice::getDistanceToSurface "
getDistanceToSurface(int cell, Point *point, int surface) -> double  

Finds the distance from a point to a particular lattice cell surface.  

Parameters
----------
* cell :  
    the cell index that the point is in.  
* point :  
    a pointer to a point being evaluated.  
* surface :  
    a surface id to get the distance to.  

Returns
-------
the distance to the lattice cell surface of interest.  
";

%feature("docstring") Lattice::getOffset "
getOffset() -> Point *  

Return a pointer to the offset for this Cell (in global coordinates).  

Returns
-------
the offset of the Cell  
";

%feature("docstring") Lattice::updateUniverse "
updateUniverse(int lat_x, int lat_y, int lat_z, Universe *universe)  

Update the Universe in a particular Lattice cell.  

This method may only be used after an array of Universes has been assigned with the
Lattice::setUniverses(...) method.  

Parameters
----------
* lat_x :  
    the Lattice cell index along x  
* lat_y :  
    the Lattice cell index along y  
* lat_z :  
    the Lattice cell index along z  
* universe :  
    the Universe to insert into the Lattice  
";

%feature("docstring") Lattice::setWidth "
setWidth(double width_x, double width_y, double width_z=std::numeric_limits< double
    >::infinity())  

Set the width of each Lattice cell.  

Parameters
----------
* width_x :  
    the width along the x-axis in centimeters  
* width_y :  
    the width along the y-axis in centimeters  
* width_z :  
    the width along the z-axis in centimeters  
";

%feature("docstring") Lattice::getUniqueUniverses "
getUniqueUniverses() -> std::map< int, Universe * >  

Aggregates a list (vector) of the IDs of all Universes within the FILL type Cells filling
this Universe.  

Note that this method only searches the first level of Cells below this Universe within
the nested Universe coordinate system.  

Returns
-------
a vector of Universe IDs  
";

%feature("docstring") Lattice::getLatticeSurface "
getLatticeSurface(int cell, Point *point) -> int  

Finds the Lattice cell surface that a point lies on. If the point is not on a surface, -1
is returned.  

The surface indices are defined in constants.h as they need to be consistent with the
surface constant definitions used in Cmfd. The index returned takes into account the cell
index and returns NUM_SURFACES*cell_index + surface_index.  

Parameters
----------
* cell :  
    the cell index that the point is in.  
* point :  
    a pointer to a point being evaluated.  

Returns
-------
the Lattice surface index.  
";

%feature("docstring") Lattice::getAllUniverses "
getAllUniverses() -> std::map< int, Universe * >  

Returns the std::map of all nested Universe IDs and Universe pointers filling this
Lattice.  

Returns
-------
std::map of Universe IDs and pointers  
";

%feature("docstring") Lattice::getNumZ "
getNumZ() const  -> int  

Return the number of Lattice cells along the z-axis.  

Returns
-------
the number of Lattice cells along z  
";

%feature("docstring") Lattice::getNumX "
getNumX() const  -> int  

Return the number of Lattice cells along the x-axis.  

Returns
-------
the number of Lattice cells along x  
";

%feature("docstring") Lattice::getNumY "
getNumY() const  -> int  

Return the number of Lattice cells along the y-axis.  

Returns
-------
the number of Lattice cells along y  
";

%feature("docstring") Lattice::getAllCells "
getAllCells() -> std::map< int, Cell * >  

Returns the std::map of Cell IDs and Cell pointers in this Lattice at all nested Universe
levels.  

Returns
-------
std::map of Cell IDs and pointers  
";

%feature("docstring") Lattice::getLatticeCell "
getLatticeCell(Point *point) -> int  

Finds the Lattice cell index that a point lies in.  

Lattice cells are numbered starting with 0 x-min/y-min/z-min corner. Lattice cell IDs then
increase monotonically from x-min to x-max, y-min to y-max, and z-min to z-max. Note that
values increase first on the x-axis, followed by the y-axis, then on the z-axis. For
example, the indices for a 4 x 4 x 1 lattice: 12 13 14 15 8 9 10 11 4 5 6 7 0 1 2 3  

Parameters
----------
* point :  
    a pointer to a point being evaluated.  

Returns
-------
the Lattice cell index.  
";

%feature("docstring") Lattice::toString "
toString() -> std::string  

Converts a Lattice's attributes to a character array representation.  

Returns
-------
character array of this Lattice's attributes  
";

%feature("docstring") Lattice::getMaxX "
getMaxX() -> double  

Returns the maximum reachable x-coordinate in the Lattice.  

Returns
-------
the maximum reachable x-coordinate  
";

%feature("docstring") Lattice::getMaxY "
getMaxY() -> double  

Returns the maximum reachable y-coordinate in the Lattice.  

Returns
-------
the maximum reachable y-coordinate  
";

%feature("docstring") Lattice::getMaxZ "
getMaxZ() -> double  

Returns the maximum reachable z-coordinate in the Lattice.  

Returns
-------
the maximum reachable z-coordinate  
";

%feature("docstring") Lattice::~Lattice "
~Lattice()  

Destructor clears memory for all of Universes pointers.  
";

// File: classLeonardPolarQuad.xml


%feature("docstring") LeonardPolarQuad "

Leonard's polar quadrature.  

C++ includes: src/Quadrature.h
";

%feature("docstring") LeonardPolarQuad::initialize "
initialize()  

Routine to initialize the polar quadrature.  

This routine uses the tabulated values for the Leonard polar angle quadrature, including
the sine thetas and weights.  
";

%feature("docstring") LeonardPolarQuad::precomputeWeights "
precomputeWeights(bool solve_3D)  

Calculates total weights for every azimuthal/polar combination based on the Leonard polar
quadrature.  

Parameters
----------
* solve_3D :  
    Boolean indicating whether this is a 3D quadrature  
";

%feature("docstring") LeonardPolarQuad::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles (4 or 6)  
";

%feature("docstring") LeonardPolarQuad::LeonardPolarQuad "
LeonardPolarQuad()  

Dummy constructor calls the parent constructor.  
";

// File: classLocalCoords.xml


%feature("docstring") LocalCoords "

The LocalCoords represents a set of local coordinates on some level of nested Universes
making up the geometry.  

C++ includes: openmoc/src/host/LocalCoords.h
";

%feature("docstring") LocalCoords::getPhi "
getPhi() const  -> double  

Returns the direction angle in radians with respect to the x-axis.  

Returns
-------
the direction angle in radians  
";

%feature("docstring") LocalCoords::~LocalCoords "
~LocalCoords()  

Destructor.  
";

%feature("docstring") LocalCoords::setPhi "
setPhi(double phi)  

Set the direction angle in radians for this LocalCoords.  

Parameters
----------
* angle :  
    the direction angle in radians  
";

%feature("docstring") LocalCoords::getLattice "
getLattice() const  -> Lattice *  

Return the Lattice within which this LocalCoords resides.  

Returns
-------
the Lattice  
";

%feature("docstring") LocalCoords::setCell "
setCell(Cell *cell)  

Set the Cell within which this LocalCoords resides.  

Parameters
----------
* cell :  
    the Cell  
";

%feature("docstring") LocalCoords::getPrev "
getPrev() const  -> LocalCoords *  

Return a pointer to the LocalCoord at the next higher nested Universe level if one exists.  

Returns
-------
pointer to the previous LocalCoord  
";

%feature("docstring") LocalCoords::getCell "
getCell() const  -> Cell *  

Return the Cell within which this LocalCoords resides.  

Returns
-------
the Cell  
";

%feature("docstring") LocalCoords::prune "
prune()  

Removes and frees memory for all LocalCoords beyond this one in the linked list.  
";

%feature("docstring") LocalCoords::setNext "
setNext(LocalCoords *next)  

Sets the pointer to the LocalCoords on the next lower nested Universe level.  

Parameters
----------
* next :  
    pointer to the next LocalCoords  
";

%feature("docstring") LocalCoords::getZ "
getZ() const  -> double  

Returns the z-coordinate for this LocalCoords location.  

Returns
-------
the z-coordinate of this LocalCoords location  
";

%feature("docstring") LocalCoords::getNext "
getNext() const  -> LocalCoords *  

Return a pointer to the LocalCoord at the next lower nested Universe level if one exists.  

Returns
-------
pointer to the next LocalCoord  
";

%feature("docstring") LocalCoords::setZ "
setZ(double z)  

Set the z-coordinate for this Localcoords.  

Parameters
----------
* z :  
    the z-coordinate  
";

%feature("docstring") LocalCoords::setX "
setX(double x)  

Set the x-coordinate for this LocalCoords.  

Parameters
----------
* x :  
    the x-coordinate  
";

%feature("docstring") LocalCoords::setY "
setY(double y)  

Set the y-coordinate for this Localcoords.  

Parameters
----------
* y :  
    the y-coordinate  
";

%feature("docstring") LocalCoords::getLowestLevel "
getLowestLevel() -> LocalCoords *  

Find and return the last LocalCoords in the linked list which represents the local
coordinates on the lowest level of a geometry of nested universes.  

Traverses a linked list of LocalCoords to find the one at the lowest nested Universe
level.  

Returns
-------
a pointer to the last LocalCoords object in the list  
";

%feature("docstring") LocalCoords::setType "
setType(coordType type)  

Set the type of LocalCoords (UNIV or LAT).  

Parameters
----------
* type :  
    the type for LocalCoords (UNIV or LAT)  
";

%feature("docstring") LocalCoords::setLatticeZ "
setLatticeZ(int lattice_z)  

Sets the z index for the Lattice cell within which this LocalCoords resides.  

Parameters
----------
* lattice_z :  
    the z Lattice cell index  
";

%feature("docstring") LocalCoords::setLatticeY "
setLatticeY(int lattice_y)  

Sets the y index for the Lattice cell within which this LocalCoords resides.  

Parameters
----------
* lattice_y :  
    the y Lattice cell index  
";

%feature("docstring") LocalCoords::getLatticeZ "
getLatticeZ() const  -> int  

Return the third index of the Lattice cell within which this LocalCoords resides.  

Returns
-------
the third Lattice cell index  
";

%feature("docstring") LocalCoords::getLatticeY "
getLatticeY() const  -> int  

Return the second index of the Lattice cell within which this LocalCoords resides.  

Returns
-------
the second Lattice cell index  
";

%feature("docstring") LocalCoords::getLatticeX "
getLatticeX() const  -> int  

Return the first index of the Lattice cell within which this LocalCoords resides.  

Returns
-------
the first Lattice cell index  
";

%feature("docstring") LocalCoords::getUniverse "
getUniverse() const  -> Universe *  

Return the Universe within which this LocalCoords resides.  

Returns
-------
the Universe  
";

%feature("docstring") LocalCoords::updateMostLocal "
updateMostLocal(Point *point)  

Update the last element in the linked list (the one at the lowest level of nested
Universes) to have the same coordinates as a given point.  

Parameters
----------
* point :  
    a pointer to a point of interest  
";

%feature("docstring") LocalCoords::toString "
toString() -> std::string  

Converts this LocalCoords's attributes to a character array representation.  

Returns
-------
a character array of the LocalCoord's attributes  
";

%feature("docstring") LocalCoords::incrementPhi "
incrementPhi(double phi)  

Increment the direction angle in radians for this LocalCoords.  

Parameters
----------
* phi :  
    the incremental direction angle in radians  
";

%feature("docstring") LocalCoords::copyCoords "
copyCoords(LocalCoords *coords)  

Copies a LocalCoords' values to this one. details Given a pointer to a LocalCoords, it
first prunes it and then creates a copy of the linked list of LocalCoords in the linked
list below this one to give to the input LocalCoords.  

Parameters
----------
* coords :  
    a pointer to the LocalCoords to give the linked list copy to  
";

%feature("docstring") LocalCoords::setLattice "
setLattice(Lattice *lattice)  

Sets the Lattice within which this LocalCoords resides.  

Parameters
----------
* lattice :  
    the Lattice  
";

%feature("docstring") LocalCoords::adjustCoords "
adjustCoords(double delta)  

Translate all of the x,y coordinates for each LocalCoords object in the linked list.  

This method will traverse the entire linked list and apply the translation to each
element.  

Parameters
----------
* delta :  
    amount we wish to move the point by  
";

%feature("docstring") LocalCoords::getHighestLevel "
getHighestLevel() -> LocalCoords *  

Find and return the first LocalCoords in the linked list which represents the local
coordinates on the highest level of a geometry of nested universes.  

Traverses a linked list of LocalCoords to find the one at the highest nested Universe
level.  

Returns
-------
a pointer to the first LocalCoords object in the list  
";

%feature("docstring") LocalCoords::setUniverse "
setUniverse(Universe *universe)  

Set the Universe within which this LocalCoords resides.  

Parameters
----------
* universe :  
    the Universe  
";

%feature("docstring") LocalCoords::getX "
getX() const  -> double  

Returns the x-coordinate for this LocalCoords location.  

Returns
-------
the x-coordinate of this LocalCoords location  
";

%feature("docstring") LocalCoords::getY "
getY() const  -> double  

Returns the y-coordinate for this LocalCoords location.  

Returns
-------
the y-coordinate of this LocalCoords location  
";

%feature("docstring") LocalCoords::getPoint "
getPoint() -> Point *  

Returns a pointer to the Point containing the coordinates for this LocalCoord.  

Returns
-------
pointer to the Point containing the x and y coordinates  
";

%feature("docstring") LocalCoords::getType "
getType() -> coordType  

Return the level (UNIV or LAT) of this LocalCoords.  

Returns
-------
the nested Universe level (UNIV or LAT)  
";

%feature("docstring") LocalCoords::LocalCoords "
LocalCoords(double x, double y, double z)  

Constructor sets the x and y coordinates.  

Parameters
----------
* x :  
    the x-coordinate  
* y :  
    the y-coordinate  
";

%feature("docstring") LocalCoords::setLatticeX "
setLatticeX(int lattice_x)  

Sets the x index for the Lattice cell within which this LocalCoords resides.  

Parameters
----------
* lattice_x :  
    the x Lattice cell index  
";

%feature("docstring") LocalCoords::setPrev "
setPrev(LocalCoords *coords)  

Sets the pointer to the LocalCoords on the next higher nested Universe level.  

Parameters
----------
* prev :  
    pointer to the previous LocalCoords  
";

// File: classMaterial.xml


%feature("docstring") Material "

The Material class represents a unique material and its relevant nuclear data (i.e.,
multigroup cross-sections) for neutron transport.  

C++ includes: src/Material.h
";

%feature("docstring") Material::setChiByGroup "
setChiByGroup(double xs, int group)  

Set the Material's chi value for some energy group.  

Parameters
----------
* xs :  
    the chi value ( $ \\Chi $)  
* group :  
    the energy group  
";

%feature("docstring") Material::getChiByGroup "
getChiByGroup(int group) -> FP_PRECISION  

Get the Material's fission spectrum for some energy group.  

Parameters
----------
* group :  
    the energy group  

Returns
-------
the fission spectrum  
";

%feature("docstring") Material::clone "
clone() -> Material *  

Create a duplicate of the Material.  

Returns
-------
a pointer to the clone  
";

%feature("docstring") Material::getSigmaFByGroup "
getSigmaFByGroup(int group) -> FP_PRECISION  

Get the Material's fission cross section for some energy group.  

Parameters
----------
* group :  
    the energy group  

Returns
-------
the fission cross section  
";

%feature("docstring") Material::setSigmaF "
setSigmaF(double *xs, int num_groups)  

Set the Material's array of fission cross-sections.  

This method is a helper function to allow OpenMOC users to assign the Material's nuclear
data in Python. A user must initialize a NumPy array of the correct size (e.g., a float64
array the length of the number of energy groups) as input to this function. This function
then fills the NumPy array with the data values for the Material's fission cross-sections.
An example of how this function might be called in Python is as follows:  


Parameters
----------
* xs :  
    the array of fission cross-sections  
* num_groups :  
    the number of energy groups  
";

%feature("docstring") Material::toString "
toString() -> std::string  

Converts this Material's attributes to a character array representation.  

The character array returned includes the user-defined ID, and each of the absorption,
total, fission, nu multiplied by fission and scattering cross-sections and chi for all
energy groups.  

Returns
-------
character array of this Material's attributes  
";

%feature("docstring") Material::buildFissionMatrix "
buildFissionMatrix()  

Builds the fission matrix from chi and the fission cross-section.  

The fission matrix is constructed as the outer product of the chi and fission cross-
section vectors. This routine is intended for internal use and is called by the Solver at
runtime.  
";

%feature("docstring") Material::getId "
getId() const  -> int  

Return the Material's user-defined ID.  

Returns
-------
the Material's user-defined ID  
";

%feature("docstring") Material::setSigmaSByGroup "
setSigmaSByGroup(double xs, int origin, int destination)  

Set the Material's scattering cross-section for some energy group.  

Parameters
----------
* xs :  
    the scattering cross-section  
* origin :  
    the column index in the scattering matrix  
* destination :  
    the row index in the scattering matrix  
";

%feature("docstring") Material::setSigmaS "
setSigmaS(double *xs, int num_groups)  

Set the Material's 2D array of scattering cross-sections.  

The array should be passed to OpenMOC as a 1D array in column-major order. This assumes
the standard convention, where column index is the origin group and the row index is the
destination group. That is, the array should be ordered as follows: 1 -> 1 1 -> 2 1 -> 3
... 2 -> 1 2 -> 2 ...  

Note that if the scattering matrix is defined in NumPy by the standard convention,
\"flat\" will put the matrix into row major order. Thus, one should transpose the matrix
before flattening.  

For cache efficiency, the transpose of the input is actually stored in OpenMOC.  

This method is a helper function to allow OpenMOC users to assign the Material's nuclear
data in Python. A user must initialize a NumPy array of the correct size (e.g., a float64
array the length of the square of the number of energy groups) as input to this function.
This function then fills the NumPy array with the data values for the Material's
scattering cross-sections. An example of how this function might be called in Python is as
follows:  


Parameters
----------
* xs :  
    the array of scattering cross-sections  
* num_groups_squared :  
    the number of energy groups squared  
";

%feature("docstring") Material::getFissionMatrix "
getFissionMatrix() -> FP_PRECISION *  

Return the array of the Material's fission matrix.  

Returns
-------
the pointer to the Material's fission matrix array  
";

%feature("docstring") Material::setSigmaT "
setSigmaT(double *xs, int num_groups)  

Set the Material's array of total cross-sections.  

This method is a helper function to allow OpenMOC users to assign the Material's nuclear
data in Python. A user must initialize a NumPy array of the correct size (e.g., a float64
array the length of the number of energy groups) as input to this function. This function
then fills the NumPy array with the data values for the Material's total cross-sections.
An example of how this function might be called in Python is as follows:  

  

         NOTE: This routine will override an zero-valued cross-sections
         (e.g., in void or gap regions) with a minimum value of 1E-10 to
         void numerical issues in the MOC solver.  

Parameters
----------
* xs :  
    the array of total cross-sections  
* num_groups :  
    the number of energy groups  
";

%feature("docstring") Material::isDataAligned "
isDataAligned() -> bool  

Returns true if the data is vector aligned, false otherwise (default).  

Returns
-------
Whether or not the Material's data is vector aligned  
";

%feature("docstring") Material::setNumEnergyGroups "
setNumEnergyGroups(const int num_groups)  

Set the number of energy groups for this Material.  

Parameters
----------
* num_groups :  
    the number of energy groups.  
";

%feature("docstring") Material::alignData "
alignData()  

Reallocates the Material's cross-section data structures along word-aligned boundaries.  

This method is used to assist with SIMD auto-vectorization of the MOC routines in the
Solver classes. Rather than using the assigned number of energy groups, this method adds
\"dummy\" energy groups such that the total number of groups is some multiple of
VEC_LENGTH (typically 4, 8, or 16). As a result, the SIMD-vectorized Solver subclasses can
break up loops over energy groups in such a way to \"expose\" the SIMD nature of the
algorithm.  
";

%feature("docstring") Material::getName "
getName() const  -> char *  

Return the user-defined name of the Material.  

Returns
-------
the Material name  
";

%feature("docstring") Material::printString "
printString()  

Prints a string representation of all of the Material's attributes to the console.  
";

%feature("docstring") Material::~Material "
~Material()  

Destructor deletes all cross-section data structures from memory.  
";

%feature("docstring") Material::getSigmaSByGroup "
getSigmaSByGroup(int origin, int destination) -> FP_PRECISION  

Get the Material's scattering cross section for some energy group.  

Parameters
----------
* origin :  
    the incoming energy group  
* destination :  
    the outgoing energy group  

Returns
-------
the scattering cross section  
";

%feature("docstring") Material::getNumInstances "
getNumInstances() -> int  

Return the number of instances of this Material in the Geometry.  

The number of instances of this Material in the Geometry is determined during track
generation.  

Returns
-------
the number of material instances  
";

%feature("docstring") Material::getNuSigmaFByGroup "
getNuSigmaFByGroup(int group) -> FP_PRECISION  

Get the Material's nu-fission cross section for some energy group.  

Parameters
----------
* group :  
    the energy group  

Returns
-------
the nu-fission cross section  
";

%feature("docstring") Material::getNumEnergyGroups "
getNumEnergyGroups() const  -> int  

Returns the number of energy groups for this Material's nuclear data.  

Returns
-------
the number of energy groups  
";

%feature("docstring") Material::isFissionable "
isFissionable() -> bool  

Returns whether or not the Material contains a fissionable (non-zero) fission cross-
section.  

Returns
-------
true if fissionable, false otherwise  
";

%feature("docstring") Material::setNumInstances "
setNumInstances(int num_instances)  

Set the number of instances of this Material.  

Parameters
----------
* num_instances :  
    the number of instances of this Material in the Geometry  
";

%feature("docstring") Material::getFissionMatrixByGroup "
getFissionMatrixByGroup(int origin, int destination) -> FP_PRECISION  

Get the Material's fission matrix for some energy group.  

Parameters
----------
* origin :  
    the incoming energy group $ E_{0} $  
* destination :  
    the outgoing energy group $ E_{1} $  

Returns
-------
the fission matrix entry $ \\nu\\Sigma_{f}(E_{0}) * \\chi(E_{1})$  
";

%feature("docstring") Material::getVolume "
getVolume() -> double  

Return the aggregate volume/area of all instances of this Material.  

The volume/area of the Material is computed from track segments which overlap this
Material during track generation.  

Returns
-------
the volume/area of the Material  
";

%feature("docstring") Material::setName "
setName(const char *name)  

Sets the name of the Material.  

Parameters
----------
* name :  
    the Material name string  
";

%feature("docstring") Material::incrementNumInstances "
incrementNumInstances()  

Increment the number of instances of this Material.  

This routine is called by the TrackGenerator during track generation and segmentation.  
";

%feature("docstring") Material::getChi "
getChi() -> FP_PRECISION *  

Return the array of the Material's chi $ \\chi $.  

Returns
-------
the pointer to the Material's array of chi $ \\chi $ values  
";

%feature("docstring") Material::setNuSigmaF "
setNuSigmaF(double *xs, int num_groups)  

Set the Material's array of fission cross-sections multiplied by $ \\nu $.  

Parameters
----------
* xs :  
    the array of fission cross-sections multiplied by nu $ \\nu $  
* num_groups :  
    the number of energy groups  
";

%feature("docstring") Material::setSigmaTByGroup "
setSigmaTByGroup(double xs, int group)  

Set the Material's total cross-section for some energy group.  

Parameters
----------
* xs :  
    the total cross-section  
* group :  
    the energy group  
";

%feature("docstring") Material::Material "
Material(int id=0, const char *name=\"\")  

Constructor sets the ID and unique ID for the Material.  

Parameters
----------
* id :  
    the user-specified optional Material ID  
* name :  
    the user-specified optional Material name  
";

%feature("docstring") Material::incrementVolume "
incrementVolume(double volume)  

Increment the volume/area of the Material by some amount.  

This routine is called by the TrackGenerator during track generation and segmentation.  

Parameters
----------
* volume :  
    the amount to increment the current volume by  
";

%feature("docstring") Material::setVolume "
setVolume(double volume)  

Set the volume/area of the Material.  

Parameters
----------
* volume :  
    the volume/area of the Material  
";

%feature("docstring") Material::getNumVectorGroups "
getNumVectorGroups() -> int  

Returns the rounded up number of energy groups to fill an integral number of vector
lengths.  

Returns
-------
The number of vector-aligned energy groups  
";

%feature("docstring") Material::getSigmaT "
getSigmaT() -> FP_PRECISION *  

Return the array of the Material's total cross-sections.  

Returns
-------
the pointer to the Material's array of total cross-sections  
";

%feature("docstring") Material::getSigmaS "
getSigmaS() -> FP_PRECISION *  

Return the array of the Material's scattering cross-section matrix.  

Returns
-------
the pointer to the Material's array of scattering cross-sections  
";

%feature("docstring") Material::transposeProductionMatrices "
transposeProductionMatrices()  

Transposes the scattering and fission matrices.  

This routine is used by the Solver when performing adjoint flux caclulations.  
";

%feature("docstring") Material::setNuSigmaFByGroup "
setNuSigmaFByGroup(double xs, int group)  

Set the Material's fission cross-section multiplied by $ \\nu $ for some energy group.  

This method is a helper function to allow OpenMOC users to assign the Material's nuclear
data in Python. A user must initialize a NumPy array of the correct size (e.g., a float64
array the length of the number of energy groups) as input to this function. This function
then fills the NumPy array with the data values for the Material's nu*fission cross-
sections. An example of how this function might be called in Python is as follows:  


Parameters
----------
* xs :  
    the fission cross-section multiplied by nu $ \\nu $  
* group :  
    the energy group  
";

%feature("docstring") Material::setSigmaFByGroup "
setSigmaFByGroup(double xs, int group)  

Set the Material's fission cross-section for some energy group.  

Parameters
----------
* xs :  
    the fission cross-section  
* group :  
    the energy group  
";

%feature("docstring") Material::setChi "
setChi(double *xs, int num_groups)  

Set the Material's array of chi $ \\chi $ values.  

This method is a helper function to allow OpenMOC users to assign the Material's nuclear
data in Python. A user must initialize a NumPy array of the correct size (e.g., a float64
array the length of the number of energy groups) as input to this function. This function
then fills the NumPy array with the data values for the Material's chi distribution. An
example of how this function might be called in Python is as follows:  


Parameters
----------
* xs :  
    the array of chi $ \\chi $ values  
* num_groups :  
    the number of energy groups  
";

%feature("docstring") Material::getSigmaF "
getSigmaF() -> FP_PRECISION *  

Return the array of the Material's fission cross-sections.  

Returns
-------
the pointer to the Material's array of fission cross-sections  
";

%feature("docstring") Material::getSigmaTByGroup "
getSigmaTByGroup(int group) -> FP_PRECISION  

Get the Material's total cross section for some energy group.  

Parameters
----------
* group :  
    the energy group  

Returns
-------
the total cross section  
";

%feature("docstring") Material::getNuSigmaF "
getNuSigmaF() -> FP_PRECISION *  

Return the array of the Material's fission cross-sections multiplied by nu $ \\nu $.  

Returns
-------
the pointer to the Material's array of fission cross-sections multiplied by nu $ \\nu $  
";

// File: classMatrix.xml


%feature("docstring") Matrix "
";

%feature("docstring") Matrix::getNumRows "
getNumRows() -> int  

Get the number of rows in the matrix.  

Returns
-------
The number of rows in the matrix.  
";

%feature("docstring") Matrix::getIA "
getIA() -> int *  

Get an array of the row indices (I) component of the CSR form of the full matrix (A).  

Returns
-------
A pointer to the I component of the CSR form of the full matrix (A).  
";

%feature("docstring") Matrix::getNNZLU "
getNNZLU() -> int  

Get the number of non-zero values in the lower + upper components of the matrix.  

Returns
-------
The number of non-zero values in the lower + upper components of the matrix.  
";

%feature("docstring") Matrix::clear "
clear()  

Clear all values in the matrix list of lists.  
";

%feature("docstring") Matrix::~Matrix "
~Matrix()  

Destructor clears list of lists and deletes the arrays used to represent the matrix in CSR
form.  
";

%feature("docstring") Matrix::getNumX "
getNumX() -> int  

Get the number of cells in the x dimension.  

Returns
-------
The number of cells in the x dimension.  
";

%feature("docstring") Matrix::getNumY "
getNumY() -> int  

Get the number of cells in the y dimension.  

Returns
-------
The number of cells in the y dimension.  
";

%feature("docstring") Matrix::getCellLocks "
getCellLocks() -> omp_lock_t *  

Return the array of cell locks for atomic cell operations.  

Returns
-------
an array of cell locks  
";

%feature("docstring") Matrix::Matrix "
Matrix(omp_lock_t *cell_locks, int num_x=1, int num_y=1, int num_groups=1)  

Constructor initializes Matrix as a list of lists and sets the matrix dimensions.  The
matrix object uses a \"lists of lists\" structure (implemented as a map of lists) to allow
for easy setting and incrementing of the values in the object. When the matrix is needed
to perform linear algebra operations, it is converted to compressed row storage (CSR) form
[1]. The matrix is ordered by cell (as opposed to by group) on the outside. Locks are used
to make the matrix thread-safe against concurrent writes the same value. One lock locks
out multiple rows of the matrix at a time reprsenting multiple groups in the same cell.  

[1] \"Sparse matrix\", Wikipedia, https://en.wikipedia.org/wiki/Sparse_matrix.  

Parameters
----------
* cell_locks :  
    Omp locks for atomic cell operations  
* num_x :  
    The number of cells in the x direction.  
* num_y :  
    The number of cells in the y direction.  
* num_groups :  
    The number of energy groups in each cell.  
";

%feature("docstring") Matrix::getDiag "
getDiag() -> FP_PRECISION *  

Get the diagonal component of the matrix object.  

Returns
-------
A pointer to the diagonal component of the matrix object.  
";

%feature("docstring") Matrix::transpose "
transpose()  

Transpose the matrix in place.  
";

%feature("docstring") Matrix::getNumGroups "
getNumGroups() -> int  

Get the number of groups in each cell.  

Returns
-------
The number of groups in each cell.  
";

%feature("docstring") Matrix::printString "
printString()  

Print the matrix object to the log file.  
";

%feature("docstring") Matrix::getILU "
getILU() -> int *  

Get an array of the row indices (I) component of the CSR form of the lower + upper (LU)
components of the matrix.  

Returns
-------
A pointer to the I component of the CSR form of the LU components of the matrix.  
";

%feature("docstring") Matrix::getA "
getA() -> FP_PRECISION *  

Get the full matrix (A) component of the CSR form of the matrix object.  

Returns
-------
A pointer to the A component of the CSR form matrix object.  
";

%feature("docstring") Matrix::getLU "
getLU() -> FP_PRECISION *  

Get the lower + upper (LU) component of the CSR form of the matrix object.  

Returns
-------
A pointer to the lower + upper (LU) component of the CSR form matrix object.  
";

%feature("docstring") Matrix::getJA "
getJA() -> int *  

Get an array of the column indices (J) component of the CSR form of the full matrix (A).  

Returns
-------
A pointer to the J component of the CSR form of the full matrix (A).  
";

%feature("docstring") Matrix::getNNZ "
getNNZ() -> int  

Get the number of non-zero values in the full matrix.  

Returns
-------
The number of non-zero values in the full matrix.  
";

%feature("docstring") Matrix::incrementValue "
incrementValue(int cell_from, int group_from, int cell_to, int group_to, FP_PRECISION val)  

Increment a value in the matrix.  This method takes a cell and group of origin (cell/group
from) and cell and group of destination (cell/group to) and floating point value. The
origin and destination are used to compute the row and column in the matrix. If a value
exists for the row/column, the value is incremented by val; otherwise, it is set to val.  

Parameters
----------
* cell_from :  
    The origin cell.  
* group_from :  
    The origin group.  
* cell_to :  
    The destination cell.  
* group_from :  
    The destination group.  
* val :  
    The value used to increment the row/column location.  
";

%feature("docstring") Matrix::getValue "
getValue(int cell_from, int group_from, int cell_to, int group_to) -> FP_PRECISION  

Get a value in the matrix.  This method takes a cell and group of origin (cell/group from)
and cell and group of destination (cell/group to). The origin and destination are used to
compute the row and column in the matrix. The value at the location specified by the
row/column is returned.  

Parameters
----------
* cell_from :  
    The origin cell.  
* group_from :  
    The origin group.  
* cell_to :  
    The destination cell.  
* group_from :  
    The destination group.  

Returns
-------
The value at the corresponding row/column location.  
";

%feature("docstring") Matrix::setValue "
setValue(int cell_from, int group_from, int cell_to, int group_to, FP_PRECISION val)  

Set a value in the matrix.  This method takes a cell and group of origin (cell/group from)
and cell and group of destination (cell/group to) and floating point value. The origin and
destination are used to compute the row and column in the matrix. The location specified
by the row/column is set to val.  

Parameters
----------
* cell_from :  
    The origin cell.  
* group_from :  
    The origin group.  
* cell_to :  
    The destination cell.  
* group_from :  
    The destination group.  
* val :  
    The value used to set the row/column location.  
";

%feature("docstring") Matrix::getJLU "
getJLU() -> int *  

Get an array of the column indices (J) component of the CSR form of the lower + upper (LU)
components of the matrix.  

Returns
-------
A pointer to the J component of the CSR form of the LU components of the matrix.  
";

// File: structmultiplyByConstant.xml


%feature("docstring") multiplyByConstant "

A functor to multiply all elements in a Thrust vector by a constant.  

Parameters
----------
* constant :  
    the constant to multiply the vector  

Attributes
----------
* constant : const T  
";

%feature("docstring") multiplyByConstant::multiplyByConstant "
multiplyByConstant(T constant)  

Constructor for the functor.  

Parameters
----------
* constant :  
    to multiply each element in a Thrust vector  
";

// File: structFixedHashMap_1_1node.xml

// File: structParallelHashMap_1_1paddedPointer.xml

// File: classParallelHashMap.xml


%feature("docstring") ParallelHashMap "

A thread-safe hash map supporting insertion and lookup operations.  

The ParallelHashMap class is built on top of the FixedHashMap class, supporting insertion
and lookup operations but not deletion as deletion is not needed in the OpenMOC
application. This hash table uses chaining for collisions, as defined in FixedHashMap. It
offers lock free lookups in O(1) time on average and fine-grained locking for insertions
in O(1) time on average as well. Resizing is conducted periodically during inserts,
although the starting table size can be chosen to limit the number of resizing operations.  

C++ includes: src/ParallelHashMap.h
";

%feature("docstring") ParallelHashMap::values "
values() -> V *  

Returns an array of the values in the underlying table.  

All buckets are scanned in order to form a list of all values present in the table and
then the list is returned. Threads announce their presence to ensure table memory is not
freed during access. WARNING: The user is responsible for freeing the allocated memory
once the array is no longer needed.  

Returns
-------
an array of values in the map whose length is the number of key/value pairs in the table.  
";

%feature("docstring") ParallelHashMap::insert_and_get_count "
insert_and_get_count(K key, V value) -> int  

Insert a given key/value pair into the parallel hash map and return the order number.  

First the underlying table is checked to determine if a resize should be conducted. Then,
the table is checked to see if it already contains the key. If so, the key/value pair is
not inserted and the function returns. Otherwise, the lock of the associated bucket is
acquired and the key/value pair is added to the bucket.  

Parameters
----------
* key :  
    key of the key/value pair to be inserted  
* value :  
    value of the key/value pair to be inserted  

Returns
-------
order number in which the key/value pair was inserted, -1 if it already exists  
";

%feature("docstring") ParallelHashMap::size "
size() -> size_t  

Returns the number of key/value pairs in the underlying table.  

Returns
-------
number of key/value pairs in the map  
";

%feature("docstring") ParallelHashMap::ParallelHashMap "
ParallelHashMap(size_t M=64, size_t L=64)  

Constructor generates initial underlying table as a fixed-sized hash map and intializes
concurrency structures.  
";

%feature("docstring") ParallelHashMap::print_buckets "
print_buckets()  

Prints the contents of each bucket to the screen.  

All buckets are scanned and the contents of the buckets are printed, which are pointers to
linked lists. If the pointer is NULL suggesting that the linked list is empty, NULL is
printed to the screen. Threads announce their presence to ensure table memory is not freed
during access.  
";

%feature("docstring") ParallelHashMap::update "
update(K key, V value)  

Updates the value associated with a key in the parallel hash map.  

The thread first acquires the lock for the bucket associated with the key is acquired,
then the linked list in the bucket is searched for the key. If the key is not found, an
exception is returned. When the key is found, the value is updated and the lock is
released.  

Parameters
----------
* key :  
    the key of the key/value pair to be updated  
* value :  
    the new value for the key/value pair  
";

%feature("docstring") ParallelHashMap::keys "
keys() -> K *  

Returns an array of the keys in the underlying table.  

All buckets are scanned in order to form a list of all keys present in the table and then
the list is returned. Threads announce their presence to ensure table memory is not freed
during access. WARNING: The user is responsible for freeing the allocated memory once the
array is no longer needed.  

Returns
-------
an array of keys in the map whose length is the number of key/value pairs in the table.  
";

%feature("docstring") ParallelHashMap::insert "
insert(K key, V value)  

Insert a given key/value pair into the parallel hash map.  

First the underlying table is checked to determine if a resize should be conducted. Then,
the table is checked to see if it already contains the key. If so, the key/value pair is
not inserted and the function returns. Otherwise, the lock of the associated bucket is
acquired and the key/value pair is added to the bucket.  

Parameters
----------
* key :  
    key of the key/value pair to be inserted  
* value :  
    value of the key/value pair to be inserted  
";

%feature("docstring") ParallelHashMap::clear "
clear()  

Clears all key/value pairs form the hash table.  
";

%feature("docstring") ParallelHashMap::~ParallelHashMap "
~ParallelHashMap()  

Destructor frees memory associated with fixed-sized hash map and concurrency structures.  
";

%feature("docstring") ParallelHashMap::num_locks "
num_locks() -> size_t  

Returns the number of locks in the parallel hash map.  

Returns
-------
number of locks in the map  
";

%feature("docstring") ParallelHashMap::contains "
contains(K key) -> bool  

Determine whether the parallel hash map contains a given key.  

First the thread accessing the table announces its presence and which table it is reading.
Then the linked list in the bucket associated with the key is searched without setting any
locks to determine whether the key is present. When the thread has finished accessing the
table, the announcement is reset to NULL. The announcement ensures that the data in the
map is not freed during a resize until all threads have finished accessing the map.  

Parameters
----------
* key :  
    key to be searched  

Returns
-------
boolean value referring to whether the key is contained in the map  
";

%feature("docstring") ParallelHashMap::bucket_count "
bucket_count() -> size_t  

Returns the number of buckets in the underlying table.  

Returns
-------
number of buckets in the map  
";

%feature("docstring") ParallelHashMap::at "
at(K key) -> V  

Determine the value associated with a given key.  

This function follows the same algorithm as <contains> except that the value associated
with the searched key is returned. First the thread accessing the table acquires the lock
corresponding with the associated bucket based on the key. Then the linked list in the
bucket is searched for the key. An exception is thrown if the key is not found. When the
thread has finished accessing the table, it releases the lock.  

Parameters
----------
* key :  
    key to be searched  

Returns
-------
value associated with the key  
";

// File: classPlane.xml


%feature("docstring") Plane "

Represents a Plane perpendicular to the xy-plane.  

C++ includes: src/Surface.h
";

%feature("docstring") Plane::getMaxZ "
getMaxZ(int halfspace) -> double  

Returns the maximum z value of INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum z value of INFINITY  
";

%feature("docstring") Plane::getMaxX "
getMaxX(int halfspace) -> double  

Returns the maximum x value of INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum x value of INFINITY  
";

%feature("docstring") Plane::getMaxY "
getMaxY(int halfspace) -> double  

Returns the maximum y value of INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum y value of INFINITY  
";

%feature("docstring") Plane::toString "
toString() -> std::string  

Converts this Plane's attributes to a character array.  

The character array returned conatins the type of Plane (ie, PLANE) and the A, B, and C
coefficients in the quadratic Surface equation.  

Returns
-------
a character array of this Plane's attributes  
";

%feature("docstring") Plane::getA "
getA() -> double  

Returns the A coefficient multiplying x in the surface equation.  

Returns
-------
the value for the A coefficient  
";

%feature("docstring") Plane::getC "
getC() -> double  

Returns the C coefficient multiplying z in the surface equation.  

Returns
-------
the value for the C coefficient  
";

%feature("docstring") Plane::getB "
getB() -> double  

Returns the B coefficient multiplying y in the surface equation.  

Returns
-------
the value for the B coefficient  
";

%feature("docstring") Plane::getD "
getD() -> double  

Returns the D constant coefficient.  

Returns
-------
the value for the D coefficient  
";

%feature("docstring") Plane::getMinX "
getMinX(int halfspace) -> double  

Returns the minimum x value of -INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum x value of -INFINITY  
";

%feature("docstring") Plane::getMinY "
getMinY(int halfspace) -> double  

Returns the minimum y value of -INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum y value of -INFINITY  
";

%feature("docstring") Plane::getMinZ "
getMinZ(int halfspace) -> double  

Returns the minimum z value of -INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum z value of -INFINITY  
";

%feature("docstring") Plane::Plane "
Plane(const double A, const double B, const double C, const double D, const int id=0,
    const char *name=\"\")  

Constructor.  

Parameters
----------
* A :  
    the first coefficient in $ A * x + B * y + C * z + D = 0 $  
* B :  
    the second coefficient in $ A * x + B * y + C * z + D = 0 $  
* C :  
    the third coefficient in $ A * x + B * y + C * z + D = 0 $  
* D :  
    the fourth coefficient in $ A * x + B * y + C * z + D = 0 $  
* id :  
    the optional Surface ID  
* name :  
    the optional name of the Surface  
";

%feature("docstring") Plane::evaluate "
evaluate(const Point *point) const  -> double  

Evaluate a Point using the Plane's quadratic Surface equation.  

Parameters
----------
* point :  
    a pointer to the Point of interest  

Returns
-------
the value of Point in the Plane's quadratic equation  
";

%feature("docstring") Plane::intersection "
intersection(Point *point, double angle, Point *points) -> int  

Finds the intersection Point with this Plane from a given Point and trajectory defined by
an angle.  

Parameters
----------
* point :  
    pointer to the Point of interest  
* angle :  
    the angle defining the trajectory in radians  
* points :  
    pointer to a Point to store the intersection Point  

Returns
-------
the number of intersection Points (0 or 1)  
";

// File: classPoint.xml


%feature("docstring") Point "

Class to represent a 3D point in space.  

C++ includes: src/Point.h
";

%feature("docstring") Point::toString "
toString() -> std::string  

Converts this Point to a character representation of its attributes.  

The character array includes the x-coordinate, y-coordinate, and z-coordinate  

Returns
-------
a character array of this Point's attributes  
";

%feature("docstring") Point::getY "
getY() const  -> double  

Returns this Point's y-coordinate.  

Returns
-------
the y-coordinate  
";

%feature("docstring") Point::getX "
getX() const  -> double  

Returns this Point's x-coordinate.  

Returns
-------
the x-coordinate  
";

%feature("docstring") Point::getZ "
getZ() const  -> double  

Returns this Point's z-coordinate.  

Returns
-------
the z-coordinate  
";

%feature("docstring") Point::Point "
Point()  

Constructor initializes an empty Point.  
";

%feature("docstring") Point::setX "
setX(const double x)  

Set the Point's x-coordinate.  

Parameters
----------
* x :  
    the new x-coordinate  
";

%feature("docstring") Point::setZ "
setZ(const double z)  

Set the Point's z-coordinate.  

Parameters
----------
* z :  
    the new z-coordinate  
";

%feature("docstring") Point::setY "
setY(const double y)  

Set the Point's y-coordinate.  

Parameters
----------
* y :  
    the new y-coordinate  
";

%feature("docstring") Point::distanceToPoint "
distanceToPoint(const Point *point) -> double  

Compute the distance from this Point to another Point of interest.  

Parameters
----------
* point :  
    a pointer to the Point of interest  

Returns
-------
distance to the Point of interest  
";

%feature("docstring") Point::~Point "
~Point()  

Destructor.  
";

%feature("docstring") Point::setCoords "
setCoords(const double x, const double y, const double z)  

Initializes a Point with two-dimensional coordinates.  

Parameters
----------
* x :  
    x-coordinate  
* y :  
    y-coordinate  
";

// File: classQuadrature.xml


%feature("docstring") Quadrature "

The arbitrary quadrature parent class.  

C++ includes: src/Quadrature.h
";

%feature("docstring") Quadrature::getPolarWeights "
getPolarWeights() -> FP_PRECISION **  

Returns a pointer to the Quadrature's array of polar weights.  

Returns
-------
a pointer to the polar weights array  
";

%feature("docstring") Quadrature::setPolarWeights "
setPolarWeights(FP_PRECISION *weights, int num_azim_times_polar)  

Set the Quadrature's array of polar weights.  

This method is a helper function to allow OpenMOC users to assign the Quadrature's polar
weights in Python. A user must initialize a NumPy array of the correct size (e.g., a
float64 array the length of the number of azimuthal times polar angles) as input to this
function. This function then fills the Quadrature's polar weights with the given values.
An example of how this function might be called in Python is as follows:  


Parameters
----------
* weights :  
    The polar weights  
* num_azim_times_polar :  
    the total number of angles in one octant (azimuthal x polar)  
";

%feature("docstring") Quadrature::setPolarWeight "
setPolarWeight(FP_PRECISION weight, int azim, int polar)  

Sets the polar weight for the given indexes.  

Parameters
----------
* weight :  
    the weight of the polar angle  
* azim :  
    the azimuthal index corresponding to the angle  
* azim :  
    the polar index corresponding to the angle  
";

%feature("docstring") Quadrature::getAzimWeights "
getAzimWeights() -> FP_PRECISION *  

Returns a pointer to the Quadrature's array of azimuthal weights.  

Returns
-------
a pointer to the azimuthal weights array  
";

%feature("docstring") Quadrature::~Quadrature "
~Quadrature()  

Destructor deletes arrray of sines of the polar angles, the weights of the polar angles
and the products of the sines and weights.  
";

%feature("docstring") Quadrature::getThetas "
getThetas() -> double **  

Returns a pointer to the Quadrature's array of polar angles $ \\theta_{p} $.  

Returns
-------
a pointer to the array of $ \\theta_{p} $  
";

%feature("docstring") Quadrature::getNumAzimAngles "
getNumAzimAngles() const  -> int  

Returns the number of azimuthal angles.  

Returns
-------
the number of azimuthal angles  
";

%feature("docstring") Quadrature::precomputeWeights "
precomputeWeights(bool solve_3D)  

This private routine computes the product of the sine thetas and weights for each angle in
the polar quadrature.  

Note that this routine must be called after populating the sine thetas and weights arrays.  
";

%feature("docstring") Quadrature::setThetas "
setThetas(double *thetas, int num_azim_times_polar)  

Sets the Quadrature's array of polar angles.  

This method is a helper function to allow OpenMOC users to assign the Quadrature's polar
angles in Python. A user must initialize a NumPy array of the correct size (e.g., a
float64 array the length of the number of azimuthal times polar angles) as input to this
function. This function then fills the Quadrature's polar angles with the given values. An
example of how this function might be called in Python is as follows:  


Parameters
----------
* thetas :  
    the array of polar angle for each azimuthal/polar angle combination  
* num_azim_times_polar :  
    the total number of angles (azimuthal x polar)  
";

%feature("docstring") Quadrature::getNumPolarAngles "
getNumPolarAngles() const  -> int  

Returns the number of polar angles.  

Returns
-------
the number of polar angles  
";

%feature("docstring") Quadrature::getAzimWeight "
getAzimWeight(int azim) -> FP_PRECISION  

Returns the azimuthal angle weight value for a particular azimuthal angle.  

Parameters
----------
* azim :  
    index of the azimuthal angle of interest  

Returns
-------
the weight for an azimuthal angle  
";

%feature("docstring") Quadrature::initialize "
initialize()  

Initialize the polar quadrature azimuthal angles.  

The parent class routine simply checks that number of polar and azimuthal angles have been
set by the user and generates the azimuthal angles if not already generated.  
";

%feature("docstring") Quadrature::getSinThetas "
getSinThetas() -> FP_PRECISION **  

Returns a pointer to the Quadrature's array of polar angle sines $ sin\\theta_{p} $.  

Returns
-------
a pointer to the array of $ sin\\theta_{p} $  
";

%feature("docstring") Quadrature::getAzimSpacings "
getAzimSpacings() -> FP_PRECISION *  

Returns an array of adjusted azimuthal spacings.  

An array of azimuthal spacings after adjustment is returned, indexed by azimuthal angle  

Returns
-------
the array of azimuthal spacings  
";

%feature("docstring") Quadrature::getPolarWeight "
getPolarWeight(int azim, int polar) -> FP_PRECISION  

Returns the polar weight for a particular azimuthal and polar angle.  

Parameters
----------
* azim :  
    index of the azimthal angle of interest  
* polar :  
    index of the polar angle of interest  

Returns
-------
the value of the polar weight for this azimuthal and polar angle  
";

%feature("docstring") Quadrature::setPhi "
setPhi(double phi, int azim)  

Sets the azimuthal angle for the given index.  

Parameters
----------
* phi :  
    the value in radians of the azimuthal angle to be set  
* azim :  
    the azimuthal index  
";

%feature("docstring") Quadrature::getPhi "
getPhi(int azim) -> double  

Returns the azimuthal angle value in radians.  

Parameters
----------
* azim :  
    index of the azimthal angle of interest  

Returns
-------
the value of the azimuthal angle  
";

%feature("docstring") Quadrature::getWeight "
getWeight(int azim, int polar) -> FP_PRECISION  

Returns the total weight for Tracks with the given azimuthal and polar indexes.  

Angular weights are multiplied by Track spcings  

Parameters
----------
* azim :  
    index of the azimuthal angle of interest  
* polar :  
    index of the polar angle of interest  

Returns
-------
the total weight of each Track with the given indexes  
";

%feature("docstring") Quadrature::setAzimSpacing "
setAzimSpacing(FP_PRECISION spacing, int azim)  

Sets the azimuthal spacing for the given index.  

Parameters
----------
* spacing :  
    the spacing (cm) in the azimuthal direction to be set  
* azim :  
    the azimuthal index  
";

%feature("docstring") Quadrature::getSinTheta "
getSinTheta(int azim, int polar) -> FP_PRECISION  

Returns the $ sin(\\theta)$ value for a particular polar angle.  

Parameters
----------
* azim :  
    index of the azimthal angle of interest  
* polar :  
    index of the polar angle of interest  

Returns
-------
the value of $ \\sin(\\theta) $ for this azimuthal and polar angle  
";

%feature("docstring") Quadrature::getPhis "
getPhis() -> double *  

Returns a pointer to the Quadrature's array of azimuthal angles $ \\phi $.  

Returns
-------
a pointer to the array of $ \\phi $  
";

%feature("docstring") Quadrature::setTheta "
setTheta(double theta, int azim, int polar)  

Sets the polar angle for the given indexes.  

Parameters
----------
* theta :  
    the value in radians of the polar angle to be set  
* azim :  
    the azimuthal index of the angle of interest  
* polar :  
    the polar index of the angle of interest  
";

%feature("docstring") Quadrature::setNumAzimAngles "
setNumAzimAngles(const int num_azim)  

Set the number of azimuthal angles to initialize.  

Parameters
----------
* num_azim :  
    the number of azimuthal angles  
";

%feature("docstring") Quadrature::toString "
toString() -> std::string  

Converts this Quadrature to a character array of its attributes.  

The character array includes the number of polar angles, the the values of the sine and
weight of each polar angle, and the product of the sine and weight of each polar angle.  

Returns
-------
a character array of the Quadrature's attributes  
";

%feature("docstring") Quadrature::getAzimSpacing "
getAzimSpacing(int azim) -> FP_PRECISION  

Returns the adjusted azimuthal spacing at the requested azimuthal angle index.  

The aziumthal spacing depends on the azimuthal angle. This function returns the azimuthal
spacing used at the desired azimuthal angle index.  

Parameters
----------
* azim :  
    the requested azimuthal angle index  

Returns
-------
the requested azimuthal spacing  
";

%feature("docstring") Quadrature::setAzimWeight "
setAzimWeight(double weight, int azim)  

Sets the azimuthal weight for the given index.  

Parameters
----------
* weight :  
    the weight of the azimuthal angle  
* azim :  
    the azimuthal index  
";

%feature("docstring") Quadrature::Quadrature "
Quadrature()  

Dummy constructor sets the default number of angles to zero.  
";

%feature("docstring") Quadrature::getQuadratureType "
getQuadratureType() -> quadratureType  

Returns the type of Quadrature created.  

Returns
-------
The quadrature type  
";

%feature("docstring") Quadrature::getWeightInline "
getWeightInline(int azim, int polar) -> FP_PRECISION  

Returns the total weight for Tracks with the given azimuthal and polar indexes without
error checking and inlined.  

Angular weights are multiplied by Track spcings  

Parameters
----------
* azim :  
    index of the azimuthal angle of interest  
* polar :  
    index of the polar angle of interest  

Returns
-------
the total weight of each Track with the given indexes  
";

%feature("docstring") Quadrature::getTheta "
getTheta(int azim, int polar) -> double  

Returns the polar angle in radians for a given azimuthal and polar angle index.  

Parameters
----------
* azim :  
    index of the azimthal angle of interest  
* polar :  
    index of the polar angle of interest  

Returns
-------
the value of the polar angle for this azimuthal and polar angle index  
";

%feature("docstring") Quadrature::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles  
";

// File: structsecond__t.xml


%feature("docstring") second_t "

A helper struct for the Universe::findCell() method.  

This is used to insert a Universe's Cells to the back of a vector of neighbor Cells in
Universe::findCell() routine. This works in symbiosis with the pair_second method template
defined below.  

C++ includes: Universe.h
";

// File: structsegment.xml


%feature("docstring") segment "

A segment represents a line segment within a single flat source region along a track.  

Attributes
----------
* _length : FP_PRECISION  
    The length of the segment (cm)  

* _material : Material *  
    A pointer to the material in which this segment resides  

* _region_id : int  
    The ID for flat source region in which this segment resides  

* _cmfd_surface_fwd : int  
    The ID for the mesh surface crossed by the Track end point  

* _cmfd_surface_bwd : int  
    The ID for the mesh surface crossed by the Track start point  

C++ includes: Track.h
";

%feature("docstring") segment::segment "
segment()  

Constructor initializes CMFD surfaces  
";

// File: classSolver.xml


%feature("docstring") Solver "

This is an abstract base class which different Solver subclasses implement for different
architectures or source iteration algorithms.  

C++ includes: src/Solver.h
";

%feature("docstring") Solver::setConvergenceThreshold "
setConvergenceThreshold(FP_PRECISION threshold)  

Sets the threshold for source/flux convergence.  

The default threshold for convergence is 1E-5.  

Parameters
----------
* source_thresh :  
    the threshold for source/flux convergence  
";

%feature("docstring") Solver::computeFlux "
computeFlux(int max_iters=1000, solverMode mode=FORWARD, bool only_fixed_source=true)  

Computes the scalar flux distribution by performing a series of transport sweeps.  

This is the main method exposed to the user through the Python interface to compute the
scalar flux distribution, e.g., for a fixed source calculation. This routine makes an
initial guess for scalar and boundary fluxes and performs transport sweep until
convergence.  

By default, this method will perform a maximum of 1000 transport sweeps with a 1E-5
threshold on the average FSR scalar flux. These values may be freely modified by the user
at runtime.  

The only_fixed_source runtime parameter may be used to control the type of source
distribution used in the calculation. By default, this paramter is true and only the fixed
sources specified by the user will be considered. Alternatively, when the parameter is
false, the source will be computed as the scattering and fission sources resulting from a
previously computed flux distribution (e.g., an eigenvalue calculation) in addition to any
user-defined fixed sources.  

This method may be called by the user to compute the scalar flux for a fixed source
distribution from Python as follows:  

  

         Alternatively, as described above, this method may be called by
         the user in Python to compute the flux from a superposition of
         fixed and / or eigenvalue sources as follows:  


Parameters
----------
* max_iters :  
    the maximum number of source iterations to allow  
* mode :  
    the solution type (FORWARD or ADJOINT)  
* only_fixed_source :  
    use only fixed sources (true by default)  
";

%feature("docstring") Solver::countFissionableFSRs "
countFissionableFSRs()  

Counts the number of fissionable flat source regions.  

This routine is used by the Solver::computeEigenvalue(...) routine which uses the number
of fissionable FSRs to normalize the residual on the fission source distribution.  
";

%feature("docstring") Solver::setExpPrecision "
setExpPrecision(FP_PRECISION precision)  

Set the precision, or maximum allowable approximation error, of the the exponential
interpolation table.  

By default, the precision is 1E-5 based on the analysis in Yamamoto's 2003 paper.  

Parameters
----------
* precision :  
    the precision of the exponential interpolation table,  
";

%feature("docstring") Solver::getTotalTime "
getTotalTime() -> double  

Returns the total time to converge the source (seconds).  

Returns
-------
the time to converge the source (seconds)  
";

%feature("docstring") Solver::getTrackGenerator "
getTrackGenerator() -> TrackGenerator *  

Returns a pointer to the TrackGenerator.  

Returns
-------
a pointer to the TrackGenerator  
";

%feature("docstring") Solver::setFixedSourceByCell "
setFixedSourceByCell(Cell *cell, int group, FP_PRECISION source)  

Assign a fixed source for a Cell and energy group.  

Parameters
----------
* cell :  
    the Cell of interest  
* group :  
    the energy group  
* source :  
    the volume-averaged source in this group  
";

%feature("docstring") Solver::computeFSRFissionSources "
computeFSRFissionSources()=0  

Computes the total fission source for each FSR and energy group.  
";

%feature("docstring") Solver::normalizeFluxes "
normalizeFluxes()=0  

Normalizes all FSR scalar fluxes and Track boundary angular fluxes to the total fission
source (times $ \\nu $).  
";

%feature("docstring") Solver::getFlux "
getFlux(int fsr_id, int group) -> FP_PRECISION  

Returns the scalar flux for some FSR and energy group.  

Parameters
----------
* fsr_id :  
    the ID for the FSR of interest  
* group :  
    the energy group of interest  

Returns
-------
the FSR scalar flux  
";

%feature("docstring") Solver::Solver "
Solver(TrackGenerator *track_generator=NULL)  

Constructor initializes an empty Solver class with array pointers set to NULL.  

Parameters
----------
* track_generator :  
    an optional pointer to a TrackGenerator object  
";

%feature("docstring") Solver::getMaxOpticalLength "
getMaxOpticalLength() -> FP_PRECISION  

Get the maximum allowable optical length for a track segment.  

Returns
-------
The max optical length  
";

%feature("docstring") Solver::computeEigenvalue "
computeEigenvalue(int max_iters=1000, solverMode mode=FORWARD, residualType
    res_type=FISSION_SOURCE)  

Computes keff by performing a series of transport sweep and source updates.  

This is the main method exposed to the user through the Python interface to perform an
eigenvalue calculation. The method makes an initial guess for the scalar and boundary
fluxes and performs transport sweeps and source updates until convergence.  

By default, this method will perform a maximum of 1000 transport sweeps with a 1E-5
threshold on the integrated FSR fission source. These values may be freely modified by the
user at runtime.  

The res_type parameter may be used to control the convergence criterion - SCALAR_FLUX,
TOTAL_SOURCE and FISSION_SOURCE (default) are all supported options in OpenMOC at this
time.  


Parameters
----------
* max_iters :  
    the maximum number of source iterations to allow  
* mode :  
    the solution type (FORWARD or ADJOINT)  
* res_type :  
    the type of residual used for the convergence criterion  
";

%feature("docstring") Solver::isUsingExponentialInterpolation "
isUsingExponentialInterpolation() -> bool  

Returns whether the Solver uses linear interpolation to compute exponentials.  

Returns
-------
true if using linear interpolation to compute exponentials  
";

%feature("docstring") Solver::setFluxes "
setFluxes(FP_PRECISION *in_fluxes, int num_fluxes)=0  
";

%feature("docstring") Solver::getFluxes "
getFluxes(FP_PRECISION *out_fluxes, int num_fluxes)=0  
";

%feature("docstring") Solver::getGeometry "
getGeometry() -> Geometry *  

Returns a pointer to the Geometry.  

Returns
-------
a pointer to the Geometry  
";

%feature("docstring") Solver::initializeSourceArrays "
initializeSourceArrays()=0  

Allocates memory for FSR source arrays.  
";

%feature("docstring") Solver::getNumPolarAngles "
getNumPolarAngles() -> int  

Returns the number of angles used for the polar quadrature.  

Returns
-------
the number of polar angles  
";

%feature("docstring") Solver::setFixedSourceByFSR "
setFixedSourceByFSR(int fsr_id, int group, FP_PRECISION source)  

Assign a fixed source for a flat source region and energy group.  

Parameters
----------
* fsr_id :  
    the flat source region ID  
* group :  
    the energy group  
* source :  
    the volume-averaged source in this group  
";

%feature("docstring") Solver::computeSource "
computeSource(int max_iters=1000, solverMode mode=FORWARD, double k_eff=1.0, residualType
    res_type=TOTAL_SOURCE)  

Computes the total source distribution by performing a series of transport sweep and
source updates.  

This is the main method exposed to the user through the Python interface to compute the
source distribution, e.g., for a fixed and/or external source calculation. This routine
makes an initial guess for the scalar and boundary fluxes and performs transport sweeps
and source updates until convergence.  

By default, this method will perform a maximum of 1000 transport sweeps with a 1E-5
threshold on the integrated FSR total source. These values may be freely modified by the
user at runtime.  

The k_eff parameter may be used for fixed source calculations with fissionable material
(e.g., start-up in a reactor from a fixed external source). In this case, the user must
\"guess\" the critical eigenvalue to be be used to scale the fission source.  

The res_type parameter may be used to control the convergence criterion - SCALAR_FLUX,
TOTAL_SOURCE (default) and FISSION_SOURCE are all supported options in OpenMOC at this
time.  

This method may be called by the user from Python as follows:  


Parameters
----------
* max_iters :  
    the maximum number of source iterations to allow  
* mode :  
    the solution type (FORWARD or ADJOINT)  
* k_eff :  
    the sub/super-critical eigenvalue (default 1.0)  
* res_type :  
    the type of residual used for the convergence criterion  
";

%feature("docstring") Solver::flattenFSRFluxes "
flattenFSRFluxes(FP_PRECISION value)=0  

Set the scalar flux for each FSR and energy group to some value.  

Parameters
----------
* value :  
    the value to assign to each FSR scalar flux  
";

%feature("docstring") Solver::getFSRVolume "
getFSRVolume(int fsr_id) -> FP_PRECISION  

Returns the calculated volume for a flat source region.  

Parameters
----------
* fsr_id :  
    the flat source region ID of interest  

Returns
-------
the flat source region volume  
";

%feature("docstring") Solver::initializeFixedSources "
initializeFixedSources()  

Assigns fixed sources assigned by Cell, Material to FSRs.  

Fixed sources assigned by Material  
";

%feature("docstring") Solver::getKeff "
getKeff() -> FP_PRECISION  

Returns the converged eigenvalue $ k_{eff} $.  

Returns
-------
the converged eigenvalue $ k_{eff} $  
";

%feature("docstring") Solver::useExponentialIntrinsic "
useExponentialIntrinsic()  

Informs the Solver to use the exponential intrinsic exp(...) function to compute the
exponential in the transport equation.  
";

%feature("docstring") Solver::computeResidual "
computeResidual(residualType res_type)=0 -> double  

Computes the residual between successive flux/source iterations.  

Parameters
----------
* res_type :  
    the residual type (SCALAR_FLUX, FISSION_SOURCE, TOTAL_SOURCE)  

Returns
-------
the total residual summed over FSRs and energy groups  
";

%feature("docstring") Solver::transportSweep "
transportSweep()=0  

This method performs one transport swep.  
";

%feature("docstring") Solver::initializeMaterials "
initializeMaterials(solverMode mode=FORWARD)  

Initializes the Material fission matrices.  

In an adjoint calculation, this routine will transpose the scattering and fission matrices
in each material.  

Parameters
----------
* mode :  
    the solution type (FORWARD or ADJOINT)  
";

%feature("docstring") Solver::isUsingDoublePrecision "
isUsingDoublePrecision() -> bool  

Returns whether the solver is using double floating point precision.  

Returns
-------
true if using double precision float point arithmetic  
";

%feature("docstring") Solver::getConvergenceThreshold "
getConvergenceThreshold() -> FP_PRECISION  

Returns the threshold for source/flux convergence.  

Returns
-------
the threshold for source/flux convergence  
";

%feature("docstring") Solver::printTimerReport "
printTimerReport()  

Prints a report of the timing statistics to the console.  
";

%feature("docstring") Solver::initializeCmfd "
initializeCmfd()  

Initializes a Cmfd object for acceleratiion prior to source iteration.  

Instantiates a dummy Cmfd object if one was not assigned to the Solver by the user and
initializes FSRs, materials, fluxes and the Mesh object. This method is for internal use
only and should not be called directly by the user.  
";

%feature("docstring") Solver::setFixedSourceByMaterial "
setFixedSourceByMaterial(Material *material, int group, FP_PRECISION source)  

Assign a fixed source for a Material and energy group.  

Parameters
----------
* material :  
    the Material of interest  
* group :  
    the energy group  
* source :  
    the volume-averaged source in this group  
";

%feature("docstring") Solver::setMaxOpticalLength "
setMaxOpticalLength(FP_PRECISION max_optical_length)  

Set the maximum allowable optical length for a track segment.  

Parameters
----------
* max_optical_length :  
    The max optical length  
";

%feature("docstring") Solver::computeKeff "
computeKeff()=0  

Compute $ k_{eff} $ from successive fission sources.  
";

%feature("docstring") Solver::getFSRSource "
getFSRSource(int fsr_id, int group) -> FP_PRECISION  

Returns the source for some energy group for a flat source region.  

This is a helper routine used by the openmoc.process module.  

Parameters
----------
* fsr_id :  
    the ID for the FSR of interest  
* group :  
    the energy group of interest  

Returns
-------
the flat source region source  
";

%feature("docstring") Solver::initializeFluxArrays "
initializeFluxArrays()=0  

Initializes Track boundary angular and FSR scalar flux arrays.  
";

%feature("docstring") Solver::scatterTransportSweep "
scatterTransportSweep()  

This method performs one transport sweep using the scatter source.  

This is a helper routine used for Krylov subspace methods.  
";

%feature("docstring") Solver::addSourceToScalarFlux "
addSourceToScalarFlux()=0  

Add the source term contribution in the transport equation to the FSR scalar flux.  
";

%feature("docstring") Solver::initializeFSRs "
initializeFSRs()  

Initializes the FSR volumes and Materials array.  

This method assigns each FSR a unique, monotonically increasing ID, sets the Material for
each FSR, and assigns a volume based on the cumulative length of all of the segments
inside the FSR.  
";

%feature("docstring") Solver::setTrackGenerator "
setTrackGenerator(TrackGenerator *track_generator)  

Sets the Solver's TrackGenerator with characteristic Tracks.  

The TrackGenerator must already have generated Tracks and have used ray tracing to
segmentize them across the Geometry. This should be initated in Python prior to assigning
the TrackGenerator to the Solver:  


Parameters
----------
* track_generator :  
    a pointer to a TrackGenerator object  
";

%feature("docstring") Solver::getNumIterations "
getNumIterations() -> int  

Returns the number of source iterations to converge the source.  

Returns
-------
the number of iterations  
";

%feature("docstring") Solver::storeFSRFluxes "
storeFSRFluxes()=0  

Stores the current scalar fluxes in the old scalar flux array.  
";

%feature("docstring") Solver::zeroTrackFluxes "
zeroTrackFluxes()=0  

Zero each Track's boundary fluxes for each energy group and polar angle in the \"forward\"
and \"reverse\" directions.  
";

%feature("docstring") Solver::fissionTransportSweep "
fissionTransportSweep()  

This method performs one transport sweep using the fission source.  

This is a helper routine used for Krylov subspace methods.  
";

%feature("docstring") Solver::~Solver "
~Solver()  

Destructor deletes arrays of boundary angular fluxes, scalar fluxes and sources for each
FSR and energy group.  

Deallocates memory for all arrays allocated for the Solver, including fluxes, sources,
quadrature weights, and exponential linear interpolation table.  
";

%feature("docstring") Solver::computeFSRSources "
computeFSRSources()=0  

Computes the total source (fission, scattering, fixed) for each FSR and energy group.  
";

%feature("docstring") Solver::setGeometry "
setGeometry(Geometry *geometry)  

Sets the Geometry for the Solver.  

This is a private setter method for the Solver and is not intended to be called by the
user.  

Parameters
----------
* geometry :  
    a pointer to a Geometry object  
";

%feature("docstring") Solver::resetMaterials "
resetMaterials(solverMode mode=FORWARD)  

Returns the Material data to its original state.  

In an adjoint calculation, the scattering and fission matrices in each material are
transposed during initialization. This routine returns both matrices to their original
(FORWARD) state at the end of a calculation.  

Parameters
----------
* mode :  
    the solution type (FORWARD or ADJOINT)  
";

%feature("docstring") Solver::useExponentialInterpolation "
useExponentialInterpolation()  

Informs the Solver to use linear interpolation to compute the exponential in the transport
equation.  
";

%feature("docstring") Solver::initializeExpEvaluator "
initializeExpEvaluator()  

Initializes new ExpEvaluator object to compute exponentials.  
";

%feature("docstring") Solver::computeFSRScatterSources "
computeFSRScatterSources()=0  

Computes the total scattering source for each FSR and energy group.  
";

%feature("docstring") Solver::computeFSRFissionRates "
computeFSRFissionRates(double *fission_rates, int num_FSRs)=0  

Computes the volume-weighted, energy integrated fission rate in each FSR and stores them
in an array indexed by FSR ID.  

This is a helper method for SWIG to allow users to retrieve FSR fission rates as a NumPy
array. An example of how this method can be called from Python is as follows:  


Parameters
----------
* fission_rates :  
    an array to store the fission rates (implicitly passed in as a NumPy array from
    Python)  
* num_FSRs :  
    the number of FSRs passed in from Python  
";

// File: structstrided__range_1_1stride__functor.xml


%feature("docstring") strided_range::stride_functor "

Attributes
----------
* stride : difference_type  
";

%feature("docstring") strided_range::stride_functor::stride_functor "
stride_functor(difference_type stride)  
";

// File: classstrided__range.xml


%feature("docstring") strided_range "
";

%feature("docstring") strided_range::end "
end(void) const  -> iterator  

Get the last element in the iterator.  

Returns
-------
the last element in the iterator  
";

%feature("docstring") strided_range::strided_range "
strided_range(Iterator first, Iterator last, difference_type stride)  

The strided iterator constructor.  
";

%feature("docstring") strided_range::begin "
begin(void) const  -> iterator  

Get the first element in the iterator.  

Returns
-------
the first element in the iterator  
";

// File: classSurface.xml


%feature("docstring") Surface "

Represents a general Surface in 3D.  

The Surface class and its subclasses are used to define the geometry for an OpenMOC
simulation using a constructive solid geometry (CSG) formalism. Surfaces are used during
ray tracing of charateristic tracks across the geometry.  

C++ includes: src/Surface.h
";

%feature("docstring") Surface::getUid "
getUid() const  -> int  

Return the Surface's unique ID.  

Returns
-------
the Surface's unique ID  
";

%feature("docstring") Surface::~Surface "
~Surface()  

Destructor.  
";

%feature("docstring") Surface::getMaxY "
getMaxY(int halfspace)=0 -> double  

Returns the maximum y value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum y value  
";

%feature("docstring") Surface::getMaxX "
getMaxX(int halfspace)=0 -> double  

Returns the maximum x value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum x value  
";

%feature("docstring") Surface::Surface "
Surface(const int id=0, const char *name=\"\")  

Constructor assigns unique ID and user-defined ID for a Surface.  

Assigns a default boundary condition for this Surface to BOUNDARY_NONE.  

Parameters
----------
* id :  
    an optional user-defined Surface ID  
* name :  
    an optional user-defined Surface name  
";

%feature("docstring") Surface::getMaxZ "
getMaxZ(int halfspace)=0 -> double  

Returns the maximum z value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the maximum z value  
";

%feature("docstring") Surface::getMinZ "
getMinZ(int halfspace)=0 -> double  

Returns the minimum z value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum z value  
";

%feature("docstring") Surface::getMinY "
getMinY(int halfspace)=0 -> double  

Returns the minimum y value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum y value  
";

%feature("docstring") Surface::getMinX "
getMinX(int halfspace)=0 -> double  

Returns the minimum x value for one of this Surface's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the Surface to consider  

Returns
-------
the minimum x value  
";

%feature("docstring") Surface::evaluate "
evaluate(const Point *point) const =0 -> double  

Evaluate a Point using the Surface's potential equation.  

This method returns the values $ f(x,y) $ for the potential function $f$ representing this
Surface.  

Parameters
----------
* point :  
    a pointer to the Soint of interest  

Returns
-------
the value of Point in the Plane's potential equation.  
";

%feature("docstring") Surface::getId "
getId() const  -> int  

Return the Surface's user-defined ID.  

Returns
-------
the Surface's user-defined ID  
";

%feature("docstring") Surface::getBoundaryType "
getBoundaryType() -> boundaryType  

Returns the type of boundary conditions for this Surface (REFLECTIVE, VACUUM or
BOUNDARY_NONE)  

Returns
-------
the type of boundary condition type for this Surface  
";

%feature("docstring") Surface::addNeighborCell "
addNeighborCell(int halfspace, Cell *cell)  

Adds a neighbor Cell to this Surface's collection of neighbors.  

Parameters
----------
* halfspace :  
    the +/-1 halfspace for the neighboring Cell  
* cell :  
    a pointer to the neighboring Cell  
";

%feature("docstring") Surface::isPointOnSurface "
isPointOnSurface(Point *point) -> bool  

Return true or false if a Point is on or off of a Surface.  

Parameters
----------
* point :  
    pointer to the Point of interest  

Returns
-------
on (true) or off (false) the Surface  
";

%feature("docstring") Surface::getMinDistance "
getMinDistance(LocalCoords *coords) -> double  

Finds the minimum distance to a Surface.  

Finds the miniumum distance to a Surface from a LocalCoords with a trajectory defined by
an angle to this Surface. If the trajectory will not intersect the Surface, returns
INFINITY.  

Parameters
----------
* coords :  
    a pointer to a localcoords object  

Returns
-------
the minimum distance to the Surface  
";

%feature("docstring") Surface::setBoundaryType "
setBoundaryType(const boundaryType boundary_type)  

Sets the boundary condition type (ie, VACUUM or REFLECTIVE) for this Surface.  

Parameters
----------
* boundary_type :  
    the boundary condition type for this Surface  
";

%feature("docstring") Surface::getSurfaceType "
getSurfaceType() -> surfaceType  

Return the type of Surface (ie, XPLANE, ZCYLINDER, etc).  

Returns
-------
the Surface type  
";

%feature("docstring") Surface::getName "
getName() const  -> char *  

Return the user-defined name of the Surface.  

Returns
-------
the Surface name  
";

%feature("docstring") Surface::isCoordOnSurface "
isCoordOnSurface(LocalCoords *coord) -> bool  

Return true or false if a LocalCoord is on or off of a Surface.  

Parameters
----------
* coord :  
    pointer to the LocalCoord of interest  

Returns
-------
on (true) or off (false) the Surface  
";

%feature("docstring") Surface::toString "
toString()=0 -> std::string  

Converts this Surface's attributes to a character array.  

The character array returned conatins the type of Surface (ie, PLANE) and the coefficients
in the potential equation.  

Returns
-------
a character array of this Surface's attributes  
";

%feature("docstring") Surface::printString "
printString()  

Prints a string representation of all of the Surface's objects to the console.  
";

%feature("docstring") Surface::setName "
setName(const char *name)  

Sets the name of the Surface.  

Parameters
----------
* name :  
    the Surface name string  
";

%feature("docstring") Surface::intersection "
intersection(Point *point, double angle, Point *points)=0 -> int  

Finds the intersection Point with this Surface from a given Point and trajectory defined
by an angle.  

Parameters
----------
* point :  
    pointer to the Point of interest  
* angle :  
    the angle defining the trajectory in radians  
* points :  
    pointer to a Point to store the intersection Point  

Returns
-------
the number of intersection Points (0 or 1)  
";

// File: structsurface__halfspace.xml


%feature("docstring") surface_halfspace "

A surface_halfspace represents a surface pointer with associated halfspace.  

Attributes
----------
* _surface : Surface *  
    A pointer to the Surface object  

* _halfspace : int  
    The halfspace associated with this surface  

C++ includes: Cell.h
";

// File: classThis.xml


%feature("docstring") This "

templated interface for a strided iterator over a Thrust device_vector on a GPU.  

This code is taken from the Thrust examples site on 1/20/2015:
https://github.com/thrust/thrust/blob/master/examples/strided_range.cu  
";

// File: classTimer.xml


%feature("docstring") Timer "

The Timer class is for timing and profiling regions of code.  

C++ includes: src/Timer.cpp
";

%feature("docstring") Timer::getTime "
getTime() -> double  

Returns the time elapsed from startTimer() to stopTimer().  

Returns
-------
the elapsed time in seconds  
";

%feature("docstring") Timer::recordSplit "
recordSplit(const char *msg)  

Records a message corresponding to a time for the current split.  

When this method is called it assumes that the Timer has been stopped and has the current
time for the process corresponding to the message.  

Parameters
----------
* msg :  
    a msg corresponding to this time split  
";

%feature("docstring") Timer::~Timer "
~Timer()  

Destructor.  
";

%feature("docstring") Timer::Get "
Get() -> Timer *  

Returns a static instance of the Timer class.  

Returns
-------
a pointer to the static Timer class  
";

%feature("docstring") Timer::printSplits "
printSplits()  

Prints the times and messages for each split to the console.  

This method will loop through all of the Timer's splits and print a formatted message
string (80 characters in length) to the console with the message and the time
corresponding to that message.  
";

%feature("docstring") Timer::Timer "
Timer()  

Constructor sets the current split elapsed time to zero.  
";

%feature("docstring") Timer::getSplit "
getSplit(const char *msg) -> double  

Returns the time associated with a particular split.  

If the split does not exist, returns 0.  

Parameters
----------
* msg :  
    the message tag for the split  

Returns
-------
the time recorded for the split (seconds)  
";

%feature("docstring") Timer::startTimer "
startTimer()  

Starts the Timer.  

This method is similar to starting a stopwatch.  
";

%feature("docstring") Timer::stopTimer "
stopTimer()  

Stops the Timer.  

This method is similar to stopping a stopwatch.  
";

%feature("docstring") Timer::clearSplit "
clearSplit(const char *msg)  

Clears the time split for this message and deletes the message's entry in the Timer's
splits log.  

Parameters
----------
* msg :  
    the message tag for the split  
";

%feature("docstring") Timer::clearSplits "
clearSplits()  

Clears all times split messages from the Timer.  
";

// File: classTrack.xml


%feature("docstring") Track "

A Track represents a characteristic line across the geometry.  

A Track has particular starting and ending points on the boundaries of the geometry and an
azimuthal angle.  

C++ includes: src/Track.h
";

%feature("docstring") Track::toString "
toString() -> std::string  

Convert this Track's attributes to a character array.  

The character array returned includes the Track's starting and ending coordinates, the
azimuthal angle and azimuthal weight.  

Returns
-------
a character array of this Track's attributes  
";

%feature("docstring") Track::setAzimAngleIndex "
setAzimAngleIndex(const int index)  

Set the index for the Track's azimuthal angle index.  

The azimuthal angle index corresponds to a an array of all azimuthal angles for $ \\theta
\\in [0, \\pi] $ owned by the TrackGenerator class.  

Parameters
----------
* index :  
    the azimuthal angle index  
";

%feature("docstring") Track::clearSegments "
clearSegments()  

Deletes each of this Track's segments.  
";

%feature("docstring") Track::getAzimAngleIndex "
getAzimAngleIndex() const  -> int  

Return the index for the Track's azimuthal angle (with respect to the x-axis).  

Returns
-------
th azimuthal angle index  
";

%feature("docstring") Track::setPeriodicTrackIndex "
setPeriodicTrackIndex(const int index)  

Set the index of a track in a periodic cycle.  

Tracks form periodic track cycles as they traverse the geometry. Tracks can be arbitrarily
decomposed into periodic track cycles and this index indicates the index in a particular
cycle.  

Parameters
----------
* index :  
    of the track in a periodic cycle  
";

%feature("docstring") Track::setUid "
setUid(int uid)  

Initializes a Track's unique ID.  

This is set by the trackgenerator to correspond to the Track's location in a 2D ragged
array of all tracks.  

Parameters
----------
* uid :  
    the Track's unique ID  
";

%feature("docstring") Track::getSegment "
getSegment(int s) -> segment *  

Returns a pointer to a segment with a given index.  

Returns a pointer to the segment or ends program if Track does not have the requested
segment.  

Parameters
----------
* segment :  
    index into the Track's segments container  

Returns
-------
a pointer to the requested segment  
";

%feature("docstring") Track::isNextOut "
isNextOut() const  -> bool  

Returns whether to give the outgoing flux to the \"forward\" (false) or \"reverse\" (true)
direction of the next Track when traveling along this Track's \"reverse\" direction.  

Returns
-------
\"forward\" (false) \"reverse\" (true) direction of outgoing Track  
";

%feature("docstring") Track::setPhi "
setPhi(const double phi)  

Set the Track's azimuthal angle.  

Parameters
----------
* phi :  
    the azimuthal angle  
";

%feature("docstring") Track::Track "
Track()  
";

%feature("docstring") Track::getPeriodicTrackIndex "
getPeriodicTrackIndex() const  -> int  

Get the index of a track in a periodic cycle.  

Returns
-------
index of the track in a periodic cycle  
";

%feature("docstring") Track::setTrackOut "
setTrackOut(Track *track_out)  

Sets the track going out along this Track's \"reverse\" direction.  

Parameters
----------
* track_out :  
    pointer to the Track going out in the \"reverse\" direction  
";

%feature("docstring") Track::getTransferFluxOut "
getTransferFluxOut() const  -> bool  

Returns a boolean to indicate whether the outgoing flux along this Track's \"reverse\"
direction should be transferred to the incoming Track.  

The bool with be false for vacuum BCs and true for all other BCs.  

Returns
-------
bool indicating whether the flux should be passed when tracking in the \"reverse\"
direction.  
";

%feature("docstring") Track::setNextOut "
setNextOut(const bool next_out)  

Sets the direction in which the flux leaving this Track along its \"reverse\" direction is
passed.  

Sets whether or not to pass the outgoing flux from this Track along its \"reverse\"
direction to the \"forward\" direction (false) or \"reverse\" direction (true) of the next
Track after intersection with the geometry boundary.  

Parameters
----------
* next_out :  
    the \"forward\" (false) or \"reverse (true) direction  
";

%feature("docstring") Track::getStart "
getStart() -> Point *  

Returns a pointer to the Track's start Point.  

Returns
-------
a pointer to the Track's start Point  
";

%feature("docstring") Track::getTransferFluxIn "
getTransferFluxIn() const  -> bool  

Returns a boolean to indicate whether the outgoing flux along this Track's \"forward\"
direction should be transferred to the outgoing Track.  

The bool with be false for vacuum BCs and true for all other BCs.  

Returns
-------
bool indicating whether the flux should be passed when tracking in the \"forward\"
direction.  
";

%feature("docstring") Track::setBCOut "
setBCOut(const boundaryType bc_out)  

Sets the boundary condition for the incoming flux along the Track's \"reverse\" direction.  

The boundaryType represents vacuum (0), reflective (1), or periodic (2) boundary
conditions.  

Parameters
----------
* bc_out :  
    boundary condition for the incoming flux in the \"reverse\" direction  
";

%feature("docstring") Track::getUid "
getUid() -> int  

Return the Track's unique ID.  

Returns
-------
the Track's unique ID  
";

%feature("docstring") Track::removeSegment "
removeSegment(int index)  

Removes a segment from this Track's list of segments.  

Parameters
----------
* index :  
    the index of the segment to remove  
";

%feature("docstring") Track::getSegments "
getSegments() -> segment *  

Returns a vector of pointers to the Track's segments.  

Returns
-------
vector of segment pointers  
";

%feature("docstring") Track::getTrackIn "
getTrackIn() const  -> Track *  

Returns the incoming Track.  

Returns
-------
a pointer to the incoming Track  
";

%feature("docstring") Track::isNextIn "
isNextIn() const  -> bool  

Returns whether to give the outgoing flux to the \"forward\" (false) or \"reverse\" (true)
direction of the next Track when traveling along this Tracks's \"forward\" direction.  

Returns
-------
\"forward\" (false) \"reverse\" (true) direction of outgoing Track  
";

%feature("docstring") Track::setNextIn "
setNextIn(const bool next_in)  

Sets the direction in which the flux leaving this Track along its \"forward\" direction is
passed.  

Sets whether or not to pass the outgoing flux from this Track along its \"forward\"
direction to the \"forward\" direction (false) or \"reverse\" direction (true) of the next
Track after intersection with the geometry boundary.  

Parameters
----------
* next_in :  
    the \"forward\" (false) or \"reverse (true) direction  
";

%feature("docstring") Track::getBCOut "
getBCOut() const  -> boundaryType  

Returns the boundary condition for the flux along the Track's \"reverse\" direction.  

Returns
-------
vacuum (0), reflective (1), or periodic (2) reflective boundary conditions  
";

%feature("docstring") Track::getReflectiveTrackIndex "
getReflectiveTrackIndex() const  -> int  

Get the index of a track in a reflective cycle.  

Returns
-------
index of the track in a reflective cycle  
";

%feature("docstring") Track::getEnd "
getEnd() -> Point *  

Returns a pointer to the Track's end Point.  

Returns
-------
a pointer to the Track's end Point  
";

%feature("docstring") Track::addSegment "
addSegment(segment *to_add)  

Adds a segment pointer to this Track's list of segments.  

This method assumes that segments are added in order of their starting location from the
Track's start point.  

Parameters
----------
* to_add :  
    a pointer to the segment to add  
";

%feature("docstring") Track::setValues "
setValues(const double start_x, const double start_y, const double start_z, const double
    end_x, const double end_y, const double end_z, const double phi)  

Set the values for the Track's start and end point and angle.  

Parameters
----------
* start_x :  
    the x-coordinate at the starting point  
* start_y :  
    the y-coordinate at the starting point  
* start_z :  
    the z-coordinate at the starting point  
* end_x :  
    the x-coordinate at the ending point  
* end_y :  
    the y-coordinate at the ending point  
* end_z :  
    the z-coordinate at the ending point  
* phi :  
    the track's azimuthal angle ( $ \\theta \\in [0, \\pi] $)  
";

%feature("docstring") Track::insertSegment "
insertSegment(int index, segment *segment)  

Inserts a segment pointer into this Track's list of segments.  

This method appends the new segment directly behind another segment in the Track. This is
a helper method for the TrackGenerator::splitTracks(...) routine.  

Parameters
----------
* index :  
    the index of the segment to insert behind in the list  
* segment :  
    a pointer to the segment to insert  
";

%feature("docstring") Track::getPhi "
getPhi() const  -> double  

Return the Track's azimuthal angle (with respect to the x-axis).  

Returns
-------
the azimuthal angle $ \\theta \\in [0, \\pi] $  
";

%feature("docstring") Track::getNumSegments "
getNumSegments() -> int  

Return the number of segments along this Track.  

Returns
-------
the number of segments  
";

%feature("docstring") Track::setReflectiveTrackIndex "
setReflectiveTrackIndex(const int index)  

Set the index of a track in a reflective cycle.  

Tracks form reflective track cycles as they traverse the geometry. Tracks can be
arbitrarily decomposed into reflective track cycles and this index indicates the index in
a particular cycle.  

Parameters
----------
* index :  
    of the track in a reflective cycle  
";

%feature("docstring") Track::~Track "
~Track()  

Destructor clears the Track segments container.  
";

%feature("docstring") Track::getTrackOut "
getTrackOut() const  -> Track *  

Returns the outgoing Track.  

Returns
-------
a pointer to the outgoing Track  
";

%feature("docstring") Track::setTrackIn "
setTrackIn(Track *track_in)  

Sets the track going out along this Track's \"forward\" direction.  

Parameters
----------
* track_in :  
    pointer to the Track going out in the \"forward\" direction  
";

%feature("docstring") Track::setBCIn "
setBCIn(const boundaryType bc_in)  

Sets the boundary condition for the incoming flux along the Track's \"forward\" direction.  

The boundaryType represents vacuum (0), reflective (1), or periodic (2) boundary
conditions.  

Parameters
----------
* bc_in :  
    boundary condition for the incoming flux in the \"forward\" direction  
";

%feature("docstring") Track::getBCIn "
getBCIn() const  -> boundaryType  

Returns the boundary condition for the flux along the Track's \"forward\" direction.  

Returns
-------
vacuum (0), reflective (1), or periodic (2) reflective boundary conditions  
";

// File: classTrackGenerator.xml


%feature("docstring") TrackGenerator "

The TrackGenerator is dedicated to generating and storing Tracks which cyclically wrap
across the Geometry.  

The TrackGenerator creates Track and initializes boundary conditions (vacuum, reflective,
or periodic) for each Track.  

C++ includes: src/TrackGenerator.h
";

%feature("docstring") TrackGenerator::retrieveTrackCoords "
retrieveTrackCoords(double *coords, int num_tracks)  

Fills an array with the x,y,z coordinates for each Track.  

This class method is intended to be called by the OpenMOC Python \"plotter\" module as a
utility to assist in plotting tracks. Although this method appears to require two
arguments, in reality it only requires on due to SWIG and would be called from within
Python as follows:  


Parameters
----------
* coords :  
    an array of coords of length NUM_VALUES_PER_RETRIEVED_TRACK times the number of Tracks  
* length_coords :  
    the total number of Tracks times NUM_VALUES_PER_RETRIEVED_TRACK  
";

%feature("docstring") TrackGenerator::setNumAzim "
setNumAzim(int num_azim)  

Set the number of azimuthal angles in $ [0, 2\\pi] $.  

Parameters
----------
* num_azim :  
    the number of azimuthal angles in $ 2\\pi $  
";

%feature("docstring") TrackGenerator::getZCoord "
getZCoord() -> double  

Returns the z-coord where the 2D Tracks should be created.  

Returns
-------
the z-coord where the 2D Tracks should be created.  
";

%feature("docstring") TrackGenerator::getQuadrature "
getQuadrature() -> Quadrature *  

Returns a pointer to the Quadrature.  

Returns
-------
a pointer to the Quadrature  
";

%feature("docstring") TrackGenerator::initializeSegments "
initializeSegments()  

Initialize track segments with pointers to FSR Materials.  

This is called by the Solver at simulation time. This initialization is necessary since
Materials in each FSR may be interchanged by the user in between different simulations.
This method links each segment and fsr_data struct with the current Material found in each
FSR.  
";

%feature("docstring") TrackGenerator::getFSRLocks "
getFSRLocks() -> omp_lock_t *  

Return the array of FSR locks for atomic FSR operations.  

Returns
-------
an array of FSR locks  
";

%feature("docstring") TrackGenerator::setZCoord "
setZCoord(double z_coord)  

Sets the z-coord where the 2D Tracks should be created.  

Parameters
----------
* z_coord :  
    the z-coord where the 2D Tracks should be created.  
";

%feature("docstring") TrackGenerator::resetFSRVolumes "
resetFSRVolumes()  

Deletes the memory associated with the FSR volumes and resets it NULL.  
";

%feature("docstring") TrackGenerator::generateFSRCentroids "
generateFSRCentroids()  

Generates the numerical centroids of the FSRs.  

This routine generates the numerical centroids of the FSRs by weighting the average x and
y values of each segment in the FSR by the segment's length and azimuthal weight. The
numerical centroid fomula can be found in R. Ferrer et. al. \"Linear Source
         Approximation in CASMO 5\", PHYSOR 2012.  
";

%feature("docstring") TrackGenerator::getFSRVolume "
getFSRVolume(int fsr_id) -> FP_PRECISION  

Computes and returns the volume of an FSR.  

Parameters
----------
* fsr_id :  
    the ID for the FSR of interest  

Returns
-------
the FSR volume  
";

%feature("docstring") TrackGenerator::retrieveSegmentCoords "
retrieveSegmentCoords(double *coords, int num_segments)  

Fills an array with the x,y,z coordinates for each Track segment.  

This class method is intended to be called by the OpenMOC Python \"plotter\" module as a
utility to assist in plotting segments. Although this method appears to require two
arguments, in reality it only requires one due to SWIG and would be called from within
Python as follows:  


Parameters
----------
* coords :  
    an array of coords of length NUM_VALUES_PER_RETRIEVED_SEGMENT times the number of
    segments  
* length_coords :  
    the total number of Track segments times NUM_VALUES_PER_RETRIEVED_SEGMENT  
";

%feature("docstring") TrackGenerator::setQuadrature "
setQuadrature(Quadrature *quadrature)  

Assign a Quadrature object to the Solver.  

This routine allows use of a Quadrature with any polar angle quadrature. Alternatively,
this routine may take in any subclass of the Quadrature parent class, including
TYPolarQuad (default), LeonardPolarQuad, GLPolarQuad, etc.  

Users may assign a Quadrature object to the Solver from Python script as follows:  


Parameters
----------
* quadrature :  
    a pointer to a Quadrature object  
";

%feature("docstring") TrackGenerator::getNumTracksByParallelGroup "
getNumTracksByParallelGroup(int group) -> int  

Return the number of tracks in a given parallel track group.  

Returns
-------
the number of tracks in a given parallel track group.  
";

%feature("docstring") TrackGenerator::getFSRVolumes "
getFSRVolumes() -> FP_PRECISION *  

Returns an array of volumes indexed by FSR.  

Returns
-------
a pointer to the array of FSR volumes  
";

%feature("docstring") TrackGenerator::getPhi "
getPhi(int azim) -> double  

Returns the azimuthal angle for a given azimuthal angle index.  

Parameters
----------
* the :  
    azimuthal angle index.  

Returns
-------
the desired azimuthal angle.  
";

%feature("docstring") TrackGenerator::splitSegments "
splitSegments(FP_PRECISION max_optical_length)  

Splits Track segments into sub-segments for a user-defined maximum optical length for the
problem.  

This routine is needed so that all segment lengths fit within the exponential
interpolation table used in the MOC transport sweep.  

Parameters
----------
* max_optical_length :  
    the maximum optical length  
";

%feature("docstring") TrackGenerator::getTracks "
getTracks() -> Track **  

Returns a 2D jagged array of the Tracks.  

The first index into the array is the azimuthal angle and the second index is the Track
number for a given azimuthal angle.  

Returns
-------
the 2D jagged array of Tracks  
";

%feature("docstring") TrackGenerator::getNumY "
getNumY(int azim) -> int  

Return the number of tracks on the y-axis for a given azimuthal angle.  

Parameters
----------
* azim :  
    An azimuthal angle index  

Returns
-------
The number of Tracks on the y-axis  
";

%feature("docstring") TrackGenerator::getNumX "
getNumX(int azim) -> int  

Return the number of tracks on the x-axis for a given azimuthal angle.  

Parameters
----------
* azim :  
    An azimuthal angle index  

Returns
-------
The number of Tracks on the x-axis  
";

%feature("docstring") TrackGenerator::getTracksByParallelGroup "
getTracksByParallelGroup() -> Track **  

Returns a 1D array of Track pointers.  

The tracks in the _tracks_by_parallel_group array are organized by parallel track group.
The index into the array is also the corresponding Track's UID.  

Returns
-------
The 1D array of Track pointers  
";

%feature("docstring") TrackGenerator::setGeometry "
setGeometry(Geometry *geometry)  

Set a pointer to the Geometry to use for track generation.  

Parameters
----------
* geometry :  
    a pointer to the Geometry  
";

%feature("docstring") TrackGenerator::getMaxOpticalLength "
getMaxOpticalLength() -> FP_PRECISION  

Finds and returns the maximum optical length amongst all segments.  

Returns
-------
the maximum optical path length  
";

%feature("docstring") TrackGenerator::printTimerReport "
printTimerReport()  

Prints a report of the timing statistics to the console.  
";

%feature("docstring") TrackGenerator::~TrackGenerator "
~TrackGenerator()  

Destructor frees memory for all Tracks.  
";

%feature("docstring") TrackGenerator::correctFSRVolume "
correctFSRVolume(int fsr_id, FP_PRECISION fsr_volume)  

Assign a correct volume for some FSR.  

This routine adjusts the length of each track segment crossing a FSR such that the
integrated volume is identical to the true volume assigned by the user.  

Parameters
----------
* fsr_id :  
    the ID of the FSR of interest  
* fsr_volume :  
    the correct FSR volume to use  
";

%feature("docstring") TrackGenerator::getDesiredAzimSpacing "
getDesiredAzimSpacing() -> double  

Return the azimuthal track spacing (cm).  

This will return the user-specified azimuthal track spacing and NOT the effective track
spacing which is computed and used to generate cyclic tracks.  

Returns
-------
the azimuthal track spacing (cm)  
";

%feature("docstring") TrackGenerator::generateTracks "
generateTracks(bool store=true, bool neighbor_cells=false)  

Generates tracks for some number of azimuthal angles and track spacing.  

Computes the effective angles and track spacing. Computes the number of Tracks for each
azimuthal angle, allocates memory for all Tracks at each angle and sets each Track's
starting and ending Points, azimuthal angle, and azimuthal angle quadrature weight.
neighbor_cells whether to use neighbor cell optimizations store whether to store the
tracks to a file for reuse  
";

%feature("docstring") TrackGenerator::getNumSegments "
getNumSegments() -> int  

Return the total number of Track segments across the Geometry.  

Returns
-------
the total number of Track segments  
";

%feature("docstring") TrackGenerator::getNumThreads "
getNumThreads() -> int  

Returns the number of shared memory OpenMP threads in use.  

Returns
-------
the number of threads  
";

%feature("docstring") TrackGenerator::getNumAzim "
getNumAzim() -> int  

Return the number of azimuthal angles in $ [0, 2\\pi] $.  

Returns
-------
the number of azimuthal angles in $ 2\\pi $  
";

%feature("docstring") TrackGenerator::setNumThreads "
setNumThreads(int num_threads)  

Sets the number of shared memory OpenMP threads to use (>0).  

Parameters
----------
* num_threads :  
    the number of threads  
";

%feature("docstring") TrackGenerator::getGeometry "
getGeometry() -> Geometry *  

Return the Geometry for this TrackGenerator if one has been set.  

Returns
-------
a pointer to the Geometry  
";

%feature("docstring") TrackGenerator::getNumParallelTrackGroups "
getNumParallelTrackGroups() -> int  

Return the number of parallel track groups.  

Returns
-------
The number of parallel track groups  
";

%feature("docstring") TrackGenerator::setDesiredAzimSpacing "
setDesiredAzimSpacing(double azim_spacing)  

Set the suggested azimuthal track spacing (cm).  

Parameters
----------
* azim_spacing :  
    the suggested azimuthal track spacing  
";

%feature("docstring") TrackGenerator::TrackGenerator "
TrackGenerator(Geometry *geometry, int num_azim, double azim_spacing)  

Constructor for the TrackGenerator assigns default values.  

Parameters
----------
* geometry :  
    a pointer to a Geometry object  
* num_azim :  
    number of azimuthal angles in $ [0, 2\\pi] $  
* azim_spacing :  
    azimuthal track spacing (cm)  
";

%feature("docstring") TrackGenerator::getNumTracks "
getNumTracks() -> int  

Return the total number of Tracks generated.  

Returns
-------
The number of Tracks generated  
";

%feature("docstring") TrackGenerator::containsTracks "
containsTracks() -> bool  

Returns whether or not the TrackGenerator contains Track that are for its current number
of azimuthal angles, track spacing and geometry.  

Returns
-------
true if the TrackGenerator conatains Tracks; false otherwise  
";

// File: classTYPolarQuad.xml


%feature("docstring") TYPolarQuad "

Tabuchi-Yamamoto's polar quadrature.  

C++ includes: src/Quadrature.h
";

%feature("docstring") TYPolarQuad::precomputeWeights "
precomputeWeights(bool solve_3D)  

Calculates total weights for every azimuthal/polar combination based on the TY quadrature.  

Parameters
----------
* solve_3D :  
    Boolean indicating whether this is a 3D quadrature  
";

%feature("docstring") TYPolarQuad::setNumPolarAngles "
setNumPolarAngles(const int num_polar)  

Set the number of polar angles to initialize.  

Parameters
----------
* num_polar :  
    the number of polar angles (maximum 6)  
";

%feature("docstring") TYPolarQuad::TYPolarQuad "
TYPolarQuad()  

Dummy constructor calls the parent constructor.  
";

%feature("docstring") TYPolarQuad::initialize "
initialize()  

Routine to initialize the polar quadrature.  

This routine uses the tabulated values for the Tabuchi-Yamamoto polar angle quadrature,
including the sine thetas and weights.  
";

// File: classUniverse.xml


%feature("docstring") Universe "

A Universe represents an unbounded space in the 3D.  

A Universe contains cell which are bounded subspaces in 3D which together form the
Universe. Universes allow for complex, repeating (i.e. lattices) geometries to be simply
represented with as few data structures as possible.  

C++ includes: src/Universe.h
";

%feature("docstring") Universe::getNumCells "
getNumCells() const  -> int  

Return the number of Cells in this Universe.  

Returns
-------
the number of Cells  
";

%feature("docstring") Universe::setType "
setType(universeType type)  

Sets the Universe type to SIMPLE or LATTICE.  

Parameters
----------
* type :  
    the Universe type  
";

%feature("docstring") Universe::getMinZ "
getMinZ() -> double  

Returns the minimum reachable z-coordinate in the Universe.  

Returns
-------
the minimum reachable z-coordinate  
";

%feature("docstring") Universe::getType "
getType() -> universeType  

Return the Universe type (SIMPLE or LATTICE).  

Returns
-------
the Universe type  
";

%feature("docstring") Universe::getMinX "
getMinX() -> double  

Aggregates a list (vector) of the IDs of all Materials within the MATERIAL type Cells
filling this Universe.  

Note that this method only searches the first level of Cells below this Universe within
the nested Universe coordinate system.  

Returns
-------
a vector of Material IDs  
";

%feature("docstring") Universe::getMinY "
getMinY() -> double  

Returns the minimum reachable y-coordinate in the Universe.  

Returns
-------
the minimum reachable y-coordinate  
";

%feature("docstring") Universe::getMaxX "
getMaxX() -> double  

Returns the maximum reachable x-coordinate in the Universe.  

Returns
-------
the maximum reachable x-coordinate  
";

%feature("docstring") Universe::getMaxY "
getMaxY() -> double  

Returns the maximum reachable y-coordinate in the Universe.  

Returns
-------
the maximum reachable y-coordinate  
";

%feature("docstring") Universe::getMaxZ "
getMaxZ() -> double  

Returns the maximum reachable z-coordinate in the Universe.  

Returns
-------
the maximum reachable z-coordinate  
";

%feature("docstring") Universe::printString "
printString()  

Prints a string representation of the Universe's attributes to the console.  
";

%feature("docstring") Universe::getMinYBoundaryType "
getMinYBoundaryType() -> boundaryType  

Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum reachable
y-coordinate in the Universe.  

Returns
-------
the boundary conditions at the minimum reachable y-coordinate  
";

%feature("docstring") Universe::buildNeighbors "
buildNeighbors()  

Builds collections of neighboring Cells for all Cells in this Universe for optimized ray
tracing.  
";

%feature("docstring") Universe::getAllCells "
getAllCells() -> std::map< int, Cell * >  

Returns the std::map of Cell IDs and Cell pointers in this Universe at all nested Universe
levels.  

Returns
-------
std::map of Cell IDs and pointers  
";

%feature("docstring") Universe::toString "
toString() -> std::string  

Convert the member attributes of this Universe to a character array.  

Returns
-------
a character array representing the Universe's attributes  
";

%feature("docstring") Universe::subdivideCells "
subdivideCells(double max_radius=INFINITY)  

Subdivides all of the Material-filled Cells within this Universe into rings and angular
sectors aligned with the z-axis.  

Parameters
----------
* max_radius :  
    the maximum allowable radius used in the subdivisions  
";

%feature("docstring") Universe::getMaxYBoundaryType "
getMaxYBoundaryType() -> boundaryType  

Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum reachable
y-coordinate in the Universe.  

Returns
-------
the boundary conditions at the maximum reachable y-coordinate  
";

%feature("docstring") Universe::Universe "
Universe(const int id=-1, const char *name=\"\")  

Constructor assigns a unique and user-specified ID for the Universe.  

Parameters
----------
* id :  
    the user-specified optional Universe ID  
* name :  
    the user-specified optional Universe ID  
";

%feature("docstring") Universe::removeCell "
removeCell(Cell *cell)  

Removes a Cell from this Universe's container of Cells.  

Parameters
----------
* cell :  
    a pointer to the Cell to remove  
";

%feature("docstring") Universe::getId "
getId() const  -> int  

Return the user-specified ID for this Universe.  

Returns
-------
the user-specified Universe ID  
";

%feature("docstring") Universe::getMaxXBoundaryType "
getMaxXBoundaryType() -> boundaryType  

Returns the boundary conditions (VACUUM or REFLECTIVE) at the maximum reachable
x-coordinate in the Universe.  

Returns
-------
the boundary conditions at the maximum reachable x-coordinate  
";

%feature("docstring") Universe::getName "
getName() const  -> char *  

Return the user-defined name of the Universe.  

Returns
-------
the Universe name  
";

%feature("docstring") Universe::~Universe "
~Universe()  

Destructor clears the Cell pointers container.  
";

%feature("docstring") Universe::getUid "
getUid() const  -> int  

Returns the Universe's unique ID.  

Returns
-------
the Universe's unique ID.  
";

%feature("docstring") Universe::getMinXBoundaryType "
getMinXBoundaryType() -> boundaryType  

Returns the boundary conditions (VACUUM or REFLECTIVE) at the minimum reachable
x-coordinate in the Universe.  

Returns
-------
the boundary conditions at the minimum reachable x-coordinate  
";

%feature("docstring") Universe::getCells "
getCells() const  -> std::map< int, Cell * >  

Return the container of Cell IDs and Cell pointers in this Universe.  

Returns
-------
std::map of Cell IDs  
";

%feature("docstring") Universe::isFissionable "
isFissionable() -> bool  

Returns true if the Universe contains a Cell filled by a fissionable Material and false
otherwise.  

This method should not be called prior to the calling of the
Geometry::computeFissionability() method.  

Returns
-------
true if contains a fissionable Material  
";

%feature("docstring") Universe::getAllMaterials "
getAllMaterials() -> std::map< int, Material * >  

Returns the std::map of all IDs and Material pointers filling this Universe.  

Returns
-------
std::map of Material IDs and pointers  
";

%feature("docstring") Universe::addCell "
addCell(Cell *cell)  

Adds a Cell to this Universe.  

Stores the user-specified Cell ID and Cell pointer in a std::map along with all of other
Cells added to this Universe.  

Parameters
----------
* cell :  
    the Cell pointer  
";

%feature("docstring") Universe::setFissionability "
setFissionability(bool fissionable)  

Sets whether or not this Universe contains a fissionable Material with a non-zero fission
cross-section.  

This method is called by the Geometry::computeFissionability() class method.  

Parameters
----------
* fissionable :  
    true if the Universe contains a fissionable Material; false otherwise  
";

%feature("docstring") Universe::findCell "
findCell(LocalCoords *coords) -> Cell *  

Finds the Cell for which a LocalCoords object resides.  

Finds the Cell that a LocalCoords object is located inside by checking each of this
Universe's Cells. Returns NULL if the LocalCoords is not in any of the Cells.  

Parameters
----------
* coords :  
    a pointer to the LocalCoords of interest  

Returns
-------
a pointer the Cell where the LocalCoords is located  
";

%feature("docstring") Universe::getCell "
getCell(int cell_id) -> Cell *  

Returns a Cell in this universe.  

Parameters
----------
* cell_id :  
    the integer the cell_id  

Returns
-------
Returns the cell pointer.  
";

%feature("docstring") Universe::clone "
clone() -> Universe *  

Clones this Universe and copy cells map.  

Returns
-------
a pointer to the Universe clone  
";

%feature("docstring") Universe::getAllUniverses "
getAllUniverses() -> std::map< int, Universe * >  

Returns the std::map of all nested Universe IDs and Universe pointers filling this
Universe.  

Returns
-------
std::map of Universe IDs and pointers  
";

%feature("docstring") Universe::setName "
setName(const char *name)  

Sets the name of the Universe.  

Parameters
----------
* name :  
    the Universe name string  
";

// File: classVector.xml


%feature("docstring") Vector "
";

%feature("docstring") Vector::getNumX "
getNumX() -> int  

Get the number of cells in the x dimension.  

Returns
-------
The number of cells in the x dimension.  
";

%feature("docstring") Vector::getNumY "
getNumY() -> int  

Get the number of cells in the y dimension.  

Returns
-------
The number of cells in the y dimension.  
";

%feature("docstring") Vector::getSum "
getSum() -> FP_PRECISION  

Get the sum of all the values in the vector.  

Returns
-------
The sum of all the values in the vector.  
";

%feature("docstring") Vector::getNumRows "
getNumRows() -> int  

Get the number of rows in the vector.  

Returns
-------
The number of rows in the vector.  
";

%feature("docstring") Vector::setValues "
setValues(int cell, int group_start, int group_end, FP_PRECISION *vals)  

Set values in the vector.  This method takes a cell, first group, last group, and floating
point value. The cell and groups are used to compute the rows in the vector. If a values
exist for the rows, the values are overwritten.  

Parameters
----------
* cell :  
    The cell location.  
* group_first :  
    The first group location to set.  
* group_last :  
    The last group location to set.  
* vals :  
    The values used to set the row locations.  
";

%feature("docstring") Vector::incrementValue "
incrementValue(int cell, int group, FP_PRECISION val)  

Increment a value in the vector.  This method takes a cell and group and floating point
value. The cell and group are used to compute the row and column in the vector. If a value
exists for the row, the value is incremented by val; otherwise, it is set to val.  

Parameters
----------
* cell :  
    The cell location.  
* group :  
    The group location.  
* val :  
    The value used to increment the row location.  
";

%feature("docstring") Vector::Vector "
Vector(omp_lock_t *cell_locks, int num_x=1, int num_y=1, int num_groups=1)  

Constructor initializes Vector object as a floating point array and sets the vector
dimensions.  The vector is ordered by cell (as opposed to by group) on the outside to be
consistent with the Matrix object. Locks are used to make the vector object thread-safe
against concurrent writes the same value. One lock locks out multiple rows of the vector
at a time representing multiple groups in the same cell.  

Parameters
----------
* cell_locks :  
    OpenMP locks for atomic cell operations.  
* num_x :  
    The number of cells in the x direction.  
* num_y :  
    The number of cells in the y direction.  
* num_groups :  
    The number of energy groups in each cell.  
";

%feature("docstring") Vector::getCellLocks "
getCellLocks() -> omp_lock_t *  

Return the array of cell locks for atomic cell operations.  

Returns
-------
an array of cell locks  
";

%feature("docstring") Vector::clear "
clear()  

Clear all values in the vector.  
";

%feature("docstring") Vector::printString "
printString()  

Print the vector object to the log file.  
";

%feature("docstring") Vector::copyTo "
copyTo(Vector *vector)  

Copy the values from the current vector to an input vector.  

Parameters
----------
* vector :  
    The vector to copy values to.  
";

%feature("docstring") Vector::getValue "
getValue(int cell, int group) -> FP_PRECISION  

Get a value at location described by a given cell and group index.  

Parameters
----------
* cell :  
    The cell location index.  
* group :  
    The group location index.  
";

%feature("docstring") Vector::scaleByValue "
scaleByValue(FP_PRECISION val)  

Scales the vector by a given value.  

Parameters
----------
* val :  
    The value to scale the vector by.  
";

%feature("docstring") Vector::getArray "
getArray() -> FP_PRECISION *  

Get the array describing the vector.  

Returns
-------
The array describing the vector.  
";

%feature("docstring") Vector::getNumGroups "
getNumGroups() -> int  

Get the number of groups in each cell.  

Returns
-------
The number of groups in each cell.  
";

%feature("docstring") Vector::setAll "
setAll(FP_PRECISION val)  
";

%feature("docstring") Vector::incrementValues "
incrementValues(int cell, int group_start, int group_end, FP_PRECISION *vals)  

Increment values in the vector.  This method takes a cell, first group, last group, and
floating point value. The cell and groups are used to compute the rows in the vector. If
values exist for the rows, the values are incremented by vals; otherwise, they are set.  

Parameters
----------
* cell :  
    The cell location.  
* group_first :  
    The first group location to increment.  
* group_last :  
    The last group location to increment.  
* vals :  
    The values used to increment the row locations.  
";

%feature("docstring") Vector::setValue "
setValue(int cell, int group, FP_PRECISION val)  

Set a value in the vector.  This method takes a cell and group and floating point value.
The cell and group are used to compute the row and column in the vector. The location of
the corresponding row is set to val.  

Parameters
----------
* cell :  
    The cell location.  
* group :  
    The group location.  
* val :  
    The value used to set the row location.  
";

%feature("docstring") Vector::~Vector "
~Vector()  

Destructor deletes the arrays used to represent the vector.  
";

// File: classVectorizedSolver.xml


%feature("docstring") VectorizedSolver "

This is a subclass of the CPUSolver class which uses memory-aligned data structures and
Intel's auto-vectorization.  

note: This class is only compiled if the Intel compiler is used when building OpenMOC. If
    building OpenMOC with the \"--with-icpc\" flag, then this class will be available in
    the \"openmoc.intel.single\" or \"openmoc.intel.double\" Python module.  

C++ includes: src/VectorizedSolver.h
";

%feature("docstring") VectorizedSolver::initializeSourceArrays "
initializeSourceArrays()  

Allocates memory for FSR source arrays.  

Deletes memory for old source arrays if they were allocated for a previous simulation.  
";

%feature("docstring") VectorizedSolver::getNumVectorWidths "
getNumVectorWidths() -> int  

Returns the number of vector lengths required to fit the number of energy groups.  

If the number of energy groups is 35 and the vector width is 4, this method will return 9
since 9*4 = 36 is the nearest integer greater than or equal to 35.  

Returns
-------
The number of vector widths  
";

%feature("docstring") VectorizedSolver::initializeFluxArrays "
initializeFluxArrays()  

Allocates memory for Track boundary angular and FSR scalar fluxes.  

Deletes memory for old flux arrays if they were allocated for a previous simulation.  
";

%feature("docstring") VectorizedSolver::VectorizedSolver "
VectorizedSolver(TrackGenerator *track_generator=NULL)  

Constructor initializes NULL arrays for source, flux, etc.  

Parameters
----------
* track_generator :  
    an optional pointer to a TrackGenerator object  
";

%feature("docstring") VectorizedSolver::initializeExpEvaluator "
initializeExpEvaluator()  

Allocates memory for the exponential linear interpolation table.  
";

%feature("docstring") VectorizedSolver::initializeMaterials "
initializeMaterials(solverMode mode=ADJOINT)  

Aligns all Material cross-section data for SIMD vector instructions.  

Parameters
----------
* mode :  
    the solution type (FORWARD or ADJOINT)  
";

%feature("docstring") VectorizedSolver::~VectorizedSolver "
~VectorizedSolver()  

Destructor deletes Track boundary angular flux and and FSR scalar flux and source arrays.  
";

%feature("docstring") VectorizedSolver::initializeFixedSources "
initializeFixedSources()  

Populates array of fixed sources assigned by FSR.  
";

%feature("docstring") VectorizedSolver::setGeometry "
setGeometry(Geometry *geometry)  

Sets the Geometry for the Solver.  

Parameters
----------
* geometry :  
    a pointer to the Geometry  
";

%feature("docstring") VectorizedSolver::initializeFSRs "
initializeFSRs()  

Initializes the FSR volumes and Materials array.  
";

%feature("docstring") VectorizedSolver::computeKeff "
computeKeff()  

Compute $ k_{eff} $ from successive fission sources.  
";

%feature("docstring") VectorizedSolver::normalizeFluxes "
normalizeFluxes()  

Normalizes all FSR scalar fluxes and Track boundary angular fluxes to the total fission
source (times $ \\nu $).  
";

%feature("docstring") VectorizedSolver::addSourceToScalarFlux "
addSourceToScalarFlux()  

Add the source term contribution in the transport equation to the FSR scalar flux.  
";

%feature("docstring") VectorizedSolver::computeFSRSources "
computeFSRSources()  

Computes the total source (fission, scattering, fixed) in each FSR.  

This method computes the total source in each FSR based on this iteration's current
approximation to the scalar flux.  
";

// File: classXPlane.xml


%feature("docstring") XPlane "

Represents a Plane perpendicular to the x-axis.  

C++ includes: src/Surface.h
";

%feature("docstring") XPlane::XPlane "
XPlane(const double x, const int id=0, const char *name=\"\")  

Constructor for a Plane perpendicular to the x-axis.  

Parameters
----------
* x :  
    the location of the Plane along the x-axis  
* id :  
    the optional Surface id  
* name :  
    the optional name of the XPlane  
";

%feature("docstring") XPlane::setX "
setX(const double x)  

Set the location of this XPlane on the x-axis.  

Parameters
----------
* x :  
    the location of the XPlane on the x-axis  
";

%feature("docstring") XPlane::toString "
toString() -> std::string  

Converts this XPlane's attributes to a character array.  

The character array returned conatins the type of Plane (ie, XPLANE) and the A, B, C, and
D coefficients in the quadratic Surface equation and the location of the Plane on the
x-axis.  

Returns
-------
a character array of this XPlane's attributes  
";

%feature("docstring") XPlane::getMinX "
getMinX(int halfspace) -> double  

Returns the minimum x value for one of this XPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the XPlane to consider  

Returns
-------
the minimum x value  
";

%feature("docstring") XPlane::getMaxX "
getMaxX(int halfspace) -> double  

Returns the maximum x value for one of this XPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the XPlane to consider  

Returns
-------
the maximum x value  
";

%feature("docstring") XPlane::getX "
getX() -> double  

Returns the location of the XPlane on the x-axis.  

Returns
-------
the location of the XPlane on the x-axis  
";

// File: classYPlane.xml


%feature("docstring") YPlane "

Represents a Plane perpendicular to the y-axis.  

C++ includes: src/Surface.h
";

%feature("docstring") YPlane::setY "
setY(const double y)  

Set the location of this YPlane on the y-axis.  

Parameters
----------
* y :  
    the location of the YPlane on the y-axis  
";

%feature("docstring") YPlane::toString "
toString() -> std::string  

Converts this yplane's attributes to a character array.  

The character array returned conatins the type of Plane (ie, YPLANE) and the A, B, C, and
D coefficients in the quadratic Surface equation and the location of the Plane on the
y-axis.  

Returns
-------
a character array of this YPlane's attributes  
";

%feature("docstring") YPlane::getY "
getY() -> double  

Returns the location of the YPlane on the y-axis.  

Returns
-------
the location of the YPlane on the y-axis  
";

%feature("docstring") YPlane::YPlane "
YPlane(const double y, const int id=0, const char *name=\"\")  

Constructor for a Plane perpendicular to the y-axis.  

Parameters
----------
* y :  
    the location of the Plane along the y-axis  
* id :  
    the optional Surface id  
* name :  
    the optional Surface name  
";

%feature("docstring") YPlane::getMaxY "
getMaxY(int halfspace) -> double  

Returns the maximum y value for one of this YPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the YPlane to consider  

Returns
-------
the maximum y value  
";

%feature("docstring") YPlane::getMinY "
getMinY(int halfspace) -> double  

Returns the minimum y value for one of this YPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the YPlane to consider  

Returns
-------
the minimum y value  
";

// File: classZCylinder.xml


%feature("docstring") ZCylinder "

Represents a Cylinder with axis parallel to the z-axis.  

C++ includes: src/Surface.h
";

%feature("docstring") ZCylinder::getMinZ "
getMinZ(int halfspace) -> double  

Returns the minimum z value of -INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the minimum z value of -INFINITY  
";

%feature("docstring") ZCylinder::getRadius "
getRadius() -> double  

Return the radius of the ZCylinder.  

Returns
-------
the radius of the ZCylinder  
";

%feature("docstring") ZCylinder::getMinX "
getMinX(int halfspace) -> double  

Returns the minimum x value for one of this ZCylinder's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the minimum x value  
";

%feature("docstring") ZCylinder::getMinY "
getMinY(int halfspace) -> double  

Returns the minimum y value for one of this ZCylinder's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the minimum y value  
";

%feature("docstring") ZCylinder::getX0 "
getX0() -> double  

Return the x-coordinate of the ZCylinder's center Point.  

Returns
-------
the x-coordinate of the ZCylinder center  
";

%feature("docstring") ZCylinder::getY0 "
getY0() -> double  

Return the y-coordinate of the ZCylinder's center Point.  

Returns
-------
the y-coordinate of the ZCylinder center  
";

%feature("docstring") ZCylinder::evaluate "
evaluate(const Point *point) const  -> double  

Evaluate a Point using the ZCylinder's quadratic Surface equation.  

Parameters
----------
* point :  
    a pointer to the Point of interest  

Returns
-------
the value of Point in the equation  
";

%feature("docstring") ZCylinder::toString "
toString() -> std::string  

Converts this ZCylinder's attributes to a character array.  

The character array returned conatins the type of Plane (ie, ZCYLINDER) and the A, B, C, D
and E coefficients in the quadratic Surface equation.  

Returns
-------
a character array of this ZCylinder's attributes  
";

%feature("docstring") ZCylinder::getMaxX "
getMaxX(int halfspace) -> double  

Returns the maximum x value for one of this ZCylinder's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the maximum x value  
";

%feature("docstring") ZCylinder::getMaxY "
getMaxY(int halfspace) -> double  

Returns the maximum y value for one of this ZCylinder's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the maximum y value  
";

%feature("docstring") ZCylinder::getMaxZ "
getMaxZ(int halfspace) -> double  

Returns the maximum z value of INFINITY.  

Parameters
----------
* halfspace :  
    the halfspace of the ZCylinder to consider  

Returns
-------
the maximum z value of INFINITY  
";

%feature("docstring") ZCylinder::ZCylinder "
ZCylinder(const double x, const double y, const double radius, const int id=0, const char
    *name=\"\")  

constructor.  

Parameters
----------
* x :  
    the x-coordinte of the ZCylinder center  
* y :  
    the y-coordinate of the ZCylinder center  
* radius :  
    the radius of the ZCylinder  
* id :  
    the optional Surface ID  
* name :  
    the optional Surface name  
";

%feature("docstring") ZCylinder::intersection "
intersection(Point *point, double angle, Point *points) -> int  

Finds the intersection Point with this zcylinder from a given Point and trajectory defined
by an angle (0, 1, or 2 points).  

Parameters
----------
* point :  
    pointer to the Point of interest  
* angle :  
    the angle defining the trajectory in radians  
* points :  
    pointer to a an array of Points to store intersection Points  
* polar :  
    the polar angle defining the trajectory in radians  

Returns
-------
the number of intersection Points (0 or 1)  
";

// File: classZPlane.xml


%feature("docstring") ZPlane "

Represents a Plane perpendicular to the z-axis.  

C++ includes: src/Surface.h
";

%feature("docstring") ZPlane::getMinZ "
getMinZ(int halfspace) -> double  

Returns the minimum z value for one of this ZPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZPlane to consider  

Returns
-------
the minimum z value  
";

%feature("docstring") ZPlane::ZPlane "
ZPlane(const double z, const int id=0, const char *name=\"\")  

Constructor for a Plane perpendicular to the z-axis.  

Parameters
----------
* z :  
    the location of the Plane along the z-axis  
* id :  
    the optional Surface ID  
* name :  
    the optional Surface name  
";

%feature("docstring") ZPlane::getMaxZ "
getMaxZ(int halfspace) -> double  

Returns the maximum z value for one of this ZPlane's halfspaces.  

Parameters
----------
* halfspace :  
    the halfspace of the ZPlane to consider  

Returns
-------
the maximum z value  
";

%feature("docstring") ZPlane::toString "
toString() -> std::string  

Converts this ZPlane's attributes to a character array.  

The character array returned conatins the type of Plane (ie, ZPLANE) and the A, B, C, and
D coefficients in the quadratic Surface equation and the location of the Plane along the
z-axis.  

Returns
-------
a character array of this ZPlane's attributes  
";

%feature("docstring") ZPlane::setZ "
setZ(const double z)  

Set the location of this ZPlane on the z-axis.  

Parameters
----------
* z :  
    the location of the ZPlane on the z-axis  
";

%feature("docstring") ZPlane::getZ "
getZ() -> double  

Returns the location of the ZPlane on the z-axis.  

Returns
-------
the location of the ZPlane on the z-axis  
";

// File: clone_8cu.xml

%feature("docstring") clone_track "
clone_track(Track *track_h, dev_track *track_d, std::map< int, int >
    &material_IDs_to_indices)  

Given a pointer to a Track on the host, a dev_track on the GPU, and the map of material
IDs to indices in the _materials array, copy all of the class attributes and segments from
the Track object on the host to the GPU.  

This routine is called by the GPUSolver::initializeTracks() private class method and is
not intended to be called directly.  

Parameters
----------
* track_h :  
    pointer to a Track on the host  
* track_d :  
    pointer to a dev_track on the GPU  
* material_IDs_to_indices :  
    map of material IDs to indices in the _materials array.  
";

%feature("docstring") clone_material "
clone_material(Material *material_h, dev_material *material_d)  

Given a pointer to a Material on the host and a dev_material on the GPU, copy all of the
properties from the Material object on the host struct to the GPU.  

This routine is called by the GPUSolver::initializeMaterials() private class method and is
not intended to be called directly.  

Parameters
----------
* material_h :  
    pointer to a Material on the host  
* material_d :  
    pointer to a dev_material on the GPU  
";

// File: clone_8h.xml

%feature("docstring") clone_track "
clone_track(Track *track_h, dev_track *track_d, std::map< int, int >
    &material_IDs_to_indices)  

Given a pointer to a Track on the host, a dev_track on the GPU, and the map of material
IDs to indices in the _materials array, copy all of the class attributes and segments from
the Track object on the host to the GPU.  

This routine is called by the GPUSolver::initializeTracks() private class method and is
not intended to be called directly.  

Parameters
----------
* track_h :  
    pointer to a Track on the host  
* track_d :  
    pointer to a dev_track on the GPU  
* material_IDs_to_indices :  
    map of material IDs to indices in the _materials array.  
";

%feature("docstring") clone_material "
clone_material(Material *material_h, dev_material *material_d)  

Given a pointer to a Material on the host and a dev_material on the GPU, copy all of the
properties from the Material object on the host struct to the GPU.  

This routine is called by the GPUSolver::initializeMaterials() private class method and is
not intended to be called directly.  

Parameters
----------
* material_h :  
    pointer to a Material on the host  
* material_d :  
    pointer to a dev_material on the GPU  
";

// File: GPUQuery_8cu.xml

%feature("docstring") get_num_threads_per_warp "
get_num_threads_per_warp() -> int  

Returns the number of threads in a CUDA warp for the attached GPU.  
";

%feature("docstring") machine_contains_gpu "
machine_contains_gpu() -> bool  

Queries a node to determine whether it contains one or more GPUs.  

Returns
-------
True if the node contains a GPU, false otherwise.  
";

%feature("docstring") print_basic_gpu_info "
print_basic_gpu_info()  

Prints the basic device info for the CUDA-enabled device in use.  

Prints the name, compute capability, # multiprocessors and the clock rate of the device.  
";

%feature("docstring") print_detailed_gpu_info "
print_detailed_gpu_info()  

Prints the detailed device info for the CUDA-enabled device in use.  

Prints the total global and constant memory, shared memory and registers per
multiprocessor, # threads per warp, maximum # threads per multiprocessor, maximum #
threads per block, maximum threadblock dimensions, and maximum grid dimensions.  
";

%feature("docstring") attach_gpu "
attach_gpu(int id)  

Sets the primary CUDA-enabled device to be the GPU with a given ID.  

Parameters
----------
* id :  
    the ID for the GPU to attache  
";

// File: GPUQuery_8h.xml

%feature("docstring") get_num_threads_per_warp "
get_num_threads_per_warp() -> int  

Returns the number of threads in a CUDA warp for the attached GPU.  
";

%feature("docstring") machine_contains_gpu "
machine_contains_gpu() -> bool  

Queries a node to determine whether it contains one or more GPUs.  

Returns
-------
True if the node contains a GPU, false otherwise.  
";

%feature("docstring") print_basic_gpu_info "
print_basic_gpu_info()  

Prints the basic device info for the CUDA-enabled device in use.  

Prints the name, compute capability, # multiprocessors and the clock rate of the device.  
";

%feature("docstring") print_detailed_gpu_info "
print_detailed_gpu_info()  

Prints the detailed device info for the CUDA-enabled device in use.  

Prints the total global and constant memory, shared memory and registers per
multiprocessor, # threads per warp, maximum # threads per multiprocessor, maximum #
threads per block, maximum threadblock dimensions, and maximum grid dimensions.  
";

%feature("docstring") attach_gpu "
attach_gpu(int id=0)  

Sets the primary CUDA-enabled device to be the GPU with a given ID.  

Parameters
----------
* id :  
    the ID for the GPU to attache  
";

// File: GPUSolver_8cu.xml

%feature("docstring") addSourceToScalarFluxOnDevice "
addSourceToScalarFluxOnDevice(FP_PRECISION *scalar_flux, FP_PRECISION *reduced_sources,
    FP_PRECISION *FSR_volumes, int *FSR_materials, dev_material *materials) -> __global__
    void  

Add the source term contribution in the transport equation to the FSR scalar flux on the
GPU.  

Parameters
----------
* scalar_flux :  
    an array of FSR scalar fluxes  
* reduced_sources :  
    an array of FSR sources / total xs  
* FSR_volumes :  
    an array of FSR volumes  
* FSR_materials :  
    an array of FSR material indices  
* materials :  
    an array of dev_material pointers  
";

%feature("docstring") computeFSRFissionRatesOnDevice "
computeFSRFissionRatesOnDevice(FP_PRECISION *FSR_volumes, int *FSR_materials, dev_material
    *materials, FP_PRECISION *scalar_flux, FP_PRECISION *fission) -> __global__ void  

Compute the total volume-intergrated fission source from all FSRs and energy groups.  

Parameters
----------
* FSR_volumes :  
    an array of the FSR volumes  
* FSR_materials :  
    an array of the FSR Material indices  
* materials :  
    an array of the dev_material pointers  
* scalar_flux :  
    an array of FSR scalar fluxes  
* fission :  
    an array of FSR nu-fission rates  
";

%feature("docstring") transferBoundaryFlux "
transferBoundaryFlux(dev_track *curr_track, int azim_index, FP_PRECISION *track_flux,
    FP_PRECISION *boundary_flux, int energy_angle_index, bool direction) -> __device__
    void  

Updates the boundary flux for a Track given boundary conditions.  

For reflective and periodic boundary conditions, the outgoing boundary flux for the Track
is given to the corresponding reflecting or periodic Track. For vacuum boundary
conditions, the outgoing flux is tallied as leakage. Note: Only one energy group is
transferred by this routine.  

Parameters
----------
* curr_track :  
    a pointer to the Track of interest  
* azim_index :  
    a pointer to the azimuthal angle index for this segment  
* track_flux :  
    an array of the outgoing Track flux  
* boundary_flux :  
    an array of all angular fluxes  
* weights :  
    an array of Quadrature weights  
* energy_angle_index :  
    the energy group index  
* direction :  
    the Track direction (forward - true, reverse - false)  
";

%feature("docstring") computeFSRFissionSourcesOnDevice "
computeFSRFissionSourcesOnDevice(int *FSR_materials, dev_material *materials, bool
    divide_sigma_t, FP_PRECISION *scalar_flux, FP_PRECISION *reduced_sources) ->
    __global__ void  

Computes the total fission source in each FSR in each energy group.  

This method is a helper routine for the openmoc.krylov submodule. This routine computes
the total fission source in each FSR. If the divide_sigma_t parameter is true then the
fission source will be divided by the total cross-section in each FSR.  

Parameters
----------
* FSR_materials :  
    an array of FSR Material indices  
* materials :  
    an array of dev_material pointers  
* divide_sigma_t :  
    a boolean indicating whether to divide by the total xs  
* scalar_flux :  
    an array of FSR scalar fluxes  
* reduced_sources :  
    an array of FSR fission sources  
";

%feature("docstring") atomicAdd "
atomicAdd(double *address, double val) -> __device__ double  

Perform an atomic addition in double precision to an array address.  

This method is straight out of CUDA C Developers Guide (cc 2013).  

Parameters
----------
* address :  
    the array memory address  
* val :  
    the value to add to the array  

Returns
-------
the atomically added array value and input value  
";

%feature("docstring") computeFissionSourcesOnDevice "
computeFissionSourcesOnDevice(FP_PRECISION *FSR_volumes, int *FSR_materials, dev_material
    *materials, FP_PRECISION *scalar_flux, FP_PRECISION *fission_sources) -> __global__
    void  

Compute the total fission source from all FSRs.  

Parameters
----------
* FSR_volumes :  
    an array of FSR volumes  
* FSR_materials :  
    an array of FSR Material indices  
* materials :  
    an array of dev_materials on the device  
* scalar_flux :  
    the scalar flux in each FSR and energy group  
* fission_sources :  
    array of fission sources in each FSR and energy group  
";

%feature("docstring") tallyScalarFlux "
tallyScalarFlux(dev_segment *curr_segment, int azim_index, int energy_group, dev_material
    *materials, FP_PRECISION *track_flux, FP_PRECISION *reduced_sources, FP_PRECISION
    *scalar_flux) -> __device__ void  

Computes the contribution to the FSR scalar flux from a Track segment in a single energy
group.  

This method integrates the angular flux for a Track segment across energy groups and polar
angles, and tallies it into the FSR scalar flux, and updates the Track's angular flux.  

Parameters
----------
* curr_segment :  
    a pointer to the Track segment of interest  
* azim_index :  
    a pointer to the azimuthal angle index for this segment  
* energy_group :  
    the energy group of interest  
* materials :  
    the array of dev_material pointers  
* track_flux :  
    a pointer to the Track's angular flux  
* reduced_sources :  
    the array of FSR sources / total xs  
* scalar_flux :  
    the array of FSR scalar fluxes  
";

%feature("docstring") computeFSRScatterSourcesOnDevice "
computeFSRScatterSourcesOnDevice(int *FSR_materials, dev_material *materials, bool
    divide_sigma_t, FP_PRECISION *scalar_flux, FP_PRECISION *reduced_sources) ->
    __global__ void  

Computes the total scattering source in each FSR and energy group.  

This method is a helper routine for the openmoc.krylov submodule. This routine computes
the total scatter source in each FSR. If the divide_sigma_t parameter is true then the
scatter source will be divided by the total cross-section in each FSR.  

Parameters
----------
* FSR_materials :  
    an array of FSR Material indices  
* materials :  
    an array of dev_material pointers  
* divide_sigma_t :  
    a boolean indicating whether to divide by the total xs  
* scalar_flux :  
    an array of FSR scalar fluxes  
* reduced_sources :  
    an array of FSR scatter sources  
";

%feature("docstring") transportSweepOnDevice "
transportSweepOnDevice(FP_PRECISION *scalar_flux, FP_PRECISION *boundary_flux,
    FP_PRECISION *reduced_sources, dev_material *materials, dev_track *tracks, int
    tid_offset, int tid_max) -> __global__ void  

This method performs one transport sweep of one halfspace of all azimuthal angles, tracks,
segments, polar angles and energy groups.  

The method integrates the flux along each track and updates the boundary fluxes for the
corresponding output Track, while updating the scalar flux in each FSR.  

Parameters
----------
* scalar_flux :  
    an array of FSR scalar fluxes  
* boundary_flux :  
    an array of Track boundary fluxes  
* reduced_sources :  
    an array of FSR sources / total xs  
* materials :  
    an array of dev_material pointers  
* tracks :  
    an array of Tracks  
* tid_offset :  
    the Track offset for azimuthal angle halfspace  
* tid_max :  
    the upper bound on the Track IDs for this azimuthal angle halfspace  
";

%feature("docstring") computeFSRSourcesOnDevice "
computeFSRSourcesOnDevice(int *FSR_materials, dev_material *materials, FP_PRECISION
    *scalar_flux, FP_PRECISION *fixed_sources, FP_PRECISION *reduced_sources, FP_PRECISION
    inverse_k_eff) -> __global__ void  

Computes the total source (fission, scattering, fixed) in each FSR.  

This method computes the total source in each region based on this iteration's current
approximation to the scalar flux.  

Parameters
----------
* FSR_materials :  
    an array of FSR Material indices  
* materials :  
    an array of dev_material pointers  
* scalar_flux :  
    an array of FSR scalar fluxes  
* fixed_sources :  
    an array of fixed (user-defined) sources  
* reduced_sources :  
    an array of FSR sources / total xs  
* inverse_k_eff :  
    the inverse of keff  
";

// File: GPUSolver_8h.xml

// File: DeviceMaterial_8h.xml

// File: DeviceTrack_8h.xml

// File: boundary__type_8h.xml

// File: Cell_8cpp.xml

%feature("docstring") maximize_cell_id "
maximize_cell_id(int cell_id)  

Maximize the auto-generated unique Cell ID counter.  

This method updates the auto-generated unique Cell ID counter if the input parameter is
greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Cell IDs do not collide with those created in OpenMC.  

Parameters
----------
* cell_id :  
    the id assigned to the auto-generated counter  
";

%feature("docstring") reset_cell_id "
reset_cell_id()  

Resets the auto-generated unique Cell ID counter to 10000.  
";

%feature("docstring") cell_id "
cell_id() -> int  

Returns an auto-generated unique Cell ID.  

This method is intended as a utility method for users writing OpenMOC input files. The
method makes use of a static Cell ID which is incremented each time the method is called
to enable unique generation of monotonically increasing IDs. The method's first ID begins
at 10000. Hence, user-defined Cell IDs greater than or equal to 10000 are prohibited.  
";

// File: Cell_8h.xml

%feature("docstring") maximize_cell_id "
maximize_cell_id(int cell_id)  

Maximize the auto-generated unique Cell ID counter.  

This method updates the auto-generated unique Cell ID counter if the input parameter is
greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Cell IDs do not collide with those created in OpenMC.  

Parameters
----------
* cell_id :  
    the id assigned to the auto-generated counter  
";

%feature("docstring") reset_cell_id "
reset_cell_id()  

Resets the auto-generated unique Cell ID counter to 10000.  
";

%feature("docstring") cell_id "
cell_id() -> int  

Returns an auto-generated unique Cell ID.  

This method is intended as a utility method for users writing OpenMOC input files. The
method makes use of a static Cell ID which is incremented each time the method is called
to enable unique generation of monotonically increasing IDs. The method's first ID begins
at 10000. Hence, user-defined Cell IDs greater than or equal to 10000 are prohibited.  
";

// File: Cmfd_8cpp.xml

// File: Cmfd_8h.xml

%feature("docstring") stencilCompare "
stencilCompare(const std::pair< int, FP_PRECISION > &firstElem, const std::pair< int,
    FP_PRECISION > &secondElem) -> bool  

Comparitor for sorting k-nearest stencil std::pair objects  
";

// File: constants_8h.xml

// File: CPUSolver_8cpp.xml

// File: CPUSolver_8h.xml

// File: ExpEvaluator_8cpp.xml

// File: ExpEvaluator_8h.xml

// File: Geometry_8cpp.xml

%feature("docstring") reset_auto_ids "
reset_auto_ids()  

Resets the auto-generated unique IDs for Materials, Surfaces, Cells and Universes/Lattices
to 10000.  
";

// File: Geometry_8h.xml

%feature("docstring") reset_auto_ids "
reset_auto_ids()  

Resets the auto-generated unique IDs for Materials, Surfaces, Cells and Universes/Lattices
to 10000.  
";

// File: linalg_8cpp.xml

%feature("docstring") matrixMultiplication "
matrixMultiplication(Matrix *A, Vector *X, Vector *B)  

Performs a matrix vector multiplication.  

This function takes in a Matrix (A), a variable Vector (X), and a solution Vector (B) and
computes the matrix vector product. The solution Vector is modified in place.  

Parameters
----------
* A :  
    a Matrix object  
* X :  
    the variable Vector object  
* B :  
    the solution Vector object  
";

%feature("docstring") eigenvalueSolve "
eigenvalueSolve(Matrix *A, Matrix *M, Vector *X, FP_PRECISION tol, FP_PRECISION
    SOR_factor) -> FP_PRECISION  

Solves a generalized eigenvalue problem using the Power method.  

This function takes in a loss + streaming Matrix (A), a fission gain Matrix (M), a flux
Vector (X), a tolerance used for both the power method and linear solve convergence (tol),
and a successive over-relaxation factor (SOR_factor) and computes the dominant eigenvalue
and eigenvector using the Power method. The eigenvalue is returned and the input X Vector
is modified in place to be the corresponding eigenvector.  

Parameters
----------
* A :  
    the loss + streaming Matrix object  
* M :  
    the fission gain Matrix object  
* X :  
    the flux Vector object  
* tol :  
    the power method and linear solve source convergence threshold  
* SOR_factor :  
    the successive over-relaxation factor  

Returns
-------
k_eff the dominant eigenvalue  
";

%feature("docstring") computeRMSE "
computeRMSE(Vector *X, Vector *Y, bool integrated) -> FP_PRECISION  

Computes the Root Mean Square Error of two Vectors.  

This function takes in two vectors (X and Y) and computes the Root Mean Square Error of
the Vector Y with respect to Vector X. The boolean integrated must also be given to
indicate whether the operation on the vector should be group-wise integrated before
performing the RMSE operation.  

Parameters
----------
* X :  
    a Vector object  
* Y :  
    a second Vector object  
* integrated :  
    a boolean indicating whether to group-wise integrate.  
";

%feature("docstring") linearSolve "
linearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, FP_PRECISION tol, FP_PRECISION
    SOR_factor)  

Solves a linear system using Red-Black Gauss Seidel with successive over-relaxation.  

This function takes in a loss + streaming Matrix (A), a fission gain Matrix (M), a flux
Vector (X), a source Vector (B), a source convergence tolerance (tol) and a successive
over-relaxation factor (SOR_factor) and computes the solution to the linear system. The
input X Vector is modified in place to be the solution vector.  

Parameters
----------
* A :  
    the loss + streaming Matrix object  
* M :  
    the fission gain Matrix object  
* X :  
    the flux Vector object  
* B :  
    the source Vector object  
* tol :  
    the power method and linear solve source convergence threshold  
* SOR_factor :  
    the successive over-relaxation factor  
";

// File: linalg_8h.xml

%feature("docstring") matrix_transpose "
matrix_transpose(T *matrix, int dim1, int dim2)  

Transpose a 2D matrix.  

Parameters
----------
* matrix :  
    array to transpose  
* dim1 :  
    first dimension length  
* dim2 :  
    second dimension length  
";

%feature("docstring") matrixMultiplication "
matrixMultiplication(Matrix *A, Vector *X, Vector *B)  

Performs a matrix vector multiplication.  

This function takes in a Matrix (A), a variable Vector (X), and a solution Vector (B) and
computes the matrix vector product. The solution Vector is modified in place.  

Parameters
----------
* A :  
    a Matrix object  
* X :  
    the variable Vector object  
* B :  
    the solution Vector object  
";

%feature("docstring") eigenvalueSolve "
eigenvalueSolve(Matrix *A, Matrix *M, Vector *X, FP_PRECISION tol, FP_PRECISION
    SOR_factor=1.5) -> FP_PRECISION  

Solves a generalized eigenvalue problem using the Power method.  

This function takes in a loss + streaming Matrix (A), a fission gain Matrix (M), a flux
Vector (X), a tolerance used for both the power method and linear solve convergence (tol),
and a successive over-relaxation factor (SOR_factor) and computes the dominant eigenvalue
and eigenvector using the Power method. The eigenvalue is returned and the input X Vector
is modified in place to be the corresponding eigenvector.  

Parameters
----------
* A :  
    the loss + streaming Matrix object  
* M :  
    the fission gain Matrix object  
* X :  
    the flux Vector object  
* tol :  
    the power method and linear solve source convergence threshold  
* SOR_factor :  
    the successive over-relaxation factor  

Returns
-------
k_eff the dominant eigenvalue  
";

%feature("docstring") computeRMSE "
computeRMSE(Vector *x, Vector *y, bool integrated) -> FP_PRECISION  

Computes the Root Mean Square Error of two Vectors.  

This function takes in two vectors (X and Y) and computes the Root Mean Square Error of
the Vector Y with respect to Vector X. The boolean integrated must also be given to
indicate whether the operation on the vector should be group-wise integrated before
performing the RMSE operation.  

Parameters
----------
* X :  
    a Vector object  
* Y :  
    a second Vector object  
* integrated :  
    a boolean indicating whether to group-wise integrate.  
";

%feature("docstring") linearSolve "
linearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, FP_PRECISION tol, FP_PRECISION
    SOR_factor=1.5)  

Solves a linear system using Red-Black Gauss Seidel with successive over-relaxation.  

This function takes in a loss + streaming Matrix (A), a fission gain Matrix (M), a flux
Vector (X), a source Vector (B), a source convergence tolerance (tol) and a successive
over-relaxation factor (SOR_factor) and computes the solution to the linear system. The
input X Vector is modified in place to be the solution vector.  

Parameters
----------
* A :  
    the loss + streaming Matrix object  
* M :  
    the fission gain Matrix object  
* X :  
    the flux Vector object  
* B :  
    the source Vector object  
* tol :  
    the power method and linear solve source convergence threshold  
* SOR_factor :  
    the successive over-relaxation factor  
";

// File: LocalCoords_8cpp.xml

// File: LocalCoords_8h.xml

// File: log_8cpp.xml

%feature("docstring") initialize_logger "
initialize_logger()  

Initializes the logger for use.  

This should be immediately called when the logger is imported into Python and before any
of its other routines are called. The routine initializes an OpenMP mutual exclusion lock
which is used to preclude race conditions from occurring when an ERROR message is reported
and program execution is terminated.  
";

%feature("docstring") get_header_character "
get_header_character() -> char  

Returns the character used to format HEADER type log messages.  

Returns
-------
the character used for HEADER type log messages  
";

%feature("docstring") set_line_length "
set_line_length(int length)  

Sets the maximum line length for log messages.  

Messages longer than this amount will be broken up into multiline messages.  

Parameters
----------
* length :  
    the maximum log message line length in characters  
";

%feature("docstring") get_log_filename "
get_log_filename() -> const char *  

Returns the log filename.  

Returns
-------
a character array for the log filename  
";

%feature("docstring") set_separator_character "
set_separator_character(char c)  

Sets the character to be used when printing SEPARATOR log messages.  

Parameters
----------
* c :  
    the character for SEPARATOR log messages  
";

%feature("docstring") set_log_level "
set_log_level(const char *new_level)  

Sets the minimum log message level which will be printed to the console and to the log
file.  

Parameters
----------
* new_level :  
    the minimum logging level as a character array  
";

%feature("docstring") get_separator_character "
get_separator_character() -> char  

Returns the character used to format SEPARATOR log messages.  

Returns
-------
the character used for SEPARATOR log messages  
";

%feature("docstring") log_printf "
log_printf(logLevel level, const char *format,...)  

Print a formatted message to the console.  

If the logging level is ERROR, this function will throw a runtime exception  

Parameters
----------
* level :  
    the logging level for this message  
* format :  
    variable list of C++ formatted arguments  
";

%feature("docstring") get_output_directory "
get_output_directory() -> const char *  

Returns the output directory for log files.  

Returns
-------
a character array for the log file directory  
";

%feature("docstring") set_log_filename "
set_log_filename(char *filename)  

Sets the name for the log file.  

Parameters
----------
* filename :  
    a character array for log filename  
";

%feature("docstring") get_log_level "
get_log_level() -> const char *  

Return the minimum level for log messages printed to the screen.  

Returns
-------
the minimum level for log messages  
";

%feature("docstring") set_output_directory "
set_output_directory(char *directory)  

Sets the output directory for log files.  

If the directory does not exist, it creates it for the user.  

Parameters
----------
* directory :  
    a character array for the log file directory  
";

%feature("docstring") set_header_character "
set_header_character(char c)  

Sets the character to be used when printing HEADER log messages.  

Parameters
----------
* c :  
    the character for HEADER log messages  
";

%feature("docstring") set_title_character "
set_title_character(char c)  

Sets the character to be used when printing TITLE log messages.  

Parameters
----------
* c :  
    the character for TITLE log messages  
";

%feature("docstring") get_title_character "
get_title_character() -> char  

Returns the character used to format TITLE log messages.  

Returns
-------
the character used for TITLE log messages  
";

%feature("docstring") create_multiline_msg "
create_multiline_msg(std::string level, std::string message) -> std::string  

Breaks up a message which is too long for a single line into a multiline message.  

This is an internal function which is called by log_printf and should not be called
directly by the user.  

Parameters
----------
* level :  
    a string containing log level prefix  
* message :  
    a string containing the log message  

Returns
-------
a string with a formatted multiline message  
";

// File: log_8h.xml

%feature("docstring") initialize_logger "
initialize_logger()  

Initializes the logger for use.  

This should be immediately called when the logger is imported into Python and before any
of its other routines are called. The routine initializes an OpenMP mutual exclusion lock
which is used to preclude race conditions from occurring when an ERROR message is reported
and program execution is terminated.  
";

%feature("docstring") get_header_character "
get_header_character() -> char  

Returns the character used to format HEADER type log messages.  

Returns
-------
the character used for HEADER type log messages  
";

%feature("docstring") set_line_length "
set_line_length(int length)  

Sets the maximum line length for log messages.  

Messages longer than this amount will be broken up into multiline messages.  

Parameters
----------
* length :  
    the maximum log message line length in characters  
";

%feature("docstring") get_log_filename "
get_log_filename() -> const char *  

Returns the log filename.  

Returns
-------
a character array for the log filename  
";

%feature("docstring") set_separator_character "
set_separator_character(char c)  

Sets the character to be used when printing SEPARATOR log messages.  

Parameters
----------
* c :  
    the character for SEPARATOR log messages  
";

%feature("docstring") set_log_level "
set_log_level(const char *new_level)  

Sets the minimum log message level which will be printed to the console and to the log
file.  

Parameters
----------
* new_level :  
    the minimum logging level as a character array  
";

%feature("docstring") get_separator_character "
get_separator_character() -> char  

Returns the character used to format SEPARATOR log messages.  

Returns
-------
the character used for SEPARATOR log messages  
";

%feature("docstring") log_printf "
log_printf(logLevel level, const char *format,...)  

Print a formatted message to the console.  

If the logging level is ERROR, this function will throw a runtime exception  

Parameters
----------
* level :  
    the logging level for this message  
* format :  
    variable list of C++ formatted arguments  
";

%feature("docstring") get_output_directory "
get_output_directory() -> const char *  

Returns the output directory for log files.  

Returns
-------
a character array for the log file directory  
";

%feature("docstring") set_log_filename "
set_log_filename(char *filename)  

Sets the name for the log file.  

Parameters
----------
* filename :  
    a character array for log filename  
";

%feature("docstring") get_log_level "
get_log_level() -> const char *  

Return the minimum level for log messages printed to the screen.  

Returns
-------
the minimum level for log messages  
";

%feature("docstring") set_output_directory "
set_output_directory(char *directory)  

Sets the output directory for log files.  

If the directory does not exist, it creates it for the user.  

Parameters
----------
* directory :  
    a character array for the log file directory  
";

%feature("docstring") set_header_character "
set_header_character(char c)  

Sets the character to be used when printing HEADER log messages.  

Parameters
----------
* c :  
    the character for HEADER log messages  
";

%feature("docstring") set_err "
set_err(const char *msg)  

A function stub used to convert C++ exceptions into Python exceptions through SWIG.  

This method is not defined in the C++ source. It is defined in the SWIG inteface files
(i.e., openmoc/openmoc.i)  

Parameters
----------
* msg :  
    a character array for the exception message  
";

%feature("docstring") set_title_character "
set_title_character(char c)  

Sets the character to be used when printing TITLE log messages.  

Parameters
----------
* c :  
    the character for TITLE log messages  
";

%feature("docstring") get_title_character "
get_title_character() -> char  

Returns the character used to format TITLE log messages.  

Returns
-------
the character used for TITLE log messages  
";

%feature("docstring") create_multiline_msg "
create_multiline_msg(std::string level, std::string message) -> std::string  

Breaks up a message which is too long for a single line into a multiline message.  

This is an internal function which is called by log_printf and should not be called
directly by the user.  

Parameters
----------
* level :  
    a string containing log level prefix  
* message :  
    a string containing the log message  

Returns
-------
a string with a formatted multiline message  
";

// File: Material_8cpp.xml

%feature("docstring") reset_material_id "
reset_material_id()  

Resets the auto-generated unique Material ID counter to 10000.  
";

%feature("docstring") maximize_material_id "
maximize_material_id(int material_id)  

Maximize the auto-generated unique Material ID counter.  

This method updates the auto-generated unique Material ID counter if the input parameter
is greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Material IDs do not collide with those created in OpenMC.  

Parameters
----------
* material_id :  
    the id assigned to the auto-generated counter  
";

%feature("docstring") material_id "
material_id() -> int  

Returns an auto-generated unique Material ID.  

This method is intended as a utility method for user's writing OpenMOC input files. The
method makes use of a static Material ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined material IDs greater than or equal to 10000 is
prohibited.  
";

// File: Material_8h.xml

%feature("docstring") reset_material_id "
reset_material_id()  

Resets the auto-generated unique Material ID counter to 10000.  
";

%feature("docstring") maximize_material_id "
maximize_material_id(int material_id)  

Maximize the auto-generated unique Material ID counter.  

This method updates the auto-generated unique Material ID counter if the input parameter
is greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Material IDs do not collide with those created in OpenMC.  

Parameters
----------
* material_id :  
    the id assigned to the auto-generated counter  
";

%feature("docstring") material_id "
material_id() -> int  

Returns an auto-generated unique Material ID.  

This method is intended as a utility method for user's writing OpenMOC input files. The
method makes use of a static Material ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined material IDs greater than or equal to 10000 is
prohibited.  
";

// File: Matrix_8cpp.xml

// File: Matrix_8h.xml

// File: pairwise__sum_8h.xml

%feature("docstring") pairwise_sum "
pairwise_sum(T *vector, int length) -> T  

Performs a pairwise sum of an array of numbers.  

This type of summation uses a divide-and-conquer algorithm which is necessary to bound the
error for summations of large sequences of numbers.  

Parameters
----------
* vector :  
    an array of numbers  
* length :  
    the length of the array  

Returns
-------
the sum of all numbers in the array  
";

// File: ParallelHashMap_8h.xml

// File: Point_8cpp.xml

// File: Point_8h.xml

// File: Quadrature_8cpp.xml

// File: Quadrature_8h.xml

// File: Solver_8cpp.xml

// File: Solver_8h.xml

// File: Surface_8cpp.xml

%feature("docstring") surface_id "
surface_id() -> int  

Returns an auto-generated unique surface ID.  

This method is intended as a utility mehtod for user's writing OpenMOC input files. The
method makes use of a static surface ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined surface IDs greater than or equal to 10000 are
prohibited.  
";

%feature("docstring") reset_surface_id "
reset_surface_id()  

Resets the auto-generated unique Surface ID counter to 10000.  
";

%feature("docstring") maximize_surface_id "
maximize_surface_id(int surface_id)  

Maximize the auto-generated unique Surface ID counter.  

This method updates the auto-generated unique Surface ID counter if the input parameter is
greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Surface IDs do not collide with those created in OpenMC.  

Parameters
----------
* surface_id :  
    the id assigned to the auto-generated counter  
";

// File: Surface_8h.xml

%feature("docstring") surface_id "
surface_id() -> int  

Returns an auto-generated unique surface ID.  

This method is intended as a utility mehtod for user's writing OpenMOC input files. The
method makes use of a static surface ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined surface IDs greater than or equal to 10000 are
prohibited.  
";

%feature("docstring") reset_surface_id "
reset_surface_id()  

Resets the auto-generated unique Surface ID counter to 10000.  
";

%feature("docstring") maximize_surface_id "
maximize_surface_id(int surface_id)  

Maximize the auto-generated unique Surface ID counter.  

This method updates the auto-generated unique Surface ID counter if the input parameter is
greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Surface IDs do not collide with those created in OpenMC.  

Parameters
----------
* surface_id :  
    the id assigned to the auto-generated counter  
";

// File: Timer_8cpp.xml

// File: Timer_8h.xml

// File: Track_8cpp.xml

// File: Track_8h.xml

// File: TrackGenerator_8cpp.xml

// File: TrackGenerator_8h.xml

// File: Universe_8cpp.xml

%feature("docstring") universe_id "
universe_id() -> int  

Returns an auto-generated unique Universe ID.  

This method is intended as a utility method for user's writing OpenMOC input files. The
method makes use of a static Universe ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined Universe IDs greater than or equal to 10000 is
prohibited.  
";

%feature("docstring") reset_universe_id "
reset_universe_id()  

Resets the auto-generated unique Universe ID counter to 10000.  
";

%feature("docstring") maximize_universe_id "
maximize_universe_id(int universe_id)  

Maximize the auto-generated unique Universe ID counter.  

This method updates the auto-generated unique Universe ID counter if the input parameter
is greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Universe IDs do not collide with those created in OpenMC.  

Parameters
----------
* universe_id :  
    the id assigned to the auto-generated counter  
";

// File: Universe_8h.xml

%feature("docstring") universe_id "
universe_id() -> int  

Returns an auto-generated unique Universe ID.  

This method is intended as a utility method for user's writing OpenMOC input files. The
method makes use of a static Universe ID which is incremented each time the method is
called to enable unique generation of monotonically increasing IDs. The method's first ID
begins at 10000. Hence, user-defined Universe IDs greater than or equal to 10000 is
prohibited.  
";

%feature("docstring") pair_second "
pair_second(const tMap &map) -> second_t< typename tMap::value_type >  

A helper routine for the Universe::findCell() method.  

This is used to insert a Universe's Cells to the back of a vector of neighbor Cells in
Universe::findCell() routine. This works in symbiosis with the second_t struct template
defined above.  

Parameters
----------
* map :  
    a std::map iterator  

Returns
-------
the second element in the iterator (e.g., map value)  
";

%feature("docstring") reset_universe_id "
reset_universe_id()  

Resets the auto-generated unique Universe ID counter to 10000.  
";

%feature("docstring") maximize_universe_id "
maximize_universe_id(int universe_id)  

Maximize the auto-generated unique Universe ID counter.  

This method updates the auto-generated unique Universe ID counter if the input parameter
is greater than the present value. This is useful for the OpenMC compatibility module to
ensure that the auto-generated Universe IDs do not collide with those created in OpenMC.  

Parameters
----------
* universe_id :  
    the id assigned to the auto-generated counter  
";

// File: Vector_8cpp.xml

// File: Vector_8h.xml

// File: VectorizedSolver_8cpp.xml

// File: VectorizedSolver_8h.xml

// File: dir_1ea949864ab1faf245facc269e7b2721.xml

// File: dir_e42498f83ad2e028b83ea18abff69fa1.xml

// File: dir_68267d1309a1af8e8297ef4c3efbcdba.xml


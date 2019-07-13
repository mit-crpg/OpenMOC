#!/usr/bin/env python

import openmoc
import opencg
import copy
import numpy as np


# A dictionary of all OpenMOC Materials created
# Keys    - Material IDs
# Values  - Materials
OPENMOC_MATERIALS = dict()

# A dictionary of all OpenCG Materials created
# Keys    - Material IDs
# Values  - Materials
OPENCG_MATERIALS = dict()

# A dictionary of all OpenMOC Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENMOC_SURFACES = dict()

# A dictionary of all OpenCG Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENCG_SURFACES = dict()

# A dictionary of all OpenMOC Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENMOC_CELLS = dict()

# A dictionary of all OpenCG Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENCG_CELLS = dict()

# A dictionary of all OpenMOC Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENMOC_UNIVERSES = dict()

# A dictionary of all OpenCG Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENCG_UNIVERSES = dict()

# A dictionary of all OpenMOC Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENMOC_LATTICES = dict()

# A dictionary of all OpenCG Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENCG_LATTICES = dict()



def get_opencg_material(openmoc_material):

  if not isinstance(openmoc_material, openmoc.Material):
    msg = 'Unable to create an OpenCG Material from {0} ' \
          'which is not an OpenMOC Material'.format(openmoc_material)
    raise ValueError(msg)

  global OPENCG_MATERIALS
  material_id = openmoc_material.getId()

  # If this Material was already created, use it
  if material_id in OPENCG_MATERIALS:
    return OPENCG_MATERIALS[material_id]

  # Create an OpenCG Material to represent this OpenMOC Material
  name = openmoc_material.getName()
  opencg_material = opencg.Material(material_id=material_id, name=name)

  # Add the OpenMOC Material to the global collection of all OpenMOC Materials
  OPENMOC_MATERIALS[material_id] = openmoc_material

  # Add the OpenCG Material to the global collection of all OpenCG Materials
  OPENCG_MATERIALS[material_id] = opencg_material

  return opencg_material


def get_openmoc_material(opencg_material):

  if not isinstance(opencg_material, opencg.Material):
    msg = 'Unable to create an OpenMOC Material from {0} ' \
          'which is not an OpenCG Material'.format(opencg_material)
    raise ValueError(msg)

  global OPENMOC_MATERIALS
  material_id = opencg_material._id

  # If this Material was already created, use it
  if material_id in OPENMOC_MATERIALS:
    return OPENMOC_MATERIALS[material_id]

  # Create an OpenMOC Material to represent this OpenCG Material
  name = opencg_material._name
  openmoc_material = openmoc.Material(id=material_id, name=name)

  # Add the OpenMOC Material to the global collection of all OpenMOC Materials
  OPENMOC_MATERIALS[material_id] = openmoc_material

  # Add the OpenCG Material to the global collection of all OpenCG Materials
  OPENCG_MATERIALS[material_id] = opencg_material

  # FIXME
  openmoc_material.thisown = 0

  return openmoc_material


def is_opencg_surface_compatible(opencg_surface):

  if not isinstance(opencg_surface, opencg.Surface):
    msg = 'Unable to check if OpenCG Surface is compatible' \
          'since {0} is not a Surface'.format(opencg_surface)
    raise ValueError(msg)

  if opencg_surface._type in ['z-squareprism']:
    return False
  else:
    return True


def get_opencg_surface(openmoc_surface):

  if not isinstance(openmoc_surface, openmoc.Surface):
    msg = 'Unable to create an OpenCG Surface from {0} ' \
          'which is not an OpenMOC Surface'.format(openmoc_surface)
    raise ValueError(msg)

  global OPENCG_SURFACES
  surface_id = openmoc_surface.getId()

  # If this Surface was already created, use it
  if surface_id in OPENCG_SURFACES:
    return OPENCG_SURFACES[surface_id]

  # Create an OpenCG Surface to represent this OpenMOC Surface
  name = openmoc_surface.getName()

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = openmoc_surface.getBoundaryType()
  if boundary == openmoc.VACUUM:
    boundary = 'vacuum'
  elif boundary == openmoc.REFLECTIVE:
    boundary = 'reflective'
  elif boundary == openmoc.BOUNDARY_NONE:
    boundary = 'interface'

  opencg_surface = None

  surface_type = openmoc_surface.getSurfaceType()

  if surface_type == openmoc.PLANE:
    openmoc_surface = openmoc.castSurfaceToPlane(openmoc_surface)
    A = openmoc_surface.getA()
    B = openmoc_surface.getB()
    C = openmoc_surface.getC()
    opencg_surface = opencg.Plane(surface_id, name, boundary, A, B, 0, C)

  elif surface_type == openmoc.XPLANE:
    openmoc_surface = openmoc.castSurfaceToXPlane(openmoc_surface)
    x0 = openmoc_surface.getX()
    opencg_surface = opencg.XPlane(surface_id, name, boundary, x0)

  elif surface_type == openmoc.YPLANE:
    openmoc_surface = openmoc.castSurfaceToYPlane(openmoc_surface)
    y0 = openmoc_surface.getY()
    opencg_surface = opencg.YPlane(surface_id, name, boundary, y0)

  elif surface_type == openmoc.ZPLANE:
    openmoc_surface = openmoc.castSurfaceToZPlane(openmoc_surface)
    z0 = openmoc_surface.getZ()
    opencg_surface = opencg.ZPlane(surface_id, name, boundary, z0)

  elif surface_type == openmoc.ZYCLINDER:
    openmoc_surface = openmoc.castSurfaceToZCylinder(openmoc_surface)
    x0 = openmoc_surface.getX0()
    y0 = openmoc_surface.getY0()
    R = openmoc_surface.getRadius()
    opencg_surface = opencg.ZCylinder(surface_id, name, boundary, x0, y0, R)

  # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = openmoc_surface

  # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
  OPENCG_SURFACES[surface_id] = opencg_surface

  return opencg_surface


def get_openmoc_surface(opencg_surface):

  if not isinstance(opencg_surface, opencg.Surface):
    msg = 'Unable to create an OpenMoC Surface from {0} which ' \
          'is not an OpenCG Surface'.format(opencg_surface)
    raise ValueError(msg)

  global OPENMOC_SURFACES
  surface_id = opencg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMOC_SURFACES:
    return OPENMOC_SURFACES[surface_id]

  # Create an OpenMOC Surface to represent this OpenCG Surface
  name = opencg_surface._name

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = opencg_surface._boundary_type
  if boundary == 'vacuum':
    boundary = openmoc.VACUUM
  elif boundary == 'reflective':
    boundary = openmoc.REFLECTIVE
  elif boundary == 'interface':
    boundary = openmoc.BOUNDARY_NONE

  if opencg_surface._type == 'plane':
    A = opencg_surface._coeffs['A']
    B = opencg_surface._coeffs['B']
    D = opencg_surface._coeffs['D']
    openmoc_surface = openmoc.Plane(A, B, D, surface_id, name)

  elif opencg_surface._type == 'x-plane':
    x0 = opencg_surface._coeffs['x0']
    openmoc_surface = openmoc.XPlane(x0, int(surface_id), name)

  elif opencg_surface._type == 'y-plane':
    y0 = opencg_surface._coeffs['y0']
    openmoc_surface = openmoc.YPlane(y0, surface_id, name)

  elif opencg_surface._type == 'z-plane':
    z0 = opencg_surface._coeffs['z0']
    openmoc_surface = openmoc.ZPlane(z0, surface_id, name)

  elif opencg_surface._type == 'z-cylinder':
    x0 = opencg_surface._coeffs['x0']
    y0 = opencg_surface._coeffs['y0']
    R = opencg_surface._coeffs['R']
    openmoc_surface = openmoc.ZCylinder(x0, y0, R, surface_id, name)

  else:
    msg = 'Unable to create an OpenMOC Surface from an OpenCG ' \
          'Surface of type {0} since it is not a compatible ' \
          'Surface type in OpenMOC'.format(opencg_surface._type)
    raise ValueError(msg)

  # Set the boundary condition for this Surface
  openmoc_surface.setBoundaryType(boundary)

  # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = openmoc_surface

  # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
  OPENCG_SURFACES[surface_id] = opencg_surface

  # FIXME
  openmoc_surface.thisown = 0

  return openmoc_surface


def get_compatible_opencg_surfaces(opencg_surface):

  if not isinstance(opencg_surface, opencg.Surface):
    msg = 'Unable to create an OpenMOC Surface from {0} which ' \
          'is not an OpenCG Surface'.format(opencg_surface)
    raise ValueError(msg)

  global OPENMOC_SURFACES
  surface_id = opencg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMOC_SURFACES:
    return OPENMOC_SURFACES[surface_id]

  # Create an OpenMOC Surface to represent this OpenCG Surface
  name = opencg_surface._name

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = opencg_surface._boundary_type

  if opencg_surface._type == 'z-squareprism':
    x0 = opencg_surface._coeffs['x0']
    y0 = opencg_surface._coeffs['y0']
    R = opencg_surface._coeffs['R']

    # Create a list of the four planes we need
    left = opencg.XPlane(x0=x0-R, name=name)
    right = opencg.XPlane(x0=x0+R, name=name)
    bottom = opencg.YPlane(y0=y0-R, name=name)
    top = opencg.YPlane(y0=y0+R, name=name)

    # Set the boundary conditions for each Surface
    left.setBoundaryType(boundary)
    right.setBoundaryType(boundary)
    bottom.setBoundaryType(boundary)
    top.setBoundaryType(boundary)

    surfaces = [left, right, bottom, top]

  elif opencg_surface._type in ['x-cylinder', 'y-cylinder',
                                 'x-squareprism', 'y-squareprism']:
    msg = 'Unable to create a compatible OpenMOC Surface from an OpenCG ' \
          'Surface of type {0} since it is not compatible with OpenMOCs 2D ' \
          'geometry formulation on the xy-plane'.format(opencg_surface._type)
    raise ValueError(msg)

  else:
    msg = 'Unable to create a compatible OpenMOC Surface from an OpenCG ' \
          'Surface of type {0} since it already a compatible ' \
          'Surface type in OpenMOC'.format(opencg_surface._type)
    raise ValueError(msg)

  # Add the OpenMOC Surface(s) to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = surfaces

  # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
  OPENCG_SURFACES[surface_id] = opencg_surface

  return surfaces


def get_opencg_cell(openmoc_cell):

  if not isinstance(openmoc_cell, openmoc.Cell):
    msg = 'Unable to create an OpenCG Cell from {0} which ' \
          'is not an OpenMOC Cell'.format(openmoc_cell)
    raise ValueError(msg)

  global OPENCG_CELLS
  cell_id = openmoc_cell.getId()

  # If this Cell was already created, use it
  if cell_id in OPENCG_CELLS:
    return OPENCG_CELLS[cell_id]

  # Create an OpenCG Cell to represent this OpenMOC Cell
  name = openmoc_cell.getName()
  opencg_cell = opencg.Cell(cell_id, name)

  fill = openmoc_cell.getFill()
  if (openmoc_cell.getType == openmoc.MATERIAL):
    opencg_cell.setFill(get_opencg_material(fill))
  elif (openmoc_cell.getType() == openmoc.FILL):
    if isinstance(fill, openmoc.Lattice):
      opencg_cell.setFill(get_opencg_lattice(fill))
    else:
      opencg_cell.setFill(get_opencg_universe(fill))

  surfaces = openmoc_cell.getSurfaces()

  for surf_id, surface_halfspace in surfaces.items():
    halfspace = surface_halfspace._halfspace
    surface = surface_halfspace._surface
    opencg_cell.addSurface(get_opencg_surface(surface), halfspace)

  # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
  OPENMOC_CELLS[cell_id] = openmoc_cell

  # Add the OpenCG Cell to the global collection of all OpenCG Cells
  OPENCG_CELLS[cell_id] = opencg_cell

  return opencg_cell


def get_compatible_opencg_cells(opencg_cell, opencg_surface, halfspace):

  if not isinstance(opencg_cell, opencg.Cell):
    msg = 'Unable to create compatible OpenMOC Cell from {0} which ' \
          'is not an OpenCG Cell'.format(opencg_cell)
    raise ValueError(msg)

  elif not isinstance(opencg_surface, opencg.Surface):
    msg = 'Unable to create compatible OpenMOC Cell since {0} is ' \
          'not an OpenCG Surface'.format(opencg_surface)
    raise ValueError(msg)

  elif not halfspace in [-1, +1]:
    msg = 'Unable to create compatible OpenMOC Cell since {0}' \
          'is not a +/-1 halfspace'.format(halfspace)
    raise ValueError(msg)

  # Initialize an empty list for the new compatible cells
  compatible_cells = list()

  # SquarePrism Surfaces
  if opencg_surface._type == 'z-squareprism':

    # Get the compatible Surfaces (XPlanes and YPlanes)
    compatible_surfaces = get_compatible_opencg_surfaces(opencg_surface)

    opencg_cell.removeSurface(opencg_surface)

    # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
    if halfspace == -1:
      opencg_cell.addSurface(compatible_surfaces[0], +1)
      opencg_cell.addSurface(compatible_surfaces[1], -1)
      opencg_cell.addSurface(compatible_surfaces[2], +1)
      opencg_cell.addSurface(compatible_surfaces[3], -1)
      compatible_cells.append(opencg_cell)

    # If Cell is outside SquarePrism, add "outside" of Surface halfspaces
    else:

      # Create 8 Cell clones to represent each of the disjoint planar
      # Surface halfspace intersections
      num_clones = 8

      for clone_id in range(num_clones):

        # Create a cloned OpenCG Cell with Surfaces compatible with OpenMOC
        clone = opencg_cell.clone()
        compatible_cells.append(clone)

        # Top left subcell - add left XPlane, top YPlane
        if clone_id == 0:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Top center subcell - add top YPlane, left/right XPlanes
        elif clone_id == 1:
          clone.addSurface(compatible_surfaces[0], +1)
          clone.addSurface(compatible_surfaces[1], -1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Top right subcell - add top YPlane, right XPlane
        elif clone_id == 2:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[3], +1)

        # Right center subcell - add right XPlane, top/bottom YPlanes
        elif clone_id == 3:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[3], -1)
          clone.addSurface(compatible_surfaces[2], +1)

        # Bottom right subcell - add right XPlane, bottom YPlane
        elif clone_id == 4:
          clone.addSurface(compatible_surfaces[1], +1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Bottom center subcell - add bottom YPlane, left/right XPlanes
        elif clone_id == 5:
          clone.addSurface(compatible_surfaces[0], +1)
          clone.addSurface(compatible_surfaces[1], -1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Bottom left subcell - add bottom YPlane, left XPlane
        elif clone_id == 6:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[2], -1)

        # Left center subcell - add left XPlane, top/bottom YPlanes
        elif clone_id == 7:
          clone.addSurface(compatible_surfaces[0], -1)
          clone.addSurface(compatible_surfaces[3], -1)
          clone.addSurface(compatible_surfaces[2], +1)

  # Remove redundant Surfaces from the Cells
  for cell in compatible_cells:
    cell.removeRedundantSurfaces()

  # Return the list of compatible OpenCG Cells
  return compatible_cells


def make_opencg_cells_compatible(opencg_universe):

  if isinstance(opencg_universe, opencg.Lattice):
    return

  elif not isinstance(opencg_universe, opencg.Universe):
    msg = 'Unable to make compatible OpenCG Cells for {0} which ' \
          'is not an OpenCG Universe'.format(opencg_universe)
    raise ValueError(msg)

  # Check all OpenCG Cells in this Universe for compatibility with OpenMOC
  opencg_cells = opencg_universe._cells

  for cell_id, opencg_cell in opencg_cells.items():

    # Check each of the OpenCG Surfaces for OpenMOC compatibility
    surfaces = opencg_cell._surfaces

    for surface_id in surfaces:
      surface = surfaces[surface_id][0]
      halfspace = surfaces[surface_id][1]

      # If this Surface is not compatible with OpenMOC, create compatible
      # OpenCG cells with a compatible version of this OpenCG Surface
      if not is_opencg_surface_compatible(surface):

        # Get the one or more OpenCG Cells that are compatible with OpenMOC
        # NOTE: This does not necessarily make the OpenCG fully compatible.
        #       It only removes the incompatible Surface and replaces it with
        #       compatible OpenCG Surface(s). The recursive call at the end
        #       of this block is necessary in the event that there are more
        #       incompatible Surfaces in this Cell that are not accounted for.
        cells = get_compatible_opencg_cells(opencg_cell, surface, halfspace)

        # Remove the non-compatible OpenCG Cell from the Universe
        opencg_universe.removeCell(opencg_cell)

        # Add the compatible OpenCG Cells to the Universe
        opencg_universe.addCells(cells)

        # Make recursive call to look at the updated state of the
        # OpenCG Universe and return
        return make_opencg_cells_compatible(opencg_universe)

  # If all OpenCG Cells in the OpenCG Universe are compatible, return
  return


def get_openmoc_cell(opencg_cell):

  if not isinstance(opencg_cell, opencg.Cell):
    msg = 'Unable to create an OpenMOC Cell from {0} which ' \
          'is not an OpenCG Cell'.format(opencg_cell)
    raise ValueError(msg)

  global OPENMOC_CELLS
  cell_id = opencg_cell._id

  # If this Cell was already created, use it
  if cell_id in OPENMOC_CELLS:
    return OPENMOC_CELLS[cell_id]

  # Create an OpenMOC Cell to represent this OpenCG Cell
  name = opencg_cell._name

  fill = opencg_cell._fill
  openmoc_cell = openmoc.Cell(cell_id, name)
  if opencg_cell._type == 'universe':
    openmoc_cell.setFill(get_openmoc_universe(fill))
  elif opencg_cell._type == 'lattice':
    openmoc_cell.setFill(get_openmoc_lattice(fill))
  else:
    openmoc_cell.setFill(get_openmoc_material(fill))

  surfaces = opencg_cell._surfaces

  for surface_id in surfaces:
    surface = surfaces[surface_id][0]
    halfspace = int(surfaces[surface_id][1])
    openmoc_cell.addSurface(halfspace, get_openmoc_surface(surface))

  # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
  OPENMOC_CELLS[cell_id] = openmoc_cell

  # Add the OpenCG Cell to the global collection of all OpenCG Cells
  OPENCG_CELLS[cell_id] = opencg_cell

  # FIXME
  openmoc_cell.thisown = 0

  return openmoc_cell


def get_opencg_universe(openmoc_universe):

  if not isinstance(openmoc_universe, openmoc.Universe):
    msg = 'Unable to create an OpenCG Universe from {0} which ' \
          'is not an OpenMOC Universe'.format(openmoc_universe)
    raise ValueError(msg)

  global OPENCG_UNIVERSES
  universe_id = openmoc_universe.getId()

  # If this Universe was already created, use it
  if universe_id in OPENCG_UNIVERSES:
    return OPENCG_UNIVERSES[universe_id]

  # Create an OpenCG Universe to represent this OpenMOC Universe
  name = openmoc_universe.getName()
  opencg_universe = opencg.Universe(universe_id, name)

  # Convert all OpenMOC Cells in this Universe to OpenCG Cells
  openmoc_cells = openmoc_universe.getCells()

  for cell_id, openmoc_cell in openmoc_cells.items():
    opencg_cell = get_opencg_cell(openmoc_cell)
    opencg_universe.addCell(opencg_cell)

  # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
  OPENMOC_UNIVERSES[universe_id] = openmoc_universe

  # Add the OpenCG Universe to the global collection of all OpenCG Universes
  OPENCG_UNIVERSES[universe_id] = opencg_universe

  return opencg_universe


def get_openmoc_universe(opencg_universe):

  if not isinstance(opencg_universe, opencg.Universe):
    msg = 'Unable to create an OpenMOC Universe from {0} which ' \
          'is not an OpenCG Universe'.format(opencg_universe)
    raise ValueError(msg)

  global OPENMOC_UNIVERSES
  universe_id = opencg_universe._id

  # If this Universe was already created, use it
  if universe_id in OPENMOC_UNIVERSES:
    return OPENMOC_UNIVERSES[universe_id]

  # Make all OpenCG Cells and Surfaces in this Universe compatible with OpenMOC
  make_opencg_cells_compatible(opencg_universe)

  # Create an OpenMOC Universe to represent this OpenCG Universe
  name = opencg_universe._name
  openmoc_universe = openmoc.Universe(universe_id, name)

  # Convert all OpenCG Cells in this Universe to OpenMOC Cells
  opencg_cells = opencg_universe._cells

  for cell_id, opencg_cell in opencg_cells.items():
    openmoc_cell = get_openmoc_cell(opencg_cell)
    openmoc_universe.addCell(openmoc_cell)

  # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
  OPENMOC_UNIVERSES[universe_id] = openmoc_universe

  # Add the OpenCG Universe to the global collection of all OpenCG Universes
  OPENCG_UNIVERSES[universe_id] = opencg_universe

  # FIXME
  openmoc_universe.thisown = 0

  return openmoc_universe


def get_opencg_lattice(openmoc_lattice):

  if not isinstance(openmoc_lattice, openmoc.Lattice):
    msg = 'Unable to create an OpenCG Lattice from {0} which ' \
          'is not an OpenMOC Lattice'.format(openmoc_lattice)
    raise ValueError(msg)

  global OPENCG_LATTICES
  lattice_id = openmoc_lattice.getId()

  # If this Lattice was already created, use it
  if lattice_id in OPENCG_LATTICES:
    return OPENCG_LATTICES[lattice_id]

  # Create an OpenCG Lattice to represent this OpenMOC Lattice
  name = openmoc_lattice.getName()
  offset = openmoc_lattice.getOffset()
  dimension = [1, openmoc_lattice.getNumY(), openmoc_lattice.getNumX()]
  width = [1, openmoc_lattice.getWidthY(), openmoc_lattice.getWidthX()]
  lower_left = [-np.inf, width[1]*dimension[1]/2. + offset.getX(),
                width[2]*dimension[2] / 2. + offset.getY()]

  # Initialize an empty array for the OpenCG nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension)[::-1]), \
                              dtype=opencg.Universe)

  # Create OpenCG Universes for each unique nested Universe in this Lattice
  unique_universes = openmoc_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_opencg_universe(universe)

  # Build the nested Universe array
  for y in range(dimension[1]):
    for x in range(dimension[0]):
      universe = openmoc_lattice.getUniverse(x, y)
      universe_id = universe.getId()
      universe_array[0][y][x] = unique_universes[universe_id]

  opencg_lattice = opencg.Lattice(lattice_id, name)
  opencg_lattice.setDimension(dimension)
  opencg_lattice.setWidth(width)
  opencg_lattice.setUniverses(universe_array)

  offset = np.array(lower_left, dtype=np.float64) - \
           ((np.array(width, dtype=np.float64) * \
             np.array(dimension, dtype=np.float64))) / -2.0
  opencg_lattice.setOffset(offset)

  # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
  OPENMOC_LATTICES[lattice_id] = openmoc_lattice

  # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
  OPENCG_LATTICES[lattice_id] = opencg_lattice

  return opencg_lattice


def get_openmoc_lattice(opencg_lattice):

  if not isinstance(opencg_lattice, opencg.Lattice):
    msg = 'Unable to create an OpenMOC Lattice from {0} which ' \
          'is not an OpenCG Lattice'.format(opencg_lattice)
    raise ValueError(msg)

  global OPENMOC_LATTICES
  lattice_id = opencg_lattice._id

  # If this Lattice was already created, use it
  if lattice_id in OPENMOC_LATTICES:
    return OPENMOC_LATTICES[lattice_id]

  name = opencg_lattice._name
  dimension = opencg_lattice._dimension
  width = opencg_lattice._width
  offset = opencg_lattice._offset
  universes = opencg_lattice._universes

  # Initialize an empty array for the OpenMOC nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension[0:2])), \
                              dtype=openmoc.Universe)

  # Create OpenMOC Universes for each unique nested Universe in this Lattice
  unique_universes = opencg_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_openmoc_universe(universe)

  # Build the nested Universe array
  for y in range(dimension[1]):
    for x in range(dimension[0]):
      universe_id = universes[0][y][x]._id
      universe_array[x][y] = unique_universes[universe_id]

  openmoc_lattice = openmoc.Lattice(lattice_id, name)
  openmoc_lattice.setWidth(width[0], width[1])
  openmoc_lattice.setUniverses(universe_array.tolist())
  openmoc_lattice.setOffset(offset[0], offset[1])

  # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
  OPENMOC_LATTICES[lattice_id] = openmoc_lattice

  # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
  OPENCG_LATTICES[lattice_id] = opencg_lattice

  # FIXME
  openmoc_lattice.thisown = 0

  return openmoc_lattice


def get_opencg_geometry(openmoc_geometry):

  if not isinstance(openmoc_geometry, openmoc.Geometry):
    msg = 'Unable to get OpenCG geometry from {0} which is ' \
          'not an OpenMOC Geometry object'.format(openmoc_geometry)
    raise ValueError(msg)

  # Clear dictionaries and auto-generated IDs
  OPENMOC_SURFACES.clear()
  OPENCG_SURFACES.clear()
  OPENMOC_CELLS.clear()
  OPENCG_CELLS.clear()
  OPENMOC_UNIVERSES.clear()
  OPENCG_UNIVERSES.clear()
  OPENMOC_LATTICES.clear()
  OPENCG_LATTICES.clear()

  openmoc_root_universe = openmoc_geometry.getRootUniverse()
  opencg_root_universe = get_opencg_universe(openmoc_root_universe)

  opencg_geometry = opencg.Geometry()
  opencg_geometry.setRootUniverse(opencg_root_universe)
  opencg_geometry.initializeCellOffsets()

  return opencg_geometry


def get_openmoc_geometry(opencg_geometry):

  if not isinstance(opencg_geometry, opencg.Geometry):
    msg = 'Unable to get OpenMOC geometry from {0} which is ' \
          'not an OpenCG Geometry object'.format(opencg_geometry)
    raise ValueError(msg)

  # Deep copy the goemetry since it may be modified to make all Surfaces
  # compatible with OpenMOC's specifications
  opencg_geometry.assignAutoIds()
  opencg_geometry = copy.deepcopy(opencg_geometry)

  # Update Cell bounding boxes in Geometry
  opencg_geometry.updateBoundingBoxes()

  # Clear dictionaries and auto-generated ID
  OPENMOC_SURFACES.clear()
  OPENCG_SURFACES.clear()
  OPENMOC_CELLS.clear()
  OPENCG_CELLS.clear()
  OPENMOC_UNIVERSES.clear()
  OPENCG_UNIVERSES.clear()
  OPENMOC_LATTICES.clear()
  OPENCG_LATTICES.clear()

  # Make the entire geometry "compatible" before assigning auto IDs
  universes = opencg_geometry.getAllUniverses()
  for universe_id, universe in universes.items():
    make_opencg_cells_compatible(universe)

  opencg_geometry.assignAutoIds()

  opencg_root_universe = opencg_geometry._root_universe
  openmoc_root_universe = get_openmoc_universe(opencg_root_universe)

  openmoc_geometry = openmoc.Geometry()
  openmoc_geometry.setRootUniverse(openmoc_root_universe)

  # FIXME
  openmoc_geometry.thisown = 0

  return openmoc_geometry

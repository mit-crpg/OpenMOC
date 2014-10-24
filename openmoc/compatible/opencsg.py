#!/usr/bin/env python

import openmoc
import opencsg
import copy
import numpy as np


# A dictionary of all OpenMOC Materials created
# Keys    - Material IDs
# Values  - Materials
OPENMOC_MATERIALS = dict()

# A dictionary of all OpenCSG Materials created
# Keys    - Material IDs
# Values  - Materials
OPENCSG_MATERIALS = dict()

# A dictionary of all OpenMOC Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENMOC_SURFACES = dict()

# A dictionary of all OpenCSG Surfaces created
# Keys    - Surface IDs
# Values  - Surfaces
OPENCSG_SURFACES = dict()

# A dictionary of all OpenMOC Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENMOC_CELLS = dict()

# A dictionary of all OpenCSG Cells created
# Keys    - Cell IDs
# Values  - Cells
OPENCSG_CELLS = dict()

# A dictionary of all OpenMOC Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENMOC_UNIVERSES = dict()

# A dictionary of all OpenCSG Universes created
# Keys    - Universes IDs
# Values  - Universes
OPENCSG_UNIVERSES = dict()

# A dictionary of all OpenMOC Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENMOC_LATTICES = dict()

# A dictionary of all OpenCSG Lattices created
# Keys    - Lattice IDs
# Values  - Lattices
OPENCSG_LATTICES = dict()



def get_opencsg_material(openmoc_material):

  if not isinstance(openmoc_material, openmoc.Material):
    msg = 'Unable to create an OpenCSG Material from {0} ' \
          'which is not an OpenMOC Material'.format(openmoc_material)
    raise ValueError(msg)

  global OPENCSG_MATERIALS
  material_id = openmoc_material.getId()

  # If this Material was already created, use it
  if material_id in OPENCSG_MATERIALS.keys():
    return OPENCSG_MATERIALS[material_id]

  # Create an OpenCSG Material to represent this OpenMOC Material
  name = openmoc_material.getName()
  opencsg_material = opencsg.Material(material_id=material_id, name=name)

  # Add the OpenMOC Material to the global collection of all OpenMOC Materials
  OPENMOC_MATERIALS[material_id] = openmoc_material

  # Add the OpenCSG Material to the global collection of all OpenCSG Materials
  OPENCSG_MATERIALS[material_id] = opencsg_material

  return opencsg_material


def get_openmoc_material(opencsg_material):

  if not isinstance(opencsg_material, opencsg.Material):
    msg = 'Unable to create an OpenMOC Material from {0} ' \
          'which is not an OpenCSG Material'.format(opencsg_material)
    raise ValueError(msg)

  global OPENMOC_MATERIALS
  material_id = opencsg_material._id

  # If this Material was already created, use it
  if material_id in OPENMOC_MATERIALS.keys():
    return OPENMOC_MATERIALS[material_id]

  # Create an OpenMOC Material to represent this OpenCSG Material
  name = opencsg_material._name
  openmoc_material = openmoc.Material(id=material_id, name=name)

  # Add the OpenMOC Material to the global collection of all OpenMOC Materials
  OPENMOC_MATERIALS[material_id] = openmoc_material

  # Add the OpenCSG Material to the global collection of all OpenCSG Materials
  OPENCSG_MATERIALS[material_id] = opencsg_material

  return openmoc_material


def is_opencsg_surface_compatible(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to check if OpenCSG Surface is compatible' \
          'since {0} is not a Surface'.format(opencsg_surface)
    raise ValueError(msg)

  if opencsg_surface._type in ['x-squareprism', 'y-squareprism',
                               'x-cylinder', 'y-cylinder']:
    return False
  else:
    return True


def get_opencsg_surface(openmoc_surface):

  if not isinstance(openmoc_surface, openmoc.Surface):
    msg = 'Unable to create an OpenCSG Surface from {0} ' \
          'which is not an OpenMOC Surface'.format(openmoc_surface)
    raise ValueError(msg)

  global OPENCSG_SURFACES
  surface_id = openmoc_surface.getId()

  # If this Surface was already created, use it
  if surface_id in OPENCSG_SURFACES.keys():
    return OPENCSG_SURFACES[surface_id]

  # Create an OpenCSG Surface to represent this OpenMOC Surface
  name = openmoc_surface.getName()

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = openmoc_surface.getBoundaryType()
  if boundary == openmoc.VACUUM:
    boundary = 'vacuum'
  elif boundary == openmoc.REFLECTIVE:
    boundary = 'reflective'
  elif boundary == openmoc.BOUNDARY_NONE:
    boundary = 'interface'

  opencsg_surface = None

  surface_type = openmoc_surface.getSurfaceType()

  if surface_type == openmoc.PLANE:
    A = openmoc_surface.getA()
    B = openmoc_surface.getB()
    C = openmoc_surface.getC()
    opencsg_surface = opencsg.Plane(surface_id, name, boundary, A, B, 0, C)

  elif surface_type == openmoc.XPLANE:
    x0 = openmoc_surface.getX()
    opencsg_surface = opencsg.XPlane(surface_id, name, boundary, x0)

  elif surface_type == openmoc.YPLANE:
    y0 = openmoc_surface.getY()
    opencsg_surface = opencsg.YPlane(surface_id, name, boundary, y0)

  elif surface_type == openmoc.ZPLANE:
    z0 = openmoc_surface.getZ()
    opencsg_surface = opencsg.ZPlane(surface_id, name, boundary, z0)

  elif surface_type == openmoc.CIRCLE:
    x0 = openmoc_surface.getX0()
    y0 = openmoc_surface.getY0()
    R = openmoc_surface.getRadius()
    opencsg_surface = opencsg.ZCylinder(surface_id, name, boundary, x0, y0, R)

  # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = openmoc_surface

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return opencsg_surface


def get_openmoc_surface(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create an OpenMoC Surface from {0} which ' \
          'is not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  global OPENMOC_SURFACES
  surface_id = opencsg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMOC_SURFACES.keys():
    return OPENMOC_SURFACES[surface_id]

  # Create an OpenMOC Surface to represent this OpenCSG Surface
  name = opencsg_surface._name

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = opencsg_surface._boundary_type
  if boundary == 'vacuum':
    boundary = openmoc.VACUUM
  elif boundary == 'reflective':
    boundary = openmoc.REFLECTIVE
  elif boundary == 'interface':
    boundary = openmoc.BOUNDARY_NONE

  if opencsg_surface._type == 'plane':
    A = opencsg_surface._coeffs['A']
    B = opencsg_surface._coeffs['B']
    D = opencsg_surface._coeffs['D']
    openmoc_surface = openmoc.Plane(A, B, D, surface_id, name)

  elif opencsg_surface._type == 'x-plane':
    x0 = opencsg_surface._coeffs['x0']
    openmoc_surface = openmoc.XPlane(x0, surface_id, name)

  elif opencsg_surface._type == 'y-plane':
    y0 = opencsg_surface._coeffs['y0']
    openmoc_surface = openmoc.ZPlane(y0, surface_id, name)

  elif opencsg_surface._type == 'z-plane':
    z0 = opencsg_surface._coeffs['z0']
    openmoc_surface = openmoc.ZPlane(z0, surface_id, name)

  elif opencsg_surface._type == 'z-cylinder':
    x0 = opencsg_surface._coeffs['x0']
    y0 = opencsg_surface._coeffs['y0']
    R = opencsg_surface._coeffs['R']
    openmoc_surface = openmoc.Circle(x0, y0, R, surface_id, name)

  else:
    msg = 'Unable to create an OpenMOC Surface from an OpenCSG ' \
          'Surface of type {0} since it is not a compatible ' \
          'Surface type in OpenMOC'.format(opencsg_surface._type)
    raise ValueError(msg)

  # Set the boundary condition for this Surface
  openmoc_surface.setBoundaryType(boundary)

  # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = openmoc_surface

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return openmoc_surface


def get_compatible_opencsg_surfaces(opencsg_surface):

  if not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create an OpenMOC Surface from {0} which ' \
          'is not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  global OPENMOC_SURFACES
  surface_id = opencsg_surface._id

  # If this Surface was already created, use it
  if surface_id in OPENMOC_SURFACES.keys():
    return OPENMOC_SURFACES[surface_id]

  # Create an OpenMC Surface to represent this OpenCSG Surface
  name = opencsg_surface._name

  # Correct for OpenMOC's syntax for Surfaces dividing Cells
  boundary = opencsg_surface._boundary_type
  if boundary == 'vacuum':
    boundary = openmoc.VACUUM
  elif boundary == 'reflective':
    boundary = openmoc.REFLECTIVE
  elif boundary == 'interface':
    boundary = openmoc.BOUNDARY_NONE

  if opencsg_surface._type == 'z-squareprism':
    x0 = opencsg_surface._coeffs['x0']
    y0 = opencsg_surface._coeffs['y0']
    R = opencsg_surface._coeffs['R']

    # Create a list of the four planes we need
    left = opencsg.XPlane(x0-R, name)
    right = opencsg.XPlane(x0+R, name)
    bottom = opencsg.YPlane(y0-R, name)
    top = opencsg.YPlane(y0+R, name)

    # Set the boundary conditions for each Surface
    left.setBoundaryType(boundary)
    right.setBoundaryType(boundary)
    bottom.setBoundaryType(boundary)
    top.setBoundaryType(boundary)

    surfaces = [left, right, bottom, top]

  elif opencsg_surface._type in ['x-cylinder', 'y-cylinder',
                                 'x-squareprism', 'y-squareprism']:
    msg = 'Unable to create a compatible OpenMOC Surface from an OpenCSG ' \
          'Surface of type {0} since it is not compatible with OpenMOCs 2D ' \
          'geometry formulation on the xy-plane'.format(opencsg_surface._type)
    raise ValueError(msg)

  else:
    msg = 'Unable to create a compatible OpenMOC Surface from an OpenCSG ' \
          'Surface of type {0} since it already a compatible ' \
          'Surface type in OpenMC'.format(opencsg_surface._type)
    raise ValueError(msg)

  # Add the OpenMOC Surface(s) to the global collection of all OpenMOC Surfaces
  OPENMOC_SURFACES[surface_id] = surfaces

  # Add the OpenCSG Surface to the global collection of all OpenCSG Surfaces
  OPENCSG_SURFACES[surface_id] = opencsg_surface

  return surfaces


def get_opencsg_cell(openmoc_cell):

  if not isinstance(openmoc_cell, openmoc.Cell):
    msg = 'Unable to create an OpenCSG Cell from {0} which ' \
          'is not an OpenMOC Cell'.format(openmoc_cell)
    raise ValueError(msg)

  global OPENCSG_CELLS
  cell_id = openmoc_cell.getId()

  # If this Cell was already created, use it
  if cell_id in OPENCSG_CELLS.keys():
    return OPENCSG_CELLS[cell_id]

  # Create an OpenCSG Cell to represent this OpenMC Cell
  name = openmoc_cell.getName()
  opencsg_cell = opencsg.Cell(cell_id, name)

  fill = openmoc_cell._fill

  if (openmoc_cell.getType == openmoc.MATERIAL):
    opencsg_cell.setFill(get_opencsg_material(fill))
  elif (openmoc_cell._type == openmoc.FILL):
    if isinstance(fill, openmoc.Lattice):
      opencsg_cell.setFill(get_opencsg_lattice(fill))
    else:
      opencsg_cell.setFill(get_opencsg_universe(fill))

  # FIXME: How will this return to Python? Use std_map.i SWIG interface
  surfaces = openmoc_cell.getSurfaces()

  for surface, halfspace in surfaces.items():
    opencsg_cell.addSurface(get_opencsg_surface(surface), halfspace)

  # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
  OPENMOC_CELLS[cell_id] = openmoc_cell

  # Add the OpenCSG Cell to the global collection of all OpenCSG Cells
  OPENCSG_CELLS[cell_id] = opencsg_cell

  return opencsg_cell


def get_compatible_opencsg_cells(opencsg_cell, opencsg_surface, halfspace):

  if not isinstance(opencsg_cell, opencsg.Cell):
    msg = 'Unable to create compatible OpenMOC Cell from {0} which ' \
          'is not an OpenCSG Cell'.format(opencsg_cell)
    raise ValueError(msg)

  elif not isinstance(opencsg_surface, opencsg.Surface):
    msg = 'Unable to create compatible OpenMOC Cell since {0} is ' \
          'not an OpenCSG Surface'.format(opencsg_surface)
    raise ValueError(msg)

  elif not halfspace in [-1, +1]:
    msg = 'Unable to create compatible OpenMOC Cell since {0}' \
          'is not a +/-1 halfspace'.format(halfspace)
    raise ValueError(msg)

  # Initialize an empty list for the new compatible cells
  compatible_cells = list()

  # SquarePrism Surfaces
  if opencsg_surface._type == 'z-squareprism':

    # Get the compatible Surfaces (XPlanes and YPlanes)
    compatible_surfaces = get_compatible_opencsg_surfaces(opencsg_surface)

    opencsg_cell.removeSurface(opencsg_surface)

    # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
    if halfspace == -1:
      opencsg_cell.addSurface(compatible_surfaces[0], +1)
      opencsg_cell.addSurface(compatible_surfaces[1], -1)
      opencsg_cell.addSurface(compatible_surfaces[2], +1)
      opencsg_cell.addSurface(compatible_surfaces[3], -1)
      compatible_cells.append(opencsg_cell)

    # If Cell is outside SquarePrism, add "outside" of Surface halfspaces
    else:

      # Create 8 Cell clones to represent each of the disjoint planar
      # Surface halfspace intersections
      num_clones = 8

      for clone_id in range(num_clones):

        # Create a cloned OpenCSG Cell with Surfaces compatible with OpenMC
        clone = opencsg_cell.clone()
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

  # Return the list of compatible OpenCSG Cells
  return compatible_cells


def make_opencsg_cells_compatible(opencsg_universe):

  if not isinstance(opencsg_universe, opencsg.Universe):
    msg = 'Unable to make compatible OpenCSG Cells for {0} which ' \
          'is not an OpenCSG Universe'.format(opencsg_universe)
    raise ValueError(msg)

  # Check all OpenCSG Cells in this Universe for compatibility with OpenMOC
  opencsg_cells = opencsg_universe._cells

  for cell_id, opencsg_cell in opencsg_cells.items():

    # Check each of the OpenCSG Surfaces for OpenMOC compatibility
    surfaces = opencsg_cell._surfaces

    for surface_id in surfaces.keys():
      surface = surfaces[surface_id][0]
      halfspace = surfaces[surface_id][1]

      # If this Surface is not compatible with OpenMOC, create compatible
      # OpenCSG cells with a compatible version of this OpenCSG Surface
      if not is_opencsg_surface_compatible(surface):

        # Get the one or more OpenCSG Cells that are compatible with OpenMOC
        # NOTE: This does not necessarily make the OpenCSG fully compatible.
        #       It only removes the incompatible Surface and replaces it with
        #       compatible OpenCSG Surface(s). The recursive call at the end
        #       of this block is necessary in the event that there are more
        #       incompatible Surfaces in this Cell that are not accounted for.
        cells = get_compatible_opencsg_cells(opencsg_cell, surface, halfspace)

        # Remove the non-compatible OpenCSG Cell from the Universe
        opencsg_universe.removeCell(opencsg_cell)

        # Add the compatible OpenCSG Cells to the Universe
        opencsg_universe.addCells(cells)

        # Make recursive call to look at the updated state of the
        # OpenCSG Universe and return
        return make_opencsg_cells_compatible(opencsg_universe)

  # If all OpenCSG Cells in the OpenCSG Universe are compatible, return
  return


def get_openmoc_cell(opencsg_cell):

  if not isinstance(opencsg_cell, opencsg.Cell):
    msg = 'Unable to create an OpenMOC Cell from {0} which ' \
          'is not an OpenCSG Cell'.format(opencsg_cell)
    raise ValueError(msg)

  global OPENMOC_CELLS
  cell_id = opencsg_cell._id

  # If this Cell was already created, use it
  if cell_id in OPENMOC_CELLS.keys():
    return OPENMOC_CELLS[cell_id]

  # Create an OpenMOC Cell to represent this OpenCSG Cell
  name = opencsg_cell._name
  openmoc_cell = openmoc.CellFill(cell_id, name)

  fill = opencsg_cell._fill
  if opencsg_cell._type == 'universe':
    openmoc_cell.setFill(get_openmoc_universe(fill))
  elif opencsg_cell._type == 'lattice':
    openmoc_cell.setFill(get_openmoc_lattice(fill))
  else:
    openmoc_cell.setFill(get_openmoc_material(fill))

  surfaces = opencsg_cell._surfaces

  for surface_id in surfaces.keys():
    surface = surfaces[surface_id][0]
    halfspace = surfaces[surface_id][1]
    openmoc_cell.addSurface(halfspace, get_openmoc_surface(surface))

  # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
  OPENMOC_CELLS[cell_id] = openmoc_cell

  # Add the OpenCSG Cell to the global collection of all OpenCSG Cells
  OPENCSG_CELLS[cell_id] = opencsg_cell

  return openmoc_cell


def get_opencsg_universe(openmoc_universe):

  if not isinstance(openmoc_universe, openmoc.Universe):
    msg = 'Unable to create an OpenCSG Universe from {0} which ' \
          'is not an OpenMOC Universe'.format(openmoc_universe)
    raise ValueError(msg)

  global OPENCSG_UNIVERSES
  universe_id = openmoc_universe.getId()

  # If this Universe was already created, use it
  if universe_id in OPENCSG_UNIVERSES.keys():
    return OPENCSG_UNIVERSES[universe_id]

  # Create an OpenCSG Universe to represent this OpenMC Universe
  name = openmoc_universe.getName()
  opencsg_universe = opencsg.Universe(universe_id, name)

  # Convert all OpenMOC Cells in this Universe to OpenCSG Cells
  # FIXME: Need to use std_map.i SWIG interface file for this
  openmoc_cells = openmoc_universe.getCells()

  for cell_id, openmoc_cell in openmoc_cells.items():
    opencsg_cell = get_opencsg_cell(openmoc_cell)
    opencsg_universe.addCell(opencsg_cell)

  # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
  OPENMOC_UNIVERSES[universe_id] = openmoc_universe

  # Add the OpenCSG Universe to the global collection of all OpenCSG Universes
  OPENCSG_UNIVERSES[universe_id] = opencsg_universe

  return opencsg_universe


def get_openmoc_universe(opencsg_universe):

  if not isinstance(opencsg_universe, opencsg.Universe):
    msg = 'Unable to create an OpenMOC Universe from {0} which ' \
          'is not an OpenCSG Universe'.format(opencsg_universe)
    raise ValueError(msg)

  global OPENMOC_UNIVERSES
  universe_id = opencsg_universe._id

  # If this Universe was already created, use it
  if universe_id in OPENMOC_UNIVERSES.keys():
    return OPENMOC_UNIVERSES[universe_id]

  # Make all OpenCSG Cells and Surfaces in this Universe compatible with OpenMOC
  make_opencsg_cells_compatible(opencsg_universe)

  # Create an OpenMOC Universe to represent this OpenCSG Universe
  name = opencsg_universe._name
  openmoc_universe = openmoc.Universe(universe_id, name)

  # Convert all OpenCSG Cells in this Universe to OpenMC Cells
  opencsg_cells = opencsg_universe._cells

  for cell_id, opencsg_cell in opencsg_cells.items():
    openmoc_cell = get_openmoc_cell(opencsg_cell)
    openmoc_universe.addCell(openmoc_cell)

  # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
  OPENMOC_UNIVERSES[universe_id] = openmoc_universe

  # Add the OpenCSG Universe to the global collection of all OpenCSG Universes
  OPENCSG_UNIVERSES[universe_id] = opencsg_universe

  return openmoc_universe


def get_opencsg_lattice(openmoc_lattice):

  if not isinstance(openmoc_lattice, openmoc.Lattice):
    msg = 'Unable to create an OpenCSG Lattice from {0} which ' \
          'is not an OpenMC Lattice'.format(openmoc_lattice)
    raise ValueError(msg)

  global OPENCSG_LATTICES
  lattice_id = openmoc_lattice.getId()

  # If this Lattice was already created, use it
  if lattice_id in OPENCSG_LATTICES.keys():
    return OPENCSG_LATTICES[lattice_id]

  # Create an OpenCSG Lattice to represent this OpenMC Lattice
  name = openmoc_lattice.getName()
  dimension = [1, openmoc_lattice.getNumY(), openmoc_lattice.getNumX()]
  width = [1, openmoc_lattice.getWidthY(), openmoc_lattice.getWidthX()]
  lower_left = [-np.inf, width[1]*dimension[1]/2., width[2]*dimension[2] / 2.]

  # Initialize an empty array for the OpenCSG nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension)[::-1]), \
                              dtype=opencsg.Universe)

  # Create OpenCSG Universes for each unique nested Universe in this Lattice
  # FIXME: Need to implement this routine for the Lattice class
  unique_universes = openmoc_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_opencsg_universe(universe)

  # Build the nested Universe array
  for y in range(dimension[1]):
    for x in range(dimension[0]):
      universe = openmoc_lattice.getUniverse(x, y)
      universe_id = universe.getId()
      universe_array[0][y][x] = unique_universes[universe_id]

  opencsg_lattice = opencsg.Lattice(lattice_id, name)
  opencsg_lattice.setDimension(dimension)
  opencsg_lattice.setWidth(width)
  opencsg_lattice.setUniverses(universe_array)

  offset = np.array(lower_left, dtype=np.float64) - \
           ((np.array(width, dtype=np.float64) * \
             np.array(dimension, dtype=np.float64))) / -2.0
  opencsg_lattice.setOffset(offset)

  # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
  OPENMOC_LATTICES[lattice_id] = openmoc_lattice

  # Add the OpenCSG Lattice to the global collection of all OpenCSG Lattices
  OPENCSG_LATTICES[lattice_id] = opencsg_lattice

  return opencsg_lattice


def get_openmoc_lattice(opencsg_lattice):

  if not isinstance(opencsg_lattice, opencsg.Lattice):
    msg = 'Unable to create an OpenMOC Lattice from {0} which ' \
          'is not an OpenCSG Lattice'.format(opencsg_lattice)
    raise ValueError(msg)

  global OPENMOC_LATTICES
  lattice_id = opencsg_lattice._id

  # If this Lattice was already created, use it
  if lattice_id in OPENMOC_LATTICES.keys():
    return OPENMOC_LATTICES[lattice_id]

  dimension = opencsg_lattice._dimension
  width = opencsg_lattice._width
  offset = opencsg_lattice._offset
  universes = opencsg_lattice._universes

  # Initialize an empty array for the OpenMOC nested Universes in this Lattice
  universe_array = np.ndarray(tuple(np.array(dimension[2:1:-1])), \
                              dtype=openmoc.Universe)

  # Create OpenMOC Universes for each unique nested Universe in this Lattice
  unique_universes = opencsg_lattice.getUniqueUniverses()

  for universe_id, universe in unique_universes.items():
    unique_universes[universe_id] = get_openmoc_universe(universe)

  # Build the nested Universe array
  for y in range(dimension[1]):
    for x in range(dimension[0]):
      universe_id = universes[0][y][x]._id
      universe_array[x][y] = unique_universes[universe_id]

  openmoc_lattice = openmoc.Lattice(lattice_id, width[2], width[1])

  # FIXME: Need to implement ability to set an array of Universe pointers for Lattice
  openmoc_lattice.setLatticeCells(universe_array)

  # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
  OPENMOC_LATTICES[lattice_id] = openmoc_lattice

  # Add the OpenCSG Lattice to the global collection of all OpenCSG Lattices
  OPENCSG_LATTICES[lattice_id] = opencsg_lattice

  return openmoc_lattice


def get_opencsg_geometry(openmoc_geometry):

  if not isinstance(openmoc_geometry, openmoc.Geometry):
    msg = 'Unable to get OpenCSG geometry from {0} which is ' \
          'not an OpenMOC Geometry object'.format(openmoc_geometry)
    raise ValueError(msg)

  # Clear dictionaries and auto-generated IDs
  OPENMOC_SURFACES.clear()
  OPENCSG_SURFACES.clear()
  OPENMOC_CELLS.clear()
  OPENCSG_CELLS.clear()
  OPENMOC_UNIVERSES.clear()
  OPENCSG_UNIVERSES.clear()
  OPENMOC_LATTICES.clear()
  OPENCSG_LATTICES.clear()

  opencsg.geometry.reset_auto_ids()

  # FIXME: Need to implement root universe in Geometry class
  openmoc_root_universe = openmoc_geometry._root_universe
  opencsg_root_universe = get_opencsg_universe(openmoc_root_universe)

  opencsg_geometry = opencsg.Geometry()
  opencsg_geometry.setRootUniverse(opencsg_root_universe)
  opencsg_geometry.initializeCellOffsets()

  return opencsg_geometry


def get_openmoc_geometry(opencsg_geometry):

  if not isinstance(opencsg_geometry, opencsg.Geometry):
    msg = 'Unable to get OpenMOC geometry from {0} which is ' \
          'not an OpenCSG Geometry object'.format(opencsg_geometry)
    raise ValueError(msg)

  # Deep copy the goemetry since it may be modified to make all Surfaces
  # compatible with OpenMOC's specifications
  opencsg_geometry = copy.deepcopy(opencsg_geometry)

  # Update Cell bounding boxes in Geometry
  opencsg_geometry.updateBoundingBoxes()

  # Clear dictionaries and auto-generated ID
  OPENMOC_SURFACES.clear()
  OPENCSG_SURFACES.clear()
  OPENMOC_CELLS.clear()
  OPENCSG_CELLS.clear()
  OPENMOC_UNIVERSES.clear()
  OPENCSG_UNIVERSES.clear()
  OPENMOC_LATTICES.clear()
  OPENCSG_LATTICES.clear()

  openmoc.reset_auto_ids()

  # FIXME: Need to implement root universe in Geometry class
  opencsg_root_universe = opencsg_geometry._root_universe
  openmoc_root_universe = get_openmoc_universe(opencsg_root_universe)

  # FIXME: Need to implement root universe in Geometry class
  openmoc_geometry = openmoc.Geometry()
  openmoc_geometry.setRootUniverse(openmoc_root_universe)

  return openmoc_geometry

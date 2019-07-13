#!/usr/bin/env python

import copy
import sys
import numpy as np

import openmoc
import opencg

# For Python 2.X.X
if sys.version_info[0] == 2:
    import checkvalue as cv
# For Python 3.X.X
else:
    import openmoc.checkvalue as cv


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
    """Return an OpenCG material corresponding to an OpenMOC material.

    Parameters
    ----------
    openmoc_material : openmoc.Material
        OpenMOC material

    Returns
    -------
    opencg_material : opencg.Material
        Equivalent OpenCG material

    """

    cv.check_type('openmoc_material', openmoc_material, openmoc.Material)

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
    """Return an OpenMOC material corresponding to an OpenCG material.

    Parameters
    ----------
    opencg_material : opencg.Material
        OpenCG material

    Returns
    -------
    openmoc_material : openmoc.Material
        Equivalent OpenMOC material

    """

    cv.check_type('opencg_material', opencg_material, opencg.Material)

    global OPENMOC_MATERIALS
    material_id = opencg_material.id

    # If this Material was already created, use it
    if material_id in OPENMOC_MATERIALS:
        return OPENMOC_MATERIALS[material_id]

    # Create an OpenMOC Material to represent this OpenCG Material
    name = str(opencg_material.name)
    openmoc_material = openmoc.Material(id=material_id, name=name)

    # Add the OpenMOC Material to the global collection of all OpenMOC Materials
    OPENMOC_MATERIALS[material_id] = openmoc_material

    # Add the OpenCG Material to the global collection of all OpenCG Materials
    OPENCG_MATERIALS[material_id] = opencg_material

    return openmoc_material


def is_opencg_surface_compatible(opencg_surface):
    """Determine whether OpenCG surface is compatible with OpenMOC geometry.

    A surface is considered compatible if there is a one-to-one correspondence
    between OpenMOC and OpenCG surface types. Note that some OpenCG surfaces,
    e.g. SquarePrism, do not have a one-to-one correspondence with OpenMOC
    surfaces but can still be converted into an equivalent collection of
    OpenMOC surfaces.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface

    Returns
    -------
    bool
        Whether OpenCG surface is compatible with OpenMOC

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    if opencg_surface.type in ['z-squareprism']:
        return False
    else:
        return True


def get_opencg_surface(openmoc_surface):
    """Return an OpenCG surface corresponding to an OpenMOC surface.

    Parameters
    ----------
    openmc_surface : openmoc.Surface
        OpenMOC surface

    Returns
    -------
    opencg_surface : opencg.Surface
        Equivalent OpenCG surface

    """

    cv.check_type('openmoc_surface', openmoc_surface, openmoc.Surface)

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
    elif boundary == openmoc.PERIODIC:
        boundary = 'periodic'
    elif boundary == openmoc.BOUNDARY_NONE:
        boundary = 'interface'

    opencg_surface = None
    surface_type = openmoc_surface.getSurfaceType()

    if surface_type == openmoc.PLANE:
        openmoc_surface = openmoc.castSurfaceToPlane(openmoc_surface)
        A = openmoc_surface.getA()
        B = openmoc_surface.getB()
        C = openmoc_surface.getC()
        D = openmoc_surface.getD()
        opencg_surface = opencg.Plane(surface_id, name, boundary, A, B, C, D)

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

    elif surface_type == openmoc.ZCYLINDER:
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
    """Return an OpenMOC surface corresponding to an OpenCG surface.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface

    Returns
    -------
    openmoc_surface : openmoc.Surface
        Equivalent OpenMOC surface

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    global OPENMOC_SURFACES
    surface_id = opencg_surface.id

    # If this Surface was already created, use it
    if surface_id in OPENMOC_SURFACES:
        return OPENMOC_SURFACES[surface_id]

    # Create an OpenMOC Surface to represent this OpenCG Surface
    name = str(opencg_surface.name)

    # Correct for OpenMOC's syntax for Surfaces dividing Cells
    boundary = opencg_surface.boundary_type
    if boundary == 'vacuum':
        boundary = openmoc.VACUUM
    elif boundary == 'reflective':
        boundary = openmoc.REFLECTIVE
    elif boundary == 'periodic':
        boundary = openmoc.PERIODIC
    elif boundary == 'interface':
        boundary = openmoc.BOUNDARY_NONE

    if opencg_surface.type == 'plane':
        A = opencg_surface.a
        B = opencg_surface.b
        C = opencg_surface.c
        D = opencg_surface.d
        openmoc_surface = openmoc.Plane(A, B, C, D, surface_id, name)

    elif opencg_surface.type == 'x-plane':
        x0 = opencg_surface.x0
        openmoc_surface = openmoc.XPlane(x0, int(surface_id), name)

    elif opencg_surface.type == 'y-plane':
        y0 = opencg_surface.y0
        openmoc_surface = openmoc.YPlane(y0, surface_id, name)

    elif opencg_surface.type == 'z-plane':
        z0 = opencg_surface.z0
        openmoc_surface = openmoc.ZPlane(z0, surface_id, name)

    elif opencg_surface.type == 'z-cylinder':
        x0 = opencg_surface.x0
        y0 = opencg_surface.y0
        R = opencg_surface.r
        openmoc_surface = openmoc.ZCylinder(x0, y0, R, surface_id, name)

    else:
        msg = 'Unable to create an OpenMOC Surface from an OpenCG ' \
              'Surface of type {0} since it is not a compatible ' \
              'Surface type in OpenMOC'.format(opencg_surface.type)
        raise ValueError(msg)

    # Set the boundary condition for this Surface
    openmoc_surface.setBoundaryType(boundary)

    # Add the OpenMOC Surface to the global collection of all OpenMOC Surfaces
    OPENMOC_SURFACES[surface_id] = openmoc_surface

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return openmoc_surface


def get_compatible_opencg_surfaces(opencg_surface):
    """Generate OpenCG surfaces that are compatible with OpenMOC equivalent to
    an OpenCG surface that is not compatible. For example, this method may be
    used to convert a ZSquarePrism OpenCG surface into a collection of
    equivalent XPlane and YPlane OpenCG surfaces.

    Parameters
    ----------
    opencg_surface : opencg.Surface
        OpenCG surface that is incompatible with OpenMOC

    Returns
    -------
    surfaces : list of opencg.Surface
        Collection of surfaces equivalent to the original one but compatible
        with OpenMOC

    """

    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)

    global OPENMOC_SURFACES
    surface_id = opencg_surface.id

    # If this Surface was already created, use it
    if surface_id in OPENMOC_SURFACES:
        return OPENMOC_SURFACES[surface_id]

    # Create an OpenMOC Surface to represent this OpenCG Surface
    name = str(opencg_surface.name)

    # Correct for OpenMOC's syntax for Surfaces dividing Cells
    boundary = opencg_surface.boundary_type

    if opencg_surface.type == 'x-squareprism':
        y0 = opencg_surface.y0
        z0 = opencg_surface.z0
        R = opencg_surface.r

        # Create a list of the four planes we need
        min_y = opencg.YPlane(y0=y0-R, name=name)
        max_y = opencg.YPlane(y0=y0+R, name=name)
        min_z = opencg.ZPlane(z0=z0-R, name=name)
        max_z = opencg.ZPlane(z0=z0+R, name=name)

        # Set the boundary conditions for each Surface
        min_y.boundary_type = boundary
        max_y.boundary_type = boundary
        min_z.boundary_type = boundary
        max_z.boundary_type = boundary

        surfaces = [min_y, max_y, min_z, max_z]

    elif opencg_surface.type == 'y-squareprism':
        x0 = opencg_surface.x0
        z0 = opencg_surface.z0
        R = opencg_surface.r

        # Create a list of the four planes we need
        min_x = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        max_x = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        min_z = opencg.ZPlane(name=name, boundary=boundary, z0=z0-R)
        max_z = opencg.ZPlane(name=name, boundary=boundary, z0=z0+R)

        # Set the boundary conditions for each Surface
        min_x.boundary_type = boundary
        max_x.boundary_type = boundary
        min_z.boundary_type = boundary
        max_z.boundary_type = boundary

        surfaces = [min_x, max_x, min_z, max_z]

    elif opencg_surface.type == 'z-squareprism':
        x0 = opencg_surface.x0
        y0 = opencg_surface.y0
        R = opencg_surface.r

        # Create a list of the four planes we need
        min_x = opencg.XPlane(name=name, boundary=boundary, x0=x0-R)
        max_x = opencg.XPlane(name=name, boundary=boundary, x0=x0+R)
        min_y = opencg.YPlane(name=name, boundary=boundary, y0=y0-R)
        max_y = opencg.YPlane(name=name, boundary=boundary, y0=y0+R)

        # Set the boundary conditions for each Surface
        min_x.boundary_type = boundary
        max_x.boundary_type = boundary
        min_y.boundary_type = boundary
        max_y.boundary_type = boundary

        surfaces = [min_x, max_x, min_y, max_y]

    else:
        msg = 'Unable to create a compatible OpenMOC Surface an OpenCG ' \
              'Surface of type "{0}" since it already a compatible ' \
              'Surface type in OpenMOC'.format(opencg_surface.type)
        raise ValueError(msg)

    # Add the OpenMOC Surface(s) to global collection of all OpenMOC Surfaces
    OPENMOC_SURFACES[surface_id] = surfaces

    # Add the OpenCG Surface to the global collection of all OpenCG Surfaces
    OPENCG_SURFACES[surface_id] = opencg_surface

    return surfaces


def get_opencg_cell(openmoc_cell):
    """Return an OpenCG cell corresponding to an OpenMOC cell.

    Parameters
    ----------
    openmoc_cell : openmoc.Cell
        OpenMOC cell

    Returns
    -------
    opencg_cell : opencg.Cell
        Equivalent OpenCG cell

    """

    cv.check_type('openmoc_cell', openmoc_cell, openmoc.Cell)

    global OPENCG_CELLS
    cell_id = openmoc_cell.getId()

    # If this Cell was already created, use it
    if cell_id in OPENCG_CELLS:
        return OPENCG_CELLS[cell_id]

    # Create an OpenCG Cell to represent this OpenMOC Cell
    name = openmoc_cell.getName()
    opencg_cell = opencg.Cell(cell_id, name)

    if (openmoc_cell.getType() == openmoc.MATERIAL):
        fill = openmoc_cell.getFillMaterial()
        opencg_cell.fill = get_opencg_material(fill)
    elif (openmoc_cell.getType() == openmoc.FILL):
        fill = openmoc_cell.getFillUniverse()
        if isinstance(fill, openmoc.Lattice):
            opencg_cell.fill = get_opencg_lattice(fill)
        else:
            opencg_cell.fill = get_opencg_universe(fill)

    if openmoc_cell.isRotated():
        rotation = openmoc_cell.getRotation(3)
        opencg_cell.rotation = rotation
    if openmoc_cell.isTranslated():
        translation = openmoc_cell.getTranslation(3)
        opencg_cell.translation = translation

    surfaces = openmoc_cell.getSurfaces()

    for surf_id, surface_halfspace in surfaces.items():
        halfspace = surface_halfspace._halfspace
        surface = surface_halfspace._surface
        opencg_cell.add_surface(get_opencg_surface(surface), halfspace)

    # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
    OPENMOC_CELLS[cell_id] = openmoc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return opencg_cell


def get_compatible_opencg_cells(opencg_cell, opencg_surface, halfspace):
    """Generate OpenCG cells that are compatible with OpenMOC equivalent to an
    OpenCG cell that is not compatible.

    Parameters
    ----------
    opencg_cell : opencg.Cell
        OpenCG cell
    opencg_surface : opencg.Surface
        OpenCG surface that causes the incompatibility, e.g. an instance of
        XSquarePrism
    halfspace : {-1, 1}
        Which halfspace defined by the surface is contained in the cell

    Returns
    -------
    compatible_cells : list of opencg.Cell
        Collection of cells equivalent to the original one but compatible with
        OpenMC

    """

    cv.check_type('opencg_cell', opencg_cell, opencg.Cell)
    cv.check_type('opencg_surface', opencg_surface, opencg.Surface)
    cv.check_value('halfspace', halfspace, (-1, +1))

    # Initialize an empty list for the new compatible cells
    compatible_cells = list()

    # SquarePrism Surfaces
    if opencg_surface.type in ['x-squareprism', 'y-squareprism',
                               'z-squareprism']:

        # Get the compatible Surfaces (XPlanes and YPlanes)
        compatible_surfaces = get_compatible_opencg_surfaces(opencg_surface)

        opencg_cell.remove_surface(opencg_surface)

        # If Cell is inside SquarePrism, add "inside" of Surface halfspaces
        if halfspace == -1:
            opencg_cell.add_surface(compatible_surfaces[0], +1)
            opencg_cell.add_surface(compatible_surfaces[1], -1)
            opencg_cell.add_surface(compatible_surfaces[2], +1)
            opencg_cell.add_surface(compatible_surfaces[3], -1)
            compatible_cells.append(opencg_cell)

        # If Cell is outside the SquarePrism (positive halfspace), add "outside"
        # of Surface halfspaces. Since OpenMOC does not have a SquarePrism
        # Surface, individual Cells are created for the 8 Cells that make up the
        # outer region of a SquarePrism.
        #                 |                    |
        #           0     |        1           |    2
        #           ______|____________________|______
        #                 |     SquarePrism    |
        #           7     |   (-)  halfspace   |    3
        #           ______|____________________|______
        #                 |                    |
        #           6     |        5           |    4
        #                 |                    |
        else:

          # Create 8 Cell clones to represent each of the disjoint planar
          # Surface halfspace intersections
          num_clones = 8

          for clone_id in range(num_clones):

              # Create a cloned OpenCG Cell with Surfaces compatible with OpenMOC
              clone = opencg_cell.clone()
              compatible_cells.append(clone)

              # Top left subcell (subcell 0)
              if clone_id == 0:
                  clone.add_surface(compatible_surfaces[0], -1)
                  clone.add_surface(compatible_surfaces[3], +1)

              # Top center subcell (subcell 1)
              elif clone_id == 1:
                  clone.add_surface(compatible_surfaces[0], +1)
                  clone.add_surface(compatible_surfaces[1], -1)
                  clone.add_surface(compatible_surfaces[3], +1)

              # Top right subcell (subcell 2)
              elif clone_id == 2:
                  clone.add_surface(compatible_surfaces[1], +1)
                  clone.add_surface(compatible_surfaces[3], +1)

              # Right center subcell (subcell 3)
              elif clone_id == 3:
                  clone.add_surface(compatible_surfaces[1], +1)
                  clone.add_surface(compatible_surfaces[3], -1)
                  clone.add_surface(compatible_surfaces[2], +1)

              # Bottom right subcell (subcell 4)
              elif clone_id == 4:
                  clone.add_surface(compatible_surfaces[1], +1)
                  clone.add_surface(compatible_surfaces[2], -1)

              # Bottom center subcell (subcell 5)
              elif clone_id == 5:
                  clone.add_surface(compatible_surfaces[0], +1)
                  clone.add_surface(compatible_surfaces[1], -1)
                  clone.add_surface(compatible_surfaces[2], -1)

              # Bottom left subcell (subcell 6)
              elif clone_id == 6:
                  clone.add_surface(compatible_surfaces[0], -1)
                  clone.add_surface(compatible_surfaces[2], -1)

              # Left center subcell (subcell 7)
              elif clone_id == 7:
                  clone.add_surface(compatible_surfaces[0], -1)
                  clone.add_surface(compatible_surfaces[3], -1)
                  clone.add_surface(compatible_surfaces[2], +1)

    # Remove redundant Surfaces from the Cells
    for cell in compatible_cells:
        cell.remove_redundant_surfaces()

    # Return the list of compatible OpenCG Cells
    return compatible_cells


def make_opencg_cells_compatible(opencg_universe):
    """Make all cells in an OpenCG universe compatible with OpenMOC.

    Parameters
    ----------
    opencg_universe : opencg.Universe
        Universe to check

    """

    if isinstance(opencg_universe, opencg.Lattice):
        return

    cv.check_type('opencg_universe', opencg_universe, opencg.Universe)

    # Check all OpenCG Cells in this Universe for compatibility with OpenMOC
    opencg_cells = opencg_universe.cells

    for cell_id, opencg_cell in opencg_cells.items():

        # Check each of the OpenCG Surfaces for OpenMOC compatibility
        surfaces = opencg_cell.surfaces

        for surface_id in surfaces:
            surface = surfaces[surface_id][0]
            halfspace = surfaces[surface_id][1]

            # If this Surface is not compatible with OpenMOC, create compatible
            # OpenCG cells with a compatible version of this OpenCG Surface
            if not is_opencg_surface_compatible(surface):

                # Get the one or more OpenCG Cells compatible with OpenMOC
                # NOTE: This does not necessarily make OpenCG fully compatible.
                # It only removes the incompatible Surface and replaces it with
                # compatible OpenCG Surface(s). The recursive call at the end
                # of this block is necessary in the event that there are more
                # incompatible Surfaces in this Cell that are not accounted for.
                cells = \
                    get_compatible_opencg_cells(opencg_cell, surface, halfspace)

                # Remove the non-compatible OpenCG Cell from the Universe
                opencg_universe.remove_cell(opencg_cell)

                # Add the compatible OpenCG Cells to the Universe
                opencg_universe.add_cells(cells)

                # Make recursive call to look at the updated state of the
                # OpenCG Universe and return
                return make_opencg_cells_compatible(opencg_universe)

    # If all OpenCG Cells in the OpenCG Universe are compatible, return
    return


def get_openmoc_cell(opencg_cell):
    """Return an OpenMOC cell corresponding to an OpenCG cell.

    Parameters
    ----------
    opencg_cell : opencg.Cell
        OpenCG cell

    Returns
    -------
    openmoc_cell : openmoc.Cell
        Equivalent OpenMOC cell

    """

    cv.check_type('openmoc_cell', opencg_cell, opencg.Cell)

    global OPENMOC_CELLS
    cell_id = opencg_cell.id

    # If this Cell was already created, use it
    if cell_id in OPENMOC_CELLS:
        return OPENMOC_CELLS[cell_id]

    # Create an OpenMOC Cell to represent this OpenCG Cell
    name = str(opencg_cell.name)
    openmoc_cell = openmoc.Cell(cell_id, name)

    fill = opencg_cell.fill
    if opencg_cell.type == 'universe':
        openmoc_cell.setFill(get_openmoc_universe(fill))
    elif opencg_cell.type == 'lattice':
        openmoc_cell.setFill(get_openmoc_lattice(fill))
    else:
        openmoc_cell.setFill(get_openmoc_material(fill))

    if opencg_cell.rotation is not None:
        rotation = np.asarray(opencg_cell.rotation, dtype=np.float64)
        openmoc_cell.setRotation(rotation)
    if opencg_cell.translation is not None:
        translation = np.asarray(opencg_cell.translation, dtype=np.float64)
        openmoc_cell.setTranslation(translation)

    surfaces = opencg_cell.surfaces

    for surface_id in surfaces:
        surface = surfaces[surface_id][0]
        halfspace = int(surfaces[surface_id][1])
        openmoc_cell.addSurface(halfspace, get_openmoc_surface(surface))

    # Add the OpenMOC Cell to the global collection of all OpenMOC Cells
    OPENMOC_CELLS[cell_id] = openmoc_cell

    # Add the OpenCG Cell to the global collection of all OpenCG Cells
    OPENCG_CELLS[cell_id] = opencg_cell

    return openmoc_cell


def get_opencg_universe(openmoc_universe):
    """Return an OpenCG universe corresponding to an OpenMOC universe.

    Parameters
    ----------
    openmoc_universe : openmoc.Universe
        OpenMOC universe

    Returns
    -------
    opencg_universe : opencg.Universe
        Equivalent OpenCG universe

    """

    cv.check_type('openmoc_universe', openmoc_universe, openmoc.Universe)

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
        opencg_universe.add_cell(opencg_cell)

    # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
    OPENMOC_UNIVERSES[universe_id] = openmoc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return opencg_universe


def get_openmoc_universe(opencg_universe):
    """Return an OpenMOC universe corresponding to an OpenCG universe.

    Parameters
    ----------
    opencg_universe : opencg.Universe
        OpenCG universe

    Returns
    -------
    openmoc_universe : openmoc.Universe
        Equivalent OpenMOC universe

    """

    cv.check_type('opencg_universe', opencg_universe, opencg.Universe)

    global OPENMOC_UNIVERSES
    universe_id = opencg_universe.id

    # If this Universe was already created, use it
    if universe_id in OPENMOC_UNIVERSES:
        return OPENMOC_UNIVERSES[universe_id]

    # Make all OpenCG Cells and Surfaces in this Universe compatible with OpenMOC
    make_opencg_cells_compatible(opencg_universe)

    # Create an OpenMOC Universe to represent this OpenCG Universe
    name = str(opencg_universe.name)
    openmoc_universe = openmoc.Universe(universe_id, name)

    # Convert all OpenCG Cells in this Universe to OpenMOC Cells
    opencg_cells = opencg_universe.cells

    for cell_id, opencg_cell in opencg_cells.items():
        openmoc_cell = get_openmoc_cell(opencg_cell)
        openmoc_universe.addCell(openmoc_cell)

    # Add the OpenMOC Universe to the global collection of all OpenMOC Universes
    OPENMOC_UNIVERSES[universe_id] = openmoc_universe

    # Add the OpenCG Universe to the global collection of all OpenCG Universes
    OPENCG_UNIVERSES[universe_id] = opencg_universe

    return openmoc_universe


def get_opencg_lattice(openmoc_lattice):
    """Return an OpenCG lattice corresponding to an OpenMOC lattice.

    Parameters
    ----------
    openmoc_lattice : openmoc.Lattice
        OpenMOC lattice

    Returns
    -------
    opencg_lattice : opencg.Lattice
        Equivalent OpenCG lattice

    """

    cv.check_type('openmoc_lattice', openmoc_lattice, openmoc.Lattice)

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
    opencg_lattice.dimension = dimension
    opencg_lattice.width = width
    opencg_lattice.universes = universe_array

    offset = np.array(lower_left, dtype=np.float64) - \
             ((np.array(width, dtype=np.float64) * \
               np.array(dimension, dtype=np.float64))) / -2.0
    opencg_lattice.offset = offset

    # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
    OPENMOC_LATTICES[lattice_id] = openmoc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return opencg_lattice


def get_openmoc_lattice(opencg_lattice):
    """Return an OpenMOC lattice corresponding to an OpenCG lattice.

    Parameters
    ----------
    opencg_lattice : opencg.Lattice
        OpenCG lattice

    Returns
    -------
    openmoc_lattice : openmoc.Lattice
        Equivalent OpenMOC lattice

    """

    cv.check_type('opencg_lattice', opencg_lattice, opencg.Lattice)

    global OPENMOC_LATTICES
    lattice_id = opencg_lattice.id

    # If this Lattice was already created, use it
    if lattice_id in OPENMOC_LATTICES:
        return OPENMOC_LATTICES[lattice_id]

    name = str(opencg_lattice.name)
    dimension = opencg_lattice.dimension
    width = opencg_lattice.width
    offset = opencg_lattice.offset
    universes = opencg_lattice.universes

    # Initialize an empty array for the OpenMOC nested Universes in this Lattice
    universe_array = np.ndarray(tuple(dimension[::-1]), dtype=openmoc.Universe)

    # Create OpenMOC Universes for each unique nested Universe in this Lattice
    unique_universes = opencg_lattice.get_unique_universes()

    for universe_id, universe in unique_universes.items():
        unique_universes[universe_id] = get_openmoc_universe(universe)

    # Build the nested Universe array
    for z in range(dimension[2]):
        for y in range(dimension[1]):
            for x in range(dimension[0]):
                universe_id = universes[z][y][x].id
                universe_array[z][dimension[1]-y-1][x] = unique_universes[universe_id]

    openmoc_lattice = openmoc.Lattice(lattice_id, name)
    openmoc_lattice.setWidth(width[0], width[1], width[2])
    openmoc_lattice.setUniverses(universe_array.tolist())
    openmoc_lattice.setOffset(offset[0], offset[1], offset[2])

    # Add the OpenMOC Lattice to the global collection of all OpenMOC Lattices
    OPENMOC_LATTICES[lattice_id] = openmoc_lattice

    # Add the OpenCG Lattice to the global collection of all OpenCG Lattices
    OPENCG_LATTICES[lattice_id] = opencg_lattice

    return openmoc_lattice


def get_opencg_geometry(openmoc_geometry):
    """Return an OpenCG geometry corresponding to an OpenMOC geometry.

    Parameters
    ----------
    openmoc_geometry : openmoc.Geometry
        OpenMOC geometry

    Returns
    -------
    opencg_geometry : opencg.Geometry
        Equivalent OpenCG geometry

    """

    cv.check_type('openmoc_geometry', openmoc_geometry, openmoc.Geometry)

    # Clear dictionaries and auto-generated IDs
    OPENMOC_MATERIALS.clear()
    OPENCG_MATERIALS.clear()
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
    opencg_geometry.root_universe = opencg_root_universe
    opencg_geometry.initialize_cell_offsets()

    return opencg_geometry


def get_openmoc_geometry(opencg_geometry):
    """Return an OpenMOC geometry corresponding to an OpenCG geometry.

    Parameters
    ----------
    opencg_geometry : opencg.Geometry
        OpenCG geometry

    Returns
    -------
    openmoc_geometry : openmoc.Geometry
        Equivalent OpenMOC geometry

    """

    cv.check_type('opencg_geometry', opencg_geometry, opencg.Geometry)

    # Deep copy the goemetry since it may be modified to make all Surfaces
    # compatible with OpenMOC's specifications
    opencg_geometry.assign_auto_ids()
    opencg_geometry = copy.deepcopy(opencg_geometry)

    # Update Cell bounding boxes in Geometry
    opencg_geometry.update_bounding_boxes()

    # Clear dictionaries and auto-generated IDs
    OPENMOC_MATERIALS.clear()
    OPENCG_MATERIALS.clear()
    OPENMOC_SURFACES.clear()
    OPENCG_SURFACES.clear()
    OPENMOC_CELLS.clear()
    OPENCG_CELLS.clear()
    OPENMOC_UNIVERSES.clear()
    OPENCG_UNIVERSES.clear()
    OPENMOC_LATTICES.clear()
    OPENCG_LATTICES.clear()

    # Make the entire geometry "compatible" before assigning auto IDs
    universes = opencg_geometry.get_all_universes()
    for universe_id, universe in universes.items():
        make_opencg_cells_compatible(universe)

    opencg_geometry.assign_auto_ids()

    opencg_root_universe = opencg_geometry.root_universe
    openmoc_root_universe = get_openmoc_universe(opencg_root_universe)

    openmoc_geometry = openmoc.Geometry()
    openmoc_geometry.setRootUniverse(openmoc_root_universe)

    # Update OpenMOC's auto-generated object IDs (e.g., Surface, Material)
    # with the maximum of those created from the OpenCG objects
    all_materials = openmoc_geometry.getAllMaterials()
    all_surfaces = openmoc_geometry.getAllSurfaces()
    all_cells = openmoc_geometry.getAllCells()
    all_universes = openmoc_geometry.getAllUniverses()

    max_material_id = max(all_materials.keys())
    max_surface_id = max(all_surfaces.keys())
    max_cell_id = max(all_cells.keys())
    max_universe_id = max(all_universes.keys())

    openmoc.maximize_material_id(max_material_id+1)
    openmoc.maximize_surface_id(max_surface_id+1)
    openmoc.maximize_cell_id(max_cell_id+1)
    openmoc.maximize_universe_id(max_universe_id+1)

    return openmoc_geometry

import sys
from abc import ABCMeta, abstractmethod

import numpy

sys.path.insert(0, 'openmoc')
import openmoc


class InputSet(object):
    """An abstract class for defining OpenMOC test geometries."""

    __metaclass__ = ABCMeta

    def __init__(self):
        self.materials = {}
        self.geometry = None

    @abstractmethod
    def create_materials(self):
        """Instantiate Materials."""
        return

    @abstractmethod
    def create_geometry(self):
        """Instantiate a Geometry."""
        return


class HomInfMedInput(InputSet):
    """A homogeneous infinite medium problem with 2-group cross sections."""

    def create_materials(self):

        sigma_f = numpy.array([0.000625, 0.135416667])
        nu_sigma_f = numpy.array([0.0015, 0.325])
        sigma_s = numpy.array([[0.1, 0.117], [0., 1.42]])
        chi = numpy.array([1.0, 0.0])
        sigma_t = numpy.array([0.2208, 1.604])

        self.materials['infinite medium'] = openmoc.Material()
        self.materials['infinite medium'].setName('2-group infinite medium')
        self.materials['infinite medium'].setNumEnergyGroups(2)
        self.materials['infinite medium'].setSigmaF(sigma_f)
        self.materials['infinite medium'].setNuSigmaF(nu_sigma_f)
        self.materials['infinite medium'].setSigmaS(sigma_s.flat)
        self.materials['infinite medium'].setChi(chi)
        self.materials['infinite medium'].setSigmaT(sigma_t)

    def create_geometry(self):
        """Instantiate an infinite medium lattice Geometry."""

        length = 2.5
        num_cells_x = 10
        num_cells_y = 10

        xmin = openmoc.XPlane(x=-length/2., name='xmin')
        xmax = openmoc.XPlane(x=+length/2., name='xmax')
        ymin = openmoc.YPlane(y=-length/2., name='ymin')
        ymax = openmoc.YPlane(y=+length/2., name='ymax')

        xmax.setBoundaryType(openmoc.REFLECTIVE)
        xmin.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        fill = openmoc.Cell(name='fill')
        fill.setFill(self.materials['infinite medium'])

        root_cell = openmoc.Cell(name='root cell')
        root_cell.addSurface(halfspace=+1, surface=xmin)
        root_cell.addSurface(halfspace=-1, surface=xmax)
        root_cell.addSurface(halfspace=+1, surface=ymin)
        root_cell.addSurface(halfspace=-1, surface=ymax)

        fill_universe = openmoc.Universe(name='homogeneous fill cell')
        fill_universe.addCell(fill)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(root_cell)

        lattice = openmoc.Lattice(name='MxN lattice')
        lattice.setWidth(width_x=length/num_cells_x, width_y=length/num_cells_y)
        lattice.setUniverses([[[fill_universe] * num_cells_x]*num_cells_y])
        root_cell.setFill(lattice)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(HomInfMedInput, self).create_geometry()


class PinCellInput(InputSet):
    """A pin cell problem from sample-input/pin-cell."""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a pin cell Geometry."""

        zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymin = openmoc.YPlane(y=-2.0, name='ymin')
        ymax = openmoc.YPlane(y=+2.0, name='ymax')

        xmin.setBoundaryType(openmoc.REFLECTIVE)
        xmax.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        fuel = openmoc.Cell(name='fuel')
        fuel.setFill(self.materials['UO2'])
        fuel.addSurface(halfspace=-1, surface=zcylinder)

        moderator = openmoc.Cell(name='moderator')
        moderator.setFill(self.materials['Water'])
        moderator.addSurface(halfspace=+1, surface=zcylinder)
        moderator.addSurface(halfspace=+1, surface=xmin)
        moderator.addSurface(halfspace=-1, surface=xmax)
        moderator.addSurface(halfspace=+1, surface=ymin)
        moderator.addSurface(halfspace=-1, surface=ymax)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(fuel)
        root_universe.addCell(moderator)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(PinCellInput, self).create_geometry()


class GridInput(InputSet):
    """Sets up a simple grid problem using explicit planes"""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a 3x3 grid Geometry."""

        # Create the planes that bound the cell and geometry
        xplanes = [None] * 4
        yplanes = [None] * 4
        for i in range(4):
            val = -3.0 + 6.0/3 * i
            xplanes[i] = openmoc.XPlane(x=val, name='xplane' + str(i))
            yplanes[i] = openmoc.YPlane(y=val, name='yplane' + str(i))

        xplanes[0].setBoundaryType(openmoc.REFLECTIVE)
        xplanes[-1].setBoundaryType(openmoc.REFLECTIVE)
        yplanes[0].setBoundaryType(openmoc.REFLECTIVE)
        yplanes[-1].setBoundaryType(openmoc.REFLECTIVE)

        # Create the cells and set the central cell to UO2, the rest as water
        cells = [[None]*3 for i in range(3)]
        for i in range(3):
            for j in range(3):
                cells[i][j] = openmoc.Cell(name='Cell ' + str(3*i+j))
                if i == 1 and j == 1:
                    cells[i][j].setFill(self.materials['UO2'])
                else:
                    cells[i][j].setFill(self.materials['Water'])
                cells[i][j].addSurface(halfspace=+1, surface=xplanes[i])
                cells[i][j].addSurface(halfspace=-1, surface=xplanes[i+1])
                cells[i][j].addSurface(halfspace=+1, surface=yplanes[j])
                cells[i][j].addSurface(halfspace=-1, surface=yplanes[j+1])

        # Add the cells to the universe
        root_universe = openmoc.Universe(name='root universe')
        for i in range(3):
            for j in range(3):
                root_universe.addCell(cells[i][j])

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(GridInput, self).create_geometry()


class LatticeGridInput(InputSet):
    """A simple lattice problem equivalent to the grid case"""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a 3x3 grid Geometry by making a lattice"""

        # Create the planes bounding the geometry
        xmin = openmoc.XPlane(x=-3.0, name='xmin')
        xmax = openmoc.XPlane(x=3.0, name='xmax')
        ymin = openmoc.YPlane(y=-3.0, name='ymin')
        ymax = openmoc.YPlane(y=3.0, name='ymax')

        xmin.setBoundaryType(openmoc.REFLECTIVE)
        xmax.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        # Create the root cell, bounded by the geometry bounds
        root_cell = openmoc.Cell(name='root cell')
        root_cell.addSurface(halfspace=+1, surface=xmin)
        root_cell.addSurface(halfspace=-1, surface=xmax)
        root_cell.addSurface(halfspace=+1, surface=ymin)
        root_cell.addSurface(halfspace=-1, surface=ymax)

        # Create UO2 and water cells
        uo2_cell = openmoc.Cell(name='UO2 Cell')
        uo2_cell.setFill(self.materials['UO2'])
        water_cell = openmoc.Cell(name='Water Cell')
        water_cell.setFill(self.materials['Water'])

        # Create universes and fill with the associated cell
        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(root_cell)
        uo2 = openmoc.Universe(name='uo2 universe')
        uo2.addCell(uo2_cell)
        water = openmoc.Universe(name='uo2 universe')
        water.addCell(water_cell)

        # Create the lattice and fill it with the appropriate universes
        lattice = openmoc.Lattice(name='3x3 lattice')
        lattice.setWidth(width_x=2.0, width_y=2.0)
        lattice.setUniverses([[[water, water, water],
                               [water, uo2,   water],
                               [water, water, water]]])
        root_cell.setFill(lattice)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(LatticeGridInput, self).create_geometry()



class SimpleLatticeInput(InputSet):
    """A 4x4 pin cell lattice problem from sample-input/simple-lattice."""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a 4x4 pin cell lattice Geometry."""

        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymax = openmoc.YPlane(y=+2.0, name='ymin')
        ymin = openmoc.YPlane(y=-2.0, name='ymax')
        boundaries = [xmin, xmax, ymin, ymax]

        large_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                            radius=0.4, name='large pin')
        medium_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                             radius=0.3, name='medium pin')
        small_zcylinder = openmoc.ZCylinder(x=0.0, y=0.0,
                                            radius=0.2, name='small pin')

        for boundary in boundaries: boundary.setBoundaryType(openmoc.REFLECTIVE)

        large_fuel = openmoc.Cell(name='large pin fuel')
        large_fuel.setNumRings(3)
        large_fuel.setNumSectors(8)
        large_fuel.setFill(self.materials['UO2'])
        large_fuel.addSurface(halfspace=-1, surface=large_zcylinder)

        large_moderator = openmoc.Cell(name='large pin moderator')
        large_moderator.setNumSectors(8)
        large_moderator.setFill(self.materials['Water'])
        large_moderator.addSurface(halfspace=+1, surface=large_zcylinder)

        medium_fuel = openmoc.Cell(name='medium pin fuel')
        medium_fuel.setNumRings(3)
        medium_fuel.setNumSectors(8)
        medium_fuel.setFill(self.materials['UO2'])
        medium_fuel.addSurface(halfspace=-1, surface=medium_zcylinder)

        medium_moderator = openmoc.Cell(name='medium pin moderator')
        medium_moderator.setNumSectors(8)
        medium_moderator.setFill(self.materials['Water'])
        medium_moderator.addSurface(halfspace=+1, surface=medium_zcylinder)

        small_fuel = openmoc.Cell(name='small pin fuel')
        small_fuel.setNumRings(3)
        small_fuel.setNumSectors(8)
        small_fuel.setFill(self.materials['UO2'])
        small_fuel.addSurface(halfspace=-1, surface=small_zcylinder)

        small_moderator = openmoc.Cell(name='small pin moderator')
        small_moderator.setNumSectors(8)
        small_moderator.setFill(self.materials['Water'])
        small_moderator.addSurface(halfspace=+1, surface=small_zcylinder)

        lattice_cell = openmoc.Cell(name='lattice cell')

        root_cell = openmoc.Cell(name='root cell')
        root_cell.addSurface(halfspace=+1, surface=boundaries[0])
        root_cell.addSurface(halfspace=-1, surface=boundaries[1])
        root_cell.addSurface(halfspace=+1, surface=boundaries[2])
        root_cell.addSurface(halfspace=-1, surface=boundaries[3])

        pin1 = openmoc.Universe(name='large pin cell')
        pin2 = openmoc.Universe(name='medium pin cell')
        pin3 = openmoc.Universe(name='small pin cell')
        assembly = openmoc.Universe(name='2x2 lattice')
        root_universe = openmoc.Universe(name='root universe')

        pin1.addCell(large_fuel)
        pin1.addCell(large_moderator)
        pin2.addCell(medium_fuel)
        pin2.addCell(medium_moderator)
        pin3.addCell(small_fuel)
        pin3.addCell(small_moderator)
        assembly.addCell(lattice_cell)
        root_universe.addCell(root_cell)

        # 2x2 assembly
        lattice = openmoc.Lattice(name='2x2 lattice')
        lattice.setWidth(width_x=1.0, width_y=1.0)
        lattice.setUniverses([[[pin1, pin2], [pin1, pin3]]])
        lattice_cell.setFill(lattice)

        # 2x2 core
        core = openmoc.Lattice(name='2x2 core')
        core.setWidth(width_x=2.0, width_y=2.0)
        core.setUniverses([[[assembly, assembly], [assembly, assembly]]])
        root_cell.setFill(core)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(SimpleLatticeInput, self).create_geometry()


class PwrAssemblyInput(InputSet):
    """A 17x17 pin cell lattice problem from sample-input/ipython-notebook."""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a 17x17 pin cell lattice Geometry."""

        # Create ZCylinder for the fuel and moderator
        fuel_radius = openmoc.ZCylinder(x=0.0, y=0.0, radius=0.54)

        # Create planes to bound the entire geometry
        xmin = openmoc.XPlane(x=-10.71, name='xmin')
        xmax = openmoc.XPlane(x=+10.71, name='xmax')
        ymin = openmoc.YPlane(y=-10.71, name='ymin')
        ymax = openmoc.YPlane(y=+10.71, name='xmax')

        xmin.setBoundaryType(openmoc.REFLECTIVE)
        xmax.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        # 4.3% MOX pin cell
        mox43_cell = openmoc.Cell()
        mox43_cell.setFill(self.materials['MOX-4.3%'])
        mox43_cell.setNumRings(3)
        mox43_cell.setNumSectors(8)
        mox43_cell.addSurface(-1, fuel_radius)

        mox43 = openmoc.Universe(name='MOX-4.3%')
        mox43.addCell(mox43_cell)

        # 7% MOX pin cell
        mox7_cell = openmoc.Cell()
        mox7_cell.setFill(self.materials['MOX-7%'])
        mox7_cell.setNumRings(3)
        mox7_cell.setNumSectors(8)
        mox7_cell.addSurface(-1, fuel_radius)

        mox7 = openmoc.Universe(name='MOX-7%')
        mox7.addCell(mox7_cell)

        # 8.7% MOX pin cell
        mox87_cell = openmoc.Cell()
        mox87_cell.setFill(self.materials['MOX-8.7%'])
        mox87_cell.setNumRings(3)
        mox87_cell.setNumSectors(8)
        mox87_cell.addSurface(-1, fuel_radius)

        mox87 = openmoc.Universe(name='MOX-8.7%')
        mox87.addCell(mox87_cell)

        # Fission chamber pin cell
        fission_chamber_cell = openmoc.Cell()
        fission_chamber_cell.setFill(self.materials['Fission Chamber'])
        fission_chamber_cell.setNumRings(3)
        fission_chamber_cell.setNumSectors(8)
        fission_chamber_cell.addSurface(-1, fuel_radius)

        fission_chamber = openmoc.Universe(name='Fission Chamber')
        fission_chamber.addCell(fission_chamber_cell)

        # Guide tube pin cell
        guide_tube_cell = openmoc.Cell()
        guide_tube_cell.setFill(self.materials['Guide Tube'])
        guide_tube_cell.setNumRings(3)
        guide_tube_cell.setNumSectors(8)
        guide_tube_cell.addSurface(-1, fuel_radius)

        guide_tube = openmoc.Universe(name='Guide Tube')
        guide_tube.addCell(guide_tube_cell)

        # Moderator cell
        moderator = openmoc.Cell()
        moderator.setFill(self.materials['Water'])
        moderator.addSurface(+1, fuel_radius)
        moderator.setNumRings(3)
        moderator.setNumSectors(8)

        # Add moderator to each pin cell
        pins = [mox43, mox7, mox87, fission_chamber, guide_tube]
        for pin in pins:
            pin.addCell(moderator)

        # CellFills for the assembly
        assembly1_cell = openmoc.Cell(name='Assembly 1')
        assembly1 = openmoc.Universe(name='Assembly 1')
        assembly1.addCell(assembly1_cell)

        # A mixed enrichment PWR MOX fuel assembly
        assembly = openmoc.Lattice(name='MOX Assembly')
        assembly.setWidth(width_x=1.26, width_y=1.26)

        # Create a template to map to pin cell types
        template = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                    [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
                    [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
                    [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
                    [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
                    [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
                    [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
                    [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
                    [1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1],
                    [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
                    [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
                    [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
                    [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
                    [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
                    [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
                    [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

        universes = {1 : mox43, 2 : mox7, 3 : mox87,
                     4 : guide_tube, 5 : fission_chamber}

        for i in range(17):
            for j in range(17):
                template[i][j] = universes[template[i][j]]

        assembly.setUniverses([template])

        # Root Cell/Universe
        root_cell = openmoc.Cell(name='Full Geometry')
        root_cell.setFill(assembly)
        root_cell.addSurface(+1, xmin)
        root_cell.addSurface(-1, xmax)
        root_cell.addSurface(+1, ymin)
        root_cell.addSurface(-1, ymax)

        root_universe = openmoc.Universe(name='Root Universe')
        root_universe.addCell(root_cell)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(PwrAssemblyInput, self).create_geometry()

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
        self.geometry.initializeFlatSourceRegions()
        return


class HomInfMedInput(InputSet):
    """A homogenized infinite medium problem with 2-group cross sections."""

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
        """Instantiate an infinite medium Geometry."""

        xmin = openmoc.XPlane(x=-100.0, name='left')
        xmax = openmoc.XPlane(x=+100.0, name='right')
        ymin = openmoc.YPlane(y=-100.0, name='bottom')
        ymax = openmoc.YPlane(y=+100.0, name='top')

        xmax.setBoundaryType(openmoc.REFLECTIVE)
        xmin.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        cell = openmoc.Cell()
        cell.setFill(self.materials['infinite_medium'])
        cell.addSurface(halfspace=+1, surface=xmin)
        cell.addSurface(halfspace=-1, surface=xmax)
        cell.addSurface(halfspace=+1, surface=ymin)
        cell.addSurface(halfspace=-1, surface=ymax)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(cell)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(PinCellInput, self).create_geometry()


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
        ymin = openmoc.YPlane(y=+2.0, name='ymin')
        ymax = openmoc.YPlane(y=-2.0, name='ymax')
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
        root_cell.addSurface(halfspace=-1, surface=boundaries[2])
        root_cell.addSurface(halfspace=+1, surface=boundaries[3])

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

import sys
import numpy as np
import openmoc
 #For Python 2.X.X
if (sys.version_info[0] == 2):
    import checkvalue as cv
# For Python 3.X.X
else:
    import openmoc.checkvalue as cv


class Subdivider(openmoc.Lattice):
    """A lattice that subdivides an arbitrary universe along a uniform mesh

    NOTE: Currently only 2D division is supported. Axial division is broken.

    Parameters
    ----------
    * division : iterable of int
        the number of divisions along each side of the mesh: (nx, ny[, nz]).
        lower_left and upper_right must match these dimensions, 2D or 3D
    * lower_left : iterable of float, cm
        the lower left corner of the universe to subdivide:  (xmin, ymin[, zmin])
    * upper_right : iterable of float, cm
        the upper right corner of the universe to subdivide: (xmax, ymax[, zmax])
        if not provided, will mirror lower_left around the origin.
    * id : int
        the user-specified optional Lattice (Universe) ID
    * name : str
        the user-specified optional Lattice (Universe) name

    Attributes:
    -----------
    * division
    * lower_left :  np.ndarray
    * upper_right : np.ndarray
    * deltas : iterable of float, cm
        the uniform mesh spacing in each dimension: (dx, dy[, dx])
    """
    
    def __init__(self, division, lower_left, upper_right=None, *args,
                 **kwargs):
        ndim = len(division)
        cv.check_iterable_type("division", division, int)
        cv.check_length("division", division, 2, 3)
        cv.check_length("lower_left", lower_left, ndim)
        if upper_right is not None:
            cv.check_length("upper_right", upper_right, ndim)
        
        self.lower_left = np.array(lower_left)
        if upper_right is None:
            self.upper_right = -self.lower_left
        else:
            self.upper_right = np.array(upper_right)
        self.division = division
        
        self.ndim = ndim
        super().__init__(*args, **kwargs)
    
    @property
    def deltas(self):
        return np.divide(self.upper_right - self.lower_left, self.division)
    
    def _subdivide2d(self, target):
        nx, ny = self.division
        dx, dy = self.deltas
        x0 = nx*dx/2
        y0 = -ny*dy/2
        universes = [[None for _ in range(nx)] for _ in range(ny)]
        # for 2D, no translation in z
        tz = 0
        for j in range(ny):
            # translation in y
            ty = y0 + dy*(j + .5)
            for i in range(nx):
                # translation in x
                tx = x0 - dx*(i + .5)
                t_vector = np.array([tx, ty, tz])
                cell = openmoc.Cell()
                cell.setFill(target)
                cell.setTranslation(t_vector)
                univ = openmoc.Universe()
                univ.addCell(cell)
                universes[j][i] = univ
        self.setUniverses([universes])
        return self._to_openmoc_universe(name=target.getName())
    
    def _subdivide3d(self, target):
        raise NotImplementedError("Axial subdivision.")
    
    def _to_openmoc_universe(self, name=""):
        if name:
            mesh_name = "x".join(np.array(self.division, dtype=str))
            name += " (subdivided " + mesh_name + ")"
        uout = openmoc.Universe(name=name)
        cout = openmoc.Cell()
        cout.setFill(self)
        uout.addCell(cout)
        return uout
    
    def get_subdivided_universe(self, target):
        """Subdivide some universe along the specified mesh.

        Parameter:
        ----------
        * target : openmoc.Universe
            the target universe to subdivide

        Returns:
        --------
        openmoc.Universe
            filled with a cell containing `target' discretized along
            the Subdivider
        """
        cv.check_type("target", target, openmoc.Universe)
        self.setWidth(*self.deltas)
        if self.ndim == 2:
            return self._subdivide2d(target)
        else:
            return self._subdivide3d(target)

/* C++ casting helper method for openmoc.process computePinPowers
 * routine and the OpenMC compatibility module */

%module casting

%inline %{

  Lattice* castUniverseToLattice(Universe* universe) {
    return dynamic_cast<Lattice*>(universe);
  }

  Universe* castLatticeToUniverse(Lattice* lattice) {
    return dynamic_cast<Universe*>(lattice);
  }

  Plane* castSurfaceToPlane(Surface* plane) {
    return dynamic_cast<Plane*>(plane);
  }

  XPlane* castSurfaceToXPlane(Surface* xplane) {
    return dynamic_cast<XPlane*>(xplane);
  }

  YPlane* castSurfaceToYPlane(Surface* yplane) {
    return dynamic_cast<YPlane*>(yplane);
  }

  ZPlane* castSurfaceToZPlane(Surface* zplane) {
    return dynamic_cast<ZPlane*>(zplane);
  }

  ZCylinder* castSurfaceToZCylinder(Surface* zcylinder) {
    return dynamic_cast<ZCylinder*>(zcylinder);
  }

%}

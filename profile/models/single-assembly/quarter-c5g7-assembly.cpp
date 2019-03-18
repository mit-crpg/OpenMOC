#include "../../../src/CPUSolver.h"
#include "../../../src/log.h"
#include <array>
#include <iostream>

int main() {

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_num_procs();
#else
  int num_threads = 1;
#endif
  double azim_spacing = 0.1;
  int num_azim = 4;
  double polar_spacing = 1.0;
  int num_polar = 6;
  double tolerance = 1e-5;
  int max_iters = 1000;
  int axial_refines = 1;

  /* Set logging information */
  set_log_level("NORMAL");
  log_printf(TITLE, "Simulating the OECD's C5G7 Benchmark Problem...");

  /* Define material properties */
  log_printf(NORMAL, "Defining material properties...");

  const size_t num_groups = 7;
  std::map<std::string, std::array<double, num_groups> > sigma_a;
  std::map<std::string, std::array<double, num_groups> > nu_sigma_f;
  std::map<std::string, std::array<double, num_groups> > sigma_f;
  std::map<std::string, std::array<double, num_groups*num_groups> > sigma_s;
  std::map<std::string, std::array<double, num_groups> > chi;
  std::map<std::string, std::array<double, num_groups> > sigma_t;

  /* Define water cross-sections */
  sigma_a["Water"] = std::array<double, num_groups> {6.0105E-4, 1.5793E-5,
      3.3716E-4, 0.0019406, 0.0057416, 0.015001, 0.037239};
  nu_sigma_f["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_f["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_s["Water"] = std::array<double, num_groups*num_groups>
      {0.0444777, 0.1134, 7.2347E-4, 3.7499E-6, 5.3184E-8, 0.0, 0.0,
      0.0, 0.282334, 0.12994, 6.234E-4, 4.8002E-5, 7.4486E-6, 1.0455E-6,
      0.0, 0.0, 0.345256, 0.22457, 0.016999, 0.0026443, 5.0344E-4,
      0.0, 0.0, 0.0, 0.0910284, 0.41551, 0.063732, 0.012139,
      0.0, 0.0, 0.0, 7.1437E-5, 0.139138, 0.51182, 0.061229,
      0.0, 0.0, 0.0, 0.0, 0.0022157, 0.699913, 0.53732,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.13244, 2.4807};
  chi["Water"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_t["Water"] = std::array<double, num_groups> {0.159206, 0.41297,
    0.59031, 0.58435, 0.718, 1.25445, 2.65038};

  /* Define UO2 cross-sections */
  sigma_a["UO2"] = std::array<double, num_groups> {0.0080248, 0.0037174,
    0.026769, 0.096236, 0.03002, 0.11126, 0.28278};
  nu_sigma_f["UO2"] = std::array<double, num_groups> {0.02005998, 0.002027303,
    0.01570599, 0.04518301, 0.04334208, 0.2020901, 0.5257105};
  sigma_f["UO2"] = std::array<double, num_groups> {0.00721206, 8.19301E-4,
    0.0064532, 0.0185648, 0.0178084, 0.0830348, 0.216004};
  sigma_s["UO2"] = std::array<double, num_groups*num_groups>
      {0.127537, 0.042378, 9.4374E-6, 5.5163E-9, 0.0, 0.0, 0.0,
      0.0, 0.324456, 0.0016314, 3.1427E-9, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.45094, 0.0026792, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.452565, 0.0055664, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.2525E-4, 0.271401, 0.010255, 1.0021E-8,
      0.0, 0.0, 0.0, 0.0, 0.0012968, 0.265802, 0.016809,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0085458, 0.27308};
  chi["UO2"] = std::array<double, num_groups> {0.58791, 0.41176, 3.3906E-4,
    1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["UO2"] = std::array<double, num_groups> {0.177949, 0.329805,
    0.480388, 0.554367, 0.311801, 0.395168, 0.564406};

  /* Define fission chamber cross-sections */
  sigma_a["Fission Chamber"] = std::array<double, num_groups> {5.1132E-4,
    7.5813E-5, 3.1643E-4, 0.0011675, 0.0033977, 0.0091886, 0.023244};
  nu_sigma_f["Fission Chamber"] = std::array<double, num_groups> {1.323401E-8,
    1.4345E-8, 1.128599E-6, 1.276299E-5, 3.538502E-7, 1.740099E-6,
    5.063302E-6};
  sigma_f["Fission Chamber"] = std::array<double, num_groups> {4.79002E-9,
    5.82564E-9, 4.63719E-7, 5.24406E-6, 1.4539E-7, 7.14972E-7, 2.08041E-6};
  sigma_s["Fission Chamber"] = std::array<double, num_groups*num_groups>
      {0.0661659, 0.05907, 2.8334E-4, 1.4622E-6, 2.0642E-8, 0.0, 0.0,
      0.0, 0.240377, 0.052435, 2.499E-4, 1.9239E-5, 2.9875E-6, 4.214E-7,
      0.0, 0.0, 0.183425, 0.092288, 0.0069365, 0.001079, 2.0543E-4,
      0.0, 0.0, 0.0, 0.0790769, 0.16999, 0.02586, 0.0049256,
      0.0, 0.0, 0.0, 3.734E-5, 0.099757, 0.20679, 0.024478,
      0.0, 0.0, 0.0, 0.0, 9.1742E-4, 0.316774, 0.23876,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.049793, 1.0991};
  chi["Fission Chamber"] = std::array<double, num_groups> {0.58791, 0.41176,
    3.3906E-4, 1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["Fission Chamber"] = std::array<double, num_groups> {0.126032,
    0.29316, 0.28425, 0.28102, 0.33446, 0.56564, 1.17214};

  /* Define guide tube cross-sections */
  sigma_a["Guide Tube"] = std::array<double, num_groups> {5.1132E-4, 7.5801E-5,
    3.1572E-4, 0.0011582, 0.0033975, 0.0091878, 0.023242};
  nu_sigma_f["Guide Tube"] = std::array<double, num_groups> {0, 0, 0, 0, 0,
    0, 0};
  sigma_f["Guide Tube"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_s["Guide Tube"] = std::array<double, num_groups*num_groups>
      {0.0661659, 0.05907, 2.8334E-4, 1.4622E-6, 2.0642E-8, 0.0, 0.0,
      0.0, 0.240377, 0.052435, 2.499E-4, 1.9239E-5, 2.9875E-6, 4.214E-7,
      0.0, 0.0, 0.183297, 0.092397, 0.0069446, 0.0010803, 2.0567E-4,
      0.0, 0.0, 0.0, 0.0788511, 0.17014, 0.025881, 0.0049297,
      0.0, 0.0, 0.0, 3.7333E-5, 0.0997372, 0.20679, 0.024478,
      0.0, 0.0, 0.0, 0.0, 9.1726E-4, 0.316765, 0.23877,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.049792, 1.09912};
  chi["Guide Tube"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_t["Guide Tube"] = std::array<double, num_groups> {0.126032, 0.29316,
    0.28424, 0.28096, 0.33444, 0.56564, 1.17215};

  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  std::map<std::string, Material*> materials;

  std::map<std::string, std::array<double, num_groups> >::iterator it;
  int id_num = 0;
  for (it = sigma_t.begin(); it != sigma_t.end(); it++) {

    std::string name = it->first;
    materials[name] = new Material(id_num, name.c_str());
    materials[name]->setNumEnergyGroups(num_groups);
    id_num++;

    materials[name]->setSigmaF(sigma_f[name].data(), num_groups);
    materials[name]->setNuSigmaF(nu_sigma_f[name].data(), num_groups);
    materials[name]->setSigmaS(sigma_s[name].data(), num_groups*num_groups);
    materials[name]->setChi(chi[name].data(), num_groups);
    materials[name]->setSigmaT(sigma_t[name].data(), num_groups);
  }

  /* Create surfaces */
  XPlane xmin(-5.355);
  XPlane xmax( 5.355);
  YPlane ymin(-5.355);
  YPlane ymax( 5.355);
  ZPlane zmin(-5.355);
  ZPlane zmax( 5.355);

  xmin.setBoundaryType(REFLECTIVE);
  xmax.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  ymax.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(REFLECTIVE);
  zmax.setBoundaryType(REFLECTIVE);

  /* Create z-cylinders for the fuel as well as to discretize the moderator
   * into rings */
  ZCylinder fuel_radius(0.0, 0.0, 0.54);
  ZCylinder moderator_inner_radius(0.0, 0.0, 0.58);
  ZCylinder moderator_outer_radius(0.0, 0.0, 0.62);

  /* Create cells and universes */
  log_printf(NORMAL, "Creating cells...");

  /* Moderator rings */
  Cell* moderator = new Cell();
  moderator->setFill(materials["Water"]);
  moderator->addSurface(+1, &fuel_radius);
  moderator->setNumRings(2);
  moderator->setNumSectors(4);

  /* UO2 pin cell */
  Cell* uo2_cell = new Cell(3, "uo2");
  uo2_cell->setNumRings(2);
  uo2_cell->setNumSectors(4);
  uo2_cell->setFill(materials["UO2"]);
  uo2_cell->addSurface(-1, &fuel_radius);

  Universe* uo2 = new Universe();
  uo2->addCell(uo2_cell);
  uo2->addCell(moderator);

  /* Fission chamber pin cell */
  Cell* fission_chamber_cell = new Cell(7, "fc");
  fission_chamber_cell->setNumRings(1);
  fission_chamber_cell->setNumSectors(4);
  fission_chamber_cell->setFill(materials["Fission Chamber"]);
  fission_chamber_cell->addSurface(-1, &fuel_radius);

  Universe* fission_chamber = new Universe();
  fission_chamber->addCell(fission_chamber_cell);
  fission_chamber->addCell(moderator);

  /* Guide tube pin cell */
  Cell* guide_tube_cell = new Cell(8, "gtc");
  guide_tube_cell->setNumRings(1);
  guide_tube_cell->setNumSectors(4);
  guide_tube_cell->setFill(materials["Guide Tube"]);
  guide_tube_cell->addSurface(-1, &fuel_radius);

  Universe* guide_tube = new Universe();
  guide_tube->addCell(guide_tube_cell);
  guide_tube->addCell(moderator);

  /* Cells */
  Cell* assembly1_cell = new Cell(10, "ac1");
  Universe* assembly1 = new Universe();
  assembly1->addCell(assembly1_cell);

  /* Root Cell* */
  Cell* root_cell = new Cell(16, "root");
  root_cell->addSurface(+1, &xmin);
  root_cell->addSurface(-1, &xmax);
  root_cell->addSurface(+1, &ymin);
  root_cell->addSurface(-1, &ymax);
  root_cell->addSurface(+1, &zmin);
  root_cell->addSurface(-1, &zmax);

  Universe* root_universe = new Universe();
  root_universe->addCell(root_cell);

  /* Create lattices */
  log_printf(NORMAL, "Creating lattices...");

  /* Top left, bottom right 17 x 17 assemblies */
  Lattice* assembly1_lattice = new Lattice();
  assembly1_lattice->setWidth(1.26, 1.26, 10.71/axial_refines);
  Universe* matrix1[9*9*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2,
                        1, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2,
                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 3};

    std::map<int, Universe*> names = {{1, uo2}, {2, guide_tube},
                                      {3, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<9*9; n++)
        matrix1[z*9*9 + n] = names[mold[n]];

    assembly1_lattice->setUniverses(axial_refines, 9, 9, matrix1);
  }
  assembly1_cell->setFill(assembly1_lattice);

  /* Fill root cell with lattice */
  root_cell->setFill(assembly1);

  /* Create CMFD mesh */
  log_printf(NORMAL, "Creating CMFD mesh...");

  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);
  cmfd->setLatticeStructure(9, 9, 3);
  std::vector<std::vector<int> > cmfd_group_structure;
  cmfd_group_structure.resize(7);
  for (int g=0; g<7; g++)
    cmfd_group_structure.at(g).push_back(g+1);
  cmfd->setGroupStructure(cmfd_group_structure);
  cmfd->setCentroidUpdateOn(false);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.setRootUniverse(root_universe);
  geometry.setCmfd(cmfd);
  geometry.initializeFlatSourceRegions();

  /* Generate tracks */
  log_printf(NORMAL, "Initializing the track generator...");
  Quadrature* quad = new EqualAnglePolarQuad();
  quad->setNumPolarAngles(num_polar);
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setQuadrature(quad);
  track_generator.setSegmentFormation(OTF_STACKS);
  std::vector<double> seg_heights {0.0};
  track_generator.setSegmentationZones(seg_heights);
  track_generator.generateTracks();

  /* Run simulation */
  CPUSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  //solver.setOTFTransport();
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  log_printf(TITLE, "Finished");
  return 0;
}

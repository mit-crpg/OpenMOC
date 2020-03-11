#include "../../../src/CPUSolver.h"
#include "../../../src/CPULSSolver.h"
#include "../../../src/log.h"
#include <array>
#include <iostream>

// FIXME This input file fails at 2D track generation time

int main(int argc, char* argv[]) {

#ifdef MPIx
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  log_set_ranks(MPI_COMM_WORLD);
#endif

  /* Define simulation parameters */
#ifdef OPENMP
  int num_threads = omp_get_max_threads();
#else
  int num_threads = 1;
#endif
  double azim_spacing = 0.1;
  int num_azim = 16;
  double polar_spacing = 1.5;
  int num_polar = 12;
  double tolerance = 1e-5;
  int max_iters = 50;
  int axial_refines = 1;

  /* Domain decomposition / modular ray tracing */
  int nx = 1;
  int ny = 1;
  int nz = 1;
#ifdef MPIx
  num_threads = std::max(1, num_threads / (nx*ny*nz));
#endif

  /* Set logging information */
  set_log_level("NORMAL");
  set_line_length(120);
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

  /* Define MOX-4.3% cross-sections */
  sigma_a["MOX-4.3%%"] = std::array<double, num_groups> {0.0084339, 0.0037577,
    0.02797, 0.10421, 0.13994, 0.40918, 0.40935};
  nu_sigma_f["MOX-4.3%%"] = std::array<double, num_groups> {0.021753,
    0.002535103, 0.01626799, 0.0654741, 0.03072409, 0.666651, 0.7139904};
  sigma_f["MOX-4.3%%"] = std::array<double, num_groups> {0.00762704,
    8.76898E-4, 0.00569835, 0.0228872, 0.0107635, 0.232757, 0.248968};
  sigma_s["MOX-4.3%%"] = std::array<double, num_groups*num_groups>
      {0.128876, 0.041413, 8.229E-6, 5.0405E-9, 0.0, 0.0, 0.0,
      0.0, 0.325452, 0.0016395, 1.5982E-9, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.453188, 0.0026142, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.457173, 0.0055394, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.6046E-4, 0.276814, 0.0093127, 9.1656E-9,
      0.0, 0.0, 0.0, 0.0, 0.0020051, 0.252962, 0.01485,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0084948, 0.265007};
  chi["MOX-4.3%%"] = std::array<double, num_groups> {0.58791, 0.41176,
    3.3906E-4, 1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["MOX-4.3%%"] = std::array<double, num_groups> {0.178731, 0.330849,
    0.483772, 0.566922, 0.426227, 0.678997, 0.68285};

  /* Define MOX-7% cross-sections */
  sigma_a["MOX-7%%"] = std::array<double, num_groups> {0.0090657, 0.0042967,
    0.032881, 0.12203, 0.18298, 0.56846, 0.58521};
  nu_sigma_f["MOX-7%%"] = std::array<double, num_groups> {0.02381395,
    0.003858689, 0.024134, 0.09436622, 0.04576988, 0.9281814, 1.0432};
  sigma_f["MOX-7%%"] = std::array<double, num_groups> {0.00825446, 0.00132565,
    0.00842156, 0.032873, 0.0159636, 0.323794, 0.362803};
  sigma_s["MOX-7%%"] = std::array<double, num_groups*num_groups>
      {0.130457, 0.041792, 8.5105E-6, 5.1329E-9, 0.0, 0.0, 0.0,
      0.0, 0.328428, 0.0016436, 2.2017E-9, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.458371, 0.0025331, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.463709, 0.0054766, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.7619E-4, 0.282313, 0.0087289, 9.0016E-9,
      0.0, 0.0, 0.0, 0.0, 0.002276, 0.249751, 0.013114,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0088645, 0.259529};
  chi["MOX-7%%"] = std::array<double, num_groups> {0.58791, 0.41176, 3.3906E-4,
    1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["MOX-7%%"] = std::array<double, num_groups> {0.181323, 0.334368,
    0.493785, 0.591216, 0.474198, 0.833601, 0.853603};

  /* Define MOX-8.7% cross-sections */
  sigma_a["MOX-8.7%%"] = std::array<double, num_groups> {0.0094862, 0.0046556,
    0.03624, 0.13272, 0.2084, 0.6587, 0.69017};
  nu_sigma_f["MOX-8.7%%"] = std::array<double, num_groups> {0.025186,
    0.004739509, 0.02947805, 0.11225, 0.05530301, 1.074999, 1.239298};
  sigma_f["MOX-8.7%%"] = std::array<double, num_groups> {0.00867209,
    0.00162426, 0.0102716, 0.0390447, 0.0192576, 0.374888, 0.430599};
  sigma_s["MOX-8.7%%"] = std::array<double, num_groups*num_groups>
      {0.131504, 0.042046, 8.6972E-6, 5.1938E-9, 0.0, 0.0, 0.0,
      0.0, 0.330403, 0.0016463, 2.6006E-9, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.461792, 0.0024749, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.468021, 0.005433, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.8597E-4, 0.285771, 0.0083973, 8.928E-9,
      0.0, 0.0, 0.0, 0.0, 0.0023916, 0.247614, 0.012322,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0089681, 0.256093};
  chi["MOX-8.7%%"] = std::array<double, num_groups> {0.58791, 0.41176,
    3.3906E-4, 1.1761E-7, 0.0, 0.0, 0.0};
  sigma_t["MOX-8.7%%"] = std::array<double, num_groups> {0.183045, 0.336705,
    0.500507, 0.606174, 0.502754, 0.921028, 0.955231};

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

  /* Define control rod cross-sections */
  sigma_a["Control Rod"] = std::array<double, num_groups> {1.70490E-03,
    8.36224E-03,8.37901E-02, 3.97797E-01, 6.98763E-01, 9.29508E-01,
    1.17836E+00};
  nu_sigma_f["Control Rod"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0,
    0};
  sigma_f["Control Rod"] = std::array<double, num_groups>  {0, 0, 0, 0, 0, 0,
    0};
  sigma_s["Control Rod"] = std::array<double, num_groups*num_groups>
      {1.70563E-01, 4.44012E-02, 9.83670E-05, 1.27786E-07, 0., 0., 0.,
        0., 4.71050E-01, 6.85480E-04, 3.91395E-10, 0., 0., 0.,
        0., 0., 8.01859E-01, 7.20132E-04, 0., 0., 0.,
        0., 0., 0., 5.70752E-01, 1.46015E-03, 0., 0.,
        0., 0., 0., 6.55562E-05, 2.07838E-01, 3.81486E-03, 3.69760E-09,
        0., 0., 0., 0., 1.02427E-03, 2.02465E-01, 4.75290E-03,
        0., 0., 0., 0., 0., 3.53043E-03, 6.58597E-01};


  chi["Control Rod"] = std::array<double, num_groups> {0, 0, 0, 0, 0, 0, 0};
  sigma_t["Control Rod"] = std::array<double, num_groups> {2.16768E-01,
    4.80098E-01, 8.86369E-01, 9.70009E-01, 9.10482E-01, 1.13775E+00,
    1.84048E+00};

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
  XPlane xmin(-42.84);
  XPlane xmax( 42.84);
  YPlane ymin(-42.84);
  YPlane ymax( 42.84);
  ZPlane zmin(-42.84);
  ZPlane zmax( 42.84);

  xmin.setBoundaryType(REFLECTIVE);
  xmax.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  ymax.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(REFLECTIVE);
  zmax.setBoundaryType(REFLECTIVE);

  /* Create z-cylinders for the fuel as well as to discretize the moderator
   * into rings */
  ZCylinder fuel_radius(0.0, 0.0, 0.54);

  /* Create cells and universes */
  log_printf(NORMAL, "Creating cells...");

  /* Moderator rings */
  Cell* moderator = new Cell(2, "mod");
  moderator->setFill(materials["Water"]);
  moderator->addSurface(+1, &fuel_radius);
  moderator->setNumSectors(8);

  /* UO2 pin cell */
  Cell* uo2_cell = new Cell(3, "uo2");
  uo2_cell->setNumSectors(4);
  uo2_cell->setFill(materials["UO2"]);
  uo2_cell->addSurface(-1, &fuel_radius);

  Universe* uo2 = new Universe();
  uo2->addCell(uo2_cell);
  uo2->addCell(moderator);

  /* 4.3% MOX pin cell */
  Cell* mox43_cell = new Cell(4, "mox43");
  mox43_cell->setNumSectors(4);
  mox43_cell->setFill(materials["MOX-4.3%%"]);
  mox43_cell->addSurface(-1, &fuel_radius);

  Universe* mox43 = new Universe();
  mox43->addCell(mox43_cell);
  mox43->addCell(moderator);

  /* 7% MOX pin cell */
  Cell* mox7_cell = new Cell(5, "mox7");
  mox7_cell->setNumSectors(4);
  mox7_cell->setFill(materials["MOX-7%%"]);
  mox7_cell->addSurface(-1, &fuel_radius);

  Universe* mox7 = new Universe();
  mox7->addCell(mox7_cell);
  mox7->addCell(moderator);

  /* 8.7% MOX pin cell */
  Cell* mox87_cell = new Cell(6, "mox87");
  mox87_cell->setNumSectors(4);
  mox87_cell->setFill(materials["MOX-8.7%%"]);
  mox87_cell->addSurface(-1, &fuel_radius);

  Universe* mox87 = new Universe();
  mox87->addCell(mox87_cell);
  mox87->addCell(moderator);

  /* Fission chamber pin cell */
  Cell* fission_chamber_cell = new Cell(7, "fc");
  fission_chamber_cell->setNumRings(5);
  fission_chamber_cell->setNumSectors(4);
  fission_chamber_cell->setFill(materials["Fission Chamber"]);
  fission_chamber_cell->addSurface(-1, &fuel_radius);

  Universe* fission_chamber = new Universe();
  fission_chamber->addCell(fission_chamber_cell);
  fission_chamber->addCell(moderator);

  /* Guide tube pin cell */
  Cell* guide_tube_cell = new Cell(8, "gtc");
  guide_tube_cell->setNumRings(5);
  guide_tube_cell->setNumSectors(4);
  guide_tube_cell->setFill(materials["Guide Tube"]);
  guide_tube_cell->addSurface(-1, &fuel_radius);

  Universe* guide_tube = new Universe();
  guide_tube->addCell(guide_tube_cell);
  guide_tube->addCell(moderator);

  /* Control rod pin cell */
  Cell* control_rod_cell = new Cell(9, "crc");
  control_rod_cell->setNumRings(5);
  control_rod_cell->setNumSectors(4);
  control_rod_cell->setFill(materials["Control Rod"]);
  control_rod_cell->addSurface(-1, &fuel_radius);

  Universe* control_rod = new Universe();
  control_rod->addCell(control_rod_cell);
  control_rod->addCell(moderator);

  /* Reflector */
  Cell* reflector_cell = new Cell(110, "rc");
  reflector_cell->setFill(materials["Water"]);

  Universe* reflector = new Universe();
  reflector->addCell(reflector_cell);

  /* Moderator pin */
  Cell* moderator_pin_cell = new Cell(21, "mpc");
  moderator_pin_cell->setNumRings(5);
  moderator_pin_cell->setNumSectors(4);
  moderator_pin_cell->setFill(materials["Water"]);
  moderator_pin_cell->addSurface(-1, &fuel_radius);

  Universe* moderator_pin = new Universe();
  moderator_pin->addCell(moderator_pin_cell);
  moderator_pin->addCell(moderator);

  /* Cells */
  Cell* assembly1_cell = new Cell(10, "ac1");
  Cell* assembly2_cell = new Cell(11, "ac2");
  Cell* refined_reflector_cell = new Cell(12, "rrc");
  Cell* right_reflector_cell = new Cell(13,"rrc2");
  Cell* corner_reflector_cell = new Cell(14, "crc");
  Cell* bottom_reflector_cell = new Cell(15, "brc");
  Cell* reflector_assembly_cell = new Cell(16, "arc");
  Cell* assembly1_cell_rodded = new Cell(17, "acr1");
  Cell* assembly2_cell_rodded = new Cell(18, "acr2");
  Cell* reflector_rodded_cell = new Cell(19, "rrc");

  Universe* assembly1 = new Universe();
  Universe* assembly2 = new Universe();
  Universe* refined_reflector = new Universe();
  Universe* right_reflector = new Universe();
  Universe* corner_reflector = new Universe();
  Universe* bottom_reflector = new Universe();
  Universe* reflector_assembly = new Universe();
  Universe* assembly1_rodded = new Universe();
  Universe* assembly2_rodded = new Universe();
  Universe* reflector_rodded = new Universe();

  assembly1->addCell(assembly1_cell);
  assembly2->addCell(assembly2_cell);
  refined_reflector->addCell(refined_reflector_cell);
  right_reflector->addCell(right_reflector_cell);
  corner_reflector->addCell(corner_reflector_cell);
  bottom_reflector->addCell(bottom_reflector_cell);
  reflector_assembly->addCell(reflector_assembly_cell);
  assembly1_rodded->addCell(assembly1_cell_rodded);
  assembly2_rodded->addCell(assembly2_cell_rodded);
  reflector_rodded->addCell(reflector_rodded_cell);

  /* Root Cell* */
  Cell* root_cell = new Cell(20, "root");
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
  assembly1_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* matrix1[17*17*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::map<int, Universe*> names = {{1, uo2}, {2, guide_tube},
                                      {3, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<17*17; n++)
        matrix1[z*17*17 + n] = names[mold[n]];

    assembly1_lattice->setUniverses(axial_refines, 17, 17, matrix1);
  }
  assembly1_cell->setFill(assembly1_lattice);

  /* Top right, bottom left 17 x 17 assemblies */
  Lattice* assembly2_lattice = new Lattice();
  assembly2_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* matrix2[17*17*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
                        1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1,
                        1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1,
                        1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1,
                        1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1,
                        1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1,
                        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::map<int, Universe*> names = {{1, mox43}, {2, mox7}, {3, mox87},
                                      {4, guide_tube}, {5, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<17*17; n++)
        matrix2[z*17*17 + n] = names[mold[n]];

    assembly2_lattice->setUniverses(axial_refines, 17, 17, matrix2);
  }
  assembly2_cell->setFill(assembly2_lattice);

  /* Top left, bottom right 17 x 17 assemblies */
  Lattice* assembly1_lattice_rodded = new Lattice();
  assembly1_lattice_rodded->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* matrix1_rodded[17*17*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::map<int, Universe*> names = {{1, uo2}, {2, control_rod},
                                      {3, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<17*17; n++)
        matrix1_rodded[z*17*17 + n] = names[mold[n]];

    assembly1_lattice_rodded->setUniverses(axial_refines, 17, 17,
                                           matrix1_rodded);
  }
  assembly1_cell_rodded->setFill(assembly1_lattice_rodded);

  /* Top right, bottom left 17 x 17 assemblies */
  Lattice* assembly2_lattice_rodded = new Lattice();
  assembly2_lattice_rodded->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* matrix2_rodded[17*17*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
                        1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1,
                        1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1,
                        1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1,
                        1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1,
                        1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1,
                        1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1,
                        1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1,
                        1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::map<int, Universe*> names = {{1, mox43}, {2, mox7}, {3, mox87},
                                      {4, control_rod}, {5, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<17*17; n++)
        matrix2_rodded[z*17*17 + n] = names[mold[n]];

    assembly2_lattice_rodded->setUniverses(axial_refines, 17, 17,
                                           matrix2_rodded);
  }
  assembly2_cell_rodded->setFill(assembly2_lattice_rodded);

  /* Top right, bottom left 17 x 17 assemblies */
  Lattice* reflector_rodded_lattice = new Lattice();
  reflector_rodded_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* matrix_ref_rodded[17*17*axial_refines];
  {
    int mold[17*17] =  {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1,
                        1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::map<int, Universe*> names = {{1, moderator_pin}, {2, control_rod},
                                      {3, fission_chamber}};
    for (int z=0; z < axial_refines; z++)
      for (int n=0; n<17*17; n++)
        matrix_ref_rodded[z*17*17 + n] = names[mold[n]];

    reflector_rodded_lattice->setUniverses(axial_refines, 17, 17,
                                           matrix_ref_rodded);
  }
  reflector_rodded_cell->setFill(reflector_rodded_lattice);

  /* Sliced up water cells - semi finely spaced */
  Lattice* refined_ref_lattice = new Lattice();
  int nr = 3;
  refined_ref_lattice->setWidth(1.26/nr, 1.26/nr, 7.14/axial_refines);
  Universe* refined_ref_matrix[nr*nr*axial_refines];
  for (int z=0; z < axial_refines; z++)
    for (int n=0; n<nr*nr; n++)
      refined_ref_matrix[z*nr*nr + n] = reflector;
  refined_ref_lattice->setUniverses(axial_refines, nr, nr, refined_ref_matrix);
  refined_reflector_cell->setFill(refined_ref_lattice);

  /* Sliced up water cells - right side of geometry */
  Lattice* right_ref_lattice = new Lattice();
  right_ref_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* right_ref_matrix[17*17*axial_refines];
  for (int z=0; z < axial_refines; z++) {
    for (int i=0; i<17; i++) {
      for (int j=0; j<17; j++) {
        int index =  17*j + i;
        if (i<11)
          right_ref_matrix[17*17*z + index] = refined_reflector;
        else
          right_ref_matrix[17*17*z + index] = reflector;
      }
    }
  }
  right_ref_lattice->setUniverses(axial_refines, 17, 17, right_ref_matrix);
  right_reflector_cell->setFill(right_ref_lattice);

  /* Sliced up water cells for bottom corner of geometry */
  Lattice* corner_ref_lattice = new Lattice();
  corner_ref_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* corner_ref_matrix[17*17*axial_refines];
  for (int z=0; z < axial_refines; z++) {
    for (int i=0; i<17; i++) {
      for (int j=0; j<17; j++) {
        int index = 17*j + i;
        if (i<11 && j<11)
          corner_ref_matrix[17*17*z + index] = refined_reflector;
        else
          corner_ref_matrix[17*17*z + index] = reflector;
      }
    }
  }
  corner_ref_lattice->setUniverses(axial_refines, 17, 17, corner_ref_matrix);
  corner_reflector_cell->setFill(corner_ref_lattice);

  /* Sliced up water cells for bottom of geometry */
  Lattice* bottom_ref_lattice = new Lattice();
  bottom_ref_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* bottom_ref_matrix[17*17*axial_refines];
  for (int z=0; z < axial_refines; z++) {
    for (int i=0; i<17; i++) {
      for (int j=0; j<17; j++) {
        int index = 17*j + i;
        if (j<11)
          bottom_ref_matrix[17*17*z + index] = refined_reflector;
        else
          bottom_ref_matrix[17*17*z + index] = reflector;
      }
    }
  }
  bottom_ref_lattice->setUniverses(axial_refines, 17, 17, bottom_ref_matrix);
  bottom_reflector_cell->setFill(bottom_ref_lattice);

  /* Reflector assembly (unrodded) */
  Lattice* assembly_ref_lattice = new Lattice();
  assembly_ref_lattice->setWidth(1.26, 1.26, 7.14/axial_refines);
  Universe* assembly_ref_matrix[17*17*axial_refines];
  for (int n=0; n < axial_refines*17*17; n++)
    assembly_ref_matrix[n] = refined_reflector;

  assembly_ref_lattice->setUniverses(axial_refines, 17, 17, assembly_ref_matrix);
  reflector_assembly_cell->setFill(assembly_ref_lattice);

  /* 3 x 3 x 9 core to represent two bundles and water */
  Lattice* full_geometry = new Lattice();
  full_geometry->setWidth(21.42, 21.42, 7.14);
/*
  Universe* universes[] = {
    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,
    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,

    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,
    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,

    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,
    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,

    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,
    assembly1_rodded,   assembly2_rodded, assembly1_rodded,   assembly2_rodded,
    assembly2_rodded,   assembly1,        assembly2_rodded,   assembly1,

    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1_rodded,   assembly2,        assembly1_rodded,   assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,

    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1,
    assembly1,          assembly2,        assembly1,          assembly2,
    assembly2,          assembly1,        assembly2,          assembly1};
*/
  Universe* universes[] = {
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,

    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1,
    assembly1,   assembly2, assembly1,   assembly2,
    assembly2,   assembly1, assembly2,   assembly1};

  full_geometry->setUniverses(12, 4, 4, universes);

  /* Fill root cell with lattice */
  root_cell->setFill(full_geometry);

  /* Create CMFD mesh */
  log_printf(NORMAL, "Creating CMFD mesh...");

  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.0);
  cmfd->setLatticeStructure(68, 68, 48);
  std::vector<std::vector<int> > cmfd_group_structure;
  cmfd_group_structure.resize(2);
  for (int g=0; g<3; g++)
    cmfd_group_structure.at(0).push_back(g+1);
  for (int g=3; g<7; g++)
    cmfd_group_structure.at(1).push_back(g+1);
  cmfd->setGroupStructure(cmfd_group_structure);
  cmfd->setKNearest(3);

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");
  Geometry geometry;
  geometry.setRootUniverse(root_universe);
#ifdef MPIx
  geometry.setDomainDecomposition(nx, ny, nz, MPI_COMM_WORLD);
#else
  geometry.setNumDomainModules(nx, ny, nz);
#endif
  geometry.setCmfd(cmfd);
  geometry.initializeFlatSourceRegions();

  /* Generate tracks */
  log_printf(NORMAL, "Initializing the track generator...");
  Quadrature* quad = new GLPolarQuad();
  quad->setNumAzimAngles(num_azim);
  quad->setNumPolarAngles(num_polar);
  TrackGenerator3D track_generator(&geometry, num_azim, num_polar, azim_spacing,
                                   polar_spacing);
  track_generator.setNumThreads(num_threads);
  track_generator.setQuadrature(quad);
  track_generator.setSegmentFormation(OTF_STACKS);
  std::vector<double> seg_zones {-42.84, -14.28, 14.28, 42.84};
  track_generator.setSegmentationZones(seg_zones);
  track_generator.generateTracks();

  /* Run simulation */
  CPULSSolver solver(&track_generator);
  solver.setNumThreads(num_threads);
  solver.setConvergenceThreshold(tolerance);
  solver.computeEigenvalue(max_iters);
  solver.printTimerReport();

  log_printf(TITLE, "Finished");
#ifdef MPIx
  MPI_Finalize();
#endif
  return 0;
}

#include <vector>
#include <iostream>

std::vector<std::vector<int> > get_group_structure(int num_groups, 
                                                   int num_cmfd_groups);

inline std::vector<std::vector<int> > get_group_structure(int num_groups, 
                                                          int num_cmfd_groups) {

  std::vector<std::vector<int> > cmfd_group_structure;
  cmfd_group_structure.resize(num_cmfd_groups);

  if (num_groups == num_cmfd_groups) {
    for (int g=0; g < num_groups; g++)
      cmfd_group_structure.at(g).push_back(g+1);
  }
  else if (num_cmfd_groups == 1) {
    for (int g=0; g < num_groups; g++)
      cmfd_group_structure.at(0).push_back(g+1);
  }
  else if (num_groups == 70) {
    if (num_cmfd_groups == 2) {
      for (int g=0; g<46; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=46; g<70; g++)
        cmfd_group_structure.at(1).push_back(g+1);
    }
    else if (num_cmfd_groups == 4) {
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<15; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=15; g<46; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=46; g<70; g++)
        cmfd_group_structure.at(3).push_back(g+1);
    }
    else if (num_cmfd_groups == 8) {
      for (int g=0; g<6; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=6; g<12; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=12; g<19; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=19; g<24; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=24; g<33; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=33; g<53; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=53; g<61; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=61; g<70; g++)
        cmfd_group_structure.at(7).push_back(g+1);
    }
    else if (num_cmfd_groups == 14) {
      for (int g=0; g<2; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=2; g<6; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=6; g<9; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=9; g<12; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=12; g<16; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=16; g<19; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=19; g<21; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=21; g<24; g++)
        cmfd_group_structure.at(7).push_back(g+1);
      for (int g=24; g<27; g++)
        cmfd_group_structure.at(8).push_back(g+1);
      for (int g=27; g<33; g++)
        cmfd_group_structure.at(9).push_back(g+1);
      for (int g=33; g<53; g++)
        cmfd_group_structure.at(10).push_back(g+1);
      for (int g=53; g<61; g++)
        cmfd_group_structure.at(11).push_back(g+1);
      for (int g=61; g<68; g++)
        cmfd_group_structure.at(12).push_back(g+1);
      for (int g=68; g<70; g++)
        cmfd_group_structure.at(13).push_back(g+1);
    }

    else {
      std::cout << "ERROR: CMFD group structure not found" << std::endl;
      exit(1);
    }
  }
  else if (num_groups == 40) {
    if (num_cmfd_groups == 2) {
      for (int g=0; g<28; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=28; g<40; g++)
        cmfd_group_structure.at(1).push_back(g+1);
    }
    else if (num_cmfd_groups == 4) {
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<9; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=9; g<28; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=28; g<40; g++)
        cmfd_group_structure.at(3).push_back(g+1);
    }
    else {
      std::cout << "ERROR: CMFD group structure not found" << std::endl;
      exit(1);
    }
  }

  else if (num_groups == 16) {
    if (num_cmfd_groups == 2) {
      for (int g=0; g<10; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=10; g<16; g++)
        cmfd_group_structure.at(1).push_back(g+1);
    }
    else if (num_cmfd_groups == 4) {
      for (int g=0; g<3; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=3; g<4; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=4; g<5; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=12; g<16; g++)
        cmfd_group_structure.at(3).push_back(g+1);
    }
    else {
      std::cout << "ERROR: CMFD group structure not found" << std::endl;
      exit(1);
    }
  }
  else if (num_groups == 8) {
    if (num_cmfd_groups == 2) {
      for (int g=0; g<3; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=3; g<8; g++)
        cmfd_group_structure.at(1).push_back(g+1);
    }
    else if (num_cmfd_groups == 4) {
      for (int g=0; g<2; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=2; g<3; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=3; g<5; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=5; g<8; g++)
        cmfd_group_structure.at(3).push_back(g+1);
    }
    else {
      std::cout << "ERROR: CMFD group structure not found" << std::endl;
      exit(1);
    }
  }
  else {
    std::cout << "ERROR: CMFD group structure not found" << std::endl;
    exit(1);
  }

  return cmfd_group_structure;
}




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
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<15; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=15; g<27; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=27; g<46; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=46; g<52; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=52; g<56; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=56; g<60; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=60; g<70; g++)
        cmfd_group_structure.at(7).push_back(g+1);
    }
    /*
    else if (num_cmfd_groups == 8) {
      for (int g=0; g<7; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=7; g<13; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=13; g<20; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=20; g<25; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=25; g<34; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=34; g<54; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=54; g<62; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=62; g<70; g++)
        cmfd_group_structure.at(7).push_back(g+1);
    }
    */

    else if (num_cmfd_groups == 11) {
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<15; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=15; g<24; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=24; g<25; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=25; g<26; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=26; g<27; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=27; g<46; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=46; g<52; g++)
        cmfd_group_structure.at(7).push_back(g+1);
      for (int g=52; g<56; g++)
        cmfd_group_structure.at(8).push_back(g+1);
      for (int g=56; g<60; g++)
        cmfd_group_structure.at(9).push_back(g+1);
      for (int g=60; g<70; g++)
        cmfd_group_structure.at(10).push_back(g+1);
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
    else if (num_cmfd_groups == 12) {
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<15; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=15; g<27; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=27; g<46; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=46; g<49; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=49; g<52; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=52; g<54; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=54; g<56; g++)
        cmfd_group_structure.at(7).push_back(g+1);
      for (int g=56; g<58; g++)
        cmfd_group_structure.at(8).push_back(g+1);
      for (int g=58; g<60; g++)
        cmfd_group_structure.at(9).push_back(g+1);
      for (int g=60; g<65; g++)
        cmfd_group_structure.at(10).push_back(g+1);
      for (int g=65; g<70; g++)
        cmfd_group_structure.at(11).push_back(g+1);
    }
    else if (num_cmfd_groups == 16) {
      for (int g=0; g<5; g++)
        cmfd_group_structure.at(0).push_back(g+1);
      for (int g=5; g<15; g++)
        cmfd_group_structure.at(1).push_back(g+1);
      for (int g=15; g<27; g++)
        cmfd_group_structure.at(2).push_back(g+1);
      for (int g=27; g<34; g++)
        cmfd_group_structure.at(3).push_back(g+1);
      for (int g=34; g<35; g++)
        cmfd_group_structure.at(4).push_back(g+1);
      for (int g=35; g<37; g++)
        cmfd_group_structure.at(5).push_back(g+1);
      for (int g=37; g<40; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=40; g<42; g++)
        cmfd_group_structure.at(7).push_back(g+1);
      for (int g=42; g<45; g++)
        cmfd_group_structure.at(8).push_back(g+1);
      for (int g=45; g<47; g++)
        cmfd_group_structure.at(9).push_back(g+1);
      for (int g=47; g<50; g++)
        cmfd_group_structure.at(10).push_back(g+1);
      for (int g=50; g<53; g++)
        cmfd_group_structure.at(11).push_back(g+1);
      for (int g=53; g<57; g++)
        cmfd_group_structure.at(12).push_back(g+1);
      for (int g=57; g<61; g++)
        cmfd_group_structure.at(13).push_back(g+1);
      for (int g=61; g<65; g++)
        cmfd_group_structure.at(14).push_back(g+1);
      for (int g=65; g<70; g++)
        cmfd_group_structure.at(15).push_back(g+1);
    }
    else if (num_cmfd_groups == 25) {
      for (int g=0; g<6; g++)
        cmfd_group_structure.at(g).push_back(g+1);
      for (int g=6; g<9; g++)
        cmfd_group_structure.at(6).push_back(g+1);
      for (int g=9; g<14; g++)
        cmfd_group_structure.at(7).push_back(g+1);
      for (int g=14; g<15; g++)
        cmfd_group_structure.at(8).push_back(g+1);
      for (int g=15; g<21; g++)
        cmfd_group_structure.at(9).push_back(g+1);
      for (int g=21; g<25; g++)
        cmfd_group_structure.at(10).push_back(g+1);
      for (int g=25; g<26; g++)
        cmfd_group_structure.at(11).push_back(g+1);
      for (int g=26; g<27; g++)
        cmfd_group_structure.at(12).push_back(g+1);
      for (int g=27; g<31; g++)
        cmfd_group_structure.at(13).push_back(g+1);
      for (int g=31; g<34; g++)
        cmfd_group_structure.at(14).push_back(g+1);
      for (int g=34; g<36; g++)
        cmfd_group_structure.at(15).push_back(g+1);
      for (int g=36; g<39; g++)
        cmfd_group_structure.at(16).push_back(g+1);
      for (int g=39; g<41; g++)
        cmfd_group_structure.at(17).push_back(g+1);
      for (int g=41; g<46; g++)
        cmfd_group_structure.at(18).push_back(g+1);
      for (int g=46; g<49; g++)
        cmfd_group_structure.at(19).push_back(g+1);
      for (int g=49; g<52; g++)
        cmfd_group_structure.at(20).push_back(g+1);
      for (int g=52; g<56; g++)
        cmfd_group_structure.at(21).push_back(g+1);
      for (int g=56; g<60; g++)
        cmfd_group_structure.at(22).push_back(g+1);
      for (int g=60; g<64; g++)
        cmfd_group_structure.at(23).push_back(g+1);
      for (int g=64; g<70; g++)
        cmfd_group_structure.at(24).push_back(g+1);
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




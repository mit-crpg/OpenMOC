/**
 * @file clone.h
 * @brief Routines to copy Material and Track objects to the GPU from CPU.
 * @author May 30, 2013
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */


#include "../DeviceMaterial.h"
#include "../DeviceTrack.h"
#include "GPUQuery.h"
#include <map>

void clone_material(Material* material_h, dev_material* material_d);
void clone_track(Track* track_h, dev_track* track_d,
                        std::map<int, int> &material_IDs_to_indices);

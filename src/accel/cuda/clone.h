#include "../DeviceMaterial.h"
#include "../DeviceTrack.h"

void cloneMaterialOnGPU(Material* material_h, dev_material* material_d);
void cloneTrackOnGPU(Track* track_h, dev_track* track_d);

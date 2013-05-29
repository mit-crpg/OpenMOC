#include "DeviceFlatSourceRegion.h"


/**
 * @brief Given a pointer to a flatsourceregion on the host and a 
 *        flatsourceregion on the device, copy all of the properties from the 
 *        flatsourceregion on the host to the device.
 * @param fsr_h pointer to a flatsourceregion on the host
 * @param fsr_d pointer to a flatsourceregion on the device
 */
void cloneOnDevice(FlatSourceRegion* fsr_h, dev_flatsourceregion* fsr_d) {

   /* Create a temporary flatsourceregion struct on the host and populate it
    * with the data from the FlatSourceRegion object on the host */
    dev_flatsourceregion new_fsr;

    new_fsr._uid = fsr_h->getUid();
    new_fsr._material_uid = fsr_h->getMaterial()->getUid();
    new_fsr._volume = fsr_h->getVolume();

    cudaMemcpy((void*)fsr_d, (void*)&new_fsr, sizeof(dev_flatsourceregion),
	       cudaMemcpyHostToDevice);

    return;
}

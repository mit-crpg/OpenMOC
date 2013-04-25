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
    dev_flatsourceregion* temp =
                   (dev_flatsourceregion*)malloc(sizeof(dev_flatsourceregion));

    int uid = fsr_h->getUid();
    int material_uid = fsr_h->getMaterial()->getUid();
    FP_PRECISION volume = fsr_h->getVolume();

    cudaMemcpy((void*)&fsr_d->_uid, (void*)&uid, sizeof(int), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&fsr_d->_material_uid, (void*)&material_uid, sizeof(int), 
	       cudaMemcpyHostToDevice);
    cudaMemcpy((void*)&fsr_d->_volume, (void*)&volume, sizeof(FP_PRECISION), 
	       cudaMemcpyHostToDevice);

    return;
}

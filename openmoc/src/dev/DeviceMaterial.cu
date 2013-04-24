#include "DeviceMaterial.h"


/**
 * @brief Given a pointer to a material on the host and a material on the 
 *        device, copy all of the properties from the material on the host 
 *        to the device.
 * @param hmaterial pointer to a material on the host
 * @param dmaterial pointer to a material on the device
 */
void cloneOnDevice(Material* hmaterial, dev_material* dmaterial) {

    /* Create a temporary material struct on the host and populate it with
     * the data from the Material objecton the host */
    dev_material temp;
    temp._uid = hmaterial->getUid();

    for (int i=0; i < 7; i++) {
        temp._sigma_a[i] = hmaterial->getSigmaA()[i];
        temp._sigma_f[i] = hmaterial->getSigmaF()[i];
        temp._nu_sigma_f[i] = hmaterial->getNuSigmaF()[i];
        temp._sigma_t[i] = hmaterial->getSigmaT()[i];
        temp._chi[i] = hmaterial->getChi()[i];

        for (int j=0; j < 7; j++)
            temp._sigma_s[i][j] = hmaterial->getSigmaS()[i*7+j];
    }

    cudaMemcpy((void*)dmaterial, (void*)&temp,
                              sizeof(dev_material), cudaMemcpyHostToDevice);

    return;
}

#include "MICQuery.h"


/**
 * @brief Queries a node to determine whether it contains one or more MICs.
 * @return True if the node contains a MIC, false otherwise.
 */
bool machineContainsMIC() {

    int num_devices = machineContainsHowManyMICs();

    if (num_devices > 0)
        return true;
    else
        return false;
}


int machineContainsHowManyMICs() {

    int num_devices = 0;

    #ifdef __INTEL_OFFLOAD
    num_devices = _Offload_number_of_devices();
    #endif

    return num_devices;
}

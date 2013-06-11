#include "DeviceQuery.h"


/**
 * @brief Queries a node to determine whether it contains one or more MICs.
 * @return True if the node contains a MIC, false otherwise.
 */
bool machineContainsMIC() {

    int count = 0;

#ifdef __INTEL_OFFLOAD
    count = _Offload_number_of_devices();
#endif

    log_printf(INFO, "%d Intel MIC devices are present", count);

    if (count > 0)
        return true;
    else
        return false;
}


/**
 * @brief Queries a node to determine whether it contains one or more MICs.
 * @return the number of MIC devices on a node
 */
int getNumMICDevices() {

    int count = 0;

#ifdef __INTEL_OFFLOAD
    count = _Offload_number_of_devices();
#endif

    return count;
}

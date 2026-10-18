#ifndef SENSORS_H
#define SENSORS_H

/* Sensor acquisition interface; hardware details belong in sensors.c. */
#include "ecu_types.h"

/* Prepare sensor peripherals. */
void Sensors_Initialize(void);

/* Start or configure the future DMA acquisition path. */
void Sensors_StartDma(void);

/* Read the newest completed sample. */
bool Sensors_Read(sensor_data_t *data);

#endif
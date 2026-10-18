#ifndef SENSORS_H
#define SENSORS_H

/* Sensor acquisition interface; hardware details belong in sensors.c. */
#include "ecu_types.h"

typedef enum
{
  /* PLACEHOLDER types. */
  SENSOR_TYPE_PRESSURE = 0,
  SENSOR_TYPE_TEMPERATURE
} sensor_type_t;

/* One sensor: its properties plus its latest reading. */
typedef struct
{
  sensor_type_t type;
  uint8_t channel;   /* ADC/SPI channel, placeholder.      */
  float scale;       /* value = raw * scale + offset.      */
  float offset;
  float min_valid;   /* Plausibility range.                */
  float max_valid;
  uint32_t raw;      /* Last raw reading.                  */
  float value;       /* Last converted value.              */
  bool valid;
} sensor_t;

/* Prepare sensor peripherals. */
void Sensors_Initialize(void);

/* Start or configure the future DMA acquisition path. */
void Sensors_StartDma(void);

/* Copy every sensor's newest value into the snapshot. */
bool Sensors_Read(sensor_data_t *data);

/* Read-only access to one sensor's properties and reading. */
const sensor_t *Sensors_Get(sensor_id_t id);

#endif

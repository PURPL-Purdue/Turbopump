#ifndef SENSORS_H
#define SENSORS_H

/* Sensor acquisition interface; hardware details belong in sensors.c. */
#include "ecu_types.h"
#include "stm32h7xx_hal.h"

typedef enum
{
  /* PLACEHOLDER types. */
  SENSOR_TYPE_PRESSURE = 0,
  SENSOR_TYPE_TEMPERATURE
} sensor_type_t;

/* One sensor: its properties, hardware binding and latest reading. */
typedef struct
{
  sensor_type_t type;

  /* Hardware binding, unassigned (NULL/0) until the pinout is final.  Add the
   * ADC/SPI handle pointer once that peripheral is enabled in CubeMX. */
  GPIO_TypeDef *cs_port; /* External ADC chip select.          */
  uint16_t cs_pin;
  uint8_t channel;   /* ADC/SPI channel, placeholder.      */
  float scale;       /* value = raw * scale + offset.      */
  float offset;
  float min_valid;   /* Plausibility range.                */
  float max_valid;
  float redline;     /* Fault limit, TBD; unset (0) for now. */
  uint32_t raw;      /* Last raw reading.                  */
  float value;       /* Last converted value.              */
  bool valid;
} sensor_t;

/* Configure each sensor's hardware from its table entry.  The sensor objects
 * themselves are static (sensors.c). */
void Sensors_Initialize(void);

/* Start or configure the future DMA acquisition path. */
void Sensors_StartDma(void);

/* Copy every sensor's newest value into the snapshot. */
bool Sensors_Read(sensor_data_t *data);

/* Read-only access to one sensor's properties and reading. */
const sensor_t *Sensors_Get(sensor_id_t id);

#endif

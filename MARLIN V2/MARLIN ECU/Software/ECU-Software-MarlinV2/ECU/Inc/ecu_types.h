#ifndef ECU_TYPES_H
#define ECU_TYPES_H

#include <stdbool.h>
#include <stdint.h>

/* Shared data contracts used by the ECU modules. */

/* Period of the control task loop. */
#define CONTROL_PERIOD_MS 1U

/* Minimum time between frames sent to telemetry. */
#define TELEMETRY_PERIOD_MS 5U

typedef enum
{
  ECU_STATE_STANDBY = 0,
  ECU_STATE_STARTUP,
  ECU_STATE_RUNNING,
  ECU_STATE_SHUTDOWN,
  ECU_STATE_FAULT
} ecu_state_t;

/*
 * Every output the ECU drives.  These names are PLACEHOLDERS until the test
 * stand's valve/actuator list is final.
 */
typedef enum
{
  //all placeholder, will put full amout of valves and solenoids later
  ACT_VALVE_0 = 0,
  ACT_VALVE_1,
  ACT_VALVE_2,
  ACT_VALVE_3,
  ACT_VALVE_4,
  ACT_VALVE_5,
  ACT_VALVE_6,
  ACT_VALVE_7,
  ACT_COUNT
} actuator_id_t;

/*
 * Every sensor the ECU reads.  These names are PLACEHOLDERS until the test
 * stand's instrumentation list is final.
 */
typedef enum
{
  SNS_PT_0 = 0,
  SNS_PT_1,
  SNS_PT_2,
  SNS_PT_3,
  SNS_TC_0,
  SNS_TC_1,
  SNS_COUNT
} sensor_id_t;

/* One value per output (0 = off/closed). */
typedef struct
{
  uint16_t out[ACT_COUNT];
} ecu_command_t;

/* Snapshot of every sensor's value for one control cycle. */
typedef struct
{
  float value[SNS_COUNT];  /* Engineering units. */
  bool valid[SNS_COUNT];
  uint32_t timestamp_ms;   /* Synced clock. */
} sensor_data_t;

/*
 * Frame received from Mission Control, delivered to the control task through
 * the inbound queue.  Contents not defined yet.
 */
typedef struct
{
  uint8_t placeholder; /* TODO: define inbound frame fields. */
} ecu_inbound_frame_t;

/*
 * Frame sent to Mission Control, delivered to the telemetry task through the
 * outbound queue.  Values only; device properties stay on the ECU.
 */
typedef struct
{
  uint32_t timestamp_ms;              /* Synced clock. */
  ecu_state_t state;
  float sensor_value[SNS_COUNT];      /* Engineering units. */
  uint16_t actuator_value[ACT_COUNT]; /* Commanded value. */
  /* TODO: further fields TBD. */
} ecu_outbound_frame_t;

#endif

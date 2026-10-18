#ifndef ECU_TYPES_H
#define ECU_TYPES_H

#include <stdbool.h>
#include <stdint.h>

/* Shared data contracts used by the ECU modules. */

/* Length of one control cycle, and so the time between two autosequencer rows. */
#define CONTROL_PERIOD_MS 1U

typedef enum
{
  ECU_STATE_STANDBY = 0, /* Outputs off, manual control from Mission Control. */
  ECU_STATE_STARTUP,     /* Main sequence running, manual control locked out. */
  ECU_STATE_RUNNING,     /* Main sequence running, entered automatically.     */
  ECU_STATE_SHUTDOWN,    /* Shutdown rows at the end of the main sequence.    */
  ECU_STATE_FAULT        /* Separate fault sequence, then back to STANDBY.    */
} ecu_state_t;

/*
 * Every output the ECU drives.  These names are PLACEHOLDERS until the test
 * stand's valve/actuator list is final.
 */
typedef enum
{
  ACT_VALVE_0 = 0,
  ACT_VALVE_1,
  ACT_VALVE_2,
  ACT_VALVE_3,
  ACT_VALVE_4,
  ACT_VALVE_5,
  ACT_IGNITER_0,
  ACT_IGNITER_1,
  ACT_COUNT
} actuator_id_t;

/* One value per output (0 = off/closed): one autosequencer row, or the manual
 * outputs in STANDBY. */
typedef struct
{
  uint16_t out[ACT_COUNT];
} ecu_command_t;

/* Messages from Mission Control, delivered to the control task by a queue. */
typedef enum
{
  MC_CMD_HEARTBEAT = 0,  /* Proves the link is alive; carries no data.       */
  MC_CMD_SET_ACTUATOR,   /* STANDBY only: out[actuator] = value.             */
  MC_CMD_START,          /* STANDBY only: begin the main sequence.           */
  MC_CMD_ABORT           /* Stop the sequence (goes to FAULT).               */
} mc_command_type_t;

typedef struct
{
  mc_command_type_t type;
  uint16_t actuator;
  uint16_t value;
} mc_command_t;

/* Each bit is one fault; zero means no fault. */
typedef uint32_t fault_flags_t;

#define FAULT_NONE         0U
#define FAULT_MC_LINK_LOST (1UL << 0) /* No MC message within the timeout. */
#define FAULT_MC_ABORT     (1UL << 1) /* Operator pressed abort.           */

typedef struct
{
  /* Generic fields only; real sensor channels can be added later. */
  uint32_t timestamp_ms;
  bool valid;
} sensor_data_t;

typedef struct
{
  /* One snapshot passed from the control task to the telemetry task. */
  uint32_t timestamp_ms;
  ecu_state_t state;
  uint32_t sequence_tick;    /* Row being executed in the active table. */
  sensor_data_t sensors;
  ecu_command_t command;     /* What was sent to the outputs.           */
  fault_flags_t faults;
} telemetry_frame_t;

#endif

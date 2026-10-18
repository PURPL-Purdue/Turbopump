#ifndef CONTROL_STATE_MACHINE_H
#define CONTROL_STATE_MACHINE_H

/*
 * ECU state machine.  Each control cycle the control task:
 *   1. passes every waiting inbound frame to HandleCommand(),
 *   2. calls Update() once with this cycle's sensors and faults,
 *   3. sends machine.output_command to the actuators.
 *
 * Control logic not defined yet.
 */

#include <stdbool.h>

#include "autosequencer.h"
#include "ecu_types.h"

typedef struct
{
  ecu_state_t state;
  uint32_t step_index;              /* Active step of the running sequence.   */

  const autosequencer_table_t *main_table;
  const autosequencer_table_t *fault_table;

  ecu_command_t output_command;     /* What goes to the actuators this cycle. */
  fault_flags_t faults;             /* Faults seen this cycle.                */
} control_state_machine_t;

/* Enter STANDBY with every output off. */
void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *main_table,
                                    const autosequencer_table_t *fault_table);

/* Take in one inbound frame from Mission Control. */
void ControlStateMachine_HandleCommand(control_state_machine_t *machine,
                                       const ecu_inbound_frame_t *frame);

/* Run one control cycle. */
void ControlStateMachine_Update(control_state_machine_t *machine,
                                const sensor_data_t *sensors,
                                fault_flags_t faults);

#endif

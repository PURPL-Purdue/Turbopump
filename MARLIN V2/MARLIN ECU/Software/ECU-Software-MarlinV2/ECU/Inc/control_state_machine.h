#ifndef CONTROL_STATE_MACHINE_H
#define CONTROL_STATE_MACHINE_H

/*
 * ECU state machine.  Each control cycle the control task:
 *   1. passes every waiting Mission Control message to HandleCommand(),
 *   2. calls Update() once with this cycle's faults,
 *   3. sends machine.output_command to the actuators.
 *
 * The state decides where the output comes from: the manual outputs in
 * STANDBY, the main autosequencer table in STARTUP/RUNNING/SHUTDOWN, and the
 * fault table in FAULT.
 */

#include <stdbool.h>

#include "autosequencer.h"
#include "ecu_types.h"

/* No MC message for this long means the link is lost (placeholder value). */
#define MC_LINK_TIMEOUT_MS 500U

typedef struct
{
  ecu_state_t state;
  uint32_t sequence_tick;           /* Next row of the active table.          */

  const autosequencer_table_t *main_table;
  const autosequencer_table_t *fault_table;

  ecu_command_t manual_command;     /* Outputs commanded by MC in STANDBY.    */
  ecu_command_t output_command;     /* What goes to the actuators this cycle. */

  bool mc_link_alive;
  uint32_t ticks_since_mc_message;
  bool start_requested;             /* Set by HandleCommand, used by Update.  */
  bool abort_requested;

  fault_flags_t faults;             /* Faults seen this cycle.                */
} control_state_machine_t;

/* Enter STANDBY with every output off. */
void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *main_table,
                                    const autosequencer_table_t *fault_table);

/* Take in one Mission Control message. */
void ControlStateMachine_HandleCommand(control_state_machine_t *machine,
                                       const mc_command_t *command);

/* Run one control cycle: state changes, then this cycle's output_command. */
void ControlStateMachine_Update(control_state_machine_t *machine,
                                fault_flags_t faults);

#endif

#ifndef CONTROL_STATE_MACHINE_H
#define CONTROL_STATE_MACHINE_H

/*
 * ECU state machine.  Called by the control task each cycle (app_tasks.c):
 *   - StartSequence() when Mission Control requests a start,
 *   - ClearFault() when Mission Control requests FAULT -> STANDBY,
 *   - Fault() when the fault logic trips (callable at any point; the FAULT
 *     state persists, so a fault raised late in one cycle holds into the next),
 *   - then Update() runs the active state and sets output_command,
 * then output_command goes to the actuators.
 *
 * Control logic not defined yet.  Only STANDBY has been designed.
 */

#include <stdbool.h>

#include "autosequencer.h"
#include "ecu_types.h"

typedef struct
{
  ecu_state_t state;
  ecu_state_t previous_state;       /* State on the last Update(); a change
                                     * means this cycle is the state's entry. */
  uint32_t step_index;              /* Active step of the sequence.           */

  const autosequencer_table_t *sequence;

  ecu_command_t output_command;     /* What goes to the actuators this cycle. */
} control_state_machine_t;

/* Start in STANDBY. */
void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *sequence);

/*
 * Requests.  Each only sets the next state; Update() runs that state's entry
 * actions on its first cycle there.
 */

/* Move to STARTUP.  Ignored outside STANDBY. */
void ControlStateMachine_StartSequence(control_state_machine_t *machine);

/* Move to FAULT.  Safe to call every cycle a fault is present: already being
 * in FAULT is not a new entry, so the fault section is not restarted. */
void ControlStateMachine_Fault(control_state_machine_t *machine);

/* Return from FAULT to STANDBY on Mission Control's command.  Ignored
 * outside FAULT. */
void ControlStateMachine_ClearFault(control_state_machine_t *machine);

/* Run the active state: its entry actions on the first cycle there, then
 * its per-cycle work (e.g. advance the autosequence) and set output_command. */
void ControlStateMachine_Update(control_state_machine_t *machine,
                                const sensor_data_t *sensors);

#endif

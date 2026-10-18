#include "control_state_machine.h"

#include <stddef.h>

/*
 * Run the main sequence: hold the active step's outputs until its exit
 * conditions are met, then advance.  Not called yet.
 */
__attribute__((unused))
static void run_main_sequence(control_state_machine_t *machine,
                              const sensor_data_t *sensors)
{
  const autosequencer_step_t *step = Autosequencer_GetStep(machine->main_table,
                                                           machine->step_index);

  if (step == NULL)
  {
    /* TODO: end of the main sequence. */
    return;
  }

  machine->output_command = step->command;

  if (Autosequencer_StepComplete(step, sensors))
  {
    machine->step_index++;
  }

  /* TODO: STARTUP/RUNNING/SHUTDOWN transitions, driven by conditions. */
}

void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *main_table,
                                    const autosequencer_table_t *fault_table)
{
  *machine = (control_state_machine_t){0};
  machine->main_table = main_table;
  machine->fault_table = fault_table;
  machine->state = ECU_STATE_STANDBY;
}

void ControlStateMachine_HandleCommand(control_state_machine_t *machine,
                                       const ecu_inbound_frame_t *frame)
{
  (void)machine;
  (void)frame;
  /* TODO: command handling. */
}

void ControlStateMachine_Update(control_state_machine_t *machine,
                                const sensor_data_t *sensors,
                                fault_flags_t faults)
{
  (void)sensors;
  machine->faults = faults;
  /* TODO: state machine logic. */
}

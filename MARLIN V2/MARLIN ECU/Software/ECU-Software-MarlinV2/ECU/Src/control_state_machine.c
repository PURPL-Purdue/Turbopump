#include "control_state_machine.h"

#include <stddef.h>

/*
 * Outputs driven on entering STANDBY (and at boot): each entry is that
 * actuator's safe value.  PLACEHOLDER (all de-energized) until the valve list
 * is final; e.g. vents may need to be commanded open here.
 */
static const ecu_command_t standby_safe_command = {0};

/*
 * Run the sequence: hold the active step's outputs until its exit conditions
 * are met, then advance.
 */
static void run_sequence(control_state_machine_t *machine,
                         const sensor_data_t *sensors)
{
  const autosequencer_step_t *step = Autosequencer_GetStep(machine->sequence,
                                                           machine->step_index);

  if (step == NULL)
  {
    /* TODO: end of the sequence. */
    return;
  }

  machine->output_command = step->command;

  if (Autosequencer_StepComplete(step, sensors))
  {
    machine->step_index++;
  }

  /* TODO: state transitions and jumps between sections, driven by conditions
   * (e.g. a completed step sets machine->state to its next state; the new
   * state's entry runs on the next Update()). */
}

void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *sequence)
{
  *machine = (control_state_machine_t){0};
  machine->sequence = sequence;
  machine->state = ECU_STATE_STANDBY; //On boot start in STANDBY
  machine->previous_state = ECU_STATE_STANDBY;
  /* Booting into STANDBY is not an entry (state == previous_state), so set
   * its safe outputs here. */
  machine->output_command = standby_safe_command;
}

void ControlStateMachine_ClockSynced(control_state_machine_t *machine)
{
  machine->clock_synced = true;
}

void ControlStateMachine_StartSequence(control_state_machine_t *machine)
{
  /* Needs: only start from STANDBY, once the clock is synced and a
   * sequence with at least one step is loaded. */
  if ((machine->state != ECU_STATE_STANDBY) || !machine->clock_synced ||
      (machine->sequence == NULL) || (machine->sequence->length == 0U))
  {
    return;
  }
  machine->state = ECU_STATE_STARTUP;
}

void ControlStateMachine_ManualCommand(control_state_machine_t *machine,
                                       const ecu_command_t *command)
{
  if (machine->state != ECU_STATE_STANDBY)
  {
    return;
  }
  machine->manual_command = *command;
  machine->manual_pending = true;
}

void ControlStateMachine_Fault(control_state_machine_t *machine)
{
  machine->state = ECU_STATE_FAULT;
}

void ControlStateMachine_ClearFault(control_state_machine_t *machine)
{
  if (machine->state != ECU_STATE_FAULT)
  {
    return;
  }
  /* TODO: only clear once the fault section has finished (and any other
   * clear conditions are met). */
  machine->state = ECU_STATE_STANDBY;
}

void ControlStateMachine_Update(control_state_machine_t *machine,
                                const sensor_data_t *sensors)
{
  /* First cycle in a state = its entry; entry actions run once.  Staying in
   * the same state (e.g. Fault() called every cycle while a fault is present)
   * is not a new entry, so the fault section is not restarted. */
  const bool entering = (machine->state != machine->previous_state);
  machine->previous_state = machine->state;

  switch (machine->state)
  {
    case ECU_STATE_STANDBY:
      if (entering)
      {
        machine->step_index = 0U;
        machine->output_command = standby_safe_command;
      }
      /* After the entry, so a manual command in the cycle STANDBY is entered
       * still applies.  Latched: output_command holds it until the next. */
      if (machine->manual_pending)
      {
        machine->output_command = machine->manual_command;
      }
      break;

    case ECU_STATE_STARTUP:
      if (entering)
      {
        machine->step_index = 0U; /* TODO: main section start index. */
      }
      run_sequence(machine, sensors);
      break;

    case ECU_STATE_RUNNING:
      /* No entry action: continues the sequence from where STARTUP ended. */
      run_sequence(machine, sensors);
      break;

    case ECU_STATE_SHUTDOWN:
      if (entering)
      {
        /* TODO: jump to the shutdown section of the sequence. */
      }
      run_sequence(machine, sensors);
      break;

    case ECU_STATE_FAULT:
      if (entering)
      {
        /* TODO: jump to the fault section of the sequence; outputs on
         * entering FAULT. */
      }
      /* TODO: run the fault section. */
      break;

    default:
      /* Corrupted state: go to FAULT; its entry runs on the next Update(). */
      machine->state = ECU_STATE_FAULT;
      break;
  }

  /* A manual command is for this cycle only: one accepted in a cycle that
   * then left STANDBY must not be applied on a later STANDBY entry. */
  machine->manual_pending = false;
}

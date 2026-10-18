#include "control_state_machine.h"

/*
 * Transitions:
 *   STANDBY  --START (no faults, MC link alive)----> STARTUP
 *   STARTUP  --reaches running_start_row-----------> RUNNING
 *   RUNNING  --reaches shutdown_start_row----------> SHUTDOWN
 *   SHUTDOWN --end of main table-------------------> STANDBY
 *   STARTUP/RUNNING --fault, abort, MC link lost---> FAULT
 *   SHUTDOWN        --fault, abort-----------------> FAULT
 *   FAULT    --end of fault table------------------> STANDBY
 */

static void enter_state(control_state_machine_t *machine, ecu_state_t state)
{
  machine->state = state;
  machine->sequence_tick = 0U;

  if (state == ECU_STATE_STANDBY)
  {
    /* Every return to STANDBY starts with all outputs off. */
    machine->manual_command = (ecu_command_t){0};
  }
}

static void run_fault(control_state_machine_t *machine)
{
  if (!Autosequencer_GetRow(machine->fault_table, machine->sequence_tick,
                            &machine->output_command))
  {
    enter_state(machine, ECU_STATE_STANDBY);
    machine->output_command = machine->manual_command;
    return;
  }
  machine->sequence_tick++;
}

static void run_main_sequence(control_state_machine_t *machine)
{
  const autosequencer_table_t *table = machine->main_table;

  if (machine->faults != FAULT_NONE)
  {
    enter_state(machine, ECU_STATE_FAULT);
    run_fault(machine);
    return;
  }

  /* Phase changes come from the row markers in the table. */
  if ((machine->state == ECU_STATE_STARTUP) &&
      (machine->sequence_tick >= table->running_start_row))
  {
    machine->state = ECU_STATE_RUNNING;
  }
  if ((machine->state == ECU_STATE_RUNNING) &&
      (machine->sequence_tick >= table->shutdown_start_row))
  {
    machine->state = ECU_STATE_SHUTDOWN;
  }

  if (!Autosequencer_GetRow(table, machine->sequence_tick, &machine->output_command))
  {
    /* End of the main sequence. */
    enter_state(machine, ECU_STATE_STANDBY);
    machine->output_command = machine->manual_command;
    return;
  }
  machine->sequence_tick++;
}

static void run_standby(control_state_machine_t *machine)
{
  /* Abort, or losing Mission Control, turns every manual output off. */
  if (machine->abort_requested || !machine->mc_link_alive)
  {
    machine->manual_command = (ecu_command_t){0};
  }

  if (machine->start_requested && !machine->abort_requested &&
      machine->mc_link_alive && (machine->faults == FAULT_NONE))
  {
    enter_state(machine, ECU_STATE_STARTUP);
    run_main_sequence(machine);
    return;
  }

  machine->output_command = machine->manual_command;
}

void ControlStateMachine_Initialize(control_state_machine_t *machine,
                                    const autosequencer_table_t *main_table,
                                    const autosequencer_table_t *fault_table)
{
  *machine = (control_state_machine_t){0};
  machine->main_table = main_table;
  machine->fault_table = fault_table;
  enter_state(machine, ECU_STATE_STANDBY);
  machine->output_command = machine->manual_command;
}

void ControlStateMachine_HandleCommand(control_state_machine_t *machine,
                                       const mc_command_t *command)
{
  /* Any message from MC proves the link is alive. */
  machine->mc_link_alive = true;
  machine->ticks_since_mc_message = 0U;

  switch (command->type)
  {
    case MC_CMD_SET_ACTUATOR:
      if ((machine->state == ECU_STATE_STANDBY) && (command->actuator < ACT_COUNT))
      {
        machine->manual_command.out[command->actuator] = command->value;
      }
      break;

    case MC_CMD_START:
      /* Checked in Update(), once this cycle's faults are known. */
      machine->start_requested = true;
      break;

    case MC_CMD_ABORT:
      machine->abort_requested = true;
      break;

    default:
      /* MC_CMD_HEARTBEAT: nothing to do beyond the link check above. */
      break;
  }
}

void ControlStateMachine_Update(control_state_machine_t *machine,
                                fault_flags_t faults)
{
  if (machine->mc_link_alive)
  {
    machine->ticks_since_mc_message++;
    if (machine->ticks_since_mc_message >= (MC_LINK_TIMEOUT_MS / CONTROL_PERIOD_MS))
    {
      machine->mc_link_alive = false;
    }
  }

  if ((machine->state == ECU_STATE_STARTUP) || (machine->state == ECU_STATE_RUNNING) ||
      (machine->state == ECU_STATE_SHUTDOWN))
  {
    if (machine->abort_requested)
    {
      faults |= FAULT_MC_ABORT;
    }
    if ((machine->state != ECU_STATE_SHUTDOWN) && !machine->mc_link_alive)
    {
      faults |= FAULT_MC_LINK_LOST;
    }
  }
  machine->faults = faults;

  switch (machine->state)
  {
    case ECU_STATE_STANDBY:
      run_standby(machine);
      break;

    case ECU_STATE_STARTUP:
    case ECU_STATE_RUNNING:
    case ECU_STATE_SHUTDOWN:
      run_main_sequence(machine);
      break;

    default: /* ECU_STATE_FAULT: already shutting down. */
      run_fault(machine);
      break;
  }

  machine->start_requested = false;
  machine->abort_requested = false;
}

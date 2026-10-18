#include "external_input.h"

#include "cmsis_os2.h"

/*
 * This module will receive non-sequence commands, such as manual valve
 * commands or an abort request, from the Mission Control interface.  The
 * transport (Ethernet) decodes packets and submits them here.
 */
static osMessageQueueId_t mc_command_queue;

void ExternalInput_Initialize(void)
{
  mc_command_queue = osMessageQueueNew(MC_COMMAND_QUEUE_LENGTH, sizeof(mc_command_t), NULL);
}

bool ExternalInput_Submit(const mc_command_t *command)
{
  return osMessageQueuePut(mc_command_queue, command, 0U, 0U) == osOK;
}

bool ExternalInput_Receive(mc_command_t *command)
{
  return osMessageQueueGet(mc_command_queue, command, NULL, 0U) == osOK;
}

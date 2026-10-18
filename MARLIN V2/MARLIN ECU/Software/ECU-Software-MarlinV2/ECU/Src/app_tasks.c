#include "app_tasks.h"

#include "cmsis_os2.h"

#include "actuators.h"
#include "autosequencer.h"
#include "control_state_machine.h"
#include "external_input.h"
#include "failure_monitor.h"
#include "sensors.h"
#include "telemetry.h"

#define INBOUND_QUEUE_LENGTH  16U
#define OUTBOUND_QUEUE_LENGTH 4U

static osMessageQueueId_t inbound_queue;  /* Mission Control -> control task. */
static osMessageQueueId_t outbound_queue; /* Control task -> telemetry task.  */

static void sync_clock(void)
{
  /* TODO: clock pulse sync at setup; frame timestamps come from this clock. */
}

void AppTasks_Initialize(void)
{
  Actuators_Initialize();

  inbound_queue = osMessageQueueNew(INBOUND_QUEUE_LENGTH, sizeof(ecu_inbound_frame_t), NULL);
  outbound_queue = osMessageQueueNew(OUTBOUND_QUEUE_LENGTH, sizeof(ecu_outbound_frame_t), NULL);

  Sensors_Initialize();
  ExternalInput_Initialize();
  Telemetry_Initialize();
  sync_clock();
  Sensors_StartDma();
}

bool AppTasks_SubmitInbound(const ecu_inbound_frame_t *frame)
{
  return osMessageQueuePut(inbound_queue, frame, 0U, 0U) == osOK;
}

void AppControlTask_Run(void *argument)
{
  control_state_machine_t state_machine;
  (void)argument;

  ControlStateMachine_Initialize(&state_machine,
                                 &autosequencer_main_sequence,
                                 &autosequencer_fault_sequence);

  for (;;)
  {
    /* TODO: control loop.  Outbound frame is filled from the synced clock,
     * state_machine.state, the sensor snapshot, output_command and the
     * actuator positions. */
    (void)osDelay(CONTROL_PERIOD_MS);
  }
}

void AppTelemetryTask_Run(void *argument)
{
  ecu_outbound_frame_t frame;

  (void)argument;
  for (;;)
  {
    if (osMessageQueueGet(outbound_queue, &frame, NULL, osWaitForever) == osOK)
    {
      Telemetry_Process(&frame);
    }
  }
}

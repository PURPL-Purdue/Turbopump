#include "app_tasks.h"

#include "cmsis_os2.h"
#include "stm32h7xx_hal.h"

#include "actuators.h"
#include "autosequencer.h"
#include "control_state_machine.h"
#include "external_input.h"
#include "failure_monitor.h"
#include "sensors.h"
#include "telemetry.h"

#define TELEMETRY_QUEUE_LENGTH 4U

static osMessageQueueId_t telemetry_queue;

void AppTasks_Initialize(void)
{
  /* Outputs first, so every valve is off as early as possible. */
  Actuators_Initialize();

  /*
   * Create the queue before the scheduler starts.  The control task will put
   * snapshots into it and the telemetry task will remove them later.
   */
  telemetry_queue = osMessageQueueNew(TELEMETRY_QUEUE_LENGTH, sizeof(telemetry_frame_t), NULL);

  /* Initialize each replaceable subsystem before the tasks are created. */
  Sensors_Initialize();
  ExternalInput_Initialize();
  Telemetry_Initialize();
  Sensors_StartDma();
}

void AppControlTask_Run(void *argument)
{
  control_state_machine_t state_machine;
  sensor_data_t sensors = {0};
  mc_command_t mc_command;
  telemetry_frame_t frame = {0};
  uint32_t next_control_release = osKernelGetTickCount();

  (void)argument;

  ControlStateMachine_Initialize(&state_machine,
                                 &autosequencer_main_sequence,
                                 &autosequencer_fault_sequence);

  for (;;)
  {
    /* Run once every CONTROL_PERIOD_MS, measured from a fixed start time. */
    next_control_release += CONTROL_PERIOD_MS;
    (void)osDelayUntil(next_control_release);

    /* Mission Control commands, at most one queue's worth per cycle. */
    for (uint32_t i = 0U; i < MC_COMMAND_QUEUE_LENGTH; i++)
    {
      if (!ExternalInput_Receive(&mc_command))
      {
        break;
      }
      ControlStateMachine_HandleCommand(&state_machine, &mc_command);
    }

    (void)Sensors_Read(&sensors);
    ControlStateMachine_Update(&state_machine,
                               FailureMonitor_Evaluate(&sensors, state_machine.state));
    Actuators_Apply(&state_machine.output_command);

    frame.timestamp_ms = osKernelGetTickCount();
    frame.state = state_machine.state;
    frame.sequence_tick = state_machine.sequence_tick;
    frame.sensors = sensors;
    frame.command = state_machine.output_command;
    frame.faults = state_machine.faults;

    /* Never wait on telemetry; if the queue is full this frame is skipped. */
    (void)osMessageQueuePut(telemetry_queue, &frame, 0U, 0U);
  }
}

void AppTelemetryTask_Run(void *argument)
{
  telemetry_frame_t frame;

  (void)argument;
  for (;;)
  {
    if (osMessageQueueGet(telemetry_queue, &frame, NULL, osWaitForever) == osOK)
    {
      Telemetry_Process(&frame);
    }
  }
}

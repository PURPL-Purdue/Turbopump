#include "app_tasks.h"

#include "FreeRTOS.h"
#include "queue.h"
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

/* Statically allocated queue memory; nothing comes from the RTOS heap. */
static StaticQueue_t inbound_queue_cb;
static uint8_t inbound_queue_storage[INBOUND_QUEUE_LENGTH * sizeof(ecu_inbound_frame_t)];
static StaticQueue_t outbound_queue_cb;
static uint8_t outbound_queue_storage[OUTBOUND_QUEUE_LENGTH * sizeof(ecu_outbound_frame_t)];

static const osMessageQueueAttr_t inbound_queue_attr = {
  .name = "inbound",
  .cb_mem = &inbound_queue_cb,
  .cb_size = sizeof(inbound_queue_cb),
  .mq_mem = inbound_queue_storage,
  .mq_size = sizeof(inbound_queue_storage),
};

static const osMessageQueueAttr_t outbound_queue_attr = {
  .name = "outbound",
  .cb_mem = &outbound_queue_cb,
  .cb_size = sizeof(outbound_queue_cb),
  .mq_mem = outbound_queue_storage,
  .mq_size = sizeof(outbound_queue_storage),
};

/* Static so it lives in .bss rather than on the control task's stack. */
static control_state_machine_t state_machine;

/* Control cycles that did not finish within CONTROL_PERIOD_MS. */
static uint32_t control_overrun_count;

/* Requests from Mission Control, decoded from this cycle's inbound frames. */
typedef struct
{
  bool clock_sync;
  bool start_sequence;
  bool clear_fault;
  /* TODO: further requests (manual valve commands, ...). */
} control_requests_t;

static void sync_clock(void)
{
  /* TODO: clock pulse sync; frame timestamps come from this clock. */
}

/* Step 1: drain the inbound queue into this cycle's requests. */
static void read_commands(control_requests_t *requests)
{
  ecu_inbound_frame_t frame;

  *requests = (control_requests_t){0};
  while (osMessageQueueGet(inbound_queue, &frame, NULL, 0U) == osOK)
  {
    /* TODO: decode frame into requests. */
  }
}

/* Step 3: copy the sensor snapshot into the outbound frame. */
static void frame_set_sensors(ecu_outbound_frame_t *frame,
                              const sensor_data_t *sensors)
{
  frame->timestamp_ms = sensors->timestamp_ms;
  for (uint32_t i = 0U; i < SNS_COUNT; i++)
  {
    frame->sensor_value[i] = sensors->value[i];
  }
}

/* Step 6: copy what is about to be commanded into the outbound frame. */
static void frame_set_actuators(ecu_outbound_frame_t *frame,
                                const control_state_machine_t *machine)
{
  frame->state = machine->state;
  for (uint32_t i = 0U; i < ACT_COUNT; i++)
  {
    frame->actuator_value[i] = machine->output_command.out[i];
  }
}

void AppTasks_Initialize(void)
{
  Actuators_Initialize();

  inbound_queue = osMessageQueueNew(INBOUND_QUEUE_LENGTH, sizeof(ecu_inbound_frame_t),
                                    &inbound_queue_attr);
  outbound_queue = osMessageQueueNew(OUTBOUND_QUEUE_LENGTH, sizeof(ecu_outbound_frame_t),
                                     &outbound_queue_attr);

  Sensors_Initialize();
  ExternalInput_Initialize();
  Telemetry_Initialize();
  Sensors_StartDma();
}

bool AppTasks_SubmitInbound(const ecu_inbound_frame_t *frame)
{
  return osMessageQueuePut(inbound_queue, frame, 0U, 0U) == osOK;
}

void AppControlTask_Run(void *argument)
{
  control_requests_t requests;
  sensor_data_t sensors = {0};
  ecu_outbound_frame_t frame = {0};

  (void)argument;

  ControlStateMachine_Initialize(&state_machine, &autosequencer_sequence);

  /* Absolute wake time: each cycle starts on its own tick, independent of
   * how long the previous cycle's work took. */
  uint32_t next_wake = osKernelGetTickCount();
  uint32_t last_telemetry_tick = next_wake;

  for (;;)
  {
    next_wake += CONTROL_PERIOD_MS;

    /* 1. Read commands from Mission Control. */
    read_commands(&requests);

    /* 2. Clock sync, only if requested and in STANDBY. */
    if (requests.clock_sync && (state_machine.state == ECU_STATE_STANDBY))
    {
      sync_clock();
    }

    /* 3. Read sensors into the local snapshot and frame.
     *    TODO: use the read result and per-sensor valid flags (see
     *    FailureMonitor_Evaluate for the bad-reading -> FAULT policy). */
    (void)Sensors_Read(&sensors);
    frame_set_sensors(&frame, &sensors);

    /* 4. Transitions: start request, clear-fault request, then fault logic.
     *    Start before clear, so FAULT cannot go to STARTUP in one cycle;
     *    fault logic last, so a fault always wins and a fault that is still
     *    present re-enters FAULT straight after a clear. */
    if (requests.start_sequence)
    {
      ControlStateMachine_StartSequence(&state_machine);
    }
    if (requests.clear_fault)
    {
      ControlStateMachine_ClearFault(&state_machine);
    }
    if (FailureMonitor_Evaluate(&sensors, state_machine.state))
    {
      ControlStateMachine_Fault(&state_machine);
    }

    /* 5. Run the active state (exit conditions checked against this cycle's
     *    snapshot) and set output_command. */
    ControlStateMachine_Update(&state_machine, &sensors);

    /* 6. Record what is about to be commanded. */
    frame_set_actuators(&frame, &state_machine);

    /* 7. Apply the actuator commands. */
    Actuators_Apply(&state_machine.output_command);

    /* 8-9. Send a frame to telemetry at most every TELEMETRY_PERIOD_MS.
     *      Never blocks; the frame is dropped if the queue is full. */
    if ((osKernelGetTickCount() - last_telemetry_tick) >= TELEMETRY_PERIOD_MS)
    {
      last_telemetry_tick = osKernelGetTickCount();
      (void)osMessageQueuePut(outbound_queue, &frame, 0U, 0U);
    }

    if ((int32_t)(osKernelGetTickCount() - next_wake) >= 0)
    {
      /* Work ran into the next cycle.  Drop the missed cycles and realign
       * rather than running back-to-back to catch up. */
      control_overrun_count++;
      /* TODO: overrun -> state machine fault (policy TBD). */
      next_wake = osKernelGetTickCount();
    }
    else
    {
      (void)osDelayUntil(next_wake);
    }
  }
}

void AppTelemetryTask_Run(void *argument)
{
  ecu_outbound_frame_t frame;

  (void)argument;
  for (;;)
  {
    /* Runs once per tick, after the higher-priority control task's cycle. */
    ExternalInput_Poll();

    while (osMessageQueueGet(outbound_queue, &frame, NULL, 0U) == osOK)
    {
      Telemetry_Process(&frame);
    }

    /* TODO: network stack timers (TCP retransmit, keepalive). */

    (void)osDelay(1U);
  }
}

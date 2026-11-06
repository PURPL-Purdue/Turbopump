#include <ESP32SPISlave.h>

ESP32SPISlave slave; //slave object

uint8_t recieve_buffer[1] = {0}; //buffer to recieve from stm32
uint8_t send_buffer[1] = {0}; //buffer send to stm32

void setup() {
  // put your setup code here, to run once:

    pinMode(22,OUTPUT);

    Serial.begin(115200);

    //must match SPI setup in Stm32, CPOL = Low, CPHA = First Edge
    slave.setDataMode(SPI_MODE0);

    slave.begin(VSPI,18,19,23,5);
}

void loop() {
  // put your main code here, to run r3epeatedly:
  
    slave.queue(send_buffer, recieve_buffer, 1);

    slave.wait();

    Serial.print("Recieved: ");

    if(recieve_buffer[0] == 1)
    {
      Serial.println(recieve_buffer[0]);

      digitalWrite(22,HIGH);
    }
    else if(recieve_buffer[0] == 0)
    {
      Serial.println(recieve_buffer[0]);

      digitalWrite(22, LOW);
    }

}

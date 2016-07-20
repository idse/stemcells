/*
   hallSensor

   Measures the analog outputs of the two halls sensors of the
   precision syringe.
*/

//Define pins for hall Sensor
#define front A0
#define back A5

//Declare variables/constants
float ana;
int readInt = 100;

void setup() {
  Serial.begin(9600);
  pinMode(A0, INPUT);
  pinMode(A1, INPUT);
}

void loop() {
  Serial.print(analogRead(front));
  Serial.print(" ");
  Serial.println(analogRead(back));
  delay(readInt);
}

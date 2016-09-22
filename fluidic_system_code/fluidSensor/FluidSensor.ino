/*
 * fluidSensor
 * 
 * Measures the resistance of a fluid while alternating the current
 * to prevent electrolysis. A measurement is taken every time the current
 * is switched. 
*/


//Define pins for fluid sensor. Please see the online notes for documentation
#define Vin1 13
#define Vin2 12
#define Vout1 A0
#define Vout2 A1

//Declare fluid sensor variables
float anaVolt;      //analog input from voltage divider
int anaRes = 1023;  //resolution of analog input
int R1 = 30;        //constant resistor 1 (in kOhms)
int R2 = 30;        //constant resistor 2 (in kOhms)
int readInt = 100; //milliseconds between each reading

void setup() {
  Serial.begin(9600);
  pinMode(Vin1, OUTPUT);
  pinMode(Vin2, OUTPUT);
  pinMode(Vout1, INPUT);
  pinMode(Vout2, INPUT);
  //analogReadResolution(10);
}

void loop() {
  //FORWARD CURRENT
  forwardCur();
  //Serial.println("Forward Current:");
  readFluidSensor(Vout2, 1);
  //BACKWARD CURRENT
  backwardCur();
  //Serial.println("Backward Current:");
  readFluidSensor(Vout1, 2);
}

void forwardCur() {
  digitalWrite(Vin2, LOW);
  delay(1);
  digitalWrite(Vin1, HIGH);
  delay(readInt);
}

void backwardCur() {
  digitalWrite(Vin1, LOW);
  delay(1);
  digitalWrite(Vin2, HIGH);
  delay(readInt);
}

int readFluidSensor(int VoutPin, int curState) {
  anaVolt = analogRead(VoutPin);
  float Rf; //resistance of fluid
  if (curState == 1)
    Rf = (anaRes / anaVolt) * R2 - R1 - R2;
  else
    Rf = (anaRes / anaVolt) * R2 - R1 - R2;
  Serial.print(anaVolt);
  Serial.println(Rf);
  //Serial.println();
  return Rf;
}


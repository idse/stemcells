/*
   MotorDrive sketch
   This codes expects a message in the format: 1,4567,890*
   The first field is an indicator of 1 or 0 for direction
   The rest of the fields can be arbitrary length
   This code requires some arbitrary character
   to indicate the end of the command. Here I have used *
*/

//Define pins for motor 1
#define stp_1 2
#define dir_1 3
#define MS1_1 4
#define MS2_1 5
#define MS3_1 6
#define EN_1  7
//Define pins for motor 2
#define stp_2 8
#define dir_2 9
#define MS1_2 10
#define MS2_2 11
#define MS3_2 12
#define EN_2  13

//Declare variables for serial reading
char ui;                        //user serial input
const int NUMBER_OF_FIELDS = 4; //expected # of fields
int fieldIndex = 0;             //to iterate through fields
int values[NUMBER_OF_FIELDS];   //store our fields here
int done = 0;                   //indicator variable

//Declare motor drive variables
long xlim;  //indicate volume (in # of steps)
long spd_t; //delay time between steps (for controlling speed)

void setup() {
  pinMode(stp_1, OUTPUT);pinMode(stp_2, OUTPUT);
  pinMode(dir_1, OUTPUT);pinMode(dir_2, OUTPUT);
  pinMode(MS1_1, OUTPUT);pinMode(MS1_2, OUTPUT);
  pinMode(MS2_1, OUTPUT);pinMode(MS2_2, OUTPUT);
  pinMode(MS3_1, OUTPUT);pinMode(MS3_2, OUTPUT);
  pinMode(EN_1, OUTPUT);pinMode(EN_2, OUTPUT);
  resetBEDPins();     //Set step, direction, microstep and enable pins to default states
  Serial.begin(9600); //Open Serial connection for debugging
  Serial.println("Begin fluid control");
  Serial.println();
  //Print control instructions
  Serial.println("Enter parameters");
  Serial.println("Usage: <motor><direction>,<microLiters>,<microLiters/sec><*>");
  Serial.println("Example: 1,1,1000,100*");
  Serial.println("Where direction is given by 1 or 0 for in or out, respectively.");
  Serial.println();
}

//Main loop
void loop() {
  //loop to acquire the fields
  while (done == 0) {
    if (Serial.available()) {
      ui = Serial.read();
      if (ui >= '0' && ui <= '9') //check ascii value
        values[fieldIndex] = (values[fieldIndex] * 10) + (ui - '0');
      else if (ui == ',') {
        if (fieldIndex < NUMBER_OF_FIELDS - 1)
          fieldIndex++;
      }
      else
      {
        //field acquisition is over
        done = 1; //indicator
      }
    }
  }
  done = 0;
  Serial.print("Motor: "); Serial.println(values[0]);
  Serial.print("Direction: "); Serial.println(values[1]);
  Serial.print("Volume: "); Serial.println(values[2]);
  Serial.print("Speed: "); Serial.println(values[3]);
  Serial.println();
  if (values[0] == 1) {
    digitalWrite(EN_1, LOW); //Pull enable pin low to set FETs active and allow motor control
    SmallStepMode1();
  }
  else if (values[0] == 2) {
    digitalWrite(EN_2, LOW); //Pull enable pin low to set FETs active and allow motor control
    SmallStepMode2();
  }
  else
    Serial.print("Not a valid motor");
  for (int i = 0; i <= fieldIndex; i++)
    values[i] = 0; //reset for next message
  fieldIndex = 0;  //reset for next message
  resetBEDPins();
}

//Reset Big Easy Driver pins to default states
void resetBEDPins()
{
  digitalWrite(stp_1, LOW); digitalWrite(stp_2, LOW);
  digitalWrite(dir_1, LOW); digitalWrite(dir_2, LOW);
  digitalWrite(MS1_1, LOW); digitalWrite(MS1_2, LOW);
  digitalWrite(MS2_1, LOW); digitalWrite(MS2_2, LOW);
  digitalWrite(MS3_1, LOW); digitalWrite(MS3_2, LOW);
  digitalWrite(EN_1, HIGH); digitalWrite(EN_2, HIGH);
}

// Drive motor 1 @ 1/16th microstep
void SmallStepMode1()
{
  digitalWrite(MS1_1, HIGH); //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
  digitalWrite(MS2_1, HIGH);
  digitalWrite(MS3_1, HIGH);
  //DIRECTION
  if (values[1] == 1)
    digitalWrite(dir_1, HIGH); //Pull direction pin high to inject
  else
    digitalWrite(dir_1, LOW); //Pull direction pin low to remove
  //VOLUME
  xlim = ceil(6.997333 * values[2]);
  Serial.println("Steps: "); Serial.println(xlim);
  //SPEED
  spd_t = ceil(71.45579 / values[3] * 1000);
  Serial.println("Delay: "); Serial.println(spd_t);
  Serial.println();
  //DRIVE
  if (spd_t <= 16383) //microsecond delay is only accurate for such values
  {
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_1, HIGH); //Trigger one step forward
      delayMicroseconds(spd_t);
      digitalWrite(stp_1, LOW);  //Pull step pin low so it can be triggered again
      delayMicroseconds(spd_t);
    }
  }
  else //must use millisecond delay
  {
    spd_t = floor(spd_t / 1000);  //convert to milliseconds
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_1, HIGH); //Trigger one step forward
      delay(spd_t);
      digitalWrite(stp_1, LOW);  //Pull step pin low so it can be triggered again
      delay(spd_t);
    }
  }
}

// Drive motor 2 @ 1/16th microstep
void SmallStepMode2()
{
  digitalWrite(MS1_2, HIGH); //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
  digitalWrite(MS2_2, HIGH);
  digitalWrite(MS3_2, HIGH);
  //DIRECTION
  if (values[1] == 1)
    digitalWrite(dir_2, HIGH); //Pull direction pin high to inject
  else
    digitalWrite(dir_2, LOW); //Pull direction pin low to remove
  //VOLUME
  xlim = ceil(6.997333 * values[2]);
  Serial.println("Steps: "); Serial.println(xlim);
  //SPEED
  spd_t = ceil(71.45579 / values[3] * 1000);
  Serial.println("Delay: "); Serial.println(spd_t);
  Serial.println();
  //DRIVE
  if (spd_t <= 16383) //microsecond delay is only accurate for such values
  {
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_2, HIGH); //Trigger one step forward
      delayMicroseconds(spd_t);
      digitalWrite(stp_2, LOW);  //Pull step pin low so it can be triggered again
      delayMicroseconds(spd_t);
    }
  }
  else //must use millisecond delay
  {
    spd_t = floor(spd_t / 1000);  //convert to milliseconds
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_2, HIGH); //Trigger one step forward
      delay(spd_t);
      digitalWrite(stp_2, LOW);  //Pull step pin low so it can be triggered again
      delay(spd_t);
    }
  }
}

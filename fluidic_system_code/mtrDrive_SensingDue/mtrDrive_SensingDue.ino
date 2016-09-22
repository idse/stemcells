/*
   mtrDrive_SensingDue sketch
   Drives syringe pumps using arduino Due, with the motors pins in parallel.
   Stops third syringe with limit switch if motor overextends either way.
   Connects second and third syringe with 3-way valve.

   This codes expects a message in the format: 1,4567,890*
   <motor>,<direction>,<microLiters>,<microLiters/sec><*>
   The fields can be of arbitrary length. They must be integers.
   This code requires some arbitrary character
   to indicate the end of the command. Here I have used "*"

   ************************************************************
   "Motor 4" actually represents the case where motor 1 and 2 are
   pumping in opposite directions, to produce a fluid flow. I didn't
   have time to implement a library like AccelStepper, which allows for
   simultaneous motor driving. This solution is temporary and should be
   replaced.
   ************************************************************
*/
//#include <Scheduler.h>

//Define how many motors we are driving
#define numMtrs 3
//Define parallel pins for all motors
#define MS1 53
#define MS2 51
#define MS3 49
//Define pins for motor 1
#define EN_1  52
#define stp_1 50
#define dir_1 48
//Define pins for motor 2
#define EN_2  46
#define stp_2 44
#define dir_2 42
//Define pins for motor 3
#define EN_3  40
#define stp_3 38
#define dir_3 36

//Define pins for limit switch
#define LSout 22
#define LSin 24

//Define pins for 3-way solenoid valve
#define valvepin 13 //actually it controls the relay

//Define pins for fluid sensor



//Declare motor pins array
int pins[numMtrs][6] = {
  {EN_1, MS1, MS2, MS3, stp_1, dir_1},
  {EN_2, MS1, MS2, MS3, stp_2, dir_2},
  {EN_3, MS1, MS2, MS3, stp_3, dir_3}
};
int dpins[6]; //the pins of the currently driven motor

//Declare serial reading variables
char ui;                        //user serial input
const int NUMBER_OF_FIELDS = 4; //expected # of fields
int fieldIndex = 0;             //to iterate through fields
int values[NUMBER_OF_FIELDS];   //store our fields here
int done = 0;                   //indicator variable

//Declare motor drive variables
long xlim;  //indicate volume (in # of steps)
long spd_t; //delay time between steps (for controlling speed)

//Declare limit switch variables
volatile boolean stopmtr = LOW;

//Declare basic syringe constants (an underscore indicates "per")
float mm_turn = 1.25;
float mcl_mm = (30.0 / 82) * 1000;
int stp_turn = 3200;
float stp_mcl = stp_turn / mcl_mm / mm_turn; // ~7
float mcs_stp = 1000000 / stp_mcl;
float errorCorrection = 1.025746;
//Declare precision syringe constants (an underscore indicates "per")
float mm_turn_PS = 5.00;
float mcl_mm_PS = (1.0 / 60) * 1000;
int stp_turn_PS = 3200;
float stp_mcl_PS = stp_turn_PS / mcl_mm_PS / mm_turn_PS; // 38.4
float mcs_stp_PS = 1000000 / stp_mcl_PS;





void setup() {
  //Set pinModes for motor pins:
  for (int i = 0; i < numMtrs; i++) {
    for (int j = 0; j < 6; j++) {
      pinMode(pins[i][j], OUTPUT);
    }
  }
  //Set pinModes for other pins:
  pinMode(LSout, OUTPUT);
  pinMode(LSin, INPUT);
  pinMode(valvepin, OUTPUT);

  //Set up interrupts:
  attachInterrupt(digitalPinToInterrupt(LSin), checkLS, CHANGE);

  //Prepare for motor driving:
  resetPins();        //Set all pins to default states
  Serial.begin(9600); //Open serial connection for displaying/debugging
  while (!Serial) {
    ; // wait for serial port to connect
  }

  //Print control instructions:
  Serial.print(
    "Begin fluid control\n\n"
    "Enter parameters\n"
    "Usage: <motor>,<direction>,<microLiters>,<microLiters/sec><*>\n"
    "Example: 1,1,1000,100*\n"
    "Where direction is given by 1 or 0 for in or out, respectively.\n\n");
  delay(1000);
}




//Main loop
void loop() {
  //Loop to acquire the serial fields
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
        //field acquisition is over
        done = 1; //indicator
    }
  }
  done = 0;

  //Diagnostic print statements
  /*
    Serial.print("Motor: "); Serial.println(values[0]);
    Serial.print("Direction: "); Serial.println(values[1]);
    Serial.print("Volume: "); Serial.println(values[2]);
    Serial.print("Speed: "); Serial.println(values[3]);
    Serial.println();
  */

  //Find motor pins
  if (values[0] != 4) { //4 is not a real motor
    for (int i = 0; i < 6; i++) {
      dpins[i] = pins[values[0] - 1][i];
    }
  }

  //Send command
  if (values[0] == 3)
    PrecisionSyringe();
  else if (values[0] == 4)
    BasicSyringeFlow();
  else
    BasicSyringe();

  //Reset for next command
  for (int i = 0; i <= fieldIndex; i++)
    values[i] = 0;
  fieldIndex = 0;
  resetPins();
  delay(100);
}




//Reset Big Easy Driver pins to default states
void resetPins()
{
  for (int i = 0; i < numMtrs; i++) {
    for (int j = 0; j < 6; j++) {
      if (j == 0)
        digitalWrite(pins[i][j], HIGH);
      else
        digitalWrite(pins[i][j], LOW);
    }
  }
  digitalWrite(LSout, LOW);
  digitalWrite(valvepin, LOW);
}




// Drive blank and waste syringes
void BasicSyringe()
{
  //Pull EN down to enable FETs and allow motor to move
  digitalWrite(dpins[0], LOW);
  //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
  digitalWrite(dpins[1], HIGH); digitalWrite(dpins[2], HIGH); digitalWrite(dpins[3], HIGH);

  //DIRECTION
  if (values[1] == 1)
    digitalWrite(dpins[5], HIGH); //Pull direction pin high to inject
  else
    digitalWrite(dpins[5], LOW); //Pull direction pin low to remove
  //VOLUME
  xlim = ceil(values[2] * stp_mcl * errorCorrection);
  Serial.print("Steps: "); Serial.println(xlim);
  //SPEED
  spd_t = ceil(mcs_stp / values[3] / 2);
  Serial.print("Delay: "); Serial.println(spd_t); Serial.println();

  //DRIVE
  if (spd_t <= 16383) //microsecond delay is only accurate for such values
  {
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(dpins[4], HIGH); //Trigger one step forward
      delayMicroseconds(spd_t);
      digitalWrite(dpins[4], LOW);  //Pull step pin low so it can be triggered again
      delayMicroseconds(spd_t);
    }
  }
  else //must use millisecond delay
  {
    spd_t = ceil(spd_t / 1000.0);  //convert to milliseconds
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(dpins[4], HIGH); //Trigger one step forward
      delay(spd_t);
      digitalWrite(dpins[4], LOW);  //Pull step pin low so it can be triggered again
      delay(spd_t);
    }
  }
}




// Drive concentrated ligand syringe
void PrecisionSyringe()
{
  //Enable limit switch and open 3-way valve
  digitalWrite(LSout, HIGH);
  digitalWrite(valvepin, HIGH); //valve is normally closed, now we open

  //Pull EN down to enable FETs and allow motor to move
  digitalWrite(dpins[0], LOW);
  //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
  digitalWrite(dpins[1], HIGH); digitalWrite(dpins[2], HIGH); digitalWrite(dpins[3], HIGH);

  //DIRECTION
  if (values[1] == 1)
    digitalWrite(dpins[5], HIGH); //Pull direction pin high to inject
  else
    digitalWrite(dpins[5], LOW); //Pull direction pin low to remove
  //VOLUME
  xlim = ceil(values[2] * stp_mcl_PS);
  Serial.print("Steps: "); Serial.println(xlim);
  //SPEED
  spd_t = ceil(mcs_stp_PS / values[3] / 2);
  Serial.print("Delay: "); Serial.println(spd_t);

  //DRIVE
  if (spd_t <= 16383) //microsecond delay is only accurate for such values
  {
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      if (stopmtr) {
        Serial.println("LIMIT REACHED");
        break;
      }
      digitalWrite(dpins[4], HIGH); //Trigger one step forward
      delayMicroseconds(spd_t);
      digitalWrite(dpins[4], LOW);  //Pull step pin low so it can be triggered again
      delayMicroseconds(spd_t);
    }
  }
  else //must use millisecond delay
  {
    spd_t = ceil(spd_t / 1000.0);  //convert to milliseconds
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      if (stopmtr) {
        Serial.println("LIMIT REACHED");
        break;
      }
      digitalWrite(dpins[4], HIGH); //Trigger one step forward
      delay(spd_t);
      digitalWrite(dpins[4], LOW);  //Pull step pin low so it can be triggered again
      delay(spd_t);
    }
  }
  Serial.println();
}




// Drive blank and waste syringes at same time
void BasicSyringeFlow()
{
  //Pull EN down to enable FETs and allow motor to move
  digitalWrite(EN_1, LOW); digitalWrite(EN_2, LOW);
  //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
  digitalWrite(MS1, HIGH); digitalWrite(MS2, HIGH); digitalWrite(MS3, HIGH);

  //DIRECTION
  if (values[1] == 1) {
    digitalWrite(dir_1, HIGH); //Pull direction pin high to inject
    digitalWrite(dir_2, LOW);  //Pull direction pin low to remove
  }
  else {
    digitalWrite(dir_1, LOW); //Pull direction pin high to inject
    digitalWrite(dir_2, HIGH);  //Pull direction pin low to remove
  }
  //VOLUME
  xlim = ceil(values[2] * stp_mcl * errorCorrection);
  Serial.print("Steps: "); Serial.println(xlim);
  //SPEED
  spd_t = ceil(mcs_stp / values[3] / 2);
  Serial.print("Delay: "); Serial.println(spd_t); Serial.println();

  //DRIVE
  if (spd_t <= 16383) //microsecond delay is only accurate for such values
  {
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_1, HIGH); //Trigger one step forward
      digitalWrite(stp_2, HIGH); //Trigger one step forward
      delayMicroseconds(spd_t);
      digitalWrite(stp_1, LOW);  //Pull step pin low so it can be triggered again
      digitalWrite(stp_2, LOW);  //Pull step pin low so it can be triggered again
      delayMicroseconds(spd_t);
    }
  }
  else //must use millisecond delay
  {
    spd_t = ceil(spd_t / 1000.0);  //convert to milliseconds
    for (int x = 1; x <= xlim; x++) //Loop the forward stepping enough times for motion to be visible
    {
      digitalWrite(stp_1, HIGH); //Trigger one step forward
      digitalWrite(stp_2, HIGH); //Trigger one step forward
      delay(spd_t);
      digitalWrite(stp_1, LOW);  //Pull step pin low so it can be triggered again
      digitalWrite(stp_2, LOW);  //Pull step pin low so it can be triggered again
      delay(spd_t);
    }
  }
}



//ISR to stop the Precision Syringe after limit switch
void checkLS()
{
  stopmtr = digitalRead(LSin);
  stopmtr = digitalRead(LSin); //must read twice. Reading once is unreliable.
}




#define valvepin 13 //actually it controls the relay



void setup() {
  pinMode(valvepin, OUTPUT);
  digitalWrite(valvepin, HIGH);
}

void loop() {
  delay(10000);
}

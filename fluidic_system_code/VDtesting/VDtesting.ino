float anavolt;

void setup() {
  Serial.begin(9600);
  pinMode(0, INPUT);
}

void loop() {
  anavolt = analogRead(0);
  Serial.println(
    String("Voltage: ") + String(analogRead(0)) + 
    String(" Resistance: ") + String(floor(60*(anavolt/1024)/(1-(anavolt/1024))))
    );
  delay(100);
}

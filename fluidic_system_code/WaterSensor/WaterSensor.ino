float anavolt;

void setup() {
  Serial.begin(9600);
  pinMode(A0, INPUT);
}

void loop() {
  analogReadResolution(12);
  anavolt = analogRead(A0);
  Serial.println(
    String("Voltage: ") + String(analogRead(A0)) + 
    String(" Resistance: ") + String(floor(60*(anavolt/4096)/(1-(anavolt/4096))))
    );
  delay(100);
}

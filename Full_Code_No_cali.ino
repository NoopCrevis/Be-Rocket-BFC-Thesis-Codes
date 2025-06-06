#include <Wire.h>
#include <SPI.h>
#include <Adafruit_MAX31865.h>
#include "SparkFun_Qwiic_Scale_NAU7802_Arduino_Library.h"
#include "HX711.h"

// --- Sensor Pins ---
const int pressurePins[4] = {0, 1, 2, 3}; // ADC pins
#define HX711_DATA 6
#define HX711_CLK  7
const int potPin = 4; // Potentiometer analog pin

// --- MAX31865 Pins ---
#define MAX31865_CS1  18
#define MAX31865_CS2  19
#define MAX31865_CS3  15
#define MAX31865_CS4  13
#define MAX31865_MOSI 12
#define MAX31865_MISO 11
#define MAX31865_SCK  10

// --- Pressure Sensor Calibration ---
const float dividerRatio = 12.0 / (12.0 + 4.7); // voltage divider ratio
const float ADC_VREF = 3.3;
const int ADC_RESOLUTION = 4095;
const float SENSOR_VMIN = 0.35;  // V
const float SENSOR_VMAX = 3.3;   // V
const float PRESSURE_MAX = 110.3; // bar (1600 psi)

// --- Globals ---
NAU7802 nau;
HX711 scale;

Adafruit_MAX31865 pt100_1 = Adafruit_MAX31865(MAX31865_CS1);
Adafruit_MAX31865 pt100_2 = Adafruit_MAX31865(MAX31865_CS2);
Adafruit_MAX31865 pt100_3 = Adafruit_MAX31865(MAX31865_CS3);
Adafruit_MAX31865 pt100_4 = Adafruit_MAX31865(MAX31865_CS4);

unsigned long lastPrintTime = 0;
const unsigned long intervalMs = 10;  // 100 Hz

void setup() {
  Serial.begin(230400); // Fast baud for logging
  Wire.begin(5, 8);      // SDA = GPIO 5, SCL = GPIO 8

  // --- SPI bus for MAX31865 ---
  SPI.begin(MAX31865_SCK, MAX31865_MISO, MAX31865_MOSI);

  pt100_1.begin(MAX31865_3WIRE);
  pt100_2.begin(MAX31865_3WIRE);
  pt100_3.begin(MAX31865_3WIRE);
  pt100_4.begin(MAX31865_3WIRE);

  // --- NAU7802 Torque Sensor ---
  if (!nau.begin()) {
    Serial.println("NAU7802 not detected");
    while (1);
  }
  nau.setGain(128);
  nau.calibrateAFE();

  // --- HX711 Load Cell ---
  scale.begin(HX711_DATA, HX711_CLK);
  if (!scale.is_ready()) {
    Serial.println("HX711 not found");
    while (1);
  }
  scale.set_scale();  // raw units
  scale.tare();

  // --- CSV Header ---
  Serial.println("Torque_raw,Load_raw,Pressure1_bar,Pressure2_bar,Pressure3_bar,Pressure4_bar,Temp1_C,Temp2_C,Temp3_C,Temp4_C,ValvePosition_V");
}

void loop() {
  unsigned long now = millis();
  if (now - lastPrintTime >= intervalMs) {
    lastPrintTime = now;

    // --- Torque Sensor ---
    long torque = nau.available() ? nau.getReading() : 0;

    // --- Load Cell ---
    long load = scale.get_units(1);  // minimal averaging for speed

    // --- Pressure Sensors ---
    float pressures[4];
    for (int i = 0; i < 4; i++) {
      int raw = analogRead(pressurePins[i]);
      float Vadc = (raw * ADC_VREF) / ADC_RESOLUTION;
      float Vsensor = Vadc / dividerRatio;
      float pressure = (Vsensor - SENSOR_VMIN) / (SENSOR_VMAX - SENSOR_VMIN) * PRESSURE_MAX;
      pressure = constrain(pressure, 0.0, PRESSURE_MAX);
      pressures[i] = pressure;
    }

    // --- Temperatures ---
    float t1 = pt100_1.temperature(100.0, 430.0);
    float t2 = pt100_2.temperature(100.0, 430.0);
    float t3 = pt100_3.temperature(100.0, 430.0);
    float t4 = pt100_4.temperature(100.0, 430.0);

    // --- Valve Position Potentiometer ---
    int potRaw = analogRead(potPin);
    float potVoltage = (potRaw * ADC_VREF) / ADC_RESOLUTION;

    // --- Output CSV ---
    Serial.print(torque); Serial.print(",");
    Serial.print(load); Serial.print(",");
    for (int i = 0; i < 4; i++) {
      Serial.print(pressures[i], 2); Serial.print(",");
    }
    Serial.print(t1, 2); Serial.print(",");
    Serial.print(t2, 2); Serial.print(",");
    Serial.print(t3, 2); Serial.print(",");
    Serial.print(t4, 2); Serial.print(",");
    Serial.println(potVoltage, 3);
  }
}

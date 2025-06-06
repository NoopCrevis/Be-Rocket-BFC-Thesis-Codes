#include <Encoder.h> /*Library for the encoder*/
#include <Servo.h> /*Library for servo, this is used for the PWM for ESC*/
#include "U8glib.h" /*Library for OLED display*/

#define ENCODER_A_PIN 2 /*Defining pin 2 as encoder pin A*/
#define ENCODER_B_PIN 3 /*Defining pin 3 as encoder pin B*/
#define ESC_PWM_PIN 9  /*Defining pin 9 as the ouput pin for the PWM signal*/

Encoder myEncoder(ENCODER_A_PIN, ENCODER_B_PIN); // "Encoder" is a class from the "Encoder.h" library, "myEncoder" is the variable of that class, and "ENCODER_A_PIN and B_PIN are the digital pins the encoder is connected to
Servo esc; // "Servo" is the class from the "Servo.h" library, "esc" is the variable of that class

long position = 0; //"long" dat type for large integers, "position" is the variable, this is also set to zero
int pwmValue = 0; //"int" also an integer but of the shorter type (no decimal point), the variable pwmValue is set to zero

// OLED setup
U8GLIB_SSD1306_128X64 u8g(U8G_I2C_OPT_NONE); //"U8GLIB_SSD1306_128X64" Class from the library for 128x63 pixel OLED, "u8g" name of the object for later use, "U8G_I2C_OPT_NONE" Constructor argument so I2C is used (SDA & SCL)

// Define RPM range for mapping
const int minPWM = 45; //Constant integer variable for the minimum PWM
const int maxPWM = 135; //Constant integer variable for maximum PWM
const int minRPM = 0; //Constant integer variable minimum RPM
const int maxRPM = 2500; //Constant integer variable maximum RPM

void drawRPM(int rpm) {
  char rpm_str[20]; //Creates a string (character array) big enough to hold the formatted text
  sprintf(rpm_str, "RPM: %d", rpm); //Formats the rpm integer into a string
  
  u8g.setFont(u8g_font_fub14);  //Sets the font to a bold, 14-point font
  u8g.drawStr(10, 25, rpm_str); //Draws the formatted string on the screen at coordinates (x=10, y=25).
}

void setup() {
  Serial.begin(9600); //For serial communication and serial monitor, ...
  esc.attach(ESC_PWM_PIN); //Attaches ESC to the defined pin (9) 
}
\
void loop() {
  // Read encoder
  position = myEncoder.read(); //Reading of the encoder position
  
  // Map encoder to PWM range
  pwmValue = map(position, -250, 250, minPWM, maxPWM); //"map" rescales the -250 to minPWM and 250 to maxPWM
  pwmValue = constrain(pwmValue, minPWM, maxPWM); //Makes sure the pwmValue will not go under or over the minPWM and maxPWM
  
  // Send to ESC
  esc.write(pwmValue); //Send the pwmValue to the ESC PWM pin

  // Estimate RPM based on PWM (linear assumption)
  int estimatedRPM = map(pwmValue, 90, maxPWM, 0, maxRPM); //"pwmValue" is current PWM signal, "90" and "maxPWM" PWM range, same for RPM, "map" rescales pwmValue into estimated RPM range
  
  // Display on OLED
  u8g.firstPage(); //Sets up the first page on OLED
  do { //Ensures drawRPM function runs to render RPM text
    drawRPM(estimatedRPM); //Text to display on OLED
  } while (u8g.nextPage()); //Prepares next page to render (avioding flickering)

  // Serial output for debugging
  Serial.print("PWM Value: ");
  Serial.print(pwmValue);
  Serial.print(" -> Estimated RPM: ");
  Serial.println(estimatedRPM);

  delay(100);  // Display refresh rate
}

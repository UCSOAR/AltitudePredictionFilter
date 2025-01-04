# AltitudePredictionFilter

Uses Madgwick Filter for the IMU which is then fused with the Barometers under Everest. The Everest output is fed into HALO which encapsulates the Unscented Kalman Filter and fuses the GPS also with the system.

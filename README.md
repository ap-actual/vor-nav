# vor-nav

## Introduction
It is no secret that GPS is vulnerable, and the United States as a whole is too reliant on its service [1]. One particular weak point is commercial aircraft’s inability to identify a spoofing attack and reject bad GPS measurements [2].

This project will build on an idea explored in a project during a previous course, “Modern Navigation Systems”. The previous project analyzed the amount of aviation navigation aids around the US and explored the efficacy of using them as an onboard consistency check against an aircraft’s GPS solution. Devices like DME’s and VOR’s are sources of terrestrial-based radio navigation scattered around the US; they are typically used standalone when GPS is not available or for use in IFR navigation [3]. 

<img width="562" height="218" alt="image" src="https://github.com/user-attachments/assets/82aacde5-83da-4087-8cb1-778d6222bed0" />

As there are many VOR’s scattered across the country with diverse geometry, these VOR’s can be used as the measurement update in a Kalman Filter coupled with an inexpensive MEMS accelerometer. Below is one example flight from BWI to ORD. The simulation displays the expected encounter lines from visible VOR’s based on the aircraft’s latitude, longitude, and altitude as well as their uncertainty cone (2.5 degrees). The location information on these VOR’s can be pulled from an FAA database. [4][5]

<img width="643" height="359" alt="BWI_to_ORD_VOR" src="https://github.com/user-attachments/assets/a1eae7c8-d62f-4c56-83bc-9793729dae46" />

## Project Overview
This project will mechanize an Extended Kalman Filter (EKF) to utilize 

## Datasets Used
### Aircraft Trajectories
Aircraft trajectories were created using data from opensky-network.org. Data was pulled from the open-source database and filtered for relevant commercial aircraft flights. 

<img width="738" height="490" alt="image" src="https://github.com/user-attachments/assets/5e56bcc9-3750-4078-a1b5-536df6115d81" />


### VOR Data
Navaid data for the entire continental United States was pulled from the FAA website at https://www.faa.gov/airtraffic/rco-vor-master-list.

<img width="788" height="461" alt="image" src="https://github.com/user-attachments/assets/abdfe919-68a2-4815-96d5-6d2fb7ecfe9b" />

### VOR Measurement Model
<img width="960" height="720" alt="VOR Measurement Model (1)" src="https://github.com/user-attachments/assets/5a6ed984-e55c-466d-9215-f78262b91298" />

> Link to draw file: https://docs.google.com/drawings/d/16Xxh9EqZsI2CmBNIY7Cf3i_a1hbbCS4FhxrkBuhYRwA/edit?usp=sharing

## References

[1]https://spacenews.com/the-race-to-back-up-vulnerable-gps/ 

[2]https://www.forbes.com/sites/erictegler/2023/09/28/someone-in-the-middle-east-is-leading-aircraft-astray-by-spoofing-gps-signals/ 

[3]https://en.wikipedia.org/wiki/VOR/DME#:~:text=In%20radio%20navigation%2C%20a%20VOR,the%20receiver%20and%20the%20station 

[4]https://adds-faa.opendata.arcgis.com/datasets/990e238991b44dd08af27d7b43e70b92_0/explore?location=5.144894%2C-1.658576%2C2.00 

[5]https://epicflightacademy.com/what-is-vor/ 

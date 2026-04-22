# vor-nav

## Introduction
It is no secret that GPS is vulnerable, and the United States as a whole is too reliant on its service [1]. One particular weak point is commercial aircraft’s inability to identify a spoofing attack and reject bad GPS measurements [2].

This project will build on an idea explored in a project during a previous course, “Modern Navigation Systems”. The previous project analyzed the amount of aviation navigation aids around the US and explored the efficacy of using them as an onboard consistency check against an aircraft’s GPS solution. Devices like DME’s and VOR’s are sources of terrestrial-based radio navigation scattered around the US; they are typically used standalone when GPS is not available or for use in IFR navigation [3]. 

<img width="562" height="218" alt="image" src="https://github.com/user-attachments/assets/82aacde5-83da-4087-8cb1-778d6222bed0" />

As there are many VOR’s scattered across the country with diverse geometry, these VOR’s can be used as the measurement update in a Kalman Filter coupled with an inexpensive MEMS accelerometer. Below is one example flight from BWI to ORD. The simulation displays the expected encounter lines from visible VOR’s based on the aircraft’s latitude, longitude, and altitude as well as their uncertainty cone (2.5 degrees). The location information on these VOR’s can be pulled from an FAA database. [4][5]

<img width="643" height="359" alt="BWI_to_ORD_VOR" src="https://github.com/user-attachments/assets/a1eae7c8-d62f-4c56-83bc-9793729dae46" />

## Project Overview
This project will
1. [Collect a dataset](#datasets-used) of representative aircraft trajectories and FAA terrestrial navigation aids (NAVAIDS)
2. [Generate simulated IMU and VOR measurements](#simulated-measurement-data) along a given flight path from the dataset generated above
3. Fuse the simulated measurements using an [Extended Kalman Filter](#the-extended-kalman-filter)

## Datasets Used
### Aircraft Trajectories
Aircraft trajectories were created using data from opensky-network.org. Data was pulled from the open-source database and filtered for relevant commercial aircraft flights. 

<img width="738" height="490" alt="image" src="https://github.com/user-attachments/assets/5e56bcc9-3750-4078-a1b5-536df6115d81" />


### VOR Data
Navaid data for the entire continental United States was pulled from the FAA website at https://www.faa.gov/airtraffic/rco-vor-master-list.

<img width="788" height="461" alt="image" src="https://github.com/user-attachments/assets/abdfe919-68a2-4815-96d5-6d2fb7ecfe9b" />

### Generation of Simulated IMU Measurements
#### Data Processing and Interpolation
This portion describes the methodology used to generate simulated IMU measurements for EKF testing. The process consists of:

1. Processing real flight trajectory data  
2. Generating "truth" IMU measurements  
3. Injecting sensor errors based on IMU specifications  

Flight trajectory data was obtained from:
- https://opensky-network.org
- Below is a picture of the data contained for each opensky-network trajectory

<img width="1592" height="696" alt="image" src="https://github.com/user-attachments/assets/74ec5894-c41d-4c07-9f3f-132008ce963f" />

Each trajectory contains data sampled at **0.1 Hz (every 10 seconds)**:
- Time (s)
- Latitude (deg)
- Longitude (deg)
- Speed (knots)
- Heading (deg)
- Vertical rate (m/s)
- Geoaltitude (m)

The following filtering steps were applied to find suitable flights:
- Removed non-flight segments (e.g., ground idle)
- Removed flights with missing data (NaNs or dropouts)
- Ensured continuous trajectory segments only

Trajectory data was then upsampled to **10 Hz** using PCHIP interpolation to prevent oscillations between sparse data points  


Velocity was converted into the North-East-Down (NED) frame which is conventionally the Nav Frame:

$$
V_h = \sqrt{V^2 - V_{\text{vertical}}^2}
$$

$$
V_N = V_h \cos(\psi)
$$

$$
V_E = V_h \sin(\psi)
$$

$$
V_D = -V_{vertical}
$$

Acceleration was computed via numerical differentiation:

$$
a_N = \frac{dV_N}{dt}, \quad
a_E = \frac{dV_E}{dt}, \quad
a_D = \frac{dV_D}{dt}
$$

These represent true acceleration in the NED/Nav frame.

To get the Euler angles from the flight data, some assumptions for large commercial aircraft were made because no orientation data was given:
- Yaw is assumed to be heading:

$$
\psi = \text{heading}
$$

- Pitch is assumed to be equal to the flight path angle:

$$
\theta = \gamma = \arctan\left(\frac{V_{vertical}}{V_h}\right)
$$

- Roll assumed coordinated no slip turns:

$$
\phi = \arctan\left(\frac{V_{heading}\dot{\psi}}{g}\right)
$$


Time rate derivatives computed via gradients:

$$
\dot{\phi}, \quad \dot{\theta}, \quad \dot{\psi}
$$

Body angular rates are computed as:

$$
\begin{pmatrix}
p \\
q \\
r
\end{pmatrix}=
\begin{pmatrix}
1 & 0 & -\sin\theta \\
0 & \cos\phi & \sin\phi \cos\theta \\
0 & -\sin\phi & \cos\phi \cos\theta
\end{pmatrix}
\begin{pmatrix}
\dot{\phi} \\
\dot{\theta} \\
\dot{\psi}
\end{pmatrix}
$$


Specific force is computed in the NED frame:

$$
\mathbf{f}_{nav} = \mathbf{a}_{nav} - \mathbf{g}_{nav}
$$

Then rotated into the body frame using a 3-2-1 Direction Cosine Matrix:

$$
\mathbf{f}_{body} = C_{b/n} \, \mathbf{f}_{nav}
$$


Gravity is computed using the WGS-84 model:
- Function: `gravitywgs84`
- Evaluated at each interpolated position

Truth IMU data is validated by integrating acceleration:

$$
\mathbf{v}(t) = \int \mathbf{a}(t)\, dt
$$

The reconstructed velocity is compared against the original trajectory to ensure consistency.

#### IMU Error Modeling
After generating truth data, IMU errors are injected based on the sensor specification sheet. Modeled error sources include: 
- Bias  
- Cross-axis sensitivity  
- White noise  

Models not included were: 
- Temperature-dependent effects (no data available)

The final IMU dataset includes:
- Body-frame specific force measurements  
- Body-frame angular rates (p, q, r)  
- Injected sensor errors  

Notes on the IMU Dataset: 
- All computations are performed at **10 Hz**
- NED is used as the navigation frame
- Body frame is derived using **3-2-1 (ZYX) rotations**

### Generation of VOR Measurements
VOR location data is loaded in from the `datasets/` routine and then processed at the top of the filter script. The `vorMeas` function then accepts aircraft positioning in either **ECEF** (Earth-Centered, Earth-Fixed) or **LLA** (Latitude, Longitude, Altitude) formats. 

To determine the spatial relationship between the aircraft and each navaid, the function utilizes two primary mathematical models:
* **Haversine Formula:** Used to calculate the great-circle distance between the aircraft and the VOR, accounting for the Earth's curvature.
* **Bearing Calculation:** Computes the forward azimuth from the VOR station to the aircraft. This represents the **true radial** the aircraft is currently occupying.

The function then applies standard FAA Service Volume constraints to determine signal "reception." A VOR is only returned if the aircraft falls within its defined range-altitude envelope:

| VOR Class | Altitude Range (ft) | Max Range (NM) |
| :--- | :--- | :--- |
| **Low (L-VOR)** | 1,000 to 14,500 | 40 |
| **High (H-VOR)** | 1,000 to 14,500 | 40 |
| | 14,500 to 18,000 | 100 |
| | 18,000 to 45,000 | 130 |
| | 45,000 to 60,000 | 100 |

The function returns an array of structs (`validVORs`), each containing the original navaid metadata appended with calculated measurement fields:
* **lat / lon**: Explicit coordinate mapping for the station.
* **distance_m**: The calculated great-circle distance in meters.
* **bearing_deg**: The true bearing from the station to the aircraft (0–359°).

### VOR Measurement Model
<img width="960" height="720" alt="VOR Measurement Model (1)" src="https://github.com/user-attachments/assets/5a6ed984-e55c-466d-9215-f78262b91298" />
> Link to draw file: https://docs.google.com/drawings/d/16Xxh9EqZsI2CmBNIY7Cf3i_a1hbbCS4FhxrkBuhYRwA/edit?usp=sharing


## The Extended Kalman Filter
An Extended Kalman Filter (EKF) was selected to fuse the generated IMU and VOR measurements into a final state estimate. An EKF was chosen as the relationship between the measurements 




## References

[1]https://spacenews.com/the-race-to-back-up-vulnerable-gps/ 

[2]https://www.forbes.com/sites/erictegler/2023/09/28/someone-in-the-middle-east-is-leading-aircraft-astray-by-spoofing-gps-signals/ 

[3]https://en.wikipedia.org/wiki/VOR/DME#:~:text=In%20radio%20navigation%2C%20a%20VOR,the%20receiver%20and%20the%20station 

[4]https://adds-faa.opendata.arcgis.com/datasets/990e238991b44dd08af27d7b43e70b92_0/explore?location=5.144894%2C-1.658576%2C2.00 

[5]https://epicflightacademy.com/what-is-vor/ 

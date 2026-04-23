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


## The Extended Kalman Filter
An Extended Kalman Filter (EKF) was selected to fuse the generated IMU and VOR measurements into a final state estimate. An EKF was chosen as the relationship between the VOR measurements and the IMU data propogation is nonlinear. 

### Time Update

The IMU measurements were propogated in the ECI frame as whole-state position, and velocity, with Euler angle orientation. The propagation equation was defined as:

$$
\begin{pmatrix}
\dot{x}_{ECI} \\
\dot{v}_{ECI} \\
\dot{\phi}
\end{pmatrix} = \begin{pmatrix}
v_{ECI} \\
T_{NED}^{ECI} * sf_{NED} + g_{ECI} \\
f(angularRates, \phi)
\end{pmatrix}
$$

With propagation of aircraft Euler angles as:

$$
f([p, q, r]_{body}, [\phi, \theta, \psi]) = \begin{pmatrix}
p + q * sin(\phi) * tan(\theta) + r * cos(\phi) * tan(\theta)\\
q * cos(\phi) - r * sin(\phi)\\
q * sin(\phi) + r * cos(\phi)
\end{pmatrix}
$$

The F matrix jacobian is calculated numerically, and used to propagate covariance in line with the Kalman filtering equations

$$
P_{n+1,n} = FP_{n,n}F^T + Q
$$

The process noise matrix is defined by the velocity random walk (VRW) and angular random walk (ARW) values defined in the IMU spec sheet, specifically:

$$
Q = \begin{pmatrix}
0_{3x3} && ... && ...\\
... && VRW{3x3} && ...\\
... && ... && ARW{3x3}
\end{pmatrix}
$$

### VOR Measurement Model
<img width="1216" height="843" alt="image" src="https://github.com/user-attachments/assets/ca992581-211e-4922-8a28-a065fa889eb5" />


### Measurement Update

The H matrix jacobian is obtained through the non-linear measurement model as described above, and fed into the **Jacob Form** of the Kalman update equations. This was done as instabilities were found to appear rapidly after a new VOR entered the plane's line of sight, and will be discussed further in the results section.

$$
H(:,1:3) = [-\frac{\Delta r_{N_n}}{\rho^2}, \frac{\Delta r_{E_n}}{\rho^2}, 0] * T_{ECI}^{NED_{VOR_n}}
$$

$$
H(:,4:9) = 0_{7xn}
$$

$$
K = \frac{P_{n,n-1}H^T}{HP_{n,n-1}H^T + R}
$$

With R defined as the VOR bearing measurement 1 sigma of 2 degrees.

$$
x_{n,n} = x_{n,n-1} + K(z - Hx_{n,n-1})
$$

With $z_n$ being the true measurement value

$$
P_{n,n} = (I_{9x9} - KH) * P_{n,n-1} * (I_{9x9} - KH)^T + K * R * K
$$

### Altitude Measurement Update

The VOR measurements only update the latitude and longitude of the aircraft, so a second measurement update was required for altitude. This was assuming a constant altimeter update within the aircraft. The innovation was the difference between the true altitude and the estimated altitude. The Jacobian of the H matrix is calculated numerically with perturbations, and an altitude measurement noise of 50 ft is applied in the Kalman filter update equations.

## Results
The trajectory used for this demonstration was a commercial flight from San Diego to Sacramento. Throughout the trajectory, the Aircraft came into range of 27 unique VORs, and had constant coverage of at least 2 VORs.

<img width="777" height="799" alt="image" src="https://github.com/user-attachments/assets/d6f8efd3-8671-4117-9b92-8dbc6401b8f7" />

The plot of the flight path with VOR stations marked as blue triangles.


<img width="806" height="585" alt="Initial Uncertainty Shrink-1(1)(1)" src="https://github.com/user-attachments/assets/84525fd8-79cd-4a85-a605-27ff7a1fd46d" />

This simulation shows the initial uncertainty reducing as the VOR measurements update the aircraft's position.


<img width="800" height="700" alt="sf_approach(1)" src="https://github.com/user-attachments/assets/4331a64d-6d58-42db-a1e1-f8e015ef3ce7" />

Ground track with uncertainty bounds throughout the end of the flight


## Challenges 

### Dynamic propagation
In our final results, the propagation model was shown to be highly inconsistent with expected reality. To combat this, the process noise was raised to increase the filter's reliability on VOR measurements.

<img width="1462" height="949" alt="ECI_estimates" src="https://github.com/user-attachments/assets/7c37fb5e-c8ac-46aa-9e5f-c0bb2933f9b5" />

The ECI position uncertainty approached a "steady state" of ~2 km, with the state estimates showing observability and general good behavior within the EKF. 


### Numerical Stability
As these trajectories were tens of thousands of indices long (multiple hours of 10hz data), the covariance matrix was found to lose symmetry due to numerical stability and diverge. This was solved with the following equation, which would enforce P symmetry at each measurment update.

```
  results.P(:,:,i+1) = IKH_alt * results.P(:,:,i+1) * IKH_alt' + K_alt * R_alt * K_alt';
  results.P(:,:,i+1) = 0.5 * (results.P(:,:,i+1) + results.P(:,:,i+1)');
```


## References

[1]https://spacenews.com/the-race-to-back-up-vulnerable-gps/ 

[2]https://www.forbes.com/sites/erictegler/2023/09/28/someone-in-the-middle-east-is-leading-aircraft-astray-by-spoofing-gps-signals/ 

[3]https://en.wikipedia.org/wiki/VOR/DME#:~:text=In%20radio%20navigation%2C%20a%20VOR,the%20receiver%20and%20the%20station 

[4]https://adds-faa.opendata.arcgis.com/datasets/990e238991b44dd08af27d7b43e70b92_0/explore?location=5.144894%2C-1.658576%2C2.00 

[5]https://epicflightacademy.com/what-is-vor/ 

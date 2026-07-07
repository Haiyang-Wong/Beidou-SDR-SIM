# Beidou-SDR-SIM
An open-source tool for generating  Beidou Navigation Satellite **B1C** Signal. It supports custom **positioning** (**static locations or dynamic trajectories**) and **timing**, and **real-time USRP transmission**.

## Requirements

1. g++
2. Cmake
3. UHD (usrp api and dev library)

#### Installation

```
git clone https://github.com/Haiyang-Wong/Beidou-SDR-SIM.git
cd Beidou-SDR-SIM
mkdir build && cd build
cmake ../
make
```

#### Simulation parameters (configured in main.cpp)

| Parameters | Description                                    | Default value                      |
| :--------- | :--------------------------------------------- | ---------------------------------- |
| eph        | Ephemeris file                                 | BRDC00IGS_R_20211700000_01D_MN.rnx |
| `dur`      | Simulation duration (seconds)                  | 120                                |
| `lla`      | Location (latitude °, longitude °, altitude m) | (50,100,450)                       |
| `rate`     | Sampling rate (Hz)                             | 10.23e6                            |
| `elvmask`  | Elevation mask (degrees)                       | 10                                 |
| `traj`     | Dynamic trajectory(time,lat,lon,alt)           | demo.csv                           |
| `usrp`     | Real-time transmission using USRP              | -                                  |

#### Execution

##### Static position

```
./build/BDS_SDR_SIM \
  --eph rinex_files/BRDC00IGS_R_20211700000_01D_MN.rnx \
  --dur 120 \
  --lla 27.9881, 86.925, 8848.86 \
  --rate 10230000 \
  --out output/bds_b1c.ishort
```

This project will generate baseband I/Q sampling data and output it to `output/` directory. 

##### Dynamic trajectory

The trajectory time is seconds from the scenario start. Positions are linearly interpolated at the simulator 0.1 s update interval.

CSV format:

```csv
time,lat,lon,alt
0.0,34.000000,108.000000,450
1.0,34.000010,108.000020,450
2.0,34.000025,108.000050,452
```

Example:

```bash
./build/BDS_SDR_SIM --dur 240 --traj traj/demo.csv --out bds_dynamic.ishort
```

##### Real-time transmission

The project automatically detects a UHD device and transmits using the default BDS B1C frequency, sample rate, and gain. 

Example:

```
./build/BDS_SDR_SIM --usrp
```

If you need to set parameters manually, you can append:

```
./build/BDS_SDR_SIM --usrp \
  --usrp-args "serial=xxxxxx" \
  --gain 30 \
  --tx-freq 1575420000 \
  --rate 10230000 \
  --tx-prebuffer 20
```

#### Evaluation

The signal can be transmitted using the SDR (USRP B210) and verified using Beitian receivers, as well as using the software receiver [CU-SDR-Collection](https://github.com/gnsscusdr/CU-SDR-Collection). Some results are shown below.

![](Beitian-BE609U-Output.png)

![](CU-SDR-Collection-Output.png)

### Future work

- Integrate **GPS L1** and **Galileo E1B** signal

## Acknowledgements

- [GPS-SDR-SIM](https://github.com/osqzss/gps-sdr-sim) and [GALILEO-SDR-SIM](https://github.com/harshadms/galileo-sdr-sim) has been an inspiration to start this project.
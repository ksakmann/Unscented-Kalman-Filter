# Unscented Kalman Filter Project 
Self-Driving Car Engineer Nanodegree Program

---
## Introduction
The unscented Kalman filter is a way to improve on the extended Kalman Filter. Unlike the EKF the UKF does not linearize the 
state equations. It relies on constructing sigma points that get propagated through the state vector model. 

Shown below are the results of this project for two datasets.

[//]: # (Image References)
[image1]: ./images/position.png
[image2]: ./images/NIS.png

![UKF prediction][image1]

The noise parameters were chosen in such a way to make the normalized innovation squared close to its statistically expected value.
The radar measurement space is three dimensional (rho, phi, rho_dot) and the chi-squared value for a 95% confidence intervall is 7.8.
The lidar measurement space is two dimensional (x,y) and the chi-squared value for a 95% confidence intervall is 6. Averaging these 
two one would expect about 5% of all predicted states to have a chi-squared value of 7 or higher. This is approximately true 
for the chosen noise parameters. 


![NIS][image2]


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Editor Settings

Use the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Generating Additional Data

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/4d0420af-0527-4c9f-a5cd-56ee0fe4f09e)
for instructions and the project rubric.

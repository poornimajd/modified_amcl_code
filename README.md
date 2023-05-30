# modified_amcl_code

Our paper on integrating external cameras with on-board lidar for robust robot localization has been accepted to the Advances in Robotics, 2023 , 6TH INTERNATIONAL CONFERENCE OF THE ROBOTICS SOCIETY.  

This is the modified version of the AMCL algorithm from the ROS2 Galactic version. The correction step has been modified that is in addition to the weights from the laser, the particles are also
assigned weights from external cameras. The camera observation model is a multivariate Gaussian with a mean and covariance matrix. 

   The code for the camera observation model will be released after the paper is published.

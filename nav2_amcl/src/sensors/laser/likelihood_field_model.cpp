/*
 *  Player - One Hell of a Robot Server
 *  Copyright (C) 2000  Brian Gerkey   &  Kasper Stoy
 *                      gerkey@usc.edu    kaspers@robotics.usc.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <math.h>
#include <assert.h>
#include <memory>
#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "nav2_amcl/sensors/laser/laser.hpp"
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"
#include "rclcpp/rclcpp.hpp"
#include <vector>
#include <cmath>
#include<bits/stdc++.h>
#include <cmath>
#include "iostream"
// #include "code.cpp"
#include <fstream>
#include <fstream>
#include <sstream>
#include "iostream"


std::ofstream fw("./amcl_poses_hardware.txt", std::ofstream::out);
namespace nav2_amcl
{
geometry_msgs::msg::PoseStamped pose_image_wt;
geometry_msgs::msg::PoseStamped pose_image_wt2;

LikelihoodFieldModel::LikelihoodFieldModel(
  double z_hit, double z_rand, double sigma_hit,
  double max_occ_dist, size_t max_beams, map_t * map)
: Laser(max_beams, map)
{
  z_hit_ = z_hit;
  z_rand_ = z_rand;
  sigma_hit_ = sigma_hit;
  map_update_cspace(map, max_occ_dist);
}

double
LikelihoodFieldModel::sensorFunction(LaserData * data, pf_sample_set_t * set)
{
  

  
  double mean_pos_x=pose_image_wt.pose.position.x;//mean x
  double mean_pos_y=pose_image_wt.pose.position.y;//mean y

  double mean_pos_theta=pose_image_wt.pose.position.z;//mean theta
  
  double cov_x=(pose_image_wt.pose.orientation.x)+0.00000001; // add a small vlaue, to avoid cov_x=0
  double cov_y=(pose_image_wt.pose.orientation.y)+0.00000001;
  double cov_z=(pose_image_wt.pose.orientation.z)+0.0000001;
  double qk[3][3]={{cov_x,0,0},{0,cov_y,0},{0,0,cov_z}};

  double mean_pos_x2=pose_image_wt2.pose.position.x;
  double mean_pos_y2=pose_image_wt2.pose.position.y;

  double mean_pos_theta2=pose_image_wt2.pose.position.z;//+0.01;//mean theta

  
  double cov_x2=(pose_image_wt2.pose.orientation.x)+0.00000001;
  double cov_y2=(pose_image_wt2.pose.orientation.y)+0.00000001;
  double cov_z2=(pose_image_wt2.pose.orientation.z)+0.0000001;
  
  double qk2[3][3]={{cov_x2,0,0},{0,cov_y2,0},{0,0,cov_z2}};

  LikelihoodFieldModel * self;
  int i, j, step;
  double z, pz;
  double p;
  double obs_range, obs_bearing;
  double total_weight;
  pf_sample_t * sample;
  pf_vector_t pose;
  pf_vector_t hit;

  self = reinterpret_cast<LikelihoodFieldModel *>(data->laser);

  total_weight = 0.0;


  double multi=pow(determinantOfMatrix(qk, 3),0.5);

  double term1=1/(pow((2*3.14),3/2.0)*multi);
  double inv[3][3];
  term1=term1+0.0;
  double wt_cam;
  std::vector< double > laserwts;
  std::vector< double > camwts;

  inverse(qk,inv);

  double multi2=pow(determinantOfMatrix(qk2, 3),0.5);

  double term12=1/(pow((2*3.14),3/2.0)*multi2);
  double inv2[3][3];
  term12=term12+0.0;
  double wt_cam2;

  std::vector< double > camwts2;

  inverse(qk2,inv2);

// to store in a file
    fw << "means"<<" "<<mean_pos_x<<" "<<mean_pos_y<<" "<<mean_pos_theta<<" "<<cov_x<<" "<<cov_y<<" "<<cov_z << "\n";
    fw << "means2"<<" "<<mean_pos_x2<<" "<<mean_pos_y2<<" "<<mean_pos_theta2<<" "<<cov_x2<<" "<<cov_y2<<" "<<cov_z2 << "\n";
   
    for (j = 0; j < set->sample_count; j++) {
    sample = set->samples + j;
    pose = sample->pose;
    wt_cam=0.0;
    wt_cam2=0.0;


    double new_pose=fmod(pose.v[2],6.28);
    if (new_pose>=3.14)
    {
      new_pose=new_pose-6.28;
    
  }
  if (new_pose<-3.14)
    {new_pose+=6.28;}

    if (mean_pos_x!=-1000.0)
    {wt_cam=gaussian_pdf(pose.v,term1,inv,mean_pos_x,mean_pos_y,mean_pos_theta,new_pose);}

    if (mean_pos_x2!=-1000)
    {
      wt_cam2=gaussian_pdf(pose.v,term12,inv2,mean_pos_x2,mean_pos_y2,mean_pos_theta2,new_pose);
    }
    
    
    pose = pf_vector_coord_add(self->laser_pose_, pose);


    p = 1.0;

    // Pre-compute a couple of things
    double z_hit_denom = 2 * self->sigma_hit_ * self->sigma_hit_;
    double z_rand_mult = 1.0 / data->range_max;

    step = (data->range_count - 1) / (self->max_beams_ - 1);

    // Step size must be at least 1
    if (step < 1) {
      step = 1;
    }
    for (i = 0; i < data->range_count; i += step) {
      obs_range = data->ranges[i][0];
      obs_bearing = data->ranges[i][1];

      // This model ignores max range readings
      if (obs_range >= data->range_max) {
        continue;
      }

      // Check for NaN
      if (obs_range != obs_range) {
        continue;
      }

      pz = 0.0;

      // Compute the endpoint of the beam
      hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
      hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

      // Convert to map grid coords.
      int mi, mj;
      mi = MAP_GXWX(self->map_, hit.v[0]);
      mj = MAP_GYWY(self->map_, hit.v[1]);

      // Part 1: Get distance from the hit to closest obstacle.
      // Off-map penalized as max distance
      if (!MAP_VALID(self->map_, mi, mj)) {
        z = self->map_->max_occ_dist;
      } else {
        z = self->map_->cells[MAP_INDEX(self->map_, mi, mj)].occ_dist;
      }
      // Gaussian model
      // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
      pz += self->z_hit_ * exp(-(z * z) / z_hit_denom);
      // Part 2: random measurements
      pz += self->z_rand_ * z_rand_mult;

      // TODO(?): outlier rejection for short readings

      assert(pz <= 1.0);
      assert(pz >= 0.0);

      p += pz * pz * pz;
      // }
    }

    laserwts.push_back(p);
    camwts.push_back(wt_cam);
    camwts2.push_back(wt_cam2);


  }

  long double maxlaser=(accumulate(laserwts.begin(),laserwts.end(), 0.0));
  long double maxcam=(accumulate(camwts.begin(),camwts.end(), 0.0));
  long double maxcam2=(accumulate(camwts2.begin(),camwts2.end(), 0.0));

  if (maxcam==0)
    maxcam=1;
  if (maxlaser==0)
    maxlaser=1;
  if (maxcam2==0)
    maxcam2=1;

  for (j = 0; j < set->sample_count; j++) {
    sample = set->samples + j;
    pose = sample->pose;

    laserwts[j]=laserwts[j]/maxlaser;

    if (maxcam==1)
    {
     
      if (camwts[j]==0)

        {
          camwts[j]=1;
        }

    }
    camwts[j]=camwts[j]/maxcam;
      
   

    if (maxcam2==1)
    {
     
      if (camwts2[j]==0)

        {
          camwts2[j]=1;
        }

    }

    camwts2[j]=camwts2[j]/maxcam2;
    

    
    sample->weight *=(laserwts[j]*camwts[j]*camwts2[j]);

    total_weight += sample->weight;
  }

  return total_weight;
  }

  

bool
LikelihoodFieldModel::sensorUpdate(pf_t * pf, LaserData * data, geometry_msgs::msg::PoseStamped & image_pose_rx_,geometry_msgs::msg::PoseStamped & image_pose_rx_2)
{


  pose_image_wt=image_pose_rx_;
  pose_image_wt2=image_pose_rx_2;


  if (max_beams_ < 2) {
    return false;
  }
  pf_update_sensor(pf, (pf_sensor_model_fn_t) sensorFunction, data);

  return true;
}

 // namespace nav2_amcl


void 
LikelihoodFieldModel::getCofactor(double mat[3][3], double temp[3][3], int p,
                 int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q)
            {
                temp[i][j++] = mat[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
 
/* Recursive function for finding determinant of matrix.
   n is current dimension of mat[][]. */
double ///corr to double  
LikelihoodFieldModel::determinantOfMatrix(double mat[3][3], int n)
{
    double D = 0.0; // Initialize result //corr double
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return mat[0][0];
 
    double temp[3][3]; // To store cofactors
 
    double sign = 1.0; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of mat[0][f]
        getCofactor(mat, temp, 0, f, n);
        D += sign * mat[0][f]
             * determinantOfMatrix(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}

void 
LikelihoodFieldModel::adjoint(double A[3][3],double adj[3][3])
{
    if (3 == 1)
    {
        adj[0][0] = 1;
        return;
    }
 
    // temp is used to store cofactors of A[][]
    double sign = 1.0;
    double temp[3][3];
 
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, 3);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1.0: -1.0;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinantOfMatrix(temp, 3-1));
        }
    }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
bool 
LikelihoodFieldModel::inverse(double A[3][3], double inverse[3][3])
{
    // Find determinant of A[][]
    double det = determinantOfMatrix(A, 3);//corr
    if (det == 0)
    {
        //cout << "Singular matrix, can't find its inverse";
        return false;
    }
 
    // Find adjoint
    double adj[3][3];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)

    inverse[0][0]=1/A[0][0];
    inverse[0][1]=0;
    inverse[0][2]=0;
    inverse[1][0]=0;
    inverse[1][1]=1/A[1][1];
    inverse[1][2]=0;
    inverse[2][0]=0;
    inverse[2][1]=0;
    inverse[2][2]=1/A[2][2];

 
    return true;
}

double 
LikelihoodFieldModel::gaussian_pdf(double posev[3],double term1,double term2[3][3],double mean_pos_x,double mean_pos_y,double mean_pos_theta,double new_pose)
  {
    double mat1[3][1];

    mat1[0][0]=posev[0]-mean_pos_x;
    mat1[1][0]=posev[1]-mean_pos_y;
    mat1[2][0]=new_pose-mean_pos_theta;

    double tmat1[1][3];
    tmat1[0][0]=posev[0]-mean_pos_x;
    tmat1[0][1]=posev[1]-mean_pos_y;
    tmat1[0][2]=new_pose-mean_pos_theta;
    double prod1[1][3];
    prod1[0][0]=0.0;
    prod1[0][1]=0.0;
    prod1[0][2]=0.0;

    for (int pp=0;pp<3;pp++)
    {
      for (int nn=0;nn<3;nn++)
      {
        prod1[0][pp]+=(tmat1[0][nn]*term2[nn][pp]);
      }
    }

    double prod2=0.0;
    for (int pp=0;pp<3;pp++)
    {
      prod2+=prod1[0][pp]*mat1[pp][0];
    }

    double term3=std::exp(-0.5*prod2);

    double cost=term1*term3;

    return cost;///pow(10,16);
  }

}//namespcae nav2
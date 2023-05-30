#include <./Eigen/Dense>
#include <./Eigen/Cholesky>
#include <./Eigen/SVD>
#include <./Eigen/Eigenvalues>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    
#include <fstream>
#include <cmath>
/*
  We need a functor that can pretend it's const,
  but to be a good random number generator 
  it needs mutable state.
*/
using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor | Eigen::AutoAlign>;
namespace Eigen {
namespace internal {
template<typename Scalar> 
struct scalar_normal_dist_op 
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

/*
  Draw nn samples from a size-dimensional normal distribution
  with a specified mean and covariance
*/
void multi_pdf(double mean_pos_x,double mean_pos_y,double mean_pos_theta,double qk[3][3],double* sample_x,double* sample_y,double* sample_theta)
{
  int size = 3; // Dimensionality (rows)
  //int nn=1;     // How many samples (columns) to draw
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  //Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng

  // Define mean and covariance of the distribution
  Eigen::VectorXd mean(size);       
  Eigen::MatrixXd covar(size,size);

  mean  <<  mean_pos_x,mean_pos_y,mean_pos_theta;
  //std::cout << qk[1];

  covar << qk[0][0],qk[0][1],qk[0][2],
           qk[1][0],qk[1][1],qk[1][2],
           qk[2][0],qk[2][1],qk[2][2];
           
  //Eigen::MatrixXd normTransform(size,size);

  //Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);
   Eigen::LLT<Matrix2D> cholSolver(covar);
  // We can only use the cholesky decomposition if 
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  Eigen::MatrixXd normTransform;
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    
    //std::cout << "Broken\n";
        Eigen::BDCSVD<Matrix2D> solver(covar, Eigen::ComputeFullV);
        normTransform = (-solver.matrixV().transpose()).array().colwise() * solver.singularValues().array().sqrt();
        // Use eigen solver
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    //normTransform = eigenSolver.eigenvectors() 
     //              * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }
 
  Matrix2D gaussianSamples = Eigen::MatrixXd::NullaryExpr(1, mean.size(), randN);

  Eigen::MatrixXd samples = (gaussianSamples * normTransform).rowwise() + mean.transpose();
  //std::cout << "samples\n" << samples << std::endl;
  //Eigen::MatrixXd samples = (normTransform 
  //                         * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise() 
  //                         + mean;
   //outdata << samples;
  //std::cout << "Mean\n" << mean << std::endl;
  //std::cout << "Covar\n" << covar << std::endl;
  //std::cout << "Samples\n" << samples(0,0)<<" "<< samples(1,0)<<" "<< samples(2,0) << std::endl;
  /*if (isnan(samples(0,0)))
  {
  std::cout << "Covar\n" << covar << std::endl;
  }
  else
  {
  std::cout << "not Covar\n" << covar << std::endl;
  }*/
  /*while (isnan(samples(0,0)))
  {
  Eigen::MatrixXd samples = (normTransform 
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise() 
                           + mean;
  std::cout << "samples\n" << samples << std::endl;
  }*/
  
  *sample_x=samples(0,0);
  *sample_y=samples(0,1);
  *sample_theta=samples(0,2);
  
  //return samples(0,0),samples(1,0),samples(2,0);
}

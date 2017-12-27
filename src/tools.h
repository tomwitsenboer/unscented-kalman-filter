#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "Eigen/SVD"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Ref;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE (root-mean-square error).
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
   * A helper method to normalize angle, i.e. make it belong to [-pi, pi].
   */
  void NormalizeAngle(double *angle);

};

#endif /* TOOLS_H_ */

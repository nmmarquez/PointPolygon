#define TMB_LIB_INIT R_init_PointPolygon
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    using namespace R_inla;
    using namespace density;
    using namespace Eigen;
    
    // Counts of observed values
    DATA_IVECTOR(yPoint);
    DATA_IVECTOR(yPoly);
    DATA_IVECTOR(idPoint);

    // Denoms
    DATA_IVECTOR(denomPoint);
    DATA_IVECTOR(denomPoly);
    DATA_IVECTOR(idPoly);

    //Covariates
    DATA_MATRIX(covs);

    // Projections
    DATA_SPARSE_MATRIX(AprojObs);
    DATA_SPARSE_MATRIX(AprojPoly);

    // SPDE objects
    DATA_SPARSE_MATRIX(M0);
    DATA_SPARSE_MATRIX(M1);
    DATA_SPARSE_MATRIX(M2);
    
    // Polygon Modeling Option
    DATA_INTEGER(moption);
    // Should I Evaluate Priors? Critical for TMBstan
    DATA_INTEGER(priors);

    // Parameters
    PARAMETER_VECTOR(beta);
    PARAMETER(log_tau);
    PARAMETER(log_kappa);
    PARAMETER_VECTOR(z);

    int Npoint = yPoint.size();
    int Npoly = yPoly.size();

    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);

    Type nll = 0.0;

    if(priors == 1){
        for(int b=0; b<beta.size(); b++){
            nll -= dnorm(beta[b], Type(0.0), Type(100.0), true);
        }
        nll -= dnorm(log_tau, Type(0.0), Type(100.0), true);
        nll -= dnorm(log_kappa, Type(0.0), Type(100.0), true);
    }

    SparseMatrix<Type> Q = spde_Q(log_kappa, log_tau, M0, M1, M2);

    nll += GMRF(Q)(z);

    vector<Type> projLatF = AprojObs * z;
    vector<Type> projCov = covs * beta;
    vector<Type> projLatObs = projLatF + projCov;
    vector<Type> projPObs = exp(projLatObs) / (Type(1.) + exp(projLatObs));
    
    for(int i=0; i<Npoint; i++){
        Type p = projPObs[idPoint[i]];
        nll -= dbinom(Type(yPoint[i]), Type(denomPoint[i]), p, true);
    }
    
    Type nllPart;
    Type p;
    
    
    if(Npoly != 0){
        // mixture model estimation
        if(moption == 0){
            for(int i=0; i<Npoly; i++){
                nllPart = Type(.0);
                for(typename SparseMatrix<Type>::InnerIterator it(AprojPoly,idPoly[i]); it; ++it){
                    // see the following link for indexing sparse matrices
                    // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title2
                    p = projPObs[it.row()];
                    nllPart += dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, false) * it.value();
                }
                nll -= log(nllPart);
            }
        }
        // Redistribute points
        if(moption == 1){
            for(int i=0; i<Npoly; i++){
                Type p = projPObs[idPoly[i]];
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
        if(moption == 3){
            // Reimann sum approximation
            SparseMatrix<Type> RAprojPoly = AprojPoly.transpose();
            vector<Type> projPoly = RAprojPoly * projPObs;
            for(int i=0; i<Npoly; i++){
                Type p = projPoly[i];
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
    }

    REPORT(z);
    return nll;
}


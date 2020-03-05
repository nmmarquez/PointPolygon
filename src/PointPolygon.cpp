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
    DATA_IVECTOR(idtPoint);

    // Denoms
    DATA_IVECTOR(denomPoint);
    DATA_IVECTOR(denomPoly);
    DATA_IVECTOR(idPoly);
    DATA_IVECTOR(idtPoly);

    //Covariates
    DATA_ARRAY(covs);

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
    PARAMETER(logit_rho);
    PARAMETER_ARRAY(z);

    int Npoint = yPoint.size();
    int Npoly = yPoly.size();

    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    Type rho = exp(logit_rho) / (Type(1.) + exp(logit_rho));

    Type nll = 0.0;

    if(priors == 1){
        for(int b=0; b<beta.size(); b++){
            nll -= dnorm(beta[b], Type(0.0), Type(100.0), true);
        }
        nll -= dnorm(log_tau, Type(0.0), Type(100.0), true);
        nll -= dnorm(log_kappa, Type(0.0), Type(100.0), true);
    }

    SparseMatrix<Type> Q = spde_Q(log_kappa, log_tau, M0, M1, M2);

    if(z.dim(1) > 1){
        nll += SEPARABLE(AR1(rho), GMRF(Q))(z);
    }
    else{
        nll += GMRF(Q)(z.matrix().col(0));
    }

    // Turn z into a mtrix here so projLatF is space by time
    matrix<Type> projLatF = AprojObs * z.matrix();
    // covs should be an array that we loop through time to have space by time
    matrix<Type> projCov(projLatF.rows(), projLatF.cols());
    for(int s=0; s<projLatF.rows(); s++){
        for(int t=0; t<projLatF.cols(); t++){
            projCov(s,t) = Type(0.0);
            for(int b=0; b<beta.size(); b++){
                projCov(s,t) += covs(s,b,t) * beta(b);
            }
        }
    }

    matrix<Type> projLatObs = projLatF + projCov;
    matrix<Type> projPObs = projLatObs.array().exp() / (Type(1.) + projLatObs.array().exp());

    for(int i=0; i<Npoint; i++){
        Type p = projPObs(idPoint[i],idtPoint[i]);
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
                    p = projPObs(it.row(),idtPoly[i]);
                    nllPart += dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, false) * it.value();
                }
                nll -= log(nllPart);
            }
        }
        // Redistribute points
        if(moption == 1){
            for(int i=0; i<Npoly; i++){
                Type p = projPObs(idPoly[i],idtPoly[i]);
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
        if(moption == 3){
            // Reimann sum approximation
            SparseMatrix<Type> RAprojPoly = AprojPoly.transpose();
            matrix<Type> projPoly = RAprojPoly * projPObs;
            for(int i=0; i<Npoly; i++){
                p = projPoly(idPoly[i],idtPoly[i]);
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
    }

    REPORT(z);
    REPORT(projLatObs);
    REPORT(projPObs);
    return nll;
}

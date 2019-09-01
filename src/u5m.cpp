
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;
template<class Type>
SparseMatrix<Type> spde2_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

template<class Type>
matrix<Type> rw2_Q(int N, Type sigma) {
    matrix<Type> Q(N,N);
    Q(0,0) = (1.) / pow(sigma, 2.);
    Q(0,1) = (-2.) / pow(sigma, 2.);
    Q(1,0) = (-2.) / pow(sigma, 2.);
    Q(1,1) = (5.) / pow(sigma, 2.);
    for (int n = 2; n < (N-1); n++) {
        Q(n,n) = (6.) / pow(sigma, 2.);
        Q(n-1,n) = (-4.) / pow(sigma, 2.);
        Q(n,n-1) = (-4.) / pow(sigma, 2.);
        Q(n-2,n) = (1.) / pow(sigma, 2.);
        Q(n,n-2) = (1.) / pow(sigma, 2.);
    }
    Q(N-2,N-2) = (5.) / pow(sigma, 2.);
    Q(N-1,N-2) = (-2.) / pow(sigma, 2.);
    Q(N-2,N-1) = (-2.) / pow(sigma, 2.);
    Q(N-1,N-1) = (1.) / pow(sigma, 2.);
    return Q;
}

template<class Type>
Type densRW2(vector<Type>x, Type sigma) {
    Type kappa_ = pow(sigma, -2.);
    int N = x.size();
    matrix<Type> xt(1, N);
    for(int n = 1; n < (N); n++){
        xt(0, n) = x[n];
    }
    matrix<Type> Q = rw2_Q(N, sigma);
    matrix<Type> temp = xt * Q;
    temp = temp * x;
    Type nll = -1. * (((Type(N) - 2.) / 2.) * log(kappa_) + (-.5 * temp(0,0)));
    return nll;
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

    // Denoms
    DATA_IVECTOR(denomPoint);
    DATA_IVECTOR(denomPoly);

    // Identifiers
    DATA_IVECTOR(idPoint); // geolocated data pixel index
    DATA_IVECTOR(idtPoint); // geolocated data time index
    DATA_IVECTOR(idaPoint); // geolocated data age index
    DATA_IVECTOR(idnPoint); // geolocated data survey index
    DATA_IVECTOR(idcPoint); // geolocated data cluster index
    DATA_IVECTOR(idPoly); // geomasked data loc index
    DATA_IVECTOR(idtPoly); // geomasked data time index
    DATA_IVECTOR(idaPoly); // geomasked data age index
    DATA_IVECTOR(idnPoly); // geomasked data survey index
    DATA_IVECTOR(idcPoly); // geomasked data cluster index

    // Covariates
    DATA_ARRAY(covs); // covariate value for entire prediction surface and time

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
    PARAMETER_VECTOR(beta); // linear beta coefficients for covariates
    PARAMETER_VECTOR(beta_age); // age effects
    PARAMETER(log_tau); // 2d spatial field parameter for z
    PARAMETER(log_kappa); // 2d spatial field parameter for z
    PARAMETER(logit_rho); // temporal autocorrelation for z
    PARAMETER_VECTOR(log_sigma_phi); // variance over time within age groups
    PARAMETER(log_sigma_nu); // sd on survey random effects
    PARAMETER(log_sigma_epsilon); // sd on time unstructured random effect
    PARAMETER(log_sigma_eta); // sd on cluster (nugget)

    // Random Effects
    PARAMETER_ARRAY(z); // Space-time random field
    PARAMETER_VECTOR(epsilon); // unstructured time effects
    PARAMETER_ARRAY(phi); // Age specific time effects
    PARAMETER_VECTOR(nu); // survey level effects
    PARAMETER_VECTOR(eta); // cluster level nugget

    std::cout << "Variables Loaded.\n";

    int Npoint = yPoint.size();
    int Npoly = yPoly.size();

    // Hyper parameter transformations
    Type tau = exp(log_tau);
    Type kappa = exp(log_kappa);
    Type rho = exp(logit_rho) / (Type(1.) + exp(logit_rho));
    vector<Type> sigma_phi = exp(log_sigma_phi);
    Type sigma_nu = exp(log_sigma_nu);
    Type sigma_epsilon = exp(log_sigma_epsilon);
    Type sigma_eta = exp(log_sigma_eta);

    //fixed effect transforms
    vector<Type> beta_agec(beta_age.size() + 1);
    beta_agec[0] = Type(0.0);
    for(int i=1; i < beta_agec.size(); i++){
        beta_agec[i] = beta_age[i-1];
    }

    std::cout << "Variables Transformed.\n";

    Type nll = 0.0;

    // Apply fixed effect and hyper-parameter priors
    if(priors == 1){
        Type tempSum;
        for(int b=0; b<beta.size(); b++){
            nll -= dnorm(beta[b], Type(0.0), Type(10.0), true);
        }
        for(int p=0; p<log_sigma_phi.size(); p++){
            nll -= dnorm(log_sigma_phi[p], Type(0.0), Type(2.), true);
        }
        for(int i=0; i<phi.dim(0); i++){
            tempSum = Type(0.0);
            for(int j=0; j<phi.dim(1); j++){
                tempSum += phi(i,j);
            }
            nll -= dnorm(tempSum, Type(0.0), Type(.001), true);
        }
        nll -= dnorm(log_tau, Type(0.0), Type(10.0), true);
        nll -= dnorm(logit_rho, Type(2.0), Type(5.0), true);
        nll -= dnorm(log_kappa, Type(0.0), Type(10.0), true);
        nll -= dnorm(log_sigma_eta, Type(0.0), Type(10.0), true);
        nll -= dnorm(log_sigma_epsilon, Type(0.0), Type(10.0), true);
        nll -= dnorm(log_sigma_nu, Type(0.0), Type(10.0), true);
    }

    std::cout << "Fixed Effect priors applied.\n";
    // apply random effect structure priors
    SparseMatrix<Type> Q = spde2_Q(log_kappa, log_tau, M0, M1, M2);

    if(z.dim(1) > 1){
        nll += SEPARABLE(AR1(rho), GMRF(Q))(z);
    }
    else{
        nll += GMRF(Q)(z.matrix().col(0));
    }

    std::cout << "SPDE priors applied.\n";

    for(int i=0; i<nu.size(); i++){
        nll -= dnorm(nu[i], Type(0.0), sigma_nu, true);
    }

    for(int i=0; i<epsilon.size(); i++){
        nll -= dnorm(epsilon[i], Type(0.0), sigma_epsilon, true);
    }

    for(int i=0; i<eta.size(); i++){
        nll -= dnorm(eta[i], Type(0.0), sigma_eta, true);
    }

    std::cout << "Basic iid priors applied.\n";

    vector<Type> aget(phi.dim(1));
    for(int i=0; i<phi.dim(0); i++){
        for(int j=0; j<phi.dim(1); j++){
            aget[j] = phi(i,j);
        }
        //nll += densRW2(aget, sigma_phi[i]);
        nll += SCALE(AR1(Type(.99)), sigma_phi[i])(aget);
    }

    std::cout << "RW2 priors applied.\n";

    // Turn z into a matrix here so projLatF is space by time
    matrix<Type> projLatF = AprojObs * z.matrix();
    std::cout << "Project space time.\n";
    // covs should be an array that we loop through time to have space by time
    matrix<Type> projCov(projLatF.rows(), projLatF.cols());
    for(int s=0; s<projLatF.rows(); s++){
        for(int t=0; t<projLatF.cols(); t++){
            projCov(s,t) = epsilon(t);
            for(int b=0; b<beta.size(); b++){
                projCov(s,t) += covs(s,t,b) * beta(b);
            }
        }
    }

    std::cout << "Betas added.\n";

    matrix<Type> projLatObs = projLatF + projCov;
    Type nllPart;
    Type p;
    Type logitp;

    for(int i=0; i<Npoint; i++){
        logitp = projLatObs(idPoint[i],idtPoint[i]) + \
            phi(idaPoint[i], idtPoint[i]) + nu[idnPoint[i]] + \
            beta_agec[idaPoint[i]] + eta[idcPoint[i]];
        p = exp(logitp) / (Type(1.) + exp(logitp));
        nll -= dbinom(Type(yPoint[i]), Type(denomPoint[i]), p, true);
    }

    std::cout << "Point data added.\n";

    if(Npoly != 0){
        // mixture model estimation
        if(moption == 0){
            for(int i=0; i<Npoly; i++){
                nllPart = Type(.0);
                for(typename SparseMatrix<Type>::InnerIterator it(AprojPoly,idPoly[i]); it; ++it){
                    // see the following link for indexing sparse matrices
                    // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title2
                    logitp = projLatObs(it.row(),idtPoly[i]) + \
                        phi(idaPoly[i], idtPoly[i]) + nu[idnPoly[i]] + \
                        beta_agec[idaPoly[i]] + eta[idcPoly[i]];
                    nllPart += dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, false) * it.value();
                }
                // std::cout << log(nllPart);
                // std::cout << "\n";
                nll -= log(nllPart);
            }
        }
        // Redistribute points
        if(moption == 1){
            for(int i=0; i<Npoly; i++){
                logitp = projLatObs(idPoly[i],idtPoly[i]) + \
                    phi(idaPoly[i], idtPoly[i]) + nu[idnPoly[i]] + \
                    beta_agec[idaPoly[i]] + eta[idcPoly[i]];
                p = exp(logitp) / (Type(1.) + exp(logitp));
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
        if(moption == 3){
            // Reimann sum approximation
            SparseMatrix<Type> RAprojPoly = AprojPoly.transpose();
            for(int i=0; i<Npoly; i++){
                vector<Type> projLatObsFull = projLatObs.col(idtPoly[i]).array() + \
                    phi(idaPoly[i], idtPoly[i]) + nu[idnPoly[i]] + \
                    beta_agec[idaPoly[i]] + eta[idcPoly[i]];
                vector<Type> projPObs = projLatObsFull.array().exp() / \
                    (Type(1.) + projLatObsFull.array().exp());
                vector<Type> projPoly = RAprojPoly * projPObs;
                p = projPoly(idPoly[i]);
                nll -= dbinom(Type(yPoly[i]), Type(denomPoly[i]), p, true);
            }
        }
    }

    REPORT(z);
    return nll;
}

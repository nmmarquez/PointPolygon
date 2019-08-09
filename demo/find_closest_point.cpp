#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
double compute_a(double lat1, double long1, double lat2, double long2) {
    
    double sin_dLat = ::sin((lat2 - lat1) / 2);
    double sin_dLon = ::sin((long2 - long1) / 2);
    
    return sin_dLat * sin_dLat + ::cos(lat1) * ::cos(lat2) * sin_dLon * sin_dLon;
}

int find_min(double lat1, double long1,
             const NumericVector& lat2,
             const NumericVector& long2,
             int current0) {
    
    int m = lat2.size();
    double lat_k, lat_min, lat_max, a, a0;
    int k, current = current0;
    
    a0 = compute_a(lat1, long1, lat2[current], long2[current]);
    // Search before current0
    lat_min = lat1 - 2 * ::asin(::sqrt(a0));
    for (k = current0 - 1; k >= 0; k--) {
        lat_k = lat2[k];
        if (lat_k > lat_min) {
            a = compute_a(lat1, long1, lat_k, long2[k]);
            if (a < a0) {
                a0 = a;
                current = k;
                lat_min = lat1 - 2 * ::asin(::sqrt(a0));
            }
        } else {
            // No need to search further
            break;
        }
    }
    // Search after current0
    lat_max = lat1 + 2 * ::asin(::sqrt(a0));
    for (k = current0 + 1; k < m; k++) {
        lat_k = lat2[k];
        if (lat_k < lat_max) {
            a = compute_a(lat1, long1, lat_k, long2[k]);
            if (a < a0) {
                a0 = a;
                current = k;
                lat_max = lat1 + 2 * ::asin(::sqrt(a0));
            }
        } else {
            // No need to search further
            break;
        }
    }
    
    return current;
} 

//' Find index of closest neighbor
//'
//' cpp code for finding closest neighbor
//'
//' @param lat1
//' @param long1
//' @param lat2
//' @param long2
//' @return integer vector indicating closest neighbor
//' @export
//' @useDynLib PointPolygon
// [[Rcpp::export]]
IntegerVector find_closest_point(const NumericVector& lat1,
                                 const NumericVector& long1,
                                 const NumericVector& lat2,
                                 const NumericVector& long2) {
    
    int n = lat1.size();
    IntegerVector res(n);
    
    int current = 0;
    for (int i = 0; i < n; i++) {
        res[i] = current = find_min(lat1[i], long1[i], lat2, long2, current);
    }
    
    return res; // need +1
}
#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <boost/math/distributions/fisher_f.hpp>
using namespace std;

#define INF 1e9

double compute_p_value(double f_stat, double dfbg, double dfwg);

struct AnovaResult {
    double f_statistic;
    double p_value;
};

AnovaResult f_oneway(const std::vector<std::vector<double>>& groups) {
    // Variables to hold sum of squares
    double ss_total = 0.0;
    double ss_between_groups = 0.0;
    double ss_within_groups = 0.0;
    
    // Total number of observations
    double n_total = 0;
    for (const auto& group : groups) {
        n_total += group.size();
    }
    
    // Overall mean
    double overall_mean = 0.0;
    for (const auto& group : groups) {
        overall_mean += std::accumulate(group.begin(), group.end(), 0.0);
    }
    overall_mean /= n_total;
    
    // Calculating sum of squares total (SST)
    for (const auto& group : groups) {
        for (double val : group) {
            ss_total += pow(val - overall_mean, 2);
        }
    }
    
    // Calculating sum of squares between groups (SSBG)
    for (const auto& group : groups) {
        double group_mean = std::accumulate(group.begin(), group.end(), 0.0) / group.size();
        ss_between_groups += group.size() * pow(group_mean - overall_mean, 2);
    }
    
    // Calculating sum of squares within groups (SSWG)
    ss_within_groups = ss_total - ss_between_groups;
    
    // Degrees of freedom
    double df_between_groups = groups.size() - 1;
    double df_within_groups = n_total - groups.size();
    
    // Mean squares
    double ms_between_groups = ss_between_groups / df_between_groups;
    double ms_within_groups = ss_within_groups / df_within_groups;
    
    // F-statistic
    if (ms_within_groups == 0) {
        cout << "Warning: Mean squares within groups is zero." << endl;
        return AnovaResult{INF, 0.0};
    }
    double f_statistic = ms_between_groups / ms_within_groups;

    double p_value = compute_p_value(f_statistic, df_between_groups, df_within_groups);
    
    return AnovaResult{f_statistic, p_value};
}

double compute_p_value(double f_stat, double dfbg, double dfwg) {
    // Create an F-distribution object with the given degrees of freedom
    boost::math::fisher_f dist(dfbg, dfwg);

    // Compute the p-value as the complement of the CDF at the F-statistic
    double p_value = boost::math::cdf(boost::math::complement(dist, f_stat));
    return p_value;
}

#include <vector>
#include <algorithm>
#include <utility> // For std::pair

std::vector<double> compute_bh_correction(const std::vector<double>& p_values) {
    size_t n = p_values.size();
    std::vector<double> bh_correction(n);
    std::vector<std::pair<double, size_t>> p_values_with_index(n);

    // Associate each p-value with its index
    for (size_t i = 0; i < n; ++i) {
        p_values_with_index[i] = {p_values[i], i};
    }

    // Sort p-values in ascending order with their indices
    std::sort(p_values_with_index.begin(), p_values_with_index.end());

    // Adjust p-values using the BH procedure
    double min_adjusted_pval = 1.0;
    for (int i = n - 1; i >= 0; --i) {
        double original_pval = p_values_with_index[i].first;
        size_t original_index = p_values_with_index[i].second;
        double adjusted_pval = std::min(min_adjusted_pval, (original_pval * n) / (i + 1));
        bh_correction[original_index] = adjusted_pval;
        min_adjusted_pval = adjusted_pval; // Update min_adjusted_pval to enforce monotonicity
    }

    return bh_correction;
}

#include <Rcpp.h>
// [[Rcpp::depends(BH)]]
#include <cmath>
#include <boost/math/distributions.hpp>
#include <vector>
#include <algorithm>

bool compareFunction(const std::pair<double, int>& a, const std::pair<double, int>& b)
{
  return a.first < b.first;
}

std::vector<std::pair<double, int>> cpugetpairvector(const std::vector<double>& vect1, const std::vector<double>& vect2)
{
  if (vect1.size() != vect2.size()) {
    Rcpp::stop("Log-rank times and event indicators must have the same length in each group.");
  }
  std::vector<std::pair<double, int>> target(vect1.size());
  for (std::size_t i = 0; i < target.size(); i++) {
    if (!std::isfinite(vect1[i]) || vect1[i] < 0) {
      Rcpp::stop("Log-rank times must be finite and non-negative.");
    }
    if (vect2[i] != 0 && vect2[i] != 1) {
      Rcpp::stop("Log-rank event indicators must be 0 or 1, with no missing values.");
    }
    target[i] = std::make_pair(vect1[i], static_cast<int>(vect2[i]));
  }
  return target;
}

// [[Rcpp::export]]
std::vector<double> logrank_instance(const std::vector<double>& groupa, const std::vector<double>& groupb, const std::vector<double>& groupacensored, const std::vector<double>& groupbcensored, SEXP onlyz){

  if (TYPEOF(onlyz) != LGLSXP || Rf_xlength(onlyz) != 1 || LOGICAL(onlyz)[0] == NA_LOGICAL) {
    Rcpp::stop("'onlyz' must be a single non-missing logical value.");
  }
  const bool return_only_z = LOGICAL(onlyz)[0];

  std::vector<std::pair<double, int>> groupaboth = cpugetpairvector(groupa, groupacensored);
  std::vector<std::pair<double, int>> groupbboth = cpugetpairvector(groupb, groupbcensored);

  std::sort(groupaboth.begin(), groupaboth.end(), compareFunction);
  std::sort(groupbboth.begin(), groupbboth.end(), compareFunction);

  unsigned long long na = groupaboth.size();
  unsigned long long nb = groupbboth.size();
  unsigned long long nasize = na;
  unsigned long long nbsize = nb;
  unsigned long long n = na + nb;
  long double oa = 0.0;
  long double ea = 0.0;
  long double var = 0.0;
  unsigned long long indexa = 0;
  unsigned long long indexb = 0;

  while (na > 0 && nb > 0) {
    double value = groupaboth[indexa].first;
    unsigned long long count = 0;
    unsigned long long countb = 0;
    unsigned long long minusna = 0;
    unsigned long long minusnb = 0;

    if (indexb < nbsize && groupbboth[indexb].first < value) {
      ++minusnb;
      value = groupbboth[indexb].first;
      if (groupbboth[indexb].second == 1) {
        ++countb;

      }
      ++indexb;
      while (indexb < nbsize && value == groupbboth[indexb].first) {
        if (groupbboth[indexb].second == 1) {
          ++countb;
        }
        ++indexb;
        ++minusnb;
      }
    }
    else if (indexb < nbsize && groupbboth[indexb].first == value) {
      ++minusna;
      if (groupaboth[indexa].second == 1) {
        ++count;

      }
      ++indexa;
      while (indexa < nasize && value == groupaboth[indexa].first) {
        if (groupaboth[indexa].second == 1) {
          ++count;
        }
        ++indexa;
        ++minusna;
      }
      while (indexb < nbsize && value == groupbboth[indexb].first) {
        if (groupbboth[indexb].second == 1) {
          ++countb;
        }
        ++indexb;
        ++minusnb;
      }
    }
    else {
      ++minusna;
      if (indexa < nasize && groupaboth[indexa].second == 1) {
        ++count;

      }
      ++indexa;
      while (indexa < nasize && value == groupaboth[indexa].first) {
        if (groupaboth[indexa].second == 1) {
          ++count;
        }
        ++indexa;
        ++minusna;
      }
    }

    const long double risk_a = static_cast<long double>(na);
    const long double risk_b = static_cast<long double>(nb);
    const long double risk_total = static_cast<long double>(n);
    const long double events = static_cast<long double>(count + countb);
    ea += risk_a * events / risk_total;
    oa += count;
    var += risk_a * risk_b * events * (risk_total - events) /
      (risk_total * risk_total * (risk_total - 1));
    na -= minusna;
    nb -= minusnb;
    n = na + nb;
  }
  if (var <= 0) {
    if (return_only_z) { return std::vector<double>{R_NaN}; }
    return std::vector<double>{R_NaN, R_NaN, R_NaN};
  }
  long double zstat = (oa-ea)/std::sqrt(var);
  if (return_only_z) { return std::vector<double>{static_cast<double>(zstat)}; }
  else {
    long double logrank = ((oa - ea) * (oa - ea)) / var;
    double pvalue = boost::math::cdf(boost::math::complement(boost::math::chi_squared(1), logrank));
    return std::vector<double>{static_cast<double>(logrank),static_cast<double>(zstat),pvalue};
  }
}

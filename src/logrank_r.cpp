#include <Rcpp.h>
// [[Rcpp::depends(BH)]]
#include <cmath>
#include <boost/math/distributions.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <thread>

std::vector<std::vector<double>> res;
std::vector<std::vector<double>> groupas;
std::vector<std::vector<int>> groupacensoreds;
std::vector<std::vector<double>> groupbs;
std::vector<std::vector<int>> groupbcensoreds;
bool onlyz;

bool compareFunction(std::pair<double, int>& a, std::pair<double, int>& b)
{
  return a.first < b.first;
}

std::vector<std::pair<double, int>> cpugetpairvector(std::vector<double>& vect1, std::vector<int>& vect2)
{
  std::vector<std::pair<double, int>> target(vect1.size());
  for (unsigned long long i = 0; i < target.size(); i++)
    target[i] = std::make_pair(vect1[i], vect2[i]);
  return target;
}

// [[Rcpp::export]]
std::vector<double> logrank_instance(std::vector<double>& groupa, std::vector<double>& groupb, std::vector<int>& groupacensored, std::vector<int>& groupbcensored, bool onlyz){

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

    ea += static_cast<long double>(na * (count + countb)) / static_cast<long double>(n);
    oa += count;
    var += (static_cast<long double>(na) * static_cast<long double>(nb * (count + countb)) * static_cast<long double>(n - (count + countb))) / (static_cast<long double>(n) * static_cast<long double>(n) * static_cast<long double>(n - 1));
    na -= minusna;
    nb -= minusnb;
    n = na + nb;
  }
  long double zstat = (oa-ea)/sqrt(var);
  if (onlyz) { return std::vector<double>(static_cast<double>(zstat)); }
  else {
    long double logrank = ((oa - ea) * (oa - ea)) / var;
    double pvalue = 1.0;
    if (logrank >= 0) {
      pvalue = (1.0 - boost::math::cdf(boost::math::chi_squared(1), logrank));
      //std::cout << pvalue << std::endl;
    }
    return std::vector<double>{static_cast<double>(logrank),static_cast<double>(zstat),pvalue};
  }
}

void startthread(unsigned long long start, unsigned long long end) {
  for (unsigned long long j = start; j < end; ++j) {
    res[j]=logrank_instance(groupas[j], groupbs[j], groupacensoreds[j], groupbcensoreds[j], onlyz);
  }
}

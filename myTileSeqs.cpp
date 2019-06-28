//[[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace std;

// Stripped down version of TileSeqs from DECIPHER
//[[Rcpp::export]]
DataFrame myTileSeqs(vector<string> identifier, vector<string> seqs,
                     size_t minLength = 26, size_t maxLength = 27) {
  if (identifier.size() != seqs.size()) {
    Rcpp::stop("identifier must have the same length as seqs.\n");
  }
  // Convert to unordered map for faster search
  std::unordered_map<std::string, std::string> seq_db;
  seq_db.reserve(identifier.size());
  for (size_t i = 0; i < identifier.size(); i++) {
    seq_db[identifier[i]] = seqs[i];
  }
  
  vector<string> id, target_site;
  vector<size_t> start, end, width;
  vector<bool> misprimes; 
  id.reserve(1400 * identifier.size());
  target_site.reserve(1400 * identifier.size());
  start.reserve(1400 * identifier.size());
  end.reserve(1400 * identifier.size());
  width.reserve(1400 * identifier.size());
  misprimes.reserve(1400 * identifier.size());
  
  size_t j = 0;
  Progress p(identifier.size(), true);
  // Generate the tiles
  for (auto k : seq_db) {
    string target = k.second;
    //if (target.length() < maxLength) {
   //   Rcpp::warning("Skipped because sequences shorter than maxLength: ", k);
   //   continue;
    //}
    size_t tl = target.length();
    size_t l = tl - maxLength + 1;
    size_t s = 0; // start
    for (size_t i = 0; i < l; i++) {
      if (j % 10000 == 0) {
        Rcpp::checkUserInterrupt();
      }
      bool misprime = false;
      string ts = target.substr(s, maxLength);
      
      // Check for misprime
      if (ts.length() < minLength || ts.length() > maxLength) {
        continue;
      }
      size_t repeats = 0, runs = 0;
      for (size_t p = 1; p < ts.length(); p++) {
        if (p > 2) {
          if ((ts[p]==ts[p - 2]) && (ts[p - 1]==ts[p - 3])) {
            repeats++;
            if (repeats > 5) {
              misprime = true;
              break;
            }
          } else {
            repeats = 0;
          }
        }
        if (ts[p] == ts[p - 1]) {
          runs++;
          if (runs > 3) {
            misprime = true;
            break;
          }
        } else {
          runs = 0;
        }
      }
      target_site.push_back(ts);
      id.push_back(k.first);
      start.push_back(s + 1);
      end.push_back(s + maxLength);
      width.push_back(tl);
      misprimes.push_back(misprime);
      s++;
      j++;
    }
    p.increment();
  }
  vector<size_t> coverage(target_site.size(), 1);
  // Fill in the data frame
  return DataFrame::create(_["start"] = start,
                           _["end"] = end,
                           _["start_aligned"] = start,
                           _["end_aligned"] = end,
                           _["misprime"] = misprimes,
                           _["width"] = width,
                           _["id"] = id,
                           _["coverage"] = coverage,
                           _["groupCoverage"] = coverage,
                           _["target_site"] = target_site);
}

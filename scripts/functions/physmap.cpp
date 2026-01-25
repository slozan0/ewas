// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/algorithm/string.hpp>
using namespace Rcpp;

namespace MosSci
{

  inline std::string IdNuc(const std::string &lChars)
  {

    if (lChars.size() == 1)
    {
      if (lChars[0] == 'A')
        return "A";
      if (lChars[0] == 'C')
        return "C";
      if (lChars[0] == 'G')
        return "G";
      if (lChars[0] == 'T')
        return "T";
      return "";
    }

    if (lChars.rfind("INS", 0) == 0)
      return "I";
    if (lChars.rfind("DEL", 0) == 0)
      return "D";

    return "";
  }

  CharacterVector GetMonoSites(CharacterVector line)
  {
    Rcpp::CharacterVector outline(9);

    outline[0] = line[0];
    outline[1] = line[1];
    outline[2] = line[2];

    if (line[2] == "A")
    {
      outline[3] = line[4];
    }
    else if (line[2] == "C")
    {
      outline[4] = line[4];
    }
    else if (line[2] == "G")
    {
      outline[5] = line[4];
    }
    else if (line[2] == "T")
    {
      outline[6] = line[4];
    }
    else
    { // do nothing
    }

    return outline;
  }

  CharacterVector GetPolySites(CharacterVector line)
  {
    CharacterVector outline(9);

    std::vector<std::string> result;

    // int nExpectedCols = 6;
    // count columns were the values are not empty
    int nColumns = line.length();
    int nNucs = 0;

    for (int j = 5; j < nColumns; j++)
    {
      if (line[j].size() > 0)
        nNucs = nNucs + 1;
    }

    NumericVector myA(nNucs);
    NumericVector myC(nNucs);
    NumericVector myG(nNucs);
    NumericVector myT(nNucs);
    NumericVector myI(nNucs);
    NumericVector myD(nNucs);

    //()row, column)
    StringMatrix nucTable(nNucs, 4);
    // keep track of errors
    int warn_count = 0;
    const int warn_max = 5;

    // create a table for each unique identifier
    for (int j = 0; j < nNucs; j++)
    {

      if (line[j + 5] == NA_STRING)
      {
        if (warn_count < warn_max)
        {
          Rcpp::warning("GetPolySites: NA in line[%d] column[%d] (j+5). Skipping.",
                        /*row?*/ -1, j + 5);
          warn_count++;
        }
        continue;
      }

      std::string cell = Rcpp::as<std::string>(line[j + 5]);
      if (cell.empty())
      {
        if (warn_count < warn_max)
        {
          Rcpp::warning("GetPolySites: empty string in column[%d]. Skipping.", j + 5);
          warn_count++;
        }
        continue;
      }

      result.clear();
      // base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads
      boost::split(result, cell, [](char c)
                   { return c == ':'; });

      if (result.size() < 2)
      {
        if (warn_count < warn_max)
        {
          Rcpp::warning("GetPolySites: malformed token '%s' in column[%d] (expected id:reads:...). Skipping.",
                        cell.c_str(), j + 5);
          warn_count++;
        }
        continue;
      }

      //"id" , "qa", "nuc" , "reads"
      nucTable(j, 0) = j;
      nucTable(j, 1) = result[0];
      nucTable(j, 2) = MosSci::IdNuc(result[0]);
      nucTable(j, 3) = result[1];
    }

    // sum across rows with the same values
    for (int i = 0; i < nNucs; i++)
    {
      if (nucTable(i, 2) == "A")
      {
        myA[i] = atoi(nucTable(i, 3));
      }
      else if (nucTable(i, 2) == "C")
      {
        myC[i] = atoi(nucTable(i, 3));
      }
      else if (nucTable(i, 2) == "G")
      {
        myG[i] = atoi(nucTable(i, 3));
      }
      else if (nucTable(i, 2) == "T")
      {
        myT[i] = atoi(nucTable(i, 3));
      }
      else if (nucTable(i, 2) == "I")
      {
        myI[i] = atoi(nucTable(i, 3));
      }
      else if (nucTable(i, 2) == "D")
      {
        myD[i] = atoi(nucTable(i, 3));
      }
      else
      {
        // do nothing
      }
    }

    int sumA = sum(myA);
    int sumC = sum(myC);
    int sumG = sum(myG);
    int sumT = sum(myT);
    int sumI = sum(myI);
    int sumD = sum(myD);

    outline[0] = line[0];
    outline[1] = line[1];
    outline[2] = line[2];

    outline[3] = std::to_string(sumA);
    outline[4] = std::to_string(sumC);
    outline[5] = std::to_string(sumG);
    outline[6] = std::to_string(sumT);
    outline[7] = std::to_string(sumI);
    outline[8] = std::to_string(sumD);

    // chr, position, a, c, g, t, i, d

    return (outline);
  }
}

// [[Rcpp::export]]
CharacterMatrix MapMonoSites(CharacterMatrix monoSites)
{

  int nRows = monoSites.nrow();
  CharacterMatrix m3(nRows, 9);

  for (int i = 0; i < nRows; i++)
  {
    m3.row(i) = MosSci::GetMonoSites(monoSites.row(i));
  }

  return m3;
}

// [[Rcpp::export]]
CharacterMatrix MapPolySites(CharacterMatrix polySites)
{

  int nRows = polySites.nrow();
  CharacterMatrix m3(nRows, 9);

  for (int i = 0; i < nRows; i++)
  {

    m3.row(i) = MosSci::GetPolySites(polySites.row(i));
  }

  return m3;
}

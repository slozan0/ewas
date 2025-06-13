// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/algorithm/string.hpp>

using namespace Rcpp;

namespace MosSci{

  String IdNuc(CharacterVector inChars) {
  
    std::string lChars = Rcpp::as<std::string>(inChars); //transform charachtervector to c++ string object
  
    String returnVar;
  
    int nChars = lChars.size();
  
    if(nChars == 1){
      if (lChars.find("A") == 0){
        returnVar = "A";
      } else if (lChars.find("C") == 0) {
        returnVar = "C";
      } else if (lChars.find("G") == 0) {
        returnVar = "G";
      } else if (lChars.find("T") == 0) {
        returnVar = "T";
      } else {
      //do nothing
      }
    } else {
      if( lChars.find("INS" ) == 0) {
        returnVar = "I";
      } else if( lChars.find("DEL") == 0) {
        returnVar = "D";
      }
    }
    return(returnVar);
  }
  CharacterVector GetMonoSites(CharacterVector line) {
    Rcpp::CharacterVector outline(9);

    outline[0] = line[0];
    outline[1] = line[1];
    outline[2] = line[2];
  
         if (line[2] == "A") { outline[3] = line[4]; } 
    else if (line[2] == "C") { outline[4] = line[4]; } 
    else if (line[2] == "G") { outline[5] = line[4]; }
    else if (line[2] == "T") { outline[6] = line[4]; }
    else { //do nothing
      }

    return outline;
  }
  
  CharacterVector GetPolySites(CharacterVector line) {
    CharacterVector outline(9);
  
    std::vector<std::string> result;
    
    //int nExpectedCols = 6;
    //count columns were the values are not empty
    int nColumns = line.length();
    int nNucs = 0;
    
    for( int j = 5; j < nColumns; j++){
      if ( line[j].size() > 0 ) nNucs = nNucs + 1;
    }
    
    NumericVector myA ( nNucs );
    NumericVector myC ( nNucs );
    NumericVector myG ( nNucs );
    NumericVector myT ( nNucs );
    NumericVector myI ( nNucs );
    NumericVector myD ( nNucs );
  
    //()row, column)
    StringMatrix nucTable( nNucs, 4);
  
    //create a table for each unique identifier
    for( int j = 0; j < nNucs; j++){
      //base:reads:strands:avg_qual:map_qual:plus_reads:minus_reads
      boost::split(result, line[j + 5], [](char c){return c == ':';});
    
      //"id" , "qa", "nuc" , "reads" 
      nucTable(j, 0) = j;
      nucTable(j, 1) = result[0];
      nucTable(j, 2) = IdNuc( result[0] );
      nucTable(j, 3) = result[1];
    }
  
  //sum across rows with the same values
    for( int i = 0; i < nNucs; i++){
      if        ( nucTable(i, 2) == "A") {
        myA[i] = atoi( nucTable(i, 3) );
      } else if ( nucTable(i, 2) == "C") {
        myC[i] = atoi( nucTable(i, 3) );
      } else if ( nucTable(i, 2) == "G") {
        myG[i] = atoi( nucTable(i, 3) );
      } else if ( nucTable(i, 2) == "T") {
        myT[i] = atoi( nucTable(i, 3) );
      } else if ( nucTable(i, 2) == "I") {
        myI[i] = atoi( nucTable(i, 3) );
      } else if ( nucTable(i, 2) == "D") {
        myD[i] = atoi( nucTable(i, 3) );
      } else {
      //do nothing
      }
    }
  
    int sumA = sum( myA );
    int sumC = sum( myC );
    int sumG = sum( myG );
    int sumT = sum( myT );
    int sumI = sum( myI );
    int sumD = sum( myD );
  
    outline[0] = line[0];
    outline[1] = line[1];
    outline[2] = line[2];
  
    outline[3] = std::to_string( sumA );
    outline[4] = std::to_string( sumC );
    outline[5] = std::to_string( sumG );
    outline[6] = std::to_string( sumT );
    outline[7] = std::to_string( sumI );
    outline[8] = std::to_string( sumD ); 
  
    //chr, pos, a, c, g, t, i, d
  
    return(outline);
  }
}

// [[Rcpp::export]]
CharacterMatrix MapMonoSites(CharacterMatrix monoSites) {
  
  int nRows = monoSites.nrow();
  CharacterMatrix m3( nRows, 9);

  for( int i = 0; i < nRows; i++) {
    m3.row(i) = MosSci::GetMonoSites( monoSites.row(i) );
  }
  
  return m3;
}
// [[Rcpp::export]]
CharacterMatrix MapPolySites(CharacterMatrix polySites) {
  
  int nRows = polySites.nrow();
  CharacterMatrix m3( nRows, 9);

  for( int i = 0; i < nRows; i++) {

   m3.row(i) = MosSci::GetPolySites( polySites.row(i) );

  }
  
  return m3;
}

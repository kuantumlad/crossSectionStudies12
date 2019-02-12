#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>

std::vector<std::string> parseLine(std::string line){

  std::vector<std::string> tokenizedLine;
  std::string buffer;

  std::istringstream tokens(line);
  while(!tokens.eof())
    {
      std::string buffer;
      tokens >> buffer;
      tokenizedLine.push_back(buffer);
    }

  return tokenizedLine;
}

std::map<int, std::vector<double> > loadECParameters(std::string inputFilename ){

  std::map<int, std::vector<double> > ParametersOut;

  std::ifstream file(inputFilename.c_str());
  
  std::string currentLine;
  int counter = 0;
  std::vector< double > fit_parameters;

  while(getline(file, currentLine)){
    // Trying to skip whitespace lines
    if (currentLine.size() > 10 && currentLine[0] != '#'){
      

      std::vector<std::string> line = parseLine(currentLine);
      
      int numberOfValues = line.size()-3;
      std::vector< double > fit_para_sect;
      //std::cout << " >> " << numberOfValues<< std::endl;
      for (int i=0; i<numberOfValues; i++){ 
	fit_parameters.push_back( atof(line[i+3].c_str()) ); 
	//std::cout << " >> " <<  atof(line[i+3].c_str()) << std::endl;
      } // Setting values	  
      
      ParametersOut[counter] = fit_parameters;
      counter++;
      //std::cout << " COUNTER " << counter << " FIT PARAMETERS VECTOR SIZE " << fit_parameters.size() << std::endl;
      fit_parameters.clear();

    }
  }
  file.close();


  return ParametersOut;
}



#include <iostream>
#include <fstream>
#include <sstream>

int NrConfigs = 9;
int nEvts = 10; 

//std::ifstream infile ("weights.txt", std::ifstream::in); 
//char test = infile.get();
//std::cout << " infile.get : " << test << std::endl;

//vector<float> readFromFile(int EvtNr){
//void readFromFile(int EvtNr, int ConfigNr){
float readFromFile(int EvtNr, int ConfigNr){
  //float weightVector[][]; //[nEvts][nrConfigs];
  //vector<float> weightVector;
  //vector<float> weightUncVector[nEvts];
  //for(int jj = 0; jj < nEvts; jj++){
  //  for(int ii = 0; ii < NrConfigs; ii++){
  //    //weightVector[jj].push_back(0);
  //    weightVector[jj][ii] = 0;
  //    //weightUncVector[jj].push_back(0);
  //  }
  //}

  std::ifstream ifs ("weights.txt", std::ifstream::in);
  std::string line;
  int evt,config,tf;
  float weight, weightUnc;
  bool weightFound = 0;
  //while( ifs >> evt >> config >> tf >> weight >> weightUnc ){
  std::cout << " Output of file.eof is : " << ifs.eof() << std::endl;
  while( std::getline(ifs,line) ){
    std::istringstream iss(line);
    if( iss >> evt >> config >> tf >> weight >> weightUnc){
      std::cout << " Looking at event number : " << evt << std::endl;
      if( evt == EvtNr && config == ConfigNr){
        //weightVector[evt-1][config-1] = weight;
        //weightUncVector[evt-1][config-1] = weightUnc;
        //std::cout << " Looking at line : " << line << std::endl;
        std::cout << " Should be storing value : " << evt-1 << " & " << config-1 << std::endl;
        std::cout << " Weight value is : " << weight << std::endl;
        return weight;
        weightFound = 1;
        continue;
        //weightVector.push_back(weight);
        //std::cout << " Stored weightVector is : " << weightVector[config-1] << std::endl;
        //std::cout << " Looking at evt : " << evt << std::endl;
      }
      else
        weight = 0;
    }
  }
  std::cout << " weight value is : " << weight << std::endl;
  std::cout << " Output of file.eof is : " << ifs.eof() << std::endl;

  ifs.close();
  std::cout << " Output of file.eof is : " << ifs.eof() << std::endl;

  if( weightFound == 0) return weight;
  //return weightVector[EvtNr];

  /*
  vector<float> weightVector[nEvts];
  vector<float> weightUncVector[nEvts];
  for(int jj = 0; jj < nEvts; jj++){
    for(int ii = 0; ii < NrConfigs; ii++){
      weightVector[jj].push_back(0);
      weightUncVector[jj].push_back(0);
    }
  }
  
  std::string line; 
  while (std::getline(infile, line)){
    std::istringstream iss(line);
    int evt,config,tf;
    float weight, weightUnc;
    if ( (iss >> evt >> config >> tf >> weight >> weightUnc)){ //Only consider the case with the numbers!
      weightVector[evt-1][config-1] = weight;
      weightUncVector[evt-1][config-1] = weightUnc;
    }
  }
  infile.close()
  return weightVector[EvtNr];   //No idea how to return the uncertainty vector then ...
  */
}

int outputTest(){

 
  value = readFromFile(9,2); 
  std::cout << " Obtained value is : " << value << std::endl;
  secondValue = readFromFile(2,2);
  std::cout << " Obtained second value is : " << secondValue << std::endl;
  //vector<float> V = readFromFile(2);  //Acces the weightsvector for a given event number!
  //for(unsigned int ii = 0; ii < V.size(); ii++){
  //  std::cout << ii << " ) " << V[ii] << std::endl;
  //}
}


//Run this before FindCSD.cpp . Run it once on each data set for which you want to find correlations. Correlations are output to RhoAndVar.txt - remember to rename the file before running the next iteration, or your original data will be overwritten. 

#include <ctime>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <list>
#include <fstream>
#include <cstring>
#include <stdlib.h>


using namespace std;


//Parameters depending on input file

const char* expDataFile = "PCG_Contr.txt"; //Name of expression data file
const char* outFile = "PCG_Controut.txt"; //Name of output data file
const int sampleSize = 24; //Number of data points per gene (normally number of individuals from which data is collected, corresponds to columns in the expression data text file)
const int numberOfGenes = 21044;// Number of distinct genes for which there is expression data (corresponds to rows in expression data)
const int subSampleSize = 2;// of subsamples for determination of variance in co-expression. 10 is a good minimum - can be increased if sampleSize is very large.                                                                      


struct Pair;


struct Pair
{
  double x;
  double y;
  double xRank;
  double yRank;
};



double spearman(list<Pair> subSample);

int main()
{
  srand(time(NULL));

  ifstream inStream;
  ofstream outStream;
  
  double expressionValues[sampleSize][numberOfGenes];
  string geneName[numberOfGenes];
  
  inStream.open(expDataFile);
  
  if(inStream.is_open())
    {
      char temp[4000000];
      int i = 0;
      while(!inStream.getline(temp, 4000000).eof() && i < numberOfGenes)
	{     
	  inStream >> geneName[i];

	  for (int k = 0; k < sampleSize; k++)
	    {
	      inStream >> expressionValues[k][i];

	    }

	  i++;
	}
      

    }
  


  int newSequence[sampleSize];
  bool alreadyTaken[sampleSize];
  for (int i = 0; i < sampleSize; i++)
    {
      alreadyTaken[i] = 0;
    }
  for (int i = 0; i < sampleSize; i++)
    {
      int r = rand() % sampleSize;
      while(alreadyTaken[r])
	{
	  r = rand() % sampleSize;
	}
      newSequence[i] = r;
      alreadyTaken[r] = 1;
    }

  double shuffledExpression[sampleSize][numberOfGenes];

  for (int i = 0; i < sampleSize; i++)
    {
      for (int j = 0; j < numberOfGenes; j++)
	{
	  shuffledExpression[i][j] = expressionValues[newSequence[i]][j];
	}
    }
  for (int i = 0; i < sampleSize; i++)
    {
      for (int j = 0; j < numberOfGenes; j++)
	{
	  expressionValues[i][j] = shuffledExpression[i][j];
	}
    }
  


  double corrcoefs[numberOfGenes][numberOfGenes];


  int selectionSequence[sampleSize][sampleSize];
  int table[sampleSize][sampleSize];

  for (int i = 0; i < sampleSize; i++)
    {
      for (int j = 0; j < sampleSize; j++)
	{
	  table[i][j] = 2;
	  selectionSequence[i][j] = 0;
	}
    }
  for (int i = 0; i < sampleSize; i++)
    {
      table[i][i] = 0;
    }
  
  int selection = 0;
  bool end = 0;

  int itr = 0; 
  int jtr = 0;
  int row;
  int col; 
  int nestlevel = 0;
  int fullSelections[numberOfGenes][numberOfGenes];
  for (int i = 0; i < numberOfGenes; i++)
    {
      for (int j = 0; j < numberOfGenes; j++)
	{
	  fullSelections[i][j] = 0;
	}
    }
  
  int prevNodesID[subSampleSize];
  int nodeCounter;
  bool nodeIsOK;
  

  double corrcoefAverage[numberOfGenes][numberOfGenes];
  double corrcoefAverageFull[numberOfGenes][numberOfGenes];



  double corrcoefAverageSquare[numberOfGenes][numberOfGenes];
  double corrcoefVar[numberOfGenes][numberOfGenes];
  double squareDev[numberOfGenes][numberOfGenes];
  double sumSquareDev[numberOfGenes][numberOfGenes];
  double meanSquareDev[numberOfGenes][numberOfGenes];
  

  for (int i = 0; i < numberOfGenes; i++)
    {
      for (int j = 0; j < numberOfGenes; j++)
	{
	  corrcoefAverage[i][j] = 0;
	  corrcoefAverageSquare[i][j] = 0;
	  corrcoefVar[i][j] = 0;
	  sumSquareDev[i][j] = 0;
	  list<Pair> completeSample;
	  for (int k = 0; k < sampleSize; k++)
	    {
	      Pair newPair;
	      newPair.x = expressionValues[k][i];
	      newPair.y = expressionValues[k][j];
	      completeSample.push_back(newPair);
	    }
	  corrcoefAverageFull[i][j] = spearman(completeSample);
	}
    }


  outStream.open(outFile);
  for (int i = 0; i < numberOfGenes; i++)
    {
      for (int j = 0; j < numberOfGenes; j++)
	{
	  outStream << geneName[i] << "\t" << geneName[j] << "\t" << corrcoefAverageFull[i][j] << "\n";
	}
    }

 outStream.close();

 outStream.open("RhoAndVarDetailed.txt");
 for (int i = 0; i < numberOfGenes; i++)
   {
     for (int j = 0; j < numberOfGenes; j++)
       {
	 outStream << geneName[i] << "\t" << geneName[j] << "\t" << corrcoefAverageFull[i][j] << "\t" << corrcoefAverage[i][j] << "\t" << fullSelections[i][j] << "\n";
       }
   }

 outStream.close();
   
}


double spearman(list<Pair> subSample)
{
  double xRank;
  double yRank;
  double xTies;
  double yTies;

  for (list<Pair>::iterator it = subSample.begin(); it != subSample.end(); it++)
    {
      xRank = 1;
      yRank = 1;
      xTies = -1;
      yTies = -1;
      for (list<Pair>::iterator jt = subSample.begin(); jt != subSample.end(); jt++)
	{
	  if ((*jt).x < (*it).x)
	    {
	      xRank++;
	    }
	  if ((*jt).x == (*it).x)
	    {
	      xTies++;
	    }

	  if ((*jt).y < (*it).y)
	    {
	      yRank++;
	    }
	  if ((*jt).y == (*it).y)
	    {
	      yTies++;
	    } 
	}
      
      (*it).xRank = xRank + xTies/2;
      (*it).yRank = yRank + yTies/2;
    }
  double averageRank = (double(subSample.size())+1)/2;
  double spearmanNum = 0;
  double spearmanDen1 = 0;
  double spearmanDen2 = 0;

  for (list<Pair>::iterator it= subSample.begin(); it != subSample.end(); it++)
    {
      spearmanNum = spearmanNum + ((*it).xRank-averageRank)*((*it).yRank-averageRank);
      spearmanDen1 = spearmanDen1 + pow(((*it).xRank-averageRank),2);
      spearmanDen2 = spearmanDen2 + pow(((*it).yRank-averageRank),2);
    }
  spearmanDen1 = sqrt(spearmanDen1);
  spearmanDen2 = sqrt(spearmanDen2);
  
  return spearmanNum/(spearmanDen1*spearmanDen2);
  
  
  
}

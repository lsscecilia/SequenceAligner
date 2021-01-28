#include <iostream>
#include <cstring>
#include <vector>

#include "alignment.h"

using namespace std;

enum edges {Up, Left, Diag, None}; 

struct matrix{
	int score; 
	edges edge; 
}; 


std::string CompressCigar(std::string uCigar){
	char previousChar=uCigar[uCigar.length()-1];
	int count=1; 
	std::string cigar="";
	for (int i=uCigar.length()-2; i>=0; i--){
		if (previousChar == uCigar[i]){
			count++;
		}
		else{
			cigar += std::to_string(count) + previousChar; 
			count = 1;
			previousChar = uCigar[i]; 
		}
	}
	cigar += std::to_string(count) + previousChar; 
	return cigar; 
};


int SemiGlobalAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin){		
		//initialise score table 
		vector<vector<matrix>> scoreTable( target_len+1 , vector<matrix> (query_len+1));   
		
		//initialise values with gap 
        for (int i=1; i<=target_len;i++){
            scoreTable[i][0]={0, None}; //initialise values
        } 
        for (int i=1;i<=query_len; i++){
            scoreTable[0][i]={gap*i, Left};
        }
		scoreTable[0][0]={0, None};
		
        //int diag, up, left, matchValue, max; //0,1,2
        
        int gMatch, gInsertion, gDeletion, matchValue, max, maxOverall=0, maxRow; 
        
        int temp; 

        for (int i=1; i<=target_len; i++){
            for (int r=1; r<=query_len; r++){
                //match value 
                if (target[i-1]==query[r-1]){
                    matchValue = match;
                }else{
                    matchValue = mismatch; 
                }

                gMatch = scoreTable[i-1][r-1].score + matchValue; 
                gDeletion = scoreTable[i-1][r].score + gap; 
                gInsertion = scoreTable[i][r-1].score+ gap; 
                
                //max
                temp = std::max(gMatch,gDeletion); 
                max = std::max(temp, gInsertion);
            
                if (max==gMatch){
					scoreTable[i][r] = {max, Diag}; 
				}
				else if (max==gInsertion) {
					scoreTable[i][r]= {max,Left}; 
				}
				else if (max==gDeletion){
					scoreTable[i][r]= {max, Up}; 
				}

				if (r==query_len){
					if (max > maxOverall){
						maxOverall = max; 
						maxRow = i; 
					}
                
				}
			}
		}
		
        //traceback + cigar 
        int row=maxRow, column=query_len; 
        std::string uCigar = ""; //uncompresssed cigar
        
        while(column!=0){
			if  (scoreTable[row][column].edge == Diag){
				if (scoreTable[row][column].score-match == scoreTable[row-1][column-1].score){
					uCigar += "M";
				}
				else{
					uCigar += "X";
					
				}
				row--;
				column--; 
			}
			else if (scoreTable[row][column].edge == Left){
				column--; 
				uCigar += "I"; 
			}
			else if (scoreTable[row][column].edge == Up){
				row--; 
				uCigar += "D"; 
			}
		}
		std::string tempCigar ; 
		if (cigar != nullptr){
			tempCigar = CompressCigar(uCigar); 
			*cigar= tempCigar.c_str();
		}
		
		
		if (target_begin != nullptr){
			unsigned int counter=1; 
			
			while (tempCigar[counter]=='D'){
				counter+=2;
			}
			*target_begin = counter;
		}
		
		return scoreTable[maxRow][query_len].score; 
};

int SmithWatermanAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin){

		//initialise score table 
		vector<vector<matrix>> scoreTable( target_len+1 , vector<matrix> (query_len+1));    
	
		//initialise values with gap 
		if (gap<0){
			for (int i=1; i<=target_len;i++){
				scoreTable[i][0]={0, None}; //initialise values
			} 
			for (int i=1;i<=query_len; i++){
				scoreTable[0][i]={0, None};
			}
		}
		else{
			for (int i=1; i<=target_len;i++){
				scoreTable[i][0]={gap*i, Up}; //initialise values
			} 
			for (int i=1;i<=query_len; i++){
				scoreTable[0][i]={gap*i, Left};
			}
		}
		scoreTable[0][0]={0, None};
        
        int gMatch, gInsertion, gDeletion, matchValue, max, temp; 
		int maxOverall=0, maxCol, maxRow; 

        for (int i=1; i<=target_len; i++){
            for (int r=1; r<=query_len; r++){
                //match value 
                if (target[i-1]==query[r-1]){
                    matchValue = match;
                }else{
                    matchValue = mismatch; 
                }

                gMatch = scoreTable[i-1][r-1].score + matchValue; 
                gDeletion = scoreTable[i-1][r].score + gap; 
                gInsertion = scoreTable[i][r-1].score+ gap; 
                
                //max
                temp = std::max(gMatch,gDeletion); 
                max = std::max(temp, gInsertion);
                
                if (max <= 0){
					scoreTable[i][r] = {0, None};
				}
                else if (max==gMatch){
					scoreTable[i][r] = {max, Diag}; 
				}
				else if (max==gInsertion) {
					scoreTable[i][r]= {max,Left}; 
				}
				else if (max==gDeletion){
					scoreTable[i][r]= {max, Up}; 
				}
				
				if (max > maxOverall) {
					maxOverall = max; 
					maxRow = i; 
					maxCol = r; 
				}
            }
        }
		
        //traceback + cigar 
        int row=maxRow, column=maxCol;  
        std::string uCigar = ""; //uncompresssed cigar
        
        
        while(scoreTable[row][column].edge!=None){
			if  (scoreTable[row][column].edge == Diag){
				if (scoreTable[row][column].score-match == scoreTable[row-1][column-1].score){
					uCigar += "M";
				}
				else{
					uCigar += "X";
					
				}
				row--;
				column--; 
			}
			else if (scoreTable[row][column].edge == Left){
				column--; 
				uCigar += "I";  
			}
			else if (scoreTable[row][column].edge == Up){
				row--; 
				uCigar += "D"; 
			}
		}
		std::string tempCigar ; 
		if (cigar != nullptr){
			tempCigar = CompressCigar(uCigar); 
			*cigar= tempCigar.c_str();
		}
		
		
		if (target_begin != nullptr){
			*target_begin =  (unsigned int) row; 
		}
		
		return maxOverall; 
};


int NeedlemanWunschAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin){
		
		//initialise score table 
		vector<vector<matrix>> scoreTable( target_len+1 , vector<matrix> (query_len+1)); 

		//initialise values with gap 
        for (int i=1; i<=target_len;i++){
            scoreTable[i][0]={gap*i, Up}; //initialise values
        } 
        for (int i=1;i<=query_len; i++){
            scoreTable[0][i]={gap*i, Left};
        }
		scoreTable[0][0]={0, None};
		
        //int diag, up, left, matchValue, max; //0,1,2
        
        int gMatch, gInsertion, gDeletion, matchValue, max; 
        
        int temp; 

		
        for (int i=1; i<=target_len; i++){
            for (int r=1; r<=query_len; r++){
                //match value 
                if (target[i-1]==query[r-1]){
                    matchValue = match;
                }else{
                    matchValue = mismatch; 
                }

                gMatch = scoreTable[i-1][r-1].score + matchValue; 
                gDeletion = scoreTable[i-1][r].score + gap; 
                gInsertion = scoreTable[i][r-1].score+ gap; 
                
                //max
                temp = std::max(gMatch,gDeletion); 
                max = std::max(temp, gInsertion);
                
                if (max==gMatch){
					scoreTable[i][r] = {max, Diag}; 
				}
				else if (max==gInsertion) {
					scoreTable[i][r]= {max,Left}; 
				}
				else if (max==gDeletion){
					scoreTable[i][r]= {max, Up}; 
				}

				
                
            }
        }
		
        //traceback + cigar 
        int row=target_len, column=query_len; 
        std::string uCigar = ""; //uncompresssed cigar
        
        
        while((row!=0)||(column!=0)){
			if  (scoreTable[row][column].edge == Diag){
				if (scoreTable[row][column].score-match == scoreTable[row-1][column-1].score){
					uCigar += "M";
				}
				else{
					uCigar += "X";
				}
				row--;
				column--; 
			}
			else if (scoreTable[row][column].edge == Left){
				column--; 
				uCigar += "I";
			}
			else if (scoreTable[row][column].edge == Up){
				row--; 
				uCigar += "D";
			}
		}
		std::string tempCigar ; 
		if (cigar != nullptr){
			tempCigar = CompressCigar(uCigar); 
			*cigar= tempCigar.c_str();
		}
		
		
		if (target_begin != nullptr){
			unsigned int counter=1; 
			while (tempCigar[counter]=='D'){	
				counter+=2;
			}
			*target_begin = counter;
		}
		
		return scoreTable[target_len][query_len].score; 
};

int Align(
	const char* query, unsigned int query_len,
	const char* target, unsigned int target_len,
	AlignmentType type,
	int match,
	int mismatch,
	int gap,
	std::string* cigar,
	unsigned int* target_begin){
	if (type == Global){
		return NeedlemanWunschAlgorithm(query,query_len,target,target_len,match,mismatch,gap,cigar,target_begin); 
	}
	else if (type == Local){
		return SmithWatermanAlgorithm(query,query_len,target,target_len,match,mismatch,gap,cigar,target_begin); 
	}
	else if (type == Semiglobal){
		return SemiGlobalAlgorithm(query,query_len,target,target_len,match,mismatch,gap,cigar,target_begin); 
	}
	else {
		return 9999; 
	}
}; 



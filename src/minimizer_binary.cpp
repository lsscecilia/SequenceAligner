#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>
#include <set>

#include "minimizer.h"
#include "utility.h"

using namespace std;

int mapNct(char n, bool seq){
	if (seq){
		switch (n){
			case 'a':
			case 'A':
				return 1; 
			case 'c':
			case 'C':
				return 0; 
			case 'g':
			case 'G':
				return 3; 
			case 't':
			case 'T':
				return 2;
			default:
				return -1; 
		}
	}
	else{
		switch (n){
			case 'a':
			case 'A':
				return 2; 
			case 'c':
			case 'C':
				return 3; 
			case 'g':
			case 'G':
				return 0; 
			case 't':
			case 'T':
				return 1;
			default:
				return -1; 
		}
	}				
}


unsigned int initFirstKmer(std::string kmer, int kmer_len, bool seq){
	unsigned int kmerValue=0;
	for (int i=0; i<kmer_len;i++){
		kmerValue = kmerValue << 2; 
		kmerValue = kmerValue | mapNct(kmer[i], seq); 
	}
	return kmerValue;  
}

unsigned int getKmer(unsigned int prevKmerValue, unsigned int mask, char nextN, bool seq){
	unsigned int kmerValue; 
	kmerValue = prevKmerValue << 2; 
	kmerValue = kmerValue | mapNct(nextN, seq); 
	kmerValue = kmerValue & mask;
	return kmerValue; 
}

//mask aft shift and or 
unsigned int getMask(int kmer_len){
	unsigned int kmerValue = 3; 
	for (int i=1; i<kmer_len; i++){
		kmerValue = kmerValue << 2; 
		kmerValue = kmerValue | 3;
	}
	return kmerValue; 
}

//seq == 0 --> reverse strand
std::vector<std::tuple<unsigned int, unsigned int, bool>> getAllKmer(const char* sequence, int sequence_len,int kmer_len, bool seq){
	std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList;
	unsigned int kmerValue, mask;  
	
	std::string stringSeq(sequence); 
	kmerValue = initFirstKmer(stringSeq.substr(0,kmer_len),kmer_len,seq); 
	mask = getMask(kmer_len); 
	//put into vector
	kmerList.emplace_back(make_tuple((unsigned int) kmerValue, 0, seq)); 
	
	for (int i=1; i<sequence_len-kmer_len+1; i++){
		kmerValue = getKmer(kmerValue,mask, stringSeq[i+kmer_len-1],seq); 
		//put into vector
		kmerList.emplace_back(make_tuple((unsigned int) kmerValue, i, seq)); 
	}
	return kmerList; 
}

void initFindMinKmer(const std::vector<std::tuple<unsigned int, unsigned int, bool>>& kmerList, 
	int window_len,std::tuple<unsigned int, unsigned int, bool>* minKmer, int start){
	int kmerValue;
	*minKmer = kmerList[start];
	for (int i=start+1; i<window_len+start; i++){
		//find min kmer
		kmerValue = get<0>(kmerList[i]);
		if (kmerValue < get<0>(*minKmer)){
			*minKmer = kmerList[i]; 
		}
	}
}
	
void findMinKmer(const std::vector<std::tuple<unsigned int, unsigned int, bool>>& kmerList,
std::tuple<unsigned int, unsigned int, bool>& nextKmer, 
int window_len,int kmer_len, std::tuple<unsigned int, unsigned int, bool>* minKmer, 
std::tuple<unsigned int, unsigned int, bool>& prevMinKmer){
	int kmerValue = get<0>(nextKmer);
	int kmerIndex = get<1>(nextKmer); 
	int substringLen = window_len + kmer_len - 1; 

	if (kmerIndex-kmer_len-1== get<1>(prevMinKmer)){
		//refind min 
		initFindMinKmer(kmerList,window_len, minKmer ,get<1>(prevMinKmer)+1);
		
	}
	else {
		if (kmerValue < get<0>(prevMinKmer)){
			*minKmer = nextKmer;
		}
		else{
			*minKmer = prevMinKmer;
		}
	}
}

				
std::vector<std::tuple<unsigned int, unsigned int, bool>> MinimizeBinary(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len){

	std::vector<std::tuple<unsigned int, unsigned int, bool>> allKmer, rAllKmer;
	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers;
	std::tuple<unsigned int, unsigned int, bool> minKmer, rMinKmer;
	//get all Kmer
	allKmer = getAllKmer(sequence, sequence_len, kmer_len, true); 
	rAllKmer = getAllKmer(sequence, sequence_len, kmer_len, false);

	//init 
	initFindMinKmer(allKmer,window_len, &minKmer, 0);
	initFindMinKmer(rAllKmer,window_len, &rMinKmer, 0);

	//put into vector 
	if (get<0>(minKmer) < get<0>(rMinKmer)){
		minimizers.emplace_back(minKmer); 
	}
	else{
		minimizers.emplace_back(rMinKmer); 
	}

	//find all other min kmer
	for (int i=1; i<sequence_len-window_len-kmer_len+2;i++){
		//cout << i << "finding minimizer.." << endl; 
		findMinKmer(allKmer,allKmer[i+window_len-1],window_len,kmer_len, &minKmer,minKmer);
		findMinKmer(rAllKmer,rAllKmer[i+window_len-1],window_len,kmer_len, &rMinKmer,rMinKmer);

		//put into vector 
		if (get<0>(minKmer) < get<0>(rMinKmer)){
			minimizers.emplace_back(minKmer); 
		}
		else{
			minimizers.emplace_back(rMinKmer); 
		}
		
	}
	
	//sort and remove duplicate in kmerList
	return removeDuplicate(minimizers); 	
}

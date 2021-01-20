#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>

#include <minimizer.h>

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
	unsigned int kmerValue;
	for (int i=0; i<kmer_len;i++){
		kmerValue = kmerValue << 2; 
		kmerValue = kmerValue | mapNct(kmer[i], seq); 
	}
	return kmerValue;
}

unsigned int getKmer(unsigned int prevKmerValue, unsigned int mask, char nextN, bool seq){
	unsigned int kmerValue; 
	kmerValue = prevKmerValue & mask;
	kmerValue = kmerValue << 2; 
	kmerValue = kmerValue | mapNct(kmer[i], seq); 
}

unsigned int getMask(int kmer_len){
	unsigned int kmerValue = 3; 
	for (int i=1; i<kmer_len-1; i++){
		kmerValue = kmerValue << 2 + 3;
	}
	return kmerValue; 
}

//seq == 0 --> reverse strand
std::vector<std::tuple<unsigned int, unsigned int, bool>> getAllKmer(const char* sequence, int sequence_len,int kmer_len, bool seq){
	//store in tuple?
	std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList;
	unsigned int kmerValue, mask;  
	
	std::string stringSeq(seqeunce); 
	kmerValue = initFirstKmer(stringSeq.substring(0,kmer_len),seq); 
	mask = getMask(kmer_len); 
	//put into vector
	kmerList.push_back(make_tuple((unsigned int) kmerValue, 0, seq)); 
	
	for (int i=1; i<sequence_len-kmer+1; i++){
		kmerValue = getKmer(kmerValue,mask, stringSeq[i],seq); 
		//put into vector
		kmerList.push_back(make_tuple((unsigned int) kmerValue, i, seq)); 
	}
	return kmerList; 
}

void initFindMinKmer(std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList, 
	int window_len,std::tuple<unsigned int, unsigned int, bool>* minKmer, int start){
	int kmerValue;
	*min = get<0>(kmerList[start]); 
	*minIndex = get<1>(kmerList[start]);
	for (int i=start+1; i<window_len+start; i++){
		//find min kmer
		kmerValue = get<0>(kmerList[i]);
		if (kmerValue < *min){
			*minKmer = kmerList[i]; 
		}
	}
}
	
void findMinKmer(std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList,
std::tuple<unsigned int, unsigned int, bool> nextKmer, 
int window_len,std::tuple<unsigned int, unsigned int, bool>* minKmer, 
std::tuple<unsigned int, unsigned int, bool> prevMinKmer){
	int kmerValue = get<0>(nextKmer);
	int kmerIndex = get<1>(nextKmer); 
	int substringLen = window_len + kmer_len - 1; 
	if (kmerIndex-substringLen == get<1>(prevMinKmer)){
		//refind min 
		initFindMinKmer(kmerList,window_len, minKmer ,get<1>(prevMinKmer)+1);
			
	}
	else {
		if (kmerValue < preMin){
			*minKmer = nextKmer;
		}
		else{
			*minKmer = prevMinKmer;
		}
	}
}

std::vector<std::tuple<unsigned int, unsigned int, bool>> removeDuplicate(
	std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList){
		std::vector<std::tuple<unsigned int, unsigned int, bool>> newKmerList; 
		std::sort(kmerList.begin(), kmerList.end()); 
		unsigned int prevSeq=get<0>(kmerList[0]), prevIndex=get<1>(kmerList[0]); 
		bool strand=get<2>(kmerList[0]); 
		newKmerList.push_back(kmerList[0]);
		
		for (int i=1; i <kmerList.size(); i++){
			if (!(get<0>(kmerList[i])==prevSeq && get<1>(kmerList[i])==prevIndex 
				&& get<2>(kmerList[i])==strand)){
					newKmerList.push_back(kmerList[i]); 
					prevSeq= get<0>(kmerList[i]); 
					prevIndex = get<1>(kmerList[i]);
					strand = get<2>(kmerList[i]); 
			}
		}
		return newKmerList; 
	}
				

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len){

	std::vector<std::tuple<unsigned int, unsigned int, bool>> allKmer, rAllKmer, minimizers;
	std::tuple<unsigned int, unsigned int, bool> minKmer, rMinKmer;
	//get all Kmer
	allKmer = getAllKmer(sequence, sequence_len, kmer_len, true); 
	rAllKmer = getAllKmer(sequence, sequence_len, kmer_len, false);

	
	//init 
	initFindMinKmer(allKmer,window_len, &minKmer, 0);
	initFindMinKmer(rAllKmer,window_len, &rMinKmer, 0);
	
	//put into vector 
	if (get<0>(minKmer) < get<0>(rMinKmer)){
		minimizers.push_back(minKmer); 
	}
	else{
		minimizers.push_back(rMinKmer); 
	}
	
	
	//find all other min kmer
	for (int i=1; i<sequence_len-kmer_len+1;i++){
		findMinKmer(allKmer,allKmer[i+window_len-1],window_len, &minKmer,minKmer);
		findMinKmer(rAllKmer,rAllKmer[i+window_len-1],window_len, &rRinKmer,rMinKmer);
		
			//put into vector 
		if (get<0>(minKmer) < get<0>(rMinKmer)){
			minimizers.push_back(minKmer); 
		}
		else{
			minimizers.push_back(rMinKmer); 
		}
	}
	
	//sort and remove duplicate in kmerList
	return removeDuplicate(minimizers); 	
}

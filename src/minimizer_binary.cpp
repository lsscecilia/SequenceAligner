#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>

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

//tested, working fine
unsigned int initFirstKmer(std::string kmer, int kmer_len, bool seq){
	unsigned int kmerValue=0;
	for (int i=0; i<kmer_len;i++){
		kmerValue = kmerValue << 2; 
		kmerValue = kmerValue | mapNct(kmer[i], seq); 
		//cout << "counter: " << i << " kmer value: " << kmerValue << endl; 
	}
	return kmerValue;  
}

//tested
unsigned int getKmer(unsigned int prevKmerValue, unsigned int mask, char nextN, bool seq){
	unsigned int kmerValue; 
	
	//cout << "prev kmer value: " << prevKmerValue << endl; 
	//cout << "mask: " << mask << endl; 
	
	kmerValue = prevKmerValue << 2; 
	//cout << " shift: "  << kmerValue << endl; 
	
	kmerValue = kmerValue | mapNct(nextN, seq); 
	//cout << "next char: " << nextN << endl; 
	//cout << "or: " << kmerValue << endl; 
	
	kmerValue = kmerValue & mask;
	//cout << "after masking: " << kmerValue << endl; 
	return kmerValue; 
}

//mask aft shift and or , tested
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
	//store in tuple?
	std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList;
	unsigned int kmerValue, mask;  
	
	std::string stringSeq(sequence); 
	kmerValue = initFirstKmer(stringSeq.substr(0,kmer_len),kmer_len,seq); 
	mask = getMask(kmer_len); 
	//put into vector
	kmerList.push_back(make_tuple((unsigned int) kmerValue, 0, seq)); 
	
	for (int i=1; i<sequence_len-kmer_len+1; i++){
		kmerValue = getKmer(kmerValue,mask, stringSeq[i+kmer_len-1],seq); 
		//put into vector
		kmerList.push_back(make_tuple((unsigned int) kmerValue, i, seq)); 
	}
	return kmerList; 
}

void initFindMinKmer(std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList, 
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
	
void findMinKmer(std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList,
std::tuple<unsigned int, unsigned int, bool> nextKmer, 
int window_len,int kmer_len, std::tuple<unsigned int, unsigned int, bool>* minKmer, 
std::tuple<unsigned int, unsigned int, bool> prevMinKmer){
	int kmerValue = get<0>(nextKmer);
	int kmerIndex = get<1>(nextKmer); 
	int substringLen = window_len + kmer_len - 1; 
	//cout << "prev min index: " << get<1>(prevMinKmer) << endl;
	//cout << " start: " << kmerIndex-kmer_len-1 << endl; 
	if (kmerIndex-kmer_len-1== get<1>(prevMinKmer)){
		//refind min 
		//cout << "init ... find min kmer" << endl;
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

	std::vector<std::tuple<unsigned int, unsigned int, bool>> allKmer, rAllKmer, minimizers;
	std::tuple<unsigned int, unsigned int, bool> minKmer, rMinKmer;
	//get all Kmer
	allKmer = getAllKmer(sequence, sequence_len, kmer_len, true); 
	rAllKmer = getAllKmer(sequence, sequence_len, kmer_len, false);

	/*
	//print kmer
	for (int i=0; i< allKmer.size();i++){
		cout <<"value: " << get<0>(allKmer[i]) << "position :" 
		<< get<1>(allKmer[i]) << "strand: " << get<2>(allKmer[i]) << endl; 
	}
	
		//print kmer
	for (int i=0; i< rAllKmer.size();i++){
		cout <<"value: " << get<0>(rAllKmer[i]) << "position :" 
		<< get<1>(rAllKmer[i]) << "strand: " << get<2>(rAllKmer[i]) << endl; 
	}*/
	
	
	
	//init 
	initFindMinKmer(allKmer,window_len, &minKmer, 0);
	initFindMinKmer(rAllKmer,window_len, &rMinKmer, 0);
	//cout << "first min kmer: " << get<0>(minKmer) << " p: "<<get<1>(minKmer)<< endl; 
	//cout << "first rMin kmer: " << get<0>(rMinKmer) << " p: "<<get<1>(rMinKmer)<< endl; 
	
	//put into vector 
	if (get<0>(minKmer) < get<0>(rMinKmer)){
		minimizers.push_back(minKmer); 
	}
	else{
		minimizers.push_back(rMinKmer); 
	}

	/*
	//print vector
	cout << "---In vector---" << endl; 
	for (int i =0; i<minimizers.size();i++){
		cout << "value: " << get<0>(minimizers[i]) << " p: " <<get<1>(minimizers[i])
		<< " strand : " << get<2>(minimizers[i]) << endl;
	}
	cout << "---------------" << endl; 	*/

	//find all other min kmer
	for (int i=1; i<sequence_len-window_len-1;i++){
		findMinKmer(allKmer,allKmer[i+window_len-1],window_len,kmer_len, &minKmer,minKmer);
		findMinKmer(rAllKmer,rAllKmer[i+window_len-1],window_len,kmer_len, &rMinKmer,rMinKmer);
		
		/*
		cout << i <<" min kmer: " << get<0>(minKmer) << " p: "<<get<1>(minKmer)<< endl; 
		cout << i <<" rMin kmer: " << get<0>(rMinKmer) << " p: "<<get<1>(rMinKmer)<< endl; 
		*/
		//put into vector 
		if (get<0>(minKmer) < get<0>(rMinKmer)){
			minimizers.push_back(minKmer); 
		}
		else{
			minimizers.push_back(rMinKmer); 
		}
		
		/*
		//print vector
		cout << "---In vector---" << endl; 
		for (int i =0; i<minimizers.size();i++){
			cout << "value: " << get<0>(minimizers[i]) << " p: " <<get<1>(minimizers[i])
			<< " strand : " << get<2>(minimizers[i]) << endl;
		}
		cout << "---------------" << endl; 	*/
	}
	
	//sort and remove duplicate in kmerList
	return removeDuplicate(minimizers); 	
}

/*
int main(){ 
	std::string seq = "TGACGTACATGGACA"; 
	unsigned int len = 15; 
	unsigned int kmer_len = 3;
	unsigned int w = 4; 
	std::vector<std::tuple<unsigned int, unsigned int, bool>> result = 
		Minimize(seq.c_str(), len,kmer_len,w); 
	
	for (int i=0; i < result.size(); i++){
		cout << get<0>(result[i]) << " " ; 
		cout << get<1>(result[i]) << " " ; 
		cout << get<2>(result[i]) << " "  << endl; 
	}
} */

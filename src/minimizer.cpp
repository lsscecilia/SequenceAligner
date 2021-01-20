#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>

#include <minimizer.h>


using namespace std; 

std::string map(const char* sequence, int seq_len){
	string mapSequence=""; 
	//0,1,2,3 CATG --> odd number bases of kmer (reverse for even base of kmer)
	for (int i=0; i<seq_len; i++){
		switch (sequence[i]){
			case 'a':
			case 'A':
				mapSequence += "1"; 
				break; 
			case 'c':
			case 'C':
				mapSequence += "0"; 
				break; 
			case 'g':
			case 'G':
				//1
				mapSequence += "3";
				break; 
			case 't':
			case 'T':
				//1
				mapSequence += "2";
				break; 
		}
	}
	//cout << mapSequence << endl; 
	return mapSequence; 
}

std::string reverseMap(const char* sequence, int seq_len){
	string rMapSequence=""; 
	//0,1,2,3 CATG --> odd number bases of kmer (reverse for even base of kmer)
	for (int i=0; i<seq_len; i++){
		switch (sequence[i]){
			case 'a':
			case 'A':
				rMapSequence += "2"; 
				break; 
			case 'c':
			case 'C':
				rMapSequence += "3"; 
				break; 
			case 'g':
			case 'G':
				//1
				rMapSequence += "0";
				break; 
			case 't':
			case 'T':
				//1
				rMapSequence += "1";
				break; 
		}
	}
	//cout << rMapSequence << endl;
	return rMapSequence; 
}
			
void initWindow(std::string substring, int substringLen, int kmer_len, 
	int* min, int* minIndex){ 
	int value; 
	*min = stoi(substring.substr(0,kmer_len)); 
	//cout << "init min ...." << *min << endl;
	for (int i=1; i<substringLen-kmer_len+1;i++){
		value = stoi(substring.substr(i,kmer_len));
		//cout << "init value ... " << value << endl; 
		if (*min > value){
			*min = value; 
			*minIndex = i; 
		}
	}
}

void compareLastKmer(std::string substring, int substringLen,int kmer_len, int preMin, 
	int preMinIndex, int* min, int* minIndex, int i){
	*min=444; 

			//cout << substring << endl; 
	int value = stoi(substring.substr(substringLen-kmer_len,kmer_len));
		//cout << "flag" << endl; 
		//cout << "value" << value << endl; 
	if (value <  preMin) {
		*min = value; 
		*minIndex = i+substringLen-kmer_len; 
		//put into tuple
		
	}
	else{
		*min = preMin;
		*minIndex = preMinIndex; 
	}
		//cout  << "min" << *min << endl; 
}

void compareKmer(std::string substring, int substringLen, int kmer_len, 
	int *min, int* minIndex, int i ){
	int value; 
	*min = 444; 
	//cout << substring << endl; 
	for (int r=0; r<1+substringLen-kmer_len; r++){
			//cout << "flag1" << endl; 
		value = stoi(substring.substr(r,kmer_len));
		//cout << "value" << value << endl; 
		if (*min > value){
			*min = value; 
			*minIndex = i+r; 
		}
		//cout  << "min" << *min << endl; 
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

	std::vector<std::tuple<unsigned int, unsigned int, bool>> kmerList;
	tuple<unsigned int, unsigned int, bool> kmer; 
	
	int minIndex,min,preMin,preMinIndex; 
	int rMinIndex, rMin, rPreMin, rPreMinIndex; 
	std::string substring, rSubstring; 
	
	//length of string = w+k-1
	int substringLen = window_len + kmer_len - 1;  
	
	//map nucleotides to number
	std::string mapSeq = map(sequence,sequence_len); 
	std::string rMapSeq = reverseMap(sequence, sequence_len); 
	
	substring = mapSeq.substr(0,substringLen);  
	initWindow(substring, substringLen,kmer_len,&preMin, &preMinIndex); 
	//cout << "init min" << preMin << endl; 
	
	rSubstring = rMapSeq.substr(0,substringLen);  
	initWindow(rSubstring, substringLen,kmer_len,&rPreMin, &rPreMinIndex); 
	//cout << "init min" << rPreMin << endl; 
	
	//kmerList.push_back(make_tuple((unsigned int) preMin, preMinIndex, true)); 
	
	if (preMin < rPreMin){
		//preMin into vetor
		kmerList.push_back(make_tuple((unsigned int) preMin, preMinIndex, true)); 
	}
	else{
		//rPreMin into vector
		kmerList.push_back(make_tuple((unsigned int) rPreMin, rPreMinIndex, false)); 
	}
	
	//cout << "all" << (sequence_len-substringLen+1 )<< endl; 
	for (int i=1; i<(sequence_len-substringLen+1); i++){
		//cout << "counter " << i <<endl; 
		substring = mapSeq.substr(i,substringLen); 
		rSubstring = rMapSeq.substr(i,substringLen);
		min = 4*kmer_len; 
		
		if (preMinIndex != i-1){
			compareLastKmer(substring, substringLen, kmer_len, preMin, 
				preMinIndex, &min, &minIndex,i); 
		}
		else{
			compareKmer(substring, substringLen, kmer_len, &min, &minIndex,i);
		}
		
				
		if (rPreMinIndex != i-1){
			compareLastKmer(rSubstring, substringLen, kmer_len, rPreMin, 
				rPreMinIndex, &rMin, &rMinIndex,i); 
		}
		else{
			compareKmer(rSubstring, substringLen, kmer_len, &rMin, &rMinIndex,i);
		}
		
		//kmerList.push_back(make_tuple((unsigned int) min, minIndex, true)); 
		
		
		if (min < rMin){
			//preMin into vetor
			kmerList.push_back(make_tuple((unsigned int) min, minIndex, true)); 
		}
		else{
			//rPreMin into vector
			kmerList.push_back(make_tuple((unsigned int) rMin, rMinIndex, false)); 
		}
	

		//save as previous
		preMinIndex = minIndex; 
		preMin = min; 
		
		rPreMinIndex = rMinIndex;
		rPreMin = rMin; 
		
	} 
	
	//sort and remove duplicate in kmerList
	return removeDuplicate(kmerList); 
	
	
	//return tuple
	//return kmerList; 
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

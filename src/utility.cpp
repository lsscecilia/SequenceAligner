#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std; 

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
				

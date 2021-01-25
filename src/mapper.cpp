#include <iostream>
#include <cstring>
#include <getopt.h>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

#include "config.h"

#include "alignment.h"
#include "minimizer_binary.h"

using namespace std;

static struct option long_options[] = {
  /* These options don’t set a flag.
	 We distinguish them by their indices. */
  {"version", no_argument,       0, 'v'},
  {"help",  no_argument,       0, 'h'},
  {"alignment_type",  required_argument, 0, 'a'},
  {"match",  required_argument, 0, 'm'},
  {"nomatch",    required_argument, 0, 'n'},
  {"gap",    required_argument, 0, 'g'},
  {0, 0, 0, 0}
};


struct Sequence {  // or any other name
	std::string name;
    std::string data;   
    std::string quality;
 public:
	 Sequence(const char* name, std::uint32_t name_len, const char* data,
           std::uint32_t data_len)
      : name(name, name_len), data(data, data_len) {}
 
	Sequence(  // required arguments
      const char* name, std::uint32_t name_len,
      const char* data, std::uint32_t data_len,
      const char* quality, std::uint32_t quality_len):
		name(name, name_len),
        data(data, data_len),
        quality(quality, quality_len) {}
}; 



void Help(){
	cout << "some help command..." << endl;
};

void ProjectVersion(){
	cout << "v" << PROJECT_VER << endl;
};

void PrintStats(const std::vector<std::unique_ptr<Sequence>>& fragments){
	int lengthAllFragments=0, minLength=0, maxLength=0, N50sum = 0, N50;
	std::string maxLengthFrag, minLengthFrag; 
	
	minLength = fragments[0]->data.length(); 
	for (auto i = 0; i != fragments.size(); i++) {
		lengthAllFragments +=  fragments[i]->data.length(); 
		if (fragments[i]->data.length() > maxLength){
			maxLength = fragments[i]->data.length(); 
			maxLengthFrag = fragments[i]->name; 
		}
		if (fragments[i]->data.length() < minLength) {
			minLength = fragments[i]->data.length(); 
			minLengthFrag = fragments[i]->name; 
		}
	}
	    
    for (int i = 0; i < fragments.size(); i++ ) {
        N50sum += fragments[i]->data.length();
        if (N50sum > 0.5 * lengthAllFragments) {
            N50 = fragments[i]->data.length();
            break;
        }
    }
    
    cout << "---------------Statistics---------------" << endl; 
	cout << "Number of sequences: " <<  fragments.size() << std::endl;
    cout << "Total length of all fragments: " << lengthAllFragments << std::endl;
    cout << "Largest fragment: " << maxLengthFrag << std::endl;
    cout << "  length: " << maxLength << std::endl;
    cout << "Smallest fragment: " << minLengthFrag << std::endl;
    cout << "   length: " << minLength << std::endl;
    cout << "Average length: " << ((double) lengthAllFragments) / fragments.size() << std::endl;
    cout << "N50 length: " << N50<< std::endl;
}; 

void PrintAlignmentResult(
		std::string query, unsigned int query_len,
		std::string target, unsigned int target_len,
		AlignmentType type,
		int match,
		int mismatch,
		int gap,
		std::string cigar,
		unsigned int target_begin,
		int alignmentScore){
			
			cout << "---------------Alignment---------------" << endl; 
			cout << "Query: " << query << endl; 
			cout << "Query len: " << query_len << endl; 
			cout << "Target: " << target << endl; 
			cout << "Target len: " << target_len << endl; 
			cout << "Alignment type: " << type << endl; 
			cout << "match: " << match << endl;
			cout << "mismatch: " << mismatch << endl; 
			cout << "gap: " << gap << endl << endl;
			cout << "---------------Results---------------" << endl; 
			cout << "alignment score: " << alignmentScore << endl;
			cout << "taget begin: " << target_begin << endl << endl; 
			//cout << "cigar: " << cigar << endl;
			
}; 
	
	
bool IsFastaFile(std::string file){
	//fasta file format extension
	//File extensions :     file.fa, file.fasta, file.fsa
	if (file.find(".fa") != std::string::npos){
		return true; 
	}
	else if (file.find(".fasta") != std::string::npos){
		return true; 
	}
	else if (file.find(".fsa") != std::string::npos){
		return true;
	}
	else {
		return false; 
	}
	
}

bool IsFastqFile(std::string file){
	//fastq file format extension 
	//File extensions: file.fastq, file.sanfastq, file.fq
	
	if (file.find(".fastq") != std::string::npos){
		
		return true; 
	}
	else if (file.find(".sanfastq") != std::string::npos){
		return true; 
	}
	else if (file.find(".fq") != std::string::npos){
		return true;
	}
	else {
		return false; 
	}
}

void getMinimizer(std::unique_ptr<Sequence> & seq, map<unsigned int,vector<tuple<unsigned int,bool>>> *minimizerIndex, int kmer_len, int window_len){
	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizer;
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	tuple<unsigned int, bool> temp; 
	vector<tuple<unsigned int, bool>> tempVector; 
	
	minimizer = MinimizeBinary(seq->data.c_str(),seq->data.length(),kmer_len,window_len);
	 

	for (int i=0; i<minimizer.size(); i++){
		tempVector.clear();
		it = minimizerIndex->find(get<0>(minimizer[i]));
		temp = make_tuple(get<1>(minimizer[i]),get<2>(minimizer[i]));
		if ( it == minimizerIndex->end() ) {
		// not found
			//cout << "not found " << endl; 
			tempVector.push_back(temp); 
			minimizerIndex->insert(pair<unsigned int,vector<tuple<unsigned int,bool>>>
			(get<0>(minimizer[i]),tempVector)); 

		} else {
		// found
			//cout << "found" << endl;
			(it->second).push_back(temp);
			//cout << it->second << endl; 
		}	
	} 
}

//wrong, parameter should be occurance
void printMinimizerIndexTable(map<unsigned int,vector<tuple<unsigned int,bool>>> minimizerIndex){
	/*
	std::map<unsigned int,tuple<unsigned int,bool>>::iterator it;
	cout << "\nMinimizer Count : \n"; 
	cout << "\tKEY\tELEMENT\n"; 
	for (it = minimizerIndex.begin(); it != minimizerIndex.end(); ++it) { 
		cout << '\t' << it->first << '\t' << get<>it->second << '\n'; 
	} 
	cout << endl;*/
}

vector<tuple<int, unsigned int>> getOccurrences(map<unsigned int,vector<tuple<unsigned int,bool>>> minimizerIndex){
	//count, minimizer
	vector<tuple<int, unsigned int>> occurrences; 
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	for (it = minimizerIndex.begin(); it != minimizerIndex.end(); ++it) { 
		occurrences.push_back(make_tuple((it->second).size(), it->first)); 
	}
	return occurrences; 
}

int getSingletonCount(vector<tuple<int, unsigned int>> occurrences){
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	int singletonCount=0;
	for (int i=0; i<occurrences.size();i++){
		if (get<0>(occurrences[i])==1){
			singletonCount++;
		}
	}
	return singletonCount; 
}

int getNumOccurrencesMostFrequentMinimizer(float f,vector<tuple<int, unsigned int>> occurrences){
	std::sort(occurrences.begin(),occurrences.end()); 
	int index = (int) occurrences.size()*f;
	return get<0>(occurrences[index]);
}

vector<tuple<unsigned int, unsigned int>> matchMinimizer(map<unsigned int,vector<tuple<unsigned int,bool>>> referenceIndex, 
map<unsigned int,vector<tuple<unsigned int,bool>>> fragmentIndex){
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator search;
	vector<tuple<unsigned int, unsigned int>> match; 
	vector<unsigned int> fragmentOrigin, fragmentReverse, refOrigin, refReverse; 

	for (it = fragmentIndex.begin(); it != fragmentIndex.end(); ++it){
		search = referenceIndex.find(it->first); 
		if ( search != referenceIndex.end() ){
			//found
			for (int i=0; i< it->second.size(); i++){
				if (get<1>(it->second[i])){
					//if origin
					fragmentOrigin.push_back(get<0>(it->second[i]));
				}
				else{
					fragmentReverse.push_back(get<0>(it->second[i])); 
				}				
			}

			for (int i=0; i<search->second.size();i++ ){
				//if origin
				if (get<1>(search->second[i])){
					refOrigin.push_back(get<0>(search->second[i]));
				}
				else{
					refReverse.push_back(get<0>(search->second[i])); 
				}
			}

			for (int i=0; i< fragmentOrigin.size(); i++){
				for (int r=0; r<refOrigin.size(); r++){
					match.push_back(make_tuple(fragmentOrigin[i], refOrigin[r])); 
				}
			}

			for (int i=0; i< fragmentReverse.size(); i++){
				for (int r=0; r<refReverse.size(); r++){
					match.push_back(make_tuple(fragmentReverse[i], refReverse[r])); 
				}
			}
		}
	}
	return match; 
}

/*
vector<tuple<unsigned int, unsigned int, vector<unsigned int>>> 
matchMinimizer(map<unsigned int,vector<tuple<unsigned int,bool>>> referenceIndex, 
map<unsigned int,vector<tuple<unsigned int,bool>>> fragmentIndex ){
	//match from fragment to genome
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	std::map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator search;
	vector<tuple<unsigned int, unsigned int, vector<unsigned int>>> matchTable; 
	//vector<tuple<unsigned int,bool>> fragmentPosition; 
	vector<unsigned int> fragmentOrigin, fragmentReverse, refOrigin, refReverse; 
	
	
	int counter =0; 
	for (it = referenceIndex.begin(); it != referenceIndex.end();++it){
		cout << it->second.size() << endl; 
		counter++; 
		if (counter >20)
			break;
	}

	cout << "flag" << endl; 
	//for each fragment
	for (it = fragmentIndex.begin(); it != fragmentIndex.end(); ++it) { 
		//cout << "in loop" << endl; 
		fragmentOrigin.clear(); 
		fragmentReverse.clear(); 
		search = referenceIndex.find(it->first); 
		//cout << "after search ... " << endl; 
		if ( search != referenceIndex.end() ) {
			//cout << "found .." << endl; 
			//cout << "fragment position size .. " << it->second.size() << endl; 
			//cout << "ref position size... " << search->second.size() << endl; 
			//found 
			//split according to strand
			for (int i=0; i< it->second.size(); i++){
				//cout << "in nested loop" << endl; 
				//split
				if (get<1>(it->second[i])){
					//if origin
					fragmentOrigin.push_back(get<0>(it->second[i]));
				}
				else{
					fragmentReverse.push_back(get<0>(it->second[i])); 
				}				
			}
			//cout << "after split" << endl; 
			//cout << "fragment origin size: " << fragmentOrigin.size() << endl; 
			//cout << "fragment reverse size: " << fragmentReverse.size() << endl; 
			if (fragmentOrigin.size()!=0 && fragmentReverse.size()!=0){
				//find reference on origin
				for (int i=0; i<search->second.size();i++ ){
					//split
					if (get<1>(search->second[i])){
						refOrigin.push_back(get<0>(search->second[i]));
					}
					else{
						refReverse.push_back(get<0>(search->second[i])); 
					}
				}
				std::sort(refReverse.begin(),refReverse.end()); 
				std::sort(refOrigin.begin(),refOrigin.end()); 
			}
			else if (fragmentOrigin.size()==0){
				//find reference on origin
				for (int i=0; i<search->second.size();i++ ){
					//split
					if (!get<1>(search->second[i])){
						refReverse.push_back(get<0>(search->second[i])); 
					}
				}
				std::sort(refReverse.begin(),refReverse.end()); 
			}
			else if (fragmentReverse.size()==0){
				for (int i=0; i<search->second.size();i++ ){
					//split
					//cout << "split ref" << endl; 
					if (get<1>(search->second[i])){
						//cout << "wth" << endl; 
						refOrigin.push_back(get<0>(search->second[i]));
						//cout << "wth is goin on here" << endl; 
					}
				}
				//cout << "sort problem" << endl; 
				std::sort(refOrigin.begin(),refOrigin.end()); 
			}

			//cout << "ref genome split" << endl;	

			//create for each minimizer with the difference position
			for (int i =0; i<fragmentOrigin.size();i++){
				
				matchTable.push_back(make_tuple(fragmentOrigin[i],it->first, refOrigin));
			}

			for (int i=0; i<fragmentReverse.size(); i++){
				matchTable.push_back(make_tuple(fragmentReverse[i],it->first,refReverse));
			}
		} 
	}
	return matchTable;
} */
// Binary search (note boundaries in the caller) 
int CeilIndex(vector<tuple<unsigned int,unsigned int>> &matches, vector<unsigned int> &T, int l, int r,
                 int key) {
  while (r - l > 1) {
    int m = l + (r - l) / 2;
    if (get<1>(matches[T[m]]) >= key)
      r = m;
    else
      l = m;
  }

  return r;
}

  
int LongestIncreasingSubsequence(vector<tuple<unsigned int, unsigned int>> &matches, int n,
                                 int &t_begin, int &t_end, int &q_begin,
                                 int &q_end) {
	if (n<=1)
		return n; 
	
	vector<unsigned int> tailIndices(n, 0);
	vector<unsigned int> prevIndices(n, -1);

	int len = 1;
	for (int i = 1; i < n; i++) {
	if (get<1>(matches[i])< get<1>(matches[tailIndices[0]])) {
		// new smallest value
		tailIndices[0] = i;
	} else if (get<1>(matches[i]) > get<1>(matches[tailIndices[len - 1]])) {
		// matches[i].second wants to extend largest subsequence
		prevIndices[i] = tailIndices[len - 1];
		tailIndices[len++] = i;
	} else {
		// matches[i].second wants to be a potential condidate of
		// future subsequence
		// It will replace ceil value in tailIndices
		int pos =
			CeilIndex(matches, tailIndices, -1, len - 1, get<1>(matches[i]));

		prevIndices[i] = tailIndices[pos - 1];
		tailIndices[pos] = i;
	}
	}
	t_begin = get<1>(matches[prevIndices[1]]);
	q_begin = get<0>(matches[prevIndices[1]]);
	t_end = get<1>(matches[tailIndices[len - 1]]);
	q_end = get<0>(matches[tailIndices[len - 1]]);

	return len;
}

int main (int argc, char **argv){
	/* getopt_long stores the option index here. */
    int option_index = 0;
	int c;
	int gap=0,match=1,mismatch=-1, alignmentType=0; 
	/*
	 * by default 
	 * gap = 0 
	 * match = 1 
	 * mismatch = -1 
	 * alignmenttype = 0 
	 * 	0 --> global
	 *  1 -> local
     *  2 -> semi-global */

	/* Detect the end of the options. */
	
	while ((c = getopt_long (argc, argv, "vh:a:m:n:g:",
				   long_options, &option_index)) != -1){
	
		switch (c) {
			case 'v':
			  ProjectVersion(); 
			  break;

			case 'h':
			  Help(); 
			  break;

			case 'a':
				alignmentType = atoi(optarg);
				/*
				 * 0 -> global
				 * 1 -> local
				 * 2 -> semi-global */
				
			  break;

			case 'm':
				match = atoi(optarg);
			  break;

			case 'n':
				mismatch = atoi(optarg);
				break;
			  
			case 'g':
				gap = atoi(optarg);
				break;

			case '?':
			  /* getopt_long already printed an error message. */
			  break;

			default:
			  //abort();
			  break;	
		}
	}
	
	/* Print any remaining command line arguments (not options). */
	if (optind < argc){		
		std::vector<std::unique_ptr<Sequence>> s1 ,s2; 
		
		//first file is always fasta file
		auto p1 = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind]);
		s1 = p1->Parse(-1);  
		
		//second file can be fasta or fastq file
		if (IsFastaFile(argv[optind+1])){
			auto p2 = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[optind+1]);
            s2 = p2->Parse(-1);
			cout << "fasta file" << endl; 
            PrintStats(s2);
		}
		else if (IsFastqFile(argv[optind+1])){
			auto p2 = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(argv[optind+1]);

			// parse in chunks
			std::vector<std::unique_ptr<Sequence>> s2;
			std::uint32_t chunk_size = 500 * 1024 * 1024;  // 500 MB
			for (auto t = p2->Parse(chunk_size); !t.empty(); t = p2->Parse(chunk_size)) {
			  s2.insert(
				  s2.end(),
				  std::make_move_iterator(t.begin()),
				  std::make_move_iterator(t.end()));
			}
			 PrintStats(s2);
		}
		else {
			//file format wrong
		}
		
		//alignment
		std::vector<std::unique_ptr	<Sequence>> shortFragments;  // fragments with length < 5000
		std::vector<std::unique_ptr	<Sequence>> longFragments;

		
		for (int i=0; i< s2.size(); i++){
			if (&s2[i]==nullptr){
				cout << "null pointer" << endl;
			} 
			if (s2[i]->data.size() < 5000) {
				//vec_a.at(0).get();
				shortFragments.push_back(std::move(s2[i]));
			}else{
				longFragments.push_back(std::move(s2[i]));
			}
		}
		
		//cout << "short fragment size.."<< shortFragments.size() << endl; 
		
		srand((unsigned int)time(NULL));
		
		int randomIndex1 = rand() % shortFragments.size();
		int randomIndex2 = rand() % shortFragments.size();
		while (randomIndex1 == randomIndex2){
			randomIndex2 = rand() % shortFragments.size();
		}
		
		//cout << "marker.." << randomIndex1 << " " << randomIndex2 << endl; 
		
		unsigned int target_begin; 
		std::string cigar = "";  
		//alignment type
		
				
		AlignmentType type; 
		switch (alignmentType){
			case 0:
				type = Global; 
				break;
			case 1:
				type = Local; 
				break; 
			case 2: 
				type = Semiglobal; 
				break; 
			default:
				break; 
		}
		
	
		const char* query = shortFragments[randomIndex1]->data.c_str(); 
		const char* target = shortFragments[randomIndex2]->data.c_str(); 
		
		int alignmentScore = Align(
			query, shortFragments[randomIndex1]->data.length(),
			target, shortFragments[randomIndex2]->data.length(),
			type,
			match,
			mismatch,
			gap,
			&cigar,
			&target_begin); 

		
		PrintAlignmentResult(
			shortFragments[randomIndex1]->name, shortFragments[randomIndex1]->data.length(),
			shortFragments[randomIndex2]->name, shortFragments[randomIndex2]->data.length(),
			type,
			match,
			mismatch,
			gap,
			cigar, 
			target_begin, 
			alignmentScore); 
			

		 
		unsigned int kmer_len = 15, window_len = 5;
		float f = 0.001; 
		
		//find distinct 

		map<unsigned int,vector<tuple<unsigned int,bool>>> referenceIndex, fragmentIndex; 
		
		//reference genome index
		getMinimizer(s1[0], &referenceIndex, kmer_len, window_len); 


		cout << "num of long fragments: " << longFragments.size() << endl;
		cout << "size of short fragment: " << shortFragments.size() << endl; 

		//fragments genome index
		/*
		for (int i=0; i<longFragments.size();i++){
			getMinimizer(longFragments[i], &fragmentIndex, kmer_len, window_len); 
		}*/

		vector<tuple<std::string, map<unsigned int,vector<tuple<unsigned int,bool>>>>> allFragmentIndex;

		for (int i=0; i<shortFragments.size();i++){
			fragmentIndex.clear(); 

			getMinimizer(shortFragments[i], &fragmentIndex, kmer_len, window_len);
			allFragmentIndex.push_back(make_tuple(shortFragments[i]->name, fragmentIndex)); 
			cout << "find minimizer for fragments ..." << shortFragments[i]->name << endl; 
		}

		cout<< "fragment index done..." << endl;
		

		//occurrences 
		vector<tuple<int, unsigned int>> occurrencesReferenceIndex, occurrencesFragmentIndex;
		occurrencesReferenceIndex = getOccurrences(referenceIndex); 
		occurrencesFragmentIndex = getOccurrences(fragmentIndex);

		//singleton count
		int referenceSingletonCount=getSingletonCount(occurrencesReferenceIndex);
		int fragmentSingletonCount = getSingletonCount(occurrencesFragmentIndex);

		//print result
		
		cout << "In reference genome: " << endl; 
		cout << "num minimizer:" << referenceIndex.size() << endl;
		cout << "num singleton: " << referenceSingletonCount << endl; 
		cout << "Singleton Fraction of refence genome: " << (float) referenceSingletonCount/referenceIndex.size() << endl;
		cout << "number of occurrences of the most frequent minimizer: " << getNumOccurrencesMostFrequentMinimizer(f,occurrencesReferenceIndex) << endl;
		

		/*
		cout << "In fragment genome: " << endl; 
		cout << "num minimizer:" << fragmentIndex.size() << endl;
		cout << "num singleton: " << fragmentSingletonCount << endl; 
		cout << "Singleton Fraction of refence genome: " << (float) fragmentSingletonCount/fragmentIndex.size() << endl;
		cout << "number of occurrences of the most frequent minimizer: " << getNumOccurrencesMostFrequentMinimizer(f,occurrencesFragmentIndex) << endl;*/
		cout << "num fragments" << allFragmentIndex.size() << endl; 

		//mapper 
		vector<tuple<unsigned int, unsigned int>> matchTable;
		int t_begin, t_end, q_begin, q_end, lenLIS; 
		//for each fragment 
		for (int i=0; i< allFragmentIndex.size(); i++){
			matchTable = matchMinimizer(referenceIndex, get<1>(allFragmentIndex[i])); 
			cout << "before longest linear chain" << endl ; 
			cout << i <<  " match table size " << matchTable.size() << endl; 

			//sort fragment position
			if (matchTable.size()>1)
				std::sort(matchTable.begin(), matchTable.end()); 
		
			//find longest linear chain
			lenLIS = LongestIncreasingSubsequence(matchTable, matchTable.size(), t_begin,
                                         t_end, q_begin, q_end);

			cout << lenLIS << " -> lens of longest increasing subsequence" << endl;

			/*
			//then do alignment (global alignment)
			string cigar; 
			unsigned int target_begin;
			cout << "REF q_begin: " << q_begin << endl; 
			cout << "REF q_end: " << q_end << endl; 		
			cout << "FRAG t_begin: " << t_begin << endl; 
			cout << "FRAG t_end: " << t_end << endl; 


			cout << "name in fragment index: " <<get<0>(allFragmentIndex[i]) << endl; 
			cout << "name in the short fragment: " <<shortFragments[i]->name << endl; 
			cout << "-------------------------------------------------------" << endl; 

			int result = 
			Align(shortFragments[i]->data.substr(q_begin,q_end).c_str(),
	(unsigned int) q_end - q_begin, s1[0]->data.substr(t_begin,t_end).c_str(),(unsigned int) t_end-t_begin,
	Global,
	match,
	mismatch,
	gap,
	&cigar,
	&target_begin);  
			cout << "alignment results: " <<  result << endl;
			cout << "#########################################" << endl; 
			//find fragment in long fragment/short fragment 
			//find substring in reference genome and then do alignment*/

			//in PAF format
			
		}

	}
	return 0;
}



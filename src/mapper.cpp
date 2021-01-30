#include <iostream>
#include <cstring>
#include <getopt.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <bits/stdc++.h>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "../3rdparty/bprinter/include/bprinter/table_printer.h"
#include "thread_pool/thread_pool.hpp"

#include "config.h"

#include "alignment.h"
#include "minimizer_binary.h"

using namespace std;

struct Sequence {  
	std::string name;
    std::string data;   
    std::string quality;
 public:
	 Sequence(const char* name, std::uint32_t name_len, const char* data,
           std::uint32_t data_len)
      : name(name, name_len), data(data, data_len) {}
 
	Sequence( 
      const char* name, std::uint32_t name_len,
      const char* data, std::uint32_t data_len,
      const char* quality, std::uint32_t quality_len):
		name(name, name_len),
        data(data, data_len),
        quality(quality, quality_len) {}
}; 



void Help(){
	cerr << "some help command..." << endl;
};

void ProjectVersion(){
	cerr << "v" << PROJECT_VER << endl;
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
    
    cerr << "---------------Statistics---------------" << endl; 
	cerr << "Number of sequences: " <<  fragments.size() << std::endl;
    cerr << "Total length of all fragments: " << lengthAllFragments << std::endl;
    cerr << "Largest fragment: " << maxLengthFrag << std::endl;
    cerr << "  length: " << maxLength << std::endl;
    cerr << "Smallest fragment: " << minLengthFrag << std::endl;
    cerr << "   length: " << minLength << std::endl;
    cerr << "Average length: " << ((double) lengthAllFragments) / fragments.size() << std::endl;
    cerr << "N50 length: " << N50<< std::endl;
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
			
			cerr << "---------------Alignment---------------" << endl; 
			cerr << "Query: " << query << endl; 
			cerr << "Query len: " << query_len << endl; 
			cerr << "Target: " << target << endl; 
			cerr << "Target len: " << target_len << endl; 
			cerr << "Alignment type: " << type << endl; 
			cerr << "match: " << match << endl;
			cerr << "mismatch: " << mismatch << endl; 
			cerr << "gap: " << gap << endl << endl;
			cerr << "---------------Results---------------" << endl; 
			cerr << "alignment score: " << alignmentScore << endl;
			cerr << "taget begin: " << target_begin << endl << endl; 
			//cerr << "cigar: " << cigar << endl;
			
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

void getMinimizer(std::unique_ptr<Sequence> & seq, std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>> *minimizerIndex, int kmer_len, int window_len){
	vector<std::tuple<unsigned int, unsigned int, bool>> minimizer;
	std:vector<std::tuple<unsigned int, unsigned int, bool>>::iterator itm; 
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	tuple<unsigned int, bool> temp; 
	vector<tuple<unsigned int, bool>> tempVector; 

	minimizer = MinimizeBinary(seq->data.c_str(),seq->data.length(),kmer_len,window_len);
	
	for (itm=minimizer.begin(); itm!=minimizer.end(); ++itm){
		tempVector.clear();
		it = minimizerIndex->find(get<0>(*itm));
		temp = make_tuple(get<1>(*itm),get<2>(*itm));
		if ( it == minimizerIndex->end() ) {
		// not found
			tempVector.emplace_back(temp); 
			minimizerIndex->insert(pair<unsigned int,vector<tuple<unsigned int,bool>>>
			(get<0>(*itm),tempVector)); 

		} else {
		// found
			(it->second).emplace_back(temp);
		}	
	} 
}

vector<tuple<int, unsigned int>> getOccurrences(std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>& minimizerIndex){
	vector<tuple<int, unsigned int>> occurrences; 
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	for (it = minimizerIndex.begin(); it != minimizerIndex.end(); ++it) { 
		occurrences.emplace_back(make_tuple((it->second).size(), it->first)); 
	}
	return occurrences; 
}

int getSingletonCount(vector<tuple<int, unsigned int>>& occurrences){
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	int singletonCount=0;
	for (int i=0; i<occurrences.size();i++){
		if (get<0>(occurrences[i])==1){
			singletonCount++;
		}
	}
	return singletonCount; 
}

int getNumOccurrencesMostFrequentMinimizer(float f,vector<tuple<int, unsigned int>>& occurrences){
	std::sort(occurrences.begin(),occurrences.end()); 
	int index = (int) occurrences.size()*f;
	return get<0>(occurrences[occurrences.size()-1-index]);
}

void ignoreTooFrequentMinimizer(float f,vector<tuple<int, unsigned int>>& occurrences, std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>& minimizerIndex){
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it; 
	std::sort(occurrences.begin(),occurrences.end()); 
	int index = (int) occurrences.size()*f;
	int ignoreFrom = occurrences.size()-index; 
	for (int i=ignoreFrom; i<occurrences.size();i++){
		minimizerIndex.erase(get<1>(occurrences[i])); 
	}
}

vector<tuple<unsigned int, unsigned int>> matchMinimizer(std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>& referenceIndex, 
std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>& fragmentIndex){
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;
	std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator search;
	vector<tuple<unsigned int, unsigned int>> match; 
	vector<unsigned int> fragmentOrigin, fragmentReverse, refOrigin, refReverse; 

	for (it = fragmentIndex.begin(); it != fragmentIndex.end(); ++it){
		fragmentOrigin.clear();
		fragmentReverse.clear();
		refOrigin.clear(); 
		refReverse.clear();
		search = referenceIndex.find(it->first); 
		if ( search != referenceIndex.end() ){
			//found
			for (int i=0; i< it->second.size(); i++){
				if (get<1>(it->second[i])){
					//if origin
					fragmentOrigin.emplace_back(get<0>(it->second[i]));
				}
				else{
					fragmentReverse.emplace_back(get<0>(it->second[i]));
				}				
			}

			for (int i=0; i<search->second.size();i++ ){
				//if origin
				if (get<1>(search->second[i])){
					refOrigin.emplace_back(get<0>(search->second[i]));
				}
				else{
					refReverse.emplace_back(get<0>(search->second[i]));  
				}
			}
			
			//sort refOrign, refReverse
			std::sort(refOrigin.begin(), refOrigin.end()); 
			std::sort(refReverse.begin(), refReverse.end()); 

			for (int i=0; i< fragmentOrigin.size(); i++){
				for (int r=0; r<refOrigin.size(); r++){
					match.emplace_back(make_tuple(fragmentOrigin[i], refOrigin[r])); 
				}
			}

			for (int i=0; i< fragmentReverse.size(); i++){
				for (int r=0; r<refReverse.size(); r++){
					match.emplace_back(make_tuple(fragmentReverse[i], refReverse[r])); 
				}
			}

			
		}
	}
	std::sort(match.begin(), match.end()); 
	return match; 
}

int LongestIncreasingSubsequence(vector<tuple<unsigned int, unsigned int>> const& matches, int &t_begin, int &t_end, int &q_begin,int &q_end) {
	int n = matches.size();
    const int INF = 1e9;
    vector<int> parent(n+1, INF);  //Tracking the predecessors/parents of elements of each subsequence.
	vector<int> increasingSub(n+1, INF); //Tracking ends of each increasing subsequence.
    int length = 0; //Length of longest subsequence.

	if (n==0)
		return 0; 

	if (n==1){
		t_begin = get<1>(matches[0]);
		q_begin = get<0>(matches[0]);
		t_end = get<1>(matches[0]);
		q_begin = get<0>(matches[0]);
		return n; 
	}

	for(int i=0; i< n; i++)
	{
		//Binary search
		int low = 1;
		int high = length;
		while(low <= high)
		{
			int mid = (int) ceil((low + high)/2);
				
			if(get<1>(matches[increasingSub[mid]]) < get<1>(matches[i]))
				low = mid + 1;
			else
				high = mid - 1;
		}
			
		int pos = low;
		//update parent/previous element for LIS
		parent[i] = increasingSub[pos-1];
		//Replace or append
		increasingSub[pos] =  i;
			
		//Update the length of the longest subsequence.
		if(pos > length)
			length=pos;
	}
		
	//Generate LIS by traversing parent array
	vector<int> refList(length, INF); 
	vector<int> fragList(length,INF); 
	int k   = increasingSub[length];
	for(int j=length-1; j>=0; j--)
	{
		refList[j] =  get<1>(matches[k]);
		fragList[j] = get<0>(matches[k]); 
		k = parent[k];
	}

	t_begin = refList[0]; 
	t_end = refList[length-1]; 
	q_begin = fragList[0]; 
	q_end = fragList[length-1]; 
    return length;
}

string generatePAFString(string queryName,int queryLen, int queryStart, int queryEnd, 
	string targetName, int targetLen, int targetStart, int targetEnd, int alignmentScore, 
	int alignmentBlockLen, string* cigar){
	string paf;
	paf = queryName + "\t"+  to_string(queryLen) + "\t" + to_string(queryStart)
		+ "\t" + to_string(queryEnd)  + "\t" + targetName + "\t" + to_string(targetStart)
		+ "\t" + to_string(targetEnd) + "\t"+ to_string(alignmentScore) + "\t"
		+ to_string(alignmentBlockLen); 
	if (cigar!=nullptr){
		paf += ("\t" + *cigar + "\n"); 
	}
	else{
		paf += "\n"; 
	}
	//paf += "-----------------------------------------------------\n";   

	return paf; 
}

void getAlignmentBlockLengthAndMatchLength(string& cigar, int *abl, int *ml){
	
	int value=0,sum=0, match =0 ; 
	bool prevIsNum = false; 
	for (int i=0; i < cigar.length(); i++ ){ 

		if ( isdigit(cigar[i]) ){
			if (!prevIsNum){
				prevIsNum = true; 
				value = atoi(cigar.substr(i,1).c_str()); 
			}
			else{
				value = value *10; 
				value += atoi(cigar.substr(i,1).c_str());
			}
		} else{
			if (prevIsNum){
				sum += value; 
				prevIsNum=false; 
			}

			if (cigar[i]=='m'){
				match += value; 
			} 
		}
	}
	*abl = sum; 
	*ml = match; 
}

void mapping(vector<tuple<std::string, std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>>>& allFragmentIndex, std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>& referenceIndex
,int match,int mismatch,int gap, std::vector<std::unique_ptr<Sequence>>& s1, std::vector<std::unique_ptr<Sequence>>& shortFragments, int i, bool cigarNeeded, int kmer_len){
	vector<tuple<unsigned int, unsigned int>> matchTable;
	int t_begin, t_end, q_begin, q_end, lenLIS, result, alignmentRefLength, alignmentQueryLength; 
	string cigar;
	unsigned int target_begin; 
	int matchLength, alignmentBlockLen, overlapLen; 
	
	//for each fragment 
	matchTable = matchMinimizer(referenceIndex, get<1>(allFragmentIndex[i])); 

	//find longest linear chain
	lenLIS = LongestIncreasingSubsequence(matchTable, t_begin, t_end, q_begin, q_end);

	alignmentRefLength = t_end-t_begin; 
	alignmentQueryLength = q_end - q_begin; 

	if (lenLIS>0 && alignmentRefLength<100000 && cigarNeeded){
		result = 
		Align(shortFragments[i]->data.substr(q_begin,q_end).c_str(),
		(unsigned int) q_end - q_begin, s1[0]->data.substr(t_begin,t_end).c_str(),(unsigned int) t_end-t_begin,
		Global,match,mismatch,gap,&cigar,&target_begin);  

		//get match length and alignment block len 
		getAlignmentBlockLengthAndMatchLength(cigar, &alignmentBlockLen, &matchLength); 

		//in PAF format
		cout << generatePAFString(shortFragments[i]->name, shortFragments[i]->data.length(), q_begin, q_end, s1[0]->name, s1[0]->data.length(), t_begin, t_end, matchLength,alignmentBlockLen, &cigar); 
	}
	else if (!cigarNeeded || alignmentRefLength>=100000  ){
		//no need alignment 
		//lenLIS*k for num of seq matches
		//num of sequence matches --> overlap length
		overlapLen = max(alignmentQueryLength, alignmentRefLength);
		cout << generatePAFString(shortFragments[i]->name, shortFragments[i]->data.length(), q_begin, q_end, s1[0]->name, s1[0]->data.length(), t_begin, t_end, lenLIS*kmer_len, overlapLen, nullptr); 
	}
	else {
		cout << generatePAFString(shortFragments[i]->name, shortFragments[i]->data.length(), -1, -1, s1[0]->name, s1[0]->data.length(), -1, -1, 0, 0, nullptr);
	}
	
}

static struct option long_options[] = {
  /* These options don’t set a flag.
	 We distinguish them by their indices. */
  {"version", no_argument,       0, 'v'},
  {"help",  no_argument,       0, 'h'},
  {"alignment_type",  required_argument, 0, 'a'},
  {"match",  required_argument, 0, 'm'},
  {"nomatch",    required_argument, 0, 'n'},
  {"gap",    required_argument, 0, 'g'},
  {"cigar", no_argument, 0, 'c'}, 
  {"kmer_len", required_argument, 0, 'k'}, 
  {"window_len", required_argument, 0, 'w'}, 
  {"thread_num", required_argument, 0, 't'}, 
  {"frequent", required_argument, 0, 'f'},
  {0, 0, 0, 0}
};

int main (int argc, char **argv){
	/* getopt_long stores the option index here. */
    int option_index = 0;
	int c=0; 
	bool cigarNeeded=false;
	int gap=0,match=1,mismatch=-1, alignmentType=0; 
	int t= 5; //number of threads
	unsigned int kmer_len = 15, window_len = 5;
	float f = 0.001; 
	/*
	 * by default 
	 * 	0 --> global
	 *  1 -> local
     *  2 -> semi-global */

	while ((c = getopt_long (argc, argv, "vhc:a:m:n:g:k:w:t:f:",
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
			
			case 'c':
				cigarNeeded=true; 
				break; 
			case 'k':
				kmer_len = (unsigned int) atoi(optarg); 
				break;
			case 'w':
				window_len = (unsigned int) atoi(optarg); 
				break; 
			case 't':
				t= atoi(optarg); 
				break; 
			case 'f': 
				f = std::stof(optarg); 
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
				cerr << "null pointer" << endl;
			} 
			if (s2[i]->data.size() < 5000) {
				shortFragments.emplace_back(std::move(s2[i]));
			}else{
				longFragments.emplace_back(std::move(s2[i]));
			}
		}
		
		
		srand((unsigned int)time(NULL));
		
		int randomIndex1 = rand() % shortFragments.size();
		int randomIndex2 = rand() % shortFragments.size();
		while (randomIndex1 == randomIndex2){
			randomIndex2 = rand() % shortFragments.size();
		}
	
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
		int alignmentScore; 
		if (cigarNeeded){
			unsigned int target_begin; 
			std::string cigar = "";
			alignmentScore = Align(query, shortFragments[randomIndex1]->data.length(),
			target, shortFragments[randomIndex2]->data.length(),
			type,match,mismatch,gap,&cigar,&target_begin); 
		}
		else{
			alignmentScore = Align(query, shortFragments[randomIndex1]->data.length(),
				target, shortFragments[randomIndex2]->data.length(),
				type,match,mismatch,gap,nullptr,nullptr); 
		}


		unsigned int target_begin; 
		std::string cigar = "";

		PrintAlignmentResult(shortFragments[randomIndex1]->name, shortFragments[randomIndex1]->data.length(),
			shortFragments[randomIndex2]->name, shortFragments[randomIndex2]->data.length(),
			type,match,mismatch,gap,cigar, target_begin, alignmentScore); 
		
	
		//find distinct 

		std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>> referenceIndex, fragmentIndex; 
		std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>::iterator it;

		//reference genome index
		getMinimizer(s1[0], &referenceIndex, kmer_len, window_len); 

		//fragments genome index
		/*
		for (int i=0; i<longFragments.size();i++){
			getMinimizer(longFragments[i], &fragmentIndex, kmer_len, window_len); 
		}*/

		vector<tuple<std::string, std::unordered_map<unsigned int,vector<tuple<unsigned int,bool>>>>> allFragmentIndex;

		for (int i=0; i<shortFragments.size();i++){
			fragmentIndex.clear(); 

			getMinimizer(shortFragments[i], &fragmentIndex, kmer_len, window_len);
			allFragmentIndex.emplace_back(make_tuple(shortFragments[i]->name, fragmentIndex)); 
		}

		//occurrences 
		vector<tuple<int, unsigned int>> occurrencesReferenceIndex, occurrencesFragmentIndex;
		occurrencesReferenceIndex = getOccurrences(referenceIndex); 

		//singleton count
		int referenceSingletonCount=getSingletonCount(occurrencesReferenceIndex);
		//print result
		
		cerr << "------------------------------------------------------------------------------" << endl; 
		cerr << "In reference genome: " << endl; 
		cerr << "num minimizer:" << referenceIndex.size() << endl;
		cerr << "num singleton: " << referenceSingletonCount << endl; 
		cerr << "Singleton Fraction of refence genome: " << (float) referenceSingletonCount/referenceIndex.size() << endl;
		cerr << "number of occurrences of the most frequent minimizer: " << getNumOccurrencesMostFrequentMinimizer(f,occurrencesReferenceIndex) << endl;
		ignoreTooFrequentMinimizer(f, occurrencesReferenceIndex, referenceIndex); 

		int fragmentSingletonCount; 
		cerr << "fragment index size: "<< allFragmentIndex.size(); 

		for (int i =0; i<allFragmentIndex.size();i++){
			cerr << "sheme..." << endl; 
			cerr << "fragment name:" <<  get<0>(allFragmentIndex[i])<< endl; 
			cerr << "num minimizer:" << get<1>(allFragmentIndex[i]).size() << endl;
			occurrencesFragmentIndex = getOccurrences(get<1>(allFragmentIndex[i]));
			cerr << "$$ occurance fragment index" << occurrencesFragmentIndex.size() << endl; 
			fragmentSingletonCount = getSingletonCount(occurrencesFragmentIndex);
			cerr << "------------------------------------------------------------------------------" << endl; 
			cerr << i << endl;
			cerr << "fragment name:" <<  get<0>(allFragmentIndex[i])<< endl; 
			cerr << "num minimizer:" << get<1>(allFragmentIndex[i]).size() << endl;
			cerr << "num singleton: " << fragmentSingletonCount << endl; 
			cerr << "Singleton Fraction of refence genome: " << (float) fragmentSingletonCount/get<1>(allFragmentIndex[i]).size()  << endl;
			cerr << "number of occurrences of the most frequent minimizer: " << getNumOccurrencesMostFrequentMinimizer(f,occurrencesFragmentIndex) << endl;

			//ignore too frequent minimizer
			ignoreTooFrequentMinimizer(f, occurrencesFragmentIndex, get<1>(allFragmentIndex[i])); 
		}

		cerr << "------------------------------------------------------------------------------" << endl; 
		
		//without multithreading 
		/*
		for (int i=0 ; i<allFragmentIndex.size();i++){
			cerr << i << " --->" << mapping(allFragmentIndex, referenceIndex ,match, mismatch, gap, s1,shortFragments, i, cigarNeeded); 
		}*/

		
		auto thread_pool = thread_pool::ThreadPool(t);
		std::vector<std::future<void>> void_futures;

		for (int i = 0; i < allFragmentIndex.size(); i++) {
			void_futures.emplace_back(thread_pool.Submit(mapping, std::ref(allFragmentIndex),std::ref(referenceIndex),
										 match, mismatch, gap, std::ref(s1),std::ref(shortFragments), i, cigarNeeded, kmer_len));
		}

		for (const auto& it : void_futures) {
			it.wait();
		}

	}
	return 0;
}



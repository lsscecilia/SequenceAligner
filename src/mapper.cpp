#include <iostream>
#include <cstring>
#include <getopt.h>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

#include "config.h"

#include "alignment.h"

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
		std::vector<std::unique_ptr<Sequence>> shortFragments;  // fragments with length < 5000
		for (int i=0; i< s2.size(); i++){
			if (s2[i]->data.size() < 5000) {
				shortFragments.push_back(std::move(s2[i]));
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
			
		 

	}
	return 0;
}



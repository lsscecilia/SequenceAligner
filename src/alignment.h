std::string CompressCigar(std::string uCigar); 

int SemiGlobalAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin);
    
int SmithWatermanAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin); 

int NeedlemanWunschAlgorithm(
	const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap, std::string* cigar, unsigned int* target_begin); 

	
enum AlignmentType {Global, Local, Semiglobal}; 


int Align(
	const char* query, unsigned int query_len,
	const char* target, unsigned int target_len,
	AlignmentType type,
	int match,
	int mismatch,
	int gap,
	std::string* cigar,
	unsigned int* target_begin); 

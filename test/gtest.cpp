#include "gtest/gtest.h"
#include "alignment.h"

TEST (AlignmentTest,GlobalAlignmentTest){
	char * query ="CTCTGTTCG"; 
	char * target = "CGTATCTTGA"; 
	
	std::string cigar; 
	unsigned int target_begin=0; 
	
	EXPECT_EQ(-5,NeedlemanWunschAlgorithm(query,9,target,10,
	AlignmentType::Global,0,-1,-1,cigar,target_begin)); 
}

int main(int argc, char **argc){
	::testing::InitGoogleTest(&argc, argv); 
	return RUN_ALL_TESTS();
}

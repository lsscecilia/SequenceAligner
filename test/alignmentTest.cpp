#include "gtest/gtest.h"

#include "alignment.h"
#include "minimizer.h"

TEST (AlignmentTest,GlobalAlignmentTest){
	const char * query ="CTCTGTTCG"; 
	const char * target = "CGTATCTTGA"; 
	
	std::string cigar; 
	unsigned int target_begin=0; 
	
	EXPECT_EQ(-5,NeedlemanWunschAlgorithm(query,9,target,10,0,-1,-1,&cigar,&target_begin)); 
}


TEST (AlignmentTest, LocalAligmentTest){
	std::string query  ="CTCTGAG"; 
	std::string target = "TGTCAGT"; 
	
	std::string cigar; 
	unsigned int target_begin = 0; 
	EXPECT_EQ(6,SmithWatermanAlgorithm(query.c_str(),query.length(),target.c_str(),target.length(),
			2,-2,-1,&cigar,&target_begin)); 
}

TEST (AlignmentTest, SemiGlobalAlignmentTest){
	const char* query = "AGCATGCAAT"; 
	const char* target = "ATCCGAACATCCAATCGAAGC"; 
	
	std::string cigar;
	unsigned int target_begin =0; 
	EXPECT_EQ(14,SemiGlobalAlgorithm(query,10,target,21,2,-1,-1,&cigar,&target_begin)); 
}
	
	


TEST (Minimizers, TestTGACGTACATGGACA) {
	std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = Minimize("TGACGTACATGGACA", 15, 3, 4);

	//std::sort(minimizers.begin(), minimizers.end());

	EXPECT_EQ(std::get<0>(minimizers.at(0)), 2);
	EXPECT_EQ(std::get<1>(minimizers.at(0)), 10);
	EXPECT_EQ(std::get<2>(minimizers.at(0)), 0);
	
	EXPECT_EQ(std::get<0>(minimizers.at(1)), 12);
	EXPECT_EQ(std::get<1>(minimizers.at(1)), 4);
	EXPECT_EQ(std::get<2>(minimizers.at(1)), 0);
	
	EXPECT_EQ(std::get<0>(minimizers.at(2)), 12);
	EXPECT_EQ(std::get<1>(minimizers.at(2)), 7);
	EXPECT_EQ(std::get<2>(minimizers.at(2)), 1);
	
	EXPECT_EQ(std::get<0>(minimizers.at(3)), 23);
	EXPECT_EQ(std::get<1>(minimizers.at(3)), 1);
	EXPECT_EQ(std::get<2>(minimizers.at(3)), 0);
	
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

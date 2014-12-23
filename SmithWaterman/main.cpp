// Implement Smith-Waterman Local Alignment Algorithm
// http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

// This does an alignment and scoring of hemoglobin amino acids.
// With enough scoring data, the taxonomic position of a species
// relative to other species can be inferred.

#include "PreCompile.h"
#include "ScorePolicy.h"
#include "SmithWaterman.h"

//---------------------------------------------------------------------------
// These sample hemoglobins were taken from ExPASy.org/SwissProt.
// Sample hemoglobins could also be taken from uniprot.org or NCBI/BLAST.

// >sp|P68871|HBB_HUMAN Hemoglobin subunit beta (Hemoglobin beta chain) (Beta-globin) - Homo sapiens (Human).
static std::string HBB_HUMAN("VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV"
                             "KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
                             "EFTPPVQAAYQKVVAGVANALAHKYH");

// >sp|P68873|HBB_PANTR Hemoglobin subunit beta (Hemoglobin beta chain) (Beta-globin) - Pan troglodytes (Chimpanzee).
static std::string HBB_PANTR("VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV"
                             "KAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGK"
                             "EFTPPVQAAYQKVVAGVANALAHKYH");

// >sp|P02088|HBB1_MOUSE Hemoglobin beta-1 subunit (Hemoglobin beta-1 chain) (Beta-1-globin) (Hemoglobin beta-major chain) - Mus musculus (Mouse).
static std::string HBB1_MOUSE("VHLTDAEKAAVSCLWGKVNSDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNAKV"
                              "KAHGKKVITAFNDGLNHLDSLKGTFASLSELHCDKLHVDPENFRLLGNMIVIVLGHHLGK"
                              "DFTPAAQAAFQKVVAGVATALAHKYH");

// >sp|P02112|HBB_CHICK Hemoglobin subunit beta (Hemoglobin beta chain) (Beta-globin) - Gallus gallus (Chicken).
static std::string HBB_CHICK("VHWTAEEKQLITGLWGKVNVAECGAEALARLLIVYPWTQRFFASFGNLSSPTAILGNPMV"
                             "RAHGKKVLTSFGDAVKNLDNIKNTFSQLSELHCDKLHVDPENFRLLGDILIIVLAAHFSK"
                             "DFTPECQAAWQKLVRVVAHALARKYH");

// >tr|Q802A3|Q802A3_FUGRU Hemoglobin beta subunit - Fugu rubripes (Japanese pufferfish) (Takifugu rubripes).
static std::string Q802A3_FUGRU("MVEWTDQERTIISNIFSTLDYEDVGSKSLIRCLIVYPWTQRYFAGFGNLYNAEAIKNNPN"
                                "IAKHGVTVLHGLDRAVKNMDNIKETYKELSELHSEKLHVDPDNFKLLSDCLTIVVATKMG"
                                "SKFTPEIQATFQKFLAVVVSALGRQYH");

// >tr|Q540F0|Q540F0_VIGUN Leghaemoglobin - Vigna unguiculata (Cowpea).
static std::string Q540F0_VIGUN("MVAFSDKQEGLVNGAYEAFKADIPKYSVVFYTTILEKAPAAKNLFSFLANGVDATNPKLT"
                                "GHAEKLFGLVRDSAAQLRASGGVVADAALGAVHSQKAVNDAQFVVVKEALVKTLKEAVGD"
                                "KWSDELGTAVELAYDELAAAIKKAY");

// >sp|P51460|INSL3_HUMAN Insulin-like 3 precursor (Leydig insulin-like peptide) (Ley-I-L) (Relaxin-like factor)
// [Contains: Insulin-like 3 B chain; Insulin-like 3 A chain] - Homo sapiens (Human).
static std::string INSL3_HUMAN("MDPRLPAWALVLLGPALVFALGPAPTPEMREKLCGHHFVRALVRVCGGPRWSTEARRPAA"
                               "GGDRELLQWLERRHLLHGLVADSNLTLGPGLQPLPQTSHHHRHHRAAATNPARYCCLSGC"
                               "TQQDLLTLCPY");

//---------------------------------------------------------------------------
int main()
{
#ifdef _DEBUG
    // Exercise the local alignment algorithm with sample vectors.
    {
        static std::string test_vector1("xxxcde");
        static std::string test_vector2("abcxdex");

        std::cout << "Aligning " << test_vector1 << " and " << test_vector2 << ":\n";
        Alignment_table table(test_vector1, test_vector2, &basic_calc_score);

        table.print_table(std::cout);

        table.print_trace_back(std::cout);
    }
#endif

    {
        static std::string sequence1("deadly");
        static std::string sequence2("ddgearlyk");

        std::cout << "\nAligning " << sequence1 << " and " << sequence2 << ":\n";
        Alignment_table table(sequence1, sequence2, &BLOSUM62_calc_score<-4>);

        table.print_table(std::cout);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, 1000);
    }

    const unsigned int num_permutations = 10000;

    {
        std::cout << "\nAligning  HBB_HUMAN and HBB_PANTR:\n";
        Alignment_table table(HBB_HUMAN, HBB_PANTR, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    {
        std::cout << "\nAligning and HBB1_MOUSE:\n";
        Alignment_table table(HBB_HUMAN, HBB1_MOUSE, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    {
        std::cout << "\nAligning HBB_HUMAN and HBB_CHICK:\n";
        Alignment_table table(HBB_HUMAN, HBB_CHICK, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    {
        std::cout << "\nAligning HBB_HUMAN and Q802A3_FUGRU:\n";
        Alignment_table table(HBB_HUMAN, Q802A3_FUGRU, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    {
        std::cout << "\nAligning HBB_HUMAN and Q540F0_VIGUN:\n";
        Alignment_table table(HBB_HUMAN, Q540F0_VIGUN, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    {
        std::cout << "\nAligning HBB_HUMAN and INSL3_HUMAN:\n";
        Alignment_table table(HBB_HUMAN, INSL3_HUMAN, &BLOSUM62_calc_score<-4>);

        table.print_trace_back(std::cout);
        table.calc_pvalue(std::cout, num_permutations);
    }

    std::cout << "Program done." << std::endl;
    return 0;
}

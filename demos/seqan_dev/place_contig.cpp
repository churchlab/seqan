#include <iostream>
#include <seqan/align.h>
#include <seqan/align_split.h>
#include <seqan/file.h>  // output of String
#include <seqan/sequence.h>
#include <seqan/score.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("place_contig");

    addOption(parser, ArgParseOption(
        "ref", "reference", "Reference genome fasta.",
        ArgParseArgument::STRING, "TEXT"));
    addOption(parser, ArgParseOption(
        "read", "read", "Read fasta.",
        ArgParseArgument::STRING, "TEXT"));
    addOption(parser, ArgParseOption(
        "v", "verbose", "Give verbose output."));
    
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    // Extract option values and print them.
    CharString ref_path;
    getOptionValue(ref_path, parser, "ref");
    CharString read_path;
    getOptionValue(read_path, parser, "read");

    bool isVerbose = isSet(parser, "verbose");

    // Read in files
    CharString ref_id;
    Dna5String ref;

    SeqFileIn refFileIn(toCString(ref_path));
    readRecord(ref_id, ref, refFileIn);

    CharString read_id;
    Dna5String read;

    SeqFileIn readFileIn(toCString(read_path));
    readRecord(read_id, read, readFileIn);

    // Prepare Gaps objects.  We need one for the left part and one for the
    // right part of the alignment.
    Align<Dna5String> alignL;
    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), ref);
    assignSource(row(alignL, 1), read);
    Align<Dna5String> alignR;
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), ref);
    assignSource(row(alignR, 1), read);

    // Define scoring scheme.
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Call split alignment function.
    splitAlignment(alignL, alignR, scoringScheme);

    // Output refSplitPosition, readSplitLPosition, and readSplitRPosition on two lines
    int refSplitPosition = toSourcePosition(row(alignL, 0), clippedEndPosition(row(alignL, 0)));
    int readSplitLPosition = toSourcePosition(row(alignL, 1), clippedEndPosition(row(alignL, 1)));
    int readSplitRPosition = toSourcePosition(row(alignR, 1), 0);
    std::cout << refSplitPosition << "\n" << readSplitLPosition << "\n" << readSplitRPosition << "\n";

    // VERBOSE OUTPUT: Print resulting alignment to stdout.
    if(isVerbose) {
        std::cout << "Resulting alignments\n"
                  << "\n"
                  << "Left\n"
                  << alignL
                  << "Right\n"
                  << alignR
                  << "\n";

        SEQAN_ASSERT_EQ(refSplitPosition, toSourcePosition(row(alignR, 0), 0));

        std::cout << "refSplitPosition   == " << refSplitPosition << "\n"
                  << "readSplitLPosition == " << readSplitLPosition << "\n"
                  << "readSplitRPosition == " << readSplitRPosition << "\n\n";

        // Print sequence parts to stdout.
        std::cout << "Reference Left  " << prefix(ref, refSplitPosition) << "\n"
                  << "Reference Right " << suffix(ref, refSplitPosition) << "\n"
                  << "\n"
                  << "Read Left       " << prefix(read, readSplitLPosition) << "\n"
                  << "Read Center     " << infix(read, readSplitLPosition, readSplitRPosition) << "\n"
                  << "Read Right      " << suffix(read, readSplitRPosition) << "\n";
    }

    return 0;
}
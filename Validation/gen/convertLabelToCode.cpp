/**
 * @file convertLabelToCode.cpp
 * @brief Clones a TTree into a new file, replacing a `std::string` column with
 *        an integer one.
 * 
 * The original tree "label" branch is a `std::string` branch which has few
 * possible options (order of 100). The code produces a new tree clone of the
 * original one, but with a 32-bit label _code_ branch of hash codes instead of
 * the string label branch. An additional tree is also produced with in a branch
 * all the observed label codes and in the other the original strings.
 * 
 * The code is extracted with the standard C++ `std::hash<std::string>`
 * keeping only the 32 least significant bits.
 * 
 * 
 * Compilation and execution
 * --------------------------
 * 
 * Compile with:
 *     
 *     g++ -Wall -Wextra -pedantic -std=c++17 -O3 \
 *       -L$(root-config --libdir --cflags) -lCore -lRIO -lTree -lTreePlayer \
 *       -o convertLabelToCode.exe convertLabelToCode.cpp
 *     
 * Run without arguments for 
 * 
 */

#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


int main(int argc, char** argv) {
  
  using namespace std::string_literals;
  
  //
  // argument parsing
  //
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " input.root treeName labelBranchName [output.root] [newBranchName]\n";
    return (argc == 1)? 0: 1;
  }

  std::string const inputFilename = argv[1];
  std::string const treeName = argv[2];
  std::string const labelBranchName = argv[3];
  std::string const outputFilename = (argc >= 5) ? argv[4]
    : (inputFilename.substr(inputFilename.length()-5) + "_withLabelCode.root");
  std::string const newBranchName = (argc >= 6 && !std::string_view(argv[5]).empty())
                                    ? std::string{ argv[5] }
                                    : (labelBranchName + "code");

  std::cout 
    <<   "Input file: '" << inputFilename << "'"
    << "\nTree: '" << treeName << "'  string branch: '" << labelBranchName
      << "'  new branch: '" << newBranchName << "'"
    << std::endl;

  //
  // clone the input tree, except for the label branch
  //
  auto Fin = std::unique_ptr<TFile>{ TFile::Open(inputFilename.c_str(), "READ") };
  auto treeIn = Fin->Get<TTree>(treeName.c_str());
  treeIn->SetBranchStatus("*", 1);
  treeIn->SetBranchStatus("label", 0);
  auto Fout = std::unique_ptr<TFile>{ TFile::Open(outputFilename.c_str(), "RECREATE") };
  Fout->cd();
  TTree* treeOut = treeIn->CloneTree(-1, "fast");
  treeIn->SetBranchStatus("label", 1);
  unsigned int const nEntries = treeIn->GetEntries();

  //
  // create the label code branch and the code-to-label map
  //
  UInt_t labelCode;
  std::unordered_map<UInt_t, std::string> codeToLabel;
  codeToLabel.reserve(100); // just a guess
  TBranch* labelBranch = treeOut->Branch(newBranchName.c_str(), &labelCode);
  TTreeReader treeReader{ treeIn };
  TTreeReaderValue<std::string> labelReader(treeReader, "label");
  unsigned int iEntry = 0;
  unsigned int pageSize = nEntries / 10;
  std::hash<std::string> hasher;
  while (treeReader.Next()) {
    if (iEntry % pageSize == 0) std::cout << (iEntry * 100 / nEntries) << "%... " << std::flush;
    ++iEntry;
    labelCode = hasher(*labelReader) & 0xFFFF'FFFF;
    labelBranch->Fill();
    codeToLabel.try_emplace(labelCode, *labelReader);
  }
  std::cout << "100%!" << std::endl;
  Fin->Close();
  treeOut->Write();
  std::cout << "Written the new tree with " << treeOut->GetEntries()
    << " entries and label code branch '" << newBranchName
    << "' instead of label branch '" << labelBranchName << "'" << std::endl;
  
  //
  // print a summary of the labels and write the code-to-label map to a tree
  //
  std::cout << "Found " << codeToLabel.size() << " labels:";
  for (auto const& [ code, label ]: codeToLabel)
    std::cout << "\n - '" << label << "'  \t(0x" << std::hex << code << ")";
  std::cout << std::endl;
  TTree* mapTree = new TTree
    ((treeOut->GetName() + "_LabelMap"s).c_str(), "Map from label codes to original labels");
  std::string label;
  mapTree->Branch("labelCode", &labelCode);
  mapTree->Branch("label", &label);
  for (const auto& [ c, l ]: codeToLabel) {
    label = l;
    labelCode = c;
    mapTree->Fill();
  }
  mapTree->Write();
  
  //
  // close and be done
  //
  std::cout << "Closing output file: '" << Fout->GetName() << "'" << std::endl;
  Fout->Close();

  return 0;
} // main()

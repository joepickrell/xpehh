#include "ehhpop.h"
#include "CmdLine.h"
using namespace std;

int main(int argc, char *argv[]){
  CCmdLine cmdline;
  string pop1in, pop2in, mapin, snpin;
  if (cmdline.SplitLine(argc, argv) < 1)
    {
      // no switches were given on the command line, abort
      cout << "No input files\n";
      exit(-1);
    }
  if (!cmdline.HasSwitch("-h")){
    cout << "ERROR: no haplotype files\n";
    exit(1);
  }
  else{
    pop1in = cmdline.GetArgument("-h", 0);
    pop2in = cmdline.GetArgument("-h", 1);
  }
  if (!cmdline.HasSwitch("-m")){
    cout << "ERROR: no map file\n";
    exit(1);
  }
  else{
    mapin = cmdline.GetArgument("-m", 0);
  }
  if (!cmdline.HasSwitch("-s")){
    cout << "ERROR: no snp file\n";
    exit(1);
  }
  else{
    snpin = cmdline.GetArgument("-s", 0);
  }
  Ehhpop pop1(pop1in.c_str(), mapin.c_str());
  Ehhpop pop2(pop2in.c_str(), mapin.c_str());
  Ehhpop combined(pop1, pop2);
  vector<string> snps;
  ifstream infile(snpin.c_str());
  string st;
  if(infile.fail()){
    cerr << "ERROR: cannot open file " << snpin << "\n";
    exit(1);
  }
  while(getline(infile, st)){
    string buf;
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    snps.push_back(line.at(0));
  }
  int i; 
  for(i = 0; i<snps.size(); i++){
    string rs = snps[i];
    int index = combined.rs2index[rs];
    int pos = combined.rs2pos[rs];
    int left = combined.findcutoff(index, 0.05, false);
    int right = combined.findcutoff(index, 0.05, true);
    if (left == -1 || right ==-1){
      continue;
    }
    double p1_left = pop1.integrate_ehh(index, left, false);
    double p1_right = pop1.integrate_ehh(index, right, true);
    double IA = p1_left+p1_right;
    double p2_left = pop2.integrate_ehh(index, left, false);
    double p2_right = pop2.integrate_ehh(index, right, true);
    double IB = p2_left+p2_right;
    double ratio = IA/IB;
    double logratio = log(ratio);
    cout << rs << " " << pos << " "<< IA << " " << IB << " " << logratio << "\n";
  }
}
    

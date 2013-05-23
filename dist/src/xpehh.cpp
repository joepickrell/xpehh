#include "ehhpop.h"
#include "CmdLine.h"
using namespace std;

bool distcorrect = true;

int main(int argc, char *argv[]){
  CCmdLine cmdline;
  string pop1in, pop2in, mapin;
  if (cmdline.SplitLine(argc, argv) < 1)
    {
      // no switches were given on the command line, abort
      cout << "No input files\n";
      cout << "-h [input file 1] [input file 2]\n";
      cout << "-m [map file]\n";
      exit(-1);
    }
  if (!cmdline.HasSwitch("-h")){
    cout << "ERROR: no haplotype files\n";
    cout << "-h [input file 1] [input file 2]\n";
    cout << "-m [map file]\n";
    exit(1);
  }
  else{
    pop1in = cmdline.GetArgument("-h", 0);
    pop2in = cmdline.GetArgument("-h", 1);
  }
  if (!cmdline.HasSwitch("-m")){
    cout << "ERROR: no map file\n";
    cout << "-h [input file 1] [input file 2]\n";
    cout << "-m [map file]\n";
    exit(1);
  }
  else{
    mapin = cmdline.GetArgument("-m", 0);
  }
  if (cmdline.HasSwitch("-nd")) distcorrect = false;
  Ehhpop pop1(pop1in.c_str(), mapin.c_str());
  Ehhpop pop2(pop2in.c_str(), mapin.c_str());
  Ehhpop combined(pop1, pop2);
  int i;
  for(i = 0; i<combined.Nsnp; i++){
    string rs = combined.index2rs[i];
    int pos = combined.rs2pos[rs];
    int left = combined.findcutoff(i, 0.05, false, distcorrect);
    int right = combined.findcutoff(i, 0.05, true, distcorrect);
    //cout << left << " "<< right << "\n";
    if (left==-1 || right ==-1){
      continue;
    }
    double p1_left = pop1.integrate_ehh(i, left, false, distcorrect);
    double p1_right = pop1.integrate_ehh(i, right, true, distcorrect);
    cout << p1_left << " "<< p1_right << " p1\n";
    double IA = p1_left+p1_right;
    double p2_left = pop2.integrate_ehh(i, left, false, distcorrect);
    double p2_right = pop2.integrate_ehh(i, right, true, distcorrect);
    double IB = p2_left+p2_right;
    cout << p2_left << " "<< p2_right << " p2\n";
    double ratio = IA/IB;
    double logratio = log(ratio);
    cout << rs << " " << pos << " "<< IA << " " << IB << " " << logratio << "\n";
  }
}


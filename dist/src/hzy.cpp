#include "ehhpop.h"
#include "CmdLine.h"
using namespace std;

int main(int argc, char *argv[]){
  CCmdLine cmdline;
  string popin, mapin;
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
    popin = cmdline.GetArgument("-h", 0);
  }
  if (!cmdline.HasSwitch("-m")){
    cout << "ERROR: no map file\n";
    exit(1);
  }
  else{
    mapin = cmdline.GetArgument("-m", 0);
  }
  Ehhpop pop(popin.c_str(), mapin.c_str());
  int i; 
  for(i = 0; i<pop.Nsnp; i++){
    string rs = pop.index2rs[i];
    int pos = pop.rs2pos[rs];
    double hzy = pop.ehh(i, i);
    hzy = 1-hzy;
    cout << rs << " " << pos << " "<< hzy << "\n";
  }
}
    

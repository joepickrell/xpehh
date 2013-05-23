#include "ehhpop.h"

Ehhpop::Ehhpop(const char hapin[],  const char mapin[]){
  readhap(hapin);
  readmap(mapin);
}

const Ehhpop& Ehhpop::operator=(const Ehhpop& pop){
  haps = pop.haps;
  rs2pos = pop.rs2pos;
  rs2gpos = pop.rs2gpos;
  rs2anc = pop.rs2anc;
  index2rs = pop.index2rs;
  rs2index = pop.rs2index;
  Nhap = pop.Nhap;
  Nsnp = pop.Nsnp;
  return *this;
}

Ehhpop::Ehhpop( Ehhpop& pop1, Ehhpop&pop2)throw(){
  haps = pop1.haps;
  rs2pos = pop1.rs2pos;
  rs2gpos = pop1.rs2gpos;
  rs2anc = pop1.rs2anc;
  index2rs = pop1.index2rs;
  rs2index = pop1.rs2index;
  Nhap = pop1.Nhap +pop2.Nhap;
  Nsnp = rs2pos.size();
  vector<vector<char> >::iterator it;
  for (it = pop2.haps.begin(); it !=pop2.haps.end(); it++){
    haps.push_back(*it);
  }
}
void Ehhpop::readhap(const char in[]){

  //
  // read the haplotype data in from a file
  //

  ifstream infile(in);
  string st;
  if (infile.fail()){
    cerr << "ERROR: cannot open file " << in << "\n";
    exit(1);
  }
  while(getline(infile, st)){
    string buf;
    stringstream ss(st);
    vector<char> line;
    while (ss >> buf){
      line.push_back(buf[0]);
    }
    haps.push_back(line);
  }
  Nhap = haps.size();
}

void Ehhpop::printhap(){
  vector< vector<char> >::iterator it;
  for (it = haps.begin(); it!= haps.end(); it++){
    vector<char>::iterator it2;
    for (it2 = it->begin(); it2!=it->end(); it2++){
      char tmp2 = *it2;
      cout << tmp2 << " ";
    }
    cout << "\n";
  }
}

void Ehhpop::printmap(){
  int i;
  for (i = 0; i<Nsnp; i++){
    string rs = index2rs[i];
    double gpos = rs2gpos[rs];
    int pos = rs2pos[rs];
    bool anc = rs2anc[rs];
    cout << rs << " " << pos << " " << gpos << " " << anc  << "\n";
  }
}
void Ehhpop::readmap(const char in[]){
  ifstream infile(in);
  string st;
  if (infile.fail()){
    cerr << "ERROR: cannot open file " << in << "\n";
    exit(1);
  }
  int i = 0;
  while(getline(infile, st)){
    string buf;
    stringstream ss(st);
    vector<string> line;
    while (ss >> buf){
      line.push_back(buf);
    }
    string rs = line.at(0);
    int ppos = atoi(line.at(1).c_str());
    double gpos = atof(line.at(2).c_str());
    bool isanc = false;
    if (line.at(3)== line.at(4)){
      isanc = true;
    }
    rs2gpos.insert(make_pair(rs, gpos));
    rs2pos.insert(make_pair(rs, ppos));
    rs2anc.insert(make_pair(rs, isanc));
    index2rs.insert(make_pair(i, rs));
    rs2index.insert(make_pair(rs, i));
    i++;
  }
  Nsnp = i;
}

double Ehhpop::ehh( vector<vector< char> > subset){
  vector<vector<char > >::iterator it;
  int ident = 0;
  int total = 0;
  for (it = subset.begin(); it!=(subset.end()-1); it++){
    vector<char> tmp = *it;
    vector<vector<char> >::iterator it2;
    for (it2 = it+1; it2!=subset.end(); it2++){
      vector<char> tmp2 = *it2;
      total++;
      if (tmp == tmp2){
	ident++;
      }
    }
  }
  double ehh = (double) ident /(double) total;
  return(ehh);
}

double Ehhpop::ehh(int start, int end){
  vector<vector<char> > tmp;
  vector<vector<char > >::iterator it;
  for (it = haps.begin(); it!=haps.end(); it++){
    int i;
    vector<char> tmp2;
    for (i = start; i<=end; i++){
      tmp2.push_back(it->at(i));
    }
    tmp.push_back(tmp2);
  }
  double out = ehh(tmp);
  return(out);
}

double Ehhpop::integrate_ehh(int core, int stop, bool right, bool distcorrect){
  double total = 0;
  double startgpos = rs2gpos[index2rs[core]];
  int startpos = rs2pos[index2rs[core]];
  double startehh = ehh(core, core);
  int i;
  if (right){
    for (i = core; i<=stop; i++){
      double nextgpos = rs2gpos[index2rs[i]];
      int nextpos = rs2pos[index2rs[i]];
      string nextrs = index2rs[i];
      double gdist = nextgpos-startgpos;
      int dist = nextpos-startpos;

      //
      // if gap is too big, return with error
      //

      if (dist>200000 and distcorrect){
		//this should never happen because of the checks for finding bounds for integration
		return(-1);
	}
      if (gdist<0){
	cerr << "ERROR: genetic distance between "<< index2rs[i-1] << " and "<< index2rs[i] << " is negative ("<< gdist << ")\n";
	exit(1);
      }
      double nextehh = ehh(core, i);
      //
      // if gap is large, correct
      //

      if (dist>20000 and distcorrect){
    	  double cfactor = 20000/(double)dist;
    	  gdist = gdist*cfactor;
      }
      double add = gdist*0.5*(startehh+nextehh);
      //cout << nextgpos << " "<< gdist<< " "<< startehh<< " "<< nextehh << " "<< add << "\n";
      total+=add;
      startgpos = nextgpos;
      startpos = nextpos;
      startehh = nextehh;
    }
  }
  else{
    for (i = core; i>=stop; i--){
      double nextgpos = rs2gpos[index2rs[i]];
      int nextpos = rs2pos[index2rs[i]];
      string nextrs = index2rs[i];
      double gdist = startgpos-nextgpos;
      int dist = startpos-nextpos;

      //
      // if gap is too big, return with error
      //

      if (dist>200000 and distcorrect){
	// this should never happen
		return(-1);
	}
      if (gdist<0){
	cerr << "ERROR: genetic distance between "<< index2rs[i] << " and "<< index2rs[i+1] << " is negative ("<< gdist << ")\n";
        exit(1);
	}
      double nextehh = ehh( i, core);
      //
      // if gap is large, correct
      //

      if (dist>20000 and distcorrect){
	double cfactor = 20000/(double)dist;
	gdist = gdist*cfactor;
      }

      double add = gdist*0.5*(startehh+nextehh);
      total+=add;
      //cout << nextgpos <<  " "<< gdist << " "<< startehh<< " "<< nextehh << " "<< add << "\n";
      //cout <<nextrs<< " "<<nextgpos << " " << startgpos << " " << startehh << " " << nextehh << " " << add  <<  " " << total<<"\n";
      startgpos = nextgpos;
      startpos = nextpos;
      startehh = nextehh;
    }
  }
  return(total);
}


int Ehhpop::findcutoff(int core, double thold, bool right, bool distcorrect){
  string corers = index2rs[core];
  int corepos = rs2pos[corers];
  int stoppos = 0;
  if (right){
    int test = core;
    int prevpos = corepos;
    while(test < Nsnp-1){
      test++;
      int testpos = rs2pos[index2rs[test]];
      if (testpos-corepos>4000000){
    	  stoppos = test;
    	  break;
      }
      if (testpos-prevpos >200000 and distcorrect){
    	  stoppos = -1;
    	  break;
      }
      double testehh = ehh(core, test);
      //cout << testehh << "\n";
      if (testehh<thold){
    	  stoppos = test;
    	  break;
      }
      prevpos = testpos;
    }
    if (test ==Nsnp-1){
      stoppos = Nsnp-1;
    }
  }
  else{
    int test = core;
    int prevpos = corepos;
    while(test >0){
      test--;
      int testpos = rs2pos[index2rs[test]];
      if (corepos-testpos>4000000){
	//cout << corepos << " "<< testpos << "\n";
    	  stoppos = test;
    	  break;
      }
      if (prevpos-testpos>200000 and distcorrect){
    	  stoppos = -1;
    	  break;
      }
      //   cout << core << " " << test << "\n";
      double testehh = ehh(test, core);
      if (testehh<thold){
    	  stoppos = test;
    	  break;
      }
      prevpos = testpos;
    }
    if (test ==0){
      stoppos = 0;
    }
  }
  //cout << stoppos << "\n";
  return(stoppos);
}


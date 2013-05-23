#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <map>
#include <sstream>
using namespace std;


class Ehhpop{
 public:
  Ehhpop(const char[], const char[]);
  Ehhpop(Ehhpop&, Ehhpop&)throw();
  void readhap(const char[]);
  void readmap(const char[]);
  const Ehhpop &operator=(const Ehhpop&);
  void printhap();
  void printmap();
  double ehh (vector< vector<char> >);
  double ehh (int, int);
  int findcutoff(int, double, bool, bool);
  double integrate_ehh(int, int, bool, bool);
  int Nhap, Nsnp;
  vector< vector<char> > haps;
  map< string, double> rs2gpos;
  map<int, string> index2rs;
  map<string, int> rs2index;
  map<string, int> rs2pos;
  map<string, bool> rs2anc; // is the 0 allele ancestral?
};

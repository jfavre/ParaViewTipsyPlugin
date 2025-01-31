
#include <iostream>
#include "tipsy_file.h"

int main(int argc, char * argv[])
{
  TipsyFile *filein = new TipsyFile(argv[1]);
  filein->read_all();
  filein->write(argv[2], true);
  return 1;
}

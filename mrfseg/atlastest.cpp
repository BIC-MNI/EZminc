#include "atlasspec.h"

int main(int argc,char** argv)
{
  AtlasSpec atlas;
  readAtlasSpec(&atlas,argv[1]);
  cout << atlas.regionnames[3] << endl;
  cout << atlas.filenames[6] << endl;
  cout << atlas.labelTypes[4].pureLabel  << endl; 
  cout << atlas.labelTypes[4].mixed[0]  << endl; 
  cout << atlas.labelTypes[4].mixed[1]  << endl; 
  return(0);

}

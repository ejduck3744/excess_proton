#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class Shift;
StChain *chain;
void readPicoDst(const Char_t *inputFile = "7p7_list_trunc.list",  Char_t *outputFile = "test_mine")
{
  ////Int_t nEvents = 999999999;//10;//483;
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

	gSystem->Load("StRefMultCorr");
    gSystem->Load("StPicoEvent");
    gSystem->Load("StPicoDstMaker");
    gSystem->Load("StEpdUtil");
    gSystem->Load("StPileupUtil");
    gSystem->Load("Shift");

    chain = new StChain();

    StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");
    Shift *anaMaker = new Shift("ana",picoMaker, outputFile);

    chain->Init();
    ////cout << "chain->Init();" << endl;
    ////int total = picoMaker->chain()->GetEntries();
	Int_t nEvents = picoMaker->chain()->GetEntries();
	
    ////cout << " Total entries = " << total << endl;
	cout << " Total entries = " << nEvents << endl;
	////if (nEvents > total) nEvents = total;
	
    for (Int_t i = 0; i < nEvents; i++)
    {
        if (i % 1000 == 0){cout << "Working on eventNumber " << i << " of " << nEvents << endl;}

        chain->Clear();
        int iret = chain->Make(i);
        if (iret)
	  {
            cout << "Bad return code!" << iret << endl;
            break;
	  }
        ////total++;
    }
    cout << "****************************************** " << endl;
    cout << "Work done... now its time to close up shop!" << endl;
    cout << "****************************************** " << endl;
    chain->Finish();
    cout << "****************************************** " << endl;
    cout << "total number of events  " << nEvents << endl;
    cout << "****************************************** " << endl;
    delete chain;
}

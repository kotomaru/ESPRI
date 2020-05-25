/*macro to draw BDC tracking for RIBF79 2019*/
/*from HIMAC 2015 online macro by terashima*/

#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TBox.h"

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"
#include "TArtCalibBDC.hh"
#include "TArtBDC.hh"


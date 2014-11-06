#ifndef DmpAlgBgoRawTrack_H
#define DmpAlgBgoRawTrack_H

#include "DmpVAlg.h"
#include "TH1D.h"
#include "DmpEvtBgoHits.h"

class DmpEvtBgoHits;
class DmpAlgBgoRawTrack : public DmpVAlg{
/*
 *  DmpAlgBgoRawTrack
 *
 */
public:
  DmpAlgBgoRawTrack();
  ~DmpAlgBgoRawTrack();

  //void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()

  void Reset();
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool Finalize();
private:
  DmpEvtBgoHits *fBgoHits;
  TH1D *Angle_xz;
  TH1D *Angle_yz; 
  double x_lc[14];
  double Ex_lc[14];

  double Energy;
  double LC_E[14];//Cluster energy

};

#endif

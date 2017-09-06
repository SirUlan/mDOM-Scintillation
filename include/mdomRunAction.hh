#ifndef mdomRunAction_h
#define mdomRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class mdomRunAction : public G4UserRunAction
{
public:
  mdomRunAction();
  ~mdomRunAction();

public:
	void BeginOfRunAction(const G4Run*);
	void EndOfRunAction(const G4Run*);

private:
  double startingtime;
};

#endif

export ROOT_INCLUDE_PATH=path-to-json-include // if not in O2Physics environment

root -l -x -b -q "HFInvMassFitter.cxx" "runMassFitter.C(\"config_massfitter.json\")"


Step 1: Generate ROOT dictionary:
rootcling -f G__HFInvMassFitter.cxx -c ../HFInvMassFitter.h ../HFInvMassFitterLinkDef.h

Step 2: Compile source code:
g++ -fPIC -I$(root-config --incdir) -I path-to-json-include -c ../HFInvMassFitter.cxx ../runMassFitter.C G__HFInvMassFitter.cxx

Step 3: Link the executable:
g++ -o runMassFitter HFInvMassFitter.o runMassFitter.o G__HFInvMassFitter.o $(root-config --libs)  -lRooFit -lRooFitCore -lEG

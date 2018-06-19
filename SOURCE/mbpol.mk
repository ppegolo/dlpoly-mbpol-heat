mbpol-2b-poly.o: mbpol-2b-poly.cpp poly-2b-v6x.h
	$(CXX) -c -o mbpol-2b-poly.o $(CXXFLAGS) mbpol-2b-poly.cpp
	
poly-2b-v6x.o: poly-2b-v6x.cpp poly-2b-v6x.h
	$(CXX) -c -o poly-2b-v6x.o $(CXXFLAGS) poly-2b-v6x.cpp
	
mbpol-3b-poly.o: mbpol-3b-poly.cpp poly-3b-v2x.h
	$(CXX) -c -o mbpol-3b-poly.o $(CXXFLAGS) mbpol-3b-poly.cpp
	
poly-3b-v2x.o: poly-3b-v2x.cpp poly-3b-v2x.h
	$(CXX) -c -o poly-3b-v2x.o $(CXXFLAGS) poly-3b-v2x.cpp
	
mbnrg-2b-h2o-f-poly.o: mbnrg-2b-h2o-f-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-f-poly.o $(CXXFLAGS) mbnrg-2b-h2o-f-poly.cpp

poly-2b-h2o-ion-v1x.o: poly-2b-h2o-ion-v1x.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o poly-2b-h2o-ion-v1x.o $(CXXFLAGS) poly-2b-h2o-ion-v1x.cpp

poly-2b-A1-A1-v1x.o: poly-2b-A1-A1-v1x.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o poly-2b-A1-A1-v1x.o $(CXXFLAGS) poly-2b-A1-A1-v1x.cpp

poly-2b-A1-B1-v1x.o: poly-2b-A1-B1-v1x.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o poly-2b-A1-B1-v1x.o $(CXXFLAGS) poly-2b-A1-B1-v1x.cpp

mbnrg-2b-h2o-cl-poly.o: mbnrg-2b-h2o-cl-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-cl-poly.o $(CXXFLAGS) mbnrg-2b-h2o-cl-poly.cpp

mbnrg-2b-h2o-br-poly.o: mbnrg-2b-h2o-br-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-br-poly.o $(CXXFLAGS) mbnrg-2b-h2o-br-poly.cpp

mbnrg-2b-h2o-i-poly.o: mbnrg-2b-h2o-i-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-i-poly.o $(CXXFLAGS) mbnrg-2b-h2o-i-poly.cpp

mbnrg-2b-h2o-li-poly.o: mbnrg-2b-h2o-li-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-li-poly.o $(CXXFLAGS) mbnrg-2b-h2o-li-poly.cpp

mbnrg-2b-h2o-na-poly.o: mbnrg-2b-h2o-na-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-na-poly.o $(CXXFLAGS) mbnrg-2b-h2o-na-poly.cpp

mbnrg-2b-h2o-k-poly.o: mbnrg-2b-h2o-k-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-k-poly.o $(CXXFLAGS) mbnrg-2b-h2o-k-poly.cpp

mbnrg-2b-h2o-rb-poly.o: mbnrg-2b-h2o-rb-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-rb-poly.o $(CXXFLAGS) mbnrg-2b-h2o-rb-poly.cpp

mbnrg-2b-h2o-cs-poly.o: mbnrg-2b-h2o-cs-poly.cpp poly-2b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-2b-h2o-cs-poly.o $(CXXFLAGS) mbnrg-2b-h2o-cs-poly.cpp

mbnrg-2b-f-f-poly.o: mbnrg-2b-f-f-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-f-f-poly.o $(CXXFLAGS) mbnrg-2b-f-f-poly.cpp

mbnrg-2b-cl-cl-poly.o: mbnrg-2b-cl-cl-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-cl-cl-poly.o $(CXXFLAGS) mbnrg-2b-cl-cl-poly.cpp

mbnrg-2b-br-br-poly.o: mbnrg-2b-br-br-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-br-br-poly.o $(CXXFLAGS) mbnrg-2b-br-br-poly.cpp

mbnrg-2b-i-i-poly.o: mbnrg-2b-i-i-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-i-i-poly.o $(CXXFLAGS) mbnrg-2b-i-i-poly.cpp

mbnrg-2b-li-li-poly.o: mbnrg-2b-li-li-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-li-li-poly.o $(CXXFLAGS) mbnrg-2b-li-li-poly.cpp

mbnrg-2b-na-na-poly.o: mbnrg-2b-na-na-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-na-na-poly.o $(CXXFLAGS) mbnrg-2b-na-na-poly.cpp

mbnrg-2b-k-k-poly.o: mbnrg-2b-k-k-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-k-k-poly.o $(CXXFLAGS) mbnrg-2b-k-k-poly.cpp

mbnrg-2b-rb-rb-poly.o: mbnrg-2b-rb-rb-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-rb-rb-poly.o $(CXXFLAGS) mbnrg-2b-rb-rb-poly.cpp

mbnrg-2b-cs-cs-poly.o: mbnrg-2b-cs-cs-poly.cpp poly-2b-A1-A1-v1x.h
	$(CXX) -c -o mbnrg-2b-cs-cs-poly.o $(CXXFLAGS) mbnrg-2b-cs-cs-poly.cpp

mbnrg-2b-li-f-poly.o: mbnrg-2b-li-f-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-li-f-poly.o $(CXXFLAGS) mbnrg-2b-li-f-poly.cpp

mbnrg-2b-li-cl-poly.o: mbnrg-2b-li-cl-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-li-cl-poly.o $(CXXFLAGS) mbnrg-2b-li-cl-poly.cpp

mbnrg-2b-li-br-poly.o: mbnrg-2b-li-br-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-li-br-poly.o $(CXXFLAGS) mbnrg-2b-li-br-poly.cpp

mbnrg-2b-li-i-poly.o: mbnrg-2b-li-i-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-li-i-poly.o $(CXXFLAGS) mbnrg-2b-li-i-poly.cpp

mbnrg-2b-na-f-poly.o: mbnrg-2b-na-f-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-na-f-poly.o $(CXXFLAGS) mbnrg-2b-na-f-poly.cpp

mbnrg-2b-na-cl-poly.o: mbnrg-2b-na-cl-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-na-cl-poly.o $(CXXFLAGS) mbnrg-2b-na-cl-poly.cpp

mbnrg-2b-na-br-poly.o: mbnrg-2b-na-br-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-na-br-poly.o $(CXXFLAGS) mbnrg-2b-na-br-poly.cpp

mbnrg-2b-na-i-poly.o: mbnrg-2b-na-i-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-na-i-poly.o $(CXXFLAGS) mbnrg-2b-na-i-poly.cpp

mbnrg-2b-k-f-poly.o: mbnrg-2b-k-f-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-k-f-poly.o $(CXXFLAGS) mbnrg-2b-k-f-poly.cpp

mbnrg-2b-k-cl-poly.o: mbnrg-2b-k-cl-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-k-cl-poly.o $(CXXFLAGS) mbnrg-2b-k-cl-poly.cpp

mbnrg-2b-k-br-poly.o: mbnrg-2b-k-br-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-k-br-poly.o $(CXXFLAGS) mbnrg-2b-k-br-poly.cpp

mbnrg-2b-k-i-poly.o: mbnrg-2b-k-i-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-k-i-poly.o $(CXXFLAGS) mbnrg-2b-k-i-poly.cpp

mbnrg-2b-rb-f-poly.o: mbnrg-2b-rb-f-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-rb-f-poly.o $(CXXFLAGS) mbnrg-2b-rb-f-poly.cpp

mbnrg-2b-rb-cl-poly.o: mbnrg-2b-rb-cl-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-rb-cl-poly.o $(CXXFLAGS) mbnrg-2b-rb-cl-poly.cpp

mbnrg-2b-rb-br-poly.o: mbnrg-2b-rb-br-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-rb-br-poly.o $(CXXFLAGS) mbnrg-2b-rb-br-poly.cpp

mbnrg-2b-rb-i-poly.o: mbnrg-2b-rb-i-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-rb-i-poly.o $(CXXFLAGS) mbnrg-2b-rb-i-poly.cpp

mbnrg-2b-cs-f-poly.o: mbnrg-2b-cs-f-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-cs-f-poly.o $(CXXFLAGS) mbnrg-2b-cs-f-poly.cpp

mbnrg-2b-cs-cl-poly.o: mbnrg-2b-cs-cl-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-cs-cl-poly.o $(CXXFLAGS) mbnrg-2b-cs-cl-poly.cpp

mbnrg-2b-cs-br-poly.o: mbnrg-2b-cs-br-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-cs-br-poly.o $(CXXFLAGS) mbnrg-2b-cs-br-poly.cpp

mbnrg-2b-cs-i-poly.o: mbnrg-2b-cs-i-poly.cpp poly-2b-A1-B1-v1x.h
	$(CXX) -c -o mbnrg-2b-cs-i-poly.o $(CXXFLAGS) mbnrg-2b-cs-i-poly.cpp

poly-3b-h2o-ion-v1x.o: poly-3b-h2o-ion-v1x.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o poly-3b-h2o-ion-v1x.o $(CXXFLAGS) poly-3b-h2o-ion-v1x.cpp

mbnrg-3b-h2o-h2o-f-poly.o: mbnrg-3b-h2o-h2o-f-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-f-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-f-poly.cpp

mbnrg-3b-h2o-h2o-cl-poly.o: mbnrg-3b-h2o-h2o-cl-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-cl-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-cl-poly.cpp

mbnrg-3b-h2o-h2o-br-poly.o: mbnrg-3b-h2o-h2o-br-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-br-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-br-poly.cpp

mbnrg-3b-h2o-h2o-i-poly.o: mbnrg-3b-h2o-h2o-i-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-i-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-i-poly.cpp

mbnrg-3b-h2o-h2o-li-poly.o: mbnrg-3b-h2o-h2o-li-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-li-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-li-poly.cpp

mbnrg-3b-h2o-h2o-na-poly.o: mbnrg-3b-h2o-h2o-na-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-na-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-na-poly.cpp

mbnrg-3b-h2o-h2o-k-poly.o: mbnrg-3b-h2o-h2o-k-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-k-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-k-poly.cpp

mbnrg-3b-h2o-h2o-rb-poly.o: mbnrg-3b-h2o-h2o-rb-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-rb-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-rb-poly.cpp

mbnrg-3b-h2o-h2o-cs-poly.o: mbnrg-3b-h2o-h2o-cs-poly.cpp poly-3b-h2o-ion-v1x.h
	$(CXX) -c -o mbnrg-3b-h2o-h2o-cs-poly.o $(CXXFLAGS) mbnrg-3b-h2o-h2o-cs-poly.cpp










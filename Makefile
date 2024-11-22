prog: Lax_Fridrihs.o
	g++ Lax_Fridrihs.o -o prog.exe && del Lax_Fridrihs.o && prog.exe
Lax_Fridrihs.o: Lax_Fridrihs.cpp
	g++ -c Lax_Fridrihs.cpp
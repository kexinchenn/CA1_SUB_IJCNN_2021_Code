#include <getCombination.h>

void getCombination (int p, bool & lesionAV, bool & lesionHD, bool & lesionPos, bool & lesionLV, bool & lesionExc, bool & lesionInh) {
	if ( (p == 0) || (p >= 6 && p <= 10) || (p >= 21 && p <= 30) || (p >= 41 && p <= 50) || (p >= 56 && p <= 60) || (p == 62) ) {
		lesionAV = true;
	}
	else {
		lesionAV = false;
	}
	if ( (p == 1) || (p == 6) || (p >= 11 && p <= 14) || (p >= 21 && p <= 24) || (p >= 31 && p <= 36) || (p >= 41 && p<= 46) 
		|| (p >= 51 && p <= 54) || (p >= 56 && p <= 59) || (p == 61) || (p == 62) ) {
		lesionHD = true;
	}
	else {
		lesionHD = false;
	}
	if ( (p == 2) || (p == 7) || (p == 11) || (p >= 15 && p <= 17) || (p == 21) || (p >= 25 && p <= 27) || (p >= 31 && p <= 33) 
		|| (p >= 37 && p <= 39) || (p >= 41 && p <= 43) || (p >= 47 && p <= 49) || (p >= 51 && p <= 53) || (p >= 55 && p <= 58) 
		|| (p >= 60 && p <= 62) ) {
		lesionPos = true;
	}
	else {
		lesionPos = false;
	}
	if ( (p == 3) || (p == 8) || (p == 12) || (p == 15) || (p == 18) || (p == 19) || (p == 22) || (p == 25) || (p == 28) || 
		(p == 29) || (p == 31) || (p == 34) || (p == 35) || (p == 37) || (p == 38) || (p == 40) || (p == 41) || (p == 44) ||
		(p == 45) || (p == 47) || (p == 48) || (p >= 50 && p <= 52) || (p >= 54 && p <= 57) || (p >= 59 && p <= 62) ) {
		lesionLV = true;
	}
	else {
		lesionLV = false;
	}
	if ( (p == 4) || (p == 9) || (p == 13) || (p == 16) || (p == 18) || (p == 20) || (p == 23) || (p == 26) || (p == 28) 
		|| (p == 30) || (p == 32) || (p == 34) || (p == 36) || (p == 37) || (p == 39) || (p == 40) || (p == 42) || (p == 44)
		|| (p == 46) || (p == 47) || (p >= 49 && p <= 51) || (p >= 53 && p <= 56) || (p >= 58 && p <= 62) ) {
		lesionExc = true;
	}
	else {
		lesionExc = false;
	}
	if ( (p == 5) || (p == 10) || (p == 14) || (p == 17) || (p == 19) || (p == 20) || (p == 24) || (p == 27) || (p == 29) ||
		(p == 30) || (p == 33) || (p == 35) || (p == 36) || (p >= 38 && p <= 40) || (p == 43) || (p == 45) || (p == 46) ||
		(p >= 48 && p <= 50) || (p >= 52 && p <= 55) || (p >= 57 && p <= 62) ) {
		lesionInh = true;
	}
	else {
		lesionInh = false;
	}
}
#pragma once
#include "map.h"
#include <iostream>
#include <string>
#include "strlib.h"
#include "vector.h"
#include "error.h"
#include "set.h"

Map<std::string, Vector<std::string>> parseElementInfo(std::string filename);
Map<std::string, std::string> symbolsMap(std::string filename);
Map<std::string, Vector<std::string>> parseDiscoveryInfo(std::string filename);
Map<std::string, Set<std::string>> parseGroupInfo(std::string filename);
Map<std::string, std::string> parseElectronConfigInfo(std::string filename);

class PeriodicTable {
public:
    float calculateMass(std::string compound); // finds the mass of a compound
    std::string generalInfo(std::string element); // lists info associated with an element
    std::string mass(std::string element); // gives the mass of an element
    std::string number(std::string element); // gives the atomic number of an element
    std::string symbol(std::string element); // gives the chemical symbol of an element
    std::string name(std::string symbol); // gives the name of a chemical symbol
    std::string discoveryDate(std::string element); // gives the date an element was discovered
    std::string scientist(std::string element); // gives the scientist(s) who discovered the element
    Set<std::string> alkalis();
    Set<std::string> alkalines();
    Set<std::string> metals();
    Set<std::string> metalloids();
    Set<std::string> halogens();
    Set<std::string> noble_gases();
    Set<std::string> lanthanides();
    Set<std::string> actinides();
    Set<std::string> radioactives();
    std::string halfLife(std::string element);
    std::string electronConfig(std::string element);
    std::string longestCommonSubsequence(std::string s1, std::string s2);
    std::string speller(std::string word);
    std::string symbolSpeller(std::string word);
    Vector<std::string> spellerElements(std::string word);
    Set<std::string> allPerfectWords(std::string filename, int size);
    void client();
};


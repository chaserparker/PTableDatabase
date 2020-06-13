#include "periodictable.h"
#include "testing/SimpleTest.h"
#include <iostream>
#include "console.h"
#include <fstream>
#include "filelib.h"
#include "strlib.h"
#include <string>
#include "float.h"
#include <algorithm>
using namespace std;

// parses through the data in "mass_name_symbol_num.txt" to obtain information about each element
Map<string, Vector<string>> parseElementInfo(string filename) {
    Vector<string> lines;
    Map<string, Vector<string>> elementAndInfo;
    Vector<string> info;

    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);

    for (int i = 0; i < lines.size(); i++) {
        lines[i] = stringReplace(lines[i], " ", "\t");
        info = stringSplit(lines[i], "\t");
        for (int n = 0; n < info.size(); n++) {
            if (info[n].size() == 0) {
                info.remove(n);
            }
        }
        string mass = info[0];
        string name = info[1];
        string symbol = info[2];
        string number = info[3];
        elementAndInfo[toLowerCase(name)] = {symbol, number, mass};
    }
    return elementAndInfo;
}

// parses through the data in "mass_name_symbol_num.txt" to map chemical symbols to their associated names
Map<string, string> symbolsMap(string filename) {
    Vector<string> lines;
    Map<string, string> symbolsAndNames;
    Vector<string> info;

    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);

    for (int i = 0; i < lines.size(); i++) {
        lines[i] = stringReplace(lines[i], " ", "\t");
        info = stringSplit(lines[i], "\t");
        for (int n = 0; n < info.size(); n++) {
            if (info[n].size() == 0) {
                info.remove(n);
            }
        }
        string name = info[1];
        string symbol = info[2];
        symbolsAndNames[toLowerCase(symbol)] = name;
    }
    return symbolsAndNames;

}

// parses through the data in "discoverydates.txt" and obtains the discovery dates + scientists of each element
Map<string, Vector<string>> parseDiscoveryInfo(string filename) {
    Vector<string> lines;
    Map<string, Vector<string>> elementsAndInfo;
    Vector<string> info;
    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);

    for (int i = 0; i < lines.size(); i++) {
        string element;
        string date;
        string scientist;
        lines[i] = stringReplace(lines[i], " ", "\t");
        info = stringSplit(lines[i], "\t");
        element = info[1];
        if (i < 13) {
            date = info[2] + info[3];
        }
        else {
            date = info[2];
            scientist = info[3] + info[4];
            if (isalpha(info[5][0])) {
                int j = 0;
                while (isalpha(info[5 + j][0])) {
                    scientist += info[5 + j];
                    j++;
                }
            }
        }
        if (stringContains(scientist, "and")) {
            scientist = stringReplace(scientist, "and", ",");
        }
        if (stringContains(scientist, ".")) {
            scientist = stringReplace(scientist, ".", ". ");
        }
        if (stringContains(scientist, "etal.")) {
            scientist = stringReplace(scientist, "etal."," et al.");
        }
        if (stringContains(scientist, ",")) {
            scientist = stringReplace(scientist, ",", ", ");
        }
        elementsAndInfo[toLowerCase(element)] = {date, scientist};
    }
    return elementsAndInfo;
}

// parses through the data in "groups.txt" and categorizes atoms into their elemental groups
Map<string, Set<string>> parseGroupInfo(string filename) {
    Vector<string> lines;
    Map<string, Set<string>> allGroups;

    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);

    Set<string> group;
    for (int i = 1; i < 7; i++) {
        group.add(lines[i]);
    }
    allGroups["Alkali metals"] = group;
    group.clear();

    for (int i = 9; i < 15; i++) {
        group.add(lines[i]);
    }
    allGroups["Alkaline metals"] = group;
    group.clear();

    for (int i = 17; i < 66; i++) {
        group.add(lines[i]);
    }
    allGroups["Transition metals"] = group;
    group.clear();

    for (int i = 68; i < 75; i++) {
        group.add(lines[i]);
    }
    allGroups["Metalloids"] = group;
    group.clear();

    for (int i = 77; i < 83; i++) {
        group.add(lines[i]);
    }
    allGroups["Halogens"] = group;
    group.clear();

    for (int i = 85; i < 92; i++) {
        group.add(lines[i]);
    }
    allGroups["Noble gases"] = group;
    group.clear();

    for (int i = 94; i < 109; i++) {
        group.add(lines[i]);
    }
    allGroups["Lanthanides"] = group;
    group.clear();

    for (int i = 111; i < 126; i++) {
        group.add(lines[i]);
    }
    allGroups["Actinides"] = group;
    group.clear();

    for (int i = 128; i < 165; i++) {
        Vector<string> info;
        lines[i] = stringReplace(lines[i], " ", "\t");
        info = stringSplit(lines[i], "\t");
        for (int n = 0; n < info.size(); n++) {
            if (info[n].size() == 0) {
                info.remove(n);
            }
        }
        group.add(info[0]);
    }

    allGroups["Radioactive"] = group;
    group.clear();
    return allGroups;
}

Set<string> PeriodicTable::alkalis() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Alkali metals"];
    return group;
}

Set<string> PeriodicTable::alkalines() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Alkaline metals"];
    return group;
}

Set<string> PeriodicTable::metals() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Transition metals"];
    return group;
}

Set<string> PeriodicTable::metalloids() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Metalloids"];
    return group;
}

Set<string> PeriodicTable::halogens() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Halogens"];
    return group;
}

Set<string> PeriodicTable::noble_gases() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Noble gases"];
    return group;
}
Set<string> PeriodicTable::lanthanides() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Lanthanides"];
    return group;
}

Set<string> PeriodicTable::actinides() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Actinides"];
    return group;
}

Set<string> PeriodicTable::radioactives() {
    Map<string, Set<string>> info = parseGroupInfo("src/res/groups");
    Set<string> group = info["Radioactive"];
    return group;
}

// parses through the data in "electron-configs.txt" and obtains the electron configurations of every element
Map<string, string> parseElectronConfigInfo(std::string filename) {
    Vector<string> lines;
    Map<string, string> configurations;

    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);
    for (int i = 0; i < lines.size(); i++) {
        Vector<string> info = stringSplit(lines[i], ' ');
        info = stringSplit(lines[i], "\t");
        for (size_t j = 0; j < info[2].size(); j++) {
            if (isalpha(info[2][j]) && isnumber(info[2][j - 1])) {
                info[2] = info[2].insert(j - 1, "_");
                j++;
            }
        }
        configurations[toLowerCase(info[1])] = info[2];
    }
    return configurations;
}

// returns the half-life of an element
string PeriodicTable::halfLife(std::string element) {
    Vector<string> lines;

    ifstream in;
    openFile(in, "src/res/groups");
    readEntireFile(in, lines);

    Map<string, string> halfLives;
    for (int i = 128; i < 165; i++) {
        lines[i] = stringReplace(lines[i], " ", "\t");
        Vector<string> info = stringSplit(lines[i], "\t");
        info = stringSplit(lines[i], "\t");

        halfLives[toLowerCase(info[0])] = info[2] + ' ' + info[3];

    }
    if (element.size() < 3) {
        element = PeriodicTable::name(element);
    }
    if (!halfLives.containsKey(toLowerCase(element))) {
        return("Element is stable");
    }
    else {
        return("Half-life of " + toLowerCase(element) + ": " + halfLives[toLowerCase(element)]);
    }
}

// returns the electron configuration of an element
string PeriodicTable::electronConfig(std::string element) {
    string config;
    Map<string, string> configurations = parseElectronConfigInfo("src/res/electron-configs");
    if (!configurations.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    config = configurations[toLowerCase(element)];
    return("Electron configuration of " + toLowerCase(element) + ": " + config);
}

// returns all general info about an element
string PeriodicTable::generalInfo(std::string element) {
    Vector<string> allInfo;
    Map<string, Vector<string>> elementAndInfo = parseElementInfo("src/res/mass_name_symbol_num");
    if (!elementAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info1 = elementAndInfo[toLowerCase(element)];
    Map<string, Vector<string>> elementsAndInfo = parseDiscoveryInfo("src/res/discoverydates");
    Vector<string> info2 = elementsAndInfo[toLowerCase(element)];
    string text = "Name: " + toLowerCase(element) + "\nSymbol: " + info1[0] + "\nAtomic number: " + info1[1] + "\nMass: " + info1[2] + " amu" + "\nYear discovered: " + info2[0] + "\nDiscoverer(s): " + info2[1];
    return text;
}

// returns the mass of an element
string PeriodicTable::mass(std::string element) {
    Map<string, Vector<string>> elementAndInfo = parseElementInfo("src/res/mass_name_symbol_num");
    if (!elementAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info = elementAndInfo[toLowerCase(element)];
    return (info[2]);
}

// returns the atomic number of an element
string PeriodicTable::number(std::string element) {
    Map<string, Vector<string>> elementAndInfo = parseElementInfo("src/res/mass_name_symbol_num");
    if (!elementAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info = elementAndInfo[toLowerCase(element)];
    return (info[1]);
}

// returns the chemical symbol of an element
string PeriodicTable::symbol(std::string element) {
    Map<string, Vector<string>> elementAndInfo = parseElementInfo("src/res/mass_name_symbol_num");
    if (!elementAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info = elementAndInfo[toLowerCase(element)];
    return (info[0]);
}

// returns the name of a chemical symbol
string PeriodicTable::name(std::string symbol) {
    Map<string, string> symbolsAndNames = symbolsMap("src/res/mass_name_symbol_num");
    if (!symbolsAndNames.containsKey(toLowerCase(symbol))) {
        return(symbol + " is not an element!");
    }
    return (symbolsAndNames[toLowerCase(symbol)]);
}

// returns the date an element was discovered
std::string discoveryDate(std::string element) {
    Map<string, Vector<string>> elementsAndInfo = parseDiscoveryInfo("src/res/discoverydates");
    if (!elementsAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info = elementsAndInfo[toLowerCase(element)];
    return (info[0]);
}

// returns who discovered an element
std::string scientist(std::string element) {
    Map<string, Vector<string>> elementsAndInfo = parseDiscoveryInfo("src/res/discoverydates");
    if (!elementsAndInfo.containsKey(toLowerCase(element))) {
        return(element + " is not an element!");
    }
    Vector<string> info = elementsAndInfo[toLowerCase(element)];
    return (info[1]);
}


// calculates the atomic mass of a compound
float PeriodicTable::calculateMass(std::string compound) {
    if (compound.size() == 1) {
        string element = PeriodicTable::name(compound);
        string total = PeriodicTable::mass(element);
        return std::stof(total);
    }
    float total = 0;
    for (size_t i = 0; i < compound.size(); i++) {
        char ch = compound[i];
        if (i < compound.size() - 1) {
            if (islower(compound[i + 1])) {
                string symbol;
                symbol = ch;
                symbol += compound[i + 1];
                string element = PeriodicTable::name(symbol);
                string mass = PeriodicTable::mass(element);
                if (i + 2 <= compound.size() && isnumber(compound[i + 2])) {
                    int j = 0;
                    while (isnumber(compound[i + 2 + j]) && i + 2 + j < compound.size()) {
                        j++;
                    }
                    float multiMass = std::stof(mass) * std::stof(compound.substr(i + 2, j + 1));
                    mass = std::to_string(multiMass);
                    total += std::stof(mass);
                    i += j + 2;
                }
                else {
                    total += std::stof(mass);
                    i++;
                }
            }
            else if (isnumber(compound[i + 1])) {
                string symbol;
                symbol = ch;
                string element = PeriodicTable::name(symbol);
                string mass = PeriodicTable::mass(element);
                int j = 0;
                while (compound[i + 1 + j] != '_' && i + 1 + j < compound.size()) {
                    j++;
                }
                float multiMass = std::stof(mass) * std::stof(compound.substr(i + 1, j + 1));
                mass = std::to_string(multiMass);
                total += std::stof(mass);
                i += j + 1;
            }
            else {
                string symbol;
                symbol = ch;
                string element = PeriodicTable::name(symbol);
                string mass = PeriodicTable::mass(element);
                total += std::stof(mass);
            }
        }
    }
    if (isalpha(compound[compound.size() - 1]) && isupper(compound[compound.size() - 1])) {
        string symbol;
        symbol = compound[compound.size() - 1];
        string element = PeriodicTable::name(symbol);
        string mass = PeriodicTable::mass(element);
        total += std::stof(mass);
    }
    return total;
}


// takes two words and returns the longest common subsequence (where order is preserved)
string PeriodicTable::longestCommonSubsequence(string s1, string s2) { // borrowed from stanford cs106b curriculum
    if (s1.length() == 0 || s2.length() == 0) {
        return "";
        } else if (_tolower(s1[0]) == _tolower(s2[0])) {
            return s1[0] + longestCommonSubsequence(s1.substr(1), s2.substr(1));

        } else {
            string choice1 = longestCommonSubsequence(s1, s2.substr(1));
            string choice2 = longestCommonSubsequence(s1.substr(1), s2);
            if (choice1.length() >= choice2.length()) {
                return choice1;
            } else {
                return choice2;
            }
        }
    }


// Takes in a single word and returns its spelling using element symbols
std::string PeriodicTable::speller(std::string word) {
    Map<string, string> symbols = symbolsMap("src/res/mass_name_symbol_num");
    Vector<string> elements;
    for (size_t i = 0; i < word.size(); i++) {
        if (i < word.size() - 1) { // Puts all letter combos into a vector (singles and adjacent pairs)

            elements.add(charToString(toupper(word[i])));
            elements.add(charToString(toupper(word[i])) + charToString(tolower(word[i + 1])));
        }
        else {
            elements.add(charToString(toupper(word[i])));
        }
    }

    for (int i = 0; i < elements.size(); i++) { // removes all non-elements
        string element = elements[i];
        if (!symbols.containsKey(toLowerCase(element))) {
            elements.remove(i);
            i--;
        }

    }
    string elementsString;
    for (string elem : elements) {
        elementsString += elem;
    }
    string longest = PeriodicTable::longestCommonSubsequence(elementsString, word);

    for (size_t i = 0; i < longest.size(); i++) {
        string next = charToString(longest[i]) + longest[i + 1];
        string prev = charToString(longest[i]) + longest[i - 1];
        if (i < longest.size() - 1 && islower(longest[i]) && islower(longest[i - 1])) {
            if (!symbols.containsKey(charToString(longest[i])) && !symbols.containsKey(toLowerCase(next)) && !symbols.containsKey(toLowerCase(prev))) {
                longest.erase(i, 1);
                i--;
            }
            else if (symbols.containsKey(charToString(longest[i])) && !symbols.containsKey(toLowerCase(next))) {
                longest[i] = toupper(longest[i]);
            }
            else if (symbols.containsKey(toLowerCase(next))) {
                longest[i] = toupper(longest[i]);
             }
            else if (symbols.containsKey(charToString(longest[longest.size() - 1]))) {
                longest[longest.size() - 1] = toupper(longest[longest.size() - 1]);
            }
            else if (!symbols.containsKey(charToString(longest[i])) && symbols.containsKey(toLowerCase(prev))) {
                longest[i - 1] = toupper(longest[i - 1]);
            }
            else if (isupper(longest[i]) && isupper(longest[i + 1]) && !symbols.containsKey(charToString(tolower(longest[i])))) {
                longest.erase(i, 1);
                i--;
            }

        }
        if (islower(longest[i]) && islower(longest[i + 1]) && isupper(longest[i - 1]) && symbols.containsKey(charToString(tolower(longest[i - 1]))) && symbols.containsKey(next)) {
            longest[i] = toupper(longest[i]);
        }
        if (i == longest.size() - 1 && isupper(longest[i]) && !symbols.containsKey(charToString(tolower(longest[i])))) {
            longest.erase(i, 1);
            i--;
        }
        if (i == longest.size() - 1 && islower(longest[i]) && !symbols.containsKey(charToString(tolower(longest[i]))) && islower(longest[i - 1])) {
            longest.erase(i, 1);
            i--;
        }
        if (i == longest.size() - 1 && islower(longest[i]) && islower(longest[i - 1]) && symbols.containsKey(charToString(tolower(longest[i])))) {
            longest[i] = toupper(longest[i]);
        }
    }
    for (size_t i = 0; i < longest.size(); i++) {
        if (isupper(longest[i]) && isupper(longest[i + 1]) && !symbols.containsKey(charToString(tolower(longest[i])))) {
            longest.erase(i, 1);
            i--;
        }
        if (islower(longest[i]) && islower(longest[i - 1]) && !symbols.containsKey(charToString(tolower(longest[i])))) {
            longest.erase(i, 1);
            i--;
        }
    }
    return longest;
}

// calls speller() for a multi-word text
string PeriodicTable::symbolSpeller(string text) {
    string newText;
    for (size_t i = 0; i < text.size(); i++) {
        if (!isalpha(text[i]) && text[i] != ' ') {
            text.erase(i, 1);
            i--;

        }
    }
    Vector<string> words = stringSplit(text, ' ');
    for (string word : words) {
        word = PeriodicTable::speller(word);
        newText += word + ' ';
    }
    return newText;
}

// given a word that's spelled with element symbols, returns a vector containing the names of each element
Vector<string> PeriodicTable::spellerElements(string word) {
    Vector<string> elements;
    Map<string, string> symbols = symbolsMap("src/res/mass_name_symbol_num");
    for (size_t i = 0; i < word.size(); i++) {
        char ch = word[i];
        if (isupper(word[i]) && isupper(word[i + 1])) {
            elements.add(symbols[charToString(tolower(ch))]);
        }
        if (isupper(word[i]) && islower(word[i + 1])) {
            string elem;
            elem += ch;
            elem += word[i + 1];
            elements.add(symbols[toLowerCase(elem)]);
        }
        if (isupper(word[i]) && i == word.size() - 1) {
            elements.add(symbols[charToString(tolower(ch))]);
        }
    }
    return elements;
}

Set<string> PeriodicTable::allPerfectWords(string filename, int size) {
    if (size < 1) {
        return {};
    }
    Set<string> perfects;
    Vector<string> lines;

    ifstream in;
    openFile(in, filename);
    readEntireFile(in, lines);

    for (int i = 0; i < lines.size(); i++) {
        PeriodicTable PeriodicTable;
        string symbolSpell;
        if (lines[i].size() == size) {
            symbolSpell = PeriodicTable.speller(lines[i]);
            if (toLowerCase(lines[i]) == toLowerCase(symbolSpell)) {
                perfects.add(lines[i]);
            }
        }
    }
    return perfects;
}


//Set<string> size(string filename, int n) {
//    int count = 0;
//    Vector<string> lines;
//    Set<string> words;
//    ifstream in;
//    openFile(in, filename);
//    readEntireFile(in, lines);
//    for (int i = 0; i < lines.size(); i++) {
//        if (lines[i].size() == n) {
//            words.add(lines[i]);
//            count++;
//        }
//    }
//    return words;
//}

//Set<string> perfectsGivenWords(int n) {
//    Set<string> words = size("src/res/EnglishWords.txt", n);
//    string symbolSpell;
//    Set<string> perfects;
//    PeriodicTable ptable;
//    for (string word : words) {
//        symbolSpell = ptable.speller(word);
//        if (toLowerCase(word) == toLowerCase(symbolSpell)) {
//            cout << word << endl;
//            perfects.add(word);
//        }
//    }
//    return perfects;
//}


string printPerfects(int n) {
    if (n < 1) {
        return "{}";
    }
    if (n > 14) {
        return("Query too big (this fix is in progress)");
    }
    Vector<string> lines;

    ifstream in;
    openFile(in, "src/res/perfectwords");
    readEntireFile(in, lines);
    string allWords = lines[n - 1];
    Vector<string> numOfWords = {"2", "47", "361", "1336", "2656", "3556", "4427", "4574", "3517", "2468", "1776", "1252", "818", "481"};
    Vector<string> percentages = {"66.67", "48.96", "37.14", "34.23", "30.52", "23.35", "19.16", "16.09", "14.14", "12.16", "11.46", "11.02", "10.45", "9.38"};
    string returnString;
    returnString = allWords + "\n" + "\n" + "\nNumber of words: " + numOfWords[n - 1] + "\n" + percentages[n - 1] + " percent of " + integerToString(n) + "-letter words";
    return returnString;
}

//Map<string, int> freqOfElemements() {
//    Map<string, int> count;
//    Map<string, string> symbols = symbolsMap("src/res/mass_name_symbol_num");
//    for (string elem : symbols.keys()) {
//        count.keys().add(elem);
//        count[elem] = 0;
//    }
//    Vector<string> lines;
//    PeriodicTable ptable;

//    ifstream in;
//    openFile(in, "src/res/perfectwords");
//    readEntireFile(in, lines);
//    for (int i = 14; i < lines.size(); i++) {
//        lines[i].erase(0, 1);
//        lines[i].erase(lines[i].size() - 1);
//        Vector<string> words = stringSplit(lines[i], ' ');
//        for (string word : words) {
//            word = ptable.speller(word);
//            cout << word << endl;
//            for (size_t j = 0; j < word.size(); j++) {
//                string next = charToString(word[j]);
//                next += word[j + 1];
//                if (isupper(word[j])) {
//                    if (isupper(word[j + 1])) {
//                        count[charToString(tolower(word[j]))]++;
//                    }
//                    if (islower(word[j + 1])) {
//                        count[toLowerCase(next)]++;
//                    }
//                    if (j == word.size() - 1) {
//                        count[charToString(tolower(word[j]))]++;
//                    }

//                }
//            }
//        }
//    }
//    return count;
//}

Map<string, int> freqOfElemements() {
    Map<string, int> count;
    Map<string, string> symbols = symbolsMap("src/res/mass_name_symbol_num");
    for (string elem : symbols.keys()) {
        count.keys().add(elem);
        count[elem] = 0;
    }
    Vector<string> lines;
    ifstream in;
    openFile(in, "src/res/allperfects.txt");
    readEntireFile(in, lines);
    for (int i = 0; i < lines.size(); i++) {
        cout << i << endl;
        for (size_t j = 0; j < lines[i].size(); j++) {
            string word = lines[i];
            if (word[j] != '/') {
                string next = charToString(word[j]);
                next += word[j + 1];
                if (isupper(word[j])) {
                    if (word[j + 1] == '/') {
                        word.erase(j + 1, 1);
                    }
                    else if (isupper(word[j + 1])) {
                        count[charToString(tolower(word[j]))]++;
                    }
                    else if (islower(word[j + 1])) {
                        count[toLowerCase(next)]++;
                    }
                    else if (j == word.size() - 1) {
                        count[charToString(tolower(word[j]))]++;
                    }
                }
            }
        }
    }
    return count;
}

void printFreqs() {
    Map<std::string, float> perfectsFreq = {{"ac",998}, {"ag",443}, {"al",2387}, {"am",777}, {"ar",2217}, {"as",1408}, {"at",2009}, {"au",298},
                                          {"b",3440}, {"ba",260}, {"be",400}, {"bh",6}, {"bi",529}, {"bk",1}, {"br",279}, {"c",6633}, {"ca",568},
                                          {"cd",4}, {"ce",666}, {"cf",0}, {"cl",246}, {"cm",0}, {"cn",0}, {"co",1615}, {"cr",383}, {"cs",204},
                                          {"cu",592}, {"db",37}, {"ds",726}, {"dy",147}, {"er",5924}, {"es",6323}, {"eu",140}, {"f",2957}, {"fe",245},
                                          {"fl",229}, {"fm",0}, {"fr",129}, {"ga",420}, {"gd",8}, {"ge",791}, {"h",3708}, {"he",652}, {"hf",5},
                                          {"hg",12}, {"ho",960}, {"hs",110}, {"i",9052}, {"in",1591}, {"ir",328}, {"k",2467}, {"kr",24}, {"la",1954},
                                          {"li",1936}, {"lr",5}, {"lu",556}, {"lv",88}, {"mc",13}, {"md",1}, {"mg",0}, {"mn",83}, {"mo",860}, {"mt",8},
                                          {"n",7931}, {"na",356}, {"nb",42}, {"nd",948}, {"ne",1037}, {"nh",32}, {"ni",1164}, {"no",910}, {"np",31},
                                          {"o",8751}, {"og",361}, {"os",701}, {"p",5727}, {"pa",342}, {"pb",18}, {"pd",13}, {"pm",5}, {"po",637},
                                          {"pr",675}, {"pt",132}, {"pu",342}, {"ra",1850}, {"rb",98}, {"re",2990}, {"rf",70}, {"rg",148}, {"rh",58},
                                          {"rn",360}, {"ru",477}, {"s",9016}, {"sb",22}, {"sc",120}, {"se",1058}, {"sg",4}, {"si",905}, {"sm",573},
                                          {"sn",25}, {"sr",2}, {"ta",1451}, {"tb",77}, {"tc",238}, {"te",2222}, {"th",805}, {"ti",3213}, {"tl",306},
                                          {"tm",25}, {"ts",1396}, {"u",4658}, {"v",2277}, {"w",1771}, {"xe",138}, {"y",1528}, {"yb",26}, {"zn",1}, {"zr",0}};

    Vector<string> elementsSorted;
    Vector<string> keys = perfectsFreq.keys();
    Vector<float> vals = perfectsFreq.values();
    Vector<float> valsSorted = vals;
    Vector<float> valsPercent;
    valsSorted.sort();
    int total = 0;
    for (int i = 0; i < vals.size(); i++) {
        total += vals[i];

    }
    for (int i = 0; i < valsSorted.size(); i++) {
        for (int j = 0; j < vals.size(); j++) {
            if (valsSorted[i] == vals[j]) {
                string elem = keys[j];
                if (!elementsSorted.contains(elem)) {
                    elementsSorted.add(elem);
                }
            }
        }
    }
    for (string element : elementsSorted) {
        element[0] = toupper(element[0]);
    }
    for (int i = 0; i < elementsSorted.size(); i++) {
        elementsSorted[i][0] = toupper(elementsSorted[i][0]);
        float percent = valsSorted[i] / total * 100;
        valsPercent.add(percent);
        cout << elementsSorted[i] << ": " << valsPercent[i] << "%, ";
    }
//    cout << elementsSorted << endl;
//    cout << valsPercent << endl;
}

//////////////////////////////////////
/////// CLIENT IS RIGHT HERE   ///////
//////////////////////////////////////


void PeriodicTable::client() {
    ifstream in;
    string decision;
    cout << "Welcome to the Periodic Table Database!" << endl;
    cout << endl;
    cout << endl;
    cout << "Informational commands:" << endl;
    cout << "ALL (1a), SYMBOL (1b), NUMBER (1c), MASS (1d), NAME (1e)" << endl;
    cout << endl;
    cout << "Computational commands:" << endl;
    cout << "CALC-MASS (2a), ELECTRON-CONFIG (2b), HALF-LIFE (2c)" << endl;
    cout << endl;
    cout << "Organizational commands:" << endl;
    cout << "ALKALIS (3a), ALKALINES (3b), METALS (3c), METALLOIDS (3d), HALOGENS (3e), NOBLE-GASES (3f), LANTHANIDES (3g), ACTINIDES (3h), RADIOACTIVES (3i)" << endl;
    cout << endl;
    cout << "Miscellaneous commands:" << endl;
    cout << "CHEMSPELL (4a), SPELLMATCH (4.1a)" << endl;
    cout << endl;

    while (toLowerCase(decision) != "quit") {

        cout << "Type a command or type '?' for a list of commands: ";
        std::getline (std::cin,decision);

    if (decision == "?") {
        cout << "Informational commands:" << endl;
        cout << "ALL (1a), SYMBOL (1b), NUMBER (1c), MASS (1d), NAME (1e)" << endl;
        cout << endl;
        cout << "Computational commands:" << endl;
        cout << "CALC-MASS (2a), ELECTRON-CONFIG (2b), HALF-LIFE (2c)" << endl;
        cout << endl;
        cout << "Organizational commands:" << endl;
        cout << "ALKALIS (3a), ALKALINES (3b), METALS (3c), METALLOIDS (3d), HALOGENS (3e), NOBLE-GASES (3f), LANTHANIDES (3g), ACTINIDES (3h), RADIOACTIVES (3i)" << endl;
        cout << endl;
        cout << "Miscellaneous commands:" << endl;
        cout << "CHEMSPELL (4a)" << endl;
    }

    if (toLowerCase(decision) == "all" || toLowerCase(decision) == "1a") {
        string element;
        cout << endl;
        cout << "Enter the element you want to learn about: ";
        std::getline (std::cin,element);
        if (element.size() < 3) {
            element = PeriodicTable::name(element);
        }
        string all = PeriodicTable::generalInfo(element);
        cout << endl;
        cout << all << endl;
    }
    else if (toLowerCase(decision) == "symbol" || toLowerCase(decision) == "1b") {
        string element;
        cout << "Enter the element whose symbol you want: ";
        std::getline (std::cin,element);
        if (element.size() < 3) {
            element = PeriodicTable::name(element);
        }
        string symbol = PeriodicTable::symbol(element);
        cout << endl;
        cout << symbol << endl;
    }
    else if (toLowerCase(decision) == "number" || toLowerCase(decision) == "1c") {
        string element;
        cout << "Enter the element whose atomic number you want: ";
        std::getline (std::cin,element);
        if (element.size() < 3) {
            element = PeriodicTable::name(element);
        }
        string number = PeriodicTable::number(element);
        cout << endl;
        cout << number << endl;
    }
    else if (toLowerCase(decision) == "mass" || toLowerCase(decision) == "1d") {
        string element;
        cout << "Enter the element whose mass you want: ";
        std::getline (std::cin,element);
        if (element.size() < 3) {
            element = PeriodicTable::name(element);
        }
        string mass = PeriodicTable::mass(element);
        cout << endl;
        cout << mass << endl;
    }
    else if (toLowerCase(decision) == "name" || toLowerCase(decision) == "1e") {
        string symbol;
        cout << "Enter the element whose name you want: ";
        std::getline (std::cin,symbol);
        string name = PeriodicTable::name(symbol);
        cout << endl;
        cout << name << endl;
    }
    else if (toLowerCase(decision) == "calc-mass" || toLowerCase(decision) == "2a") {
        string compound;
        cout << endl;
        cout << "Examples of compound formatting: CO2, H2_O, C6_H12_O6, CaBr" << endl;
        cout << "WARNING: If format is incorrect, the program will crash" << endl;
        cout << endl;
        cout << "Enter the compound whose mass you want: ";
        std::getline (std::cin,compound);
        float calc_mass = PeriodicTable::calculateMass(compound);
        cout << endl;
        cout << calc_mass << " amu" << endl;
    }
    else if (toLowerCase(decision) == "electron-config" || toLowerCase(decision) == "2b") {
        string element;
        cout << endl;
        cout << "Enter the element whose electron configuration you want:" << endl;
        std::getline (std::cin,element);
        if (element.size() < 3) {
            element = PeriodicTable::name(element);
        }
        string config = PeriodicTable::electronConfig(element);
        cout << endl;
        cout << config << endl;
    }
    else if (toLowerCase(decision) == "half-life" || toLowerCase(decision) == "2c") {
        string element;
        cout << endl;
        cout << "Radioactive elements: " << endl;
        cout << PeriodicTable::radioactives() << endl;
        cout << endl;
        cout << "Enter the radioactive element whose half-life you want:" << endl;
        std::getline (std::cin,element);
        string halfLife = PeriodicTable::halfLife(element);
        cout << endl;
        cout << halfLife << endl;
    }
    else if (toLowerCase(decision) == "alkalis" || toLowerCase(decision) == "3a") {
        cout << endl;
        cout << PeriodicTable::alkalis() << endl;
    }
    else if (toLowerCase(decision) == "alkalines" || toLowerCase(decision) == "3b") {
        cout << endl;
        cout << PeriodicTable::alkalines() << endl;
    }
    else if (toLowerCase(decision) == "metals" || toLowerCase(decision) == "3c") {
        cout << endl;
        cout << PeriodicTable::metals() << endl;
    }
    else if (toLowerCase(decision) == "metalloids" || toLowerCase(decision) == "3d") {
        cout << endl;
        cout << PeriodicTable::metalloids() << endl;
    }
    else if (toLowerCase(decision) == "halogens" || toLowerCase(decision) == "3e") {
        cout << endl;
        cout << PeriodicTable::halogens() << endl;
    }
    else if (toLowerCase(decision) == "noble-gases" || toLowerCase(decision) == "3f") {
        cout << endl;
        cout << PeriodicTable::noble_gases() << endl;
    }
    else if (toLowerCase(decision) == "lanthanides" || toLowerCase(decision) == "3g") {
        cout << endl;
        cout << PeriodicTable::lanthanides() << endl;
    }
    else if (toLowerCase(decision) == "actinides" || toLowerCase(decision) == "3h") {
        cout << endl;
        cout << PeriodicTable::actinides() << endl;
    }
    else if (toLowerCase(decision) == "radioactives" || toLowerCase(decision) == "3i") {
        cout << endl;
        cout << PeriodicTable::radioactives() << endl;
    }
    else if (toLowerCase(decision) == "chemspell" || toLowerCase(decision) == "4a") {
        cout << endl;
        cout << "Welcome to Chemspell! Enter a message to see it spelled using the periodic table!" << endl;
        cout << "*******Keep individual words less than 15 letters (to remedy this, insert spaces in the word)" << endl;
        while (toLowerCase(decision) != "q") {
            string text;
            cout << endl;
            cout << "Message:" << endl;
            std::getline (std::cin,text);
            cout << endl;
            string newText = PeriodicTable::symbolSpeller(text);
            cout << newText << endl;
            cout << endl;
            cout << "Type 'SEE' to look at the elements involved" << endl;
            cout << "Type 'C' to do another entry" << endl;
            cout << "Type 'Q' to leave" << endl;
            cout << endl;
            std::getline (std::cin,decision);
            if (toLowerCase(decision) == "see") {
                Vector<string> words = stringSplit(newText, ' ');
                for (string word : words) {
                    cout << word << ": " << PeriodicTable::spellerElements(word) << endl;
                    cout << endl;
                }
                cout << endl;
                decision = "c";
            }
            if (toLowerCase(decision) != "c" && toLowerCase(decision) != "q") {
                while (toLowerCase(decision) != "c" && toLowerCase(decision) != "q") {
                    cout << "Invalid command" << endl;
                    std::getline (std::cin,decision);
                }
            }
        }
    }

    else if (toLowerCase(decision) == "spellmatch" || toLowerCase(decision) == "4.1a") {
        cout << endl;
        cout << "Welcome to Spellmatch!" << endl;
        cout << "This program shows every English word that can be made with Chemspell" << endl;
        cout << endl;
        while (toLowerCase(decision != "q")) {
            string strNum;
            cout << endl;
            cout << "Enter a number to see all words of that length that can be made with Chemspell" << endl;
            cout << endl;
            cout << "OR" << endl;
            cout << endl;
            cout << "Type 'F' to see the frequency of each element" << endl;
            std::getline (std::cin,strNum);
            if (toLowerCase(strNum) != "f") {
                int num = stringToInteger(strNum);
                string perfects = printPerfects(num);
                cout << perfects << endl;
                cout << endl;
            }
            if (toLowerCase(strNum) == "f") {
                printFreqs();
                cout << endl;
                cout << endl;
            }
            cout << "Type 'C' to do another entry" << endl;
            cout << "Type 'Q' to leave" << endl;
            cout << endl;
            std::getline (std::cin,decision);
            if (toLowerCase(decision) != "c" && toLowerCase(decision) != "q") {
                while (toLowerCase(decision) != "c" && toLowerCase(decision) != "q") {
                    cout << "Invalid command" << endl;
                    std::getline (std::cin,decision);
                }
            }
        }

    }

    else if (toLowerCase(decision) == "continue") {
        cout << endl;
    }
    else if (toLowerCase(decision) != "?"){
        cout << "Invalid command" << endl;
        cout << endl;
    }
    cout << endl;
    cout << endl;
    }

}


//PROVIDED_TEST("symbolSpeller() test") {
//    string word = "chaser";
//    PeriodicTable ptable;
//    string longest = ptable.symbolSpeller(word);
//    cout << "TEST CASE: " << longest << endl;
//}

//PROVIDED_TEST("count english words") {
//    string file = "src/res/EnglishWords.txt";
//    Set<string> a = size(file, 14);
//    cout << "SIZE: " << a << endl;
//}

//PROVIDED_TEST("perfectsGivenWords() test") {
//    Set<string> a = perfectsGivenWords(14);
//    cout << a << endl;
//    cout << "SIZE: " << a.size() << endl;
//}

//PROVIDED_TEST("printPerfects() test") {
//    cout << printPerfects(5) << endl;
//}

PROVIDED_TEST("freqOfElements() test") {
    cout << freqOfElemements() << endl;
}

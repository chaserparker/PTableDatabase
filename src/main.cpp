#include <iostream>
#include "console.h"
#include "periodictable.h"
#include "testing/SimpleTest.h"
using namespace std;

int main() {
//        if (runSimpleTests(SELECTED_TESTS)) {
//            return 0;
//        }

    PeriodicTable periodictable;
    periodictable.client();
    return 0;
}

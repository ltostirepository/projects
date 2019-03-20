#include <iostream>
#include "myclass.h"

using namespace std;

int main(int argc, char* argv[])
{
    // suppress warnings
    (void)argc;
    (void)argv;
    
    MyClass a; // no longer produces an error, because MyClass is defined
    a.print();
    
    a.value(5);
    a.Get_value();
    
    cout<< "\nHello World!" <<endl;
    return 0;
}

#include <iostream>
#include "myclass.h"

using namespace std;

void MyClass::foo(){
}

void MyClass::print(){
    cout<<"\nprinting....\n";
    cout<<"private value a is ->>>> "<< a <<endl;
}

void MyClass::value(int n){
    
    cout<<"\nmodifing private value from "<< a <<" to "<< n <<"\n";
    a=n;
}

void MyClass::Get_value(){
    print();
}

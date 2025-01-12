#include <stdio.h>

void testf(int* a) {
    *a = 1;
}

int main() {

    int a = 0; 
    testf(&a);

    printf("Value of a = %i\n", a);

    return 0;
}
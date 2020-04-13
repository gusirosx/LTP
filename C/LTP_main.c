#include <global.h>
int main() 
{
    readfile();
    if(sheme_flag) 
        Scheme_routine();
    if(rms_flag) 
        ERROR_routine();
    return(0);
}

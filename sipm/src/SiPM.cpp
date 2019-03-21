#include "general_lib.h"
#include "Peaks.h"
#include "functions_sipm_settings.h"
#include "functions_sipm_smoothing_and_DLED.h"
#include "functions_sipm_Peaks_charge.h"
#include "functions_sipm_Peaks_research.h"
#include "functions_sipm_print_and_trials.h"
#include "drawnolaser.h"

using namespace std;

int main(int argc, char* argv[])
{
    // suppress warnings
    (void)argc;
    (void)argv;
    
    cout<<"........INIZIALIZING PROGRAM........\n";
    drawnolaser("f_in_only_tgraphs.root",0, false,-999,-999,false ,8e-10);
    cout<< ".......COMPLETED RUN.......\n" <<endl;
    return 0;
}


omniarcanum@vbox:~/tau_generator$ root -l
root [0] .x generator_tautau.cpp
In file included from input_line_8:1:
/home/omniarcanum/tau_generator/generator_tautau.cpp:85:1: error: expected expression
@@ -111,40 +113,52 @@ void makeSample(Long64_t N, int spinMode, const char *outName)
^
/home/omniarcanum/tau_generator/generator_tautau.cpp:89:5: error: use of undeclared identifier 'tree'
    tree.Write();
    ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:90:5: error: use of undeclared identifier 'hCos'
    hCos.Write();
    ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:91:5: error: use of undeclared identifier 'fout'
    fout.Close();
    ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:92:30: error: use of undeclared identifier 'accepted'
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
                             ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:92:56: error: use of undeclared identifier 'outName'
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
                                                     ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:92:80: error: use of undeclared identifier 'spinMode'
    std::cout<<"Generated "<<accepted<<" events → "<<outName<<" (spinMode="<<spinMode<<")\n";
                                                                             ^
/home/omniarcanum/tau_generator/generator_tautau.cpp:93:1: error: extraneous closing brace ('}')
}

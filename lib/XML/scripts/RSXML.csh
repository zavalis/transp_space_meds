if(`test -n "-L/opt/R/arm64/lib -lxml2 -lz -L/opt/R/arm64/lib -llzma -lpthread -liconv -lm"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-L/opt/R/arm64/lib -lxml2 -lz -L/opt/R/arm64/lib -llzma -lpthread -liconv -lm
else
   setenv LD_LIBRARY_PATH -L/opt/R/arm64/lib -lxml2 -lz -L/opt/R/arm64/lib -llzma -lpthread -liconv -lm
endif

endif

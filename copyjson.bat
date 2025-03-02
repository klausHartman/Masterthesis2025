setlocal
set "source=C:\Users\User\Downloads\PerceptLFPs\Percept-Pat-LFP\PD071PH_sub-20210622PStn"
set "target=C:\Users\User\MATLAB Drive\masterthesis\Patient2"
for /r "%source%" %%f in (*.jsoncopy) do (
   echo kopiere %%f nach %target%
   copy "%%f" "%target%"
)
echo Alle dateien wurden kopiert!
endlocal

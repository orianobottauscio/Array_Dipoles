aus=version;
kk=strfind(aus,'(');
versioneMatlab = aus(kk+2:kk+6);

if strcmp(versioneMatlab(1:4),'2020')
   mcc -m Array_Dipoles.m -a dipoles_sub.m util.m monte_carlo.m read_labels.m toml.m -o Array_Dipoles_2020
else
   mcc -m Array_Dipoles.m -a dipoles_sub.m util.m monte_carlo.m read_labels.m toml.m -o Array_Dipoles
end
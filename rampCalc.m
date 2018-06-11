Vinitial = 250;
Vfinal = 300;
Cfinal = 8;
Nsteps = 8;
%Cstock = 50*10^3;

%------------

dC = Cfinal/Nsteps;
Ctreat = Cfinal*Vfinal/(Vfinal - Vinitial);

%------------ summary

disp('------------------------------------------');
disp('------------------------------------------');
disp('INPUT:');
disp(['initial volume: ' num2str(Vinitial) ' ul']);
disp(['final volume: ' num2str(Vfinal) ' ul']);
disp(['final concentration: ' num2str(Cfinal) ' ng/ml']);
%disp(['stock concentration: ' num2str(Cstock/10^3) ' ug/ml']);
disp(['number of steps: ' num2str(Nsteps)]);
disp('----------------------------------');
disp('RESULT:');
disp(['treatment concentration: ' num2str(Ctreat) ' ng/ml']);
%disp(['treatment dilution: ' num2str(1) ':' num2str(round(Cstock/Ctreat))]);

disp('--');
%------------

V = Vinitial;   % total volume
L = 0;          % total amount of ligand

for i = 1:Nsteps
   
    Ctarget = dC*i;
    
    dV = (L - Ctarget*V)/(Ctarget-Ctreat);
    disp(['treatment ' num2str(i) ' volume: ' num2str(dV,3)]);
    
    V = V + dV;
    L = L + dV*Ctreat;
end
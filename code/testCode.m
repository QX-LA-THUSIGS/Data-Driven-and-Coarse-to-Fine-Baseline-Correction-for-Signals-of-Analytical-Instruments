%**************************************************************************
%test program
%*************************************************************************

%% mass spectral data test
close all;

for i = 1:1:length(ms_data)
    
    Sig = ms_data(i).signal;
    [DBSig] = DD_CF_v1(Sig);
    figure(3+i);
    t1 = 1:1:length(Sig);t2 = 1:1:length(DBSig);
    zero_line = zeros(1,length(DBSig));
    plot(t1,Sig,t2,DBSig,'r',t2,zero_line,'k');xlabel('Variable');ylabel('Amp');
    legend('Original signal','Corrected signal');
    close Figure 1
    
end

close Figure 2
%% Chromatograph data test
close all;

for i = 1:1:length(chro_data)
    
    Sig = chro_data(i).signal;
    [DBSig] = DD_CF_v1(Sig);
    figure(3+i);
    t1 = 1:1:length(Sig);t2 = 1:1:length(DBSig);
    zero_line = zeros(1,length(DBSig));
    plot(t1,Sig,t2,DBSig,'r',t2,zero_line,'k');xlabel('Variable');ylabel('Amp');axis tight;
    legend('Original signal','Corrected signal');
    close Figure 1
    
end

close Figure 2

%% Migration spectrum data test
close all;

for i = 1:1:length(migration_spectrum_data)
    
    Sig = migration_spectrum_data(i).signal;
    [DBSig] = DD_CF_v1(Sig);
    figure(3+i);
    t1 = 1:1:length(Sig);t2 = 1:1:length(DBSig);
    zero_line = zeros(1,length(DBSig));
    plot(t1,Sig,t2,DBSig,'r',t2,zero_line,'k');xlabel('Variable');ylabel('Amp');axis tight;
    legend('Original signal','Corrected signal');
    close Figure 1
    
end

close Figure 2
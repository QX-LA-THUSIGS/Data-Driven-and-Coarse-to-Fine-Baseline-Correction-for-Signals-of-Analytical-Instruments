%*************************************************************************************************************************************************************
%DD-CF algorithm(Data-driven and coarse-to-fine baseline correction algorithm)
%
%Function: 
%You can use this program to realize the baseline adaptive correction of signals from a variety of analytical instruments, 
%including but not limited to mass spectrometers, ion mobility spectrometers, and chromatographs. 
%The feature of this program is that in the process of different kinds of chemical signal baseline calibration, human intervention is generally not required.
%
%The copyright belongs to the Analytical Instrument Group of Shenzhen International Graduate School of Tsinghua University.
%This program is designed and maintained by the QX-LA-THUSIGS team.
%
%References:(If you need to cite or improve this algorithm, please declare the original document of this algorithm in your work.)
%Xiangchun Xu, Xiang Qian et.al, Data-Driven and Coarse-to-Fine Baseline Correction for Signals of Analytical Instruments[J],
%Analytica Chimica Acta,2021,1157(338386):1873-4324.��DOI:10.1016/j.aca.2021.338386��
%
% How to use?
%You can use the DD_CF_v1 function "[DBSig,baseline] = DD_CF_v1(Sig)" to reslize the baseline adaptive
%correction.where, DBSig is the signal which represents signal after
%baseline correction;baseline represents the fitting baseline;Sig
%represents signal needed to be removed baseline.
% 
%*************************************************************************************************************************************************************


function [DBSig,baseline] = DD_CF_v1(Sig)
    %The FliplrLength is the proportion of the number of extreme points contained in the signal segment that needs to be flipped (percentage, this parameter does not need to be modified)
    FliplrLength = 100; %��Ҫ��ת�źŶΰ����ļ�ֵ�����ռ�ȣ��ٷֱȣ��ò�������Ҫ�޸ģ�

    SigButter = Sig;

    [tmax,tmin,SigButtermax,SigButtermin]  = fink_local_peaks(SigButter);   %Find all local maximum and local minimum points of the signal/�ҵ��ź����еľֲ�����ֵ����ֲ���Сֵ��

    %%%Step1:Location of high-amplitude peaks/�߷�ֵ�׷�Ķ�λ%%%
    xdata = tmax;
    ydata = SigButtermax;
    [peak_number,zcrnumber,gof,x,y,col] = polyfit_ployval_peaknum_zcrnumber(xdata,ydata,0.95,SigButter);    
  
    zcrnumber = zcrnumber+2;      %Add 2 to prevent extreme situations where the actual baseline inflection point appears at the end point/��2����ֹ����ʵ�ʻ��߹յ�����ڶ˵�ļ������

    %%%Step2:Truncate each peak, find the length of each segment, and judge
    %%%whether there is interference between them
    %%%��ÿ���׷���нضϣ����ÿһ�εĳ��ȣ����ж�֮���Ƿ�������%%%
    [subsection1,subsection_min] = delete_ms_peaks(col,tmax,tmin,SigButtermax,SigButtermin);
    
    %%%Cubic Spline Interpolation Method to Complement Spectral Peaks/����������ֵ�����׷�%%%
    [SigButterRIS] = cubic_splin_interpolation(subsection1,subsection_min,SigButtermax,SigButtermin);
    
    

    %%%Flip and extend both ends of the signal from which the spectral peak is removed/��ȥ����Ъ��ֵ���źŽ������˷�ת%%%
    [CompSigButter,Lsection] = reduce_end_effect(SigButterRIS,FliplrLength);

    %%%Step3:EMD/����ģ̬�ֽ�%%%
    N = length(CompSigButter);
    Nstd = 0.2;
    NR = 500;
    MaxIter = 5000;
    ecg = CompSigButter;
    % emd(ecg,'Interpolation','pchip', 'MaxNumIMF',30);
    [modes,res,info]=emd(ecg, 'MaxNumIMF',30);%,'Interpolation','pchip');
    modes = [modes res];
    modes = modes';
    t=1:length(ecg);

    %%%all IMF figure/�Ը���ģ̬���л�ͼ%%%
    [a,b]=size(modes);
    t = 1:1:length(CompSigButter);
    figure(2);
    subplot(a+1,1,1);
    plot(t,ecg,'color',[65 105 225]./255,'linewidth',1.6);        % the  signal is in the first row of the subplot
    ylabel('OriSig')
    set(gca,'xtick',[])
    axis tight;
    for i=2:a
        subplot(a+1,1,i);
        plot(t,modes(i-1,:),'color',[65 105 225]./255,'linewidth',1.6);
        ylabel (['IMF ' num2str(i-1)]);
        set(gca,'xtick',[])
        axis tight;
    end
    subplot(a+1,1,a+1)
    plot(t,modes(a,:),'color',[65 105 225]./255,'linewidth',1.6)
    ylabel(['IMF ' num2str(a)]);
    axis tight;xlabel('Variables','FontSize',20);
    
    %%%Step4:��IMFs���з���%%%
    %%%����㷨
    [m,n] = size(modes);
    zcrnum = [];
    for ii = 1:m
        [num] = scr2(modes(ii,:));      
        zcrnum(ii) = num;
    end

    m = [];
    m = find(zcrnum<=zcrnumber);        %Select the imf dominated by the baseline according to the number of zero crossings/���ݹ����ĸ���ѡ���ɻ���������imf
    imfbegin = m(1,1);
    imfend = m(1,end);

    
    %%%Find baseline%%%
    if length(Lsection)+length(SigButter) >= n   %Determine whether the index exceeds the matrix dimension/�ж������Ƿ񳬳�����ά��
        temp_array = ones(m,length(Lsection)+length(SigButter)-n+1);
        temp_array = modes(:,end).*temp_array;
        modes = [modes temp_array];
        
    end
    
    if imfbegin == imfend   %Determine whether the number of IMFs dominated by the baseline is 1
        baseline = modes(imfend,length(Lsection)+1:length(Lsection)+length(SigButter));
        DBSig = SigButter - baseline;
    else
        baseline = sum((modes(imfbegin:imfend,length(Lsection)+1:length(Lsection)+length(SigButter))));
        DBSig = SigButter - baseline;   
    end
        

    %**************************************************************************
    %Find all local extreme points of the signal/���ź����оֲ���ֵ��
    %[tmax,tmin,SigButtermax,SigButtermin] = fink_local_peaks(sig)
    %*************************************************************************
    function [tmax,tmin,SigButtermax,SigButtermin]  = fink_local_peaks(sig)
        %%%�Ե�֮ͨ����ź���ȡ��ֵ��%%%
        t1 = 1:length(sig);
        %%%���źŵ����м�ֵ���λ��
        Lmax = diff(sign(diff(sig)))== -2;
        Lmin = diff(sign(diff(sig)))== 2; 
        %%%���źŲ����ԭ�е����ݳ���
        Lmax = [false, Lmax, false];        
        Lmin =  [false, Lmin, false];
        %%%��ע��ֵ��������λ��
        tmax = t1(Lmax);
        tmin = t1(Lmin); 
        %%%��ע��ֵ���������λ��
        SigButtermax = sig(Lmax);
        SigButtermin = sig(Lmin); 
    end 
    
    
    %**************************************************************************
    %Truncate each peak, find the length of each segment, and judge whether there is interference between them
    %��ÿ���׷���нضϣ����ÿһ�εĳ��ȣ����ж�֮���Ƿ�������
    %[subsection1,subsection_min] = delete_ms_peaks(col,tmax,tmin,SigButtermax,SigButtermin)
    %*************************************************************************
    function [subsection1,subsection_min] = delete_ms_peaks(col,tmax,tmin,SigButtermax,SigButtermin)
        
        subsection = zeros(3,length(col));      %Store each start and end point of the truncation into the subsection matrix/���ضϵ�ÿһ����ʼ��ֹ�����subsection������
        for  i2 = 1:length(col)
            n = find(xdata == col(i2)) - 2;
            while SigButtermax(n) > SigButtermax(find(xdata == col(i2)))
                n = n-1;
            end
            subsection(1,i2) = tmax(n);
            m = find(xdata == col(i2)) + 2;
            while SigButtermax(m) > SigButtermax(find(xdata == col(i2)))
                m = m + 1;
            end
            subsection(2,i2) = tmax(m);
        end
        %%%����subsection�����е�ÿ�γ��ȴ��������
        for i2 = 1:length(col)
            subsection(3,i2) = (subsection(2,i2)-1) - (subsection(1,i2)+1)+1;
        end
        %%%��һ�θ����жϣ��ж϶˵�֮���Ƿ�������
        LogicM1 = ones(1,length(col));
        for i2 = 2:length(col)
            if subsection(1,i2) - subsection(2,i2-1) < 0
                LogicM1(i2) = 0;
            end
        end
        %%%����ɾ��subsection�����������Ƭ�Σ�Ҳ����LogicM1���ҵ��������Ϊ0�Ķ�   
        i2 = 1;
        while i2 <length(LogicM1)
            if LogicM1(1,i2) ==0 && LogicM1(1,i2+1) == 0
               LogicM1(:,i2) = [];
               subsection(:,i2) = [];
            else
                 i2 = i2 + 1;
            end
        end
        %%��һ�θ����ж�֮��ϲ����������λ��
        subsection1 = zeros(size(subsection));
        for i2 = 1:length(LogicM1)
            if LogicM1(i2) ~= 0
                subsection1(:,i2) = subsection(:,i2);
            else
                subsection1(2,i2-1) = subsection(2,i2);
            end
        end
        subsection1(:,all(subsection1==0)) = [];        %ɾ��ȫ����
        %����subsection1�����е�ÿ�γ��ȴ��������
        [m,n] = size(subsection1);
        for i2 = 1:n
            subsection1(3,i2) = (subsection1(2,i2)-1) - (subsection1(1,i2)+1)+1;
        end
        %%%Ѱ��subsection1�еļ�Сֵ���Ƭ��subsection_min
        Sigmax = SigButtermax;Sigmin = SigButtermin;
        [m,n] = size(subsection1);
        subsection_min = [];
        for i2 = 1:n
            zR = find(tmin > subsection1(2,i2));
            zR1 = zR(1,1);
            zR2 = zR(1,2);
            subsection_min(3,i2) = tmin(1,zR1);
            subsection_min(4,i2) = tmin(1,zR2);
            zL = find(tmin < subsection1(1,i2));
            zL1 = zL(1,end);
            zL2 = zL(1,end-1);
            subsection_min(1,i2) = tmin(1,zL1);  
            subsection_min(2,i2) = tmin(1,zL2);  
        end
        for i2  =1:n
            subsection_min(6,i2) = (subsection_min(3,i2)-1) - (subsection_min(1,i2)+1) + 1;
        end
        clear m 

    end

    
    %**************************************************************************
    %Cubic Spline Interpolation Method to Complement Spectral Peaks/����������ֵ�������׷�
    %[SigButterRIS] = cubic_splin_interpolation(subsection,subsection_min)
    %*************************************************************************
    function [SigButterRIS] = cubic_splin_interpolation(subsection,subsection_min,SigButtermax,SigButtermin)
        subsection1 = subsection;
        trunction_sig = cell(2*n,1);
        ttmax = cell(1,n);      
        ttmin = cell(1,n);
        
        %the minPoint is used to store the ordinate of the minimum point, the first line is the minimum point on the left side of the spectrum peak, 
        %and the second line is the minimum point on the right side of the spectrum peak
        %minPoint�������漫Сֵ���������,��һ�����׷���༫Сֵ�㣬�ڶ������׷��Ҳ༫Сֵ��
        minPoint = [];      
        for i =1:n
            %Cubic Spline Difference of Maximum Point Pair/����ֵ����������ֵ
            z = find(abs(tmax - subsection1(1,i)) < 1);
            Lysigmax1 = SigButtermax(1,z);
            z = find(abs(tmax - subsection1(2,i)) < 1);
            Rysigmax1 = SigButtermax(1,z);
            x = [subsection1(1,i),subsection1(2,i)];
            y = [Lysigmax1,Rysigmax1];
            ttmax{1,i} = linspace(subsection_min(1,i)+1,subsection_min(3,i)-1,subsection_min(6,i));
            trunction_sig{(i-1)*2+1} = spline(x,y, ttmax{1,i});
            %Cubic spline difference of minimum point pairs/��Сֵ����������ֵ
            z = find(abs(tmin - subsection_min(1,i)) < 1);
            Lysigmin1 = SigButtermin(1,z);
            minPoint(1,i) = Lysigmin1;
            z = find(abs(tmin - subsection_min(3,i)) < 1);            
            Rysigmin1 = SigButtermin(1,z);
            minPoint(2,i) = Rysigmin1;
            x = [ subsection_min(1,i),subsection_min(3,i)];
            y = [Lysigmin1,Rysigmin1];
            ttmin{1,i} = linspace(subsection_min(1,i)+1,subsection_min(3,i)-1,subsection_min(6,i));
            trunction_sig{(i-1)*2+2} = spline(x,y,ttmin{1,i});
        end
        %%%Calculate the mean value of the difference between the maximum and minimum values of each segment of the cubic spline
        %%%��ÿһ�μ���ֵ��Сֵ����������ֵ�ľ�ֵ
        trunction_mean = cell(n,1);
        for i  = 1:n
            trunction_mean{i,1} = (trunction_sig{(i-1)*2+1}+trunction_sig{(i-1)*2+2})./2;
            trunction_amp(1,i) = (mean(trunction_sig{(i-1)*2+1})-mean(trunction_sig{(i-1)*2+2}))./2;
        end
        %%%���뵽ԭʼ�ź�֮�ϣ��۲�Ч��
        SigButterRIS = SigButter;
        for i = 1:n
            z = trunction_mean{i,1};
            LengthRIS = subsection1(2,i)-subsection1(1,i);
            SigButterRIS(1,subsection1(1,i):subsection1(2,i)) = z(1,subsection1(1,i)-subsection_min(1,i):LengthRIS+subsection1(1,i)-subsection_min(1,i));
        end
    end


    %**************************************************************************
    %Flip both ends to remove the end effect of emd/���˷�תȥ���˵�ЧӦ
    %[CompSigButter] = reduce_end_effect(SigButterRIS,FliplrLength)
    %*************************************************************************
    function [CompSigButter,Lsection] = reduce_end_effect(SigButterRIS,FliplrLength)
        %Cut off the signal from the two ends at the outermost maximum point/���źŴ������������ļ���ֵ�㴦�ض�
        TruSigButterRIS = SigButterRIS(1,tmax(1,1):tmax(1,end));        %TruSigButterRISȥ����Ъ�źź�֮����ź�
        tmaxL = round(length(tmax)/(FliplrLength));      %tmaxL-�ӵ�tmaxL������ֵ��ض�
        if tmaxL == 1
            tmaxL = tmaxL + 2;
        end
        Lendpoint = tmax(1,tmaxL);
        Rendpoint = tmax(1,end-tmaxL);
        Lsection = fliplr(SigButterRIS(1,tmax(1,1)+1:Lendpoint));
        Rsection = fliplr(SigButterRIS(1,Rendpoint:tmax(1,end)-1));
        CompSigButter = [Lsection,TruSigButterRIS,Rsection];
    end
    
    
    %**************************************************************************
    %Single-step method to find the number of zero crossings/���������������
    %[scrnum] = scr2(sig)
    %*************************************************************************
    function [scrnum] = scr2(sig)
        sig1 = sig(1,1:end-1);
        sig2 = sig(1,2:end);
        N = length(sig1);
        scrnum = 0;
        i = 1;
        while i <= N
            if sig1(1,i)*sig2(1,i) < 0
                scrnum = scrnum + 1;
            end
            i = i + 1;
        end
    end


    %**************************************************************************
    %Polynomial fitting to determine the number of peaks/����ʽ���ȷ����ֵ����
    %[peak_number,zcrnumber,gof,x,y,col] = polyfit_ployval_peaknum_zcrnumber(xdata,ydata,alpha,SigButter)
    %*************************************************************************
    function [peak_number,zcrnumber,gof,x,y,col] = polyfit_ployval_peaknum_zcrnumber(xdata,ydata,alpha,SigButter)
        [xData, yData] = prepareCurveData(xdata,ydata);

        % Set up fittype and options.
        ft = fittype( 'poly9' );
        opts = fitoptions( 'Method', 'LinearLeastSquares' );
        opts.Robust = 'LAR';

        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );


        figure(1)
        sigbutter_x = 1:1:length(SigButter);
        plot(sigbutter_x,SigButter,'linewidth',2.0,'color',[190 190 190]./255);
        hold on;
        h = plot( fitresult, xData, yData, 'predobs', alpha);set(h,'linewidth',2.0);
        set(gca,'FontSize',16,'Fontname', 'Times New Roman');
        xlabel('Variables','FontSize',20);ylabel('Intensity','FontSize',20);
        h=findobj(gca,'Type','Line');   % Extract curve data object/��ȡ�������ݶ���
        x = get(h,'xdata');                                                            
        y = get(h,'ydata');
        legend('original signal', 'maximal point',  'maximal point fit line','prediction bounds');

        y_No1 = [];
        y_length = length(y);
        for i = 1:y_length-2
            y_No1(1,i) = y{i,1}(1,1);       
        end
        y_max_loc = find(y_No1 == max(y_No1));
        up_confidence_line_y = y{y_max_loc,1};      %Output the upper bound coordinates of the confidence interval/����ͼ����ȡ�����ļ����߶��������������������ߣ��������������ߵ�y��ֵ��x��ֵ
        up_confidence_line_x = x{y_max_loc,1};
        %Output the number of peaks removed peak_number, and store the abscissa of each peak to be removed/���ȥ����ֵ����peak_number,������ÿ����Ҫȥ����ֵ�ĺ�����
        peak_number = 0;
        col = [];       %the col is used to store the abscissa of the extreme point higher than the upper limit of the confidence interval/col����������������������޵ļ�ֵ������꣬��ʵ�����׷������
        i = 1;
        while i <= length(xdata)
            if ydata(i) > up_confidence_line_y(find(up_confidence_line_x == xdata(i)))
                peak_number = peak_number+ 1;
                col = [col xdata(i)];
            end
            i = i + 1;
        end
        %Output the number of zero-crossing points or extreme points of the rough baseline/����ֲڻ��ߵĹ�������/��ֵ�����
        %The number of zero crossing points of rough baseline/�ֲڻ��߹�������
%         up_confidence_line_y_mean = up_confidence_line_y - mean(up_confidence_line_y);
%         [up_confidence_line_y_mean_scrnum] = scr2(up_confidence_line_y_mean);

        %The number of extreme points/��ֵ�����
        [tmax_baseline,tmin_baseline,~,~]  = fink_local_peaks(up_confidence_line_y);
        up_confidence_line_y_mean_scrnum = length(tmax_baseline) + length(tmin_baseline);
        
        display(up_confidence_line_y_mean_scrnum);
        zcrnumber = up_confidence_line_y_mean_scrnum;     %�������
    end
end

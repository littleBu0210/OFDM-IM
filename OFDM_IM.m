clc;clear;close all;
% Matlab Version 2023b
%% main loop parameter
M=4;                                    % QAM order                  
EbN0dB=(0:5:25)';                       % Eb/N0 'dB'
N_iter=1e3;                             % Number of iterations for each EbN0
Ngroup=4096;                            % Number of subcarrier group per iteration
channel = "Rayleigh";                   % "Rayleigh" or "Gaussian"
%% variable
Nbps=log2(M);                           % Bit per symbols
BER = zeros(length(EbN0dB),1);          % Bit error rate
SAPER = zeros(length(EbN0dB),1);        % subcarrier activation pattern error rate
QAMBER = zeros(length(EbN0dB),1);       % qam bit error rate
%% OFDM-IM parameter
NsubCarry = 4;                                      % All carry numbers per group
NactiveCarry = 2;                                   % active carry numbers per group
p1 = floor(log2(nchoosek(NsubCarry,NactiveCarry))); % index bit length per group
p2 = NactiveCarry*Nbps;                             % information bit length per group
p=p1+p2;                                            % total number of bits per group
NbitPerIter = p*Ngroup;                             % total bit per iteration 
SNRdB = EbN0dB+10*log10(p/NsubCarry);               % EbN0 to SNRdB
SNR = 10.^(SNRdB/10);                               % SNRdB to SNR
N_T = 1./SNR;                                       % Time domain noise power,Assume the signal average power is 1 per symbol
N_F = NactiveCarry/NsubCarry.*N_T;                  % Frequency domain noise power
% look up table
if(NactiveCarry==2&&NsubCarry==4)
    allPattern = [1 0;2 0;3 1;3 2;3 0;2 1];          % optimal combination for this case (maximum distance)
    % allPattern = Combin_Md(NsubCarry,NactiveCarry);
else
    allPattern = Combin_Md(NsubCarry,NactiveCarry);
end
index_allz= fliplr(uint16(allPattern+1))';
% index select offset
IndexMat = uint16(reshape(repmat((0:1:Ngroup-1)*NsubCarry,NactiveCarry,1),[],1));
% reference symbol 
s_ref = qammod(0:M-1,M,'gray',UnitAveragePower=true)';
Ntb=NbitPerIter*N_iter;                             % Initialize the number of total bits
NtSAP=Ngroup*N_iter;                                % Initialize the number of error subcarrier activation patton, SAP
NtQam=p2*Ngroup*N_iter;                             % Initialize the number of error QAM bit
%% main loop
tic;
for i=1:length(EbN0dB) 
    Neb=0;   % Initialize the number of error bits
    NeSAP=0; % Initialize the number of error subcarrier activation patton, SAP
    NeQam=0; % Initialize the number of error QAM bit
    for j = 1:N_iter
        %all bits
        data = randi([0 1], NbitPerIter,1,'uint16');

        %index bits
        DataIndex = data(1:p1*Ngroup);

        %information bit
        DataQam = data(p1*Ngroup+1:end);


        %Index bit processing,using the offset array to comfirm the final index
        DecimalIndex = bit2int(DataIndex,p1)+1;
        FinIndex = IndexMat+reshape(index_allz(:,DecimalIndex),[],1);


        %QAM modulation,unitAvgPower is needed
        ConstellData = qammod(DataQam,M,'gray','InputType','bit',UnitAveragePower=true);

        %Transmit Symbol
        TransmitSymbol = zeros(NsubCarry,Ngroup,'like',1i);
        TransmitSymbol(FinIndex) = ConstellData;

        %channel select
        if channel=="Rayleigh"
            h = 1/sqrt(2)*(randn(NsubCarry,Ngroup)+1i*randn(NsubCarry,Ngroup));
        elseif channel=="Gaussian"
            h = ones(NsubCarry,Ngroup);
        end
        
        %noise
        noise = sqrt(N_F(i)/2)*(randn(NsubCarry,Ngroup)+1i*randn(NsubCarry,Ngroup));

        %receive signal
        yF = TransmitSymbol.*h+noise;

        %Using LLR Recieve
        LLR = LLRCalculate(yF,NsubCarry,NactiveCarry,h,s_ref,N_F(i));

        %Decode index bits
        [~,Order] = sort(LLR,'descend');
        %find the max likely position
        SortOrder = sort(uint16(Order(1:NactiveCarry,:)),1,"ascend");
        DecodeIndex = zeros(size(DecimalIndex),'uint16');
        %match
        for k =1:length(index_allz)
           DecodeIndex(sum(SortOrder==index_allz(:,k),1)==NactiveCarry) = k;
        end
        ErrorSAP = Ngroup - length(find((DecodeIndex-DecimalIndex)==0));
        NeSAP = NeSAP+ErrorSAP;

        %Statistical pattern error bits
        Neb=Neb+sum(xor(int2bit(DecodeIndex-1,p1),DataIndex));
        %information bit index
        QamDecodeIndex = IndexMat+reshape(SortOrder,[],1);
        %balanced
        yB = yF(QamDecodeIndex).*conj(h(QamDecodeIndex));
        %QAM demodulation
        DataQamDecode = qamdemod(yB,M,'gray',UnitAveragePower=true,NoiseVariance=N_F(i));
        %Convert to binary
        DataQamDecodeBit = int2bit(DataQamDecode,Nbps);
        %Number of information bit errors
        ErrorQam = sum(xor(DataQamDecodeBit,DataQam));
        Neb=Neb+ErrorQam;
        NeQam=NeQam+ErrorQam;
    end
    fprintf('EbN0=%.2f[dB]\n', EbN0dB(i))
    SAPER(i) = NeSAP/NtSAP;
    fprintf('SAP_ER=%4d/%8d =%11.3e\n',NeSAP,NtSAP,SAPER(i))
    BER(i) = Neb/Ntb;
    fprintf('BER=%4d/%8d =%11.3e\n', Neb,Ntb,BER(i))
    QAMBER(i) = NeQam/NtQam;
    fprintf('QAM_ER=%4d/%8d =%11.3e\n', NeQam,NtQam,QAMBER(i))
end
toc;
%% Plot
semilogy(EbN0dB,SAPER,'b--s',"DisplayName","SPA error");  
hold on
semilogy(EbN0dB,QAMBER,'g--d',"DisplayName","QamDecode error");  
hold on
semilogy(EbN0dB,BER,'r--o',"DisplayName",sprintf("BER n=%d t=%d",NsubCarry,NactiveCarry));
hold on
%Theoretical OFDM Curve Plotting
ber_Rayleigh = ber_QAM(EbN0dB,M,'Rayleigh');
semilogy(EbN0dB,ber_Rayleigh,'r-',"DisplayName","tradition M-ary OFDM");
hold on

grid on
axis([EbN0dB(1) EbN0dB(end) 1e-5 1])
xlabel("Eb/N0")
ylabel("BER")
legend()
%% save
isSave = 1;
if isSave
    folder='./res/';
    if exist(folder,"dir")==0
        mkdir(folder)
    end
    TimeString = string(datetime("now"),'yyyy-MM-dd');
    FileName = strcat(folder,TimeString,'.mat');
    save(FileName,'BER','EbN0dB','N_iter',"Ngroup","channel",'M',"NsubCarry","NactiveCarry");
end
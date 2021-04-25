% Examples of relevant channel models, 20MHz BW (labels with x,y,z,etc contain configurable parameters):

% LTE (https://www.etsi.org/deliver/etsi_ts/136100_136199/136104/10.02.00_60/ts_136104v100200p.pdf)
% 'EPA_5Hz' -> Extended Pedestrian A model (EPA) with 5Hz doppler frequency (table B.2-1, page 85)
% 'EVA_70Hz' -> Extended Vehicular A model (EPA) with 70Hz doppler frequency (table B.2-1)
% 'ETU_300Hz' -> Extended Typical Urban model (ETU) with 300Hz doppler frequency (table B.2-3)

% NR (section 7.7, page 65, https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/15.00.00_60/tr_138901v150000p.pdf)
% 'TDL-x_yHz_zns' -> Tapped Delay Line (CDL) models (simplified case, non-MIMO evaluations)
% parameters: x = A,B,...,E power delay profiles (pdp), y = doppler frequency, z = desired delay spread (table 7.7.3-1, from 10ns to 1000ns)
function H = EPA_Uplink_TDD

ChannelModel = 'EPA_5Hz';%'TDL-A_300Hz_1000ns'%

% auxiliary parameters
SF = 5;      % Number of subframes (see http://niviuk.free.fr/lte_resource_grid.html or http://niviuk.free.fr/nr_frame.php)
SCS = 15;    % subcarrier space in KHz
NDLRB = 100; % Number of resource blocks [6,15,25,50,100 for 1.4MHz,3MHz,5Mhz,10Mhz,20MHz]
NULRB = 100;

if strcmp(ChannelModel(1),'E') % LTE
    % LTE parameters
    enb.NDLRB = NDLRB;              % Number of resource blocks
    enb.NULRB = NULRB;              % Number of resource blocks
    enb.CellRefP = 1;               % One transmit antenna port
    enb.NCellID = 0 ;               % Cell ID
    enb.CyclicPrefixUL = 'Normal';  % 'Normal' or 'Extended' cyclic prefix
    enb.DuplexMode = 'TDD';         % FDD/TDD
    enb.TotSubframes = SF;          % Number of subframes
    enb.NTxAnts = 1;
    
    % OFDM parameters
    info = lteOFDMInfo(enb);
    
    % LTE channel
    channel.Seed = randi(2^31); % seed for multipath generation
    channel.DelayProfile = ChannelModel(1:3); %EVA/EPA/ETU
    channel.DopplerFreq = str2num(ChannelModel(5:(find(ChannelModel=='H')-1))); % doppler frequency in Hz
    channel.ModelType = 'GMEDS';
    channel.MIMOCorrelation = 'UplinkMedium';
    channel.NRxAnts = 64;
    channel.InitTime = 0;
    channel.InitPhase = 'Random';
    channel.NormalizePathGains = 'On';
    channel.NormalizeTxAnts = 'On';
    channel.SamplingRate = info.SamplingRate;
    H=single(lteULPerfectChannelEstimate(enb,channel));
elseif strcmp(ChannelModel(1),'T') || strcmp(ChannelModel(1),'C') % NR
    SR = 30.72e6;                   % Sample rate (e.g. SR = 30.72e6 for 20MHz BW)

    % NR channel
    tdl = nrTDLChannel;
    tdl.NumReceiveAntennas = 1;
    tdl.DelayProfile = ChannelModel(1:5);
    tdl.MaximumDopplerShift = str2num(ChannelModel(7:(find(ChannelModel=='H')-1)));
    tdl.DelaySpread = 1e-9*str2num(ChannelModel((find(ChannelModel=='H')+3):(find(ChannelModel=='n')-1)));
    tdl.SampleRate = SR;

    % waveform parameters
    T = SR*1e-3*SF;                 % create a waveform with the duration of SF subframes (SF*1ms)
    tdlInfo = info(tdl);
    Nt = tdlInfo.NumTransmitAntennas;
    in = complex(randn(T,Nt),randn(T,Nt));
    pathFilters = getPathFilters(tdl);
    [~,pathGains] = tdl(in);

    nSlot = 0; % initial slot

    H = single(nrPerfectChannelEstimate(pathGains,pathFilters,NDLRB,SCS,nSlot));
end

H = H(1:64,1:64,:);
% figure
% surf(abs(H))
% t=title(ChannelModel);set(t,'Interpreter','none');
% xlabel('symbols')
% ylabel('subcarriers')
% zlabel('Magnitude')

end

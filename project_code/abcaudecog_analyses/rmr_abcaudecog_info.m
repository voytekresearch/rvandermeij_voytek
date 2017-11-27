function info = rmr_abcaudecog_info


% get hostname for paths
[dum hostname] = system('hostname');
if ~isempty(strfind(hostname,'sdsc.edu'))
  info.ontscc = true;
else
  info.ontscc = false;
end

% paths
if info.ontscc
  pathprefix = '/projects/ps-voyteklab/';
else
  pathprefix = '/Volumes/voyteklab/';
end
info.datapath = [pathprefix 'common/data2/BidetCaulet_ECoG/'];
info.savepath = [pathprefix 'roevdmei/abcaudecog/'];


% subject names
info.subj = {'GP15','GP22','GP28','GP35','JH2','ST1','ST6','ST8','JH13'}; % GP13 and JH11 removed due to bad performance/lack of button press. GP14 removed due to anatomical issues underneath electrodes.

% session names
%info.session.GP14   = {'GP14_B8','GP14_B9','GP14_B13','GP14_B15'}; % GP14_B7/GP14_B14 discarded by ABC
info.session.GP15   = {'GP15_B12','GP15_B14','GP15_B19','GP15_B20','GP15_B21'}; % GP15_B13 discarded by ABC
info.session.GP22   = {'GP22_B40','GP22_B41','GP22_B43','GP22_B44','GP22_B45','GP22_B46','GP22_B47'}; % GP22_B42 discarded by ABC
info.session.GP28   = {'GP28_B36','GP28_B38','GP28_B40','GP28_B41','GP28_B42','GP28_B43'}; % GP28_B39 discarded by ABC
info.session.GP35   = {'GP35_B29','GP35_B30','GP35_B31','GP35_B34','GP35_B35','GP35_B36'}; % GP35_B32/GP35_B33 discard by ABC
info.session.JH2    = {'Beta_Edfimport_EEG_6507'};
info.session.JH13   = {'Edfimport_EEG_10001'};
info.session.ST1    = {'stanford01_B23','stanford01_B24','stanford01_B25','stanford01_B26','stanford01_B27','stanford01_B28','stanford01_B29'};
info.session.ST6    = {'ST06_09','ST06_10','ST06_11','ST06_12','ST06_13','ST06_14','ST06_15'};
info.session.ST8    = {'ANP027','ANP028','ANP029','ANP030','ANP031','ANP032','ANP033'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandstop frequencies found using rmr_findlinespectra
% For those with sampling at ~3k, checked till 540, with the intention of really sharp low pass at 540. 
% %info.GP14.bsfreq     = [60-2 60+2; 120-2 120+2; 180-3 180+3; 256-4 256+4]; % checked
info.GP15.bsfilt.peak      = [60.0   98.1   98.7  100.2  119.9  148.2  180.0  196.2  197.3  200.2  239.9  251.9  294.3  296.0  300.3  340.0  359.8  392.5 394.6  400.5  419.9  460.0  479.5  480.0  490.6  493.2  500.7  539.5  539.9];
info.GP15.bsfilt.halfbandw = [0.25    0.25    0.25    0.30    0.25    0.25    0.50    0.60    0.60    0.80    0.50    0.25    0.80    0.80    1.30    0.25    0.50    1.10   1.10    1.90    0.70    0.25    0.25    0.25    1.40    1.30    2.    0.30    0.25];
info.GP22.bsfilt.peak      = [60.00  119.9  123.3  124.3  180.0  203.7  239.8  246.6  248.2  249.2  256.1  284.6  299.9  325.3  359.8  369.9  372.3  373.8 376.1  392.1  419.7  447.2  479.8  493.3  496.3  498.4  512.1  528.5  539.6  365.9];
info.GP22.bsfilt.halfbandw = [0.4    0.4    0.2500    0.9    1.9    0.2500    0.8    0.4    1.3    0.4    0.3    0.2500    1.4    0.2500    0.9    0.6    1.9    0.6  0.2500    0.2500    1.0    0.2500    0.9    0.9    2.6    0.8    0.3    0.2500    0.8    0.2500];
info.GP28.bsfilt.peak      = [60.0   89.5  176.1  176.3  176.5  177.5  178.8  180.1  256.0  264.4  267.5  300.1  352.5  356.7  420.1  440.6  445.9  512.1 527.2  528.7  535.0];
info.GP28.bsfilt.halfbandw = [0.25    0.70    0.25    0.25    0.25    0.50    1.70    0.30    0.25    1.10    4.    0.40    1.70    5.60    0.30    1.70    6.40    0.25  0.25    2.60    8.95];
info.GP35.bsfilt.peak      = [60.0   85.4   86.0   87.3   88.8   92.7   93.0   93.6   94.3  171.7  174.6  177.4  180.0  185.7  187.3  188.5  257.5  261.8 266.2  280.4  300.0  343.3  349.1  354.8  373.9  426.7  429.5  436.4  443.6  464.1  468.0  471.3  514.9  523.6  532.3];
info.GP35.bsfilt.halfbandw = [0.25    0.25    0.70    0.50    0.50    0.25    0.25    0.40    0.30    2.10    1.30    1.60    0.25    1.10    1.10    0.80    3.65    2.20  3.25    6.85    0.50    5.35    3.10    4.35    9.05    0.60    5.25    3.10    4.45    3.85    2.60    2.75    7.35    4.75    5.25];
info.JH2.bsfilt.peak       = [60.0  180.0  299.9  340.2  359.6  419.9  460.2]; % noise floor is probably around 200
info.JH2.bsfilt.halfbandw  = [0.25    0.30    0.40    0.25    0.25    0.40    0.25];
info.JH13.bsfilt.peak      = [60 119.9 180 239.8 300.1];
info.JH13.bsfilt.halfbandw = [0.25 0.25 0.50 0.25 0.50];
info.ST1.bsfilt.peak       = [51.8   60.0   83.2   83.5   84.5   84.7   85.1   85.5   85.9   86.9   87.2  120.1  166.8  170.3  173.7  180.1  200.1  201.8 202.4  203.2  203.7  204.2  204.7  205.8  206.4  206.9  207.2  207.5  207.8  208.1  208.4  210.7  211.5  212.2  240.1  250.3 257.0  300.1  325.5  333.8  342.6  360.1  417.0  420.8  423.8  427.3  430.1  433.7  436.0  500.5  507.3  514.4  521.0]; % TERRIBLE FOR SLOPE FITTING!
info.ST1.bsfilt.halfbandw  = [0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    1.20    4.50    1.80    0.25    0.50    0.50 0.50    0.25    0.30    0.25    0.50    0.25    0.50    0.25    0.25    0.25    0.25    0.25    0.50    0.50    0.50    0.50    0.25    1.90   10.35    0.25    0.55    2.90   13.85    0.25    3.25    1.80    3.60    2.50    2.40    3.30    1.    4.15    7.    6.70    6.05];
info.ST6.bsfilt.peak       = [60.0  120.0  180.0  239.7  240.7  241.6  284.6  300.0  325.3  360.1  378.9  400.5  420.0  480.0  487.8  539.9];
info.ST6.bsfilt.halfbandw  = [0.80    0.25    0.25    1.50    0.25    0.60    0.25    0.40    0.25    0.30    0.40    6.40    0.40    0.25    0.25    0.25];
info.ST8.bsfilt.peak       = [60.1   97.7  100.0  120.0  180.0  195.3  240.0  248.0  249.0  289.7  290.1  291.0  292.8  300.0  360.0  384.7  394.3  420.0  480.0  482.9  483.6  487.2  500.0  539.9];
info.ST8.bsfilt.halfbandw  = [2.50    0.25    0.25    0.25    0.40    0.25    0.40    0.50    0.25    0.25    0.25    0.60    2.50    0.90    0.30    0.25    0.80    0.40 0.25    0.25    0.30    5.60    0.25    0.25];
%%%%%%%%%%%%%%%%%%%%%
% Length of edge artifact longest filter above
info.GP15.bsfilt.edgeartlen = 3.8332;
info.GP22.bsfilt.edgeartlen = 3.8332;
info.GP28.bsfilt.edgeartlen = 3.8332;
info.GP35.bsfilt.edgeartlen = 3.8332;
info.JH2.bsfilt.edgeartlen  = 3.8360;
info.JH13.bsfilt.edgeartlen = 3.8360;
info.ST1.bsfilt.edgeartlen  = 3.8358;
info.ST6.bsfilt.edgeartlen  = 3.8325;
info.ST8.bsfilt.edgeartlen  = 3.8325;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD: Bandstop frequencies for eliminating line noise (not doing it session specific unless needed)
% For Checked beyond 250+, WITH INTENTION of low pass at 250 (some without it)
% %info.GP14.bsfreq     = [60-2 60+2; 120-2 120+2; 180-3 180+3; 256-4 256+4]; % checked
% info.GP15.bsfreq     = [60-.25 60+.25; 98-1 98+1; 100-.25 100.5+.25; 120-.25 120+.25; 148-.25 148+.25; 179.8-.25 179.8+.25; 180-.25 180+.25; 195-.25 198+.25; 199.5-.25 200.5+.25; 216-.25 220+.25; 239-.25 240+.25; 251.8-.25 251.8+.25; 293.5-.25 296.3+.25; 299-.25 301+.25; 310-.25 313+.25; 325.2-.25 325.2+.25; 340-.25 340+.25; 359.6-.25 360.5+.25; 392-1 395.5+.25; 399-.25 402+.25; 419-.25 420+.25; 434.5-.25 439.25; 460-.25 460+.25; 479.5-.25 480.5+.25; 489-.25 494+.25; 499-.25 502+.25; 503.7-.25 503.7+.25; 539 541; 547 549]; % checked
% info.GP22.bsfreq     = [60-.25 60+.25; 120-.25 120+.25; 180-1 180+1; 240-.25 240+.25; 247-.25 247+.25; 256-.25 256+.25]; % checked
% info.GP28.bsfreq     = [60-.25 60+.25; 88-.25 88+.25; 89.3-.25 89.3+.25; 176-.25 176+.25; 177-.25 177+.25; 178-.25 178+.25; 178.5-.25 178.5+.25; 179-.25 179+.25; 180.05-.25 180.05+.25; 256-.25 256+.25; 264-.25 264+.25]; % checked; 
% info.GP35.bsfreq     = [60-.25 60+.25; 85-.25 89+.25; 85-.25 89+.25; 93.2-.25 93.6+.25; 170.3-.25 172.5+.25; 173-.25 178+.25; 180-.25 180+.25; 185-.25 188.5+.25; 255.8-.25 255.8+.25; 256-.25 256+.25; 261-.25 261+.25]; % checked; 
% info.JH2.bsfreq      = [60-.25 60+.25; 180-.25 180+.25]; % checked; 
% info.JH13.bsfreq     = [60-.25 60+.25; 120-.25 120+.25; 180-.25 180+.25]; % checked;
% info.ST1.bsfreq      = [51.8-.25 51.8+.25; 60-.25 60+.25; 83.1-.25 83.1+.25; 103.5-.25 103.5+.25; 120-.25 120+.25; 165-.25 174.6+.25; 240-.25 240+.25; 248-.25 261+.25]; % checked
% info.ST6.bsfreq      = [60-.25 60+.25; 120-.25 120+.25; 180-.25 180+.25; 240-.25 240+.25]; % checked 
% info.ST8.bsfreq      = [60-.25 60+.25; 97.3-.25 97.3+.25; 100-.25 100+.25; 120-.25 120+.25; 180-.25 180+.25;]; % checked; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD - narrowed bands down, selectively, based on aliassing/narrowness of line noise and 2 order filters
% Bandstop frequencies for eliminating line noise (not doing it session specific unless needed)
% NOT checked beyond 250+, WITH INTENTION of low pass at 250
% %info.GP14.bsfreq     = [60-2 60+2; 120-2 120+2; 180-3 180+3; 256-4 256+4]; % checked
% info.GP15.bsfreq     = [60-2 60+2; 100-3 100+3; 120-3 120+3; 180-3 180+3; 196-3 196+3; 200-4 200+4; 240-4 240+4]; % checked
% info.GP22.bsfreq     = [60-2 60+2; 120-3 120+3; 180-3 180+3; 240-4 240+4; 247-4 247+4; 256-4 256+4]; % checked
% info.GP28.bsfreq     = [60-2 60+2; 120-3 120+3; 176-3 176+3; 180-3 180+3; 264-4 264+4]; % checked; 120 probably not needed
% info.GP35.bsfreq     = [60-2 60+2; 120-3 120+3; 171-3 188+3; 256-4 256+4]; % checked; 120 probably not needed. Collapsed 171, 174, 180, 185 and 188 into one
% info.JH2.bsfreq      = [60-2 60+2; 120-3 120+3; 180-3 180+3]; % checked; 120 probably not needed
% info.JH13.bsfreq     = [60-2 60+2; 120-3 120+3; 180-3 180+3]; % checked; 120 probably not needed
% info.ST1.bsfreq      = [83-2 83+2; 166-3 166+3; 171-3 256+4]; % checked combined several high freq because there was a continous 'ribbon' riding many channels
% info.ST6.bsfreq      = [60-2 60+2; 120-3 120+3; 180-3 180+3]; % checked 
% info.ST8.bsfreq      = [60-2 60+2; 120-3 120+3; 180-3 180+3]; % checked; 120 probably not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad channels (not doing it session specific unless absolutely needed)
%info.GP14.badchan     = {'chan13','chan45','chan46','chan47','chan48','chan61','chan62','chan63','chan64'}; % checked; first 5 from ABC (65-84 empty). Missed by ABC: 61-64 are bad in last 3 blocks.
%info.GP14.badchan     = [info.GP14.badchan, cellstr([repmat('chan',[20 1]) strtrim(num2str((65:84)'))])']; % adding the empties
info.GP15.badchan     = {'chan8','chan46'}; % checked; all by ABC
info.GP15.badchan     = [info.GP15.badchan, cellstr([repmat('chan',[20 1]) strtrim(num2str((65:84)'))])']; % adding the empties
info.GP22.badchan     = {'chan33','chan34','chan35','chan36','chan44','chan48','chan17'}; % checked; 33-36,44,48 by ABC. ABC also removed 1,2 due to flat, and 12,13 due to beta bursts. However, I don't see them...; Missed by ABC: 17 flat
info.GP22.badchan     = [info.GP22.badchan, cellstr([repmat('chan',[35 1]) strtrim(num2str((65:99)'))])' cellstr([repmat('chan',[29 1]) strtrim(num2str((100:128)'))])']; % adding the empties
info.GP28.badchan     = {'chan1','chan13','chan14','chan15','chan16','chan43','chan49','chan34'}; % checked; all by ABC. 34 by RMR. Behaves "funny"
info.GP28.badchan     = [info.GP28.badchan, cellstr([repmat('chan',[35 1]) strtrim(num2str((65:99)'))])' cellstr([repmat('chan',[29 1]) strtrim(num2str((100:128)'))])']; % adding the empties, strictly only necessary for session B39
info.GP35.badchan     = {'chan1','chan7','chan8','chan9','chan10','chan16','chan61','chan62','chan64','chan63','chan6','chan13','chan14','chan49'}; % checked; Most from ABC. 7-10,19 because beta buz. Missed by ABC: 63,6. 6 very similar to 7-9. 13/14/49 flat after reref
info.JH2.badchan      = {'chan15','chan47','chan73','chan107','chan108','chan109','chan40','chan43','chan44','chan45','chan48','chan52','chan58','chan86','chan87','chan96','chan104','chan12','chan20'}; % checked; ABC marked 94 as noise, appears fine. ABC also mentioned analogs, but they're separate on disk.; RMR: 102/103 appear silent, but they're not fully dead. After 109, all by RMR, "epilepsyish"
info.JH13.badchan     = {'chan109','chan110'}; % from ABC. Quite a few silent looking channels.
info.JH13.badchan     = [info.JH13.badchan, {'chan83','chan84','chan85','chan86','chan87','chan88','chan91','chan92','chan93','chan94','chan95'}]; % add channels that are very very silent, but not completely empty
info.ST1.badchan      = {'chan5','chan11','chan16','chan41','chan57','chan58'}; % from ABC
info.ST1.badchan      = [info.ST1.badchan, cellstr([repmat('chan',[12 1]) strtrim(num2str((73:84)'))])']; % 73 to 80 mostly noise (from ABC) 81-84 zeros
info.ST6.badchan      = {'chan20','chan31','chan36','chan61','chan62','chan75','chan76','chan33','chan34','chan38','chan39','chan47','chan48','chan65','chan68'}; % from ABC. ABC also threw out 1-19 due to noise, but I don't see anything. They're partially very quit, but nothing else. RMR: 33/34/38/39/48/47/65/68/ epileptic??
info.ST6.badchan      = [info.ST6.badchan, cellstr([repmat('chan',[23 1]) strtrim(num2str((77:99)'))])' cellstr([repmat('chan',[19 1]) strtrim(num2str((100:118)'))])']; % adding 77,78 EKG 79-96 EEG, rest empty 
info.ST8.badchan      = {'chan8','chan9','chan10','chan15','chan16','chan17','chan18','chan32','chan42','chan49','chan50','chan51','chan52','chan61','chan62','chan63','chan64','chan7','chan14','chan48','chan41'}; % from ABC. ABC also marked 8,32,42,49,50,51 as noisy, but can't see anything, they are pretty silent though. RMR last four added because mostly silent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presence of sample-wise artifact specification files
%info.GP14.trlartfctflg     = 0; % still needs commentsand 
info.GP15.trlartfctflg     = 1; % pretty good, some channels might be a bit inter-ictally, but overall looks clean, only a couple of rejected events (conservatively rejected)
info.GP22.trlartfctflg     = 1; % not seeing signs of epileptic activity, might be awesome data; sess1, through out some high var slow freq bits. Some of the lower num channels. sess2 same for very few. sess3 same. sess4 same. sess 5 same. sess6 same. sess7 same
info.GP28.trlartfctflg     = 1; % not seeing much epileptic stuff, could be nice and clean; some outlier events with high amplitudes, some HFO muscles probably in some of the sessions
info.GP35.trlartfctflg     = 1; % might be quite some epileptic stuff, could be baaaaaaaad data
info.JH2.trlartfctflg      = 1; % quite some epileptic stuff, "not cleain'ish", could be bad data
info.JH13.trlartfctflg     = 1; % looks like some epileptic stuff maybe, "cleain'ish";. Also, kinda weird, intermittent very silent channels, but not empty..., many channels but flattish
info.ST1.trlartfctflg      = 1; % looks pretty clean epilepsy wise; yumyum!
info.ST6.trlartfctflg      = 1; % many epileptic waves maybe?, Could be bad data; bleghhh
info.ST8.trlartfctflg      = 1; % might have many epileptic waves; could be bad data
% Summary Good/Bad:
% GOOD: GP22/GP28/GP35/ST1/GP15     BAD: JH2/JH13/ST6/ST8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visual artifact definition in samples
info.GP15.artfctdef.visual.GP15_B12    = []';
info.GP15.artfctdef.visual.GP15_B14    = []';
info.GP15.artfctdef.visual.GP15_B19    = []';
info.GP15.artfctdef.visual.GP15_B20    = [272581;...
                                          274680]';
info.GP15.artfctdef.visual.GP15_B21    = [270506;...
                                          273271]';
info.GP22.artfctdef.visual.GP22_B40    = [ 98135   128708   131883   161759   165112   199006;...
                                          100547   131096   134061   164737   167172   201271]';
info.GP22.artfctdef.visual.GP22_B41    = [119029;...
                                          119949]';
info.GP22.artfctdef.visual.GP22_B43    = [3900    17008    35254    39677    42729   131237   167861;...
                                          6066    18312    39659    41566    44584   133190   169740]';
info.GP22.artfctdef.visual.GP22_B44    = [   1    52558    77423   101193   108400   192277   248583   253317   267239;...
                                          2171    54722    79352   103708   109686   194339   250264   254997   268576]';
info.GP22.artfctdef.visual.GP22_B45    = [34215    48833    89081   110386   138766   171241;...
                                          36624    50444    90823   112841   140392   173157]';
info.GP22.artfctdef.visual.GP22_B46    = [1239    62351    67145    82438;...
                                          2754    64092    68723    83880]';
info.GP22.artfctdef.visual.GP22_B47    = [51263    79675   104733;...
                                          52766    80945   106608]';
info.GP28.artfctdef.visual.GP28_B36    = [161825;...
                                          164325]';
info.GP28.artfctdef.visual.GP28_B38    = [15261    43605    46981    79784   165071   166782;...
                                          16228    44879    47826    81298   166120   167775]';
info.GP28.artfctdef.visual.GP28_B40    = [53767   134289;...
                                          54936   135356]';
info.GP28.artfctdef.visual.GP28_B41    = [10285   119177   122103   150624   164809   174819   177017;...
                                          14258   122080   134073   152600   165685   176698   177934]';
info.GP28.artfctdef.visual.GP28_B42    = [   1    42729    49884    81618   105292   109927   113705   120629   125478;...
                                          2156    44015    50554    82772   109872   111807   115970   122358   128184]';
info.GP28.artfctdef.visual.GP28_B43    = [   1    56714    92866   153939   188142   207537   221100;...
                                          1147    58491    94081   157582   189756   208733   222967]';
info.GP35.artfctdef.visual.GP35_B29    = [29574    70308   132189   187904;...
                                          33193    73248   137359   191397]';
info.GP35.artfctdef.visual.GP35_B30    = [5428    84013   123862;...
                                          8338    85456   125795]';
info.GP35.artfctdef.visual.GP35_B31    = [7517;...
                                          11656]';
info.GP35.artfctdef.visual.GP35_B34    = [4838    33229    47188;...
                                          7571    36624    48832]';
info.GP35.artfctdef.visual.GP35_B35    = [91561   203924;...
                                          96419   207029]';
info.GP35.artfctdef.visual.GP35_B36    = [39829   101423;...
                                          42728   106248]';
info.JH2.artfctdef.visual.Beta_Edfimport_EEG_6507 = [     1  3001    27001    33321    44171    45001    50399    54144    74001    75405    76112    77411   109256,...
                                                     112605   132695   135001   138651   166787   168030   218001   224400   229292   247001   248223   259541   317684,...
                                                     356001   390727   391925   394001;...
                                                        426  3972    27519    33882    44922    45372    50841    54536    74367    76000    76653    78000   109754,...
                                                     112900   133387   135416   139303   167000   168577   218581   224888   230000   247664   248426   260185   318199,...
                                                     356514   391260   392565   394551]';
info.JH13.artfctdef.visual.Edfimport_EEG_10001    = [200001;...
                                                     200001]';
info.ST1.artfctdef.visual.stanford01_B23    = [14511   36625   56636   59979;...
                                               16870   37435   58721   61040]';
info.ST1.artfctdef.visual.stanford01_B24    = [158705;...
                                               160590]';
info.ST1.artfctdef.visual.stanford01_B25    = [1757  128469  140468;...
                                               3592  129721  143365]';
info.ST1.artfctdef.visual.stanford01_B26    = [5502   99569  148319  175368
                                               7641  102474  150601  176838]';
info.ST1.artfctdef.visual.stanford01_B27    = [178   72010  153123;...
                                               1611   73033  154535]';
info.ST1.artfctdef.visual.stanford01_B28    = [39140  105562  147384;...
                                               40470  106888  149249]';
info.ST1.artfctdef.visual.stanford01_B29    = [11294   31806   33623   67145   70169  130823;...
                                               12208   32824   34678   68201   72143  131690]';
info.ST6.artfctdef.visual.ST06_09    = []';
info.ST6.artfctdef.visual.ST06_10    = []';
info.ST6.artfctdef.visual.ST06_11    = []';
info.ST6.artfctdef.visual.ST06_12    = []';
info.ST6.artfctdef.visual.ST06_13    = [115120;...
                                        118050]';
info.ST6.artfctdef.visual.ST06_14    = []';
info.ST6.artfctdef.visual.ST06_15    = [12985;...
                                        16981]';
info.ST8.artfctdef.visual.ANP027    = [19885  105546  134479;...
                                       20694  106043  135136]';
info.ST8.artfctdef.visual.ANP028    = [153332  171833;...
                                       154309  172716]';
info.ST8.artfctdef.visual.ANP029    = [171714;...
                                       172511]';
info.ST8.artfctdef.visual.ANP030    = [37295   63474   80519;...
                                       38501   64328   81123]';
info.ST8.artfctdef.visual.ANP031    = [1;...
                                       885]';
info.ST8.artfctdef.visual.ANP032    = [89674  135303  167579;...
                                       90725  136211  168457]';
info.ST8.artfctdef.visual.ANP033    = [138363;...
                                       139077]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left or right dominance of electrodes (used for selection which trials were (un)attended)
%info.GP14.leftrightelec     =  
info.GP15.leftrightelec     = 'left';
info.GP22.leftrightelec     = 'right';
info.GP28.leftrightelec     = 'right';
info.GP35.leftrightelec     = 'left';
info.JH2.leftrightelec      = 'left';
info.JH13.leftrightelec     = 'left';
info.ST1.leftrightelec      = 'right';
info.ST6.leftrightelec      = 'left';
info.ST8.leftrightelec      = 'left';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low or high pitch dominance in attended/ignored trials, contralateral stimuli, see below
%info.GP14.leftrightelec     =  
info.GP15.lohipitch     = {'high'};         % 
info.GP22.lohipitch     = {'low'};           % 
info.GP28.lohipitch     = {'low','high'};   % 
info.GP35.lohipitch     = {'high'};         % 
info.JH2.lohipitch      = {'low','high'};   % 
info.JH13.lohipitch     = {'low','high'};   % 
info.ST1.lohipitch      = {'low','high'};   % 
info.ST6.lohipitch      = {'low','high'};   % 
info.ST8.lohipitch      = {'low','high'};   % 
%%% CONTRA:
% GP15  nattlo0  nignlo38  natthi37  nignhi35
% GP15  nbattlo36  nbatthi37
% GP22  nattlo68  nignlo61  natthi0  nignhi0
% GP22  nbattlo67  nbatthi35
% GP28  nattlo38  nignlo33  natthi30  nignhi30
% GP28  nbattlo30  nbatthi27
% GP35  nattlo34  nignlo0  natthi33  nignhi30
% GP35  nbattlo29  nbatthi63
% JH2  nattlo35  nignlo30  natthi33  nignhi34
% JH2  nbattlo31  nbatthi35
% ST1  nattlo32  nignlo69  natthi28  nignhi31
% ST1  nbattlo33  nbatthi33
% ST6  nattlo38  nignlo33  natthi51  nignhi35
% ST6  nbattlo24  nbatthi11
% ST8  nattlo37  nignlo34  natthi38  nignhi63
% ST8  nbattlo34  nbatthi31
% JH13  nattlo38  nignlo39  natthi32  nignhi34
% JH13  nbattlo34  nbatthi37
%%% IPSI
% GP15  nattlo39  nignlo37  natthi35  nignhi0
% GP15  nbattlo37  nbatthi38
% GP22  nattlo0  nignlo0  natthi70  nignhi61
% GP22  nbattlo32  nbatthi63
% GP28  nattlo30  nignlo33  natthi29  nignhi34
% GP28  nbattlo31  nbatthi24
% GP35  nattlo36  nignlo33  natthi0  nignhi36
% GP35  nbattlo63  nbatthi34
% JH2  nattlo38  nignlo32  natthi30  nignhi33
% JH2  nbattlo37  nbatthi35
% ST1  nattlo30  nignlo32  natthi71  nignhi31
% ST1  nbattlo34  nbatthi33
% ST6  nattlo39  nignlo58  natthi36  nignhi36
% ST6  nbattlo22  nbatthi5
% ST8  nattlo67  nignlo31  natthi36  nignhi36
% ST8  nbattlo31  nbatthi32
% JH13  nattlo38  nignlo37  natthi38  nignhi37
% JH13  nbattlo37  nbatthi39
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















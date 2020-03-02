
% ========================================================================
% ========================================================================
    % Dissertation Worflow Demonstration
    % Groundwater Recoverability
    % Justin C. Thompson
% ========================================================================
% ========================================================================

% Conversion Constants
    % 1 Horsepower = 745.7 watts = 3,960 gallons/minute/foot
    % 1 Gallon = 231 cubic inches = 0.133681 cubic feet
    % 1 Acre = 6,272,640 square inches = 43,560 square feet
    % 1 Acre-Foot = 43,560 cubic feet = 325,851 gallons

clc
clear
format shortg
    
% ------------------------------------------------------------------------
% MODEL INPUTS
% ------------------------------------------------------------------------

% Agricultural Characteristics ...........................................

% Harvest Value [$/acre/year]
    % Dollar value gained per year per acre of irrigated farmland
        % *** Assumed net of expenditures unrelated to pumping!***
            % JCT ref: where IrArea = 130 and IrDepth = .25
            % JCT ref: 99.4 to economic limit in unconfined
            % JCT ref: 314.5 to economic limit in confined
HarvestValue = 99.4;

% Irrigation Area [acres]
    % ***CAUTION***
        % Model assumes irrigated area is served by a single well!
        % Combines with irrigation depth to derive pumping demand!
IrArea = 130;

% Irrigation Depth [inches/acre/day]
    % ***CAUTION***
        % Model assumes flood irrigation!
        % Combines with irrigation area to derive pumping demand!    
IrDepth = 0.25;

% Growing Season [weeks/year]
    % Texas Almanac = 259 days or 37 weeks
IrWk = 37;
    
% Irrigation Days [days/week]
IrDay = 3;

% Well Characteristics ...................................................

% Power Cost Rate required to run the well [$/kilowatt-hour]
PowerCost = 0.07;

% Radius of the well in the screened interval [inches]
rWell = 3;

% Length of the well screen interval [feet]
ScreenInterval = 100;

% Percentage of well screen open area [percent]
    % Percentage of the screen open area over the total screen area
    % Open area derived from slot size, see manufacturer
ScreenOpen = 10;

% Nonlinear well loss coefficient [seconds^2/foot^5]
    % Walton (1962): 5-10 mild, >10 severe
C = 5;

% Aquifer Characteristics ................................................

% Setting: Confined or Unconfined?
    % Unconfined: Setting = 0
    % Confined: Setting = 1
Setting = 0;

% Initial equilibrium depth-to-water from land surface [feet]
H = 350;

% Depth of aquifer bottom from land surface [feet]
Bottom = 700;

% Depth of aquifer top from land surface[feet]
    % Uncofined: Top = 0
Top = 1650;

% Storage Coefficients
    % Specific Yield [dimensionless]: unconfined storage coefficient
        % Carrizo-Wilcox ranges from 0.05 to 0.3
            % Dutton et.al. (2003)
Sy = 0.15;
    % Storativity [dimensionless]: confined storage coefficient
        % Carrizo-Wilcox ranges from 10^-6 to 10^-1
            % Dutton et.al. (2003)
St = 10^(-3.5);

% Hydraulic Conductivity [feet/day]
    % Carrizo-Wilcox ranges from 50 to 1
        % Dutton et.al. (2003)
K = 7;

% ========================================================================
% MODEL CALCULATIONS
% ========================================================================

% Simple Storage & Recovery Estimations ..................................

% Total Storage [acre-ft]
if Setting == 0
    TS = IrArea*(Bottom-H)*Sy;
elseif Setting == 1
    TS = (IrArea*(Bottom-Top)*Sy)+(IrArea*(Top-H)*St);
    TS_con = (IrArea*(Top-H)*St);
    TS_uncon = (IrArea*(Bottom-Top)*Sy);
end

% Total Estimated Recoverable Storage (TERS) [acre-ft]
TERS_25 = TS*0.25;
TERS_75 = TS*0.75;

fprintf('SIMPLE STORAGE ESTIMATES \n');
fprintf(' Total Storage: %.0f acre-feet. \n', TS);
if Setting == 1
    fprintf('  - Confined: %.0f acre-feet. \n', TS_con);
    fprintf('  - Unconfined: %.0f acre-feet. \n', TS_uncon);
end
fprintf(' TERS @ 25%%: %.0f acre-feet. \n', TERS_25);
fprintf(' TERS @ 75%%: %.0f acre-feet. \n', TERS_75);
fprintf('\n')

% Simulating Potential Changes in Depth-to-Water .........................

% Initial saturated thickness [feet]
if Setting == 0
    SatThk = Bottom-H;
elseif Setting == 1
    SatThk = Bottom-Top;
end

% Array of all possible saturated thickness [feet]
    % Simulates hypothetical dewatering over time in one-foot increments
B = SatThk-1:-1:0;

% Array of all possible transmissivity [square feet/day]
    % Transmissivity = Hydraulic Conductivity x Saturated Thickness
    % Simulates hypothetical dewatering over time in one-foot increments
if Setting == 0
    T = K.*B;
elseif Setting == 1
    T = zeros(1,Bottom-H);
    T(1:Top-H) = K*B(1)+1;
    T(Top-H+1:end) = K.*B;
end

% Array of all possible depth-to-water [feet]
    % Simulates hypothetical drawdown over time in one-foot increments
Z = H+1:1:Bottom;
    
% Calculating Demand Volumes .............................................

% Demand volume from irrigation area and depth [acre-feet/day]
DailyDemand = IrArea*(IrDepth/12);

% Demand volume from irrigation area and depth [acre-feet/year]
AnnualDemand = DailyDemand*IrDay*IrWk;

% Minimum pumping rate [gallons per minute] necessary to meet demand
    % Maximum 24-hour pumping period to meet irrigation demand!
MinGPM = (DailyDemand*325851)/1440;
MinGPMline = zeros(1,length(Z));
MinGPMline(1:end) = MinGPM;
    % Conversion [cubic feet/day]
Qmin = (MinGPM*1440)/7.481;

% Time to depletion of storage [years]
    % *** Assumes constant demand and zero contribution from recharge!
Depletion = TS/AnnualDemand;

fprintf('DEMAND ESTIMATES \n');
fprintf(' Daily Demand Volume: %.2f acre-feet. \n', DailyDemand);
fprintf(' Annual Demand Volume: %.0f acre-feet. \n', AnnualDemand);
fprintf(' Minimum Pumping Rate to Meet Demand: %.0f GPM. \n', MinGPM);
fprintf(' Total Storage / Annual Demand: %.0f years. \n', Depletion);

% Well Screen Entry Speed Check ..........................................

% Maximum Allowable Well Entry Speed [feet/second]
    % ***Commonly assumed 0.1 ft/second!***
    % Efficiecny losses (turbulent) if well entry speed > 0.1 ft/second
MaxEntryVelocity = 0.1;

% Well Entry Speed [feet/second]
MinWellEntryVelocity = MinGPM/(235*rWell*ScreenInterval*(ScreenOpen/100));

% Minimum well screen interval to meet well entry speed maximum at demand
MinWellScreen = (MinGPM/MaxEntryVelocity)/(235*rWell*(ScreenOpen/100));

if MinWellEntryVelocity > MaxEntryVelocity
    fprintf('WARNING! \n')
    fprintf('Well entry speed exceeds 0.1 feet/second! \n')
    fprintf('Minimum well interval @ given demand: %.2f feet. \n', ...
        MinWellScreen);
    fprintf('\n')
else 
    fprintf(' Minimum well screen interval: %.2f feet. \n', ...
        MinWellScreen);
    fprintf('\n')
end

% ------------------------------------------------------------------------
% Demand-Capacity Constraints
% ------------------------------------------------------------------------

% Pumping Time [days]
    % t_pump = the maximum duration of one pumping session
    % Assumes a 24-hr maximum pumping period
t_pump = 1;

% Array of Maximum Possible Drawdown under Pumping (s) [feet]
    % Difference between the depth-to-water and the top of the well screen
        % *** Well/formation losses not explicitly considered
if Setting == 0
    s_max = B-ScreenInterval;
    s_max(s_max < 0) = 0;
elseif Setting == 1
    s_max = fliplr(Z-H-1-ScreenInterval);
    s_max(s_max < 0) = 0;
end

% Array of maximum change in depth to aquifer bottom
    % For theoretical lifting cost calculations
if Setting == 0
    s_max_hypo = fliplr(((Bottom-B)-H)-1);
end
%s_max_hypo = fliplr(Z-H-1-Bottom);

% Specific Capacity Check ................................................

% Array of all possible Specific Capacity [square feet/day]
    % Specific Capacity = (4*pi*T)/[ln((2.25*T*t)/(r^2*S))]
        % Developed from Theis (1935) non-equilibrium by Mace etal (2000)
        % Changes in Specific Capacity from changes in Transmissivity
if Setting == 0
    SC1 = (4*pi*T);
    SC2 = log((2.25*T*t_pump)/((rWell^2)*Sy));
    SC = SC1./SC2;
elseif Setting == 1
    SC1 = (4*pi*T);
    SC2A = log((2.25*T*t_pump)/((rWell^2)*St));
    SC2B = log((2.25*T*t_pump)/((rWell^2)*Sy));
    SC = zeros(1,Bottom-H);
    SC(1:Top-H) = SC1(1:Top-H)./SC2A(1:Top-H);
    SC(Top-H+1:end) = SC1(Top-H+1:end)./SC2B(Top-H+1:end);
end

% Specific Capacity applied to maximum possible drawdown
maxcubicfeetperday = SC.*s_max;

% Maximum pumping rate from Specific Capacity and maximum drawdown
MaxGPM = (maxcubicfeetperday*7.48052)/1440;
MaxGPM(MaxGPM == 0) = NaN;
    % Conversion [cubic feet/day]
Qmax = (MaxGPM*1440)/7.481;

% Maximum depth supported by Maximum Pumping Rate by Specific Capacity
DifGPM = MaxGPM - MinGPM;
DifGPM(DifGPM < 0) = 0;
[y,x] = min(DifGPM);
MaxH_SC = H+x;

% Maximum pumping rate available without turblent flow
    % At commmon limit [0.1 feet/second]:
WellVelocity = MaxGPM./(235*rWell*ScreenInterval*(ScreenOpen/100));
Dif_Velocity1 = WellVelocity - MaxEntryVelocity;
Dif_Velocity1(Dif_Velocity1 < 0) = 0;
[y,x] = min(Dif_Velocity1);
        % Maximum pumping rate [gallons/minute] at velocity limit
CommonVelGPM = MaxGPM(x);
    % At American Water Works Association limit of 1.5 [feet/second]
AWWA = 1.5;
Dif_Velocity2 = WellVelocity - AWWA;
Dif_Velocity2(Dif_Velocity2 < 0) = 0;
[y,x] = min(Dif_Velocity2);
        % Maximum pumping rate [gallons/minute] at velocity limit
AWWAVelGPM = zeros(1,length(Z));
AWWAVelGPM(1:end) = MaxGPM(x);

% Pumping Drawdown Check .................................................

% ************ COOPER-JACOB (1946) simplified approximation **************

% Drawdown observation distance [feet]
    % 1 foot to simulate the outter radius of the gravel pack
r = 1;

% *SOLUTION 1* for drawdown at Minimum Pumping Rate
    % Formation loss (aka "s" or "BQ") solution
if Setting == 0
    CJ1 = log(((r^2)*Sy)./(4*T*t_pump));
    CJ2 = (Qmin./(4*pi*T)).*(-0.5772 - CJ1);
    fLossmin = CJ2;
elseif Setting == 1
    CJ1 = log(((r^2)*St)./(4*T*t_pump));
    CJ2 = (Qmin./(4*pi*T)).*(-0.5772 - CJ1);
    fLossmin = CJ2;
end
    % Well loss (aka "CQ^2") solution
Qmin_sec = Qmin/(1440*60);
wLossmin = C.*Qmin_sec.^2;
    % Total losses (i.e. drawdown, or "BQ+CQ^2") at the well
tLossmin = fLossmin+wLossmin;
    % Difference in Cooper-Jacob (1946) drawdown and available drawdown
sdiff = s_max - tLossmin;
sdiff(sdiff < 0) = 0;
[y,x] = min(sdiff);
maxtru = x - 1;
MaxH_CJ_tru = maxtru+H;

% *SOLUTION 2* for drawdown at Maximum Pumping Rate by Specific Capacity
    % Formation loss (aka "s" or "BQ") solution
if Setting == 0
    CJ1 = log(((r^2)*Sy)./(4*T*t_pump));
    CJ2 = (Qmax./(4*pi*T)).*(-0.5772 - CJ1);
    fLossmax = CJ2;
elseif Setting == 1
    CJ1 = log(((r^2)*St)./(4*T*t_pump));
    CJ2 = (Qmax./(4*pi*T)).*(-0.5772 - CJ1);
    fLossmax = CJ2;
end
    % Well loss (aka "CQ^2") solution
Qmax_sec = Qmax./(1440*60);
wLossmax = C.*Qmax_sec.^2;
    % Total losses (i.e. drawdown) at the well
tLossmax = fLossmax+wLossmax;

% *SOLUTION 3* for Maximum Pumping Rate from drawdown
    % At given range of Qmax from Specific Capacity
    % Formation loss (aka "s" or "BQ") solution
B_CJ = (log(r/(rWell/12)))./(2*pi*T);
BQ = zeros(length(T));
for i = 1:length(T)
    BQ(i,1:end) = B_CJ(i)*Qmax;
end
    % Well loss (aka "CQ^2") solution at given range of Q
CQ = C.*Qmax_sec.^2;
    % Total losses (i.e. drawdown, or "BQ+CQ^2") at the well
tLossCJ = bsxfun(@plus,BQ,CQ);
tLossCJ(end,:) = 0;
    % Identifying difference in drawdown over available thickness
s_diff = zeros(length(T));
s_avail = s_max.';
for i = 1:length(T)
    diffcheck = s_avail - tLossCJ(:,i);
    s_diff(:,i) = diffcheck;
end
s_diff(s_diff<0) = 0;
    % Indexing of first non-zero to locate maximum pumping rate
[r,c] = find(s_diff);
index1  = accumarray(r,c,[size(s_diff,1),1],@min,NaN);
index1 = index1.';
[y,x] = max(index1);
Qmax_CJ = zeros(1,length(T));
for i = 1:x
    ck1 = index1(i);
    ck2 = Qmax(ck1);
    Qmax_CJ(i) = ck2;
end
MaxGPM_CJ = (Qmax_CJ*7.48052)/1440;
MaxGPM_CJ(MaxGPM_CJ == 0) = NaN;
    % Determination of maximum pumping depth
pumpdiff = MaxGPM_CJ - MinGPM;
pumpdiff(pumpdiff < 0) = 0;
[y,x] = min(pumpdiff);
MaxH_CJ = x+H;

fprintf('DEMAND-CAPACITY CONSTRAINTS \n');
fprintf(' * Cooper-Jacob (1946) simplification * \n');
fprintf(' Maximum Depth by Maximum Pumping Rate: %.0f feet. \n', MaxH_CJ);
fprintf(' Maximum Depth by Specific Capacity: %.0f feet. \n', MaxH_SC);
fprintf(' Maximum Depth by Pumping Drawdown: %.0f feet. \n', MaxH_CJ_tru);
fprintf(' *See demand-capacity constraint plot! \n');
fprintf('\n')

% Plotting Demand-Capacity Constraints ...................................

% Formation and well loss in the unconfined setting
if Setting == 0
    figure
    plot(MaxGPM,fLossmax,'b-',MaxGPM,wLossmax,'r-',MaxGPM,tLossmax,'k-')
    line([(CommonVelGPM) (CommonVelGPM)],ylim,'Color','black',...
        'LineStyle',':') 
    set(gca,'XLim',[0 MaxGPM(1)])
    line([MinGPM MinGPM],ylim,'Color','black','LineStyle','-.')
    title('Cooper-Jacob (1946) Formation and Well Loss')
    xlabel('Maximum Pumping Rate [gallons/minute] by Specific Capacity')
    ylabel('Drawdown [feet]')
    legend('Formation Loss: "BQ" or "s"','Well Loss: "CQ^2"',...
        'Formation Loss + Well Loss',...
        'Turbulence Limit @ 0.1 [feet/second]',...
        'Minimum Pumping Rate Demand')
    legend('location','northwest')
    grid on
end

if Setting == 0
    figure
    plot(Z,MaxGPM,Z,MaxGPM_CJ,Z,MinGPMline,'k-')
    set(gca,'XLim',[H Bottom])
    line([(MaxH_SC) (MaxH_SC)],ylim,'Color','black','LineStyle',':');
    line([(MaxH_CJ) (MaxH_CJ)],ylim,'Color','black','LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    title('Demand-Capacity Constraints')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Pumping Rate [gallons/minute]')
    legend('Maximum Pumping Rate by Specific Capacity',...
        'Maximum Pumping Rate by Cooper-Jacob (1946)',...
        'Minimum Pumping Rate Demand',...
        'MaxH from Specific Capacity','MaxH from Cooper-Jacob (1946)',...
        'Top of Well Screen')
    legend('location','northeast')
    grid on
elseif Setting == 1
    figure
    plot(Z,MaxGPM,Z,MaxGPM_CJ,Z,MinGPMline,'k-')
    set(gca,'XLim',[H Bottom])
    line([(Top) (Top)],ylim,'Color','m','LineStyle','--');
    line([(MaxH_SC) (MaxH_SC)],ylim,'Color','black','LineStyle',':');
    line([(MaxH_CJ) (MaxH_CJ)],ylim,'Color','black','LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    title('Demand-Capacity Constraints')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Pumping Rate [gallons/minute]')
    legend('Maximum Pumping Rate by Specific Capacity',...
        'Maximum Pumping Rate by Cooper-Jacob (1946)',...
        'Minimum Pumping Rate Demand','Aquifer Top',...
        'MaxH from Specific Capacity','MaxH from Cooper-Jacob (1946)',...
        'Top of Well Screen')
    legend('location','northeast')
    grid on
end

% ------------------------------------------------------------------------
% ECONOMIC CONSTRAINTS
% ------------------------------------------------------------------------

% Water horsepower at depth-to-water plus drawdown under pumping
HplusS = tLossmin+Z;
WaterHP = (HplusS.*MinGPM)./3960;
WaterHP(WaterHP == 0) = NaN;

% Converted Power Cost Rate [$/watt-minute]
WPwrCst = (PowerCost/60000);

% Pumping Costs ..........................................................

% Pumping cost rate [$/acre-foot]
PmpCost_Rate = (((WaterHP.*745.7).*WPwrCst)./MinGPM).*325851;

% Daily lifting costs [$] at daily demand
PmpCost_Daily = PmpCost_Rate.*DailyDemand;
PmpCost_Daily(end) = PmpCost_Daily(end)*(-1);

% Annual lifting costs [$] 
PmpCost_Annual = PmpCost_Rate.*AnnualDemand;
PmpCost_Annual(end) = PmpCost_Annual(end)*(-1);

% Profit Function ........................................................

% Total annual harvest value [$/year]
TotalAnnualHV = HarvestValue*IrArea;

% Estimated value of water pumped [$/acre-foot]
    % *** Harvest Value assumed net of costs unrelated to pumping!
EstWaterVal = TotalAnnualHV / AnnualDemand;

% Profit function [$/acre-foot]
Profit = EstWaterVal - PmpCost_Rate;
Profit(end) = Profit(end)*(-1);

% Depth-to-water where Profit - Pumping Cost = 0 [feet]
ProfitMod = Profit;
ProfitMod(ProfitMod < 0) = 0;
[y,x] = min(ProfitMod);
MaxH_NoProfit = (x-1)+H;

fprintf('ECONOMIC CONSTRAINTS \n');
fprintf(' Minimum Daily Pumping Cost: $%.2f. \n', min(PmpCost_Daily));
fprintf(' Minimum Annual Pumping Cost: $%.2f. \n', min(PmpCost_Annual));
fprintf(' Daily Pumping Cost @ Demand-Capacity: $%.2f. \n',...
    PmpCost_Daily((MaxH_CJ_tru - H)));
fprintf(' Annual Pumping Cost @ Demand-Capacity: $%.2f. \n',...
    PmpCost_Annual((MaxH_CJ_tru - H)));
fprintf(' Estimated value of water pumped: $%.2f /acre-foot. \n',...
    EstWaterVal);
fprintf(' Depth-to-water where Profit = 0: %.0f feet. \n', MaxH_NoProfit);
fprintf(' *See pumping cost rate plot! \n');
fprintf('\n')

% Plotting Economic Constraints ..........................................

if Setting == 0
    figure
    plot(Z,PmpCost_Rate)
    line([(MaxH_CJ_tru) (MaxH_CJ_tru)],ylim,'Color','black',...
        'LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    set(gca,'YLim',[0 100])
    set(gca,'XLim',[H Bottom])
    title('Pumping Costs with Changes in Transmissivity')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Pumping Cost [$/acre-foot]')
    legend('Pumping Cost Rate [$/acre-foot]','Demand-Capacity Constraint',...
        'Top of Well Screen')
    legend('location','northwest')
    grid on
elseif Setting == 1
    figure
    plot(Z,PmpCost_Rate)
    line([(MaxH_CJ_tru) (MaxH_CJ_tru)],ylim,'Color','black',...
        'LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    line([(Top) (Top)],ylim,'Color','m','LineStyle','--');
    set(gca,'YLim',[0 200])
    set(gca,'XLim',[H Bottom])
    title('Pumping Costs with Changes in Transmissivity')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Pumping Cost [$/acre-foot]')
    legend('Pumping Cost Rate [$/acre-foot]','Demand-Capacity Constraint',...
        'Top of Well Screen','Aquifer Top')
    legend('location','northwest')
    grid on
end

% ------------------------------------------------------------------------
% BINDING CONSTRAINT AND RECOVERABLE YIELD
% ------------------------------------------------------------------------

% Identifying binding constraint as most restrictive depth-to-water
if MaxH_CJ_tru < MaxH_NoProfit
    Constraint = 0;
    MaxH = MaxH_CJ_tru;
elseif MaxH_CJ_tru >= MaxH_NoProfit
    Constraint = 1;
    MaxH = MaxH_NoProfit;
end

% Recoverable storage yield calculation [acre-feet]
if Setting == 0
    RS = IrArea*(MaxH-H)*Sy;
elseif Setting == 1
    if MaxH > Top
        RS = (IrArea*(MaxH-Top)*Sy)+(IrArea*(Top-H)*St);
        RS_con = (IrArea*(Top-H)*St);
        RS_uncon = (IrArea*(MaxH-Top)*Sy);
    elseif MaxH <= Top
        RS = (IrArea*(Top-H)*St);
    end
end
    
fprintf('MAXIMUM ECONOMICALLY RECOVERABLE STORAGE \n');
if Constraint == 0
    fprintf(' Binding constraint: demand-capacity limit \n');
    fprintf(' Maximum depth-to-water: %.0f feet. \n', MaxH_CJ_tru);
elseif Constraint == 1
    fprintf(' Binding constraint: economic limit \n');
    fprintf(' Maximum depth-to-water: %.0f feet. \n', MaxH_NoProfit);
end
if Setting == 0
    fprintf(' Recoverable storage: %.0f acre-feet. \n', RS);
elseif Setting == 1
    fprintf(' Recoverable storage: %.0f acre-feet. \n', RS);
    if MaxH > Top
        fprintf('  - Confined recoverable storage: %.0f acre-feet. \n',...
            RS_con);
        fprintf('  - Unconfined recoverable storage: %.0f acre-feet \n',...
            RS_uncon);
    end
end
fprintf(' Recoverable Storage / Total Storage: %.0f%%. \n', ((RS/TS)*100))
fprintf(' Recoverable Storage / Annual Demand: %.0f years. \n',...
    (RS/AnnualDemand))        
fprintf(' *See maximum recoverable depth-to-water plot! \n')

% Plotting Recoverability ................................................

if Setting == 0
    figure
    plot(Z,Profit)
    line([(MaxH_CJ_tru) (MaxH_CJ_tru)],ylim,'Color','black',...
        'LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    set(gca,'YLim',[0 max(Profit)+(0.1*max(Profit))])
    set(gca,'XLim',[H Bottom])
    title('Maximum Recoverable Depth-to-Water')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Profit [$/acre-foot]')
    legend('Profit Function [$/acre-foot]',...
        'Demand-Capacity Constraint','Top of Well Screen')
elseif Setting == 1
    figure
    plot(Z,Profit)
    line([(MaxH_CJ_tru) (MaxH_CJ_tru)],ylim,'Color','black',...
        'LineStyle','-.');
    line([(Bottom-ScreenInterval) (Bottom-ScreenInterval)],ylim, ...
        'Color','magenta','LineStyle','-');
    line([(Top) (Top)],ylim,'Color','m','LineStyle','--');
    set(gca,'YLim',[0 max(Profit)+(0.1*max(Profit))])
    set(gca,'XLim',[H Bottom])
    title('Maximum Recoverable Depth-to-Water')
    xlabel('Depth-to-Water: Current Equilibrium to Aquifer Bottom [feet]')
    ylabel('Profit [$/acre-foot]')
    legend('Profit Function [$/acre-foot]',...
        'Demand-Capacity Constraint','Top of Well Screen','Aquifer Top')
end

% ========================================================================
    % SIMULATION END
% ========================================================================

function DP()
tic
clc
warning off
DP_input_data
GMIN            = gen_data(:, 2);
GMAX            = gen_data(:, 3);
GINC            = gen_data(:, 4);
GNLC            = gen_data(:, 5);
GSC             = gen_data(:, 6);
GFC             = gen_data(:, 7);
GMINUP          = gen_data(:, 8);
GMINDOWN        = gen_data(:, 9);
GSTATINI        = gen_data(:,10);
GSH             = gen_data(:,11);
GCSTIME         = gen_data(:,12);
GRAMPUP         = gen_data(:,13);
GRAMPDOWN       = gen_data(:,14);
COEF_A          = gen_data(:,15);
COEF_B          = gen_data(:,16);
COEF_C          = gen_data(:,17);
GSDC            = gen_data(:,18);
TAU             = gen_data(:,19);
NG              = size(gen_data,1);
NT              = size(DEMAND,1);
if (DISPATCH_METHOD == 2 | DISPATCH_METHOD == 3) & (any(isnan(GNLC)) | any(isnan(GFC)) | any(isnan(GINC)))
    STR = ['To use linear cost model, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
elseif DISPATCH_METHOD == 1 & (any(isnan(COEF_A)) | any(isnan(COEF_B)) | any(isnan(COEF_C)))
    STR = ['To use quadratic cost model, you must provide data for the cost coefficients:,'...
        'COEFF_A (£), COEFF_B (£/MWh) and COEFF_C (£/MW^2h).']; 
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if MIN_UP_DOWN_TIME_FLAG == 1 & (any(isnan(GMINUP)) | any(isnan(GMINDOWN)))
    STR = ['To use minimum up and down time constraints, you must provide data for GMINUP'...
        ' and GMINDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if RAMP_UP_DOWN_FLAG == 1 & (any(isnan(GRAMPUP)) | any(isnan(GRAMPDOWN)))
    STR = ['To use rump constraints, you must provide data for GRAMPUP '...
        'and GRAMPDOWN.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if COMPLETE_ENUMERATION_FLAG == 0 & (any(isnan(GNLC)) | any(isnan(GFC)) | any(isnan(GINC)))
    STR = ['To use priority list, you must provide data for NO LOAD COSTS,'...
        'FUEL COSTS and INCREMENTAL COSTS.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if START_UP_COST_METHOD == 2 & (any(isnan(GSH)) | any(isnan(GCSTIME)) )
    STR = ['To use cold/hot start up cost method, you must provide data for GSH '...
        'and GCSTIME.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
elseif START_UP_COST_METHOD == 3 & (any(isnan(GSH)) | any(isnan(TAU)) )
    STR = ['To use exponential start up cost method, you must provide data for GSH '...
        'and TAU.'];
    msgbox(STR,'DATA AVAILABILITY CHECK FAILURE!','warn');
    return
end
if MIN_UP_DOWN_TIME_FLAG == 0
    GMINUP(:)   = 1;
    GMINDOWN(:) = 1;
end
if RAMP_UP_DOWN_FLAG == 0              
    GRAMPUP(:)   = Inf;                
    GRAMPDOWN(:) = Inf;                
end
if COMPLETE_ENUMERATION_FLAG == 0        
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG);
else
    [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG);
end
if RESERVE_FLAG == 1
    if exist('RES_UP','var') ~= 1 | isempty(RES_UP) 
        RES_UP = K_RES_UP * DEMAND;                 
    end
    if exist('RES_DN','var') ~= 1 | isempty(RES_DN) 
        RES_DN = K_RES_DN * DEMAND;                 
    end
else
    RES_UP = zeros(size(DEMAND));                   
    RES_DN = zeros(size(DEMAND));                   
end
if START_UP_COST_METHOD == 3                        
    ALPHA = GSH;                                    
    BETA  = GSC;                                    
else
    ALPHA = NaN*ones(NG,1);                         
    BETA  = NaN*ones(NG,1);                         
end
INI_STATE = (GSTATINI > 0);
[I, INI_STATE_NUM]= ismember(INI_STATE',LIST_STATES','rows');
for HOUR = 1:NT
    fprintf('Currently processing hour: %2d \n',HOUR)
    if HOUR == 1
        PREV_STATES_NUM = INI_STATE_NUM;           
        X_PREV  = GSTATINI;                        
        PRODUCTION_PREV = zeros(size(X_PREV));     
        TR_PREV = PREV_STATES_NUM;                 
        FC_PREV = 0;                               
    else
        X_PREV = X_CURR;                           
        PRODUCTION_PREV = PRODUCTION_CURR;
        TR_PREV = TR;                     
        FC_PREV = FC;                     
        PREV_STATES_NUM = TR_PREV(1:COUNTER,end);
    end
    [FEASIBLE_STATES_NUM,SUCCESS] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN);
    if SUCCESS == 0
        return
    end
    
    MN = min(length(PREV_STATES_NUM),N_PRED);
    X_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));
    PRODUCTION_CURR = zeros(NG,MN*length(FEASIBLE_STATES_NUM));
    TR = zeros(MN*length(FEASIBLE_STATES_NUM),HOUR+1);         
    FC = zeros(MN*length(FEASIBLE_STATES_NUM),1);              
    COUNTER = 0;
    for J = 1: length(FEASIBLE_STATES_NUM)
        GEN_START_SHUT_COST = zeros(NG,1);
        TOTAL_COST = zeros(1,length(PREV_STATES_NUM));
        X_TEMP = zeros(NG,length(PREV_STATES_NUM));
        PRODUCTION_TEMP = zeros(NG,length(PREV_STATES_NUM));
        CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM(J));
        for K = 1: length(PREV_STATES_NUM)
            if HOUR == 1;
                PREVIOUS_STATE = INI_STATE;
            else
                PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM(K));
            end
            [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV(:,K),GMINUP,GMINDOWN,NG);
            if SUCCESS==0
                GEN_START_SHUT_COST(:,K) = Inf;
                PROD_COST = ones(NG,1)*Inf;
                GEN_PRODUCTION = ones(NG,1)*NaN;
            else                           
                STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
                if START_UP_COST_METHOD == 1
                    GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* GSC;
                elseif START_UP_COST_METHOD == 2
                    GEN_START_SHUT_COST(:,K) =                            ((STATE_DIFF > 0) & (-X_PREV(:,K) >= (GMINDOWN + GCSTIME))) .* GSC;  % cold start-up cost
                    GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + ((STATE_DIFF > 0) & (-X_PREV(:,K) <  (GMINDOWN + GCSTIME))) .* GSH;  % hot start-up cost
                else
                    GEN_START_SHUT_COST(:,K) = (STATE_DIFF > 0) .* (ALPHA + BETA .* (1-exp(X_PREV(:,K) ./ TAU)));   % exponential start-up costs
                end
                GEN_START_SHUT_COST(:,K) = GEN_START_SHUT_COST(:,K) + (STATE_DIFF  < 0) .* GSDC;   % shut down cost

                [GEN_PRODUCTION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV(:,K),GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD);
            end
            X_TEMP(:,K) = X;
            PRODUCTION_TEMP(:,K) = GEN_PRODUCTION;
            if HOUR == 1
                TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K));
            else
                TOTAL_COST(K) = sum(PROD_COST) + sum(GEN_START_SHUT_COST(:,K)) + FC_PREV(K);
            end
        end
        [MM,II] = sort(TOTAL_COST(TOTAL_COST ~= 0));
        for K = 1:MN
            if MM(K) ~= Inf
                COUNTER = COUNTER +1;
                FC(COUNTER,1) = MM(K);
                TR(COUNTER,1:size(TR_PREV,2)) = TR_PREV(II(K),:);
                TR(COUNTER,end) = FEASIBLE_STATES_NUM(J);
                X_CURR(:,COUNTER) = X_TEMP(:,II(K));
                PRODUCTION_CURR(:,COUNTER) = PRODUCTION_TEMP(:,II(K));
            end
        end
    end
    if COUNTER == 0;                                                        
        STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];
        msgbox(STR,'NO FEASIBLE STATES','warn');                                     
        return
    end
end
[M,I]=min(FC(1:COUNTER));
BEST_PATH = TR(I,:).'; 
evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,DEMAND,GEN_ORDER,GNLC,GFC,GINC,GSC,INI_STATE,NG,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,DETAIL_PRINT_FLAG,GSDC,GSTATINI,...
    GMINUP,GMINDOWN,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU)
warning on
t=toc;
fprintf('\n Elapsed time: %15.4f sec.\n\n',t)
end
function [GMAXcum,GMINcum,LIST_STATES,LIST_INDEX] = priority_list(GNLC,GFC,GMAX,GMIN,GINC,NG)
GFULLAVECOST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX;  
[M,LIST_INDEX] = sort(GFULLAVECOST);        
LIST_STATES = triu(ones(NG));                        
LIST_STATES(LIST_INDEX,:) = LIST_STATES(1:NG,:);     
LIST_STATES = logical(LIST_STATES);               
GMAXcum = cumsum(GMAX(LIST_INDEX));            
GMINcum = cumsum(GMIN(LIST_INDEX));           
prints_states(NG,GMINcum,GMAXcum,LIST_STATES)
end
function [GMAXlst,GMINlst,LIST_STATES,GEN_ORDER] = complete_enumeration(GNLC,GFC,GMAX,GMIN,GINC,NG)
GFULLAVECOST = (0*GNLC + GFC.*GMAX.*GINC/1000)./GMAX; 
[M,GEN_ORDER] = sort(GFULLAVECOST);                   
LIST_STATES=dec2bin(0:2^NG-1)';                       
LIST_STATES = logical(sscanf(LIST_STATES,'%1d',size(LIST_STATES)));
GMINlst = LIST_STATES.' * GMIN; 
GMAXlst = LIST_STATES.' * GMAX; 
[GMAXlst,INDEX]=sort(GMAXlst);  
GMINlst = GMINlst(INDEX);       
LIST_STATES = LIST_STATES(:,INDEX);
prints_states(NG,GMINlst,GMAXlst,LIST_STATES)
end
function [FEASIBLE_STATES_NUM,SUCCESS] = find_feasible_states(GMINlst,GMAXlst,DEMAND,HOUR,RES_UP,RES_DN)
FEASIBLE_STATES_NUM = find((GMINlst <= DEMAND(HOUR)-RES_DN(HOUR)) & (DEMAND(HOUR)+RES_UP(HOUR) <= GMAXlst));
if isempty(FEASIBLE_STATES_NUM)   
    SUCCESS = 0; 
    STR = ['NO FEASIBLE STATES FOR HOUR ',num2str(HOUR),'! PROGRAM TERMINATES!'];
    msgbox(STR,'NO FEASIBLE STATES','warn');
    return
else
    SUCCESS = 1;
end
end
function [GENERATION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD)
if DISPATCH_METHOD == 3
    lb = zeros(size(GMIN));
    ub = zeros(size(GMAX));
    if HOUR ==1
        lb = GMIN .* CURRENT_STATE;
        ub = GMAX .* CURRENT_STATE;
    else
        lb(CURRENT_STATE == 1) = max([GMIN(CURRENT_STATE == 1),PRODUCTION_PREV(CURRENT_STATE == 1)-GRAMPDOWN(CURRENT_STATE == 1)],[],2)...
            .* CURRENT_STATE(CURRENT_STATE == 1);

        ub(PREVIOUS_STATE == 0) = min([GMAX(PREVIOUS_STATE == 0),max([GRAMPUP(PREVIOUS_STATE == 0),GMIN(PREVIOUS_STATE == 0)],[],2)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 0);
        ub(PREVIOUS_STATE == 1) = min([GMAX(PREVIOUS_STATE == 1),PRODUCTION_PREV(PREVIOUS_STATE == 1)+GRAMPUP(PREVIOUS_STATE == 1)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 1);
    end
    if (sum(lb) > DEMAND(HOUR)) | (sum(ub) < DEMAND(HOUR)) | any((ub-lb) < 0)
        GENERATION = ones(NG,1)*NaN;
        PROD_COST = ones(NG,1)*Inf;
    else
        GENERATION = dispatch(CURRENT_STATE,lb,ub,DEMAND,HOUR,GEN_ORDER);
        PROD_COST =  GNLC .* CURRENT_STATE + GFC .* GINC .* GENERATION .* CURRENT_STATE / 1000; % and calculate their costs
    end
    return
else
    Aeq = double(CURRENT_STATE.');
    beq = DEMAND(HOUR);
    lb = zeros(size(GMIN));
    ub = zeros(size(GMAX));
    if HOUR ==1
        lb = GMIN .* CURRENT_STATE;
        ub = GMAX .* CURRENT_STATE;
    else
        lb(CURRENT_STATE == 1) = max([GMIN(CURRENT_STATE == 1),PRODUCTION_PREV(CURRENT_STATE == 1)-GRAMPDOWN(CURRENT_STATE == 1)],[],2)...
            .* CURRENT_STATE(CURRENT_STATE == 1);
        ub(PREVIOUS_STATE == 0) = min([GMAX(PREVIOUS_STATE == 0),max([GRAMPUP(PREVIOUS_STATE == 0),GMIN(PREVIOUS_STATE == 0)],[],2)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 0);
        ub(PREVIOUS_STATE == 1) = min([GMAX(PREVIOUS_STATE == 1),PRODUCTION_PREV(PREVIOUS_STATE == 1)+GRAMPUP(PREVIOUS_STATE == 1)],[],2)...
            .* CURRENT_STATE(PREVIOUS_STATE == 1);
    end
    options = optimset('Display','Off');        
    if DISPATCH_METHOD == 2                     
        f = GFC .* GINC .* CURRENT_STATE / 1000;
        [GENERATION,FVAL,EXITFLAG] = linprog(f,[],[],Aeq,beq,lb,ub,[],options);
        if EXITFLAG > 0
            PROD_COST =  GNLC .* CURRENT_STATE + GFC .* GINC .* GENERATION .* CURRENT_STATE / 1000; % and calculate their costs
        else
            GENERATION = ones(NG,1)*NaN;
            PROD_COST = ones(NG,1)*Inf;
        end
    end
    if DISPATCH_METHOD == 1
        GENERATION = zeros(NG,1);
        X0 = [];
        H = 2*diag(COEF_C(CURRENT_STATE));
        f = COEF_B(CURRENT_STATE);
        Aeq = Aeq(:,CURRENT_STATE);
        lb = lb(CURRENT_STATE);
        ub = ub(CURRENT_STATE);
        [GENERATION1,FVAL,EXITFLAG] = quadprog(H,f,[],[],Aeq,beq,lb,ub,X0,options);      % calculate the optimal production for each generator
        if EXITFLAG > 0
            GENERATION(CURRENT_STATE) = GENERATION1;
            PROD_COST =  (COEF_A.*CURRENT_STATE) + (COEF_B.*GENERATION.*CURRENT_STATE) + (COEF_C.*GENERATION.^2.*CURRENT_STATE); % and calculate their costs
        else
            GENERATION = ones(NG,1)*NaN;
            PROD_COST = ones(NG,1)*Inf;
        end
    end
end
end
function prints_states(NG,GMINcum,GMAXcum,LIST_STATES)
fprintf('   State No.      MW min        MW max                     Units\n')
fprintf('%s',repmat(' ',1,23))
fprintf(['               ',repmat('    %5d ', 1, NG)],1:NG)
fprintf('\n %s \n',repmat('-',1,80'))
for I=1:size(LIST_STATES,2)
    fprintf('      %2d       %8.1f      %8.1f ',I,GMINcum(I),GMAXcum(I))
    fprintf([repmat('       %2d ', 1, size(LIST_STATES,1)) '\n'], LIST_STATES(:,I));
end
end
function evaluate_solution(NT,BEST_PATH,LIST_STATES,GMIN,GMAX,DEMAND,GEN_ORDER,GNLC,GFC,GINC,GSC,INI_STATE,NG,...
    GRAMPUP,GRAMPDOWN,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD,DETAIL_PRINT_FLAG,GSDC,GSTATINI,...
    GMINUP,GMINDOWN,START_UP_COST_METHOD,GCSTIME,GSH,ALPHA,BETA,TAU)
GEN_START_SHUT_COST1 = zeros(NG,NT);
GEN_PRODUCTION1      = zeros(NG,NT);
PROD_COST1           = zeros(NG,NT);
FCOST1               = zeros(NT,1);
GENERATING_COST1     = zeros(NT,1);
GEN_PRODUCTION       = zeros(NG,1);
X  = GSTATINI;
for HOUR = 1:NT
    PREV_STATES_NUM = BEST_PATH(HOUR);
    FEASIBLE_STATES_NUM = BEST_PATH(HOUR+1);
    X_PREV = X;
    if HOUR==1 & PREV_STATES_NUM == 0
        PREVIOUS_STATE = INI_STATE;
    else
        PREVIOUS_STATE = LIST_STATES(:,PREV_STATES_NUM);
    end
    CURRENT_STATE  = LIST_STATES(:,FEASIBLE_STATES_NUM);
    PRODUCTION_PREV = GEN_PRODUCTION;
    [GEN_PRODUCTION,PROD_COST] = production(CURRENT_STATE,PREVIOUS_STATE,GMIN,GMAX,DEMAND,HOUR,GNLC,GFC,GINC,NG,GRAMPUP,GRAMPDOWN,PRODUCTION_PREV,GEN_ORDER,COEF_A,COEF_B,COEF_C,DISPATCH_METHOD);

    STATE_DIFF = CURRENT_STATE - PREVIOUS_STATE;
    [X,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG);
    if START_UP_COST_METHOD == 1 
        GEN_START_SHUT_COST = (STATE_DIFF > 0) .* GSC;
    elseif START_UP_COST_METHOD == 2
        GEN_START_SHUT_COST =                       ((STATE_DIFF > 0) & (-X_PREV >= (GMINDOWN + GCSTIME))) .* GSC;  % hot start-up cost
        GEN_START_SHUT_COST = GEN_START_SHUT_COST + ((STATE_DIFF > 0) & (-X_PREV <  (GMINDOWN + GCSTIME))) .* GSH;  % cold start-up cost
    else
        GEN_START_SHUT_COST = (STATE_DIFF > 0) .* (ALPHA + BETA .* (1-exp(X_PREV ./ TAU)));
    end
    GEN_START_SHUT_COST = GEN_START_SHUT_COST + (STATE_DIFF < 0 ) .* GSDC;
    if HOUR == 1
        TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST);
    else
        TOTAL_COST = sum(PROD_COST) + sum(GEN_START_SHUT_COST) + FCOST1(HOUR-1);
    end 
    FCOST1(HOUR) = TOTAL_COST;
    GENERATING_COST1(HOUR) = sum(PROD_COST);
    GEN_PRODUCTION1(:,HOUR) = GEN_PRODUCTION;
    PROD_COST1(:,HOUR) = PROD_COST;
    GEN_START_SHUT_COST1(:,HOUR) = GEN_START_SHUT_COST;
end
GEN_START_SHUT_COST_TOTAL = sum(GEN_START_SHUT_COST1).';
print_results(BEST_PATH,LIST_STATES,INI_STATE,NT,NG,GMIN,GMAX,DEMAND,FCOST1,GENERATING_COST1,GEN_PRODUCTION1,PROD_COST1,GEN_START_SHUT_COST1,DETAIL_PRINT_FLAG)
end
function print_results(BEST_PATH,LIST_STATES,INI_STATE,NT,NG,GMIN,GMAX,DEMAND,FCOST1,GENERATING_COST1,GENERATION,PROD_COST,START_COST,DETAIL_PRINT_FLAG)
global UnitStates;
if DETAIL_PRINT_FLAG == 0
    S = ['Hour         '
        'Demand       '
        'Tot.Gen      '
        'Min MW       '
        'Max MW       '
        'ST-UP Cost   '
        'Prod.Cost    '
        'F-Cost       '
        'State        '
        'Units ON/OFF '];
    fprintf('\n%s',repmat('=',1,150'))
    fprintf('\n       HOURLY RESULTS:')
    fprintf('\n%s \n',repmat('=',1,150'))
    fprintf([repmat('%12s ', 1, size(S,1))], S');
    fprintf('\n%s\n',repmat('-',1,150'))
else
    S = ['UNITS          '
        'ON/OFF         '
        'GENERATION     '
        'MIN MW         '
        'MAX MW         '
        'ST-UP Cost     '
        'PROD.COST      '];
end
if BEST_PATH(1) == 0
    LIST_STATES = [LIST_STATES,INI_STATE];
    BEST_PATH(1) = size(LIST_STATES,2);
end
for HOUR = 1:length(BEST_PATH)-1
    CURRENT_STATES_NUM  = BEST_PATH(HOUR+1);
    CURRENT_STATE   = LIST_STATES(:,CURRENT_STATES_NUM);
    MIN_MW = CURRENT_STATE.*GMIN;
    MAX_MW = CURRENT_STATE.*GMAX;
    if HOUR ==1 & DETAIL_PRINT_FLAG == 0
        fprintf('%3d  %12s  %12s %12.0f %12.0f %12.0f  %12.0f %12.0f %10.0f ',HOUR-1, '-','-',sum(LIST_STATES(:,BEST_PATH(HOUR)).*GMIN),sum(LIST_STATES(:,BEST_PATH(HOUR)).*GMAX),0,0,0,BEST_PATH(HOUR))
        fprintf(['       ',repmat('%2d', 1, size(LIST_STATES(:,BEST_PATH(HOUR)),1)),'\n'], LIST_STATES(:,BEST_PATH(HOUR)));
    end

    if DETAIL_PRINT_FLAG == 0;
        fprintf('%3d  %12.0f  %12.0f %12.0f %12.0f ',HOUR, DEMAND(HOUR), sum(GENERATION(:,HOUR)), sum(MIN_MW), sum(MAX_MW))
        fprintf('%12.0f  %12.0f %12.0f %10d ',sum(START_COST(:,HOUR)),sum(PROD_COST(:,HOUR)),FCOST1(HOUR),CURRENT_STATES_NUM)
        fprintf(['       ',repmat('%2d', 1, size(CURRENT_STATE,1)),'\n'], CURRENT_STATE);
    else
        TEMP = [(1:NG).',CURRENT_STATE,GENERATION(:,HOUR),MIN_MW,MAX_MW,START_COST(:,HOUR),PROD_COST(:,HOUR)];
        fprintf('\n\n\nHOUR: %2d             DEMAND:%7.1f MW           F-COST: %6.1f £',HOUR,DEMAND(HOUR),FCOST1(HOUR))
        fprintf('\n%s \n',repmat('-',1,120'))
        fprintf([repmat('%15s ', 1, size(S,1)) '\n\n'], S');fprintf('\n');
        fprintf(['%3d %15d ',repmat('%15.1f', 1, size(TEMP,2)-2) '\n'], TEMP.');
        fprintf('%s \n',repmat('-',1,120'))
        fprintf('TOTAL: %12d  %14.1f %14.1f %14.1f %14.1f %14.1f\n',sum(CURRENT_STATE),sum(GENERATION(:,HOUR)), sum(MIN_MW), sum(MAX_MW),sum(START_COST(:,HOUR)),sum(PROD_COST(:,HOUR)))
    end
UnitStates(:,HOUR)=CURRENT_STATE;
end
UnitStates
end
function [X_CURR,SUCCESS] = check_up_down_time(CURRENT_STATE,PREVIOUS_STATE,X_PREV,GMINUP,GMINDOWN,NG)
X_CURR = zeros(NG,1 );
SUCCESS = 1;
if all((X_PREV - GMINUP).*(PREVIOUS_STATE - CURRENT_STATE) >=0 & (-X_PREV - GMINDOWN).*(CURRENT_STATE - PREVIOUS_STATE) >=0)
    for I=1:NG
        if (X_PREV(I) >= 1) & (CURRENT_STATE(I) == 1)
            X_CURR(I) = X_PREV(I) + 1;
        elseif (X_PREV(I) <= -GMINDOWN(I)) & (CURRENT_STATE(I) == 1)
            X_CURR(I) = 1;
        elseif (X_PREV(I) <= -1) & (CURRENT_STATE(I) == 0)
            X_CURR(I) = X_PREV(I) - 1;
        elseif (X_PREV(I) >= GMINUP(I)) & (CURRENT_STATE(I) == 0)
            X_CURR(I) = -1;
        end
    end
else
    SUCCESS = 0;
    X_CURR = ones(NG,1 )*NaN;
    return
end
end
function GENERATION = dispatch(CURRENT_STATE,GMIN,GMAX,DEMAND,HOUR,GEN_ORDER)
GENERATION = GMIN.*CURRENT_STATE;
LOAD = DEMAND(HOUR) - sum(GENERATION);
for K = 1:length(CURRENT_STATE);
    L = GEN_ORDER(K);
    GENERATION(L) = GENERATION(L) + min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L);
    LOAD = LOAD - min(GMAX(L)-GMIN(L),LOAD)*CURRENT_STATE(L);
end
end
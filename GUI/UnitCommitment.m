function varargout = UnitCommitment(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @UnitCommitment_OpeningFcn, ...
    'gui_OutputFcn',  @UnitCommitment_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function UnitCommitment_OpeningFcn(hObject, eventdata, handles, varargin)
global y;
handles.output = hObject;
guidata(hObject, handles);
set(hObject,'toolbar','figure');
axes(handles.maingridpic)
imshow('grid.jpg');
axes(handles.solarpic)
imshow('PV.jpg');
axes(handles.batterypic)
imshow('battery.jpg');
axes(handles.windpic)
imshow('wind.jpg');
axes(handles.forecastaxes);
y = [900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000;900;1060;1200;1080;800;560;580;1000];
plot(y,'-y','LineWidth',2)
whitebg('k')
grid on
title 'Hourly Load Curve of Grid Station'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Peak Demand (MW)'
legend('(MW)','Location','best');
function varargout = UnitCommitment_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
function mainchk_Callback(hObject, eventdata, handles)
extraconnections(handles)
function solarchk_Callback(hObject, eventdata, handles)
extraconnections(handles)
function batterychk_Callback(hObject, eventdata, handles)
extraconnections(handles)
function windchk_Callback(hObject, eventdata, handles)
extraconnections(handles)
function prioritychk_Callback(hObject, eventdata, handles)
% hObject    handle to prioritychk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global UnitStates y1;
UnitStates=zeros(8,24);
gen_data = [...
    1        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN
    2        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN
    3        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN
    4        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN
    5        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN
    6        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN
    7        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN
    8        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN
    ];
if (get(hObject,'Value')==1)
    set(handles.dpchk,'Value',0);
    PL
    for i=1:24
        unitdata(:,1)=UnitStates(:,i);
        hourval=0;
        for j=1:(length(gen_data(:,1)))
            if unitdata(j,1)==1
                unitdata(j,2)=gen_data(j,3);
                hourval=hourval+gen_data(j,3);
            else
                unitdata(j,2)=0;
            end
        end
        set(handles.hourtxt,'String',i);
        set(handles.demandtxt,'String',y1(i));
        set(handles.generationtxt,'String',hourval);
        
        t=uitable(handles.uc,'Data',unitdata,'Position',[475 50 185 205]);
        t.ColumnName = {'Status','Output'};
        pause(2)
    end
end
UnitStates
function dpchk_Callback(hObject, eventdata, handles)
% hObject    handle to dpchk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global UnitStates y1;
UnitStates=zeros(8,24);
gen_data = [...
    1        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN
    2        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN
    3        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN
    4        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN
    5        25     80      10440           213.00       350            2.00          4             2           -5          150                4                50             75          NaN         NaN           NaN               0               NaN
    6        60    250       9000           585.62       400            2.00          5             3           +8          170                5                80            120          NaN         NaN           NaN               0               NaN
    7        75    300       8730           684.74      1100            2.00          5             4           +8          500                5               100            150          NaN         NaN           NaN               0               NaN
    8        20     60      11900           252.00      0.02            2.00          1             1           -6            0                0                80            120          NaN         NaN           NaN               0               NaN
    ];
if (get(hObject,'Value')==1)
    set(handles.prioritychk,'Value',0);
    DP
    for i=1:24
        unitdata(:,1)=UnitStates(:,i);
        hourval=0;
        for j=1:(length(gen_data(:,1)))
            if unitdata(j,1)==1
                unitdata(j,2)=gen_data(j,3);
                hourval=hourval+gen_data(j,3);
            else
                unitdata(j,2)=0;
            end
        end
        set(handles.hourtxt,'String',i);
        set(handles.demandtxt,'String',y1(i));
        set(handles.generationtxt,'String',hourval);
        t=uitable(handles.uc,'Data',unitdata,'Position',[475 50 185 205]);
        t.ColumnName = {'Status','Output'};
        pause(2)
    end
end
UnitStates
function allchk_Callback(hObject, eventdata, handles)
% hObject    handle to allchk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(hObject,'Value')==1)
    set(handles.mainchk,'Value',1);
    set(handles.solarchk,'Value',1);
    set(handles.batterychk,'Value',1);
    set(handles.windchk,'Value',1);
else
    set(handles.mainchk,'Value',0);
    set(handles.solarchk,'Value',0);
    set(handles.batterychk,'Value',0);
    set(handles.windchk,'Value',0);
end
extraconnections(handles)
function extraconnections(handles)
global y y1 ysolar ybattery ywind;
axes(handles.solaraxes);
if get(handles.solarchk,'Value')==1
    ysolar = 200*[0.0008
        0.0004
        0.0053
        0.0035
        0.0032
        0.0345
        0.0909
        0.1314
        0.2034
        0.3605
        0.6993
        0.8269
        0.8168
        0.8069
        0.7280
        0.4920
        0.2317
        0.1155
        0.0841
        0.0672
        0.0062
        0.0024
        0.0081
        0.0097
        ];
else
    ysolar = zeros(24,1);
end
plot(ysolar,'-y','LineWidth',1)
whitebg('k')
grid on
axis([0 24 0 200]);
title 'Power Curve of Solar'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');



axes(handles.batteryaxes);
if get(handles.batterychk,'Value')==1
    ybattery = 50*[0.6692
        0.1904
        0.3689
        0.4607
        0.9816
        0.1564
        0.8555
        0.6448
        0.3763
        0.1909
        0.4283
        0.4820
        0.1206
        0.5895
        0.2262
        0.3846
        0.5830
        0.2518
        0.2904
        0.6171
        0.2653
        0.8244
        0.9827
        0.7302
        ];
else
    ybattery = zeros(24,1);
end
plot(ybattery,'-y','LineWidth',1)
whitebg('k')
grid on
title 'Battery Storage'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');
axis([0 24 0 50]);




axes(handles.windaxes);
if get(handles.windchk,'Value')==1
    ywind = 100*[0.3395
        0.9516
        0.9203
        0.0527
        0.7379
        0.2691
        0.4228
        0.5479
        0.9427
        0.4177
        0.9831
        0.3015
        0.7011
        0.6663
        0.5391
        0.6981
        0.6665
        0.1781
        0.1280
        0.9991
        0.1711
        0.0326
        0.5612
        0.8819
        ];
else
    ywind = zeros(24,1);
end
plot(ywind,'-y','LineWidth',1)
whitebg('k')
grid on
axis([0 24 0 100]);
title 'Power Curve of Wind'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');



y1 = y-ysolar-ybattery-ywind;
axes(handles.maingridaxes);
plot(y1,'-y','LineWidth',1)
whitebg('k')
grid on
axis([0 24 0 1500]);
title 'Load Curve for Generators'
xlabel 'Time (Hour of Day)'
ylabel 'Hourly Power (MW)'
legend('(MW)','Location','best');



function hourtxt_Callback(hObject, eventdata, handles)
function hourtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function demandtxt_Callback(hObject, eventdata, handles)
function demandtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function generationtxt_Callback(hObject, eventdata, handles)
function generationtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

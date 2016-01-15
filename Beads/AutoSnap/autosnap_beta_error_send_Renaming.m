function [ output_args ] = autosnap_beta1(pads, positions, padnames, colors, outputfolder,fname_title,auto_print)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% inut parameters:
% pads: number of pads
% positions: numbero of positions per pads
% padnames: name of pads
% colors: in the form 'ryp'
% outputfolder: name of output folder

handles.outputfolder= outputfolder;
mkdir(handles.outputfolder);
handles.current_directory=cd;
for i=1:pads
    cd([handles.current_directory]);
    cd([handles.current_directory, '\',handles.outputfolder]);
    fname=[padnames{i} '_' num2str(str3(i))];
    fname_array{i}=[padnames{i} '_' num2str(str3(i))];
    redseg_better_del_no_cells_auto_no_schmutz(fname, 'naturalback',1);
    num_col=length(colors);
    analysesnaphist_auto(fname,colors,fname_title{i});
end
%Generating print file
figure;
j=1;
k=1;
for i=1:pads
    if j>9
        j=1;
        saveas(gcf,['output_',num2str(k),'.jpg']);
        saveas(gcf,['output_',num2str(k),'.fig']);
        figure;
        k=k+1;
    end
    subplot_tight(3,3,j);
    if exist(['S',fname_array{i},'-1.jpg'])~=0
        image(imread(['S',fname_array{i},'-1.jpg']));
        axis off;
        j=j+1;
    end
end
saveas(gcf,['output_',num2str(k),'.jpg']);
saveas(gcf,['output_',num2str(k),'.fig']);

if auto_print==1
    h=hgload('output.fig');
    set(gcf,'Units','centimeters');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperType','A4');
    set(gcf, 'PaperPosition', [0.2 0.5 29 19.5 ]);
    print(h, '-PSouth Wing Colour Laserjet');
end


function send_error
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Defining Account details
myaddress = 'nikon.bacillus.1@gmail.com';
mypassword = 'P455w0rd!';

%Setting Prefences
setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

%Setting Properties
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

%sending Email
sendmail(myaddress, 'Acquisition Error', 'Come and check me. I am having a problems!');
end
end



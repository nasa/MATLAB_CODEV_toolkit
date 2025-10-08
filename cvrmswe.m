function [rmsval,xshift,yshift,strehl,output] = cvrmswe(nrd,f,w,z,mode)
%CVRMSWE  gets RMS wavefront information from the RMSWE macro
%
%   function [rmsval,xshift,yshift,strehl,all] = cvrmswe(nrd,f,w,z,focus_mode)
%
%   INPUTS: nrd = NxN grid of rays across pupil
%           f = field number
%           w = wavelength number
%           z = zoom position
%           focus_mode = nominal (1=default) or best focus (0)
%
%   OUTPUT: Indivual outputs are vectors of length(max_fields+1), where the
%               last number is the composite value.  If only one field
%               exists, then no composite value is returned.
%           The "all" output is an array as follows:
%               -For nominal focus: 5x(N+1), N=# of fields
%               Columns = RMS WFE,Focus,X-Shift,Y-Shift,Strehl
%               -For best focus:   10x(N+1), N=# of fields
%               RMS WFE (Comp),Focus (Comp),X-Shift (Comp),Y-Shift (Comp),
%               Strehl (Comp),RMS WFE (Individual),Focus shift (Individual),
%               X-Shift (Individual),Y-Shift (Individual),Strehl (Individual)
%
%   See also: CVMACRO

[nums,numf,numw,numz] = cvnum;
if nargin<1, nrd = 20; end
if nargin<2, f = 1:numf; end
if nargin<3, w = 1; end
if nargin<4, z = 1; end
if nargin<5, mode = 1; end

[nums,numf,numw,numz] = cvnum;
rms_vals=zeros(10,numf+1);

cvmacro('cvrmswe.seq',[nrd, w, z, mode]);

field_lab = [num2cell(1:numf),{'Comp'}];
out2(2:numf+2,1)=field_lab';
if mode==1,
    output=cvbuf(1,1:numf+1,1:5);
    labels = {'Field','RMS WFE','Focus','X-Shift','Y-Shift','Strehl'};
    out2(1,1:6)=labels;
    out2(2:numf+2,2:6)=num2cell(output);
else
    output=cvbuf(1,1:numf+1,1:10);
    labels = {'Field','RMS WFE (Comp)','Focus (Comp)','X-Shift (Comp)','Y-Shift (Comp)','Strehl (Comp)',...
          'RMS WFE (Ind)','Focus (Ind)','X-Shift (Ind)','Y-Shift (Ind)','Strehl (Ind)'};
      out2(1,1:11)=labels;
    out2(2:numf+2,2:11)=num2cell(output);
end

rmsval = output(f,1);
xshift = output(f,3);   
yshift = output(f,4);
strehl = output(f,5);

if nargout<1
    disp(out2);
    if numf>1
        figure,
        plot([1,numf],[output(numf+1,1),output(numf+1,1)],'k');
        hold on
        plot(1:numf,output(1:numf,1),'.-'); hold off
        title('RMS WFE for all fields')
        xlabel('Field Number')
        ylabel('RMS WFE (waves)')
        legend('Composite RMS')
    end
end
    
% Copyright © 2004-2005 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration.  No copyright is claimed in the United States 
% under Title 17, U.S. Code. All Other Rights Reserved.
% 
% Authors: Joseph M. Howard, Blair L. Unger, Mark E. Wilson, NASA
% Revision Date: 2007.08.22     
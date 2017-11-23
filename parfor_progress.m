function percent = parfor_progress(N)
%   Progress monitor (progress bar) that works with parfor.
%   parfor_progress works by creating a file called parfor_progress.txt in
%   matlab folder, and then keeping track of the parfor loop's
%   progress within that file.
%
%   parfor_progress(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   parfor_progress updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   parfor_progress(0) deletes parfor_progress.txt and finalizes progress
%   bar.

narginchk(0, 1);

if nargin < 1
    N = -1;
end

percent = 0;

bar_width = 50;  

if N > 0
    
    f = fopen('parfor_progress.txt', 'w');
    fprintf(f, '%d\n', N);                 % Save N at the top of progress.txt
    fclose(f);
    
%     if nargout == 0
%         disp(['  0%[>', repmat(' ', 1, bar_width), ']']);
%     end
elseif N == 0
    delete('parfor_progress.txt');
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (bar_width+9)), char(10), '100%[', repmat('=', 1, bar_width+1), ']']);
    end
else
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run parfor_progress(N) before parfor_progress to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); 
        disp([repmat(char(8), 1, (bar_width+9)), char(10), perc, '[', repmat('=', 1, round(percent*bar_width/100)), '>', repmat(' ', 1, bar_width - round(percent*bar_width/100)), ']']);
    end
end

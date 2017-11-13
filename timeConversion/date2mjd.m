function [mjd]=date2mjd(date)
    
    mjd2000 = date2mjd2000(date);
    [mjd]=mjd20002mjd(mjd2000);
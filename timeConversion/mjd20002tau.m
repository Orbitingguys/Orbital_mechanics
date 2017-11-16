function [tau] = mjd20002tau(mjd2000)
%
%       - tau:      represents te value of the time at the Modified Julian Date since 2000 Jan
%                   12:00:00 (mjd2000) following the criteria from Matlab to measure time.

    Date = mjd20002date(mjd2000);
    tau = datenum(datetime(Date));

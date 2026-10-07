%% sat2ijk.m
% Converts satellite data into IJK (geocentric equatorial / ECI) position
% and velocity vectors.
%
% Input options (set inputMode below):
%   'TLE' - two-line element set (e.g., from CelesTrak or Space-Track)
%   'COE' - classical orbital elements (a, e, i, RAAN, argp, nu)
%
% Method: COEs -> perifocal (PQW) frame -> rotate to IJK
%   [IJK <- PQW] = R3(-RAAN) * R1(-i) * R3(-argp)
%
% Units: km, km/s, degrees (inputs)
% Requires MATLAB R2016b or newer (local functions in scripts).

clear; clc; close all;

%% Constants
mu = 398600.4418;   % Earth gravitational parameter [km^3/s^2]
Re = 6378.137;      % Earth equatorial radius [km]

%% ---------------- USER INPUT ----------------
inputMode = 'TLE';  % 'TLE' or 'COE'

switch upper(inputMode)
    case 'TLE'
        % Example TLE (ISS format) - replace with a current TLE
        line1 = '1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9005';
        line2 = '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391428322';
        coe = tle2coe(line1, line2, mu);

    case 'COE'
        coe.a    = 7000;    % semi-major axis [km]
        coe.e    = 0.01;    % eccentricity [-]
        coe.i    = 45;      % inclination [deg]
        coe.RAAN = 30;      % right ascension of ascending node [deg]
        coe.argp = 60;      % argument of perigee [deg]
        coe.nu   = 90;      % true anomaly [deg]
        % Special-case angles (only used if orbit is circular/equatorial):
        % coe.u      = ...; % argument of latitude [deg] (circular inclined)
        % coe.lambda = ...; % true longitude [deg]      (circular equatorial)
        % coe.wtrue  = ...; % longitude of periapsis [deg] (elliptical equatorial)

    otherwise
        error('inputMode must be ''TLE'' or ''COE''.');
end
%% --------------------------------------------

%% Convert to IJK
[r_ijk, v_ijk] = coe2ijk(coe, mu);

%% Display results
fprintf('===== Orbital Elements =====\n');
fprintf('a    = %12.4f km\n',  coe.a);
fprintf('e    = %12.7f\n',     coe.e);
fprintf('i    = %12.4f deg\n', coe.i);
fprintf('RAAN = %12.4f deg\n', coe.RAAN);
fprintf('argp = %12.4f deg\n', coe.argp);
fprintf('nu   = %12.4f deg\n', coe.nu);
if isfield(coe, 'epoch')
    fprintf('Epoch: %s UTC\n', datestr(coe.epoch, 'yyyy-mm-dd HH:MM:SS.FFF'));
end

fprintf('\n===== IJK State Vector =====\n');
fprintf('r = %12.4f I + %12.4f J + %12.4f K   km   (|r| = %.4f km)\n', r_ijk, norm(r_ijk));
fprintf('v = %12.6f I + %12.6f J + %12.6f K   km/s (|v| = %.6f km/s)\n', v_ijk, norm(v_ijk));
fprintf('Altitude = %.4f km\n', norm(r_ijk) - Re);

%% Plot the orbit in IJK
nuSweep = linspace(0, 360, 361);
orbit = zeros(3, numel(nuSweep));
tmp = coe;
for k = 1:numel(nuSweep)
    tmp.nu = nuSweep(k);
    tmp.u = coe.argp + nuSweep(k);  % keeps special cases consistent
    tmp.lambda = coe.RAAN + coe.argp + nuSweep(k);
    orbit(:,k) = coe2ijk(tmp, mu);
end

figure; hold on; grid on; axis equal;
[xs, ys, zs] = sphere(50);
surf(Re*xs, Re*ys, Re*zs, 'FaceColor', [0.3 0.5 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
plot3(orbit(1,:), orbit(2,:), orbit(3,:), 'k-', 'LineWidth', 1.5);
plot3(r_ijk(1), r_ijk(2), r_ijk(3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
quiver3(0,0,0, 1.5*Re,0,0, 'r', 'LineWidth', 1.5);  text(1.6*Re,0,0,'I');
quiver3(0,0,0, 0,1.5*Re,0, 'g', 'LineWidth', 1.5);  text(0,1.6*Re,0,'J');
quiver3(0,0,0, 0,0,1.5*Re, 'b', 'LineWidth', 1.5);  text(0,0,1.6*Re,'K');
xlabel('I [km]'); ylabel('J [km]'); zlabel('K [km]');
title('Satellite Orbit in IJK Frame'); view(3);

%% ======================== LOCAL FUNCTIONS ========================

function [r_ijk, v_ijk] = coe2ijk(coe, mu)
% COE2IJK  Classical orbital elements -> IJK position/velocity.
    tol  = 1e-8;
    a    = coe.a;
    e    = coe.e;
    i    = deg2rad(coe.i);
    RAAN = deg2rad(coe.RAAN);
    argp = deg2rad(coe.argp);
    nu   = deg2rad(coe.nu);

    % Special cases (Vallado convention)
    if e < tol && i < tol                       % circular equatorial
        argp = 0; RAAN = 0;
        if isfield(coe, 'lambda'), nu = deg2rad(coe.lambda); end
    elseif e < tol                              % circular inclined
        argp = 0;
        if isfield(coe, 'u'), nu = deg2rad(coe.u); end
    elseif i < tol                              % elliptical equatorial
        RAAN = 0;
        if isfield(coe, 'wtrue'), argp = deg2rad(coe.wtrue); end
    end

    p = a*(1 - e^2);  % semi-latus rectum [km]

    % Perifocal (PQW) frame
    r_pqw = [p*cos(nu)/(1 + e*cos(nu));
             p*sin(nu)/(1 + e*cos(nu));
             0];
    v_pqw = [-sqrt(mu/p)*sin(nu);
              sqrt(mu/p)*(e + cos(nu));
              0];

    % Rotation matrix PQW -> IJK = R3(-RAAN) * R1(-i) * R3(-argp)
    cO = cos(RAAN); sO = sin(RAAN);
    ci = cos(i);    si = sin(i);
    cw = cos(argp); sw = sin(argp);

    Q = [ cO*cw - sO*sw*ci,  -cO*sw - sO*cw*ci,   sO*si;
          sO*cw + cO*sw*ci,  -sO*sw + cO*cw*ci,  -cO*si;
          sw*si,              cw*si,              ci   ];

    r_ijk = Q*r_pqw;
    v_ijk = Q*v_pqw;
end

function coe = tle2coe(line1, line2, mu)
% TLE2COE  Parse a two-line element set into classical orbital elements.
    line1 = char(line1); line2 = char(line2);

    % Line 1: epoch
    coe.satnum = strtrim(line1(3:7));
    yy  = str2double(line1(19:20));
    doy = str2double(line1(21:32));
    if yy < 57, yr = 2000 + yy; else, yr = 1900 + yy; end
    coe.epoch = datetime(yr,1,1,0,0,0) + days(doy - 1);

    % Line 2: elements
    coe.i    = str2double(line2(9:16));                    % deg
    coe.RAAN = str2double(line2(18:25));                   % deg
    coe.e    = str2double(['0.' strtrim(line2(27:33))]);   % implied decimal
    coe.argp = str2double(line2(35:42));                   % deg
    M        = deg2rad(str2double(line2(44:51)));          % mean anomaly [rad]
    n_rev    = str2double(line2(53:63));                   % rev/day

    % Semi-major axis from mean motion
    n = n_rev*2*pi/86400;          % rad/s
    coe.a = (mu/n^2)^(1/3);        % km

    % Mean anomaly -> eccentric anomaly -> true anomaly
    E  = solveKepler(M, coe.e);
    nu = 2*atan2(sqrt(1 + coe.e)*sin(E/2), sqrt(1 - coe.e)*cos(E/2));
    coe.nu = mod(rad2deg(nu), 360);
end

function E = solveKepler(M, e)
% SOLVEKEPLER  Newton-Raphson solution of M = E - e*sin(E).
    if e < 0.8, E = M; else, E = pi; end
    for k = 1:100
        dE = (M - E + e*sin(E)) / (1 - e*cos(E));
        E = E + dE;
        if abs(dE) < 1e-12, break; end
    end
end

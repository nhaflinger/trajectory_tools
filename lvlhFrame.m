function varargout = lvlhFrame(arg1, arg2, varargin)
%LVLHFRAME  Transform position/velocity between ECI and LVLH (Hill / RSW) frame.
%
%   LVLH frame definition (RSW):
%     x (R): radial outward (along position vector)
%     y (S): along-track (perpendicular to R in orbital plane, in direction of motion)
%     z (W): normal to orbital plane (= h_hat = R x V direction)
%
%   Usage:
%     [r_lvlh, v_lvlh] = lvlhFrame(r_eci, v_eci)
%       ECI -> LVLH transformation using r_eci as the reference.
%
%     DCM = lvlhFrame(r_eci, v_eci, 'dcm')
%       Return the 3x3 rotation matrix (DCM) from ECI to LVLH only.
%
%     [r_eci_out, v_eci_out] = lvlhFrame(r_lvlh, v_lvlh, r_ref_eci, v_ref_eci, 'inverse')
%       LVLH -> ECI using reference orbit r_ref_eci, v_ref_eci.
%
%     [dr_lvlh, dv_lvlh] = lvlhFrame(r_chaser_eci, v_chaser_eci, r_target_eci, v_target_eci, 'relative')
%       Relative position/velocity of chaser w.r.t. target, in target LVLH frame.
%
%   Inputs:
%     r_eci, v_eci       - ECI position (km) and velocity (km/s), 3x1
%     r_ref_eci, v_ref_eci - reference orbit for inverse/relative transforms, 3x1
%
%   Outputs:
%     r_lvlh, v_lvlh     - position (km) and velocity (km/s) in LVLH frame, 3x1
%     DCM                - 3x3 rotation matrix from ECI to LVLH

r1 = arg1(:);
r2 = arg2(:);

%% ── Parse mode from varargin ─────────────────────────────────────────────────
mode = 'forward';   % default: ECI -> LVLH
if ~isempty(varargin)
    last_arg = varargin{end};
    if ischar(last_arg) || isstring(last_arg)
        mode = lower(strtrim(last_arg));
    end
end

%% ── Build DCM for reference orbit ────────────────────────────────────────────
switch mode
    case {'forward', 'dcm'}
        r_ref = r1;
        v_ref = r2;

    case {'inverse', 'relative'}
        if numel(varargin) < 2
            error('lvlhFrame: ''%s'' mode requires r_ref_eci and v_ref_eci as 3rd and 4th arguments', mode);
        end
        r_ref = varargin{1}(:);
        v_ref = varargin{2}(:);

    otherwise
        error('lvlhFrame: unknown mode ''%s''. Valid: (none), ''dcm'', ''inverse'', ''relative''', mode);
end

% Compute LVLH basis vectors from reference orbit
R_hat = r_ref / norm(r_ref);
h_vec = cross(r_ref, v_ref);
W_hat = h_vec / norm(h_vec);
S_hat = cross(W_hat, R_hat);

% DCM: rows are the LVLH basis vectors expressed in ECI
DCM_eci2lvlh = [R_hat'; S_hat'; W_hat'];   % 3x3

%% ── Apply requested transformation ──────────────────────────────────────────
switch mode
    case 'forward'
        varargout{1} = DCM_eci2lvlh * r1;
        varargout{2} = DCM_eci2lvlh * r2;

    case 'dcm'
        varargout{1} = DCM_eci2lvlh;

    case 'inverse'
        % LVLH -> ECI: DCM transpose (orthogonal matrix)
        DCM_lvlh2eci = DCM_eci2lvlh';
        varargout{1} = DCM_lvlh2eci * r1;
        varargout{2} = DCM_lvlh2eci * r2;

    case 'relative'
        % Chaser relative to target in target LVLH frame
        % r1/r2 are chaser ECI position/velocity
        % r_ref/v_ref are target ECI position/velocity
        dr_eci = r1 - r_ref;
        dv_eci = r2 - v_ref;
        varargout{1} = DCM_eci2lvlh * dr_eci;
        varargout{2} = DCM_eci2lvlh * dv_eci;
end
end

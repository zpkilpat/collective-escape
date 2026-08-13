function st = paper_style()
% PAPER_STYLE  Single source of truth for figure styling.
st.LW   = 5;    st.LWt = 3;    st.FSZ = 32;   st.FIG = [22 20];
st.col_H1   = [0.80 0.15 0.10];   % red    threat present / detection
st.col_H0   = [0.18 0.45 0.70];   % blue   no threat
st.col_naiv = [0.45 0.45 0.45];   % grey   naive rule
st.col_bay  = [0.13 0.55 0.13];   % green  Bayesian / heuristic
st.col_ref  = [0.55 0.55 0.55];   % grey   reference lines
st.green_ramp = @(n) [0.70 0.86 0.74].*(1-linspace(0,1,n)') + ...
                     [0.04 0.32 0.17].*linspace(0,1,n)';
st.blue_ramp  = @(n) [0.62 0.80 0.98].*(1-linspace(0,1,n)') + ...
                     [0.04 0.16 0.42].*linspace(0,1,n)';
end

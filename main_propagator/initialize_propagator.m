%% before_loop 
load('Sat_params.mat');
Aux_fileread;
AuxParam.n      = 10;
AuxParam.m      = 10;
AuxParam.n_a    = 10;
AuxParam.m_a    = 10;
AuxParam.n_G    = 10;
AuxParam.m_G    = 10;

radius_limit = 1;
options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
